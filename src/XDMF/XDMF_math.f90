#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

! This module provides general routines used by the XDMF system
module XDMF_math

    use Tools_log
    use SFC_data_types
    use XDMF_data_types

    use, intrinsic :: iso_fortran_env

    implicit none

    contains

    ! This function calculates a cells offset in the refinement quadtree using its position
    ! Essentially, this is a form of convenient hash on the cell position and type
    subroutine xdmf_hash_element(element, grid_scale, bq_hdf5_offs)
        type(t_element_base), intent(in)					            :: element
        integer(INT64), intent(in)	                                    :: grid_scale
        integer(INT64), intent(out)                                     :: bq_hdf5_offs
        
        integer(INT64), dimension(2, 3)			                        :: pos_grid
        integer(INT64), dimension(2)                                    :: pos_grid_tresh, tresh_cmp, bq_pos
        integer(INT64)                                                  :: type_flag, bq_size_half, bq_size_cur_offs, q_path_i, i, q_depth
    
        do i = 1, 3
            pos_grid(:, i) = grid_scale * element%nodes(i)%ptr%position
        end do
    
        q_depth = ishft(element%cell%geometry%i_depth, -1)
        bq_size_half = ishft(grid_scale, -(q_depth + 1))
        type_flag = abs(element%cell%geometry%i_plotter_type)
    
        ! Compute the position of this cell normalized to quad bounds
        bq_pos = (/ minval(pos_grid(1, :)), minval(pos_grid(2, :)) /)
        if(type_flag.eq.1) then
            bq_pos(2) = bq_pos(2) - bq_size_half
        else if (type_flag.eq.7) then
            bq_pos(1) = bq_pos(1) - bq_size_half
        end if
    
        ! Do a binary search in two dimensions to find the fitting quad for this cell
        bq_hdf5_offs = 0
        pos_grid_tresh = (/ ishft(grid_scale, -1), ishft(grid_scale, -1) /)
        do i = 1, q_depth
            bq_size_cur_offs = ishft(grid_scale, -(i + 1))
            tresh_cmp = pos_grid_tresh - bq_pos
            if(tresh_cmp(1).le.0) then
                pos_grid_tresh(1) = pos_grid_tresh(1) + bq_size_cur_offs
                if(tresh_cmp(2).le.0) then
                    pos_grid_tresh(2) = pos_grid_tresh(2) + bq_size_cur_offs
                    q_path_i = 1
                else
                    pos_grid_tresh(2) = pos_grid_tresh(2) - bq_size_cur_offs
                    q_path_i = 4
                end if
            else
                pos_grid_tresh(1) = pos_grid_tresh(1) - bq_size_cur_offs
                if(tresh_cmp(2).le.0) then
                    pos_grid_tresh(2) = pos_grid_tresh(2) + bq_size_cur_offs
                    q_path_i = 2
                else
                    pos_grid_tresh(2) = pos_grid_tresh(2) - bq_size_cur_offs
                    q_path_i = 3
                end if
            end if
            bq_hdf5_offs = ishft(bq_hdf5_offs, 2) + q_path_i
        end do

        ! Mix the cell type into the hash
        bq_hdf5_offs = ishft(bq_hdf5_offs, 3) + (type_flag - 1) + 1
    end subroutine

    ! This subroutine checks whether a given odd integer is prime
    subroutine is_prime(n, isp)
        integer(INT64), intent(in)      :: n
        logical, intent(out)            :: isp

        integer(INT64)                  :: i
        integer(INT64)    	            :: end

        ! Only test integers up to sqrt(n)
        end = int(sqrt(real(n)), INT64) + 1
        isp = .false.

        ! Only test odd integers
        do i = 3, end, 2
            if (mod(n, i).eq.0) then
                return
            end if
        end do

        isp = .true.
    end subroutine

    ! This subroutine return the closest prime in either direction to an even integer
    subroutine close_prime(n, bigger, i)
        integer(INT64), intent(in)      :: n
        logical, intent(in)             :: bigger
        integer(INT64), intent(out)     :: i

        logical                         :: isp

        isp = .false.
        if (bigger) then
            ! Make sure new integer is bigger and odd
            i = n + 1
            ! Make sure it is prime
            call is_prime(i, isp)
            do while (.not.isp)
                i = i + 2
                call is_prime(i, isp)
            end do
        else
            ! See above
            i = n - 1
            call is_prime(i, isp)
            do while (.not.isp)
                i = i - 2
                call is_prime(i, isp)
            end do
        end if
    end subroutine

    ! This subroutine calculates the prime size of a hashtable 
    ! and the prime used for the secondary hash function
    subroutine xdmf_hashtable_params(num_cells, table_size, hash2_prime)
        integer(INT64), intent(in)      :: num_cells
        integer(INT64), intent(out)     :: table_size, hash2_prime

        ! Compute table size, set a load factor of <= 0.5
        select case (num_cells)
            case (0)
                table_size = 0
                hash2_prime = 1
            case (1)
                table_size = 3
                hash2_prime = 2
            case (2)
                table_size = 5
                hash2_prime = 3
            case default
                call close_prime(int(real(num_cells, REAL64) * 2.0_REAL64, INT64), .true., table_size)
                ! Compute secondary hash function base
                call close_prime(table_size - 1, .false., hash2_prime)
        end select
    end subroutine

    ! This subroutine hashes an element using double hashing
    subroutine xdmf_hashtable_hash(element, table_size, hash2_prime, try, key)
        integer(INT64), intent(in)      :: element, table_size, hash2_prime, try
        integer(INT64), intent(out)     :: key

        key = mod(mod(element, table_size) + (try * (hash2_prime - mod(element, hash2_prime))), table_size)
    end subroutine

end module