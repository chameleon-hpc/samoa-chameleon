#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

! This module provides the XDMF core API for the output traversal
module XDMF_output_base
        
    use SFC_data_types
    use Tools_openmp
    use Tools_mpi

    use XDMF_data_types
    use XDMF_output_base_data_types
    use XDMF_math
    use XDMF_hdf5
    use XDMF_config

    use, intrinsic :: iso_fortran_env
    use, intrinsic :: iso_c_binding

    implicit none

    type(t_xdmf_section_buffer), save, target                           :: sect_store_data_cells, sect_store_data_patches

    contains

#   if defined (_XDMF_PATCH)    
        subroutine xdmf_base_pre_traversal_grid_op(base, sections_ptr, grid, param_cells, param_patches)
#   else
        subroutine xdmf_base_pre_traversal_grid_op(base, sections_ptr, grid, param_cells)
#   endif
        type(t_xdmf_base_output_traversal), intent(inout)				:: base
        type(t_xdmf_base_output_traversal_ptr), dimension(:)            :: sections_ptr 
        type(t_grid), intent(inout)							            :: grid
        type(t_xdmf_parameter), intent(in)                              :: param_cells
#       if defined (_XDMF_PATCH)
            type(t_xdmf_parameter), intent(in)                          :: param_patches
#       endif

        integer                                                         :: file_basename_index, i, j, sect_j
        logical                                                         :: write_cp, no_filter
        type(t_xdmf_layout_descriptor)                                  :: layout_desc_bku
        integer(XDMF_GRID_SI)                                           :: num_cells_complete
        real                                                            :: filter_amount

        ! Compute several meta infos about the domain
        base%grid_scale = ishft(1_XDMF_GRID_SI, ishft(cfg%i_max_depth + 1, -1))
        file_basename_index = index(base%s_file_stamp, "/", .true.)
        base%s_file_stamp_base = trim(base%s_file_stamp(file_basename_index+1:))
        if ((base%xdmf_remove_lines .eq. 0) .and. cfg%xdmf%l_xdmfcheckpoint) then
            base%xdmf_remove_lines = 1 + (cfg%xdmf%i_xdmfoutput_iteration - cfg%xdmf%i_xdmfcheckpoint_iteration)
        else
            base%xdmf_remove_lines = 1
        end if

        ! Compute whether to output tree (checkpoint) data and to read from filter
        if(cfg%xdmf%i_xdmfcpint.eq.0) then
            write_cp = .false.
        else
            write_cp = mod(int(base%output_iteration, GRID_SI), int(cfg%xdmf%i_xdmfcpint, GRID_SI)).eq.0
        end if
        no_filter = cfg%xdmf%i_xdmffilter_index .eq. 0 .or. write_cp

        ! Compute number of cells
        call xdmf_compute_num_cells(size(sections_ptr), grid, num_cells_complete, .false.)
        if (no_filter) then
            base%num_cells = num_cells_complete
            filter_amount = 100
        else
            call xdmf_compute_num_cells(size(sections_ptr), grid, base%num_cells, .true.)
            filter_amount = (real(base%num_cells) / real(num_cells_complete)) * 100
        end if
        
        ! Calculate section offsets in file
        call xdmf_allocate_compute_sections_descs(base, size(sections_ptr), grid, .not. no_filter)

        if (rank_MPI .eq. 0) then
            if (no_filter) then
                _log_write(1, '(A, I0, A, I0)') " XDMF: Output step: ", base%output_iteration, &
                    ", simulation step: ", base%i_sim_iteration
            else
                _log_write(1, '(A, I0, A, I0, A, I0, A, I0, A, F0.2, A)') " XDMF: Output step: ", base%output_iteration, &
                    ", simulation step: ", base%i_sim_iteration, ", filter: ", base%num_cells, &
                    " of ", num_cells_complete, " cells (", filter_amount, "%)"
            end if
        end if

        ! Allocate buffers for cell data for the all sections of this rank
        if ((iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_cells) .ne. 0) .or. write_cp) then
            call sect_store_data_cells%allocate(base%root_layout_desc%ranks(rank_MPI + 1)%num_cells, param_cells)
        end if
#       if defined (_XDMF_PATCH)
            if ((iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_patches) .ne. 0) .or. write_cp) then
                call sect_store_data_patches%allocate(base%root_layout_desc%ranks(rank_MPI + 1)%num_cells, param_patches)
            end if
#       endif
        
        ! Scatter computated data across all sections
        do i = 1, size(sections_ptr)
            call base%root_layout_desc%scatter_to(sections_ptr(i)%ptr%root_layout_desc)
            sections_ptr(i)%ptr%sect_store_cells%ptr => sect_store_data_cells
            sections_ptr(i)%ptr%sect_store_patches%ptr => sect_store_data_patches
            sections_ptr(i)%ptr%s_file_stamp_base = base%s_file_stamp_base
            sections_ptr(i)%ptr%s_file_stamp = base%s_file_stamp
            sections_ptr(i)%ptr%grid_scale = base%grid_scale
            sections_ptr(i)%ptr%output_iteration = base%output_iteration
            sections_ptr(i)%ptr%num_cells = base%num_cells
            sections_ptr(i)%ptr%xdmf_remove_lines = base%xdmf_remove_lines
        end do
    end subroutine

#   if defined (_XDMF_PATCH)
        subroutine xdmf_base_post_traversal_grid_op(base, sections_ptr, grid, param_cells, param_patches)
#   else
        subroutine xdmf_base_post_traversal_grid_op(base, sections_ptr, grid, param_cells)
#   endif
        type(t_xdmf_base_output_traversal), intent(inout)				:: base
        type(t_xdmf_base_output_traversal_ptr), dimension(:)            :: sections_ptr   
        type(t_grid), intent(inout)							            :: grid
        type(t_xdmf_parameter), intent(in)                              :: param_cells
#       if defined (_XDMF_PATCH)   
            type(t_xdmf_parameter), intent(in)                          :: param_patches
#       endif

        integer                                                         :: error, result, num_cells, r
        integer (HSIZE_T)                                               :: hdf5_vals_length_cells, hdf5_tree_length_cells, hdf5_vals_offset_cells
#       if defined(_XDMF_PATCH)
            integer (HSIZE_T)                                           :: hdf5_vals_length_patches, hdf5_vals_offset_patches
#       endif
        integer (HSIZE_T)                                               :: i, j

        integer (XDMF_GRID_SI)                                          :: output_meta_iteration
        logical                                                         :: hdf5_step_exists, new_file, write_cp
        character(len = 256)					                        :: file_name_h5

        integer (INT64)                                                 :: htbl_size, htbl_hash2, htbl_try, htbl_key, htbl_num_cells
        integer (INT32)                                                 :: htbl_element, hdf5_tree_index
        integer (HSIZE_T)                                               :: htbl_size_local, htbl_size_local_spec

#       if defined(_MPI)
            integer (MPI_ADDRESS_KIND)                                  :: htbl_lookup_offset
            integer                                                     :: htbl_mpi_win, htbl_lookup_rank
            integer (INT32), dimension(1)                               :: htbl_switch_origin, htbl_switch_comp, htbl_switch_result, htbl_source
            type (c_ptr)                                                :: htbl_mpi_win_baseptr
            integer(INT32), dimension(:, :), pointer                    :: htbl_alloc
#       else
            integer(INT32), dimension(:, :), allocatable                :: htbl_alloc
#       endif

        ! Every n steps, a new HDF5 file will be created
        output_meta_iteration = int(real(base%output_iteration, XDMF_GRID_SR) / cfg%xdmf%i_xdmfspf, XDMF_GRID_SI)
        new_file = mod(int(base%output_iteration, XDMF_GRID_SI), int(cfg%xdmf%i_xdmfspf, XDMF_GRID_SI)).eq.0

        ! Compute whether to output tree (checkpoint) data
        if(cfg%xdmf%i_xdmfcpint.eq.0) then
            write_cp = .false.
        else
            write_cp = mod(int(base%output_iteration, XDMF_GRID_SI), int(cfg%xdmf%i_xdmfcpint, XDMF_GRID_SI)).eq.0
        end if

        ! Open or create file and init meta hdf5 ids
        write (file_name_h5, "(A, A, I0, A, A)") trim(base%s_file_stamp), "_", output_meta_iteration, "_xdmf.h5", char(0)
        call base%root_desc%hdf5_meta_ids%createopen(file_name_h5, base%output_iteration, .not.new_file, .true., .true., .true., hdf5_step_exists, result)
        if (result.ne.0) then
            _log_write(1, '(A, I0, A)') " XDMF: Error: Cannot open or create file: ", trim(base%s_file_stamp)
            assert(.false.)
            call exit()
        end if
        if (hdf5_step_exists .and. (rank_MPI .eq. 0)) then
            _log_write(1, '(A, I0)') " XDMF: WARNING: Overwriting output step: ", base%output_iteration
        end if

        ! Get neccessary meta infos for this rank (how much data to write where)
        hdf5_vals_length_cells = base%root_layout_desc%ranks(rank_MPI + 1)%num_cells
        hdf5_vals_offset_cells = base%root_layout_desc%ranks(rank_MPI + 1)%offset_cells
        hdf5_tree_length_cells = hdf5_vals_length_cells
#       if defined (_XDMF_PATCH)  
            hdf5_tree_length_cells = hdf5_tree_length_cells / _XDMF_PATCH_ORDER_SQUARE
            hdf5_vals_length_patches = base%root_layout_desc%ranks(rank_MPI + 1)%num_patches
            hdf5_vals_offset_patches = base%root_layout_desc%ranks(rank_MPI + 1)%offset_patches
#       endif

        ! Compute tree hashtable parameters
        htbl_size = 0
        if (write_cp) then
            call xdmf_hashtable_params(int(base%num_cells, HSIZE_T), htbl_size, htbl_hash2)
        end if

        ! Create hdf5 datasets
#       if defined (_XDMF_PATCH)
            num_cells = base%num_cells * _XDMF_PATCH_ORDER_SQUARE
#       else
            num_cells = base%num_cells
#       endif
        if ((iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_cells) .ne. 0) .or. write_cp) then
            call base%root_desc%hdf5_ids_cells%create(param_cells, base%root_desc%hdf5_meta_ids%step_group_cells_id, num_cells, htbl_size)
        end if
#       if defined (_XDMF_PATCH)
            if ((iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_patches) .ne. 0) .or. write_cp) then
                call base%root_desc%hdf5_ids_patches%create(param_patches, base%root_desc%hdf5_meta_ids%step_group_patches_id, base%num_cells, 0_XDMF_GRID_DI)
            end if
#       endif

        ! Check for empty domain, otherwise MPI-IO might fail
        if (base%num_cells .gt. 0) then
            if (write_cp) then
                if (rank_MPI.eq.0) then
                    _log_write(1, '(A, I0)') " XDMF: Output unfiltered checkpoint: ", base%output_iteration
                end if
                ! Compute and allocate tree hash table allocation map distribution
#               if defined(_MPI)
                    htbl_size_local = max(htbl_size / size_MPI, 1)
                    ! Last rank takes the remainder
                    htbl_size_local_spec = htbl_size_local
                    if ((rank_MPI .eq. (size_MPI - 1)) .and. ((htbl_size_local * size_MPI) .le. htbl_size)) then
                        htbl_size_local_spec = htbl_size_local + (htbl_size - (htbl_size_local * size_MPI))
                    end if
                    ! Allocate MPI shared memeory, field size = int32 = 4 bytes * num fields
                    call mpi_win_allocate(htbl_size_local_spec * 4 * hdf5_tree_width, 4, MPI_INFO_NULL, MPI_COMM_WORLD, &
                        htbl_mpi_win_baseptr, htbl_mpi_win, error); assert_eq(error, 0)
                    ! Make memory accessible to fortran
                    call c_f_pointer(htbl_mpi_win_baseptr, htbl_alloc, (/ hdf5_tree_width, htbl_size_local_spec /))
                    ! Zero memory
                    htbl_alloc = 0
                    ! Enable shared access
#               else
                    htbl_size_local = htbl_size
                    htbl_size_local_spec = htbl_size
                    allocate(htbl_alloc(hdf5_tree_width, htbl_size), stat = error); assert_eq(error, 0)
                    htbl_alloc = 0
#               endif

                ! Generate, compress and store tree information into distributed hash table
                if(hdf5_vals_length_cells.gt.0) then
                    ! Build the offsets of tree data in hash table in memory
                    htbl_num_cells = 0
                    do i = 1, hdf5_tree_length_cells
#                       if defined (_XDMF_PATCH)                    
                            hdf5_tree_index = 1 + hdf5_vals_offset_cells + ((i - 1) * _XDMF_PATCH_ORDER_SQUARE)
#                       else
                            hdf5_tree_index = 1 + hdf5_vals_offset_cells + (i - 1)
#                       endif
                        ! Hash the tree offset to compute hash table offset
                        htbl_element = sect_store_data_cells%tree(i)
                        htbl_try = 0
                        do
                            ! Calculate hash table position
                            call xdmf_hashtable_hash(int(htbl_element, INT64), htbl_size, htbl_hash2, htbl_try, htbl_key)
                            ! Check if rehash is needed
#                           if defined(_MPI)
                                ! Compute target rank and offset for lookup in allocation table
                                if (htbl_key .lt. (htbl_size_local * size_MPI)) then
                                    htbl_lookup_rank = (htbl_key / htbl_size_local)
                                    htbl_lookup_offset = mod(htbl_key, htbl_size_local)
                                else
                                    htbl_lookup_rank = size_MPI - 1
                                    htbl_lookup_offset = htbl_key - (htbl_size_local * (size_MPI - 1))
                                end if
                                ! Atomically read and switch in shared memory
                                htbl_switch_origin(1) = htbl_element + 1
                                htbl_switch_comp(1) = 0
                                call mpi_win_lock(MPI_LOCK_EXCLUSIVE, htbl_lookup_rank, 0, htbl_mpi_win, error); assert_eq(error, 0)
                                call mpi_compare_and_swap(htbl_switch_origin, htbl_switch_comp, htbl_switch_result, MPI_INTEGER4, &
                                    htbl_lookup_rank, htbl_lookup_offset * hdf5_tree_width, htbl_mpi_win, error); assert_eq(error, 0)
                                call mpi_win_unlock(htbl_lookup_rank, htbl_mpi_win, error); assert_eq(error, 0)
                                ! Either it did switch, and result (prev val) = 0, 
                                ! or cell (bucket) is allocated. If with this quads hash, thats fine too, add data
                                if ((htbl_switch_result(1) .eq. htbl_switch_comp(1)) .or. &
                                    (htbl_switch_result(1) .eq. htbl_switch_origin(1))) then
                                        htbl_source(1) = hdf5_tree_index
                                        call mpi_win_lock(MPI_LOCK_EXCLUSIVE, htbl_lookup_rank, 0, htbl_mpi_win, error); assert_eq(error, 0)
                                        call mpi_put(htbl_source, 1, MPI_INTEGER4, htbl_lookup_rank, &
                                            (htbl_lookup_offset * hdf5_tree_width) + 1, &
                                            1, MPI_INTEGER4, htbl_mpi_win, error); assert_eq(error, 0)
                                        call mpi_win_unlock(htbl_lookup_rank, htbl_mpi_win, error); assert_eq(error, 0)
                                        exit
                                else
                                    ! Hash collision, rehash
                                    htbl_try = htbl_try + 1
                                end if
#                           else
                                ! See above
                                if ((htbl_alloc(1, htbl_key + 1) .eq. 0) .or. &
                                    (htbl_alloc(1, htbl_key + 1) .eq. (htbl_element + 1))) then
                                    htbl_alloc(:, htbl_key + 1) = (/ htbl_element + 1, hdf5_tree_index /)
                                    exit
                                else
                                    htbl_try = htbl_try + 1
                                end if
#                           endif
                        end do
                    end do
                end if

                ! Write the hashtable chunk into the HDF5 dataset
                if ((rank_MPI * htbl_size_local) .lt. htbl_size) then
                    call hdf5_write_chunk_int(base%root_desc%hdf5_ids_cells%tree%dset_id, &
                        base%root_desc%hdf5_ids_cells%tree%dspace_id, &
                        (/ 0_HSIZE_T, (rank_MPI * htbl_size_local) /), &
                        (/ hdf5_tree_width, htbl_size_local_spec /), &
                        hdf5_rank, hdf5_subset_stride, hdf5_subset_block, &
                        htbl_alloc, base%root_desc%hdf5_meta_ids%access_dset_id)
                else
                    call hdf5_write_chunk_int(base%root_desc%hdf5_ids_cells%tree%dset_id, &
                        base%root_desc%hdf5_ids_cells%tree%dspace_id, &
                        (/ 0_HSIZE_T, 0_HSIZE_T /), (/ 0_HSIZE_T, 0_HSIZE_T /), &
                        hdf5_rank, hdf5_subset_stride, hdf5_subset_block, &
                        htbl_alloc, base%root_desc%hdf5_meta_ids%access_dset_id)
                end if

                ! Free memory for tree hash table allocation table
#               if defined(_MPI)
                    call mpi_win_free(htbl_mpi_win, error); assert_eq(error, 0)
#               else
                    deallocate(htbl_alloc, stat = error); assert_eq(error, 0)
#               endif
            end if

            ! Write the sections cell attribute buffers to the HDF5 file
            if ((iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_cells) .ne. 0) .or. write_cp) then
                call xdmf_write_data_buffers(base, base%root_desc%hdf5_ids_cells, sect_store_data_cells, &
                    hdf5_vals_offset_cells, hdf5_vals_length_cells, param_cells)
            end if
#           if defined (_XDMF_PATCH)
                if ((iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_patches) .ne. 0) .or. write_cp) then
                    call xdmf_write_data_buffers(base, base%root_desc%hdf5_ids_patches, sect_store_data_patches, &
                        hdf5_vals_offset_patches, hdf5_vals_length_patches, param_patches)
                end if
#           endif
        else
            if (rank_MPI.eq.0) then
                _log_write(1, '(A, I0)') " XDMF: WARNING: No output data at step: ", base%output_iteration
            end if
        end if

        ! Deallocate per section arrays
        do i = 1, size(sections_ptr)
            call sections_ptr(i)%ptr%root_layout_desc%deallocate()
        end do
        call base%root_layout_desc%deallocate()

        ! Free section buffers
        if ((iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_cells) .ne. 0) .or. write_cp) then
            call sect_store_data_cells%deallocate()
        end if
#       if defined (_XDMF_PATCH)
            if ((iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_patches) .ne. 0) .or. write_cp) then
                call sect_store_data_patches%deallocate()
            end if
#       endif

        ! Close hdf5 ids
        if ((iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_cells) .ne. 0) .or. write_cp) then
            call base%root_desc%hdf5_ids_cells%close(param_cells)
        end if
#       if defined (_XDMF_PATCH)
            if ((iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_patches) .ne. 0) .or. write_cp) then
                call base%root_desc%hdf5_ids_patches%close(param_patches)
            end if
#       endif
        call base%root_desc%hdf5_meta_ids%close()
    end subroutine

    ! This routine writes out the actual buffers into the hdf file
    subroutine xdmf_write_data_buffers(base, hdf5_ids, sect_store_data, hdf5_vals_offset, hdf5_vals_length, param)
        type(t_xdmf_base_output_traversal), intent(in)				    :: base
        type(t_xdmf_hdf5), intent(in)				                    :: hdf5_ids
        type(t_xdmf_section_buffer)                                     :: sect_store_data
        integer (HSIZE_T)                                               :: hdf5_vals_length, hdf5_vals_offset
        type(t_xdmf_parameter), intent(in)                              :: param

        integer (HSIZE_T)                                               :: i, j
        integer                                                         :: error
        integer (INT32), dimension(:, :), allocatable                   :: hdf5_sect_topo_buffer

        do i = 1, param%hdf5_valsi_width
            call hdf5_write_chunk_int(hdf5_ids%batch_valsi(i)%dset_id, &
                hdf5_ids%batch_valsi(i)%dspace_id, &
                (/ 0_HSIZE_T, hdf5_vals_offset /), (/ 1_HSIZE_T, hdf5_vals_length /), &
                hdf5_rank, hdf5_subset_stride, hdf5_subset_block, &
                sect_store_data%valsi(:, i), base%root_desc%hdf5_meta_ids%access_dset_id)
        end do
        do i = 1, param%hdf5_valsr_width
            call hdf5_write_chunk_real(hdf5_ids%batch_valsr(i)%dset_id, &
                hdf5_ids%batch_valsr(i)%dspace_id, &
                (/ 0_HSIZE_T, hdf5_vals_offset /), (/ param%hdf5_attr_width, hdf5_vals_length /), &
                hdf5_rank, hdf5_subset_stride, hdf5_subset_block, &
                sect_store_data%valsr(:, :, i), base%root_desc%hdf5_meta_ids%access_dset_id)
        end do
        do i = 1, param%hdf5_valsuv_width
            call hdf5_write_chunk_real(hdf5_ids%batch_valsuv(i)%dset_id, &
                hdf5_ids%batch_valsuv(i)%dspace_id, &
                (/ 0_HSIZE_T, 0_HSIZE_T, hdf5_vals_offset /), (/ hdf5_vector_width, param%hdf5_attr_width, &
                hdf5_vals_length /), hdf5_rank_v, hdf5_subset_stride_v, hdf5_subset_block_v, &
                sect_store_data%valsuv(:, :, :, i), base%root_desc%hdf5_meta_ids%access_dset_id)
        end do
        call hdf5_write_chunk_real(hdf5_ids%valsg%dset_id, &
            hdf5_ids%valsg%dspace_id, &
            (/ 0_HSIZE_T, hdf5_vals_offset * param%hdf5_valst_width /), &
            (/ param%hdf5_valsg_width, hdf5_vals_length * param%hdf5_valst_width /), &
            hdf5_rank, hdf5_subset_stride, hdf5_subset_block, &
            sect_store_data%valsg, base%root_desc%hdf5_meta_ids%access_dset_id)

        ! Generate and write topology indices for this section
        allocate(hdf5_sect_topo_buffer(param%hdf5_valst_width, hdf5_vals_length), stat = error); assert_eq(error, 0)
        do i = 0, hdf5_vals_length - 1
            do j = 0, param%hdf5_valst_width - 1
                hdf5_sect_topo_buffer(j + 1 , i + 1) = ((hdf5_vals_offset + i) * param%hdf5_valst_width) + j
            end do
        end do
        call hdf5_write_chunk_int(hdf5_ids%valst%dset_id, &
            hdf5_ids%valst%dspace_id, &
            (/ 0_HSIZE_T, hdf5_vals_offset /), (/ param%hdf5_valst_width, hdf5_vals_length /), &
            hdf5_rank, hdf5_subset_stride, hdf5_subset_block, &
            hdf5_sect_topo_buffer, base%root_desc%hdf5_meta_ids%access_dset_id)
        deallocate(hdf5_sect_topo_buffer, stat = error); assert_eq(error, 0)
    end subroutine

    ! This routine computes how many cells the complete domain contains
    subroutine xdmf_compute_num_cells(sections_size, grid, num_cells, from_filter)
        integer(XDMF_GRID_SI), intent(in)                               :: sections_size 
        type(t_grid), intent(inout)							            :: grid
        integer(XDMF_GRID_SI), intent(out)							    :: num_cells
        logical, intent(in)                                             :: from_filter

        integer                                                         :: i, error
        type(t_section_info)                                            :: section_info

#       if defined(_MPI)
            integer(XDMF_GRID_SI), dimension(:), allocatable            :: section_sizes
#       endif

        num_cells = 0
        ! Sum over the number of cells in each section
#       if defined(_MPI)
            allocate(section_sizes(sections_size), stat = error); assert_eq(error, 0)
            do i = 1, sections_size
                if (.not. from_filter) then
                    section_info = grid%sections%elements(i)%get_info()
                    section_sizes(i) = section_info%i_cells
                else
                    section_sizes(i) = grid%sections%elements(i)%xdmf_filter_count_cells
                end if
            end do
            ! Using MPI, gather the number of cells per rank and sum it
            call reduce(num_cells, section_sizes, MPI_SUM, .true.)
            deallocate(section_sizes, stat = error); assert_eq(error, 0)
#       else
            do i = 1, sections_size
                if (.not. from_filter) then
                    section_info = grid%sections%elements(i)%get_info()
                    num_cells = num_cells + section_info%i_cells
                else
                    num_cells = num_cells + grid%sections%elements(i)%xdmf_filter_count_cells
                end if
            end do
#       endif
    end subroutine

    ! This routine computes the length and offset in the HDF5 file for each section across all ranks
    subroutine xdmf_allocate_compute_sections_descs(base, sections_size, grid, from_filter)
        type(t_xdmf_base_output_traversal), intent(inout)		        :: base
        integer(XDMF_GRID_SI), intent(in)                               :: sections_size  
        type(t_grid), intent(inout)							            :: grid
        logical, intent(in)                                             :: from_filter

        integer                                                         :: i, j, error
        type(t_section_info)                                            :: section_info
        integer(XDMF_GRID_SI)                                           :: offset_cells_current, offset_cells_current_buffer
        integer(XDMF_GRID_SI)                                           :: num_cells_section, num_cells_section_buffer
#       if defined(_XDMF_PATCH)
            integer(XDMF_GRID_SI)                                       :: offset_patches_current, offset_patches_current_buffer
            integer(XDMF_GRID_SI)                                       :: num_patches_section, num_patches_section_buffer
#       endif
#       if defined(_MPI)
            integer(INT32), dimension(:), allocatable                   :: num_sections, num_cells_local, sect_index_local, sect_index_local_inv
            integer(INT32), dimension(:, :), allocatable                :: num_cells, sect_index
            integer(INT32), dimension(1)                                :: num_sections_local
            integer                                                     :: max_num_sections
#           if defined(_XDMF_PATCH)
                integer(INT32), dimension(:), allocatable               :: num_patches_local
                integer(INT32), dimension(:, :), allocatable            :: num_patches
#           endif
#       else
            integer                                                     :: section_index
#       endif
        
#       if defined(_MPI)
            ! Compute the amount of sections of a rank from every rank
            num_sections_local(1) = sections_size
            ! Prepare for incoming data from other ranks
            allocate(num_sections(size_MPI), stat = error); assert_eq(error, 0)
            ! And distribute this information across all ranks
            call mpi_allgather(num_sections_local, 1, MPI_INTEGER4, num_sections, 1, MPI_INTEGER4, MPI_COMM_WORLD, error); assert_eq(error, 0)
            ! Then, find the global maximum
            max_num_sections = maxval(num_sections)

            ! Gather the amount of cells in all sections of this rank
            ! If this rank has less sections than max_num_sections, the remaining fields will be -1
            allocate(num_cells_local(max_num_sections), stat = error); assert_eq(error, 0)
            num_cells_local = -1
#           if defined(_XDMF_PATCH)
                allocate(num_patches_local(max_num_sections), stat = error); assert_eq(error, 0)
                num_patches_local = -1
#           endif
            if (.not. from_filter) then
                do j = 1, sections_size
                    section_info = grid%sections%elements(j)%get_info()
                    num_cells_local(j) = section_info%i_cells
#                   if defined(_XDMF_PATCH)
                        num_cells_local(j) = num_cells_local(j) * _XDMF_PATCH_ORDER_SQUARE
                        num_patches_local(j) = section_info%i_cells
#                   endif
                end do
            else
                do j = 1, sections_size
                    num_cells_local(j) = grid%sections%elements(j)%xdmf_filter_count_cells
#                   if defined(_XDMF_PATCH)
                        num_patches_local(j) = grid%sections%elements(j)%xdmf_filter_count_patches
#                   endif
                end do
            end if
            ! Prepare for incoming data from other ranks
            allocate(num_cells(max_num_sections, size_MPI), stat = error); assert_eq(error, 0)
            ! And distribute this information across all ranks
            call mpi_allgather(num_cells_local, max_num_sections, MPI_INTEGER4, num_cells, max_num_sections, MPI_INTEGER4, MPI_COMM_WORLD, error); assert_eq(error, 0)
#           if defined(_XDMF_PATCH)
                allocate(num_patches(max_num_sections, size_MPI), stat = error); assert_eq(error, 0)
                call mpi_allgather(num_patches_local, max_num_sections, MPI_INTEGER4, num_patches, max_num_sections, MPI_INTEGER4, MPI_COMM_WORLD, error); assert_eq(error, 0)
#           endif

            ! Gather the index of each section of this rank
            ! Again, if this rank has less sections than max_num_sections, the remaining fields will be -1
            allocate(sect_index_local(max_num_sections), stat = error); assert_eq(error, 0)
            sect_index_local = -1
            do j = 1, sections_size
                sect_index_local(j) = grid%sections%elements(j)%index
            end do
            ! Again, Prepare for incoming data from other ranks
            allocate(sect_index(max_num_sections, size_MPI), stat = error); assert_eq(error, 0)
            ! And, again, distribute this information across all ranks
            call mpi_allgather(sect_index_local, max_num_sections, MPI_INTEGER4, sect_index, max_num_sections, MPI_INTEGER4, MPI_COMM_WORLD, error); assert_eq(error, 0)

            ! Allocate the global section descriptor table
            call base%root_layout_desc%allocate(size_MPI)
            ! Fill the global section descriptor table by giving each section a offset and length in the big data base
            ! Also, the sections will be reordered according to their index given by the scheduler
            do i = 1, size_MPI
                call base%root_layout_desc%ranks(i)%allocate(num_sections(i))
                do j = 1, num_sections(i)
                    base%root_layout_desc%ranks(i)%sections(sect_index(j, i))%num_cells = num_cells(j, i)
#                   if defined(_XDMF_PATCH)
                        base%root_layout_desc%ranks(i)%sections(sect_index(j, i))%num_patches = num_patches(j, i)
#                   endif
                end do
            end do

            ! Based on the section lengths, compute their offsets in the data set
            offset_cells_current = 0
#           if defined(_XDMF_PATCH)
                offset_patches_current = 0
#           endif
            do i = 1, size_MPI
                offset_cells_current_buffer = 0
                num_cells_section_buffer = 0
                base%root_layout_desc%ranks(i)%offset_cells = offset_cells_current
#               if defined(_XDMF_PATCH)
                    offset_patches_current_buffer = 0
                    num_patches_section_buffer = 0
                    base%root_layout_desc%ranks(i)%offset_patches = offset_patches_current
#               endif
                do j = 1, num_sections(i)
                    num_cells_section = base%root_layout_desc%ranks(i)%sections(j)%num_cells
                    num_cells_section_buffer = num_cells_section_buffer + num_cells_section
                    if (num_cells_section .ne. 0) then
                        base%root_layout_desc%ranks(i)%sections(j)%offset_cells = offset_cells_current
                        base%root_layout_desc%ranks(i)%sections(j)%offset_cells_buffer = offset_cells_current_buffer
                        offset_cells_current = offset_cells_current + num_cells_section
                        offset_cells_current_buffer = offset_cells_current_buffer + num_cells_section
                    end if
#                   if defined(_XDMF_PATCH)
                        num_patches_section = base%root_layout_desc%ranks(i)%sections(j)%num_patches
                        num_patches_section_buffer = num_patches_section_buffer + num_patches_section
                        if (num_patches_section .ne. 0) then
                            base%root_layout_desc%ranks(i)%sections(j)%offset_patches = offset_patches_current
                            base%root_layout_desc%ranks(i)%sections(j)%offset_patches_buffer = offset_patches_current_buffer
                            offset_patches_current = offset_patches_current + num_patches_section
                            offset_patches_current_buffer = offset_patches_current_buffer + num_patches_section
                        end if
#                   endif
                end do
                base%root_layout_desc%ranks(i)%num_cells = num_cells_section_buffer
#               if defined(_XDMF_PATCH)
                    base%root_layout_desc%ranks(i)%num_patches = num_patches_section_buffer
#               endif
            end do

            ! Free memory from local MPI exchange arrays
            deallocate(num_cells, stat = error); assert_eq(error, 0)
            deallocate(num_cells_local, stat = error); assert_eq(error, 0)
#           if defined(_XDMF_PATCH)
                deallocate(num_patches, stat = error); assert_eq(error, 0)
                deallocate(num_patches_local, stat = error); assert_eq(error, 0)
#           endif
            deallocate(sect_index_local, stat = error); assert_eq(error, 0)
            deallocate(sect_index, stat = error); assert_eq(error, 0)
            deallocate(num_sections, stat = error); assert_eq(error, 0)
#       else
            ! No MPI, only one rank
            call base%root_layout_desc%allocate(1)
            call base%root_layout_desc%ranks(1)%allocate(sections_size)
            offset_cells_current = 0
            num_cells_section_buffer = 0
#           if defined(_XDMF_PATCH)
                offset_patches_current = 0
                num_patches_section_buffer = 0
#           endif
            base%root_layout_desc%ranks(1)%offset_cells = offset_cells_current
            do j = 1, sections_size
                section_index = grid%sections%elements(j)%index
                if (.not. from_filter) then
                    section_info = grid%sections%elements(j)%get_info()
                    num_cells_section = section_info%i_cells
#                   if defined(_XDMF_PATCH)
                        num_cells_section = num_cells_section * _XDMF_PATCH_ORDER_SQUARE
                        num_patches_section = section_info%i_cells
#                   endif
                else
                    num_cells_section = grid%sections%elements(j)%xdmf_filter_count_cells
#                   if defined(_XDMF_PATCH)
                        num_patches_section = grid%sections%elements(j)%xdmf_filter_count_patches
#                   endif
                end if
                ! Set section length and compute offset
                num_cells_section_buffer = num_cells_section_buffer + num_cells_section
                base%root_layout_desc%ranks(1)%sections(section_index)%num_cells = num_cells_section
                base%root_layout_desc%ranks(1)%sections(section_index)%offset_cells = offset_cells_current
                base%root_layout_desc%ranks(1)%sections(section_index)%offset_cells_buffer = offset_cells_current
                offset_cells_current = offset_cells_current + num_cells_section
#               if defined(_XDMF_PATCH)
                    num_patches_section_buffer = num_patches_section_buffer + num_patches_section
                    base%root_layout_desc%ranks(1)%sections(section_index)%num_patches = num_patches_section
                    base%root_layout_desc%ranks(1)%sections(section_index)%offset_patches = offset_patches_current
                    base%root_layout_desc%ranks(1)%sections(section_index)%offset_patches_buffer = offset_patches_current
                    offset_patches_current = offset_patches_current + num_patches_section
#               endif
            end do
            base%root_layout_desc%ranks(1)%num_cells = num_cells_section_buffer
#           if defined(_XDMF_PATCH)
                base%root_layout_desc%ranks(1)%num_patches = num_patches_section_buffer
#           endif
#       endif
    end subroutine

end module
