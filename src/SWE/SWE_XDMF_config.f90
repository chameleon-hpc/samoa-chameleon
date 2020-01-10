#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

! This module defines the XDMF configuration for SWE
#if defined(_SWE)
    module SWE_XDMF_Config

        use Samoa_swe
        use XDMF_data_types
        use XDMF_config

        use, intrinsic :: iso_fortran_env

        implicit none

        ! Integer value names
        character(len = 1), parameter			:: swe_hdf5_attr_depth_dname_nz = "d"
        character(len = 1), parameter			:: swe_hdf5_attr_rank_dname_nz = "o"
        character(len = 1), parameter			:: swe_hdf5_attr_plotter_dname_nz = "l"
        character(len = 1), parameter			:: swe_hdf5_attr_section_dname_nz = "s"
        character(len = 2), parameter			:: swe_hdf5_attr_depth_dname = swe_hdf5_attr_depth_dname_nz//char(0)
        character(len = 2), parameter			:: swe_hdf5_attr_rank_dname = swe_hdf5_attr_rank_dname_nz//char(0)
        character(len = 2), parameter			:: swe_hdf5_attr_plotter_dname = swe_hdf5_attr_plotter_dname_nz//char(0)
        character(len = 2), parameter			:: swe_hdf5_attr_section_dname = swe_hdf5_attr_section_dname_nz//char(0)

        ! Real value names
        character(len = 1), parameter			:: swe_hdf5_attr_b_dname_nz = "b"
        character(len = 1), parameter			:: swe_hdf5_attr_h_dname_nz = "h"
        character(len = 1), parameter			:: swe_hdf5_attr_bh_dname_nz = "k"
        character(len = 2), parameter			:: swe_hdf5_attr_b_dname = swe_hdf5_attr_b_dname_nz//char(0)
        character(len = 2), parameter			:: swe_hdf5_attr_h_dname = swe_hdf5_attr_h_dname_nz//char(0)
        character(len = 2), parameter			:: swe_hdf5_attr_bh_dname = swe_hdf5_attr_bh_dname_nz//char(0)

        ! 2D vector real value names
        character(len = 1), parameter			:: swe_hdf5_attr_f_dname_nz = "f"
        character(len = 1), parameter			:: swe_hdf5_attr_g_dname_nz = "g"
        character(len = 2), parameter			:: swe_hdf5_attr_f_dname = swe_hdf5_attr_f_dname_nz//char(0)
        character(len = 2), parameter			:: swe_hdf5_attr_g_dname = swe_hdf5_attr_g_dname_nz//char(0)
       
        ! Convenient offsets
        integer, parameter                      :: swe_hdf5_valsi_plotter_offset = 3
        integer, parameter                      :: swe_hdf5_valsr_b_offset = 1, swe_hdf5_valsr_bh_offset = 3
        integer, parameter                      :: swe_hdf5_valsuv_f_offset = 1, swe_hdf5_valsuv_g_offset = 2

        ! Parameter for the XDMF core API
        type(t_xdmf_parameter), save            :: swe_xdmf_param = t_xdmf_parameter( &
            hdf5_valsg_width = 2, &     ! 2 geometry data fields: dimensions X and Y
            hdf5_valst_width = 3, &     ! 3 geometry entries per triangle: corners
            hdf5_attr_width = 3, &      ! 3 attribute structure instances per triangle: one per corner 
            hdf5_valsi_width = 4, &     ! 4 int32 values per cell: d, o, l, s
            hdf5_valsuv_width = 2, &    ! 2 2D vector real32 values in attribute structure: f, g
            hdf5_valsr_width = 3 &      ! 3 real32 values in attribute structure: b, h, k
        )

        ! Max amount of probes
        integer, parameter                      :: xdmf_filter_max_probes = (xdmf_filter_params_width - 1) / 2

        contains

        ! This routine loads the XDMF config for SWE
        subroutine SWE_xdmf_config_load()
            call swe_xdmf_param%allocate()
          
            swe_xdmf_param%hdf5_valsi_dnames = &
                (/ swe_hdf5_attr_depth_dname_nz, swe_hdf5_attr_rank_dname_nz, swe_hdf5_attr_plotter_dname_nz, swe_hdf5_attr_section_dname_nz /)
            swe_xdmf_param%hdf5_valsr_dnames = &
                (/ swe_hdf5_attr_b_dname_nz, swe_hdf5_attr_h_dname_nz, swe_hdf5_attr_bh_dname_nz /)
            swe_xdmf_param%hdf5_valsuv_dnames = &
                (/ swe_hdf5_attr_f_dname_nz, swe_hdf5_attr_g_dname_nz /)
        end subroutine

        ! Subsection filter routines
        subroutine SWE_xdmf_filter(element, index, argc, argv, result)
            type(t_element_base), intent(inout) :: element
            integer, intent(in) 	            :: index, argc
            real, dimension(xdmf_filter_params_width), intent(in) :: argv
            logical, intent(out)                :: result

            select case (index)
                case (1)
                    call SWE_xdmf_filter_rect(element, argv(1), argv(2), argv(3), argv(4), result)
                case (2)
                    call SWE_xdmf_filter_h_treshold(element, argv(1), argv(2), result)
                case (3)
                    call SWE_xdmf_filter_probes(element, argv(1), (argc - 1) / 2, &
                        reshape(argv(2:), (/ 2, xdmf_filter_max_probes /), (/ 0.0 /)), result)
                case default
                    result = .true.
            end select
        end subroutine

        ! Rectangular selection filter
        subroutine SWE_xdmf_filter_rect(element, xmin, xmax, ymin, ymax, result)
            type(t_element_base), intent(inout) :: element
            real, intent(in)                    :: xmin, xmax, ymin, ymax
            logical, intent(out)                :: result

            integer                             :: i
            real(XDMF_GRID_SR), dimension(2)	:: coord

            do i = 1, 3
                ! Compute the cells cartesian position
                coord = real((cfg%scaling * element%nodes(i)%ptr%position) + cfg%offset(:), XDMF_GRID_SR)
                ! Compare each vertex to selection
                if ((coord(1) .ge. xmin .and. coord(1) .le. xmax) .and. &
                    (coord(2) .ge. ymin .and. coord(2) .le. ymax)) then
                        result = .true.
                        return
                end if
            end do
            result = .false.
        end subroutine

        ! Water height treshold selection filter
        subroutine SWE_xdmf_filter_h_treshold(element, min, max, result)
            type(t_element_base), intent(inout) :: element
            real, intent(in)                    :: min, max
            logical, intent(out)                :: result

            integer                             :: i
            real(XDMF_GRID_SR)		            :: h

            do i = 1, 3
                ! Get cells water height
                h = real(element%cell%data_pers%Q(i)%h, XDMF_GRID_SR)
                ! Compare each vertex to selection
                if (h .ge. min .and. h .le. max) then
                    result = .true.
                    return
                end if
            end do
            result = .false.
        end subroutine

        ! Include probes output filter
        subroutine SWE_xdmf_filter_probes(element, radius, num_probes, probes, result)
            type(t_element_base), intent(inout) :: element
            real, intent(in)                    :: radius
            integer, intent(in)                 :: num_probes
            real, dimension(2, xdmf_filter_max_probes), intent(in) :: probes
            logical, intent(out)                :: result

            integer                             :: i, j
            real(XDMF_GRID_SR), dimension(2)	:: coord
            real                                :: radius_sq

            radius_sq = radius ** 2

            do i = 1, 3
                ! Compute the cells cartesian position
                coord = real((cfg%scaling * element%nodes(i)%ptr%position) + cfg%offset(:), XDMF_GRID_SR)
                ! Compare each vertex to probe locations
                do j = 1, num_probes
                    if((((coord(1) - probes(1, j)) ** 2) + &
                        ((coord(2) - probes(2, j)) ** 2)) .le. radius_sq) then
                            result = .true.
                            return
                    end  if
                end do
            end do
            result = .false.
        end subroutine

    end module
#endif