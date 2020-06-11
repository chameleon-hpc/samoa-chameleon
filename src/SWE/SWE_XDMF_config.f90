#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

! This module defines the XDMF configuration for SWE
#if defined(_SWE)
    module SWE_XDMF_Config

        use Samoa_swe
        use XDMF_data_types
        use XDMF_config
#       if defined(_SWE_DG)
            use SWE_DG_Limiter
#       endif

        use, intrinsic :: iso_fortran_env

        implicit none

        ! Integer value names
        character(len = 1), parameter			:: swe_hdf5_attr_depth_dname_nz = "d"
        character(len = 1), parameter			:: swe_hdf5_attr_rank_dname_nz = "o"
        character(len = 1), parameter			:: swe_hdf5_attr_plotter_dname_nz = "l"
        character(len = 1), parameter			:: swe_hdf5_attr_section_dname_nz = "s"
#       if defined(_SWE_DG)
            character(len = 1), parameter		:: swe_hdf5_attr_troubled_dname_nz = "t"
#       endif
        character(len = 2), parameter			:: swe_hdf5_attr_depth_dname = swe_hdf5_attr_depth_dname_nz//char(0)
        character(len = 2), parameter			:: swe_hdf5_attr_rank_dname = swe_hdf5_attr_rank_dname_nz//char(0)
        character(len = 2), parameter			:: swe_hdf5_attr_plotter_dname = swe_hdf5_attr_plotter_dname_nz//char(0)
        character(len = 2), parameter			:: swe_hdf5_attr_section_dname = swe_hdf5_attr_section_dname_nz//char(0)
#       if defined(_SWE_DG)
            character(len = 2), parameter	    :: swe_hdf5_attr_troubled_dname = swe_hdf5_attr_troubled_dname_nz//char(0)
#       endif

        ! Real value names
        character(len = 1), parameter			:: swe_hdf5_attr_b_dname_nz = "b"
        character(len = 1), parameter			:: swe_hdf5_attr_h_dname_nz = "h"
        character(len = 1), parameter			:: swe_hdf5_attr_bh_dname_nz = "k"
        character(len = 2), parameter			:: swe_hdf5_attr_b_dname = swe_hdf5_attr_b_dname_nz//char(0)
        character(len = 2), parameter			:: swe_hdf5_attr_h_dname = swe_hdf5_attr_h_dname_nz//char(0)
        character(len = 2), parameter			:: swe_hdf5_attr_bh_dname = swe_hdf5_attr_bh_dname_nz//char(0)

        ! 2D vector real value names
        character(len = 1), parameter			:: swe_hdf5_attr_f_dname_nz = "f"
        character(len = 2), parameter			:: swe_hdf5_attr_f_dname = swe_hdf5_attr_f_dname_nz//char(0)
       
        ! Convenient offsets
        integer, parameter                      :: swe_hdf5_valsi_plotter_offset = 3
        integer, parameter                      :: swe_hdf5_valsr_b_offset = 1, swe_hdf5_valsr_h_offset = 2, swe_hdf5_valsr_bh_offset = 3
        integer, parameter                      :: swe_hdf5_valsuv_f_offset = 1

        ! Parameter for the XDMF core API
        type(t_xdmf_parameter), save            :: swe_xdmf_param_cells = t_xdmf_parameter( &
            hdf5_valsg_width = 2, &     ! 2 geometry data fields: dimensions X and Y
            hdf5_valst_width = 3, &     ! 3 geometry entries per triangle: corners
            hdf5_attr_width = 3, &      ! 1 attribute structure instances per triangle: in the middle. Duplicate to all corners for compatibility with dg elements
#           if defined(_SWE_DG)
                hdf5_valsi_width = 5, &     ! 5 int32 values per cell: d, o, l, s, t
#           else
                hdf5_valsi_width = 4, &     ! 4 int32 values per cell: d, o, l, s
#           endif
            hdf5_valsuv_width = 1, &    ! 1 2D vector real32 values in attribute structure: f
            hdf5_valsr_width = 3 &      ! 3 real32 values in attribute structure: b, h, k
        )

#       if defined (_SWE_PATCH)
            type(t_xdmf_parameter), save        :: swe_xdmf_param_patches = t_xdmf_parameter( &
                hdf5_valsg_width = 2, &     ! 2 geometry data fields: dimensions X and Y
                hdf5_valst_width = 3, &     ! 3 geometry entries per triangle: corners
#           if defined(_SWE_DG)
                hdf5_attr_width = _SWE_DG_DOFS, &      ! _SWE_DG_DOFS attribute structure instances per triangle: at the sampling nodes
#           else
                hdf5_attr_width = 1, &      ! fallback 1 attribute structure instances per triangle: in the middle
#           endif
                hdf5_valsi_width = 5, &     ! 5 int32 values per cell: d, o, l, s, t
                hdf5_valsuv_width = 1, &    ! 1 2D vector real32 values in attribute structure: f
                hdf5_valsr_width = 3 &      ! 3 real32 values in attribute structure: b, h, k
            )
#       endif

        ! Max amount of probes
        integer, parameter                      :: xdmf_filter_max_probes = (xdmf_filter_params_width - 1) / 2

        contains

        ! This routine loads the XDMF config for SWE
        subroutine SWE_xdmf_config_load()
            call swe_xdmf_param_cells%allocate()
            call assign_names(swe_xdmf_param_cells)
#           if defined (_SWE_PATCH)
                call swe_xdmf_param_patches%allocate()
                call assign_names(swe_xdmf_param_patches)
#           endif
        end subroutine

        subroutine assign_names(param)
            type(t_xdmf_parameter), intent(inout)   :: param

            param%hdf5_valsi_dnames = &
                (/ swe_hdf5_attr_depth_dname_nz, swe_hdf5_attr_rank_dname_nz, swe_hdf5_attr_plotter_dname_nz, swe_hdf5_attr_section_dname_nz &
#           if defined(_SWE_DG)
                    , swe_hdf5_attr_troubled_dname_nz & 
#           endif
                /)
            param%hdf5_valsr_dnames = &
                (/ swe_hdf5_attr_b_dname_nz, swe_hdf5_attr_h_dname_nz, swe_hdf5_attr_bh_dname_nz /)
            param%hdf5_valsuv_dnames = &
                (/ swe_hdf5_attr_f_dname_nz /)
        end subroutine

        ! Subsection filter routines
        function SWE_xdmf_filter(element, index, argc, argv, is_cell_layer)
            type(t_element_base), intent(in) :: element
            integer, intent(in) 	            :: index, argc
            real, dimension(xdmf_filter_params_width), intent(in) :: argv
            logical, intent(in)                 :: is_cell_layer
            logical                             :: SWE_xdmf_filter, res

            ! Evaluate hybrid compositing filter
            if(iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_all) .eq. xdmf_output_mode_all) then
                res = SWE_xdmf_filter_hybrid_selector(element)
                if((is_cell_layer .and. .not. res) .or. (.not. is_cell_layer .and. res)) then
                    SWE_xdmf_filter = .false.
                    return
                end if
            end if

            ! Evaluate regular filter
            select case (index)
                case (1)
                    SWE_xdmf_filter = SWE_xdmf_filter_rect(element, argv(1), argv(2), argv(3), argv(4))
                case (2)
                    SWE_xdmf_filter = SWE_xdmf_filter_h_treshold(element, argv(1), argv(2))
                case (3)
                    SWE_xdmf_filter = SWE_xdmf_filter_probes(element, argv(1), (argc - 1) / 2, &
                        reshape(argv(2:), (/ 2, xdmf_filter_max_probes /), (/ 0.0 /)))
                case default
                    SWE_xdmf_filter = .true.
            end select
        end function

        ! Hybrid layer selector
        function SWE_xdmf_filter_hybrid_selector(element)
            type(t_element_base), intent(in) :: element
            logical                             :: SWE_xdmf_filter_hybrid_selector

            SWE_xdmf_filter_hybrid_selector = element%cell%data_pers%troubled .eq. COAST
            ! SWE_xdmf_filter_hybrid_selector = SWE_xdmf_filter_probes(element, 1.0, 1, &
            !     reshape((/ 0.0, 2.0 /), (/ 2, xdmf_filter_max_probes /), (/ 0.0 /)))
        end function

        ! Rectangular selection filter
        function SWE_xdmf_filter_rect(element, xmin, xmax, ymin, ymax)
            type(t_element_base), intent(in) :: element
            real, intent(in)                    :: xmin, xmax, ymin, ymax
            logical                             :: SWE_xdmf_filter_rect

            integer                             :: i
            real(XDMF_GRID_SR), dimension(2)	:: coord

            do i = 1, 3
                ! Compute the cells cartesian position
                coord = real((cfg%scaling * element%nodes(i)%ptr%position) + cfg%offset(:), XDMF_GRID_SR)
                ! Compare each vertex to selection
                if (coord(1) .ge. xmin .and. coord(1) .le. xmax .and. &
                    coord(2) .ge. ymin .and. coord(2) .le. ymax) then
                        SWE_xdmf_filter_rect = .true.
                        return
                end if
            end do
            SWE_xdmf_filter_rect = .false.
        end function

        ! Water height treshold selection filter
        function SWE_xdmf_filter_h_treshold(element, min, max)
            type(t_element_base), intent(in) :: element
            real, intent(in)                    :: min, max
            logical                             :: SWE_xdmf_filter_h_treshold

            integer                             :: i
            real(XDMF_GRID_SR)		            :: h

            do i = 1, 3
                ! Get cells water height
                h = real(element%cell%data_pers%Q(i)%h, XDMF_GRID_SR)
                ! Compare each vertex to selection
                if (h .ge. min .and. h .le. max) then
                    SWE_xdmf_filter_h_treshold = .true.
                    return
                end if
            end do
            SWE_xdmf_filter_h_treshold = .false.
        end function

        ! Include probes output filter
        function SWE_xdmf_filter_probes(element, radius, num_probes, probes)
            type(t_element_base), intent(in) :: element
            real, intent(in)                    :: radius
            integer, intent(in)                 :: num_probes
            real, dimension(2, xdmf_filter_max_probes), intent(in) :: probes
            logical                             :: SWE_xdmf_filter_probes

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
                            SWE_xdmf_filter_probes = .true.
                            return
                    end  if
                end do
            end do
            SWE_xdmf_filter_probes = .false.
        end function

    end module
#endif