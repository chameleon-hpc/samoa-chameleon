#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

#if defined(_SWE)
    module SWE_XDMF_output_filter
        
        use SFC_edge_traversal
        use Samoa_swe
        use Tools_openmp

        use SWE_XDMF_config
        use XDMF_output_base

#       if defined(_SWE_PATCH)
            use SWE_PATCH
#       endif
#       if defined(_SWE_DG)
            use SWE_dg_matrices
            use SWE_data_types
            use SWE_dg_solver
#       endif   

        implicit none

        type num_traversal_data
            type(t_xdmf_base_output_filter_traversal)   :: base
        end type

#       define _GT_EDGES
#		define _GT_NAME								    t_swe_xdmf_output_filter_traversal

#		define _GT_PRE_TRAVERSAL_OP					    pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP				    post_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP				pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP				post_traversal_grid_op

#		define _GT_ELEMENT_OP						    element_op

#       define _GT_CELL_TO_EDGE_OP				        cell_to_edge_op_dg

#		include "SFC_generic_traversal_ringbuffer.f90"

        subroutine pre_traversal_grid_op(traversal, grid)
            type(t_swe_xdmf_output_filter_traversal), intent(inout)		:: traversal
            type(t_grid), intent(inout)							            :: grid

            integer                                                         :: i

            do i = 1, size(traversal%sections)
                traversal%sections(i)%base%output_iteration = traversal%base%output_iteration
            end do
        end subroutine

        subroutine post_traversal_grid_op(traversal, grid)
            type(t_swe_xdmf_output_filter_traversal), intent(inout)		:: traversal
            type(t_grid), intent(inout)							            :: grid

            traversal%base%output_iteration = traversal%base%output_iteration + 1
        end subroutine

        subroutine pre_traversal_op(traversal, section)
            type(t_swe_xdmf_output_filter_traversal), intent(inout)		:: traversal
            type(t_grid_section), intent(inout)							    :: section

            section%xdmf_filter_count_cells = 0
#           if defined(_SWE_PATCH)
                section%xdmf_filter_count_patches = 0
#           endif
        end subroutine

        subroutine post_traversal_op(traversal, section)
            type(t_swe_xdmf_output_filter_traversal), intent(inout)		:: traversal
            type(t_grid_section), intent(inout)							    :: section

        end subroutine

        subroutine element_op(traversal, section, element)
            type(t_swe_xdmf_output_filter_traversal), intent(inout)		:: traversal
            type(t_grid_section), intent(inout)							    :: section
            type(t_element_base), intent(inout)					            :: element

            logical                                                         :: write_cp, filter_result_cells
#           if defined(_SWE_PATCH)
                logical                                                     :: filter_result_patches
#           endif

            filter_result_cells = .true.
#           if defined(_SWE_PATCH)
                filter_result_patches = .true.
#           endif

            ! Compute whether to output tree (checkpoint) data
            if(cfg%xdmf%i_xdmfcpint.eq.0) then
                write_cp = .false.
            else
                write_cp = mod(int(traversal%base%output_iteration, GRID_SI), int(cfg%xdmf%i_xdmfcpint, GRID_SI)).eq.0
            end if

            ! Evaluate filter if this step is not a checkpoint
            if ((.not. write_cp) .and. ((cfg%xdmf%i_xdmffilter_index .ne. 0) .or. &
                (iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_all) .eq. xdmf_output_mode_all))) then
                if (iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_cells) .eq. xdmf_output_mode_cells) then
                    filter_result_cells = SWE_xdmf_filter(element, cfg%xdmf%i_xdmffilter_index, cfg%xdmf%i_xmdffilter_params_count, &
                        cfg%xdmf%r_xdmffilter_params_vector, .true.)
                end if
#               if defined(_SWE_PATCH)
                    if (iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_patches) .eq. xdmf_output_mode_patches) then
                        filter_result_patches = SWE_xdmf_filter(element, cfg%xdmf%i_xdmffilter_index, cfg%xdmf%i_xmdffilter_params_count, &
                            cfg%xdmf%r_xdmffilter_params_vector, .false.)
                    end if
#               endif
            end if

            ! Update cell count
            if (write_cp .or. filter_result_cells .or. filter_result_patches) then
                if (write_cp .or. filter_result_cells) then
#                   if defined(_SWE_PATCH)
                        section%xdmf_filter_count_cells = section%xdmf_filter_count_cells + _SWE_PATCH_ORDER_SQUARE
#                   else
                        section%xdmf_filter_count_cells = section%xdmf_filter_count_cells + 1
#                   endif
                end if
#               if defined(_SWE_PATCH)
                    if (write_cp .or. filter_result_patches) then
                        section%xdmf_filter_count_patches = section%xdmf_filter_count_patches + 1
                    end if
#               endif
            end if
        end subroutine

    end module
#endif
