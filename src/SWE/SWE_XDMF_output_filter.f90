#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

#if defined(_FLASH)
    module FLASH_XDMF_output_filter
        
        use SFC_edge_traversal
        use Samoa_flash
        use FLASH_heun1_timestep
        use Tools_openmp

        use FLASH_XDMF_config
        use XDMF_output_base

        implicit none

        type num_traversal_data
            type(t_xdmf_base_output_filter_traversal)   :: base
        end type

#       define _GT_EDGES
#       if defined(LIMITER_BJ_EDGE)
#           define _GT_EDGE_MPI_TYPE
#       elif defined(LIMITER_BJ_VERTEX)
#           define _GT_NODE_MPI_TYPE
#       endif  

#		define _GT_NAME								    t_flash_xdmf_output_filter_traversal

#		define _GT_PRE_TRAVERSAL_OP					    pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP				    post_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP				pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP				post_traversal_grid_op

#		define _GT_ELEMENT_OP						    element_op

#       define _GT_CELL_TO_EDGE_OP				        cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

        subroutine pre_traversal_grid_op(traversal, grid)
            type(t_flash_xdmf_output_filter_traversal), intent(inout)		:: traversal
            type(t_grid), intent(inout)							            :: grid

            integer                                                         :: i

            do i = 1, size(traversal%sections)
                traversal%sections(i)%base%output_iteration = traversal%base%output_iteration
            end do
        end subroutine

        subroutine post_traversal_grid_op(traversal, grid)
            type(t_flash_xdmf_output_filter_traversal), intent(inout)		:: traversal
            type(t_grid), intent(inout)							            :: grid

            traversal%base%output_iteration = traversal%base%output_iteration + 1
        end subroutine

        subroutine pre_traversal_op(traversal, section)
            type(t_flash_xdmf_output_filter_traversal), intent(inout)		:: traversal
            type(t_grid_section), intent(inout)							    :: section

            section%xdmf_filter_count = 0
        end subroutine

        subroutine post_traversal_op(traversal, section)
            type(t_flash_xdmf_output_filter_traversal), intent(inout)		:: traversal
            type(t_grid_section), intent(inout)							    :: section

        end subroutine

        subroutine element_op(traversal, section, element)
            type(t_flash_xdmf_output_filter_traversal), intent(inout)		:: traversal
            type(t_grid_section), intent(inout)							    :: section
            type(t_element_base), intent(inout)					            :: element

            logical                                                         :: write_cp, filter_result = .true.

            ! Compute whether to output tree (checkpoint) data
            if(cfg%xdmf%i_xdmfcpint.eq.0) then
                write_cp = .false.
            else
                write_cp = mod(int(traversal%base%output_iteration, GRID_SI), int(cfg%xdmf%i_xdmfcpint, GRID_SI)).eq.0
            end if

            ! Evaluate filter if this step is not a checkpoint
            if ((.not. write_cp) .and. (cfg%xdmf%i_xdmffilter_index .ne. 0)) then
                call flash_xdmf_filter(element, cfg%xdmf%i_xdmffilter_index, cfg%xdmf%i_xmdffilter_params_count, &
                    cfg%xdmf%r_xdmffilter_params_vector, filter_result)
            end if

            ! Update cell count
            if (write_cp .or. filter_result) then
                section%xdmf_filter_count = section%xdmf_filter_count + 1
            end if
        end subroutine

    end module
#endif
