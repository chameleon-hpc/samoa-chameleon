#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

! This module defines the XDMF input adaption traversal
! It is essentially empty and exists only to trigger the splitting of cells
#if defined(_SWE)
    module SWE_XDMF_Adapt
        use SFC_edge_traversal
        use Conformity

        use Samoa_swe
        use Tools_noise
        use SWE_heun1_timestep

#       if defined(_SWE_PATCH)
            use SWE_PATCH
#       endif
#       if defined(_SWE_DG)
            use SWE_dg_matrices
            use SWE_data_types
#       endif   
    
        implicit none

        type num_traversal_data
        end type
        
#       define _GT_EDGES
#		define _GT_NAME							t_swe_xdmf_adaption_traversal

#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP			post_traversal_op

#		define _GT_TRANSFER_OP					transfer_op
#		define _GT_REFINE_OP					refine_op
#		define _GT_COARSEN_OP					coarsen_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op

#		include "SFC_generic_adaptive_traversal.f90"

        subroutine pre_traversal_grid_op(traversal, grid)
            type(t_swe_xdmf_adaption_traversal), intent(inout)		        :: traversal
            type(t_grid), intent(inout)							            :: grid
        end subroutine

        subroutine post_traversal_grid_op(traversal, grid)
            type(t_swe_xdmf_adaption_traversal), intent(inout)				:: traversal
            type(t_grid), intent(inout)							            :: grid
        end subroutine

        subroutine pre_traversal_op(traversal, section)
            type(t_swe_xdmf_adaption_traversal), intent(inout)				:: traversal
            type(t_grid_section), intent(inout)							    :: section
        end subroutine

        subroutine post_traversal_op(traversal, section)
            type(t_swe_xdmf_adaption_traversal), intent(inout)	            :: traversal
            type(t_grid_section), intent(inout)							    :: section
        end subroutine

        subroutine transfer_op(traversal, section, src_element, dest_element)
             type(t_swe_xdmf_adaption_traversal), intent(inout)	             :: traversal
            type(t_grid_section), intent(inout)								 :: section
            type(t_traversal_element), intent(inout)						 :: src_element
            type(t_traversal_element), intent(inout)						 :: dest_element
        end subroutine

        subroutine refine_op(traversal, section, src_element, dest_element, refinement_path)
            type(t_swe_xdmf_adaption_traversal), intent(inout)	             :: traversal
            type(t_grid_section), intent(inout)								 :: section
            type(t_traversal_element), intent(inout)						 :: src_element
            type(t_traversal_element), intent(inout)						 :: dest_element
            integer, dimension(:), intent(in)								 :: refinement_path
        end subroutine

        subroutine coarsen_op(traversal, section, src_element, dest_element, refinement_path)
            type(t_swe_xdmf_adaption_traversal), intent(inout)	             :: traversal
            type(t_grid_section), intent(inout)								 :: section
            type(t_traversal_element), intent(inout)						 :: src_element
            type(t_traversal_element), intent(inout)						 :: dest_element
            integer, dimension(:), intent(in)								 :: refinement_path
        end subroutine

    end module
#endif

