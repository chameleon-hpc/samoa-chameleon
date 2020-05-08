! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
MODULE SWE_Displace
  use SFC_edge_traversal
  use SWE_initialize_bathymetry
  use SWE_PATCH
  use SWE_DG_solver
  use SWE_DG_Limiter
  use Samoa_swe

  implicit none

  type num_traversal_data
     integer (kind = GRID_DI)			:: i_refinements_issued
  end type num_traversal_data

  type(t_gv_Q)							:: gv_Q
  type(t_lfs_flux)						:: lfs_flux

#		define _GT_NAME							t_swe_displace_traversal
#		define _GT_EDGES

#		define _GT_ELEMENT_OP				element_op
#		define _GT_CELL_TO_EDGE_OP  cell_to_edge_op_dg

#		include "SFC_generic_traversal_ringbuffer.f90"

  
  !******************
  !Geometry operators
  !******************
  
  subroutine element_op(traversal, section, element)
    type(t_swe_displace_traversal), intent(inout)	:: traversal
    type(t_grid_section), intent(inout)					  :: section
    type(t_element_base), intent(inout)					  :: element
    call alpha_volume_op(traversal, section, element)
  end subroutine element_op

  !*******************************
  !Volume and DoF operators
  !*******************************

  subroutine alpha_volume_op(traversal, section, element)
    type(t_swe_displace_traversal), intent(inout)		:: traversal
    type(t_grid_section), intent(inout)							:: section
    type(t_element_base), intent(inout)						  :: element
    
    !local variables
    real (kind = GRID_SR), dimension(_SWE_PATCH_ORDER_SQUARE) :: b_new

    !evaluate initial function values at dof positions and compute DoFs
#if defined(_ASAGI)
    if(.not.isCoast(element%cell%data_pers%troubled)) then
       element%cell%data_pers%Q%B = get_bathymetry_at_dg_patch(section, element, section%r_time)
    else
       b_new = get_bathymetry_at_patch(section, element, section%r_time)
       element%cell%data_pers%H   =  element%cell%data_pers%H - element%cell%data_pers%B + b_new
       element%cell%data_pers%B   =  b_new
    end if
#endif
    !no coarsening while the earthquake takes place
    !element%cell%geometry%refinement = max(0, element%cell%geometry%refinement)

  end subroutine alpha_volume_op
END MODULE SWE_Displace
#endif
