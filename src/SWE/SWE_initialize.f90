
! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
MODULE SWE_Initialize_Bathymetry
  use iso_c_binding
  use SFC_edge_traversal

  use Samoa_swe
  use SWE_PATCH
  use SWE_DG_Limiter
  
  ! No ASAGI -> Artificial scenario selector 
#       if !defined(_ASAGI)
  use SWE_Scenario
#       endif 

  implicit none

  type num_traversal_data
  end type num_traversal_data

  type(t_gv_Q)							:: gv_Q
  type(t_lfs_flux)						:: lfs_flux

  PUBLIC get_bathymetry_at_element, get_bathymetry_at_position
  PUBLIC get_bathymetry_at_patch
  PUBLIC get_bathymetry_at_dg_patch

#define _GT_NAME                    t_swe_init_b_traversal
#define _GT_EDGES                    
#define _GT_EDGES_TEMP

#define _GT_PRE_TRAVERSAL_OP        pre_traversal_op
#define _GT_PRE_TRAVERSAL_GRID_OP   pre_traversal_grid_op
#define _GT_POST_TRAVERSAL_GRID_OP  post_traversal_grid_op
#define _GT_ELEMENT_OP              element_op

#include "SFC_generic_traversal_ringbuffer.f90"

  subroutine pre_traversal_grid_op(traversal, grid)
    type(t_swe_init_b_traversal), intent(inout)   :: traversal
    type(t_grid), intent(inout)			  :: grid

    grid%r_time = 0.0_GRID_SR
    grid%r_dt = 0.0_GRID_SR
    grid%r_dt_new = 0.0_GRID_SR

    call scatter(grid%r_time, grid%sections%elements_alloc%r_time)
  end subroutine pre_traversal_grid_op

  subroutine post_traversal_grid_op(traversal, grid)
    type(t_swe_init_b_traversal), intent(inout)   :: traversal
    type(t_grid), intent(inout)	                  :: grid

    call reduce(grid%r_dt_new, grid%sections%elements_alloc%r_dt_new, MPI_MIN, .true.)
    grid%r_dt = grid%r_dt_new

#if defined(_ASAGI)
    !----- TODO: refactor this -----!
    ! get new minimal an maximal bathymetry
    call reduce(grid%b_min, grid%sections%elements_alloc%b_min_new, MPI_MIN, .true.)
    call reduce(grid%b_max, grid%sections%elements_alloc%b_max_new, MPI_MAX, .true.)

    ! scatter to sections
    call scatter(grid%b_min, grid%sections%elements_alloc%b_min)
    call scatter(grid%b_max, grid%sections%elements_alloc%b_max)
    !-------------------------------!
#endif !_ASAGI
  end subroutine post_traversal_grid_op

  subroutine pre_traversal_op(traversal, section)
    type(t_swe_init_b_traversal), intent(inout)  :: traversal
    type(t_grid_section), intent(inout)		 :: section
    
    section%r_dt_new = huge(1.0_GRID_SR)
  end subroutine pre_traversal_op

  !******************
  !Geometry operators
  !******************

  subroutine element_op(traversal, section, element)
    type(t_swe_init_b_traversal), intent(inout)               :: traversal
    type(t_grid_section), intent(inout)                       :: section
    type(t_element_base), intent(inout)                       :: element
    type(t_state), dimension(_SWE_PATCH_ORDER_SQUARE)	      :: Q

    call alpha_volume_op(traversal, section, element, Q)

    element%cell%data_pers%B        = Q(:)%b
    element%cell%data_pers%troubled = WET_DRY_INTERFACE

    
#if defined(_ASAGI)
    section%b_min_new=min( section%b_min_new, minval(element%cell%data_pers%B))
    section%b_max_new=max( section%b_max_new, maxval(element%cell%data_pers%B))    
#endif !_ASAGI
  end subroutine element_op

  !*******************************
  !Volume and DoF operators
  !*******************************

  subroutine alpha_volume_op(traversal, section, element, Q)
    type(t_swe_init_b_traversal), intent(inout)		            :: traversal
    type(t_grid_section), intent(inout)				    :: section
    type(t_element_base), intent(inout)				    :: element
    type(t_state), dimension(_SWE_PATCH_ORDER_SQUARE), intent(out)  :: Q

    Q%b = get_bathymetry_at_patch(section, element, section%r_time)
  end subroutine alpha_volume_op

  recursive function refine_2D_recursive(section, x1, x2, x3, t, depth) result(bath)
    type(t_grid_section), intent(inout)     :: section
    real (kind = GRID_SR), intent(in)       :: x1(:), x2(:), x3(:), t
    real (kind = GRID_SR)                   :: bath
    integer, intent(in)                     :: depth

    if (depth > 0) then
       bath = 0.5_SR * (  refine_2D_recursive(section, x1, 0.5_SR * (x1 + x3), x2, t, depth - 1) &
            + refine_2D_recursive(section, x2, 0.5_SR * (x1 + x3), x3, t, depth - 1))
    else
       bath = get_bathymetry_at_position(section, (x1 + x2 + x3) / 3.0_SR, t)
    end if
  end function refine_2D_recursive

  function get_bathymetry_at_element(section, element, t) result(bathymetry)
    type(t_grid_section), intent(inout)     :: section
    type(t_element_base), intent(inout)     :: element
    real (kind = GRID_SR), intent(in)       :: t
    real (kind = GRID_SR)                   :: bathymetry

    real (kind = GRID_SR)   :: x(2), x1(2), x2(2), x3(2)
    integer                 :: ddepth, data_depth

#if defined(_ADAPT_INTEGRATE)
    !limit to 16 refinement levels and maximum depth + 3
#if defined(_ASAGI)
    data_depth = nint(log(cfg%scaling ** 2 / (asagi_grid_delta(cfg%afh_bathymetry, 0) * asagi_grid_delta(cfg%afh_bathymetry, 1))) / log(2.0_SR))
    ddepth = min(16, min(cfg%i_max_depth + 3, data_depth + 3) - element%cell%geometry%i_depth)
#else
    ddepth = min(16, cfg%i_max_depth - element%cell%geometry%i_depth)
#endif

    x1 = samoa_barycentric_to_world_point(element%transform_data, [1.0_SR, 0.0_SR])
    x2 = samoa_barycentric_to_world_point(element%transform_data, [0.0_SR, 0.0_SR])
    x3 = samoa_barycentric_to_world_point(element%transform_data, [0.0_SR, 1.0_SR])

    bathymetry = refine_2D_recursive(section, x1, x2, x3, t, ddepth)
#elif defined(_ADAPT_SAMPLE)
    x = samoa_barycentric_to_world_point(element%transform_data, [1.0_SR / 3.0_SR, 1.0_SR / 3.0_SR])
    bathymetry = get_bathymetry_at_position(section, x, t)
#endif
  end function get_bathymetry_at_element


  function get_bathymetry_at_patch(section, element, t) result(bathymetry)
    type(t_grid_section), intent(inout)     :: section
    type(t_element_base), intent(inout)     :: element
    real (kind = GRID_SR), intent(in)       :: t
    real (kind = GRID_SR), dimension(_SWE_PATCH_ORDER_SQUARE) :: bathymetry
    integer (kind = GRID_SI)  :: i, j, row, col, cell_id

    real (kind = GRID_SR)     :: x(2), x1(2), x2(2), x3(2)
    integer                   :: ddepth, data_depth

    !iterate through cells in patch
    row = 1
    col = 1
    do i=1, _SWE_PATCH_ORDER_SQUARE
       ! if orientation is backwards, the plotter uses a transformation that mirrors the cell...
       ! this simple change solves the problem :)
       if (element%cell%geometry%i_plotter_type > 0) then 
          cell_id = i
       else
          cell_id = (row-1)*(row-1) + 2 * row - col
       end if

#if defined(_ADAPT_INTEGRATE)
       !limit to 16 refinement levels and maximum depth + 3
#if defined(_ASAGI)
       data_depth = nint(log(cfg%scaling ** 2 / (asagi_grid_delta(cfg%afh_bathymetry, 0) * asagi_grid_delta(cfg%afh_bathymetry, 1))) / log(2.0_SR))
       ddepth = min(16, min(cfg%i_max_depth + 3, data_depth + 3) - element%cell%geometry%i_depth)
#else  !_ASAGI
       ddepth = min(16, cfg%i_max_depth - element%cell%geometry%i_depth)
#endif !_ASAGI

       if (mod(col,2) == 1) then 
          x1 = samoa_barycentric_to_world_point(element%transform_data, SWE_PATCH_geometry%coords(:,1,cell_id) )
          x2 = samoa_barycentric_to_world_point(element%transform_data, SWE_PATCH_geometry%coords(:,2,cell_id) )
          x3 = samoa_barycentric_to_world_point(element%transform_data, SWE_PATCH_geometry%coords(:,3,cell_id) )
       else
          x3 = samoa_barycentric_to_world_point(element%transform_data, SWE_PATCH_geometry%coords(:,1,cell_id) )
          x2 = samoa_barycentric_to_world_point(element%transform_data, SWE_PATCH_geometry%coords(:,2,cell_id) )
          x1 = samoa_barycentric_to_world_point(element%transform_data, SWE_PATCH_geometry%coords(:,3,cell_id) )
       end if

       bathymetry(i) = refine_2D_recursive(section, x1, x2, x3, t, ddepth)
!_ADAPT       
#elif defined(_ADAPT_SAMPLE) 
       ! x = average of the 3 vertices
       x = 0
       do j=1,3
          x = x  + SWE_PATCH_geometry%coords(:,j,cell_id) 
       end do
       x = x / 3
       ! transform x to world coordinates
       x = samoa_barycentric_to_world_point(element%transform_data, x)
       bathymetry(i) = get_bathymetry_at_position(section, x, t)
#endif !_ADAPT
       col = col + 1
       if (col == 2*row) then
          col = 1
          row = row + 1
       end if
    end do
  end function get_bathymetry_at_patch

  function get_bathymetry_at_dg_patch(section, element, t) result(bathymetry)
    type(t_grid_section), intent(inout)            :: section
    type(t_element_base), intent(inout)            :: element
    real (kind = GRID_SR), intent(in)              :: t
    real (kind = GRID_SR), dimension(_SWE_DG_DOFS) :: bathymetry
    integer (kind = GRID_SI)                       :: i,index
    real (kind = GRID_SR)                          :: x(2)


    !------- TODO: refactor me -------!    
# if (_SWE_DG_ORDER == 1)
    integer (kind = GRID_SI), parameter, dimension(_SWE_DG_DOFS) :: mirrored_coords = [1, 3, 2]
# elif (_SWE_DG_ORDER == 2)
    integer (kind = GRID_SI), parameter, dimension(_SWE_DG_DOFS) :: mirrored_coords = [1, 4, 6, 2, 5, 3]
# elif (_SWE_DG_ORDER == 3)
    integer (kind = GRID_SI), parameter, dimension(_SWE_DG_DOFS) :: mirrored_coords = [1, 5, 8, 10, 2, 6, 9, 3, 7, 4]
# elif (_SWE_DG_ORDER == 4)
    integer (kind = GRID_SI), parameter, dimension(_SWE_DG_DOFS) :: mirrored_coords = [ 1, 6, 10, 13, 15, 2, 7, 11, 14, 3, 8, 12, 4, 9, 5]
#endif
    real (kind = GRID_SR), parameter, dimension(2, _SWE_DG_DOFS) :: coords = nodes
    !------------------------- -------!    

    !iterate through cells in patch
    
    do i=1, _SWE_DG_DOFS
       if (element%cell%geometry%i_plotter_type > 0) then 
          index = i
       else
          index = mirrored_coords(i)
       end if
       x = samoa_barycentric_to_world_point(element%transform_data, coords(:,index))
       bathymetry(i) = get_bathymetry_at_position(section, x, t)
   end do

  end function get_bathymetry_at_dg_patch

  function get_bathymetry_at_position(section, x, t) result(bathymetry)
    type(t_grid_section), intent(inout)              :: section
    real (kind = GRID_SR), dimension(:), intent(in)  :: x
    real (kind = GRID_SR), intent(in)                :: t	
    real (kind = GRID_SR)                            :: bathymetry

    real (kind = c_double)            :: xs(3)
    real (kind=GRID_SR)               :: inner_height=1.20_GRID_SR, outer_height=1.10_GRID_SR
    real (kind = GRID_SR), parameter  :: dam_radius=0.1
    real (kind=GRID_SR),Dimension(2)  :: dam_center=[0.5,0.5]

    xs(1:2) = real(cfg%scaling * x + cfg%offset, c_double)
#if !defined(_ASAGI)
    bathymetry = SWE_Scenario_get_bathymetry(real(xs, GRID_SR))
#else !_ASAGI
    
#if defined(_ASAGI_TIMING)
    call section%stats%start_time(asagi_time)
#endif !_ASAGI_TIMING

    if (asagi_grid_min(cfg%afh_bathymetry, 0) <= xs(1) .and. asagi_grid_min(cfg%afh_bathymetry, 1) <= xs(2) &
         .and. xs(1) <= asagi_grid_max(cfg%afh_bathymetry, 0) .and. xs(2) <= asagi_grid_max(cfg%afh_bathymetry, 1)) then

       xs(3) = 0.0
       bathymetry = asagi_grid_get_float(cfg%afh_bathymetry, xs, 0)
    else
       bathymetry = -5000.0_SR
    end if

    if (asagi_grid_min(cfg%afh_displacement, 0) <= xs(1) .and. asagi_grid_min(cfg%afh_displacement, 1) <= xs(2) &
         .and. xs(1) <= asagi_grid_max(cfg%afh_displacement, 0) .and. xs(2) <= asagi_grid_max(cfg%afh_displacement, 1) &
         .and. t > cfg%t_min_eq) then

       xs(3) = real(min(t, cfg%t_max_eq), c_double)
       bathymetry = bathymetry + asagi_grid_get_float(cfg%afh_displacement, xs, 0)
    end if

#if defined(_ASAGI_TIMING)
    call section%stats%stop_time(asagi_time)
#endif !_ASAGI_TIMING
#endif !_ASAGI
  end function get_bathymetry_at_position
END MODULE SWE_Initialize_Bathymetry

MODULE SWE_Initialize_Dofs
  use Tools_noise

  use iso_c_binding
  use SFC_edge_traversal
  use SWE_Initialize_Bathymetry
  use Samoa_swe
#       if defined(_SWE_PATCH)
  use SWE_PATCH
#       endif

# if defined(_SWE_DG)           
  use SWE_DG_predictor
# endif          

  ! No ASAGI -> Artificial scenario selector 
#       if !defined(_ASAGI)
  use SWE_Scenario
#       endif 

  implicit none

  type num_traversal_data
     integer (kind = GRID_DI)			:: i_refinements_issued
  end type num_traversal_data

  type(t_gv_Q)							:: gv_Q
  type(t_lfs_flux)						:: lfs_flux

#		define _GT_NAME							t_swe_init_dofs_traversal

#		define _GT_EDGES
#		define _GT_EDGES_TEMP

#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_ELEMENT_OP					element_op

#		include "SFC_generic_traversal_ringbuffer.f90"

  subroutine pre_traversal_grid_op(traversal, grid)
    type(t_swe_init_dofs_traversal), intent(inout)		        :: traversal
    type(t_grid), intent(inout)							    :: grid

    grid%r_time = 0.0_GRID_SR
    grid%r_dt = 0.0_GRID_SR
    grid%r_dt_new = 0.0_GRID_SR


    call scatter(grid%r_time, grid%sections%elements_alloc%r_time)
  end subroutine pre_traversal_grid_op

  subroutine post_traversal_grid_op(traversal, grid)
    type(t_swe_init_dofs_traversal), intent(inout)		        :: traversal
    type(t_grid), intent(inout)							    :: grid

    call reduce(traversal%i_refinements_issued, traversal%sections%i_refinements_issued, MPI_SUM, .true.)
    call reduce(grid%r_dt_new, grid%sections%elements_alloc%r_dt_new, MPI_MIN, .true.)

    grid%r_dt = cfg%courant_number * grid%r_dt_new

#if defined(_ASAGI)    
#if defined(_SWE_DG)
    ! get new minimal an maximal bathymetry
    call reduce(grid%b_min, grid%sections%elements_alloc%b_min_new, MPI_MIN, .true.)
    call reduce(grid%b_max, grid%sections%elements_alloc%b_max_new, MPI_MAX, .true.)
    
    ! scatter to sections
    call scatter(grid%b_min, grid%sections%elements_alloc%b_min)
    call scatter(grid%b_max, grid%sections%elements_alloc%b_max)
    print*,"B_min: ", grid%b_min
    print*,"B_max: ", grid%b_max
#endif
#endif
    
  end subroutine post_traversal_grid_op

  subroutine pre_traversal_op(traversal, section)
    type(t_swe_init_dofs_traversal), intent(inout)				    :: traversal
    type(t_grid_section), intent(inout)							:: section

    !this variable will be incremented for each cell with a refinement request
    traversal%i_refinements_issued = 0
    section%r_dt_new  = huge(1.0_GRID_SR)
#if defined(_ASAGI)    
    section%b_min_new = huge(1.0_GRID_SR)
    section%b_max_new = tiny(1.0_GRID_SR)
#endif    
  end subroutine pre_traversal_op

  !******************
  !Geometry operators
  !******************

  subroutine element_op(traversal, section, element)
    type(t_swe_init_dofs_traversal), intent(inout)    :: traversal
    type(t_grid_section), intent(inout)               :: section
    type(t_element_base), intent(inout)               :: element
    type(t_state), dimension(_SWE_PATCH_ORDER_SQUARE)   :: Q
    integer :: i,j,count
    real (kind = GRID_SR)   :: x(2)
    type(t_dof_state) :: QS
    type(t_state), dimension(_SWE_DG_DOFS)   :: Q_DG

    if(element%cell%data_pers%troubled .ge.1)then
       element%cell%data_pers%Q%B = get_bathymetry_at_dg_patch(section, element, section%r_time)
    end if

    Q_DG(:)%b = element%cell%data_pers%Q%B
    Q(:)%b = element%cell%data_pers%B    

#if defined(_ASAGI)    
    section%b_min_new=min( section%b_min_new, minval(element%cell%data_pers%B))
    section%b_max_new=max( section%b_max_new, maxval(element%cell%data_pers%B))    
#endif

    call alpha_volume_op_dg(traversal, section, element, Q_DG)
    element%cell%data_pers%Q%H    = Q_DG(:)%h-Q_DG(:)%b    
    element%cell%data_pers%Q%p(1) = Q_DG(:)%p(1)
    element%cell%data_pers%Q%p(2) = Q_DG(:)%p(2)


    !------------- set initial state -------------------!
    element%cell%data_pers%troubled = DG    
    if(isWetDryInterface(element%cell%data_pers%Q%H))then
       element%cell%data_pers%troubled = WET_DRY_INTERFACE
    end if

    if(checkIfCellIsDry(element%cell%data_pers%Q%H)) then
       element%cell%data_pers%troubled = DRY
    end if

    if(any(abs(element%cell%data_pers%Q%B) < cfg%coast_height)) then
       element%cell%data_pers%troubled = COAST
    end if

    !----------------------------------------------------!

    !----------- set initial dofs for FV cells ----------!
    if(isFV(element%cell%data_pers%troubled))then
       if(isCoast(element%cell%data_pers%troubled)) then
          element%cell%data_pers%B = get_bathymetry_at_patch(section, element, section%r_time)
       else
          call apply_phi(element%cell%data_pers%Q%B,element%cell%data_pers%B)
       end if
       call alpha_volume_op(traversal, section, element, Q)
       element%cell%data_pers%H    = Q(:)%h
       element%cell%data_pers%HU   = Q(:)%p(1)
       element%cell%data_pers%HV   = Q(:)%p(2)

    end if
    
  end subroutine element_op


  !Volume and DoF operators
  !*******************************

  subroutine alpha_volume_op_dg(traversal, section, element, Q)
    type(t_swe_init_dofs_traversal), intent(inout)				    :: traversal
    type(t_grid_section), intent(inout)							:: section
    type(t_element_base), intent(inout)						    :: element
    type(t_state), dimension(_SWE_DG_DOFS), intent(out)  :: Q
    real (kind = GRID_SR), dimension(_SWE_DG_DOFS)       :: max_wave_speed
# if (_SWE_DG_ORDER == 1)
    integer (kind = GRID_SI), parameter, dimension(_SWE_DG_DOFS) :: mirrored_coords = [1, 3, 2]
# elif (_SWE_DG_ORDER == 2)
    integer (kind = GRID_SI), parameter, dimension(_SWE_DG_DOFS) :: mirrored_coords = [1, 4, 6, 2, 5, 3]
# elif (_SWE_DG_ORDER == 3)
    integer (kind = GRID_SI), parameter, dimension(_SWE_DG_DOFS) :: mirrored_coords = [1, 5, 8, 10, 2, 6, 9, 3, 7, 4]
# elif (_SWE_DG_ORDER == 4)
    integer (kind = GRID_SI), parameter, dimension(_SWE_DG_DOFS) :: mirrored_coords = [ 1, 6, 10, 13, 15, 2, 7, 11, 14, 3, 8, 12, 4, 9, 5]
#endif
    real (kind = GRID_SR), parameter, dimension(2, _SWE_DG_DOFS)	:: coords = nodes
    real (kind = GRID_SR), dimension(2)	:: x
    integer                             :: i, index
    integer                             :: bathy_depth
    real (kind = GRID_SR)               :: dQ_norm

    
    !evaluate initial DoFs (not the bathymetry!)
    element%cell%geometry%refinement = 0

    do i=1,_SWE_DG_DOFS
       if (element%cell%geometry%i_plotter_type > 0) then 
          index = i
       else
          index = mirrored_coords(i)
       end if
       x = samoa_barycentric_to_world_point(element%transform_data, coords(:,index))
       Q(i)%t_dof_state = get_initial_dof_state_at_position(section,x)
    end do

    if (element%cell%geometry%i_depth < cfg%i_min_depth) then
       !refine if the minimum depth is not met
       element%cell%geometry%refinement = 1
       traversal%i_refinements_issued = traversal%i_refinements_issued + 1
    end if

#if defined (_ASAGI)
    if (element%cell%geometry%i_depth < cfg%i_max_depth) then

       dQ_norm = maxval(abs(get_bathymetry_at_dg_patch(section, element, real(cfg%t_max_eq + 1.0, GRID_SR) ) - element%cell%data_pers%Q%B))

       if (dQ_norm > 2.0_SR) then
          element%cell%data_pers%troubled = WET_DRY_INTERFACE
       end if
    end if
    
#endif !_ASAGI
    

    !estimate initial time step
    where (Q%h - Q%b > 0.0_GRID_SR)
       max_wave_speed = sqrt(g * (Q%h - Q%b)) + sqrt((Q%p(1) * Q%p(1) + Q%p(2) * Q%p(2)) / ((Q%h - Q%b) * (Q%h - Q%b)))
    elsewhere
       max_wave_speed = 0.0_GRID_SR
    end where

    !This will cause a division by zero if the wave speeds are 0.
    !Bue to the min operator, the error will not affect the time step.

    if (maxval(max_wave_speed) > 0.0_GRID_SR) then
       section%r_dt_new = min(section%r_dt_new, element%cell%geometry%get_volume() / (sum(element%cell%geometry%get_edge_sizes()) * maxval(max_wave_speed) * _SWE_PATCH_ORDER))
    end if
    
  end subroutine alpha_volume_op_dg


  
  !*******************************
  !Volume and DoF operators
  !*******************************

  subroutine alpha_volume_op(traversal, section, element, Q)
    type(t_swe_init_dofs_traversal), intent(inout)				    :: traversal
    type(t_grid_section), intent(inout)							:: section
    type(t_element_base), intent(inout)						    :: element
#           if defined(_SWE_PATCH)
    type(t_state), dimension(_SWE_PATCH_ORDER_SQUARE), intent(out)  :: Q
    real (kind = GRID_SR), dimension(_SWE_PATCH_ORDER_SQUARE)       :: max_wave_speed
    real (kind = GRID_SR), DIMENSION(2)                             :: r_coords     !< cell coords within patch
    integer (kind = GRID_SI)                                        :: j, row, col, cell_id
    real (kind = GRID_SR)       :: b_min,b_max                
#           else
    type(t_state), dimension(_SWE_CELL_SIZE), intent(out)   :: Q
    real (kind = GRID_SR)               :: max_wave_speed(_SWE_CELL_SIZE)
#           endif

    real (kind = GRID_SR)               :: x(2)
    real (kind = GRID_SR)               :: dQ_norm
    real (kind = GRID_SR), parameter    :: probes(2, 3) = reshape([1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [2, 3])
    type(t_dof_state)	                :: Q_test(3)
    integer                             :: i,bathy_depth



    !evaluate initial DoFs (not the bathymetry!)

#           if defined(_SWE_PATCH)
    row = 1
    col = 1
    do i=1, _SWE_PATCH_ORDER_SQUARE

       ! if orientation is backwards, the plotter uses a transformation that mirrors the cell...
       ! this simple change solves the problem :)
       if (element%cell%geometry%i_plotter_type > 0) then 
          cell_id = i
       else
          cell_id = (row-1)*(row-1) + 2 * row - col
       end if

       r_coords = [0_GRID_SR, 0_GRID_SR]
       do j=1,3
          r_coords(:) = r_coords(:) + SWE_PATCH_geometry%coords(:,j,cell_id) 
       end do
       r_coords = r_coords / 3
       Q(i)%t_dof_state = get_initial_dof_state_at_position(section, samoa_barycentric_to_world_point(element%transform_data, r_coords))

       col = col + 1
       if (col == 2*row) then
          col = 1
          row = row + 1
       end if
    end do

#           else			
    Q%t_dof_state = get_initial_dof_state_at_element(section, element)

    ! dry cells
    if (Q(1)%h < Q(1)%b + cfg%dry_tolerance) then
       Q%h  = Q%b
       Q(1)%p(:) = [0.0_GRID_SR, 0.0_GRID_SR]
    end if
#           endif

#if defined(_SWE_DG)
    traversal%i_refinements_issued = traversal%i_refinements_issued - element%cell%geometry%refinement
#endif

    element%cell%geometry%refinement = 0
!    print*,traversal%i_refinements_issued

    if (element%cell%geometry%i_depth < cfg%i_min_depth) then
       !refine if the minimum depth is not met
       element%cell%geometry%refinement = 1
       traversal%i_refinements_issued = traversal%i_refinements_issued + 1
    else if (element%cell%geometry%i_depth < cfg%i_max_depth) then
#               if defined (_ASAGI)
       !refine the displacements

       x = samoa_barycentric_to_world_point(element%transform_data, [1.0_SR/3.0_SR, 1.0_SR/3.0_SR])


# if defined(_SWE_DG)


# if defined (_SWE_PATCH)
       dQ_norm = maxval(abs(get_bathymetry_at_patch(section, element, real(cfg%t_max_eq + 1.0, GRID_SR) ) - element%cell%data_pers%B))
#                   else
       dQ_norm = abs(get_bathymetry_at_element(section, element, real(cfg%t_max_eq + 1.0, GRID_SR) ) - Q(1)%b)
#                   endif
       
       if (dQ_norm > 2.0_SR) then
          element%cell%geometry%refinement = 1
          traversal%i_refinements_issued = traversal%i_refinements_issued + 1
       end if
#endif
#               else
       !refine any slopes in the initial state
       do i = 1, 3
          x = samoa_barycentric_to_world_point(element%transform_data, probes(:, i))
          Q_test(i) = get_initial_dof_state_at_position(section, x)
       end do

       if (maxval(Q_test%h) > minval(Q_test%h)) then
          element%cell%geometry%refinement = 1
          traversal%i_refinements_issued = traversal%i_refinements_issued + 1
       end if
#               endif

    end if

    !estimate initial time step

    where (Q%h - Q%b > 0.0_GRID_SR)
       max_wave_speed = sqrt(g * (Q%h - Q%b)) + sqrt((Q%p(1) * Q%p(1) + Q%p(2) * Q%p(2)) / ((Q%h - Q%b) * (Q%h - Q%b)))
    elsewhere
       max_wave_speed = 0.0_GRID_SR
    end where

    !This will cause a division by zero if the wave speeds are 0.
    !Bue to the min operator, the error will not affect the time step.
#           if defined _SWE_PATCH
    if(maxval(max_wave_speed) > 0.0_GRID_SR)then
    section%r_dt_new = min(section%r_dt_new, element%cell%geometry%get_volume() / (sum(element%cell%geometry%get_edge_sizes()) * maxval(max_wave_speed) * _SWE_PATCH_ORDER))
#           else
    section%r_dt_new = min(section%r_dt_new, element%cell%geometry%get_volume() / (sum(element%cell%geometry%get_edge_sizes()) * maxval(max_wave_speed)))
#           endif

    end if

  end subroutine alpha_volume_op

  function get_initial_dof_state_at_element(section, element) result(Q)
    type(t_grid_section), intent(inout)		:: section
    type(t_element_base), intent(inout)     :: element
    type(t_dof_state)						:: Q
    real (kind = GRID_SR)		            :: x(2)

    x = samoa_barycentric_to_world_point(element%transform_data, [1.0_SR / 3.0_SR, 1.0_SR / 3.0_SR])
    Q = get_initial_dof_state_at_position(section, x)
  end function get_initial_dof_state_at_element


  function get_initial_dof_state_at_position(section, x) result(Q)
    type(t_grid_section), intent(inout)					:: section
    real (kind = GRID_SR), intent(in)		            :: x(:)            !< position in world coordinates
    type(t_dof_state)							        :: Q
    real (kind = GRID_SR), parameter		            :: hL = 1.20_GRID_SR, hR = 1.30_GRID_SR
    real (kind = GRID_SR), parameter		            :: dam_radius=0.2
    real(kind=GRID_SR),Dimension(2)                             :: dam_center=[0.5,0.5]
    real (kind = GRID_SR)                                       :: xs(2)

    xs = cfg%scaling * x + cfg%offset


#			if defined(_ASAGI)
    Q%h = 0.0_GRID_SR                                    
    Q%p = 0.0_GRID_SR            
#			else

    Q = SWE_Scenario_get_initial_Q(xs)
#			endif

  end function get_initial_dof_state_at_position


END MODULE SWE_Initialize_Dofs
#endif

