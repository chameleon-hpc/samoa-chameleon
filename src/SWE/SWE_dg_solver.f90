! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2017 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE
! Author: Leonhard Rannabauer rannabau (at) in.tum.de

#if defined(_SWE_DG)
#include "Compilation_control.f90"

MODULE SWE_DG_solver
  use SFC_edge_traversal
  !use SWE_euler_timestep
  use Samoa_swe
  use SWE_PDE
  use SWE_FV_solver

#if defined(CHAMELEON)
   use Chameleon_lib
#define _GT_USE_CHAMELEON
#endif
#if defined(CHAMELEON_CALL)
#  define _GT_USE_CHAMELEON_CALL
#endif
   
  use SWE_DG_Matrices
  use SWE_DG_Limiter  
  use SWE_initialize_bathymetry

#if defined(_HLLE_FLUX)
   use SWE_HLLE
#endif
#if defined(_BOUNDARY_FUNC)
   use SWE_Scenario
#endif
#if defined(_BOUNDARY_FILE)
   use Tools_boundary_file
#endif


  implicit none
  type num_traversal_data
     integer (kind = GRID_DI)			:: i_refinements_issued=0
  end type num_traversal_data

  interface skeleton_op_dg
     module procedure skeleton_array_op_dg
     module procedure skeleton_scalar_op_dg
  end interface skeleton_op_dg

  interface bnd_skeleton_op_dg
     module procedure bnd_skeleton_array_op_dg
     module procedure bnd_skeleton_scalar_op_dg
  end interface bnd_skeleton_op_dg

#		define _GT_NAME	           		t_swe_dg_timestep_traversal

#		define _GT_EDGES
!#   define _GT_EDGE_MPI_TYPE

#		define _GT_PRE_TRAVERSAL_OP      	pre_traversal_op_dg
#		define _GT_PRE_TRAVERSAL_GRID_OP	pre_traversal_grid_op_dg
#		define _GT_POST_TRAVERSAL_GRID_OP	post_traversal_grid_op_dg

#		define _GT_SKELETON_OP		       	skeleton_op_dg
#		define _GT_BND_SKELETON_OP	        bnd_skeleton_op_dg
#		define _GT_ELEMENT_OP                   element_op_dg

#		define _GT_CELL_UPDATE_OP		cell_update_op_dg
#		define _GT_CELL_TO_EDGE_OP	cell_to_edge_op_dg

#if defined(_SWE_DG_LIMITER_ALL) || defined(_SWE_DG_LIMITER_HEIGHT)
#   define _GT_NODES
#   define _GT_ELEMENT_OP               element_op
#   define _GT_NODE_FIRST_TOUCH_OP      node_first_touch_op
#   define _GT_NODE_MERGE_OP            node_merge_op
#endif

  public cell_to_edge_op_dg

#		include "SFC_generic_traversal_ringbuffer.f90"

  subroutine post_traversal_grid_op_dg(traversal, grid)
    type(t_swe_dg_timestep_traversal), intent(inout)		:: traversal
    type(t_grid), intent(inout)					:: grid

    grid%r_time = grid%r_time + grid%r_dt
    
    call reduce(traversal%i_refinements_issued, traversal%sections%i_refinements_issued, MPI_SUM, .true.)
    call reduce(grid%r_dt_new, grid%sections%elements_alloc%r_dt_new, MPI_MIN, .true.)

    grid%min_courant = min(grid%r_dt_new / grid%r_dt , grid%min_courant)
    
    grid%r_dt_new = cfg%courant_number * grid%r_dt_new

    if (rank_MPI == 0) then
       if (cfg%courant_number > grid%r_dt_new / grid%r_dt) then
          _log_write(1, '("WARNING! Time step size was too big. dt (old): ", ES10.3, ", dt (CFL): ", ES10.3, ", maximum courant number: ", F0.3)') grid%r_dt, grid%r_dt_new / cfg%courant_number, grid%r_dt_new / grid%r_dt
       end if
    end if
    
    grid%r_dt = grid%r_dt_new
    grid%r_dt = min(cfg%r_max_time-grid%r_time, grid%r_dt)    
    call scatter(grid%r_time, grid%sections%elements_alloc%r_time)

    ! reduce and scatter refinement errors
    call reduce(grid%min_error, grid%sections%elements_alloc%min_error_new, MPI_MIN, .true.)
    call reduce(grid%max_error, grid%sections%elements_alloc%max_error_new, MPI_MAX, .true.)
    
    call scatter(grid%min_error, grid%sections%elements_alloc%min_error)
    call scatter(grid%max_error, grid%sections%elements_alloc%max_error)
    
  end subroutine post_traversal_grid_op_dg

  subroutine pre_traversal_grid_op_dg(traversal, grid)
    type(t_swe_dg_timestep_traversal), intent(inout)		:: traversal
    type(t_grid), intent(inout)							    :: grid
  end subroutine pre_traversal_grid_op_dg


  subroutine pre_traversal_op_dg(traversal, section)
    type(t_swe_dg_timestep_traversal), intent(inout)				:: traversal
    type(t_grid_section), intent(inout)							:: section
    !this variable will be incremented for each cell with a refinement request
    traversal%i_refinements_issued = 0_GRID_DI
    section%r_dt_new      =  huge(1.0_GRID_SR)
    section%min_error_new =  huge(1.0_GRID_SR)
    section%max_error_new = -huge(1.0_GRID_SR)

  end subroutine pre_traversal_op_dg

  subroutine element_op_dg(traversal, section, element)
    type(t_swe_dg_timestep_traversal), intent(inout)				:: traversal
    type(t_grid_section), intent(inout)							:: section
    type(t_element_base), intent(inout)						:: element
  end subroutine element_op_dg

  subroutine writeFVBoundaryFields(data,H,HU,HV,B,edge,is_coast)    
    class(num_cell_data_pers), intent(in) :: data
    real(KIND=GRID_SR), intent(out)       :: H ( _SWE_PATCH_ORDER)
    real(KIND=GRID_SR), intent(out)       :: HU( _SWE_PATCH_ORDER)
    real(KIND=GRID_SR), intent(out)       :: HV( _SWE_PATCH_ORDER)
    real(KIND=GRID_SR), intent(out)       :: B ( _SWE_PATCH_ORDER)    
    integer, intent(in)                   :: edge
    logical, intent(in)                   :: is_coast
    integer                               :: i

    ! if(.not.is_coast) then
    !    if(edge == 1) then
    !       H  = matmul(phi_r(:,:),data%Q%h)
    !       HU = matmul(phi_r(:,:),data%Q%p(1))
    !       HV = matmul(phi_r(:,:),data%Q%p(2))
    !       B  = matmul(phi_r(:,:),data%Q%b)
    !    end if
    !    if(edge == 2) then
    !       H  = matmul(phi_m(:,:),data%Q%h)
    !       HU = matmul(phi_m(:,:),data%Q%p(1))
    !       HV = matmul(phi_m(:,:),data%Q%p(2))
    !       B  = matmul(phi_m(:,:),data%Q%b)
    !    end if
    !    if(edge == 3) then
    !       H  = matmul(phi_l(:,:),data%Q%h)
    !       HU = matmul(phi_l(:,:),data%Q%p(1))
    !       HV = matmul(phi_l(:,:),data%Q%p(2))
    !       B  = matmul(phi_l(:,:),data%Q%b)
    !    end if
    !      H  = H + B
    ! else
    do i=1,_SWE_PATCH_ORDER         
       H (i) = data%H (idx_project_FV(edge,i))
       HU(i) = data%HU(idx_project_FV(edge,i))
       HV(i) = data%HV(idx_project_FV(edge,i))
       B (i) = data%B (idx_project_FV(edge,i))
    end do
!    end if
  end subroutine writeFVBoundaryFields

  function cell_to_edge_op_dg(element, edge) result(rep)
    type(t_element_base), intent(in)			    :: element
    type(t_edge_data), intent(in)			    :: edge
    
    type(num_cell_rep)				            :: rep
    integer                               :: edge_type !1=right, 2=hypotenuse, 3=left
    
    integer :: i,j,k,indx,indx_mir
    integer(kind=BYTE) :: orientation
    type(num_cell_rep) ::rep_fv

#if defined (_PREFETCH) 
    call mm_prefetch(element%cell%data_pers%Q_DG_UPDATE,0)
#endif

#if defined (_DEBUG)
    rep%debug_flag = element%cell%data_pers%debug_flag
#endif
    edge_type = get_edge_type(abs(element%cell%geometry%i_plotter_type),edge%transform_data%normal)
    
    rep%QP(:,:)   = element%cell%data_pers%QP(edge_type  ,:,:)
    rep%FP(:,:,:) = element%cell%data_pers%FP(edge_type,:,:,:)

    call writeFVBoundaryFields(element%cell%data_pers,rep%H,rep%HU,rep%HV,rep%B,edge_type,&
    isCoast(element%cell%data_pers%troubled) .or.&
    element%cell%data_pers%troubled .eq. NEIGHBOUR_WAS_TROUBLED .or.&
    element%cell%data_pers%troubled .eq. WET_DRY_INTERFACE)
    !call writeFVBoundaryFields(element%cell%data_pers,rep%H,rep%HU,rep%HV,rep%B,edge_type,isCoast(element%cell%data_pers%troubled))

    rep%troubled=element%cell%data_pers%troubled
    
  end function cell_to_edge_op_dg


function get_edge_type(cell_type,norm)
  real(KIND=GRID_SR) :: norm(2)
  integer :: cell_type
  integer :: get_edge_type

  get_edge_type = norm_to_edge_type(cell_type,nint(norm(1) * 2.0_GRID_SR),nint(norm(2) * 2.0_GRID_SR))
end function
  
  
subroutine skeleton_array_op_dg(traversal, grid, edges, rep1, rep2, update1, update2)
type(t_swe_dg_timestep_traversal), intent(in)	            :: traversal
type(t_grid_section), intent(in)		            	    :: grid
type(t_edge_data), intent(in)				    :: edges(:)
type(num_cell_rep), intent(in)				    :: rep1(:), rep2(:)
type(num_cell_update), intent(out)				    :: update1(:), update2(:)
integer (kind = GRID_SI)                                       :: i

do i = 1, size(edges)
   call skeleton_scalar_op_dg(traversal, grid, edges(i), rep1(i), rep2(i), update1(i), update2(i))
end do

end subroutine skeleton_array_op_dg

subroutine general_dg_riemannsolver(edge,rep1,rep2,update1,update2)
  type(t_edge_data), intent(in)	                          :: edge
  type(num_cell_rep), intent(in)                          :: rep1   , rep2
  type(num_cell_update), intent(out)                      :: update1, update2
  real(kind=GRID_SR),dimension((_SWE_DG_ORDER+1)**2,4)    :: QLp ,QRp
  real(kind=GRID_SR),dimension(_SWE_DG_ORDER+1,4)    :: QL_s ,QR_s
  real(kind=GRID_SR),dimension((_SWE_DG_ORDER+1),4)       :: QL ,QR
  real(kind=GRID_SR),dimension((_SWE_DG_ORDER+1),3)    :: FLn_s,FRn_s
  real(kind=GRID_SR),dimension(2,(_SWE_DG_ORDER+1)**2,3)    :: FLg,FRg
  real(kind=GRID_SR),dimension((_SWE_DG_ORDER+1)**2,3)    :: FLng,FRng
  real(kind=GRID_SR),dimension((_SWE_DG_ORDER+1),3)    :: FLn,FRn
  real(kind=GRID_SR),dimension((_SWE_DG_ORDER+1),3)    :: FRn_temp
  real(kind=GRID_SR),dimension(2,(_SWE_DG_ORDER+1),3)  :: FL ,FR
  
  real(kind=GRID_SR)        nF1(_SWE_DG_ORDER+1,_SWE_DG_ORDER+1,4), nF2(_SWE_DG_ORDER+1,_SWE_DG_ORDER+1,4)
  real(kind=GRID_SR)        nF3(_SWE_DG_ORDER+1,_SWE_DG_ORDER+1,3), nF4(_SWE_DG_ORDER+1,_SWE_DG_ORDER+1,3)

  type(t_update), dimension((_SWE_DG_ORDER+1)**2)         :: flux_temp
  real(kind = GRID_SR)                   	          :: normal(2)
  integer                                                 :: i,j
! real(kind = GRID_SR),Parameter :: error = 1.0d-10

  normal = edge%transform_data%normal/NORM2(edge%transform_data%normal)
   !        ______
   ! |6 \   \3 2 1|
   ! |4 5 \   \5 4|
   ! |1_2_3\   \ 6|
   
  QL = rep1%QP
  FL = rep1%FP
  
  QR = rep2%QP
  FR = rep2%FP

  FLn = FL(1,:,:) * normal(1) + FL(2,:,:) * normal(2)
  FRn = FR(1,:,:) * normal(1) + FR(2,:,:) * normal(2)
  
  call compute_flux_pred(normal,QL,QR,FLn,FRn)

  update1%flux%h    = FLn(:,1)
  update1%flux%p(1) = FLn(:,2)
  update1%flux%p(2) = FLn(:,3)
  
  update2%flux%h    = FRn(:,1)
  update2%flux%p(1) = FRn(:,2)
  update2%flux%p(2) = FRn(:,3)
  
end subroutine

subroutine skeleton_scalar_op_dg(traversal, grid, edge, rep1, rep2, update1, update2)
type(t_swe_dg_timestep_traversal), intent(in)		:: traversal
type(t_grid_section), intent(in)			:: grid
type(t_edge_data), intent(in)				:: edge
type(num_cell_rep), intent(in)				:: rep1, rep2
type(num_cell_update), intent(out)			:: update1, update2
integer                                                 :: i,j, edge_type
real(kind = GRID_SR)	          :: normal(2)
type(num_cell_rep)				:: rep2_rev
type(num_cell_update)			:: update2_rev



update1%troubled      =rep2%troubled
update2%troubled      =rep1%troubled

if(isDG(rep1%troubled) .and. isDG(rep2%troubled)) then
   rep2_rev%troubled = rep2%troubled
   do i = 1,_SWE_DG_ORDER+1
      rep2_rev%QP(  i,:) = rep2%QP(  _SWE_DG_ORDER+2-i,:)
      rep2_rev%FP(1,i,:) = rep2%FP(1,_SWE_DG_ORDER+2-i,:)
      rep2_rev%FP(2,i,:) = rep2%FP(2,_SWE_DG_ORDER+2-i,:)
   end do

   call general_dg_riemannsolver(edge,rep1,rep2_rev,update1,update2_rev)
   
   do i = 1,_SWE_DG_ORDER+1
      update2%flux(i) = update2_rev%flux(_SWE_DG_ORDER+2-i)
   end do
   
end if


update2%H  = rep1%H
update2%HU = rep1%HU
update2%HV = rep1%HV
update2%B  = rep1%B

update1%H  = rep2%H
update1%HU = rep2%HU
update1%HV = rep2%HV
update1%B  = rep2%B

! update1%minObservables = rep2%minObservables
! update1%maxObservables = rep2%maxObservables
! update2%minObservables = rep1%minObservables
! update2%maxObservables = rep1%maxObservables
end subroutine skeleton_scalar_op_dg

subroutine bnd_skeleton_array_op_dg(traversal, grid, edges, rep, update)
type(t_swe_dg_timestep_traversal), intent(in)				:: traversal
type(t_grid_section), intent(in)						:: grid
type(t_edge_data), intent(in)						:: edges(:)
type(num_cell_rep), intent(in)						:: rep(:)
type(num_cell_update), intent(out)						:: update(:)
integer :: i
do i = 1, size(edges)
   call bnd_skeleton_scalar_op_dg(traversal, grid, edges(i), rep(i), update(i))
end do
end subroutine bnd_skeleton_array_op_dg

subroutine bnd_skeleton_scalar_op_dg(traversal, grid, edge, rep, update)
type(t_swe_dg_timestep_traversal), intent(in)       :: traversal
type(t_grid_section), intent(in)		    :: grid
type(t_edge_data), intent(in)		            :: edge
type(num_cell_rep), intent(in)		            :: rep
type(num_cell_rep)        		            :: rep_bnd
type(num_cell_update), intent(out)		    :: update
type(num_cell_update)	                	    :: update_bnd
real(kind=GRID_SR)                                  :: normal(2)
real(kind=GRID_SR)                                  :: length_flux
integer                                             :: i,j
real(kind=GRID_SR), dimension(_DMP_NUM_OBSERVABLES)             :: observables
real(kind=GRID_SR)                                  :: h_psi,hu_psi,hv_psi

normal=(edge%transform_data%normal)/NORM2(edge%transform_data%normal)
update%troubled=rep%troubled
rep_bnd%troubled = rep%troubled
if(isDG(rep%troubled)) then
   update_bnd = update
   rep_bnd%QP(:,1) = rep%QP(:,1)
   rep_bnd%QP(:,4) = rep%QP(:,4)
#  if defined(_BOUNDARY)
   if(get_edge_boundary_align(normal)) then
      rep_bnd%FP(:,:,:) = 0.0_GRID_SR
      rep_bnd%QP(:,1:3) = 0.0_GRID_SR
      rep_bnd%QP(:,4) = rep%QP(:,4)
      h_psi  = 0.0_GRID_SR
      hu_psi = 0.0_GRID_SR
      hv_psi = 0.0_GRID_SR   
         ! Time-dependent condition, generate virtual wave from data or function
         do i=1, _SWE_DG_ORDER+1
            do j = 1, _SWE_DG_ORDER+1
               call get_edge_boundary_Q(&
               grid%r_time + (gll_nodes(j) * grid%r_dt),&
               h_psi, &
               hu_psi, &
               hv_psi, &
               rep_bnd%QP(i,4), .true.)
               
               rep_bnd%QP(i,1) = rep_bnd%QP(i,1) + h_psi *gll_weights(j)
               rep_bnd%QP(i,2) = rep_bnd%QP(i,2) + hu_psi*gll_weights(j)
               rep_bnd%QP(i,3) = rep_bnd%QP(i,3) + hv_psi*gll_weights(j)
               
               rep_bnd%FP(1,i,:)=rep_bnd%FP(1,i,:)+flux_1_b((/ h_psi,hu_psi,hv_psi /))*gll_weights(j)
               rep_bnd%FP(2,i,:)=rep_bnd%FP(2,i,:)+flux_2_b((/ h_psi,hu_psi,hv_psi /))*gll_weights(j)

            end do
         end do
      else
#  endif
         ! Generate mirrored wave to reflect out incoming wave
         do i=1,(_SWE_DG_ORDER+1)
            length_flux = dot_product(rep%QP(i,2:3), normal)
            rep_bnd%QP(i,2) = rep%QP(i,2)-2.0_GRID_SR*length_flux*normal(1)
            rep_bnd%QP(i,3) = rep%QP(i,3)-2.0_GRID_SR*length_flux*normal(2)
         end do

         rep_bnd%FP(1,:,1) = -rep%FP(1,:,1) 
         rep_bnd%FP(1,:,2) =  rep%FP(1,:,2) 
         rep_bnd%FP(1,:,3) = -rep%FP(1,:,3) 

         rep_bnd%FP(2,:,1) = -rep%FP(2,:,1) 
         rep_bnd%FP(2,:,2) = -rep%FP(2,:,2) 
         rep_bnd%FP(2,:,3) =  rep%FP(2,:,3) 

#  if defined(_BOUNDARY)
      end if
#  endif
   call general_dg_riemannsolver(edge,rep,rep_bnd,update,update_bnd)  
end if
#  if defined(_BOUNDARY)
if(get_edge_boundary_align(normal)) then
   ! Time-dependent condition, generate virtual wave from data or function
   do i=1, _SWE_PATCH_ORDER
      update%H(i) = rep%H(i)
      update%B(i) = rep%B(i)
      call get_edge_boundary_Q(grid%r_time, update%H(i), &
           update%HU(i), update%HV(i), update%B(i), .false.)
   end do
! else if ((normal(1) .eq. 1.0_GRID_SR) .and. (normal(2) .eq. 0.0_GRID_SR))then
!    if(isDG(update%troubled))then
!       update%troubled=1
!    end if
!    update%H=rep%H
!    update%B=rep%B
!    update%HU=rep%HU
!    update%HV=rep%HV
else
#  endif
   ! Generate mirrored wave to reflect out incoming wave
   do i=1, _SWE_PATCH_ORDER    
      length_flux = dot_product([ rep%HU(i), rep%HV(i) ], normal)
      update%H(_SWE_PATCH_ORDER + 1 - i)=rep%H(i)
      update%B(_SWE_PATCH_ORDER + 1 - i)=rep%B(i)
      update%HU(_SWE_PATCH_ORDER + 1 - i)=rep%HU(i) - 2.0_GRID_SR*length_flux*normal(1)
      update%HV(_SWE_PATCH_ORDER + 1 - i)=rep%HV(i) - 2.0_GRID_SR*length_flux*normal(2)
   end do
   !compute new extrema
   do i = 1,_SWE_PATCH_ORDER
      call getObservables(update%H(i),update%HU(i),update%HV(i),update%B(i),observables)
      do j = 1,_DMP_NUM_OBSERVABLES
         update%minObservables(j) = min(update%minObservables(j),observables(j))
         update%maxObservables(j) = max(update%maxObservables(j),observables(j))
      end do
   end do

#  if defined(_BOUNDARY)
end if
#  endif

end subroutine bnd_skeleton_scalar_op_dg

#if defined(_BOUNDARY)
   function get_edge_boundary_align(normal) result(align)
      real(GRID_SR), dimension(2), intent(in)   :: normal
      logical                                   :: align

      select case(cfg%i_boundary_side)
         case(0)
            align = (normal(1) .eq. 0.0_GRID_SR) .and. (normal(2) .eq. 1.0_GRID_SR)
         case(1)
            align = (normal(1) .eq. 1.0_GRID_SR) .and. (normal(2) .eq. 0.0_GRID_SR)
         case(2)
            align = (normal(1) .eq. 0.0_GRID_SR) .and. (normal(2) .eq. -1.0_GRID_SR)
         case(3)
            align = (normal(1) .eq. -1.0_GRID_SR) .and. (normal(2) .eq. 0.0_GRID_SR)
      end select
   end function get_edge_boundary_align

   subroutine get_edge_boundary_Q(t, H, HU, HV, B, isDG)
      real(GRID_SR), intent(in)                 :: t
      real(GRID_SR), intent(out)              :: H, HU, HV
      real(GRID_SR), intent(in)                 :: B
      logical, intent(in)                       :: isDG

      real(GRID_SR)                             :: velocity, h_init

#     if defined(_BOUNDARY_FUNC)
         H = max(0.0_GRID_SR, SWE_Scenario_get_boundary_height(B, t))
         h_init = max(0.0_GRID_SR, SWE_Scenario_get_boundary_height(B, 0.0_GRID_SR))
#     elif defined(_BOUNDARY_FILE)
         H = max(0.0_GRID_SR, boundary_file_get(t) - B)
         h_init = max(0.0_GRID_SR, boundary_file_get(0.0_GRID_SR) - B)
#     endif
      velocity = 2.0_GRID_SR * (sqrt(_GRAV_CONSTANT * H) - sqrt(_GRAV_CONSTANT * h_init))
      select case(cfg%i_boundary_side)
         case(0)
            HU = 0.0_GRID_SR
            HV = -velocity
         case(1)
            HU = -velocity
            HV = 0.0_GRID_SR
         case(2)
            HU = 0.0_GRID_SR
            HV = velocity
         case(3)
            HU = velocity
            HV = 0.0_GRID_SR
      end select
      HU = HU * H
      HV = HV * H
      if(.not. isDG) then
         ! FV stores H + B, DG stores H
         H = H + B
      end if
   end subroutine
#endif

subroutine cell_update_op_dg(traversal, section, element, update1, update2, update3)
type(t_swe_dg_timestep_traversal), intent(inout)		:: traversal
type(t_grid_section), intent(inout)			:: section
type(t_element_base), intent(inout)			:: element
type(num_cell_update), intent(inout)			:: update1, update2, update3
type(num_cell_update) :: tmp
type(t_update),  DIMENSION((_SWE_DG_ORDER+1)**2)		:: flux1, flux2, flux3,flux2_temp
type(t_state), dimension((_SWE_DG_ORDER+1)**2)             :: edge_l, edge_m,edge_r
real(kind= GRID_SR) :: dt,dx,delta(3),max_neighbour(3),min_neighbour(3),data_max_val(3),data_min_val(3)
integer :: indx,indx_mir,k,i,j
real (kind=GRID_SR) :: max_wave_speed, wave_speed=0,dQ_norm,b_min,b_max,maxWaveSpeedFV
real (kind=GRID_SR) :: refinement_threshold = 10.25_GRID_SR/_SWE_DG_ORDER
real (kind = GRID_SR), dimension (_SWE_PATCH_ORDER_SQUARE) :: H_old, HU_old, HV_old, B_old
real (kind = GRID_SR), dimension (_SWE_PATCH_ORDER_SQUARE) :: H, HU, HV
logical :: refine,coarsen
integer :: i_depth,bathy_depth
real(kind=GRID_SR),Dimension(3) :: edge_sizes
real(kind=GRID_SR) :: error
real(kind=GRID_SR) :: refine_factor = 10.0d-1, coarsen_factor=10.0d-4

!> For DMP
real(kind=GRID_SR) :: minVals(_DMP_NUM_OBSERVABLES)
real(kind=GRID_SR) :: maxVals(_DMP_NUM_OBSERVABLES)

associate(data => element%cell%data_pers)

!----If a cell is dry and all neighbours are we can skip the solver step---!
 if(isDry(data%troubled)) then
    if(allNeighboursDry(update1,update2,update3))then
#if defined(_CELL_METRICS)
      data%dt(:)    = -1.0_GRID_SR
      data%dt_fv(:) = -1.0_GRID_SR
#endif      
      return         
   end if
end if
!--------------------------------------------------------------------------!


!--------------For a positive orientation update are switched--------------!
if (element%cell%geometry%i_plotter_type > 0) then ! 
   tmp=update1
   update1=update3
   update3=tmp

end if
!--------------------------------------------------------------------------!

do i = 1,_SWE_DG_ORDER+1
   tmp%flux(i) = update3%flux(_SWE_DG_ORDER+2-i)
end do
update3%flux = tmp%flux

!------- Update cell status and compute next timestep size --------!
call updateCellStatusPre(data,update1,update2,update3)

edge_sizes=element%cell%geometry%get_edge_sizes()

if(isDG(data%troubled)) then
   data%troubled = DG
   H_old  = data%H
   HU_old = data%HU
   HV_old = data%HV
   B_old  = data%B

   call getObservableLimits(element%cell%data_pers,minVals,maxVals)

   call dg_solver(element,update1%flux,update2%flux,update3%flux,section%r_dt)

   call apply_phi_cons(data%Q(:)%h,data%Q(:)%p(1),data%Q(:)%p(2),&
                       data%H(:),data%HU(:),data%HV(:))             
   call apply_phi(data%Q(:)%b    ,data%b)
   
   data%h = data%h + data%b
   
   if(isWetDryInterface(data%Q%H)) then
      data%troubled = WET_DRY_INTERFACE
   else if(checkDMP(element%cell%data_pers,minVals,maxVals,&
        element%nodes(1)%ptr%data_pers,&
        element%nodes(2)%ptr%data_pers,&
        element%nodes(3)%ptr%data_pers)) then
      data%troubled = TROUBLED
   end if
   
   if(isFV(data%troubled)) then
      !--if troubled or drying perform rollback--!
      data%H  = H_old
      data%HU = HU_old
      data%HV = HV_old
      data%B  = B_old

      ! call apply_phi_cons(data%Q(:)%h,data%Q(:)%p(1),data%Q(:)%p(2),&
      !                     data%H(:),data%HU(:),data%HV(:))             
      ! call apply_phi(data%Q(:)%b    ,data%b)
      !data%H = data%H + data%b
   end if
end if

if(isFV(data%troubled)) then
   !-----Call FV patch solver----!    
   call fv_patch_solver(element, update1, update2, update3, section%r_dt, maxWaveSpeedFV)
end if
!------ Update cell status------!
call updateCellStatus(data)
!------------------------------------!


!-------compute next timestep size --------!
dx = cfg%scaling *  edge_sizes(1)
if(.not.isDry(data%troubled)) then
   if(isDG(data%troubled))then
      do i=1,_SWE_DG_DOFS
         section%r_dt_new = min(section%r_dt_new, get_next_time_step_size(data%Q(i)%h,data%Q(i)%p(1),data%Q(i)%p(2),dx,cfg%dry_tolerance))
      end do
   else
      dt = dx  / (maxWaveSpeedFV * (_SWE_DG_ORDER*2.0_GRID_SR+2.0_GRID_SR))
      ! if(isCoast(data%troubled)) then
      !    !limit by courant condition on coast
      !    dt = max( dt , section%r_dt * cfg%courant_number)
      ! end if
      section%r_dt_new = min(section%r_dt_new, dt)
   endif
else
   section%r_dt_new = min(section%r_dt_new,1000.0_GRID_SR)
endif

#if defined(_CELL_METRICS)
if(.not.isDry(data%troubled)) then
   do i=1,_SWE_DG_DOFS
      if(data%Q(i)%h > cfg%dry_tolerance) then
         dt = get_next_time_step_size(data%Q(i)%h,data%Q(i)%p(1),data%Q(i)%p(2),dx,cfg%dry_tolerance)
      else
         dt = section%r_dt / cfg%courant_number
      end if
      data%dt(i)    = dt * cfg%courant_number
   end do
   dt = dx  / (maxWaveSpeedFV * (_SWE_DG_ORDER*2.0_GRID_SR+2.0_GRID_SR))
   data%dt_fv(:) =  dt * cfg%courant_number
else
   data%dt(:)    = 0.0_GRID_SR
   data%dt_fv(:) = 0.0_GRID_SR
endif
#endif  
!------------------------------------------------------------------!


!----------- Refinement ------------!
element%cell%geometry%refinement = 0

!refine      
!consider wave
if (.not.isCoast(element%cell%data_pers%troubled))then
   error = get_error_estimate(data%Q)      
   if(error > section%max_error * refine_factor) then
      refine = .true.
   endif
   if(error < section%max_error * coarsen_factor) then
      coarsen = .true.
   endif
   
   i_depth = element%cell%geometry%i_depth
   
   !----- Compute and update error estimate -----!
   section%min_error_new = min(error,section%min_error_new)
   section%max_error_new = max(error,section%max_error_new)
   !---------------------------------------------!
   
   if (i_depth < cfg%i_max_depth .and. refine) then
      element%cell%geometry%refinement = 1
      traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
   else if (i_depth > cfg%i_min_depth .and. coarsen) then
      element%cell%geometry%refinement = -1         
   end if
end if
!-----------------------------------!

end associate

end subroutine cell_update_op_dg


subroutine compute_flux_pred(normal,QL, QR, FLn, FRn)

real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,4), intent(in)    :: QL
real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,4), intent(in)    :: QR
real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,3), intent(inout) :: FLn
real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,3), intent(inout) :: FRn
real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1)                  :: VelL
real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1)                  :: VelR
real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1)                  :: hRoe,bm,bml,bmr
real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1)                  :: Deta
real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1)                  :: Djump
real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,3)                :: Fn_avg
real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,3)                :: Q_rus

real(kind = GRID_SR), intent(in)  :: normal(2)

real(kind = GRID_SR)	          :: epsilon
real(kind=GRID_SR)	          :: vL, vR, alpha
integer                           :: i

Fn_avg = 0.5_GRID_SR * ( FRn - FLn )

VelL(:) = QL(:,2)/QL(:,1) * normal(1) +  QL(:,3)/QL(:,1) * normal(2)
VelR(:) = QR(:,2)/QR(:,1) * normal(1) +  QR(:,3)/QR(:,1) * normal(2)

alpha=max(maxval(abs(VelL) + sqrt(g * QL(:,1))),&
     maxval(abs(VelR) + sqrt(g * QR(:,1))))

hRoe  = 0.5_GRID_SR * (QR(:,1) + QL(:,1))
bm    = max(QR(:,4),QL(:,4))
bmr   = merge(0.0_GRID_SR, QR(:,4) - bm , abs(QR(:,4) - bm) < 1.0e-16)
bml   = merge(0.0_GRID_SR, QL(:,4) - bm , abs(QL(:,4) - bm) < 1.0e-16)
Deta  = max(QR(:,1) + bmr, 0.0) - max(QL(:,1) + bml, 0.0)
Deta  = merge(Deta, 0.0_GRID_SR, abs(Deta) > 1.0e-16)

Djump = 0.5_GRID_SR * g * hRoe * Deta

Q_rus(:,  1) = 0.5_GRID_SR * alpha * Deta
Q_rus(:,2:3) = 0.5_GRID_SR * alpha * (QR(:,2:3) - QL(:,2:3))


FLn = Fn_avg - Q_rus
FRn = Fn_avg + Q_rus

FLn(:,1) = FLn(:,1) 
FLn(:,2) = FLn(:,2) + Djump * normal(1)
FLn(:,3) = FLn(:,3) + Djump * normal(2)

FRn(:,1) = FRn(:,1)    
FRn(:,2) = FRn(:,2) + Djump * normal(1)
FRn(:,3) = FRn(:,3) + Djump * normal(2)

end subroutine compute_flux_pred

subroutine dg_solver(element,update1,update2,update3,dt)
  type(t_element_base), intent(in)				:: element
  real(kind=GRID_SR), intent(in)                                :: dt
  !type(t_update), dimension((_SWE_DG_ORDER+1)**2),intent(in)  :: update1, update2, update3
  type(t_update), dimension(_SWE_DG_ORDER+1),intent(in)  :: update1, update2, update3
  
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3)                 :: bnd_flux_l,bnd_flux_m,bnd_flux_r
  real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,_SWE_DG_DOFS,3) :: bnd_flux_l_st,bnd_flux_m_st,bnd_flux_r_st
  
  real(kind=GRID_SR),Dimension(2,2) :: jacobian,jacobian_inv
  real(kind=GRID_SR),Dimension(3)   :: edge_sizes
  real(kind=GRID_SR)                :: dx
  integer :: i
  real(kind=GRID_SR) :: q(_SWE_DG_DOFS,3)

  associate(data        => element%cell%data_pers, &
            Q_DG        => element%cell%data_pers%Q)

    call data%get_dofs_dg(q)
    
    edge_sizes=element%cell%geometry%get_edge_sizes()
    dx=edge_sizes(1)*cfg%scaling
    
    !---------compute boundary integrals------------------!
    bnd_flux_l(:,1) = matmul(s_m_inv,matmul(s_b_1_l,update1(:)%h   ))
    bnd_flux_l(:,2) = matmul(s_m_inv,matmul(s_b_1_l,update1(:)%p(1)))
    bnd_flux_l(:,3) = matmul(s_m_inv,matmul(s_b_1_l,update1(:)%p(2)))
    bnd_flux_m(:,1) = matmul(s_m_inv,matmul(s_b_2_m,update2(:)%h   ))
    bnd_flux_m(:,2) = matmul(s_m_inv,matmul(s_b_2_m,update2(:)%p(1)))
    bnd_flux_m(:,3) = matmul(s_m_inv,matmul(s_b_2_m,update2(:)%p(2)))
    bnd_flux_r(:,1) = matmul(s_m_inv,matmul(s_b_3_r,update3(:)%h   ))
    bnd_flux_r(:,2) = matmul(s_m_inv,matmul(s_b_3_r,update3(:)%p(1)))
    bnd_flux_r(:,3) = matmul(s_m_inv,matmul(s_b_3_r,update3(:)%p(2)))
    !-----------------------------------------------------!

#if defined(_DEBUG)
    do i=1,_SWE_DG_DOFS
       if(isnan(bnd_flux_l(i,1)).or.isnan(bnd_flux_m(i,1)).or.isnan(bnd_flux_r(i,1)))then
          !print*,"nan bnd"
          exit
       endif
    end do

    do i=1,_SWE_DG_DOFS
       if(isnan(q(i,1))) then
          !print*,"nan q"
          !print*,q
       end if
    end do
#endif
    !!----update dofs----!!
    q=q+(data%Q_DG_UPDATE - bnd_flux_l - bnd_flux_m - bnd_flux_r)* dt/dx
    !!-------------------!!
    call data%set_dofs_dg(q)
  end associate
  
end subroutine dg_solver

        function get_next_time_step_size(h,hu,hv,dx,dry_tolerance) result(dt)
          real(kind=GRID_SR), INTENT(in)  :: h,hu,hv
          real(kind=GRID_SR), INTENT(in)  :: dx
          real(kind=GRID_SR), INTENT(in)  :: dry_tolerance
          real(kind=GRID_SR)              :: dt
          
          real(kind=GRID_SR)              :: waveSpeed

          dt = huge(1.0_GRID_SR)
          if(h > dry_tolerance) then
             waveSpeed =  sqrt(g * (h)) + max(abs(hu/h),abs(hv/h))
             dt = min(dx  / (waveSpeed * (_SWE_DG_ORDER*2.0_GRID_SR+2.0_GRID_SR)),dt)
          end if

        end function get_next_time_step_size


#if defined(_SWE_DG_LIMITER_ALL) || defined(_SWE_DG_LIMITER_HEIGHT)
        subroutine element_op(traversal, section, element)
          type(t_swe_dg_timestep_traversal), intent(inout)    :: traversal
          type(t_grid_section), intent(inout)							    :: section
          type(t_element_base), intent(inout)                 :: element
          integer                                             :: i,j
          real(kind=GRID_SR), dimension(_DMP_NUM_OBSERVABLES) :: minVals
          real(kind=GRID_SR), dimension(_DMP_NUM_OBSERVABLES) :: maxVals
          
          
          call getObservableLimits(element%cell%data_pers,minVals,maxVals)

          do i=1,3
             associate( node => element%nodes(i)%ptr%data_pers)
               do j = 1,_DMP_NUM_OBSERVABLES
                  node%minObservables(j) = min(node%minObservables(j),minVals(j))
                  node%maxObservables(j) = max(node%maxObservables(j),maxVals(j))
               end do
             end associate
          end do
               
        end subroutine element_op

        elemental subroutine node_merge_op(local_node, neighbor_node)
          type(t_node_data), intent(inout)		:: local_node
          type(t_node_data), intent(in)		    :: neighbor_node
          integer :: i
          
          do i = 1,_DMP_NUM_OBSERVABLES
             local_node%data_pers%minObservables(i) = &
                  min(local_node%data_pers%minObservables(i),&
                  neighbor_node%data_pers%maxObservables(i))
             local_node%data_pers%maxObservables(i) = &
                  max(local_node%data_pers%maxObservables(i),&
                  neighbor_node%data_pers%maxObservables(i))
          end do
          
        end subroutine node_merge_op
        
        elemental subroutine node_first_touch_op(traversal, section, node)
          type(t_swe_dg_timestep_traversal), intent(inout)   :: traversal
          type(t_grid_section), intent(in)                   :: section
          type(t_node_data), intent(inout)                   :: node

          node%data_pers%minObservables(:) =  Huge(1.0_GRID_SR)
          node%data_pers%maxObservables(:) = -Huge(1.0_GRID_SR)
          
        end subroutine node_first_touch_op

#endif
        
END MODULE SWE_DG_solver
!_SWE_DG
#endif
