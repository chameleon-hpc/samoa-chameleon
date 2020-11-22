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
  use SWE_Riemann_solver_wrapper
  use c_bind_riemannsolvers
  use SWE_PDE

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
  use SWE_PATCH
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
#		define _GT_CELL_TO_EDGE_OP		cell_to_edge_op_dg

  public cell_to_edge_op_dg

#		include "SFC_generic_traversal_ringbuffer.f90"

  subroutine post_traversal_grid_op_dg(traversal, grid)
    type(t_swe_dg_timestep_traversal), intent(inout)		:: traversal
    type(t_grid), intent(inout)					:: grid

    grid%r_time = grid%r_time + grid%r_dt
    
    call reduce(traversal%i_refinements_issued, traversal%sections%i_refinements_issued, MPI_SUM, .true.)
    call reduce(grid%r_dt_new, grid%sections%elements_alloc%r_dt_new, MPI_MIN, .true.)
    
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

  function cell_to_edge_op_dg(element, edge) result(rep)
    type(t_element_base), intent(in)			    :: element
    type(t_edge_data), intent(in)			    :: edge
    type(num_cell_rep)				            :: rep
    integer                               :: edge_type !1=right, 2=hypotenuse, 3=left
    
    integer :: i,j,k,indx,indx_mir
    integer(kind=BYTE) :: orientation
    type(num_cell_rep) ::rep_fv


#if defined (_DEBUG)
    rep%debug_flag = element%cell%data_pers%debug_flag
#endif
    ! associate(cell_edge => ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%edges)
    !   if(((edge%transform_data%normal(1) ==  cell_edge(2)%normal(1)) .and.&
    !       (edge%transform_data%normal(2) ==  cell_edge(2)%normal(2))).or.&
    !      ((edge%transform_data%normal(1) == -cell_edge(2)%normal(1)) .and.&
    !       (edge%transform_data%normal(2) == -cell_edge(2)%normal(2)))) then
    !      edge_type = 2
    !   else if(((edge%transform_data%normal(1) ==  cell_edge(1)%normal(1)) .and.&
    !            (edge%transform_data%normal(2) ==  cell_edge(1)%normal(2))).or.&
    !           ((edge%transform_data%normal(1) == -cell_edge(1)%normal(1)) .and.&
    !            (edge%transform_data%normal(2) == -cell_edge(1)%normal(2)))) then
    !      edge_type = 1
    !   else
    !      edge_type = 3
    !   endif
    ! end associate

!    call writeFVBoundaryFields(element)

    edge_type = get_edge_type(abs(element%cell%geometry%i_plotter_type),edge%transform_data%normal)
    
    rep%QP(:,:)   = element%cell%data_pers%QP(edge_type  ,:,:)
    rep%FP(:,:,:) = element%cell%data_pers%FP(edge_type,:,:,:)
    
    rep%H  = element%cell%data_pers%QFV(edge_type  ,:,1)
    rep%HU = element%cell%data_pers%QFV(edge_type  ,:,2)
    rep%HV = element%cell%data_pers%QFV(edge_type  ,:,3)
    rep%B  = element%cell%data_pers%QFV(edge_type  ,:,4)
    
    rep%troubled=element%cell%data_pers%troubled
    call getObservableLimits(element%cell%data_pers%Q,rep%minObservables,rep%maxObservables)
    
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

update1%minObservables = rep2%minObservables
update1%maxObservables = rep2%maxObservables
update2%minObservables = rep1%minObservables
update2%maxObservables = rep1%maxObservables
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
integer                                             :: i

normal=(edge%transform_data%normal)/NORM2(edge%transform_data%normal)
update%troubled=rep%troubled
rep_bnd%troubled = rep%troubled
update%minObservables = rep%minObservables
update%maxObservables = rep%maxObservables

if(isDG(rep%troubled)) then
   update_bnd = update
   rep_bnd%QP(:,1) = rep%QP(:,1)
   rep_bnd%QP(:,4) = rep%QP(:,4)
#  if defined(_BOUNDARY)
      if(get_edge_boundary_align(normal)) then
         ! Time-dependent condition, generate virtual wave from data or function
         do i=1, _SWE_DG_ORDER+1
            call get_edge_boundary_Q(grid%r_time, rep_bnd%QP(i,1), &
               rep_bnd%QP(i,2), rep_bnd%QP(i,3), rep_bnd%QP(i,4), .true.)
         end do

         rep_bnd%FP(1, :, :) = flux_1(rep_bnd%QP(:, 1:3), _SWE_DG_ORDER+1)
         rep_bnd%FP(2, :, :) = flux_2(rep_bnd%QP(:, 1:3), _SWE_DG_ORDER+1)
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
else if(isFV(rep%troubled)) then
#  if defined(_BOUNDARY)
      if(get_edge_boundary_align(normal)) then
         ! Time-dependent condition, generate virtual wave from data or function
         do i=1, _SWE_PATCH_ORDER
            update%H(i) = rep%H(i)
            update%B(i) = rep%B(i)
            call get_edge_boundary_Q(grid%r_time, update%H(i), &
               update%HU(i), update%HV(i), update%B(i), .false.)
         end do
      else
#  endif
         ! Generate mirrored wave to reflect out incoming wave
         update%H=rep%H
         update%B=rep%B
         do i=1, _SWE_PATCH_ORDER
            length_flux = dot_product([ rep%HU(i), rep%HV(i) ], normal)
            update%HU(i)=rep%HU(i) - 2.0_GRID_SR*length_flux*normal(1)
            update%HV(i)=rep%HV(i) - 2.0_GRID_SR*length_flux*normal(2)
         end do
#  if defined(_BOUNDARY)
      end if
#  endif
end if   
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
      real(GRID_SR), intent(inout)              :: H, HU, HV
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
real (kind=GRID_SR) :: max_wave_speed, wave_speed=0,dQ_norm,b_min,b_max
real (kind=GRID_SR) :: refinement_threshold = 10.25_GRID_SR/_SWE_DG_ORDER
real (kind = GRID_SR), dimension (_SWE_DG_DOFS) :: H_old, HU_old, HV_old, B_old
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

!------- Update cell status and compute next timestep size --------!
call updateCellStatusPre(data,update1,update2,update3)

edge_sizes=element%cell%geometry%get_edge_sizes()

if(isDG(data%troubled)) then
   data%troubled = DG
   H_old  = data%Q%H
   B_old  = data%Q%B
   HU_old = data%Q%p(1)
   HV_old = data%Q%p(2)

   call getObservableLimits(element%cell%data_pers%Q,minVals,maxVals)

   ! !print*,"update1"
   ! !print*,update1%flux(:)%h
   ! !print*,update1%flux(:)%p(1)
   ! !print*,update1%flux(:)%p(2)
   ! !print*,"update2"
   ! !print*,update2%flux(:)%h
   ! !print*,update2%flux(:)%p(1)
   ! !print*,update2%flux(:)%p(2)
   ! !print*,"update3"
   ! !print*,update3%flux(:)%h
   ! !print*,update3%flux(:)%p(1)
   ! !print*,update3%flux(:)%p(2)
   
   ! if(any(abs(update1%flux%p(1)+update2%flux%p(1)+update3%flux%p(1)) > 1.5d-4))then
   !    stop
   ! endif

   do i = 1,_SWE_DG_ORDER+1
      tmp%flux(i) = update3%flux(_SWE_DG_ORDER+2-i)
   end do
   update3%flux = tmp%flux


   call dg_solver(element,update1%flux,update2%flux,update3%flux,section%r_dt)
   
   if(isWetDryInterface(data%Q%H)) then
      data%troubled = WET_DRY_INTERFACE
   else if(checkDMP(element%cell%data_pers%Q,minVals,maxVals,update1,update2,update3)) then
      data%troubled = TROUBLED
   end if
   
   if(isFV(data%troubled)) then
      !--if troubled or drying perform rollback--!
      data%Q%H   = H_old
      data%Q%p(1)= HU_old
      data%Q%p(2)= HV_old
      data%Q%B   = B_old

      call apply_phi(data%Q(:)%h+data%Q(:)%b,data%h)
      call apply_phi(data%Q(:)%p(1)         ,data%hu)
      call apply_phi(data%Q(:)%p(2)         ,data%hv)
      call apply_phi(data%Q(:)%b            ,data%b)
   end if
end if

if(isFV(data%troubled)) then
   !-----Call FV patch solver----!    
   call fv_patch_solver(traversal, section, element, update1, update2, update3)
end if


!------- Update cell status and compute next timestep size --------!
call updateCellStatus(data)

dx = cfg%scaling *  edge_sizes(1)
if(isDG(data%troubled)) then
   do i=1,_SWE_DG_DOFS
      ! if(get_next_time_step_size(data%Q(i)%h,data%Q(i)%p(1),data%Q(i)%p(2),dx,cfg%dry_tolerance) < section%r_dt_new)then
      !    print*,"DG"
      !    print*,data%troubled
      !    print*,data%Q(i)%h,data%Q(i)%p(1),data%Q(i)%p(2)
      !    print*,get_next_time_step_size(data%Q(i)%h,data%Q(i)%p(1),data%Q(i)%p(2),dx,cfg%dry_tolerance)
      ! endif
      section%r_dt_new = min(section%r_dt_new, get_next_time_step_size(data%Q(i)%h,data%Q(i)%p(1),data%Q(i)%p(2),dx,cfg%dry_tolerance))
   end do
else
   do i=1,_SWE_PATCH_ORDER_SQUARE
      ! if(get_next_time_step_size(data%h(i)-data%b(i),data%hu(i),data%hv(i),dx,cfg%dry_tolerance) < section%r_dt_new)then
      !    print*,"FV"
      !    print*,data%h(i)-data%b(i),data%hu(i),data%hv(i)
      !    print*,get_next_time_step_size(data%h(i)-data%b(i),data%hu(i),data%hv(i),dx,cfg%dry_tolerance)
      ! endif
      section%r_dt_new = min(section%r_dt_new, get_next_time_step_size(data%h(i)-data%b(i),data%hu(i),data%hv(i),dx,cfg%dry_tolerance))
end do
endif
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
   ! if(error > 10.d-8) then
   !    !print*,"error"
   !    !print*,error
   !    !print*,section%max_error
   ! end if
   
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
real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1)                  :: hRoe,bm
real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1)                  :: Deta
real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1)                  :: Djump
real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,3)                :: Fn_avg
real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,3)                :: Q_rus

real(kind = GRID_SR), intent(in)  :: normal(2)

real(kind = GRID_SR)	          :: epsilon
real(kind=GRID_SR)	          :: vL, vR, alpha
integer                           :: i

Fn_avg = 0.5_GRID_SR * ( FRn - FLn )
! print*,"Fn_avg"
! print*,Fn_avg(:,1)
! print*,Fn_avg(:,2)
! print*,Fn_avg(:,3)

VelL(:) = QL(:,2)/QL(:,1) * normal(1) +  QL(:,3)/QL(:,1) * normal(2)
VelR(:) = QR(:,2)/QR(:,1) * normal(1) +  QR(:,3)/QR(:,1) * normal(2)

alpha=max(maxval(abs(VelL) + sqrt(g * QL(:,1))),&
     maxval(abs(VelR) + sqrt(g * QR(:,1))))

hRoe  = 0.5_GRID_SR * (QR(:,1) + QL(:,1))
bm    = max(QR(:,4),QL(:,4))
Deta  = max(QR(:,1) + QR(:,4) - bm, 0.0) - max(QL(:,1) + QL(:,4) - bm, 0.0)

Djump = 0.5_GRID_SR * g * hRoe * Deta
!print*,"Djump"
!print*,Djump

Q_rus(:,  1) = 0.5_GRID_SR * alpha * Deta
Q_rus(:,2:3) = 0.5_GRID_SR * alpha * (QR(:,2:3) - QL(:,2:3))

FLn    =  Fn_avg - Q_rus
FRn    =  Fn_avg + Q_rus

FLn(:,1) = FLn(:,1) 
FLn(:,2) = FLn(:,2) + Djump * normal(1)
FLn(:,3) = FLn(:,3) + Djump * normal(2)

FRn(:,1) = FRn(:,1)    
FRn(:,2) = FRn(:,2) - Djump * normal(1)
FRn(:,3) = FRn(:,3) - Djump * normal(2)

! if (any(abs(FLn(:,2)) > 10.0e-5)) then
! print*,"FLn1"   
! print*,FLn
! print*,"FRn2"   
! print*,FRn
!    !print*,"QR"   
!    !print*,QR(:,1)
!    !print*,QR(:,2)
!    !print*,QR(:,3)
!    !print*,QR(:,4)
!    !print*,"QL"   
!    !print*,QL(:,1)
!    !print*,QL(:,2)
!    !print*,QL(:,3)
!    !print*,QL(:,4)
!    !print*,"Djump"
!    !print*,Djump
!    !print*,"Qrus"
!    !print*,Q_rus
!    !print*,"MARK"
! end if

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
          !print*,bnd_flux_l(:,:)
          !print*,bnd_flux_m(:,:)
          !print*,bnd_flux_r(:,:)
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

subroutine fv_patch_solver(traversal, section, element, update1, update2, update3)
            class(t_traversal), intent(inout)		  :: traversal
            type(t_grid_section), intent(inout)	  :: section
            type(t_element_base), intent(inout)	  :: element
            type(num_cell_update), intent(inout)  :: update1, update2, update3
            integer                               :: i, j, ind
            type(num_cell_update)                 :: tmp !> ghost cells in correct order 
            real(kind = GRID_SR)                  :: volume, edge_lengths(3), maxWaveSpeed, dQ_max_norm, dt_div_volume, maxWaveSpeedLocal, maxWaveSpeed_node


            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_ORDER_SQUARE)                :: dQ_H, dQ_HU, dQ_HV !> deltaQ, used to compute cell updates
            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)           :: hL, huL, hvL, bL
            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)           :: hR, huR, hvR, bR
            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)           :: upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR

#if !defined (_SWE_USE_PATCH_SOLVER)
            type(t_state), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE)      :: edges_a, edges_b
            type(t_update)                                              :: update_a, update_b
            real(kind = GRID_SR), dimension(2,3)                        :: normals
#endif

            type(t_update),  DIMENSION((_SWE_DG_ORDER+1)**2)		:: flux1, flux2, flux3
            type(t_state), dimension((_SWE_DG_ORDER+1)**2)              :: edge_l, edge_m,edge_r
            real(kind= GRID_SR):: dt,dx,delta(3),max_neighbour(3),min_neighbour(3),data_max_val(3),data_min_val(3)
            integer:: indx,indx_mir,k,edgesLeft
            real (kind = GRID_SR), dimension (_SWE_PATCH_ORDER_SQUARE):: H_old, HU_old, HV_old
            logical :: drying,troubled,coarsen,refine
            real (kind=GRID_SR) :: refinement_threshold = 0.50_GRID_SR
            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)                :: normals_x, normals_y


#if defined(_OPT_KERNELS)            
            !DIR$ ASSUME_ALIGNED hL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED hR: ALIGNMENT
            !DIR$ ASSUME_ALIGNED huL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED huR: ALIGNMENT
            !DIR$ ASSUME_ALIGNED hvL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED hvR: ALIGNMENT
            !DIR$ ASSUME_ALIGNED bL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED bR: ALIGNMENT
            !DIR$ ASSUME_ALIGNED upd_hL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED upd_hR: ALIGNMENT
            !DIR$ ASSUME_ALIGNED upd_huL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED upd_huR: ALIGNMENT
            !DIR$ ASSUME_ALIGNED upd_hvL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED upd_hvR: ALIGNMENT
#endif
            upd_hL = 0.0_GRID_SR
            upd_hR = 0.0_GRID_SR
            upd_huL= 0.0_GRID_SR
            upd_huR= 0.0_GRID_SR
            upd_hvL= 0.0_GRID_SR
            upd_hvR= 0.0_GRID_SR

            
            ! using patches, but applying geoclaw solvers on single edges            
            ! the normals are only needed in this case.            
            ! copy/compute normal vectors
            ! normal for type 2 edges is equal to the 2nd edge's normal
            normals(:,2) = element%edges(2)%ptr%transform_data%normal
            ! normal for types 1 and 3 depend on cell orientation.
            ! notice that normal for type 1 points inwards. That's why it is multiplied by -1.
            if (element%cell%geometry%i_plotter_type < 0) then ! orientation = backward
               normals(:,1) = - element%edges(1)%ptr%transform_data%normal
               normals(:,3) = element%edges(3)%ptr%transform_data%normal
            else ! orientation = forward, so reverse edge order
               normals(:,1) = - element%edges(3)%ptr%transform_data%normal
               normals(:,3) = element%edges(1)%ptr%transform_data%normal
            end if

            ! init some variables
            dQ_H = 0.0_GRID_SR
            dQ_HU = 0.0_GRID_SR
            dQ_HV = 0.0_GRID_SR
            maxWaveSpeed = 0.0_GRID_SR
            volume = cfg%scaling * cfg%scaling * element%cell%geometry%get_volume() / (_SWE_PATCH_ORDER_SQUARE)
            dt_div_volume = section%r_dt / volume
            edge_lengths = cfg%scaling * element%cell%geometry%get_edge_sizes() / _SWE_PATCH_ORDER

            associate(data => element%cell%data_pers, geom => SWE_PATCH_geometry)

              ! if(element%cell%data_pers%troubled .eq. 2) then
              !    print*,"patch start"
              !    print*,data%h 
              !    print*,data%hu
              !    print*,data%hv
              !    print*,data%b 
              ! endif


              ! copy cell values to arrays edges_a and edges_b
              ! obs: cells with id > number of cells are actually ghost cells and come from edges "updates"
              ! see t_SWE_PATCH_geometry for an explanation about ghost cell ordering
              do i=1, _SWE_PATCH_NUM_EDGES, _SWE_PATCH_SOLVER_CHUNK_SIZE ! i -> init of chunk

                 ! if this is the last chunk and it is not a full chunk, it is necessary to set everything to 0 to avoid using last iteration values
                 if (i + _SWE_PATCH_SOLVER_CHUNK_SIZE -1 > _SWE_PATCH_NUM_EDGES) then
                    hL = 0.0_GRID_SR
                    huL = 0.0_GRID_SR
                    hvL = 0.0_GRID_SR
                    bL = 0.0_GRID_SR
                    hR = 0.0_GRID_SR
                    huR = 0.0_GRID_SR
                    hvR = 0.0_GRID_SR
                    bR = 0.0_GRID_SR
                 end if
                 !number of edges that still haven't been processed
                 edgesLeft = _SWE_PATCH_NUM_EDGES - i + 1

                 do j=1, min(edgesLeft,_SWE_PATCH_SOLVER_CHUNK_SIZE) !j -> position inside chunk
                    ind = i + j - 1 ! actual index
                    
                    normals_x(j) = normals(1,geom%edges_orientation(ind))
                    normals_y(j) = normals(2,geom%edges_orientation(ind))

                    ! left
                    if (geom%edges_a(ind) <= _SWE_PATCH_ORDER_SQUARE) then
                       hL(j) = data%H(geom%edges_a(ind))
                       huL(j) = data%HU(geom%edges_a(ind))
                       hvL(j) = data%HV(geom%edges_a(ind))
                       bL(j) = data%B(geom%edges_a(ind))
                    else if (geom%edges_a(ind) <= _SWE_PATCH_ORDER_SQUARE + _SWE_PATCH_ORDER) then
                       hL(j) = update1%H(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE)
                       huL(j) = update1%HU(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE)
                       hvL(j) = update1%HV(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE)
                       bL(j) = update1%B(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE)
                    else if (geom%edges_a(ind) <= _SWE_PATCH_ORDER_SQUARE + 2*_SWE_PATCH_ORDER) then
                       hL(j) = update2%H(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                       huL(j) = update2%HU(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                       hvL(j) = update2%HV(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                       bL(j) = update2%B(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                    else 
                       hL(j) = update3%H(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                       huL(j) = update3%HU(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                       hvL(j) = update3%HV(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                       bL(j) = update3%B(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                    end if

                    ! right
                    if (geom%edges_b(ind) <= _SWE_PATCH_ORDER_SQUARE) then
                       hR(j) = data%H(geom%edges_b(ind))
                       huR(j) = data%HU(geom%edges_b(ind))
                       hvR(j) = data%HV(geom%edges_b(ind))
                       bR(j) = data%B(geom%edges_b(ind))
                    else if (geom%edges_b(ind) <= _SWE_PATCH_ORDER_SQUARE + _SWE_PATCH_ORDER) then
                       hR(j) = update1%H(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE)
                       huR(j) = update1%HU(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE)
                       hvR(j) = update1%HV(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE)
                       bR(j) = update1%B(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE)
                    else if (geom%edges_b(ind) <= _SWE_PATCH_ORDER_SQUARE + 2*_SWE_PATCH_ORDER) then
                       hR(j) = update2%H(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                       huR(j) = update2%HU(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                       hvR(j) = update2%HV(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                       bR(j) = update2%B(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                    else 
                       hR(j) = update3%H(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                       huR(j) = update3%HU(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                       hvR(j) = update3%HV(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                       bR(j) = update3%B(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                    end if
                 end do

                 ! compute net_updates -> solve Riemann problems within chunk
#if defined(_OPT_KERNELS)
                 ! Vectorization! (Requires OpenMP 4.0 or later)
                 !$OMP SIMD PRIVATE(maxWaveSpeedLocal) REDUCTION(max: maxWaveSpeed)
                 do j=1, min(edgesLeft, _SWE_PATCH_SOLVER_CHUNK_SIZE)
! #if defined(_PATCH_VEC_INLINE)
!                     ! Warning: inlining this subroutine into an OMP SIMD loop may cause
!                     ! bugs and incorrect calculations depending on the compiler version
                    !                     ! Check the 
!                     ! Recommended compiler: ifort 17.0 or 19.0
!                     !DIR$ FORCEINLINE
! #endif !_PATCH_VEC_INLINE
                    call SWE_compute_fluxes(normals_x(j), normals_y(j), hL(j), hR(j), huL(j), huR(j), hvL(j), hvR(j), bL(j), bR(j), upd_hL(j), upd_hR(j), upd_huL(j), upd_huR(j), upd_hvL(j), upd_hvR(j), maxWaveSpeedLocal)
                    maxWaveSpeed = max(maxWaveSpeed, maxWaveSpeedLocal)
                 end do


#else 
                 ! using original geoclaw solver
                                     
                 do j=1, min(edgesLeft,_SWE_PATCH_SOLVER_CHUNK_SIZE) !j -> position inside chunk
                    ind = i + j - 1 ! actual index

                    edges_a(j)%h = hL(j)
                    edges_a(j)%p(1) = huL(j)
                    edges_a(j)%p(2) = hvL(j)
                    edges_a(j)%b = bL(j)
                    edges_b(j)%h = hR(j)
                    edges_b(j)%p(1) = huR(j)
                    edges_b(j)%p(2) = hvR(j)
                    edges_b(j)%b = bR(j)
                    maxWaveSpeedLocal = 0.0_GRID_SR
                    call compute_geoclaw_flux(normals(:,geom%edges_orientation(ind)), edges_a(j), edges_b(j), update_a, update_b, maxWaveSpeedLocal)
                    maxWaveSpeed=max(maxWaveSpeed,maxWaveSpeedLocal)
                    

                    upd_hL(j)  = update_a%h
                    upd_huL(j) = update_a%p(1)
                    upd_hvL(j) = update_a%p(2)
                    upd_hR(j)  = update_b%h
                    upd_huR(j) = update_b%p(1)
                    upd_hvR(j) = update_b%p(2)
                 end do
#                       endif
                  ! compute dQ
                  do j=1, min(edgesLeft,_SWE_PATCH_SOLVER_CHUNK_SIZE)
                     ind = i + j - 1 ! actual index

                     if (geom%edges_a(ind) <= _SWE_PATCH_ORDER_SQUARE) then !ignore ghost cells
                        dQ_H(geom%edges_a(ind)) = dQ_H(geom%edges_a(ind)) + upd_hL(j) * edge_lengths(geom%edges_orientation(ind))
                        dQ_HU(geom%edges_a(ind)) = dQ_HU(geom%edges_a(ind)) + upd_huL(j) * edge_lengths(geom%edges_orientation(ind))
                        dQ_HV(geom%edges_a(ind)) = dQ_HV(geom%edges_a(ind)) + upd_hvL(j) * edge_lengths(geom%edges_orientation(ind))
                     end if
                     if (geom%edges_b(ind) <= _SWE_PATCH_ORDER_SQUARE) then
                        dQ_H(geom%edges_b(ind)) = dQ_H(geom%edges_b(ind)) + upd_hR(j) * edge_lengths(geom%edges_orientation(ind))
                        dQ_HU(geom%edges_b(ind)) = dQ_HU(geom%edges_b(ind)) + upd_huR(j) * edge_lengths(geom%edges_orientation(ind))
                        dQ_HV(geom%edges_b(ind)) = dQ_HV(geom%edges_b(ind)) + upd_hvR(j) * edge_lengths(geom%edges_orientation(ind))
                     end if
                  end do
               end do

               !set refinement condition -> Here I am using the same ones as in the original no-patches implementation, but considering only the max value.
               element%cell%geometry%refinement = 0
               dQ_max_norm = maxval(abs(dQ_H))

               ! ! if the water level falls below the dry tolerance, set water surface to 0 and velocity to 0
                              
               dQ_H  = dQ_H * (-dt_div_volume)
               dQ_HU = dQ_HU * (-dt_div_volume)
               dQ_HV = dQ_HV * (-dt_div_volume)

            !    if(data%troubled == 4)then
            !    if(update1%troubled == -2)then
            !       print*,"data"
            !       print*,data%H 
            !       print*,data%HU
            !       print*,data%HV
            !       print*,"update1"
            !       print*,update1%H(:)
            !       print*,update1%HU(:)
            !       print*,update1%HV(:)
            !       print*,update1%B(:)
            !       print*,"dQ"
            !       print*,dQ_H
            !       print*,dQ_HU
            !       print*,dQ_HV
            !    end if
            !    if(update2%troubled == -2)then
            !       print*,"data"
            !       print*,data%H 
            !       print*,data%HU
            !       print*,data%HV
            !       print*,"update2"
            !       print*,update2%H(:)
            !       print*,update2%HU(:)
            !       print*,update2%HV(:)
            !       print*,update2%B(:)
            !       print*,"dQ"
            !       print*,dQ_H
            !       print*,dQ_HU
            !       print*,dQ_HV
            !    end if
            !    if(update3%troubled == -2)then
            !       print*,"data"
            !       print*,data%H 
            !       print*,data%HU
            !       print*,data%HV
            !       print*,"update3"
            !       print*,update3%H(:)
            !       print*,update3%HU(:)
            !       print*,update3%HV(:)
            !       print*,update3%B(:)
            !       print*,"dQ"
            !       print*,dQ_H
            !       print*,dQ_HU
            !       print*,dQ_HV                  
            !    end if
            ! end if

            ! if(any(abs(dQ_HU) > 10.0e-8))then
            !    stop
            ! end if

               ! if land is flooded, init water height to dry tolerance and
               ! velocity to zero
               where (data%H < data%B + cfg%dry_tolerance .and. dQ_H > 0.0_GRID_SR)
#if defined (_ASAGI)                  
                  data%H = min(0.0_GRID_SR,data%B)
#else
                  !                  data%H = data%B + cfg%dry_tolerance
                  data%H = data%B + cfg%dry_tolerance 
#endif                  
                  data%HU = 0.0_GRID_SR
                  data%HV = 0.0_GRID_SR
               end where

               data%H  = data%H + dQ_H
               data%HU = data%HU + dQ_HU
               data%HV = data%HV + dQ_HV
               
               ! if the water level falls below the dry tolerance, set water level to 0 and velocity to 0          
               where (data%H < data%B + cfg%dry_tolerance)
#if defined (_ASAGI)                  
                  data%H = min(0.0_GRID_SR,data%B)
#else
                 data%H = data%B 
!                  data%H = data%B

#endif                  
                  data%HU = 0.0_GRID_SR
                  data%HV = 0.0_GRID_SR
               end where


               ! if(all(data%h - data%b < cfg%dry_tolerance * 10.0d4)) then
               !    where (data%HU > cfg%dry_tolerance * 10.0d-2) 
               !       data%HU = 0.0_GRID_SR
               !    end where
               !    where (data%HV > cfg%dry_tolerance * 10.0d-2) 
               !       data%HV = 0.0_GRID_SR
               !    end where
               ! end if
               
               if(.not.isCoast(element%cell%data_pers%troubled)) then
                  call apply_mue(data%h ,data%Q%h)
                  call apply_mue(data%hu,data%Q%p(1))
                  call apply_mue(data%hv,data%Q%p(2))
!                  call apply_mue(data%b ,data%Q%b)                                  
                  data%Q%h=data%Q%h-data%Q%b
               end if
          end associate

        end subroutine fv_patch_solver

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

        subroutine compute_geoclaw_flux(normal, QL, QR, fluxL, fluxR,max_wave_speed)
          type(t_state), intent(in)           :: QL, QR
          type(t_update), intent(out)         :: fluxL, fluxR
          real(kind = GRID_SR), intent(in)    :: normal(2)
          real(kind = GRID_SR)                :: max_wave_speed
          
          real(kind = GRID_SR)				:: transform_matrix(2, 2)
          real(kind = GRID_SR)			    :: net_updatesL(3), net_updatesR(3)
          real(kind = GRID_SR)                :: pL(2), pR(2), hL, hR, bL, bR
          
          transform_matrix(1, :) = normal
          transform_matrix(2, :) = [-normal(2), normal(1)]
          
          pL = matmul(transform_matrix, QL%p)
          pR = matmul(transform_matrix, QR%p)
          hL = QL%h - QL%b
          hR = QR%h - QR%b
          bL = QL%b
          bR = QR%b
          
#           if defined(_FWAVE_FLUX)
          call c_bind_geoclaw_solver(GEOCLAW_FWAVE, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_AUG_RIEMANN_FLUX)
          call c_bind_geoclaw_solver(GEOCLAW_AUG_RIEMANN, 10, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_HLLE_FLUX)
          call compute_updates_hlle_single(hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, net_updatesL, net_updatesR, max_wave_speed)
#           endif
          
          fluxL%h = net_updatesL(1)
          fluxL%p = matmul(net_updatesL(2:3), transform_matrix)
          fluxR%h = net_updatesR(1)
          fluxR%p = matmul(net_updatesR(2:3), transform_matrix)
        end subroutine compute_geoclaw_flux
        
END MODULE SWE_DG_solver
!_SWE_DG
#endif
