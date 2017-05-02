! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2017 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE
! Author: Leonhard Rannabauer rannabau (at) in.tum.de

#if defined(_SWE_DG)
#include "Compilation_control.f90"

MODULE SWE_DG_solver
  use SFC_edge_traversal

  use Samoa_swe
  use c_bind_riemannsolvers

  use SWE_DG_Matrices
  use SWE_DG_Predictor
  use SWE_Euler_Timestep
  use SWE_PATCH
#               if defined(_HLLE_FLUX)
  use SWE_HLLE
#               endif


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

  interface edge_first_touch_op_dg
     module procedure edge_first_touch_array_op_dg
     module procedure edge_first_touch_scalar_op_dg
  end interface edge_first_touch_op_dg


#		define _GT_NAME	           		t_swe_dg_timestep_traversal

#		define _GT_EDGES

#		define _GT_PRE_TRAVERSAL_OP      	pre_traversal_op_dg
#		define _GT_PRE_TRAVERSAL_GRID_OP	pre_traversal_grid_op_dg
#		define _GT_POST_TRAVERSAL_GRID_OP	post_traversal_grid_op_dg

#		define _GT_SKELETON_OP			skeleton_op_dg
#		define _GT_BND_SKELETON_OP		bnd_skeleton_op_dg
#		define _GT_ELEMENT_OP                   element_op_dg

#		define _GT_CELL_UPDATE_OP		cell_update_op_dg
  
#		define _GT_CELL_TO_EDGE_OP		cell_to_edge_op_dg

#               define _GT_EDGE_FIRST_TOUCH_OP          edge_first_touch_op_dg

  public cell_to_edge_op_dg

#		include "SFC_generic_traversal_ringbuffer.f90"

!  use SWE_Euler_Timestep

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
    grid%r_dt = min(cfg%r_max_time, grid%r_dt)    
    call scatter(grid%r_time, grid%sections%elements_alloc%r_time)
  end subroutine post_traversal_grid_op_dg


  
  subroutine edge_first_touch_scalar_op_dg(traversal, section, edge)
    type(t_swe_dg_timestep_traversal)  ::traversal
    type(t_grid_section), intent(inout)::section
    type(t_edge_data), intent(inout)   ::edge
  end subroutine edge_first_touch_scalar_op_dg

  subroutine edge_first_touch_array_op_dg(traversal, section, edge)
    type(t_swe_dg_timestep_traversal)  ::traversal
    type(t_grid_section), intent(inout)::section
    type(t_edge_data), intent(inout)   ::edge(:)
  end subroutine edge_first_touch_array_op_dg


  subroutine pre_traversal_grid_op_dg(traversal, grid)
    type(t_swe_dg_timestep_traversal), intent(inout)		:: traversal
    type(t_grid), intent(inout)							    :: grid
  end subroutine pre_traversal_grid_op_dg


  subroutine pre_traversal_op_dg(traversal, section)
    type(t_swe_dg_timestep_traversal), intent(inout)				:: traversal
    type(t_grid_section), intent(inout)							:: section
    !this variable will be incremented for each cell with a refinement request
    traversal%i_refinements_issued = 0_GRID_DI
    section%r_dt_new = huge(1.0_SR)

  end subroutine pre_traversal_op_dg

  subroutine element_op_dg(traversal, section, element)
    type(t_swe_dg_timestep_traversal), intent(inout)				:: traversal
    type(t_grid_section), intent(inout)							:: section
    type(t_element_base), intent(inout)						:: element
  end subroutine element_op_dg

  function cell_to_edge_op_dg(element, edge) result(rep)
    type(t_element_base), intent(in)						:: element
    type(t_edge_data), intent(in)						    :: edge
    type(num_cell_rep)										:: rep
    integer                                             :: edge_type=-5 !1=right, 2=hypotenuse, 3=left
    integer :: i,k,indx,indx_mir
    integer(kind=BYTE) :: orientation
    type(num_cell_rep) ::rep_fv

    associate(cell_edge => ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%edges, normal => edge%transform_data%normal)
      do i=1,3
         if ( normal(1) == cell_edge(i)%normal(1) .and. normal(2) == cell_edge(i)%normal(2) )   then
            edge_type = i
         else if( normal(1) == -cell_edge(i)%normal(1) .and. normal(2) == -cell_edge(i)%normal(2) ) then
            edge_type = -i
         endif
      end do

    end associate

    if(element%cell%data_pers%troubled.le.0) then

     do k=0,_SWE_DG_ORDER
        do i=0,_SWE_DG_ORDER
           indx_mir=i+1+k*(_SWE_DG_ORDER+1)
           select case (edge_type)
           case(-1 ,1) !right
              indx=_SWE_DG_DOFS-(_SWE_DG_ORDER-i+1)*(_SWE_DG_ORDER-i+2)/2 +1+k*_SWE_DG_DOFS
              rep%Q_DG_P(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
              rep%Q_DG_P(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
              rep%Q_DG_P(indx_mir)%b = element%cell%data_pers%Q_DG(indx-k*_SWE_DG_DOFS)%b
           case(-2,2) !mid
              indx=_SWE_DG_DOFS-(_SWE_DG_ORDER-i)*(_SWE_DG_ORDER-i+1)/2 + k*_SWE_DG_DOFS
              rep%Q_DG_P(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
              rep%Q_DG_P(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
              rep%Q_DG_P(indx_mir)%b = element%cell%data_pers%Q_DG(indx-k*_SWE_DG_DOFS)%b
           case(-3 ,3) !left
              indx=1+i+k*_SWE_DG_DOFS
              rep%Q_DG_P(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
              rep%Q_DG_P(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
              rep%Q_DG_P(indx_mir)%b = element%cell%data_pers%Q_DG(indx-k*_SWE_DG_DOFS)%b
           case default
              !print*,"ERROR: In cell_to_edge_op edge type not known: ",edge_type
              stop
              
           end select
        end do
     end do
  end if

  rep%Q_DG = element%cell%data_pers%Q_DG
  rep%troubled=element%cell%data_pers%troubled
  
end function cell_to_edge_op_dg

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

subroutine skeleton_scalar_op_dg(traversal, grid, edge, rep1, rep2, update1, update2)
type(t_swe_dg_timestep_traversal), intent(in)		:: traversal
type(t_grid_section), intent(in)			:: grid
type(t_edge_data), intent(in)				:: edge
type(num_cell_rep), intent(in)				:: rep1, rep2
type(num_cell_update), intent(out)			:: update1, update2
type(t_update),Allocatable                              :: flux_temp(:)
type(t_state),Allocatable                               :: rep_temp(:)
integer                                                 :: i,j


if(rep1%troubled.ne.1 .and. rep2%troubled.ne.2) then
   
   !        ______
   ! |6 \   \3 2 1|
   ! |4 5 \   \5 4|
   ! |1_2_3\   \ 6|
   
   !---- Dofs for hypothenuse need to be permuted ----!
   if(edge%transform_data%index.eq.2) then
      Allocate(rep_temp((_SWE_DG_ORDER+1)**2))
      do i=0,_SWE_DG_ORDER
         do j=1,_SWE_DG_ORDER+1
            rep_temp(i*(_SWE_DG_ORDER+1)+j) = rep2%Q_DG_P((1+i)*(_SWE_DG_ORDER+1)-j+1)
         end do
      end do
   else
      rep_temp=rep2%Q_DG_P
   end if

   
   call compute_flux_pred(edge%transform_data%normal,rep1%Q_DG_P,rep_temp,update1%flux,update2%flux)

   ! if(edge%transform_data%is_hypothenuse) then
   !       do i=0,_SWE_DG_ORDER
   !          do j=1,_SWE_DG_ORDER+1
   !             flux_temp(i*(_SWE_DG_ORDER+1)+j) = update1%flux((1+i)*(_SWE_DG_ORDER+1)-j+1)
   !          end do
   !       end do
   !       update1%flux=flux_temp
   ! end if

   !      if(rep2%inverted) then
   
   !---- Result for hypothenuse needs to be permuted back----!
   if(edge%transform_data%index.eq.2) then
      Deallocate(rep_temp)
      Allocate(flux_temp((_SWE_DG_ORDER+1)**2))
      do i=0,_SWE_DG_ORDER
         do j=1,_SWE_DG_ORDER+1
            flux_temp(i*(_SWE_DG_ORDER+1)+j) = update2%flux((1+i)*(_SWE_DG_ORDER+1)-j+1)
         end do
      end do
      update2%flux=flux_temp
      Deallocate(flux_temp)
   end if

end if   

update1%Q_DG          =rep2%Q_DG
update2%Q_DG          =rep1%Q_DG

update1%troubled      =rep2%troubled
update2%troubled      =rep1%troubled


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
type(t_swe_dg_timestep_traversal), intent(in)  :: traversal
type(t_grid_section), intent(in)		    :: grid
type(t_edge_data), intent(in)		    :: edge
type(num_cell_rep), intent(in)		    :: rep
type(num_cell_update), intent(out)		    :: update
type(num_cell_update)		    :: temp_update
real(kind=GRID_SR)                             :: normal_normed(2)
real(kind=GRID_SR)                             :: length_flux
integer                                        :: i
type(t_state),Dimension((_SWE_DG_ORDER+1)**2)                      :: temp_Q_DG_P

normal_normed=(edge%transform_data%normal)/NORM2(edge%transform_data%normal)


if(rep%troubled.le.0) then
temp_Q_DG_P(:)%h = rep%Q_DG_P(:)%h
temp_Q_DG_P(:)%b = rep%Q_DG_P(:)%b

update%troubled=rep%troubled

!---Predictor velocities---!
!---TODO: Scons variable or config for boundary---!
do i=1,(_SWE_DG_ORDER+1)**2
   !                          reflecting boundary
   length_flux = dot_product(rep%Q_DG_P(i)%p, normal_normed)

   temp_Q_DG_P(i)%p(1) = rep%Q_DG_P(i)%p(1)-2.0_GRID_SR*length_flux*normal_normed(1)
   temp_Q_DG_P(i)%p(2) = rep%Q_DG_P(i)%p(2)-2.0_GRID_SR*length_flux*normal_normed(2)

   !simple outflowing boundary
   !temp_Q_DG_P(i)%p(1) =  rep%Q_DG_P(i)%p(1)
   !temp_Q_DG_P(i)%p(2) =  rep%Q_DG_P(i)%p(2)
   !                            end if

   !zero vel bnd
   !                          update%Q_DG_P(i)%p(1) =  0
   !                          update%Q_DG_P(i)%p(2) =  0
end do

call compute_flux_pred(edge%transform_data%normal,rep%Q_DG_P,temp_Q_DG_P,update%flux,temp_update%flux)


end if

update%Q_DG(:)%h = rep%Q_DG(:)%h
update%Q_DG(:)%b = rep%Q_DG(:)%b



!---DG velocities---!
do i=1,(_SWE_DG_DOFS)
   
   length_flux = dot_product(rep%Q_DG(i)%p, normal_normed)
   ! reflecting
   update%Q_DG(i)%p(1) = rep%Q_DG(i)%p(1)-2.0_GRID_SR*length_flux*normal_normed(1)
   update%Q_DG(i)%p(2) = rep%Q_DG(i)%p(2)-2.0_GRID_SR*length_flux*normal_normed(2)
   !outflow
   !update%Q_DG(i)%p(1) = rep%Q_DG(i)%p(1)
   !update%Q_DG(i)%p(2) = rep%Q_DG(i)%p(2)
end do


end subroutine bnd_skeleton_scalar_op_dg

!cup_time
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
real (kind=GRID_SR) :: max_wave_speed, wave_speed=0,dQ_norm
real (kind=GRID_SR) :: refinement_threshold=0.5_GRID_SR
real (kind = GRID_SR), dimension (_SWE_DG_DOFS) :: H_old, HU_old, HV_old, B_old
real (kind = GRID_SR), dimension (_SWE_PATCH_ORDER_SQUARE) :: H, HU, HV
logical :: drying,troubled,neighbours_troubled
integer :: i_depth



if (element%cell%geometry%i_plotter_type > 0) then ! if orientation = forward, reverse updates
   tmp=update1
   update1=update3
   update3=tmp
end if


neighbours_troubled=(update1%troubled.ge.1).or.(update2%troubled.ge.1).or.(update3%troubled.ge.1)

!!!print*,"Troubled solver "
!!!print*,update1%troubled
!!!print*,update2%troubled
!!!print*,update3%troubled

associate(data => element%cell%data_pers)
  
  
if(neighbours_troubled) then
   data%troubled=merge(data%troubled,2,data%troubled.ge.1)
!   call element%cell%data_pers%convert_dg_to_fv_bathymetry()
end if

if(data%troubled.le.0) then
   data%troubled=0

!!!!!print*,"updatesl"
!!!!!print*,edge_l(:)%H
!!!!!print*,edge_l(:)%p(1)
!!!!!print*,edge_l(:)%p(2)
!!!!!print*,edge_l(:)%b
!!!!!print*,update1%Q_DG_P(:)%H
!!!!!print*,update1%Q_DG_P(:)%p(1)
!!!!!print*,update1%Q_DG_P(:)%p(2)
!!!!!print*,update1%Q_DG_P(:)%b
!!!!!print*,flux1(:)%H
!!!!!print*,flux1(:)%p(1)
!!!!!print*,flux1(:)%p(2)

   H_old =data%Q_DG%H
   B_old =data%Q_DG%B
   HU_old=data%Q_DG%p(1)
   HV_old=data%Q_DG%p(2)

   data_max_val(1)=maxval(H_old)
   data_max_val(2)=maxval(HU_old)
   data_max_val(3)=maxval(HV_old)

   data_min_val(1)=minval(H_old)
   data_min_val(2)=minval(HU_old)
   data_min_val(3)=minval(HV_old)

   max_neighbour(1)=max(maxval(update1%Q_DG(:)%H),maxval(update2%Q_DG(:)%H),maxval(update3%Q_DG(:)%H),data_max_val(1))
   max_neighbour(2)=max(maxval(update1%Q_DG(:)%p(1)),maxval(update2%Q_DG(:)%p(1)),maxval(update3%Q_DG(:)%p(1)),data_max_val(1))
   max_neighbour(3)=max(maxval(update1%Q_DG(:)%p(2)),maxval(update2%Q_DG(:)%p(2)),maxval(update3%Q_DG(:)%p(2)),data_max_val(1))

   min_neighbour(1)=min(minval(update1%Q_DG(:)%H),minval(update2%Q_DG(:)%H),minval(update3%Q_DG(:)%H),data_min_val(1))
   min_neighbour(2)=min(minval(update1%Q_DG(:)%p(1)),minval(update2%Q_DG(:)%p(1)),minval(update3%Q_DG(:)%p(1)),data_min_val(1))
   min_neighbour(3)=min(minval(update1%Q_DG(:)%p(2)),minval(update2%Q_DG(:)%p(2)),minval(update3%Q_DG(:)%p(2)),data_min_val(1))

   delta(1) = max(0.1_GRID_SR,max_neighbour(1)-min_neighbour(1))
   delta(2) = max(0.1_GRID_SR,max_neighbour(2)-min_neighbour(2))
   delta(3) = max(0.1_GRID_SR,max_neighbour(3)-min_neighbour(3))

   delta=delta*1.0e-3_GRID_SR
   !       !!!!!!!!!!print*,"solver"

    ! !print*,element%cell%geometry%i_plotter_type
    ! !print*,update1%flux(:)%H
    ! !print*,update1%flux(:)%p(1)
    ! !print*,update1%flux(:)%p(2)
    ! !print*,update2%flux(:)%H
    ! !print*,update2%flux(:)%p(1)
    ! !print*,update2%flux(:)%p(2)
    ! !print*,update3%flux(:)%H
    ! !print*,update3%flux(:)%p(1)
    ! !print*,update3%flux(:)%p(2)
   
   
   call dg_solver(element,update1%flux,update2%flux,update3%flux,section%r_dt)

   !       call element%cell%data_pers%convert_dg_to_fv_bathymetry()
   !call data%convert_dg_to_fv()

   ! H =data%H
   ! HU=data%HU
   ! HV=data%HV


#if defined(_SWE_DG_LIMITER_UNLIMITED)
   troubled=.false.
#elif defined(_SWE_DG_LIMITER_HEIGHT)
   troubled=.not.all(data%Q_DG%H< max_neighbour(1)+delta(1).and.data%Q_DG%H >min_neighbour(1)-delta(1))
#elif defined(_SWE_DG_LIMITER_ALL)
   troubled=&
        .not.all(data%Q_DG%H< max_neighbour(1)+delta(1).and.data%Q_DG%H >min_neighbour(1)-delta(1)) .or.&
        .not.all(data%Q_DG%p(1)< max_neighbour(2)+delta(2).and.data%Q_DG%p(1)>min_neighbour(2)-delta(2)) .or.&
        .not.all(data%Q_DG%p(2)< max_neighbour(3)+delta(3).and.data%Q_DG%p(2)>min_neighbour(3)-delta(3))
#endif
   
   !   drying=.not.all(data%H - data%B > cfg%dry_tolerance*50.0).or..not.all(data%Q_DG%H > cfg%dry_tolerance*50.0)
   drying=.not.all(data%Q_DG%H > cfg%dry_tolerance*50.0)

   if(troubled.or.drying) then
      !--if troubled or drying perform rollback--!
      data%Q_DG%H=H_old
      data%Q_DG%p(1)=HU_old
      data%Q_DG%p(2)=HV_old
      data%Q_DG%B=B_old
      data%troubled=merge(1,3,drying)
   else
      do i=1,_SWE_DG_DOFS
         wave_speed =  sqrt(g * (data%Q_DG(i)%h)) + maxval(abs(data%Q_DG(i)%p/data%Q_DG(i)%h))
         section%r_dt_new = min(section%r_dt_new,cfg%scaling*  element%transform_data%custom_data%scaling  / (wave_speed* (_SWE_DG_ORDER*4.0_GRID_SR +2.0_GRID_SR)))
         max_wave_speed = max(max_wave_speed,wave_speed)
      end do
      i_depth = element%cell%geometry%i_depth

      !--- Refinement TODO ---!
      
! #define  _SWE_DG_REFINEMENT_COAST_HEIGHT 200
! #define  _SWE_DG_REFINEMENT_COAST_DEPTH 12
! #define  _SWE_DG_REFINEMENT_WAVE_DEPTH 10
! #define  _SWE_DG_REFINEMENT_LAKE_AT_REST_DEPTH 4
      
      dQ_norm = maxval(data%Q_DG%H-H_old)
      
      if (i_depth < cfg%i_max_depth .and. dQ_norm > refinement_threshold * cfg%scaling * get_edge_size(cfg%i_max_depth)/_SWE_PATCH_ORDER) then
         element%cell%geometry%refinement = 1
         traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
      else if (i_depth > cfg%i_min_depth .and. dQ_norm < refinement_threshold * cfg%scaling * get_edge_size(cfg%i_max_depth)/(_SWE_PATCH_ORDER*8_GRID_SR)) then
         element%cell%geometry%refinement = -1
      endif
      
!       if (minval(H-data%B).le._SWE_DG_REFINEMENT_COAST_HEIGHT*cfg%dry_tolerance .and. i_depth < min(_SWE_DG_REFINEMENT_COAST_DEPTH,cfg%i_max_depth)) then
!          element%cell%geometry%refinement = 1
!          traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
!       else if(max_wave_speed > 1e-15 .and. i_depth < min(_SWE_DG_REFINEMENT_WAVE_DEPTH,cfg%i_max_depth)) then
!          element%cell%geometry%refinement = 1
!          traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
!       else if(max_wave_speed < 1e-15 .and. i_depth > max(_SWE_DG_REFINEMENT_LAKE_AT_REST_DEPTH,cfg%i_min_depth)) then
!          element%cell%geometry%refinement = -1
!       end if
      
   end if

end if   

!if cell is troubled, compute fv solution and mark all edges as troubled

if(data%troubled.ge.1) then


   call get_fv_update(update1,1)
   call get_fv_update(update2,2)
   call get_fv_update(update3,3)


   call data%convert_dg_to_fv()

   !-----Call FV patch solver----!
   call fv_patch_solver(traversal, section, element, update1, update2, update3)
   
   ! do i=1,size(element%edges)
   !    element%edges(i)%ptr%data_pers%troubled= .true.
   ! end do
   
end if
end associate

end subroutine cell_update_op_dg


subroutine compute_flux_pred(normal,QL, QR,flux_l,flux_r)

type(t_state),Dimension((_SWE_DG_ORDER+1)**2), intent(in)   :: QL
type(t_state),Dimension((_SWE_DG_ORDER+1)**2), intent(in)   :: QR
type(t_update),Dimension((_SWE_DG_ORDER+1)**2), intent(out) :: flux_l,flux_r

#if !defined(_SWE_DG_NODAL)
! type(t_update),Dimension(size(bnd_gl_node_vals,1)):: flux_l2
! type(t_update),Dimension(size(bnd_gl_node_vals,1)):: flux_r_l2
! type(t_state),Dimension(size(bnd_gl_node_vals,1)) :: QL_l2
! type(t_state),Dimension(size(bnd_gl_node_vals,1)) :: QR_l2
#endif

real(kind=GRID_SR),Dimension(1,3) :: f1r_s,f1l_s,f2r_s,f2l_s,f,ql_v,qr_v
real(kind = GRID_SR), intent(in)  :: normal(2)
real(kind = GRID_SR)	      :: normal_normed(2)
real(kind = GRID_SR)	      :: epsilon

integer ::i

#if defined(_LLF_FLUX_DG)
real(kind=GRID_SR)								:: vL, vR, alpha
real(kind=GRID_SR) :: max_wave_speedR

#endif

normal_normed=normal/NORM2(normal)

flux_l%h=0
flux_l%p(1)=0
flux_l%p(2)=0
flux_r%h=0
flux_r%p(1)=0
flux_r%p(2)=0

#if !defined(_SWE_DG_NODAL)
! QL_l2(:)%h=matmul(bnd_gl_node_vals,QL%h)
! QL_l2(:)%p(1)=matmul(bnd_gl_node_vals,QL%p(1))
! QL_l2(:)%p(2)=matmul(bnd_gl_node_vals,QL%p(2))
! QL_l2(:)%b=matmul(bnd_gl_node_vals,QL%b)

! QR_l2(:)%h=matmul(bnd_gl_node_vals,QR%h)
! QR_l2(:)%p(1)=matmul(bnd_gl_node_vals,QR%p(1))
! QR_l2(:)%p(2)=matmul(bnd_gl_node_vals,QR%p(2))
! QR_l2(:)%b=matmul(bnd_gl_node_vals,QR%b)
#endif

#if defined(_LLF_FLUX_DG)
#if defined(_SWE_DG_NODAL)
alpha=0

do i=1,_SWE_DG_ORDER+1

 if(.not.(QL(i)%h < cfg%dry_tolerance) ) then 
    vL = DOT_PRODUCT(normal_normed, QL(i)%p/QL(i)%h)
    flux_l(i)%max_wave_speed = sqrt(g * (QL(i)%h)) + abs(vL)
 else
    flux_l(i)%max_wave_speed = 0
 end if

 if(.not.(QR(i)%h < cfg%dry_tolerance) ) then 
    vR = DOT_PRODUCT(normal_normed, QR(i)%p/QR(i)%h)
    flux_r%max_wave_speed = sqrt(g * (QR(i)%h)) + abs(vR)
 else
    flux_r%wave_speedR = 0
 end if

 alpha = max(alpha,flux_l(i)%max_wave_speed, flux_r(i)%max_wave_speed)

end do

do i=1,(_SWE_DG_ORDER+1)**2
 ql_v(1,1)= QL(i)%h
 ql_v(1,2:3)= QL(i)%p

 qr_v(1,1)= QR(i)%h
 qr_v(1,2:3)= QR(i)%p

 f1r_s=f1(qr_v,1)
 f2r_s=f2(qr_v,1)
 f1l_s=f1(ql_v,1)
 f2l_s=f2(ql_v,1)
 f=0.5_GRID_SR*((f1r_s-f1l_s) * normal_normed(1) + (f2r_s-f2l_s) *normal_normed(2)+alpha*(ql_v-qr_v))
 flux_l(i)%h=f(1,1)
 flux_l(i)%p=f(1,2:3)
 flux_r(i)%h=-f(1,1)
 flux_r(i)%p=-f(1,2:3)

 ! end if
end do
#else
! alpha=0

! do i=1,size(bnd_gl_node_vals,1)

!  if(.not.(QL_L2(i)%h < cfg%dry_tolerance) ) then 
!     vL = DOT_PRODUCT(normal_normed, QL_L2(i)%p/QL_L2(i)%h)
!     flux_L2(i)%max_wave_speed = sqrt(g * (QL_L2(i)%h)) + abs(vL)
!  else
!     flux_L2(i)%max_wave_speed = 0
!  end if

!  if(.not.(QR_L2(i)%h < cfg%dry_tolerance) ) then 
!     vR = DOT_PRODUCT(normal_normed, QR_L2(i)%p/QR_L2(i)%h)
!     max_wave_speedR = sqrt(g * (QR_L2(i)%h)) + abs(vR)
!  else
!     max_wave_speedR = 0
!  end if

!  alpha = max(alpha,flux_L2(i)%max_wave_speed, max_wave_speedR)

! end do

! do i=1,size(bnd_gl_node_vals,1)
!  ql_v(1,1)= QL_L2(i)%h
!  ql_v(1,2:3)= QL_L2(i)%p

!  qr_v(1,1)= QR_L2(i)%h
!  qr_v(1,2:3)= QR_L2(i)%p

!  f1r_s=f1(qr_v,1)
!  f2r_s=f2(qr_v,1)
!  f1l_s=f1(ql_v,1)
!  f2l_s=f2(ql_v,1)

!  f=0.5_GRID_SR*((f1r_s-f1l_s) * normal_normed(1) + (f2r_s-f2l_s) *normal_normed(2)+alpha*(ql_v-qr_v))
!  flux_l2(i)%h=f(1,1)
!  flux_l2(i)%p=f(1,2:3)

!  ! end if
! end do
#endif
#elif                 (defined(_FWAVE_FLUX)|defined(_AUG_RIEMANN_FLUX)|defined(_HLLE_FLUX))

#if defined(_SWE_DG_NODAL)

flux_l%h=0
flux_l%p(1)=0
flux_l%p(2)=0
flux_r%h=0
flux_r%p(1)=0
flux_r%p(2)=0

do i=1,(_SWE_DG_ORDER+1)**2
   !!!print*,"in"   
   !!!print*,QL(i)%h
   !!!print*,QL(i)%p(1)
   !!!print*,QL(i)%p(2)
   !!!print*,QR(i)%h
   !!!print*,QR(i)%p(1)
   !!!print*,QR(i)%p(2)
   
   
   call compute_geoclaw_flux_pred(normal_normed, QL(i), QR(i), flux_l(i), flux_r(i))
   !!!print*,"out"
   !!!print*,flux_l(i)%h
   !!!print*,flux_r(i)%h

end do

#else
! do i=1,size(bnd_gl_node_vals,1)
!  call compute_geoclaw_flux_pred(normal_normed, QL_l2(i), QR_l2(i), flux_l2(i), flux_r_l2(i))
! end do
#endif

#endif

#if !defined(_SWE_DG_NODAL)
! flux%h=matmul(bnd_gl_weights,flux_l2%h)
! flux%p(1)=matmul(bnd_gl_weights,flux_l2%p(1))
! flux%p(2)=matmul(bnd_gl_weights,flux_l2%p(2))
#endif
end subroutine compute_flux_pred

subroutine dg_solver(element,update1,update2,update3,dt)
  type(t_element_base), intent(in)				:: element
  real(kind=GRID_SR), intent(in)                                :: dt
  type(t_update), dimension((_SWE_DG_ORDER+1)**2),intent(in)    :: update1, update2, update3
  
  real(kind=GRID_SR)                                            :: q(_SWE_DG_DOFS,3),&
       q_p((_SWE_DG_ORDER+1)*_SWE_DG_DOFS,3),q_temp(_SWE_DG_DOFS,3),&
       nF1(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3),&
       nF2(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3),&
       nF3(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3)
  
  real(kind=GRID_SR) :: dx
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3) :: f1_ref,f1,f2_ref,f2
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3) :: volume_flux1,volume_flux2,bnd_flux_contrib
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3) :: source,source_ref
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3) :: bnd_source_contrib,bnd_source_contrib_ref
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3) :: bnd_flux_l,bnd_flux_m,bnd_flux_r
  
#if defined(_SWE_DG_NODAL)
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1)) :: H_x,H_y,H_x_temp,H_y_temp,b
#else
  ! real(kind=GRID_SR),Dimension(size(st_gl_node_vals,1)) :: b_x,b_y
  ! real(kind=GRID_SR),Dimension(size(st_gl_node_vals,1),3) :: q_v
  ! integer :: st_nodes = size(st_gl_node_vals,1)
  ! integer :: s_nodes = size(s_der_x_gl_node_vals,1)
#endif
  
  integer :: i,j,k,inx_l,inx_m,inx_r
  real(kind=GRID_SR),Dimension(2,2) :: jacobian,jacobian_inv
  real(kind=GRID_SR),Dimension(3) :: edge_sizes     

  jacobian=ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%jacobian_normalized
  jacobian_inv=ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%jacobian_inv_normalized

  associate(data => element%cell%data_pers, Q_DG_P => element%cell%data_pers%Q_DG_P,  Q_DG => element%cell%data_pers%Q_DG)
    
    call data%get_dofs_dg(q)
    call data%get_dofs_pred(q_p)
    
    edge_sizes=element%cell%geometry%get_edge_sizes()
    dx=edge_sizes(1)*cfg%scaling

#if defined(_SWE_DG_NODAL)
    
    !---TODO: replace boundary matrices to work without zeros---!
    nF1=0
    nF2=0
    nF3=0
    do k=0,_SWE_DG_ORDER
       do i=0,_SWE_DG_ORDER
          j=i+1+k*(_SWE_DG_ORDER+1)
          inx_l=i+1+k*_SWE_DG_DOFS
          inx_m=_SWE_DG_DOFS-(_SWE_DG_ORDER-i)*(_SWE_DG_ORDER-i+1)/2+k*_SWE_DG_DOFS
          inx_r=_SWE_DG_DOFS-(_SWE_DG_ORDER-i+1)*(_SWE_DG_ORDER-i+2)/2+1+k*_SWE_DG_DOFS
          
          nF1(inx_l,1)=update1(j)%h
          nF1(inx_l,2)=update1(j)%p(1)
          nF1(inx_l,3)=update1(j)%p(2)
          
          NF2(inx_m,1)=update2(j)%h
          NF2(inx_m,2)=update2(j)%p(1)
          NF2(inx_m,3)=update2(j)%p(2)
          
          NF3(inx_r,1)=update3(j)%h
          NF3(inx_r,2)=update3(j)%p(1)
          NF3(inx_r,3)=update3(j)%p(2)
       end do
    end do
    
    bnd_flux_l=matmul(b_m_1,nF1)
    bnd_flux_m=matmul(b_m_2,nF2)
    bnd_flux_r=matmul(b_m_3,nF3)

#else

! bnd_flux_l=0
! bnd_flux_m=0
! bnd_flux_r=0

! do j=0,_SWE_DG_ORDER
!  do i=0,_SWE_DG_ORDER
!     inx_l=i+1
!     inx_m=_SWE_DG_DOFS-(_SWE_DG_ORDER-i)*(_SWE_DG_ORDER-i+1)/2
!     inx_r=_SWE_DG_DOFS-(_SWE_DG_ORDER-i+1)*(_SWE_DG_ORDER-i+2)/2+1

!     bnd_flux_l(inx_l,1)=bnd_flux_l(inx_l,1)+update1(1+j*(_SWE_DG_ORDER+1)+i)%h
!     bnd_flux_m(inx_m,1)=bnd_flux_m(inx_m,1)+update2(1+j*(_SWE_DG_ORDER+1)+i)%h * sqrt(2.0_GRID_SR)
!     bnd_flux_r(inx_r,1)=bnd_flux_r(inx_r,1)+update3(1+j*(_SWE_DG_ORDER+1)+i)%h

!     bnd_flux_l(inx_l,2)=bnd_flux_l(inx_l,2)+update1(1+j*(_SWE_DG_ORDER+1)+i)%p(1)
!     bnd_flux_m(inx_m,2)=bnd_flux_m(inx_m,2)+update2(1+j*(_SWE_DG_ORDER+1)+i)%p(1) * sqrt(2.0_GRID_SR)
!     bnd_flux_r(inx_r,2)=bnd_flux_r(inx_r,2)+update3(1+j*(_SWE_DG_ORDER+1)+1)%p(1)

!     bnd_flux_l(inx_l,3)=bnd_flux_l(inx_l,3)+update1(1+j*(_SWE_DG_ORDER+1)+1)%p(2)
!     bnd_flux_m(inx_m,3)=bnd_flux_m(inx_m,3)+update2(1+j*(_SWE_DG_ORDER+1)+1)%p(2) * sqrt(2.0_GRID_SR)

!     bnd_flux_r(inx_r,3)=bnd_flux_r(inx_r,3)+update3(1+j*(_SWE_DG_ORDER+1)+1)%p(2)
!  end do
! end do

#endif

#if defined(_SWE_DG_NODAL)

    do i=0,_SWE_DG_ORDER
       !print*,i*_SWE_DG_DOFS+1
       !print*,(i+1)*_SWE_DG_DOFS
       !print*,size(b)
       !print*,size(Q_DG(:)%b)
       b(i*_SWE_DG_DOFS+1:(i+1)*_SWE_DG_DOFS) = Q_DG(:)%b
    end do

    H_x=matmul(basis_der_st_x,Q_DG_P(:)%h+B)
    H_y=matmul(basis_der_st_y,Q_DG_P(:)%h+B)
    
    source_ref=0
    source_ref(:,2)=-matmul(s_m,(g*Q_DG_P%H*H_x))-matmul(s_k_x,(0.5_GRID_SR*g*Q_DG_P%H**2))
    source_ref(:,3)=-matmul(s_m,(g*Q_DG_P%H*H_y))-matmul(s_k_y,(0.5_GRID_SR*g*Q_DG_P%H**2))
    
    source=0
    source(:,2)=jacobian(1,1) * source_ref(:,2) + jacobian(1,2) * source_ref(:,3)
    source(:,3)=jacobian(2,1) * source_ref(:,2) + jacobian(2,2) * source_ref(:,3)
    
    bnd_source_contrib_ref=0
    bnd_source_contrib_ref(:,2)=-matmul(b_m_3,(0.5_GRID_SR*g*Q_DG_P%H**2))+matmul(b_m_2,0.5_GRID_SR*g*Q_DG_P%H**2)/sqrt(2.0_GRID_SR)
    bnd_source_contrib_ref(:,3)=-matmul(b_m_1,(0.5_GRID_SR*g*Q_DG_P%H**2))+matmul(b_m_2,0.5_GRID_SR*g*Q_DG_P%H**2)/sqrt(2.0_GRID_SR)
    
    bnd_source_contrib=0
    bnd_source_contrib(:,2)= bnd_source_contrib_ref(:,2)*jacobian(1,1)  + bnd_source_contrib_ref(:,3)*jacobian(1,2)
    bnd_source_contrib(:,3)= bnd_source_contrib_ref(:,2)*jacobian(2,1)  + bnd_source_contrib_ref(:,3)*jacobian(2,2)  
    
!    b_x=Q_DG_P%b_x
!    b_y=Q_DG_P%b_y

    if(any(q_p(:,1) <= 0))then
       print*,"Negative water height in dg solver"
       !!!print*,"PRED"
       !!!print*,q_p(:,1)
       !!!print*,"FV"
       !!!print*,element%cell%data_pers%H(:)
       stop
    end if
    
    f1_ref=flux_1(q_p,_SWE_DG_DOFS*(_SWE_DG_ORDER+1))
    f2_ref=flux_2(q_p,_SWE_DG_DOFS*(_SWE_DG_ORDER+1))
#else
    
    ! q_v=matmul(st_gl_node_vals,q_p)
    
    ! f1_hat=matmul(st_gl_weights,f1(q_v,st_nodes))
    ! f2_hat=matmul(st_gl_weights,f2(q_v,st_nodes))
    
    ! b_x=element%cell%data_pers%b_x
    ! b_y=element%cell%data_pers%b_y
    
    ! S_s =matmul(st_gl_weights,S(q_v(:,1),b_x,b_y,st_nodes))
    
    ! source=0
    
    ! do i=0,_SWE_DG_ORDER
    !    source=source+S_s(1+_SWE_DG_DOFS*i:_SWE_DG_DOFS*(1+i),:)
    ! end do
    
    ! call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f1_hat(:,1))
    ! call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f2_hat(:,1))
    
    ! call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f1_hat(:,2))
    ! call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f2_hat(:,2))
    
    ! call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f1_hat(:,3))
    ! call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f2_hat(:,3))
    
#endif                        
    
    f1=(jacobian_inv(1,1)*f1_ref+jacobian_inv(1,2)*f2_ref)
    f2=(jacobian_inv(2,1)*f1_ref+jacobian_inv(2,2)*f2_ref)
    
    volume_flux1=matmul(s_k_x,f1)
    volume_flux2=matmul(s_k_y,f2)
    
    bnd_flux_contrib=-matmul(b_m_1,f2)-matmul(b_m_3,f1)+matmul(b_m_2,f1+f2)/sqrt(2.0_GRID_SR)
    
    !!print*,"solver"
    !!print*, "orig ",q
    !!print*,"orientation",element%cell%geometry%i_plotter_type
    !!print*,"jacobian_inv",jacobian_inv
    !!print*, "dx ",dx
    !!print*, "dt ",dt
    !!print*
    !!print*, "elt flux 1 "
    !!print*,volume_flux1(:,1)
    !!print*,volume_flux1(:,2)
    !!print*,volume_flux1(:,3)
    !!print*
    !!print*, "elt flux 2 "
    !!print*,volume_flux2(:,1)
    !!print*,volume_flux2(:,2)
    !!print*,volume_flux2(:,3)
    !!print*
    !!print*, "elt flux sum"
    !!print*,volume_flux2(:,1)+volume_flux1(:,1)
    !!print*,volume_flux2(:,2)+volume_flux1(:,2)
    !!print*,volume_flux2(:,3)+volume_flux1(:,3)
    !!print*
    !!print*, "bnd_flux_contribution "
    !!print*,bnd_flux_contrib(:,1)
    !!print*,bnd_flux_contrib(:,2)
    !!print*,bnd_flux_contrib(:,3)
    !!print*
    !!print*, "bnd_flux "
    !!print*,bnd_flux_l(:,1)+bnd_flux_r(:,1)+bnd_flux_m(:,1)
    !!print*,bnd_flux_l(:,2)+bnd_flux_r(:,2)+bnd_flux_m(:,2)
    !!print*,bnd_flux_l(:,3)+bnd_flux_r(:,3)+bnd_flux_m(:,3)
    !!print*,                        
    !!print*, "volume_flux comp"
    !!print*, volume_flux1(:,1)+volume_flux2(:,1)-bnd_flux_contrib(:,1)
    !!print*, volume_flux1(:,2)+volume_flux2(:,2)-bnd_flux_contrib(:,2)
    !!print*, volume_flux1(:,3)+volume_flux2(:,3)-bnd_flux_contrib(:,3)
    !!print*
    !!print*, "source"
    !!print*,source(:,1)
    !!print*,source(:,2)
    !!print*,source(:,3)
    !!print*, "bnd source"
    !!print*, bnd_source_contrib(:,1)
    !!print*, bnd_source_contrib(:,2)
    !!print*, bnd_source_contrib(:,3)
    !!print*,"source complete"
    !!print*,-source(:,1) - bnd_source_contrib(:,1)
    !!print*,-source(:,2) - bnd_source_contrib(:,2)
    !!print*,-source(:,3) - bnd_source_contrib(:,3)
    !!print*
    !!print*,
    !!print*,"comp"
    !!print*,(bnd_flux_l(:,1)+bnd_flux_r(:,1)+bnd_flux_m(:,1)-(volume_flux1(:,1)+volume_flux2(:,1)-bnd_flux_contrib(:,1)+bnd_source_contrib(:,1)+source(:,1)))*dt/dx
    !!print*,(bnd_flux_l(:,2)+bnd_flux_r(:,2)+bnd_flux_m(:,2)-(volume_flux1(:,2)+volume_flux2(:,2)-bnd_flux_contrib(:,2)+bnd_source_contrib(:,2)+source(:,2)))*dt/dx
    !!print*,(bnd_flux_l(:,3)+bnd_flux_r(:,3)+bnd_flux_m(:,3)-(volume_flux1(:,3)+volume_flux2(:,3)-bnd_flux_contrib(:,3)+bnd_source_contrib(:,3)+source(:,3)))*dt/dx
    !!print*
    !!print*, "bnd_flux_1_dof h"
    !!print*,update1%h
    !!print*,update1%p(1)
    !!print*,update1%p(2)
    !!print*
    !!print*, "bnd_flux_2_dof"
    !!print*,update2%h
    !!print*,update2%p(1)
    !!print*,update2%p(2)
    !!print*
    !!print*, "bnd_flux_3_dof"
    !!print*,update3%h
    !!print*,update3%p(1)
    !!print*,update3%p(2)
    !!print*
    !!print*, "volume_flux1 ",volume_flux1
    !!print*,                           
    !!print*, "volume_flux2 ",volume_flux2
    !!print*,
     !!print*, "f1"
     !!print*,f1(:,1)
     !!print*,f1(:,2)
     !!print*,f1(:,3)
     !!print*,
     !!print*, "f2"
     !!print*,f2(:,1)
     !!print*,f2(:,2)
     !!print*,f2(:,3)
     !!print*,
    !!print*, "f1_ref"
    !!print*,f1_ref(:,1)
    !!print*,f1_ref(:,2)
    !!print*,f1_ref(:,3)
    !!print*,
    !!print*, "f2_ref"
    !!print*,f2_ref(:,1)
    !!print*,f2_ref(:,2)
    !!print*,f2_ref(:,3)
    !!Print*, "elt_source"
    !!print*,source_ref(:,1)
    !!print*,source_ref(:,2)
    !!print*,source_ref(:,3)
    !!print*
    !!print*,"tot Height x"
    !!print*,H_x
    !!print*,"tot Height y"
    !!print*,H_y
    
    !!print*,"orig"
    !!print*,q
    !!print*
    !!print*,"Q_DG_P"
    !!print*,q_p(:,1)
    !!print*,q_p(:,2)
    !!print*,q_p(:,3)
    !!print*,b
    !!print*
    !!print*,"pred height"
    !!print*,element%cell%data_pers%Q_DG_P(:)%h
    !!print*
    !!print*,"pred bathy"
    !!print*,b
    !!print*

    q_temp = ( bnd_flux_l &
         + bnd_flux_m &
         + bnd_flux_r &
         - volume_flux1 &
         - volume_flux2 &
         - source & 
         - bnd_source_contrib &
         + bnd_flux_contrib &
         )* dt/dx

    ! !print*,"dg temp"
    ! !print*,q_temp(:,1)
    ! !print*,q_temp(:,2)
    ! !print*,q_temp(:,3)
    ! !print*

    

    call lusolve(s_m_lu, _SWE_DG_DOFS, s_m_lu_pivot,q_temp(:,1))
    call lusolve(s_m_lu, _SWE_DG_DOFS, s_m_lu_pivot,q_temp(:,2))
    call lusolve(s_m_lu, _SWE_DG_DOFS, s_m_lu_pivot,q_temp(:,3))
    
    q=q-q_temp

    ! !print*,"dg temp lu"
    ! !print*,q_temp(:,1)
    ! !print*,q_temp(:,2)
    ! !print*,q_temp(:,3)
    ! !print*
    
    call data%set_dofs_dg(q)
    
  end associate
  
end subroutine dg_solver

subroutine compute_geoclaw_flux_pred(normal, QL, QR, fluxL, fluxR)
type(t_state), intent(in)           :: QL, QR
type(t_update), intent(out)         :: fluxL, fluxR
real(kind = GRID_SR), intent(in)    :: normal(2)
real(kind = GRID_SR)		      :: transform_matrix(2, 2)
real(kind = GRID_SR)		      :: net_updatesL(3), net_updatesR(3), max_wave_speed
real(kind = GRID_SR)                :: pL(2), pR(2), hL, hR, bL, bR

transform_matrix(1, :) = normal
transform_matrix(2, :) = [-normal(2), normal(1)]

pL = matmul(transform_matrix, QL%p)
pR = matmul(transform_matrix, QR%p)

hL = QL%h
hR = QR%h

bL = QL%b
bR = QR%b

#           if defined(_FWAVE_FLUX)
call c_bind_geoclaw_solver(GEOCLAW_FWAVE, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_AUG_RIEMANN_FLUX)
call c_bind_geoclaw_solver(GEOCLAW_AUG_RIEMANN, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_HLLE_FLUX)
call compute_updates_hlle_single(hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, net_updatesL, net_updatesR, max_wave_speed)
#           endif

fluxL%h = net_updatesL(1)
fluxL%p = matmul(net_updatesL(2:3), transform_matrix)
fluxL%max_wave_speed = max_wave_speed

fluxR%h = net_updatesR(1)
fluxR%p = matmul(net_updatesR(2:3), transform_matrix)
fluxR%max_wave_speed = max_wave_speed
end subroutine compute_geoclaw_flux_pred


subroutine get_fv_update(update,leg)
  type(num_cell_update), intent(inout) :: update
  integer,intent(in) :: leg
  real(kind=GRID_SR),Dimension(_SWE_PATCH_ORDER_SQUARE) :: H, HU, HV ,B
  integer :: i,j

  call apply_phi(update%Q_DG(:)%h+update%Q_DG(:)%b,H)
  call apply_phi(update%Q_DG(:)%p(1),HU)
  call apply_phi(update%Q_DG(:)%p(2),HV)
  call apply_phi(update%Q_DG(:)%b,B)

  
  !----FV updates need to be in reverse order----!
  select case (leg)
  case (3) !cells with id i*i+1 (left leg)
     do i=0, _SWE_PATCH_ORDER - 1
        j=_SWE_PATCH_ORDER-1-i
        update%H(i+1) = H(j*j + 1)
        update%HU(i+1)= HU(j*j + 1)
        update%HV(i+1)= HV(j*j + 1)
        update%B(i+1) = B(j*j + 1)
     end do
  case (2) ! hypotenuse
     do i=1, _SWE_PATCH_ORDER
        j=_SWE_PATCH_ORDER+1-i
        update%H(i) = H ((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
        update%HU(i)= HU((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
        update%HV(i)= HV((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
        update%B(i) = B ((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
     end do
  case (1) !cells with id i*i (right leg)
     do i=1, _SWE_PATCH_ORDER
!        j=_SWE_PATCH_ORDER+1-i
        update%H(i)  = H(i*i)
        update%HU(i) = HU(i*i)
        update%HV(i) = HV(i*i)
        update%B(i)  = B(i*i)
     end do
  end select


  where (update%H < update%B + cfg%dry_tolerance) 
     update%H = update%B
     update%HU = 0.0_GRID_SR
     update%HV = 0.0_GRID_SR
  end where

end subroutine get_fv_update

subroutine fv_patch_solver(traversal, section, element, update1, update2, update3)
            class(t_traversal), intent(inout)		                            :: traversal
            type(t_grid_section), intent(inout)				            :: section
            type(t_element_base), intent(inout)				            :: element
            type(num_cell_update), intent(inout)			            :: update1, update2, update3
            integer                                                                 :: i, j, ind
            type(num_cell_update)                                                   :: tmp !> ghost cells in correct order 
            real(kind = GRID_SR)                                                    :: volume, edge_lengths(3), maxWaveSpeed, dQ_max_norm, dt_div_volume

            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_ORDER_SQUARE)                :: dQ_H, dQ_HU, dQ_HV !> deltaQ, used to compute cell updates
            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)           :: hL, huL, hvL, bL
            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)           :: hR, huR, hvR, bR
            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)           :: upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR
            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE,2,2)       :: transf
            real(kind= GRID_SR) :: maxWaveSpeeds (_SWE_PATCH_ORDER_SQUARE)

#if !defined (_SWE_USE_PATCH_SOLVER)
            type(t_state), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE)      :: edges_a, edges_b
            type(t_update)                                              :: update_a, update_b
            real(kind = GRID_SR), dimension(2,3)                        :: normals
#endif
#if defined(_SWE_DG)
            type(t_update),  DIMENSION((_SWE_DG_ORDER+1)**2)		:: flux1, flux2, flux3
            type(t_state), dimension((_SWE_DG_ORDER+1)**2)              :: edge_l, edge_m,edge_r
            real(kind= GRID_SR):: dt,dx,delta(3),max_neighbour(3),min_neighbour(3),data_max_val(3),data_min_val(3)
            integer:: indx,indx_mir,k
            real (kind = GRID_SR), dimension (_SWE_PATCH_ORDER_SQUARE):: H_old, HU_old, HV_old
            logical :: drying,troubled
#endif

#if !defined(_SWE_USE_PATCH_SOLVER)
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
#endif

            ! init some variables
            dQ_H = 0.0_GRID_SR
            dQ_HU = 0.0_GRID_SR
            dQ_HV = 0.0_GRID_SR
            maxWaveSpeed = 0.0_GRID_SR

            volume = cfg%scaling * cfg%scaling * element%cell%geometry%get_volume() / (_SWE_PATCH_ORDER_SQUARE)
            dt_div_volume = section%r_dt / volume
            edge_lengths = cfg%scaling * element%cell%geometry%get_edge_sizes() / _SWE_PATCH_ORDER

            associate(data => element%cell%data_pers, geom => SWE_PATCH_geometry)

               ! copy cell values to arrays edges_a and edges_b
               ! obs: cells with id > number of cells are actually ghost cells and come from edges "updates"
               ! see t_SWE_PATCH_geometry for an explanation about ghost cell ordering
               do i=1, _SWE_PATCH_NUM_EDGES, _SWE_PATCH_SOLVER_CHUNK_SIZE ! i -> init of chunk
               ! if this is the last chunk and it is not a full chunk, 
               ! it is necessary to set everything to 0 to avoid using last iteration values
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

                  do j=1, _SWE_PATCH_SOLVER_CHUNK_SIZE ! j -> position inside chunk
                     ind = i + j - 1 ! actual index

                     ! don't go outside array limits
                     if (ind > _SWE_PATCH_NUM_EDGES) then
                        exit
                     end if

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

                     !copy transformation matrices
                     transf(j,:,:) = geom%transform_matrices(geom%edges_orientation(ind),:,:,element%cell%geometry%i_plotter_type)
                  end do

                  ! compute net_updates -> solve Riemann problems within chunk
#                       if defined (_SWE_USE_PATCH_SOLVER)
#                       if defined(_FWAVE_FLUX) || defined(_AUG_RIEMANN_FLUX)
                  call compute_updates_simd(transf, hL, huL, hvL, bL, hR, huR, hvR, bR, upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR, maxWaveSpeed)
#                       elif defined(_HLLE_FLUX)
                  call compute_updates_hlle_simd(transf, hL, huL, hvL, bL, hR, huR, hvR, bR, upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR, maxWaveSpeed)
#                       else
                  ! this should never happen -> SCons rules should avoid this before compiling
#                       error "No valid SWE solver for patches/simd implementation has been defined!"
#                       endif
#                       else                    
                  ! using original geoclaw solver
                  maxWaveSpeeds=0
                  do j=1, _SWE_PATCH_SOLVER_CHUNK_SIZE
                     ind = i + j - 1 ! actual index

                     ! don't go outside array limits
                     if (ind > _SWE_PATCH_NUM_EDGES) then
                        exit
                     end if

                     edges_a(j)%h = hL(j)
                     edges_a(j)%p(1) = huL(j)
                     edges_a(j)%p(2) = hvL(j)
                     edges_a(j)%b = bL(j)
                     edges_b(j)%h = hR(j)
                     edges_b(j)%p(1) = huR(j)
                     edges_b(j)%p(2) = hvR(j)
                     edges_b(j)%b = bR(j)
                     call compute_geoclaw_flux(normals(:,geom%edges_orientation(ind)), edges_a(j), edges_b(j), update_a, update_b)
                     upd_hL(j) = update_a%h
                     upd_huL(j) = update_a%p(1)
                     upd_hvL(j) = update_a%p(2)
                     upd_hR(j) = update_b%h
                     upd_huR(j) = update_b%p(1)
                     upd_hvR(j) = update_b%p(2)
                     maxWaveSpeeds(j)=update_a%max_wave_speed
                  end do
#                       endif

                  ! compute dQ
                  do j=1, _SWE_PATCH_SOLVER_CHUNK_SIZE
                     ind = i + j - 1 ! actual index

                     ! don't go outside array limits
                     if (ind > _SWE_PATCH_NUM_EDGES) then
                        exit
                     end if

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

               if (element%cell%geometry%i_depth < cfg%i_max_depth .and. dQ_max_norm > 5.0_GRID_SR * cfg%scaling * get_edge_size(cfg%i_max_depth) / _SWE_PATCH_ORDER ) then
                  element%cell%geometry%refinement = 1
                  traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
               else if (element%cell%geometry%i_depth > cfg%i_min_depth .and. dQ_max_norm < 5.0_GRID_SR * cfg%scaling * get_edge_size(cfg%i_max_depth) / (_SWE_PATCH_ORDER * 8.0_SR) ) then
                  element%cell%geometry%refinement = -1
               endif

               dQ_H = dQ_H * (-dt_div_volume)
               dQ_HU = dQ_HU * (-dt_div_volume)
               dQ_HV = dQ_HV * (-dt_div_volume)


               ! if the water level falls below the dry tolerance, set water level to 0 and velocity to 0
               where (data%H < data%B + cfg%dry_tolerance) 
                  data%H = data%B
                  data%HU = 0.0_GRID_SR
                  data%HV = 0.0_GRID_SR
               end where


               ! if land is flooded, init water height to dry tolerance and
               ! velocity to zero
               where (data%H < data%B + cfg%dry_tolerance .and. dQ_H > 0.0_GRID_SR)
                  data%H = data%B + cfg%dry_tolerance
                  data%HU = 0.0_GRID_SR
                  data%HV = 0.0_GRID_SR
               end where

               ! update unknowns
               data%H = data%H + dQ_H
               data%HU = data%HU + dQ_HU
               data%HV = data%HV + dQ_HV
               
               ! if the water level falls below the dry tolerance, set water level to 0 and velocity to 0
               where (data%H < data%B + cfg%dry_tolerance) 
                  data%H = data%B
                  data%HU = 0.0_GRID_SR
                  data%HV = 0.0_GRID_SR
               end where

!               !print*,data%H
!               !print*,data%HU
!               !print*,data%HV
               
               maxWaveSpeed=maxval(maxWaveSpeeds)
               section%r_dt_new = min(section%r_dt_new, volume / (edge_lengths(2) * maxWaveSpeed) )
               maxWaveSpeed=0


               call data%convert_fv_to_dg()
               if(all(data%H - data%B > cfg%dry_tolerance*50.0)) then
                  if(all(data%Q_DG%H > cfg%dry_tolerance*50.0)) then
                     data%troubled = -data%troubled
                  else
                     data%troubled = 1
                  end if
               else
                  data%troubled = 1
               end if
               
            !consider new dg cells in r_dt
               if(data%troubled .le. 0) then
                  do i=1,_SWE_DG_DOFS
                     maxWaveSpeed =  sqrt(g * (data%Q_DG(i)%h)) + maxval(abs(data%Q_DG(i)%p/data%Q_DG(i)%h))
                     section%r_dt_new = min(section%r_dt_new,cfg%scaling*  element%transform_data%custom_data%scaling  / (maxWaveSpeed* (_SWE_DG_ORDER*4.0_GRID_SR +2.0_GRID_SR)))
                     maxWaveSpeed = 0.0
                  end do
               end if
               
          end associate

        end subroutine

        
END MODULE SWE_DG_solver
!_SWE_DG
#endif
