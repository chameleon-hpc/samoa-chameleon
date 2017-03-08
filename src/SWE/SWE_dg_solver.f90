! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2017 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE
! Author: Leonhard Rannabauer rannabau (at) in.tum.de

#include "Compilation_control.f90"

MODULE SWE_DG_solver
  use SFC_edge_traversal

  use Samoa_swe
  use c_bind_riemannsolvers

  use SWE_DG_Matrices
  use SWE_dg_predictor

#               if defined(_HLLE_FLUX)
  use SWE_HLLE
#               endif


  implicit none
  type num_traversal_data
     integer (kind = GRID_DI)			:: i_refinements_issued
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

#		define _GT_SKELETON_OP			skeleton_op_dg
#		define _GT_BND_SKELETON_OP		bnd_skeleton_op_dg
#		define _GT_ELEMENT_OP                   element_op_dg

#		define _GT_CELL_UPDATE_OP		cell_update_op_dg
#		define _GT_CELL_TO_EDGE_OP		cell_to_edge_op_dg

#               define _GT_EDGE_FIRST_TOUCH_OP          edge_first_touch_op_dg

  public cell_to_edge_op_dg

#		include "SFC_generic_traversal_ringbuffer.f90"

  subroutine edge_first_touch_scalar_op_dg(traversal, section, edge)
    type(t_swe_dg_timestep_traversal)  ::traversal
    type(t_grid_section), intent(inout)::section
    type(t_edge_data), intent(inout)   ::edge
    edge%data_pers%troubled=.false.
  end subroutine edge_first_touch_scalar_op_dg

  subroutine edge_first_touch_array_op_dg(traversal, section, edge)
    type(t_swe_dg_timestep_traversal)  ::traversal
    type(t_grid_section), intent(inout)::section
    type(t_edge_data), intent(inout)   ::edge(:)
    edge%data_pers%troubled=.false.
  end subroutine edge_first_touch_array_op_dg


  subroutine pre_traversal_grid_op_dg(traversal, grid)
    type(t_swe_dg_timestep_traversal), intent(inout)		:: traversal
    type(t_grid), intent(inout)							    :: grid
  end subroutine pre_traversal_grid_op_dg


  subroutine pre_traversal_op_dg(traversal, section)
    type(t_swe_dg_timestep_traversal), intent(inout)				:: traversal
    type(t_grid_section), intent(inout)							:: section
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
    integer                                             :: edge_type !1=left, 2=hypotenuse, 3=right
    integer :: i,k,indx,indx_mir
    type(num_cell_rep) ::rep_fv

#if defined(_SWE_DG)		

    associate(cell_edge => ref_plotter_data(-abs(element%cell%geometry%i_plotter_type))%edges, normal => edge%transform_data%normal)
    do i=1,3
       if ((normal(1) == cell_edge(i)%normal(1) .and. normal(2) == cell_edge(i)%normal(2))   &
            .or. (normal(1) == -cell_edge(i)%normal(1) .and. normal(2) == -cell_edge(i)%normal(2))) then
          edge_type = i
       end if
    end do
  end associate

  if(element%cell%data_pers%troubled.le.0) then
     do k=0,_SWE_DG_ORDER
        do i=0,_SWE_DG_ORDER
           indx_mir=i+1+k*(_SWE_DG_ORDER+1)
           select case (edge_type)
           case(1) !left
              indx=1+i+k*_SWE_DG_DOFS
              rep%Q_DG_P(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
              rep%Q_DG_P(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
              rep%Q_DG_P(indx_mir)%b = element%cell%data_pers%Q_DG_P(indx)%b
           case(2) !mid
              indx=_SWE_DG_DOFS-(i)*(i+1)/2+k*_SWE_DG_DOFS
              rep%Q_DG_P(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
              rep%Q_DG_P(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
              rep%Q_DG_P(indx_mir)%b = element%cell%data_pers%Q_DG_P(indx)%b
           case(3) !right
              indx=_SWE_DG_DOFS-(_SWE_DG_ORDER-i+1)*(_SWE_DG_ORDER-i+2)/2 +1+k*_SWE_DG_DOFS
              rep%Q_DG_P(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
              rep%Q_DG_P(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
              rep%Q_DG_P(indx_mir)%b = element%cell%data_pers%Q_DG_P(indx)%b
           end select
        end do
     end do
  end if


  !     call element%cell%data_pers%convert_dg_to_fv()

  associate(H => element%cell%data_pers%H, HU => element%cell%data_pers%HU, HV => element%cell%data_pers%HV, B => element%cell%data_pers%B)
  select case (edge_type)
  case (1) !cells with id i*i+1 (left leg)
     do i=0, _SWE_PATCH_ORDER - 1
        rep%H(i+1) = H(i*i + 1)
        rep%HU(i+1) = HU(i*i + 1)
        rep%HV(i+1) = HV(i*i + 1)
        rep%B(i+1) = B(i*i + 1)
     end do
  case (2) ! hypotenuse
     do i=1, _SWE_PATCH_ORDER
        rep%H(i) = H((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
        rep%HU(i) = HU((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
        rep%HV(i) = HV((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
        rep%B(i) = B((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
     end do
  case (3) !cells with id i*i (right leg)
     do i=1, _SWE_PATCH_ORDER
        rep%H(_SWE_PATCH_ORDER + 1 - i) = H(i*i)
        rep%HU(_SWE_PATCH_ORDER + 1 - i) = HU(i*i)
        rep%HV(_SWE_PATCH_ORDER + 1 - i) = HV(i*i)
        rep%B(_SWE_PATCH_ORDER + 1 - i) = B(i*i)
     end do
  end select

end associate


rep%troubled=element%cell%data_pers%troubled_old  

rep%max_val(1)=maxval(element%cell%data_pers%H)
rep%max_val(2)=maxval(element%cell%data_pers%HU)
rep%max_val(3)=maxval(element%cell%data_pers%HV)

rep%min_val(1)=minval(element%cell%data_pers%H)
rep%min_val(2)=minval(element%cell%data_pers%HU)
rep%min_val(3)=minval(element%cell%data_pers%HV)

#endif				

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
type(t_grid_section), intent(in)				:: grid
type(t_edge_data), intent(in)				:: edge
type(num_cell_rep), intent(in)				:: rep1, rep2
type(num_cell_update), intent(out)				:: update1, update2

update1%Q_DG_P(:)%h   =rep2%Q_DG_P(:)%h
update1%Q_DG_P(:)%p(1)=rep2%Q_DG_P(:)%p(1)
update1%Q_DG_P(:)%p(2)=rep2%Q_DG_P(:)%p(2)
update1%Q_DG_P(:)%b   =rep2%Q_DG_P(:)%b
update1%troubled      =rep2%troubled
update1%max_val       =rep2%max_val
update1%min_val       =rep2%min_val

update2%Q_DG_P(:)%h   =rep1%Q_DG_P(:)%h
update2%Q_DG_P(:)%p(1)=rep1%Q_DG_P(:)%p(1)
update2%Q_DG_P(:)%p(2)=rep1%Q_DG_P(:)%p(2)
update2%Q_DG_P(:)%b   =rep1%Q_DG_P(:)%b
update2%troubled      =rep1%troubled
update2%max_val       =rep1%max_val
update2%min_val       =rep1%min_val

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
real(kind=GRID_SR)                             :: normal_normed(2)
real(kind=GRID_SR)                             :: length_flux
integer                                        :: i

normal_normed=(edge%transform_data%normal)/NORM2(edge%transform_data%normal)

update%Q_DG_P(:)%h = rep%Q_DG_P(:)%h
update%Q_DG_P(:)%b = rep%Q_DG_P(:)%b

update%max_val=rep%max_val
update%min_val=rep%min_val

update%troubled=rep%troubled

do i=1,(_SWE_DG_ORDER+1)**2
   !                          reflecting boundary
   length_flux = dot_product(rep%Q_DG_P(i)%p, normal_normed)

   update%Q_DG_P(i)%p(1) = rep%Q_DG_P(i)%p(1)-2.0_GRID_SR*length_flux*normal_normed(1)
   update%Q_DG_P(i)%p(2) = rep%Q_DG_P(i)%p(2)-2.0_GRID_SR*length_flux*normal_normed(2)

   !simple outflowing boundary
   !        update%Q_DG_P(i)%p(1) =  rep%Q_DG_P(i)%p(1)
   !        update%Q_DG_P(i)%p(2) =  rep%Q_DG_P(i)%p(2)
   !                            end if

   !zero vel bnd
   !                          update%Q_DG_P(i)%p(1) =  0
   !                          update%Q_DG_P(i)%p(2) =  0
end do

end subroutine bnd_skeleton_scalar_op_dg

!cup_time
subroutine cell_update_op_dg(traversal, section, element, update1, update2, update3)
type(t_swe_dg_timestep_traversal), intent(inout)		:: traversal
type(t_grid_section), intent(inout)			:: section
type(t_element_base), intent(inout)			:: element
type(num_cell_update), intent(inout)			:: update1, update2, update3
type(num_cell_update) :: tmp
type(t_update),  DIMENSION((_SWE_DG_ORDER+1)**2)		:: flux1, flux2, flux3
type(t_state), dimension((_SWE_DG_ORDER+1)**2)             :: edge_l, edge_m,edge_r
real(kind= GRID_SR) :: dt,dx,delta(3),max_neighbour(3),min_neighbour(3),data_max_val(3),data_min_val(3)
integer :: indx,indx_mir,k,i
real (kind=GRID_SR) :: max_wave_speed, wave_speed=0
real (kind = GRID_SR), dimension (_SWE_PATCH_ORDER_SQUARE) :: H_old, HU_old, HV_old, H, HU, HV, B_old
logical :: drying,troubled,neighbours_troubled
integer :: i_depth


if (element%cell%geometry%i_plotter_type > 0) then ! if orientation = forward, reverse updates
   tmp=update1
   update1=update3
   update3=tmp
end if

neighbours_troubled=(update1%troubled.ge.1).or.(update2%troubled.ge.1).or.(update3%troubled.ge.1)

associate(data => element%cell%data_pers)

if(neighbours_troubled) then
   data%troubled=merge(data%troubled,2,data%troubled.ge.1)
   !        call element%cell%data_pers%convert_dg_to_fv_bathymetry()
end if

if(data%troubled.le.0) then
   data%troubled=0
   do k=0,_SWE_DG_ORDER
      do i=0,_SWE_DG_ORDER
         indx_mir=i+1+k*(_SWE_DG_ORDER+1)

         indx=i+1+k*_SWE_DG_DOFS
         edge_l(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
         edge_l(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
         edge_l(indx_mir)%b = element%cell%data_pers%Q_DG_P(indx)%b

         indx=_SWE_DG_DOFS-(_SWE_DG_ORDER-i)*(_SWE_DG_ORDER-i+1)/2+k*_SWE_DG_DOFS
         edge_m(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
         edge_m(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
         edge_m(indx_mir)%b = element%cell%data_pers%Q_DG_P(indx)%b

         indx=_SWE_DG_DOFS-(_SWE_DG_ORDER-i+1)*(_SWE_DG_ORDER-i+2)/2+1+k*_SWE_DG_DOFS
         edge_r(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
         edge_r(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
         edge_r(indx_mir)%b = element%cell%data_pers%Q_DG_P(indx)%b
      end do
   end do

   call compute_flux_pred(ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%edges(3)%normal,edge_l,update1%Q_DG_P,flux1)
   call compute_flux_pred(ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%edges(2)%normal,edge_m,update2%Q_DG_P,flux2)
   call compute_flux_pred(ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%edges(1)%normal,edge_r,update3%Q_DG_P,flux3)

!!!print*,"updatesl"
!!!print*,edge_l(:)%H
!!!print*,edge_l(:)%p(1)
!!!print*,edge_l(:)%p(2)
!!!print*,edge_l(:)%b
!!!print*,update1%Q_DG_P(:)%H
!!!print*,update1%Q_DG_P(:)%p(1)
!!!print*,update1%Q_DG_P(:)%p(2)
!!!print*,update1%Q_DG_P(:)%b
!!!print*,flux1(:)%H
!!!print*,flux1(:)%p(1)
!!!print*,flux1(:)%p(2)



   H_old =data%H
   B_old =data%B
   HU_old=data%HU
   HV_old=data%HV

   data_max_val(1)=maxval(H_old)
   data_max_val(2)=maxval(HU_old)
   data_max_val(3)=maxval(HV_old)

   data_min_val(1)=minval(H_old)
   data_min_val(2)=minval(HU_old)
   data_min_val(3)=minval(HV_old)

   max_neighbour=max(update1%max_val,update2%max_val,update3%max_val,data_max_val)
   min_neighbour=min(update1%min_val,update2%min_val,update3%min_val,data_min_val)

   delta(1) = max(0.1_GRID_SR,max_neighbour(1)-min_neighbour(1))
   delta(2) = max(0.1_GRID_SR,max_neighbour(2)-min_neighbour(2))
   delta(3) = max(0.1_GRID_SR,max_neighbour(3)-min_neighbour(3))

   delta=delta*1.0e-3_GRID_SR
   !       !!!!!!!!print*,"solver"
   call dg_solver(element,flux1,flux2,flux3,section%r_dt)

   !       call element%cell%data_pers%convert_dg_to_fv_bathymetry()
   call data%convert_dg_to_fv()

   H =data%H
   HU=data%HU
   HV=data%HV


#if defined(_SWE_DG_LIMITER_UNLIMITED)
   troubled=.false.
#elif defined(_SWE_DG_LIMITER_HEIGHT)
   troubled=.not.all(data%H< max_neighbour(1)+delta(1).and.data%H >min_neighbour(1)-delta(1))
#elif defined(_SWE_DG_LIMITER_ALL)
   troubled=&
        .not.all(data%H< max_neighbour(1)+delta(1).and.data%H >min_neighbour(1)-delta(1)) .or.&
        .not.all(data%HU< max_neighbour(2)+delta(2).and.data%HU>min_neighbour(2)-delta(2)) .or.&
        .not.all(data%HV< max_neighbour(3)+delta(3).and.data%HV>min_neighbour(3)-delta(3))
#endif       
   drying=.not.all(data%H - data%B > cfg%dry_tolerance*50.0).or..not.all(data%Q_DG%H > cfg%dry_tolerance*50.0)

   if(troubled.or.drying) then
      data%H=H_old
      data%HU=HU_old
      data%HV=HV_old
      data%B=B_old
      data%troubled=merge(1,3,drying)
   else

      do i=1,_SWE_DG_DOFS
         wave_speed =  sqrt(g * (data%Q_DG(i)%h)) + maxval(abs(data%Q_DG(i)%p/data%Q_DG(i)%h))
         section%r_dt_new = min(section%r_dt_new,cfg%scaling*  element%transform_data%custom_data%scaling  / (wave_speed* (_SWE_DG_ORDER*4.0_GRID_SR +2.0_GRID_SR)))
         max_wave_speed = max(max_wave_speed,wave_speed)
      end do
      i_depth = element%cell%geometry%i_depth

#define  _SWE_DG_REFINEMENT_COAST_HEIGHT 200
#define  _SWE_DG_REFINEMENT_COAST_DEPTH 12
#define  _SWE_DG_REFINEMENT_WAVE_DEPTH 10
#define  _SWE_DG_REFINEMENT_LAKE_AT_REST_DEPTH 4

      if (minval(H-data%B).le._SWE_DG_REFINEMENT_COAST_HEIGHT*cfg%dry_tolerance .and. i_depth < min(_SWE_DG_REFINEMENT_COAST_DEPTH,cfg%i_max_depth)) then
         element%cell%geometry%refinement = 1
         traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
      else if(max_wave_speed > 1e-15 .and. i_depth < min(_SWE_DG_REFINEMENT_WAVE_DEPTH,cfg%i_max_depth)) then
         element%cell%geometry%refinement = 1
         traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
      else if(max_wave_speed < 1e-15 .and. i_depth > max(_SWE_DG_REFINEMENT_LAKE_AT_REST_DEPTH,cfg%i_min_depth)) then
         element%cell%geometry%refinement = -1
      end if
   end if
end if

!if cell is troubled mark all edges as troubled
if(data%troubled.ge.1) then
   do i=1,size(element%edges)
      element%edges(i)%ptr%data_pers%troubled= .true.
   end do
end if
end associate

end subroutine cell_update_op_dg


subroutine compute_flux_pred(normal,QL, QR,flux)

type(t_state),Dimension((_SWE_DG_ORDER+1)**2), intent(inout)	:: QL
type(t_state),Dimension((_SWE_DG_ORDER+1)**2), intent(inout)   :: QR
type(t_update),Dimension((_SWE_DG_ORDER+1)**2), intent(out) :: flux
type(t_update),Dimension((_SWE_DG_ORDER+1)**2) 		:: flux_r
#if !defined(_SWE_DG_NODAL)
type(t_update),Dimension(size(bnd_gl_node_vals,1)):: flux_l2
type(t_update),Dimension(size(bnd_gl_node_vals,1)):: flux_r_l2
type(t_state),Dimension(size(bnd_gl_node_vals,1)) :: QL_l2
type(t_state),Dimension(size(bnd_gl_node_vals,1)) :: QR_l2
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

flux%h=0
flux%p(1)=0
flux%p(2)=0
flux_r%h=0
flux_r%p(1)=0
flux_r%p(2)=0

#if !defined(_SWE_DG_NODAL)
QL_l2(:)%h=matmul(bnd_gl_node_vals,QL%h)
QL_l2(:)%p(1)=matmul(bnd_gl_node_vals,QL%p(1))
QL_l2(:)%p(2)=matmul(bnd_gl_node_vals,QL%p(2))
QL_l2(:)%b=matmul(bnd_gl_node_vals,QL%b)

QR_l2(:)%h=matmul(bnd_gl_node_vals,QR%h)
QR_l2(:)%p(1)=matmul(bnd_gl_node_vals,QR%p(1))
QR_l2(:)%p(2)=matmul(bnd_gl_node_vals,QR%p(2))
QR_l2(:)%b=matmul(bnd_gl_node_vals,QR%b)
#endif

#if defined(_LLF_FLUX_DG)
#if defined(_SWE_DG_NODAL)
alpha=0

do i=1,_SWE_DG_ORDER+1

 if(.not.(QL(i)%h < cfg%dry_tolerance) ) then 
    vL = DOT_PRODUCT(normal_normed, QL(i)%p/QL(i)%h)
    flux(i)%max_wave_speed = sqrt(g * (QL(i)%h)) + abs(vL)
 else
    flux(i)%max_wave_speed = 0
 end if

 if(.not.(QR(i)%h < cfg%dry_tolerance) ) then 
    vR = DOT_PRODUCT(normal_normed, QR(i)%p/QR(i)%h)
    max_wave_speedR = sqrt(g * (QR(i)%h)) + abs(vR)
 else
    max_wave_speedR = 0
 end if

 alpha = max(alpha,flux(i)%max_wave_speed, max_wave_speedR)

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
 flux(i)%h=f(1,1)
 flux(i)%p=f(1,2:3)

 ! end if
end do
#else
alpha=0

do i=1,size(bnd_gl_node_vals,1)

 if(.not.(QL_L2(i)%h < cfg%dry_tolerance) ) then 
    vL = DOT_PRODUCT(normal_normed, QL_L2(i)%p/QL_L2(i)%h)
    flux_L2(i)%max_wave_speed = sqrt(g * (QL_L2(i)%h)) + abs(vL)
 else
    flux_L2(i)%max_wave_speed = 0
 end if

 if(.not.(QR_L2(i)%h < cfg%dry_tolerance) ) then 
    vR = DOT_PRODUCT(normal_normed, QR_L2(i)%p/QR_L2(i)%h)
    max_wave_speedR = sqrt(g * (QR_L2(i)%h)) + abs(vR)
 else
    max_wave_speedR = 0
 end if

 alpha = max(alpha,flux_L2(i)%max_wave_speed, max_wave_speedR)

end do

do i=1,size(bnd_gl_node_vals,1)
 ql_v(1,1)= QL_L2(i)%h
 ql_v(1,2:3)= QL_L2(i)%p

 qr_v(1,1)= QR_L2(i)%h
 qr_v(1,2:3)= QR_L2(i)%p

 f1r_s=f1(qr_v,1)
 f2r_s=f2(qr_v,1)
 f1l_s=f1(ql_v,1)
 f2l_s=f2(ql_v,1)

 f=0.5_GRID_SR*((f1r_s-f1l_s) * normal_normed(1) + (f2r_s-f2l_s) *normal_normed(2)+alpha*(ql_v-qr_v))
 flux_l2(i)%h=f(1,1)
 flux_l2(i)%p=f(1,2:3)

 ! end if
end do
#endif
#elif                 (defined(_FWAVE_FLUX)|defined(_AUG_RIEMANN_FLUX)|defined(_HLLE_FLUX))

#if defined(_SWE_DG_NODAL)

flux%h=0
flux%p(1)=0
flux%p(2)=0
flux_r%h=0
flux_r%p(1)=0
flux_r%p(2)=0

do i=1,(_SWE_DG_ORDER+1)**2
   !    epsilon=1.0e-14
   !    if(abs(QL(i)%p(1))<epsilon)then
   !       QL(i)%p(1)=0.0_GRID_SR
   !    end if
   !    if(abs(QR(i)%p(1))<epsilon)then
   !       QR(i)%p(1)=0.0_GRID_SR
   !    end if
   !    if(abs(QL(i)%p(2))<epsilon)then
   !       QL(i)%p(2)=0.0_GRID_SR
   !    end if
   !    if(abs(QR(i)%p(2))<epsilon)then
   !       QR(i)%p(2)=0.0_GRID_SR
   !    end if
   !    if(QL(i)%b.ne.0) then
   !       if(abs(QL(i)%b-QR(i)%b)/QL(i)%b < epsilon)then
   !          QL(i)%b=QR(i)%b
   !       end if
   !    else if(abs((QL(i)%b-QR(i)%b)) < epsilon)then
   !          QL(i)%b=QR(i)%b
   !    end if

   !    if(QL(i)%h.ne.0) then
   !       if(abs(QL(i)%h-QR(i)%h)/QL(i)%b < epsilon)then
   !          QL(i)%h=QR(i)%h
   !       end if
   !    else if(abs(QL(i)%h-QR(i)%h) < epsilon)then
   !       QL(i)%h=QR(i)%h
   !    end if



 call compute_geoclaw_flux_pred(normal_normed, QL(i), QR(i), flux(i), flux_r(i))
end do

#else
do i=1,size(bnd_gl_node_vals,1)
 call compute_geoclaw_flux_pred(normal_normed, QL_l2(i), QR_l2(i), flux_l2(i), flux_r_l2(i))
end do
#endif

#endif

#if !defined(_SWE_DG_NODAL)

flux%h=matmul(bnd_gl_weights,flux_l2%h)
flux%p(1)=matmul(bnd_gl_weights,flux_l2%p(1))
flux%p(2)=matmul(bnd_gl_weights,flux_l2%p(2))
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
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1)) :: H_x,H_y,H_x_temp,H_y_temp,b_x,b_y
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

  associate(Q_DG => element%cell%data_pers,Q_DG_P => element%cell%data_pers%Q_DG_P)
    
    call Q_DG%get_dofs_dg(q)
    call Q_DG%get_dofs_pred(q_p)
    
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

    H_x=matmul(basis_der_st_x,Q_DG_P(:)%h+Q_DG_P(:)%B)
    H_y=matmul(basis_der_st_y,Q_DG_P(:)%h+Q_DG_P(:)%B)
    
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
    
    b_x=Q_DG_P%b_x
    b_y=Q_DG_P%b_y
      
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
    
    q_temp = ( bnd_flux_l &
         + bnd_flux_m &
         + bnd_flux_r &
         - volume_flux1 &
         - volume_flux2 &
         - source & 
         - bnd_source_contrib &
         + bnd_flux_contrib &
         )* dt/dx  
    
    call lusolve(s_m_lu, _SWE_DG_DOFS, s_m_lu_pivot,q_temp(:,1))
    call lusolve(s_m_lu, _SWE_DG_DOFS, s_m_lu_pivot,q_temp(:,2))
    call lusolve(s_m_lu, _SWE_DG_DOFS, s_m_lu_pivot,q_temp(:,3))
    
    q=q-q_temp
    
    call Q_DG%set_dofs_dg(q)
    
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

END MODULE SWE_DG_solver
