! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2017 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE
! Author: Leonhard Rannabauer rannabau (at) in.tum.de

#if defined(_SWE_DG)
#include "Compilation_control.f90"

MODULE SWE_DG_solver
  use SFC_edge_traversal
  use SWE_euler_timestep
  use Samoa_swe
  use c_bind_riemannsolvers


#       if defined(CHAMELEON)
             use Chameleon_lib
#            define _GT_USE_CHAMELEON
#       endif
#       if defined(CHAMELEON_CALL)
#            define _GT_USE_CHAMELEON_CALL
#       endif
  use SWE_DG_Matrices
  use SWE_DG_Predictor
  use SWE_initialize_bathymetry
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
    integer :: i,j,k,indx,indx_mir
    integer(kind=BYTE) :: orientation
    type(num_cell_rep) ::rep_fv

#if defined (_DEBUG)
    rep%debug_flag = element%cell%data_pers%debug_flag
#endif    

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
    !---- Project time averaged predictor and fluxes ----!
    do i=0,_SWE_DG_ORDER
       select case (edge_type)
       case(-1 ,1) !right
          indx=_SWE_DG_DOFS-(_SWE_DG_ORDER-i+1)*(_SWE_DG_ORDER-i+2)/2 +1
       case(-2,2) !mid
          indx=_SWE_DG_DOFS-(_SWE_DG_ORDER-i)*(_SWE_DG_ORDER-i+1)/2
       case(-3 ,3) !left
          indx=1+i
       case default
          stop
       end select
       rep%Q(i+1)%h                        = element%cell%data_pers%QP(indx,1)
       rep%Q(i+1)%p(1)                     = element%cell%data_pers%QP(indx,2)
       rep%Q(i+1)%p(2)                     = element%cell%data_pers%QP(indx,3)
       rep%Q(i+1)%b                        = element%cell%data_pers%QP(indx,4)
       
       rep%Q((_SWE_DG_ORDER+1)+i+1)%h      = element%cell%data_pers%FP(1,indx,1)
       rep%Q((_SWE_DG_ORDER+1)+i+1)%p(1)   = element%cell%data_pers%FP(1,indx,2)
       rep%Q((_SWE_DG_ORDER+1)+i+1)%p(2)   = element%cell%data_pers%FP(1,indx,3)
       rep%Q((_SWE_DG_ORDER+1)+i+1)%b      = 0.0_GRID_SR
       
       rep%Q(2*(_SWE_DG_ORDER+1)+i+1)%h    = element%cell%data_pers%FP(2,indx,1)
       rep%Q(2*(_SWE_DG_ORDER+1)+i+1)%p(1) = element%cell%data_pers%FP(2,indx,2)
       rep%Q(2*(_SWE_DG_ORDER+1)+i+1)%p(2) = element%cell%data_pers%FP(2,indx,3)
       rep%Q(2*(_SWE_DG_ORDER+1)+i+1)%b    = 0.0_GRID_SR
    end do
    !-----------------------------------------------------!
    
    !------------------Project FV dofs------------------!
    select case (edge_type)
    case (-3,3) !cells with id i*i+1 (left leg)
       rep%H  = matmul(phi_l,element%cell%data_pers%Q_DG%h) 
       rep%HU = matmul(phi_l,element%cell%data_pers%Q_DG%p(1))
       rep%HV = matmul(phi_l,element%cell%data_pers%Q_DG%p(2)) 
       rep%B  = matmul(phi_l,element%cell%data_pers%Q_DG%b)
    case (-2,2) ! hypotenuse
       rep%H  = matmul(phi_m,element%cell%data_pers%Q_DG%h)
       rep%HU = matmul(phi_m,element%cell%data_pers%Q_DG%p(1))
       rep%HV = matmul(phi_m,element%cell%data_pers%Q_DG%p(2))
       rep%B  = matmul(phi_m,element%cell%data_pers%Q_DG%b)
    case (-1,1) !cells with id i*i (right leg)
       rep%H  = matmul(phi_r,element%cell%data_pers%Q_DG%h)
       rep%HU = matmul(phi_r,element%cell%data_pers%Q_DG%p(1))
       rep%HV = matmul(phi_r,element%cell%data_pers%Q_DG%p(2))
       rep%B  = matmul(phi_r,element%cell%data_pers%Q_DG%b)
    end select
    !---scale by element size---!
    rep%H = (rep%H + rep%B) * _REF_TRIANGLE_SIZE_INV
    rep%HU = rep%HU * _REF_TRIANGLE_SIZE_INV
    rep%HV = rep%HV * _REF_TRIANGLE_SIZE_INV
    rep%B  = rep%B  * _REF_TRIANGLE_SIZE_INV
    !-----------------------------------------------------!
 else
    select case (edge_type)
    case (-3,3) !cells with id i*i+1 (left leg)
       do i=0, _SWE_PATCH_ORDER - 1
          j=_SWE_PATCH_ORDER-1-i
          rep%H(i+1) = element%cell%data_pers%H(j*j + 1)
          rep%HU(i+1)= element%cell%data_pers%HU(j*j + 1)
          rep%HV(i+1)= element%cell%data_pers%HV(j*j + 1)
          rep%B(i+1) = element%cell%data_pers%B(j*j + 1)
       end do
    case (-2,2) ! hypotenuse
       do i=1, _SWE_PATCH_ORDER
          j=_SWE_PATCH_ORDER+1-i
          rep%H(i) = element%cell%data_pers%H ((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
          rep%HU(i)= element%cell%data_pers%HU((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
          rep%HV(i)= element%cell%data_pers%HV((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
          rep%B(i) = element%cell%data_pers%B ((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
       end do
    case (-1,1) !cells with id i*i (right leg)
       do i=1, _SWE_PATCH_ORDER
          !        j=_SWE_PATCH_ORDER+1-i
          rep%H(i)  = element%cell%data_pers%H(i*i)
          rep%HU(i) = element%cell%data_pers%HU(i*i)
          rep%HV(i) = element%cell%data_pers%HV(i*i)
          rep%B(i)  = element%cell%data_pers%B(i*i)
       end do
    end select
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

  normal = edge%transform_data%normal/NORM2(edge%transform_data%normal)
   !        ______
   ! |6 \   \3 2 1|
   ! |4 5 \   \5 4|
   ! |1_2_3\   \ 6|
   
  QL(:,1) = rep1%Q(1:(_SWE_DG_ORDER+1))%h   
  QL(:,2) = rep1%Q(1:(_SWE_DG_ORDER+1))%p(1)
  QL(:,3) = rep1%Q(1:(_SWE_DG_ORDER+1))%p(2)
  QL(:,4) = rep1%Q(1:(_SWE_DG_ORDER+1))%b

  FL(1,:,1) = rep1%Q((_SWE_DG_ORDER+1)+1:(_SWE_DG_ORDER+1)*2)%h   
  FL(1,:,2) = rep1%Q((_SWE_DG_ORDER+1)+1:(_SWE_DG_ORDER+1)*2)%p(1)
  FL(1,:,3) = rep1%Q((_SWE_DG_ORDER+1)+1:(_SWE_DG_ORDER+1)*2)%p(2)

  FL(2,:,1) = rep1%Q(2*(_SWE_DG_ORDER+1)+1:(_SWE_DG_ORDER+1)*3)%h   
  FL(2,:,2) = rep1%Q(2*(_SWE_DG_ORDER+1)+1:(_SWE_DG_ORDER+1)*3)%p(1)
  FL(2,:,3) = rep1%Q(2*(_SWE_DG_ORDER+1)+1:(_SWE_DG_ORDER+1)*3)%p(2)
  
  !---- Dofs for hypothenuse need to be permuted ----!
  if(edge%transform_data%index.eq.2) then
     do j=0,_SWE_DG_ORDER
        QR(j+1,1) = rep2%Q((_SWE_DG_ORDER+1)-j)%h   
        QR(j+1,2) = rep2%Q((_SWE_DG_ORDER+1)-j)%p(1)
        QR(j+1,3) = rep2%Q((_SWE_DG_ORDER+1)-j)%p(2)
        QR(j+1,4) = rep2%Q((_SWE_DG_ORDER+1)-j)%b
        
        FR(1,j+1,1) = rep2%Q((_SWE_DG_ORDER+1)*2-j)%h   
        FR(1,j+1,2) = rep2%Q((_SWE_DG_ORDER+1)*2-j)%p(1)
        FR(1,j+1,3) = rep2%Q((_SWE_DG_ORDER+1)*2-j)%p(2)
        
        FR(2,j+1,1) = rep2%Q((_SWE_DG_ORDER+1)*3-j)%h   
        FR(2,j+1,2) = rep2%Q((_SWE_DG_ORDER+1)*3-j)%p(1)
        FR(2,j+1,3) = rep2%Q((_SWE_DG_ORDER+1)*3-j)%p(2)
     end do
  else
     QR(:,1) = rep2%Q(1:(_SWE_DG_ORDER+1))%h   
     QR(:,2) = rep2%Q(1:(_SWE_DG_ORDER+1))%p(1)
     QR(:,3) = rep2%Q(1:(_SWE_DG_ORDER+1))%p(2)
     QR(:,4) = rep2%Q(1:(_SWE_DG_ORDER+1))%b
     
     FR(1,:,1) = rep2%Q((_SWE_DG_ORDER+1)+1:(_SWE_DG_ORDER+1)*2)%h   
     FR(1,:,2) = rep2%Q((_SWE_DG_ORDER+1)+1:(_SWE_DG_ORDER+1)*2)%p(1)
     FR(1,:,3) = rep2%Q((_SWE_DG_ORDER+1)+1:(_SWE_DG_ORDER+1)*2)%p(2)
     
     FR(2,:,1) = rep2%Q(2*(_SWE_DG_ORDER+1)+1:(_SWE_DG_ORDER+1)*3)%h   
     FR(2,:,2) = rep2%Q(2*(_SWE_DG_ORDER+1)+1:(_SWE_DG_ORDER+1)*3)%p(1)
     FR(2,:,3) = rep2%Q(2*(_SWE_DG_ORDER+1)+1:(_SWE_DG_ORDER+1)*3)%p(2)
  end if

  FLn = FL(1,:,:) * normal(1) + FL(2,:,:) * normal(2)
  FRn = FR(1,:,:) * normal(1) + FR(2,:,:) * normal(2)
  
  call compute_flux_pred(normal,QL,QR,FLn,FRn)

   !---- Result for hypothenuse needs to be permuted back----!
   if(edge%transform_data%index.eq.2) then
         do j=1,_SWE_DG_ORDER+1
            FRn_temp(j,:) = FRn((_SWE_DG_ORDER+1)-j+1,:)
         end do
   else
      FRn_temp = FRn
   end if

   update1%flux%h    = FLn(:,1)
   update1%flux%p(1) = FLn(:,2)
   update1%flux%p(2) = FLn(:,3)
   update2%flux%h    = FRn_temp(:,1)
   update2%flux%p(1) = FRn_temp(:,2)
   update2%flux%p(2) = FRn_temp(:,3)
end subroutine

subroutine skeleton_scalar_op_dg(traversal, grid, edge, rep1, rep2, update1, update2)
type(t_swe_dg_timestep_traversal), intent(in)		:: traversal
type(t_grid_section), intent(in)			:: grid
type(t_edge_data), intent(in)				:: edge
type(num_cell_rep), intent(in)				:: rep1, rep2
type(num_cell_update), intent(out)			:: update1, update2
integer                                                 :: i,j, edge_type
real(kind = GRID_SR)	          :: normal(2)


update1%Q_DG          =rep2%Q_DG
update2%Q_DG          =rep1%Q_DG

update1%troubled      =rep2%troubled
update2%troubled      =rep1%troubled

if(rep1%troubled.le.0 .and. rep2%troubled.le.0) then
   call general_dg_riemannsolver(edge,rep1,rep2,update1,update2)
end if

update2%H  = rep1%H
update2%HU = rep1%HU
update2%HV = rep1%HV
update2%B  = rep1%B

update1%H  = rep2%H
update1%HU = rep2%HU
update1%HV = rep2%HV
update1%B  = rep2%B
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

if(rep%troubled.le.0) then
   update_bnd = update
   
   rep_bnd%Q(:)%h = rep%Q(:)%h
   rep_bnd%Q(:)%b = rep%Q(:)%b

   do i=1,(_SWE_DG_ORDER+1)
      length_flux = dot_product(rep%Q(i)%p, normal)
      rep_bnd%Q(i)%p(1) = rep%Q(i)%p(1)-2.0_GRID_SR*length_flux*normal(1)
      rep_bnd%Q(i)%p(2) = rep%Q(i)%p(2)-2.0_GRID_SR*length_flux*normal(2)
   end do
   normal = abs(normal)
   do i=1+(_SWE_DG_ORDER+1)  ,2*(_SWE_DG_ORDER+1)
      rep_bnd%Q(i)%h    = rep%Q(i)%h    - 2.0_GRID_SR * rep%Q(i-(_SWE_DG_ORDER+1))%p(1) * normal(1)
      rep_bnd%Q(i)%p(1) = rep%Q(i)%p(1) - 2.0_GRID_SR * rep%Q(i)%p(1) * normal(2)
      rep_bnd%Q(i)%p(2) = rep%Q(i)%p(2) - 2.0_GRID_SR * rep%Q(i)%p(2) * normal(1)
   end do

   do i=1+(_SWE_DG_ORDER+1)*2,3*(_SWE_DG_ORDER+1)
      rep_bnd%Q(i)%h    = rep%Q(i)%h    - 2.0_GRID_SR * rep%Q(i-(_SWE_DG_ORDER+1)*2)%p(2) * normal(2)
      rep_bnd%Q(i)%p(1) = rep%Q(i)%p(1) - 2.0_GRID_SR * rep%Q(i)%p(1) * normal(2)
      rep_bnd%Q(i)%p(2) = rep%Q(i)%p(2) - 2.0_GRID_SR * rep%Q(i)%p(2) * normal(1)
   end do
   
   call general_dg_riemannsolver(edge,rep,rep_bnd,update,update_bnd)   
end if
   
update%H=rep%H
update%HU=rep%HU
update%HV=rep%HV
update%B=rep%B

update%Q_DG(:)%h = rep%Q_DG(:)%h
update%Q_DG(:)%b = rep%Q_DG(:)%b

normal=(edge%transform_data%normal)/NORM2(edge%transform_data%normal)
!---DG velocities---!
do i=1,(_SWE_DG_DOFS)
   length_flux = dot_product(rep%Q_DG(i)%p, normal)
   ! reflecting
   update%Q_DG(i)%p(1) = rep%Q_DG(i)%p(1)-2.0_GRID_SR*length_flux*normal(1)
   update%Q_DG(i)%p(2) = rep%Q_DG(i)%p(2)-2.0_GRID_SR*length_flux*normal(2)
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
real (kind=GRID_SR) :: max_wave_speed, wave_speed=0,dQ_norm,b_min,b_max
real (kind=GRID_SR) :: refinement_threshold = 0.25_GRID_SR/_SWE_DG_ORDER
real (kind = GRID_SR), dimension (_SWE_DG_DOFS) :: H_old, HU_old, HV_old, B_old
real (kind = GRID_SR), dimension (_SWE_PATCH_ORDER_SQUARE) :: H, HU, HV
logical :: drying,troubled,neighbours_troubled,refine,coarsen
integer :: i_depth,bathy_depth
real(kind=GRID_SR),Dimension(3) :: edge_sizes     

if (element%cell%geometry%i_plotter_type > 0) then ! if orientation = forward, reverse updates
   tmp=update1
   update1=update3
   update3=tmp
end if

neighbours_troubled=(update1%troubled.ge.1).or.(update2%troubled.ge.1).or.(update3%troubled.ge.1)

associate(data => element%cell%data_pers)
  
if(neighbours_troubled) then
   data%troubled=merge(data%troubled,2,data%troubled.ge.1)
end if

if(data%troubled.eq.2) then

   call apply_phi(data%Q_DG(:)%h+data%Q_DG(:)%b,data%h)
   call apply_phi(data%Q_DG(:)%p(1),data%hu)
   call apply_phi(data%Q_DG(:)%p(2),data%hv)

   data%b=get_bathymetry_at_patch(section, element, section%r_time)

end if

if(data%troubled.le.0) then
   data%troubled=0

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

   max_neighbour(1)=max(maxval(update1%Q_DG(:)%H)   ,maxval(update2%Q_DG(:)%H)   ,maxval(update3%Q_DG(:)%H)   ,data_max_val(1))
   max_neighbour(2)=max(maxval(update1%Q_DG(:)%p(1)),maxval(update2%Q_DG(:)%p(1)),maxval(update3%Q_DG(:)%p(1)),data_max_val(1))
   max_neighbour(3)=max(maxval(update1%Q_DG(:)%p(2)),maxval(update2%Q_DG(:)%p(2)),maxval(update3%Q_DG(:)%p(2)),data_max_val(1))

   min_neighbour(1)=min(minval(update1%Q_DG(:)%H)   ,minval(update2%Q_DG(:)%H)   ,minval(update3%Q_DG(:)%H)   ,data_min_val(1))
   min_neighbour(2)=min(minval(update1%Q_DG(:)%p(1)),minval(update2%Q_DG(:)%p(1)),minval(update3%Q_DG(:)%p(1)),data_min_val(1))
   min_neighbour(3)=min(minval(update1%Q_DG(:)%p(2)),minval(update2%Q_DG(:)%p(2)),minval(update3%Q_DG(:)%p(2)),data_min_val(1))

   delta(1) = max(0.1_GRID_SR,max_neighbour(1)-min_neighbour(1))
   delta(2) = max(0.1_GRID_SR,max_neighbour(2)-min_neighbour(2))
   delta(3) = max(0.1_GRID_SR,max_neighbour(3)-min_neighbour(3))

   delta=delta*1.0e-3_GRID_SR
   

   call dg_solver(element,update1%flux,update2%flux,update3%flux,section%r_dt)


#if defined(_SWE_DG_LIMITER_UNLIMITED)
   troubled=.false.
#elif defined(_SWE_DG_LIMITER_HEIGHT)
   troubled=.not.all(data%Q_DG%H < max_neighbour(1)+delta(1).and.data%Q_DG%H > min_neighbour(1)-delta(1))
#elif defined(_SWE_DG_LIMITER_ALL)
   troubled=&
        .not.all(data%Q_DG%H   < max_neighbour(1)+delta(1).and.data%Q_DG%H >min_neighbour(1)-delta(1)) .or.&
        .not.all(data%Q_DG%p(1)< max_neighbour(2)+delta(2).and.data%Q_DG%p(1)>min_neighbour(2)-delta(2)) .or.&
        .not.all(data%Q_DG%p(2)< max_neighbour(3)+delta(3).and.data%Q_DG%p(2)>min_neighbour(3)-delta(3))
#endif
   
   !   drying=.not.all(data%H - data%B > cfg%dry_tolerance*50.0).or..not.all(data%Q_DG%H > cfg%dry_tolerance*50.0)
   drying=.not.all(data%Q_DG%H > cfg%coast_height)

   if(troubled.or.drying) then
      !--if troubled or drying perform rollback--!
      data%Q_DG%H   = H_old
      data%Q_DG%p(1)= HU_old
      data%Q_DG%p(2)= HV_old
      data%Q_DG%B   = B_old
      data%troubled=merge(1,3,drying)

      call apply_phi(data%Q_DG(:)%h+data%Q_DG(:)%b,data%h)
      call apply_phi(data%Q_DG(:)%p(1),data%hu)
      call apply_phi(data%Q_DG(:)%p(2),data%hv)
   !   call apply_phi(data%Q_DG(:)%b,data%b)
      data%b=get_bathymetry_at_patch(section, element, section%r_time)

   else
      do i=1,_SWE_DG_DOFS
         wave_speed =  sqrt(g * (data%Q_DG(i)%h)) + maxval(abs(data%Q_DG(i)%p/data%Q_DG(i)%h))

         section%r_dt_new = min(section%r_dt_new,cfg%scaling*  element%transform_data%custom_data%scaling  / (wave_speed* (_SWE_DG_ORDER*4.0_GRID_SR +2.0_GRID_SR)))
         max_wave_speed = max(max_wave_speed,wave_speed)
      end do
      i_depth = element%cell%geometry%i_depth
      
      dQ_norm = maxval(abs(data%Q_DG%H-H_old))
      
      edge_sizes=element%cell%geometry%get_edge_sizes()

      !refine      
      !consider wave
      refine = (dQ_norm/section%r_dt  * edge_sizes(1) * edge_sizes(1)) > refinement_threshold *get_edge_size(cfg%i_max_depth)**2
      
      !coarsen
      !consider wave
      
#if defined (_SWE_SCENARIO_OSCILLATING_LAKE)
      coarsen=.false.
#else      
      coarsen=(( dQ_norm/section%r_dt  * edge_sizes(1) * edge_sizes(1)) < refinement_threshold/8_GRID_SR *get_edge_size(cfg%i_max_depth)**2)
#endif
      

#if defined(_ASAGI)
      ! consider bathymetrie
      if(coarsen) then
         section%b_min=-0.1_GRID_SR
         section%b_max=0.1_GRID_SR

         if((minval(element%cell%data_pers%Q_DG%B) * section%b_min) > 0) then
            
            bathy_depth= cfg%i_max_depth &
            -floor(abs(minval(element%cell%data_pers%Q_DG%B)/abs(section%b_min)) &
            *(cfg%i_max_depth-cfg%i_min_depth))
            
         else if((maxval(element%cell%data_pers%Q_DG%B) * section%b_max) > 0) then

            bathy_depth= cfg%i_max_depth &
            -floor(abs(maxval(element%cell%data_pers%Q_DG%B)/abs(section%b_max)) &
            *(cfg%i_max_depth-cfg%i_min_depth))
            
         end if
         
         bathy_depth = max(min(bathy_depth,cfg%i_max_depth),cfg%i_min_depth)
      end if

      !!      coarsen=coarsen.and. (cfg%i_max_depth-floor(min(element%cell%data_pers%Q_DG%b)/(section%section%b_min/(cfg%i_max_depth-cfg%i_min_depth))) < i_depth)
      coarsen=coarsen.and. (i_depth > bathy_depth)
      
!      refine=.false.
!      coarsen=.false.      
#endif      
      
      if (i_depth < cfg%i_max_depth .and. refine) then
         element%cell%geometry%refinement = 1
         traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
      else if (i_depth > cfg%i_min_depth .and. coarsen) then
         element%cell%geometry%refinement = -1         
      endif
      
   end if

end if   

!if cell is troubled, compute fv solution and mark all edges as troubled

if(data%troubled.ge.1) then
   !-----Call FV patch solver----!
   call fv_patch_solver(traversal, section, element, update1, update2, update3)
end if
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

VelL(:) = QL(:,2)/QL(:,1) * normal(1) +  QL(:,3)/QL(:,1) * normal(2)
VelR(:) = QR(:,2)/QR(:,1) * normal(1) +  QR(:,3)/QR(:,1) * normal(2)

alpha=max(maxval(abs(VelL) + sqrt(g * QL(:,1))),&
     maxval(abs(VelR) + sqrt(g * QR(:,1))))

hRoe  = 0.5_GRID_SR * (QR(:,1) + QL(:,1))
bm    = max(QR(:,4),QL(:,4))
Deta  = max(QR(:,1) + QR(:,4) - bm, 0.0) - max(QL(:,1) + QL(:,4) - bm, 0.0)
Djump = 0.5_GRID_SR * g * hRoe * Deta

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
            Q_DG        => element%cell%data_pers%Q_DG, &
            Q_DG_UPDATE => element%cell%data_pers%Q_DG_UPDATE)

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
    !!----update dofs----!!
    q=q-((bnd_flux_l + bnd_flux_m + bnd_flux_r) - Q_DG_UPDATE) * dt/dx
    !!-------------------!!

    call data%set_dofs_dg(q)
    
  end associate
  
end subroutine dg_solver

subroutine get_fv_update(update,leg)
  type(num_cell_update), intent(inout) :: update
  integer,intent(in) :: leg
  real(kind=GRID_SR),Dimension(_SWE_PATCH_ORDER_SQUARE) :: H, HU, HV ,B
  integer :: i,j
  
  call apply_phi(update%Q_DG(:)%h,H)
  call apply_phi(update%Q_DG(:)%p(1),HU)
  call apply_phi(update%Q_DG(:)%p(2),HV)
  call apply_phi(update%Q_DG(:)%b,B)

  H=H+B

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
  case (2) ! hypothenuse
     do i=1, _SWE_PATCH_ORDER
        j=_SWE_PATCH_ORDER+1-i
        update%H(i) = H ((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
        update%HU(i)= HU((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
        update%HV(i)= HV((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
        update%B(i) = B ((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
     end do
  case (1) !cells with id i*i (right leg)
     do i=1, _SWE_PATCH_ORDER
        update%H(i)  = H(i*i)
        update%HU(i) = HU(i*i)
        update%HV(i) = HV(i*i)
        update%B(i)  = B(i*i)
     end do
  end select

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
            logical :: drying,troubled,coarsen,refine
            real (kind=GRID_SR) :: refinement_threshold = 0.50_GRID_SR

            !print*,"updates"
            !print*,element%cell%geometry%i_plotter_type
            !print*,update1%H
            !print*,update1%HU
            !print*,update1%HV
            !print*,update1%B
            !print*
            !print*,update2%H
            !print*,update2%HU
            !print*,update2%HV
            !print*,update2%B
            !print*                 
            !print*,update3%H
            !print*,update3%HU
            !print*,update3%HV
            !print*,update3%B
            !print*                 

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
#               endif


            ! if (element%cell%geometry%i_plotter_type > 0) then ! if orientation = forward, reverse updates
            !    tmp=update1
            !    update1=update3
            !    update3=tmp
            ! end if

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
#                           if defined(_FWAVE_FLUX) || defined(_AUG_RIEMANN_FLUX)
                 call compute_updates_simd(transf, hL, huL, hvL, bL, hR, huR, hvR, bR, upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR, maxWaveSpeed)
#                           elif defined(_HLLE_FLUX)
                 call compute_updates_hlle_simd(transf, hL, huL, hvL, bL, hR, huR, hvR, bR, upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR, maxWaveSpeed)
#                           else
                 ! this should never happen -> SCons rules should avoid this before compiling
#                               error "No valid SWE solver for patches/simd implementation has been defined!"
                 print *, 
#                           endif
#                       else                    
                 ! using original geoclaw solver
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
!                    maxWaveSpeed = max(maxWaveSpeed, update_a%max_wave_speed)
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

               ! ! if the water level falls below the dry tolerance, set water surface to 0 and velocity to 0
               dQ_H  = dQ_H * (-dt_div_volume)
               dQ_HU = dQ_HU * (-dt_div_volume)
               dQ_HV = dQ_HV * (-dt_div_volume)

               ! if land is flooded, init water height to dry tolerance and
               ! velocity to zero
               where (data%H < data%B + cfg%dry_tolerance .and. dQ_H > 0.0_GRID_SR)
#if defined (_ASAGI)                  
                  data%H = min(0.0_GRID_SR,data%B)
#else
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
                  data%H = data%B + cfg%dry_tolerance
#endif                  
                  data%HU = 0.0_GRID_SR
                  data%HV = 0.0_GRID_SR
               end where

               ! maxWaveSpeed=maxval(maxWaveSpeeds)
               if(maxWaveSpeed > 0.0_GRID_SR) then
                  section%r_dt_new = min(section%r_dt_new, volume / (edge_lengths(2) * maxWaveSpeed) )
                  maxWaveSpeed=0
               end if


               call apply_mue(data%h          ,data%Q_DG%h)
               call apply_mue(data%hu         ,data%Q_DG%p(1))
               call apply_mue(data%hv         ,data%Q_DG%p(2))

               data%Q_DG%b=get_bathymetry_at_dg_patch(section, element, section%r_time)
               call bathymetry_derivatives(element%cell%data_pers,ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%jacobian)
               
               data%Q_DG%h=data%Q_DG%h-data%Q_DG%b
               if(all(data%H - data%B > cfg%coast_height)) then
                  if(all(data%Q_DG%H > cfg%coast_height)) then
                     data%troubled = 5+data%troubled
                     
                  else
                     data%troubled = 1
                  end if
               else
                  data%troubled = 1
               end if

               coarsen = all(data%H - data%B < cfg%dry_tolerance)
               refine  = .not.all(data%H - data%B < cfg%dry_tolerance)


#if defined(_ASAGI)
               coarsen=.false.
               refine =.false.
#endif               
               
               if(data%troubled.eq.1) then
                  if ( element%cell%geometry%i_depth < cfg%i_max_depth .and. refine) then
                     element%cell%geometry%refinement = cfg%i_max_depth - element%cell%geometry%i_depth
                     traversal%i_refinements_issued = traversal%i_refinements_issued+ cfg%i_max_depth - element%cell%geometry%i_depth
                  else if ( element%cell%geometry%i_depth > cfg%i_min_depth .and. coarsen) then
                     element%cell%geometry%refinement = -1
                  endif
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
