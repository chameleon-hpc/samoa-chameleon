! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2017 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE
! Author: Leonhard Rannabauer rannabau (at) in.tum.de


#include "Compilation_control.f90"

#if defined(_SWE_DG)

MODULE SWE_dg_predictor
  use SFC_edge_traversal
  use Samoa_swe
  use c_bind_riemannsolvers

#		if defined(_SWE_PATCH)
!  use SWE_PATCH
!  use SWE_PATCH_Solvers
#               if defined(_SWE_DG)
  use SWE_DG_Matrices
#		endif
#		endif

#               if defined(_HLLE_FLUX)
  use SWE_HLLE
#               endif

  implicit none
  type num_traversal_data
  end type num_traversal_data

  interface skeleton_op
     module procedure skeleton_array_op
     module procedure skeleton_scalar_op
  end interface skeleton_op
  
  interface bnd_skeleton_op
     module procedure bnd_skeleton_array_op
     module procedure bnd_skeleton_scalar_op
  end interface bnd_skeleton_op
  
  
  public dg_predictor,flux_1,flux_2
  
#		define _GT_NAME				t_swe_dg_predictor_traversal
#		define _GT_EDGES
  
#		define _GT_PRE_TRAVERSAL_OP		pre_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP	pre_traversal_grid_op

#		define _GT_CELL_TO_EDGE_OP	        cell_to_edge_op
  
#		define _GT_SKELETON_OP			skeleton_op
#		define _GT_BND_SKELETON_OP		bnd_skeleton_op
#		define _GT_ELEMENT_OP                   element_op
#		define _GT_CELL_UPDATE_OP		cell_update_op

#		include "SFC_generic_traversal_ringbuffer.f90"


  subroutine pre_traversal_grid_op(traversal, grid)
    type(t_swe_dg_predictor_traversal), intent(inout)		:: traversal
    type(t_grid), intent(inout)							    :: grid

    if (cfg%r_max_time > 0.0_SR) then
       grid%r_dt = min(cfg%r_max_time-grid%r_time, grid%r_dt) 
       if(mod(grid%r_time,cfg%r_output_time_step) > 0.0_GRID_SR) then
          grid%r_dt = min(cfg%r_output_time_step-mod(grid%r_time,cfg%r_output_time_step), grid%r_dt)
       end if
    end if


    if (cfg%r_output_time_step > 0.0_SR) then
       grid%r_dt = min(cfg%r_output_time_step, grid%r_dt)
    end if

#           if defined(_ASAGI)
    !if we are in the earthquake phase, limit the simulation time step by the earthquake time step
    if (grid%r_time < cfg%t_max_eq) then
       grid%r_dt = min(grid%r_dt, cfg%dt_eq)
    end if
#           endif
    call scatter(grid%r_dt, grid%sections%elements_alloc%r_dt)

  end subroutine pre_traversal_grid_op


  subroutine pre_traversal_op(traversal, section)
    type(t_swe_dg_predictor_traversal), intent(inout)	:: traversal
    type(t_grid_section), intent(inout)		:: section
  end subroutine pre_traversal_op
  
  subroutine dg_predictor(cell,dt)
    class(t_cell_data_ptr), intent(in)		:: cell
    real(kind=GRID_SR) ,intent(in)             :: dt
    real(kind=GRID_SR)                         :: dx
    
    real(kind=GRID_SR)                         :: q_0(_SWE_DG_DOFS,3)
    real(kind=GRID_SR)                         :: q_i(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3)
    real(kind=GRID_SR)                         :: q_temp(_SWE_DG_DOFS*_SWE_DG_ORDER,3)
    
    !--local variables--!
    integer                                    :: i,j,offset
    real(kind=GRID_SR),Dimension(2,2)          :: jacobian,jacobian_inv
    real(kind=GRID_SR),Dimension(3)            :: edge_sizes
    real(kind=GRID_SR)                         :: epsilon=1.0_GRID_SR
    
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3)  :: f1,f1_ref,f2,f2_ref
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*_SWE_DG_ORDER,3)      :: source,source_ref
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*_SWE_DG_ORDER,3)      :: volume_flux1,volume_flux2
    
    
#    if defined(_SWE_DG_NODAL)
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1))    :: H_x,H_y
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS)    :: b,H
#    else
    !     real(kind=GRID_SR),Dimension(size(st_gl_node_vals,1))   :: b_x,b_y
    !     real(kind=GRID_SR),Dimension(size(st_gl_node_vals,1),3) :: q_v
    !     integer :: st_nodes = size(st_gl_node_vals,1)
    !     integer :: s_nodes = size(s_der_x_gl_node_vals,1)
#    endif
    
    associate(Q_DG =>cell%data_pers%Q_DG ,Q_DG_P =>cell%data_pers%Q_DG_P)
      
      !--Normalize jacobian--!
      !TODO:do this only one time
      jacobian=ref_plotter_data(abs(cell%geometry%i_plotter_type))%jacobian_normalized
      
      jacobian_inv=ref_plotter_data(abs(cell%geometry%i_plotter_type))%jacobian_inv_normalized
      
      edge_sizes=cell%geometry%get_edge_sizes()
      dx=edge_sizes(1)*cfg%scaling
      
      call cell%data_pers%get_dofs_dg(q_0)
      
      !--Set initial conditions for discrete picard iteration--!
      q_i(1:_SWE_DG_DOFS,:) = q_0
      b=Q_DG(:)%B
      
      do i=1,_SWE_DG_ORDER
         q_i(1+i*_SWE_DG_DOFS:(i+1)*_SWE_DG_DOFS,:) = q_0
      end do
      
      i=0
      do while(epsilon > 1.0e-14_GRID_SR)
         i=i+1
#if defined(_SWE_DG_NODAL)

#if defined(_DEBUG)  
         if (any(q_i(:,1).le.0)) then
            write(*), "ERROR: Waterheigth less than zero in predictor"
            exit
         end if
#endif

         !--- Flux on ref
         f1_ref = flux_1(q_i,(_SWE_DG_ORDER+1)*_SWE_DG_DOFS)
         f2_ref = flux_2(q_i,(_SWE_DG_ORDER+1)*_SWE_DG_DOFS)
         
         !--space derivations needed for total waterheigth needed for well-balanced scheme--!
         do j=0,(_SWE_DG_ORDER)
            offset=_SWE_DG_DOFS * j
            H=q_i(offset+1:offset+_SWE_DG_DOFS,1)+b
            H_x(offset+1:offset+_SWE_DG_DOFS)= matmul(basis_der_x,H)
            H_y(offset+1:offset+_SWE_DG_DOFS)= matmul(basis_der_y,H)
         end do
         
         
         !---Source on reference element---!
         source_ref=0
         source_ref(:,2)=-matmul(st_m,(g*q_i(:,1)*H_x))
         source_ref(:,3)=-matmul(st_m,(g*q_i(:,1)*H_y))
         source_ref(:,2)=source_ref(:,2)+matmul(st_k_x,(0.5_GRID_SR*g*q_i(:,1)**2))
         source_ref(:,3)=source_ref(:,3)+matmul(st_k_y,(0.5_GRID_SR*g*q_i(:,1)**2))
         
         !---Transform back to element ---!
         source=0
         source(:,2) = jacobian(1,1) * source_ref(:,2) + jacobian(1,2) * source_ref(:,3)
         source(:,3) = jacobian(2,1) * source_ref(:,2) + jacobian(2,2) * source_ref(:,3)
         
#else
         ! q_v=matmul(st_gl_node_vals,q_i)
         ! f1_hat=matmul(st_gl_weights,f1(q_v,st_nodes))
         ! f2_hat=matmul(st_gl_weights,f2(q_v,st_nodes))
         
         ! b_x=cell%data_pers%b_x
         ! b_y=cell%data_pers%b_y
         
         ! S_s=matmul(st_gl_weights,S(q_v(:,1),b_x,b_y,st_nodes))
         
         ! call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f1_hat(:,1))
         ! call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f1_hat(:,2))
         ! call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f1_hat(:,3))
         ! call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f2_hat(:,1))
         ! call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f2_hat(:,2))
         ! call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f2_hat(:,3))
         
         ! source = S_s(_SWE_DG_DOFS+1:,:)
#endif
         
         !---Transform flux back to element---!
         f1=(jacobian_inv(1,1)*f1_ref+jacobian_inv(1,2)*f2_ref)
         f2=(jacobian_inv(2,1)*f1_ref+jacobian_inv(2,2)*f2_ref)
         
         source  = source*dt/dx
         volume_flux1 =matmul(st_k_x,f1)*dt/dx
         volume_flux2 =matmul(st_k_y,f2)*dt/dx
         
         !TODO only calculate new dofs
         q_temp= source &
              - matmul(st_w_k_t_1_0,q_i(1:_SWE_DG_DOFS,:))&
              - volume_flux1 &
              - volume_flux2
         
         call lusolve(st_w_k_t_1_1_lu,_SWE_DG_ORDER*_SWE_DG_DOFS,st_w_k_t_1_1_lu_pivot ,q_temp(:,1))
         call lusolve(st_w_k_t_1_1_lu,_SWE_DG_ORDER*_SWE_DG_DOFS,st_w_k_t_1_1_lu_pivot ,q_temp(:,2))
         call lusolve(st_w_k_t_1_1_lu,_SWE_DG_ORDER*_SWE_DG_DOFS,st_w_k_t_1_1_lu_pivot ,q_temp(:,3))
         
         epsilon=0.0_GRID_SR                          
         do j=1,(_SWE_DG_ORDER)*_SWE_DG_DOFS
            if(.not.all(q_i(_SWE_DG_DOFS+j,:).eq.0)) then
               epsilon=max(epsilon,NORM2(q_temp(j,:)-q_i(_SWE_DG_DOFS+j,:))/NORM2(q_i(_SWE_DG_DOFS+j,:)))
            else if(.not.all(q_temp(j,:).eq. 0)) then
               epsilon=max(epsilon,1.0_GRID_SR)
            end if
         end do
         
         if(i > 200) then                           
            print*,"predictor not converging"
            exit
         end if
         
         q_i(1+_SWE_DG_DOFS:(_SWE_DG_ORDER+1)*_SWE_DG_DOFS,:) = q_temp
         
      end do
      
      call cell%data_pers%set_dofs_pred(q_i)
      
    end  associate
  end subroutine dg_predictor
  

  subroutine element_op(traversal, section, element) 
  type(t_swe_dg_predictor_traversal), intent(inout)				:: traversal
  type(t_grid_section), intent(inout)						:: section
  type(t_element_base), intent(inout)						:: element
  integer :: i


  if(element%cell%data_pers%troubled.le.0) then
!!!print*,"predict"
     call dg_predictor(element%cell,section%r_dt)
!!!print*,"predict done"
     ! else
     !    element%cell%data_pers%Q_DG_P(:)%H=0
     !    element%cell%data_pers%Q_DG_P(:)%p(1)=0
     !    element%cell%data_pers%Q_DG_P(:)%p(2)=0
     !    element%cell%data_pers%Q_DG_P(1:_SWE_DG_DOFS) = element%cell%data_pers%Q_DG
  end if
end subroutine element_op

function flux_1(q,N)
  real(kind=GRID_SR)             ::flux_1(N,3)
  real(kind=GRID_SR) ,intent(in) ::q(N,3)
  integer :: N
!  integer                        ::i
  flux_1(:,1) = q(:,2)
  flux_1(:,2) = q(:,2)**2/q(:,1) + 0.5_GRID_SR * g * q(:,1)**2
  flux_1(:,3) = q(:,2)*q(:,3)/q(:,1)
  
  ! do i=1,N
  !    if(q(i,1).eq.0) then
  !       f1(i,1) = 0
  !       f1(i,2) = 0
  !       f1(i,3) = 0
  !    else
  !    end if
  ! end do

end function flux_1

function flux_2(q,N)
  real(kind=GRID_SR)             ::flux_2(N,3)
  real(kind=GRID_SR) ,intent(in) ::q(N,3)
  integer :: N
  !  integer                        ::i
  flux_2(:,1) = q(:,3)
  flux_2(:,2) = q(:,2)*q(:,3)/q(:,1)
  flux_2(:,3) = q(:,3)**2/q(:,1) + 0.5_GRID_SR * g * q(:,1)**2
  ! do i=1,N
  !    if(q(i,1).eq.0) then
  !       f2(i,1) = 0
  !       f2(i,2) = 0
  !       f2(i,3) = 0
  !    else
  !    end if
  ! end do
end function flux_2

! !source
! function S(q,b_x,b_y,N)
!   integer ::N,i
!   real(kind=GRID_SR) ,intent(in) ::q(N,3), b_x(N), b_y(N)
!   real(kind=GRID_SR)             ::S(N,3)

!   S(:,1) = 0
!   do i=1,N
!      S(i,2)=-q(i,1)*(b_x(i))*g
!      S(i,3)=-q(i,1)*(b_y(i))*g 

!   end do

! end function S

function cell_to_edge_op(element, edge) result(rep)
  type(t_element_base), intent(in)						:: element
  type(t_edge_data), intent(in)						    :: edge
  type(num_cell_rep)										:: rep

end function cell_to_edge_op

subroutine skeleton_array_op(traversal, grid, edges, rep1, rep2, update1, update2)
  type(t_swe_dg_predictor_traversal), intent(in)				:: traversal
  type(t_grid_section), intent(in)							    :: grid
  type(t_edge_data), intent(in)								    :: edges(:)
  type(num_cell_rep), intent(in)									:: rep1(:), rep2(:)
  type(num_cell_update), intent(out)								:: update1(:), update2(:)
  integer (kind = GRID_SI)                                        :: i

  do i = 1, size(edges)
     call skeleton_scalar_op(traversal, grid, edges(i), rep1(i), rep2(i), update1(i), update2(i))
  end do

end subroutine skeleton_array_op

subroutine skeleton_scalar_op(traversal, grid, edge, rep1, rep2, update1, update2)
  type(t_swe_dg_predictor_traversal), intent(in)				:: traversal
  type(t_grid_section), intent(in)							    :: grid
  type(t_edge_data), intent(in)								    :: edge
  type(num_cell_rep), intent(in)									:: rep1, rep2
  type(num_cell_update), intent(out)								:: update1, update2
!!!!!!!!!!!print*,rep1%H

end subroutine skeleton_scalar_op

subroutine bnd_skeleton_array_op(traversal, grid, edges, rep, update)
  type(t_swe_dg_predictor_traversal), intent(in)				:: traversal
  type(t_grid_section), intent(in)							    :: grid
  type(t_edge_data), intent(in)								    :: edges(:)
  type(num_cell_rep), intent(in)									:: rep(:)
  type(num_cell_update), intent(out)								:: update(:)
  integer :: i

  do i = 1, size(edges)
     call bnd_skeleton_scalar_op(traversal, grid, edges(i), rep(i), update(i))
  end do


end subroutine bnd_skeleton_array_op

subroutine bnd_skeleton_scalar_op(traversal, grid, edge, rep, update)
  type(t_swe_dg_predictor_traversal), intent(in)				:: traversal
  type(t_grid_section), intent(in)							    :: grid
  type(t_edge_data), intent(in)								    :: edge
  type(num_cell_rep), intent(in)									:: rep
  type(num_cell_update), intent(out)								:: update

end subroutine bnd_skeleton_scalar_op


subroutine cell_update_op(traversal, section, element, update1, update2, update3)
  type(t_swe_dg_predictor_traversal), intent(inout)				:: traversal
  type(t_grid_section), intent(inout)							:: section
  type(t_element_base), intent(inout)						:: element
  type(num_cell_update), intent(inout)						:: update1, update2, update3
  integer :: i

  element%cell%data_pers%troubled_old = element%cell%data_pers%troubled
end subroutine cell_update_op

END MODULE SWE_DG_predictor

#endif
