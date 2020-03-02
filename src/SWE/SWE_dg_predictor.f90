! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2017 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE
! Author: Leonhard Rannabauer rannabau (at) in.tum.de

#include "Compilation_control.f90"

#if defined(_SWE_DG)

MODULE SWE_dg_predictor
  use SFC_edge_traversal
  use Samoa_swe
  use SWE_DG_Limiter
  implicit none
  
  public dg_predictor,flux_1,flux_2
  
contains

    subroutine dg_predictor(element,dt)
    class(t_element_base), intent(inout)       :: element
    real(kind=GRID_SR) ,intent(in)             :: dt
    real(kind=GRID_SR)                         :: dx
    real(kind=GRID_SR)                         :: q_0(_SWE_DG_DOFS,3)
    real(kind=GRID_SR)                         :: q_i(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3)
    
    real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,_SWE_DG_DOFS,3)  :: q_i_st
    real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,_SWE_DG_DOFS,3)  :: source_st,source_ref_st
    real(kind=GRID_SR),Dimension(_SWE_DG_ORDER,_SWE_DG_DOFS,3)    :: source_st_red
    real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,_SWE_DG_DOFS,3)  :: volume_flux
    real(kind=GRID_SR),Dimension(_SWE_DG_ORDER,_SWE_DG_DOFS,3)    :: volume_flux_red
    
    real(kind=GRID_SR),Dimension(2,_SWE_DG_ORDER+1,_SWE_DG_DOFS,3) :: f,f_ref

    real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1, _SWE_DG_DOFS, 3, 2) :: s_ref
!    real(kind=GRID_SR),Dimension(2, _SWE_DG_ORDER+1, _SWE_DG_DOFS,3) :: s
    
    real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,_SWE_DG_DOFS) :: H_x_st,H_y_st
    real(kind=GRID_SR),Dimension(_SWE_DG_ORDER,_SWE_DG_DOFS,3) :: q_temp_st
    
    !--local variables--!
    integer                                    :: iteration,i,j,offset,edge_type,indx
    real(kind=GRID_SR),Dimension(2,2)          :: jacobian,jacobian_inv
    real(kind=GRID_SR),Dimension(3)            :: edge_sizes
    real(kind=GRID_SR)                         :: epsilon=1.0_GRID_SR

    real(kind=GRID_SR), DIMENSION(_SWE_DG_DOFS,3)    :: Q_DG_UPDATE
    real(kind=GRID_SR), DIMENSION(2,_SWE_DG_DOFS,3)  :: FP
    real(kind=GRID_SR), DIMENSION(  _SWE_DG_DOFS,4)  :: QP


    associate(Q_DG        => element%cell%data_pers%Q,&
              cell  => element%cell,&
              edges => element%edges)
!              FP => cell%data_pers%FP,&
!              QP => element%cell%data_pers%QP,&      

      ! TODO make this a precompiled matrix

      
      !--Normalize jacobian--!
      jacobian=ref_plotter_data(abs(cell%geometry%i_plotter_type))%jacobian_normalized
      jacobian_inv=ref_plotter_data(abs(cell%geometry%i_plotter_type))%jacobian_inv_normalized
      
      edge_sizes=cell%geometry%get_edge_sizes()
      dx=edge_sizes(1)*cfg%scaling
      
      call cell%data_pers%get_dofs_dg(q_0)
      
      !!--Set initial conditions for discrete picard iteration--!!
      !--span dofs at t=0 over time basis--!
      do i=1,_SWE_DG_ORDER+1
         q_i_st(i,:,:) = q_0
      end do
      
      iteration=0
      epsilon=1.0_GRID_SR
      !!---------------------------!!
      do while(epsilon > 1.0e-14_GRID_SR)
         iteration=iteration+1
         
         if (any(q_i_st(:,:,1).le.0)) then
! #if defined(_DEBUG)  
!             print*,"PRED"
!             print*,q_i_st
!             print*,"FV"
!             print*,cell%data_pers%H
!             write(*,*), "ERROR: Waterheight less than zero in predictor"
! #endif
            exit
         end if

         !!--- Compute F S1 and S2---!!
         s_ref = 0
         do i=1,_SWE_DG_ORDER+1
            s_ref(i,:,2,1) = ( g * q_i_st(i,:,1)    * matmul(basis_der_x,q_i_st(i,:,1) + Q_DG(:)%B) )
            s_ref(i,:,3,1) = ( g * q_i_st(i,:,1)    * matmul(basis_der_y,q_i_st(i,:,1) + Q_DG(:)%B) )
            s_ref(i,:,2,2) = ( g * q_i_st(i,:,1)**2 * 0.5_GRID_SR )
            s_ref(i,:,3,2) = ( g * q_i_st(i,:,1)**2 * 0.5_GRID_SR )
         end do
         
         f_ref = 0
         do i=1,_SWE_DG_ORDER+1
            f_ref(:,i,:,:) = flux(q_i_st(i,:,:),_SWE_DG_DOFS)
         end do
        
         !!--------Run kernels-------!!
         
         !!--------flux terms--------!!
         volume_flux = 0

         do i=1,_SWE_DG_ORDER+1
            volume_flux(i,:,:) = matmul(jacobian_inv(1,1) * s_m_inv_s_k_x_t+&
                                        jacobian_inv(2,1) * s_m_inv_s_k_y_t, f_ref(1,i,:,:)) +&
                                 matmul(jacobian_inv(1,2) * s_m_inv_s_k_x_t+&
                                        jacobian_inv(2,2) * s_m_inv_s_k_y_t, f_ref(2,i,:,:))
         end do

         do j=1,_SWE_DG_DOFS
            q_temp_st(:,j,:) = matmul(-t_k_t_11_inv_t_m_1, volume_flux(:,j,:))
         end do
         !!----- end flux terms -----!!
         
         !!------- source terms ------!!
         
         !!------------S1-------------!!
         source_st = 0
         source_st(:,:,2) = s_ref(:,:,2,1) * jacobian(1,1) + s_ref(:,:,3,1) * jacobian(1,2)
         source_st(:,:,3) = s_ref(:,:,2,1) * jacobian(2,1) + s_ref(:,:,3,1) * jacobian(2,2)

         do j=1,_SWE_DG_DOFS
            q_temp_st(:,j,:) = q_temp_st(:,i,:) + matmul(t_k_t_11_inv_t_m_1, source_st(:,j,:))
         end do

         !!------------S2-------------!!
         source_ref_st = 0
         do i=1,_SWE_DG_ORDER + 1
            source_ref_st(i,:,2) =  matmul(s_m_inv_s_k_x_t, s_ref(i,:,2,2))
            source_ref_st(i,:,3) =  matmul(s_m_inv_s_k_y_t, s_ref(i,:,3,2))
         end do

         source_st = 0
         source_st(:,:,2) = source_ref_st(:,:,2) * jacobian(1,1) + source_ref_st(:,:,3) * jacobian(1,2)
         source_st(:,:,3) = source_ref_st(:,:,2) * jacobian(2,1) + source_ref_st(:,:,3) * jacobian(2,2)

         do j=1,_SWE_DG_DOFS
            q_temp_st(:,j,:) = q_temp_st(:,i,:) + matmul(t_k_t_11_inv_t_m_1, source_st(:,j,:))
         end do
         !!---- end source terms ----!!

         q_temp_st = q_temp_st * dt/dx
         
         do i=1,_SWE_DG_ORDER
            do j=1,_SWE_DG_DOFS         
               q_temp_st(i,j,:) = q_temp_st(i,j,:) - t_k_t_11_inv_x_t_k_t_10(i,1) * q_0(j,:)
            end do
         end do

         !------ compute error ------!
         epsilon=0.0_GRID_SR
         do i=1,_SWE_DG_ORDER
            do j=1,_SWE_DG_DOFS
               if(.not.all(q_i_st(i+1,j,:).eq.0)) then
                  epsilon=max(epsilon, &
                       NORM2(q_temp_st(i,j,:)-q_i_st(i+1,j,:)) / &
                       NORM2(q_i_st(i+1,j,:)))
               else if(.not.all(q_temp_st(i,j,:) .eq. 0)) then
                  epsilon=max(epsilon,1.0_GRID_SR)
               end if
            end do
         end do
         !--------------------------!

         !--------------Update predictor-------------!
         do i=2,_SWE_DG_ORDER+1
            q_i_st(i,:,:) = q_temp_st(i-1,:,:)
         end do
         !-------------------------------------------!

         !------Guard for diverging Picard Loop------!
         if(iteration > cfg%max_picard_iterations) then                           
            !print*,"predictor not converging"
            !------ if predictor diverges continue with fv cell ------!
            call apply_phi(Q_DG(:)%h+Q_DG(:)%b     ,cell%data_pers%H)
            call apply_phi(Q_DG(:)%p(1)            ,cell%data_pers%HU)
            call apply_phi(Q_DG(:)%p(2)            ,cell%data_pers%HV)
            call apply_phi(Q_DG(:)%b               ,cell%data_pers%B)
            cell%data_pers%troubled=PREDICTOR_DIVERGED
            !---------------------------------------------------------!
            exit
         end if
         !-------------------------------------------!
      end do

      do i=0,_SWE_DG_ORDER
         offset   = _SWE_DG_DOFS * i
         q_i(offset+1:offset+_SWE_DG_DOFS,:) = q_i_st(i+1,:,:)
      end do

      !!--- compute volume update----!!
      !-- Question: is it better to recompte the flux then to store it --!
      f(1,:,:,:) = jacobian_inv(1,1) * f_ref(1,:,:,:) + jacobian_inv(1,2) * f_ref(2,:,:,:)
      f(2,:,:,:) = jacobian_inv(2,1) * f_ref(1,:,:,:) + jacobian_inv(2,2) * f_ref(2,:,:,:)

      do i=1,_SWE_DG_ORDER+1
         volume_flux(i,:,:) = matmul(s_m_inv, &
                              matmul(s_k_x_s_b_3_s_b_2, f(1,i,:,:)) + &
                              matmul(s_k_y_s_b_1_s_b_2, f(2,i,:,:)))
      end do
      
      do i=1,_SWE_DG_ORDER+1
         source_ref_st(i,:,2) = matmul(s_m_inv, -matmul(s_m, s_ref(i,:,2,1))  - matmul(s_k_x_s_b_3_s_b_2, s_ref(i,:,2,2)))
         source_ref_st(i,:,3) = matmul(s_m_inv, -matmul(s_m, s_ref(i,:,3,1))  - matmul(s_k_y_s_b_1_s_b_2, s_ref(i,:,3,2)))
      end do

      source_st(:,:,2) = source_ref_st(:,:,2) * jacobian(1,1) + source_ref_st(:,:,3) * jacobian(1,2)
      source_st(:,:,3) = source_ref_st(:,:,2) * jacobian(2,1) + source_ref_st(:,:,3) * jacobian(2,2)

      volume_flux = volume_flux + source_st
      
      Q_DG_UPDATE(:,1) = reshape(matmul(t_a,volume_flux(:,:,1)),(/_SWE_DG_DOFS/))
      Q_DG_UPDATE(:,2) = reshape(matmul(t_a,volume_flux(:,:,2)),(/_SWE_DG_DOFS/))
      Q_DG_UPDATE(:,3) = reshape(matmul(t_a,volume_flux(:,:,3)),(/_SWE_DG_DOFS/))
      !!------------------------------!!
      
      !!---- set values for riemannsolve and project on edges----!!
      do i = 1,_SWE_DG_DOFS
         QP(  i,1:3) = reshape( matmul(t_a,q_i_st(:,i,:)) ,(/ 3 /))
         FP(1,i, : ) = reshape( matmul(t_a,f_ref(1,:,i,:)),(/ 3 /))
         FP(2,i, : ) = reshape( matmul(t_a,f_ref(2,:,i,:)),(/ 3 /))
      end do
      QP(:,4) = Q_DG(:)%B
      
      do j = 1,3
         if(cell%geometry%i_plotter_type < 0)then
            edge_type = -j
         else
            edge_type =  j
         end if
         do i = 1,_SWE_DG_ORDER+1
            select case(edge_type)
            case(-3 ,1) !right
               indx=_SWE_DG_DOFS-(_SWE_DG_ORDER-i+2)*(_SWE_DG_ORDER-i+3)/2 +1
            case(-2, 2) !mid
               indx=_SWE_DG_DOFS-(_SWE_DG_ORDER-i)*(_SWE_DG_ORDER-i+1)/2
            case(-1 ,3) !left
               indx=1+i
            case default
               stop
            end select
            edges(j)%ptr%data_pers%QP(i,:)     = QP(indx,:)
            edges(j)%ptr%data_pers%FP(1,i,:)   = FP(1,indx,:)
            edges(j)%ptr%data_pers%FP(2,i,:)   = FP(2,indx,:)        
         end do
         select case(edge_type)
         case (-1,3) !cells with id i*i+1 (left leg)
            edges(j)%ptr%data_pers%H  = matmul(phi_l,Q_DG%h) 
            edges(j)%ptr%data_pers%HU = matmul(phi_l,Q_DG%p(1))
            edges(j)%ptr%data_pers%HV = matmul(phi_l,Q_DG%p(2)) 
            edges(j)%ptr%data_pers%B  = matmul(phi_l,Q_DG%b)
         case (-2,2) ! hypotenuse
            edges(j)%ptr%data_pers%H  = matmul(phi_m,Q_DG%h)
            edges(j)%ptr%data_pers%HU = matmul(phi_m,Q_DG%p(1))
            edges(j)%ptr%data_pers%HV = matmul(phi_m,Q_DG%p(2))
            edges(j)%ptr%data_pers%B  = matmul(phi_m,Q_DG%b)
         case (-3,1) !cells with id i*i (right leg)
            edges(j)%ptr%data_pers%H  = matmul(phi_r,Q_DG%h)
            edges(j)%ptr%data_pers%HU = matmul(phi_r,Q_DG%p(1))
            edges(j)%ptr%data_pers%HV = matmul(phi_r,Q_DG%p(2))
            edges(j)%ptr%data_pers%B  = matmul(phi_r,Q_DG%b)
         end select
      end do

      !!---- updateQ ----!!
      Q_DG%h    = Q_DG%h    + Q_DG_UPDATE(:,1) * dt/dx
      Q_DG%p(1) = Q_DG%p(1) + Q_DG_UPDATE(:,2) * dt/dx
      Q_DG%p(2) = Q_DG%p(2) + Q_DG_UPDATE(:,3) * dt/dx
      
    end associate
  end subroutine dg_predictor  

function flux(q,N)
  real(kind=GRID_SR)             :: flux(2,N,3)
  real(kind=GRID_SR), intent(in) :: q(N,3)
  integer                        :: N
  
  flux(1,:,1) = q(:,2)
  flux(1,:,2) = q(:,2)**2/q(:,1) + 0.5_GRID_SR * g * q(:,1)**2
  flux(1,:,3) = q(:,2)*q(:,3)/q(:,1)
  
  flux(2,:,1) = q(:,3)
  flux(2,:,2) = q(:,2)*q(:,3)/q(:,1)
  flux(2,:,3) = q(:,3)**2/q(:,1) + 0.5_GRID_SR * g * q(:,1)**2
  
end function flux
  
function flux_1(q,N)
  real(kind=GRID_SR)             ::flux_1(N,3)
  real(kind=GRID_SR) ,intent(in) ::q(N,3)
  integer :: N
  flux_1(:,1) = q(:,2)
  flux_1(:,2) = q(:,2)**2/q(:,1) + 0.5_GRID_SR * g * q(:,1)**2
  flux_1(:,3) = q(:,2)*q(:,3)/q(:,1)
end function flux_1

function flux_2(q,N)
  real(kind=GRID_SR)             ::flux_2(N,3)
  real(kind=GRID_SR) ,intent(in) ::q(N,3)
  integer :: N
  flux_2(:,1) = q(:,3)
  flux_2(:,2) = q(:,2)*q(:,3)/q(:,1)
  flux_2(:,3) = q(:,3)**2/q(:,1) + 0.5_GRID_SR * g * q(:,1)**2
end function flux_2

END MODULE SWE_DG_predictor

#endif
