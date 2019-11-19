! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2017 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE
! Author: Leonhard Rannabauer rannabau (at) in.tum.de

#include "Compilation_control.f90"

#if defined(_SWE_DG)

MODULE SWE_dg_predictor
  use SFC_edge_traversal
  use Samoa_swe
  implicit none
  
  public dg_predictor,flux_1,flux_2
  
contains
  
  subroutine dg_predictor_old(cell,dt)
    class(t_cell_data_ptr), intent(in)         :: cell
    real(kind=GRID_SR) ,intent(in)             :: dt
    real(kind=GRID_SR)                         :: dx
    real(kind=GRID_SR)                         :: q_0(_SWE_DG_DOFS,3)
    real(kind=GRID_SR)                         :: q_i(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3)
    real(kind=GRID_SR)                         :: q_temp(_SWE_DG_DOFS*_SWE_DG_ORDER,3)
    
    !--local variables--!
    integer                                    :: iteration,i,j,offset
    real(kind=GRID_SR),Dimension(2,2)          :: jacobian,jacobian_inv
    real(kind=GRID_SR),Dimension(3)            :: edge_sizes
    real(kind=GRID_SR)                         :: epsilon=1.0_GRID_SR
    
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3)  :: f1, f1_ref, f2, f2_ref
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*_SWE_DG_ORDER,3)      :: source,source_ref
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*_SWE_DG_ORDER,3)      :: volume_flux1,volume_flux2,volume_flux

    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1))    :: H_x,H_y
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS)    :: b,H
    associate(Q_DG =>cell%data_pers%Q_DG ,Q_DG_P =>cell%data_pers%Q_DG_P)
      
      !--Normalize jacobian--!
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

      print*,"q_i"
      do j=1,(_SWE_DG_ORDER)*_SWE_DG_DOFS
         print*,q_i(j,:)
      end do
      
      i=0
      epsilon=1.0_GRID_SR
      do while(epsilon > 1.0e-14_GRID_SR)
         i=i+1
         
         if (any(q_i(:,1).le.0)) then
#if defined(_DEBUG)  
            print*,"PRED"
            print*,q_i
            print*,"FV"
            print*,cell%data_pers%H
            write(*,*), "ERROR: Waterheigth less than zero in predictor"
#endif
            exit
         end if

         !--- Flux on ref
         f1_ref = flux_1(q_i,(_SWE_DG_ORDER+1)*_SWE_DG_DOFS)
         f2_ref = flux_2(q_i,(_SWE_DG_ORDER+1)*_SWE_DG_DOFS)
         
         !--spatial derivations needed for well-balanced scheme--!
         do j=0,(_SWE_DG_ORDER)
            offset=_SWE_DG_DOFS * j
            H=q_i(offset+1:offset+_SWE_DG_DOFS,1) + b
            H_x(offset+1:offset+_SWE_DG_DOFS)= matmul(basis_der_x,H)
            H_y(offset+1:offset+_SWE_DG_DOFS)= matmul(basis_der_y,H)
         end do
         
         !---Source on reference element---!
         source_ref=0
         source_ref(:,2)=matmul(st_k_x,(0.5_GRID_SR*g*q_i(:,1)**2))
         source_ref(:,3)=matmul(st_k_y,(0.5_GRID_SR*g*q_i(:,1)**2))

         source_ref(:,2)=source_ref(:,2)-matmul(st_m,(g*q_i(:,1)*H_x))
         source_ref(:,3)=source_ref(:,3)-matmul(st_m,(g*q_i(:,1)*H_y))
         
         !---Transform back to element ---!
         source=0
         source(:,2) = jacobian(1,1) * source_ref(:,2) + jacobian(1,2) * source_ref(:,3)
         source(:,3) = jacobian(2,1) * source_ref(:,2) + jacobian(2,2) * source_ref(:,3)
         
         !---Transform flux back to element---!
         f1=(jacobian_inv(1,1)*f1_ref+jacobian_inv(1,2)*f2_ref)
         f2=(jacobian_inv(2,1)*f1_ref+jacobian_inv(2,2)*f2_ref)

         source  = source*dt/dx
         volume_flux1 =matmul(st_k_x,f1)*dt/dx
         volume_flux2 =matmul(st_k_y,f2)*dt/dx

         volume_flux = volume_flux1 + volume_flux2
         
         print*,"volume_flux 1"
         do j=1,(_SWE_DG_ORDER)*_SWE_DG_DOFS
            print*,volume_flux(j,:)
         end do


         print*,"volume_flux"
         volume_flux = matmul(s_m,volume_flux1 + volume_flux2)
         
         do j=1,(_SWE_DG_ORDER)*_SWE_DG_DOFS
            print*,volume_flux(j,:)
         end do
         
         !TODO only calculate new dofs
         q_temp= source &
              - matmul(st_w_k_t_1_0,q_i(1:_SWE_DG_DOFS,:))&
              - volume_flux1 &
              - volume_flux2

         q_temp(:,1)=matmul(st_w_k_t_1_1_inv,q_temp(:,1))
         q_temp(:,2)=matmul(st_w_k_t_1_1_inv,q_temp(:,2))
         q_temp(:,3)=matmul(st_w_k_t_1_1_inv,q_temp(:,3))

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
            cell%data_pers%troubled=4
            exit
         end if
         
         q_i(1+_SWE_DG_DOFS:(_SWE_DG_ORDER+1)*_SWE_DG_DOFS,:) = q_temp
         
         print*,"q_i"
         do j=1,(_SWE_DG_ORDER)*_SWE_DG_DOFS
            print*,q_i(j,:)
         end do
      end do

      call cell%data_pers%set_dofs_pred(q_i)
      
    end associate
  end subroutine dg_predictor_old

    subroutine dg_predictor(cell,dt)
    class(t_cell_data_ptr), intent(in)         :: cell
    real(kind=GRID_SR) ,intent(in)             :: dt
    real(kind=GRID_SR)                         :: dx
    real(kind=GRID_SR)                         :: q_0(_SWE_DG_DOFS,3)
    real(kind=GRID_SR)                         :: q_i(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3)
    real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,_SWE_DG_DOFS,3)  :: q_i_st
    real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,_SWE_DG_DOFS,3)  :: source_st,source_ref_st
    real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,_SWE_DG_DOFS,3)  :: volume_flux
    real(kind=GRID_SR),Dimension(2,_SWE_DG_ORDER+1,_SWE_DG_DOFS,3):: f,f_ref
    real(kind=GRID_SR),Dimension(_SWE_DG_ORDER+1,_SWE_DG_DOFS)    :: H_x_st,H_y_st
    real(kind=GRID_SR),Dimension(_SWE_DG_ORDER,_SWE_DG_DOFS,3)    :: q_temp_st
    
    !--local variables--!
    integer                                    :: iteration,i,j,offset
    real(kind=GRID_SR),Dimension(2,2)          :: jacobian,jacobian_inv
    real(kind=GRID_SR),Dimension(3)            :: edge_sizes
    real(kind=GRID_SR)                         :: epsilon=1.0_GRID_SR
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS) :: b

    real(kind=GRID_SR),Dimension(_SWE_DG_ORDER,1) :: t_k_t_11_inv_x_t_k_t_10 

    associate(Q_DG =>cell%data_pers%Q_DG ,Q_DG_P =>cell%data_pers%Q_DG_P)

      ! TODO make this a precompiled matrix
      t_k_t_11_inv_x_t_k_t_10=matmul(t_k_t_11_inv,t_k_t_10)
      
      !--Normalize jacobian--!
      jacobian=ref_plotter_data(abs(cell%geometry%i_plotter_type))%jacobian_normalized
      jacobian_inv=ref_plotter_data(abs(cell%geometry%i_plotter_type))%jacobian_inv_normalized
      
      edge_sizes=cell%geometry%get_edge_sizes()
      dx=edge_sizes(1)*cfg%scaling
      
      call cell%data_pers%get_dofs_dg(q_0)
      
      !--Set initial conditions for discrete picard iteration--!
      b=Q_DG(:)%B

      !span over time basis
      do i=1,_SWE_DG_ORDER+1
         q_i_st(i,:,:) = q_0
      end do

      print*,"q_i"
      do i=1,(_SWE_DG_ORDER)
         do j=1,_SWE_DG_DOFS
            print*,q_i_st(i,j,:)
         end do
      end do
      
      iteration=0
      epsilon=1.0_GRID_SR
      do while(epsilon > 1.0e-14_GRID_SR)
         iteration=iteration+1
         
         if (any(q_i_st(:,:,1).le.0)) then
#if defined(_DEBUG)  
            print*,"PRED"
            print*,q_i_st
            print*,"FV"
            print*,cell%data_pers%H
            write(*,*), "ERROR: Waterheigth less than zero in predictor"
#endif
            exit
         end if
        
         !!--spatial derivations needed for well-balanced scheme--!!
         do i=1,_SWE_DG_ORDER+1
            H_x_st(i,:) = matmul(basis_der_x,q_i_st(i,:,1) + b)
            H_y_st(i,:) = matmul(basis_der_y,q_i_st(i,:,1) + b)
         end do
         !!---------------------------!!

         !!------- source terms ------!!
         source_ref_st = 0
         do j=1,_SWE_DG_DOFS
            source_ref_st(:,j,2) = matmul(t_m, (0.5_GRID_SR * g * q_i_st(:,j,1)**2))
            source_ref_st(:,j,3) = matmul(t_m, (0.5_GRID_SR * g * q_i_st(:,j,1)**2))
         end do
         
         do i=1,_SWE_DG_ORDER+1
            source_ref_st(i,:,2) = matmul(s_m_inv, matmul(transpose(s_k_x), source_ref_st(i,:,2)))
            source_ref_st(i,:,3) = matmul(s_m_inv, matmul(transpose(s_k_y), source_ref_st(i,:,3)))
         end do

         do j=1,_SWE_DG_DOFS
            source_ref_st(:,j,2) = source_ref_st(:,j,2) - matmul(t_m,(g * q_i_st(:,j,1) * H_x_st(:,j)))
            source_ref_st(:,j,3) = source_ref_st(:,j,3) - matmul(t_m,(g * q_i_st(:,j,1) * H_y_st(:,j)))
         end do

         source_st(:,:,1) = 0
         source_st(:,:,2) = jacobian(1,1) * source_ref_st(:,:,2) + jacobian(1,2) * source_ref_st(:,:,3)
         source_st(:,:,3) = jacobian(2,1) * source_ref_st(:,:,2) + jacobian(2,2) * source_ref_st(:,:,3)
         !!---- end source terms ----!!
         
         !!--------flux terms--------!!
         !--- Flux on ref element---!
         f_ref = 0
         do i=1,_SWE_DG_ORDER+1
            f_ref(:,i,:,:) = flux(q_i_st(i,:,:),_SWE_DG_DOFS)
         end do

         f = 0
         !---Transform flux back to local element---!
         f(1,:,:,:) = jacobian_inv(1,1) * f_ref(1,:,:,:) + jacobian_inv(1,2) * f_ref(2,:,:,:)
         f(2,:,:,:) = jacobian_inv(2,1) * f_ref(1,:,:,:) + jacobian_inv(2,2) * f_ref(2,:,:,:)

         volume_flux = 0

         !---evaluate spatial derivative---!

         do i=1,_SWE_DG_ORDER+1
            volume_flux(i,:,:) = matmul(transpose(s_k_x),f(1,i,:,:)) &
                               + matmul(transpose(s_k_y),f(2,i,:,:))
         end do
   
         do j=1,_SWE_DG_DOFS
            volume_flux(:,j,:) =  matmul(t_m, volume_flux(:,j,:))
         end do

         print*,"volume_flux1"
         do i=1,_SWE_DG_ORDER+1
            do j=1,_SWE_DG_DOFS
               print*,volume_flux(i,j,:)*dt/dx
            end do
         end do

         do i=1,_SWE_DG_ORDER+1
            volume_flux(i,:,:) = matmul(s_m_inv,volume_flux(i,:,:))
         end do

         print*,"volume_flux"
         do i=1,_SWE_DG_ORDER+1
            do j=1,_SWE_DG_DOFS
               print*,volume_flux(i,j,:)*dt/dx
            end do
         end do

         !!----- end flux terms -----!!
         

         !---- add flux ----!
         !LR: get rid of that one colume thats to much
         do i=2,_SWE_DG_ORDER+1
            q_temp_st(i-1,:,:) = source_st(i,:,:) - volume_flux(i,:,:)
         end do
         
         do i=1,_SWE_DG_DOFS         
            q_temp_st(:,i,:) = matmul(t_k_t_11_inv,q_temp_st(:,i,:))
         end do

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

         !-----Guard for diverging Picard Loop-------!
         if(iteration > 200) then                           
            print*,"predictor not converging"
            cell%data_pers%troubled=4
            exit
         end if

         !--------------Update predictor-------------!
         do i=2,_SWE_DG_ORDER+1
            q_i_st(i,:,:) = q_temp_st(i-1,:,:)
         end do
         !-------------------------------------------!
         print*,"q_i"
         do i=1,(_SWE_DG_ORDER)
            do j=1,_SWE_DG_DOFS
               print*,q_i_st(i,j,:)
            end do
         end do
         
      end do

      do i=0,_SWE_DG_ORDER
         offset   = _SWE_DG_DOFS * i
         q_i(offset+1:offset+_SWE_DG_DOFS,:) = q_i_st(i+1,:,:)
      end do
      
      call cell%data_pers%set_dofs_pred(q_i)
      
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
