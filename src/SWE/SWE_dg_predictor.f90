! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2017 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE
! Author: Leonhard Rannabauer rannabau (at) in.tum.de

#include "Compilation_control.f90"

#if defined(_SWE_DG)

MODULE SWE_DG_predictor
  use SFC_edge_traversal
  use Samoa_swe
  use SWE_DG_Limiter
  use SWE_dg_solver
#if defined(_OPT_KERNELS)  
  use yateto_interface
#endif

#if defined(CHAMELEON)
   use Chameleon_lib
#define _GT_USE_CHAMELEON
#endif
#if defined(CHAMELEON_CALL)
#  define _GT_USE_CHAMELEON_CALL
#endif

  implicit none
  
  type num_traversal_data     
  end type num_traversal_data

# define _GT_NAME	              		t_swe_dg_predictor_traversal
# define _GT_EDGES
!# define _GT_EDGE_MPI_TYPE

# define _GT_ELEMENT_OP             element_op
# define _GT_CELL_TO_EDGE_OP        cell_to_edge_op_dg  
  
  public dg_predictor!,flux_1,flux_2
  
  public writeFVBoundaryFields
#		include "SFC_generic_traversal_ringbuffer.f90"
  
  subroutine element_op(traversal, section, element)
    type(t_swe_dg_predictor_traversal), intent(inout) :: traversal
    type(t_grid_section), intent(inout)		      :: section
    type(t_element_base), intent(inout)		      :: element
    integer(kind = selected_int_kind(16))             :: iterations

    iterations = 0
    !-------Compute predictor --------!
    if(isDG(element%cell%data_pers%troubled)) then
       call dg_predictor    (element,section%r_dt,iterations)
    end if

    call writeFVBoundaryFields(element)
    call traversal%stats%add_counter(predictor_iterations, iterations)
    !---------------------------------!
  end subroutine element_op

  subroutine dg_predictor(element,dt,iterations)
    class(t_element_base), intent(inout)         :: element
    real(kind=GRID_SR) ,intent(in)               :: dt
    real(kind=GRID_SR)                           :: dx
    integer(kind = selected_int_kind(16)),intent(inout) :: iterations

    !local variables
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3) :: q_0
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,3)  :: q_i_st

    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,2,3) :: f_ref
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,  2) :: s_ref
    
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER,3) :: q_temp_st
    real(kind=GRID_SR)   , DIMENSION(_SWE_DG_DOFS,3)           :: q_dg_update_t

#if defined(_OPT_KERNELS)
    !dir$ attributes align:ALIGNMENT :: q_0
    !dir$ attributes align:ALIGNMENT :: q_i_st
    !dir$ attributes align:ALIGNMENT :: f_ref
    !dir$ attributes align:ALIGNMENT :: s_ref
    !dir$ attributes align:ALIGNMENT :: q_temp_st
    !dir$ attributes align:ALIGNMENT :: dtdx
    !dir$ attributes align:ALIGNMENT :: q_dg_update_t
#endif    
    
    integer                                    :: cell_type
    integer                                    :: i,j
    real(kind=GRID_SR),Dimension(1)            :: dtdx
    
    !--local variables--!
    integer                                    :: iteration
    real(kind=GRID_SR)                         :: epsilon


    associate(Q_DG        => element%cell%data_pers%Q,&
              Q_DG_UPDATE => element%cell%data_pers%Q_DG_UPDATE,&
              cell  => element%cell,&
              edges => element%edges)

      cell_type = abs(cell%geometry%i_plotter_type)

      !!--Set initial conditions for discrete picard iteration--!!
      call initialise_predictor(Q_DG, cell, dt, dtdx, iteration, epsilon , q_0 , q_i_st)
      !!---------------------------!!
      
      do while(epsilon > cfg%max_picard_error .and.&
           (.not.(cell%data_pers%troubled.eq.PREDICTOR_DIVERGED)))

         iteration  = iteration+1
         iterations = iteration
         
         !!--- Compute F S1 and S2---!!
         call compute_sref(s_ref,q_i_st(:,:,1),Q_DG%B)
         call compute_fref(f_ref,q_i_st)
         
#if defined(_DEBUG)
         do i=1,_SWE_DG_ORDER+1
            do j=1,_SWE_DG_DOFS
               if(isnan(f_ref(j,i,1,1)).or.isnan(f_ref(j,i,2,1))) then
                  print*,"fref"
                  print*,f_ref
                  exit
               end if
            end do
         end do
#endif


#if defined(_OPT_KERNELS)
         q_temp_st = 0
         call yateto_predictor_execute(q_temp_st, f_ref, q_0, s_ref, dtdx , cell_type-1)
#else         
         call update_predictor(q_temp_st, q_0, dtdx, s_ref, f_ref, cell_type)
#endif         
         !------ compute error ------!
         ! print*,q_i_st
         ! print*,q_temp_st
         epsilon =  compute_epsilon(q_i_st,q_temp_st)
         !--------------------------!

         !--------------Update predictor-------------!
         do i=2,_SWE_DG_ORDER+1
            q_i_st(:,i,:) = q_temp_st(:,i-1,:)
         end do
         !-------------------------------------------!

         !------Guard for diverging Picard Loop------!
         if(iteration > cfg%max_picard_iterations) then                           
            cell%data_pers%troubled=PREDICTOR_DIVERGED
         end if
         !-------------------------------------------!

          !------Guard for negative water height------!
         if (any(q_i_st(:,:,1).le.0)) then
            cell%data_pers%troubled=PREDICTOR_DIVERGED
         end if
         !-------------------------------------------!

      end do

      if(.not.(cell%data_pers%troubled.eq.PREDICTOR_DIVERGED)) then
#if defined(_OPT_KERNELS)
         call yateto_volume_execute(q_dg_update_t, f_ref, s_ref, cell_type - 1)
         Q_DG_UPDATE = q_dg_update_t
#else
         call compute_volume_update(Q_DG_UPDATE, s_ref, f_ref, cell_type)
#endif         

         call initialise_riemann_arguments(cell,q_i_st,Q_DG)

      end if
    end associate
  end subroutine dg_predictor
  
  subroutine writeFVBoundaryFields(element)
    class(t_element_base), intent(inout)       :: element
    integer                                    :: i,j,edge
    if(.not.isCoast(element%cell%data_pers%troubled)) then
       associate(Q_DG => element%cell%data_pers%Q,&
            QFV  => element%cell%data_pers%QFV)
         do edge = 1,3
            select case(edge)
            case (1) !right
               QFV(edge,:,1) = matmul(phi_r,Q_DG%h+Q_DG%b)
               QFV(edge,:,2) = matmul(phi_r,Q_DG%p(1))
               QFV(edge,:,3) = matmul(phi_r,Q_DG%p(2))
               QFV(edge,:,4) = matmul(phi_r,Q_DG%b)
            case (2) !mid
               QFV(edge,:,1) = matmul(phi_m,Q_DG%h+Q_DG%b)
               QFV(edge,:,2) = matmul(phi_m,Q_DG%p(1))
               QFV(edge,:,3) = matmul(phi_m,Q_DG%p(2))
               QFV(edge,:,4) = matmul(phi_m,Q_DG%b)
            case (3) !left
               QFV(edge,:,1) = matmul(phi_l,Q_DG%h+Q_DG%b) 
               QFV(edge,:,2) = matmul(phi_l,Q_DG%p(1))    
               QFV(edge,:,3) = matmul(phi_l,Q_DG%p(2))    
               QFV(edge,:,4) = matmul(phi_l,Q_DG%b)
            end select
         end do
       end associate

    else
       associate(data => element%cell%data_pers,&
                 QFV  => element%cell%data_pers%QFV)
         do edge = 1,3
            select case (edge)       
            case (3) !cells with id i*i+1 (left leg)
               do i=0, _SWE_PATCH_ORDER - 1
                  j=_SWE_PATCH_ORDER-1-i
                  QFV(edge,i+1,1) = data%H(j*j + 1)
                  QFV(edge,i+1,2) = data%HU(j*j + 1)
                  QFV(edge,i+1,3) = data%HV(j*j + 1)
                  QFV(edge,i+1,4) = data%B(j*j + 1)
               end do
            case (2) ! hypotenuse
               do i=1, _SWE_PATCH_ORDER
                  j=_SWE_PATCH_ORDER+1-i
                  QFV(edge,i,1) = data%H ((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
                  QFV(edge,i,2) = data%HU((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
                  QFV(edge,i,3) = data%HV((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
                  QFV(edge,i,4) = data%B ((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*j - 1)
               end do
            case (1) !cells with id i*i (right leg)
               do i=1, _SWE_PATCH_ORDER
                  QFV(edge,i,1) = data%H(i*i)
                  QFV(edge,i,2) = data%HU(i*i)
                  QFV(edge,i,3) = data%HV(i*i)
                  QFV(edge,i,4) = data%B(i*i)
               end do
            end select
         end do
       end associate
    end if
  end subroutine writeFVBoundaryFields

function flux(q,N)
  real(kind=GRID_SR)             :: flux(N,2,3)
  real(kind=GRID_SR), intent(in) :: q(N,3)
  integer                        :: N
  
  flux(:,1,1) = q(:,2)
  flux(:,1,2) = q(:,2)**2/q(:,1) + 0.5_GRID_SR * g * q(:,1)**2
  flux(:,1,3) = q(:,2)*q(:,3)/q(:,1)
  
  flux(:,2,1) = q(:,3)
  flux(:,2,2) = q(:,2)*q(:,3)/q(:,1)
  flux(:,2,3) = q(:,3)**2/q(:,1) + 0.5_GRID_SR * g * q(:,1)**2
  
end function flux
  
function flux_1(q,N)
  real(kind=GRID_SR)             ::flux_1(N,3)
  real(kind=GRID_SR) ,intent(in) ::q(N,3)
  integer :: N
  flux_1(:,1) = q(:,2)
  flux_1(:,2) = q(:,2)**2/q(:,1) + 0.5_GRID_SR * g * q(:,1)**2
  !flux_1(:,2) = q(:,2)**2/q(:,1)
  flux_1(:,3) = q(:,2)*q(:,3)/q(:,1)
end function flux_1

function flux_2(q,N)
  real(kind=GRID_SR)             ::flux_2(N,3)
  real(kind=GRID_SR) ,intent(in) ::q(N,3)
  integer :: N
  flux_2(:,1) = q(:,3)
  flux_2(:,2) = q(:,2)*q(:,3)/q(:,1)
  flux_2(:,3) = q(:,3)**2/q(:,1) + 0.5_GRID_SR * g * q(:,1)**2
  !flux_2(:,3) = q(:,3)**2/q(:,1)
end function flux_2

function flux_no_grav(q,N)
  real(kind=GRID_SR)             ::flux_no_grav(N,2,3)
  real(kind=GRID_SR) ,intent(in) ::q(N,3)
  integer :: N

  flux_no_grav(:,1,1) = q(:,2)
  flux_no_grav(:,1,2) = q(:,2)**2/q(:,1)
  flux_no_grav(:,1,3) = q(:,2)*q(:,3)/q(:,1)

  flux_no_grav(:,2,1) = q(:,3)
  flux_no_grav(:,2,2) = q(:,2)*q(:,3)/q(:,1)
  flux_no_grav(:,2,3) = q(:,3)**2/q(:,1)

end function flux_no_grav

function compute_epsilon(q_i_st,q_temp_st) result(epsilon)
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,3)  :: q_i_st
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER,3)    :: q_temp_st
  real(kind=GRID_SR)                         :: epsilon
  integer i,j
  
  epsilon=0.0_GRID_SR
  do j=1,_SWE_DG_DOFS
     do i=1,_SWE_DG_ORDER
        if(.not.all(q_i_st(j,i+1,:).eq.0)) then
           epsilon=max(epsilon, &
                NORM2(q_temp_st(j,i,:)-q_i_st(j,i+1,:)) / &
                NORM2(q_i_st(j,i+1,:)))
        else if(.not.all(q_temp_st(j,i,:) .eq. 0)) then
           epsilon=max(epsilon,1.0_GRID_SR)
        end if
     end do
  end do
end function compute_epsilon

subroutine initialise_predictor(Q_DG, cell, dt, dtdx, iteration, epsilon, q_0 , q_i_st)
  type(t_state)        , DIMENSION(_SWE_DG_DOFS),intent(in) :: Q_DG
  type(t_cell_data_ptr),intent(inout)          :: cell
  real(kind=GRID_SR),intent(in)                :: dt
  real(kind=GRID_SR),intent(out),dimension(1)  :: dtdx
  integer,intent(inout)                        :: iteration
  real(kind=GRID_SR),intent(inout)             :: epsilon
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3) :: q_0
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,3)  :: q_i_st

  !local variables
  real(kind=GRID_SR),Dimension(3)             :: edge_sizes
  real(kind=GRID_SR)                          :: dx
  integer :: i

  !--transform Q_DG to tensor--!
  !--TODO: Should be default representation--!
  q_0(:,1) = Q_DG%h
  q_0(:,2) = Q_DG%p(1)
  q_0(:,3) = Q_DG%p(2)
  !----!

  !------Guard for negative water height------!
  if (any(q_0(:,1).le.0)) then
     cell%data_pers%troubled=PREDICTOR_DIVERGED
  end if
  !-------------------------------------------!
  
  !--span dofs at t=0 over time basis--!
  do i=1,_SWE_DG_ORDER+1
     q_i_st(:,i,:) = q_0
  end do
  
  edge_sizes=cell%geometry%get_edge_sizes()
  dx=edge_sizes(1)*cfg%scaling
  dtdx(1) = dt/dx
  
  iteration=0
  epsilon=1.0_GRID_SR

end subroutine initialise_predictor

subroutine compute_sref(s_ref,h,b)
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,  2),intent(out) :: s_ref
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1),intent(in)      :: h
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS),intent(in)                      :: b
  !local variables
  integer :: i,j,k

#if defined(_OPT_KERNELS)
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS, _SWE_DG_ORDER+1, 2) :: h2
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS, _SWE_DG_ORDER+1)    :: w
  !dir$ attributes align:32 :: w
#endif
  
  s_ref = 0
  
#if defined(_OPT_KERNELS)

  do i=1,_SWE_DG_ORDER+1
     w(:,i) = h(:,i) + b
  end do
  
  call yateto_compute_source_execute(s_ref, w)
  
  do i=1,2
     do j=1,_SWE_DG_ORDER+1
        !$omp simd
        do k=1,_SWE_DG_DOFS
           s_ref(k,j,i) = g* h(k,j)*s_ref(k,j,i)
        end do
     end do
  end do
  
#else  

  do i=1,_SWE_DG_ORDER+1
     s_ref(:,i,1) = ( g * h(:,i) * matmul(basis_der_x,h(:,i) + b) )
     s_ref(:,i,2) = ( g * h(:,i) * matmul(basis_der_y,h(:,i) + b) )
  end do
  
#endif

  
end subroutine compute_sref

subroutine compute_fref(f_ref,q_i_st)
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,2,3),intent(out) :: f_ref
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,3),intent(in)    :: q_i_st
  integer :: i,j
  
  f_ref = 0
  do i=1,_SWE_DG_ORDER+1
     f_ref(:,i,:,:) = flux_no_grav(q_i_st(:,i,:),_SWE_DG_DOFS)
  end do
  
#if defined(_DEBUG)
  do i=1,_SWE_DG_ORDER+1
     do j=1,_SWE_DG_DOFS
        if(isnan(f_ref(j,i,1,1)).or.isnan(f_ref(j,i,2,1))) then
           print*,"fref"
           print*,f_ref
           exit
        end if
     end do
  end do
#endif
  
end subroutine compute_fref

subroutine update_predictor(q_temp_st, q_0 , dtdx, s_ref, f_ref, cell_type)
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER,3),intent(out)     :: q_temp_st
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3),intent(in)                    :: q_0
  real(kind=GRID_SR),dimension(1),intent(in)                                 :: dtdx
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,  2),intent(in)  :: s_ref
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,2,3),intent(in)  :: f_ref
  integer,intent(in)                                                         :: cell_type

  
  !local variables
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,3)  :: volume_flux
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,3)  :: source_st
  integer :: i,j

  volume_flux = 0
  
  do i=1,_SWE_DG_ORDER+1
     volume_flux(:,i,:) = matmul(jacobian_inv(cell_type,1,1) * s_m_inv_s_k_x_t+&
                                 jacobian_inv(cell_type,2,1) * s_m_inv_s_k_y_t, f_ref(:,i,1,:)) +&
                          matmul(jacobian_inv(cell_type,1,2) * s_m_inv_s_k_x_t+&
                                 jacobian_inv(cell_type,2,2) * s_m_inv_s_k_y_t, f_ref(:,i,2,:))
  end do
  
  do j=1,_SWE_DG_DOFS
     q_temp_st(j,:,:) = matmul(-t_k_t_11_inv_t_m_1, volume_flux(j,:,:))
  end do

  !!----- end flux terms -----!!
  
  !!------- source terms ------!!
  
  !!------------S1-------------!!
  source_st = 0
  source_st(:,:,2) = s_ref(:,:,1) * jacobian(1,1,cell_type) + s_ref(:,:,2) * jacobian(1,2,cell_type)
  source_st(:,:,3) = s_ref(:,:,1) * jacobian(2,1,cell_type) + s_ref(:,:,2) * jacobian(2,2,cell_type)
  
  do j=1,_SWE_DG_DOFS
     q_temp_st(j,:,:) = q_temp_st(j,:,:) + matmul(t_k_t_11_inv_t_m_1, source_st(j,:,:))
  end do
  !!------- end source term -------!!
  
  q_temp_st = q_temp_st * dtdx(1)
  
  do j=1,_SWE_DG_DOFS                  
     do i=1,_SWE_DG_ORDER
        q_temp_st(j,i,:) = q_temp_st(j,i,:) - t_k_t_11_inv_x_t_k_t_10(i,1) * q_0(j,:)
     end do
  end do
  
end subroutine update_predictor

subroutine compute_volume_update(Q_DG_UPDATE, s_ref, f_ref, cell_type)
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3)                               :: Q_DG_UPDATE
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,  2),intent(in)  :: s_ref
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,2,3),intent(in)  :: f_ref
  integer,intent(in)                                                         :: cell_type

  !local variables
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,2  ) :: s_ref_a
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,2,3) :: f_ref_a

  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3)   :: source_st, source_ref_st
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3)   :: volume_flux
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,2,3) :: f
  integer :: i,j

  volume_flux = 0
  source_st = 0
  f = 0
  
  do j = 1,_SWE_DG_DOFS
     f_ref_a(j,1,:) = reshape(matmul(t_a,f_ref(j,:,1,:)),(/ 3 /))
     f_ref_a(j,2,:) = reshape(matmul(t_a,f_ref(j,:,2,:)),(/ 3 /))
  end do

  do j = 1,_SWE_DG_DOFS
     s_ref_a(j,:) = reshape(matmul(t_a,s_ref(j,:,:)),(/ 2 /))
  end do

  
  f(:,1,:) = jacobian_inv(1,1,cell_type) * f_ref_a(:,1,:) +&
               jacobian_inv(1,2,cell_type) * f_ref_a(:,2,:)
  f(:,2,:) = jacobian_inv(2,1,cell_type) * f_ref_a(:,1,:) +&
               jacobian_inv(2,2,cell_type) * f_ref_a(:,2,:)
  
  volume_flux(:,:) = matmul(s_m_inv, &
       matmul(s_k_x_s_b_3_s_b_2, f(:,1,:)) + &
       matmul(s_k_y_s_b_1_s_b_2, f(:,2,:)))
  
  source_st(:,2) = -s_ref_a(:,1) * jacobian(1,1,cell_type)&
                   -s_ref_a(:,2) * jacobian(1,2,cell_type)
  source_st(:,3) = -s_ref_a(:,1) * jacobian(2,1,cell_type)&
                   -s_ref_a(:,2) * jacobian(2,2,cell_type)
  
#if defined(_DEBUG)
  do j=1,_SWE_DG_DOFS
     if(isnan(volume_flux(j,1)).or.isnan(source_st(j,1))) then
        print*,"vol"
        print*,volume_flux
        print*,"src"
        print*,source_st
        exit
     end if
  end do
#endif
  
  Q_DG_UPDATE = volume_flux + source_st

  !!------------------------------!!  
end subroutine compute_volume_update

subroutine initialise_riemann_arguments(cell, q_i_st, Q_DG)
  type(t_cell_data_ptr),intent(inout)                                     :: cell
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,3),intent(in) :: q_i_st
  type(t_state), DIMENSION(_SWE_DG_DOFS),intent(in)                       :: Q_DG
  
  !local variables
  real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,2,3)  :: f_ref
  real(kind=GRID_SR), DIMENSION(2,_SWE_DG_DOFS,3)  :: FP
  real(kind=GRID_SR), DIMENSION(  _SWE_DG_DOFS,4)  :: QP
  integer :: i,j,indx,edge_type 

  FP = 0
  QP = 0
  f_ref = 0
  
  !!---- set values for riemannsolve and project on edges----!!
  if(isDG(cell%data_pers%troubled)) then
     do i=1,_SWE_DG_ORDER+1
        f_ref(:,i,:,:) = flux(q_i_st(:,i,:),_SWE_DG_DOFS)
     end do
  endif
  
  do j = 1,_SWE_DG_DOFS
     QP(  j,1:3) = reshape( matmul(t_a,q_i_st(j,:,:))  ,(/ 3 /))
     FP(1,j, : ) = reshape( matmul(t_a,f_ref (j,:,1,:)),(/ 3 /))
     FP(2,j, : ) = reshape( matmul(t_a,f_ref (j,:,2,:)),(/ 3 /))
  end do
  QP(:,4) = Q_DG(:)%B

  do j = 1,3
     edge_type =  j
     do i = 1,_SWE_DG_ORDER+1
        select case(edge_type)
        case(1) !right
           indx=_SWE_DG_DOFS-(i+1)*(i)/2 +1
        case(2) !mid
           indx=_SWE_DG_DOFS-(_SWE_DG_ORDER-i+2)*(_SWE_DG_ORDER-i+1)/2
        case(3) !left
           indx=i
        case default
           stop
        end select
        cell%data_pers%QP(j,  i,:) = QP(indx,:)
        cell%data_pers%FP(j,1,i,:) = FP(1,indx,:)               
        cell%data_pers%FP(j,2,i,:) = FP(2,indx,:)
     end do
  end do
end subroutine initialise_riemann_arguments


! subroutine compute_volume_update_o(Q_DG_UPDATE, s_ref, f_ref, cell_type)
!   real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3)                               :: Q_DG_UPDATE
!   real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,  2),intent(in)  :: s_ref
!   real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,2,3),intent(in)  :: f_ref
!   integer,intent(in)                                                         :: cell_type

!   !local variables
!   real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,3)   :: source_st, source_ref_st
!   real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,3)   :: volume_flux
!   real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,2,3) :: f
!   integer :: i,j

!   volume_flux = 0
!   source_st = 0
!   f = 0
  
!   f(:,:,1,:) = jacobian_inv(1,1,cell_type) * f_ref(:,:,1,:) +&
!                jacobian_inv(1,2,cell_type) * f_ref(:,:,2,:)
!   f(:,:,2,:) = jacobian_inv(2,1,cell_type) * f_ref(:,:,1,:) +&
!                jacobian_inv(2,2,cell_type) * f_ref(:,:,2,:)
  
!   do i=1,_SWE_DG_ORDER+1
!      volume_flux(:,i,:) = matmul(s_m_inv, &
!           matmul(s_k_x_s_b_3_s_b_2, f(:,i,1,:)) + &
!           matmul(s_k_y_s_b_1_s_b_2, f(:,i,2,:)))
!   end do
  
!   do i=1,_SWE_DG_ORDER+1
!      source_ref_st(:,i,2) = matmul(s_m_inv, -matmul(s_m, s_ref(:,i,1)))
!      source_ref_st(:,i,3) = matmul(s_m_inv, -matmul(s_m, s_ref(:,i,2)))
!   end do
  
!   source_st(:,:,2) = source_ref_st(:,:,2) * jacobian(1,1,cell_type) +&
!                      source_ref_st(:,:,3) * jacobian(1,2,cell_type)
!   source_st(:,:,3) = source_ref_st(:,:,2) * jacobian(2,1,cell_type) +&
!                      source_ref_st(:,:,3) * jacobian(2,2,cell_type)
  
! #if defined(_DEBUG)
!   do j=1,_SWE_DG_DOFS
!      do i=1,_SWE_DG_ORDER+1
!         if(isnan(volume_flux(j,i,1)).or.isnan(source_st(j,i,1))) then
!            print*,epsilon
!            print*,iteration
!            print*,cell%data_pers%troubled
!            print*,"vol"
!            print*,volume_flux
!            print*,"src"
!            print*,source_st
!            exit
!         end if
!      end do
!   end do
! #endif
  
!   volume_flux = volume_flux + source_st
!   do j = 1,_SWE_DG_DOFS
!      Q_DG_UPDATE(j,:) = reshape(matmul(t_a,volume_flux(j,:,:)),(/ 3 /))
!   end do

! #if defined(_DEBUG)
!          do i=1,_SWE_DG_DOFS
!             if(isnan(Q_DG_UPDATE(i,1))) then
!                print*,"nan Q_DG_update"
!                print*,Q_DG_UPDATE
!                print*,"q_i_st"
!                print*,q_i_st
!                exit
!             end if
!          end do
! #endif

!   !!------------------------------!!  
!        end subroutine compute_volume_update_o
  
END MODULE SWE_DG_predictor

#endif
