module gen_matrices
use elementary
use gl_integration
use LU
implicit none


interface matrixToFile
   module procedure matrixToFile_real
   module procedure vectorToFile_int
   module procedure vectorToFile_real
end interface


contains

  subroutine settings(N,basis)
    integer :: N
    character(len=*):: basis
    
    call initGLIntegration(N)
    call set_basis(basis,N)
  end subroutine settings

  subroutine generate_matrices(N)
    integer :: N
    call generate_limiter_matrices(N)
    if (nodal) then
       call generate_nodal_matrices(N)
    else
       print*, "Modal Basis not yet implemented"
    end if
  end subroutine

  !------generate Common (Modal and Nodal) matrices ------!
  !---limiter matrices---!  
  subroutine generate_limiter_matrices(N)
    integer,intent(in)            :: N
    real(kind=GRID_SR),allocatable:: phi(:,:)
    real(kind=GRID_SR),allocatable:: mue_lu(:,:)
    integer           ,allocatable:: mue_lu_pivot(:)

    allocate( phi((2*N+1)*(2*N+1),(N+2)*(N+1)/2))
    allocate( mue_lu((N+2)*(N+1)/2+1,(N+2)*(N+1)/2+1))
    allocate( mue_lu_pivot((N+2)*(N+1)/2+1))
    
    phi=0
    mue_lu=0
    mue_lu_pivot=0
    
    call generatePhi(phi,N)
    call generateMueLU(phi,N,mue_lu,mue_lu_pivot)
    
    !Transformation matrices
    call matrixToFile(N,"phi",phi,"_SWE_PATCH_ORDER*_SWE_PATCH_ORDER,_SWE_DG_DOFS")
    call matrixToFile(N,"mue_lu",mue_lu,"_SWE_DG_DOFS+1,_SWE_DG_DOFS+1")
    call matrixToFile(N,"mue_lu_pivot",mue_lu_pivot,"_SWE_DG_DOFS+1")
  end subroutine generate_limiter_matrices
  
  !------generate Nodal matrices ------!
  subroutine  generate_nodal_matrices(N)

    integer :: N
    
    !---predictor matrices---!
    real(kind=GRID_SR),allocatable:: st_m  (:,:)

    real(kind=GRID_SR),allocatable:: st_m_lu  (:,:)
    integer           ,allocatable:: st_m_lu_pivot  (:)

    real(kind=GRID_SR),allocatable:: st_k_x(:,:)
    real(kind=GRID_SR),allocatable:: st_k_y(:,:)
    real(kind=GRID_SR),allocatable:: st_k_t(:,:)

    real(kind=GRID_SR),allocatable:: st_w_k_t         (:,:)
    real(kind=GRID_SR),allocatable:: st_w_k_t_1_0     (:,:)

    real(kind=GRID_SR),allocatable:: st_w_k_t_1_1_lu  (:,:)
    integer           ,allocatable:: st_w_k_t_1_1_lu_pivot(:)

  
    !---solver matrices---!

    real(kind=GRID_SR),allocatable:: tau(:,:)

    integer           ,allocatable:: s_m_lu_pivot(:)
    real(kind=GRID_SR),allocatable:: ss_m(:,:)
    real(kind=GRID_SR),allocatable:: ss_k_x(:,:)
    real(kind=GRID_SR),allocatable:: ss_k_y(:,:)

    real(kind=GRID_SR),allocatable:: s_m(:,:)
    real(kind=GRID_SR),allocatable:: s_m_lu(:,:)
    real(kind=GRID_SR),allocatable:: s_k_x(:,:)
    real(kind=GRID_SR),allocatable:: s_k_y(:,:)

    !---boundary matrices---!
    real(kind=GRID_SR),allocatable:: b_m_1(:,:)
    real(kind=GRID_SR),allocatable:: b_m_2(:,:)
    real(kind=GRID_SR),allocatable:: b_m_3(:,:)

    !---local variables---!
    integer :: code,info,M
    
    !---boundary matrices---!
    allocate( b_m_1((N+2)*(N+1)/2,(N+2)*(N+1)/2*(N+1)))
    allocate( b_m_2((N+2)*(N+1)/2,(N+2)*(N+1)/2*(N+1)))
    allocate( b_m_3((N+2)*(N+1)/2,(N+2)*(N+1)/2*(N+1)))
    
    b_m_1 = 0
    b_m_2 = 0
    b_m_3 = 0

    call generate_b_m_1(b_m_1,N)
    call generate_b_m_2(b_m_2,N)
    call generate_b_m_3(b_m_3,N)

    call matrixToFile(N,"b_m_1",b_m_1,"_SWE_DG_DOFS,_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
    call matrixToFile(N,"b_m_2",b_m_2,"_SWE_DG_DOFS,_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
    call matrixToFile(N,"b_m_3",b_m_3,"_SWE_DG_DOFS,_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
    !-----------------------!

    !---solver only space matrices---!
    allocate( s_m((N+2)*(N+1)/2,(N+2)*(N+1)/2))
    allocate( s_m_lu((N+2)*(N+1)/2,(N+2)*(N+1)/2))
    allocate( s_m_lu_pivot((N+2)*(N+1)/2))
    allocate( s_k_x((N+2)*(N+1)/2,(N+2)*(N+1)/2))
    allocate( s_k_y((N+2)*(N+1)/2,(N+2)*(N+1)/2))

    s_m=0
    s_m_lu=0
    s_m_lu_pivot=0
    s_k_x=0
    s_k_y=0

    call generate_s_k_x(s_k_x,N)
    call generate_s_k_y(s_k_y,N)

    !---solver matrices: space time test functions---!
    
    allocate( ss_m((N+2)*(N+1)/2,(N+2)*(N+1)/2*(N+1)))
    allocate( ss_k_x((N+2)*(N+1)/2,(N+2)*(N+1)/2*(N+1)))
    allocate( ss_k_y((N+2)*(N+1)/2,(N+2)*(N+1)/2*(N+1)))
    
    call generate_s_m(s_m,N)
    s_m_lu = s_m
    call ludcmp(s_m_lu,(N+2)*(N+1)/2, s_m_lu_pivot, code,info)

    ss_m=0
    ss_k_x=0
    ss_k_y=0
    
    call generate_ss_k(ss_k_x,s_k_x,N)
    call generate_ss_k(ss_k_y,s_k_y,N)
    call generate_ss_k(ss_m,s_m,N)

    call matrixToFile(N,"s_m"  ,ss_m  ,"_SWE_DG_DOFS,_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
    call matrixToFile(N,"s_k_x",ss_k_x,"_SWE_DG_DOFS,_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
    call matrixToFile(N,"s_k_y",ss_k_y,"_SWE_DG_DOFS,_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
    call matrixToFile(N,"s_m_lu",s_m_lu,"_SWE_DG_DOFS,_SWE_DG_DOFS")
    call matrixToFile(N,"s_m_lu_pivot",s_m_lu_pivot,"_SWE_DG_DOFS")

    allocate( tau(N+1,N+1))
    allocate( st_m  ((N+2)*(N+1)/2*(N+1),(N+2)*(N+1)/2*(N+1)))
    allocate( st_m_lu  ((N+2)*(N+1)/2*(N+1),(N+2)*(N+1)/2*(N+1)))
    allocate( st_m_lu_pivot  ((N+2)*(N+1)/2*(N+1)))

    allocate( st_k_x((N+2)*(N+1)/2*(N+1),(N+2)*(N+1)/2*(N+1)))
    allocate( st_k_y((N+2)*(N+1)/2*(N+1),(N+2)*(N+1)/2*(N+1)))
    allocate( st_k_t((N+2)*(N+1)/2*(N+1),(N+2)*(N+1)/2*(N+1)))

    tau=0
    st_m=0
    st_k_x=0
    st_k_y=0
    st_k_t=0


    call generate_tau(tau,N)
    call generate_st_m(st_m,s_m,tau,N)
    call generate_st_k_t(st_k_t,s_m,N)
    call generate_st_k(st_k_x,s_k_x,tau,N)
    call generate_st_k(st_k_y,s_k_y,tau,N)

    M=(N+2)*(N+1)/2
    
    call matrixToFile(N,"st_m",st_m(M+1:M*(N+1),:)    ,"_SWE_DG_DOFS*(_SWE_DG_ORDER),_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
    call matrixToFile(N,"st_k_x",st_k_x(M+1:M*(N+1),:),"_SWE_DG_DOFS*(_SWE_DG_ORDER),_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
    call matrixToFile(N,"st_k_y",st_k_y(M+1:M*(N+1),:),"_SWE_DG_DOFS*(_SWE_DG_ORDER),_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
    
    st_m_lu = st_m
    st_m_lu_pivot=0
    
    call ludcmp(st_m_lu,(N+1)*(M), st_m_lu_pivot, code,info)
    call matrixToFile(N,"st_m_lu",st_m_lu,"_SWE_DG_DOFS*(_SWE_DG_ORDER+1),_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
    call matrixToFile(N,"st_m_lu_pivot",st_m_lu_pivot,"_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")

    allocate(st_w_k_t( (N+2)*(N+1)/2*(N+1) , (N+2)*(N+1)/2*(N+1) ))
    st_w_k_t=st_k_t
        
    allocate(st_w_k_t_1_0((N+2)*(N+1)/2*N,(N+2)*(N+1)/2))
    st_w_k_t_1_0   =st_w_k_t( M+1 : M*(N+1),1:M)
    call matrixToFile(N,"st_w_k_t_1_0",st_w_k_t_1_0,"_SWE_DG_DOFS*(_SWE_DG_ORDER),(_SWE_DG_DOFS)")
    
    allocate( st_w_k_t_1_1_lu((N+2)*(N+1)/2*N,(N+2)*(N+1)/2*N))
    allocate( st_w_k_t_1_1_lu_pivot((N+2)*(N+1)/2*N))

    st_w_k_t_1_1_lu=st_w_k_t(M+1:M*(N+1),M+1:M*(N+1))
!    st_w_k_t_1_1_lu=0
    call ludcmp(st_w_k_t_1_1_lu,(N+2)*(N+1)/2*N, st_w_k_t_1_1_lu_pivot, code,info)
    call matrixToFile(N,"st_w_k_t_1_1_lu",st_w_k_t_1_1_lu,"_SWE_DG_DOFS*(_SWE_DG_ORDER),_SWE_DG_DOFS*(_SWE_DG_ORDER)")
    call matrixToFile(N,"st_w_k_t_1_1_lu_pivot",st_w_k_t_1_1_lu_pivot ,"_SWE_DG_DOFS*(_SWE_DG_ORDER)")


    

    call generate_nodal_derivation_matrices(N)

  end subroutine generate_nodal_matrices


!------generate Modal matrices ------!
  subroutine  generate_modal_matrices(N)
    integer,intent(in)            :: N
    !     !---predictor matrices---!
    !     real(kind=GRID_SR),allocatable:: st_gl_node_vals(:,:)
    !     real(kind=GRID_SR),allocatable:: st_gl_weights(:,:)

    !     !---boundary matrices---!
    !     real(kind=GRID_SR),allocatable:: bnd_gl_node_vals(:,:)
    !     real(kind=GRID_SR),allocatable:: bnd_gl_weights(:,:)

    !     real(kind=GRID_SR),allocatable:: temp(:,:)

    ! real(kind=GRID_SR),allocatable:: s_der_x_gl_node_vals(:,:)
    ! real(kind=GRID_SR),allocatable:: s_der_y_gl_node_vals(:,:)
    ! real(kind=GRID_SR),allocatable:: st_der_x_gl_node_vals(:,:)
    ! real(kind=GRID_SR),allocatable:: st_der_y_gl_node_vals(:,:)


    ! allocate( st_gl_node_vals(N_nodes**3,((N+2)*(N+1))/2 * (N+1)))
    ! allocate( st_gl_weights(((N+2)*(N+1))/2 * (N+1),N_nodes**3))
    
    ! allocate( s_der_x_gl_node_vals(N_nodes**2,((N+2)*(N+1))/2 ))
    ! allocate( s_der_y_gl_node_vals(N_nodes**2,((N+2)*(N+1))/2 ))
    
    ! allocate( st_der_x_gl_node_vals(N_nodes**3,((N+2)*(N+1))/2 *(N+1)))
    ! allocate( st_der_y_gl_node_vals(N_nodes**3,((N+2)*(N+1))/2 *(N+1)))
    
    ! allocate(bnd_gl_node_vals(N_nodes**2,(N+1)**2))
    ! allocate(bnd_gl_weights((N+1)**2,N_nodes**2))
    
    ! allocate( temp(((N+2)*(N+1))/2 ,M * (N+1)))

    ! write(n_nodes_string,"(I0,A)"),N_nodes**3,",_SWE_DG_DOFS*(_SWE_DG_ORDER+1)"
    ! call generate_l2_projection(st_gl_node_vals,st_gl_weights,N)
    
    ! call matrixToFile(N,"st_gl_node_vals",st_gl_node_vals,trim(n_nodes_string))
    
    ! write(n_nodes_string,"(A,I0)"),"_SWE_DG_DOFS*(_SWE_DG_ORDER+1),",N_nodes**3
    ! call matrixToFile(N,"st_gl_weights",st_gl_weights,trim(n_nodes_string))
    
    ! s_der_x_gl_node_vals=0
    ! s_der_y_gl_node_vals=0
    ! st_der_x_gl_node_vals=0
    ! st_der_y_gl_node_vals=0
    
    ! call generate_l2_projection_der(s_der_x_gl_node_vals,s_der_y_gl_node_vals,st_der_x_gl_node_vals,st_der_y_gl_node_vals,N)
    
    ! write(n_nodes_string,"(I0,A)"),N_nodes**3,",_SWE_DG_DOFS*(_SWE_DG_ORDER+1)"
    ! call matrixToFile(N,"st_der_x_gl_node_vals",st_der_x_gl_node_vals,trim(n_nodes_string))
    ! call matrixToFile(N,"st_der_y_gl_node_vals",st_der_y_gl_node_vals,trim(n_nodes_string))
    
    ! write(n_nodes_string,"(I0,A)"),N_nodes**2,",_SWE_DG_DOFS"
    ! call matrixToFile(N,"s_der_x_gl_node_vals",s_der_x_gl_node_vals,trim(n_nodes_string))
    ! call matrixToFile(N,"s_der_y_gl_node_vals",s_der_y_gl_node_vals,trim(n_nodes_string))
    
    ! call generate_l2_projection_bnd(bnd_gl_node_vals,bnd_gl_weights,N)
    
    ! write(n_nodes_string,"(I0,A)"),N_nodes**2,",(_SWE_DG_ORDER+1)*(_SWE_DG_ORDER+1)"
    ! call matrixToFile(N,"bnd_gl_node_vals",bnd_gl_node_vals,trim(n_nodes_string))
    ! write(n_nodes_string,"(A,I0)"),"(_SWE_DG_ORDER+1)*(_SWE_DG_ORDER+1),",N_nodes**2
    ! call matrixToFile(N,"bnd_gl_weights",bnd_gl_weights,trim(n_nodes_string))
    
    
  end subroutine

  subroutine generate_nodal_derivation_matrices(N)
    integer,intent(in)            :: N
    !---Matrices used for derivations---!
    real(kind=GRID_SR),allocatable:: basis_der_x(:,:)
    real(kind=GRID_SR),allocatable:: basis_der_y(:,:)
    real(kind=GRID_SR),allocatable:: basis_der_st_x(:,:)
    real(kind=GRID_SR),allocatable:: basis_der_st_y(:,:)

    allocate( basis_der_x((N+1)*(N+2)/2,(N+1)*(N+2)/2))
    allocate( basis_der_y((N+1)*(N+2)/2,(N+1)*(N+2)/2))
    allocate( basis_der_st_x((N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2*(N+1)))
    allocate( basis_der_st_y((N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2*(N+1)))

    call generate_basis_st_der_y(basis_der_st_y,N)
    call generate_basis_st_der_x(basis_der_st_x,N)
    call generate_basis_der_y(basis_der_y,N)
    call generate_basis_der_x(basis_der_x,N)
    
    call matrixToFile(N,"basis_der_st_y",basis_der_st_y,"_SWE_DG_DOFS*(_SWE_DG_ORDER+1),_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
    call matrixToFile(N,"basis_der_st_x",basis_der_st_x,"_SWE_DG_DOFS*(_SWE_DG_ORDER+1),_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
    call matrixToFile(N,"basis_der_y",basis_der_y,"_SWE_DG_DOFS,_SWE_DG_DOFS")
    call matrixToFile(N,"basis_der_x",basis_der_x,"_SWE_DG_DOFS,_SWE_DG_DOFS")
  end subroutine generate_nodal_derivation_matrices

  !------Single matrices-----!

  !--- space mass matrix ---!
  subroutine generate_s_m(s_m,N)
    integer :: i_1,j_1,i_2, j_2, l_1,l_2,N
    real(kind=GRID_SR), intent (inout) :: s_m((N+1)*(N+2)/2,(N+1)*(N+2)/2)
    do j_1=0,N
       do i_1=0,N-j_1
          l_1 = i_1+1 + (N+1)*(N+2)/2 - (N+1-j_1)*(N+2-j_1)/2 
          do j_2=0,N
             do i_2=0,N-j_2
                l_2 = i_2+1 + (N+1)*(N+2)/2 - (N+1-j_2)*(N+2-j_2)/2 
                s_m (l_1,l_2)= gl_quadr_2d(basis_polynomial_array(get_gl_nodes_2d((/ 0.0_GRID_SR,0.0_GRID_SR /),(/ 1.0_GRID_SR,1.0_GRID_SR /)),N,i_1,j_1) * basis_polynomial_array(get_gl_nodes_2d((/ 0.0_GRID_SR,0.0_GRID_SR /),(/ 1.0_GRID_SR,1.0_GRID_SR /)),N,i_2,j_2),0.0_GRID_SR,1.0_GRID_SR)
             end do
          end do
       end do
    end do

  end subroutine

   !--- space stiffness matrix x dimension ---!
   subroutine generate_s_k_x(s_k_x,N)
    integer :: i_1,j_1,i_2, j_2, l_1,l_2,N
    real(kind=GRID_SR), intent (inout) :: s_k_x((N+1)*(N+2)/2,(N+1)*(N+2)/2)

    do j_1=0,N
       do i_1=0,N-j_1
          l_1 = i_1+1 + (N+1)*(N+2)/2 - (N+1-j_1)*(N+2-j_1)/2 
          do j_2=0,N
             do i_2=0,N-j_2
                l_2 = i_2+1 + (N+1)*(N+2)/2 - (N+1-j_2)*(N+2-j_2)/2 
                s_k_x (l_1,l_2)= gl_quadr_2d( &
basis_polynomial_der_x_array(get_gl_nodes_2d((/ 0.0_GRID_SR,0.0_GRID_SR /),(/ 1.0_GRID_SR,1.0_GRID_SR /)),N,i_1,j_1) * &
basis_polynomial_array(get_gl_nodes_2d((/ 0.0_GRID_SR,0.0_GRID_SR /),(/ 1.0_GRID_SR,1.0_GRID_SR /)),N,i_2,j_2)&
,0.0_GRID_SR,1.0_GRID_SR)
             end do
          end do
       end do
    end do
  
   end subroutine generate_s_k_x

   !--- space stiffness matrix y dimension ---!
   subroutine generate_s_k_y(s_k_y,N)
    integer :: i_1,j_1,i_2, j_2, l_1,l_2,N
    real(kind=GRID_SR), intent (inout) :: s_k_y((N+1)*(N+2)/2,(N+1)*(N+2)/2)

    do j_1=0,N
       do i_1=0,N-j_1
          l_1 = i_1+1 + (N+1)*(N+2)/2 - (N+1-j_1)*(N+2-j_1)/2 
          do j_2=0,N
             do i_2=0,N-j_2
                l_2 = i_2+1 + (N+1)*(N+2)/2 - (N+1-j_2)*(N+2-j_2)/2 
                s_k_y (l_1,l_2)= gl_quadr_2d( &
 basis_polynomial_der_y_array(get_gl_nodes_2d((/ 0.0_GRID_SR,0.0_GRID_SR /),(/ 1.0_GRID_SR,1.0_GRID_SR /)) ,N,i_1,j_1) * &
basis_polynomial_array(get_gl_nodes_2d((/ 0.0_GRID_SR,0.0_GRID_SR /),(/ 1.0_GRID_SR,1.0_GRID_SR /)),N,i_2,j_2)&
,0.0_GRID_SR,1.0_GRID_SR)
             end do
          end do
       end do
    end do
  
 
   end subroutine generate_s_k_y

   !---predictor stiffness matrices with only spacial test functions ---!
   subroutine generate_ss_k(ss_k,s_k,N)
    integer :: i,N
    real(kind=GRID_SR), intent (in) :: s_k((N+1)*(N+2)/2,(N+1)*(N+2)/2)
    real(kind=GRID_SR) :: t(1,N+1)
    real(kind=GRID_SR), intent (inout) :: ss_k((N+1)*(N+2)/2,(N+1)*(N+2)/2*(N+1))

    do i=0,N
       t(1,i+1)=gl_quadr_1d(basis_polynomial_1d_array(gl_nodes,N,i),0.0_GRID_SR,1.0_GRID_SR)
    end do
    
    ss_k = kroneckerMultiplication(t,s_k)
   end subroutine generate_ss_k

   !---space time stiffness matrix in time ---!
   subroutine generate_st_k_t(st_k_t,s_m,N)
     integer :: N,i,j
     real(kind=GRID_SR), intent (inout) :: st_k_t((N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2*(N+1))     
     real(kind=GRID_SR), intent (in) :: s_m((N+1)*(N+2)/2,(N+1)*(N+2)/2)     
     real(kind=GRID_SR)              :: tau_dt(N+1,N+1)
     
     do i=0,N
        do j=0,N
           tau_dt(i+1,j+1)=gl_quadr_1d(basis_polynomial_1d_array(gl_nodes,N,i)*basis_polynomial_1d_der_array(gl_nodes,N,j),0.0_GRID_SR,1.0_GRID_SR)
        end do
     end do
     
     st_k_t=kroneckerMultiplication(tau_dt,s_m)
   end subroutine

   !---space time mass matrix along time boundary ---!
   !--- TODO: here only for nodes on time boundray !!! ---!
   subroutine generate_st_w(w,s_m,N)
     real(kind=GRID_SR), intent (inout) :: w((N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2*(N+1))     
     real(kind=GRID_SR), intent (in) :: s_m((N+1)*(N+2)/2,(N+1)*(N+2)/2)     
     integer,intent(in) :: N

     integer ::i

     w(1:(N+1)*(N+2)/2,1:(N+1)*(N+2)/2)=-s_m
     w((N+1)*(N+2)/2*N+1:(N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2*N+1:(N+1)*(N+2)/2*(N+1))=s_m

   end subroutine generate_st_w

   subroutine generate_st_w_lu(w_lu,w_lu_pivot,s_m,N)
     real(kind=GRID_SR), intent (inout) :: w_lu((N+1)*(N+2)/2*(N),(N+1)*(N+2)/2*(N))
     integer, intent(inout)            :: w_lu_pivot((N+1)*(N+2)/2*(N))     

     real(kind=GRID_SR), intent (in) :: s_m((N+1)*(N+2)/2,(N+1)*(N+2)/2)     

     real(kind=GRID_SR)  ::w_sm((N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2*(N+1))
     integer,intent(in)  :: N
     integer :: code,info

     call generate_st_w(w_sm,s_m,N)
     
     w_lu = w_sm((N+1)*(N+2)/2+1:(N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2+1:(N+1)*(N+2)/2*(N+1))

     call ludcmp(w_lu, (N+1)*(N+2)/2*N ,w_lu_pivot ,code ,info)

   end subroutine generate_st_w_lu

   !----- 1d Polynomials integrated in time, used for space time matrices -----!
   subroutine generate_tau(tau,N)
     real(kind=GRID_SR),intent(inout) :: tau(N+1,N+1)
     integer ,intent(in) :: N
     integer  :: i,j
     
     do i=0,N
        do j=0,N
           tau(i+1,j+1)=gl_quadr_1d(basis_polynomial_1d_array(gl_nodes,N,i)*basis_polynomial_1d_array(gl_nodes,N,j),0.0_GRID_SR,1.0_GRID_SR)
        end do
     end do
     
   end subroutine generate_tau

   !----- space time stiffness matrices  -----!
    subroutine generate_st_k(st_k,s_k,tau,N)
      integer,intent(in) :: N
      real(kind=GRID_SR),intent(in) :: tau(N+1,N+1),s_k((N+1)*(N+2)/2,(N+1)*(N+2)/2)
      real(kind=GRID_SR),intent(inout) :: st_k((N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2*(N+1))

      st_k=kroneckerMultiplication(tau,transpose(s_k))
            
    end subroutine generate_st_k

   !----- space time mass matrix  -----!    
   subroutine generate_st_m(st_m,s_m ,tau,N)      
     integer,intent(in) :: N
     real(kind=GRID_SR),intent(in) :: tau(N+1,N+1),s_m((N+1)*(N+2)/2,(N+1)*(N+2)/2)
     real(kind=GRID_SR),intent(inout) :: st_m((N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2*(N+1))
     st_m=kroneckerMultiplication(tau,s_m)
   end subroutine generate_st_m
   
   !-----deriving matrices-----!   
   subroutine generate_basis_der_x(basis_der_x,N)
     integer :: i_1,j_1,i_2,j_2,l_1,l_2,N

     real(kind=GRID_SR) :: basis_der_x((N+1)*(N+2)/2,(N+1)*(N+2)/2)
     real(kind=GRID_SR) ::x,y

     basis_der_x=0
     do j_1=0,N
        do i_1=0,N-j_1
           l_1 = i_1+1 + (N+1)*(N+2)/2 - (N+1-j_1)*(N+2-j_1)/2 
           x=Real(i_1,kind=GRID_SR)/Real(N,kind=GRID_SR)
           y=Real(j_1,kind=GRID_SR)/Real(N,kind=GRID_SR)
           do j_2=0,N
              do i_2=0,N-j_2
                 l_2 = i_2+1 + (N+1)*(N+2)/2 - (N+1-j_2)*(N+2-j_2)/2 
                 basis_der_x(l_1,l_2)= basis_polynomial_der_x(x,y,N,i_2,j_2)
              end do
           end do
        end do
     end do

   end subroutine generate_basis_der_x

   subroutine generate_basis_der_y(basis_der_y,N)
     integer :: i_1,j_1,i_2,j_2,l_1,l_2,N
     real(kind=GRID_SR) :: basis_der_y((N+1)*(N+2)/2,(N+1)*(N+2)/2)
     Real(kind=GRID_SR) ::x,y
     basis_der_y=0
     do j_1=0,N
        do i_1=0,N-j_1
           l_1 = i_1+1 + (N+1)*(N+2)/2 - (N+1-j_1)*(N+2-j_1)/2 
           x=Real(i_1,kind=GRID_SR)/Real(N,kind=GRID_SR)
           y=Real(j_1,kind=GRID_SR)/Real(N,kind=GRID_SR)
           do j_2=0,N
              do i_2=0,N-j_2
                 l_2 = i_2+1 + (N+1)*(N+2)/2 - (N+1-j_2)*(N+2-j_2)/2 
                 basis_der_y(l_1,l_2)= basis_polynomial_der_y(x,y,N,i_2,j_2)
              end do
           end do
        end do
     end do

   end subroutine generate_basis_der_y

   subroutine generate_basis_st_der_x(basis_st_der_x,N)
     integer :: i_1,j_1,i_2,j_2,l_1,l_2,k_1,k_2,N

     real(kind=GRID_SR) :: basis_st_der_x((N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2*(N+1))
     real(kind=GRID_SR) ::x,y,t

     basis_st_der_x=0
     do k_1=0,N
        t=Real(k_1,kind=GRID_SR)/Real(N,kind=GRID_SR)
        do k_2=0,N
           do j_1=0,N
              do i_1=0,N-j_1
                 l_1 = i_1+1 + (N+1)*(N+2)/2 - (N+1-j_1)*(N+2-j_1)/2 + (N+1)*(N+2)/2 * k_1
                 x=Real(i_1,kind=GRID_SR)/Real(N,kind=GRID_SR)
                 y=Real(j_1,kind=GRID_SR)/Real(N,kind=GRID_SR)
                 do j_2=0,N
                    do i_2=0,N-j_2
                       l_2 = i_2+1 + (N+1)*(N+2)/2 - (N+1-j_2)*(N+2-j_2)/2  + (N+1)*(N+2)/2 * k_2
                       basis_st_der_x(l_1,l_2)= basis_polynomial_der_x(x,y,N,i_2,j_2)*basis_polynomial_1d(t,N,k_2)
                    end do
                 end do
              end do
           end do
        end do
     end do

   end subroutine generate_basis_st_der_x

   subroutine generate_basis_st_der_y(basis_st_der_y,N)
     integer :: i_1,j_1,i_2,j_2,l_1,l_2,k_1,k_2,N

     real(kind=GRID_SR) :: basis_st_der_y((N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2*(N+1))
     real(kind=GRID_SR) ::x,y,t

     basis_st_der_y=0
     do k_1=0,N
        t=Real(k_1,kind=GRID_SR)/Real(N,kind=GRID_SR)
        do k_2=0,N
           do j_1=0,N
              do i_1=0,N-j_1
                 l_1 = i_1+1 + (N+1)*(N+2)/2 - (N+1-j_1)*(N+2-j_1)/2 + (N+1)*(N+2)/2 * k_1
                 x=Real(i_1,kind=GRID_SR)/Real(N,kind=GRID_SR)
                 y=Real(j_1,kind=GRID_SR)/Real(N,kind=GRID_SR)
                 do j_2=0,N
                    do i_2=0,N-j_2
                       l_2 = i_2+1 + (N+1)*(N+2)/2 - (N+1-j_2)*(N+2-j_2)/2  + (N+1)*(N+2)/2 * k_2
                       basis_st_der_y(l_1,l_2)= basis_polynomial_der_y(x,y,N,i_2,j_2)*basis_polynomial_1d(t,N,k_2)
                    end do
                 end do
              end do
           end do
        end do
     end do

   end subroutine generate_basis_st_der_y

   !-----L2 Projection matrices -----!
   subroutine generate_l2_projection(basis_gl_node_vals,basis_gl_weights,N)
     integer,intent(in) :: N

     integer            :: i,j,k,l,N_gl_nodes

     !  real(kind=GRID_SR) :: basis_gl_node_vals(N_nodes**3,((N+2)*(N+1))/2 * (N+1))
     !  real(kind=GRID_SR) :: basis_gl_weights(((N+2)*(N+1))/2 * (N+1),N_nodes**3)

     real(kind=GRID_SR) :: basis_gl_node_vals(:,:)
     real(kind=GRID_SR) :: basis_gl_weights(:,:)

     !  real(kind=GRID_SR) :: vals_space(N_nodes**2,(N+2)*(N+1)/2)
     !  real(kind=GRID_SR) :: vals_time(N_nodes,N+1)

     real(kind=GRID_SR),allocatable :: vals_space(:,:)
     real(kind=GRID_SR),allocatable :: vals_time(:,:)
     real(kind=GRID_SR),allocatable :: gl_nodes_2d(:,:)

     !  real(kind=GRID_SR) :: gl_nodes_2d(2,N_nodes**2)

     N_gl_nodes=get_num_gl_nodes()

     allocate(vals_space(N_gl_nodes**2,(N+2)*(N+1)/2))
     allocate(vals_time(N_gl_nodes,N+1))
     allocate(gl_nodes_2d(2,N_gl_nodes**2))


     basis_gl_node_vals=0
     vals_space=0
     vals_time=0

     gl_nodes_2d = get_gl_nodes_2d((/0.0_GRID_SR,0.0_GRID_SR/),(/1.0_GRID_SR,1.0_GRID_SR/))

     do j=0,N
        do i=0,N-j
           l = i+1 + (N+1)*(N+2)/2 - (N+1-j)*(N+2-j)/2 
           vals_space(:,l) = basis_polynomial_array(gl_nodes_2d,N,i,j)
        end do
     end do

     do j=0,N
        vals_time(:,1+j)=basis_polynomial_1d_array(gl_nodes,N,j)
     end do

     basis_gl_node_vals= kroneckerMultiplication(vals_time,vals_space)

     do j=0,N_nodes-1
        do i=0,N_nodes-1
           vals_space((i+1)+N_nodes*j,:)= vals_space((i+1)+N_nodes*j,:)*gl_weights(i+1)*gl_weights(j+1)*gl_nodes(N_nodes-j)
        end do
     end do

     do i=1,N_nodes
        vals_time(i,:)=vals_time(i,:)*gl_weights(i)
     end do

     basis_gl_weights=Transpose(kroneckerMultiplication(vals_time,vals_space))

   end subroutine generate_l2_projection

   subroutine generate_l2_projection_der(s_der_x_gl_node_vals,s_der_y_gl_node_vals,st_der_x_gl_node_vals,st_der_y_gl_node_vals,N)
     real(kind=GRID_SR),intent(inout) :: s_der_y_gl_node_vals(:,:)
     real(kind=GRID_SR),intent(inout) :: s_der_x_gl_node_vals(:,:)
     real(kind=GRID_SR),intent(inout) :: st_der_y_gl_node_vals(:,:)
     real(kind=GRID_SR),intent(inout) :: st_der_x_gl_node_vals(:,:)
     real(kind=GRID_SR),allocatable :: gl_nodes_2d(:,:)
     real(kind=GRID_SR),allocatable :: vals_space_x(:,:)
     real(kind=GRID_SR),allocatable :: vals_space_y(:,:)
     real(kind=GRID_SR),allocatable :: vals_time(:,:)
     integer :: i,j,l,N_gl_nodes
     integer,intent(in) :: N

     N_gl_nodes=get_num_gl_nodes()

     allocate(vals_space_x(N_gl_nodes**2,(N+2)*(N+1)/2))
     allocate(vals_space_y(N_gl_nodes**2,(N+2)*(N+1)/2))
     allocate(vals_time(N_gl_nodes,N+1))

     allocate(gl_nodes_2d(2,N_gl_nodes**2))


     gl_nodes_2d = get_gl_nodes_2d((/0.0_GRID_SR,0.0_GRID_SR/),(/1.0_GRID_SR,1.0_GRID_SR/))

     s_der_x_gl_node_vals=0
     s_der_y_gl_node_vals=0

     l=0
     do j=0,N
        do i=0,N-j
           l=l+1
           vals_space_x(:,l)= basis_polynomial_der_x_array(gl_nodes_2d,N,i,j)
        end do
     end do

     l=0
     do j=0,N
        do i=0,N-j
           l=l+1
           vals_space_y(:,l)= basis_polynomial_der_y_array(gl_nodes_2d,N,i,j)
        end do
     end do

     do j=0,N
        vals_time(:,1+j)=basis_polynomial_1d_array(gl_nodes,N,j)
     end do

     s_der_x_gl_node_vals= vals_space_x
     s_der_y_gl_node_vals= vals_space_y

     st_der_x_gl_node_vals= kroneckerMultiplication(vals_time,vals_space_x)
     st_der_y_gl_node_vals= kroneckerMultiplication(vals_time,vals_space_y)

   end subroutine generate_l2_projection_der

   subroutine generate_l2_projection_bnd(bnd_gl_node_vals,bnd_gl_weights,N)
     integer,intent(in) :: N
     real(kind=GRID_SR),intent(inout) ::bnd_gl_weights(:,:),bnd_gl_node_vals(:,:)
     real(kind=GRID_SR),allocatable :: vals_space(:,:)
     real(kind=GRID_SR),allocatable :: vals_time(:,:)
     integer :: i,j,l,N_gl_nodes  
     N_gl_nodes=get_num_gl_nodes()  

     allocate(vals_space(N_gl_nodes,N+1))
     allocate(vals_time(N_gl_nodes,N+1))

     do i=0,N
        l=i
        vals_space(:,i+1)=basis_polynomial_1d_array(gl_nodes,N,l)
     end do

     do i=0,N
        l=i
        vals_time(:,i+1)=basis_polynomial_1d_array(gl_nodes,N,l)
     end do

     bnd_gl_node_vals=kroneckerMultiplication(vals_time,vals_space)

     do i=0,N
        vals_space(:,i+1)=vals_space(:,i+1)*gl_weights
        vals_time(:,i+1)=vals_time(:,i+1)*gl_weights
     end do

     bnd_gl_weights=Transpose(kroneckerMultiplication(vals_time,vals_space))

   end subroutine generate_l2_projection_bnd



   !-----Nodal Boundary Matrices-----!
   subroutine generate_b_m_1(b_m_1,N) 
     integer,intent(in) :: N
     real(kind=GRID_SR),intent(inout) :: b_m_1((N+1)*(N+2)/2,(N+1)*(N+2)/2*(N+1))
     real(kind=GRId_SR) :: nodes(2,N_nodes)

     real(kind=GRID_SR) :: time_int
     integer :: i_1,i_2,j_1,j_2,l_1,l_2,k,j,M,M_st

     nodes(1,:)=gl_nodes
     nodes(2,:)=0

     do j_1=0,N
        do i_1=0,N-j_1
           l_1 = i_1+1 + (N+1)*(N+2)/2 - (N+1-j_1)*(N+2-j_1)/2 

           do k=0,N
              time_int= gl_quadr_1d(basis_polynomial_1d_array(gl_nodes,N,k),0.0_GRID_SR,1.0_GRID_SR)
              do j_2=0,N
                 do i_2=0,N-j_2
                    l_2 = i_2+1 + (N+1)*(N+2)/2 - (N+1-j_2)*(N+2-j_2)/2 + (N+1)*(N+2)/2 * k

                    b_m_1(l_1,l_2)=gl_quadr_1d(basis_polynomial_array(nodes,N,i_1,j_1)*&
                         basis_polynomial_array(nodes,N,i_2,j_2),0.0_GRID_SR,1.0_GRID_SR)*&
                         time_int

                 end do
              end do
           end do
        end do
     end do
   end subroutine generate_b_m_1


   subroutine generate_b_m_2(b_m_2,N)
     integer,intent(in) :: N
     real(kind=GRID_SR),intent(inout) :: b_m_2((N+1)*(N+2)/2,(N+1)*(N+2)/2*(N+1))
     real(kind=GRID_SR) :: nodes(2,N_nodes)
     real(kind=GRID_SR) :: time_int
     integer :: i_1,i_2,j_1,j_2,l_1,l_2,k,j,M,M_st

     nodes(1,:)=gl_nodes
     nodes(2,:)=gl_nodes(N_nodes:1:-1)


     do j_1=0,N
        do i_1=0,N-j_1
           l_1 = i_1+1 + (N+1)*(N+2)/2 - (N+1-j_1)*(N+2-j_1)/2 
           do k=0,N
              time_int= gl_quadr_1d(basis_polynomial_1d_array(gl_nodes,N,k),0.0_GRID_SR,1.0_GRID_SR)
              do j_2=0,N
                 do i_2=0,N-j_2
                    l_2 = i_2+1 + (N+1)*(N+2)/2 - (N+1-j_2)*(N+2-j_2)/2 + (N+1)*(N+2)/2 * k

                    b_m_2(l_1,l_2)=gl_quadr_1d(basis_polynomial_array(nodes,N,i_1,j_1)*&
                         basis_polynomial_array(nodes,N,i_2,j_2),0.0_GRID_SR,sqrt(2.0_GRID_SR))*&
                         time_int
                 end do
              end do
           end do
        end do
     end do
   end subroutine generate_b_m_2


   subroutine generate_b_m_3(b_m_3,N)
     integer,intent(in) :: N
     real(kind=GRID_SR),intent(inout) :: b_m_3((N+1)*(N+2)/2,(N+1)*(N+2)/2*(N+1))
     real(kind=GRId_SR) :: nodes(2,N_nodes)
     real(kind=GRID_SR) :: time_int

     integer :: i_1,i_2,j_1,j_2,l_1,l_2,k,j,M,M_st

     nodes(1,:)=0
     nodes(2,:)=gl_nodes

     do j_1=0,N
        do i_1=0,N-j_1
           l_1 = i_1+1 + (N+1)*(N+2)/2 - (N+1-j_1)*(N+2-j_1)/2 
           do k=0,N
              time_int= gl_quadr_1d(basis_polynomial_1d_array(gl_nodes,N,k),0.0_GRID_SR,1.0_GRID_SR)
              do j_2=0,N
                 do i_2=0,N-j_2
                    l_2 = i_2+1 + (N+1)*(N+2)/2 - (N+1-j_2)*(N+2-j_2)/2 + (N+1)*(N+2)/2 * k
                    b_m_3(l_1,l_2)=gl_quadr_1d(basis_polynomial_array(nodes,N,i_1,j_1)*&
                         basis_polynomial_array(nodes,N,i_2,j_2),0.0_GRID_SR,1.0_GRID_SR)*&
                         time_int
                 end do
              end do
           end do
        end do
     end do
   end subroutine generate_b_m_3

   !-----Limiter Matrices-----!

   function generateTriangles(N) result(triangles)
     integer :: N_s,N_s_squared,M,i,j,l,pnt,k,step, count=1 ,N

     real(kind=GRID_SR) :: triangles(2,3*(N*2 + 1)*(N*2 + 1))
     real(kind=GRID_SR),Allocatable :: points(:,:)

     N_s = N*2 + 1
     N_s_squared = N_s*N_s
     M = (N_s+1)*(N_s+2)/2

     step =2
     allocate(points(2,M))
     count = 1
     do j=0,N_s
        do i=0,j
           points(1,count) = real(j-i,kind=GRID_SR)/real(N_s,kind=GRID_SR)
           points(2,count) = real(i,kind=GRID_SR)/real(N_s,kind=GRID_SR)
           count = count +1
        end do
     end do

     !dedicate vertices to triangles
     !lower
     do i=1,N_s
        do k = 1,i
           j = (i-1)**2 + k * 2 -1
           pnt = (i-1)*i/2 + k
           triangles(:,(j-1)*3+1)=points(:,pnt)
           pnt = (i+1)*i/2 + k
           triangles(:,(j-1)*3+2)=points(:,pnt)
           triangles(:,((j-1))*3+3)=points(:,pnt+1)
        end do
     end do

     !upper
     do i=2,N_s
        do k = 1,i-1
           j = (i-1)**2 + k * 2
           pnt = (i-1)*(i)/2 + k
           triangles(1:2,(j-1)*3+1)=points(1:2,pnt)
           triangles(:,(j-1)*3+2)=points(:,pnt+1)
           pnt = (i+1)*(i)/2 +k+1
           triangles(:,(j-1)*3+3)=points(:,pnt)

        end do
     end do

   end function generateTriangles

   subroutine generatePhi(phi,N)
     integer,intent(in)                 :: N
     real(kind=GRID_SR),intent(inout)   :: phi((N*2 + 1)*(N*2 + 1), (N+1)*(N+2)/2)

     !--local Variables--!
     real(kind=GRID_SR)                 :: nodes_2d(2,N_nodes*N_nodes),a(2),b(2)
     real(kind=GRID_SR)                 :: triangles(2,(N*2 + 1)*(N*2 + 1)*3)
     integer                            :: N_s,N_s_squared,M,i,j,k,l
     logical                            :: lower

     triangles = generateTriangles(N)
     N_s = N*2 + 1
     N_s_squared = N_s*N_s
     M = (N+1)*(N+2)/2

     !1. generate subcell vertices
     do i = 1,N_s_squared

        if(triangles(2,(i-1)*3+2)-triangles(2,(i-1)*3+1) >  1/real(2*N_s+1,kind=GRID_SR) ) then
           lower=.false.
           a(1)=triangles(1,(i-1)*3+2)
           a(2)=triangles(2,(i-1)*3+1)
           b(1)=triangles(1,(i-1)*3+3)
           b(2)=triangles(2,(i-1)*3+3)
        else
           lower=.true.
           a(1)=triangles(1,(i-1)*3+1)
           a(2)=triangles(2,(i-1)*3+1)
           b(1)=triangles(1,(i-1)*3+2)
           b(2)=triangles(2,(i-1)*3+3)
        end if

        nodes_2d=get_gl_nodes_2d(a,b,lower)

        do k = 0,N
           do j = 0,N-k
              l=M-(N+1-k)*(N+2-k)/2+j+1
              phi(i,l)=gl_quadr_2d(basis_polynomial_array(nodes_2d,N,j,k),a(1),b(1))
           end do
        end do
     end do

   end subroutine generatePhi

   subroutine generateMueLU(phi,N,mue_lu,ipiv)
     integer,intent(in)                :: N
     real(kind=GRID_SR),intent(in)     :: phi((N*2+1)*(N*2+1),((N+1)*(N+2)/2))

     integer,intent(inout)             :: ipiv(((N+1)*(N+2)/2)+1)
     real(kind=GRID_SR),intent(inout)  ::mue_lu(((N+1)*(N+2)/2)+1,(N+1)*(N+2)/2+1)

     !--local Variables--!
     integer                           :: M,j,k,l,code, info

     M = (N+1)*(N+2)/2
     mue_lu(1:M,1:M) = matmul(transpose(phi),phi) * 2.0_GRID_SR
     mue_lu(M+1,M+1) = 0.0_GRID_SR

     do k = 0,N
        do j = 0,N-k
           l=M-(N+1-k)*(N+2-k)/2+j+1
           mue_lu(M+1,l) = gl_quadr_2d(basis_polynomial_array(get_gl_nodes_2d((/ 0.0_GRID_SR,0.0_GRID_SR /),(/ 1.0_GRID_SR,1.0_GRID_SR /)),N,j,k),0.0_GRID_SR,1.0_GRID_SR)
           mue_lu(l,M+1) = mue_lu(M+1,l)
        end do
     end do

     call ludcmp(mue_lu,M+1,ipiv,code,info)

   end subroutine generateMueLU


   !-----Matrix Output Routines-----!   
   subroutine matrixToFile_real(N,name,matrix,dims)
     real(kind=GRID_SR),intent(in) :: matrix(:,:)
     real(kind=GRID_SR) ::tolerance=1q-34
     character(len=*), intent(in)       :: name,dims
     integer,intent(in)           :: N
     character(len=1024)      :: file_name
     integer :: i,j
     character(len=12)       :: format_string

     if(GRID_SR.eq.kind(1.0q0)) then
        format_string="(A,E38.31,A)"
     else if(GRID_SR.eq.kind(1.0d0)) then
        format_string="(A,E24.17,A)"
     end if

     write(file_name,'(A,A,A,I0,A)') "output/",name,"_",N,".incl"
     file_name=trim(file_name)
     write(*,*), "!Printing ",name," to ",trim(file_name)
     open(11, FILE=trim(file_name), ACTION="write" , STATUS="replace")
     write(11,"(A,A,A,A,A)") "real(kind=GRID_SR), Parameter :: ", name ,"(",dims,")= reshape((/&"

     j=1
     i=1

     if (abs(matrix(1,1)) < tolerance) then
        write(11, format_string), " ",0.0_GRID_SR, "_GRID_SR&"
     else
        write(11, format_string), " ",matrix(1,1), "_GRID_SR&"
     end if

     do i = 2,size(matrix,1)
        if(abs(matrix(i,j)) < tolerance) then
           write(11, format_string), ",", 0.0_GRID_SR, "_GRID_SR&"
        else
           write(11, format_string), ",", matrix(i,j), "_GRID_SR&"
        end if
     end do


     do j=2,size(matrix,2)
        do i = 1,size(matrix,1)

           if(abs(matrix(i,j)) < tolerance) then
              write(11, format_string), ",", 0.0_GRID_SR, "_GRID_SR&"
           else
              write(11, format_string), ",", matrix(i,j), "_GRID_SR&"
           end if

        end do
     end do
     write(11, "(A,A,A)") , "/),(/",dims,"/))"

     close(11)

   end subroutine matrixToFile_real

   subroutine vectorToFile_int(N,name,vector,dims)
     integer,intent(in) :: vector(:)
     character(len=*), intent(in)       :: name,dims
     integer,intent(in)           :: N
     character(len=1024)      :: file_name
     integer :: i

     write(file_name,'(A,A,A,I0,A)') "output/",name,"_",N,".incl"
     file_name=trim(file_name)
     write(*,*), "!Printing ",name," to ",trim(file_name)
     open(11, FILE=trim(file_name), ACTION="write" , STATUS="replace")

     write(11,"(A,A,A,A,A)") "integer, Parameter :: ", name ,"(",dims,")= (/&"



     i=1

     write(11,"(I10)" ,advance='no' ),  vector(i)

     do i = 2,size(vector)
        write(11 , "(A,I10)" ,advance='no'), ",", vector(i)
     end do
     write(11 , "(A)" ,advance='no'), "/)"

     close(11)

   end subroutine vectorToFile_int

   subroutine vectorToFile_real(N,name,vector,dims)
     real(kind=GRID_SR),intent(in) :: vector(:)
     character(len=*), intent(in)       :: name,dims
     integer,intent(in)           :: N
     character(len=1024)      :: file_name
     integer :: i

     write(file_name,'(A,A,A,I0,A)') "output/",name,"_",N,".incl"
     file_name=trim(file_name)
     write(*,*), "!Printing ",name," to ",trim(file_name)
     open(11, FILE=trim(file_name), ACTION="write" , STATUS="replace")

     write(11,"(A,A,A,A,A)") "real(kind=GRID_SR), Parameter :: ", name ,"(",dims,")= (/&"

     i=1

     write(11,"(E38.31,A)" ,advance='no' ),  vector(i),"_GRID_SR"

     do i = 2,size(vector)
        write(11 , "(A,E38.31,A)" ,advance='no'), ",", vector(i),"_GRID_SR"
     end do
     write(11 , "(A)" ,advance='no'), "/)"

     close(11)

   end subroutine vectorToFile_real

end module gen_matrices


