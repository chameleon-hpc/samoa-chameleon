program matrices
 use gen_matrices
integer :: N

character(len=32) :: arg
character(len=1024) :: n_nodes_string


integer :: M, N_s, N_s_s, code, info ,i,j,N_gl_nodes
!real(kind=GRID_SR),allocatable :: identity(:,:)

call get_command_argument(1,arg)
write(*,*), "using order: ",arg

read(arg,'(I1)') N

M=(N+2)*(N+1)/2
N_s=2*N+1
N_s_s=N_s*N_s

call initGLIntegration(N)

N_gl_nodes=get_num_gl_nodes()

! allocate( identity((M)*(N+1),(M)*(N+1)))
! identity=0
! do i=1,(M)*(N+1)
!    identity(i,i) =1
! end do
!Print*,identity

allocate( phi((2*N+1)*(2*N+1),(N+2)*(N+1)/2))
allocate( mue_lu((N+2)*(N+1)/2+1,(N+2)*(N+1)/2+1))
allocate( mue_lu_pivot((N+2)*(N+1)/2+1))

allocate( s_m((N+2)*(N+1)/2,(N+2)*(N+1)/2))

allocate( s_m_lu((N+2)*(N+1)/2,(N+2)*(N+1)/2))
allocate( s_m_lu_pivot((N+2)*(N+1)/2))
allocate( s_k_x((N+2)*(N+1)/2,(N+2)*(N+1)/2))
allocate( s_k_y((N+2)*(N+1)/2,(N+2)*(N+1)/2))
allocate( ss_m((N+2)*(N+1)/2,(N+2)*(N+1)/2*(N+1)))
allocate( ss_k_x((N+2)*(N+1)/2,(N+2)*(N+1)/2*(N+1)))
allocate( ss_k_y((N+2)*(N+1)/2,(N+2)*(N+1)/2*(N+1)))

allocate( tau(N+1,N+1))
allocate( st_m  ((N+2)*(N+1)/2*(N+1),(N+2)*(N+1)/2*(N+1)))
allocate( st_m_lu  ((N+2)*(N+1)/2*(N+1),(N+2)*(N+1)/2*(N+1)))
allocate( st_m_lu_pivot  ((N+2)*(N+1)/2*(N+1)))
allocate( st_w  ((N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2*(N+1)))
allocate( st_k_x((N+2)*(N+1)/2*(N+1),(N+2)*(N+1)/2*(N+1)))
allocate( st_k_y((N+2)*(N+1)/2*(N+1),(N+2)*(N+1)/2*(N+1)))
allocate( st_k_t((N+2)*(N+1)/2*(N+1),(N+2)*(N+1)/2*(N+1)))

allocate( st_w_k_t         ( (N+2)*(N+1)/2*(N+1) , (N+2)*(N+1)/2*(N+1) ))
allocate( st_w_k_t_1_0     ( (N+2)*(N+1)/2*N     , (N+2)*(N+1)/2))
allocate( st_w_k_t_1_1_lu  ( (N+2)*(N+1)/2*N     , (N+2)*(N+1)/2*N))
allocate( st_w_k_t_1_1_lu_pivot( (N+2)*(N+1)/2*N ))

allocate( b_m_1((N+2)*(N+1)/2,(N+2)*(N+1)/2*(N+1)))
allocate( b_m_2((N+2)*(N+1)/2,(N+2)*(N+1)/2*(N+1)))
allocate( b_m_3((N+2)*(N+1)/2,(N+2)*(N+1)/2*(N+1)))

allocate( st_gl_node_vals(N_nodes**3,((N+2)*(N+1))/2 * (N+1)))
allocate( st_gl_weights(((N+2)*(N+1))/2 * (N+1),N_nodes**3))

allocate( s_der_x_gl_node_vals(N_nodes**2,((N+2)*(N+1))/2 ))
allocate( s_der_y_gl_node_vals(N_nodes**2,((N+2)*(N+1))/2 ))

allocate( st_der_x_gl_node_vals(N_nodes**3,((N+2)*(N+1))/2 *(N+1)))
allocate( st_der_y_gl_node_vals(N_nodes**3,((N+2)*(N+1))/2 *(N+1)))

allocate(bnd_gl_node_vals(N_nodes**2,(N+1)**2))
allocate(bnd_gl_weights((N+1)**2,N_nodes**2))

allocate( temp(((N+2)*(N+1))/2 ,M * (N+1)))
allocate( basis_der_x((N+1)*(N+2)/2,(N+1)*(N+2)/2))
allocate( basis_der_y((N+1)*(N+2)/2,(N+1)*(N+2)/2))

allocate( basis_der_st_x((N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2*(N+1)))
allocate( basis_der_st_y((N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2*(N+1)))




phi=0
mue_lu=0
mue_lu_pivot=0
phi=generatePhi(N)

call generateMueLU(phi,N,mue_lu,mue_lu_pivot)

!Transformation matrices
call matrixToFile(N,"phi",phi,"_SWE_PATCH_ORDER*_SWE_PATCH_ORDER,_SWE_DG_DOFS")
call matrixToFile(N,"mue_lu",mue_lu,"_SWE_DG_DOFS+1,_SWE_DG_DOFS+1")
call matrixToFile(N,"mue_lu_pivot",mue_lu_pivot,"_SWE_DG_DOFS+1")

!Space stiffness an mass matrices

s_m=0
s_k_x=0
s_k_y=0
s_m_lu=0
s_m_lu_pivot=0

call generate_s_m(s_m,N)
s_m_lu = s_m
call ludcmp(s_m_lu,(N+2)*(N+1)/2, s_m_lu_pivot, code,info)
call generate_s_k_x(s_k_x,N)
call generate_s_k_y(s_k_y,N)

call generate_ss_k(ss_k_x,s_k_x,N)
call generate_ss_k(ss_k_y,s_k_y,N)
call generate_ss_k(ss_m,s_m,N)


! do i=0,N
!    do j=0,N
!    ss_k_x(i+1,i+1+j*M)=0
!    ss_k_y(i+1,i+1+j*M)=0
!    end do
! end do

call matrixToFile(N,"s_m"  ,ss_m  ,"_SWE_DG_DOFS,_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
call matrixToFile(N,"s_k_x",ss_k_x,"_SWE_DG_DOFS,_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
call matrixToFile(N,"s_k_y",ss_k_y,"_SWE_DG_DOFS,_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
call matrixToFile(N,"s_m_lu",s_m_lu,"_SWE_DG_DOFS,_SWE_DG_DOFS")
call matrixToFile(N,"s_m_lu_pivot",s_m_lu_pivot,"_SWE_DG_DOFS")

!predictor stiffness an mass matrices

tau=0
st_m=0
st_k_x=0
st_k_y=0
st_k_t=0

call generate_tau(tau,N)
call generate_st_m(st_m,s_m,tau,N)
call generate_st_k(st_k_x,s_k_x,tau,N)
call generate_st_k(st_k_y,s_k_y,tau,N)
call generate_st_k_t(st_k_t,s_m,N)
!call generate_st_w(st_w,s_m,N)

! st_w_k_t=st_w-st_k_t
st_w_k_t=st_k_t

 st_w_k_t_1_1_lu=st_w_k_t( M+1 : M*(N+1) , M+1 : M*(N+1) )
 st_w_k_t_1_0   =st_w_k_t( M+1 : M*(N+1) , 1 : M)

call ludcmp(st_w_k_t_1_1_lu,(N+2)*(N+1)/2*N, st_w_k_t_1_1_lu_pivot, code,info)

st_m_lu = st_m
call ludcmp(st_m_lu,(N+1)*(M), st_m_lu_pivot, code,info)


call matrixToFile(N,"st_m",st_m(M+1:M*(N+1),:)    ,"_SWE_DG_DOFS*(_SWE_DG_ORDER),_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
call matrixToFile(N,"st_k_x",st_k_x(M+1:M*(N+1),:),"_SWE_DG_DOFS*(_SWE_DG_ORDER),_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
call matrixToFile(N,"st_k_y",st_k_y(M+1:M*(N+1),:),"_SWE_DG_DOFS*(_SWE_DG_ORDER),_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")

 call matrixToFile(N,"st_w_k_t_1_0",st_w_k_t_1_0,"_SWE_DG_DOFS*(_SWE_DG_ORDER),(_SWE_DG_DOFS)")
 call matrixToFile(N,"st_w_k_t_1_1_lu",st_w_k_t_1_1_lu,"_SWE_DG_DOFS*(_SWE_DG_ORDER),_SWE_DG_DOFS*(_SWE_DG_ORDER)")
 call matrixToFile(N,"st_w_k_t_1_1_lu_pivot",st_w_k_t_1_1_lu_pivot ,"_SWE_DG_DOFS*(_SWE_DG_ORDER)")

b_m_1 = 0
b_m_2 = 0
b_m_3 = 0


call generate_b_m_1(b_m_1,N)
call generate_b_m_2(b_m_2,N)
call generate_b_m_3(b_m_3,N)

call matrixToFile(N,"b_m_1",b_m_1,"_SWE_DG_DOFS,_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
call matrixToFile(N,"b_m_2",b_m_2,"_SWE_DG_DOFS,_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
call matrixToFile(N,"b_m_3",b_m_3,"_SWE_DG_DOFS,_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")

write(n_nodes_string,"(I0,A)"),N_nodes**3,",_SWE_DG_DOFS*(_SWE_DG_ORDER+1)"
call generate_l2_projection(st_gl_node_vals,st_gl_weights,N)

call matrixToFile(N,"st_gl_node_vals",st_gl_node_vals,trim(n_nodes_string))

write(n_nodes_string,"(A,I0)"),"_SWE_DG_DOFS*(_SWE_DG_ORDER+1),",N_nodes**3
call matrixToFile(N,"st_gl_weights",st_gl_weights,trim(n_nodes_string))

s_der_x_gl_node_vals=0
s_der_y_gl_node_vals=0
st_der_x_gl_node_vals=0
st_der_y_gl_node_vals=0

call generate_l2_projection_der(s_der_x_gl_node_vals,s_der_y_gl_node_vals,st_der_x_gl_node_vals,st_der_y_gl_node_vals,N)

write(n_nodes_string,"(I0,A)"),N_nodes**3,",_SWE_DG_DOFS*(_SWE_DG_ORDER+1)"
call matrixToFile(N,"st_der_x_gl_node_vals",st_der_x_gl_node_vals,trim(n_nodes_string))
call matrixToFile(N,"st_der_y_gl_node_vals",st_der_y_gl_node_vals,trim(n_nodes_string))

write(n_nodes_string,"(I0,A)"),N_nodes**2,",_SWE_DG_DOFS"
call matrixToFile(N,"s_der_x_gl_node_vals",s_der_x_gl_node_vals,trim(n_nodes_string))
call matrixToFile(N,"s_der_y_gl_node_vals",s_der_y_gl_node_vals,trim(n_nodes_string))

call generate_l2_projection_bnd(bnd_gl_node_vals,bnd_gl_weights,N)

write(n_nodes_string,"(I0,A)"),N_nodes**2,",(_SWE_DG_ORDER+1)*(_SWE_DG_ORDER+1)"
call matrixToFile(N,"bnd_gl_node_vals",bnd_gl_node_vals,trim(n_nodes_string))
write(n_nodes_string,"(A,I0)"),"(_SWE_DG_ORDER+1)*(_SWE_DG_ORDER+1),",N_nodes**2
call matrixToFile(N,"bnd_gl_weights",bnd_gl_weights,trim(n_nodes_string))

call matrixToFile(N,"st_m_lu",st_m_lu,"_SWE_DG_DOFS*(_SWE_DG_ORDER+1),_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
call matrixToFile(N,"st_m_lu_pivot",st_m_lu_pivot,"_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")


call generate_basis_st_der_y(basis_der_st_y,N)
call generate_basis_st_der_x(basis_der_st_x,N)
call generate_basis_der_y(basis_der_y,N)
call generate_basis_der_x(basis_der_x,N)

call matrixToFile(N,"basis_der_st_y",basis_der_st_y,"_SWE_DG_DOFS*(_SWE_DG_ORDER+1),_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
call matrixToFile(N,"basis_der_st_x",basis_der_st_x,"_SWE_DG_DOFS*(_SWE_DG_ORDER+1),_SWE_DG_DOFS*(_SWE_DG_ORDER+1)")
call matrixToFile(N,"basis_der_y",basis_der_y,"_SWE_DG_DOFS,_SWE_DG_DOFS")
call matrixToFile(N,"basis_der_x",basis_der_x,"_SWE_DG_DOFS,_SWE_DG_DOFS")

end program matrices
