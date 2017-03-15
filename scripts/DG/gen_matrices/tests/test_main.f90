Program test_main
use gl_integration
use asserts
use elementary
use gen_matrices
use test_defs

call set_basis('BERNSTEIN')
call testFactorial()
call testkronecker()
call testBernsteinpolynomial()
call testGL()
call test_triangles()
call test_phi_and_mue_lu()
call test_space_dg_matrices()
call test_space_time_predictor_matrices()
call test_boundary_matrices()
call test_basis_der_matrices()
call endTests()

contains


  ! TEST GL_Quadrature

  subroutine testGL()
    implicit none
    integer :: i,j,count
    real(kind=GRID_SR)              ::a(2),b(2)
    real(kind=GRID_SR),dimension(5) ::testvalues_1d
    real(kind=GRID_SR),dimension(5*5) ::testvalues_2d
    real(kind=GRID_SR)                ::nodes_2d(2,25)

   call initGLIntegration(4)

    !test GL1D
    testvalues_1d=(gl_nodes*6-3)**2*(-1)
    call assertEqual(-18.0_GRID_SR,gl_quadr_1d(testvalues_1d,-3.0_GRID_SR,3.0_GRID_SR),"Gl quadarture 1d f(x)=-x^2 on -3 to 3")
    testvalues_1d= (/(0,I=1,5)/)
    call assertEqual(-0.0_GRID_SR,gl_quadr_1d(testvalues_1d,-3.0_GRID_SR,3.0_GRID_SR),"Gl quadarture 1d f(x)=0 on -3 to 3")

    !test GL2D
    testvalues_2d= (/(1,I=1,5*5)/)
    call assertEqual(0.5_GRID_SR,gl_quadr_2d(testvalues_2d,0.0_GRID_SR,1.0_GRID_SR)  ,"Gl quadrature 2d f(x)=1 on Triangle (0,0) (0,1) (1,0)")
    call assertEqual(0.125_GRID_SR,gl_quadr_2d(testvalues_2d,0.0_GRID_SR,0.5_GRID_SR),"Gl quadratur 2d f(x)=1 on Triangle (0,0) (0,0.5) (0.5,0)")

    a=(/ (1.0_GRID_SR/2.0_GRID_SR),(1.0_GRID_SR/2.0_GRID_SR) /)
    b=(/ (3.0_GRID_SR/4.0_GRID_SR),(3.0_GRID_SR/4.0_GRID_SR) /)
    nodes_2d=get_gl_nodes_2d(a,b)
    testvalues_2d=nodes_2d(1,:)*nodes_2d(2,:)
    
    call assertTolerance(0.0105794270833333333333333333333333333_GRID_SR,gl_quadr_2d(testvalues_2d,0.5_GRID_SR,0.75_GRID_SR),&
"Gl quadratur 2d f(x,y)=x*y on Triangle (0.5,0.5) (0.5,0.75) (0.75,0.5)")

    a=(/ (1.0_GRID_SR/9.0_GRID_SR),(1.0_GRID_SR/3.0_GRID_SR) /)
    b=(/ (2.0_GRID_SR/9.0_GRID_SR),(4.0_GRID_SR/9.0_GRID_SR) /)

    nodes_2d=get_gl_nodes_2d(a,b)
    testvalues_2d = basis_polynomial_array(nodes_2d,4,2,0)

    Call assertTolerance(0.000191052879498069086377101754161484216_GRID_SR &
,gl_quadr_2d(testvalues_2d,1.0_GRID_SR/9.0_GRID_SR,2.0_GRID_SR/9.0_GRID_SR), &
"Gl quadratur 2d f(x)=bernstein(3) for N=4 on Triangle (1/9,4/9) (2/9,4/9) (1/9,1/3)")

    nodes_2d=get_gl_nodes_2d(a,b,.false.) 
    testvalues_2d = basis_polynomial_array(nodes_2d,4,2,0)
    call assertTolerance(0.00021231582307976489080318103671589759415_GRID_SR&
         ,gl_quadr_2d(testvalues_2d,1.0_GRID_SR/9.0_GRID_SR,2.0_GRID_SR/9.0_GRID_SR), &
         "Gl quadratur 2d f(x)=bernstein(3) for N=4 on Triangle (1/9,4/9) (2/9,4/9) (2/9,1/3)")
    

  end subroutine testGL

  subroutine testkronecker()
    real(kind=GRID_SR) :: a(2,2),b(2,2),c(4,4),c_ref(16)
    
    c_ref=(/1.0_GRID_SR,1.0_GRID_SR,3.0_GRID_SR,3.0_GRID_SR,1.0_GRID_SR,1.0_GRID_SR,3.0_GRID_SR,3.0_GRID_SR,2.0_GRID_SR,2.0_GRID_SR,4.0_GRID_SR,4.0_GRID_SR,2.0_GRID_SR,2.0_GRID_SR,4.0_GRID_SR,4.0_GRID_SR/)

    a=reshape((/1,3,2,4/),(/2,2/))
    b=reshape((/1,1,1,1/),(/2,2/))
    c = kroneckerMultiplication(a,b)
    call assertEqualArray(Real(c_ref,kind=GRID_SR),reshape(c,(/16/)),"Test Kronecker Multiplication")
  end subroutine testkronecker
! TEST Factorial
  subroutine testFactorial
    call assertEqual(3628800,faculty(10),"Test Factorial 10")
    call assertEqual(1,faculty(0),"Test Factorial 0")
    call assertEqual(-1,faculty(-100),"Test Factorial negative")
  end subroutine testFactorial

  subroutine testBernsteinpolynomial
    implicit none

    call assertEqual(.1318359375_GRID_SR,bernstein_polynomial_1d(0.25_GRID_SR,6,3),"1D Bernsteinpolynomial for N=6, l=3 at x=0,25")

    call assertEqual(0.0_GRID_SR,bernstein_polynomial(0.0_GRID_SR,1.0_GRID_SR,4,1,1),"Bernsteinpolynomial for N=4, l=7 at x=0.25, y=0.25")
    call assertEqual(0.125_GRID_SR,bernstein_polynomial(0.25_GRID_SR,0.25_GRID_SR,4,0,1),"Bernsteinpolynomial for N=4, l=6 at x=0.25, y=0.25")
    call assertTolerance(0.04305422068841828748647865654489403932459_GRID_SR,&
bernstein_polynomial(0.8523456_GRID_SR,0.1237648_GRID_SR,10,7,2),"Bernsteinpolynomial for N=10, l=29 at x=0.8523456, y=0.1237648")

    call assertTolerance(0.585e-2_GRID_SR, bernstein_polynomial_der_x(0.25_GRID_SR,0.1_GRID_SR,6,1,3),"Derivative x Bernsteinpolynomial x=0.25 y=0.1 N=6 j=1 k=3")
    call assertTolerance(-0.2535e-1_GRID_SR, bernstein_polynomial_der_x(0.25_GRID_SR,0.1_GRID_SR,6,0,3),"Derivative x Bernsteinpolynomial x=0.25 y=0.1 N=6 j=0 k=3")

    call assertTolerance(1.70625e-1_GRID_SR, bernstein_polynomial_der_y(0.25_GRID_SR,0.1_GRID_SR,6,1,3),"Derivative y Bernsteinpolynomial x=0.25 y=0.1 N=6 j=1 k=3")

  end subroutine testBernsteinpolynomial

  subroutine test_phi_and_mue_lu
    implicit none
    real(kind=GRID_SR),Allocatable:: phi(:,:)
    real(kind=GRID_SR),Allocatable:: mue_lu(:,:),b_16(:)
    integer,Allocatable :: ipiv(:)
    allocate(phi(9,3))
    call generatePhi(phi,1)
    call assertEqualArray(Real(phi_ref_1,kind=GRID_SR),reshape(transpose(phi) ,(/9*3/)),"Compare Phi for N=1")
    deallocate(phi)
    allocate(phi(81,15))
    call generatePhi(phi,4)
    call assertEqualArray(Real(phi_ref_4,kind=GRID_SR),reshape(transpose(phi) ,(/81*15/)),"Compare Phi for N=4")

    allocate(mue_lu(16,16))
    allocate(ipiv(16))

    call generateMueLU(phi,4,mue_lu,ipiv)

    allocate(b_16(16))

    b_16(1:15)= 2* matmul(transpose(phi),c_81)
    b_16(16)=sum(c_81)
    b_16=b_16/(2.0_GRID_SR*81.0_GRID_SR)
    call lubksb(mue_lu,16,ipiv,b_16)

    call assertEqualArray(Real(ref_dofs_4,GRID_SR),b_16(1:15),"Test Mue LU for N=4")
    

  end subroutine test_phi_and_mue_lu
  
  subroutine test_space_dg_matrices()
    real(kind=GRID_SR) :: s_m(15,15), s_k_y(15,15), s_k_x(15,15),ss_k_x(15,75), ss_k_y(15,75)
    
    call generate_s_m(s_m,4)
    call generate_s_k_x(s_k_x,4)
    call generate_s_k_y(s_k_y,4)

    call assertEqualArray(Real(s_m_ref_4,kind=GRID_SR),reshape(s_m,(/15*15/)),"Test Mass matrix N=4")
    call assertEqualArray(Real(s_k_y_ref_4,kind=GRID_SR),reshape(s_k_y,(/15*15/)),"Test Stiffness y matrix N=4")
    call assertEqualArray(Real(s_k_x_ref_4,kind=GRID_SR),reshape(s_k_x,(/15*15/)),"Test Stiffness x matrix N=4")

    call generate_ss_k(ss_k_x,s_k_x,4)
    call generate_ss_k(ss_k_y,s_k_y,4)

    call assertEqualArray(Real(ss_k_x_ref_4,kind=GRID_SR),reshape(ss_k_x,(/15*75/)),"Test Stiffness x matrix N=4")
    call assertEqualArray(Real(ss_k_y_ref_4,kind=GRID_SR),reshape(ss_k_y,(/15*75/)),"Test Stiffness y matrix N=4")

    
  end subroutine test_space_dg_matrices


  
  subroutine test_space_time_predictor_matrices()
    real(kind=GRID_SR) :: s_m(15,15), s_k_y(15,15), s_k_x(15,15), tau(5,5)
    real(kind=GRID_SR)   :: w(75,75),st_k_x(75,75),st_k_y(75,75),st_m(75,75),st_k_t(75,75)

!    call initGLIntegration(4)

     call generate_s_m(s_m,4)
     call generate_st_w(w,s_m,4)

     call assertEqualArray(Real(w_ref_4,kind=GRID_SR),reshape(w,(/75*75/)),"Test W-Matrix N=4")

     call generate_s_k_x(s_k_x,4)
     call generate_tau(tau,4)

     call generate_st_k(st_k_x,s_k_x,tau,4)

     call assertEqualArray(Real(st_k_x_ref_4,kind=GRID_SR),reshape(st_k_x,(/75*75/)),"Test Space-Time K-x-Matrix N=4")

     call generate_s_k_y(s_k_y,4)
     call generate_st_k(st_k_y,s_k_y,tau,4)

     call assertEqualArray(Real(st_k_y_ref_4,kind=GRID_SR),reshape(st_k_y,(/75*75/)),"Test Space-Time K-y-Matrix N=4")

     call generate_st_m(st_m,s_m,tau,4)

     call assertEqualArray(Real(st_m_ref_4,kind=GRID_SR),reshape(st_m,(/75*75/)),"Test Space-Time Mass-Matrix N=4")

     call generate_st_k_t(st_k_t,s_m,4)

     call assertEqualArray(Real(st_k_t_ref_4,kind=GRID_SR),reshape(transpose(st_k_t),(/75*75/)),"Test Space-Time K-t-Matrix N=4")
     
  end subroutine test_space_time_predictor_matrices

  subroutine test_boundary_matrices()
    real(kind=kind(1_GRID_SR)) :: b_m_1(15,75)
    real(kind=kind(1_GRID_SR)) :: b_m_2(15,75)
    real(kind=kind(1_GRID_SR)) :: b_m_3(15,75)
    call generate_b_m_1(b_m_1,4)
    call generate_b_m_2(b_m_2,4)
    call generate_b_m_3(b_m_3,4)

    call assertEqualArray(Real(b_m_1_ref_4,kind=GRID_SR),reshape(b_m_1,(/15*75/)),"Test Boundary Mass Matrix 1 N=4")
    call assertEqualArray(Real(b_m_2_ref_4,kind=GRID_SR),reshape(b_m_2,(/15*75/)),"Test Boundary Mass Matrix 2 N=4")
    call assertEqualArray(Real(b_m_3_ref_4,kind=GRID_SR),reshape(b_m_3,(/15*75/)),"Test Boundary Mass Matrix 3 N=4")
  end subroutine

  subroutine test_triangles()
    real(kind=GRID_SR) :: triangles(2,81*3)
    triangles = generateTriangles(4)
    call assertEqualArray(Real(triangle_ref,kind=GRID_SR),reshape(triangles,(/2*81*3/)),"Test triangles N=4")
  end subroutine test_triangles

  subroutine test_basis_der_matrices()
    real(kind=GRID_SR) :: basis_der_x(15,15)
    real(kind=GRID_SR) :: basis_der_y(15,15)
    
    call generate_basis_der_x(basis_der_x,4)
    call assertEqualArray(Real(basis_der_x_ref_4,kind=GRID_SR),reshape(basis_der_x,(/15*15/)),"Test basis derivative-matrix x")
    call generate_basis_der_y(basis_der_y,4)
    call assertEqualArray(Real(basis_der_y_ref_4,kind=GRID_SR),reshape(basis_der_y,(/15*15/)),"Test basis derivative-matrix y")

  end subroutine

end program test_main
