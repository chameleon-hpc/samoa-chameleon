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

  !------generate Common (Modal and Nodal) matrices ------!
  
  subroutine generate_limiter_matrices(N)
    !---limiter matrices---!
    real(kind=GRID_SR),allocatable:: phi(:,:)
    real(kind=GRID_SR),allocatable:: mue_lu(:,:)
    integer           ,allocatable:: mue_lu_pivot(:)

    phi=0
    mue_lu=0
    mue_lu_pivot=0
    phi=generatePhi(N)
    
    call generateMueLU(phi,N,mue_lu,mue_lu_pivot)
    
    !Transformation matrices
    call matrixToFile(N,"phi",phi,"_SWE_PATCH_ORDER*_SWE_PATCH_ORDER,_SWE_DG_DOFS")
    call matrixToFile(N,"mue_lu",mue_lu,"_SWE_DG_DOFS+1,_SWE_DG_DOFS+1")
    call matrixToFile(N,"mue_lu_pivot",mue_lu_pivot,"_SWE_DG_DOFS+1")
  end subroutine generate_common_matrices
  
  !------generate Nodal matrices ------!
  subroutine  generate_nodal_matrices(N)
    !---predictor matrices---!
    real(kind=GRID_SR),allocatable:: basis_der_x(:,:)
    real(kind=GRID_SR),allocatable:: basis_der_y(:,:)
    real(kind=GRID_SR),allocatable:: basis_der_st_x(:,:)
    real(kind=GRID_SR),allocatable:: basis_der_st_y(:,:)


    real(kind=GRID_SR),allocatable:: st_m  (:,:)
    real(kind=GRID_SR),allocatable:: st_k_x(:,:)
    real(kind=GRID_SR),allocatable:: st_k_y(:,:)

    real(kind=GRID_SR),allocatable:: st_w_k_t         (:,:)
    real(kind=GRID_SR),allocatable:: st_w_k_t_1_0     (:,:)
    
    real(kind=GRID_SR),allocatable:: st_w_k_t_1_1_lu  (:,:)
    integer           ,allocatable:: st_w_k_t_1_1_lu_pivot(:)


    !---solver matrices---!
    real(kind=GRID_SR),allocatable:: s_m(:,:)
    real(kind=GRID_SR),allocatable:: s_m_lu(:,:)
    real(kind=GRID_SR),allocatable:: s_k_x(:,:)
    real(kind=GRID_SR),allocatable:: s_k_y(:,:)

    !---boundary matrices---!
    real(kind=GRID_SR),allocatable:: b_m_1(:,:)
    real(kind=GRID_SR),allocatable:: b_m_2(:,:)
    real(kind=GRID_SR),allocatable:: b_m_3(:,:)


  end subroutine generate_nodal_matrices
  
real(kind=GRID_SR),allocatable:: temp(:,:)

integer               ,allocatable:: s_m_lu_pivot(:)
real(kind=GRID_SR),allocatable:: ss_m(:,:)
real(kind=GRID_SR),allocatable:: ss_k_x(:,:)
real(kind=GRID_SR),allocatable:: ss_k_y(:,:)

real(kind=GRID_SR),allocatable:: tau(:,:)

real(kind=GRID_SR),allocatable:: st_m_lu  (:,:)
integer               ,allocatable:: st_m_lu_pivot  (:)
real(kind=GRID_SR),allocatable:: st_w  (:,:)
real(kind=GRID_SR),allocatable:: st_k_t(:,:)

real(kind=GRID_SR),allocatable:: s_der_x_gl_node_vals(:,:)
real(kind=GRID_SR),allocatable:: s_der_y_gl_node_vals(:,:)
real(kind=GRID_SR),allocatable:: st_der_x_gl_node_vals(:,:)
real(kind=GRID_SR),allocatable:: st_der_y_gl_node_vals(:,:)


  !------generate Modal matrices ------!

  subroutine  generate_modal_matrices(N)

    !---predictor matrices---!
    real(kind=GRID_SR),allocatable:: st_gl_node_vals(:,:)
    real(kind=GRID_SR),allocatable:: st_gl_weights(:,:)

    !---boundary matrices---!
    real(kind=GRID_SR),allocatable:: bnd_gl_node_vals(:,:)
    real(kind=GRID_SR),allocatable:: bnd_gl_weights(:,:)


    

  end subroutine




  

  !------Single matrices-----!

  ! space mass matrix
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

   ! space stiffness matrix x
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

   ! space stiffness matrix y
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

   subroutine generate_ss_k(ss_k,s_k,N)
    integer :: i,N
    real(kind=GRID_SR), intent (in) :: s_k((N+1)*(N+2)/2,(N+1)*(N+2)/2)
    real(kind=GRID_SR) :: t(1,N+1)
    real(kind=GRID_SR), intent (inout) :: ss_k((N+1)*(N+2)/2,(N+1)*(N+2)/2*(N+1))

    do i=0,N
       t(1,i+1)=gl_quadr_1d(basis_polynomial_1d_array(gl_nodes,N,i),0.0_GRID_SR,1.0_GRID_SR)
       !print*,t(1,i+1)
    end do

    ss_k = kroneckerMultiplication(t,s_k)
      
   end subroutine generate_ss_k

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

 !routines to generate predictor matrices

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

    subroutine generate_tau(tau,N)
      real(kind=GRID_SR),intent(inout) :: tau(N+1,N+1)
      integer ,intent(in) :: N
      integer  :: i,j

      do i=0,N
        do j=0,N
           tau(i+1,j+1)=gl_quadr_1d(basis_polynomial_1d_array(gl_nodes,N,i)*basis_polynomial_1d_array(gl_nodes,N,j),0.0_GRID_SR,1.0_GRID_SR)
!           tau(i,j)= gl_quadr_1d(basis_polynomial_1d_array(nodes(1,:),N,i)*basis_polynomial_1d_array(nodes(2,:),N,j),0.0_GRID_SR,1q0)
         end do
      end do

    end subroutine generate_tau

    subroutine generate_st_k(st_k,s_k,tau,N)
      integer,intent(in) :: N
      real(kind=GRID_SR),intent(in) :: tau(N+1,N+1),s_k((N+1)*(N+2)/2,(N+1)*(N+2)/2)
      real(kind=GRID_SR),intent(inout) :: st_k((N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2*(N+1))

      st_k=kroneckerMultiplication(tau,transpose(s_k))
            
    end subroutine generate_st_k

   subroutine generate_st_m(st_m,s_m ,tau,N)      
     integer,intent(in) :: N
     real(kind=GRID_SR),intent(in) :: tau(N+1,N+1),s_m((N+1)*(N+2)/2,(N+1)*(N+2)/2)
     real(kind=GRID_SR),intent(inout) :: st_m((N+1)*(N+2)/2*(N+1),(N+1)*(N+2)/2*(N+1))
     st_m=kroneckerMultiplication(tau,s_m)

   end subroutine generate_st_m


!generate ref triangle vertices
  subroutine generateMueLU(phi,N,mue_lu,ipiv)
    integer:: N,M,j,k,l,code
    real(kind=GRID_SR) :: phi((N*2 + 1)*(N*2 + 1),((N+1)*(N+2)/2)),b((N+1)*(N+2)/2+1),c((N*2 + 1)*(N*2 + 1))
    real(kind=GRID_SR) :: mue(((N+1)*(N+2)/2)+1,(N+1)*(N+2)/2+1),mue_lu(((N+1)*(N+2)/2)+1,(N+1)*(N+2)/2+1)
    integer            :: ipiv(((N+1)*(N+2)/2)+1), info

    M = (N+1)*(N+2)/2
    mue(1:M,1:M) = matmul(transpose(phi),phi) * 2.0_GRID_SR
    mue(M+1,M+1) = 0.0_GRID_SR

    do k = 0,N
       do j = 0,N-k
          l=M-(N+1-k)*(N+2-k)/2+j+1
          mue(M+1,l) = gl_quadr_2d(basis_polynomial_array(get_gl_nodes_2d((/ 0.0_GRID_SR,0.0_GRID_SR /),(/ 1.0_GRID_SR,1.0_GRID_SR /)),N,j,k),0.0_GRID_SR,1.0_GRID_SR)
          mue(l,M+1) = mue(M+1,l)
       end do
    end do
    mue_lu=mue
    
    call ludcmp(mue_lu,M+1,ipiv,code,info)

  end subroutine generateMueLU


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

  function generatePhi(N) result(phi)
    real(kind=GRID_SR):: nodes_2d(2,N_nodes*N_nodes),a(2),b(2)
    integer :: N_s,N_s_squared,M,i,j,k,N,l
    real(kind=GRID_SR) :: triangles(2,(N*2 + 1)*(N*2 + 1)*3)
    real(kind=GRID_SR) :: phi((N*2 + 1)*(N*2 + 1), (N+1)*(N+2)/2)
    logical :: lower

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

  end function generatePhi

!integer,Parameter  ::GRID_SR=GRID_SR

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
                                            basis_polynomial_array(nodes,N,i_2,j_2),0.0_GRID_SR,1.0_GRID_SR)*&
                                            time_int* &
                                            sqrt(2.0_GRID_SR)

              end do
           end do
        end do
     end do
  end do


end subroutine generate_b_m_2

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

end subroutine

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

end subroutine

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

end subroutine

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

end subroutine



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

end module gen_matrices


