module elementary
  use gl_integration
  use dubiner
implicit none

integer :: basis
logical :: nodal

contains


  subroutine set_basis(basis_in,N)
    character(len=*) :: basis_in
    integer :: basis_int
    integer :: N
    
    select case(basis_in)
    case('BERNSTEIN')
       basis=1
       nodal=.TRUE.
    case('EQUIDISTAND')
       call gen_vandermonde_inv_equidistand(N)
       basis=2
       nodal=.TRUE.
    case('LGL')
       call gen_vandermonde_inv_lgl(N)
       basis=3
       nodal=.TRUE.
       print*, "Not ready yet"
       stop
    case default
       print*,"Basis doesn't exists, possible are:"
       print*,'BERNSTEIN'
       print*,'EQUIDISTAND'
       print*,'LGL'
    end select
    
  end subroutine set_basis


  function basis_polynomial_1d(t,N,l)
    implicit none
    real(kind=GRID_SR),intent(in) :: t
    real(kind=GRID_SR),allocatable :: basis_polynomial_1d
    integer,intent(in) :: N,l
    integer :: i

    select case(basis)
    case (1)
       basis_polynomial_1d=bernstein_polynomial_1d(t,N,l)
    case(2)
       basis_polynomial_1d=lagrange_polynomial_1d(t,N,l)
    end select


  end function basis_polynomial_1d
    

  function basis_polynomial(x,y,N,j,k)
    real(kind=GRID_SR) :: basis_polynomial
    real(kind=GRID_SR) ::x,y 
    integer ::N,j,k,i

    select case(basis)
    case (1)
       basis_polynomial = bernstein_polynomial(x,y,N,j,k)
    case (2)
       basis_polynomial = equ_nodalBasis(x,y,N,j,k)
    end select

  end function

  function basis_polynomial_1d_der(t,N,l)
    real(kind=GRID_SR) :: basis_polynomial_1d_der
    real(kind=GRID_SR) :: t
    integer :: N,l,i
    
    select case(basis)
    case (1)
       basis_polynomial_1d_der = bernstein_polynomial_1d_der(t,N,l)
    case (2)
       basis_polynomial_1d_der = lagrange_polynomial_1d_der(t,N,l)
    end select

  end function

  
  function basis_polynomial_der_x(x,y,N,j,k)
    real(kind=GRID_SR) :: basis_polynomial_der_x
    real(kind=GRID_SR) :: x,y 
    integer :: N,j,k,i
    select case(basis)
    case (1)
       basis_polynomial_der_x = bernstein_polynomial_der_x(x,y,N,j,k)
    case (2)
       basis_polynomial_der_x = equ_nodal_basis_der_x(x,y,N,j,k)
    end select
    
  end function basis_polynomial_der_x

 function basis_polynomial_der_y(x,y,N,j,k)
    real(kind=GRID_SR) ::basis_polynomial_der_y
    real(kind=GRID_SR) ::x,y 
    integer ::N,j,k,i

    select case(basis)
    case (1)
       basis_polynomial_der_y = bernstein_polynomial_der_y(x,y,N,j,k)
    case (2)
       basis_polynomial_der_y = equ_nodal_basis_der_y(x,y,N,j,k)
    end select
  end function basis_polynomial_der_y

  
  function lagrange_polynomial_1d(t,N,l)
    real(kind=GRID_SR) :: Nodes(0:N)
    real(kind=GRID_SR) :: t
    real(kind=GRID_SR) :: lagrange_polynomial_1d
    integer :: N,l,i,ind
    

    ! equidistand nodes !
    do i=0,N
       Nodes(i) = (i)/Real(N)
    end do
        
    lagrange_polynomial_1d=1

    do i=0,N
       if(i.eq.l) then
          cycle
       end if
       lagrange_polynomial_1d=lagrange_polynomial_1d*(t-Nodes(i))/(Nodes(l)-Nodes(i))
    end do   

  end function

  function lagrange_polynomial(x,y,N,j,k)
    real(kind=GRID_SR) :: lagrange_polynomial
    integer N,j,k
    real(kind=GRID_SR) :: x,y

    lagrange_polynomial=(lagrange_polynomial_1d(x,N,j)*lagrange_polynomial_1d(y,N,k))

  end function

  function lagrange_polynomial_1d_der(t,N,l)
    real(kind=GRID_SR) :: lagrange_polynomial_1d_der
    real(kind=GRID_SR) :: Nodes(0:N)
    integer N,l,i,j
    real(kind=GRID_SR) :: t
    real(kind=GRID_SR) :: temp

    ! equidistand nodes !
    do i=0,N
       Nodes(i) = (i)/Real(N)
    end do

    
    lagrange_polynomial_1d_der=0
    do i=0,N
       if(i.eq.l) then
          cycle
       end if
       temp=1.0_GRID_SR
       do j=0,N

          if(j.eq.l .or. j.eq.i) then
             cycle
          end if
          temp=temp*(t-Nodes(j))/(Nodes(l)-Nodes(j))
       end do

       lagrange_polynomial_1d_der =  lagrange_polynomial_1d_der + (temp)* (1.0_GRID_SR/(Nodes(l)-Nodes(i)))
    end do
    

  end function

  function lagrange_polynomial_der_x(x,y,N,j,k)
    real(kind=GRID_SR) :: Nodes(N+1)
    real(kind=GRID_SR) :: lagrange_polynomial_der_x
    integer N,j,k,i
    real(kind=GRID_SR) :: x,y
    
    lagrange_polynomial_der_x=lagrange_polynomial_1d_der(x,N,j)*lagrange_polynomial_1d(y,N,k)

  end function

  function lagrange_polynomial_der_y(x,y,N,j,k)
    real(kind=GRID_SR) :: lagrange_polynomial_der_y
    real(kind=GRID_SR) :: Nodes(N+1)
    integer N,j,k,i
    real(kind=GRID_SR) :: x,y
    
    lagrange_polynomial_der_y=lagrange_polynomial_1d_der(y,N,k)*lagrange_polynomial_1d(x,N,j)

  end function


  function bernstein_polynomial(x,y,N,j,k)
    implicit none
    real(kind=GRID_SR), intent(in) :: x,y
    real(kind=GRID_SR)             :: bernstein_polynomial,facp,pnt_data
    integer,intent(in)             :: N,k,j

    facp=REAL(faculty(N),GRID_SR)/REAL((faculty(k)*faculty(j)*faculty(N-k-j)),GRID_SR)

    pnt_data=x**j * y**k * (1.0_GRID_SR-x-y) ** (N-k-j)

    bernstein_polynomial= facp * pnt_data

    return
  end function bernstein_polynomial

!parital derivative in x direction
  function bernstein_polynomial_der_x(x,y,N,j,k)
    implicit none
    real(kind=GRID_SR), intent(in) :: x,y
    real(kind=GRID_SR)             :: bernstein_polynomial_der_x
    integer,intent(in)             :: N,k,j

    if(j+k.eq.N) then
       if(j.eq.0) then
          bernstein_polynomial_der_x=0
       else
          bernstein_polynomial_der_x= x**(j-1) * y**k * j * faculty(N)/REAL(faculty(k)*faculty(j),GRID_SR)
       end if
    else if(j.ne.0) then
       bernstein_polynomial_der_x=-x**(j-1)*y**k*(1-x-y)**(N-j-k-1)*(-j+j*y+x*N-x*k)*faculty(N)/REAL(faculty(N-j-k)*faculty(k)*faculty(j),GRID_SR)
    else if (k.ne.N) then
       bernstein_polynomial_der_x=-N*bernstein_polynomial(x,y,N-1,0,k)
    else 
       bernstein_polynomial_der_x=0
    end if
  end function bernstein_polynomial_der_x

  function bernstein_polynomial_der_y(x,y,N,j,k)
    real(kind=GRID_SR) :: bernstein_polynomial_der_y
    real(kind=GRID_SR) :: x,y
    integer :: N,j,k
    bernstein_polynomial_der_y=bernstein_polynomial_der_x(y,x,N,k,j)
  end function bernstein_polynomial_der_y


  function bernstein_polynomial_1d(t,N,l)
    implicit none
    real(kind=GRID_SR),intent(in) :: t
    real(kind=GRID_SR) :: bernstein_polynomial_1d
    integer,intent(in) :: N,l
    bernstein_polynomial_1d=(1.0_GRID_SR - t) **(N-l) * t ** (l) * REAL(faculty(N)/REAL(faculty(l)*faculty(N-l),GRID_SR),GRID_SR)

  end function bernstein_polynomial_1d

  function bernstein_polynomial_1d_der(t,N,l)
    implicit none
    real(kind=GRID_SR),intent(in) :: t
    real(kind=GRID_SR) :: bernstein_polynomial_1d_der
    integer,intent(in) :: N,l

    if (l.eq.0) then
       bernstein_polynomial_1d_der = (1.0_GRID_SR - t) **(N-1) * Real(-N,Grid_SR)
    else if (l.eq.N) then
       bernstein_polynomial_1d_der = t ** (N-1) * Real(N,Grid_SR)
    else
       bernstein_polynomial_1d_der = (1.0_GRID_SR - t) **(N-l-1) * t ** (l-1) * REAL(faculty(N)/REAL(faculty(l)*faculty(N-l),GRID_SR),GRID_SR)*(Real(l,GRID_SR)-Real(N,Grid_SR)*t)
    end if
  end function bernstein_polynomial_1d_der
  
  

  function basis_polynomial_1d_array(t,N,l)
    implicit none
    real(kind=GRID_SR),intent(in) :: t(:)
    real(kind=GRID_SR),allocatable :: basis_polynomial_1d_array(:)
    integer,intent(in) :: N,l
    integer :: i

    allocate(basis_polynomial_1d_array(size(t)))
    do i=1,size(t)
       basis_polynomial_1d_array(i)=basis_polynomial_1d(t(i),N,l)
    end do

  end function basis_polynomial_1d_array


  function basis_polynomial_1d_der_array(t,N,l)
    implicit none
    real(kind=GRID_SR),intent(in) :: t(:)
    real(kind=GRID_SR),allocatable :: basis_polynomial_1d_der_array(:)
    integer,intent(in) :: N,l
    integer :: i

    allocate(basis_polynomial_1d_der_array(size(t)))

    do i=1,size(t)
       basis_polynomial_1d_der_array(i)=basis_polynomial_1d_der(t(i),N,l)
    end do

  end function basis_polynomial_1d_der_array



  function basis_polynomial_array(nodes,N,j,k)
    integer,intent(in)               :: N,k,j
    real(kind=GRID_SR),intent(in)    :: nodes(:,:)
    real(kind=GRID_SR),allocatable                 :: basis_polynomial_array(:)
    integer                          :: i

    allocate(basis_polynomial_array(1:size(nodes,2)))

    do i=1,size(nodes,2)
       basis_polynomial_array(i) = basis_polynomial(nodes(1,i),nodes(2,i),N,j,k)
    end do

  end function basis_polynomial_array


  function basis_polynomial_der_x_array(nodes,N,j,k)
    integer,intent(in)               :: N,k,j
    real(kind=GRID_SR),intent(in)    :: nodes(:,:)
    real(kind=GRID_SR),allocatable   :: basis_polynomial_der_x_array(:)
    integer                          :: i

    allocate(basis_polynomial_der_x_array(1:size(nodes,2)))

    do i=1,size(nodes,2)
       basis_polynomial_der_x_array(i) = basis_polynomial_der_x(nodes(1,i),nodes(2,i),N,j,k)
    end do
  end function basis_polynomial_der_x_array

  function basis_polynomial_der_y_array(nodes,N,j,k)
    integer,intent(in)               :: N,k,j
    real(kind=GRID_SR),intent(in)    :: nodes(:,:)
    real(kind=GRID_SR),allocatable   :: basis_polynomial_der_y_array(:)
    integer                          :: i

    allocate(basis_polynomial_der_y_array(1:size(nodes,2)))

    do i=1,size(nodes,2)
       basis_polynomial_der_y_array(i) = basis_polynomial_der_y(nodes(1,i),nodes(2,i),N,j,k)
    end do
  end function basis_polynomial_der_y_array


  recursive function faculty(I,temp) result(fac)
    !
    implicit none
    integer, intent(in) :: I
    integer, optional   :: temp
    integer             :: fac,t

    if(I.lt.0) then
       fac=-1
       return
    end if
    if (.not.present(temp)) then
       t = 1
    else
       t=temp
    end if

    if (I.eq.0) then
       fac=t
    else
       fac=faculty(I-1,t*I)
    end if
    
    return
  end function faculty

  function kroneckerMultiplication(a,b) result(c)
    real(kind=GRID_SR) :: a(:,:), b(:,:)
    real(kind=GRID_SR),Allocatable :: c(:,:)
    integer :: i,j

    allocate(c(size(a,1)*size(b,1),size(a,2)*size(b,2)))

    do i=1,size(a,1)
       do j=1,size(a,2)
       c((i-1)*size(b,1)+1:i*size(b,1),(j-1)*size(b,2)+1:j*size(b,2))=a(i,j)*b
    end do
    end do

  end function kroneckerMultiplication
end module elementary
