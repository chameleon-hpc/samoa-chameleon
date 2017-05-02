module gl_integration
!
implicit none

integer,Parameter :: GRID_SR= kind(1.0d0)
real(kind=GRID_SR),allocatable ::gl_nodes(:),gl_weights(:)
integer              :: N_nodes

contains

!nodeValues ordering
!lower

!...
!7---8---9\
!4----5----6\  
!1-----2-----3\
!a             b

!upper
!a            b
!\1-----2-----3|
!  \4----5----6|  
!    \7---8---9|
!            ...
!


  !Takes an array of values of f(x) on the integration nodes and in interval [a,b]
  !and returns int_a^b f(x) dx

  function gl_quadr_2d(nodeValues,a,b) result(int)
    real(kind=GRID_SR),dimension(N_nodes*N_nodes) :: nodeValues
    real(kind=GRID_SR),intent(in)   :: a,b
    real(kind=GRID_SR)              :: int,a_temp
    integer                      :: i

    int = 0
    do i = 1,N_nodes
       a_temp = a+gl_nodes(i)*(b-a)
       int = int+gl_quadr_1d(nodeValues(1+(i-1)*N_nodes:i*N_nodes),a_temp,b)*gl_weights(i)
    end do
    int = int*(b-a)
  end function gl_quadr_2d

  !Takes an array of values of f(x) on the integration nodes and in interval [a,b]
  !and returns int_a^b f(x) dx
  function gl_quadr_1d(nodeValues,a,b) result(int)
    real(kind=GRID_SR),dimension(N_nodes) :: nodeValues
    real(kind=GRID_SR),intent(in)   :: a,b
    real(kind=GRID_SR)              :: int
    integer                         :: i
    int=sum(nodeValues*gl_weights)*(b-a)
  end function gl_quadr_1d
  
  ! ! subroutine gl_iniPhi(nodes,weights)
    
  ! ! end subroutine gl_iniPhi

   subroutine initGLIntegration(N)
     integer,intent(in):: N

!     N_nodes=max(N+1,2)
     N_nodes=5
     
     allocate(gl_nodes(N_nodes))
     allocate(gl_weights(N_nodes))
     !-- TODO: is it possible to generate nodes and weights when needed ? --!
     select case(N_nodes) 
     case(2)
        gl_nodes = (/.21132486540518711774542560974902127217619912436493q0,&
                     .78867513459481288225457439025097872782380087563507q0/)
        gl_weights= (/ 0.5q0,0.5q0/)
     case(3)
        gl_nodes = (/.11270166537925831148207346002176003891670782947084q0,&
                     .50000000000000000000000000000000000000000000000000q0,& 
                     .88729833462074168851792653997823996108329217052916q0/)
        gl_weights = (/.27777777777777777777777777777777777777778q0,&
                       .44444444444444444444444444444444444444444q0,&
                       .27777777777777777777777777777777777777778q0/)
     case(4)
        gl_nodes = (/0.6943184420297371238802675555359524745213731018515q-1,&
                      .33000947820757186759866712044837765639971206511455q0,&
                      .66999052179242813240133287955162234360028793488545q0,&
                      .93056815579702628761197324444640475254786268981485q0/)

        gl_weights = (/.17392742256872692868653197461099970361767q0,&
                       .32607257743127307131346802538900029638233q0,&
                       .32607257743127307131346802538900029638233q0,&
                       .17392742256872692868653197461099970361767q0/)
     case(5)
         gl_nodes = (/ .04691007703066800360118656085030351743716q0,&
                       .23076534494715845448184278964989559751635q0,&
                       .50000000000000000000000000000000000000000q0,&
                       .76923465505284154551815721035010440248365q0,&
                       .95308992296933199639881343914969648256284q0 /)

         gl_weights=(/ .11846344252809454375713202035995868132163q0,&
                       .23931433524968323402064575741781909645615q0,&
                       .28444444444444444444444444444444444444444q0,&
                       .23931433524968323402064575741781909645615q0,&
                       .11846344252809454375713202035995868132163q0/)
     end select
  !   !
  !   implicit none
  !   integer :: stat1, stat2 ,i
  !   character(len=128) ::weights_file, nodes_file


  !   write(weights_file, "(A23,I1,A4)")"../matrices/gl_weights_",N_nodes,".mtx"
  !   write(nodes_file, "(A21,I1,A4)")"../matrices/gl_nodes_",N_nodes,".mtx"

  !   print*,"'Reading weights from : "
  !   print*,weights_file

  !   print*,"Reading nodes from : "
  !   print*,nodes_file

  !   open(10,file=weights_file,status='old',iostat=stat1)
  !   open(11,file=nodes_file,status='old',iostat=stat2)

  !   allocate(weights(N_nodes))
  !   allocate(gl_nodes(N_nodes))

  !   print*,N_nodes
  !   do i=1,N_nodes
  !      read(unit=10,fmt="(E41.40)"),weights(i)
  !      read(unit=11,fmt="(E41.40)"),gl_nodes(i)
  !   end do
   
  !   close(10,iostat=stat1)
  !   close(11,iostat=stat2)

  !   return
   end subroutine initGLIntegration
  
  ! Returns 2 dimensional nodes on the lower triangel with vertices a,b

  !nodeValues ordering
  !           (b(1),b(2))
  !      | \ --|
  !      |* \ -|  
  !      |** \ |
  !(a(1),a(2))
  
  function get_gl_nodes(N)
    integer :: N
    real(kind=GRID_SR) :: get_gl_nodes(N)

     select case(N) 
     case(2)
        get_gl_nodes = (/.21132486540518711774542560974902127217619912436493q0,&
                     .78867513459481288225457439025097872782380087563507q0/)

     case(3)
        get_gl_nodes = (/.11270166537925831148207346002176003891670782947084q0,&
                     .50000000000000000000000000000000000000000000000000q0,& 
                     .88729833462074168851792653997823996108329217052916q0/)

     case(4)
        get_gl_nodes = (/0.6943184420297371238802675555359524745213731018515q-1,&
                      .33000947820757186759866712044837765639971206511455q0,&
                      .66999052179242813240133287955162234360028793488545q0,&
                      .93056815579702628761197324444640475254786268981485q0/)

     case(5)
         get_gl_nodes = (/ .04691007703066800360118656085030351743716q0,&
                       .23076534494715845448184278964989559751635q0,&
                       .50000000000000000000000000000000000000000q0,&
                       .76923465505284154551815721035010440248365q0,&
                       .95308992296933199639881343914969648256284q0 /)

     end select

 end function get_gl_nodes


  function get_gl_nodes_2d(a,b,lower_arg) result(nodes_2d)
    real(kind=GRID_SR),dimension(2,N_nodes*N_nodes)  :: nodes_2d
    real(kind=GRID_SR),dimension(N_nodes)      :: nodes_scaled
    real(kind=GRID_SR),dimension(2)      :: a,b
    integer                           :: i,count
    logical                           :: lower
    logical,optional                  :: lower_arg

    if(present(lower_arg)) then
       lower = lower_arg
    else 
       lower = .true.
    end if
    nodes_scaled=gl_nodes*(b(2)-a(2))

    count = 0
    if(lower) then
       do i = 1,N_nodes
          nodes_2d(1,count+1:count+N_nodes) = nodes_scaled*(1.0q0-gl_nodes(i)) + a(1)

          nodes_2d(2,count+1:count+N_nodes) = nodes_scaled(i) + a(2)
          count= count+N_nodes
       end do
    else
       do i = 1,N_nodes
          nodes_2d(1,count+1:count+N_nodes) = b(1)-nodes_scaled(N_nodes:1:-1)*(1.0q0-gl_nodes(i))
          nodes_2d(2,count+1:count+N_nodes) = b(2)-nodes_scaled(i)
          count= count+N_nodes
       end do
    end if
    return
  end function get_gl_nodes_2d

  function get_num_gl_nodes() 
    integer::get_num_gl_nodes
    get_num_gl_nodes=N_nodes
  end function get_num_gl_nodes

end module gl_integration
