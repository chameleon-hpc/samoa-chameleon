! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

MODULE SWE_DG_Limiter
  use SWE_data_types
  use Samoa_swe
  
  !*************************!
  !States for Troubled cells!
  !*************************!
  public DG,TROUBLED,NEIGHBOUR_TROUBLED,NEIGHBOUR_WAS_TROUBLED,COAST,DRY
  public node_fist_touch, node_merge

  enum,bind( C )
     enumerator :: NEIGHBOUR_WAS_TROUBLED = -2 ,&
                   DG = 0                 ,&
                   WET_DRY_INTERFACE  = 1 ,&
                   NEIGHBOUR_TROUBLED = 2 ,&
                   TROUBLED           = 3 ,&
                   COAST              = 4 ,&
                   DRY                = 5 ,&
                   PREDICTOR_DIVERGED = 6
  end enum

contains

  function isDG(troubled)
    integer :: troubled
    logical :: isDG
    isDG = troubled.le.0
  end function isDG
  
  function isFV(troubled)
    integer :: troubled
    logical :: isFV
    isFV = troubled.ge.1
  end function isFV

  function isCoast(troubled)
    integer :: troubled
    logical :: isCoast
    isCoast = troubled.eq.COAST
  end function isCoast

  function isDry(troubled)
    integer :: troubled
    logical :: isDry
    isDry = troubled.eq.DRY
  end function isDry

  function checkIfCellIsDry(H) result(isDry)
    real(kind=GRID_SR),intent(in) :: H(:)
    logical                       :: isDry
    if (all(H < cfg%dry_tolerance)) then
       isDry = .True.
    else
       isDry = .False.
    end if
  end function checkIfCellIsDry

  subroutine updateCellStatusPre(data,update1,update2,update3)
    type(num_cell_data_pers),intent(inout) :: data
    type(num_cell_update) , intent(in) :: update1
    type(num_cell_update) , intent(in) :: update2
    type(num_cell_update) , intent(in) :: update3

    !----If any neighbour is troubled use FV scheme ---!
    if(isDG(data%troubled)) then
       if(neighbourTroubled(update1,update2,update3)) then
           data%troubled = NEIGHBOUR_TROUBLED
       end if
    end if
    !--------------------------------------------------!
    
  end subroutine updateCellStatusPre


  subroutine updateCellStatus(data)
    type(num_cell_data_pers),intent(inout) :: data

    !---- Coast stays coast ----!
    if(data%troubled == COAST) then
       return
    end if

    ! if(data%troubled == NEIGHBOUR_TROUBLED) then
    !    return
    ! end if
    !---------------------------!

    if (data%troubled .le. DG) then
       data%troubled = DG
    else
       data%troubled = -data%troubled
    end if
    
    if(isWetDryInterface(data%Q%H))then
       data%troubled = WET_DRY_INTERFACE
    end if
    
    if(checkIfCellIsDry(data%Q%H)) then
       data%troubled = DRY
    end if
    
  end subroutine updateCellStatus

  function isWetDryInterface(H)
    real(kind=GRID_SR),intent(in) :: H(:)
    logical                       :: isWetDryInterface
    if (.not.all(H > cfg%dry_dg_guard)) then
       isWetDryInterface = .True.
    else
       isWetDryInterface = .False.
    end if
  end function isWetDryInterface

  function allNeighboursDry(update1,update2,update3) result(neighbours_dry)
    type(num_cell_update) , intent(in) :: update1
    type(num_cell_update) , intent(in) :: update2
    type(num_cell_update) , intent(in) :: update3
    logical                            :: neighbours_dry
    neighbours_dry = (update1%troubled.eq.DRY) .and. &
                     (update2%troubled.eq.DRY) .and. &
                     (update3%troubled.eq.DRY)
  end function allNeighboursDry


  function neighbourTroubled(update1,update2,update3) result(neighbour_troubled)
    type(num_cell_update) , intent(in) :: update1
    type(num_cell_update) , intent(in) :: update2
    type(num_cell_update) , intent(in) :: update3
    logical                            :: neighbour_troubled
    neighbour_troubled = (update1%troubled.ge.1) .or. &
                         (update2%troubled.ge.1) .or. &
                         (update3%troubled.ge.1)
  end function neighbourTroubled


  subroutine getObservables(H,HU,HV,B,observables)
    real(kind=GRID_SR), intent(in)	:: H,HU,HV,B
    real(kind=GRID_SR), intent(out) :: observables(_DMP_NUM_OBSERVABLES)
#if defined(_SWE_DG_LIMITER_UNLIMITED)
    return
#elif defined(_SWE_DG_LIMITER_HEIGHT)
    observables(1) = H - B
#elif defined(_SWE_DG_LIMITER_ALL)
    observables(1) = H - B
    observables(2) = HU
    observables(3) = HV
#endif
  end subroutine getObservables

  subroutine getObservableLimits(data,minVals,maxVals)
    class (num_cell_data_pers),intent(in) :: data
    real(kind=GRID_SR), dimension(_DMP_NUM_OBSERVABLES),intent(out) :: minVals
    real(kind=GRID_SR), dimension(_DMP_NUM_OBSERVABLES),intent(out) :: maxVals
    real(kind=GRID_SR), dimension(_DMP_NUM_OBSERVABLES)             :: observables
    integer                                                         :: i,j

    minVals(:) =  Huge(1.0_GRID_SR)
    maxVals(:) = -Huge(1.0_GRID_SR)

    do i = 1,_SWE_PATCH_ORDER_SQUARE
       call getObservables(data%H(i),data%HU(i),data%HV(i),data%B(i),observables)
       do j = 1,_DMP_NUM_OBSERVABLES
          minVals(j) = min(minVals(j),observables(j))
          maxVals(j) = max(maxVals(j),observables(j))
       end do
    end do
  end subroutine getObservableLimits

  function checkDMP(data,minVals,maxVals,node1,node2,node3) result(troubled)
    class(num_cell_data_pers),intent(in) :: data
    ! type(num_cell_update) , intent(in) :: update1
    ! type(num_cell_update) , intent(in) :: update2
    ! type(num_cell_update) , intent(in) :: update3
    type(num_node_data_pers) , intent(in) :: node1
    type(num_node_data_pers) , intent(in) :: node2
    type(num_node_data_pers) , intent(in) :: node3
    real(kind=GRID_SR)    , intent(in)    :: maxVals(_DMP_NUM_OBSERVABLES)
    real(kind=GRID_SR)    , intent(in)    :: minVals(_DMP_NUM_OBSERVABLES)
    real(kind=GRID_SR)                    :: maxNeighbour(_DMP_NUM_OBSERVABLES)
    real(kind=GRID_SR)                    :: minNeighbour(_DMP_NUM_OBSERVABLES)
    real(kind=GRID_SR)                    :: delta(_DMP_NUM_OBSERVABLES)
    real(kind=GRID_SR)                    :: observables(_DMP_NUM_OBSERVABLES)
    integer                               :: i,j
    logical                               :: troubled

    troubled = .FALSE.
    
    do i = 1,_DMP_NUM_OBSERVABLES
       maxNeighbour(i) = maxVals(i)
       maxNeighbour(i) = max(maxNeighbour(i),node1%maxObservables(i))
       maxNeighbour(i) = max(maxNeighbour(i),node2%maxObservables(i))
       maxNeighbour(i) = max(maxNeighbour(i),node3%maxObservables(i))

       minNeighbour(i) = minVals(i)
       minNeighbour(i) = min(minNeighbour(i),node1%minObservables(i))
       minNeighbour(i) = min(minNeighbour(i),node2%minObservables(i))
       minNeighbour(i) = min(minNeighbour(i),node3%minObservables(i))
    end do

    do i = 1,_DMP_NUM_OBSERVABLES
       delta(i) = max(0.01_GRID_SR,maxNeighbour(i) - minNeighbour(i)) * cfg%limiter_buffer
    end do

    do i = 1,_SWE_PATCH_ORDER_SQUARE
       call getObservables(data%H(i),data%HU(i),data%HV(i),data%B(i),observables)
       do j = 1,_DMP_NUM_OBSERVABLES
          troubled = troubled .or. observables(j) > maxNeighbour(j)+delta(j) &
                              .or. observables(j) < minNeighbour(j)-delta(j)
       end do
    end do

  end function checkDMP

  function get_error_estimate(Q) result(error)
    type(t_state)        , DIMENSION(_SWE_DG_DOFS) :: Q
    real(kind=GRID_SR) :: error
    real(kind=GRID_SR) :: min_h
    real(kind=GRID_SR) :: max_h
    integer :: i
    max_h= -huge(1.0_GRID_SR)
    min_h=  huge(1.0_GRID_SR)
    do i=1,_SWE_DG_DOFS
       min_h = min(Q(i)%h+Q(i)%b,min_h)
       max_h = max(Q(i)%h+Q(i)%b,max_h)
    end do
    error = abs(max_h-min_h)
  end function get_error_estimate

  function plotFV(troubled)
    integer :: troubled
    logical :: plotFV
    
    !plotFV = isCoast(troubled) .or. (troubled .eq. -NEIGHBOUR_TROUBLED)
    !plotFV = isCoast(troubled)
    plotFV = isFV(troubled)
  end function plotFV

  subroutine node_first_touch(node)
    type(t_node_data), intent(inout)                   :: node
    
  end subroutine node_first_touch

  
  subroutine node_merge(local_node, neighbor_node)
    type(t_node_data), intent(inout)		:: local_node
    type(t_node_data), intent(in)		    :: neighbor_node

  end subroutine node_merge

END MODULE SWE_DG_Limiter
