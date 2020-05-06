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
  public DG,TROUBLED,NEIGHBOUR_TROUBLED,COAST,DRY

  enum,bind( C )
     enumerator :: DG = 0                 ,&
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



  subroutine updateCellStatus(data)
    type(num_cell_data_pers),intent(inout) :: data

    !---- Coast stays coast ----!
    if(data%troubled == COAST) then
       return
    end if
    !---------------------------!

    if (data%troubled == DG) then
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


  subroutine getObservables(Q,observables)
    type(t_state)      :: Q
    real(kind=GRID_SR) :: observables(_DMP_NUM_OBSERVABLES)
#if defined(_SWE_DG_LIMITER_UNLIMITED)
    return
#elif defined(_SWE_DG_LIMITER_HEIGHT)
    observables(1) = Q%h
#elif defined(_SWE_DG_LIMITER_ALL)
    observables(1) = Q%h
    observables(2) = Q%p(1)
    observables(3) = Q%p(2)
#endif
  end subroutine getObservables

  subroutine getObservableLimits(Q,minVals,maxVals)
    type(t_state)     , dimension(_SWE_DG_DOFS),intent(in)          :: Q
    real(kind=GRID_SR), dimension(_DMP_NUM_OBSERVABLES),intent(out) :: minVals
    real(kind=GRID_SR), dimension(_DMP_NUM_OBSERVABLES),intent(out) :: maxVals
    real(kind=GRID_SR), dimension(_DMP_NUM_OBSERVABLES)             :: observables
    integer                                                         :: i,j

    minVals(:) =  Huge(1.0_GRID_SR)
    maxVals(:) = -Huge(1.0_GRID_SR)
    do i = 1,_SWE_DG_DOFS
       observables = 0
       call getObservables(Q(i),observables)
       do j = 1,_DMP_NUM_OBSERVABLES
          minVals(j) = min(minVals(j),observables(j))
          maxVals(j) = max(minVals(j),observables(j))
       end do
    end do
  end subroutine getObservableLimits

  function checkDMP(Q,minVals,maxVals,update1,update2,update3) result(troubled)
    type(t_state)         , intent(in) :: Q(_SWE_DG_DOFS)
    type(num_cell_update) , intent(in) :: update1
    type(num_cell_update) , intent(in) :: update2
    type(num_cell_update) , intent(in) :: update3
    real(kind=GRID_SR)    , intent(in) :: maxVals(_DMP_NUM_OBSERVABLES)
    real(kind=GRID_SR)    , intent(in) :: minVals(_DMP_NUM_OBSERVABLES)
    real(kind=GRID_SR)                 :: maxNeighbour(_DMP_NUM_OBSERVABLES)
    real(kind=GRID_SR)                 :: minNeighbour(_DMP_NUM_OBSERVABLES)
    real(kind=GRID_SR)                 :: delta(_DMP_NUM_OBSERVABLES)
    real(kind=GRID_SR)                 :: observables(_DMP_NUM_OBSERVABLES)
    integer                            :: i,j
    logical                            :: troubled

    troubled = .FALSE.
    
    do i = 1,_DMP_NUM_OBSERVABLES
       maxNeighbour(i) = maxVals(i)
       maxNeighbour(i) = max(maxNeighbour(i),update1%maxObservables(i))
       maxNeighbour(i) = max(maxNeighbour(i),update2%maxObservables(i))
       maxNeighbour(i) = max(maxNeighbour(i),update3%maxObservables(i))

       minNeighbour(i) = minVals(i)
       minNeighbour(i) = min(minNeighbour(i),update1%minObservables(i))
       minNeighbour(i) = min(minNeighbour(i),update2%minObservables(i))
       minNeighbour(i) = min(minNeighbour(i),update3%minObservables(i))
    end do

    do i = 1,_DMP_NUM_OBSERVABLES
       delta(i) = max(0.1_GRID_SR,maxNeighbour(i)-minNeighbour(i))
    end do

    do i = 1,_SWE_DG_DOFS
       call getObservables(Q(i),observables)
       do j = 1,_DMP_NUM_OBSERVABLES
          troubled = troubled .or. observables(i) > maxNeighbour(i)+delta(i) &
                              .or. observables(i) < minNeighbour(i)-delta(i)
       end do
    end do
    if(troubled) then
       print*,troubled
    endif
  end function checkDMP


END MODULE SWE_DG_Limiter
