module asserts
  implicit none
  integer,Parameter :: KIND_q=kind(1.0q0)
  integer           :: tests=0
  integer           :: failed=0

     
  interface assertEqual
     module procedure assertEqualReal_quad, assertEqualInt,assertEqualReal_double
  end interface assertEqual

  interface assertEqualArray
     module procedure assertEqualRealArray_quad
     module procedure assertEqualRealArray_double
  end interface assertEqualArray


  interface assertTolerance
     module procedure assertToleranceReal_quad
     module procedure assertToleranceReal_double
  end interface assertTolerance

contains

  subroutine endTests()
    character(100) :: resultMessage
    
    if(failed.eq.0) then
       resultMessage="All tests SUCCEEDED"
    else
       resultMessage="At least one test failed"
    endif
    write(*,*),"[",tests,"|",tests-failed,"] ",resultMessage
  end subroutine endTests

  subroutine assertEqualReal_quad(expected,actual,testMessage)
    real(kind=kind(1.0q0)),intent(in)   ::expected,actual
    logical                             ::isTrue
    character(len=*) ,intent(in)   ::testMessage

    call testHeader("Assert Equal",testMessage)

    isTrue=((expected-actual) == 0.0q0)
    write(*,*), "Expecting : ",char(9),expected
    write(*,*), "Actual: ",char(9),actual

    call testTail(isTrue)
  end subroutine assertEqualReal_quad

  subroutine assertEqualReal_double(expected,actual,testMessage)
    real(kind=kind(1.0d0)),intent(in)   ::expected,actual
    logical                             ::isTrue
    character(len=*) ,intent(in)   ::testMessage

    call testHeader("Assert Equal",testMessage)

    isTrue=((expected-actual) == 0.0d0)
    write(*,*), "Expecting : ",char(9),expected
    write(*,*), "Actual: ",char(9),actual

    call testTail(isTrue)
  end subroutine assertEqualReal_double


  subroutine assertToleranceReal_quad(expected,actual,testMessage)
    real(kind=kind(1.0q0)),intent(in)   ::expected,actual
    logical                             ::isTrue
    character(len=*) ,intent(in)   ::testMessage
    real(kind=kind(1.0q0))         ::tolerance=.2q-33

    call testHeader("Assert Tolerance",testMessage)

    isTrue=(abs(expected-actual).lt.tolerance)

    write(*,*), "Expecting : ",char(9),expected
    write(*,*), "Actual: ",char(9),actual
    write(*,*), "Difference: ",char(9),abs(actual-expected)

    call testTail(isTrue)
  end subroutine assertToleranceReal_quad

  subroutine assertToleranceReal_double(expected,actual,testMessage)
    real(kind=kind(1.0d0)),intent(in)   ::expected,actual
    logical                             ::isTrue
    character(len=*) ,intent(in)   ::testMessage
    real(kind=kind(1.0d0))         ::tolerance=.2d-33

    call testHeader("Assert Tolerance",testMessage)

    isTrue=(abs(expected-actual).lt.tolerance)

    write(*,*), "Expecting : ",char(9),expected
    write(*,*), "Actual: ",char(9),actual
    write(*,*), "Difference: ",char(9),abs(actual-expected)

    call testTail(isTrue)
  end subroutine assertToleranceReal_double


  subroutine assertEqualRealArray_quad(expected,actual,testMessage,output)
    real(kind=kind(1.0q0)),intent(in),dimension(:)   ::expected,actual
    logical,optional                             ::output
    logical                             ::isTrue,doPrint
    character(len=*) ,intent(in)   ::testMessage
    integer                         :: i
    real(kind=kind(1.0q0))         ::tolerance=.2q-29,maxError

    maxError=0.0q0
  
    if (.not.present(output)) then
       doPrint=.false.
    else 
       doPrint=output
    end if

    call testHeader("Assert Equal Array",testMessage)
!    write(*,'(A)') 'Using tolerance:'
!    write(*,'(E45.34)') tolerance
    isTrue=.true.

    if(size(expected).ne.size(actual)) then
       isTrue=.false.
    endif
    
    do i=1,size(expected)
       if(abs(actual(i)).ge.tolerance) then
          isTrue=(abs(1-expected(i)/actual(i)).lt.tolerance).and.isTrue
          maxError=max(maxError,abs(1-expected(i)/actual(i)))
       else 
          isTrue=(abs(expected(i)).le.tolerance).and.isTrue
          maxError=max(maxError,abs(expected(i))-abs(actual(i)))
       end if
    end do
    
    if (doPrint) then
       write(*,'(A)'), '[Expecting | Actual]'
       do i=1,size(expected)
          write(*,'(I4)',advance='no'), i
          write(*,'(A)',advance='no'), ':['
          write(*,'(E45.34)',advance='no'), expected(i)
          write(*,'(A)',advance='no'), "|"
          write(*,'(E45.34)',advance='no'), actual(i)
          if (abs(actual(i)) >= tolerance) then
             if (abs(1-expected(i)/actual(i)).lt.tolerance) then
                write(*,'(A)'), ']'
!                write(*,'(A)'), '1'
             else
                write(*,'(A)'), ']x'
!                write(*,'(A)'), '2'
             end if
          else
             if (abs(expected(i)) <= tolerance) then
                write(*,'(A)'), ']'
!                write(*,'(A)'), '3'
             else
                write(*,'(A)'), ']x'
!                write(*,'(A)'), '4'
             end if

          endif
       end do
    endif
    write(*,'(A,E45.38)'),' Maximal Error: ',maxError       
    call testTail(isTrue)
  end subroutine assertEqualRealArray_quad

  subroutine assertEqualRealArray_double(expected,actual,testMessage,output)
    real(kind=kind(1.0d0)),intent(in),dimension(:)   ::expected,actual
    logical,optional                             ::output
    logical                             ::isTrue,doPrint
    character(len=*) ,intent(in)   ::testMessage
    integer                         :: i
    real(kind=kind(1.0d0))         ::tolerance=.2d-12,maxError

    maxError=0.0q0
  
    if (.not.present(output)) then
       doPrint=.false.
    else 
       doPrint=output
    end if

    call testHeader("Assert Equal Array",testMessage)
!    write(*,'(A)') 'Using tolerance:'
!    write(*,'(E45.34)') tolerance
    isTrue=.true.

    if(size(expected).ne.size(actual)) then
       isTrue=.false.
    endif
    
    do i=1,size(expected)
       if(abs(actual(i)).ge.tolerance) then
          isTrue=(abs(1-expected(i)/actual(i)).lt.tolerance).and.isTrue
          maxError=max(maxError,abs(1-expected(i)/actual(i)))
       else 
          isTrue=(abs(expected(i)).le.tolerance).and.isTrue
          maxError=max(maxError,abs(expected(i))-abs(actual(i)))
       end if
    end do
    
    if (doPrint) then
       write(*,'(A)'), '[Expecting | Actual]'
       do i=1,size(expected)
          write(*,'(I4)',advance='no'), i
          write(*,'(A)',advance='no'), ':['
          write(*,'(E24.34)',advance='no'), expected(i)
          write(*,'(A)',advance='no'), "|"
          write(*,'(E24.34)',advance='no'), actual(i)
          if (abs(actual(i)) >= tolerance) then
             if (abs(1-expected(i)/actual(i)).lt.tolerance) then
                write(*,'(A)'), ']'
!                write(*,'(A)'), '1'
             else
                write(*,'(A)'), ']x'
!                write(*,'(A)'), '2'
             end if
          else
             if (abs(expected(i)) <= tolerance) then
                write(*,'(A)'), ']'
!                write(*,'(A)'), '3'
             else
                write(*,'(A)'), ']x'
!                write(*,'(A)'), '4'
             end if

          endif
       end do
    endif
    write(*,'(A,E24.17)'),' Maximal Error: ',maxError       
    call testTail(isTrue)
  end subroutine assertEqualRealArray_double


  ! subroutine assertEqualRealArray_double(expected,actual,testMessage,output)
  !   real(kind=kind(1.0d0)),intent(in),dimension(:)   ::expected,actual
  !   logical,optional                             ::output
  !   logical                             ::isTrue,doPrint
  !   character(len=*) ,intent(in)   ::testMessage
  !   integer                         :: i
  !   real(kind=kind(1.0d0))         ::tolerance=.3d-15, maxError=.0d0


  !   if (.not.present(output)) then
  !      doPrint=.false.
  !   else 
  !      doPrint=.true.
  !   end if

  !   call testHeader("Assert Equal Array",testMessage)
  !   isTrue=.true.

  !   if(size(expected).ne.size(actual)) then
  !      print *, "Sizes differ ",size(expected)," to ",size(actual)
  !      isTrue=.false.
  !      return
  !   endif
    
  !   do i=1,size(expected)
  !      isTrue=(abs(1-expected(i)/actual(i)).lt.tolerance).and.isTrue
  !      maxError=max(maxError,abs(1-expected(i)/actual(i)))
  !   end do
    
  !   if (doPrint) then
  !      write(*,'(A)'), '[Expecting | Actual]'
  !      do i=1,size(expected)
  !         write(*,'(I4)',advance='no'), i
  !         write(*,'(A)',advance='no'), ':['
  !         write(*,'(E20.15)',advance='no'), expected(i)
  !         write(*,'(A)',advance='no'), "|"
  !         write(*,'(E20.15)',advance='no'), actual(i)
  !         if (abs(1-expected(i)/actual(i)).lt.tolerance) then
  !         write(*,'(A)'), ']'
  !         else
  !         write(*,'(A)'), ']x'
  !         end if
  !      end do
  !   endif
  !   write(*,'(A,E20.15)'),'Maximal Error: ',maxError
       
  !   call testTail(isTrue)
  ! end subroutine assertEqualRealArray_double


  subroutine assertEqualInt(Expected,actual,testMessage)
    implicit none
    integer,intent(in)        ::expected,actual
    logical                   ::isTrue
    character(len=*) ,intent(in)   ::testMessage

    call testHeader("Assert Equal",testMessage)

    isTrue=((expected-actual) == 0)

    write(*,*), "Expecting : ",char(9),expected
    write(*,*), "Actual: ",char(9),actual

    call testTail(isTrue)
  end subroutine assertEqualInt

  subroutine testHeader(testName,testMessage)
    implicit none
    character(len=*) ::testName,testMessage

    tests = tests+1
    write(*,*),"_________________________"
    write(*,*),testName
    write(*,*),testMessage
  end subroutine testHeader

  subroutine testTail(isTrue)
    implicit none
    logical, intent(in) ::isTrue

    if (.not.isTrue) then
       write(*,*), "Test FAILED"
       failed = failed+1
    else
       write(*,*), "Test SUCCEEDED"
    endif
    write(*,*),"_________________________"
    write(*,*),""
  end subroutine testTail

end module asserts
