! Simple self-contained solver for system of linear equations: 
! Source: https://github.com/cfinch/Shocksolution_Examples/tree/master/FORTRAN/LinearEquationSolver

!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************

    !  ***************************************************************
    !  * Given an N x N matrix A, this routine replaces it by the LU *
    !  * decomposition of a rowwise permutation of itself. A and N   *
    !  * are input. INDX is an output vector which records the row   *
    !  * permutation effected by the partial pivoting; D is output   *
    !  * as -1 or 1, depending on whether the number of row inter-   *
    !  * changes was even or odd, respectively. This routine is used *
    !  * in combination with LUBKSB to solve linear equations or to  *
    !  * invert a matrix. Return code is 1, if matrix is singular.   *
    !  ***************************************************************
     Subroutine LUDCMP(A,N,INDX,D,CODE)
     IMPLICIT NONE
     !integer, parameter :: nmax = 100
     real, parameter :: tiny = 1.5D-16

     real*8, intent(inout):: A(6,6)
     integer, intent(in) :: N
     integer, intent(out) :: D, CODE
     integer, intent(out) :: INDX(6)
     !f2py depend(N) A, indx

     REAL*8  :: AMAX, DUM, SUMM, VV(10)
     INTEGER :: i, j, k, imax

     D=1; CODE=0

     DO I=1,N
       AMAX=0.d0
       DO J=1,N
         IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
       END DO ! j loop
       IF(AMAX.LT.TINY) THEN
         CODE = 1
         !RETURN
       END IF
       VV(I) = 1.d0 / AMAX
     END DO ! i loop

     DO J=1,N
       DO I=1,J-1
         SUMM = A(I,J)
         DO K=1,I-1
           SUMM = SUMM - A(I,K)*A(K,J) 
         END DO ! k loop
         A(I,J) = SUMM
       END DO ! i loop
       AMAX = 0.d0
       DO I=J,N
         SUMM = A(I,J)
         DO K=1,J-1
           SUMM = SUMM - A(I,K)*A(K,J) 
         END DO ! k loop
         A(I,J) = SUMM
         DUM = VV(I)*DABS(SUMM)
         IF(DUM.GE.AMAX) THEN
           IMAX = I
           AMAX = DUM
         END IF
       END DO ! i loop  
       
       IF(J.NE.IMAX) THEN
         DO K=1,N
           DUM = A(IMAX,K)
           A(IMAX,K) = A(J,K)
           A(J,K) = DUM
         END DO ! k loop
         D = -D
         VV(IMAX) = VV(J)
       END IF

       INDX(J) = IMAX
       IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

       IF(J.NE.N) THEN
         DUM = 1.d0 / A(J,J)
         DO I=J+1,N
           A(I,J) = A(I,J)*DUM
         END DO ! i loop
       END IF 
     END DO ! j loop

     RETURN
 END subroutine LUDCMP


!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
 Subroutine LUBKSB(A, N, INDX, B)
 integer, intent(in) :: N 
 real*8, intent(in) :: A(6,6)
 integer, intent(in) :: INDX(6)
 real*8, intent(inout) :: B(6)
 !f2py depend(N) A, INDX, B

 REAL*8  SUMM
 integer II, LL, I, J

 II = 0

 DO I=1,N
   LL = INDX(I)
   SUMM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUMM = SUMM - A(I,J)*B(J)
     END DO ! j loop
   ELSE IF(SUMM.NE.0.d0) THEN
     II = I
   END IF
   B(I) = SUMM
 END DO ! i loop

 
 DO I=N,1,-1
   SUMM = B(I)
   IF(I < N) THEN
     DO J=I+1,N
       SUMM = SUMM - A(I,J)*B(J)
     END DO ! j loop
   END IF
   B(I) = SUMM / A(I,I)
 END DO ! i loop

 RETURN
 END subroutine LUBKSB
 
 
 
!  LU() = LUDCMP() + LUBKSB()
subroutine LU(A, N, B)
    integer, intent(in) :: N 
    real*8, intent(inout) :: A(6,6)
    real*8, intent(inout) :: B(6)
    !f2py depend(N) A, INDX, B

    integer II, LL, I, J, K
    !integer, parameter :: nmax = 100
    real, parameter :: tiny = 1.5D-16
    integer :: INDX(6)
    REAL*8  :: AMAX, DUM, SUMM, VV(10)
    INTEGER :: imax
    
    !LUDCMP
!     D=1; CODE=0

    DO I=1,N
        AMAX=0.d0
        DO J=1,N
            IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
        END DO ! j loop
!         IF(AAX.LT.TINY) THEN
!             CODE = 1
!             !RETURN
!         END IF
        VV(I) = 1.d0 / AMAX
    END DO ! i loop

    DO J=1,N
        DO I=1,J-1
            SUMM = A(I,J)
            DO K=1,I-1
                SUMM = SUMM - A(I,K)*A(K,J) 
            END DO ! k loop
            A(I,J) = SUMM
        END DO ! i loop
        AMAX = 0.d0
        DO I=J,N
            SUMM = A(I,J)
            DO K=1,J-1
                SUMM = SUMM - A(I,K)*A(K,J) 
            END DO ! k loop
            A(I,J) = SUMM
            DUM = VV(I)*DABS(SUMM)
            IF(DUM.GE.AMAX) THEN
                IMAX = I
                AMAX = DUM
            END IF
        END DO ! i loop  
    
        IF(J.NE.IMAX) THEN
            DO K=1,N
                DUM = A(IMAX,K)
                A(IMAX,K) = A(J,K)
                A(J,K) = DUM
            END DO ! k loop
!             D = -D
            VV(IMAX) = VV(J)
        END IF

        INDX(J) = IMAX
        IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

        IF(J.NE.N) THEN
            DUM = 1.d0 / A(J,J)
            DO I=J+1,N
                A(I,J) = A(I,J)*DUM
            END DO ! i loop
        END IF 
    END DO ! j loop
    
    
    !LUBKSB
    II = 0

    DO I=1,N
        LL = INDX(I)
        SUMM = B(LL)
        B(LL) = B(I)
        IF(II.NE.0) THEN
            DO J=II,I-1
                SUMM = SUMM - A(I,J)*B(J)
            END DO ! j loop
        ELSE IF(SUMM.NE.0.d0) THEN
            II = I
        END IF
        B(I) = SUMM
    END DO ! i loop
    
    DO I=N,1,-1
        SUMM = B(I)
        IF(I < N) THEN
            DO J=I+1,N
            SUMM = SUMM - A(I,J)*B(J)
            END DO ! j loop
        END IF
        B(I) = SUMM / A(I,I)
    END DO ! i loop
    
end subroutine


subroutine LU_scalars(  A11, A12, A13, A14, A15, A16, &
                        A21, A22, A23, A24, A25, A26, &
                        A31, A32, A33, A34, A35, A36, &
                        A41, A42, A43, A44, A45, A46, &
                        A51, A52, A53, A54, A55, A56, &
                        A61, A62, A63, A64, A65, A66, &
                        N, B1, B2, B3, B4, B5, B6)
                   
    integer, intent(in) :: N 
    real*8, intent(inout) ::    A11, A12, A13, A14, A15, A16, &
                                A21, A22, A23, A24, A25, A26, &
                                A31, A32, A33, A34, A35, A36, &
                                A41, A42, A43, A44, A45, A46, &
                                A51, A52, A53, A54, A55, A56, &
                                A61, A62, A63, A64, A65, A66
                                
    real*8, intent(inout) ::    B1, B2, B3, B4, B5, B6
    !f2py depend(N) A, INDX, B

    integer II, LL, I, J
    !integer, parameter :: nmax = 100
    real, parameter :: tiny = 1.5D-16
    integer :: INDX1, INDX2, INDX3, INDX4, INDX5, INDX6
    REAL*8  :: AMAX, DUM, SUMM, VV1, VV2, VV3, VV4, VV5, VV6 
    INTEGER :: imax
    

    ! This is a HUGE code where all matrices were replaced with scalar variables and all loops were unrolled.
    ! The code was generated with generate_lu_scalars.cpp 
    
#   include "lu_scalars.f90"
    
end subroutine
