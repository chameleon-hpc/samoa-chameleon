!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!*                        (www.jpmoreau.fr)            *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************
MODULE LU
use gl_integration

CONTAINS

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
 subroutine ludcmp(a,n,indx,d,code)
 integer ::  n,ll,i,imax,j,k
 integer,parameter :: nmax=100
 real(kind=GRID_SR),parameter :: tiny=1.5e-33_GRID_SR
 real(kind=GRID_SR) :: amax, dum, sum, a(n,n), vv(nmax)
 integer code, d, indx(n)

 d=1; code=0

 do i=1,n
   amax=0.0_GRID_SR
   do j=1,n
     if (abs(a(i,j)).gt.amax) amax=abs(a(i,j))
   end do ! j loop
   if(amax.lt.tiny) then
     code = 1
     return
   end if
   vv(i) = 1.0_GRID_SR / amax
 end do ! i loop

 do j=1,n
   do i=1,j-1
     sum = a(i,j)
     do k=1,i-1
       sum = sum - a(i,k)*a(k,j) 
     end do ! k loop
     a(i,j) = sum
   end do ! i loop
   amax = 0.0_GRID_SR
   do i=j,n
     sum = a(i,j)
     do k=1,j-1
       sum = sum - a(i,k)*a(k,j) 
     end do ! k loop
     a(i,j) = sum
     dum = vv(i)*abs(sum)
     if(dum.gt.amax) then
       imax = i
       amax = dum
     end if
   end do ! i loop  

   
   if(j.ne.imax) then
     do k=1,n
       dum = a(imax,k)
       a(imax,k) = a(j,k)
       a(j,k) = dum
     end do ! k loop
     d = -d
     vv(imax) = vv(j)
   end if

   indx(j) = imax
   if(abs(a(j,j)) < tiny) a(j,j) = tiny

   if(j.ne.n) then
     dum = 1.0_GRID_SR / a(j,j)
     do i=j+1,n
       a(i,j) = a(i,j)*dum
     end do ! i loop
   end if 
 end do ! j loop

 return
 end subroutine ludcmp


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
 Subroutine LUBKSB(A,N,INDX,B)
 INTEGER :: N,INDX(N),II,I,J,LL
 real(kind=GRID_SR) ::  SUM, A(N,N),B(N)

 II = 0

 DO I=1,N
   LL = INDX(I)
   SUM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   ELSE IF(SUM.NE.0.0_GRID_SR) THEN
     II = I
   END IF
   B(I) = SUM
 END DO ! i loop

 DO I=N,1,-1
   SUM = B(I)
   IF(I < N) THEN
     DO J=I+1,N
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   END IF
   B(I) = SUM / A(I,I)
 END DO ! i loop

 RETURN
 END subroutine LUBKSB

END MODULE LU

! end of file lu.f90
