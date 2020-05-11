#include "Compilation_control.f90"

#if defined(_SWE_DG)
#define _REF_TRIANGLE_SIZE_INV  (2.0q0 * real(_SWE_PATCH_ORDER_SQUARE,kind=kind(1.0q0)))
MODULE SWE_dg_matrices
	implicit none
	
#       if defined(_SINGLE_PRECISION)
            integer, PARAMETER :: GRID_SR = kind(1.0e0)
#       elif defined(_DOUBLE_PRECISION)
            integer, PARAMETER :: GRID_SR = kind(1.0d0)
#       elif defined(_QUAD_PRECISION)
            integer, PARAMETER :: GRID_SR = kind(1.0q0)
#       else
#           error "No floating point precision is chosen!"
#       endif

! #if defined (_SWE_DG_EXTRAPOLATE_BND)
! #define _SWE_DG_BND_INTEGRATION_POINTS 0
! #else
! #define _SWE_DG_BND_INTEGRATION_POINTS _SWE_DG_ORDER+1
! #endif
! #define 

!Time matrices
#include "dg_matrices/t_m_1_3.incl"
#include "dg_matrices/t_a_3.incl"
#include "dg_matrices/t_k_t_10_3.incl"
#include "dg_matrices/t_k_t_11_inv_3.incl"
    	    
!Space matrices
#include "dg_matrices/s_m_3.incl"
#include "dg_matrices/s_m_inv_3.incl"
#include "dg_matrices/s_k_x_3.incl"
#include "dg_matrices/s_k_y_3.incl"

!Boundary matrices
#include "dg_matrices/s_b_1_3.incl"
#include "dg_matrices/s_b_2_3.incl"
#include "dg_matrices/s_b_3_3.incl"

#include "dg_matrices/s_b_1_l_3.incl"
#include "dg_matrices/s_b_2_m_3.incl"
#include "dg_matrices/s_b_3_r_3.incl"

!FV projection
#include "dg_matrices/phi_l_3.incl"
#include "dg_matrices/phi_m_3.incl"
#include "dg_matrices/phi_r_3.incl"
#include "dg_matrices/phi_3.incl"
#include "dg_matrices/mue_inv_3.incl"

!Coordinates	    
#include "dg_matrices/nodes_3.incl"
#include "dg_matrices/mirrored_coords_3.incl"

!Basis derivatives
#include "dg_matrices/basis_der_x_3.incl"
#include "dg_matrices/basis_der_y_3.incl"
	    
!Refinement matrices
#include "dg_matrices/ref1_3.incl"
#include "dg_matrices/ref2_3.incl"
#include "dg_matrices/coarsen1_3.incl"
#include "dg_matrices/coarsen2_3.incl"

real(kind=GRID_SR), Parameter :: s_m_inv_s_k_x_t   (_SWE_DG_DOFS,_SWE_DG_DOFS)     = matmul(s_m_inv,transpose(s_k_x))
real(kind=GRID_SR), Parameter :: s_m_inv_s_k_y_t   (_SWE_DG_DOFS,_SWE_DG_DOFS)     = matmul(s_m_inv,transpose(s_k_y))
real(kind=GRID_SR), Parameter :: t_k_t_11_inv_t_m_1(_SWE_DG_ORDER,_SWE_DG_ORDER+1) = matmul(t_k_t_11_inv,t_m_1)
real(kind=GRID_SR), Parameter :: s_k_x_s_b_3_s_b_2 (_SWE_DG_DOFS,_SWE_DG_DOFS)     = s_k_x + s_b_3 - s_b_2
real(kind=GRID_SR), Parameter :: s_k_y_s_b_1_s_b_2 (_SWE_DG_DOFS,_SWE_DG_DOFS)     = s_k_y + s_b_1 - s_b_2
real(kind=GRID_SR), Parameter :: t_k_t_11_inv_x_t_k_t_10(_SWE_DG_ORDER,1)          = matmul(t_k_t_11_inv,t_k_t_10)

real(kind=GRID_SR),Parameter :: s_m_inv_ref1(_SWE_DG_DOFS,_SWE_DG_DOFS)  = matmul(s_m_inv,ref1)
real(kind=GRID_SR),Parameter :: s_m_inv_ref2(_SWE_DG_DOFS,_SWE_DG_DOFS)  = matmul(s_m_inv,ref2)
real(kind=GRID_SR),Parameter :: s_m_inv_coarsen1(_SWE_DG_DOFS,_SWE_DG_DOFS) =&
			     matmul(s_m_inv,coarsen1)
real(kind=GRID_SR),Parameter :: s_m_inv_coarsen2(_SWE_DG_DOFS,_SWE_DG_DOFS) =&
			     matmul(s_m_inv,coarsen2)


contains 

subroutine lusolve(mat,n,pivot,b)
 integer :: n,ii,i,j,ll
 integer, intent(in) :: pivot(n)

 real(kind=GRID_SR) ::  sum
 real(kind=GRID_SR),intent(inout) ::  b(n)
 real(kind=GRID_SR),intent(in) ::  mat(n,n)

 ii = 0

 do i=1,n
   ll = pivot(i)
   sum = b(ll)
   b(ll) = b(i)
   if(ii.ne.0) then
     do j=ii,i-1
       sum = sum - mat(i,j)*b(j)
     end do ! j loop
   else if(sum.ne.0.0_GRID_SR) then
     ii = i
   end if
   b(i) = sum
 end do ! i loop

 do i=n,1,-1
   sum = b(i)
   if(i < n) then
     do j=i+1,n
       sum = sum - mat(i,j)*b(j)
     end do ! j loop
   end if
   b(i) = sum / mat(i,i)
 end do ! i loop

 return
end subroutine lusolve




END MODULE SWE_dg_matrices
#endif

