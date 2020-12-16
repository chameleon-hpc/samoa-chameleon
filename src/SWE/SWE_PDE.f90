! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2017 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE
! Author: Leonhard Rannabauer rannabau (at) in.tum.de

#include "Compilation_control.f90"

#if defined(_SWE_DG)

MODULE SWE_PDE
  use SWE_data_types 

  implicit none

contains


  function flux(q)
    real(kind=GRID_SR)             :: flux(_SWE_DG_DOFS,2,3)
    real(kind=GRID_SR), intent(in) :: q(_SWE_DG_DOFS,3)
    integer                        :: k
    real(kind=GRID_SR)             :: u(_SWE_DG_DOFS), v(_SWE_DG_DOFS)

    !$omp simd
    do k=1,_SWE_DG_DOFS
       u(k) = q(k,2)/q(k,1)
       v(k) = q(k,3)/q(k,1)
    end do

    !$omp simd
    do k=1,_SWE_DG_DOFS       
       flux(k,1,1) = q(k,2)
       flux(k,1,2) = u(k)*q(k,2) + 0.5_GRID_SR * g * q(k,1)**2
       flux(k,1,3) = u(k)*q(k,3)
       
       flux(k,2,1) = q(k,3)
       flux(k,2,2) = v(k)*q(k,2)
       flux(k,2,3) = v(k)*q(k,3) + 0.5_GRID_SR * g * q(k,1)**2
    end do
       
  end function flux

    function flux_(q)
    real(kind=GRID_SR)             :: flux_(_SWE_DG_ORDER+1,2,3)
    real(kind=GRID_SR), intent(in) :: q(_SWE_DG_ORDER+1,3)
    integer                        :: k
    real(kind=GRID_SR)             :: u(_SWE_DG_ORDER+1), v(_SWE_DG_ORDER+1)

    !$omp simd
    do k=1,_SWE_DG_ORDER+1
       u(k) = q(k,2)/q(k,1)
       v(k) = q(k,3)/q(k,1)
    end do

    !$omp simd
    do k=1,_SWE_DG_ORDER+1       
       flux_(k,1,1) = q(k,2)
       flux_(k,1,2) = u(k)*q(k,2)
       flux_(k,1,3) = u(k)*q(k,3)
       
       flux_(k,2,1) = q(k,3)
       flux_(k,2,2) = v(k)*q(k,2)
       flux_(k,2,3) = v(k)*q(k,3)
    end do
       
  end function flux_


  function flux_no_grav(q)
    real(kind=GRID_SR)             ::flux_no_grav(_SWE_DG_DOFS,2,3)
    real(kind=GRID_SR) ,intent(in) ::q(_SWE_DG_DOFS,3)
    integer :: k
    real(kind=GRID_SR)             :: u(_SWE_DG_DOFS), v(_SWE_DG_DOFS)

    !$omp simd
    do k=1,_SWE_DG_DOFS
       u(k) = q(k,2)/q(k,1)
       v(k) = q(k,3)/q(k,1)
    end do

    !$omp simd
    do k=1,_SWE_DG_DOFS       
       flux_no_grav(k,1,1) = q(k,2)
       flux_no_grav(k,1,2) = u(k)*q(k,2)
       flux_no_grav(k,1,3) = u(k)*q(k,3)
       
       flux_no_grav(k,2,1) = q(k,3)
       flux_no_grav(k,2,2) = v(k)*q(k,2)
       flux_no_grav(k,2,3) = v(k)*q(k,3) 
    end do

  end function flux_no_grav

  function flux_1(q)
    real(kind=GRID_SR), intent(in) :: q(_SWE_DG_ORDER+1,3)
    real(kind=GRID_SR)             :: flux_1(_SWE_DG_ORDER+1,3)
    real(kind=GRID_SR)             :: u(_SWE_DG_ORDER+1), v(_SWE_DG_ORDER+1)
    integer                        :: k
    
    !$omp simd
    do k=1,_SWE_DG_ORDER+1
       u(k) = q(k,2)/q(k,1)
       v(k) = q(k,3)/q(k,1)
    end do
    
    !$omp simd
    do k=1,_SWE_DG_ORDER+1       
       flux_1(k,1) = q(k,2)
       flux_1(k,2) = u(k)*q(k,2) + 0.5_GRID_SR * g * q(k,1)**2
       flux_1(k,3) = u(k)*q(k,3)
    end do
  end function flux_1

  function flux_2(q)
    real(kind=GRID_SR), intent(in) :: q(_SWE_DG_ORDER+1,3)
    real(kind=GRID_SR)             :: flux_2(_SWE_DG_ORDER+1,3)
    real(kind=GRID_SR)             :: u(_SWE_DG_ORDER+1), v(_SWE_DG_ORDER+1)
    integer                        :: k
    
    !$omp simd
    do k=1,_SWE_DG_ORDER+1
       u(k) = q(k,2)/q(k,1)
       v(k) = q(k,3)/q(k,1)
    end do
    
    !$omp simd
    do k=1,_SWE_DG_ORDER+1
       flux_2(k,1) = q(k,3)
       flux_2(k,2) = v(k)*q(k,2)
       flux_2(k,3) = v(k)*q(k,3) + 0.5_GRID_SR * g * q(k,1)**2
    end do
  end function flux_2


end module SWE_PDE

#endif
