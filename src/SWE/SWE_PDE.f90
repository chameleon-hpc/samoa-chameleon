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


  function flux(q,N)
    real(kind=GRID_SR)             :: flux(N,2,3)
    real(kind=GRID_SR), intent(in) :: q(N,3)
    integer                        :: N

    flux(:,1,1) = q(:,2)
    flux(:,1,2) = q(:,2)**2/q(:,1) + 0.5_GRID_SR * g * q(:,1)**2
    flux(:,1,3) = q(:,2)*q(:,3)/q(:,1)

    flux(:,2,1) = q(:,3)
    flux(:,2,2) = q(:,2)*q(:,3)/q(:,1)
    flux(:,2,3) = q(:,3)**2/q(:,1) + 0.5_GRID_SR * g * q(:,1)**2

  end function flux

  function flux_1(q,N)
    real(kind=GRID_SR)             ::flux_1(N,3)
    real(kind=GRID_SR) ,intent(in) ::q(N,3)
    integer :: N
    flux_1(:,1) = q(:,2)
    flux_1(:,2) = q(:,2)**2/q(:,1) + 0.5_GRID_SR * g * q(:,1)**2
    !flux_1(:,2) = q(:,2)**2/q(:,1)
    flux_1(:,3) = q(:,2)*q(:,3)/q(:,1)
  end function flux_1

  function flux_2(q,N)
    real(kind=GRID_SR)             ::flux_2(N,3)
    real(kind=GRID_SR) ,intent(in) ::q(N,3)
    integer :: N
    flux_2(:,1) = q(:,3)
    flux_2(:,2) = q(:,2)*q(:,3)/q(:,1)
    flux_2(:,3) = q(:,3)**2/q(:,1) + 0.5_GRID_SR * g * q(:,1)**2
    !flux_2(:,3) = q(:,3)**2/q(:,1)
  end function flux_2

  function flux_no_grav(q,N)
    real(kind=GRID_SR)             ::flux_no_grav(N,2,3)
    real(kind=GRID_SR) ,intent(in) ::q(N,3)
    integer :: N

    flux_no_grav(:,1,1) = q(:,2)
    flux_no_grav(:,1,2) = q(:,2)**2/q(:,1)
    flux_no_grav(:,1,3) = q(:,2)*q(:,3)/q(:,1)

    flux_no_grav(:,2,1) = q(:,3)
    flux_no_grav(:,2,2) = q(:,2)*q(:,3)/q(:,1)
    flux_no_grav(:,2,3) = q(:,3)**2/q(:,1)

  end function flux_no_grav

end module SWE_PDE

#endif
