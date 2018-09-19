! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE

! SWE2L: THIS FILE WAS MODIFIED FOR THE TWO-LAYER CODE

! Artificial scenario selector for the SWE scenario. 
! To add a new scenario:
! 1) create a new module according to the template below
! 2) add an option for it in the Sconstruct file (swe_scenario argument) and define a macro if it is chosen
! 3) "USE" its module in the SWE_Scenario module (at the end of this file)
! 4) Don't forget to use #if defined(_NEW_MACRO) to avoid conflicts!

#include "Compilation_control.f90"
#if defined _SWE2L

!*******************
!* MODULE TEMPLATE *
!*******************
#if 0
MODULE SWE2L_Scenario_template
    use iso_c_binding
    use Samoa_swe2l
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 10.0_GRID_SR
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = SWE_Scenario_get_scaling() * [-0.5_GRID_SR, -0.5_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = c_double), intent(in) :: x(3)
        real (kind = GRID_SR) :: bathymetry
        
        bathymetry = 0.0_GRID_SR
    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = c_double), intent(in) :: x(3)
        type(t_dof_state) :: Q
        
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        Q%h = 0.0_GRID_SR
        
        Q%p2 = [0.0_GRID_SR, 0.0_GRID_SR]
        Q%h2 = 0.0_GRID_SR
    end function

END MODULE 
#endif



!***************
!* Bowl radial *
!****************

! Scenario with two layers of water based on the GeoClaw scenario:
! https://github.com/clawpack/geoclaw/tree/master/examples/multi-layer/bowl-radial
! https://github.com/clawpack/geoclaw/blob/master/examples/multi-layer/bowl-radial/maketopo.py

!#if defined (_SWE_SCENARIO_RADIAL_BOWL_RADIAL)
MODULE SWE2L_Scenario_bowl_radial
    use iso_c_binding
    use Samoa_swe2l
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 200.0_GRID_SR
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = SWE_Scenario_get_scaling() * [-0.5_GRID_SR, -0.5_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = c_double), intent(in) :: x(3)
        real (kind = GRID_SR) :: bathymetry, zmin
        
        zmin = 80.0_GRID_SR
        bathymetry = 1.0e-2_GRID_SR * (x(1)*x(1) + x(2)*x(2)) - zmin
    end function
    
    function SWE_Scenario_get_initial_Q(x) result (Q)
        real (kind = c_double), intent(in) :: x(3)
        type(t_dof_state) :: Q
        real (kind = GRID_SR) :: b, z
        
        ! All velocities/momenta are zero
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        Q%p2 = [0.0_GRID_SR, 0.0_GRID_SR]
        
        ! Bottom layer surface is -20
        Q%h2 = -20.0_GRID_SR
        
        ! Top layer surface is given by qinit_func in GeoClaw script maketopo.py (hump.xyz)
        ! (from https://github.com/clawpack/geoclaw/blob/master/examples/multi-layer/bowl-radial/maketopo.py )
        z = -(x(1)*x(1) + x(2)*x(2))*0.01_GRID_SR
        if (z > -10_GRID_SR) then
            z = 4.0_GRID_SR * exp(z)
        else
            z = 0.0_GRID_SR
        end if
        
        if (z >= 0.05) then
            Q%h = z 
        else
            Q%h = 0.0_GRID_SR
        end if
    end function

END MODULE 
!#endif


MODULE SWE2L_Scenario

    ! Currently there is only one scenario for SWE2L (bowl_radial)
    ! Thus, the current implementation igores the SCons option swe_scenario and associated #defines.
    
    ! If you want to create additional scenarios, try to use the same organization
    ! as we used for the single-layer SWEs (check src/SWE/SWE_Scenario.f90)
    
    USE SWE2L_Scenario_bowl_radial

END MODULE

#endif
