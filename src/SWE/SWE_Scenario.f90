! Artificial scenario selector for the SWE scenario. 
! To add a new scenario:
! 1) create a new module according to the template below
! 2) add an option for it in the Sconstruct file (swe_scenario argument) and define a macro if it is chosen
! 3) "USE" its module in the SWE_Scenario module (at the end of this file)
! 4) Don't forget to use #if defined(_NEW_MACRO) to avoid conflicts!

#include "Compilation_control.f90"


!*******************
!* MODULE TEMPLATE *
!*******************
#if 0
MODULE SWE_Scenario_template
    use Samoa_swe
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
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry
        
        bathymetry = 0.0_GRID_SR
    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q
        
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        Q%h = 0.0_GRID_SR
    end function

END MODULE SWE_Scenario_template
#endif

#if defined (_SWE_SCENARIO_ALL_RAREFACTION)
MODULE SWE_Scenario_all_rarefaction
    use Samoa_swe
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 4.0_GRID_SR
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = SWE_Scenario_get_scaling() * [-0.5_GRID_SR, -0.5_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry
        
        !        bathymetry = -0.5_GRID_SR
        bathymetry = 0.0_GRID_SR        
    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q


!        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]        
        ! if(x(1) > -1.5 .and. x(1) < 1.5) then
        !    if(x(2) > -1.5 .and. x(2) < 1.5) then
           
              if(x(2) < -x(1))then
                 Q%p = [-0.5_GRID_SR, -0.5_GRID_SR]/sqrt(2.0_GRID_SR)
              else
                 Q%p = [0.5_GRID_SR, 0.5_GRID_SR]/sqrt(2.0_GRID_SR)
              end if
        !    end if
        ! end if
        Q%h = 1.0
    end function

END MODULE SWE_Scenario_all_rarefaction
#endif



#if defined (_SWE_SCENARIO_GAUSSIAN_CURVE)
MODULE SWE_Scenario_gaussian_curve
    use Samoa_swe
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 20.0_GRID_SR
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = SWE_Scenario_get_scaling() * [-0.5_GRID_SR, -0.5_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry
        
        bathymetry = 0.0_GRID_SR
    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: height_offset=4.0_GRID_SR
        real (kind = GRID_SR) :: curve_height=2.0_GRID_SR
        type(t_dof_state) :: Q
        
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]

        if(NORM2(x)<8.0_GRID_SR) then
           Q%h = 1/6.28 * exp(-0.5_GRID_SR*(x(1)**2+x(2)**2)) * 8.0_GRID_SR *curve_height + height_offset
        else
           Q%h=height_offset
        end if
    end function

  END MODULE SWE_Scenario_gaussian_curve
#endif

#if defined (_SWE_SCENARIO_SPLASHING_POOL)
MODULE SWE_Scenario_splashing_pool
    use Samoa_swe
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 20.0_GRID_SR
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = SWE_Scenario_get_scaling() * [-0.5_GRID_SR, -0.5_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry
        real (kind = GRID_SR) :: height_offset=4.0_GRID_SR
!        bathymetry = 0.0
        bathymetry = 0.5+x(1) / cfg%scaling 
    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: height_offset=4.0_GRID_SR
        real (kind = GRID_SR) :: curve_height=1.0_GRID_SR
        type(t_dof_state) :: Q
        
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]

        Q%h = 0.5+x(1) / cfg%scaling *curve_height + height_offset
    end function

  END MODULE SWE_Scenario_splashing_pool
#endif



!********************
!* Resting lake: center island *
!********************
#if defined(_SWE_SCENARIO_RESTING_LAKE)
MODULE SWE_Scenario_resting_lake
    use SWE_data_types
    use iso_c_binding
!    ! use Samoa_swe
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
#   if defined(_BOUNDARY_FUNC)
      public SWE_Scenario_get_boundary_height
#   endif
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 1.0_GRID_SR
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = SWE_Scenario_get_scaling() * [0.0_GRID_SR, 0.0_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry
        
        bathymetry = x(1) *4 + 16 * x(2) - 20.0_GRID_SR
            
    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q
        
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        Q%h = max(0.0_GRID_SR, 0.0_GRID_SR)
    end function

#   if defined(_BOUNDARY_FUNC)
      function SWE_Scenario_get_boundary_height(b, t) result(height)
        real (kind = GRID_SR), intent(in)   :: b
        real (kind = GRID_SR), intent(in)   :: t
        real (kind = GRID_SR)               :: height

        height = max(0.0_GRID_SR, (0.1_GRID_SR + (sin(t * 10.0_GRID_SR) * 0.05_GRID_SR)) - b)
      end function
#   endif

  END MODULE SWE_Scenario_resting_lake
#endif


!********************
!* Resting lake 2: overlapping *
!********************
#if defined(_SWE_SCENARIO_RESTING_LAKE2)
MODULE SWE_Scenario_resting_lake2
    use SWE_data_types
    use iso_c_binding
!    ! use Samoa_swe
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 1.0_GRID_SR
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = SWE_Scenario_get_scaling() * [0.0_GRID_SR, 0.0_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry
        
        logical               :: om_1, om_2, om_3, om_4

        om_1 = norm2(x - (/ 0.35_GRID_SR, 0.65_GRID_SR /)) .lt. 0.1_GRID_SR
        om_2 = norm2(x - (/ 0.55_GRID_SR, 0.45_GRID_SR /)) .lt. 0.1_GRID_SR
        om_3 = (abs(x(1) - 0.47_GRID_SI) .lt. 0.25_GRID_SR) .and. (abs(x(2) - 0.55_GRID_SI) .lt. 0.25_GRID_SR)
        om_4 = norm2(x - (/ 0.5_GRID_SR, 0.5_GRID_SR /)) .lt. 0.45_GRID_SR

        if (om_1) then
          bathymetry = 0.12_GRID_SR
        else if (om_2) then
          bathymetry = 0.05_GRID_SR
        else if (om_3 .and. .not. (om_1 .or. om_2)) then
          bathymetry = 0.07_GRID_SR
        else if (om_4 .and. .not. om_3) then
          bathymetry = 0.03_GRID_SR
        else
          bathymetry = 0.0_GRID_SR
       end if

       bathymetry = bathymetry - 0.1_GRID_SR
    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q
        
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        !        Q%h = max(0.0_GRID_SR, - SWE_Scenario_get_bathymetry(x))
        Q%h = 0.0_GRID_SR
    end function

  END MODULE SWE_Scenario_resting_lake2
#endif


!********************
!* Linear beach *
!********************
#if defined(_SWE_SCENARIO_LINEAR_BEACH)
MODULE SWE_Scenario_linear_beach
    use SWE_data_types
    use iso_c_binding
!    ! use Samoa_swe
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 50400.0_GRID_SR                                                                ! 10 * L
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = [-400.0_GRID_SR, 0.0_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry

        bathymetry = -0.1_GRID_SR * x(1)                                                        ! Linear beach: L - alpha * x
    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q
        
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        Q%h = max((((0.006_GRID_SR * exp(-0.4444_GRID_SR    * (((x(1) / 5000.0_GRID_SR) - 4.1209_GRID_SR) ** 2))) - & ! Two gaussian terms of the submarine landslide n-wave
                    (0.018_GRID_SR * exp(-4.0_GRID_SR       * (((x(1) / 5000.0_GRID_SR) - 1.6384_GRID_SR) ** 2)))) &  ! x / L transforms x to non-dimensional form
                  * 500.0_GRID_SR) &                                                                     ! Backtransform eta: x *= alpha * L
                - SWE_Scenario_get_bathymetry(x), 0.0_GRID_SR)                                                     ! Water level to water height
    end function

  END MODULE SWE_Scenario_linear_beach
#endif


!********************
!* Long wave in basin *
!********************
#if defined(_SWE_SCENARIO_LONGWAVE_BASIN)
MODULE SWE_Scenario_longwave_basin
    use SWE_data_types
    use iso_c_binding
!    ! use Samoa_swe
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 8000.0_GRID_SR
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = [-4000.0_GRID_SR, -4000.0_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry

        bathymetry = ((x(1) ** 2) + (x(2) ** 2)) / (2500.0_GRID_SR ** 2)
    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q

        real (kind = GRID_SR) :: A
        
        A = ((2500.0_GRID_SR ** 4) - (2000.0_GRID_SR ** 4)) / &
            ((2500.0_GRID_SR ** 4) + (2000.0_GRID_SR ** 4))

        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        Q%h = max(( &
                    (sqrt(1.0_GRID_SR - (A ** 2))) / &
                    (1.0_GRID_SR      - A) &
                  ) - ( &
                    (((x(1) ** 2) + (x(2) ** 2))  * (1.0_GRID_SR - (A ** 2))) / &
                    ((2500.0_GRID_SR ** 2)          * ((1.0_GRID_SR - A) ** 2)) &
                  ), 0.0_GRID_SR)
    end function

  END MODULE SWE_Scenario_longwave_basin
#endif


!********************
!* Concical island*
!********************
#if defined(_SWE_SCENARIO_CONICAL_ISLAND_A) || defined(_SWE_SCENARIO_CONICAL_ISLAND_C)
MODULE SWE_Scenario_conical_island
    use SWE_data_types
    use iso_c_binding
!    ! use Samoa_swe
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
#   if defined(_BOUNDARY_FUNC)
      public SWE_Scenario_get_boundary_height
#   endif
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 28.2_GRID_SR
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = [0.0_GRID_SR, 0.0_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry

        real(GRID_SR)           :: r

        r = sqrt(((x(1) - 12.96_GRID_SR) ** 2) + ((x(2) - 13.8_GRID_SR) ** 2))
        
        if(r .le. 1.1_GRID_SR) then
            bathymetry = 0.625_GRID_SR
        else if(r .ge. 1.1_GRID_SR .and. r .le. 3.6_GRID_SR) then
            bathymetry = (3.6_GRID_SR - r) / 4.0_GRID_SR
        else
            bathymetry = 0
        end if
    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q
        
        real(GRID_SR)               :: h0
        
        h0 = 0.32_GRID_SR

        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        Q%h = max(0.0_GRID_SR, h0 - SWE_Scenario_get_bathymetry(x))
    end function

#   if defined(_BOUNDARY_FUNC)
      function SWE_Scenario_get_boundary_height(b, t) result(height)
        real (kind = GRID_SR), intent(in)   :: b
        real (kind = GRID_SR), intent(in)   :: t
        real (kind = GRID_SR)               :: height

        real(GRID_SR)                       :: h0, K, c, a, bigT, x0

        h0 = 0.32_GRID_SR
#       if defined(_SWE_SCENARIO_CONICAL_ISLAND_A)
            a = 0.014_GRID_SR
            bigT = 8.85_GRID_SR
            x0 = 5.76_GRID_SR
#       elif defined(_SWE_SCENARIO_CONICAL_ISLAND_C)
            a = 0.057_GRID_SR
            bigT = 7.77_GRID_SR
            x0 = 7.56_GRID_SR
#       endif

        K = sqrt((3.0_GRID_SR * a) / (4.0_GRID_SR * (h0 ** 3)))
        c = sqrt(_GRAV_CONSTANT * h0) * (1.0_GRID_SR + (a / (2.0_GRID_SR * h0)))
        height = h0 + (a * (( &
            1.0_GRID_SR / cosh(K * ((c * bigT) - (c * t) - x0))) ** 2))
      end function
#   endif

  END MODULE SWE_Scenario_conical_island
#endif


!********************
!* Radial Dam Break *
!********************
#if defined (_SWE_SCENARIO_RADIAL_DAM_BREAK)
MODULE SWE_Scenario_radial_dam_break
    use Samoa_swe
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
#   if defined(_BOUNDARY_FUNC)
        public SWE_Scenario_get_boundary_height
#   endif
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 20.0_GRID_SR
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = SWE_Scenario_get_scaling() * [-0.5_GRID_SR, -0.5_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry
        
        if (x(1)*x(1) + x(2)*x(2) < 9.0_GRID_SR) then
            bathymetry = -20.0_GRID_SR
        else 
            bathymetry = -22.0_GRID_SR
        end if
    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q
        
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        
        if (x(1)*x(1) + x(2)*x(2) < 4.5_GRID_SR) then
            Q%h = 20.0_GRID_SR
        else 
            Q%h = 0.0_GRID_SR
        end if
    end function

#   if defined(_BOUNDARY_FUNC)
        function SWE_Scenario_get_boundary_height(b, t) result(height)
            real (kind = GRID_SR), intent(in)   :: b
            real (kind = GRID_SR), intent(in)   :: t
            real (kind = GRID_SR)               :: height

            height = (sin(t * 10.0_GRID_SR) + 1.0_GRID_SR) * 7.0_GRID_SR
        end function
#   endif

END MODULE SWE_Scenario_radial_dam_break
#endif


!********************
!* Linear Dam Break *
!********************
#if defined (_SWE_SCENARIO_LINEAR_DAM_BREAK)
MODULE SWE_Scenario_linear_dam_break
    use Samoa_swe
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 100.0_GRID_SR
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = SWE_Scenario_get_scaling() * [-0.5_GRID_SR, -0.5_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry
        
        bathymetry = 0.0_GRID_SR
    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q
        
        real (kind = GRID_SR), parameter :: hL = 0.005_GRID_SR, hR = 0.001_GRID_SR
        
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        
!        if (x(1) < 0.0_GRID_SR) then
        if ( x(1) < 5.0_GRID_SR ) then
            Q%h = hL
        else 
            Q%h = hR
        end if
    end function

END MODULE SWE_Scenario_linear_dam_break
#endif

!********************
!* Oscillating Lake *
!********************
! Proposed in: 
! [1] Gallardo et al., 2007. On a well-balanced high-order finite volume scheme for shallow water equations with topography and dry areas. Journal of Computational Physics, volume 227
! [2] Meister & Ortleb, 2014. On unconditionally positive implicit time integration for the dg scheme applied to shallow water flows. International Journal for Numerical Methods in Fluids, volume 76, 69-94.
!
! Domain: [-2, 2]�
! Bathymetry: b(x,y) = 0.1 (x� + y�)
!
! Analytic solution:
! H(x,y,t) = max{ 0.0, 0.05 * (2x cos(wt) + 2y sin(wt) + 0.075 - b(x,y) ) }
! u(x,y,t) = -0.5w sin(wt)
! v(x,y,t) =  0.5w cos(wt)
! --> with w = sqrt( 0.2 g )

#if defined (_SWE_SCENARIO_OSCILLATING_LAKE)
MODULE SWE_Scenario_oscillating_lake
    use Samoa_swe
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 4.0_GRID_SR
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = SWE_Scenario_get_scaling() * [-0.5_GRID_SR, -0.5_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR):: x_scal(2)
        real (kind = GRID_SR) :: bathymetry

        bathymetry = 0.1_GRID_SR * (x(1)*x(1) + x(2)*x(2))

    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
      real (kind = GRID_SR), intent(in) :: x(2)
      type(t_dof_state) :: Q
      double precision :: w, t, sinwt, coswt, b
      double precision :: x_scal(2)
        
      b = SWE_Scenario_get_bathymetry(x)

        
        t = 0.0
        w = sqrt(0.2_GRID_SR * g)
        sinwt = 0.0_GRID_SR ! t = 0
        coswt = 1.0_GRID_SR ! t = 0

        Q%h = max( 0.0_GRID_SR, 0.05_GRID_SR * (2.0_GRID_SR*x(1)*coswt + 2.0_GRID_SR*x(2)*sinwt) + 0.075_GRID_SR  - b )
        
        Q%p(1) = -0.5_GRID_SR*w*sinwt * (Q%h) 
        Q%p(2) =  0.5_GRID_SR*w*coswt * (Q%h)
        
        Q%h = Q%h + b
        
    end function

END MODULE SWE_Scenario_oscillating_lake
#endif


#if defined(_SWE_SCENARIO_PARABOLIC_ISLE)
MODULE SWE_Scenario_resting_isle
    use Samoa_swe
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
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry
        
        if(NORM2(x) <= 3.0) then
        !if(x(1) <= 0.0) then        
           bathymetry = -0.5_GRID_SR + 1.0_GRID_SR * (1.0_GRID_SR- (NORM2(x)/3.0_GRID_SR)**(2))
        !bathymetry = -10.0_GRID_SR + cos(3.1415 * 0.5 * NORM2(x)/5.0)**2 * 20.0
        else
           bathymetry = -0.5_GRID_SR
        end if

        ! if(norm2(x)<2.0)then
        !    bathymetry = -10.0 + (4-(x(1)**2 + x(2)**2))*15.0
        ! end if

    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
      real (kind = GRID_SR), intent(in) :: x(2)
      real (kind = GRID_SR) :: b
      type(t_dof_state) :: Q
        
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        Q%h = 0.0_GRID_SR
        ! if (SWE_Scenario_get_bathymetry(x) > 0 ) then
        b=SWE_Scenario_get_bathymetry(x)
        ! end if

!         if(x(1) >= 0.0) then
! !        if(NORM2(x) <= 3.0) then
!            Q%h = -10.0_GRID_SR + 20.0 * (1.0_GRID_SR- (NORM2(x)/3.0_GRID_SR)**(2))
!            !bathymetry = -10.0_GRID_SR + cos(3.1415 * 0.5 * NORM2(x)/5.0)**2 * 20.0
!            !Q%h = x(1) * 0.5                      
!         elseo
!            Q%h = 0.0_GRID_SR
!         end if

        if(x(1)<-2.0 .and. x(1) > -2.5)then
!         !    Q%h = 8.0_GRID_SR
           Q%h = 0.2_GRID_SR* (0.25**2-(x(1)+2.25)**2)
        else
           Q%h = 0.0_GRID_SR
        end if
        
!        Q%h = max(Q%h,b)

!        Q%h = 0.0_GRID_SR
    end function

END MODULE SWE_Scenario_resting_isle
#endif


#if defined(_SWE_SCENARIO_SINGLE_WAVE)
MODULE SWE_Scenario_single_wave_on_the_beach
    use Samoa_swe
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 50400.0_GRID_SR                                                                ! 10 * L
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = [-400.0_GRID_SR, 0.0_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry

        bathymetry = -0.1_GRID_SR * x(1)                                                        ! Linear beach: L - alpha * x
    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q
        
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        Q%h = max((((0.006_GRID_SR * exp(-0.4444_GRID_SR    * (((x(1) / 5000.0_GRID_SR) - 4.1209_GRID_SR) ** 2))) - & ! Two gaussian terms of the submarine landslide n-wave
                    (0.018_GRID_SR * exp(-4.0_GRID_SR       * (((x(1) / 5000.0_GRID_SR) - 1.6384_GRID_SR) ** 2)))) &  ! x / L transforms x to non-dimensional form
                  * 500.0_GRID_SR) &                                                                     ! Backtransform eta: x *= alpha * L
                    , SWE_Scenario_get_bathymetry(x))                                                     ! Water level to water height
    end function SWE_Scenario_get_initial_Q

END MODULE SWE_Scenario_single_wave_on_the_beach
#endif


#if defined(_SWE_SCENARIO_CONVERGENCE_TEST)
MODULE  SWE_Convergence_test
    use Samoa_swe
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 6.0_GRID_SR
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = -[0.5_GRID_SR, 0.5_GRID_SR] * SWE_Scenario_get_scaling()
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry
        
        bathymetry=sin(8*atan(1.0)*x(1)) + cos(8*atan(1.0)*x(2))
                
    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q
        real(kind=GRID_SR) :: d = 20.0
        real(kind=GRID_SR) :: H = 2.0
        real(kind=GRID_SR) :: x_s(2) = [-25,-25]
        real(kind=GRID_SR) :: sech
        real(kind=GRID_SR) :: gamma
        real(kind=GRID_SR) :: z
        
        Q%h= 10.0_GRID_SR + exp(sin(8*atan(1.0))*x(1))*cos(8*atan(1.0)*x(2)) + SWE_Scenario_get_bathymetry(x)

        Q%p(1) = sin(cos(8*atan(1.0)*x(1)))*sin(8*atan(1.0)*x(2))
        Q%p(2) = cos(sin(8*atan(1.0)*x(1)))*cos(8*atan(1.0)*x(1))
    end function

END MODULE SWE_Convergence_test
#endif

#if defined (_SWE_SCENARIO_SMOOTH_WAVE)
MODULE SWE_Scenario_smooth_wave
    use Samoa_swe
    public SWE_Scenario_get_scaling, SWE_Scenario_get_offset, SWE_Scenario_get_bathymetry, SWE_Scenario_get_initial_Q
    contains

    function SWE_Scenario_get_scaling() result(scaling)
        real (kind = GRID_SR) :: scaling
        
        scaling = 4.0_GRID_SR
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = SWE_Scenario_get_scaling() * [0.5_GRID_SR, 0.5_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
      real (kind = GRID_SR), intent(in) :: x(2)
      real (kind = GRID_SR):: x_temp      
      real (kind = GRID_SR) :: bathymetry

      x_temp=(x(1)+x(2))/sqrt(2.0_GRID_SR)
      x_temp=x(1)

      bathymetry = -(0.5*x_temp**2/g + g/x_temp)
    end function SWE_Scenario_get_bathymetry
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q
        double precision :: w, t, sinwt, coswt, b
        double precision :: x_scal(2)
        real (kind = GRID_SR):: x_temp      

        b = SWE_Scenario_get_bathymetry(x)        
        t = 0.0

        x_temp=(x(1)+x(2))*sqrt(2.0_GRID_SR)/2.0_GRID_SR
        x_temp=x(1)
        Q%h = (1.0_GRID_SR/x_temp + 1.0_GRID_SR)*g

        Q%p(1) = g+x_temp*g
        Q%p(2) = x_temp * Q%h
        Q%p(2) = 0.0_GRID_SR
        
        Q%h = Q%h + b
        
    end function

END MODULE SWE_Scenario_smooth_wave
#endif




MODULE SWE_Scenario

#   if defined(_SWE_SCENARIO_RADIAL_DAM_BREAK)
        USE SWE_Scenario_radial_dam_break
#   elif defined(_SWE_SCENARIO_LINEAR_DAM_BREAK)
        USE SWE_Scenario_linear_dam_break
#   elif defined(_SWE_SCENARIO_OSCILLATING_LAKE)
        USE SWE_Scenario_oscillating_lake
#   elif defined(_SWE_SCENARIO_RESTING_LAKE)
        USE SWE_Scenario_resting_lake
#   elif defined(_SWE_SCENARIO_RESTING_LAKE2)
        USE SWE_Scenario_resting_lake2
#   elif defined(_SWE_SCENARIO_LINEAR_BEACH)
        USE SWE_Scenario_linear_beach
#   elif defined(_SWE_SCENARIO_LONGWAVE_BASIN)
        USE SWE_Scenario_longwave_basin
#   elif defined(_SWE_SCENARIO_CONICAL_ISLAND_A) || defined(_SWE_SCENARIO_CONICAL_ISLAND_C)
        USE SWE_Scenario_conical_island
#   elif defined(_SWE_SCENARIO_GAUSSIAN_CURVE)
        USE SWE_Scenario_gaussian_curve
#   elif defined(_SWE_SCENARIO_SPLASHING_POOL)
        USE SWE_Scenario_splashing_pool
#   elif defined(_SWE_SCENARIO_PARABOLIC_ISLE)
        USE SWE_Scenario_resting_isle
#   elif defined(_SWE_SCENARIO_SINGLE_WAVE)
        USE SWE_Scenario_single_wave_on_the_beach
#   elif defined(_SWE_SCENARIO_CONVERGENCE_TEST)
        USE SWE_Convergence_test
#   elif defined(_SWE_SCENARIO_ALL_RAREFACTION)
        USE SWE_Scenario_all_rarefaction
#   elif defined(_SWE_SCENARIO_SMOOTH_WAVE)
        USE SWE_Scenario_smooth_wave
#   elif defined(_SWE_SCENARIO_ASAGI)
        ! USE SWE_Scenario_asagi
#   endif

END MODULE SWE_Scenario
