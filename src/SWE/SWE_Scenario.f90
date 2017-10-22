
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


#if defined (_SWE_SCENARIO_RESTING_LAKE)
MODULE SWE_Scenario_resting_lake
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
        
        !bathymetry = -4
        bathymetry = -1.5 + 1/(x(1)+100_GRID_SR)
    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
        real (kind = GRID_SR), intent(in) :: x(2)
        type(t_dof_state) :: Q
        
        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
        Q%h = 5.0_GRID_SR
    end function

END MODULE SWE_Scenario_resting_lake
#endif


!********************
!* Radial Dam Break *
!********************
#if defined (_SWE_SCENARIO_RADIAL_DAM_BREAK)
MODULE SWE_Scenario_radial_dam_break
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
            Q%h = 10.0_GRID_SR
        else 
            Q%h = 0.0_GRID_SR
        end if
    end function

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
! Domain: [-2, 2]²
! Bathymetry: b(x,y) = 0.1 (x² + y²)
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

        bathymetry = 0.1 * (x(1)*x(1) + x(2)*x(2))

    end function
    
    function SWE_Scenario_get_initial_Q(x) result(Q)
      real (kind = GRID_SR), intent(in) :: x(2)
      type(t_dof_state) :: Q
      double precision :: w, t, sinwt, coswt, b
      double precision :: x_scal(2)
        
      b = SWE_Scenario_get_bathymetry(x)

        
        t = 0.0
        w = sqrt(0.2 * g)
        sinwt = 0.0_GRID_SR ! t = 0
        coswt = 1.0_GRID_SR ! t = 0

        Q%h = max( 0.0_GRID_SR, 0.05 * (2*x(1)*coswt + 2*x(2)*sinwt) + 0.075  - b )
        
        Q%p(1) = -0.5*w*sinwt * (Q%h) 
        Q%p(2) =  0.5*w*coswt * (Q%h)
        
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
        
        scaling = 200.0_GRID_SR
    end function

    function SWE_Scenario_get_offset() result(offset)
        real (kind = GRID_SR) :: offset(2)
        
        offset = SWE_Scenario_get_scaling() * [-0.5_GRID_SR, -0.5_GRID_SR]
    end function
    
    function SWE_Scenario_get_bathymetry(x) result(bathymetry)
        real (kind = GRID_SR), intent(in) :: x(2)
        real (kind = GRID_SR) :: bathymetry
        real (kind = GRID_SR) :: x_0(2)
        real (kind = GRID_SR) :: desc=0.5
        real (kind = GRID_SR) :: depth=-25.0
                
        x_0 =[0.0,0.0]

        if(x(1)< x_0(1)) then
           bathymetry = depth
        else
           bathymetry = depth+x(1)*desc
        end if
        
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
        
        gamma=sqrt(3*H/(4*d))
        z = gamma*(x(1)-x_s(1))/d
        sech=2/(exp(z)+exp(-z))

        Q%h=sech**2 * (H) + d

        Q%p = [0.0_GRID_SR, 0.0_GRID_SR]
    end function

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
#   endif

END MODULE SWE_Scenario
