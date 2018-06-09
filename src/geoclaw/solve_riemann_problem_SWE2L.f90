subroutine solve_rpn2(meqn,mwaves, n_index,t_index, h_l,h_r, hu_l,hu_r, hv_l,hv_r, b_l,b_r, h_hat_l,h_hat_r, fw_priv, sw_priv)


    use amr_module, only: mcapa

    use geoclaw_module, only: g => grav, rho, pi, earth_radius

    use multilayer_module, only: num_layers, eigen_func
    use multilayer_module, only: dry_tolerance, aux_layer_index, r
    use multilayer_module, only: eigen_method, inundation_method
    use multilayer_module, only: eigen_func, inundation_eigen_func
    use statistics

    implicit none
    
    ! Input arguments
    integer, intent(in) :: meqn,mwaves, n_index,t_index
    real(kind=8), dimension(2), intent(inout) :: h_l, h_r, hu_l, hu_r, hv_l, hv_r
    real(kind=8), intent(inout) :: b_l, b_r
    real(kind=8), dimension(2), intent(inout) :: h_hat_l, h_hat_r

    ! Output arguments
    real(kind=8), dimension(6, 6), intent(out) :: fw_priv
    real(kind=8), dimension(6), intent(out) :: sw_priv
    
    ! Counters
    integer :: i, j, mw, info, info2
    integer :: layer_index
    
    ! Physics
    real(kind=8), dimension(2) :: u_l, u_r, v_l, v_r
    real(kind=8), dimension(2) :: h_ave, momentum_transfer
    real(kind=8) :: flux_transfer_l, flux_transfer_r, lambda(6)
    
    ! Input for single-layer solver
    real(kind=8) :: hL, hR, huL, huR, hvL, hvR

    ! real(kind=8) :: advected_speed, eta_l, eta_r, gamma_l, gamma_r, kappa_l, kappa_r, w_normal, w_transverse

    ! Solver variables
    integer :: num_dry_states
    real(kind=8), dimension(2) :: eigen_h_l, eigen_h_r
    real(kind=8), dimension(2) :: eigen_u_l, eigen_u_r
    real(kind=8), dimension(2) :: eigen_v_l, eigen_v_r
    real(kind=8), dimension(2) :: flux_h_l, flux_h_r
    real(kind=8), dimension(2) :: flux_hu_l, flux_hu_r
    real(kind=8), dimension(2) :: flux_hv_l, flux_hv_r
    real(kind=8), dimension(2) :: flux_u_l, flux_u_r
    real(kind=8), dimension(2) :: flux_v_l, flux_v_r
    real(kind=8), dimension(6) :: delta, flux_r, flux_l
    real(kind=8) :: delta1,delta2,delta3,delta4,delta5,delta6
    integer, dimension(6) :: pivot
    real(kind=8) :: eig_vec(6,6),   A11, A12, A13, A14, A15, A16, &
                                    A21, A22, A23, A24, A25, A26, &
                                    A31, A32, A33, A34, A35, A36, &
                                    A41, A42, A43, A44, A45, A46, &
                                    A51, A52, A53, A54, A55, A56, &
                                    A61, A62, A63, A64, A65, A66
    real(kind=8) :: beta(6), alpha(4), fw(3, 3), sw(3)
    logical, dimension(2) :: dry_state_l, dry_state_r
    logical :: inundation

    external dgesv
    
    dry_state_l = .false.
    dry_state_r = .false.
    inundation = .false.
    
    do j=1,2
        h_ave(:) = 0.5d0 * (h_l(:) + h_r(:))
        
        ! Check for dry states
        if (h_l(j) < dry_tolerance(j)) then
            dry_state_l(j) = .true.
            hu_l(j) = 0.d0
            hv_l(j) = 0.d0
            u_l(j) = 0.d0
            v_l(j) = 0.d0
        else
            u_l(j) = hu_l(j) / h_l(j)
            v_l(j) = hv_l(j) / h_l(j)
        endif
        if (h_r(j) < dry_tolerance(j)) then
            dry_state_r(j) = .true.
            hu_r(j) = 0.d0
            hv_r(j) = 0.d0
            u_r(j) = 0.d0
            v_r(j) = 0.d0
        else
            u_r(j) = hu_r(j) / h_r(j)
            v_r(j) = hv_r(j) / h_r(j)
        endif
    end do
    
    ! For ease of checking below, count up number of dry states
    num_dry_states = 0
    do mw=1,2
        if (dry_state_l(mw)) then
            num_dry_states = num_dry_states + 1
        end if
        if (dry_state_r(mw)) then
            num_dry_states = num_dry_states + 1
        end if
    end do

    ! ===============================
    !  Completely dry cell - (T T T T)
    ! ===============================
    if (num_dry_states == 4) then
        sw_priv = 0.d0
        fw_priv = 0.d0

    ! ===============================================
    !  Single-layer problem - 
    !   (T F T F) or (F T F T) or (F T T T) or
    !   (T F T T) or (T T F T) or (T T T F)
    ! ===============================================
    !  Here check explicitly for the completely wet single-layer cases and 
    !  rely on the count to determine the other cases
    else if ((      dry_state_l(1) .and. .not. dry_state_l(2) .and.        &  ! T F T F
                    dry_state_r(1) .and. .not. dry_state_r(2))       .or.  &
                (.not. dry_state_l(1) .and.       dry_state_l(2) .and.        &  ! F T F T
                .not. dry_state_r(1) .and.       dry_state_r(2))       .or.  &
                num_dry_states == 3) then

        ! Set wet layer index so that we can handle both the bottom and top
        ! layers being dry while the other is wet
        if (.not.dry_state_l(1) .or. .not.dry_state_r(1)) then
            layer_index = 1
            hL = h_l(1)
            hR = h_r(1)
            huL = hu_l(1)
            huR = hu_r(1)
            hvR = hv_r(1)
            hvL = hv_l(1)
        else if (.not.dry_state_l(2) .or. .not. dry_state_r(2)) then
            layer_index = 2
            hL = h_l(2)
            hR = h_r(2)
            huL = hu_l(2)
            huR = hu_r(2)
            hvL = hv_l(2)
            hvR = hv_r(2)
        else
!             print *, "Invalid dry layer state reached."
!             print *, "dry states: ", dry_state_l, dry_state_r
!             print *, "        left            |             right"
!             print *, "====================================================="
!             print "(2d16.8)", h_l(1), h_r(1)
!             print "(2d16.8)", hu_l(1), hu_r(1)
!             print "(2d16.8)", hv_l(1), hv_r(1)
!             print "(2d16.8)", h_l(2), h_r(2)
!             print "(2d16.8)", hu_l(2), hu_r(2)
!             print "(2d16.8)", hv_l(2), hv_r(2)
!             print "(2d16.8)", b_l, b_r
!             stop
        end if

        !DIR$ FORCEINLINE
        call solve_single_layer_rp(layer_index, hL, hR, huL, huR,      &
                                              hvL, hvR, b_l, b_r,      &
                                              fw, sw)
          

        ! Update speeds and waves
        ! Note that we represent all the waves in the first three arrays
        ! so it does not directly correspond to the two-layer case's wave
        ! structure
        sw_priv = 0.d0
        fw_priv = 0.d0
        
        if (layer_index == 1) then
            sw_priv(1:3) = sw(:)
            fw_priv(1, 1:3) = fw(1, :) * rho(layer_index)
            fw_priv(n_index, 1:3) = fw(2, :) * rho(layer_index)
            fw_priv(t_index, 1:3) = fw(3, :) * rho(layer_index)
        else
            sw_priv(4:6) = sw(:)
            fw_priv(1, 4:6) = fw(1, :) * rho(layer_index)
            fw_priv(n_index, 4:6) = fw(2, :) * rho(layer_index)
            fw_priv(t_index, 4:6) = fw(3, :) * rho(layer_index)
        end if

        ! Go on to next cell, lat-long and fluctuation calculations are 
        ! outside of this loop

    ! ======================================================================
    !  Multi-layer system must be solved
    !   In each case a special eigen system is solved and then the flux
    !     difference is evaluated.
    !   Note that the parameter *eigen_method* controls the method being 
    !     used if the cells are completely wet or there exists a wall 
    !     boundary problem in the bottom layer.  Otherwise the parameter 
    !     *inundation_method* is used.
    ! ======================================================================
    else
        ! By default fill in the eigen and flux evaluation states with their
        ! side values
        if (eigen_method == 1) then
            eigen_h_l = h_hat_l
            eigen_h_r = h_hat_r
        else
            eigen_h_l = h_l
            eigen_h_r = h_r
        end if
        eigen_u_l = u_l
        eigen_u_r = u_r
        eigen_v_l = v_l
        eigen_v_r = v_r

        flux_h_l = h_l
        flux_h_r = h_r
        flux_hu_l = hu_l
        flux_hu_r = hu_r
        flux_hv_l = hv_l
        flux_hv_r = hv_r
        flux_u_l = u_l
        flux_u_r = u_r
        flux_v_l = v_l
        flux_v_r = v_r

        ! Also intialize other flux evaluation stuff
        flux_transfer_r = 0.d0
        flux_transfer_l = 0.d0
        momentum_transfer = 0.d0

        ! ==================================================================
        !  Right state is completely dry - (F F T T)
        if (.not. dry_state_l(1) .and. .not. dry_state_l(2) .and.          &
                    dry_state_r(1) .and.       dry_state_r(2)) then

            ! Inundation occurs
            inundation = sum(h_l) + b_l > b_r
            if (inundation) then
!                 print *, "Inundation in this case not yet handled."
!                 print *, "  dry_state = ", dry_state_l, dry_state_r
!                 stop 

            ! Wall boundary
            else
                ! Wall state - Mirror left state onto right
                if (eigen_method /= 1) eigen_h_r = h_l
                eigen_u_r = -u_l
                eigen_v_r =  v_r

                ! Flux evaluation
                flux_h_r = h_l
                flux_hu_r = -hu_l
                flux_u_r = -u_l
                flux_hv_r = hv_l
                flux_v_r = v_l

                flux_transfer_r = 0.d0
                flux_transfer_l = 0.d0
                momentum_transfer(1) = 0.d0
                momentum_transfer(2) = 0.d0
            endif

        ! ==================================================================
        !  Left state is completely dry - (T T F F)
        else if (      dry_state_l(1) .and.       dry_state_l(2) .and.     &
                    .not. dry_state_r(1) .and. .not. dry_state_r(2)) then

            ! Inundation
            inundation = sum(h_r) + b_r > b_l
            if (inundation) then
!                 print *, "Inundation in this case not yet handled."
!                 print *, "  dry_state = ", dry_state_l, dry_state_r
!                 stop 

            ! Wall
            else
                if (eigen_method /= 1) eigen_h_l = h_r
                eigen_u_l = -u_r
                eigen_v_l = v_l
            
                ! Flux evaluation
                flux_h_l = h_r
                flux_hu_l = -hu_r
                flux_u_l = -u_r
                flux_hv_l = hv_r
                flux_v_l = v_r

                flux_transfer_r = 0.d0
                flux_transfer_l = 0.d0
                momentum_transfer(1) = 0.d0
                momentum_transfer(2) = 0.d0

            end if

        ! ==================================================================
        !  Right bottom state is dry - (F F F T)
        else if (.not. dry_state_l(1) .and. .not. dry_state_l(2) .and.     &
                    .not. dry_state_r(1) .and.       dry_state_r(2)) then

            ! Inundation
            inundation = h_l(2) + b_l > b_r
            if (inundation) then
                if (inundation_method == 1 .or.                            &
                    inundation_method == 4) then
                    eigen_h_r = [h_r(1), 0.d0]
                else if (inundation_method == 2 .or.                       &
                            inundation_method == 3 .or.                       &
                            inundation_method == 5) then
                    eigen_h_r = [h_r(1), 0.d0]
                end if

                ! Flux evaluation
                momentum_transfer(1) =  g * rho(1) * h_ave(1) * (h_r(2) - h_l(2) + b_r - b_l)
                momentum_transfer(2) = -g * rho(1) * h_ave(1) * (h_r(2) - h_l(2)) + g * rho(2) * h_ave(2) * (b_r - b_l)
                flux_transfer_r = g * rho(1) * h_r(1) * h_r(2)
                flux_transfer_l = g * rho(1) * h_l(1) * h_l(2)

            ! Wall
            else
                if (eigen_method /= 1) eigen_h_r = [h_r(1), 0.d0]
                eigen_u_r = [u_r(1), -u_l(2)]
                eigen_v_r = [v_r(1),  v_l(2)]

                ! Flux evaluation
                flux_h_r(2) = h_l(2)
                flux_hu_r(2) = -hu_l(2)
                flux_u_r(2) = -u_l(2)
                flux_hv_r(2) = hv_l(2)
                flux_v_r(2) = v_l(2)
            
                flux_transfer_r = 0.d0
                flux_transfer_l = 0.d0
                momentum_transfer(1) = g * rho(1) * h_ave(1) * (b_r - flux_h_l(2) - b_l)
                momentum_transfer(2) = 0.d0

            end if

        ! ==================================================================
        !  Left bottom state is dry - (F T F F)
        else if (.not. dry_state_l(1) .and.       dry_state_l(2) .and.     &
                    .not. dry_state_r(1) .and. .not. dry_state_r(2)) then

            ! Inundation
            inundation = (h_r(2) + b_r > b_l)
            if (inundation) then
                if (inundation_method == 1 .or. inundation_method == 5) then
                    eigen_h_l = [h_l(1), 0.d0]
                else if (inundation_method == 2 .or.                       &
                            inundation_method == 3 .or.                       &
                            inundation_method == 4) then
                    eigen_h_l = [h_l(1), dry_tolerance(1)]
                end if

                ! Flux evaluation
                momentum_transfer(1) =  g * rho(1) * h_ave(1) * (h_r(2) - h_l(2) + b_r - b_l)
                momentum_transfer(2) = -g * rho(1) * h_ave(1) * (h_r(2) - h_l(2)) + g * rho(2) * h_ave(2) * (b_r - b_l)
                flux_transfer_r = g * rho(1) * h_r(1) * h_r(2)
                flux_transfer_l = g * rho(1) * h_l(1) * h_l(2)

            ! Wall
            else
                if (eigen_method /= 1) eigen_h_l = [h_l(1), 0.d0]
                eigen_u_l = [u_l(1), -u_r(2)]
                eigen_v_l = [v_l(1),  v_r(2)]

                ! Flux evaluation
                flux_h_l(2) = h_r(2)
                flux_hu_l(2) = -hu_r(2)
                flux_u_l(2) = -u_r(2)
                flux_hv_l(2) = hv_r(2)
                flux_v_l(2) = v_r(2)
            
                flux_transfer_r = 0.d0
                flux_transfer_l = 0.d0
                momentum_transfer(1) = g * rho(1) * h_ave(1) * (b_r + flux_h_r(2) - b_l)
                momentum_transfer(2) = 0.d0

            end if

        ! ==================================================================
        !  All-states are wet - (F F F F)
!             else if (.not. dry_state_l(1) .and. .not. dry_state_l(2) .and.     &
!                      .not. dry_state_r(1) .and. .not. dry_state_r(2)) then
        else if (num_dry_states == 0) then

            ! Nothing to do for eigenspace evaluation

            ! Flux evaulation
            momentum_transfer(1) =  g * rho(1) * h_ave(1) * (h_r(2) - h_l(2) + b_r - b_l)
            momentum_transfer(2) = -g * rho(1) * h_ave(1) * (h_r(2) - h_l(2)) + g * rho(2) * h_ave(2) * (b_r - b_l)
            flux_transfer_r = g * rho(1) * h_r(1) * h_r(2)
            flux_transfer_l = g * rho(1) * h_l(1) * h_l(2)

        ! ==================================================================
        !  We do not yet handle this case - F F F F and F F F F 
        else
!             print *, "Unhandled dry-state condition reached."
!             print *, "dry states: ", dry_state_l, dry_state_r
!             print *, "        left            |             right"
!             print *, "====================================================="
!             print "(2d16.8)", h_l(1), h_r(1)
!             print "(2d16.8)", hu_l(1), hu_r(1)
!             print "(2d16.8)", hv_l(1), hv_r(1)
!             print "(2d16.8)", h_l(2), h_r(2)
!             print "(2d16.8)", hu_l(2), hu_r(2)
!             print "(2d16.8)", hv_l(2), hv_r(2)
!             print "(2d16.8)", b_l, b_r
!             stop
        end if

        ! ==================================================================
        !  Compute eigen space
        ! ==================================================================
        if (inundation) then
            !DIR$ FORCEINLINE
            call eigen(eigen_h_l, eigen_h_r,               &
                                        eigen_u_l, eigen_u_r,               &
                                        eigen_v_l, eigen_v_r,               &
                                        n_index, t_index,                   &
                                        lambda, eig_vec)
                            
            ! Internal wave corrections
            if (inundation_method == 1) then
                ! Left bottom state dry
                if (.not. dry_state_l(1) .and.       dry_state_l(2) .and.  &
                    .not. dry_state_r(1) .and. .not. dry_state_r(2)) then

                    sw_priv(2) = u_r(2) - 2.d0 * sqrt(g*(1.d0-r)*h_r(2))
                    alpha(2) = r * g * h_r(2) / ((sw_priv(2) - u_r(2))**2 - g*h_r(2))
                    eig_vec(1,2) = 1.d0
                    eig_vec(n_index,2) = sw_priv(2) 
                    eig_vec(t_index,2) = v_l(1)
                    eig_vec(4,2) = alpha(2)
                    eig_vec(n_index,2) = alpha(2)*sw_priv(2)
                    eig_vec(t_index,2) = alpha(2)*v_l(2)
                ! Right bottom state dry
                else if (.not. dry_state_l(1) .and. .not. dry_state_l(2) .and. &
                            .not. dry_state_r(1) .and.       dry_state_r(2)) then

                    sw_priv(5) = u_l(2) + 2.d0 * sqrt(g*(1.d0-r)*h_l(2))
                    alpha(3) = r * g * h_l(2) / ((sw_priv(5) - u_l(2))**2 - g * h_l(2))
                    
                    eig_vec(1,5) = 1.d0
                    eig_vec(n_index,5) = sw_priv(i)
                    eig_vec(t_index,5) = v_r(1)
                    eig_vec(4,5) = alpha(3)
                    eig_vec(n_index,5) = sw_priv(i) * alpha(3)
                    eig_vec(t_index,6) = v_r(2) * alpha(3)
                end if
            end if
            if (inundation_method /= 5) then
                ! Left bottom state dry
                if (.not. dry_state_l(1) .and.       dry_state_l(2) .and.  &
                    .not. dry_state_r(1) .and. .not. dry_state_r(2)) then

                    sw_priv(1) = u_l(1) - sqrt(g*h_l(1))
                    eig_vec(1,1) = 1.d0
                    eig_vec(n_index,1) = sw_priv(1)
                    eig_vec(t_index,1) = v_l(1)
                    eig_vec(4:6,1) = 0.d0

                ! Right bottom state dry
                else if (.not. dry_state_l(1) .and. .not. dry_state_l(2) .and. &
                            .not. dry_state_r(1) .and.       dry_state_r(2)) then

                    sw_priv(6) = u_r(1) + sqrt(g*h_r(1))
                    eig_vec(1,6) = 1.d0
                    eig_vec(n_index,6) = sw_priv(6)
                    eig_vec(t_index,6) = v_r(1)
                    eig_vec(4:6,6) = 0.d0

                end if

            end if
        else
            !DIR$ FORCEINLINE
            call eigen(eigen_h_l, eigen_h_r,                    &
                            eigen_u_l, eigen_u_r,                    &
                            eigen_v_l, eigen_v_r,                    &
                            n_index, t_index,                        &
                            lambda, eig_vec)

        end if

        sw_priv(:) = lambda

        ! ========================num_layers==============================================
        !  Compute flux differences
        ! ======================================================================
        do j=1,2
            layer_index = 3*(j-1)
            flux_r(layer_index+1) = rho(j) * flux_hu_r(j)
            flux_r(layer_index+n_index) = rho(j) * (flux_h_r(j) * flux_u_r(j)**2 + 0.5d0 * g * flux_h_r(j)**2)
            flux_r(layer_index+t_index) = rho(j) * flux_h_r(j) * flux_u_r(j) * flux_v_r(j)
            
            flux_l(layer_index+1) = rho(j) * flux_hu_l(j)
            flux_l(layer_index+n_index) = rho(j) * (flux_h_l(j) * flux_u_l(j)**2 + 0.5d0 * g * flux_h_l(j)**2)
            flux_l(layer_index+t_index) = rho(j) * flux_h_l(j) * flux_u_l(j) * flux_v_l(j)
        enddo
        ! Add extra flux terms
        flux_r(3 + n_index) = flux_r(3 + n_index) + flux_transfer_r
        flux_l(3 + n_index) = flux_l(3 + n_index) + flux_transfer_l
        
        delta = flux_r - flux_l
            
        ! Momentum transfer and bathy terms
        delta(n_index) = delta(n_index) + momentum_transfer(1)
        delta(n_index+3) = delta(n_index+3) + momentum_transfer(2)

        ! ======================================================================
        ! Project jump in fluxes - Use LAPACK's dgesv routine
        !    N - (int) - Number of linear equations (6)
        !    NRHS - (int) - Number of right hand sides (1)
        !    A - (dp(num_layers6,6)) - Coefficient matrix, in this case eig_vec
        !    LDA - (int) - Leading dimension of A (6)
        !    IPIV - (int(N)) - Pivot indices
        !    B - (dp(LDB,NRHS)) - RHS of equations (delta)
        !    LDB - (int) - Leading dimension of B (6)
        !    INFO - (int) - Status of result
        !  Note that the solution (betas) are in delta after the call
        ! ======================================================================
!         A = eig_vec ! We need to do this as the return matrix is modified and
!                     ! we have to use eig_vec again to compute fwaves

        delta1 = delta(1)
        delta2 = delta(2)
        delta3 = delta(3)
        delta4 = delta(4)
        delta5 = delta(5)
        delta6 = delta(6)

        A11 = eig_vec(1,1)
        A12 = eig_vec(1,2)
        A13 = eig_vec(1,3)
        A14 = eig_vec(1,4)
        A15 = eig_vec(1,5)
        A16 = eig_vec(1,6)
        !===
        A21 = eig_vec(2,1)
        A22 = eig_vec(2,2)
        A23 = eig_vec(2,3)
        A24 = eig_vec(2,4)
        A25 = eig_vec(2,5)
        A26 = eig_vec(2,6)
        !===
        A31 = eig_vec(3,1)
        A32 = eig_vec(3,2)
        A33 = eig_vec(3,3)
        A34 = eig_vec(3,4)
        A35 = eig_vec(3,5)
        A36 = eig_vec(3,6)
        !===
        A41 = eig_vec(4,1)
        A42 = eig_vec(4,2)
        A43 = eig_vec(4,3)
        A44 = eig_vec(4,4)
        A45 = eig_vec(4,5)
        A46 = eig_vec(4,6)
        !===
        A51 = eig_vec(5,1)
        A52 = eig_vec(5,2)
        A53 = eig_vec(5,3)
        A54 = eig_vec(5,4)
        A55 = eig_vec(5,5)
        A56 = eig_vec(5,6)
        !===
        A61 = eig_vec(6,1)
        A62 = eig_vec(6,2)
        A63 = eig_vec(6,3)
        A64 = eig_vec(6,4)
        A65 = eig_vec(6,5)
        A66 = eig_vec(6,6)

        ! instead of using a library, use a self-contained code -> vectorization
        !DIR$ FORCEINLINE
        call LU_scalars(A11, A12, A13, A14, A15, A16, &
                A21, A22, A23, A24, A25, A26, &
                A31, A32, A33, A34, A35, A36, &
                A41, A42, A43, A44, A45, A46, &
                A51, A52, A53, A54, A55, A56, &
                A61, A62, A63, A64, A65, A66, &
                6,delta1,delta2,delta3,delta4,delta5,delta6)
        delta(1) = delta1
        delta(2) = delta2
        delta(3) = delta3
        delta(4) = delta4
        delta(5) = delta5
        delta(6) = delta6

            
        !call dgesv(6,1,A,6,pivot,delta,6,info)
!             if (.not.(info == 0)) then
!                 print *, "dry states: ", dry_state_l, dry_state_r
!                 print *, "        left            |             right"
!                 print *, "====================================================="
!                 print "(2d16.8)", h_l(1), h_r(1)
!                 print "(2d16.8)", hu_l(1), hu_r(1)
!                 print "(2d16.8)", hv_l(1), hv_r(1)
!                 print "(2d16.8)", h_l(2), h_r(2)
!                 print "(2d16.8)", hu_l(2), hu_r(2)
!                 print "(2d16.8)", hv_l(2), hv_r(2)
!                 print "(2d16.8)", b_l, b_r
!                 print *,""
!                 print "(a,i2)","In normal solver: ixy=",ixy
!                 print "(a,i3)","  Error solving R beta = delta,",info
!                 print "(a,6d16.8)","  Eigenspeeds: ",s(:,i)
!                 print "(a)","  Eigenvectors:"
!                 do j=1,6
!                     print "(a,6d16.8)","  ",(eig_vec(j,mw),mw=1,6)
!                 enddo
!                 stop
!             endif
        beta = delta

        ! ======================================================================
        ! Compute fwaves
        forall(mw=1:mwaves)
            fw_priv(:,mw) = eig_vec(:,mw) * beta(mw)
        end forall
    end if

end subroutine solve_rpn2


! ==================================
!  Solve the single-layer equations
! ==================================
! subroutine solve_sinlge_layer_rp(layer_index, h_l, h_r,                        &
!                                               hu_l, hu_r,                      &
!                                               hv_l, hv_r,                      &
!                                               u_l, u_r,                        &
!                                               v_l, v_r,                        &
!                                               b_l, b_r,                        &
!                                               fw, sw)

!     use geoclaw_module, only: g => grav
!     use multilayer_module, only: dry_tolerance

!     implicit none

!     ! Input
!     integer, intent(in) :: layer_index
!     real(kind=8), intent(in out), dimension(2) :: h_l, h_r
!     real(kind=8), intent(in out), dimension(2) :: hu_l, hu_r
!     real(kind=8), intent(in out), dimension(2) :: hv_l, hv_r
!     real(kind=8), intent(in out), dimension(2) :: u_r, u_l
!     real(kind=8), intent(in out), dimension(2) :: v_r, v_l
!     real(kind=8), intent(in out) :: b_l, b_r

!     ! Output
!     real(kind=8), intent(out) :: fw(3, 3), sw(3)

!     ! Local storage
!     integer :: m, mw
!     logical :: rare(2)
!     real(kind=8) :: wall(3), h_star, h_star_test, sm(2)
!     real(kind=8) :: phi_l, phi_r, s_l, s_r, u_hat, c_hat, s_roe(2), s_E(2)

!     ! Algorithm parameters
!     integer, parameter :: MAX_ITERATIONS = 1

!     wall = 1.d0
    
!     ! Calculate momentum fluxes
!     phi_l = 0.5d0 * g * h_l(layer_index)**2      &
!                         + h_l(layer_index) * u_l(layer_index)**2
!     phi_r = 0.5d0 * g * h_r(layer_index)**2      &
!                         + h_r(layer_index) * u_r(layer_index)**2
     
!     ! Check for dry state to right
!     if (h_r(layer_index) < dry_tolerance(layer_index)) then
!         call riemanntype(h_l(layer_index), h_l(layer_index),   &
!                          u_l(layer_index),-u_l(layer_index),   &
!                          h_star, sm(1), sm(2), rare(1), rare(2), 1,    &
!                          dry_tolerance(layer_index), g)
!         h_star_test = max(h_l(layer_index), h_star)
!         ! Right state should become ghost values that mirror left for wall problem
!         if (h_star_test + b_l < b_r) then 
!             wall(2:3) = 0.d0
!             h_r(layer_index) = h_l(layer_index)
!             hu_r(layer_index) = -hu_l(layer_index)
!             b_r = b_l
!             phi_r = phi_l
!             u_r(layer_index) = -u_l(layer_index)
!             v_r(layer_index) = v_l(layer_index)
!         else if (h_l(layer_index) + b_l < b_r) then
!             b_r = h_l(layer_index) + b_l
!         endif
!     ! Check for drystate to left, i.e right surface is lower than left topo
!     else if (h_l(layer_index) < dry_tolerance(layer_index)) then 
!         call riemanntype(h_r(layer_index), h_r(layer_index),   &
!                         -u_r(layer_index), u_r(layer_index),   &
!                          h_star, sm(1), sm(2), rare(1), rare(2), 1,    &
!                          dry_tolerance(layer_index), g)
!         h_star_test = max(h_r(layer_index), h_star)
!         ! Left state should become ghost values that mirror right
!         if (h_star_test + b_r < b_l) then  
!            wall(1:2) = 0.d0
!            h_l(layer_index) = h_r(layer_index)
!            hu_l(layer_index) = -hu_r(layer_index)
!            b_l = b_r
!            phi_l = phi_r
!            u_l(layer_index) = -u_r(layer_index)
!            v_l(layer_index) = v_r(layer_index)
!         else if (h_r(layer_index) + b_r < b_l) then
!            b_l = h_r(layer_index) + b_r
!         endif
!     endif

!     ! Determine wave speeds
!     ! 1 wave speed of left state
!     s_l = u_l(layer_index) - sqrt(g * h_l(layer_index))
!     ! 2 wave speed of right state
!     s_r = u_r(layer_index) + sqrt(g * h_r(layer_index))
    
!     ! Roe average
!     u_hat = (sqrt(g * h_l(layer_index)) * u_l(layer_index)     &
!            + sqrt(g * h_r(layer_index)) * u_r(layer_index))    &
!            / (sqrt(g * h_r(layer_index)) + sqrt(g * h_l(layer_index))) 
!     c_hat = sqrt(g * 0.5d0 * (h_r(layer_index) + h_l(layer_index))) 
!     s_roe(1) = u_hat - c_hat ! Roe wave speed 1 wave
!     s_roe(2) = u_hat + c_hat ! Roe wave speed 2 wave
!     s_E(1) = min(s_l, s_roe(1)) ! Eindfeldt speed 1 wave
!     s_E(2) = max(s_r, s_roe(2)) ! Eindfeldt speed 2 wave
    
!     ! Solve Riemann problem
!     call riemann_aug_JCP(MAX_ITERATIONS, 3, 3,                         &
!                          h_l(layer_index), h_r(layer_index),   &
!                          hu_l(layer_index), hu_r(layer_index), &
!                          hv_l(layer_index), hv_r(layer_index), &
!                          b_l, b_r,                                     &
!                          u_l(layer_index), u_r(layer_index),   &
!                          v_l(layer_index), v_r(layer_index),   &
!                          phi_l, phi_r, s_E(1), s_E(2),           &
!                          dry_tolerance(layer_index), g, sw, fw)
    
!     ! Eliminate ghost fluxes for wall
!     do mw=1,3
!         sw(mw) = sw(mw) * wall(mw)
!         do m=1,3
!            fw(m, mw) = fw(m, mw) * wall(mw)
!         enddo
!     enddo

! end subroutine solve_sinlge_layer_rp

subroutine solve_single_layer_rp(layer_index, hL, hR, huL, huR, hvL, hvR, bL, bR, fw, sw)

    use geoclaw_module, only: g => grav
    use multilayer_module, only: dry_tolerance

    implicit none

    ! Input
    integer, intent(in) :: layer_index
    !real(kind=8), intent(in), dimension(2) :: h_l, h_r, hu_l, hu_r, hv_l, hv_r
    !real(kind=8), intent(in) :: b_l, b_r
    real(kind=8), intent(in out) :: hL, hR, huL, huR, hvL, hvR, bL, bR

    ! Output
    real(kind=8), intent(in out) :: fw(3, 3), sw(3)

    ! Locals
    integer :: mw
    real(kind=8) :: uL, uR, vL, vR, pL, pR
    real(kind=8) :: phiL, phiR, wall(3), drytol
    real(kind=8) :: hstar, hstartest, s1m, s2m, rare1, rare2, sL, sR, uhat, chat, sRoe1, sRoe2, sE1, sE2

    ! Parameters (should be anyway)
    integer :: maxiter

    ! Density used in pressure gradient calculation
    ! This needs to be realistic here as compared to air density so we cannot
    ! use what is stored in the multilayer_module unless real densities are 
    ! being used
    ! TODO - fix this limitation 
    real(kind=8), parameter :: rho = 1025.d0

    drytol = dry_tolerance(layer_index)

!     hL = h_l(layer_index)
!     hR = h_r(layer_index)
!     huL = hu_l(layer_index)
!     huR = hu_r(layer_index)
!     hvL = hv_l(layer_index)
!     hvR = hv_r(layer_index)
!     bL = b_l
!     bR = b_r
    pL = 0.d0
    pR = 0.d0

    ! ========================================
    !  Begin Snipped Code From rpn2_geoclaw.f
    ! ========================================
    !check for wet/dry boundary
         if (hR.gt.drytol) then
            uR=huR/hR
            vR=hvR/hR
            phiR = 0.5d0*g*hR**2 + huR**2/hR
         else
            hR = 0.d0
            huR = 0.d0
            hvR = 0.d0
            uR = 0.d0
            vR = 0.d0
            phiR = 0.d0
         endif

         if (hL.gt.drytol) then
            uL=huL/hL
            vL=hvL/hL
            phiL = 0.5d0*g*hL**2 + huL**2/hL
         else
            hL=0.d0
            huL=0.d0
            hvL=0.d0
            uL=0.d0
            vL=0.d0
            phiL = 0.d0
         endif

         wall(1) = 1.d0
         wall(2) = 1.d0
         wall(3) = 1.d0
         if (hR.le.drytol) then
            !DIR$ FORCEINLINE
            call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m,rare1,rare2,1,drytol,g)
            hstartest=max(hL,hstar)
            if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
!                bR=hstartest+bL
               wall(2)=0.d0
               wall(3)=0.d0
               hR=hL
               huR=-huL
               bR=bL
               phiR=phiL
               uR=-uL
               vR=vL
            elseif (hL+bL.lt.bR) then
               bR=hL+bL
            endif
         elseif (hL.le.drytol) then ! right surface is lower than left topo
            !DIR$ FORCEINLINE
            call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m,rare1,rare2,1,drytol,g)
            hstartest=max(hR,hstar)
            if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
!               bL=hstartest+bR
               wall(1)=0.d0
               wall(2)=0.d0
               hL=hR
               huL=-huR
               bL=bR
               phiL=phiR
               uL=-uR
               vL=vR
            elseif (hR+bR.lt.bL) then
               bL=hR+bR
            endif
         endif

         !determine wave speeds
         sL=uL-sqrt(g*hL) ! 1 wave speed of left state
         sR=uR+sqrt(g*hR) ! 2 wave speed of right state

         uhat=(sqrt(g*hL)*uL + sqrt(g*hR)*uR)/(sqrt(g*hR)+sqrt(g*hL)) ! Roe average
         chat=sqrt(g*0.5d0*(hR+hL)) ! Roe average
         sRoe1=uhat-chat ! Roe wave speed 1 wave
         sRoe2=uhat+chat ! Roe wave speed 2 wave

         sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
         sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

         !--------------------end initializing...finally----------
         !solve Riemann problem.

         maxiter = 1

         !DIR$ FORCEINLINE
         call riemann_aug_JCP(maxiter,3,3,hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR, &
                                          vL,vR,phiL,phiR,pL,pR,sE1,sE2,     &
                                          drytol,g,rho,sw,fw)

!         call riemann_ssqfwave(maxiter,meqn,mwaves,hL,hR,huL,huR,
!     &     hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,sw,fw)

!          call riemann_fwave(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,
!     &      bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,sw,fw)

!        !eliminate ghost fluxes for wall
         do mw=1,3
            sw(mw)=sw(mw)*wall(mw)

               fw(1,mw)=fw(1,mw)*wall(mw) 
               fw(2,mw)=fw(2,mw)*wall(mw)
               fw(3,mw)=fw(3,mw)*wall(mw)
         enddo
         
end subroutine solve_single_layer_rp




! ==========================================================================
!  Code from linearized_eigen(), multilayer_module
! ==========================================================================
subroutine eigen(h_l,h_r,u_l,u_r,v_l,v_r,n_index,t_index,s,eig_vec)
    use multilayer_module
    use geoclaw_module, only: grav

    implicit none
    
    ! Input
    double precision, dimension(2), intent(in) :: h_l,h_r,u_l,u_r,v_l,v_r
    integer, intent(in) :: n_index,t_index
    
    ! Output
    double precision, intent(inout) :: s(6),eig_vec(6,6)
        
    ! Local
    double precision :: gamma_l,gamma_r,alpha(4),g
    
    g = grav
        
    ! Calculate relevant quantities
    gamma_l = h_l(2) / h_l(1)
    gamma_r = h_r(2) / h_r(1)

    alpha(1) = 0.5d0*(gamma_l-1.d0+sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
    alpha(2) = 0.5d0*(gamma_l-1.d0-sqrt((gamma_l-1.d0)**2+4.d0*r*gamma_l))
    alpha(3) = 0.5d0*(gamma_r-1.d0-sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))
    alpha(4) = 0.5d0*(gamma_r-1.d0+sqrt((gamma_r-1.d0)**2+4.d0*r*gamma_r))

    s(1) = -sqrt(g*h_l(1)*(1+alpha(1)))
    s(2) = -sqrt(g*h_l(1)*(1+alpha(2)))
    s(3:4) = 0.5d0 * (u_l(:) + u_r(:))
    s(5) = sqrt(g*h_r(1)*(1+alpha(3)))
    s(6) = sqrt(g*h_r(1)*(1+alpha(4)))

    ! Compute eigenspace exactly based on eigenvalues provided
    eig_vec(1,:) = [1.d0,1.d0,0.d0,0.d0,1.d0,1.d0]
    
    eig_vec(n_index,:) = [s(1),s(2),0.d0,0.d0,s(5),s(6)]
    
    eig_vec(t_index,:) = [v_l(1),v_l(1),1.d0,0.d0,v_r(1),v_r(1)]

    eig_vec(4,1:2) = alpha(1:2)
    eig_vec(4,3:4) = 0.d0
    eig_vec(4,5:6) = alpha(3:4)
    
    eig_vec(n_index+3,:) = s * eig_vec(4,:)
    
    eig_vec(t_index+3,1:2) = v_l(2) * alpha(1:2)
    eig_vec(t_index+3,3:4) = [0.d0,1.d0]
    eig_vec(t_index+3,5:6) = v_r(2) * alpha(3:4)

end subroutine eigen
    


