#include "Compilation_control.f90"

#if defined(_SWE_PATCH)

MODULE SWE_PATCH_Solvers
	use Samoa_swe
	use SWE_PATCH
	
	contains

    subroutine compute_updates_simd(transform_matrices, hL, huL, hvL, bL, hR, huR, hvR, bR,upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR, maxWaveSpeed)
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,2,2), intent(in) :: transform_matrices
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE), intent(inout)  :: hL, hR, huL, huR, hvL, hvR, bL, bR
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE), intent(out)    :: upd_hL, upd_hR, upd_huL, upd_huR, upd_hvL, upd_hvR
        real(kind = GRID_SR), intent(inout)                                 :: maxWaveSpeed
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE)             :: uL, uR, vL, vR
        
        !local
        integer                             :: i, j
        real(kind = GRID_SR)                                        :: hstar, s1m, s2m
        logical                                                     :: rare1, rare2
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,3)       :: waveSpeeds
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,3,3) :: fwaves
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,3)       :: wall
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE)     :: delphi
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE)     :: sL, sR, uhat, chat, sRoe1, sRoe2, sE1, sE2

        !DIR$ ASSUME_ALIGNED transform_matrices: 64
        !DIR$ ASSUME_ALIGNED hL:64, hR:64, huL:64, huR:64, hvL:64, hvR:64, bL:64, bR: 64
        !DIR$ ASSUME_ALIGNED upd_hL:64, upd_hR:64, upd_huL:64, upd_huR:64, upd_hvL:64, upd_hvR: 64
        !DIR$ ASSUME_ALIGNED uL:64, uR:64, vL:64, vR: 64
        
        !DIR$ ASSUME_ALIGNED waveSpeeds: 64
        !DIR$ ASSUME_ALIGNED fwaves: 64
        !DIR$ ASSUME_ALIGNED wall: 64
        !DIR$ ASSUME_ALIGNED delphi: 64
        !DIR$ ASSUME_ALIGNED sL:64, sR:64, uhat:64, chat:64, sRoe1:64, sRoe2:64, sE1:64, sE2: 64



        ! *** F-Wave/AugRie solvers *** (based on geoclaw implementation)
        
        !samoa considers bathymetry included in h, the solver doesn't
        hL = hL - bL
        hR = hR - bR

        ! change base so hu/hv become ortogonal/perperdicular to edge
        call apply_transformations_before(transform_matrices, huL, hvL)
        call apply_transformations_before(transform_matrices, huR, hvR)
        
        ! initialize Riemann problem for grid interfaces
        waveSpeeds=0.0_GRID_SR
        fWaves=0.0_GRID_SR
        
        ! check for wet/dry boundary
        where (hR > cfg%dry_tolerance) 
            uR = huR / hR
            vR = hvR / hR
        elsewhere
            hR = 0.0_GRID_SR
            huR = 0.0_GRID_SR
            hvR = 0.0_GRID_SR
            uR = 0.0_GRID_SR
            vR = 0.0_GRID_SR
        end where
        
        where (hL > cfg%dry_tolerance)
            uL = huL / hL
            vL = hvL / hL
        elsewhere
            hL = 0.0_GRID_SR
            huL = 0.0_GRID_SR
            hvL = 0.0_GRID_SR
            uL = 0.0_GRID_SR
            vL = 0.0_GRID_SR
        end where
        
        ! per default there is no wall
        wall = 1.0_GRID_SR
        do i=1,_SWE_PATCH_SOLVER_CHUNK_SIZE
            if (hR(i) <= cfg%dry_tolerance) then
#               if defined(_SINGLE_PRECISION)
                    call riemanntype_sp(hL(i), hL(i), uL(i), -uL(i), hstar, s1m, s2m, rare1, rare2, 1, cfg%dry_tolerance, g)
#               elif defined(_DOUBLE_PRECISION)
                    call riemanntype(hL(i), hL(i), uL(i), -uL(i), hstar, s1m, s2m, rare1, rare2, 1, cfg%dry_tolerance, g)
#               endif
                hstar = max(hL(i), hstar)
                if (hstar + bL(i) < bR(i)) then !right state should become ghost values that mirror left for wall problem
                    wall(i,2) = 0.0_GRID_SR
                    wall(i,3) = 0.0_GRID_SR
                    hR(i) = hL(i)
                    huR(i) = - huL(i)
                    bR(i) = bL(i)
                    uR(i) = -uL(i)
                    vR(i) = vL(i)
                else if (hL(i) + bL(i) < bR(i)) then
                    bR(i) = hL(i)+bL(i)
                end if
            else if (hL(i) <= cfg%dry_tolerance) then ! right surface is lower than left topo
#               if defined(_SINGLE_PRECISION)
                    call riemanntype_sp(hR(i), hR(i), -uR(i), uR(i), hstar, s1m, s2m, rare1, rare2, 1, cfg%dry_tolerance, g)
#               elif defined(_DOUBLE_PRECISION)
                    call riemanntype(hR(i), hR(i), -uR(i), uR(i), hstar, s1m, s2m, rare1, rare2, 1, cfg%dry_tolerance, g)
#               endif
                hstar = max (hR(i), hstar)
                if (hstar + bR(i) < bL(i)) then !left state should become ghost values that mirror right
                    wall(i,1) = 0.0_GRID_SR
                    wall(i,2) = 0.0_GRID_SR
                    hL(i) = hR(i)
                    huL(i) = -huR(i)
                    bL(i) = bR(i)
                    uL(i) = -uR(i)
                    vL(i) = vR(i)
                else if (hR(i) + bR(i) < bL(i)) then
                    bL(i) = hR(i) + bR(i)
                end if
            end if
        end do
        
        ! BUGFIX:
        ! Problem: loss of significance may occur in phiR-phiL, causing divergence of the steady state.
        ! Action:  Compute delphi=phiR-phiL explicitly. delphi is arithmetically equivalent to phiR-phiL, but with far smaller numerical loss.
        delphi = (huR - huL)*(uL + uR) - uL*uR*(hR-hL) + (0.5_GRID_SR * g *(bR +hR - bL - hL)*(hR + hL)) - 0.5_GRID_SR*g*(hR + hL)*(bR - bL)
        
        ! determine wave speeds
        sL=uL-sqrt(g*hL) ! 1 wave speed of left state
        sR=uR+sqrt(g*hR) ! 2 wave speed of right state

        uhat=(sqrt(g*hL)*uL + sqrt(g*hR)*uR)/(sqrt(g*hR)+sqrt(g*hL)) ! Roe average
        chat=sqrt(g*0.5_GRID_SR*(hR+hL)) ! Roe average
        sRoe1=uhat-chat ! Roe wave speed 1 wave
        sRoe2=uhat+chat ! Roe wave speed 2 wave

        sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
        sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave
        
        !*******************
        !* call the solver *
        !*******************

#       if defined(_SWE_FWAVE)
            call riemann_fwave_simd(hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vL,vR,delphi,sE1,sE2,waveSpeeds,fWaves)
#       elif defined(_SWE_AUG_RIEMANN)
            call riemann_augrie_simd(1,hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vL,vR,delphi,sE1,sE2,waveSpeeds,fWaves)
#       endif
        !*****************
        !* end of solver *
        !*****************

        ! eliminate ghost fluxes for wall
        do i=1, 3 !waveNumber
            waveSpeeds(:,i) = waveSpeeds(:,i) * wall(:,i)
            do j=1,3 !equationNumber
                fwaves(:,j,i) = fwaves(:,j,i) * wall(:,i)
            end do
        end do
        
        ! compute net updates
        upd_hL = 0.0_GRID_SR
        upd_huL = 0.0_GRID_SR
        upd_hvL = 0.0_GRID_SR
        upd_hR = 0.0_GRID_SR
        upd_huR = 0.0_GRID_SR
        upd_hvR = 0.0_GRID_SR
        
        do i=1,3 ! waveNumber
            where (waveSpeeds(:,i) < 0)
                upd_hL = upd_hL + fwaves(:,1,i)
                upd_huL = upd_huL + fwaves(:,2,i)
                upd_hvL = upd_hvL + fwaves(:,3,i)
            elsewhere (waveSpeeds(:,i) > 0)
                upd_hR = upd_hR + fwaves(:,1,i)
                upd_huR = upd_huR + fwaves(:,2,i)
                upd_hvR = upd_hvR + fwaves(:,3,i)
            elsewhere
                upd_hL = upd_hL + 0.5_GRID_SR * fwaves(:,1,i)
                upd_huL = upd_huL + 0.5_GRID_SR * fwaves(:,2,i)
                upd_hvL = upd_hvL + 0.5_GRID_SR * fwaves(:,3,i)
                upd_hR = upd_hR + 0.5_GRID_SR * fwaves(:,1,i)
                upd_huR = upd_huR + 0.5_GRID_SR * fwaves(:,2,i)
                upd_hvR = upd_hvR + 0.5_GRID_SR * fwaves(:,3,i)
            end where
        end do
        
        ! no net updates in DryDry case -> a non-vectorized solver would have stopped the computation
        ! right after checking wet/dry boundary
        where (hL < cfg%dry_tolerance .and. hR < cfg%dry_tolerance)
            waveSpeeds(:,1) = 0.0_GRID_SR
            waveSpeeds(:,2) = 0.0_GRID_SR
            waveSpeeds(:,3) = 0.0_GRID_SR
            upd_hL = 0.0_GRID_SR
            upd_huL = 0.0_GRID_SR
            upd_hvL = 0.0_GRID_SR
            upd_hR = 0.0_GRID_SR
            upd_huR = 0.0_GRID_SR
            upd_hvR = 0.0_GRID_SR
        end where

        ! compute maximum wave speed
        maxWaveSpeed = max(maxWaveSpeed, maxVal(abs(waveSpeeds(:,:))))
        
        
        ! inverse transformations
        call apply_transformations_after(transform_matrices, upd_huL, upd_hvL)
        call apply_transformations_after(transform_matrices, upd_huR, upd_hvR)
        
    end subroutine
	
	
    subroutine riemann_fwave_simd(hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vL,vR,delphi,s1,s2,sw,fw)
        
        ! --> implementation with vectorization on patches
      
        ! solve shallow water equations given single left and right states
        ! solution has two waves.
        ! flux - source is decomposed.
        
        implicit none

        !input
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE), intent(inout) :: hL,hR,huL,huR,bL,bR,uL,uR,delphi,s1,s2
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE), intent(inout) :: hvL,hvR,vL,vR

        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,3), intent(inout) ::  sw
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,3,3), intent(inout) ::  fw

        !local
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE) :: delh, delhu, delb, deldelphi, delphidecomp, beta1, beta2


        !DIR$ ASSUME_ALIGNED hL:64, hR:64, huL:64, huR:64, bL:64, bR:64, uL:64, uR:64, delphi:64, s1:64, s2:64, hvL:64, hvR:64, vL:64, vR:64
        !DIR$ ASSUME_ALIGNED sw:64, fw:64
        !DIR$ ASSUME_ALIGNED delh:64, delhu:64, delb:64, deldelphi:64, delphidecomp:64, beta1:64, beta2:64
        
        !determine del vectors
        delh = hR-hL
        delhu = huR-huL
        delb = bR-bL

        deldelphi = -g*0.5d0*(hR+hL)*delb
        delphidecomp = delphi - deldelphi

        !flux decomposition
        beta1 = (s2*delhu - delphidecomp)/(s2-s1)
        beta2 = (delphidecomp - s1*delhu)/(s2-s1)

        sw(:,1)=s1
        sw(:,2)=0.5d0*(s1+s2)
        sw(:,3)=s2
        ! 1st nonlinear wave
        fw(:,1,1) = beta1
        fw(:,2,1) = beta1*s1
        fw(:,3,1) = beta1*vL
        ! 2nd nonlinear wave
        fw(:,1,3) = beta2
        fw(:,2,3) = beta2*s2
        fw(:,3,3) = beta2*vR
        ! advection of transverse wave
        fw(:,1,2) = 0.d0
        fw(:,2,2) = 0.d0
        fw(:,3,2) = huR*vR - huL*vL -fw(:,3,1)-fw(:,3,3)
    end subroutine
	
	subroutine riemann_augrie_simd(maxiter,hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vL,vR,delphi,sE1,sE2,sw,fw)
        
        ! --> implementation with vectorization on patches
      
        ! solve shallow water equations given single left and right states
        ! This solver is described in J. Comput. Phys. (6): 3089-3113, March 2008
        ! Augmented Riemann Solvers for the Shallow Equations with Steady States and Inundation

        ! To use the original solver call with maxiter=1.

        ! This solver allows iteration when maxiter > 1. The iteration seems to help with
        ! instabilities that arise (with any solver) as flow becomes transcritical over variable topo
        ! due to loss of hyperbolicity.



        implicit none

        !input
        integer maxiter
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,3,3), intent(inout) :: fw
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,3), intent(inout) :: sw
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE), intent(inout) :: hL,hR,huL,huR,bL,bR,uL,uR,delphi,sE1,sE2
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE), intent(inout) :: hvL,hvR,vL,vR

        !local
        integer, parameter :: mwaves = 3, meqn = 3
        integer :: m,mw,k,iter, i
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,3,3) :: A, r
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,3) :: lambda, del, beta

        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE) :: delh,delhu,delb,delnorm
        !real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE) :: rare1st,rare2st,sdelta,raremin,raremax
        real(kind = GRID_SR)                                            :: criticaltol,convergencetol,raretol
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE) :: s1s2bar,s1s2tilde,hbar,hLstar,hRstar,hustar
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE) :: huRstar,huLstar,uRstar,uLstar,hstarHLL
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE) :: deldelh,deldelphi
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE) :: s1m,s2m,hm
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE) :: det1,det2,det3,determinant

        logical, dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE) :: rare1,rare2,sonic
        !logical, dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE) :: rarecorrector,rarecorrectortest
        
        !DIR$ ASSUME_ALIGNED fw:64, sw:64
        !DIR$ ASSUME_ALIGNED hL:64, hR:64, huL:64, huR:64, bL:64, bR:64, uL:64, uR:64, delphi:64, sE1:64,sE2:64
        !DIR$ ASSUME_ALIGNED hL:64, hvL:64, hvR:64, vL:64, vR:64
        !DIR$ ASSUME_ALIGNED A:64, r:64, lambda:64, del:64, beta:64
        !DIR$ ASSUME_ALIGNED delh:64, delhu:64, delb:64, delnorm:64
        !---DIR$ ASSUME_ALIGNED rare1st:64, rare2st:64, sdelta:64, raremin:64, raremax:64
        !DIR$ ASSUME_ALIGNED s1s2bar:64, s1s2tilde:64, hbar:64, hLstar:64, hRstar:64, hustar:64
        !DIR$ ASSUME_ALIGNED huRstar:64, huLstar:64, uRstar:64, uLstar:64, hstarHLL:64
        !DIR$ ASSUME_ALIGNED deldelh:64, deldelphi:64
        !DIR$ ASSUME_ALIGNED s1m:64, s2m:64, hm:64
        !DIR$ ASSUME_ALIGNED det1:64, det2:64, det3:64, determinant:64
        !DIR$ ASSUME_ALIGNED rare1:64, rare2:64, sonic:64
        !---DIR$ ASSUME_ALIGNED rarecorrector:64, rarecorrectortest:64


        !determine del vectors
        delh = hR-hL
        delhu = huR-huL
        delb = bR-bL
        delnorm = delh**2 + delphi**2

        do i=1,_SWE_PATCH_SOLVER_CHUNK_SIZE
#           if defined(_SINGLE_PRECISION)
                call riemanntype_sp(hL(i), hR(i), uL(i), uR(i), hm(i), s1m(i), s2m(i), rare1(i), rare2(i), 1, cfg%dry_tolerance, g)
#           elif defined(_DOUBLE_PRECISION)
                call riemanntype(hL(i), hR(i), uL(i), uR(i), hm(i), s1m(i), s2m(i), rare1(i), rare2(i), 1, cfg%dry_tolerance, g)
#           endif
        end do

        lambda(:,1)= min(sE1,s2m) !Modified Einfeldt speed
        lambda(:,3)= max(sE2,s1m) !Modified Eindfeldt speed
        sE1=lambda(:,1)
        sE2=lambda(:,3)
        hstarHLL = max((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.0_GRID_SR) ! middle state in an HLL solve

        !determine the middle entropy corrector wave------------------------
        ! commenting this because rarecorrectortest and rarecorrector are always false!
!         rarecorrectortest=.false.
!         rarecorrector=.false.
!         where (rarecorrectortest) 
!             sdelta=lambda(:,3)-lambda(:,1)
!             raremin = 0.5_GRID_SR
!             raremax = 0.9_GRID_SR
!             where (rare1.and.sE1*s1m.lt.0.d0) raremin=0.2_GRID_SR
!             where (rare2.and.sE2*s2m.lt.0.d0) raremin=0.2_GRID_SR
!             where (rare1.or.rare2)
!                 !see which rarefaction is larger
!                 rare1st=3.0_GRID_SR*(sqrt(g*hL)-sqrt(g*hm))
!                 rare2st=3.0_GRID_SR*(sqrt(g*hR)-sqrt(g*hm))
!                 where (max(rare1st,rare2st).gt.raremin*sdelta .and. max(rare1st,rare2st).lt.raremax*sdelta)
!                     rarecorrector=.true.
!                 where (rare1st.gt.rare2st)
!                     lambda(:,2)=s1m
!                 elsewhere (rare2st.gt.rare1st)
!                     lambda(:,2)=s2m
!                 elsewhere
!                     lambda(:,2)=0.5_GRID_SR*(s1m+s2m)
!                 end where
!                 end where
!             end where
!             where (hstarHLL.lt.min(hL,hR)*0.2_GRID_SR) rarecorrector=.false.
!         end where
        
        where (abs(lambda(:,2)) .lt. 1.e-20_GRID_SR) lambda(:,2) = 0.0_GRID_SR
        
        do mw=1,mwaves
            r(:,1,mw)=1.0_GRID_SR
            r(:,2,mw)=lambda(:,mw)
            r(:,3,mw)=(lambda(:,mw))**2
        enddo

        ! commenting this because rarecorrectortest and rarecorrector are always false!
!         where (.not.rarecorrector)
            lambda(:,2) = 0.5_GRID_SR*(lambda(:,1)+lambda(:,3))
            lambda(:,2) = max(min(0.5_GRID_SR*(s1m+s2m),sE2),sE1)
            where (abs(lambda(:,2)) .lt. 1.d-20) lambda(:,2) = 0.0_GRID_SR
            r(:,1,2)=0.0_GRID_SR
            r(:,2,2)=0.0_GRID_SR
            r(:,3,2)=1.0_GRID_SR
!         end where
        !---------------------------------------------------

        !determine the steady state wave -------------------
        criticaltol = 1.e-6_GRID_SR
        deldelh = -delb
        deldelphi = -g*0.5_GRID_SR*(hR+hL)*delb

        !determine a few quanitites needed for steady state wave if iterated
        hLstar=hL
        hRstar=hR
        uLstar=uL
        uRstar=uR
        huLstar=uLstar*hLstar
        huRstar=uRstar*hRstar

        !iterate to better determine the steady state wave
        convergencetol=1.e-6_GRID_SR
        do iter=1,maxiter
            !determine steady state wave (this will be subtracted from the delta vectors)
            ! commenting this because rarecorrectortest and rarecorrector are always false!
!             where (min(hLstar,hRstar).lt.cfg%dry_tolerance.and.rarecorrector)
!                 rarecorrector=.false.
!                 hLstar=hL
!                 hRstar=hR
!                 uLstar=uL
!                 uRstar=uR
!                 huLstar=uLstar*hLstar
!                 huRstar=uRstar*hRstar
!                 lambda(:,2) = 0.5_GRID_SR*(lambda(:,1)+lambda(:,3))
!                 lambda(:,2) = max(min(0.5_GRID_SR*(s1m+s2m),sE2),sE1)
!                 where (abs(lambda(:,2)) .lt. 1.0_GRID_SR) lambda(:,2) = 0.0_GRID_SR
!                 r(:,1,2)=0.0_GRID_SR
!                 r(:,2,2)=0.0_GRID_SR
!                 r(:,3,2)=1.0_GRID_SR
!             end where

            hbar =  max(0.5_GRID_SR*(hLstar+hRstar),0.0_GRID_SR)
            s1s2bar = 0.25_GRID_SR*(uLstar+uRstar)**2 - g*hbar
            s1s2tilde= max(0.0_GRID_SR,uLstar*uRstar) - g*hbar

            !find if sonic problem
            sonic = abs(s1s2bar).le.criticaltol .OR. s1s2bar*s1s2tilde.le.criticaltol .OR. s1s2bar*sE1*sE2.le.criticaltol .OR. min(abs(sE1),abs(sE2)).lt.criticaltol
            sonic = sonic .OR. sE1.lt.0.0_GRID_SR.and.s1m.gt.0.0_GRID_SR .OR. sE2.gt.0.0_GRID_SR.and.s2m.lt.0.0_GRID_SR .OR. (uL+sqrt(g*hL))*(uR+sqrt(g*hR)).lt.0.0_GRID_SR .OR. (uL-sqrt(g*hL))*(uR-sqrt(g*hR)).lt.0.0_GRID_SR

            !find jump in h, deldelh
            where (sonic) 
                deldelh =  -delb
            elsewhere
                deldelh = delb*g*hbar/s1s2bar
            end where


            
            !find bounds in case of critical state resonance, or negative states
            where (sE1.lt.-criticaltol.and.sE2.gt.criticaltol) 
                deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE2)
                deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE1)
            elsewhere (sE1.ge.criticaltol)
                deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE1)
                deldelh = max(deldelh,-hL)
            elsewhere (sE2.le.-criticaltol)
                deldelh = min(deldelh,hR)
                deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE2)
            end where

            !find jump in phi, deldelphi
            where (sonic) 
                deldelphi = -g*hbar*delb
            elsewhere
                deldelphi = -delb*g*hbar*s1s2tilde/s1s2bar
            end where
         
            !find bounds in case of critical state resonance, or negative states
            deldelphi=min(deldelphi,g*max(-hLstar*delb,-hRstar*delb))
            deldelphi=max(deldelphi,g*min(-hLstar*delb,-hRstar*delb))

            del(:,1)=delh-deldelh
            del(:,2)=delhu
            del(:,3)=delphi-deldelphi

            !Determine determinant of eigenvector matrix========
            det1=r(:,1,1)*(r(:,2,2)*r(:,3,3)-r(:,2,3)*r(:,3,2))
            det2=r(:,1,2)*(r(:,2,1)*r(:,3,3)-r(:,2,3)*r(:,3,1))
            det3=r(:,1,3)*(r(:,2,1)*r(:,3,2)-r(:,2,2)*r(:,3,1))
            determinant=det1-det2+det3
            


            !solve for beta(k) using Cramers Rule=================
            do k=1,3
                A = r
                A(:,:,k) = del(:,:)

                det1=A(:,1,1)*(A(:,2,2)*A(:,3,3)-A(:,2,3)*A(:,3,2))
                det2=A(:,1,2)*(A(:,2,1)*A(:,3,3)-A(:,2,3)*A(:,3,1))
                det3=A(:,1,3)*(A(:,2,1)*A(:,3,2)-A(:,2,2)*A(:,3,1))
                beta(:,k)=(det1-det2+det3)/determinant
            enddo

            !exit if things aren't changing --> not anymore, with vectorization we can't leave early
            !if (abs(del(1)**2+del(3)**2-delnorm).lt.convergencetol) exit
            
            delnorm = del(:,1)**2+del(:,3)**2
            !find new states qLstar and qRstar on either side of interface
            hLstar=hL
            hRstar=hR
            uLstar=uL
            uRstar=uR
            huLstar=uLstar*hLstar
            huRstar=uRstar*hRstar
            do mw=1,mwaves
                where (lambda(:,mw).lt.0.0_GRID_SR)
                hLstar= hLstar + beta(:,mw)*r(:,1,mw)
                huLstar= huLstar + beta(:,mw)*r(:,2,mw)
                end where
            enddo
            do mw=mwaves,1,-1
                where (lambda(:,mw).gt.0.0_GRID_SR)
                hRstar= hRstar - beta(:,mw)*r(:,1,mw)
                huRstar= huRstar - beta(:,mw)*r(:,2,mw)
                end where
            enddo

            where (hLstar.gt.cfg%dry_tolerance) 
                uLstar=huLstar/hLstar
            else where
                hLstar=max(hLstar,0.0_GRID_SR)
                uLstar=0.0_GRID_SR
            end where
            where (hRstar.gt.cfg%dry_tolerance) 
                uRstar=huRstar/hRstar
            else where
                hRstar=max(hRstar,0.0_GRID_SR)
                uRstar=0.0_GRID_SR
            end where

        enddo ! end iteration on Riemann problem

        do mw=1,mwaves
            sw(:,mw)=lambda(:,mw)
            fw(:,1,mw)=beta(:,mw)*r(:,2,mw)
            fw(:,2,mw)=beta(:,mw)*r(:,3,mw)
            fw(:,3,mw)=beta(:,mw)*r(:,2,mw)
        enddo
        !find transverse components (ie huv jumps).
        fw(:,3,1)=fw(:,3,1)*vL
        fw(:,3,3)=fw(:,3,3)*vR
        fw(:,3,2)= hR*uR*vR - hL*uL*vL - fw(:,3,1)- fw(:,3,3)
    end subroutine
	
	subroutine compute_updates_hlle_simd(transform_matrices, hL, huL, hvL, bL, hR, huR, hvR, bR, upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR, maxWaveSpeed)
		real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,2,2),intent(in)	:: transform_matrices
		real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE), intent(inout)	:: hL, hR, huL, huR, hvL, hvR, bL, bR
		real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE), intent(out)	:: upd_hL, upd_hR, upd_huL, upd_huR, upd_hvL, upd_hvR
		real(kind = GRID_SR), intent(inout)									:: maxWaveSpeed
		
		!local
		integer														:: i, j
		real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE)		:: uL, uR
		real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE)		:: sqrt_hL, sqrt_hR, sqrt_ghL, sqrt_ghR
		real(kind = GRID_SR)										:: half_g, sqrt_g
		integer(kind = BYTE), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE)		:: wetDryState
		
		real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,2)		:: characteristicSpeeds, roeSpeeds, extEinfeldtSpeeds, steadyStateWave
		real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE)		:: hRoe, uRoe, sqrt_g_hRoe, hLLMiddleHeight, inverseDiff
		real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,3)		:: eigenValues, rightHandSide, beta
		real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,3,3)	:: eigenVectors
		real(kind = GRID_SR), parameter								:: r_eps = 1e-7
		
		real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,2,3)	:: fWaves
		real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,3)		:: waveSpeeds
		
		enum, bind(c) !constants to classify wet-dry-state of pairs of cells
			enumerator :: DryDry = 0
			enumerator :: WetWet = 1
			enumerator :: WetDryInundation = 2
			enumerator :: WetDryWall = 3
			enumerator :: WetDryWallInundation = 4
			enumerator :: DryWetInundation = 5
			enumerator :: DryWetWall = 6
			enumerator :: DryWetWallInundation = 7
		end enum
		
		enum, bind(c) ! constants to classify Riemann state of a pair of cells:
			enumerator :: DrySingleRarefaction = 0
			enumerator :: SingleRarefactionDry = 1
			enumerator :: ShockShock = 2
			enumerator :: ShockRarefaction = 3
			enumerator :: RarefactionShock = 4
			enumerator :: RarefactionRarefaction = 5
		end enum

		!DIR$ ASSUME_ALIGNED hL: 64
		!DIR$ ASSUME_ALIGNED hR: 64
		!DIR$ ASSUME_ALIGNED huL: 64
		!DIR$ ASSUME_ALIGNED huR: 64
		!DIR$ ASSUME_ALIGNED hvL: 64
		!DIR$ ASSUME_ALIGNED hvR: 64
		!DIR$ ASSUME_ALIGNED bL: 64
		!DIR$ ASSUME_ALIGNED bR: 64

		!DIR$ ASSUME_ALIGNED transform_matrices: 64
		!DIR$ ASSUME_ALIGNED upd_hL: 64
		!DIR$ ASSUME_ALIGNED upd_hR: 64
		!DIR$ ASSUME_ALIGNED upd_huL: 64
		!DIR$ ASSUME_ALIGNED upd_huR: 64
		!DIR$ ASSUME_ALIGNED upd_hvL: 64
		!DIR$ ASSUME_ALIGNED upd_hvR: 64
		!DIR$ ASSUME_ALIGNED uL: 64
		!DIR$ ASSUME_ALIGNED uR: 64
		
		!DIR$ ASSUME_ALIGNED sqrt_hL: 64
		!DIR$ ASSUME_ALIGNED sqrt_hR: 64
		!DIR$ ASSUME_ALIGNED sqrt_ghL: 64
		!DIR$ ASSUME_ALIGNED sqrt_ghR: 64
		!DIR$ ASSUME_ALIGNED wetDryState: 64
		
		!DIR$ ASSUME_ALIGNED characteristicSpeeds: 64
		!DIR$ ASSUME_ALIGNED roeSpeeds: 64
		!DIR$ ASSUME_ALIGNED extEinfeldtSpeeds: 64
		!DIR$ ASSUME_ALIGNED steadyStateWave: 64
		!DIR$ ASSUME_ALIGNED hRoe: 64
		!DIR$ ASSUME_ALIGNED uRoe: 64
		!DIR$ ASSUME_ALIGNED sqrt_g_hRoe: 64
		!DIR$ ASSUME_ALIGNED hLLMiddleHeight: 64
		!DIR$ ASSUME_ALIGNED inverseDiff: 64
		!DIR$ ASSUME_ALIGNED eigenValues: 64
		!DIR$ ASSUME_ALIGNED rightHandSide: 64
		!DIR$ ASSUME_ALIGNED eigenVectors: 64
		!DIR$ ASSUME_ALIGNED beta: 64
		
		!DIR$ ASSUME_ALIGNED fWaves: 64
		!DIR$ ASSUME_ALIGNED waveSpeeds: 64

		!samoa considers bathymetry included in h, the solver doesn't
		hL = hL - bL
		hR = hR - bR

		! change base so hu/hv become ortogonal/perperdicular to edge
		call apply_transformations_before(transform_matrices, huL, hvL)
		call apply_transformations_before(transform_matrices, huR, hvR)
			
			
		! initialize uL and uR (will be computed later, if h>dry_tolerance	
		uL = 0.0_GRID_SR
		uR = 0.0_GRID_SR
		
		!!!!! actual solver
		! reset net updates
		upd_hL = 0.0_GRID_SR
		upd_hR = 0.0_GRID_SR
		upd_huL = 0.0_GRID_SR
		upd_huR = 0.0_GRID_SR
		upd_hvL = 0.0_GRID_SR
		upd_hvR = 0.0_GRID_SR
		
		! declare variables which are used over and over again
		half_g = 0.5_GRID_SR * g
		sqrt_g = sqrt(g)
		
		!************************************************************************
		!* Determine Wet Dry State Begin
		!* (determine the wet/dry state and compute local variables correspondly)
		!*************************************************************************
		
		! compute speeds or set them to zero (dry cells)
		where (hL >= cfg%dry_tolerance)
			uL = huL / hL
		elsewhere
			bL = bL + hL
			hL = 0.0_GRID_SR
			huL = 0.0_GRID_SR
			!uL = 0.0_GRID_SR  !is already zero
		end where
		
		where (hR >= cfg%dry_tolerance)
			uR = huR / hR
		elsewhere
			bR = bR + hR
			hR = 0.0_GRID_SR
			huR = 0.0_GRID_SR
			!uR = 0_GRID_SR  !is already zero
		end where				
		
		! MB: determine wet/dry-state - try to start with most frequent case
		where (hL >= cfg%dry_tolerance)
			where (hR >= cfg%dry_tolerance)
				! simple wet/wet case - expected as most frequently executed branch
				wetDryState = WetWet
			elsewhere ! hL >= cfg%dry_tolerance .and. hR < cfg%dry_tolerance
				! we have a shoreline: left cell wet, right cell dry
				! => check for simple inundation problems
				where (hL + bL > bR)
					! => dry cell lies lower than wet cell
					wetDryState = WetDryInundation
				elsewhere ! hL >= cfg%dry_tolerance .and. hR < cfg%dry_tolerance .and. hL + bL <= bR
					! => dry cell (right) lies higher than the wet cell (left)
					! Removed middle state computation -> assuming momentum is never high enough to overcome boundary on its own
					hR = hL
					uR = -uL
					huR = -huL
					bR = 0.0_GRID_SR
					bL = 0.0_GRID_SR
					wetDryState = WetDryWall
				end where
			end where
		elsewhere ! hL < cfg%dry_tolerance
			where (hR >= cfg%dry_tolerance)
				! we have a shoreline: left cell dry, right cell wet
				! => check for simple inundation problems
				where (hR + bR > bL)
					! => dry cell lies lower than wet cell
					wetDryState = DryWetInundation
				elsewhere ! hL < cfg%dry_tolerance .and. hR >= cfg%dry_tolerance .and. hR + bR <= bL
					! => dry cell (left) lies higher than the wet cell (right)
					! Removed middle state computation -> assuming momentum is never high enough to overcome boundary on its own
					hL = hR
					uL = -uR
					huL = -huR
					bL = 0.0_GRID_SR
					bR = 0.0_GRID_SR
					wetDryState = DryWetWall
				end where
			elsewhere ! hL < cfg%dry_tolerance .and. hR < cfg%dry_tolerance
				wetDryState = DryDry
				! nothing to do for dry/dry case, all netUpdates and maxWaveSpeed are 0
			end where
		end where

		!************************************************************************
		!* Determine Wet Dry State End
		!*************************************************************************

		! precompute some terms which are fixed during
		! the computation after some specific point
		sqrt_hL = sqrt(hL)
		sqrt_hR = sqrt(hR)
		sqrt_ghR = sqrt_g * sqrt_hR
		sqrt_ghL = sqrt_g * sqrt_hL
		
		!************************************************************************
		!* Compute Wave Decomposition Begin
		!************************************************************************
		
		! compute eigenvalues of the jacobian matrices in states Q_{i-1} and Q_{i} (char. speeds)
		characteristicSpeeds(:,1) = uL - sqrt_ghL
		characteristicSpeeds(:,2) = uR + sqrt_ghR
		
		! compute "Roe Speeds"
		hRoe = 0.5_GRID_SR * (hR + hL)
		uRoe = (uL * sqrt_hL + uR * sqrt_hR) / (sqrt_hL + sqrt_hR)
		
		! optimization for dumb compilers
		sqrt_g_hRoe = sqrt_g * sqrt(hRoe)
		roeSpeeds(:,1) = uRoe - sqrt_g_hRoe
		roeSpeeds(:,2) = uRoe + sqrt_g_hRoe
		
		! middle state computation removed, using basic Einfeldt speeds
		
		extEinfeldtSpeeds = 0_GRID_SR
		where (wetDryState == WetWet .or. wetDryState == WetDryWall .or. wetDryState == DryWetWall)
			extEinfeldtSpeeds(:,1) = min(characteristicSpeeds(:,1), roeSpeeds(:,1))
			extEinfeldtSpeeds(:,2) = max(characteristicSpeeds(:,2), roeSpeeds(:,2))
		elsewhere (hL < cfg%dry_tolerance) ! MB: !!! cases DryWetInundation, DryWetWallInundation
			!ignore undefined speeds
			extEinfeldtSpeeds(:,1) = roeSpeeds(:,1)
			extEinfeldtSpeeds(:,2) = max(characteristicSpeeds(:,2), roeSpeeds(:,2))
		elsewhere (hR < cfg%dry_tolerance) ! MB: !!! cases WetDryInundation, WetDryWallInuncation
			! ignore undefined speeds		
			extEinfeldtSpeeds(:,1) = min(characteristicSpeeds(:,1), roeSpeeds(:,1))
			extEinfeldtSpeeds(:,2) = roeSpeeds(:,2)
		end where
		
		! HLL middle state
		!  \cite[theorem 3.1]{george2006finite}, \cite[ch. 4.1]{george2008augmented}
		hLLMiddleHeight = max(0.0_GRID_SR, (huL - huR + extEinfeldtSpeeds(:,2) * hR - extEinfeldtSpeeds(:,1) * hL) / (extEinfeldtSpeeds(:,2) - extEinfeldtSpeeds(:,1)) )
		
		! define eigenvalues
		eigenvalues(:,1) = extEinfeldtSpeeds(:,1)
		eigenvalues(:,2) = 0.5_GRID_SR * (extEinfeldtSpeeds(:,1) + extEinfeldtSpeeds(:,2))
		eigenvalues(:,3) = extEinfeldtSpeeds(:,2)
		
		! define eigenvectors
		! MB: no longer used as system matrix
		!     -> but still used to compute f-waves
		eigenVectors(:,1,1) = 1.0_GRID_SR
		eigenVectors(:,2,1) = 0.0_GRID_SR
		eigenVectors(:,3,1) = 1.0_GRID_SR
		
		eigenVectors(:,1,2) = eigenValues(:,1)
		eigenVectors(:,2,2) = 0.0_GRID_SR
		eigenVectors(:,3,2) = eigenValues(:,3)
		
		eigenVectors(:,1,3) = eigenValues(:,1) * eigenValues(:,1)
		eigenVectors(:,2,3) = 1.0_GRID_SR
		eigenVectors(:,3,3) = eigenValues(:,3) * eigenValues(:,3)
		
		! compute the jump in state
		rightHandSide(:,1) = hR - hL
		rightHandSide(:,2) = huR - huL
		rightHandSide(:,3) = (huR * uR + half_g * hR*hR) - (huL * uL + half_g * hL*hL)
		
		! compute steady state wave
		steadyStateWave(:,1) = -(bR - bL)
		steadyStateWave(:,2) = -half_g * (hL + hR) * (bR - bL)
		
		! preserve depth-positivity
		!  \cite[ch. 6.5.2]{george2006finite}, \cite[ch. 4.2.3]{george2008augmented}
		where (eigenvalues(:,1) < -r_eps .and. eigenvalues(:,3) > r_eps)
			! subsonic
			steadyStateWave(:,1) = max(steadyStateWave(:,1), hLLMiddleHeight * (eigenValues(:,3) - eigenValues(:,1)) / eigenValues(:,1))
			steadyStateWave(:,1) = min(steadyStateWave(:,1), hLLMiddleHeight * (eigenValues(:,3) - eigenValues(:,1)) / eigenValues(:,3))
		elsewhere (eigenValues(:,1) > r_eps)
			! supersonic right TODO: motivation?
			steadyStateWave(:,1) = max(steadyStateWave(:,1), -hL)
			steadyStateWave(:,1) = min(steadyStateWave(:,1), hLLMiddleHeight * (eigenValues(:,3) - eigenValues(:,1)) / eigenValues(:,1))
		elsewhere (eigenvalues(:,3) < -r_eps)
			! supersonic left TODO: motivation?
			steadyStateWave(:,1) = max(steadyStateWave(:,1), hLLMiddleHeight * (eigenValues(:,3) - eigenValues(:,1)) / eigenValues(:,3))
			steadyStateWave(:,1) = min(steadyStateWave(:,1), hR)
		end where
		
		! Limit the effect of the source term
		!   \cite[ch. 6.4.2]{george2006finite}
		steadyStateWave(:,2) = min(steadyStateWave(:,2), g* max(-hL * (bR - bL), -hR * (bR - bL) ) ) !TODO: use extra array for precomputing delta_B?
		steadyStateWave(:,2) = max(steadyStateWave(:,2), g* min(-hL * (bR - bL), -hR * (bR - bL) ) )
		
		rightHandSide(:,1) = rightHandSide(:,1) - steadyStateWave(:,1)
		!rightHandSide(:,2): no source term
		rightHandSide(:,3) = rightHandSide(:,3) - steadyStateWave(:,2)
		
		! everything is ready, solve the equations!
		!************************************************************************
		!* Solve linear equation begin
		!************************************************************************
		
		! direct solution of specific 3x3 system:
		! (       1           0       1           ) ( beta[0] ) = ( rightHandSide[0] )
		! ( eigenValues[0]    0  eigenValues[2]   ) ( beta[1] )   ( rightHandSide[1] )
		! ( eigenValues[0]^2  1  eigenValues[2]^2 ) ( beta[2] )   ( rightHandSide[2] )
		
		! step 1: solve the following 2x2 system (1st and 3rd column):
		! (       1              1           ) ( beta[0] ) = ( rightHandSide[0] )
		! ( eigenValues[0]  eigenValues[2]   ) ( beta[2] )   ( rightHandSide[1] )

		! compute the inverse of the wave speed difference:
		inverseDiff = 1.0_GRID_SR / (eigenValues(:,3) - eigenValues(:,1))
		! compute f-waves:
		beta(:,1) = (  eigenValues(:,3) * rightHandSide(:,1) - rightHandSide(:,2) ) * inverseDiff
		beta(:,3) = ( -eigenValues(:,1) * rightHandSide(:,1) + rightHandSide(:,2) ) * inverseDiff
		
		! step 2: solve 3rd row for beta[1]:
		beta(:,2) = rightHandSide(:,3) - eigenValues(:,1)*eigenValues(:,1) * beta(:,1) - eigenValues(:,3)*eigenValues(:,3) * beta(:,3)
		!************************************************************************
		!* Solve linear equation end
		!************************************************************************
		
		! initialize all fWaves and waveSpeeds to zero, so we don't need to set some of them to zero afterwards
		! (see commented code below)
		fWaves = 0.0_GRID_SR
		waveSpeeds = 0.0_GRID_SR
		
		! compute f-waves and wave-speeds
		where (wetDryState == WetDryWall)
			! zero ghost updates (wall boundary)
			! care about the left going wave (0) only
			fWaves(:,1,1) = beta(:,1) * eigenVectors(:,1,2)
			fWaves(:,2,1) = beta(:,1) * eigenVectors(:,1,3)
			
			! set the rest to zero
			!fWaves(:,1,2) = 0.0_GRID_SR   !is already zero
			!fWaves(:,2,2) = 0.0_GRID_SR   !is already zero
			!fwaves(:,1,3) = 0.0_GRID_SR   !is already zero
			!fwaves(:,2,3) = 0.0_GRID_SR   !is already zero
			
			waveSpeeds(:,1) = eigenValues(:,1)
			!waveSpeeds(:,2) = 0.0_GRID_SR   !is already zero
			!waveSpeeds(:,3) = 0.0_GRID_SR   !is already zero
		
		elsewhere (wetDryState == DryWetWall)
			! zero ghost updates (wall boundary)
			! care about the right going wave (2) only
			fWaves(:,1,3) = beta(:,3) * eigenVectors(:,3,2)
			fWaves(:,2,3) = beta(:,3) * eigenVectors(:,3,3)
			
			! set the rest to zero
			!fWaves(:,1,1) = 0.0_GRID_SR   !is already zero
			!fWaves(:,2,1) = 0.0_GRID_SR   !is already zero
			!fwaves(:,1,2) = 0.0_GRID_SR   !is already zero
			!fwaves(:,2,2) = 0.0_GRID_SR   !is already zero
			
			!waveSpeeds(:,1) = 0.0_GRID_SR   !is already zero
			!waveSpeeds(:,2) = 0.0_GRID_SR   !is already zero
			waveSpeeds(:,3) = eigenValues(:,3)
			
		elsewhere
			! compute f-waves (default)
			!do i=1,3
			!	fWaves(:,1,i) = beta(:,i) * eigenVectors(:,i,2)
			!	fWaves(:,2,i) = beta(:,i) * eigenVectors(:,i,3)
			!end do
			! loop unrolled:
			fWaves(:,1,1) = beta(:,1) * eigenVectors(:,1,2)
			fWaves(:,2,1) = beta(:,1) * eigenVectors(:,1,3)

			fWaves(:,1,2) = beta(:,2) * eigenVectors(:,2,2)
			fWaves(:,2,2) = beta(:,2) * eigenVectors(:,2,3)

			fWaves(:,1,3) = beta(:,3) * eigenVectors(:,3,2)
			fWaves(:,2,3) = beta(:,3) * eigenVectors(:,3,3)
			
			waveSpeeds(:,1) = eigenValues(:,1)
			waveSpeeds(:,2) = eigenValues(:,2)
			waveSpeeds(:,3) = eigenValues(:,3)
		end where
		
		!************************************************************************
		!* Compute Wave Decomposition End
		!************************************************************************
		! compute the updates from the three propagating waves
		! A^-\delta Q = \sum{s[i]<0} \beta[i] * r[i] = A^-\delta Q = \sum{s[i]<0} Z^i
		! A^+\delta Q = \sum{s[i]>0} \beta[i] * r[i] = A^-\delta Q = \sum{s[i]<0} Z^i
		do i=1,3
			where (waveSpeeds(:,i) < -r_eps)
				! left going
				upd_hL = upd_hL + fWaves(:,1,i)
				upd_huL = upd_huL + fWaves(:,2,i)
			elsewhere (waveSpeeds(:,i) > r_eps)
				! right going
				upd_hR = upd_hR + fWaves(:,1,i)
				upd_huR = upd_huR + fWaves(:,2,i)
			elsewhere
				! TODO: this case should not happen mathematically, but it does. Where is the bug? Machine accuracy only?
				! MB: >>> waveSpeeds set to 0 for cases WetDryWall and DryWetWall 
				upd_hL = upd_hL + 0.5_GRID_SR * fWaves(:,1,i)
				upd_huL = upd_huL + 0.5_GRID_SR * fWaves(:,2,i)
				
				upd_hR = upd_hR + 0.5_GRID_SR * fWaves(:,1,i)
				upd_huR = upd_huR + 0.5_GRID_SR * fWaves(:,2,i)
			end where
		end do
		
		! no net updates in DryDry case -> a non-vectorized solver would have stopped the computation
		! right after determining the Wet-Dry-State
		where (wetDryState == DryDry)
			waveSpeeds(:,1) = 0.0_GRID_SR
			waveSpeeds(:,2) = 0.0_GRID_SR
			waveSpeeds(:,3) = 0.0_GRID_SR
			upd_hL = 0.0_GRID_SR
			upd_huL = 0.0_GRID_SR
			upd_hR = 0.0_GRID_SR
			upd_huR = 0.0_GRID_SR
		end where
		
		! compute maximum wave speed (-> CFL-condition)
        maxWaveSpeed = max(maxWaveSpeed, maxVal(abs(waveSpeeds(:,:))))
		
		! inverse transformations
		call apply_transformations_after(transform_matrices, upd_huL, upd_hvL)
		call apply_transformations_after(transform_matrices, upd_huR, upd_hvR)
		
	end subroutine
	
    ! change base so hu/hv become ortogonal/perperdicular to edge
    subroutine apply_transformations_before(transform_matrices, hu, hv)
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,2,2), intent(in) :: transform_matrices
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE), intent(inout)      :: hu, hv           
        
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE)                 :: temp
        !DIR$ ASSUME_ALIGNED transform_matrices: 64
        !DIR$ ASSUME_ALIGNED hu: 64
        !DIR$ ASSUME_ALIGNED hv: 64
        !DIR$ ASSUME_ALIGNED temp: 64
        
        temp = hu
        hu = transform_matrices(:,1,1) * hu + transform_matrices(:,1,2) * hv
        hv = transform_matrices(:,2,1) * temp + transform_matrices(:,2,2) * hv
    end subroutine
    ! transform back to original base
    subroutine apply_transformations_after(transform_matrices, hu, hv)
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE,2,2), intent(in) :: transform_matrices
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE), intent(inout)      :: hu, hv           
        
        real(kind = GRID_SR), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE)                 :: temp
        !DIR$ ASSUME_ALIGNED transform_matrices: 64
        !DIR$ ASSUME_ALIGNED hu: 64
        !DIR$ ASSUME_ALIGNED hv: 64
        !DIR$ ASSUME_ALIGNED temp: 64
        
        temp = hu
        hu = transform_matrices(:,1,1) * hu + transform_matrices(:,2,1) * hv
        hv = transform_matrices(:,1,2) * temp + transform_matrices(:,2,2) * hv
    end subroutine  
	
END MODULE
#endif


#if defined(_SWE_HLLE) 

MODULE SWE_HLLE
    use Samoa_swe

    contains

    ! HLLE not vectorized - solves a single Riemann Problem
    subroutine compute_updates_hlle_single(hL, hR, huL, huR, hvL, hVR, bL, bR, net_updatesL, net_updatesR, maxWaveSpeed)
        real(kind = GRID_SR), intent(inout)  :: hL, hR, huL, huR, hvL, hvR, bL, bR
        real(kind = GRID_SR), intent(inout)  :: net_updatesL(3), net_updatesR(3)
        real(kind = GRID_SR), intent(inout)  :: maxWaveSpeed
        
        !local
        integer                :: i, j
        real(kind = GRID_SR)   :: uL, uR
        real(kind = GRID_SR)   :: sqrt_hL, sqrt_hR, sqrt_ghL, sqrt_ghR
        real(kind = GRID_SR)   :: half_g, sqrt_g
        integer(kind = BYTE)   :: wetDryState
        
        real(kind = GRID_SR), dimension(2)       :: characteristicSpeeds, roeSpeeds, extEinfeldtSpeeds, steadyStateWave
        real(kind = GRID_SR)                     :: hRoe, uRoe, sqrt_g_hRoe, hLLMiddleHeight, inverseDiff
        real(kind = GRID_SR), dimension(3)       :: eigenValues, rightHandSide, beta
        real(kind = GRID_SR), dimension(3,3)     :: eigenVectors
        real(kind = GRID_SR), parameter          :: r_eps = 1e-7
        
        real(kind = GRID_SR), dimension(2,3)     :: fWaves
        real(kind = GRID_SR), dimension(3)       :: waveSpeeds
        
        enum, bind(c) !constants to classify wet-dry-state of pairs of cells
            enumerator :: DryDry = 0
            enumerator :: WetWet = 1
            enumerator :: WetDryInundation = 2
            enumerator :: WetDryWall = 3
            enumerator :: WetDryWallInundation = 4
            enumerator :: DryWetInundation = 5
            enumerator :: DryWetWall = 6
            enumerator :: DryWetWallInundation = 7
        end enum
        
        enum, bind(c) ! constants to classify Riemann state of a pair of cells:
            enumerator :: DrySingleRarefaction = 0
            enumerator :: SingleRarefactionDry = 1
            enumerator :: ShockShock = 2
            enumerator :: ShockRarefaction = 3
            enumerator :: RarefactionShock = 4
            enumerator :: RarefactionRarefaction = 5
        end enum

     
        ! reset net updates
        net_updatesL = 0.0_GRID_SR
        net_updatesR = 0.0_GRID_SR
        maxWaveSpeed = 0.0_GRID_SR
        
        ! declare variables which are used over and over again
        half_g = 0.5_GRID_SR * g
        sqrt_g = sqrt(g)
        
        !************************************************************************
        !* Determine Wet Dry State Begin
        !* (determine the wet/dry state and compute local variables correspondly)
        !*************************************************************************
        
        ! compute speeds or set them to zero (dry cells)
        if (hL >= cfg%dry_tolerance) then
            uL = huL / hL
        else
            bL = bL + hL
            hL = 0.0_GRID_SR
            huL = 0.0_GRID_SR
            uL = 0.0_GRID_SR  !is already zero
        end if
        
        if (hR >= cfg%dry_tolerance) then
            uR = huR / hR
        else
            bR = bR + hR
            hR = 0.0_GRID_SR
            huR = 0.0_GRID_SR
            uR = 0_GRID_SR  !is already zero
        end if
        
        
        ! MB: determine wet/dry-state - try to start with most frequent case
        if (hL >= cfg%dry_tolerance) then
            if (hR >= cfg%dry_tolerance) then
                ! simple wet/wet case - expected as most frequently executed branch
                wetDryState = WetWet
            else ! hL >= cfg%dry_tolerance .and. hR < cfg%dry_tolerance
                ! we have a shoreline: left cell wet, right cell dry
                ! => check for simple inundation problems
                if (hL + bL > bR) then
                    ! => dry cell lies lower than wet cell
                    wetDryState = WetDryInundation
                else ! hL >= cfg%dry_tolerance .and. hR < cfg%dry_tolerance .and. hL + bL <= bR
                    ! => dry cell (right) lies higher than the wet cell (left)
                    ! Removed middle state computation -> assuming momentum is never high enough to overcome boundary on its own
                    hR = hL
                    uR = -uL
                    huR = -huL
                    bR = 0.0_GRID_SR
                    bL = 0.0_GRID_SR
                    wetDryState = WetDryWall
                end if
            end if
        else ! hL < cfg%dry_tolerance
            if (hR >= cfg%dry_tolerance) then
                ! we have a shoreline: left cell dry, right cell wet
                ! => check for simple inundation problems
                if (hR + bR > bL) then
                    ! => dry cell lies lower than wet cell
                    wetDryState = DryWetInundation
                else ! hL < cfg%dry_tolerance .and. hR >= cfg%dry_tolerance .and. hR + bR <= bL
                    ! => dry cell (left) lies higher than the wet cell (right)
                    ! Removed middle state computation -> assuming momentum is never high enough to overcome boundary on its own
                    hL = hR
                    uL = -uR
                    huL = -huR
                    bL = 0.0_GRID_SR
                    bR = 0.0_GRID_SR
                    wetDryState = DryWetWall
                end if
            else ! hL < cfg%dry_tolerance .and. hR < cfg%dry_tolerance
                wetDryState = DryDry
                ! nothing to do for dry/dry case, all netUpdates and maxWaveSpeed are 0
                return;
            end if
        end if

        !************************************************************************
        !* Determine Wet Dry State End
        !*************************************************************************
        

        ! precompute some terms which are fixed during
        ! the computation after some specific point
        sqrt_hL = sqrt(hL)
        sqrt_hR = sqrt(hR)
        sqrt_ghR = sqrt_g * sqrt_hR
        sqrt_ghL = sqrt_g * sqrt_hL
        
        !************************************************************************
        !* Compute Wave Decomposition Begin
        !************************************************************************
        
        ! compute eigenvalues of the jacobian matrices in states Q_{i-1} and Q_{i} (char. speeds)
        characteristicSpeeds(1) = uL - sqrt_ghL
        characteristicSpeeds(2) = uR + sqrt_ghR
        
        ! compute "Roe Speeds"
        hRoe = 0.5_GRID_SR * (hR + hL)
        uRoe = (uL * sqrt_hL + uR * sqrt_hR) / (sqrt_hL + sqrt_hR)
        
        ! optimization for dumb compilers
        sqrt_g_hRoe = sqrt_g * sqrt(hRoe)
        roeSpeeds(1) = uRoe - sqrt_g_hRoe
        roeSpeeds(2) = uRoe + sqrt_g_hRoe

        
        ! middle state computation removed, using basic Einfeldt speeds
        
        extEinfeldtSpeeds = 0_GRID_SR
        if (wetDryState == WetWet .or. wetDryState == WetDryWall .or. wetDryState == DryWetWall) then
            extEinfeldtSpeeds(1) = min(characteristicSpeeds(1), roeSpeeds(1))
            extEinfeldtSpeeds(2) = max(characteristicSpeeds(2), roeSpeeds(2))
        else if (hL < cfg%dry_tolerance) then ! MB: !!! cases DryWetInundation, DryWetWallInundation
            !ignore undefined speeds
            extEinfeldtSpeeds(1) = roeSpeeds(1)
            extEinfeldtSpeeds(2) = max(characteristicSpeeds(2), roeSpeeds(2))
        else if (hR < cfg%dry_tolerance) then ! MB: !!! cases WetDryInundation, WetDryWallInuncation
            ! ignore undefined speeds       
            extEinfeldtSpeeds(1) = min(characteristicSpeeds(1), roeSpeeds(1))
            extEinfeldtSpeeds(2) = roeSpeeds(2)
        end if

        
        ! HLL middle state
        !  \cite[theorem 3.1]{george2006finite}, \cite[ch. 4.1]{george2008augmented}
        hLLMiddleHeight = max(0.0_GRID_SR, (huL - huR + extEinfeldtSpeeds(2) * hR - extEinfeldtSpeeds(1) * hL) / (extEinfeldtSpeeds(2) - extEinfeldtSpeeds(1)) )
        
        ! define eigenvalues
        eigenvalues(1) = extEinfeldtSpeeds(1)
        eigenvalues(2) = 0.5_GRID_SR * (extEinfeldtSpeeds(1) + extEinfeldtSpeeds(2))
        eigenvalues(3) = extEinfeldtSpeeds(2)
        
        ! define eigenvectors
        ! MB: no longer used as system matrix
        !     -> but still used to compute f-waves
        eigenVectors(1,1) = 1.0_GRID_SR
        eigenVectors(2,1) = 0.0_GRID_SR
        eigenVectors(3,1) = 1.0_GRID_SR
        
        eigenVectors(1,2) = eigenValues(1)
        eigenVectors(2,2) = 0.0_GRID_SR
        eigenVectors(3,2) = eigenValues(3)
        
        eigenVectors(1,3) = eigenValues(1) * eigenValues(1)
        eigenVectors(2,3) = 1.0_GRID_SR
        eigenVectors(3,3) = eigenValues(3) * eigenValues(3)
        
        ! compute the jump in state
        rightHandSide(1) = hR - hL
        rightHandSide(2) = huR - huL
        rightHandSide(3) = (huR * uR + half_g * hR*hR) - (huL * uL + half_g * hL*hL)
        
        ! compute steady state wave
        steadyStateWave(1) = -(bR - bL)
        steadyStateWave(2) = -half_g * (hL + hR) * (bR - bL)

        
        ! preserve depth-positivity
        !  \cite[ch. 6.5.2]{george2006finite}, \cite[ch. 4.2.3]{george2008augmented}
        if (eigenvalues(1) < -r_eps .and. eigenvalues(3) > r_eps) then
            ! subsonic
            steadyStateWave(1) = max(steadyStateWave(1), hLLMiddleHeight * (eigenValues(3) - eigenValues(1)) / eigenValues(1))
            steadyStateWave(1) = min(steadyStateWave(1), hLLMiddleHeight * (eigenValues(3) - eigenValues(1)) / eigenValues(3))
        else if (eigenValues(1) > r_eps) then
            ! supersonic right TODO: motivation?
            steadyStateWave(1) = max(steadyStateWave(1), -hL)
            steadyStateWave(1) = min(steadyStateWave(1), hLLMiddleHeight * (eigenValues(3) - eigenValues(1)) / eigenValues(1))
        else if (eigenvalues(3) < -r_eps) then
            ! supersonic left TODO: motivation?
            steadyStateWave(1) = max(steadyStateWave(1), hLLMiddleHeight * (eigenValues(3) - eigenValues(1)) / eigenValues(3))
            steadyStateWave(1) = min(steadyStateWave(1), hR)
        end if
        
        ! Limit the effect of the source term
        !   \cite[ch. 6.4.2]{george2006finite}
        steadyStateWave(2) = min(steadyStateWave(2), g* max(-hL * (bR - bL), -hR * (bR - bL) ) ) !TODO: use extra array for precomputing delta_B?
        steadyStateWave(2) = max(steadyStateWave(2), g* min(-hL * (bR - bL), -hR * (bR - bL) ) )
        
        rightHandSide(1) = rightHandSide(1) - steadyStateWave(1)
        !rightHandSide(:,2): no source term
        rightHandSide(3) = rightHandSide(3) - steadyStateWave(2)

        
        ! everything is ready, solve the equations!
        !************************************************************************
        !* Solve linear equation begin
        !************************************************************************
        
        ! direct solution of specific 3x3 system:
        ! (       1           0       1           ) ( beta[0] ) = ( rightHandSide[0] )
        ! ( eigenValues[0]    0  eigenValues[2]   ) ( beta[1] )   ( rightHandSide[1] )
        ! ( eigenValues[0]^2  1  eigenValues[2]^2 ) ( beta[2] )   ( rightHandSide[2] )
        
        ! step 1: solve the following 2x2 system (1st and 3rd column):
        ! (       1              1           ) ( beta[0] ) = ( rightHandSide[0] )
        ! ( eigenValues[0]  eigenValues[2]   ) ( beta[2] )   ( rightHandSide[1] )

        ! compute the inverse of the wave speed difference:
        inverseDiff = 1.0_GRID_SR / (eigenValues(3) - eigenValues(1))
        ! compute f-waves:
        beta(1) = (  eigenValues(3) * rightHandSide(1) - rightHandSide(2) ) * inverseDiff
        beta(3) = ( -eigenValues(1) * rightHandSide(1) + rightHandSide(2) ) * inverseDiff
        
        ! step 2: solve 3rd row for beta[1]:
        beta(2) = rightHandSide(3) - eigenValues(1)*eigenValues(1) * beta(1) - eigenValues(3)*eigenValues(3) * beta(3)
        !************************************************************************
        !* Solve linear equation end
        !************************************************************************
        
        ! initialize all fWaves and waveSpeeds to zero, so we don't need to set some of them to zero afterwards
        ! (see commented code below)
        fWaves = 0.0_GRID_SR
        waveSpeeds = 0.0_GRID_SR
        
        ! compute f-waves and wave-speeds
        if (wetDryState == WetDryWall) then
            ! zero ghost updates (wall boundary)
            ! care about the left going wave (0) only
            fWaves(1,1) = beta(1) * eigenVectors(1,2)
            fWaves(2,1) = beta(1) * eigenVectors(1,3)
            
            ! set the rest to zero
            !fWaves(1,2) = 0.0_GRID_SR   !is already zero
            !fWaves(2,2) = 0.0_GRID_SR   !is already zero
            !fwaves(1,3) = 0.0_GRID_SR   !is already zero
            !fwaves(2,3) = 0.0_GRID_SR   !is already zero
            
            waveSpeeds(1) = eigenValues(1)
            !waveSpeeds(2) = 0.0_GRID_SR   !is already zero
            !waveSpeeds(3) = 0.0_GRID_SR   !is already zero
        
        else if (wetDryState == DryWetWall) then
            ! zero ghost updates (wall boundary)
            ! care about the right going wave (2) only
            fWaves(1,3) = beta(3) * eigenVectors(3,2)
            fWaves(2,3) = beta(3) * eigenVectors(3,3)
            
            ! set the rest to zero
            !fWaves(1,1) = 0.0_GRID_SR   !is already zero
            !fWaves(2,1) = 0.0_GRID_SR   !is already zero
            !fwaves(1,2) = 0.0_GRID_SR   !is already zero
            !fwaves(2,2) = 0.0_GRID_SR   !is already zero
            
            !waveSpeeds(1) = 0.0_GRID_SR   !is already zero
            !waveSpeeds(2) = 0.0_GRID_SR   !is already zero
            waveSpeeds(3) = eigenValues(3)
            
        else 
            ! compute f-waves (default)
            !do i=1,3
            !   fWaves(1,i) = beta(i) * eigenVectors(i,2)
            !   fWaves(2,i) = beta(i) * eigenVectors(i,3)
            !end do
            ! loop unrolled:
            fWaves(1,1) = beta(1) * eigenVectors(1,2)
            fWaves(2,1) = beta(1) * eigenVectors(1,3)

            fWaves(1,2) = beta(2) * eigenVectors(2,2)
            fWaves(2,2) = beta(2) * eigenVectors(2,3)

            fWaves(1,3) = beta(3) * eigenVectors(3,2)
            fWaves(2,3) = beta(3) * eigenVectors(3,3)
            
            waveSpeeds(1) = eigenValues(1)
            waveSpeeds(2) = eigenValues(2)
            waveSpeeds(3) = eigenValues(3)
        end if
        
        !************************************************************************
        !* Compute Wave Decomposition End
        !************************************************************************
        

        ! compute the updates from the three propagating waves
        ! A^-\delta Q = \sum{s[i]<0} \beta[i] * r[i] = A^-\delta Q = \sum{s[i]<0} Z^i
        ! A^+\delta Q = \sum{s[i]>0} \beta[i] * r[i] = A^-\delta Q = \sum{s[i]<0} Z^i
        do i=1,3
            if (waveSpeeds(i) < -r_eps) then
                ! left going
                net_updatesL(1) = net_updatesL(1) + fWaves(1,i)
                net_updatesL(2) = net_updatesL(2) + fWaves(2,i)
            else if (waveSpeeds(i) > r_eps) then
                ! right going
                net_updatesR(1) = net_updatesR(1) + fWaves(1,i)
                net_updatesR(2) = net_updatesR(2) + fWaves(2,i)
            else
                ! TODO: this case should not happen mathematically, but it does. Where is the bug? Machine accuracy only?
                ! MB: >>> waveSpeeds set to 0 for cases WetDryWall and DryWetWall 
                net_updatesL(1) = net_updatesL(1) + 0.5_GRID_SR * fWaves(1,i)
                net_updatesL(2) = net_updatesL(2) + 0.5_GRID_SR * fWaves(2,i)
                
                net_updatesR(1) = net_updatesR(1) + 0.5_GRID_SR * fWaves(1,i)
                net_updatesR(2) = net_updatesR(2) + 0.5_GRID_SR * fWaves(2,i)
            end if
        end do
        
        ! compute maximum wave speed (-> CFL-condition)
        maxWaveSpeed = maxVal(abs(waveSpeeds(:)))

      
    end subroutine

end module
    
#endif
