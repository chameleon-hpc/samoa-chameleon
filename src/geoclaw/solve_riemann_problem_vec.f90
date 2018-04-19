subroutine solve_riemann_problem_vec(hL, hR, huL, huR, hvL, hvR, bL, bR, pL, pR, drytol, g, sw, fw)      
    implicit none

    !input
    double precision hL, hR, huL, huR, hvL, hvR, bL, bR, pL, pR, g, drytol

    !output
    double precision sw(3), fw(3,3)

    !local only
    integer m,i,mw,maxiter
    double precision wall(3)
    double precision uR,uL,vR,vL,phiR,phiL,delphi
    double precision sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
    double precision s1m,s2m
    double precision hstar,hstartest,hstarHLL,sLtest,sRtest
    double precision tw,dxdc      
    double precision sqrt_ghL, sqrt_ghR
    logical rare1,rare2

    ! For completely dry states, do not skip problem (hinders
    ! vectorization), but rather solve artificial 0-valued problem.
    if (hL < drytol .and. hR < drytol) then
        hL = 0
        hR = 0
    endif

    !check for wet/dry boundary
    if (hR.gt.drytol) then
        uR=huR/hR
        vR=hvR/hR
        !phiR = 0.5d0*g*hR**2 + huR**2/hR
    else
        hR = 0.d0
        huR = 0.d0
        hvR = 0.d0
        uR = 0.d0
        vR = 0.d0
        !phiR = 0.d0
    endif

    if (hL.gt.drytol) then
        uL=huL/hL
        vL=hvL/hL
        !phiL = 0.5d0*g*hL**2 + huL**2/hL
    else
        hL=0.d0
        huL=0.d0
        hvL=0.d0
        uL=0.d0
        vL=0.d0
        !phiL = 0.d0
    endif

    wall(1) = 1.d0
    wall(2) = 1.d0
    wall(3) = 1.d0
    
    if (hR.le.drytol) then
        !DIR$ FORCEINLINE
        call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m, &
                                rare1,rare2,1,drytol,g)
        hstartest=max(hL,hstar)
        if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
            !bR=hstartest+bL
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
    else if (hL.le.drytol) then ! right surface is lower than left topo
        !DIR$ FORCEINLINE
        call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m, &
                                rare1,rare2,1,drytol,g)
        hstartest=max(hR,hstar)
        if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
            !bL=hstartest+bR
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

    ! pre-compute square roots
    sqrt_ghL = sqrt(g*hL)
    sqrt_ghR = sqrt(g*hR)

    !BUGFIX:
    !Problem: loss of significance may occur in phiR-phiL, causing divergence of the steady state.
    !Action:  Compute delphi=phiR-phiL explicitly. delphi is arithmetically equivalent to phiR-phiL, but with far smaller numerical loss.
    delphi = (huR - huL)*(uL + uR) - uL*uR*(hR - hL) + (0.5d0*g*(bR + hR - bL - hL)*(hR + hL)) - 0.5d0*g*(hR + hL)*(bR - bL)


    !determine wave speeds
    sL=uL-sqrt_ghL ! 1 wave speed of left state
    sR=uR+sqrt_ghR ! 2 wave speed of right state

    uhat=(sqrt_ghL*uL + sqrt_ghR*uR)/(sqrt_ghR+sqrt_ghL) ! Roe average
    chat=sqrt(g*0.5d0*(hR+hL)) ! Roe average
    sRoe1=uhat-chat ! Roe wave speed 1 wave
    sRoe2=uhat+chat ! Roe wave speed 2 wave

    sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
    sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave

    !--------------------end initializing...finally----------
    !solve Riemann problem.

    maxiter = 1

    !-DIR$ FORCEINLINE
    call riemann_aug_JCP(1,3,3, hL, hR, huL, huR, hvL, hvR, bL, bR, uL, uR, vL, vR, delphi, sE1, sE2, dryTol, g, sw, fw)
!      print *, "sw", sw
!      print *, "fw", fw
!      print *, "h_end", hL, hR
!      print *, "======"
!     call riemann_aug_JCP(maxiter,3,3,hL,hR,huL, &
!         huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,pL,pR,sE1,sE2, &
!         drytol,g,rho,sw1,sw2,sw3,fw11,fw12,fw13,fw21,fw22,fw23,fw31,fw32,fw33)

    ! For completely dry states, waves should be always zero!
    ! (these problems could have been skipped, but to allow 
    ! vectorization they weren't)
    if (hL < drytol .and. hR < drytol) then
        sw(:) = 0.0d0
        fw(:,:) = 0.0d0   
    endif


    !eliminate ghost fluxes for wall
    do mw=1,3
        sw(mw)=sw(mw)*wall(mw)

        fw(1,mw)=fw(1,mw)*wall(mw) 
        fw(2,mw)=fw(2,mw)*wall(mw)
        fw(3,mw)=fw(3,mw)*wall(mw)
    enddo

end subroutine

