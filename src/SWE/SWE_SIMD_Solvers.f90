#if defined(_SWE_SIMD)

#include "Compilation_control.f90"

MODULE SWE_SIMD_Solvers
	use Samoa_swe
	use SWE_SIMD
	
	contains

	subroutine compute_updates_fwave_simd(normals, hL, huL, hvL, bL, hR, huR, hvR, bR, maxWaveSpeed)
		real(kind = GRID_SR), intent(in)    	:: normals(2,3)
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT), intent(inout)	:: hL, hR, huL, huR, hvL, hvR, bL, bR
		real(kind = GRID_SR), intent(inout)									:: maxWaveSpeed
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT)				:: uL, uR, vL, vR
		
		!local
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT,2,2)	:: transform_matrices
		integer								:: i, j
		real(kind = GRID_SR)										:: hstar, s1m, s2m
		logical														:: rare1, rare2
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT)		:: upd_hL, upd_hR, upd_huL, upd_huR, upd_hvL, upd_hvR
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT,3)		:: waveSpeeds
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT,3,3)	:: fwaves
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT,3)		:: wall
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT)		:: delphi
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT)		:: sL, sR, uhat, chat, sRoe1, sRoe2, sE1, sE2
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT)		:: delh, delhu, delb, deldelphi, delphidecomp, beta1, beta2
		!DIR$ ASSUME_ALIGNED transform_matrices: 64
		!DIR$ ASSUME_ALIGNED hL: 64
		!DIR$ ASSUME_ALIGNED hR: 64
		!DIR$ ASSUME_ALIGNED huL: 64
		!DIR$ ASSUME_ALIGNED huR: 64
		!DIR$ ASSUME_ALIGNED hvL: 64
		!DIR$ ASSUME_ALIGNED hvR: 64
		!DIR$ ASSUME_ALIGNED uL: 64
		!DIR$ ASSUME_ALIGNED uR: 64
		!DIR$ ASSUME_ALIGNED vL: 64
		!DIR$ ASSUME_ALIGNED vR: 64
		!DIR$ ASSUME_ALIGNED bL: 64
		!DIR$ ASSUME_ALIGNED bR: 64
		
		!DIR$ ASSUME_ALIGNED upd_hL: 64
		!DIR$ ASSUME_ALIGNED upd_hR: 64
		!DIR$ ASSUME_ALIGNED upd_huL: 64
		!DIR$ ASSUME_ALIGNED upd_huR: 64
		!DIR$ ASSUME_ALIGNED upd_hvL: 64
		!DIR$ ASSUME_ALIGNED upd_hvR: 64

		!DIR$ ASSUME_ALIGNED waveSpeeds: 64
		!DIR$ ASSUME_ALIGNED fwaves: 64
		!DIR$ ASSUME_ALIGNED wall: 64
		!DIR$ ASSUME_ALIGNED delphi: 64
		!DIR$ ASSUME_ALIGNED sL: 64
		!DIR$ ASSUME_ALIGNED sR: 64
		!DIR$ ASSUME_ALIGNED uhat: 64
		!DIR$ ASSUME_ALIGNED chat: 64
		!DIR$ ASSUME_ALIGNED sRoe1: 64
		!DIR$ ASSUME_ALIGNED sRoe2: 64
		!DIR$ ASSUME_ALIGNED sE1: 64
		!DIR$ ASSUME_ALIGNED sE2: 64
		!DIR$ ASSUME_ALIGNED delh: 64
		!DIR$ ASSUME_ALIGNED delhu: 64
		!DIR$ ASSUME_ALIGNED delb: 64
		!DIR$ ASSUME_ALIGNED deldelphi: 64
		!DIR$ ASSUME_ALIGNED delphidecomp: 64
		!DIR$ ASSUME_ALIGNED beta1: 64
		!DIR$ ASSUME_ALIGNED beta2: 64

		! compute transformations matrices
		associate(geom => SWE_SIMD_geometry)
			do i=1,_SWE_SIMD_NUM_EDGES 
				transform_matrices(i,1,:) = normals(:,geom%edges_orientation(i))
				transform_matrices(i,2,:) = [ - normals(2,geom%edges_orientation(i)), normals(1,geom%edges_orientation(i)) ]
			end do
		end associate


		! *** F-Wave solver *** (based on geoclaw implementation)
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
		do i=1,_SWE_SIMD_NUM_EDGES_ALIGNMENT
			if (hR(i) <= cfg%dry_tolerance) then
				call riemanntype(hL(i), hL(i), uL(i), -uL(i), hstar, s1m, s2m, rare1, rare2, 1, cfg%dry_tolerance, g)
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
				call riemanntype(hR(i), hR(i), -uR(i), uR(i), hstar, s1m, s2m, rare1, rare2, 1, cfg%dry_tolerance, g)
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
			
			!determine del vectors
			delh = hR - hL
			delhu = huR - huL
			delb = bR - bL
			
			deldelphi = -g * 0.5_GRID_SR * (hR + hL) * delb
			delphidecomp = delphi - deldelphi
			
			!flux decomposition
			beta1 = (sE2*delhu - delphidecomp) / (sE2 - sE1)
			beta2 = (delphidecomp - sE1*delhu) / (sE2 - sE1)
			
			waveSpeeds(:,1) = sE1
			waveSpeeds(:,2) = 0.5_GRID_SR * (sE1+sE2)
			waveSpeeds(:,3) = sE2
			! 1st nonlinear wave
			fwaves(:,1,1) = beta1
			fwaves(:,2,1) = beta1*sE1
			fwaves(:,3,1) = beta1*vL
			! 2nd nonlinear wave
			fwaves(:,1,3) = beta2
			fwaves(:,2,3) = beta2*sE2
			fwaves(:,3,3) = beta2*vR
			! advection of transverse wave
			fwaves(:,1,2) = 0.0_GRID_SR
			fwaves(:,2,2) = 0.0_GRID_SR
			fwaves(:,3,2) = huR*vR - huL*vL - fwaves(:,3,1)-fwaves(:,3,3)
			
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
		
		! compute maximum wave speed
		maxWaveSpeed = maxVal(abs(waveSpeeds(1:_SWE_SIMD_NUM_EDGES,:)))
		
		
		! inverse transformations
		call apply_transformations_after(transform_matrices, upd_huL, upd_hvL)
		call apply_transformations_after(transform_matrices, upd_huR, upd_hvR)
		
		! input variables are used as output for updates
		hL = upd_hL
		huL = upd_huL
		hvL = upd_hvL
		hR = upd_hR
		huR = upd_huR
		hvR = upd_hvR
				
	end subroutine
	
	subroutine compute_updates_hlle_simd(normals, hL, huL, hvL, bL, hR, huR, hvR, bR, maxWaveSpeed)
		real(kind = GRID_SR), intent(in)    	:: normals(2,3)
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT), intent(inout)	:: hL, hR, huL, huR, hvL, hvR, bL, bR
		real(kind = GRID_SR), intent(inout)									:: maxWaveSpeed
		
		!local
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT,2,2)	:: transform_matrices
		integer														:: i, j
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT)		:: upd_hL, upd_hR, upd_huL, upd_huR, upd_hvL, upd_hvR
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT)		:: uL, uR
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT)		:: sqrt_hL, sqrt_hR, sqrt_ghL, sqrt_ghR
		real(kind = GRID_SR)										:: half_g, sqrt_g
		integer(kind = BYTE), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT)		:: wetDryState
		
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT,2)		:: characteristicSpeeds, roeSpeeds, extEinfeldtSpeeds, steadyStateWave
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT)		:: hRoe, uRoe, sqrt_g_hRoe, hLLMiddleHeight, inverseDiff
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT,3)		:: eigenValues, rightHandSide, beta
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT,3,3)	:: eigenVectors
		real(kind = GRID_SR), parameter								:: r_eps = 1e-7
		
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT,2,3)	:: fWaves
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT,3)		:: waveSpeeds
		
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

		associate(geom => SWE_SIMD_geometry)
			! compute transformations matrices
			do i=1,_SWE_SIMD_NUM_EDGES
				transform_matrices(i,1,:) = normals(:,geom%edges_orientation(i))
				transform_matrices(i,2,:) = [ - normals(2,geom%edges_orientation(i)), normals(1,geom%edges_orientation(i)) ]
			end do
		end associate

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
		maxWaveSpeed = maxVal(abs(waveSpeeds(1:_SWE_SIMD_NUM_EDGES,:)))
		
		! inverse transformations
		call apply_transformations_after(transform_matrices, upd_huL, upd_hvL)
		call apply_transformations_after(transform_matrices, upd_huR, upd_hvR)
		
		! input variables are used as output for updates
		hL = upd_hL
		huL = upd_huL
		hvL = upd_hvL
		hR = upd_hR
		huR = upd_huR
		hvR = upd_hvR
		
			
	end subroutine
	
	! change base so hu/hv become ortogonal/perperdicular to edge
	subroutine apply_transformations_before(transform_matrices, hu, hv)
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT,2,2), intent(in)	:: transform_matrices
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT), intent(inout)		:: hu, hv			
		
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT)					:: temp
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
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT,2,2), intent(in)	:: transform_matrices
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT), intent(inout)		:: hu, hv			
		
		real(kind = GRID_SR), dimension(_SWE_SIMD_NUM_EDGES_ALIGNMENT)					:: temp
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
