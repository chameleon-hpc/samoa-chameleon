#include "Compilation_control.f90"

MODULE SWE_FV_solver
  use SFC_data_types
  use SWE_data_types
  use SWE_DG_limiter
  use SWE_Patch
  use SWE_Riemann_solver_wrapper
  use c_bind_riemannsolvers

  public fv_patch_solver
  
contains
subroutine fv_patch_solver(element, update1, update2, update3, dt, maxWaveSpeed)
            type(t_element_base) , intent(inout)	:: element
            type(num_cell_update), intent(inout)  :: update1, update2, update3
            real(kind = GRID_SR),intent(in)       :: dt
            real(kind = GRID_SR),intent(inout)    :: maxWaveSpeed
            
            integer                               :: i, j, ind, index_bnd
            type(num_cell_update)                 :: tmp !> ghost cells in correct order 
            real(kind = GRID_SR)                  :: volume, edge_lengths(3), dQ_max_norm, dt_div_volume, maxWaveSpeedLocal


            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_ORDER_SQUARE)                :: dQ_H, dQ_HU, dQ_HV !> deltaQ, used to compute cell updates
            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)           :: hL, huL, hvL, bL
            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)           :: hR, huR, hvR, bR
            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)           :: upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR

#if !defined (_SWE_USE_PATCH_SOLVER)
            type(t_state), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE)      :: edges_a, edges_b
            type(t_update)                                              :: update_a, update_b
            real(kind = GRID_SR), dimension(2,3)                        :: normals
#endif

            type(t_update),  DIMENSION((_SWE_DG_ORDER+1)**2)		:: flux1, flux2, flux3
            type(t_state), dimension((_SWE_DG_ORDER+1)**2)              :: edge_l, edge_m,edge_r
            real(kind= GRID_SR):: dx,delta(3),max_neighbour(3),min_neighbour(3),data_max_val(3),data_min_val(3)
            integer:: indx,indx_mir,k,edgesLeft
            real (kind = GRID_SR), dimension (_SWE_PATCH_ORDER_SQUARE):: H_old, HU_old, HV_old
            logical :: drying,troubled,coarsen,refine
            real (kind=GRID_SR) :: refinement_threshold = 0.50_GRID_SR
            real(kind=GRID_SR),Dimension(_SWE_PATCH_ORDER_SQUARE,3) :: bnd_flux_l,bnd_flux_m,bnd_flux_r
            real(kind=GRID_SR),Dimension(3) :: edge_sizes

            call compute_dg_boundary_fluxes(update1,update2,update3,&
                                           bnd_flux_l,bnd_flux_m,bnd_flux_r)
#if defined(_OPT_KERNELS)            
            !DIR$ ASSUME_ALIGNED hL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED hR: ALIGNMENT
            !DIR$ ASSUME_ALIGNED huL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED huR: ALIGNMENT
            !DIR$ ASSUME_ALIGNED hvL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED hvR: ALIGNMENT
            !DIR$ ASSUME_ALIGNED bL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED bR: ALIGNMENT
            !DIR$ ASSUME_ALIGNED upd_hL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED upd_hR: ALIGNMENT
            !DIR$ ASSUME_ALIGNED upd_huL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED upd_huR: ALIGNMENT
            !DIR$ ASSUME_ALIGNED upd_hvL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED upd_hvR: ALIGNMENT
#endif
            
            ! using patches, but applying geoclaw solvers on single edges            
            ! the normals are only needed in this case.            
            ! copy/compute normal vectors
            ! normal for type 2 edges is equal to the 2nd edge's normal
            normals(:,2) = element%edges(2)%ptr%transform_data%normal
            ! normal for types 1 and 3 depend on cell orientation.
            ! notice that normal for type 1 points inwards. That's why it is multiplied by -1.
            if (element%cell%geometry%i_plotter_type < 0) then ! orientation = backward
               normals(:,1) = - element%edges(1)%ptr%transform_data%normal
               normals(:,3) = element%edges(3)%ptr%transform_data%normal
            else ! orientation = forward, so reverse edge order
               normals(:,1) = - element%edges(3)%ptr%transform_data%normal
               normals(:,3) = element%edges(1)%ptr%transform_data%normal
            end if

            ! init some variables

            dQ_H = 0.0_GRID_SR
            dQ_HU = 0.0_GRID_SR
            dQ_HV = 0.0_GRID_SR
            
            maxWaveSpeed = 0.0_GRID_SR
            volume = cfg%scaling * cfg%scaling * element%cell%geometry%get_volume() / (_SWE_PATCH_ORDER_SQUARE)
            dt_div_volume = dt / volume
            edge_lengths = cfg%scaling * element%cell%geometry%get_edge_sizes() / _SWE_PATCH_ORDER
            edge_sizes=element%cell%geometry%get_edge_sizes()
            dx=edge_sizes(1)*cfg%scaling

            associate(data => element%cell%data_pers, geom => SWE_PATCH_geometry)
              !iterate over chunks
              do i=1, _SWE_PATCH_NUM_EDGES, _SWE_PATCH_SOLVER_CHUNK_SIZE ! i -> init of chunk

                 edgesLeft = _SWE_PATCH_NUM_EDGES - i + 1
                 
                 call init_chunk_data(i,edgesLeft,element,update1,update2,update3,&
                                      hL,huL,hvL,bL,hR,huR,hvR,bR)

                 ! compute net_updates -> solve Riemann problems within chunk
                 call compute_net_updates(i,edgesLeft,maxWaveSpeed,normals,&
                                       hL,huL,hvL,bL,hR,huR,hvR,bR,&
                                       upd_hL,upd_huL,upd_hvL,&
                                       upd_hR,upd_huR,upd_hvR)
                 ! compute dQ
                 call compute_dQ(i,edgesLeft,&
                      element,&
                      edge_lengths,&
                      update3%troubled.eq.DG,&
                      update2%troubled.eq.DG,&
                      update1%troubled.eq.DG,&
                      upd_hL,upd_huL,upd_hvL, upd_hR, upd_huR, upd_hvR,&
                      dQ_H,dQ_HU,dQ_HV)
              end do

               !set refinement condition -> Here I am using the same ones as in the original no-patches implementation, but considering only the max value.
               element%cell%geometry%refinement = 0
               dQ_max_norm = maxval(abs(dQ_H))

               dQ_H  = dQ_H * (-dt_div_volume)
               dQ_HU = dQ_HU * (-dt_div_volume)
               dQ_HV = dQ_HV * (-dt_div_volume)

               ! if land is flooded, init water height to dry tolerance and
               ! velocity to zero
               where (data%H < data%B + cfg%dry_tolerance .and. dQ_H > 0.0_GRID_SR)
                  !                  data%H  = data%B + cfg%dry_tolerance
                  data%HU = 0.0_GRID_SR
                  data%HV = 0.0_GRID_SR
               end where

                data%H  = data%H  + dQ_H
                data%HU = data%HU + dQ_HU
                data%HV = data%HV + dQ_HV

                call apply_dg_boundary_fluxes(element,&
                update1%troubled.eq.DG,&
                update2%troubled.eq.DG,&
                update3%troubled.eq.DG,&
                bnd_flux_l,bnd_flux_m,bnd_flux_r,dt,dx)
                
               ! if the water level falls below the dry tolerance, set water level to 0 and velocity to 0          
                where (data%H < data%B + cfg%dry_tolerance)
 !                  data%H = data%B 
                   data%HU = 0.0_GRID_SR
                   data%HV = 0.0_GRID_SR
                end where

               call update_dg_solution(element)
          end associate

        end subroutine fv_patch_solver


        subroutine compute_geoclaw_flux(normal, QL, QR, fluxL, fluxR,max_wave_speed)
          type(t_state), intent(in)           :: QL, QR
          type(t_update), intent(out)         :: fluxL, fluxR
          real(kind = GRID_SR), intent(in)    :: normal(2)
          real(kind = GRID_SR)                :: max_wave_speed
          
          real(kind = GRID_SR)				:: transform_matrix(2, 2)
          real(kind = GRID_SR)			    :: net_updatesL(3), net_updatesR(3)
          real(kind = GRID_SR)                :: pL(2), pR(2), hL, hR, bL, bR
          
          transform_matrix(1, :) = normal
          transform_matrix(2, :) = [-normal(2), normal(1)]
          
          pL = matmul(transform_matrix, QL%p)
          pR = matmul(transform_matrix, QR%p)
          hL = QL%h - QL%b
          hR = QR%h - QR%b
          bL = QL%b
          bR = QR%b
          
#           if defined(_FWAVE_FLUX)
          call c_bind_geoclaw_solver(GEOCLAW_FWAVE, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_AUG_RIEMANN_FLUX)
          call c_bind_geoclaw_solver(GEOCLAW_AUG_RIEMANN, 10, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_HLLE_FLUX)
          call compute_updates_hlle_single(hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, net_updatesL, net_updatesR, max_wave_speed)
#           endif
          
          fluxL%h = net_updatesL(1)
          fluxL%p = matmul(net_updatesL(2:3), transform_matrix)
          fluxR%h = net_updatesR(1)
          fluxR%p = matmul(net_updatesR(2:3), transform_matrix)
        end subroutine compute_geoclaw_flux

        subroutine update_dg_solution(element)
          type(t_element_base), intent(inout) :: element
          if(.not.isCoast(element%cell%data_pers%troubled)) then
             call apply_mue_h(element%cell%data_pers%h,&
                              element%cell%data_pers%b,&
                              element%cell%data_pers%Q%h,&
                              element%cell%data_pers%Q%b)
             call apply_mue(element%cell%data_pers%hu,element%cell%data_pers%Q%p(1))
             call apply_mue(element%cell%data_pers%hv,element%cell%data_pers%Q%p(2))
          endif
        end subroutine update_dg_solution


        subroutine init_chunk_data(i,edgesLeft,element,update1,update2,update3,hL,huL,hvL,bL,hR,huR,hvR,bR)
          integer, intent(in) :: i,edgesLeft
          type(t_element_base) , intent(in)	:: element
          type(num_cell_update), intent(in) :: update1, update2, update3
          real(kind=GRID_SR),DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE),intent(out)::hL,huL,hvL,bL
          real(kind=GRID_SR),DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE),intent(out)::hR,huR,hvR,bR
          !local variables
          integer :: j,ind

#if defined(_OPT_KERNELS)            
            !DIR$ ASSUME_ALIGNED hL:  ALIGNMENT
            !DIR$ ASSUME_ALIGNED hR:  ALIGNMENT
            !DIR$ ASSUME_ALIGNED huL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED huR: ALIGNMENT
            !DIR$ ASSUME_ALIGNED hvL: ALIGNMENT
            !DIR$ ASSUME_ALIGNED hvR: ALIGNMENT
            !DIR$ ASSUME_ALIGNED bL:  ALIGNMENT
            !DIR$ ASSUME_ALIGNED bR:  ALIGNMENT
#endif
          
          associate(data => element%cell%data_pers, geom => SWE_PATCH_geometry)
            ! if this is the last chunk and it is not a full chunk, it is necessary to set everything to 0 to avoid using last iteration values
            if (i + _SWE_PATCH_SOLVER_CHUNK_SIZE -1 > _SWE_PATCH_NUM_EDGES) then
               hL = 0.0_GRID_SR
               huL = 0.0_GRID_SR
               hvL = 0.0_GRID_SR
               bL = 0.0_GRID_SR
               hR = 0.0_GRID_SR
               huR = 0.0_GRID_SR
               hvR = 0.0_GRID_SR
               bR = 0.0_GRID_SR
            end if
            
            do j=1, min(edgesLeft,_SWE_PATCH_SOLVER_CHUNK_SIZE) !j -> position inside chunk
               ind = i + j - 1 ! actual index
               
               ! left
               if (geom%edges_a(ind) <= _SWE_PATCH_ORDER_SQUARE) then
                  hL(j) = data%H(geom%edges_a(ind))
                  huL(j) = data%HU(geom%edges_a(ind))
                  hvL(j) = data%HV(geom%edges_a(ind))
                  bL(j) = data%B(geom%edges_a(ind))
               else if (geom%edges_a(ind) <= _SWE_PATCH_ORDER_SQUARE + _SWE_PATCH_ORDER) then
                  hL(j) = update1%H(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE)
                  huL(j) = update1%HU(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE)
                  hvL(j) = update1%HV(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE)
                  bL(j) = update1%B(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE)
               else if (geom%edges_a(ind) <= _SWE_PATCH_ORDER_SQUARE + 2*_SWE_PATCH_ORDER)then
                  hL(j) =update2%H(geom%edges_a(ind)-_SWE_PATCH_ORDER_SQUARE-_SWE_PATCH_ORDER)
                  huL(j)=update2%HU(geom%edges_a(ind)-_SWE_PATCH_ORDER_SQUARE-_SWE_PATCH_ORDER)
                  hvL(j)=update2%HV(geom%edges_a(ind)-_SWE_PATCH_ORDER_SQUARE-_SWE_PATCH_ORDER)
                  bL(j) =update2%B(geom%edges_a(ind)-_SWE_PATCH_ORDER_SQUARE-_SWE_PATCH_ORDER)
               else
                  hL(j) = update3%H(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                  huL(j) = update3%HU(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                  hvL(j) = update3%HV(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                  bL(j) = update3%B(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
               end if
               
               ! right
               if (geom%edges_b(ind) <= _SWE_PATCH_ORDER_SQUARE) then
                  hR(j) = data%H(geom%edges_b(ind))
                  huR(j) = data%HU(geom%edges_b(ind))
                  hvR(j) = data%HV(geom%edges_b(ind))
                  bR(j) = data%B(geom%edges_b(ind))
               else if (geom%edges_b(ind) <= _SWE_PATCH_ORDER_SQUARE + _SWE_PATCH_ORDER) then
                  hR(j) = update1%H(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE)
                  huR(j) = update1%HU(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE)
                  hvR(j) = update1%HV(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE)
                  bR(j) = update1%B(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE)
               else if (geom%edges_b(ind) <= _SWE_PATCH_ORDER_SQUARE + 2*_SWE_PATCH_ORDER) then
                  hR(j) = update2%H(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                  huR(j) = update2%HU(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                  hvR(j) = update2%HV(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                  bR(j) = update2%B(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
               else
                  hR(j) = update3%H(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                  huR(j) = update3%HU(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                  hvR(j) = update3%HV(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                  bR(j) = update3%B(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
               end if
            end do
          end associate
        end subroutine init_chunk_data

        subroutine compute_net_updates(i,edgesLeft,maxWaveSpeed,&
                                       normals,&
                                       hL,huL,hvL,bL,hR,huR,hvR,bR,&
                                       upd_hL,upd_huL,upd_hvL,&
                                       upd_hR,upd_huR,upd_hvR)
          integer, intent(in) :: i,edgesLeft
          real(kind = GRID_SR),intent(inout) :: maxWaveSpeed
          real(kind = GRID_SR), dimension(2,3),intent(in) :: normals
          real(kind=GRID_SR),DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE),intent(inout)::hL,huL,hvL,bL
          real(kind=GRID_SR),DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE),intent(inout)::hR,huR,hvR,bR
          real(kind=GRID_SR),DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE),intent(out)  ::upd_hL,upd_huL,upd_hvL, upd_hR, upd_huR, upd_hvR
          !local variables
          integer :: j,ind
          real(kind=GRID_SR) :: maxWaveSpeedLocal
          real(kind=GRID_SR),DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE) :: normals_x, normals_y
#if !defined(_OPT_KERNELS)          
          type(t_state), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE)      :: edges_a, edges_b
          type(t_update)                                              :: update_a, update_b
#endif          

          upd_hL = 0.0_GRID_SR
          upd_hR = 0.0_GRID_SR
          upd_huL= 0.0_GRID_SR
          upd_huR= 0.0_GRID_SR
          upd_hvL= 0.0_GRID_SR
          upd_hvR= 0.0_GRID_SR
          associate(geom => SWE_PATCH_geometry)
#if defined(_OPT_KERNELS)

          do j=1, min(edgesLeft,_SWE_PATCH_SOLVER_CHUNK_SIZE) !j -> position inside chunk
             ind = i + j - 1 ! actual index
             normals_x(j) = normals(1,geom%edges_orientation(ind))
             normals_y(j) = normals(2,geom%edges_orientation(ind))
          end do
          
          !$OMP SIMD PRIVATE(maxWaveSpeedLocal) REDUCTION(max: maxWaveSpeed)
          do j=1, min(edgesLeft, _SWE_PATCH_SOLVER_CHUNK_SIZE)
             call SWE_compute_fluxes(normals_x(j), normals_y(j),&
             hL(j), hR(j),&
             huL(j), huR(j),&
             hvL(j), hvR(j),&
             bL(j), bR(j), upd_hL(j),&
             upd_hR(j), upd_huL(j),&
             upd_huR(j), upd_hvL(j),&
             upd_hvR(j), maxWaveSpeedLocal)
             maxWaveSpeed = max(maxWaveSpeed, maxWaveSpeedLocal)
          end do
#else 
          do j=1, min(edgesLeft,_SWE_PATCH_SOLVER_CHUNK_SIZE) !j -> position inside chunk
             ind = i + j - 1 ! actual index
             edges_a(j)%h    = hL(j)
             edges_a(j)%p(1) = huL(j)
             edges_a(j)%p(2) = hvL(j)
             edges_a(j)%b    = bL(j)
             edges_b(j)%h    = hR(j)
             edges_b(j)%p(1) = huR(j)
             edges_b(j)%p(2) = hvR(j)
             edges_b(j)%b = bR(j)
             maxWaveSpeedLocal = 0.0_GRID_SR
             call compute_geoclaw_flux(normals(:,geom%edges_orientation(ind)), edges_a(j), edges_b(j), update_a, update_b, maxWaveSpeedLocal)
             maxWaveSpeed=max(maxWaveSpeed,maxWaveSpeedLocal)
             
             upd_hL(j)  = update_a%h
             upd_huL(j) = update_a%p(1)
             upd_hvL(j) = update_a%p(2)
             upd_hR(j)  = update_b%h
             upd_huR(j) = update_b%p(1)
             upd_hvR(j) = update_b%p(2)
          end do
#                       endif

        end associate
      end subroutine compute_net_updates


      subroutine compute_dg_boundary_fluxes(update1,update2,update3,&
           bnd_flux_l,bnd_flux_m,bnd_flux_r)

        type(num_cell_update), intent(in)  :: update1, update2, update3
        real(kind=GRID_SR),Dimension(_SWE_PATCH_ORDER_SQUARE,3),intent(out) :: bnd_flux_l,bnd_flux_m,bnd_flux_r
        
        bnd_flux_l=0.0_GRID_SR
        bnd_flux_m=0.0_GRID_SR
        bnd_flux_r=0.0_GRID_SR
        !---------compute boundary integrals------------------!
        if(update1%troubled.eq.DG) then
           call apply_phi(matmul(s_m_inv,matmul(s_b_1_l,update1%flux(:)%h   )),bnd_flux_l(:,1))
           call apply_phi(matmul(s_m_inv,matmul(s_b_1_l,update1%flux(:)%p(1))),bnd_flux_l(:,2))
           call apply_phi(matmul(s_m_inv,matmul(s_b_1_l,update1%flux(:)%p(2))),bnd_flux_l(:,3))
        end if
        if(update2%troubled.eq.DG) then                                              
           call apply_phi(matmul(s_m_inv,matmul(s_b_2_m,update2%flux(:)%h   )),bnd_flux_m(:,1))
           call apply_phi(matmul(s_m_inv,matmul(s_b_2_m,update2%flux(:)%p(1))),bnd_flux_m(:,2))
           call apply_phi(matmul(s_m_inv,matmul(s_b_2_m,update2%flux(:)%p(2))),bnd_flux_m(:,3))
        end if
        if(update3%troubled.eq.DG) then
           call apply_phi(matmul(s_m_inv,matmul(s_b_3_r,update3%flux(:)%h   )),bnd_flux_r(:,1))
           call apply_phi(matmul(s_m_inv,matmul(s_b_3_r,update3%flux(:)%p(1))),bnd_flux_r(:,2))
           call apply_phi(matmul(s_m_inv,matmul(s_b_3_r,update3%flux(:)%p(2))),bnd_flux_r(:,3))
        end if
        !-----------------------------------------------------!
      end subroutine compute_dg_boundary_fluxes


      subroutine apply_dg_boundary_fluxes(element,&
                                          apply_l,apply_m,apply_r,&
                                          bnd_flux_l,bnd_flux_m,bnd_flux_r,&
                                          dt,dx)
        type(t_element_base),intent(inout)	:: element
        logical,intent(in) :: apply_l,apply_m,apply_r
        real(kind=GRID_SR),Dimension(_SWE_PATCH_ORDER_SQUARE,3),intent(in) :: bnd_flux_l,bnd_flux_m,bnd_flux_r
        real(kind = GRID_SR),intent(in) :: dt, dx
        
        associate(data => element%cell%data_pers, geom => SWE_PATCH_geometry)
          if(data%troubled.eq.NEIGHBOUR_TROUBLED) then
             if(apply_l) then
                data%H =data%H +bnd_flux_l(:,1)*(-dt/dx)
                data%HU=data%HU+bnd_flux_l(:,2)*(-dt/dx)
                data%HV=data%HV+bnd_flux_l(:,3)*(-dt/dx)
             end if
             if(apply_m) then
                data%H =data%H +bnd_flux_m(:,1)*(-dt/dx) 
                data%HU=data%HU+bnd_flux_m(:,2)*(-dt/dx)
                data%HV=data%HV+bnd_flux_m(:,3)*(-dt/dx)
             end if
             if(apply_r) then
                data%H =data%H +bnd_flux_r(:,1)*(-dt/dx) 
                data%HU=data%HU+bnd_flux_r(:,2)*(-dt/dx)
                data%HV=data%HV+bnd_flux_r(:,3)*(-dt/dx)
             end if
          end if
        end associate
      end subroutine apply_dg_boundary_fluxes
      
      subroutine compute_dQ(i,edgesLeft,&
           element,&
           edge_lengths,&
           left_dg,mid_dg,right_dg,&
           upd_hL,upd_huL,upd_hvL, upd_hR, upd_huR, upd_hvR,&
           dQ_H,dQ_HU,dQ_HV)
        integer :: i, edgesLeft
        type(t_element_base) :: element
        real(kind = GRID_SR),intent(in) :: edge_lengths(3)
        real(kind=GRID_SR),DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE),intent(inout) :: upd_hL,upd_huL,upd_hvL, upd_hR, upd_huR, upd_hvR
        logical,intent(in) :: left_dg,mid_dg,right_dg
        real(kind=GRID_SR),DIMENSION(_SWE_PATCH_ORDER_SQUARE),intent(out) :: dQ_H, dQ_HU, dQ_HV
        
        !local variables
        integer :: j,ind
        
        associate(data => element%cell%data_pers, geom => SWE_PATCH_geometry)
          
        do j=1, min(edgesLeft,_SWE_PATCH_SOLVER_CHUNK_SIZE)
           ind = i + j - 1 ! actual index
           if(data%troubled.eq.NEIGHBOUR_TROUBLED) then
              if(left_dg) then
                 if(any(ind == geom%right_bnd))then
                    upd_hL(j) = 0.0
                    upd_hR(j) = 0.0
                    upd_huL(j) = 0.0
                    upd_huR(j) = 0.0
                    upd_hvL(j) = 0.0
                    upd_hvR(j) = 0.0
                 end if
              end if
              if(mid_dg) then
                 if(any(ind == geom%mid_bnd))then
                    upd_hL(j) = 0.0
                    upd_hR(j) = 0.0
                    upd_huL(j) = 0.0
                    upd_huR(j) = 0.0
                    upd_hvL(j) = 0.0
                    upd_hvR(j) = 0.0
                 end if
              end if
              if(right_dg) then
                 if(any(ind == geom%left_bnd))then
                    upd_hL(j) = 0.0
                    upd_hR(j) = 0.0
                    upd_huL(j) = 0.0
                    upd_huR(j) = 0.0
                    upd_hvL(j) = 0.0
                    upd_hvR(j) = 0.0
                 end if
              end if
           end if
           
           if (geom%edges_a(ind) <= _SWE_PATCH_ORDER_SQUARE) then !ignore ghost cells
              dQ_H(geom%edges_a(ind)) = dQ_H(geom%edges_a(ind)) + upd_hL(j) * edge_lengths(geom%edges_orientation(ind))
              dQ_HU(geom%edges_a(ind)) = dQ_HU(geom%edges_a(ind)) + upd_huL(j) * edge_lengths(geom%edges_orientation(ind))
              dQ_HV(geom%edges_a(ind)) = dQ_HV(geom%edges_a(ind)) + upd_hvL(j) * edge_lengths(geom%edges_orientation(ind))
           end if
           if (geom%edges_b(ind) <= _SWE_PATCH_ORDER_SQUARE) then
              dQ_H(geom%edges_b(ind)) = dQ_H(geom%edges_b(ind)) + upd_hR(j) * edge_lengths(geom%edges_orientation(ind))
              dQ_HU(geom%edges_b(ind)) = dQ_HU(geom%edges_b(ind)) + upd_huR(j) * edge_lengths(geom%edges_orientation(ind))
              dQ_HV(geom%edges_b(ind)) = dQ_HV(geom%edges_b(ind)) + upd_hvR(j) * edge_lengths(geom%edges_orientation(ind))
           end if
        end do
      end associate
      end subroutine compute_dQ
        
      END MODULE SWE_FV_solver
