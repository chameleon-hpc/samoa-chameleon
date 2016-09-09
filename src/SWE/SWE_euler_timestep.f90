! Sa(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_Euler_Timestep
		use SFC_edge_traversal

		use Samoa_swe
		use c_bind_riemannsolvers

#		if defined(_SWE_PATCH)
			use SWE_PATCH
			use SWE_PATCH_Solvers
#               if defined(_SWE_DG)
                        use SWE_DG_Matrices
#		endif
#		endif
#               if defined(_HLLE_FLUX)
                        use SWE_HLLE
#               endif

                        implicit none

        type num_traversal_data
            integer (kind = GRID_DI)			:: i_refinements_issued
        end type

        interface skeleton_op
            module procedure skeleton_array_op
            module procedure skeleton_scalar_op
        end interface

        interface bnd_skeleton_op
            module procedure bnd_skeleton_array_op
            module procedure bnd_skeleton_scalar_op
        end interface

		PUBLIC cell_to_edge_op
#if defined(_SWE_DG)
                PUBLIC dg_predictor
#endif
		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

#		define _GT_NAME							t_swe_euler_timestep_traversal


#		define _GT_EDGES

#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op

#		define _GT_SKELETON_OP					skeleton_op
#		define _GT_BND_SKELETON_OP				bnd_skeleton_op
#		define _GT_CELL_UPDATE_OP				cell_update_op
#		define _GT_CELL_LAST_TOUCH_OP			cell_last_touch_op
#		define _GT_NODE_WRITE_OP			    node_write_op
#		define _GT_EDGE_WRITE_OP			    edge_write_op

#		define _GT_NODE_MPI_TYPE

#		include "SFC_generic_traversal_ringbuffer.f90"

        subroutine create_node_mpi_type(mpi_node_type)
            integer, intent(out)            :: mpi_node_type

#           if defined(_MPI)
                type(t_node_data)                       :: node
                integer                                 :: blocklengths(2), types(2), disps(2), type_size, i_error
                integer (kind = MPI_ADDRESS_KIND)       :: lb, ub

                blocklengths(1) = 1
                blocklengths(2) = 1

                disps(1) = 0
                disps(2) = sizeof(node)

                types(1) = MPI_LB
                types(2) = MPI_UB

                call MPI_Type_struct(2, blocklengths, disps, types, mpi_node_type, i_error); assert_eq(i_error, 0)
                call MPI_Type_commit(mpi_node_type, i_error); assert_eq(i_error, 0)

                call MPI_Type_size(mpi_node_type, type_size, i_error); assert_eq(i_error, 0)
                call MPI_Type_get_extent(mpi_node_type, lb, ub, i_error); assert_eq(i_error, 0)

                assert_eq(0, lb)
                assert_eq(0, type_size)
                assert_eq(sizeof(node), ub)
#           endif
        end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		subroutine pre_traversal_grid_op(traversal, grid)
			type(t_swe_euler_timestep_traversal), intent(inout)		:: traversal
   type(t_grid), intent(inout)							    :: grid
            !! !!!!!!!print*,grid%r_dt                        
			if (cfg%r_max_time > 0.0_SR) then
                grid%r_dt = min(cfg%r_max_time, grid%r_dt)
            end if

			if (cfg%r_output_time_step > 0.0_SR) then
                grid%r_dt = min(cfg%r_output_time_step, grid%r_dt)
            end if

            !! !!!!!!!print*,grid%r_dt                        


#           if defined(_ASAGI)
                !if we are in the earthquake phase, limit the simulation time step by the earthquake time step
                if (grid%r_time < cfg%t_max_eq) then
                    grid%r_dt = min(grid%r_dt, cfg%dt_eq)
                end if
#           endif

#if defined(_SWE_DG)
 				grid%r_dt = grid%r_dt / (_SWE_DG_ORDER*4.0_GRID_SR +2.0_GRID_SR)
#else

#			if defined(_SWE_PATCH)
				grid%r_dt = grid%r_dt / _SWE_PATCH_ORDER
#			endif

#endif
                !!!!!!print*,grid%r_dt                        
			call scatter(grid%r_dt, grid%sections%elements_alloc%r_dt)
   
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_swe_euler_timestep_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

			grid%r_time = grid%r_time + grid%r_dt

            call reduce(traversal%i_refinements_issued, traversal%sections%i_refinements_issued, MPI_SUM, .true.)
            call reduce(grid%r_dt_new, grid%sections%elements_alloc%r_dt_new, MPI_MIN, .true.)

            grid%r_dt_new = cfg%courant_number * grid%r_dt_new

            if (rank_MPI == 0) then
                if (cfg%courant_number > grid%r_dt_new / grid%r_dt) then
                    _log_write(1, '("WARNING! Time step size was too big. dt (old): ", ES10.3, ", dt (CFL): ", ES10.3, ", maximum courant number: ", F0.3)') grid%r_dt, grid%r_dt_new / cfg%courant_number, grid%r_dt_new / grid%r_dt
                end if
            end if

            grid%r_dt = grid%r_dt_new
            !!! !!!!!!!print*,grid%r_dt
			call scatter(grid%r_time, grid%sections%elements_alloc%r_time)
		end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_swe_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section

			!this variable will be incremented for each cell with a refinement request
			traversal%i_refinements_issued = 0_GRID_DI
			section%r_dt_new = huge(1.0_SR)
		end subroutine

		function cell_to_edge_op(element, edge) result(rep)
			type(t_element_base), intent(in)						:: element
			type(t_edge_data), intent(in)						    :: edge
			type(num_cell_rep)										:: rep
			type(t_state), dimension(_SWE_CELL_SIZE)				:: Q
			integer(kind = GRID_SI)									:: i, j, i_edge
			real(kind = GRID_SR), dimension(2, _SWE_EDGE_SIZE)		:: dof_pos
			real(kind = GRID_SR), dimension(2, 3), parameter		:: edge_offsets = reshape([0.0, 0.0, 0.0, 1.0, 1.0, 0.0], [2, 3])
			real(kind = GRID_SR), dimension(2, 3), parameter		:: edge_vectors = reshape([0.0, 1.0, 1.0, -1.0, -1.0, 0.0], [2, 3])
#           if defined(_SWE_PATCH)
                integer                                             :: edge_type !1=left, 2=hypotenuse, 3=right
#if defined (_SWE_DG)                
                integer:: indx,indx_mir,k
#endif
#           else
                call gv_Q%read(element, Q)
                _log_write(6, '(3X, A)') "swe cell to edge op:"
                _log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q in: ", Q
                _log_write(6, '(4X, A, F0.3, 1X, F0.3)') "normal in : ", edge%transform_data%normal
#           endif


#           if (_SWE_CELL_SIZE > 1)
                i_edge = edge%transform_data%index
                _log_write(6, '(4X, A, I0)') "edge ", i_edge

				forall (i = 1 : _SWE_EDGE_SIZE)
					dof_pos(:, i) = edge_offsets(:, i_edge) + t_basis_flux_get_dof_coords(i) * edge_vectors(:, i_edge)
				end forall

				call lfs_flux%transform(edge%transform_data, dof_pos(1, :))
				call lfs_flux%transform(edge%transform_data, dof_pos(2, :))

				forall (i = 1 : _SWE_EDGE_SIZE)
					rep%Q(i)%h = t_basis_Q_eval(dof_pos(:, i), Q%h)
					rep%Q(i)%p(1) = t_basis_Q_eval(dof_pos(:, i), Q%p(1))
					rep%Q(i)%p(2) = t_basis_Q_eval(dof_pos(:, i), Q%p(2))
				end forall
#           elif defined(_SWE_PATCH)

                !find out which edge it is comparing its normal with cell normals
                ! obs.: negative plotter_types describe triangles in desired order: left, hypotenuse, right
                associate(cell_edge => ref_plotter_data(- abs(element%cell%geometry%i_plotter_type))%edges, normal => edge%transform_data%normal)
                    do i=1,3
                        if ((normal(1) == cell_edge(i)%normal(1) .and. normal(2) == cell_edge(i)%normal(2))   &
                                .or. (normal(1) == -cell_edge(i)%normal(1) .and. normal(2) == -cell_edge(i)%normal(2))) then
                            edge_type = i
                        end if
                    end do
                end associate
                
                associate(H => element%cell%data_pers%H, HU => element%cell%data_pers%HU, HV => element%cell%data_pers%HV, B => element%cell%data_pers%B)
                    ! copy boundary values to respective edges
                    ! left leg cells go to edge 1
                    ! hypotenuse cells go to edge 2
                    ! right leg cells go to edge 3
                    select case (edge_type)
                    case (1) !cells with id i*i+1 (left leg)
                        do i=0, _SWE_PATCH_ORDER - 1
                            rep%H(i+1) = H(i*i + 1)
                            rep%HU(i+1) = HU(i*i + 1)
                            rep%HV(i+1) = HV(i*i + 1)
                            rep%B(i+1) = B(i*i + 1)
                        end do
                    case (2) ! hypotenuse
                        do i=1, _SWE_PATCH_ORDER
                            rep%H(i) = H((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
                            rep%HU(i) = HU((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
                            rep%HV(i) = HV((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
                            rep%B(i) = B((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
                        end do
                    case (3) !cells with id i*i (right leg)
                        do i=1, _SWE_PATCH_ORDER
                            rep%H(_SWE_PATCH_ORDER + 1 - i) = H(i*i)
                            rep%HU(_SWE_PATCH_ORDER + 1 - i) = HU(i*i)
                            rep%HV(_SWE_PATCH_ORDER + 1 - i) = HV(i*i)
                            rep%B(_SWE_PATCH_ORDER + 1 - i) = B(i*i)
                        end do
                    end select

                end associate

#if defined(_SWE_DG)			


!permute dofs for edge update
                    do k=0,_SWE_DG_ORDER
                       do j=0,_SWE_DG_ORDER
                          do i=0,_SWE_DG_ORDER-j
                             indx=(i+1)+_SWE_DG_DOFS-(_SWE_DG_ORDER+1-j)*(_SWE_DG_ORDER+2-j)/2+k*_SWE_DG_DOFS
                             indx_mir=(j+1)+_SWE_DG_DOFS-(_SWE_DG_ORDER+1-i)*(_SWE_DG_ORDER+2-i)/2+k*_SWE_DG_DOFS
                             rep%Q_DG_P(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
                             rep%Q_DG_P(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
                             rep%Q_DG_P(indx_mir)%b = element%cell%data_pers%Q_DG_P(indx)%b
                          end do
                       end do
                    end do
!                    print*
!                    print*,element%cell%data_pers%troubled
!                    print*
                    rep%troubled=element%cell%data_pers%troubled
#endif				

                
#           else
				rep%Q(1) = Q(1)
			_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q out: ", rep%Q
#           endif
		end function

		subroutine skeleton_array_op(traversal, grid, edges, rep1, rep2, update1, update2)
			type(t_swe_euler_timestep_traversal), intent(in)				:: traversal
			type(t_grid_section), intent(in)							    :: grid
			type(t_edge_data), intent(in)								    :: edges(:)
			type(num_cell_rep), intent(in)									:: rep1(:), rep2(:)
			type(num_cell_update), intent(out)								:: update1(:), update2(:)

            integer (kind = GRID_SI)                                        :: i

            do i = 1, size(edges)
                call skeleton_scalar_op(traversal, grid, edges(i), rep1(i), rep2(i), update1(i), update2(i))
            end do
		end subroutine

		subroutine skeleton_scalar_op(traversal, grid, edge, rep1, rep2, update1, update2)
			type(t_swe_euler_timestep_traversal), intent(in)				:: traversal
			type(t_grid_section), intent(in)							    :: grid
			type(t_edge_data), intent(in)								    :: edge
			type(num_cell_rep), intent(in)									:: rep1, rep2
			type(num_cell_update), intent(out)								:: update1, update2
			integer                                                         :: i

			_log_write(6, '(3X, A)') "swe skeleton op:"

#if defined(_SWE_DG)

                         update1%Q_DG_P(:)%h   =rep2%Q_DG_P(:)%h
                         update1%Q_DG_P(:)%p(1)=rep2%Q_DG_P(:)%p(1)
                         update1%Q_DG_P(:)%p(2)=rep2%Q_DG_P(:)%p(2)
                         update1%Q_DG_P(:)%b   =rep2%Q_DG_P(:)%b
                         update1%troubled      =rep2%troubled
                         update2%Q_DG_P(:)%h   =rep1%Q_DG_P(:)%h
                         update2%Q_DG_P(:)%p(1)=rep1%Q_DG_P(:)%p(1)
                         update2%Q_DG_P(:)%p(2)=rep1%Q_DG_P(:)%p(2)
                         update2%Q_DG_P(:)%b   =rep1%Q_DG_P(:)%b
                         update2%troubled      =rep1%troubled
                         ! print*
                         ! print*, update1%troubled
                         ! print*, update2%troubled                
                         ! print*

#endif


#			if defined (_SWE_LF) || defined (_SWE_LF_BATH) || defined (_SWE_LLF) || defined (_SWE_LLF_BATH) &&  !defined (_SWE_DG)
			_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q 1 in: ", rep1%Q
			_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q 2 in: ", rep2%Q
   call compute_lf_flux(edge%transform_data%normal, rep1%Q(1), rep2%Q(1), update1%flux(1), update2%flux(1))
#           elif defined (_SWE_PATCH)
                ! invert values in edges
                ! cells are copied in inverse order because the neighbor 
                ! ghost cells will have a mirrored numbering! See a (poorly-drawn) example:
                !         ___
                !  /\3   1\  |
                ! /  \2   2\ |
                !/____\1   3\|
                
                do i=1, _SWE_PATCH_ORDER
                    update1%H(i) = rep2%H(_SWE_PATCH_ORDER + 1 - i)
                    update1%HU(i) = rep2%HU(_SWE_PATCH_ORDER + 1 - i)
                    update1%HV(i) = rep2%HV(_SWE_PATCH_ORDER + 1 - i)
                    update1%B(i) = rep2%B(_SWE_PATCH_ORDER + 1 - i)
                    
                    
                    update2%H(i) = rep1%H(_SWE_PATCH_ORDER + 1 - i)
                    update2%HU(i) = rep1%HU(_SWE_PATCH_ORDER + 1 - i)
                    update2%HV(i) = rep1%HV(_SWE_PATCH_ORDER + 1 - i)
                    update2%B(i) = rep1%B(_SWE_PATCH_ORDER + 1 - i)
                end do
#           else
				call compute_geoclaw_flux(edge%transform_data%normal, rep1%Q(1), rep2%Q(1), update1%flux(1), update2%flux(1))
#			endif

			_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "flux 1 out: ", update1%flux
			_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "flux 2 out: ", update2%flux

		end subroutine

		subroutine bnd_skeleton_array_op(traversal, grid, edges, rep, update)
			type(t_swe_euler_timestep_traversal), intent(in)				:: traversal
			type(t_grid_section), intent(in)							    :: grid
			type(t_edge_data), intent(in)								    :: edges(:)
			type(num_cell_rep), intent(in)									:: rep(:)
			type(num_cell_update), intent(out)								:: update(:)

            integer (kind = GRID_SI)                                        :: i

            do i = 1, size(edges)
                call bnd_skeleton_scalar_op(traversal, grid, edges(i), rep(i), update(i))
            end do
		end subroutine

		subroutine bnd_skeleton_scalar_op(traversal, grid, edge, rep, update)
			type(t_swe_euler_timestep_traversal), intent(in)				:: traversal
			type(t_grid_section), intent(in)							    :: grid
			type(t_edge_data), intent(in)								    :: edge
			type(num_cell_rep), intent(in)									:: rep
			type(num_cell_update), intent(out)								:: update
			type(num_cell_update) 								:: temp_update
 
			type(t_state)													:: bnd_rep
                        type(t_update)													:: bnd_flux
			integer 														:: i,j,k,indx,indx_mir
                        real(kind=GRID_SR) :: length_flux
            !SLIP: reflect momentum at normal
			!bnd_rep = t_state(rep%Q(1)%h, rep%Q(1)%p - dot_product(rep%Q(1)%p, edge%transform_data%normal) * edge%transform_data%normal, rep%Q(1)%b)

            !NOSLIP: invert momentum (stable)
			!bnd_rep = t_state(rep%Q(1)%h, -rep%Q(1)%p, rep%Q(1)%b)

#                       ifdef _SWE_DG


                       temp_update%Q_DG_P(:)%h = rep%Q_DG_P(:)%h
                       temp_update%Q_DG_P(:)%b = rep%Q_DG_P(:)%b
                       update%troubled=rep%troubled
                       update%bnd=.true.

                       do i=1,(_SWE_DG_ORDER+1)*(_SWE_DG_DOFS)
                          !reflecting boundary
                          length_flux = dot_product(rep%Q_DG_P(i)%p, edge%transform_data%normal)/NORM2(edge%transform_data%normal)

                            temp_update%Q_DG_P(i)%p(1) = rep%Q_DG_P(i)%p(1)-2.0_GRID_SR*length_flux*edge%transform_data%normal(1)/NORM2(edge%transform_data%normal)
                            temp_update%Q_DG_P(i)%p(2) = rep%Q_DG_P(i)%p(2)-2.0_GRID_SR*length_flux*edge%transform_data%normal(2)/NORM2(edge%transform_data%normal)

                            !simple outflowing boundary
!                            temp_update%Q_DG_P(i)%p(1) =  rep%Q_DG_P(i)%p(1)
!                            temp_update%Q_DG_P(i)%p(2) =  rep%Q_DG_P(i)%p(2)

 
                           !zero vel bnd
!                            update%Q_DG_P(i)%p(1) =  0
!                            update%Q_DG_P(i)%p(2) =  0


                         end do

!                         if(edge%transform_data%normal(1).eq.-1.0_GRID_SR .and. edge%transform_data%normal(2).eq.0) then
!                            temp_update%Q_DG_P(:)%h = 8.0_GRID_SR
!                         end if

                         do k=0,_SWE_DG_ORDER
                            do j=0,_SWE_DG_ORDER
                               do i=0,_SWE_DG_ORDER-j
                                  indx=(i+1)+_SWE_DG_DOFS-(_SWE_DG_ORDER+1-j)*(_SWE_DG_ORDER+2-j)/2+k*_SWE_DG_DOFS
                                  indx_mir=(j+1)+_SWE_DG_DOFS-(_SWE_DG_ORDER+1-i)*(_SWE_DG_ORDER+2-i)/2+k*_SWE_DG_DOFS
                                  update%Q_DG_P(indx_mir)%h = temp_update%Q_DG_P(indx)%h
                                  update%Q_DG_P(indx_mir)%p = temp_update%Q_DG_P(indx)%p
                                  update%Q_DG_P(indx_mir)%b = temp_update%Q_DG_P(indx)%b
                               end do
                            end do
                         end do

#                       endif
			!OUTFLOW: copy values
			bnd_rep = rep%Q(1)


#			if defined (_LF_FLUX) || defined (_LF_BATH_FLUX) || defined (_LLF_FLUX) || defined (_LLF_BATH_FLUX)
				call compute_lf_flux(edge%transform_data%normal, rep%Q(1), bnd_rep, update%flux(1), bnd_flux)
#           elif defined (_SWE_PATCH)
                ! boundary conditions on ghost cells
                update%H = rep%H
                update%HU = rep%HU
                update%HV = rep%HV
                update%B = rep%B
#           else

                call compute_geoclaw_flux(edge%transform_data%normal, rep%Q(1), bnd_rep, update%flux(1), bnd_flux)
#			endif

		end subroutine

		subroutine cell_update_op(traversal, section, element, update1, update2, update3)
			type(t_swe_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)						:: element
			type(num_cell_update), intent(inout)						:: update1, update2, update3
#           if defined (_SWE_PATCH)
                integer                                                         :: i, j, ind
                type(num_cell_update)                                           :: tmp !> ghost cells in correct order 
                real(kind = GRID_SR)                                            :: volume, edge_lengths(3), maxWaveSpeed, dQ_max_norm, dt_div_volume
                real(kind = GRID_SR), DIMENSION(_SWE_PATCH_ORDER_SQUARE)                :: dQ_H, dQ_HU, dQ_HV !> deltaQ, used to compute cell updates
                real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)           :: hL, huL, hvL, bL
                real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)           :: hR, huR, hvR, bR
                real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)           :: upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR
                real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE,2,2)       :: transf

#               if !defined (_SWE_USE_PATCH_SOLVER)
                    type(t_state), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE)              :: edges_a, edges_b
                    type(t_update)                                              :: update_a, update_b
                    real(kind = GRID_SR), dimension(2,3)                        :: normals
#               endif
#if defined(_SWE_DG)
                         type(t_update),  DIMENSION(_SWE_DG_DOFS*(_SWE_DG_ORDER+1))			:: flux1, flux2, flux3
                         real(kind= GRID_SR) :: dt,dx
#endif

                !DIR$ ASSUME_ALIGNED hL: 64
                !DIR$ ASSUME_ALIGNED hR: 64
                !DIR$ ASSUME_ALIGNED huL: 64
                !DIR$ ASSUME_ALIGNED huR: 64
                !DIR$ ASSUME_ALIGNED hvL: 64
                !DIR$ ASSUME_ALIGNED hvR: 64
                !DIR$ ASSUME_ALIGNED bL: 64
                !DIR$ ASSUME_ALIGNED bR: 64
                !DIR$ ASSUME_ALIGNED upd_hL: 64
                !DIR$ ASSUME_ALIGNED upd_hR: 64
                !DIR$ ASSUME_ALIGNED upd_huL: 64
                !DIR$ ASSUME_ALIGNED upd_huR: 64
                !DIR$ ASSUME_ALIGNED upd_hvL: 64
                !DIR$ ASSUME_ALIGNED upd_hvR: 64
                !DIR$ ASSUME_ALIGNED transf: 64
                
#               if !defined(_SWE_USE_PATCH_SOLVER)
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
#               endif

                
                if (element%cell%geometry%i_plotter_type > 0) then ! if orientation = forward, reverse updates
                    tmp=update1
                    update1=update3
                    update3=tmp
                end if
                
                ! init some variables
                dQ_H = 0.0_GRID_SR
                dQ_HU = 0.0_GRID_SR
                dQ_HV = 0.0_GRID_SR
                maxWaveSpeed = 0.0_GRID_SR
                volume = cfg%scaling * cfg%scaling * element%cell%geometry%get_volume() / (_SWE_PATCH_ORDER_SQUARE)
                dt_div_volume = section%r_dt / volume
                edge_lengths = cfg%scaling * element%cell%geometry%get_edge_sizes() / _SWE_PATCH_ORDER
                
                associate(data => element%cell%data_pers, geom => SWE_PATCH_geometry)

#                               if defined (_SWE_DG)

                      ! print*
                      ! print*,data%troubled
                      ! print*,update1%troubled
                      ! print*,update2%troubled
                      ! print*,update3%troubled

                  if(data%troubled.le.0) then
                     if(.not.(update1%troubled.le.0 .and. update2%troubled.le.0 .and. update3%troubled.le.0)) then
                       data%troubled=2
                    else

                                call compute_flux_pred(ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%edges(3)%normal,data%Q_DG_P,update1%Q_DG_P,flux1)
                                maxWaveSpeed = max(maxWaveSpeed, maxval(flux1%max_wave_speed))
                                call compute_flux_pred(ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%edges(2)%normal,data%Q_DG_P,update2%Q_DG_P,flux2)
                                maxWaveSpeed = max(maxWaveSpeed, maxval(flux2%max_wave_speed))
                                call compute_flux_pred(ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%edges(1)%normal,data%Q_DG_P,update3%Q_DG_P,flux3)
                                maxWaveSpeed = max(maxWaveSpeed, maxval(flux3%max_wave_speed))

                                call dg_solver(element,flux1,flux2,flux3,section%r_dt)
                             end if
                          end if
#                               endif                
                    ! copy cell values to arrays edges_a and edges_b
                    ! obs: cells with id > number of cells are actually ghost cells and come from edges "updates"
                    ! see t_SWE_PATCH_geometry for an explanation about ghost cell ordering
                    do i=1, _SWE_PATCH_NUM_EDGES, _SWE_PATCH_SOLVER_CHUNK_SIZE ! i -> init of chunk
                        
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
                        
                        do j=1, _SWE_PATCH_SOLVER_CHUNK_SIZE ! j -> position inside chunk
                            ind = i + j - 1 ! actual index
                           
                            ! don't go outside array limits
                            if (ind > _SWE_PATCH_NUM_EDGES) then
                                exit
                            end if
                            
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
                            else if (geom%edges_a(ind) <= _SWE_PATCH_ORDER_SQUARE + 2*_SWE_PATCH_ORDER) then
                                hL(j) = update2%H(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                                huL(j) = update2%HU(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                                hvL(j) = update2%HV(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                                bL(j) = update2%B(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
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
                            
                            !copy transformation matrices
                            transf(j,:,:) = geom%transform_matrices(geom%edges_orientation(ind),:,:,element%cell%geometry%i_plotter_type)
                        end do

                        ! compute net_updates -> solve Riemann problems within chunk
#                       if defined (_SWE_USE_PATCH_SOLVER)
#                           if defined(_FWAVE_FLUX) || defined(_AUG_RIEMANN_FLUX)
                                call compute_updates_simd(transf, hL, huL, hvL, bL, hR, huR, hvR, bR, upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR, maxWaveSpeed)
#                           elif defined(_HLLE_FLUX)
                                call compute_updates_hlle_simd(transf, hL, huL, hvL, bL, hR, huR, hvR, bR, upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR, maxWaveSpeed)
#                           else
                                ! this should never happen -> SCons rules should avoid this before compiling
#                               error "No valid SWE solver for patches/simd implementation has been defined!"
                                    ! ! !!! !!!!!!!print *, 
#                           endif
#                       else                    
                            ! using original geoclaw solver
                            do j=1, _SWE_PATCH_SOLVER_CHUNK_SIZE
                                ind = i + j - 1 ! actual index
                                
                                ! don't go outside array limits
                                if (ind > _SWE_PATCH_NUM_EDGES) then
                                    exit
                                end if
                                
                                edges_a(j)%h = hL(j)
                                edges_a(j)%p(1) = huL(j)
                                edges_a(j)%p(2) = hvL(j)
                                edges_a(j)%b = bL(j)
                                edges_b(j)%h = hR(j)
                                edges_b(j)%p(1) = huR(j)
                                edges_b(j)%p(2) = hvR(j)
                                edges_b(j)%b = bR(j)
                                call compute_geoclaw_flux(normals(:,geom%edges_orientation(ind)), edges_a(j), edges_b(j), update_a, update_b)
                                upd_hL(j) = update_a%h
                                upd_huL(j) = update_a%p(1)
                                upd_hvL(j) = update_a%p(2)
                                upd_hR(j) = update_b%h
                                upd_huR(j) = update_b%p(1)
                                upd_hvR(j) = update_b%p(2)
                                maxWaveSpeed = max(maxWaveSpeed, update_a%max_wave_speed)
                            end do
#                       endif

                        ! compute dQ
                        do j=1, _SWE_PATCH_SOLVER_CHUNK_SIZE
                            ind = i + j - 1 ! actual index
                            
                            ! don't go outside array limits
                            if (ind > _SWE_PATCH_NUM_EDGES) then
                                exit
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
                    end do 

                    ! if land is flooded, init water height to dry tolerance and
                    ! velocity to zero
                    where (data%H < data%B + cfg%dry_tolerance .and. data%H + dQ_H > data%B + cfg%dry_tolerance)
                        data%H = data%B + cfg%dry_tolerance
                        data%HU = 0.0_GRID_SR
                        data%HV = 0.0_GRID_SR
                    end where

                    ! update unknowns
                    data%H = data%H + dQ_H * (-dt_div_volume)
                    data%HU = data%HU + dQ_HU * (-dt_div_volume)
                    data%HV = data%HV + dQ_HV * (-dt_div_volume)
                    
                    ! if the water level falls below the dry tolerance, set water surface to 0 and velocity to 0
                    where (data%H < data%B + cfg%dry_tolerance) 
!                        data%H = min(data%B, 0.0_GRID_SR)
                        data%H = data%B
                        data%HU = 0.0_GRID_SR
                        data%HV = 0.0_GRID_SR
                    end where
                    
                    !set refinement condition -> Here I am using the same ones as in the original no-patches implementation, but considering only the max value.
                    element%cell%geometry%refinement = 0
                    dQ_max_norm = maxval(abs(dQ_H))

                    if (element%cell%geometry%i_depth < cfg%i_max_depth .and. dQ_max_norm > 5.0_GRID_SR * cfg%scaling * get_edge_size(cfg%i_max_depth)) then
                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
                    else if (element%cell%geometry%i_depth > cfg%i_min_depth .and. dQ_max_norm < 5.0_GRID_SR * cfg%scaling * get_edge_size(cfg%i_max_depth) / 8.0_SR) then
                        element%cell%geometry%refinement = -1
                    endif

                    ! compute next dt
                    !!! !!!!!!!print*
                    !! !!!!!!!print*,maxWaveSpeed
                    section%r_dt_new = min(section%r_dt_new, volume / (edge_lengths(2) * maxWaveSpeed) )
                    maxWaveSpeed=0
                    !!! !!!!!!!print*,section%r_dt_new
#if defined (_SWE_DG)
                
                !detect troubled cell
                print*,data%troubled
                call data%set_troubled(cfg%dry_tolerance)
                
                !! !!!!!!!print*,"predict"
                if(data%troubled.le.0) then
                   call dg_predictor(element,section%r_dt)
                end if
                !! !!!!!!!print*,"done"
#endif

                end associate

#           else
			!local variables

			type(t_state)   :: dQ(_SWE_CELL_SIZE)

			call volume_op(element%cell%geometry, traversal%i_refinements_issued, element%cell%geometry%i_depth, &
                element%cell%geometry%refinement, section%r_dt_new, dQ, [update1%flux, update2%flux, update3%flux], section%r_dt)

			!if land is flooded, init water height to dry tolerance and velocity to 0
			if (element%cell%data_pers%Q(1)%h < element%cell%data_pers%Q(1)%b + cfg%dry_tolerance .and. dQ(1)%h > 0.0_GRID_SR) then
                element%cell%data_pers%Q(1)%h = element%cell%data_pers%Q(1)%b + cfg%dry_tolerance
                element%cell%data_pers%Q(1)%p = [0.0_GRID_SR, 0.0_GRID_SR]
                !! ! !!! !!!!!!!print '("Wetting:", 2(X, F0.0))', cfg%scaling * element%transform_data%custom_data%offset + cfg%offset
            end if

            call gv_Q%add(element, dQ)

			!if the water level falls below the dry tolerance, set water surface to 0 and velocity to 0
			if (element%cell%data_pers%Q(1)%h < element%cell%data_pers%Q(1)%b + cfg%dry_tolerance) then
                element%cell%data_pers%Q(1)%h = min(element%cell%data_pers%Q(1)%b, 0.0_GRID_SR)
                element%cell%data_pers%Q(1)%p = [0.0_GRID_SR, 0.0_GRID_SR]
           end if
#          endif

		end subroutine

		subroutine cell_last_touch_op(traversal, section, cell)
			type(t_swe_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_cell_data_ptr), intent(inout)				:: cell
			real(kind = GRID_SR)							    :: b_norm

			b_norm = minval(abs(cell%data_pers%Q%h - cell%data_pers%Q%b))

			!refine also on the coasts
			if (cell%geometry%i_depth < cfg%i_max_depth .and. b_norm < 20.0_GRID_SR) then
				!cell%geometry%refinement = 1
				!traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
			else if (b_norm < 100.0_GRID_SR) then
				!cell%geometry%refinement = max(cell%geometry%refinement, 0)
			endif
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************
#if defined(_SWE_DG)
                subroutine dg_solver(element,update1,update2,update3,dt)
			type(t_element_base), intent(in)						:: element
                        type(t_update), dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1)),intent(inout)  	:: update1, update2, update3
                        real(kind=GRID_SR)                                                              :: q(_SWE_DG_DOFS,3), q_p((_SWE_DG_ORDER+1)*_SWE_DG_DOFS,3),q_temp(_SWE_DG_DOFS,3) ,nF1(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3),nF2(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3) ,nF3(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3)
                        real(kind=GRID_SR)                                                              :: dt
                        integer :: i,j,k,inx_l,inx_m
                        real(kind=GRID_SR) :: dx
                        real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3)                                                              :: f1_hat,f1_s,f2_hat,f2_s,S_s
                        real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3)                                                              :: volume_flux1,volume_flux2,orig_flux1,orig_flux2,bnd_flux_contrib,source

                        associate(Q_DG => element%cell%data_pers, jacobian_inv =>  ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%jacobian_inv) 

                        call Q_DG%dofs_to_vec_dg(q)
                        call Q_DG%dofs_to_vec_dg_p(q_p)

                        dx=element%transform_data%custom_data%scaling
                        do i=1,_SWE_DG_DOFS*(_SWE_DG_ORDER+1)
                           nF1(i,1)=update1(i)%h
                           nF1(i,2)=update1(i)%p(1)
                           nF1(i,3)=update1(i)%p(2)
                           
                           NF2(i,1)=update2(i)%h
                           NF2(i,2)=update2(i)%p(1)
                           NF2(i,3)=update2(i)%p(2)
                           
                           NF3(i,1)=update3(i)%h
                           NF3(i,2)=update3(i)%p(1)
                           NF3(i,3)=update3(i)%p(2)
                        end do

                        f1_hat = f1(q_p,(_SWE_DG_ORDER+1)*_SWE_DG_DOFS)
                        f2_hat = f2(q_p,(_SWE_DG_ORDER+1)*_SWE_DG_DOFS)
                        
                        f1_s=(jacobian_inv(1,1)*f1_hat+jacobian_inv(1,2)*f2_hat)/NORM2(jacobian_inv(:,1))
                        f2_s=(jacobian_inv(2,1)*f1_hat+jacobian_inv(2,2)*f2_hat)/NORM2(jacobian_inv(:,1))
                        S_s = S(q_p,Q_DG%Q_DG_P(:)%b_x_g,Q_DG%Q_DG_P(:)%b_y_g,(_SWE_DG_ORDER+1)*_SWE_DG_DOFS)                                                  
                        volume_flux1=matmul(s_k_x,f1_s)
                        volume_flux2=matmul(s_k_y,f2_s)
                        source= matmul(s_m,S_s)

                        ! do j=0,_SWE_DG_ORDER
                        !    source(j+1,:)=0
                        !    source(_SWE_DG_DOFS-(j+1)*(j)/2,:)=0
                        !    source(_SWE_DG_DOFS-(j+2)*(j+1)/2+1,:)=0
                        ! end do

                        bnd_flux_contrib=matmul(b_m_1,f2_s)+matmul(b_m_3,f1_s)-matmul(b_m_2,(f1_s+f2_s))/sqrt(2.0_GRID_SR)

                        q_temp = ( matmul(b_m_1,nF1) &
                             + matmul(b_m_2,nF2)  &
                             + matmul(b_m_3,nF3)  &
                             - volume_flux1          &
                             - volume_flux2          &
                             - bnd_flux_contrib)* dt/dx   &
                             - source * dt/dx

                            !!print*, "orig ",q
                            !!print*, "dx ",dx
                            !!print*, "dt ",dt
                            !!print*, "bnd_flux ",matmul(b_m_1,nF1)+matmul(b_m_2,nF2) + matmul(b_m_3,nF3)  
                        
                            !!print*, "elt flux with boundary" ,volume_flux1 + volume_flux2 +bnd_flux_contrib
                        
                            !!print*, "elt flux ",volume_flux1 + volume_flux2 
                            
                            !!print*, "bnd flux contrib" , bnd_flux_contrib
                            !!print*,
                            !!print*, "Source ", matmul(s_m,S_s) * dt
                            !!print*,
                            !!print*, "Without Source",( matmul(b_m_1,nF1)+matmul(b_m_2,nF2) + matmul(b_m_3,nF3)  - volume_flux1  -volume_flux2  -bnd_flux_contrib) * dt/dx
                            !!print*,
                            !!print*,"b_x_g"
                            !!print*,Q_DG%Q_DG_P(:)%b_x_g
                            !!print*,"b_y_g"
                            !!print*,Q_DG%Q_DG_P(:)%b_y_g

                            !!print*,
                              !!print*, "bnd_flux1 ",matmul(b_m_1,nF1) 
                              !!print*,
                              !!print*, "bnd_flux2 ", matmul(b_m_2,nF2) 
                              !!print*,
                              !!print*, "bnd_flux3 ",matmul(b_m_3,nF3)  
                              !!print*,
                              !!print*, "num_flux1 ",nF1 
                              !!print*,
                              !!print*, "num_flux2 ",nF2 
                              !!print*,
                              !!print*, "num_flux3 ",nF3
                              !!print*,
                              !!print*, "elt_flux1 ",f1_s
                              !!print*,
                              !!print*, "elt_flux2 ",f2_s
                              !!print*,
                             !!print*, "elt_source",S_s
                             !!print*,
                             !!print*,"q_temp" , q_temp
                             !!print*
                             
                        call lusolve(s_m_lu, _SWE_DG_DOFS, s_m_lu_pivot,q_temp(:,1))
                        call lusolve(s_m_lu, _SWE_DG_DOFS, s_m_lu_pivot,q_temp(:,2))
                        call lusolve(s_m_lu, _SWE_DG_DOFS, s_m_lu_pivot,q_temp(:,3))
                        ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !!! !!!!!!!!print*, q_temp

                        q = q - q_temp
                        call Q_DG%vec_to_dofs_dg(q)

                      end associate

                end subroutine dg_solver
               
                subroutine dg_predictor(element,dt)
			class(t_element_base), intent(in)						:: element
                        integer                                                                         :: i,j
                        real(kind=GRID_SR)                                                              :: dt
                        real(kind=GRID_SR)                                                              :: q_0(_SWE_DG_DOFS,3)
                        real(kind=GRID_SR)                                                              :: q_i(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3)
                        real(kind=GRID_SR)                                                              :: q_temp(_SWE_DG_DOFS*_SWE_DG_ORDER,3)
                        real(kind=GRID_SR)                                                              :: q1_temp(_SWE_DG_DOFS*_SWE_DG_ORDER,3),q2_temp(_SWE_DG_DOFS*_SWE_DG_ORDER,3),q3_temp(_SWE_DG_DOFS*_SWE_DG_ORDER,3),q4_temp(_SWE_DG_DOFS*_SWE_DG_ORDER,3)
                        real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3)                                                              :: f1_hat,f1_s,f2_hat,f2_s,S_s
                        real(kind=GRID_SR) :: dx
                        real(kind=GRID_SR) :: epsilon
                        ! ! ! ! ! ! ! ! ! ! ! ! ! ! !!! !!!!!!!!print*,"predict"
                        associate(Q_DG =>element%cell%data_pers%Q_DG ,Q_DG_P =>element%cell%data_pers%Q_DG_P)

                        call element%cell%data_pers%dofs_to_vec_dg(q_0)

                        q_i(1:_SWE_DG_DOFS,:) = q_0
                        q_i(1+_SWE_DG_DOFS:(_SWE_DG_ORDER+1)*_SWE_DG_DOFS,:) = 0
                        
                        dx=element%transform_data%custom_data%scaling
                        associate(jacobian_inv =>  ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%jacobian_inv)
                        epsilon=1.0_GRID_SR                          
                        i=0

                        do while(epsilon > 1.0e-12_GRID_SR)
                           i=i+1
                           ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !!! !!!!!!!!print*,i
                           !                        do i = 1,_SWE_DG_DOFS*_SWE_DG_DOFS

                           f1_hat = f1_not_drying(q_i,(_SWE_DG_ORDER+1)*_SWE_DG_DOFS)
                           f2_hat = f2_not_drying(q_i,(_SWE_DG_ORDER+1)*_SWE_DG_DOFS)
                           
                           f1_s=(jacobian_inv(1,1)*f1_hat+jacobian_inv(1,2)*f2_hat)/NORM2(jacobian_inv(:,1))
                           f2_s=(jacobian_inv(2,1)*f1_hat+jacobian_inv(2,2)*f2_hat)/NORM2(jacobian_inv(:,1))

                           S_s =S(q_i,Q_DG_P(:)%b_x_g,Q_DG_P(:)%b_y_g,(_SWE_DG_ORDER+1)*_SWE_DG_DOFS)*dt/dx
                           
                           f1_s = f1_s*dt/dx
                           f2_s = f2_s*dt/dx

                           q_temp=  matmul(st_m,S_s) &
                                  - matmul(st_w_k_t_1_0,q_i(1:_SWE_DG_DOFS,:))&
                                  - matmul(st_k_x,f1_s)&
                                  - matmul(st_k_y,f2_s)
                           ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !!! !!!!!!!!print*,q_temp
                           call lusolve(st_w_k_t_1_1_lu,_SWE_DG_ORDER*_SWE_DG_DOFS,st_w_k_t_1_1_lu_pivot ,q_temp(:,1))
                           call lusolve(st_w_k_t_1_1_lu,_SWE_DG_ORDER*_SWE_DG_DOFS,st_w_k_t_1_1_lu_pivot ,q_temp(:,2))
                           call lusolve(st_w_k_t_1_1_lu,_SWE_DG_ORDER*_SWE_DG_DOFS,st_w_k_t_1_1_lu_pivot ,q_temp(:,3))
                           
                           if(.not.(all(q_i(1+_SWE_DG_DOFS:(_SWE_DG_ORDER+1)*_SWE_DG_DOFS,:) .eq. 0))) then
                              epsilon=NORM2(q_temp-q_i(1+_SWE_DG_DOFS:(_SWE_DG_ORDER+1)*_SWE_DG_DOFS,:))/NORM2(q_i(1+_SWE_DG_DOFS:(_SWE_DG_ORDER+1)*_SWE_DG_DOFS,:))
                           else
                              if(all(q_temp .eq. 0)) then
                              epsilon=0
                              else 
                              epsilon=1
                              end if
                           end if


                        !     if(i > (_SWE_DG_ORDER+1)*_SWE_DG_DOFS) then
                           !print*,i
                           !print*
                           !print*,"epsilon " , epsilon
                           !print*,"jacobian ", jacobian_inv
                           !print*,"scal ", dt/dx
                           !print*, dt
                           !print*, dx
                           !print*
                           !print*,"q_0 ", q_0
                           
                           !print*                           
                           !print*,"st_k_w_t ", matmul(st_w_k_t_1_0,q_i(1:_SWE_DG_DOFS,:))
                           !print*,"st_k_x ",matmul(st_k_x,f1_s)
                           !print*,"st_k_y ",matmul(st_k_y,f2_s)
                           !print*,"st_m ",matmul(st_m,S_s)

                           !print*,"Steady ",matmul(st_k_x,f1_s)+matmul(st_k_y,f2_s)-matmul(st_m,S_s)

                           !print*,"q_temp " , matmul(st_m,S_s)- matmul(st_w_k_t_1_0,q_i(1:_SWE_DG_DOFS,:)) - matmul(st_k_x,f1_s) - matmul(st_k_y,f2_s)
                           

                           !print*
                           
                           !print*
                           !print*,Q_DG(:)%b
                           !print*,Q_DG(:)%b_x_g
                           !print*,Q_DG(:)%b_y_g
                           !print*
                           !print*,"S_s 0 ", S_s(1:_SWE_DG_DOFS,:)
                           !print*,"f1_s 0", f1_s(1:_SWE_DG_DOFS,:)
                           !print*,"f2_s 0", f2_s(1:_SWE_DG_DOFS,:)


                           !print*,"S_s ", S_s(1+_SWE_DG_DOFS:(_SWE_DG_ORDER+1)*_SWE_DG_DOFS,:)
                           !print*,"f1_s ", f1_s(1+_SWE_DG_DOFS:(_SWE_DG_ORDER+1)*_SWE_DG_DOFS,:)
                           !print*,"f2_s ", f2_s(1+_SWE_DG_DOFS:(_SWE_DG_ORDER+1)*_SWE_DG_DOFS,:)
                           !print*,"q_i ", q_i(1+_SWE_DG_DOFS:(_SWE_DG_ORDER+1)*_SWE_DG_DOFS,:)
                           !print*
                           !print*,"q_temp lu", q_temp
                           
                           
                        !                         end if
                          
                            q_i(1+_SWE_DG_DOFS:(_SWE_DG_ORDER+1)*_SWE_DG_DOFS,:) = q_temp
                         end do


                        end associate
                                                   
                        call element%cell%data_pers%vec_to_dofs_dg_p(q_i)

                        end  associate
                        ! ! ! ! ! ! ! ! ! ! ! ! ! ! !!! !!!!!!!!print*,"predict done"
                end subroutine dg_predictor

                function f1_not_drying(q,N)
                  integer,intent(in) :: N
                  real(kind=GRID_SR) ,intent(in) ::q(N,3)
                  real(kind=GRID_SR)             ::f1_not_drying(N,3)
                  integer                        ::i
                  
                  f1_not_drying(:,1) = q(:,2)

                  do i=1,N
                     if(q(i,1)==0.0_GRID_SR) then
                        f1_not_drying(i,1) = 0
                        f1_not_drying(i,2) = 0
                        f1_not_drying(i,3) = 0
                     else
                        f1_not_drying(i,2) = q(i,2)**2/q(i,1) + 0.5_GRID_SR * g * q(i,1)**2
                        f1_not_drying(i,3) = q(i,2)*q(i,3)/q(i,1)
                     end if
                  end do
                end function f1_not_drying

                function f2_not_drying(q,N)
                  integer,intent(in) :: N
                  real(kind=GRID_SR) ,intent(in) ::q(N,3)
                  real(kind=GRID_SR)             ::f2_not_drying(N,3)
                  integer                        ::i
                  
                  f2_not_drying(:,1) = q(:,3)

                  do i=1,N
                     if(q(i,1)==0.0_GRID_SR) then
                        f2_not_drying(i,1) = 0
                        f2_not_drying(i,2) = 0
                        f2_not_drying(i,3) = 0
                     else
                        f2_not_drying(i,3) = q(i,3)**2/q(i,1) + 0.5_GRID_SR * g * q(i,1)**2
                        f2_not_drying(i,2) = q(i,2)*q(i,3)/q(i,1)
                     end if
                  end do
                end function f2_not_drying

                function f1(q,N)
                  integer :: N
                  real(kind=GRID_SR) ,intent(in) ::q(N,3)
                  real(kind=GRID_SR)             ::f1(N,3)
                  integer                        ::i

                  f1(:,1) = q(:,2)
                  do i=1,N
                     if(q(i,1)<=cfg%dry_tolerance) then
                        f1(i,1) = 0
                        f1(i,2) = 0
                        f1(i,3) = 0
                     else
                        f1(i,2) = q(i,2)**2/q(i,1) + 0.5_GRID_SR * g * q(i,1)**2
                        f1(i,3) = q(i,2)*q(i,3)/q(i,1)
                     end if
                  end do

                end function f1
                
                function f2(q,N)
                  integer :: N
                  real(kind=GRID_SR) ,intent(in) ::q(N,3)
                  real(kind=GRID_SR)             ::f2(N,3)
                  integer                        ::i

                  f2(:,1) = q(:,3)
                  do i=1,N
                     if(q(i,1)<=cfg%dry_tolerance) then
                        f2(i,1) = 0
                        f2(i,2) = 0
                        f2(i,3) = 0
                     else
                        f2(i,2) = q(i,2)*q(i,3)/q(i,1)
                        f2(i,3) = q(i,3)**2/q(i,1) + 0.5_GRID_SR * g * q(i,1)**2
                     end if
                  end do

                end function f2

                function S(q,b_x_g,b_y_g,N)
                  integer ::N,i
                  real(kind=GRID_SR) ,intent(in) ::q(N,3), b_x_g(N), b_y_g(N)
                  real(kind=GRID_SR)             ::S(N,3)

                     S(:,1) = 0
                  do i=1,N
                     S(i,2)=-q(i,1)*b_x_g(i)
                     S(i,3)=-q(i,1)*b_y_g(i) 
                 end do
                  
                end function S



#endif              

		subroutine volume_op(cell, i_refinements_issued, i_depth, i_refinement, r_dt_new, dQ, fluxes, r_dt)

			type(fine_triangle), intent(in)				                        :: cell
			integer (kind = GRID_DI), intent(inout)							    :: i_refinements_issued
			integer (kind = BYTE), intent(in)							        :: i_depth
			integer (kind = BYTE), intent(out)							        :: i_refinement
			real(kind = GRID_SR), intent(inout)								    :: r_dt_new
			type(t_state), dimension(:), intent(out)						    :: dQ
			type(t_update), dimension(:), intent(in)						    :: fluxes
			real(kind = GRID_SR), intent(in)								    :: r_dt

			real(kind = GRID_SR)											    :: volume, dQ_norm, edge_lengths(3)
			integer (kind = BYTE)												:: i
			real (kind = GRID_SR), parameter                                    :: refinement_threshold = 5.0_SR

			_log_write(6, '(3X, A)') "swe cell update op:"
			_log_write(6, '(4X, A, 4(X, F0.3))') "edge 1 flux in:", fluxes(1)
			_log_write(6, '(4X, A, 4(X, F0.3))') "edge 2 flux in:", fluxes(2)
			_log_write(6, '(4X, A, 4(X, F0.3))') "edge 3 flux in:", fluxes(3)

			volume = cfg%scaling * cfg%scaling * cell%get_volume()
			edge_lengths = cfg%scaling * cell%get_edge_sizes()

			dQ%h = dot_product(edge_lengths, fluxes%h)
			dQ%p(1) = dot_product(edge_lengths, fluxes%p(1))
			dQ%p(2) = dot_product(edge_lengths, fluxes%p(2))
			dQ%b = 0.0_GRID_SR

			!set refinement condition

			i_refinement = 0
			dQ_norm = abs(dQ(1)%h)

			if (i_depth < cfg%i_max_depth .and. dQ_norm > refinement_threshold * cfg%scaling * get_edge_size(cfg%i_max_depth)) then
				i_refinement = 1
				i_refinements_issued = i_refinements_issued + 1_GRID_DI
			else if (i_depth > cfg%i_min_depth .and. dQ_norm < refinement_threshold * cfg%scaling * get_edge_size(cfg%i_max_depth) / 8.0_SR) then
				i_refinement = -1
			endif

            !This will cause a division by zero if the wave speeds are 0.
            !Bue to the min operator, the error will not affect the time step.
            !r_dt_new = min(r_dt_new, volume / dot_product(edge_lengths, fluxes%max_wave_speed))

            r_dt_new = min(r_dt_new, volume / (edge_lengths(2) * maxval(fluxes%max_wave_speed)))


            do i = 1, _SWE_CELL_SIZE
                dQ(i)%t_dof_state = dQ(i)%t_dof_state * (-r_dt / volume)
            end do

			_log_write(6, '(4X, A, 4(X, F0.3))') "dQ out: ", dQ
		end subroutine volume_op

		!> Lax Friedrichs flux. Depending on compiler flags, the function implements
		!> the global or local variant with or without bathymetry
		subroutine compute_lf_flux(normal, QL, QR, fluxL, fluxR)
			type(t_state), intent(in)							:: QL, QR
			type(t_update), intent(out) 						:: fluxL, fluxR
			real(kind = GRID_SR), intent(in)		            :: normal(2)

			real(kind = GRID_SR)								:: vL, vR, hL, hR, alpha

#           if defined(_LF_BATH_FLUX) || defined(_LLF_BATH_FLUX)
                if (QL%h - QL%b < cfg%dry_tolerance .or. QR%h - QR%b < cfg%dry_tolerance) then
                    hL = 0.0_SR; hR = 0.0_SR
                    vL = 0.0_SR; vR = 0.0_SR

                    fluxL%max_wave_speed = 0.0_SR; fluxR%max_wave_speed = 0.0_SR
                    fluxL%h = 0.0_SR; fluxR%h = 0.0_SR
                    fluxL%p = 0.0_SR; fluxR%p = 0.0_SR

                    !This boundary treatment assumes a wall condition.
                    !For the mass flux, we choose pR := -hL * vL, vR = 0 (walls are immovable), hence hR must be infinite.
                    !For the momentum flux we choose hR := 0, vR := 0 (implying there is no hydrostatic pressure), bR := bL + hL (there is a wall to the right)

                    if (QL%h - QL%b < cfg%dry_tolerance .and. QR%h - QR%b < cfg%dry_tolerance) then
                    else if (QL%h - QL%b < cfg%dry_tolerance) then
                        hR = max(QR%h - QR%b, 0.0_SR)
                        vR = dot_product(normal, QR%p / (QR%h - QR%b))
                        fluxR%max_wave_speed = sqrt(g * hR) + abs(vR)
                        fluxR%p = -0.5_SR * vR * QR%p - 0.5_GRID_SR * g * hR * hR * normal + 0.5_SR * fluxR%max_wave_speed * QR%p
                    else if (QR%h - QR%b < cfg%dry_tolerance) then
                        hL = max(QL%h - QL%b, 0.0_SR)
                        vL = dot_product(normal, QL%p / (QL%h - QL%b))
                        fluxL%max_wave_speed = sqrt(g * hL) + abs(vL)
                        fluxL%p = 0.5_SR * vL * QL%p + 0.5_GRID_SR * g * hL * hL * normal + 0.5_SR * fluxL%max_wave_speed * QL%p
                    end if

                    return
                end if

                hL = max(QL%h - QL%b, 0.0_SR)
                hR = max(QR%h - QR%b, 0.0_SR)

                vL = dot_product(normal, QL%p / (QL%h - QL%b))
                vR = dot_product(normal, QR%p / (QR%h - QR%b))

                fluxL%max_wave_speed = sqrt(g * hL) + abs(vL)
                fluxR%max_wave_speed = sqrt(g * hR) + abs(vR)

#               if defined(_LLF_BATH_FLUX)
                    alpha = max(fluxL%max_wave_speed, fluxR%max_wave_speed)
#               else
                    alpha = 100.0_GRID_SR
#               endif

                !Except for the diffusion term, the mass flux is the standard LF flux
                fluxL%h = 0.5_GRID_SR * (hL * vL + hR * vR + alpha * (QL%h - QR%h))
                fluxR%h = -fluxL%h

                !The base momentum flux is similar to the standard LF flux.
                fluxL%p = 0.5_GRID_SR * (QL%p * vL + QR%p * vR + 0.5_GRID_SR * g * (hL * hL + hR * hR) * normal + alpha * (QL%p - QR%p))
                fluxR%p = -fluxL%p

                !The source term $\Delta x \ \Psi = $-1/2 g \ \frac{1}{2} \ (h_l + h_r) \ (b_r - b_l)$ [LeVeque] is added on both sides with a weight of 1/2.
                !This factor ensures that the method is well-balanced.
                fluxL%p = fluxL%p + 0.25_SR * g * (hL + hR) * (QR%b - QL%b) * normal
                fluxR%p = fluxR%p + 0.25_SR * g * (hL + hR) * (QR%b - QL%b) * normal
#           elif defined(_LF_FLUX) || defined(_LLF_FLUX)
                real(kind = GRID_SR), parameter					:: b = 0.0_GRID_SR     !default constant bathymetry

                hL = max(QL%h - b, 0.0_SR)
                hR = max(QR%h - b, 0.0_SR)

                vL = dot_product(normal, QL%p / (QL%h - b))
                vR = dot_product(normal, QR%p / (QR%h - b))

                fluxL%max_wave_speed = sqrt(g * (QL%h - b)) + abs(vL)
                fluxR%max_wave_speed = sqrt(g * (QR%h - b)) + abs(vR)

#               if defined(_LLF_FLUX)
                    alpha = max(fluxL%max_wave_speed, fluxR%max_wave_speed)
#               else
                    alpha = 100.0_GRID_SR
#               endif


                fluxL%h = 0.5_GRID_SR * (vL * hL + vR * hR + alpha * (QL%h - QR%h))
                fluxR%h = -fluxL%h

                fluxL%p = 0.5_GRID_SR * (vL * vL * hL + vR * vR * hR + 0.5_GRID_SR * g * (hL * hL + hR * hR) * normal + alpha * (QL%p - QR%p))
                fluxR%p = -fluxL%p

#           endif
		end subroutine

#if defined(_SWE_DG)
		subroutine compute_flux_pred(normal,QL, QR,flux)

			type(t_state),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1)), intent(in)	:: QL
                        type(t_state),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1)), intent(in)    :: QR
			type(t_update),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1)), intent(out) 						:: flux
                        type(t_update),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1)) 						:: flux_r

                        real(kind=GRID_SR),Dimension(1,3) :: f1r_s,f1l_s,f2r_s,f2l_s,f,ql_v,qr_v
			real(kind = GRID_SR), intent(in)		            :: normal(2)

                        integer ::i


#if                    !(defined(_FWAVE_FLUX)|defined(_AUG_RIEMANN_FLUX)|defined(_HLLE_FLUX))

			real(kind=GRID_SR)								:: vL, vR, alpha
                        real(kind=GRID_SR) :: max_wave_speedR


                        alpha=0

                        do i=1,_SWE_DG_DOFS*(_SWE_DG_ORDER+1) 

                        if(.not.(QL(i)%h < cfg%dry_tolerance) ) then 
                           vL = DOT_PRODUCT(normal, QL(i)%p/QL(i)%h)
                           flux(i)%max_wave_speed = sqrt(g * (QL(i)%h)) + abs(vL)
                        else
                           flux(i)%max_wave_speed = 0
                        end if

                        if(.not.(QR(i)%h < cfg%dry_tolerance) ) then 
                           vR = DOT_PRODUCT(normal, QR(i)%p/QR(i)%h)
                           max_wave_speedR = sqrt(g * (QR(i)%h)) + abs(vR)
                        else
                           max_wave_speedR = 0
                        end if

                        alpha = max(alpha,flux(i)%max_wave_speed, max_wave_speedR)

                        end do

                        do i=1,_SWE_DG_DOFS*(_SWE_DG_ORDER+1)

                           ql_v(1,1)= QL(i)%h
                           ql_v(1,2:3)= QL(i)%p
                           
                           qr_v(1,1)= QR(i)%h
                           qr_v(1,2:3)= QR(i)%p
                           
                           f1r_s=f1(qr_v,1)
                           f2r_s=f2(qr_v,1)
                           f1l_s=f1(ql_v,1)
                           f2l_s=f2(ql_v,1)
                           
                           f=0.5_GRID_SR*((f1r_s-f1l_s) * normal(1) + (f2r_s-f2l_s) *normal(2)+alpha*(ql_v-qr_v))
                           flux(i)%h=f(1,1)
                           flux(i)%p=f(1,2:3)

                           ! end if
                        end do

#else
                        real(kind=GRID_SR) :: s_s(1,3)
                        flux%h=0
                        flux%p(1)=0
                        flux%p(2)=0
                        flux_r%h=0
                        flux_r%p(1)=0
                        flux_r%p(2)=0

                        do i=1,_SWE_DG_DOFS*(_SWE_DG_ORDER+1)
                           call compute_geoclaw_flux(normal, QL(i), QR(i), flux(i), flux_r(i),.false.)
                           ql_v(1,1)= QL(i)%h
                           ql_v(1,2:3)= QL(i)%p
                           f1l_s=f1(ql_v,1)
                           f2l_s=f2(ql_v,1)
                        end do

#endif
                      end subroutine compute_flux_pred
#           endif



		!> Augmented Riemann solver from geoclaw
		subroutine compute_geoclaw_flux(normal, QL, QR, fluxL, fluxR,with_bathymetry_in)
			type(t_state), intent(in)           :: QL, QR
			type(t_update), intent(out)         :: fluxL, fluxR
			real(kind = GRID_SR), intent(in)    :: normal(2)
                        logical,optional,intent(in)         :: with_bathymetry_in

			real(kind = GRID_SR)				:: transform_matrix(2, 2)
			real(kind = GRID_SR)			    :: net_updatesL(3), net_updatesR(3), max_wave_speed
                        logical         :: with_bathymetry_local
			real(kind = GRID_SR)                :: pL(2), pR(2), hL, hR, bL, bR
                        
                        if(present(with_bathymetry_in))then
                           with_bathymetry_local=with_bathymetry_in
                        else
                           with_bathymetry_local=.true.
                        end if


			transform_matrix(1, :) = normal
			transform_matrix(2, :) = [-normal(2), normal(1)]

			pL = matmul(transform_matrix, QL%p)
			pR = matmul(transform_matrix, QR%p)
    
                        if(with_bathymetry_local) then
                           hL = QL%h-QL%b
                           hR = QR%h-QR%b
                        else
                           hL = QL%h
                           hR = QR%h
                        end if

			bL = QL%b
			bR = QR%b
			
#           if defined(_FWAVE_FLUX)
                call c_bind_geoclaw_solver(GEOCLAW_FWAVE, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_AUG_RIEMANN_FLUX)
                call c_bind_geoclaw_solver(GEOCLAW_AUG_RIEMANN, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_HLLE_FLUX)
                call compute_updates_hlle_single(hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, net_updatesL, net_updatesR, max_wave_speed)
#           endif
                
                fluxL%h = net_updatesL(1)
                fluxL%p = matmul(net_updatesL(2:3), transform_matrix)
                fluxL%max_wave_speed = max_wave_speed

                fluxR%h = net_updatesR(1)
                fluxR%p = matmul(net_updatesR(2:3), transform_matrix)
                fluxR%max_wave_speed = max_wave_speed
                      
                ! if(.not.(bL - bR)<1.0e-9)then
                        !!!!!print*
                        !!!!!print*,"hL"
                        !!!!!print*,hL
                        !!!!!print*,"hR"
                        !!!!!print*,hR
                        !!!!!print*
                        !!!!!print*,"pL"
                        !!!!!print*,pL
                        !!!!!print*,"pR"
                        !!!!!print*,pR
                        !!!!!print*
                        !!!!!print*,"bL"
                        !!!!!print*,bL
                        !!!!!print*,"bR"
                        !!!!!print*,bR
                        !!!!!print*
                        !!!!!print*,"Bathymetry diff"
                        !!!!!print*,bL- bR
                        !!!!!print*
                        !!!!!print*,"flux left h"
                        !!!!!print*,fluxL%h
                        !!!!!print*,"flux left p"
                        !!!!!print*,fluxL%p
                        !!!!!print*
                        !!!!!print*,"flux right h"
                        !!!!!print*,fluxR%h
                        !!!!!print*,"flux right p"
                        !!!!!print*,fluxR%p
                        !!!!!print*

                ! end if

		end subroutine

        pure subroutine node_write_op(local_node, neighbor_node)
            type(t_node_data), intent(inout)			    :: local_node
            type(t_node_data), intent(in)				    :: neighbor_node

            !do nothing
        end subroutine


        pure subroutine edge_write_op(local_node, neighbor_node)
            type(t_edge_data), intent(inout)			    :: local_node
            type(t_edge_data), intent(in)				    :: neighbor_node

            !do nothing
        end subroutine
	END MODULE
#endif
