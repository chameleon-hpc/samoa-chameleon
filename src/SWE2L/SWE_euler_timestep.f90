! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE2L)
	MODULE SWE2L_Euler_Timestep
		use SFC_edge_traversal

		use Samoa_swe2l
		use c_bind_riemannsolvers
#       if defined(_SWE_PATCH)
            use SWE_PATCH
#       endif
#       if defined(_HLLE_FLUX)
            use SWE_HLLE
#       endif
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

			if (cfg%r_max_time > 0.0_SR) then
                grid%r_dt = min(cfg%r_max_time, grid%r_dt)
            end if

			if (cfg%r_output_time_step > 0.0_SR) then
                grid%r_dt = min(cfg%r_output_time_step, grid%r_dt)
            end if

#           if defined(_ASAGI)
                !if we are in the earthquake phase, limit the simulation time step by the earthquake time step
                if (grid%r_time < cfg%t_max_eq) then
                    grid%r_dt = min(grid%r_dt, cfg%dt_eq)
                end if
#           endif

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
                
                associate(data => element%cell%data_pers)
                    ! copy boundary values to respective edges
                    ! left leg cells go to edge 1
                    ! hypotenuse cells go to edge 2
                    ! right leg cells go to edge 3
                    select case (edge_type)
                    case (1) !cells with id i*i+1 (left leg)
                        do i=0, _SWE_PATCH_ORDER - 1
                            rep%H(i+1) = data%H(i*i + 1)
                            rep%HU(i+1) = data%HU(i*i + 1)
                            rep%HV(i+1) = data%HV(i*i + 1)
                            rep%B(i+1) = data%B(i*i + 1)
                            rep%H2(i+1) = data%H2(i*i + 1)
                            rep%HU2(i+1) = data%HU2(i*i + 1)
                            rep%HV2(i+1) = data%HV2(i*i + 1)
                        end do
                    case (2) ! hypotenuse
                        do i=1, _SWE_PATCH_ORDER
                            rep%H(i) = data%H((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
                            rep%HU(i) = data%HU((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
                            rep%HV(i) = data%HV((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
                            rep%B(i) = data%B((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
                            rep%H2(i) = data%H2((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
                            rep%HU2(i) = data%HU2((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
                            rep%HV2(i) = data%HV2((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
                        end do
                    case (3) !cells with id i*i (right leg)
                        do i=1, _SWE_PATCH_ORDER
                            rep%H(_SWE_PATCH_ORDER + 1 - i) = data%H(i*i)
                            rep%HU(_SWE_PATCH_ORDER + 1 - i) = data%HU(i*i)
                            rep%HV(_SWE_PATCH_ORDER + 1 - i) = data%HV(i*i)
                            rep%B(_SWE_PATCH_ORDER + 1 - i) = data%B(i*i)
                            rep%H2(_SWE_PATCH_ORDER + 1 - i) = data%H2(i*i)
                            rep%HU2(_SWE_PATCH_ORDER + 1 - i) = data%HU2(i*i)
                            rep%HV2(_SWE_PATCH_ORDER + 1 - i) = data%HV2(i*i)
                        end do
                    end select
                end associate
                
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
			_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q 1 in: ", rep1%Q
			_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q 2 in: ", rep2%Q

#			if defined (_LF_FLUX) || defined (_LF_BATH_FLUX) || defined (_LLF_FLUX) || defined (_LLF_BATH_FLUX)
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
                    update1%H2(i) = rep2%H2(_SWE_PATCH_ORDER + 1 - i)
                    update1%HU2(i) = rep2%HU2(_SWE_PATCH_ORDER + 1 - i)
                    update1%HV2(i) = rep2%HV2(_SWE_PATCH_ORDER + 1 - i)
                    
                    update2%H(i) = rep1%H(_SWE_PATCH_ORDER + 1 - i)
                    update2%HU(i) = rep1%HU(_SWE_PATCH_ORDER + 1 - i)
                    update2%HV(i) = rep1%HV(_SWE_PATCH_ORDER + 1 - i)
                    update2%B(i) = rep1%B(_SWE_PATCH_ORDER + 1 - i)
                    update2%H2(i) = rep1%H2(_SWE_PATCH_ORDER + 1 - i)
                    update2%HU2(i) = rep1%HU2(_SWE_PATCH_ORDER + 1 - i)
                    update2%HV2(i) = rep1%HV2(_SWE_PATCH_ORDER + 1 - i)
                end do
#           else
                call compute_geoclaw_flux_wrapper(edge, rep1, rep2, update1, update2)
#	endif

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

			type(num_cell_rep)												:: bnd_rep
			type(num_cell_update)											:: bnd_flux
			integer                                                         :: i

            !SLIP: reflect momentum at normal
			!bnd_rep = t_state(rep%Q(1)%h, rep%Q(1)%p - dot_product(rep%Q(1)%p, edge%transform_data%normal) * edge%transform_data%normal, rep%Q(1)%b)

            !NOSLIP: invert momentum (stable)
			!bnd_rep = t_state(rep%Q(1)%h, -rep%Q(1)%p, rep%Q(1)%b)

			!OUTFLOW: copy values
			bnd_rep = rep

#			if defined (_LF_FLUX) || defined (_LF_BATH_FLUX) || defined (_LLF_FLUX) || defined (_LLF_BATH_FLUX)
				call compute_lf_flux(edge%transform_data%normal, rep%Q(1), bnd_rep, update%flux(1), bnd_flux)
#           elif defined (_SWE_PATCH)
                ! boundary conditions on ghost cells
                update%H = rep%H
                update%HU = rep%HU
                update%HV = rep%HV
                update%B = rep%B
                update%H2 = rep%H2
                update%HU2 = rep%HU2
                update%HV2 = rep%HV2
#           else
                call compute_geoclaw_flux_wrapper(edge, rep, bnd_rep, update, bnd_flux)
#			endif
		end subroutine
		
		subroutine compute_geoclaw_flux_wrapper(edge, rep1, rep2, update1, update2)
            type(t_edge_data), intent(in)                                   :: edge
            type(num_cell_rep), intent(in)                                  :: rep1, rep2
            type(num_cell_update), intent(out)                              :: update1, update2
            
            real(kind = GRID_SR) :: normal_x, normal_y
            real(kind = GRID_SR) :: hL, hR, huL, huR, hvL, hvR, bL, bR
            real(kind = GRID_SR) :: hL2, hR2, huL2, huR2, hvL2, hvR2
            real(kind = GRID_SR) :: upd_hL, upd_hR, upd_huL, upd_huR, upd_hvL, upd_hvR
            real(kind = GRID_SR) :: upd_hL2, upd_hR2, upd_huL2, upd_huR2, upd_hvL2, upd_hvR2
            real(kind = GRID_SR) :: maxWaveSpeed
            
                
            normal_x = edge%transform_data%normal(1)
            normal_y = edge%transform_data%normal(2)
            
            hL   = rep1%Q(1)%h
            huL  = rep1%Q(1)%p(1)
            hvL  = rep1%Q(1)%p(2)
            bL   = rep1%Q(1)%b
            hL2  = rep1%Q(1)%h2
            huL2 = rep1%Q(1)%p2(1)
            hvL2 = rep1%Q(1)%p2(2)
            
            hR   = rep2%Q(1)%h
            huR  = rep2%Q(1)%p(1)
            hvR  = rep2%Q(1)%p(2)
            bR   = rep2%Q(1)%b
            hR2  = rep2%Q(1)%h2
            huR2 = rep2%Q(1)%p2(1)
            hvR2 = rep2%Q(1)%p2(2)         

            call compute_geoclaw_flux_in_patch(normal_x, normal_y, &
                                                hL, hR, huL, huR, hvL, hvR, bL, bR, &
                                                hL2, hR2, huL2, huR2, hvL2, hvR2, &
                                                upd_hL, upd_hR, upd_huL, upd_huR, upd_hvL, upd_hvR, &
                                                upd_hL2, upd_hR2, upd_huL2, upd_huR2, upd_hvL2, upd_hvR2, &
                                                maxWaveSpeed)
            
            update1%flux(1)%h     = upd_hL
            update1%flux(1)%p(1)  = upd_huL
            update1%flux(1)%p(2)  = upd_hvL
            update1%flux(1)%h2    = upd_hL2
            update1%flux(1)%p2(1) = upd_huL2
            update1%flux(1)%p2(2)  = upd_hvL2
            update1%flux(1)%max_wave_speed = maxWaveSpeed
            
            update2%flux(1)%h     = upd_hR
            update2%flux(1)%p(1)  = upd_huR
            update2%flux(1)%p(2)  = upd_hvR
            update2%flux(1)%h2    = upd_hR2
            update2%flux(1)%p2(1) = upd_huR2
            update2%flux(1)%p2(2)  = upd_hvR2
            update2%flux(1)%max_wave_speed = maxWaveSpeed            
		end subroutine 

		subroutine cell_update_op(traversal, section, element, update1, update2, update3)
			type(t_swe_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)						:: element
			type(num_cell_update), intent(inout)						:: update1, update2, update3
#           if defined (_SWE_PATCH)
                integer                                                         :: i
                type(num_cell_update)                                           :: tmp !> ghost cells in correct order 
                real(kind = GRID_SR)                                            :: volume, edge_lengths(3), maxWaveSpeed, maxWaveSpeedLocal, dQ_max_norm, dt_div_volume
                real(kind = GRID_SR), DIMENSION(_SWE_PATCH_ORDER_SQUARE)        :: dQ_H, dQ_HU, dQ_HV !> deltaQ, used to compute cell updates
                real(kind = GRID_SR), DIMENSION(_SWE_PATCH_ORDER_SQUARE)        :: dQ_H2, dQ_HU2, dQ_HV2
                real(kind = GRID_SR), DIMENSION(_SWE_PATCH_NUM_EDGES)           :: hL, huL, hvL, bL
                real(kind = GRID_SR), DIMENSION(_SWE_PATCH_NUM_EDGES)           :: hR, huR, hvR, bR
                real(kind = GRID_SR), DIMENSION(_SWE_PATCH_NUM_EDGES)           :: hL2, huL2, hvL2
                real(kind = GRID_SR), DIMENSION(_SWE_PATCH_NUM_EDGES)           :: hR2, huR2, hvR2
                real(kind = GRID_SR), DIMENSION(_SWE_PATCH_NUM_EDGES)           :: upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR
                real(kind = GRID_SR), DIMENSION(_SWE_PATCH_NUM_EDGES)           :: upd_hL2, upd_huL2, upd_hvL2, upd_hR2, upd_huR2, upd_hvR2
                real(kind = GRID_SR), dimension(2,3)                            :: normals ! 1st index = x or y, 2nd index = edge orientation in patch (1, 2 or 3, see SWE_PATCH)
                real(kind = GRID_SR)                                            :: normal_x, normal_y
                
                !DIR$ ASSUME_ALIGNED hL: 64
                !DIR$ ASSUME_ALIGNED hR: 64
                !DIR$ ASSUME_ALIGNED huL: 64
                !DIR$ ASSUME_ALIGNED huR: 64
                !DIR$ ASSUME_ALIGNED hvL: 64
                !DIR$ ASSUME_ALIGNED hvR: 64
                !DIR$ ASSUME_ALIGNED bL: 64
                !DIR$ ASSUME_ALIGNED bR: 64
                !DIR$ ASSUME_ALIGNED hL2: 64
                !DIR$ ASSUME_ALIGNED hR2: 64
                !DIR$ ASSUME_ALIGNED huL2: 64
                !DIR$ ASSUME_ALIGNED huR2: 64
                !DIR$ ASSUME_ALIGNED hvL2: 64
                !DIR$ ASSUME_ALIGNED hvR2: 64
                !DIR$ ASSUME_ALIGNED upd_hL2: 64
                !DIR$ ASSUME_ALIGNED upd_hR2: 64
                !DIR$ ASSUME_ALIGNED upd_huL2: 64
                !DIR$ ASSUME_ALIGNED upd_huR2: 64
                !DIR$ ASSUME_ALIGNED upd_hvL2: 64
                !DIR$ ASSUME_ALIGNED upd_hvR2: 64
                
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

                
                if (element%cell%geometry%i_plotter_type > 0) then ! if orientation = forward, reverse updates
                    tmp=update1
                    update1=update3
                    update3=tmp
                end if
                
                ! init some variables
                dQ_H = 0.0_GRID_SR
                dQ_HU = 0.0_GRID_SR
                dQ_HV = 0.0_GRID_SR
                dQ_H2 = 0.0_GRID_SR
                dQ_HU2 = 0.0_GRID_SR
                dQ_HV2 = 0.0_GRID_SR
                maxWaveSpeed = 0.0_GRID_SR
                
                volume = cfg%scaling * cfg%scaling * element%cell%geometry%get_volume() / (_SWE_PATCH_ORDER_SQUARE)
                dt_div_volume = section%r_dt / volume
                edge_lengths = cfg%scaling * element%cell%geometry%get_edge_sizes() / _SWE_PATCH_ORDER
                
                associate(data => element%cell%data_pers, geom => SWE_PATCH_geometry)
                
                    ! copy cell values to arrays edges_a and edges_b
                    ! obs: cells with id > number of cells are actually ghost cells and come from edges "updates"
                    ! see t_SWE_PATCH_geometry for an explanation about ghost cell ordering
                    do i=1, _SWE_PATCH_NUM_EDGES
                            
                        ! cells left to the edges
                        if (geom%edges_a(i) <= _SWE_PATCH_ORDER_SQUARE) then
                            ! data for internal cells come from cell data
                            hL(i) = data%H(geom%edges_a(i))
                            huL(i) = data%HU(geom%edges_a(i))
                            hvL(i) = data%HV(geom%edges_a(i))
                            bL(i) = data%B(geom%edges_a(i))
                            hL2(i) = data%H2(geom%edges_a(i))
                            huL2(i) = data%HU2(geom%edges_a(i))
                            hvL2(i) = data%HV2(geom%edges_a(i))
                        else if (geom%edges_a(i) <= _SWE_PATCH_ORDER_SQUARE + _SWE_PATCH_ORDER) then
                            ! data for ghost cells come from edge (update) date
                            hL(i) = update1%H(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE)
                            huL(i) = update1%HU(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE)
                            hvL(i) = update1%HV(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE)
                            bL(i) = update1%B(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE)
                            hL2(i) = update1%H2(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE)
                            huL2(i) = update1%HU2(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE)
                            hvL2(i) = update1%HV2(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE)
                        else if (geom%edges_a(i) <= _SWE_PATCH_ORDER_SQUARE + 2*_SWE_PATCH_ORDER) then
                            ! data for ghost cells come from edge (update) date
                            hL(i) = update2%H(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                            huL(i) = update2%HU(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                            hvL(i) = update2%HV(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                            bL(i) = update2%B(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                            hL2(i) = update2%H2(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                            huL2(i) = update2%HU2(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                            hvL2(i) = update2%HV2(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                        else 
                            ! data for ghost cells come from edge (update) date
                            hL(i) = update3%H(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                            huL(i) = update3%HU(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                            hvL(i) = update3%HV(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                            bL(i) = update3%B(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                            hL2(i) = update3%H2(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                            huL2(i) = update3%HU2(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                            hvL2(i) = update3%HV2(geom%edges_a(i) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                        end if
                        
                        ! cells right to the edges
                        if (geom%edges_b(i) <= _SWE_PATCH_ORDER_SQUARE) then
                            ! data for internal cells come from cell data
                            hR(i) = data%H(geom%edges_b(i))
                            huR(i) = data%HU(geom%edges_b(i))
                            hvR(i) = data%HV(geom%edges_b(i))
                            bR(i) = data%B(geom%edges_b(i))
                            hR2(i) = data%H2(geom%edges_b(i))
                            huR2(i) = data%HU2(geom%edges_b(i))
                            hvR2(i) = data%HV2(geom%edges_b(i))
                        else if (geom%edges_b(i) <= _SWE_PATCH_ORDER_SQUARE + _SWE_PATCH_ORDER) then
                            ! data for ghost cells come from edge (update) date
                            hR(i) = update1%H(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE)
                            huR(i) = update1%HU(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE)
                            hvR(i) = update1%HV(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE)
                            bR(i) = update1%B(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE)
                            hR2(i) = update1%H2(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE)
                            huR2(i) = update1%HU2(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE)
                            hvR2(i) = update1%HV2(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE)
                        else if (geom%edges_b(i) <= _SWE_PATCH_ORDER_SQUARE + 2*_SWE_PATCH_ORDER) then
                            ! data for ghost cells come from edge (update) date
                            hR(i) = update2%H(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                            huR(i) = update2%HU(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                            hvR(i) = update2%HV(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                            bR(i) = update2%B(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                            hR2(i) = update2%H2(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                            huR2(i) = update2%HU2(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                            hvR2(i) = update2%HV2(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                        else 
                            ! data for ghost cells come from edge (update) date
                            hR(i) = update3%H(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                            huR(i) = update3%HU(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                            hvR(i) = update3%HV(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                            bR(i) = update3%B(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                            hR2(i) = update3%H2(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                            huR2(i) = update3%HU2(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                            hvR2(i) = update3%HV2(geom%edges_b(i) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                        end if
                    end do

#                   if defined(_SWE_PATCH_VEC_SIMD) || defined(_SWE_PATCH_VEC_INLINE)
                        ! Vectorization! (Requires OpenMP 4.0 or later)
                        !$OMP SIMD PRIVATE(maxWaveSpeedLocal,normal_x,normal_y) REDUCTION(max: maxWaveSpeed)
#                   endif
                    do i=1, _SWE_PATCH_NUM_EDGES
                    
                        normal_x = normals(1,geom%edges_orientation(i))
                        normal_y = normals(2,geom%edges_orientation(i))

#                       if defined(_SWE_PATCH_VEC_INLINE)
                            ! Warning: inlining this subroutine into an OMP SIMD loop may cause
                            ! bugs and incorrect calculations depending on the compiler version
                            ! Check the results!
                            
                            ! Recommended compiler: ifort 17.0
                            
                            !DIR$ FORCEINLINE
#                       endif
                        call compute_geoclaw_flux_in_patch(normal_x, normal_y, &
                                                            hL(i), hR(i), huL(i), huR(i), hvL(i), hvR(i), bL(i), bR(i), &
                                                            hL2(i), hR2(i), huL2(i), huR2(i), hvL2(i), hvR2(i), &
                                                            upd_hL(i), upd_hR(i), upd_huL(i), upd_huR(i), upd_hvL(i), upd_hvR(i), &
                                                            upd_hL2(i), upd_hR2(i), upd_huL2(i), upd_huR2(i), upd_hvL2(i), upd_hvR2(i), &
                                                            maxWaveSpeedLocal)
                        maxWaveSpeed = max(maxWaveSpeed, maxWaveSpeedLocal)
                    end do
                    
                    ! compute dQ
                    do i=1, _SWE_PATCH_NUM_EDGES
                        if (geom%edges_a(i) <= _SWE_PATCH_ORDER_SQUARE) then !ignore ghost cells
                            dQ_H(geom%edges_a(i)) = dQ_H(geom%edges_a(i)) + upd_hL(i) * edge_lengths(geom%edges_orientation(i))
                            dQ_HU(geom%edges_a(i)) = dQ_HU(geom%edges_a(i)) + upd_huL(i) * edge_lengths(geom%edges_orientation(i))
                            dQ_HV(geom%edges_a(i)) = dQ_HV(geom%edges_a(i)) + upd_hvL(i) * edge_lengths(geom%edges_orientation(i))
                            dQ_H2(geom%edges_a(i)) = dQ_H2(geom%edges_a(i)) + upd_hL2(i) * edge_lengths(geom%edges_orientation(i))
                            dQ_HU2(geom%edges_a(i)) = dQ_HU2(geom%edges_a(i)) + upd_huL2(i) * edge_lengths(geom%edges_orientation(i))
                            dQ_HV2(geom%edges_a(i)) = dQ_HV2(geom%edges_a(i)) + upd_hvL2(i) * edge_lengths(geom%edges_orientation(i))
                        end if
                        if (geom%edges_b(i) <= _SWE_PATCH_ORDER_SQUARE) then
                            dQ_H(geom%edges_b(i)) = dQ_H(geom%edges_b(i)) + upd_hR(i) * edge_lengths(geom%edges_orientation(i))
                            dQ_HU(geom%edges_b(i)) = dQ_HU(geom%edges_b(i)) + upd_huR(i) * edge_lengths(geom%edges_orientation(i))
                            dQ_HV(geom%edges_b(i)) = dQ_HV(geom%edges_b(i)) + upd_hvR(i) * edge_lengths(geom%edges_orientation(i))
                            dQ_H2(geom%edges_b(i)) = dQ_H2(geom%edges_b(i)) + upd_hR2(i) * edge_lengths(geom%edges_orientation(i))
                            dQ_HU2(geom%edges_b(i)) = dQ_HU2(geom%edges_b(i)) + upd_huR2(i) * edge_lengths(geom%edges_orientation(i))
                            dQ_HV2(geom%edges_b(i)) = dQ_HV2(geom%edges_b(i)) + upd_hvR2(i) * edge_lengths(geom%edges_orientation(i))
                        end if
                    end do
                    
                    !set refinement condition -> Here I am using the same ones as in the original no-patches implementation, but considering only the max value.
                    element%cell%geometry%refinement = 0
                    dQ_max_norm = maxval(max(abs(dQ_H),abs(dQ_H2)))

                    if (element%cell%geometry%i_depth < cfg%i_max_depth .and. dQ_max_norm > 5.0_GRID_SR * cfg%scaling * get_edge_size(cfg%i_max_depth) / _SWE_PATCH_ORDER ) then
                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
                    else if (element%cell%geometry%i_depth > cfg%i_min_depth .and. dQ_max_norm < 5.0_GRID_SR * cfg%scaling * get_edge_size(cfg%i_max_depth) / (_SWE_PATCH_ORDER * 8.0_SR) ) then
                        element%cell%geometry%refinement = -1
                    endif

                    dQ_H = dQ_H * (-dt_div_volume)
                    dQ_HU = dQ_HU * (-dt_div_volume)
                    dQ_HV = dQ_HV * (-dt_div_volume)
                    dQ_H2 = dQ_H2 * (-dt_div_volume)
                    dQ_HU2 = dQ_HU2 * (-dt_div_volume)
                    dQ_HV2 = dQ_HV2 * (-dt_div_volume)
                    
                    ! if land is flooded, init water height to dry tolerance and
                    ! velocity to zero
!                     where (data%H < data%H2 + cfg%dry_tolerance .and. dQ_H > 0.0_GRID_SR)
!                         data%H = data%H2 + cfg%dry_tolerance
!                         data%HU = 0.0_GRID_SR
!                         data%HV = 0.0_GRID_SR
!                     end where
!                     where (data%H2 < data%B + cfg%dry_tolerance .and. dQ_H2 > 0.0_GRID_SR)
!                         data%H2 = data%B + cfg%dry_tolerance
!                         data%HU2 = 0.0_GRID_SR
!                         data%HV2 = 0.0_GRID_SR
!                     end where
! 

                    ! update unknowns
                    data%H = data%H + dQ_H + dQ_H2
                    data%HU = data%HU + dQ_HU
                    data%HV = data%HV + dQ_HV
                    data%H2 = data%H2 + dQ_H2
                    data%HU2 = data%HU2 + dQ_HU2
                    data%HV2 = data%HV2 + dQ_HV2
                    
                    !if the water level falls below the dry tolerance, set water level to 0 and velocity to 0
                    where (data%H2 < data%B + cfg%dry_tolerance) 
                        data%H2 = data%B
                        data%HU2 = 0.0_GRID_SR
                        data%HV2 = 0.0_GRID_SR
                    end where
                    where (data%H < data%H2 + cfg%dry_tolerance) 
                        data%H = data%H2
                        data%HU = 0.0_GRID_SR
                        data%HV = 0.0_GRID_SR
                    end where
                    
                    ! compute next dt
                    section%r_dt_new = min(section%r_dt_new, volume / (edge_lengths(2) * maxWaveSpeed) )

                end associate
#           else
			!local variables

			type(t_state)   :: dQ(_SWE_CELL_SIZE)

			call volume_op(element%cell%geometry, traversal%i_refinements_issued, element%cell%geometry%i_depth, &
                element%cell%geometry%refinement, section%r_dt_new, dQ, [update1%flux, update2%flux, update3%flux], section%r_dt)

! 			!if land is flooded, init water height to dry tolerance and velocity to 0
! 			if (element%cell%data_pers%Q(1)%h < element%cell%data_pers%Q(1)%b + cfg%dry_tolerance .and. dQ(1)%h > 0.0_GRID_SR) then
!                 element%cell%data_pers%Q(1)%h = element%cell%data_pers%Q(1)%b + cfg%dry_tolerance
!                 element%cell%data_pers%Q(1)%p = [0.0_GRID_SR, 0.0_GRID_SR]
!                 !print '("Wetting:", 2(X, F0.0))', cfg%scaling * element%transform_data%custom_data%offset + cfg%offset
!             end if

            call gv_Q%add(element, dQ)

			!if the water level falls below the dry tolerance, set water level to 0 and velocity to 0
           if (element%cell%data_pers%Q(1)%h < element%cell%data_pers%Q(1)%h2 + cfg%dry_tolerance) then
                element%cell%data_pers%Q(1)%h = element%cell%data_pers%Q(1)%h2
                element%cell%data_pers%Q(1)%p = [0.0_GRID_SR, 0.0_GRID_SR]
           end if
           if (element%cell%data_pers%Q(1)%h2 < element%cell%data_pers%Q(1)%b + cfg%dry_tolerance) then
                element%cell%data_pers%Q(1)%h2 = element%cell%data_pers%Q(1)%b
                element%cell%data_pers%Q(1)%p2 = [0.0_GRID_SR, 0.0_GRID_SR]
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
			dQ%h2= dot_product(edge_lengths, fluxes%h2)
			dQ%p(1) = dot_product(edge_lengths, fluxes%p(1))
			dQ%p(2) = dot_product(edge_lengths, fluxes%p(2))
			dQ%p2(1) = dot_product(edge_lengths, fluxes%p2(1))
			dQ%p2(2) = dot_product(edge_lengths, fluxes%p2(2))
			dQ%b = 0.0_GRID_SR
			
			!set refinement condition

			i_refinement = 0
			dQ_norm = max(abs(dQ(1)%h),abs(dQ(1)%h2))

			if (i_depth < cfg%i_max_depth .and. dQ_norm > refinement_threshold * cfg%scaling * get_edge_size(cfg%i_max_depth)) then
				i_refinement = 1
				i_refinements_issued = i_refinements_issued + 1_GRID_DI
			else if (i_depth > cfg%i_min_depth .and. dQ_norm < refinement_threshold * cfg%scaling * get_edge_size(cfg%i_max_depth) / 8.0_SR) then
				i_refinement = -1
			endif

            !This will cause a division by zero if the wave speeds are 0.
            !Bue to the min operator, the error will not affect the time step.
            !r_dt_new = min(r_dt_new, volume / dot_product(edge_lengths, fluxes%max_wave_speed))
#           if defined(_DEBUG)
                if ( maxval(fluxes%max_wave_speed) .ne. 0 ) then
                       r_dt_new = min(r_dt_new, volume / (edge_lengths(2) * maxval(fluxes%max_wave_speed)))
                endif
#           else
            r_dt_new = min(r_dt_new, volume / (edge_lengths(2) * maxval(fluxes%max_wave_speed)))
#           endif   

            do i = 1, _SWE_CELL_SIZE
                dQ(i)%t_dof_state = dQ(i)%t_dof_state * (-r_dt / volume)
            end do

			_log_write(6, '(4X, A, 4(X, F0.3))') "dQ out: ", dQ
		end subroutine

		!> Lax Friedrichs flux. Depending on compiler flags, the function implements
		!> the global or local variant with or without bathymetry
		pure subroutine compute_lf_flux(normal, QL, QR, fluxL, fluxR)
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

		!> Riemann solvers from geoclaw
		subroutine compute_geoclaw_flux(normal, QL, QR, fluxL, fluxR)
			type(t_state), intent(in)           :: QL, QR
			type(t_update), intent(out)         :: fluxL, fluxR
			real(kind = GRID_SR), intent(in)    :: normal(2)

			real(kind = GRID_SR)				:: transform_matrix(2, 2)
			real(kind = GRID_SR)			    :: net_updatesL(3), net_updatesR(3), max_wave_speed
			real(kind = GRID_SR)                :: pL(2), pR(2), pL2(2), pR2(2), hL, hR, hL2, hR2, bL, bR

! 			transform_matrix(1, :) = normal
! 			transform_matrix(2, :) = [-normal(2), normal(1)]
! 
! 			pL = matmul(transform_matrix, QL%p)
! 			pR = matmul(transform_matrix, QR%p)
! 			pL2 = matmul(transform_matrix, QL%p2)
! 			pR2 = matmul(transform_matrix, QR%p2)
! 			hL2 = QL%h2 - QL%b
! 			hR2 = QR%h2 - QR%b
! 			hL = QL%h - hL2
! 			hR = QR%h - hR2
! 			bL = QL%b
! 			bR = QR%b
! 			
! 			call solve_riemann_problem_SWE2L(6,6,2,3, h_l,h_r, hu_l,hu_r, hv_l,hv_r, b_l,b_r, h_hat_l,h_hat_r, &
!                                                     fWaves, waveSpeeds, [real(cfg%dry_tolerance,kind=8),real(cfg%dry_tolerance,kind=8)], g) 

			
! #           if defined(_FWAVE_FLUX)
!                 call c_bind_geoclaw_solver(GEOCLAW_FWAVE, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
! #           elif defined(_AUG_RIEMANN_FLUX)
!                 call c_bind_geoclaw_solver(GEOCLAW_AUG_RIEMANN, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
! #           endif
! 
! 			fluxL%h = net_updatesL(1)
! 			fluxL%p = matmul(net_updatesL(2:3), transform_matrix)
! 			fluxL%max_wave_speed = max_wave_speed
! 
! 			fluxR%h = net_updatesR(1)
! 			fluxR%p = matmul(net_updatesR(2:3), transform_matrix)
! 			fluxR%max_wave_speed = max_wave_speed
	end subroutine

!#       if defined(_SWE_PATCH)
            ! version for SWE patches (vectorizable)
            subroutine compute_geoclaw_flux_in_patch(normal_x, normal_y, &
                                                        hL, hR, huL, huR, hvL, hvR, bL, bR, &
                                                        hL2, hR2, huL2, huR2, hvL2, hvR2, &
                                                        upd_hL, upd_hR, upd_huL, upd_huR, upd_hvL, upd_hvR, &
                                                        upd_hL2, upd_hR2, upd_huL2, upd_huR2, upd_hvL2, upd_hvR2, &
                                                        maxWaveSpeed)
#               if defined(_SWE_PATCH_VEC_SIMD)
                    !$OMP DECLARE SIMD(compute_geoclaw_flux_in_patch)
#               endif
                real(kind = GRID_SR), intent(in)    :: normal_x, normal_y
                real(kind = GRID_SR), intent(inout)    :: hL, hR, huL, huR, hvL, hvR, bL, bR
                real(kind = GRID_SR), intent(inout)    :: hL2, hR2, huL2, huR2, hvL2, hvR2
                real(kind = GRID_SR), intent(out)   :: upd_hL, upd_hR, upd_huL, upd_huR, upd_hvL, upd_hvR
                real(kind = GRID_SR), intent(out)   :: upd_hL2, upd_hR2, upd_huL2, upd_huR2, upd_hvL2, upd_hvR2
                real(kind = GRID_SR), intent(out)   :: maxWaveSpeed

                real(kind = GRID_SR)                :: transform_matrix(2, 2)
                real(kind = GRID_SR)                :: pL, pR ! pressure forcing, not considered here
                real(kind = GRID_SR)                :: waveSpeeds(6) ! output of Riemann solver: sw
                real(kind = GRID_SR)                :: fWaves(6,6) ! output of Riemann solver: fw
                integer                             :: equationNumber, waveNumber, i, j
                
                ! Used to call the Riemann solver:
                real(kind=8), dimension(2) :: h_l, h_r, hu_l, hu_r, hv_l, hv_r
                real(kind=8):: b_l, b_r
                real(kind=8), dimension(2) :: h_hat_l, h_hat_r

                transform_matrix(1, :) = [ normal_x, normal_y ]
                transform_matrix(2, :) = [-normal_y, normal_x]

                ! TODO: I think these can now be replace with simple matrix operations
                call apply_transformations_before(transform_matrix, huL, hvL)
                call apply_transformations_before(transform_matrix, huR, hvR)
                call apply_transformations_before(transform_matrix, huL2, hvL2)
                call apply_transformations_before(transform_matrix, huR2, hvR2)
                
                ! in this implementation we consider h1 and h2 the water elevation relative to some mean level,
                ! but the Riemann solvers considers that these variables represent the actual height of the water column:
                hL = hL - hL2
                hL2 = hL2 - bL
                hR = hR - hR2
                hR2 = hR2 - bR
                
                ! pressure forcing is not considered in these solvers
                pL = 0.0_GRID_SR
                pR = 0.0_GRID_SR
                
                ! init output from solve_riemann_problem
                waveSpeeds(:) = 0.0_GRID_SR
                fWaves(:,:) = 0.0_GRID_SR
                
                ! copy data to input vars for Riemann solver
                h_l  = [hL,hL2]
                hu_l = [huL,huL2]
                hv_l = [hvL,hvL2]
                h_r  = [hR,hR2]
                hu_r = [huR,huR2]
                hv_r = [hvR,hvR2]
                b_l = bL
                b_r = bR
                ! These are not considered in this sam(oa)^2 implementation:
                h_hat_l = 0.0_GRID_SR
                h_hat_r = 0.0_GRID_SR

#               if defined(_SWE_PATCH_VEC_INLINE)
                    !DIR$ FORCEINLINE
#               endif
                call solve_riemann_problem_SWE2L(6,6,2,3, h_l,h_r, hu_l,hu_r, hv_l,hv_r, b_l,b_r, h_hat_l,h_hat_r, &
                                                    fWaves, waveSpeeds, [real(cfg%dry_tolerance,kind=8),real(cfg%dry_tolerance,kind=8)], g) 

                upd_hL = 0.0_GRID_SR
                upd_huL = 0.0_GRID_SR
                upd_hvL = 0.0_GRID_SR
                upd_hR = 0.0_GRID_SR
                upd_huR = 0.0_GRID_SR
                upd_hvR = 0.0_GRID_SR
                upd_hL2 = 0.0_GRID_SR
                upd_huL2 = 0.0_GRID_SR
                upd_hvL2 = 0.0_GRID_SR
                upd_hR2 = 0.0_GRID_SR
                upd_huR2 = 0.0_GRID_SR
                upd_hvR2 = 0.0_GRID_SR
                ! use Riemann solution to compute net updates
                do  waveNumber=1,6
                    if (waveSpeeds(waveNumber) <= 0.d0) then
                        upd_hL = upd_hL + fWaves(1,waveNumber)
                        upd_huL = upd_huL + fWaves(2,waveNumber)
                        upd_hvL = upd_hvL + fWaves(3,waveNumber)
                        upd_hL2 = upd_hL2 + fWaves(4,waveNumber)
                        upd_huL2 = upd_huL2 + fWaves(5,waveNumber)
                        upd_hvL2 = upd_hvL2 + fWaves(6,waveNumber)
                    else !if (waveSpeeds(waveNumber) > 0.d0) then
                        upd_hR = upd_hR + fWaves(1,waveNumber)
                        upd_huR = upd_huR + fWaves(2,waveNumber)
                        upd_hvR = upd_hvR + fWaves(3,waveNumber)                    
                        upd_hR2 = upd_hR2 + fWaves(4,waveNumber)
                        upd_huR2 = upd_huR2 + fWaves(5,waveNumber)
                        upd_hvR2 = upd_hvR2 + fWaves(6,waveNumber)  
                    endif
                enddo

                !compute maximum wave speed
                waveSpeeds = abs(waveSpeeds)
                maxWaveSpeed = maxVal(waveSpeeds)
                
                ! inverse transformations
                call apply_transformations_after(transform_matrix, upd_huL, upd_hvL)
                call apply_transformations_after(transform_matrix, upd_huR, upd_hvR)
                call apply_transformations_after(transform_matrix, upd_huL2, upd_hvL2)
                call apply_transformations_after(transform_matrix, upd_huR2, upd_hvR2)

            end subroutine
            
            ! change base so hu/hv become ortogonal/perperdicular to edge
            subroutine apply_transformations_before(transform_matrix, hu, hv)
#               if defined(_SWE_PATCH_VEC_SIMD)
                    !$OMP DECLARE SIMD(apply_transformations_before)
#               endif
                real(kind = GRID_SR), intent(in)         :: transform_matrix(2,2)
                real(kind = GRID_SR), intent(inout)      :: hu, hv           
                
                real(kind = GRID_SR)                     :: temp
                
                temp = hu
                hu = transform_matrix(1,1) * hu + transform_matrix(1,2) * hv
                hv = transform_matrix(2,1) * temp + transform_matrix(2,2) * hv
            end subroutine
            
            ! transform back to original base
            subroutine apply_transformations_after(transform_matrix, hu, hv)
#               if defined(_SWE_PATCH_VEC_SIMD)
                    !$OMP DECLARE SIMD(apply_transformations_after)
#               endif
                real(kind = GRID_SR), intent(in)         :: transform_matrix(2,2)
                real(kind = GRID_SR), intent(inout)      :: hu, hv           
                
                real(kind = GRID_SR)                     :: temp
                
                temp = hu
                hu = transform_matrix(1,1) * hu + transform_matrix(2,1) * hv
                hv = transform_matrix(1,2) * temp + transform_matrix(2,2) * hv
            end subroutine
!#       endif

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
