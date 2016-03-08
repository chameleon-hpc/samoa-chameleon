! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_Euler_Timestep
		use SFC_edge_traversal

		use Samoa_swe
		use c_bind_riemannsolvers
#		if defined(_SWE_SIMD)
			use SWE_SIMD
			use SWE_SIMD_Solvers
#		endif
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
		
		
#if defined (_SWE_SIMD)
		! arrays used to compute updates within each triangular patch.
		! for each edge, the left/right values are copied to these arrays 
		! before computation takes place. 
		type(t_state), dimension(_SWE_SIMD_NUM_EDGES)								:: edges_a, edges_b
		!$omp threadprivate(edges_a, edges_b)
		!..DIR$ ASSUME_ALIGNED edges_a: 64
		!..DIR$ ASSUME_ALIGNED edges_b: 64
#endif

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

#		define _GT_NODE_MPI_TYPE

#		include "SFC_generic_traversal_ringbuffer.f90"

        subroutine create_node_mpi_type(mpi_node_type)
            integer, intent(out)            :: mpi_node_type

            type(t_node_data)               :: node
            integer                         :: blocklengths(2), types(2), disps(2), i_error, extent

#           if defined(_MPI)
                blocklengths(1) = 1
                blocklengths(2) = 1

                disps(1) = 0
                disps(2) = sizeof(node)

                types(1) = MPI_LB
                types(2) = MPI_UB

                call MPI_Type_struct(2, blocklengths, disps, types, mpi_node_type, i_error); assert_eq(i_error, 0)
                call MPI_Type_commit(mpi_node_type, i_error); assert_eq(i_error, 0)

                call MPI_Type_extent(mpi_node_type, extent, i_error); assert_eq(i_error, 0)
                assert_eq(sizeof(node), extent)

                call MPI_Type_size(mpi_node_type, extent, i_error); assert_eq(i_error, 0)
                assert_eq(0, extent)
#           endif
        end subroutine

		!*******************************
		!Geometry operators
		!*******************************

		subroutine pre_traversal_grid_op(traversal, grid)
			type(t_swe_euler_timestep_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

			grid%r_dt = cfg%courant_number * cfg%scaling * get_edge_size(grid%d_max) / ((2.0_GRID_SR + sqrt(2.0_GRID_SR)) * grid%u_max)
#			if defined(_SWE_SIMD)
				grid%r_dt = grid%r_dt / _SWE_SIMD_ORDER
#			endif
		
#           if defined(_ASAGI)
                if (grid%r_time < asagi_grid_max(cfg%afh_displacement,2)) then
                    grid%r_dt = min(grid%r_dt, asagi_grid_delta(cfg%afh_displacement,2))
                end if
#           endif

			call scatter(grid%r_dt, grid%sections%elements_alloc%r_dt)

			grid%u_max = 0.0_GRID_SR
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_swe_euler_timestep_traversal), intent(inout)		:: traversal
			type(t_grid), intent(inout)							    :: grid

			real (kind = GRID_SR)       :: r_dt_cfl, r_courant_cfl

            call reduce(traversal%i_refinements_issued, traversal%children%i_refinements_issued, MPI_SUM, .true.)
            call reduce(grid%u_max, grid%sections%elements_alloc%u_max, MPI_MAX, .true.)
			grid%r_time = grid%r_time + grid%r_dt

            r_dt_cfl = cfg%scaling * get_edge_size(grid%d_max) / ((2.0_GRID_SR + sqrt(2.0_GRID_SR)) * grid%u_max)
#			if defined(_SWE_SIMD)
				r_dt_cfl = r_dt_cfl / _SWE_SIMD_ORDER
#			endif

            if (grid%r_dt > r_dt_cfl) then
                r_courant_cfl = r_dt_cfl * cfg%courant_number / grid%r_dt

                if (rank_MPI == 0) then
                    _log_write(1, '("WARNING! Time step size was too big. dt (used): ", ES10.3, ", dt (CFL): ", ES10.3, ", correct courant number: ", F0.3)') grid%r_dt, r_dt_cfl, r_courant_cfl
                end if
            end if

			call scatter(grid%r_time, grid%sections%elements_alloc%r_time)
		end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_swe_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section

			
			!this variable will be incremented for each cell with a refinement request
			traversal%i_refinements_issued = 0_GRID_DI
			section%u_max = 0.0_GRID_SR
			
		end subroutine

		function cell_to_edge_op(element, edge) result(rep)
			type(t_element_base), intent(in)						:: element
			type(t_edge_data), intent(in)						    :: edge
			type(num_cell_rep)										:: rep
			integer(kind = GRID_SI)									:: i, j, i_edge
			real(kind = GRID_SR), dimension(2, _SWE_EDGE_SIZE)		:: dof_pos
			real(kind = GRID_SR), dimension(2, 3), parameter		:: edge_offsets = reshape([0.0, 0.0, 0.0, 1.0, 1.0, 0.0], [2, 3])
			real(kind = GRID_SR), dimension(2, 3), parameter		:: edge_vectors = reshape([0.0, 1.0, 1.0, -1.0, -1.0, 0.0], [2, 3])
#			if defined(_SWE_SIMD)
				integer												:: edge_type !1=left, 2=hypotenuse, 3=right
				
#			else
				type(t_state), dimension(_SWE_CELL_SIZE)			:: Q
				call gv_Q%read(element, Q)
				_log_write(6, '(3X, A)') "swe cell to edge op:"
				_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q in: ", Q
				_log_write(6, '(4X, A, F0.3, 1X, F0.3)') "normal in : ", edge%transform_data%normal
#			endif

			
			i_edge = edge%transform_data%index
#           if (_SWE_CELL_SIZE > 1)
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
#           elif defined(_SWE_SIMD)
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
						do i=0, _SWE_SIMD_ORDER - 1
							rep%H(i+1) = H(i*i + 1)
							rep%HU(i+1) = HU(i*i + 1)
							rep%HV(i+1) = HV(i*i + 1)
							rep%B(i+1) = B(i*i + 1)
						end do
					case (2) ! hypotenuse
						do i=1, _SWE_SIMD_ORDER
							rep%H(i) = H((_SWE_SIMD_ORDER-1)*(_SWE_SIMD_ORDER-1) + 2*i - 1)
							rep%HU(i) = HU((_SWE_SIMD_ORDER-1)*(_SWE_SIMD_ORDER-1) + 2*i - 1)
							rep%HV(i) = HV((_SWE_SIMD_ORDER-1)*(_SWE_SIMD_ORDER-1) + 2*i - 1)
							rep%B(i) = B((_SWE_SIMD_ORDER-1)*(_SWE_SIMD_ORDER-1) + 2*i - 1)
						end do
					case (3) !cells with id i*i (right leg)
						do i=1, _SWE_SIMD_ORDER
							rep%H(_SWE_SIMD_ORDER + 1 - i) = H(i*i)
							rep%HU(_SWE_SIMD_ORDER + 1 - i) = HU(i*i)
							rep%HV(_SWE_SIMD_ORDER + 1 - i) = HV(i*i)
							rep%B(_SWE_SIMD_ORDER + 1 - i) = B(i*i)
						end do
					end select
				end associate
				
#			else
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
			integer 														:: i

			_log_write(6, '(3X, A)') "swe skeleton op:"

#			if defined (_SWE_LF) || defined (_SWE_LF_BATH) || defined (_SWE_LLF) || defined (_SWE_LLF_BATH)
				call compute_lf_flux(edge%transform_data%normal, rep1%Q(1), rep2%Q(1), update1%flux(1), update2%flux(1))
#			elif defined (_SWE_SIMD)
				! invert values in edges
				! cells are copied in inverse order because the neighbor 
				! ghost cells will have a mirrored numbering! See a (poorly-drawn) example:
				!         ___
				!  /\3   1\  |
				! /  \2   2\ |
				!/____\1   3\|
				
				do i=1, _SWE_SIMD_ORDER
					update1%H(i) = rep2%H(_SWE_SIMD_ORDER + 1 - i)
					update1%HU(i) = rep2%HU(_SWE_SIMD_ORDER + 1 - i)
					update1%HV(i) = rep2%HV(_SWE_SIMD_ORDER + 1 - i)
					update1%B(i) = rep2%B(_SWE_SIMD_ORDER + 1 - i)
					
					
					update2%H(i) = rep1%H(_SWE_SIMD_ORDER + 1 - i)
					update2%HU(i) = rep1%HU(_SWE_SIMD_ORDER + 1 - i)
					update2%HV(i) = rep1%HV(_SWE_SIMD_ORDER + 1 - i)
					update2%B(i) = rep1%B(_SWE_SIMD_ORDER + 1 - i)
				end do
#			else
				_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q 1 in: ", rep1%Q
				_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q 2 in: ", rep2%Q
				call compute_geoclaw_flux(edge%transform_data%normal, rep1%Q(1), rep2%Q(1), update1%flux(1), update2%flux(1))
				_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "flux 1 out: ", update1%flux
				_log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "flux 2 out: ", update2%flux
#			endif


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

			type(t_state)													:: bnd_rep
			type(t_update)													:: bnd_flux
			integer 														:: i
            !SLIP: reflect momentum at normal
			!bnd_rep = t_state(rep%Q(1)%h, rep%Q(1)%p - dot_product(rep%Q(1)%p, edge%transform_data%normal) * edge%transform_data%normal, rep%Q(1)%b)

#			ifndef _SWE_SIMD
            !NOSLIP: invert momentum (stable)
			bnd_rep = t_state(rep%Q(1)%h, -rep%Q(1)%p, rep%Q(1)%b)
#			endif			

			!OUTFLOW: copy values
			!bnd_rep = rep%Q(1)

#			if defined (_SWE_LF) || defined (_SWE_LF_BATH) || defined (_SWE_LLF) || defined (_SWE_LLF_BATH)
				call compute_lf_flux(edge%transform_data%normal, rep%Q(1), bnd_rep, update%flux(1), bnd_flux)
#			elif defined (_SWE_SIMD)
				! boundary conditions on ghost cells
				update%H = rep%H
				update%HU = - rep%HU
				update%HV = - rep%HV
				update%B = rep%B
#			else
				call compute_geoclaw_flux(edge%transform_data%normal, rep%Q(1), bnd_rep, update%flux(1), bnd_flux)
#			endif
		end subroutine

		subroutine cell_update_op(traversal, section, element, update1, update2, update3)
			type(t_swe_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)								:: section
			type(t_element_base), intent(inout)								:: element
			type(num_cell_update), intent(inout)							:: update1, update2, update3
		
#			if defined (_SWE_SIMD)
				integer 														:: i
				type(num_cell_update)											:: tmp !> ghost cells in correct order 
				real(kind = GRID_SR)									 	    :: normals(2,3)
				type(t_update)													:: update_a, update_b
				real(kind = GRID_SR)										    :: volume, edge_lengths(3), dt_div_volume, maxWaveSpeed
				real(kind = GRID_SR), DIMENSION(_SWE_SIMD_NUM_EDGES_ALIGNMENT)			:: hL, huL, hvL, bL
				real(kind = GRID_SR), DIMENSION(_SWE_SIMD_NUM_EDGES_ALIGNMENT)			:: hR, huR, hvR, bR
				!DIR$ ASSUME_ALIGNED hL: 64
				!DIR$ ASSUME_ALIGNED hR: 64
				!DIR$ ASSUME_ALIGNED huL: 64
				!DIR$ ASSUME_ALIGNED huR: 64
				!DIR$ ASSUME_ALIGNED hvL: 64
				!DIR$ ASSUME_ALIGNED hvR: 64
				!DIR$ ASSUME_ALIGNED bL: 64
				!DIR$ ASSUME_ALIGNED bR: 64
				
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
				
				associate(data => element%cell%data_pers, geom => SWE_SIMD_geometry)
				
					! copy cell values to arrays edges_a and edges_b
					! obs: cells with id > number of cells are actually ghost cells and come from edges "updates"
					! see t_SWE_SIMD_geometry for an explanation about ghost cell ordering
					do i=1, _SWE_SIMD_NUM_EDGES
						! left
						if (geom%edges_a(i) <= _SWE_SIMD_ORDER_SQUARE) then
							hL(i) = data%H(geom%edges_a(i))
							huL(i) = data%HU(geom%edges_a(i))
							hvL(i) = data%HV(geom%edges_a(i))
							bL(i) = data%B(geom%edges_a(i))
						else if (geom%edges_a(i) <= _SWE_SIMD_ORDER_SQUARE + _SWE_SIMD_ORDER) then
							hL(i) = update1%H(geom%edges_a(i) - _SWE_SIMD_ORDER_SQUARE)
							huL(i) = update1%HU(geom%edges_a(i) - _SWE_SIMD_ORDER_SQUARE)
							hvL(i) = update1%HV(geom%edges_a(i) - _SWE_SIMD_ORDER_SQUARE)
							bL(i) = update1%B(geom%edges_a(i) - _SWE_SIMD_ORDER_SQUARE)
						else if (geom%edges_a(i) <= _SWE_SIMD_ORDER_SQUARE + 2*_SWE_SIMD_ORDER) then
							hL(i) = update2%H(geom%edges_a(i) - _SWE_SIMD_ORDER_SQUARE - _SWE_SIMD_ORDER)
							huL(i) = update2%HU(geom%edges_a(i) - _SWE_SIMD_ORDER_SQUARE - _SWE_SIMD_ORDER)
							hvL(i) = update2%HV(geom%edges_a(i) - _SWE_SIMD_ORDER_SQUARE - _SWE_SIMD_ORDER)
							bL(i) = update2%B(geom%edges_a(i) - _SWE_SIMD_ORDER_SQUARE - _SWE_SIMD_ORDER)
						else 
							hL(i) = update3%H(geom%edges_a(i) - _SWE_SIMD_ORDER_SQUARE - 2*_SWE_SIMD_ORDER)
							huL(i) = update3%HU(geom%edges_a(i) - _SWE_SIMD_ORDER_SQUARE - 2*_SWE_SIMD_ORDER)
							hvL(i) = update3%HV(geom%edges_a(i) - _SWE_SIMD_ORDER_SQUARE - 2*_SWE_SIMD_ORDER)
							bL(i) = update3%B(geom%edges_a(i) - _SWE_SIMD_ORDER_SQUARE - 2*_SWE_SIMD_ORDER)
						end if
						
						! right
						if (geom%edges_b(i) <= _SWE_SIMD_ORDER_SQUARE) then
							hR(i) = data%H(geom%edges_b(i))
							huR(i) = data%HU(geom%edges_b(i))
							hvR(i) = data%HV(geom%edges_b(i))
							bR(i) = data%B(geom%edges_b(i))
						else if (geom%edges_b(i) <= _SWE_SIMD_ORDER_SQUARE + _SWE_SIMD_ORDER) then
							hR(i) = update1%H(geom%edges_b(i) - _SWE_SIMD_ORDER_SQUARE)
							huR(i) = update1%HU(geom%edges_b(i) - _SWE_SIMD_ORDER_SQUARE)
							hvR(i) = update1%HV(geom%edges_b(i) - _SWE_SIMD_ORDER_SQUARE)
							bR(i) = update1%B(geom%edges_b(i) - _SWE_SIMD_ORDER_SQUARE)
						else if (geom%edges_b(i) <= _SWE_SIMD_ORDER_SQUARE + 2*_SWE_SIMD_ORDER) then
							hR(i) = update2%H(geom%edges_b(i) - _SWE_SIMD_ORDER_SQUARE - _SWE_SIMD_ORDER)
							huR(i) = update2%HU(geom%edges_b(i) - _SWE_SIMD_ORDER_SQUARE - _SWE_SIMD_ORDER)
							hvR(i) = update2%HV(geom%edges_b(i) - _SWE_SIMD_ORDER_SQUARE - _SWE_SIMD_ORDER)
							bR(i) = update2%B(geom%edges_b(i) - _SWE_SIMD_ORDER_SQUARE - _SWE_SIMD_ORDER)
						else 
							hR(i) = update3%H(geom%edges_b(i) - _SWE_SIMD_ORDER_SQUARE - 2*_SWE_SIMD_ORDER)
							huR(i) = update3%HU(geom%edges_b(i) - _SWE_SIMD_ORDER_SQUARE - 2*_SWE_SIMD_ORDER)
							hvR(i) = update3%HV(geom%edges_b(i) - _SWE_SIMD_ORDER_SQUARE - 2*_SWE_SIMD_ORDER)
							bR(i) = update3%B(geom%edges_b(i) - _SWE_SIMD_ORDER_SQUARE - 2*_SWE_SIMD_ORDER)
						end if
					end do
					
					! compute net_updates
#					if defined (_USE_SIMD)
#						if defined(_SWE_FWAVE)
							call compute_updates_fwave_simd(normals, hL, huL, hvL, bL, hR, huR, hvR, bR, maxWaveSpeed)
#						else
							call compute_updates_hlle_simd(normals, hL, huL, hvL, bL, hR, huR, hvR, bR, maxWaveSpeed)
#						endif
						section%u_max = max(section%u_max, maxWaveSpeed)
#					else					
					! NO-SIMD version
 					do i=1, _SWE_SIMD_NUM_EDGES
						edges_a(i)%h = hL(i)
						edges_a(i)%p(1) = huL(i)
						edges_a(i)%p(2) = hvL(i)
						edges_a(i)%b = bL(i)
						edges_b(i)%h = hR(i)
						edges_b(i)%p(1) = huR(i)
						edges_b(i)%p(2) = hvR(i)
						edges_b(i)%b = bR(i)
 						call compute_geoclaw_flux(normals(:,geom%edges_orientation(i)), edges_a(i), edges_b(i), update_a, update_b)
 						hL(i) = update_a%h
 						huL(i) = update_a%p(1)
 						hvL(i) = update_a%p(2)
 						hR(i) = update_b%h
 						huR(i) = update_b%p(1)
 						hvR(i) = update_b%p(2)
 						section%u_max = max(section%u_max, update_a%max_wave_speed)
 					end do
#					endif
 					
					! if land is flooded, init water height to dry tolerance and
					! velocity to zero
					do i=1,_SWE_SIMD_NUM_EDGES
						if (geom%edges_a(i) <= _SWE_SIMD_ORDER_SQUARE .and. hL(i) > 0 .and. data%H(geom%edges_a(i)) < data%B(geom%edges_a(i)) + cfg%dry_tolerance) then
							data%H(geom%edges_a(i)) = data%B(geom%edges_a(i)) + cfg%dry_tolerance
							data%HU(geom%edges_a(i)) = 0.0_GRID_SR
							data%HV(geom%edges_a(i)) = 0.0_GRID_SR
						endif
					
						if (geom%edges_b(i) <= _SWE_SIMD_ORDER_SQUARE .and. hR(i) > 0 .and. data%H(geom%edges_b(i)) < data%B(geom%edges_b(i)) + cfg%dry_tolerance) then
							data%H(geom%edges_b(i)) = data%B(geom%edges_a(i)) + cfg%dry_tolerance
							data%HU(geom%edges_b(i)) = 0.0_GRID_SR
							data%HV(geom%edges_b(i)) = 0.0_GRID_SR
						endif
					end do

					! update unknowns
					volume = cfg%scaling * cfg%scaling * element%cell%geometry%get_volume() / (_SWE_SIMD_ORDER_SQUARE)
					edge_lengths = cfg%scaling * element%cell%geometry%get_edge_sizes() / _SWE_SIMD_ORDER
					dt_div_volume = section%r_dt / volume
					do i=1, _SWE_SIMD_NUM_EDGES
						if (geom%edges_a(i) <= _SWE_SIMD_ORDER_SQUARE) then !ignore ghost cells
							data%H(geom%edges_a(i)) = data%H(geom%edges_a(i)) - hL(i) * edge_lengths(geom%edges_orientation(i)) * dt_div_volume
							data%HU(geom%edges_a(i)) = data%HU(geom%edges_a(i)) - huL(i) * edge_lengths(geom%edges_orientation(i)) * dt_div_volume
							data%HV(geom%edges_a(i)) = data%HV(geom%edges_a(i)) - hvL(i) * edge_lengths(geom%edges_orientation(i)) * dt_div_volume
						end if
						if (geom%edges_b(i) <= _SWE_SIMD_ORDER_SQUARE) then
							data%H(geom%edges_b(i)) = data%H(geom%edges_b(i)) - hR(i) * edge_lengths(geom%edges_orientation(i)) * dt_div_volume
							data%HU(geom%edges_b(i)) = data%HU(geom%edges_b(i)) - huR(i) * edge_lengths(geom%edges_orientation(i)) * dt_div_volume
							data%HV(geom%edges_b(i)) = data%HV(geom%edges_b(i)) - hvR(i) * edge_lengths(geom%edges_orientation(i)) * dt_div_volume
						end if
					end do
					
					! if the water level falls below the dry tolerance, set water surface to 0 and velocity to 0
					where (data%H < data%B + cfg%dry_tolerance) 
						data%H = min(data%B, 0.0_GRID_SR)
						data%HU = 0.0_GRID_SR
						data%HV = 0.0_GRID_SR
					end where

				end associate
#			else
				!local variables

				type(t_state)   :: dQ(_SWE_CELL_SIZE)

				call volume_op(element%cell%geometry, traversal%i_refinements_issued, element%cell%geometry%i_depth, &
					element%cell%geometry%refinement, section%u_max, dQ, [update1%flux, update2%flux, update3%flux], section%r_dt)

				!if land is flooded, init water height to dry tolerance and velocity to 0
				if (element%cell%data_pers%Q(1)%h < element%cell%data_pers%Q(1)%b + cfg%dry_tolerance .and. dQ(1)%h > 0.0_GRID_SR) then
					element%cell%data_pers%Q(1)%h = element%cell%data_pers%Q(1)%b + cfg%dry_tolerance
					element%cell%data_pers%Q(1)%p = [0.0_GRID_SR, 0.0_GRID_SR]
					!print '("Wetting:", 2(X, F0.0))', cfg%scaling * element%transform_data%custom_data%offset + cfg%offset
				end if

				call gv_Q%add(element, dQ)

				!if the water level falls below the dry tolerance, set water surface to 0 and velocity to 0
				if (element%cell%data_pers%Q(1)%h < element%cell%data_pers%Q(1)%b + cfg%dry_tolerance) then
					element%cell%data_pers%Q(1)%h = min(element%cell%data_pers%Q(1)%b, 0.0_GRID_SR)
					element%cell%data_pers%Q(1)%p = [0.0_GRID_SR, 0.0_GRID_SR]
				end if
#			endif
		end subroutine

		subroutine cell_last_touch_op(traversal, section, cell)
			type(t_swe_euler_timestep_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_cell_data_ptr), intent(inout)				:: cell
			real(kind = GRID_SR)							    :: b_norm

#			if defined (_SWE_SIMD)
				b_norm = minval(abs(cell%data_pers%H - cell%data_pers%B))
#			else
				b_norm = minval(abs(cell%data_pers%Q%h - cell%data_pers%Q%b))
#			endif

			!refine also on the coasts
			if (cell%geometry%i_depth < cfg%i_max_depth .and. b_norm < 100.0_GRID_SR) then
				cell%geometry%refinement = 1
				traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
			else if (b_norm < 300.0_GRID_SR) then
				cell%geometry%refinement = max(cell%geometry%refinement, 0)
			endif
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine volume_op(cell, i_refinements_issued, i_depth, i_refinement, u_max, dQ, fluxes, r_dt)
			type(fine_triangle), intent(in)				                        :: cell
			integer (kind = GRID_DI), intent(inout)							    :: i_refinements_issued
			integer (kind = BYTE), intent(in)							        :: i_depth
			integer (kind = BYTE), intent(out)							        :: i_refinement
			real(kind = GRID_SR), intent(inout)								    :: u_max
			type(t_state), dimension(:), intent(out)						    :: dQ
			type(t_update), dimension(:), intent(in)						    :: fluxes
			real(kind = GRID_SR), intent(in)								    :: r_dt

			real(kind = GRID_SR)											    :: volume, dQ_norm, edge_lengths(3)
			integer (kind = BYTE)												:: i

			_log_write(6, '(3X, A)') "swe cell update op:"
			_log_write(6, '(4X, A, 4(X, F0.3))') "edge 1 flux in:", fluxes(1)
			_log_write(6, '(4X, A, 4(X, F0.3))') "edge 2 flux in:", fluxes(2)
			_log_write(6, '(4X, A, 4(X, F0.3))') "edge 3 flux in:", fluxes(3)

			volume = cfg%scaling * cfg%scaling * cell%get_volume()
			edge_lengths = cfg%scaling * cell%get_edge_sizes()

			dQ%h = sum(edge_lengths * fluxes%h)
			dQ%p(1) = sum(edge_lengths * fluxes%p(1))
			dQ%p(2) = sum(edge_lengths * fluxes%p(2))
			dQ%b = 0.0_GRID_SR

			!set refinement condition

			i_refinement = 0
			dQ_norm = dot_product(dQ(1)%p, dQ(1)%p)

			if (i_depth < cfg%i_max_depth .and. dQ_norm > (cfg%scaling * 2.0_GRID_SR) ** 2) then
				i_refinement = 1
				i_refinements_issued = i_refinements_issued + 1_GRID_DI
			else if (i_depth > cfg%i_min_depth .and. dQ_norm < (cfg%scaling * 1.0_GRID_SR) ** 2) then
				i_refinement = -1
			endif

			u_max = max(u_max, maxval(fluxes%max_wave_speed))

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

			!real(kind = GRID_SR), parameter						:: dry_tol = 0.01_GRID_SR       --- replaced by flag parameter cfg%dry_tolerance
			real(kind = GRID_SR)								:: vL, vR, alpha

#           if defined(_SWE_LF_BATH) || defined(_SWE_LLF_BATH)
                if (QL%h - QL%b < cfg%dry_tolerance) then
                    vL = 0.0_GRID_SR
                    fluxL%max_wave_speed = 0.0_GRID_SR
                else
                    vL = DOT_PRODUCT(normal, QL%p / (QL%h - QL%b))
                    fluxL%max_wave_speed = sqrt(g * (QL%h - QL%b)) + abs(vL)
                end if

                if (QR%h - QR%b < cfg%dry_tolerance) then
                    vR = 0.0_GRID_SR
                    fluxR%max_wave_speed = 0.0_GRID_SR
                else
                    vR = DOT_PRODUCT(normal, QR%p / (QR%h - QR%b))
                    fluxR%max_wave_speed = sqrt(g * (QR%h - QR%b)) + abs(vR)
                end if

#               if defined(_SWE_LLF_BATH)
                    alpha = max(fluxL%max_wave_speed, fluxR%max_wave_speed)
#               else
                    alpha = 100.0_GRID_SR
#               endif

                fluxL%h = 0.5_GRID_SR * (vL * (QL%h - QL%b) + vR * (QR%h - QR%b) + alpha * (QL%h - QR%h))
                fluxR%h = -fluxL%h

                fluxL%p = 0.5_GRID_SR * ((vL + alpha) * QL%p + (vR - alpha) * QR%p) + 0.5_GRID_SR * g * (max(QR%h - QL%b, 0.0_GRID_SR) ** 2) * normal
                fluxR%p = -0.5_GRID_SR * ((vL + alpha) * QL%p + (vR - alpha) * QR%p) - 0.5_GRID_SR * g * (max(QL%h - QR%b, 0.0_GRID_SR) ** 2) * normal
#           else
                real(kind = GRID_SR), parameter					:: b = -1000.0_GRID_SR     !default constant bathymetry

                !use the height of the water pillars for computation
                vL = DOT_PRODUCT(normal, QL%p / (QL%h - b))
                vR = DOT_PRODUCT(normal, QR%p / (QR%h - b))

                fluxL%max_wave_speed = sqrt(g * (QL%h - b)) + abs(vL)
                fluxR%max_wave_speed = sqrt(g * (QR%h - b)) + abs(vR)

#               if defined(_SWE_LLF)
                    alpha = max(fluxL%max_wave_speed, fluxR%max_wave_speed)
#               else
                    alpha = 100.0_GRID_SR
#               endif

                fluxL%h = 0.5_GRID_SR * (vL * (QL%h - b) + vR * (QR%h - b) + alpha * (QL%h - QR%h))
                fluxR%h = -fluxL%h

                fluxL%p = 0.5_GRID_SR * ((vL + alpha) * QL%p + (vR - alpha) * QR%p) + 0.5_GRID_SR * g * (QR%h - b) * (QR%h - b) * normal
                fluxR%p = -0.5_GRID_SR * ((vL + alpha) * QL%p + (vR - alpha) * QR%p) - 0.5_GRID_SR * g * (QL%h - b) * (QL%h - b) * normal
#           endif
		end subroutine

		!> Augmented Riemann solver from geoclaw
		subroutine compute_geoclaw_flux(normal, QL, QR, fluxL, fluxR)
			type(t_state), intent(in)           :: QL, QR
			type(t_update), intent(out)         :: fluxL, fluxR
			real(kind = GRID_SR), intent(in)    :: normal(2)

			real(kind = GRID_SR)				:: transform_matrix(2, 2)
			real(kind = GRID_SR)			    :: net_updatesL(3), net_updatesR(3), max_wave_speed
			real(kind = GRID_SR)                :: pL(2), pR(2), hL, hR, bL, bR

			transform_matrix(1, :) = normal
			transform_matrix(2, :) = [-normal(2), normal(1)]

			pL = matmul(transform_matrix, QL%p)
			pR = matmul(transform_matrix, QR%p)
			hL = QL%h - QL%b
			hR = QR%h - QR%b
			bL = QL%b
			bR = QR%b

#           if defined(_SWE_FWAVE)
                call c_bind_geoclaw_solver(GEOCLAW_FWAVE, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_SWE_SSQ_FWAVE)
                call c_bind_geoclaw_solver(GEOCLAW_SSQ_FWAVE, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_SWE_AUG_RIEMANN)
                call c_bind_geoclaw_solver(GEOCLAW_AUG_RIEMANN, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
#           endif

			fluxL%h = net_updatesL(1)
			fluxL%p = matmul(net_updatesL(2:3), transform_matrix)
			fluxL%max_wave_speed = max_wave_speed

			fluxR%h = net_updatesR(1)
			fluxR%p = matmul(net_updatesR(2:3), transform_matrix)
			fluxR%max_wave_speed = max_wave_speed
		end subroutine

	END MODULE
#endif
