! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_Adapt
		use SFC_edge_traversal
		use Conformity

		use Samoa_swe
		use Tools_noise
		use SWE_initialize
		use SWE_euler_timestep
#		if defined (_SWE_SIMD)
			use SWE_SIMD
#		endif
		implicit none

        type num_traversal_data
#			if defined (_SWE_SIMD)
				type(t_state), dimension(_SWE_SIMD_ORDER_SQUARE, 2)					:: Q_in
#			else
				type(t_state), dimension(_SWE_CELL_SIZE, 2)							:: Q_in
#			endif
        end type

		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

#		define _GT_NAME							t_swe_adaption_traversal

#		define _GT_EDGES

#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP			post_traversal_op

#		define _GT_TRANSFER_OP					transfer_op
#		define _GT_REFINE_OP					refine_op
#		define _GT_COARSEN_OP					coarsen_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op
#		define _GT_CELL_LAST_TOUCH_OP			cell_last_touch_op

#		define _GT_NODE_MPI_TYPE

#		include "SFC_generic_adaptive_traversal.f90"

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

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_swe_adaption_traversal), intent(inout)				:: traversal
			type(t_grid), intent(inout)							        :: grid

			call reduce(grid%d_max, grid%sections%elements_alloc%d_max, MPI_MAX, .true.)
		end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_swe_adaption_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section

			section%d_max = 0
		end subroutine

		subroutine post_traversal_op(traversal, section)
			type(t_swe_adaption_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section
		end subroutine
		!******************
		!Adaption operators
		!******************

		subroutine transfer_op(traversal, section, src_element, dest_element)
 			type(t_swe_adaption_traversal), intent(inout)	                            :: traversal
			type(t_grid_section), intent(inout)											:: section
			type(t_traversal_element), intent(inout)									:: src_element
			type(t_traversal_element), intent(inout)									:: dest_element
#			if defined(_SWE_SIMD)
				type(t_state), dimension(_SWE_SIMD_ORDER_SQUARE)						:: Q
				real (kind = GRID_SR), dimension(_SWE_SIMD_ORDER_SQUARE)				:: H, HU, HV, B
#			else
				type(t_state), dimension(_SWE_CELL_SIZE)								:: Q
#			endif

#			if defined (_SWE_SIMD)
				dest_element%cell%data_pers%H = src_element%cell%data_pers%H
				dest_element%cell%data_pers%HU = src_element%cell%data_pers%HU
				dest_element%cell%data_pers%HV = src_element%cell%data_pers%HV
				dest_element%cell%data_pers%B = src_element%cell%data_pers%B
#			else
				call gv_Q%read( src_element%t_element_base, Q)
				call gv_Q%write( dest_element%t_element_base, Q)
#			endif
		end subroutine

		subroutine refine_op(traversal, section, src_element, dest_element, refinement_path)
 			type(t_swe_adaption_traversal), intent(inout)	                            :: traversal
			type(t_grid_section), intent(inout)										    :: section
			type(t_traversal_element), intent(inout)									:: src_element
			type(t_traversal_element), intent(inout)									:: dest_element
			integer, dimension(:), intent(in)											:: refinement_path
#			if defined(_SWE_SIMD)
				real (kind=GRID_SR), dimension(_SWE_SIMD_ORDER_SQUARE)						:: H_in, HU_in, HV_in
				real (kind=GRID_SR), dimension(_SWE_SIMD_ORDER_SQUARE)						:: H_out, HU_out, HV_out, B_out
				integer																		:: i_plotter_type, j, row, col, cell_id
				integer, DIMENSION(_SWE_SIMD_ORDER_SQUARE,2)								:: child
				real (kind = GRID_SR), DIMENSION(2)											:: r_coords		!< cell coords within patch
#			else
				type(t_state), dimension(_SWE_CELL_SIZE)									:: Q_in
				type(t_state), dimension(_SWE_CELL_SIZE, 2)									:: Q_out
#			endif

			integer																		:: i

#			if defined (_SWE_SIMD)
				H_in = src_element%cell%data_pers%H
				HU_in = src_element%cell%data_pers%HU
				HV_in = src_element%cell%data_pers%HV
				i_plotter_type = src_element%cell%geometry%i_plotter_type
				
				do i=1, size(refinement_path)
					! compute 1st or 2nd child depending on orientation and refinement_path(i)
					! and store it in output arrays.
					! Store it in input arrays for next refinement.
					
					! decide which child is being computed (first = left to hypot, second = right to hypot.)
					if ( (refinement_path(i) == 1 .and. i_plotter_type>0) .or. (refinement_path(i) == 2 .and. i_plotter_type<0)) then
						child = SWE_SIMD_GEOMETRY%first_child
					else
						child = SWE_SIMD_GEOMETRY%second_child
					end if
					
					! the child patch values are always averages of two values from the parent patch
					do j=1, _SWE_SIMD_ORDER_SQUARE
						H_out(j) = 0.5*( H_in(child(j,1)) + H_in(child(j,2)) )
						HU_out(j) = 0.5*( HU_in(child(j,1)) + HU_in(child(j,2)) )
						HV_out(j) = 0.5*( HV_in(child(j,1)) + HV_in(child(j,2)) )
					end do
					
					! refined patches will have inverse orientation
					i_plotter_type = -1 * i_plotter_type
					
					! copy to input arrays for next iteration
					H_in = H_out
					HU_in = HU_out
					HV_in = HV_out
				end do
				
				! store results in dest_element
				dest_element%cell%data_pers%H = H_out
				dest_element%cell%data_pers%HU = HU_out
				dest_element%cell%data_pers%HV = HV_out
				
				! get bathymetry
				row = 1
				col = 1
				do i=1, _SWE_SIMD_ORDER_SQUARE
				
					! if orientation is backwards, the plotter uses a transformation that mirrors the cell...
					! this simple change solves the problem :)
					if (dest_element%cell%geometry%i_plotter_type > 0) then 
						cell_id = i
					else
						cell_id = (row-1)*(row-1) + 2 * row - col
					end if
				
					r_coords = [0_GRID_SR, 0_GRID_SR]
					do j=1,3
						r_coords(:) = r_coords(:) + SWE_SIMD_geometry%coords(:,j,cell_id) 
					end do
					r_coords = r_coords / 3
					dest_element%cell%data_pers%B(i) = get_bathymetry(section, samoa_barycentric_to_world_point(dest_element%transform_data, r_coords), section%r_time, dest_element%cell%geometry%i_depth / 2_GRID_SI)

					col = col + 1
					if (col == 2*row) then
						col = 1
						row = row + 1
					end if
				end do

				
#			else
			
			

			!state vector
			call gv_Q%read( src_element%t_element_base, Q_in)

            !convert momentum to velocity
			!Q_in(1)%p = 1.0_GRID_SR / (Q_in(1)%h - Q_in(1)%b) * Q_in(1)%p

			do i = 1, size(refinement_path)
				call t_basis_Q_split(Q_in%h, 	Q_out(:, 1)%h, 		Q_out(:, 2)%h)
				call t_basis_Q_split(Q_in%p(1),	Q_out(:, 1)%p(1),	Q_out(:, 2)%p(1))
				call t_basis_Q_split(Q_in%p(2),	Q_out(:, 1)%p(2),	Q_out(:, 2)%p(2))

				Q_in = Q_out(:, refinement_path(i))
			end do

			Q_in%b = get_bathymetry(section, samoa_barycentric_to_world_point(dest_element%transform_data, [1.0_GRID_SR / 3.0_GRID_SR, 1.0_GRID_SR / 3.0_GRID_SR]), section%r_time, dest_element%cell%geometry%i_depth / 2_GRID_SI)

            !convert velocity back to momentum
			!Q_in(1)%p = (Q_in(1)%h - Q_in(1)%b) * Q_in(1)%p

			call gv_Q%write( dest_element%t_element_base, Q_in)
#endif
		end subroutine

		subroutine coarsen_op(traversal, section, src_element, dest_element, refinement_path)
 			type(t_swe_adaption_traversal), intent(inout)	                			:: traversal
			type(t_grid_section), intent(inout)											:: section
			type(t_traversal_element), intent(inout)									:: src_element
			type(t_traversal_element), intent(inout)									:: dest_element
			integer, dimension(:), intent(in)											:: refinement_path
			integer																		:: i
#			if defined(_SWE_SIMD)
				type(t_state), dimension(_SWE_SIMD_ORDER_SQUARE)						:: Q_out
#			else
				type(t_state), dimension(_SWE_CELL_SIZE)								:: Q_out
#			endif


			!state vector

			i = refinement_path(1)
#			if defined (_SWE_SIMD)
				traversal%Q_in(:, i)%h = src_element%cell%data_pers%H
				traversal%Q_in(:, i)%p(1) = src_element%cell%data_pers%HU
				traversal%Q_in(:, i)%p(2) = src_element%cell%data_pers%HV
				traversal%Q_in(:, i)%b = src_element%cell%data_pers%B
#			else
				call gv_Q%read( src_element%t_element_base, traversal%Q_in(:, i))
#			endif

            !convert momentum to velocity
			!traversal%Q_in(1, i)%p = 1.0_GRID_SR / (traversal%Q_in(1, i)%h - traversal%Q_in(1, i)%b) * traversal%Q_in(1, i)%p

			if (i > 1) then
				call t_basis_Q_merge(traversal%Q_in(:, 1)%h,		traversal%Q_in(:, 2)%h,		Q_out%h)
				call t_basis_Q_merge(traversal%Q_in(:, 1)%p(1),	    traversal%Q_in(:, 2)%p(1),	Q_out%p(1))
				call t_basis_Q_merge(traversal%Q_in(:, 1)%p(2),	    traversal%Q_in(:, 2)%p(2),	Q_out%p(2))
				call t_basis_Q_merge(traversal%Q_in(:, 1)%b,		traversal%Q_in(:, 2)%b,		Q_out%b)

                !convert velocity back to momentum
                !Q_out(1)%p = (Q_out(1)%h - Q_out(1)%b) * Q_out(1)%p
#			if defined (_SWE_SIMD)
				dest_element%cell%data_pers%H = Q_out(:)%h
				dest_element%cell%data_pers%B = Q_out(:)%b
				dest_element%cell%data_pers%HU = Q_out(:)%p(1)
				dest_element%cell%data_pers%HV = Q_out(:)%p(2)
#			else
				call gv_Q%write( dest_element%t_element_base, Q_out)
#			endif
			end if
		end subroutine

		subroutine cell_last_touch_op(traversal, section, cell)
 			type(t_swe_adaption_traversal), intent(inout)	                :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_cell_data_ptr), intent(inout)				:: cell

			!set maximum depth
			section%d_max = max(section%d_max, cell%geometry%i_depth)
		end subroutine
	END MODULE
#endif
