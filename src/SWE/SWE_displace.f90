! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_Displace
		use SFC_edge_traversal
		use SWE_euler_timestep
		use SWE_initialize
#		if defined(_SWE_SIMD)
			use SWE_SIMD
#		endif

		use Samoa_swe

		implicit none

        type num_traversal_data
            integer (kind = GRID_DI)			:: i_refinements_issued
        end type

		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

#		define _GT_NAME							t_swe_displace_traversal

#		define _GT_EDGES

#		define _GT_ELEMENT_OP					element_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op

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

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
			type(t_swe_displace_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)					:: element

			type(t_state), dimension(_SWE_CELL_SIZE)			:: Q

			call gv_Q%read(element, Q)

			call alpha_volume_op(traversal, section, element, Q)

			call gv_Q%write(element, Q)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(traversal, section, element, Q)
			type(t_swe_displace_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)						    :: element
			type(t_state), dimension(_SWE_CELL_SIZE), intent(inout)	    :: Q
			integer (kind = GRID_SI)								    :: i
#			if defined(_SWE_SIMD)
				real (kind = GRID_SR), dimension(_SWE_SIMD_ORDER_SQUARE)        :: db
				real (kind = GRID_SR), DIMENSION(2)								:: r_coords		!< cell coords within patch
				integer (kind = GRID_SI)										:: j, row, col, cell_id
#			else
				real (kind = GRID_SR)		                                :: db
#			endif

			!evaluate initial function values at dof positions and compute DoFs
#           if defined(_ASAGI)
#				if defined(_SWE_SIMD)
					row = 1
					col = 1
					do i=1, _SWE_SIMD_ORDER_SQUARE
					
						! if orientation is backwards, the plotter uses a transformation that mirrors the cell...
						! this simple change solves the problem :)
						if (element%cell%geometry%i_plotter_type > 0) then 
							cell_id = i
						else
							cell_id = (row-1)*(row-1) + 2 * row - col
						end if
					
						r_coords = [0_GRID_SR, 0_GRID_SR]
						do j=1,3
							r_coords(:) = r_coords(:) + SWE_SIMD_geometry%coords(:,j,cell_id) 
						end do
						r_coords = r_coords / 3
						db(i) = -element%cell%data_pers%B(i) + get_bathymetry(section, samoa_barycentric_to_world_point(element%transform_data, r_coords), section%r_time, element%cell%geometry%i_depth / 2_GRID_SI)

						col = col + 1
						if (col == 2*row) then
							col = 1
							row = row + 1
						end if
					end do
					element%cell%data_pers%H = element%cell%data_pers%H + db
					element%cell%data_pers%B = element%cell%data_pers%B + db
#				else
					do i = 1, _SWE_CELL_SIZE
						db = -Q(i)%b + get_bathymetry(section, samoa_barycentric_to_world_point(element%transform_data, t_basis_Q_get_dof_coords(i)), section%r_time, element%cell%geometry%i_depth / 2_GRID_SI)
						Q(i)%h = Q(i)%h + db
						Q(i)%b = Q(i)%b + db
					end do
#				endif
#           endif

            !no coarsening while the earthquake takes place
			element%cell%geometry%refinement = max(0, element%cell%geometry%refinement)
		end subroutine
	END MODULE
#endif
