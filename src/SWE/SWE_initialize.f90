! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_Initialize
 		use Tools_noise

		use SFC_edge_traversal
		use SWE_euler_timestep

		use Samoa_swe
#		if defined(_SWE_PATCH)
			use SWE_PATCH
#		endif
		implicit none

        type num_traversal_data
            integer (kind = GRID_DI)			:: i_refinements_issued
        end type

		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

		PUBLIC get_bathymetry

#		define _GT_NAME							t_swe_init_traversal

#		define _GT_EDGES
#		define _GT_EDGES_TEMP

#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_ELEMENT_OP					element_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
			type(t_swe_init_traversal), intent(inout)		        :: traversal
			type(t_grid), intent(inout)							    :: grid

			grid%r_time = 0.0_GRID_SR
			grid%r_dt = 0.0_GRID_SR
			grid%d_max = cfg%i_max_depth
			grid%u_max = sqrt(g)

            call scatter(grid%r_time, grid%sections%elements_alloc%r_time)
		end subroutine

		subroutine post_traversal_grid_op(traversal, grid)
			type(t_swe_init_traversal), intent(inout)		        :: traversal
			type(t_grid), intent(inout)							    :: grid

            call reduce(traversal%i_refinements_issued, traversal%children%i_refinements_issued, MPI_SUM, .true.)
            call reduce(grid%u_max, grid%sections%elements_alloc%u_max, MPI_MAX, .true.)

            grid%r_time = grid%r_time + grid%r_dt
		end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_swe_init_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section

			!this variable will be incremented for each cell with a refinement request
			traversal%i_refinements_issued = 0
			section%u_max = 0.0_GRID_SR
		end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
			type(t_swe_init_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout)						:: element
			
#			if defined(_SWE_PATCH)
				type(t_state), dimension(_SWE_PATCH_ORDER_SQUARE)	:: Q
#			else
				type(t_state), dimension(_SWE_CELL_SIZE)			:: Q
#			endif

			call alpha_volume_op(traversal, section, element, Q)

# 			if defined(_SWE_PATCH)
				element%cell%data_pers%H = Q(:)%h
				element%cell%data_pers%HU = Q(:)%p(1)
				element%cell%data_pers%HV = Q(:)%p(2)
				element%cell%data_pers%B = Q(:)%b
#			else
				call gv_Q%write(element, Q)
#			endif

		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(traversal, section, element, Q)
			type(t_swe_init_traversal), intent(inout)				    :: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)							:: element
#			if defined(_SWE_PATCH)
				type(t_state), dimension(_SWE_PATCH_ORDER_SQUARE), intent(out)	:: Q
				real (kind = GRID_SR), dimension(_SWE_PATCH_ORDER_SQUARE)		:: lambda
				real (kind = GRID_SR), DIMENSION(2)								:: r_coords		!< cell coords within patch
				integer (kind = GRID_SI)										:: j, row, col, cell_id
#			else
				type(t_state), dimension(_SWE_CELL_SIZE), intent(out)	:: Q
				real (kind = GRID_SR), dimension(_SWE_CELL_SIZE)		:: lambda
#			endif

			real (kind = GRID_SR), dimension(2)							:: pos
			integer (kind = GRID_SI)									:: i
			real (kind = GRID_SR), parameter		                	:: r_test_points(2, 3) = reshape([1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [2, 3])
			real (kind = GRID_SR)                                   	:: centroid_square(2), centroid_triangle(2)
			type(t_state), dimension(3)									:: Q_test
			

			!evaluate initial function values at dof positions and compute DoFs

			do i = 1, _SWE_CELL_SIZE
				Q(i) = get_initial_state(section, samoa_barycentric_to_world_point(element%transform_data, t_basis_Q_get_dof_coords(i)), element%cell%geometry%i_depth / 2_GRID_SI)
			end do
		
#			if defined(_SWE_PATCH)
				row = 1
				col = 1
				do i=1, _SWE_PATCH_ORDER_SQUARE
				
					! if orientation is backwards, the plotter uses a transformation that mirrors the cell...
					! this simple change solves the problem :)
					if (element%cell%geometry%i_plotter_type > 0) then 
						cell_id = i
					else
						cell_id = (row-1)*(row-1) + 2 * row - col
					end if
				
					r_coords = [0_GRID_SR, 0_GRID_SR]
					do j=1,3
						r_coords(:) = r_coords(:) + SWE_PATCH_geometry%coords(:,j,cell_id) 
					end do
					r_coords = r_coords / 3
					Q(i) = get_initial_state(section, samoa_barycentric_to_world_point(element%transform_data, r_coords), element%cell%geometry%i_depth / 2_GRID_SI)

					col = col + 1
					if (col == 2*row) then
						col = 1
						row = row + 1
					end if
				end do
#			endif
			element%cell%geometry%refinement = 0

			if (element%cell%geometry%i_depth < cfg%i_min_depth) then
                !refine if the minimum depth is not met

 				element%cell%geometry%refinement = 1
				traversal%i_refinements_issued = traversal%i_refinements_issued + 1
            else if (element%cell%geometry%i_depth < cfg%i_max_depth) then
                do i = 1, 3
                    Q_test(i) = get_initial_state(section, samoa_barycentric_to_world_point(element%transform_data, r_test_points(:, i)), element%cell%geometry%i_depth / 2_GRID_SI)
                end do

#               if defined (_ASAGI)
                    centroid_square = 0.5_GRID_SR * [asagi_grid_min(cfg%afh_displacement,0) + asagi_grid_max(cfg%afh_displacement,0), asagi_grid_min(cfg%afh_displacement,1) + asagi_grid_max(cfg%afh_displacement,1)]
                    centroid_square = 1.0_GRID_SR / cfg%scaling * (centroid_square - cfg%offset)
                    centroid_square = samoa_world_to_barycentric_point(element%transform_data, centroid_square)

                    centroid_triangle = [1.0_GRID_SR/3.0_GRID_SR, 1.0_GRID_SR/3.0_GRID_SR]
                    centroid_triangle = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, centroid_triangle) + cfg%offset
                    centroid_triangle = [ &
                        (centroid_triangle(1) - asagi_grid_min(cfg%afh_displacement,0)) / (asagi_grid_max(cfg%afh_displacement,0) - asagi_grid_min(cfg%afh_displacement,0)), &
                        (centroid_triangle(2) - asagi_grid_min(cfg%afh_displacement,1)) / (asagi_grid_max(cfg%afh_displacement,1) - asagi_grid_min(cfg%afh_displacement,1)) &
                    ]

                    if (maxval(Q_test%h - Q_test%b) > 0.0 .and. minval(Q_test%h - Q_test%b) <= 0.0) then
                        !refine coast lines

                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1
                    elseif (centroid_square(1) >= 0.0 .and. centroid_square(2) >= 0.0 .and. centroid_square(1) + centroid_square(2) <= 1.0) then
                        !refine the triangle if it contains the centroid of the initial condition

                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1
                    elseif (centroid_triangle(1) >= 0.0 .and. centroid_triangle(2) >= 0.0 .and. centroid_triangle(1) <= 1.0 .and. centroid_triangle(2) <= 1.0) then
                        !refine the triangle if its centroid is contained in the initial condition

                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1
                    end if
#               else
                    if (maxval(Q_test%h - Q_test%b) > 0.0 .and. minval(Q_test%h - Q_test%b) <= 0.0 .or. any(Q_test%h .ne. 0.0)) then
                        !refine coast lines and initial displacement

                        element%cell%geometry%refinement = 1
                        traversal%i_refinements_issued = traversal%i_refinements_issued + 1
                    end if
#               endif
			end if

			!estimate initial u_max

			where (Q%h - Q%b > 0)
				lambda = sqrt(g * (Q%h - Q%b)) + sqrt((Q%p(1) * Q%p(1) + Q%p(2) * Q%p(2)) / ((Q%h - Q%b) * (Q%h - Q%b)))
			elsewhere
				lambda = 0.0_GRID_SR
			end where

			section%u_max = max(section%u_max, maxval(lambda))
		end subroutine

		function get_initial_state(section, x, lod) result(Q)
			type(t_grid_section), intent(inout)					:: section
			real (kind = GRID_SR), dimension(:), intent(in)		:: x
			integer (kind = GRID_SI), intent(in)				:: lod						!< level of detail
			type(t_state)										:: Q

			!init height, momentum and bathymetry
			Q%t_dof_state = get_initial_dof_state(section, x, lod)
			Q%b = get_bathymetry(section, x, 0.0_GRID_SR, lod)
		end function

		function get_initial_dof_state(section, x, lod) result(Q)
			type(t_grid_section), intent(inout)							:: section
			real (kind = GRID_SR), dimension(:), intent(in)		:: x						!< position in world coordinates
			integer (kind = GRID_SI), intent(in)				:: lod						!< level of detail
			type(t_dof_state)									:: Q						!< initial state

			real (kind = GRID_SR), dimension(2), parameter		:: dam_center = [0.5, 0.5]
			real (kind = GRID_SR), parameter					:: dam_radius = 0.2
			real (kind = GRID_SR), parameter					:: outer_height = 0.0
			real (kind = GRID_SR), parameter					:: inner_height = 10.0
            real (kind = GRID_SR)                               :: xs(2)

            xs = cfg%scaling * x + cfg%offset

#			if defined(_ASAGI)
				Q%h = 0.0_GRID_SR
#			else
				Q%h = 0.5_GRID_SR * (inner_height + outer_height) + (inner_height - outer_height) * sign(0.5_GRID_SR, (dam_radius ** 2) - dot_product(xs - dam_center, xs - dam_center))
#			endif

			Q%p = 0.0_GRID_SR
		end function

		function get_bathymetry(section, x, t, lod) result(bathymetry)
			type(t_grid_section), intent(inout)					:: section
			real (kind = GRID_SR), dimension(:), intent(in)		:: x						!< position in world coordinates
			real (kind = GRID_SR), intent(in)		            :: t						!< simulation time
			integer (kind = GRID_SI), intent(in)				:: lod						!< level of detail
            real (kind = GRID_SR)								:: bathymetry				!< bathymetry

            double precision			                               :: xs(3), ts



#			if defined(_ASAGI)
                xs(1:2) = cfg%scaling * x + cfg%offset

#               if defined(_ASAGI_TIMING)
                    section%stats%r_asagi_time = section%stats%r_asagi_time - get_wtime()
#               endif

		if (asagi_grid_min(cfg%afh_bathymetry,0) <= xs(1) .and. asagi_grid_min(cfg%afh_bathymetry,1) <= xs(2) &
                        .and. xs(1) <= asagi_grid_max(cfg%afh_bathymetry,0) .and. xs(2) <= asagi_grid_max(cfg%afh_bathymetry,1)) then

                    bathymetry = asagi_get_float(cfg%afh_bathymetry, xs, 0)
                else
                    bathymetry = -5000.0 !we assume that the sea floor is constant here
                end if

                if (asagi_grid_min(cfg%afh_displacement,0) <= xs(1) .and. asagi_grid_min(cfg%afh_displacement,1) <= xs(2) &
                        .and. xs(1) <= asagi_grid_max(cfg%afh_displacement,0) .and. xs(2) <= asagi_grid_max(cfg%afh_displacement,1) &
                        .and. asagi_grid_min(cfg%afh_displacement,2) < t) then

                    ts = min(t, asagi_grid_max(cfg%afh_displacement,2))
                    
                    xs(3) = ts

                    bathymetry = bathymetry + asagi_get_float(cfg%afh_displacement, xs, 0)
                end if

#               if defined(_ASAGI_TIMING)
                    section%stats%r_asagi_time = section%stats%r_asagi_time + get_wtime()
#               endif
#			else
                real (kind = GRID_SR), dimension(2), parameter		:: dam_center = [0.5, 0.5]
                real (kind = GRID_SR), parameter					:: dam_radius = 0.1
                real (kind = GRID_SR), parameter					:: outer_height = -100.0
                real (kind = GRID_SR), parameter					:: inner_height = -5.0

                xs(1:2) = cfg%scaling * x + cfg%offset
				bathymetry = 0.5_GRID_SR * (inner_height + outer_height) + (inner_height - outer_height) * sign(0.5_GRID_SR, (dam_radius ** 2) - dot_product(xs(1:2) - dam_center, xs(1:2) - dam_center))
#			endif
		end function
	END MODULE
#endif

