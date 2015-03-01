! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_initialize_pressure
		use SFC_edge_traversal

        use iso_c_binding
		use Samoa_darcy
		use Tools_noise

		implicit none

        type num_traversal_data
        end type

		type(darcy_gv_p)						:: gv_p
		type(darcy_gv_saturation)				:: gv_saturation

		public get_base_permeability, get_porosity

#		define	_GT_NAME						t_darcy_init_pressure_traversal

#		define _GT_NODES

#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op
#		define _GT_INNER_NODE_FIRST_TOUCH_OP	inner_node_first_touch_op

#		define	_GT_ELEMENT_OP					element_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine pre_traversal_grid_op(traversal, grid)
 			type(t_darcy_init_pressure_traversal), intent(inout)      	:: traversal
 			type(t_grid), intent(inout)							        :: grid

			call scatter(grid%r_time, grid%sections%elements_alloc%r_time)
		end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
 			type(t_darcy_init_pressure_traversal)                       :: traversal
 			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)											:: element

			call alpha_volume_op(section, element, element%cell%data_pers%base_permeability, element%cell%data_pers%porosity)
		end subroutine

		subroutine node_first_touch_op(traversal, section, nodes)
 			type(t_darcy_init_pressure_traversal), intent(in)   :: traversal
 			type(t_grid_section), intent(inout)				    :: section
			type(t_node_data), intent(inout)			        :: nodes(:)

            integer :: i

            do i = 1, size(nodes)
                call inner_node_first_touch_op(traversal, section, nodes(i))
            end do
		end subroutine

		subroutine inner_node_first_touch_op(traversal, section, node)
 			type(t_darcy_init_pressure_traversal), intent(in)                       :: traversal
 			type(t_grid_section), intent(inout)		    :: section
			type(t_node_data), intent(inout)			:: node
			integer										:: i

			call pressure_pre_dof_op(real(cfg%r_p_prod, GRID_SR), node%data_pers%p, node%data_pers%r, node%data_pers%d, node%data_pers%A_d)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(section, element, base_permeability, porosity)
 			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(in)									:: element

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR), intent(out)									:: base_permeability(:, :)
                real (kind = GRID_SR), intent(out)									:: porosity(:)

                real (kind = GRID_SR)       :: x(3)
                integer (kind = GRID_SI)    :: i

                x(1:2) = samoa_barycentric_to_world_point(element%transform_data, [1.0_SR / 3.0_SR, 1.0_SR / 3.0_SR])

                do i = 1, _DARCY_LAYERS
                    x(3) = (i - 0.5_SR) / real(_DARCY_LAYERS, SR)

                    !horizontal permeability is assumed to be isotropic
                    base_permeability(i, :) = get_base_permeability(section, x, element%cell%geometry%i_depth / 2_GRID_SI)
                    porosity(i) = get_porosity(section, x)
                end do
#           else
                real (kind = GRID_SR), intent(out)									:: base_permeability
                real (kind = GRID_SR), intent(out)									:: porosity

                real (kind = GRID_SR)   :: x(3)

                x(1:2) = samoa_barycentric_to_world_point(element%transform_data, [1.0_SR / 3.0_SR, 1.0_SR / 3.0_SR])
                x(3) = 1.0e-5_SR

                base_permeability = get_base_permeability(section, x, element%cell%geometry%i_depth / 2_GRID_SI)
                porosity = get_porosity(section, x)
#           endif
 		end subroutine

		function get_base_permeability(section, x, lod) result(r_base_permeability)
 			type(t_grid_section), intent(inout)					:: section
			real (kind = GRID_SR), intent(in)		            :: x(:)						!< position in world coordinates
			integer (kind = GRID_SI), intent(in)				:: lod						!< level of detail

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR)							:: r_base_permeability(2)	!< permeability tensor
#           else
                real (kind = GRID_SR)							:: r_base_permeability		!< permeability tensor
#           endif

            real (kind = GRID_SR)                               :: xs(3)
            real (kind = GRID_SR)						        :: buffer(3)

            xs(1:2) = cfg%scaling * x(1:2) + cfg%offset
            xs(3) = cfg%scaling * max(1, _DARCY_LAYERS) * cfg%dz * x(3)

            assert_ge(x(1), 0.0); assert_ge(x(2), 0.0)
            assert_le(x(1), 1.0); assert_le(x(2), 1.0)

#           if defined(_ASAGI_TIMING)
                section%stats%r_asagi_time = section%stats%r_asagi_time - get_wtime()
#           endif

#           if defined(_ASAGI)
                if (grid_min_x(cfg%afh_permeability_X) <= xs(1) .and. grid_min_y(cfg%afh_permeability_X) <= xs(2) &
                        .and. xs(1) <= grid_max_x(cfg%afh_permeability_X) .and. xs(2) <= grid_max_y(cfg%afh_permeability_X)) then

#                   if (_DARCY_LAYERS > 0)
                        buffer(1) = grid_get_double_3d(cfg%afh_permeability_X, dble(xs(1)), dble(xs(2)), dble(xs(3)), 0)
                        !buffer(2) = grid_get_double_3d(cfg%afh_permeability_Y, dble(xs(1)), dble(xs(2)), dble(xs(3)), 0)
                        buffer(3) = grid_get_double_3d(cfg%afh_permeability_Z, dble(xs(1)), dble(xs(2)), dble(xs(3)), 0)

                        !assume horizontally isotropic permeability
                        !assert(abs(buffer(1) - buffer(2)) < epsilon(1.0_SR))

                        r_base_permeability(1) = buffer(1)
                        r_base_permeability(2) = buffer(3)
#                   else
                        buffer(1) = grid_get_double_3d(cfg%afh_permeability_X, dble(xs(1)), dble(xs(2)), dble(xs(3)), 0)
                        !buffer(2) = grid_get_double_3d(cfg%afh_permeability_Y, dble(xs(1)), dble(xs(2)), dble(xs(3)), 0)
                        !buffer(3) = grid_get_double_3d(cfg%afh_permeability_Z, dble(xs(1)), dble(xs(2)), dble(xs(3)), 0)

                        r_base_permeability(1) = buffer(1)
#                   endif

                    !convert from mD to m^2 to um^2
                    r_base_permeability = r_base_permeability * 9.869233e-16_SR / (cfg%scaling ** 2)
                else
                    r_base_permeability = 0.0_SR
                end if
#           else
                r_base_permeability = 5000.0_SR * 9.869233e-16_SR / (cfg%scaling ** 2)
#           endif

#           if defined(_ASAGI_TIMING)
                section%stats%r_asagi_time = section%stats%r_asagi_time + get_wtime()
#           endif
		end function

		function get_porosity(section, x) result(porosity)
 			type(t_grid_section), intent(inout)					:: section
			real (kind = GRID_SR), dimension(:), intent(in)		:: x						!< position in world coordinates
			real (kind = GRID_SR)								:: porosity		        !< porosity

            real (kind = GRID_SR)                               :: xs(3)

            xs(1:2) = cfg%scaling * x(1:2) + cfg%offset
            xs(3) = cfg%scaling * max(1, _DARCY_LAYERS) * cfg%dz * x(3)

            assert_ge(x(1), 0.0); assert_ge(x(2), 0.0)
            assert_le(x(1), 1.0); assert_le(x(2), 1.0)

#           if defined(_ASAGI_TIMING)
                section%stats%r_asagi_time = section%stats%r_asagi_time - get_wtime()
#           endif

#           if defined(_ASAGI)
                if (grid_min_x(cfg%afh_porosity) <= xs(1) .and. grid_min_y(cfg%afh_porosity) <= xs(2) &
                        .and. xs(1) <= grid_max_x(cfg%afh_porosity) .and. xs(2) <= grid_max_y(cfg%afh_porosity)) then

                    porosity = grid_get_float_3d(cfg%afh_porosity, dble(xs(1)), dble(xs(2)), dble(xs(3)), 0)
                else
                    porosity = 1.0e-8_SR
                end if
#           else
                porosity = 0.2_SR
#           endif

#           if defined(_ASAGI_TIMING)
                section%stats%r_asagi_time = section%stats%r_asagi_time + get_wtime()
#           endif
		end function

		elemental subroutine pressure_pre_dof_op(p_prod, p, r, d, A_d)
 			real (kind = GRID_SR), intent(in)					:: p_prod
			real (kind = GRID_SR), intent(out)					:: p
			real (kind = GRID_SR), intent(out)					:: r
			real (kind = GRID_SR), intent(out)					:: d
			real (kind = GRID_SR), intent(out)					:: A_d

			p = p_prod
			r = 0.0_SR
			d = 0.0_SR
			A_d = 0.0_SR
		end subroutine
	END MODULE

	MODULE Darcy_initialize_saturation
		use SFC_edge_traversal
		use Samoa_darcy

		implicit none

        type num_traversal_data
            integer (kind = GRID_DI)			:: i_refinements_issued
        end type

		type(darcy_gv_saturation)				:: gv_saturation
		type(darcy_gv_p)						:: gv_p
		type(darcy_gv_rhs)						:: gv_rhs
		type(darcy_gv_is_dirichlet_boundary)    :: gv_is_dirichlet
		type(darcy_gv_volume)					:: gv_inflow

		public initialize_rhs

#		define	_GT_NAME						t_darcy_init_saturation_traversal

#		define _GT_NODES

#		define _GT_NODE_FIRST_TOUCH_OP		    node_first_touch_op
#		define _GT_NODE_LAST_TOUCH_OP		    node_last_touch_op

#		define	_GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define	_GT_PRE_TRAVERSAL_OP			pre_traversal_op
#		define	_GT_ELEMENT_OP					element_op

#		define _GT_NODE_MERGE_OP		        node_merge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

		subroutine post_traversal_grid_op(traversal, grid)
 			type(t_darcy_init_saturation_traversal), intent(inout)      :: traversal
 			type(t_grid), intent(inout)							        :: grid

			call reduce(traversal%i_refinements_issued, traversal%children%i_refinements_issued, MPI_SUM, .true.)
		end subroutine

		subroutine pre_traversal_op(traversal, section)
 			type(t_darcy_init_saturation_traversal), intent(inout)      :: traversal
 			type(t_grid_section), intent(inout)							:: section

			!this variable will be incremented for each cell with a refinement request
			traversal%i_refinements_issued = 0
		end subroutine

		!******************
		!Geometry operators
		!******************

		elemental subroutine node_first_touch_op(traversal, section, node)
 			type(t_darcy_init_saturation_traversal), intent(in)                     :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			call flow_pre_dof_op(node%data_pers%saturation, node%data_pers%rhs, node%data_temp%is_dirichlet_boundary, node%data_temp%volume)
		end subroutine

		subroutine element_op(traversal, section, element)
 			type(t_darcy_init_saturation_traversal)                 :: traversal
 			type(t_grid_section), intent(inout)						:: section
			type(t_element_base), intent(inout)				        :: element

			real (kind = GRID_SR)   :: saturation(_DARCY_LAYERS + 1, 3)
			real (kind = GRID_SR)   :: p(_DARCY_LAYERS + 1, 3)
			real (kind = GRID_SR)   :: rhs(_DARCY_LAYERS + 1, 3)

			call gv_saturation%read_from_element(element, saturation)
			call gv_p%read_from_element(element, p)

			!call element operator
			call alpha_volume_op(traversal, element, saturation, p, rhs, element%cell%data_pers%base_permeability)

			call gv_rhs%add_to_element(element, rhs)
		end subroutine

		elemental subroutine node_merge_op(local_node, neighbor_node)
			type(t_node_data), intent(inout)		:: local_node
 			type(t_node_data), intent(in)		    :: neighbor_node

            where (neighbor_node%data_temp%is_dirichlet_boundary)
                local_node%data_pers%p = neighbor_node%data_pers%p
            end where

			local_node%data_pers%saturation = local_node%data_pers%saturation + neighbor_node%data_pers%saturation
			local_node%data_temp%is_dirichlet_boundary = local_node%data_temp%is_dirichlet_boundary .or. neighbor_node%data_temp%is_dirichlet_boundary
			local_node%data_pers%rhs = local_node%data_pers%rhs + neighbor_node%data_pers%rhs
			local_node%data_temp%volume = local_node%data_temp%volume + neighbor_node%data_temp%volume
		end subroutine

		elemental subroutine node_last_touch_op(traversal, section, node)
 			type(t_darcy_init_saturation_traversal), intent(in)                     :: traversal
 			type(t_grid_section), intent(in)							:: section
			type(t_node_data), intent(inout)			:: node

			real (kind = GRID_SR) :: total_inflow

            total_inflow = tiny(1.0_SR) + sum(node%data_temp%volume)

			call flow_post_dof_op(node%data_pers%saturation, node%data_pers%rhs, node%data_temp%is_dirichlet_boundary, node%data_temp%volume, total_inflow)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine flow_pre_dof_op(saturation, rhs, is_dirichlet, inflow)
			real (kind = GRID_SR), intent(out)					:: saturation
			real (kind = GRID_SR), intent(out)					:: rhs
			logical, intent(out)			                    :: is_dirichlet
			real (kind = GRID_SR), intent(out)					:: inflow

			saturation = 0.0_SR
			rhs = 0.0_SR
			is_dirichlet = .false.
			inflow = 0.0_SR
		end subroutine

		elemental subroutine flow_post_dof_op(saturation, rhs, is_dirichlet, inflow, total_inflow)
			real (kind = GRID_SR), intent(inout)				:: saturation
			real (kind = GRID_SR), intent(inout)				:: rhs
			real (kind = GRID_SR), intent(in)				    :: inflow
			real (kind = GRID_SR), intent(in)				    :: total_inflow
			logical, intent(in)			                        :: is_dirichlet

            saturation = min(1.0_SR, saturation)

            rhs = rhs + cfg%r_inflow * inflow / total_inflow

            if (is_dirichlet) then
                rhs = 0.0_SR
            end if
		end subroutine

		subroutine alpha_volume_op(traversal, element, saturation, p, rhs, base_permeability)
 			type(t_darcy_init_saturation_traversal)                             :: traversal
			type(t_element_base), intent(inout)									:: element
            real (kind = GRID_SR), intent(inout)		            :: saturation(:, :)
            real (kind = GRID_SR), intent(inout)			        :: p(:, :)
            real (kind = GRID_SR), intent(out)				        :: rhs(:, :)

			real (kind = GRID_SR)					                :: pos_prod(2), pos_in(2), radius
			integer (kind = GRID_SI)	                            :: i_depth
			logical 								                :: l_refine_initial, l_refine_solution, l_relevant

#           if (_DARCY_LAYERS > 0)
                real (kind = GRID_SR), intent(inout)                :: base_permeability(:, :)

                l_relevant = any(base_permeability > 0.0_GRID_SR)
#           else
                real (kind = GRID_SR), intent(in)                   :: base_permeability

                l_relevant = (base_permeability > 0.0_GRID_SR)
#           endif

            !set boundary conditions and source terms

            radius = cfg%r_well_radius / element%transform_data%custom_data%scaling

            l_refine_initial = .false.


            !saturation = 1.0_SR
            !call gv_saturation%write_to_element(element, saturation)

            pos_in = samoa_world_to_barycentric_point(element%transform_data, cfg%r_pos_in)
            if (norm2(pos_in) < 1.0_SR + radius) then
                if (norm2(pos_in - 0.5_SR) < sqrt(0.5_SR) + radius) then
                    l_refine_initial = .true.
                end if
            end if

            pos_prod = sign(cfg%r_pos_prod - 0.5_SR, element%transform_data%custom_data%offset - 0.5_SR) + 0.5_SR
            pos_prod = samoa_world_to_barycentric_point(element%transform_data, pos_prod)
            if (norm2(pos_prod) < 1.0_SR + radius) then
                if (norm2(pos_prod - 0.5_SR) < sqrt(0.5_SR) + radius) then
                    l_refine_initial = .true.
                end if
            end if

#           if !defined(_ASAGI)
                pos_in = samoa_barycentric_to_world_point(element%transform_data, [1.0_SR / 3.0_SR, 1.0_SR / 3.0_SR])

                if (pos_in(1) < 0.5_SR) then
                    saturation = 1.0_SR
                    call gv_saturation%add_to_element(element, saturation)
                end if

                if (pos_in(1) < 0.5_SR + 1.5_SR * element%transform_data%custom_data%scaling .and. pos_in(1) > 0.5_SR - 1.5_SR * element%transform_data%custom_data%scaling) then
                    l_refine_initial = .true.
                end if
#           endif

#           if (_DARCY_LAYERS > 0)
                call initialize_rhs(element, saturation, p, rhs, base_permeability)
#           else
                call initialize_rhs(element, saturation(1, :), p(1, :), rhs(1, :), base_permeability)
#           endif

			!check refinement condition

			l_refine_solution = max(maxval(abs(saturation(:, 3) - saturation(:, 2))), maxval(abs(saturation(:, 1) - saturation(:, 2)))) > 0.05_SR
			l_refine_solution = l_refine_solution .or. max(maxval(abs(p(:, 3) - p(:, 2))), maxval(abs(p(:, 1) - p(:, 2)))) > 0.003_SR * cfg%r_p_prod

			!refine the cell if necessary (no coarsening in the initialization!)

			i_depth = element%cell%geometry%i_depth

			if (i_depth < cfg%i_max_depth .and. ((l_relevant .and. (i_depth < cfg%i_min_depth .or. l_refine_solution)) .or. l_refine_initial)) then
				element%cell%geometry%refinement = 1
				traversal%i_refinements_issued = traversal%i_refinements_issued + 1_DI
			else
				element%cell%geometry%refinement = 0
			end if
		end subroutine

#       if (_DARCY_LAYERS > 0)
            subroutine initialize_rhs(element, saturation, p, rhs, base_permeability)
                type(t_element_base), intent(inout)				                    :: element
                real (kind = GRID_SR), intent(inout)		                        :: saturation(:, :)
                real (kind = GRID_SR), intent(inout)			                    :: p(:, :)
                real (kind = GRID_SR), intent(out)									:: rhs(:, :)
                real (kind = GRID_SR), intent(inout)                                :: base_permeability(:, :)

                real (kind = GRID_SR)					            :: g_local(3), pos_prod(2), pos_in(2), radius, weights(3), edge_length, surface, dz, permeability_sum
                real (kind = GRID_SR)                               :: lambda_w(_DARCY_LAYERS + 1, 3), lambda_n(_DARCY_LAYERS + 1, 3), lambda_t(_DARCY_LAYERS, 7)
                real (kind = GRID_SR)                               :: inflow(_DARCY_LAYERS + 1, 3)
                integer                                             :: i, layer
                logical		                                        :: is_dirichlet(_DARCY_LAYERS + 1, 3)

                rhs = 0.0_SR

                radius = cfg%r_well_radius / element%transform_data%custom_data%scaling

                pos_in = samoa_world_to_barycentric_point(element%transform_data, cfg%r_pos_in)

                lambda_w = (saturation * saturation) / cfg%r_nu_w
                lambda_n = (1.0_SR - saturation) * (1.0_SR - saturation) / cfg%r_nu_n

                if (norm2(pos_in) < 1.0_SR + epsilon(1.0_SR)) then
                    if (norm2(pos_in - 0.5_SR) < sqrt(0.5_SR) + epsilon(1.0_SR)) then
                        !injection well:
                        !set an inflow pressure condition and a constant saturation condition

                        !assume the whole well is filled with water
                        weights = [pos_in(1), 1.0_SR - (pos_in(1) + pos_in(2)), pos_in(2)]
                        saturation(:, 1) = max(0.0_SR, min(1.0_SR, saturation(:, 1) + weights(1)))
                        saturation(:, 2) = max(0.0_SR, min(1.0_SR, saturation(:, 2) + weights(2)))
                        saturation(:, 3) = max(0.0_SR, min(1.0_SR, saturation(:, 3) + weights(3)))
                        call gv_saturation%add_to_element(element, spread(weights, 1, _DARCY_LAYERS + 1))

                        lambda_w = (saturation * saturation) / cfg%r_nu_w
                        lambda_n = (1.0_SR - saturation) * (1.0_SR - saturation) / cfg%r_nu_n

                        !The inflow condition is given in um^3 / s
                        !If we devide this by the number of vertical layers, we obtain the 3D inflow for a vertical dual cell column
                        !Split the inflow over all primary cells in each layer that share the dual cell column

!#                       define _DARCY_INJ_TOP_INFLOW
#                       define _DARCY_INJ_ALL_INFLOW
!#                       define _DARCY_INJ_TOP_PRESSURE
!#                       define _DARCY_INJ_ALL_PRESSURE
#                       if defined(_DARCY_INJ_TOP_INFLOW)
                            weights = cfg%r_inflow / 8.0_SR * [pos_in(1), 2.0_SR - 2.0_SR * (pos_in(1) + pos_in(2)), pos_in(2)]

                            !set inflow condition on the top layer
                            rhs(_DARCY_LAYERS + 1, :) = rhs(_DARCY_LAYERS + 1, :) + weights

                            !set a free flow permeability in the well
                            base_permeability(:, 2) = 1.0e6 * 9.869233e-16_SR / (cfg%scaling ** 2)
#                       elif defined(_DARCY_INJ_ALL_INFLOW)
                            !Using Peaceman's well model we consider the well as an internal boundary and assume that
                            !near the well the following conditions hold:
                            !1) K_z is huge => p_z = rho_w g because lambda_w(S) K_z (-p_z + rho_w g) < Q / (r^2 pi) << K_z
                            !2) p is radially symmetric and the derivative p_r is constant.
                            !
                            !Thus the inflow q is q(r,phi,z) = lambda_w(S) K_r(r,phi,z) (-p_r)
                            !With \integral_{well boundary} q(r,phi,z) * dS = Q we obtain
                            ! \integral_{well boundary} lambda_w(S) K_r(r,phi,z) (-p_r) dOmega = Q
                            ! p_r = - Q / integral_(well boundary) lambda_w(S) K_r(r,phi,z) dOmega
                            ! q(r,phi,z) = Q K_r(r,phi,z) / integral_(well boundary) K_r(r,phi,z) dOmega

                            weights = 0.25_SR * [pos_in(1), 2.0_SR - 2.0_SR * (pos_in(1) + pos_in(2)), pos_in(2)]

                            inflow = 0.0_SR

                            !split local inflows in primary layer and assign half contributions to dual layers
                            do i = 1, 3
                                inflow(1:_DARCY_LAYERS, i) = inflow(1:_DARCY_LAYERS, i) + 0.5_SR * base_permeability(:, 1) * weights(i)
                                inflow(2:_DARCY_LAYERS + 1, i) = inflow(2:_DARCY_LAYERS + 1, i) + 0.5_SR * base_permeability(:, 1) * weights(i)
                            end do

                            call gv_inflow%add_to_element(element, inflow)
#                       elif defined(_DARCY_INJ_TOP_PRESSURE)
                            weights = [pos_in(1), 1.0_SR - (pos_in(1) + pos_in(2)), pos_in(2)]

                            is_dirichlet(_DARCY_LAYERS + 1, 1) = weights(1) > epsilon(1.0_SR)
                            is_dirichlet(_DARCY_LAYERS + 1, 2) = weights(2) > epsilon(1.0_SR)
                            is_dirichlet(_DARCY_LAYERS + 1, 3) = weights(3) > epsilon(1.0_SR)

                            do i = 1, 3
                                if (is_dirichlet(_DARCY_LAYERS + 1, i)) then
                                    p(_DARCY_LAYERS + 1, i) = cfg%r_p_in
                                end if
                            end do

                            !set a free flow permeability in the well
                            base_permeability(:, 2) = 1.0e6 * 9.869233e-16_SR / (cfg%scaling ** 2)

                            call gv_p%write_to_element(element, p)
                            call gv_is_dirichlet%add_to_element(element, is_dirichlet)
#                       elif defined(_DARCY_INJ_ALL_PRESSURE)
                            weights = [pos_in(1), 1.0_SR - (pos_in(1) + pos_in(2)), pos_in(2)]

                            is_dirichlet(:, 1) = weights(1) > epsilon(1.0_SR)
                            is_dirichlet(:, 2) = weights(2) > epsilon(1.0_SR)
                            is_dirichlet(:, 3) = weights(3) > epsilon(1.0_SR)

                            do i = 1, 3
                                if (is_dirichlet(_DARCY_LAYERS + 1, i)) then
                                    p(_DARCY_LAYERS + 1, i) = cfg%r_p_in

                                    do layer = _DARCY_LAYERS, 1, -1
                                        p(layer, i) = p(layer + 1, i) - &
                                            (lambda_w(layer, i) * cfg%r_rho_w + lambda_n(layer, i) * cfg%r_rho_n) &
                                            / (lambda_w(layer, i) + lambda_n(layer, i)) * cfg%dz * g(3)
                                    end do
                                end if
                            end do

                            call gv_p%write_to_element(element, p)
                            call gv_is_dirichlet%add_to_element(element, is_dirichlet)
#                       else
#                           error Injection condition must be defined!
#                       endif
                    end if
                end if

                pos_prod = sign(cfg%r_pos_prod - 0.5_SR, element%transform_data%custom_data%offset - 0.5_SR) + 0.5_SR
                pos_prod = samoa_world_to_barycentric_point(element%transform_data, pos_prod)

                if (norm2(pos_prod) < 1.0_SR + radius) then
                    if (norm2(pos_prod - 0.5_SR) < sqrt(0.5_SR) + radius) then
                        !production well:
                        !set a constant pressure condition and an outflow saturation condition

                        weights = [pos_prod(1), 1.0_SR - (pos_prod(1) + pos_prod(2)), pos_prod(2)]

!#                       define _DARCY_PROD_TOP_PRESSURE
#                       define _DARCY_PROD_ALL_PRESSURE
#                       if defined(_DARCY_PROD_TOP_PRESSURE)
                            is_dirichlet(_DARCY_LAYERS + 1, 1) = weights(1) > epsilon(1.0_SR)
                            is_dirichlet(_DARCY_LAYERS + 1, 2) = weights(2) > epsilon(1.0_SR)
                            is_dirichlet(_DARCY_LAYERS + 1, 3) = weights(3) > epsilon(1.0_SR)

                            do i = 1, 3
                                if (is_dirichlet(_DARCY_LAYERS + 1, i)) then
                                    p(_DARCY_LAYERS + 1, i) = cfg%r_p_prod
                                end if
                            end do

                            !set a free flow permeability in the well
                            base_permeability(:, 2) = 1.0e6 * 9.869233e-16_SR / (cfg%scaling ** 2)
#                       elif defined(_DARCY_PROD_ALL_PRESSURE)
                            is_dirichlet(:, 1) = weights(1) > epsilon(1.0_SR)
                            is_dirichlet(:, 2) = weights(2) > epsilon(1.0_SR)
                            is_dirichlet(:, 3) = weights(3) > epsilon(1.0_SR)

                            do i = 1, 3
                                if (is_dirichlet(_DARCY_LAYERS + 1, i)) then
                                    p(_DARCY_LAYERS + 1, i) = cfg%r_p_prod

                                    do layer = _DARCY_LAYERS, 1, -1
                                        p(layer, i) = p(layer + 1, i) - &
                                            (lambda_w(layer, i) * cfg%r_rho_w + lambda_n(layer, i) * cfg%r_rho_n) &
                                            / (lambda_w(layer, i) + lambda_n(layer, i)) * cfg%dz * g(3)
                                    end do
                                end if
                            end do
#                       else
#                           error Production condition must be defined!
#                       endif

                        call gv_p%write_to_element(element, p)
                        call gv_is_dirichlet%add_to_element(element, is_dirichlet)
                    end if
                end if

                !rotate g so it points in the right direction (no scaling!)
                g_local = g
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))
                g_local(1:2) = g_local(1:2) / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))

                edge_length = element%cell%geometry%get_leg_size()
                surface = element%cell%geometry%get_volume()
                dz = cfg%dz

                do i = 1, _DARCY_LAYERS
                    call compute_rhs_1D(edge_length, 0.25_SR * edge_length * dz, base_permeability(i, 1), p(i, 2), p(i, 1), lambda_w(i, 2), lambda_w(i, 1), lambda_n(i, 2), lambda_n(i, 1), g_local(1), rhs(i, 2), rhs(i, 1), lambda_t(i, 1))
                    call compute_rhs_1D(edge_length, 0.25_SR * edge_length * dz, base_permeability(i, 1), p(i, 2), p(i, 3), lambda_w(i, 2), lambda_w(i, 3), lambda_n(i, 2), lambda_n(i, 3), g_local(2), rhs(i, 2), rhs(i, 3), lambda_t(i, 2))

                    call compute_rhs_1D(dz, 0.25_SR * surface, base_permeability(i, 2), p(i, 1), p(i + 1, 1), lambda_w(i, 1), lambda_w(i + 1, 1), lambda_n(i, 1), lambda_n(i + 1, 1), g_local(3), rhs(i, 1), rhs(i + 1, 1), lambda_t(i, 3))
                    call compute_rhs_1D(dz, 0.50_SR * surface, base_permeability(i, 2), p(i, 2), p(i + 1, 2), lambda_w(i, 2), lambda_w(i + 1, 2), lambda_n(i, 2), lambda_n(i + 1, 2), g_local(3), rhs(i, 2), rhs(i + 1, 2), lambda_t(i, 4))
                    call compute_rhs_1D(dz, 0.25_SR * surface, base_permeability(i, 2), p(i, 3), p(i + 1, 3), lambda_w(i, 3), lambda_w(i + 1, 3), lambda_n(i, 3), lambda_n(i + 1, 3), g_local(3), rhs(i, 3), rhs(i + 1, 3), lambda_t(i, 5))

                    call compute_rhs_1D(edge_length, 0.25_SR * edge_length * dz, base_permeability(i, 1), p(i + 1, 2), p(i + 1, 1), lambda_w(i + 1, 2), lambda_w(i + 1, 1), lambda_n(i + 1, 2), lambda_n(i + 1, 1), g_local(1), rhs(i + 1, 2), rhs(i + 1, 1), lambda_t(i, 6))
                    call compute_rhs_1D(edge_length, 0.25_SR * edge_length * dz, base_permeability(i, 1), p(i + 1, 2), p(i + 1, 3), lambda_w(i + 1, 2), lambda_w(i + 1, 3), lambda_n(i + 1, 2), lambda_n(i + 1, 3), g_local(2), rhs(i + 1, 2), rhs(i + 1, 3), lambda_t(i, 7))
                end do

                if (element%transform_data%plotter_data%orientation > 0) then
                    element%cell%data_pers%lambda_t = lambda_t
                else
                    element%cell%data_pers%lambda_t(:, 1) = lambda_t(:, 2)
                    element%cell%data_pers%lambda_t(:, 2) = lambda_t(:, 1)
                    element%cell%data_pers%lambda_t(:, 3) = lambda_t(:, 5)
                    element%cell%data_pers%lambda_t(:, 4) = lambda_t(:, 4)
                    element%cell%data_pers%lambda_t(:, 5) = lambda_t(:, 3)
                    element%cell%data_pers%lambda_t(:, 6) = lambda_t(:, 7)
                    element%cell%data_pers%lambda_t(:, 7) = lambda_t(:, 6)
                end if
            end subroutine
#       else
            subroutine initialize_rhs(element, saturation, p, rhs, base_permeability)
                type(t_element_base), intent(inout)				                    :: element
                real (kind = GRID_SR), intent(inout)		                        :: saturation(:)
                real (kind = GRID_SR), intent(inout)			                    :: p(:)
                real (kind = GRID_SR), intent(out)									:: rhs(:)
                real (kind = GRID_SR), intent(in)                                   :: base_permeability

                real (kind = GRID_SR)					            :: pos_prod(2), pos_in(2), radius, inflow, edge_length
                real (kind = GRID_SR)                               :: g_local(2), lambda_w(3), lambda_n(3), lambda_t(2)
                logical		                                        :: is_dirichlet(3)

                rhs = 0.0_SR

                radius = cfg%r_well_radius / element%transform_data%custom_data%scaling

                pos_in = samoa_world_to_barycentric_point(element%transform_data, cfg%r_pos_in)

                if (norm2(pos_in) < 1.0_SR + radius) then
                    if (norm2(pos_in - 0.5_SR) < sqrt(0.5_SR) + radius) then
                        !injection well:
                        !set an inflow pressure condition and a constant saturation condition

                        saturation = max(0.0_SR, min(1.0_SR, saturation + [pos_in(1), 1.0_SR - (pos_in(1) + pos_in(2)), pos_in(2)]))
                        call gv_saturation%add_to_element(element, [pos_in(1), 1.0_SR - (pos_in(1) + pos_in(2)), pos_in(2)])

                        if (base_permeability > 0) then
                            !The inflow condition is given in um^3 / s
                            !If we devide this by the height of the domain cfg%dz, we obtain the 2D inflow in um^2/s
                            !Split the inflow over all primary cells that share the dual cell

                            inflow = cfg%r_inflow / cfg%dz

                            rhs = rhs + inflow / 8.0_SR * [pos_in(1), 2.0_SR - 2.0_SR * (pos_in(1) + pos_in(2)), pos_in(2)]
                        end if
                    end if
                end if

                pos_prod = sign(cfg%r_pos_prod - 0.5_SR, element%transform_data%custom_data%offset - 0.5_SR) + 0.5_SR
                pos_prod = samoa_world_to_barycentric_point(element%transform_data, pos_prod)

                if (norm2(pos_prod) < 1.0_SR + radius) then
                    if (norm2(pos_prod - 0.5_SR) < sqrt(0.5_SR) + radius) then
                        !production well:
                        !set a constant pressure condition and an outflow saturation condition

                        is_dirichlet = [pos_prod(1) > epsilon(1.0_SR), 1.0_SR - (pos_prod(1) + pos_prod(2)) > epsilon(1.0_SR), pos_prod(2) > epsilon(1.0_SR)]

                        where (is_dirichlet)
                            p = cfg%r_p_prod
                        end where

                        call gv_p%write_to_element(element, p)
                        call gv_is_dirichlet%add_to_element(element, is_dirichlet)

                        rhs = 0.0_SR
                    end if
                end if

                lambda_w = (saturation * saturation) / cfg%r_nu_w
                lambda_n = (1.0_SR - saturation) * (1.0_SR - saturation) / cfg%r_nu_n

                !rotate g so it points in the right direction (no scaling!)
                g_local = g
                g_local(1:2) = samoa_world_to_barycentric_normal(element%transform_data, g_local(1:2))
                g_local(1:2) = g_local(1:2) / (element%transform_data%custom_data%scaling * sqrt(abs(element%transform_data%plotter_data%det_jacobian)))

                edge_length = element%cell%geometry%get_leg_size()

                call compute_rhs_1D(edge_length, 0.5_SR * edge_length, base_permeability, p(2), p(1), lambda_w(2), lambda_w(1), lambda_n(2), lambda_n(1), g_local(1), rhs(2), rhs(1), lambda_t(1))
                call compute_rhs_1D(edge_length, 0.5_SR * edge_length, base_permeability, p(2), p(3), lambda_w(2), lambda_w(3), lambda_n(2), lambda_n(3), g_local(2), rhs(2), rhs(3), lambda_t(2))

                if (element%transform_data%plotter_data%orientation > 0) then
                    element%cell%data_pers%lambda_t = lambda_t
                else
                    element%cell%data_pers%lambda_t(1) = lambda_t(2)
                    element%cell%data_pers%lambda_t(2) = lambda_t(1)
                end if
            end subroutine
#       endif

        subroutine compute_rhs_1D(dx, area, base_permeability, pL, pR, lambda_wL, lambda_wR, lambda_nL, lambda_nR, g_local, rhsL, rhsR, lambda_t)
            real (kind = GRID_SR), intent(in)       :: area, dx, base_permeability, pL, pR, lambda_wL, lambda_wR, lambda_nL, lambda_nR, g_local
            real (kind = GRID_SR), intent(inout)    :: rhsL, rhsR, lambda_t

            real (kind = SR)    :: u_w, u_n, lambda_w_local, lambda_n_local, rhs_local

            u_w = area * base_permeability * (-(pR - pL) / dx + cfg%r_rho_w * g_local)
            u_n = area * base_permeability * (-(pR - pL) / dx + cfg%r_rho_n * g_local)

            if (u_w > 0) then
                lambda_w_local = lambda_wL
            else
                lambda_w_local = lambda_wR
            end if

            if (u_n > 0) then
                lambda_n_local = lambda_nL
            else
                lambda_n_local = lambda_nR
            end if

            rhs_local = area * base_permeability * (lambda_w_local * cfg%r_rho_w + lambda_n_local * cfg%r_rho_n) * g_local

            rhsL = rhsL - rhs_local
            rhsR = rhsR + rhs_local
            lambda_t = area / dx * base_permeability * (lambda_w_local + lambda_n_local)
        end subroutine
	END MODULE
#endif
