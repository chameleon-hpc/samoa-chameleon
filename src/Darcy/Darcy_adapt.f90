! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_Adapt
		use SFC_edge_traversal
		use Conformity

		use Samoa_darcy
		use Tools_noise
		use Darcy_initialize_pressure

        implicit none

        type num_traversal_data
        end type

		type(darcy_gv_p)							:: gv_p
		type(darcy_gv_saturation)					:: gv_saturation
		type(darcy_gv_volume)						:: gv_volume
		type(darcy_gv_mat_diagonal)					:: gv_p_volume

#		define	_GT_NAME							t_darcy_adaption_traversal

#		define	_GT_NODES

#		define	_GT_TRANSFER_OP						transfer_op
#		define	_GT_REFINE_OP						refine_op
#		define	_GT_COARSEN_OP						coarsen_op

#		define _GT_CELL_FIRST_TOUCH_OP				cell_first_touch_op
#		define _GT_NODE_FIRST_TOUCH_OP				node_first_touch_op
#		define _GT_NODE_LAST_TOUCH_OP				node_last_touch_op
#		define _GT_INNER_NODE_LAST_TOUCH_OP			inner_node_last_touch_op

#		define _GT_NODE_MERGE_OP				    node_merge_op
#		define _GT_EDGE_MPI_TYPE

#		include "SFC_generic_adaptive_traversal.f90"

        subroutine create_edge_mpi_type(mpi_edge_type)
            integer, intent(out)            :: mpi_edge_type

#           if defined(_MPI)
                type(t_edge_data)                       :: edge
                integer                                 :: blocklengths(2), types(2), disps(2), type_size, i_error
                integer (kind = MPI_ADDRESS_KIND)       :: lb, ub

                blocklengths(1) = 1
                blocklengths(2) = 1

                disps(1) = 0
                disps(2) = sizeof(edge)

                types(1) = MPI_LB
                types(2) = MPI_UB

                call MPI_Type_struct(2, blocklengths, disps, types, mpi_edge_type, i_error); assert_eq(i_error, 0)
                call MPI_Type_commit(mpi_edge_type, i_error); assert_eq(i_error, 0)

                call MPI_Type_size(mpi_edge_type, type_size, i_error); assert_eq(i_error, 0)
                call MPI_Type_get_extent(mpi_edge_type, lb, ub, i_error); assert_eq(i_error, 0)

                assert_eq(0, lb)
                assert_eq(0, type_size)
                assert_eq(sizeof(edge), ub)
#           endif
        end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine transfer_op(traversal, grid, src_element, dest_element)
 			type(t_darcy_adaption_traversal), intent(inout)							:: traversal
 			type(t_grid_section), intent(inout)							            :: grid
			type(t_traversal_element), intent(inout)									:: src_element
			type(t_traversal_element), intent(inout)									:: dest_element

			real (kind = GRID_SR)   :: p_in(_DARCY_LAYERS + 1, 3), saturation_in(_DARCY_LAYERS + 1, 3), porosity(_DARCY_LAYERS + 1), volume(3)
			integer					:: i, level

            !make sure, the effective volume never turns out to be 0
            porosity = epsilon(1.0_SR)

#           if (_DARCY_LAYERS > 0)
                porosity(1 : _DARCY_LAYERS) = porosity(1 : _DARCY_LAYERS) + 0.5_SR * src_element%cell%data_pers%porosity
                porosity(2 : _DARCY_LAYERS + 1) = porosity(2 : _DARCY_LAYERS + 1) + 0.5_SR * src_element%cell%data_pers%porosity
#           else
                porosity = porosity + src_element%cell%data_pers%porosity
#           endif

			call gv_p%read_from_element( src_element%t_element_base, p_in)
			call gv_saturation%read_from_element( src_element%t_element_base, saturation_in)

            volume(1) = src_element%cell%geometry%get_volume() * 0.25_SR
            volume(2) = src_element%cell%geometry%get_volume() * 0.5_SR
            volume(3) = src_element%cell%geometry%get_volume() * 0.25_SR

			call gv_p%add_to_element( dest_element%t_element_base, [p_in(:, 1) * volume(1), p_in(:, 2) * volume(2), p_in(:, 3) * volume(3)])
			call gv_p_volume%add_to_element( dest_element%t_element_base, [spread(volume(1), 1, _DARCY_LAYERS + 1), spread(volume(2), 1, _DARCY_LAYERS + 1), spread(volume(3), 1, _DARCY_LAYERS + 1)])
			call gv_saturation%add_to_element( dest_element%t_element_base, [saturation_in(:, 1) * volume(1) * porosity, saturation_in(:, 2) * volume(2) * porosity, saturation_in(:, 3) * volume(3) * porosity])
			call gv_volume%add_to_element( dest_element%t_element_base, [volume(1) * porosity, volume(2) * porosity, volume(3) * porosity])

			!permeability

			dest_element%cell%data_pers%base_permeability = src_element%cell%data_pers%base_permeability
			dest_element%cell%data_pers%porosity = src_element%cell%data_pers%porosity
		end subroutine

		subroutine refine_op(traversal, grid, src_element, dest_element, refinement_path)
 			type(t_darcy_adaption_traversal), intent(inout)							:: traversal
 			type(t_grid_section), intent(inout)										:: grid
			type(t_traversal_element), intent(inout)								:: src_element
			type(t_traversal_element), intent(inout)								:: dest_element
			integer, intent(in)										                :: refinement_path(:)

			real (kind = GRID_SR)   :: p_in(_DARCY_LAYERS + 1, 3), saturation_in(_DARCY_LAYERS + 1, 3), porosity(_DARCY_LAYERS + 1), volume(3), weights(3, 3), r_cells(4)
            real (kind = GRID_SR)   :: x(3), p(2), v1(2), v2(2), buffer(2), alpha, beta
            integer					:: i, j, k, level, nx, nz, k_start, k_end, no_samples(_DARCY_LAYERS)

            !make sure, the effective volume never turns out to be 0
            porosity = epsilon(1.0_SR)

#           if (_DARCY_LAYERS > 0)
                porosity(1 : _DARCY_LAYERS) = porosity(1 : _DARCY_LAYERS) + 0.5_SR * src_element%cell%data_pers%porosity
                porosity(2 : _DARCY_LAYERS + 1) = porosity(2 : _DARCY_LAYERS + 1) + 0.5_SR * src_element%cell%data_pers%porosity
#           else
                porosity = porosity + src_element%cell%data_pers%porosity
#           endif

			call gv_p%read_from_element( src_element%t_element_base, p_in)
			call gv_saturation%read_from_element( src_element%t_element_base, saturation_in)

			r_cells = [0.0_SR, 0.25_SR, 0.75_SR, 1.0_SR]

            do level = 1, size(refinement_path)
                if (refinement_path(level) > 1) then
                    r_cells = 0.5_SR + r_cells / 2.0_SR
                else
                    r_cells = r_cells / 2.0_SR
                end if
            end do

            weights = transpose(reshape([max(0.0_SR, min(0.25_SR, r_cells(2:4)) - max(0.0_SR, r_cells(1:3))), &
                        max(0.0_SR, min(0.75_SR, r_cells(2:4)) - max(0.25_SR, r_cells(1:3))), &
                        max(0.0_SR, min(1.0_SR, r_cells(2:4)) - max(0.75_SR, r_cells(1:3)))], [3, 3]))

            weights = src_element%cell%geometry%get_volume() * weights

            forall (level = 1:_DARCY_LAYERS + 1, i = 1:3)
                p_in(level, i) = dot_product(weights(:, i), p_in(level, :))
                saturation_in(level, i) = dot_product(weights(:, i), saturation_in(level, :))
            end forall

            volume(1) = sum(weights(:, 1))
            volume(2) = sum(weights(:, 2))
            volume(3) = sum(weights(:, 3))

 			call gv_p%add_to_element(dest_element%t_element_base, p_in)
			call gv_p_volume%add_to_element( dest_element%t_element_base, [spread(volume(1), 1, _DARCY_LAYERS + 1), spread(volume(2), 1, _DARCY_LAYERS + 1), spread(volume(3), 1, _DARCY_LAYERS + 1)])
			call gv_saturation%add_to_element( dest_element%t_element_base, [saturation_in(:, 1) * porosity, saturation_in(:, 2) * porosity, saturation_in(:, 3) * porosity])
			call gv_volume%add_to_element( dest_element%t_element_base, [volume(1) * porosity, volume(2) * porosity, volume(3) * porosity])

			!permeability & porosity

#           if defined(_ADAPT_INTEGRATE)
#               if (_DARCY_LAYERS > 0)
#                   if defined(_ASAGI)
                        nx = max(1, int(cfg%scaling * get_edge_size(dest_element%cell%geometry%i_depth) / asagi_grid_delta(cfg%afh_permeability_X, 1)))
                        nz = max(1, int(cfg%scaling * cfg%dz * real(_DARCY_LAYERS, SR) / asagi_grid_delta(cfg%afh_permeability_X, 2)))
#                   else
                        nx = max(1, int(512.0_SR * get_edge_size(dest_element%cell%geometry%i_depth)))
                        nz = 1
#                   endif

                    p = samoa_barycentric_to_world_point(dest_element%transform_data, [0.0_SR, 0.0_SR])
                    v1 = samoa_barycentric_to_world_vector(dest_element%transform_data, [0.5_SR, 0.5_SR])
                    v2 = samoa_barycentric_to_world_vector(dest_element%transform_data, [0.5_SR, -0.5_SR])

                    dest_element%cell%data_pers%base_permeability = 0.0_SR
                    dest_element%cell%data_pers%porosity = 0.0_SR

                    no_samples = 0

                    do i = 0, nx - 1
                        alpha = (i + 0.5_SR) / real(nx, SR)

                        do j = -i, i
                            beta = real(j, SR) / real(nx, SR)
                            x(1:2) = p + alpha * v1 + beta * v2

                            do level = 1, _DARCY_LAYERS
                                !find the k-range in the source data that we have to read from
                                !and ensure that we read at least one value from the source data

                                k_start = ((level - 1) * nz) / _DARCY_LAYERS + 1
                                k_end = max(k_start, (level * nz) / _DARCY_LAYERS)

                                do k = k_start, k_end
                                    x(3) = (real(k, SR) - 0.5_SR) / real(nz, SR)
                                    buffer = get_base_permeability(grid, x, dest_element%cell%geometry%i_depth / 2_GRID_SI)

                                    dest_element%cell%data_pers%base_permeability(level, :) = dest_element%cell%data_pers%base_permeability(level, :) + mean_transform(buffer)
                                    dest_element%cell%data_pers%porosity(level) = dest_element%cell%data_pers%porosity(level) + get_porosity(grid, x)
                                end do

                                no_samples(level) = no_samples(level) + k_end - k_start + 1
                            end do
                        end do
                    end do

                    dest_element%cell%data_pers%base_permeability(:, 1) = mean_invert(dest_element%cell%data_pers%base_permeability(:, 1) / no_samples)
                    dest_element%cell%data_pers%base_permeability(:, 2) = mean_invert(dest_element%cell%data_pers%base_permeability(:, 2) / no_samples)

                    dest_element%cell%data_pers%porosity = dest_element%cell%data_pers%porosity / no_samples
#               else
#                   if defined(_ASAGI)
                        nx = max(1, int(cfg%scaling * dest_element%transform_data%custom_data%scaling / asagi_grid_delta(cfg%afh_permeability_X, 1)))
                        nz = max(1, int(cfg%scaling * cfg%dz / asagi_grid_delta(cfg%afh_permeability_X, 2)))
#                   else
                        nx = max(1, 512 * int(dest_element%transform_data%custom_data%scaling))
                        nz = 1
#                   endif

                    p = samoa_barycentric_to_world_point(dest_element%transform_data, [0.0_SR, 0.0_SR])
                    v1 = samoa_barycentric_to_world_vector(dest_element%transform_data, [0.5_SR, 0.5_SR])
                    v2 = samoa_barycentric_to_world_vector(dest_element%transform_data, [0.5_SR, -0.5_SR])

                    dest_element%cell%data_pers%base_permeability = 0.0_SR
                    dest_element%cell%data_pers%porosity = 0.0_SR

                    do i = 0, nx - 1
                        alpha = (i + 0.5_SR) / real(nx, SR)

                        do j = -i, i
                            beta = real(j, SR) / real(nx, SR)
                            x(1:2) = p + alpha * v1 + beta * v2

                            !find the k-range in the source data that we have to read from
                            !and ensure that we read at least one value from the source data

                            k_start = 1
                            k_end = 1

                            do k = k_start, k_end
                                x(3) = (real(k, SR) - 0.5_SR) / real(nz, SR)
                                buffer = get_base_permeability(grid, x, dest_element%cell%geometry%i_depth / 2_GRID_SI)

                                dest_element%cell%data_pers%base_permeability = dest_element%cell%data_pers%base_permeability + mean_transform(buffer(1))
                                dest_element%cell%data_pers%porosity = dest_element%cell%data_pers%porosity + get_porosity(grid, x)
                            end do

                            no_samples(1) = no_samples(1) + k_end - k_start + 1
                        end do
                    end do

                    dest_element%cell%data_pers%base_permeability = mean_invert(dest_element%cell%data_pers%base_permeability / no_samples(1))
                    dest_element%cell%data_pers%porosity = dest_element%cell%data_pers%porosity / no_samples(1)
#               endif
#           elif defined(_ADAPT_SAMPLE)
                x(1:2) = samoa_barycentric_to_world_point(dest_element%transform_data, [1.0_SR / 3.0_SR, 1.0_SR / 3.0_SR])

#               if (_DARCY_LAYERS > 0)
                    do level = 1, _DARCY_LAYERS
                        x(3) = (level - 0.5_SR) / real(_DARCY_LAYERS, SR)

                        dest_element%cell%data_pers%base_permeability(level, :) = get_base_permeability(grid, x, dest_element%cell%geometry%i_depth / 2_GRID_SI)
                        dest_element%cell%data_pers%porosity(level) = get_porosity(grid, x)
                    end do
#               else
                    x(3) = 1.0e-5_SR

                    dest_element%cell%data_pers%base_permeability = get_base_permeability(grid, x, dest_element%cell%geometry%i_depth / 2_GRID_SI)
                    dest_element%cell%data_pers%porosity = get_porosity(grid, x)
#               endif
#           endif
		end subroutine

        !> converts the value p to a transformed space where the respective mean is additive
		elemental function mean_transform(p) result(p_trans)
            real (kind = SR), intent(in)    :: p
            real (kind = SR)                :: p_trans

#           if defined(_PERM_MEAN_ARITHMETIC)
                p_trans = p
#           elif defined(_PERM_MEAN_HARMONIC)
                !this will cause an intended floating point overflow if p = 0. Backtranformation will return the value 0.
                p_trans = 1.0_SR / p
#           elif defined(_PERM_MEAN_GEOMETRIC)
                !this will cause an intended floating point overflow if p = 0. Backtranformation will return the value 0.
                p_trans = log(p)
#           endif
		end function

        !> converts the value p_trans back to normal space
        elemental function mean_invert(p_trans) result(p)
            real (kind = SR), intent(in)    :: p_trans
            real (kind = SR)                :: p

#           if defined(_PERM_MEAN_ARITHMETIC)
                p = p_trans
#           elif defined(_PERM_MEAN_HARMONIC)
                p = 1.0_SR / p_trans
#           elif defined(_PERM_MEAN_GEOMETRIC)
                p = exp(p_trans)
#           endif
		end function

		subroutine coarsen_op(traversal, grid, src_element, dest_element, coarsening_path)
  			type(t_darcy_adaption_traversal), intent(inout)							    :: traversal
			type(t_grid_section), intent(inout)										    :: grid
			type(t_traversal_element), intent(inout)									:: src_element
			type(t_traversal_element), intent(inout)									:: dest_element
			integer, dimension(:), intent(in)											:: coarsening_path

			real (kind = GRID_SR)   :: p_in(_DARCY_LAYERS + 1, 3), saturation_in(_DARCY_LAYERS + 1, 3), porosity(_DARCY_LAYERS + 1), volume(3), weights(3, 3), r_cells(4)
			integer					:: i, level

            !make sure, the effective volume never reaches 0
            porosity = epsilon(1.0_SR)

#           if (_DARCY_LAYERS > 0)
                porosity(1 : _DARCY_LAYERS) = porosity(1 : _DARCY_LAYERS) + 0.5_SR * src_element%cell%data_pers%porosity
                porosity(2 : _DARCY_LAYERS + 1) = porosity(2 : _DARCY_LAYERS + 1) + 0.5_SR * src_element%cell%data_pers%porosity
#           else
                porosity = porosity + src_element%cell%data_pers%porosity
#           endif

			call gv_p%read_from_element( src_element%t_element_base, p_in)
			call gv_saturation%read_from_element( src_element%t_element_base, saturation_in)

			r_cells = [0.0_SR, 0.25_SR, 0.75_SR, 1.0_SR]

            do level = 1, size(coarsening_path)
                if (coarsening_path(level) > 1) then
                    r_cells = 0.5_SR + r_cells / 2.0_SR
                else
                    r_cells = r_cells / 2.0_SR
                end if
            end do

            weights = transpose(reshape([max(0.0_SR, min(0.25_SR, r_cells(2:4)) - max(0.0_SR, r_cells(1:3))), &
                        max(0.0_SR, min(0.75_SR, r_cells(2:4)) - max(0.25_SR, r_cells(1:3))), &
                        max(0.0_SR, min(1.0_SR, r_cells(2:4)) - max(0.75_SR, r_cells(1:3)))], [3, 3]))

            weights = dest_element%cell%geometry%get_volume() * weights

            forall (level = 1:_DARCY_LAYERS + 1, i = 1:3)
                p_in(level, i) = dot_product(weights(:, i), p_in(level, :))
                saturation_in(level, i) = dot_product(weights(:, i), saturation_in(level, :))
            end forall

            volume(1) = sum(weights(:, 1))
            volume(2) = sum(weights(:, 2))
            volume(3) = sum(weights(:, 3))

 			call gv_p%add_to_element(dest_element%t_element_base, p_in)
			call gv_p_volume%add_to_element( dest_element%t_element_base, [spread(volume(1), 1, _DARCY_LAYERS + 1), spread(volume(2), 1, _DARCY_LAYERS + 1), spread(volume(3), 1, _DARCY_LAYERS + 1)])
			call gv_saturation%add_to_element( dest_element%t_element_base, [saturation_in(:, 1) * porosity, saturation_in(:, 2) * porosity, saturation_in(:, 3) * porosity])
			call gv_volume%add_to_element( dest_element%t_element_base, [volume(1) * porosity, volume(2) * porosity, volume(3) * porosity])

            !permeability and porosity: compute the average

            dest_element%cell%data_pers%base_permeability = mean_invert(mean_transform(dest_element%cell%data_pers%base_permeability) + (0.5_SR ** size(coarsening_path) * mean_transform(src_element%cell%data_pers%base_permeability)))
            dest_element%cell%data_pers%porosity = dest_element%cell%data_pers%porosity + (0.5 ** size(coarsening_path)) * src_element%cell%data_pers%porosity
		end subroutine

		!first touch ops
		elemental subroutine cell_first_touch_op(traversal, grid, cell)
 			type(t_darcy_adaption_traversal), intent(in)							:: traversal
 			type(t_grid_section), intent(in)							:: grid
			type(t_cell_data_ptr), intent(inout)				:: cell

            cell%data_pers%base_permeability = mean_invert(0.0_SR)
			cell%data_pers%porosity = 0.0_SR
		end subroutine

		elemental subroutine node_first_touch_op(traversal, grid, node)
 			type(t_darcy_adaption_traversal), intent(in)							:: traversal
 			type(t_grid_section), intent(in)							:: grid
			type(t_node_data), intent(inout)				:: node

			call pre_dof_op(node%data_pers%saturation, node%data_temp%volume, node%data_pers%p, node%data_temp%mat_diagonal, node%data_pers%d, node%data_pers%A_d, node%data_pers%rhs)
		end subroutine

		!merge ops

		elemental subroutine node_merge_op(local_node, neighbor_node)
 			type(t_node_data), intent(inout)			    :: local_node
			type(t_node_data), intent(in)				    :: neighbor_node

            local_node%data_pers%p = local_node%data_pers%p + neighbor_node%data_pers%p
            local_node%data_pers%saturation = local_node%data_pers%saturation + neighbor_node%data_pers%saturation
            local_node%data_temp%volume = local_node%data_temp%volume + neighbor_node%data_temp%volume
            local_node%data_temp%mat_diagonal = local_node%data_temp%mat_diagonal + neighbor_node%data_temp%mat_diagonal
        end subroutine

		!last touch ops

		subroutine node_last_touch_op(traversal, section, nodes)
 			type(t_darcy_adaption_traversal), intent(in)	    :: traversal
 			type(t_grid_section), intent(inout)				    :: section
			type(t_node_data), intent(inout)					:: nodes(:)

            integer :: i

            do i = 1, size(nodes)
                call inner_node_last_touch_op(traversal, section, nodes(i))
            end do
		end subroutine

		subroutine inner_node_last_touch_op(traversal, section, node)
 			type(t_darcy_adaption_traversal), intent(in)    :: traversal
 			type(t_grid_section), intent(inout)			    :: section
			type(t_node_data), intent(inout)			    :: node

			call post_dof_op(node%data_pers%saturation, node%data_pers%p, node%data_temp%volume, node%data_temp%mat_diagonal)
		end subroutine
		!*******************************
		!Volume and DoF operators
		!*******************************

		elemental subroutine pre_dof_op(saturation, volume, p, mat_diagonal, d, A_d, rhs)
			real (kind = GRID_SR), intent(out)		:: saturation
			real (kind = GRID_SR), intent(out)		:: volume
			real (kind = GRID_SR), intent(out)		:: p
			real (kind = GRID_SR), intent(out)		:: mat_diagonal
			real (kind = GRID_SR), intent(out)		:: d
			real (kind = GRID_SR), intent(out)		:: A_d
			real (kind = GRID_SR), intent(out)		:: rhs

			saturation = 0.0_GRID_SR
			volume = 0.0_GRID_SR
			p = 0.0_GRID_SR
			mat_diagonal = 0.0_GRID_SR
			d = 0.0_GRID_SR
			A_d = 0.0_GRID_SR
			rhs = 0.0_SR
		end subroutine

		elemental subroutine post_dof_op(saturation, p, volume, p_volume)
 			real (kind = GRID_SR), intent(inout)	:: saturation
 			real (kind = GRID_SR), intent(inout)	:: p
			real (kind = GRID_SR), intent(in)		:: volume
			real (kind = GRID_SR), intent(in)	    :: p_volume

            !assert_pure(volume > 0.0_SR)

            saturation = saturation / volume
            p = p / p_volume

            !assert_pure(saturation <= 1.0_SR)
		end subroutine
	END MODULE
#endif
