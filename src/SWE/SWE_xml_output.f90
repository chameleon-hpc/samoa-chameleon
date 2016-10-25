! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_xml_output
		use LIB_VTK_IO

		use SFC_edge_traversal

		use Samoa_swe
		use SWE_euler_timestep

#       if defined(_SWE_PATCH)
            use SWE_PATCH
#       endif
		
		implicit none

		!> Output point data
		type t_output_point_data
			type(t_state)											:: Q
			real (kind = GRID_SR), dimension(2)						:: coords		!< position
		end type

		!> Output cell data
		type t_output_cell_data
			type(t_state)											:: Q
			integer (kind = GRID_SI)								:: rank
			integer (kind = GRID_SI)								:: section_index
			integer (kind = BYTE)									:: depth
			integer (kind = BYTE)									:: refinement
                        integer (kind=BYTE) :: troubled

            integer (kind = BYTE)                                   :: i_plotter_type
#           if defined(_SWE_PATCH)          
                integer (kind = BYTE)                               :: id_in_patch
#           endif
		end type

        type num_traversal_data
            type(t_output_point_data), allocatable		            :: point_data(:)
            type(t_output_cell_data), allocatable			        :: cell_data(:)
            character(len=256)							            :: s_file_stamp

            integer (kind = GRID_SI)								:: i_output_iteration = 0
            integer (kind = GRID_SI)								:: i_point_data_index
            integer (kind = GRID_SI)								:: i_cell_data_index
        end type

		integer, parameter											:: i_element_order = 0

		type(t_gv_Q)												:: gv_Q

#		define _GT_NAME								t_swe_xml_output_traversal

#		define _GT_EDGES

#		define _GT_PRE_TRAVERSAL_OP					pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP				post_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP				pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP				post_traversal_grid_op

#		define _GT_ELEMENT_OP						element_op

!#		define _GT_CELL_TO_EDGE_OP				    cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

        subroutine pre_traversal_grid_op(traversal, grid)
			type(t_swe_xml_output_traversal), intent(inout)				:: traversal
			type(t_grid), intent(inout)							        :: grid

            if (rank_MPI == 0) then
                _log_write(1, '(A, I0)') " SWE: output step: ", traversal%i_output_iteration
            end if


            ! if(traversal%i_output_iteration .eq. 3) then
            !    stop
            ! end if
            call scatter(traversal%s_file_stamp, traversal%sections%s_file_stamp)
            call scatter(traversal%i_output_iteration, traversal%sections%i_output_iteration)

		end subroutine

        subroutine post_traversal_grid_op(traversal, grid)
			type(t_swe_xml_output_traversal), intent(inout)				:: traversal
			type(t_grid), intent(inout)							        :: grid

            character (len = 256)							:: s_file_name
            integer                                         :: i_error
			integer(4)										:: i_rank, i_section, e_io
			logical                                         :: l_exists
            type(t_vtk_writer)                              :: vtk
#           if defined(_MPI)
                call mpi_barrier(MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#           endif

#           if defined(_QUAD_PRECISION)
#               warning VTK output does not work for quad precision
#           else
                if (rank_MPI == 0) then
                    write (s_file_name, "(A, A, I0, A)") TRIM(traversal%s_file_stamp), "_", traversal%i_output_iteration, ".pvtu"
                    e_io = vtk%VTK_INI_XML('ascii', s_file_name, 'PUnstructuredGrid')

                    e_io = vtk%VTK_DAT_XML('pfield', 'OPEN')
                        e_io = vtk%VTK_VAR_XML('time', 1.0_GRID_SR, 1)
                    e_io = vtk%VTK_DAT_XML('pfield', 'CLOSE')

                    e_io = vtk%VTK_DAT_XML('pnode', 'OPEN')
                        if (i_element_order > 0) then
                            e_io = vtk%VTK_VAR_XML('water height', 1.0_GRID_SR, 1)
                            e_io = vtk%VTK_VAR_XML('bathymetry', 1.0_GRID_SR, 1)
                            e_io = vtk%VTK_VAR_XML('momentum', 1.0_GRID_SR, 3)
                        end if
                    e_io = vtk%VTK_DAT_XML('pnode', 'CLOSE')

                    e_io = vtk%VTK_DAT_XML('pcell', 'OPEN')
                        if (i_element_order == 0) then
                            e_io = vtk%VTK_VAR_XML('water height', 1.0_GRID_SR, 1)
                            e_io = vtk%VTK_VAR_XML('bathymetry', 1.0_GRID_SR, 1)
                            e_io = vtk%VTK_VAR_XML('momentum', 1.0_GRID_SR, 3)
                        end if

                        e_io = vtk%VTK_VAR_XML('rank', 1_GRID_SI, 1)
                        e_io = vtk%VTK_VAR_XML('section index', 1_GRID_SI, 1)
                        e_io = vtk%VTK_VAR_XML('depth', 1_1, 1)
                        e_io = vtk%VTK_VAR_XML('refinement flag', 1_1, 1)
                        e_io = vtk%VTK_VAR_XML('i_plotter_type', 1_1, 1)
                        e_io = vtk%VTK_VAR_XML('troubled', 1_1, 1)
#                       if defined(_SWE_PATCH)
                            e_io = vtk%VTK_VAR_XML('id_in_patch', 1_1, 1)
#                       endif
                    e_io = vtk%VTK_DAT_XML('pcell', 'CLOSE')

                    e_io = vtk%VTK_GEO_XML(1.0_GRID_SR)

                    do i_rank = 0, size_MPI
                        do i_section = 1, huge(1)
                            write (s_file_name, "(A, A, I0, A, I0, A, I0, A)") trim(traversal%s_file_stamp), "_", traversal%i_output_iteration, "_r", i_rank, "_s", i_section, ".vtu"
                            inquire(file = s_file_name, exist = l_exists)

                            if (l_exists) then
                                write(s_file_name, "(A)") trim(s_file_name(scan(s_file_name, "/\&
" , .true.) + 1 : len(s_file_name)))
                                e_io = vtk%VTK_GEO_XML(s_file_name)
                            else
                                exit
                            end if
                        end do
                    end do

                    e_io = vtk%VTK_END_XML()
                end if
#           endif

            traversal%i_output_iteration = traversal%i_output_iteration + 1
        end subroutine

		subroutine pre_traversal_op(traversal, section)
			type(t_swe_xml_output_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section

			type(t_section_info)                                           :: grid_info
			integer (kind = GRID_SI)	:: i_error, i_cells, i_points

            grid_info = section%get_info()
#           if defined (_SWE_PATCH)
                i_cells = grid_info%i_cells * _SWE_PATCH_ORDER_SQUARE
#           else
                i_cells = grid_info%i_cells
#           endif

			if (i_element_order > 1) then
				i_points = 6 * i_cells
			else
				i_points = 3 * i_cells
			end if

			allocate(traversal%cell_data(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(traversal%point_data(i_points), stat = i_error); assert_eq(i_error, 0)

			traversal%i_cell_data_index = 1
			traversal%i_point_data_index = 1
		end subroutine

		subroutine post_traversal_op(traversal, section)
			type(t_swe_xml_output_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section

			integer (kind = GRID_SI), dimension(:), allocatable			:: i_offsets
			integer (1), dimension(:), allocatable						:: i_types
			integer (kind = GRID_SI), dimension(:), allocatable			:: i_connectivity
			real (kind = GRID_SR), dimension(:), allocatable			:: r_empty
            type(t_vtk_writer)                                          :: vtk

			type(t_section_info)                                        :: grid_info
			integer (kind = GRID_SI)									:: i_error, i_cells, i_points, i
			integer(4)													:: e_io
            character (len = 256)							            :: s_file_name

            grid_info = section%get_info()
#           if defined (_SWE_PATCH)
                i_cells = grid_info%i_cells * _SWE_PATCH_ORDER_SQUARE
#           else
                i_cells = grid_info%i_cells
#           endif
			if (i_element_order > 1) then
				i_points = 6 * i_cells
				allocate(i_connectivity(6 * i_cells), stat = i_error); assert_eq(i_error, 0)
			else
				i_points = 3 * i_cells
				allocate(i_connectivity(3 * i_cells), stat = i_error); assert_eq(i_error, 0)
			end if

			allocate(i_offsets(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(i_types(i_cells), stat = i_error); assert_eq(i_error, 0)
			allocate(r_empty(max(i_cells, i_points)), stat = i_error); assert_eq(i_error, 0)

			r_empty = 0.0_GRID_SR

			if (i_element_order > 1) then
				i_types = 22_1

				forall (i = 1 : i_cells)
					i_offsets(i) = 6_GRID_SI * i
					i_connectivity(6 * i - 5 : 6 * i) = [6 * i - 6, 6 * i - 5, 6 * i - 4, 6 * i - 3, 6 * i - 2, 6 * i - 1]
				end forall
			else
				i_types = 5_1

				forall (i = 1 : i_cells)
					i_offsets(i) = 3_GRID_SI * i
					i_connectivity(3 * i - 2 : 3 * i) = [3 * i - 3, 3 * i - 2, 3 * i - 1]
				end forall
			end if

			write (s_file_name, "(A, A, I0, A, I0, A, I0, A)") trim(traversal%s_file_stamp), "_", traversal%i_output_iteration, "_r", rank_MPI, "_s", section%index, ".vtu"

#           if defined(_QUAD_PRECISION)
#               warning VTK output does not work for quad precision
#           else
                e_io = vtk%VTK_INI_XML('binary', s_file_name, 'UnstructuredGrid')
                    e_io = vtk%VTK_DAT_XML('field', 'OPEN')
                        e_io = vtk%VTK_VAR_XML(1, 'time', [section%r_time])
                    e_io = vtk%VTK_DAT_XML('field', 'CLOSE')

                    e_io = vtk%VTK_GEO_XML(i_points, i_cells, traversal%point_data%coords(1), traversal%point_data%coords(2), r_empty(1:i_points))

                    e_io = vtk%VTK_CON_XML(i_cells, i_connectivity, i_offsets, i_types)

                    e_io = vtk%VTK_DAT_XML('node', 'OPEN')
                        if (i_element_order > 0) then
                            e_io = vtk%VTK_VAR_XML(i_points, 'water height', traversal%point_data%Q%h)
                            e_io = vtk%VTK_VAR_XML(i_points, 'bathymetry', traversal%point_data%Q%b)
                            e_io = vtk%VTK_VAR_XML(i_points, 'momentum',  traversal%point_data%Q%p(1), traversal%point_data%Q%p(2), r_empty(1:i_points))
                        end if
                    e_io = vtk%VTK_DAT_XML('node', 'CLOSE')

                    e_io = vtk%VTK_DAT_XML('cell', 'OPEN')
                        if (i_element_order == 0) then
                            e_io = vtk%VTK_VAR_XML(i_cells, 'water height', traversal%cell_data%Q%h)
                            e_io = vtk%VTK_VAR_XML(i_cells, 'bathymetry', traversal%cell_data%Q%b)
                            e_io = vtk%VTK_VAR_XML(i_cells, 'momentum', traversal%cell_data%Q%p(1), traversal%cell_data%Q%p(2), r_empty(1:i_cells))
                        end if

                        e_io = vtk%VTK_VAR_XML(i_cells, 'rank', traversal%cell_data%rank)
                        e_io = vtk%VTK_VAR_XML(i_cells, 'section index', traversal%cell_data%section_index)
                        e_io = vtk%VTK_VAR_XML(i_cells, 'depth', traversal%cell_data%depth)
                        e_io = vtk%VTK_VAR_XML(i_cells, 'refinement flag', traversal%cell_data%refinement)
                        e_io = vtk%VTK_VAR_XML(i_cells, 'i_plotter_type', traversal%cell_data%i_plotter_type)
                        e_io = vtk%VTK_VAR_XML(i_cells, 'troubled', traversal%cell_data%troubled)
#                       if defined(_SWE_PATCH)                        
                            e_io = vtk%VTK_VAR_XML(i_cells, 'id_in_patch', traversal%cell_data%id_in_patch)
#                       endif
                    e_io = vtk%VTK_DAT_XML('cell', 'CLOSE')

                    e_io = vtk%VTK_GEO_XML()
                e_io = vtk%VTK_END_XML()
#           endif

			deallocate(i_offsets, stat = i_error); assert_eq(i_error, 0)
			deallocate(i_types, stat = i_error); assert_eq(i_error, 0)
			deallocate(i_connectivity, stat = i_error); assert_eq(i_error, 0)
			deallocate(r_empty, stat = i_error); assert_eq(i_error, 0)

			deallocate(traversal%cell_data, stat = i_error); assert_eq(i_error, 0)
			deallocate(traversal%point_data, stat = i_error); assert_eq(i_error, 0)

			traversal%i_output_iteration = traversal%i_output_iteration + 1
		end subroutine


		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
			type(t_swe_xml_output_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)					:: element

			!local variables

			integer (kind = GRID_SI)							:: i
			real (kind = GRID_SR), parameter, dimension(2, 6)	:: r_test_points_forward = reshape([1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.5, 0.0, 0.5, 0.5 ], [2, 6 ])
			real (kind = GRID_SR), parameter, dimension(2, 6)	:: r_test_points_backward = reshape([0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.5, 0.5 ], [2, 6 ])
			real (kind = GRID_SR), parameter, dimension(2)		:: r_test_point0 = [1.0_GRID_SR/3.0_GRID_SR, 1.0_GRID_SR/3.0_GRID_SR]
#			if defined(_SWE_PATCH)
				type(t_state), dimension(_SWE_PATCH_ORDER_SQUARE):: Q
				integer											:: j, row, col, cell_id
#if defined(_SWE_DG)    
    type(num_cell_data_pers) :: data_temp
    if(element%cell%data_pers%troubled .le. 0) then
       call element%cell%data_pers%convert_dg_to_fv()
    end if
#endif
               
                row=1
                col=1

                do j=1,_SWE_PATCH_ORDER * _SWE_PATCH_ORDER
                    traversal%cell_data(traversal%i_cell_data_index)%rank = rank_MPI
                    traversal%cell_data(traversal%i_cell_data_index)%section_index = section%index
                    traversal%cell_data(traversal%i_cell_data_index)%depth = element%cell%geometry%i_depth
                    traversal%cell_data(traversal%i_cell_data_index)%i_plotter_type = element%cell%geometry%i_plotter_type
                    traversal%cell_data(traversal%i_cell_data_index)%refinement = element%cell%geometry%refinement
                    traversal%cell_data(traversal%i_cell_data_index)%troubled = element%cell%data_pers%troubled

                    ! traversal%cell_data(traversal%i_cell_data_index)%troubled = 0
                    ! do i=1,size(element%edges)
                    !    if(element%edges(i)%ptr%data_pers%troubled)then
                    !       traversal%cell_data(traversal%i_cell_data_index)%troubled = traversal%cell_data(traversal%i_cell_data_index)%troubled +1
                    !    end if
                    ! end do

                    do i=1,3
                        traversal%point_data(traversal%i_point_data_index + i - 1)%coords = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, SWE_PATCH_geometry%coords(:,i, j)) + cfg%offset
                    end do

                    traversal%i_point_data_index = traversal%i_point_data_index + 3
                    
                    ! if orientation is backwards, the plotter uses a transformation that mirrors the cell...
                    ! this simple change solves the problem :)
                    if (element%cell%geometry%i_plotter_type > 0) then 
                        cell_id = j
                    else
                        cell_id = (row-1)*(row-1) + 2 * row - col
                    end if
                    traversal%cell_data(traversal%i_cell_data_index)%id_in_patch = cell_id

                    traversal%cell_data(traversal%i_cell_data_index)%Q%h = element%cell%data_pers%H(cell_id)
                    traversal%cell_data(traversal%i_cell_data_index)%Q%b = element%cell%data_pers%B(cell_id)
                    traversal%cell_data(traversal%i_cell_data_index)%Q%p(1) = element%cell%data_pers%HU(cell_id)
                    traversal%cell_data(traversal%i_cell_data_index)%Q%p(2) = element%cell%data_pers%HV(cell_id)

                    
                    ! prepare for next cell
                    traversal%i_cell_data_index = traversal%i_cell_data_index + 1
                    col = col + 1
                    if (col == 2*row) then
                        col = 1
                        row = row + 1
                    end if
                end do

#           else
			type(t_state), dimension(_SWE_CELL_SIZE)			:: Q
			type(t_state), dimension(6)							:: Q_test

			call gv_Q%read(element, Q)

                        traversal%cell_data(traversal%i_cell_data_index)%rank = rank_MPI
                        traversal%cell_data(traversal%i_cell_data_index)%section_index = section%index
			traversal%cell_data(traversal%i_cell_data_index)%depth = element%cell%geometry%i_depth
                        traversal%cell_data(traversal%i_cell_data_index)%i_plotter_type = element%cell%geometry%i_plotter_type
                        traversal%cell_data(traversal%i_cell_data_index)%refinement = element%cell%geometry%refinement
                        traversal%cell_data(traversal%i_cell_data_index)%troubled =  element%cell%data_pers%troubled

			select case (i_element_order)
				case (2)
                    if (element%transform_data%plotter_data%orientation > 0) then
                        forall (i = 1 : 6)
                            traversal%point_data(traversal%i_point_data_index + i - 1)%coords = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, r_test_points_forward(:, i)) + cfg%offset
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%h = t_basis_Q_eval(r_test_points_forward(:, i), Q%h)
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%b = t_basis_Q_eval(r_test_points_forward(:, i), Q%b)
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(1) = t_basis_Q_eval(r_test_points_forward(:, i), Q%p(1))
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(2) = t_basis_Q_eval(r_test_points_forward(:, i), Q%p(2))
                        end forall
                    else
                        forall (i = 1 : 6)
                            traversal%point_data(traversal%i_point_data_index + i - 1)%coords = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, r_test_points_backward(:, i)) + cfg%offset
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%h = t_basis_Q_eval(r_test_points_backward(:, i), Q%h)
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%b = t_basis_Q_eval(r_test_points_backward(:, i), Q%b)
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(1) = t_basis_Q_eval(r_test_points_backward(:, i), Q%p(1))
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(2) = t_basis_Q_eval(r_test_points_backward(:, i), Q%p(2))
                        end forall
                    end if

					traversal%i_point_data_index = traversal%i_point_data_index + 6
				case (1)
                    if (element%transform_data%plotter_data%orientation > 0) then
                        forall (i = 1 : 3)
                            traversal%point_data(traversal%i_point_data_index + i - 1)%coords = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, r_test_points_forward(:, i)) + cfg%offset
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%h = t_basis_Q_eval(r_test_points_forward(:, i), Q%h)
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%b = t_basis_Q_eval(r_test_points_forward(:, i), Q%b)
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(1) = t_basis_Q_eval(r_test_points_forward(:, i), Q%p(1))
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(2) = t_basis_Q_eval(r_test_points_forward(:, i), Q%p(2))
                        end forall
                    else
                        forall (i = 1 : 3)
                            traversal%point_data(traversal%i_point_data_index + i - 1)%coords = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, r_test_points_backward(:, i)) + cfg%offset
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%h = t_basis_Q_eval(r_test_points_backward(:, i), Q%h)
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%b = t_basis_Q_eval(r_test_points_backward(:, i), Q%b)
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(1) = t_basis_Q_eval(r_test_points_backward(:, i), Q%p(1))
                            traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(2) = t_basis_Q_eval(r_test_points_backward(:, i), Q%p(2))
                        end forall
                    end if

					traversal%i_point_data_index = traversal%i_point_data_index + 3
				case (0)
                    if (element%transform_data%plotter_data%orientation > 0) then
                        forall (i = 1 : 3)
                            traversal%point_data(traversal%i_point_data_index + i - 1)%coords = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, r_test_points_forward(:, i)) + cfg%offset
                        end forall
                    else
                        forall (i = 1 : 3)
                            traversal%point_data(traversal%i_point_data_index + i - 1)%coords = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, r_test_points_backward(:, i)) + cfg%offset
                        end forall
                    end if

					traversal%i_point_data_index = traversal%i_point_data_index + 3

					traversal%cell_data(traversal%i_cell_data_index)%Q%h = t_basis_Q_eval(r_test_point0, Q%h)
					traversal%cell_data(traversal%i_cell_data_index)%Q%b = t_basis_Q_eval(r_test_point0, Q%b)
					traversal%cell_data(traversal%i_cell_data_index)%Q%p(1) = t_basis_Q_eval(r_test_point0, Q%p(1))
					traversal%cell_data(traversal%i_cell_data_index)%Q%p(2) = t_basis_Q_eval(r_test_point0, Q%p(2))
			end select

			traversal%i_cell_data_index = traversal%i_cell_data_index + 1
#           endif

                ! set water height of dry cells to zero

                where (traversal%cell_data(:)%Q%h < traversal%cell_data(:)%Q%b + cfg%dry_tolerance )
                    traversal%cell_data(:)%Q%h = min(0.0_GRID_SR, traversal%cell_data(:)%Q%b)
                end where
                
		end subroutine
	END MODULE
#endif
