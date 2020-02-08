! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE

#include "Compilation_control.f90"

#if defined(_SWE)
MODULE SWE_xml_point_output
  use LIB_VTK_IO

  use SFC_edge_traversal
  use SWE_DG_limiter

  use Samoa_swe

#       if defined(_SWE_PATCH)
  use SWE_PATCH
#       endif
#if defined(_SWE_DG)
  use SWE_dg_matrices
  use SWE_data_types
#endif            

  implicit none

  !> Output point data
  type t_output_point_data
     type(t_state)											:: Q
     real (kind = GRID_SR), dimension(2)						:: coords		!< position
     integer (kind = GRID_SI)								:: rank
     integer (kind = GRID_SI)								:: section_index
     integer (kind = BYTE)									:: depth
     integer (kind = BYTE)									:: refinement

#if defined(_SWE_DG)
     integer (kind=BYTE) :: troubled
#endif

     integer (kind = BYTE)                                   :: i_plotter_type
#           if defined(_SWE_PATCH)          
     integer (kind = BYTE)                               :: id_in_patch
#           endif
  end type t_output_point_data

  !> Output cell data
  type t_output_cell_data
     integer (kind = GRID_SI), dimension(3) :: connectivity 
     integer (kind=BYTE) :: troubled
  end type t_output_cell_data

  type num_traversal_data
     type(t_output_point_data), allocatable		            :: point_data(:)
     type(t_output_cell_data), allocatable			        :: cell_data(:)
     ! We need this because the number of cells depends on used solver for each element)
     integer (kind=BYTE), allocatable :: troubled(:)
     character(len=256)							            :: s_file_stamp

     integer (kind = GRID_SI)								:: i_output_iteration = 0
     integer (kind = GRID_SI)								:: i_element_data_index
     integer (kind = GRID_SI)								:: i_point_data_index
     integer (kind = GRID_SI)								:: i_cell_data_index
  end type num_traversal_data

  integer, parameter											:: i_element_order = 0

  type(t_gv_Q)												:: gv_Q

#		define _GT_NAME								t_swe_xml_point_output_traversal

#		define _GT_EDGES

#		define _GT_PRE_TRAVERSAL_OP					pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP				post_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP				pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP				post_traversal_grid_op

#		define _GT_ELEMENT_OP						element_op

  !#		define _GT_CELL_TO_EDGE_OP				    cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

  subroutine pre_traversal_grid_op(traversal, grid)
    type(t_swe_xml_point_output_traversal), intent(inout)				:: traversal
    type(t_grid), intent(inout)							        :: grid

    if (rank_MPI == 0) then
       _log_write(1, '(A, I0)') " SWE: output step: ", traversal%i_output_iteration
    end if


    call scatter(traversal%s_file_stamp, traversal%sections%s_file_stamp)
    call scatter(traversal%i_output_iteration, traversal%sections%i_output_iteration)

  end subroutine pre_traversal_grid_op

  subroutine post_traversal_grid_op(traversal, grid)
    type(t_swe_xml_point_output_traversal), intent(inout)				:: traversal
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
       e_io = vtk%VTK_VAR_XML('water height', 1.0_GRID_SR, 1)
       e_io = vtk%VTK_VAR_XML('bathymetry', 1.0_GRID_SR, 1)
       e_io = vtk%VTK_VAR_XML('momentum', 1.0_GRID_SR, 3)

       e_io = vtk%VTK_VAR_XML('rank', 1_GRID_SI, 1)
       e_io = vtk%VTK_VAR_XML('section index', 1_GRID_SI, 1)
       e_io = vtk%VTK_VAR_XML('depth', 1_1, 1)
       e_io = vtk%VTK_VAR_XML('refinement flag', 1_1, 1)
       e_io = vtk%VTK_VAR_XML('i_plotter_type', 1_1, 1)
#if defined(_SWE_DG)
       e_io = vtk%VTK_VAR_XML('troubled', 1_1, 1)
#endif
#                       if defined(_SWE_PATCH)
       e_io = vtk%VTK_VAR_XML('id_in_patch', 1_1, 1)
#                       endif

       e_io = vtk%VTK_DAT_XML('pnode', 'CLOSE')

       e_io = vtk%VTK_DAT_XML('pcell', 'OPEN')
       ! No cell data
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
  end subroutine post_traversal_grid_op

  subroutine pre_traversal_op(traversal, section)
    type(t_swe_xml_point_output_traversal), intent(inout)				:: traversal
    type(t_grid_section), intent(inout)							:: section

    type(t_section_info)                                           :: grid_info
    integer (kind = GRID_SI)	:: i_error, i_cells, i_points

    grid_info = section%get_info()
    ! Upper bound of needed cells and points
    i_cells = grid_info%i_cells * max(_SWE_DG_ORDER * _SWE_DG_ORDER, _SWE_PATCH_ORDER_SQUARE)
    i_points = i_cells * 3

    allocate(traversal%troubled(grid_info%i_cells), stat = i_error); assert_eq(i_error, 0)
    allocate(traversal%cell_data(i_cells), stat = i_error); assert_eq(i_error, 0)
    allocate(traversal%point_data(i_points), stat = i_error); assert_eq(i_error, 0)

    traversal%i_element_data_index = 1
    traversal%i_cell_data_index = 1
    traversal%i_point_data_index = 1
  end subroutine pre_traversal_op

  subroutine post_traversal_op(traversal, section)
    type(t_swe_xml_point_output_traversal), intent(inout)				:: traversal
    type(t_grid_section), intent(inout)							:: section

    integer (kind = GRID_SI), dimension(:), allocatable			:: i_offsets
    integer (1), dimension(:), allocatable						:: i_types
    integer (kind = GRID_SI), dimension(:), allocatable			:: i_connectivity
    real (kind = GRID_SR), dimension(:), allocatable			:: r_empty
    logical, dimension(:), allocatable :: l_point_mask, l_cell_mask
    type(t_vtk_writer)                                          :: vtk

    type(t_section_info)                                        :: grid_info
    integer (kind = GRID_SI)									:: i_error, i_cells, i_points, i, j, i_real_points, i_real_cells
    integer(4)													:: e_io
    character (len = 256)							            :: s_file_name
    integer (kind = GRID_SI), parameter :: i_cells_per_vf = _SWE_PATCH_ORDER_SQUARE
    integer (kind = GRID_SI), parameter :: i_cells_per_dg = _SWE_DG_ORDER * _SWE_DG_ORDER 

    ! These are all cells that are saved. Some of them are invalid.
    grid_info = section%get_info()
    i_cells = grid_info%i_cells * max(i_cells_per_vf, i_cells_per_dg)
    i_points = 3 * i_cells

    ! Because we have a hybrid method, each element contains a varying number of points/cells.
    i_real_cells = 0
    i_real_points = 0
    allocate(l_point_mask(i_points), stat = i_error); assert_eq(i_error, 0)
    allocate(l_cell_mask(i_cells), stat = i_error); assert_eq(i_error, 0)

    ! All points/cells that have a mask value of true contain valid data.
    ! This makes writing only correct cells easy.
    l_point_mask(:) = .false.
    l_cell_mask(:) = .false.

    ! Iterate over all elements to count stuff correctly.
    do i=1, grid_info%i_cells
       if (.not.isCoast(traversal%troubled(i))) then
       !if (traversal%troubled(i).le.0) then       
          l_point_mask(i_real_points + 1 : i_real_points + _SWE_DG_DOFS) = .true.
          l_cell_mask(i_real_cells + 1: i_real_cells + i_cells_per_dg) = .true.
          i_real_cells = i_real_cells + i_cells_per_dg
          i_real_points = i_real_points + _SWE_DG_DOFS
       else
          l_point_mask(i_real_points + 1 : i_real_points + i_cells_per_vf * 3) = .true.
          l_cell_mask(i_real_cells + 1: i_real_cells + i_cells_per_vf) = .true.
          i_real_cells = i_real_cells + i_cells_per_vf
          i_real_points = i_real_points + i_cells_per_vf * 3
       end if
    end do

    ! We now know how many cells/points we are going to write.
    allocate(i_connectivity(i_real_cells * 3), stat = i_error); assert_eq(i_error, 0)
    allocate(i_offsets(i_real_cells), stat = i_error); assert_eq(i_error, 0)
    allocate(i_types(i_real_cells), stat = i_error); assert_eq(i_error, 0)
    allocate(r_empty(max(i_real_cells, i_real_points)), stat = i_error); assert_eq(i_error, 0)
    r_empty = 0.0_GRID_SR
    i_types = 5_1 ! Use triangle output type for every cells.

    ! Find out which cells consist of which points. This allows us to avoid writing duplicate points.
    ! Reuse as counters, they have the same value after this loop again.
    i_real_cells = 0
    do i=1,i_cells
       if (l_cell_mask(i)) then
          i_connectivity(i_real_cells * 3 + 1 : i_real_cells * 3 + 3) = traversal%cell_data(i)%connectivity(:)
          i_real_cells = i_real_cells + 1
       end if
    end do
    ! Offset is easier, all cells are triangles and are defined by exactly three points.
    i_offsets(1 : i_real_cells) = (/ ( (i) * 3, i=1,i_real_cells) /)

    write (s_file_name, "(A, A, I0, A, I0, A, I0, A)") trim(traversal%s_file_stamp), "_", traversal%i_output_iteration, "_r", rank_MPI, "_s", section%index, ".vtu"

#           if defined(_QUAD_PRECISION)
#               warning VTK output does not work for quad precision
#           else
    e_io = vtk%VTK_INI_XML('binary', s_file_name, 'UnstructuredGrid')
    e_io = vtk%VTK_DAT_XML('field', 'OPEN')
    e_io = vtk%VTK_VAR_XML(1, 'time', [section%r_time])
    e_io = vtk%VTK_DAT_XML('field', 'CLOSE')

    e_io = vtk%VTK_GEO_XML(i_real_points, i_real_cells, pack(traversal%point_data%coords(1), l_point_mask), pack(traversal%point_data%coords(2), l_point_mask), r_empty(1:i_real_points))

    e_io = vtk%VTK_CON_XML(i_real_cells, i_connectivity, i_offsets, i_types)

    e_io = vtk%VTK_DAT_XML('node', 'OPEN')
    e_io = vtk%VTK_VAR_XML(i_real_points, 'water height', pack(traversal%point_data%Q%h, l_point_mask))
    e_io = vtk%VTK_VAR_XML(i_real_points, 'bathymetry', pack(traversal%point_data%Q%b, l_point_mask))
    e_io = vtk%VTK_VAR_XML(i_real_points, 'momentum',  pack(traversal%point_data%Q%p(1), l_point_mask), pack(traversal%point_data%Q%p(2), l_point_mask), r_empty(1:i_real_points))

    e_io = vtk%VTK_VAR_XML(i_real_points, 'rank', pack(traversal%point_data%rank, l_point_mask))
    e_io = vtk%VTK_VAR_XML(i_real_points, 'section index', pack(traversal%point_data%section_index, l_point_mask))
    e_io = vtk%VTK_VAR_XML(i_real_points, 'depth', pack(traversal%point_data%depth, l_point_mask))
    e_io = vtk%VTK_VAR_XML(i_real_points, 'refinement flag', pack(traversal%point_data%refinement, l_point_mask))
    e_io = vtk%VTK_VAR_XML(i_real_points, 'i_plotter_type', pack(traversal%point_data%i_plotter_type, l_point_mask))
#if defined(_SWE_DG)
    e_io = vtk%VTK_VAR_XML(i_real_points, 'troubled', pack(traversal%point_data%troubled, l_point_mask))
#endif
#                       if defined(_SWE_PATCH)                        
    e_io = vtk%VTK_VAR_XML(i_real_points, 'id_in_patch', pack(traversal%point_data%id_in_patch, l_point_mask))
#                       endif

    e_io = vtk%VTK_DAT_XML('node', 'CLOSE')

    e_io = vtk%VTK_DAT_XML('cell', 'OPEN')
    ! No cell data
    e_io = vtk%VTK_DAT_XML('cell', 'CLOSE')

    e_io = vtk%VTK_GEO_XML()
    e_io = vtk%VTK_END_XML()
#           endif

    deallocate(i_offsets, stat = i_error); assert_eq(i_error, 0)
    deallocate(i_types, stat = i_error); assert_eq(i_error, 0)
    deallocate(i_connectivity, stat = i_error); assert_eq(i_error, 0)
    deallocate(r_empty, stat = i_error); assert_eq(i_error, 0)
    deallocate(l_point_mask, stat = i_error); assert_eq(i_error, 0)
    deallocate(l_cell_mask, stat = i_error); assert_eq(i_error, 0)

    deallocate(traversal%cell_data, stat = i_error); assert_eq(i_error, 0)
    deallocate(traversal%point_data, stat = i_error); assert_eq(i_error, 0)
    deallocate(traversal%troubled, stat = i_error); assert_eq(i_error, 0)    

    traversal%i_output_iteration = traversal%i_output_iteration + 1
  end subroutine post_traversal_op


  !******************
  !Geometry operators
  !******************

  subroutine element_op(traversal, section, element)
    type(t_swe_xml_point_output_traversal), intent(inout)				:: traversal
    type(t_grid_section), intent(inout)							:: section
    type(t_element_base), intent(inout)					:: element

    !local variables

    ! Only defined for orders 1..4
    ! We need a different number of cells and points for different orders.
    ! The array mirrored_coords maps non-mirrored coordinates to mirrored coordinates.
    integer, parameter :: num_cells = _SWE_DG_ORDER * _SWE_DG_ORDER
# if (_SWE_DG_ORDER == 1)
!    real (kind = GRID_SR), parameter, dimension(2, _SWE_DG_DOFS) :: coords = reshape([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], [2, _SWE_DG_DOFS ])
    integer (kind = GRID_SI), parameter, dimension(_SWE_DG_DOFS) :: mirrored_coords = [1, 3, 2]
# elif (_SWE_DG_ORDER == 2)
!    real (kind = GRID_SR), parameter, dimension(2, _SWE_DG_DOFS) :: coords = reshape([0.0, 0.0, 0.5, 0.0, 1.0, 0.0, &
!         0.0, 0.5, 0.5, 0.5, 0.0, 1.0], [2, _SWE_DG_DOFS ])
    integer (kind = GRID_SI), parameter, dimension(_SWE_DG_DOFS) :: mirrored_coords = [1, 4, 6, 2, 5, 3]
# elif (_SWE_DG_ORDER == 3)
    ! real (kind = GRID_SR), parameter, dimension(2, _SWE_DG_DOFS) :: coords = reshape([0.0, 0.0, 1.0/3.0, 0.0, 2.0/3.0, 0.0, 1.0, 0.0, &
    !      0.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 2.0/3.0, 1.0/3.0, &
    !      0.0, 2.0/3.0, 1.0/3.0, 2.0/3.0, &
    !      0.0, 1.0], [2, _SWE_DG_DOFS ])
    integer (kind = GRID_SI), parameter, dimension(_SWE_DG_DOFS) :: mirrored_coords = [1, 5, 8, 10, 2, 6, 9, 3, 7, 4]
# elif (_SWE_DG_ORDER == 4)
    ! real (kind = GRID_SR), parameter, dimension(2, 15)	:: coords = reshape([0.0, 0.0, 0.25, 0.0, 0.5, 0.0, 0.75, 0.0, 1.0, 0.0, &
    !      0.0, 0.25, 0.25, 0.25, 0.5, 0.25, 0.75, 0.25, &
    !      0.0, 0.5, 0.25, 0.5, 0.5, 0.5, &
    !      0.0, 0.75, 0.25, 0.75, 0.0, 1.0], [2, _SWE_DG_DOFS])
    integer (kind = GRID_SI), parameter, dimension(_SWE_DG_DOFS) :: mirrored_coords = [ 1, 6, 10, 13, 15, 2, 7, 11, 14, 3, 8, 12, 4, 9, 5]
#endif
    real (kind = GRID_SR), parameter, dimension(2, _SWE_DG_DOFS)	:: coords = nodes

    ! Only implemented for patches and _SWE_DG!
# if defined(_SWE_PATCH)
    type(t_state), dimension(_SWE_PATCH_ORDER_SQUARE):: Q
    integer	:: i, j, row, col, cell_id

    ! We need this to count the number of valid cells in the post-traversal op.
    traversal%troubled(traversal%i_element_data_index) = element%cell%data_pers%troubled
    if (.not.isCoast(element%cell%data_pers%troubled)) then
!    if (element%cell%data_pers%troubled.le.0 .or. element%cell%data_pers%troubled.ge.5) then
!    if (element%cell%data_pers%troubled.le.0) then    
       ! Calulated by DG.

       traversal%cell_data(traversal%i_cell_data_index : traversal%i_cell_data_index + num_cells - 1)%troubled = element%cell%data_pers%troubled
       ! The number of cells per element depends on the order of the dg method. This is needed for a pleasing output!
       ! Beware of off-by-one errors: the connectivity data needs the offset of -2 because vtk counts cells from indices and we count both triangle vertices and points
       ! starting with 1! 
# if (_SWE_DG_ORDER == 1)
          traversal%cell_data(traversal%i_cell_data_index +  0)%connectivity(:) = [ 1, 2, 3] + traversal%i_point_data_index - 2
# elif (_SWE_DG_ORDER == 2)
          traversal%cell_data(traversal%i_cell_data_index +  0)%connectivity(:) = [ 1, 2, 4] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  1)%connectivity(:) = [ 2, 3, 5] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  2)%connectivity(:) = [ 4, 5, 2] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  3)%connectivity(:) = [ 4, 5, 6] + traversal%i_point_data_index - 2
# elif (_SWE_DG_ORDER == 3)
          traversal%cell_data(traversal%i_cell_data_index +  0)%connectivity(:) = [ 1, 2, 5] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  1)%connectivity(:) = [ 2, 3, 6] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  2)%connectivity(:) = [ 3, 4, 7] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  3)%connectivity(:) = [ 5, 6, 2] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  4)%connectivity(:) = [ 5, 6, 8] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  5)%connectivity(:) = [ 6, 7, 3] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  6)%connectivity(:) = [ 6, 7, 9] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  7)%connectivity(:) = [ 8, 9, 6] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  8)%connectivity(:) = [ 8, 9,10] + traversal%i_point_data_index - 2
# elif (_SWE_DG_ORDER == 4)
          traversal%cell_data(traversal%i_cell_data_index +  0)%connectivity(:) = [ 1, 2, 6] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  1)%connectivity(:) = [ 2, 3, 7] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  2)%connectivity(:) = [ 3, 4, 8] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  3)%connectivity(:) = [ 4, 5, 9] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  4)%connectivity(:) = [ 6, 7, 2] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  5)%connectivity(:) = [ 7, 8, 3] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  6)%connectivity(:) = [ 8, 9, 4] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  7)%connectivity(:) = [ 6, 7,10] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  8)%connectivity(:) = [ 7, 8,11] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index +  9)%connectivity(:) = [ 8, 9,12] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index + 10)%connectivity(:) = [10,11, 7] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index + 11)%connectivity(:) = [11,12, 8] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index + 12)%connectivity(:) = [10,11,13] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index + 13)%connectivity(:) = [11,12,14] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index + 14)%connectivity(:) = [13,14,11] + traversal%i_point_data_index - 2
          traversal%cell_data(traversal%i_cell_data_index + 15)%connectivity(:) = [13,14,15] + traversal%i_point_data_index - 2
# endif
          do i=1, _SWE_DG_DOFS
             ! Ever second output all cells are mirrored, this fixes it.
             if (element%cell%geometry%i_plotter_type > 0) then 
                cell_id = i
             else
                cell_id = mirrored_coords(i)
             end if
            
             traversal%point_data(traversal%i_point_data_index + i - 1)%coords = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, coords(:, i)) + cfg%offset
             traversal%point_data(traversal%i_point_data_index + i - 1)%rank = rank_MPI
             traversal%point_data(traversal%i_point_data_index + i - 1)%section_index = section%index
             traversal%point_data(traversal%i_point_data_index + i - 1)%depth = element%cell%geometry%i_depth
             traversal%point_data(traversal%i_point_data_index + i - 1)%id_in_patch = -1 ! Has no patches!
             traversal%point_data(traversal%i_point_data_index + i - 1)%i_plotter_type = element%cell%geometry%i_plotter_type
             traversal%point_data(traversal%i_point_data_index + i - 1)%refinement = element%cell%geometry%refinement
             traversal%point_data(traversal%i_point_data_index + i - 1)%troubled = element%cell%data_pers%troubled
             traversal%point_data(traversal%i_point_data_index + i - 1)%Q%b = element%cell%data_pers%Q_DG(cell_id)%b
             traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(1) = element%cell%data_pers%Q_DG(cell_id)%p(1)
             traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(2) = element%cell%data_pers%Q_DG(cell_id)%p(2)
             traversal%point_data(traversal%i_point_data_index + i - 1)%Q%h = element%cell%data_pers%Q_DG(cell_id)%h
          end do

          ! prepare for next cell
          traversal%i_cell_data_index = traversal%i_cell_data_index + num_cells
          traversal%i_point_data_index = traversal%i_point_data_index + _SWE_DG_DOFS
    else
       ! Calculated by FV.
       row=1
       col=1

       do j=1,_SWE_PATCH_ORDER * _SWE_PATCH_ORDER
          do i=1,3
             traversal%point_data(traversal%i_point_data_index + i - 1)%coords = cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, SWE_PATCH_geometry%coords(:,i, j)) + cfg%offset
          end do

          traversal%cell_data(traversal%i_cell_data_index)%troubled = element%cell%data_pers%troubled
          traversal%cell_data(traversal%i_cell_data_index)%connectivity(1:3) = (/ (traversal%i_point_data_index + j - 1, j=0,2)  /)

          ! if orientation is backwards, the plotter uses a transformation that mirrors the cell...
          ! this simple change solves the problem :)
          if (element%cell%geometry%i_plotter_type > 0) then 
             cell_id = j
          else
             cell_id = (row-1)*(row-1) + 2 * row - col
          end if

          do i=1,3
             traversal%point_data(traversal%i_point_data_index + i - 1)%rank = rank_MPI
             traversal%point_data(traversal%i_point_data_index + i - 1)%section_index = section%index
             traversal%point_data(traversal%i_point_data_index + i - 1)%depth = element%cell%geometry%i_depth
             traversal%point_data(traversal%i_point_data_index + i - 1)%id_in_patch = cell_id
             traversal%point_data(traversal%i_point_data_index + i - 1)%i_plotter_type = element%cell%geometry%i_plotter_type
             traversal%point_data(traversal%i_point_data_index + i - 1)%refinement = element%cell%geometry%refinement
             traversal%point_data(traversal%i_point_data_index + i - 1)%troubled = element%cell%data_pers%troubled
             traversal%point_data(traversal%i_point_data_index + i - 1)%Q%b = element%cell%data_pers%B(cell_id)
             traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(1) = element%cell%data_pers%HU(cell_id)
             traversal%point_data(traversal%i_point_data_index + i - 1)%Q%p(2) = element%cell%data_pers%HV(cell_id)
             ! FV calculates H + B, DG only H -> we adjust FV output to match DG-Output.
             traversal%point_data(traversal%i_point_data_index + i - 1)%Q%h = element%cell%data_pers%H(cell_id) - element%cell%data_pers%B(cell_id)
          end do

          ! prepare for next cell
          traversal%i_cell_data_index = traversal%i_cell_data_index + 1
          traversal%i_point_data_index = traversal%i_point_data_index + 3
          col = col + 1
          if (col == 2*row) then
             col = 1
             row = row + 1
          end if
       end do
    end if
    ! Prepare for next element.
    traversal%i_element_data_index = traversal%i_element_data_index + 1
# endif ! _SWE_PATCHES

  end subroutine element_op
END MODULE SWE_xml_point_output
#endif ! _SWE
