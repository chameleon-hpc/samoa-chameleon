! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


!> Generic traversal template
!> Warning: this a template module body and requires preprocessor commands before inclusion.
!>
!> The resulting method is defined as _GT_NAME
!> @author Oliver Meister

!multiple levels of indirection are necessary to properly resolve the names
#define _GT						_GT_NAME

!if no dedicated inner operators exists, use the default operators
#if defined(_GT_ELEMENT_OP) && !defined(_GT_INNER_ELEMENT_OP)
#	define _GT_INNER_ELEMENT_OP		        _GT_ELEMENT_OP
#endif

#if defined(_GT_EDGE_FIRST_TOUCH_OP) && !defined(_GT_INNER_EDGE_FIRST_TOUCH_OP)
#	define _GT_INNER_EDGE_FIRST_TOUCH_OP	_GT_EDGE_FIRST_TOUCH_OP
#endif

#if defined(_GT_EDGE_LAST_TOUCH_OP) && !defined(_GT_INNER_EDGE_LAST_TOUCH_OP)
#	define _GT_INNER_EDGE_LAST_TOUCH_OP		_GT_EDGE_LAST_TOUCH_OP
#endif

#if defined(_GT_EDGE_REDUCE_OP) && !defined(_GT_INNER_EDGE_REDUCE_OP)
#	define _GT_INNER_EDGE_REDUCE_OP		    _GT_EDGE_REDUCE_OP
#endif

#if defined(_GT_NODE_FIRST_TOUCH_OP) && !defined(_GT_INNER_NODE_FIRST_TOUCH_OP)
#	define _GT_INNER_NODE_FIRST_TOUCH_OP	_GT_NODE_FIRST_TOUCH_OP
#endif

#if defined(_GT_NODE_LAST_TOUCH_OP) && !defined(_GT_INNER_NODE_LAST_TOUCH_OP)
#	define _GT_INNER_NODE_LAST_TOUCH_OP		_GT_NODE_LAST_TOUCH_OP
#endif

#if defined(_GT_NODE_REDUCE_OP) && !defined(_GT_INNER_NODE_REDUCE_OP)
#	define _GT_INNER_NODE_REDUCE_OP		    _GT_NODE_REDUCE_OP
#endif

PRIVATE
PUBLIC _GT

!Module types:

#ifdef _GT_USE_CHAMELEON
type, bind(C)   :: t_section_metadata 
     logical :: forward
     integer(kind=c_int) :: size_cells
     integer(kind=c_int) :: size_crossed_edges_in
     integer(kind=c_int) :: size_crossed_edges_out   
     integer(kind=c_int) :: size_color_edges_in
     integer(kind=c_int) :: size_color_edges_out 
     integer(kind=c_int) :: size_nodes_in
     integer(kind=c_int) :: size_nodes_out
     integer(kind=c_int) :: size_boundary_edges_red
     integer(kind=c_int) :: size_boundary_edges_green
     integer(kind=c_int) :: size_boundary_nodes_red
     integer(kind=c_int) :: size_boundary_nodes_green

end type
#endif

!> Traversal element ring buffer structure that provides local storage for some of the grid data
type, extends(t_element_base) :: t_traversal_element
	type(num_cell_data_temp)							    :: cell_data_temp							!< cell temporary data

#	if defined(_GT_EDGES)
		type(t_edge_data)				                    :: next_edge_data            				!< next crossed edge temp + local data
#	endif

	type(t_traversal_element), pointer			   		    :: previous, next							!< pointer to previous and next traversal element in the ringbuffer
end type

type, extends(num_traversal_data) :: t_traversal
    type(t_statistics)                                      :: stats

    contains

    procedure, private, pass :: assign_traversal => t_traversal_assign

    generic, public :: assignment(=) => assign_traversal
end type

type, extends(t_traversal) :: t_thread_traversal
    type(t_traversal_element)                               :: elements(8)                                                      !< Element ring buffer (must c
end type

#define t_section_traversal _GT

type, extends(t_traversal) :: _GT
    type(t_section_traversal), pointer                      :: sections(:) => null()                   !< section data
    type(t_thread_traversal), pointer                       :: threads(:) => null()             !< thread data
    integer                                                 :: mpi_node_type, mpi_edge_type

    contains

    procedure, private, pass :: assign_gt => gt_assign

    procedure, public, pass :: create
    procedure, public, pass :: destroy
    procedure, public, pass :: traverse
    procedure, public, pass :: reduce_stats
    procedure, public, pass :: clear_stats

    generic, public :: assignment(=) => assign_gt
end type

contains

subroutine t_traversal_assign(t1, t2)
    class(t_traversal), intent(inout) :: t1
    type(t_traversal), intent(in)  :: t2

    t1%num_traversal_data = t2%num_traversal_data
    t1%stats = t2%stats
end subroutine

subroutine gt_assign(gt1, gt2)
    class(_GT), intent(inout)   :: gt1
    type(_GT), intent(in)       :: gt2

    gt1%t_traversal = gt2%t_traversal
    gt1%mpi_node_type = gt2%mpi_node_type
    gt1%mpi_edge_type = gt2%mpi_edge_type

    !pointers are not copied, they must be set manually
end subroutine

subroutine reduce_stats(traversal, mpi_op, global)
    class(_GT)              :: traversal
    integer, intent(in)     :: mpi_op
    logical, intent(in)     :: global

    if (associated(traversal%threads)) then
        call traversal%stats%reduce(traversal%threads(:)%stats, mpi_op)
    end if

    if (global) then
        call traversal%stats%reduce(mpi_op)
    end if
end subroutine

subroutine clear_stats(traversal)
    class(_GT)              :: traversal
    integer                 :: i

    if (associated(traversal%threads)) then
        do i = 1, size(traversal%threads)
            call traversal%threads(i)%stats%clear()
        end do
    end if

    call traversal%stats%clear()
end subroutine

subroutine create(traversal)
    class(_GT)      :: traversal
	integer         :: i_error

#    if defined(_GT_NODE_MPI_TYPE)
        call create_node_mpi_type(traversal%mpi_node_type)
#    endif

#    if defined(_GT_EDGE_MPI_TYPE)
        call create_edge_mpi_type(traversal%mpi_edge_type)
#    endif
end subroutine

subroutine destroy(traversal)
    class(_GT)      :: traversal
	integer         :: i_error

#    if defined(_GT_NODE_MPI_TYPE) && defined(_MPI)
        call MPI_Type_free(traversal%mpi_node_type, i_error); assert_eq(i_error, 0)
#    endif

#    if defined(_GT_EDGE_MPI_TYPE) && defined(_MPI)
        call MPI_Type_free(traversal%mpi_edge_type, i_error); assert_eq(i_error, 0)
#    endif

    if (associated(traversal%sections)) then
        deallocate(traversal%sections, stat = i_error); assert_eq(i_error, 0)
    end if

    if (associated(traversal%threads)) then
        deallocate(traversal%threads, stat = i_error); assert_eq(i_error, 0)
    end if
end subroutine

function edge_merge_wrapper_op(local_edges, neighbor_edges) result(l_conform)
    type(t_edge_data), intent(inout)    :: local_edges
    type(t_edge_data), intent(in)       :: neighbor_edges
    logical                             :: l_conform

    assert_eq(local_edges%min_distance, neighbor_edges%min_distance)
    assert(local_edges%owned_locally)

#   if defined(_GT_EDGE_MERGE_OP)
        call _GT_EDGE_MERGE_OP(local_edges, neighbor_edges)
#   endif

    l_conform = .true.
end function

function node_merge_wrapper_op(local_nodes, neighbor_nodes) result(l_conform)
    type(t_node_data), intent(inout)    :: local_nodes
    type(t_node_data), intent(in)       :: neighbor_nodes
    logical                             :: l_conform

    assert_eq(local_nodes%distance, neighbor_nodes%distance)
    assert_veq(local_nodes%position, neighbor_nodes%position)
    assert(local_nodes%owned_locally)

#   if defined(_GT_NODE_MERGE_OP)
        call _GT_NODE_MERGE_OP(local_nodes, neighbor_nodes)
#   endif

    l_conform = .true.
end function

function edge_write_wrapper_op(local_edges, neighbor_edges) result(l_conform)
    type(t_edge_data), intent(inout)    :: local_edges
    type(t_edge_data), intent(in)       :: neighbor_edges
    logical                             :: l_conform

    assert_eq(local_edges%min_distance, neighbor_edges%min_distance)
    assert(.not. local_edges%owned_locally)
    assert(neighbor_edges%owned_locally)

#   if defined(_GT_EDGE_WRITE_OP)
        call _GT_EDGE_WRITE_OP(local_edges, neighbor_edges)
#   else
        local_edges%data_pers = neighbor_edges%data_pers
        local_edges%data_temp = neighbor_edges%data_temp
#   endif

    l_conform = .true.
end function

function node_write_wrapper_op(local_nodes, neighbor_nodes) result(l_conform)
    type(t_node_data), intent(inout)    :: local_nodes
    type(t_node_data), intent(in)       :: neighbor_nodes
    logical                             :: l_conform

    assert_eq(local_nodes%distance, neighbor_nodes%distance)
    assert_veq(local_nodes%position, neighbor_nodes%position)
    assert(.not. local_nodes%owned_locally)
    assert(neighbor_nodes%owned_locally)

#   if defined(_GT_NODE_WRITE_OP)
        call _GT_NODE_WRITE_OP(local_nodes, neighbor_nodes)
#   else
        local_nodes%data_pers = neighbor_nodes%data_pers
        local_nodes%data_temp = neighbor_nodes%data_temp
#   endif

    l_conform = .true.
end function

!*****************************************************************
! Generic traversal
!*****************************************************************

!> Generic iterative traversal subroutine
!> @author Oliver Meister
subroutine traverse(traversal, grid)
	class(_GT), intent(inout)	                        :: traversal
	type(t_grid), intent(inout)					        :: grid

	integer (kind = GRID_SI)                            :: i_section, i_first_local_section, i_last_local_section
	integer (kind = GRID_SI)                            :: i_error

	type(t_statistics)                                  :: thread_stats

#ifdef _GT_USE_CHAMELEON
       type(t_section_metadata), dimension(:), allocatable,target    :: section_metadata
       type(map_entry), dimension(:,:), allocatable, target :: map_entries
       
       !procedure(traverse_section_wrapper_chameleon), pointer:: funPointer => traverse_section_wrapper_chameleon

       type(c_ptr) :: anno_ptr, task_c_ptr
       integer :: result_val
#endif 

    if (.not. associated(traversal%sections) .or. size(traversal%sections) .ne. grid%sections%get_size()) then
        !$omp barrier

        !$omp single
        if (associated(traversal%sections)) then
            deallocate(traversal%sections, stat = i_error); assert_eq(i_error, 0)
        end if

        allocate(traversal%sections(grid%sections%get_size()), stat = i_error); assert_eq(i_error, 0)
        !$omp end single
    end if

    if (.not. associated(traversal%threads) .or. size(traversal%threads) .ne. cfg%i_threads) then
        !$omp barrier

        !$omp single
        if (associated(traversal%threads)) then
            deallocate(traversal%threads, stat = i_error); assert_eq(i_error, 0)
        end if

        allocate(traversal%threads(cfg%i_threads), stat = i_error); assert_eq(i_error, 0)
        !$omp end single

        call create_ringbuffer(traversal%threads(i_thread)%elements)
    end if

    assert_eq(size(traversal%threads), omp_get_max_threads())

    call grid%get_local_sections(i_first_local_section, i_last_local_section)

    do i_section = i_first_local_section, i_last_local_section
        call traversal%sections(i_section)%stats%clear()
    end do

#   if defined(_ASAGI_TIMING)
        call grid%sections%elements_alloc(i_first_local_section : i_last_local_section)%stats%clear_time(asagi_time)
#   endif

    call thread_stats%start_time(traversal_time)

    call thread_stats%start_time(barrier_time)

    select type (traversal)
        type is (_GT)
            !$omp single
            call pre_traversal_grid(traversal, grid)
            !$omp end single nowait
        class default
            assert(.false.)
    end select

    call thread_stats%stop_time(barrier_time)

#   if defined(_GT_SKELETON_OP)
        do i_section = i_first_local_section, i_last_local_section
            assert_eq(i_section, grid%sections%elements_alloc(i_section)%index)

            call traversal%sections(i_section)%stats%start_time(pre_compute_time)
            call boundary_skeleton(traversal%sections(i_section), grid%sections%elements_alloc(i_section))
            call traversal%sections(i_section)%stats%stop_time(pre_compute_time)
        end do
#   endif

    !$omp barrier

    do i_section = i_first_local_section, i_last_local_section
#       if defined(_OPENMP_TASKS)
            !$omp task default(shared) firstprivate(i_section) mergeable
#       endif

        call traversal%sections(i_section)%stats%start_time(pre_compute_time)
        call pre_traversal_wrapper(traversal%sections(i_section), grid%sections%elements_alloc(i_section))
        call traversal%sections(i_section)%stats%stop_time(pre_compute_time)

        !WARNING: Do not use thread_stats to track sync times here,
        !a race condition can occur otherwise.
        call traversal%sections(i_section)%stats%start_time(sync_time)
#       if !defined(_GT_NODE_MPI_TYPE) && !defined(_GT_EDGE_MPI_TYPE)
            call recv_mpi_boundary(grid%sections%elements_alloc(i_section))
#       elif defined(_GT_NODE_MPI_TYPE) && !defined(_GT_EDGE_MPI_TYPE)
            call recv_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_node_type_optional=traversal%mpi_node_type)
#       elif defined(_GT_EDGE_MPI_TYPE) && !defined(_GT_NODE_MPI_TYPE)
            call recv_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_edge_type_optional=traversal%mpi_edge_type)
#       else
            call recv_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_edge_type_optional=traversal%mpi_edge_type, mpi_node_type_optional=traversal%mpi_node_type)
#       endif
        call traversal%sections(i_section)%stats%stop_time(sync_time)

#       if defined(_OPENMP_TASKS)
            !$omp end task
#       endif
    end do

#   if defined(_OPENMP_TASKS)
        !$omp taskwait
#   endif

    call thread_stats%start_time(sync_time)
    call duplicate_boundary_data(grid, edge_write_wrapper_op, node_write_wrapper_op)

    !wait until all boundary data has been copied
    !$omp barrier
    call thread_stats%stop_time(sync_time)

#if defined _GT_USE_CHAMELEON
    allocate(section_metadata(i_first_local_section:i_last_local_section))
    allocate(map_entries(14,i_first_local_section:i_last_local_section))
    ! write(*,*)  'Executing Packing/Unpacking version'
    do i_section = i_first_local_section, i_last_local_section
        call traversal%sections(i_section)%stats%start_time(packing_time)
        
        map_entries(1, i_section)%valptr = c_loc(traversal%sections(i_section))
        map_entries(1, i_section)%size = sizeof(traversal%sections(i_section))
        map_entries(1, i_section)%type = 3

        map_entries(2, i_section)%valptr = c_loc(section_metadata(i_section))
        map_entries(2, i_section)%size = sizeof(section_metadata(i_section))
        map_entries(2, i_section)%type = 3

        map_entries(3, i_section)%valptr = c_loc(grid%sections%elements_alloc(i_section)%t_global_data)
        map_entries(3, i_section)%size = sizeof(grid%sections%elements_alloc(i_section)%t_global_data)
        map_entries(3, i_section)%type = 3

        map_entries(4, i_section)%valptr = c_loc(grid%sections%elements_alloc(i_section)%cells%get_c_pointer())
        map_entries(4, i_section)%size = sizeof(grid%sections%elements_alloc(i_section)%cells%elements)
        map_entries(4, i_section)%type = 3
        
        !write(*,*)  'size elements', size(grid%sections%elements_alloc(i_section)%cells%elements), 'sizeof', sizeof(grid%sections%elements_alloc(i_section)%cells%elements)

        map_entries(5, i_section)%valptr = c_loc(grid%sections%elements_alloc(i_section)%crossed_edges_in%get_c_pointer())
        map_entries(5, i_section)%size = sizeof(grid%sections%elements_alloc(i_section)%crossed_edges_in%elements)
        map_entries(5, i_section)%type = 1

        map_entries(6, i_section)%valptr = c_loc(grid%sections%elements_alloc(i_section)%crossed_edges_out%get_c_pointer())
        map_entries(6, i_section)%size = sizeof(grid%sections%elements_alloc(i_section)%crossed_edges_out%elements)
        map_entries(6, i_section)%type = 2

        map_entries(7, i_section)%valptr = c_loc(grid%sections%elements_alloc(i_section)%color_edges_in%get_c_pointer())
        map_entries(7, i_section)%size = sizeof(grid%sections%elements_alloc(i_section)%color_edges_in%elements)
        map_entries(7, i_section)%type = 1

        map_entries(8, i_section)%valptr = c_loc(grid%sections%elements_alloc(i_section)%color_edges_out%get_c_pointer())
        map_entries(8, i_section)%size = sizeof(grid%sections%elements_alloc(i_section)%color_edges_out%elements)
        map_entries(8, i_section)%type = 2

        map_entries(9, i_section)%valptr = c_loc(grid%sections%elements_alloc(i_section)%nodes_in%get_c_pointer())
        map_entries(9, i_section)%size = sizeof(grid%sections%elements_alloc(i_section)%nodes_in%elements)
        map_entries(9, i_section)%type = 1

        map_entries(10, i_section)%valptr = c_loc(grid%sections%elements_alloc(i_section)%nodes_out%get_c_pointer())
        map_entries(10, i_section)%size = sizeof(grid%sections%elements_alloc(i_section)%nodes_out%elements)
        map_entries(10, i_section)%type = 2

        map_entries(11, i_section)%valptr = c_loc(grid%sections%elements_alloc(i_section)%boundary_edges(RED)%get_c_pointer())
        map_entries(11, i_section)%size = sizeof(grid%sections%elements_alloc(i_section)%boundary_edges(RED)%elements)
        map_entries(11, i_section)%type = 3

        map_entries(12, i_section)%valptr = c_loc(grid%sections%elements_alloc(i_section)%boundary_edges(GREEN)%get_c_pointer())
        map_entries(12, i_section)%size = sizeof(grid%sections%elements_alloc(i_section)%boundary_edges(GREEN)%elements)
        map_entries(12, i_section)%type = 3

        map_entries(13, i_section)%valptr = c_loc(grid%sections%elements_alloc(i_section)%boundary_nodes(RED)%get_c_pointer())
        map_entries(13, i_section)%size = sizeof(grid%sections%elements_alloc(i_section)%boundary_nodes(RED)%elements)
        map_entries(13, i_section)%type = 3

        map_entries(14, i_section)%valptr = c_loc(grid%sections%elements_alloc(i_section)%boundary_nodes(GREEN)%get_c_pointer())
        map_entries(14, i_section)%size = sizeof(grid%sections%elements_alloc(i_section)%boundary_nodes(GREEN)%elements)
        map_entries(14, i_section)%type = 3

        section_metadata(i_section)%forward = grid%sections%elements_alloc(i_section)%cells%forward
        section_metadata(i_section)%size_cells = size(grid%sections%elements_alloc(i_section)%cells%elements)
        section_metadata(i_section)%size_crossed_edges_in = size(grid%sections%elements_alloc(i_section)%crossed_edges_in%elements)
        section_metadata(i_section)%size_crossed_edges_out = size(grid%sections%elements_alloc(i_section)%crossed_edges_out%elements)
        section_metadata(i_section)%size_color_edges_in = size(grid%sections%elements_alloc(i_section)%color_edges_in%elements)
        section_metadata(i_section)%size_color_edges_out = size(grid%sections%elements_alloc(i_section)%color_edges_out%elements)
        section_metadata(i_section)%size_nodes_in = size(grid%sections%elements_alloc(i_section)%nodes_in%elements)
        section_metadata(i_section)%size_nodes_out = size(grid%sections%elements_alloc(i_section)%nodes_out%elements)
        section_metadata(i_section)%size_boundary_edges_red = size(grid%sections%elements_alloc(i_section)%boundary_edges(RED)%elements)
        section_metadata(i_section)%size_boundary_edges_green = size(grid%sections%elements_alloc(i_section)%boundary_edges(GREEN)%elements)
        section_metadata(i_section)%size_boundary_nodes_red = size(grid%sections%elements_alloc(i_section)%boundary_nodes(RED)%elements)
        section_metadata(i_section)%size_boundary_nodes_green = size(grid%sections%elements_alloc(i_section)%boundary_nodes(GREEN)%elements)
        call traversal%sections(i_section)%stats%stop_time(packing_time)

#if defined _GT_USE_CHAMELEON_CALL
    !    write(*,*)  'Executing Chameleon Tasks'

        task_c_ptr = chameleon_create_task(c_funloc(traverse_section_wrapper_chameleon), 14, map_entries(:,i_section))
        i_error = chameleon_add_task_fortran(task_c_ptr)
        ! anno_ptr    = chameleon_create_annotation_container_fortran()
        ! i_error     = chameleon_set_annotation_int_fortran(anno_ptr, grid%sections%elements_alloc(i_section)%cells%get_size())
        ! i_error     = chameleon_add_task_manual_fortran_w_annotations(c_funloc(traverse_section_wrapper_chameleon), 14, c_loc(map_entries(:,i_section)), anno_ptr)
        ! result_val  = chameleon_get_annotation_int_fortran(anno_ptr)
        ! if (result_val .ne. grid%sections%elements_alloc(i_section)%cells%get_size()) then
        !     write(*,*)  'Annotation not matching', grid%sections%elements_alloc(i_section)%cells%get_size(), result_val
        ! end if 
#else
#       if defined(_OPENMP_TASKS)
        !$omp task default(shared) firstprivate(i_section) mergeable
#       endif
        !call traversal%sections(i_section)%stats%start_time(inner_compute_time)
        call traverse_section_wrapper_chameleon(traversal%sections(i_section),\
                                              section_metadata(i_section),\
                                              grid%sections%elements_alloc(i_section)%t_global_data,\
                                              c_loc(grid%sections%elements_alloc(i_section)%cells%get_c_pointer()),\
                                              c_loc(grid%sections%elements_alloc(i_section)%crossed_edges_in%get_c_pointer()),\
                                              c_loc(grid%sections%elements_alloc(i_section)%crossed_edges_out%get_c_pointer()),\
                                              c_loc(grid%sections%elements_alloc(i_section)%color_edges_in%get_c_pointer()),\
                                              c_loc(grid%sections%elements_alloc(i_section)%color_edges_in%get_c_pointer()),\
                                              c_loc(grid%sections%elements_alloc(i_section)%nodes_in%get_c_pointer()),\
                                              c_loc(grid%sections%elements_alloc(i_section)%nodes_out%get_c_pointer()),\
                                              c_loc(grid%sections%elements_alloc(i_section)%boundary_edges(RED)%get_c_pointer()),\
                                              c_loc(grid%sections%elements_alloc(i_section)%boundary_edges(GREEN)%get_c_pointer()),\
                                              c_loc(grid%sections%elements_alloc(i_section)%boundary_nodes(RED)%get_c_pointer()),\
                                              c_loc(grid%sections%elements_alloc(i_section)%boundary_nodes(GREEN)%get_c_pointer()))
        !call traversal%sections(i_section)%stats%stop_time(inner_compute_time)
#       if defined(_OPENMP_TASKS)
        !$omp end task
#       endif
#endif
    end do

#if defined _GT_USE_CHAMELEON_CALL
    !call thread_stats%start_time(inner_compute_time)
    i_error = chameleon_distributed_taskwait(0)
    !call thread_stats%stop_time(inner_compute_time)
#else
#   if defined(_OPENMP_TASKS)
        !$omp taskwait
#   endif
#endif
    deallocate(map_entries)
    deallocate(section_metadata)

    do i_section = i_first_local_section, i_last_local_section
        !WARNING: Do not use thread_stats to track sync times here,
        !a race condition can occur otherwise.
        call traversal%sections(i_section)%stats%start_time(sync_time)
#       if !defined(_GT_NODE_MPI_TYPE) && !defined(_GT_EDGE_MPI_TYPE)
            call send_mpi_boundary(grid%sections%elements_alloc(i_section))
#       elif defined(_GT_NODE_MPI_TYPE) && !defined(_GT_EDGE_MPI_TYPE)
            call send_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_node_type_optional=traversal%mpi_node_type)
#       elif defined(_GT_EDGE_MPI_TYPE) && !defined(_GT_NODE_MPI_TYPE)
            call send_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_edge_type_optional=traversal%mpi_edge_type)
#       else
            call send_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_edge_type_optional=traversal%mpi_edge_type, mpi_node_type_optional=traversal%mpi_node_type)
#       endif
        call traversal%sections(i_section)%stats%stop_time(sync_time)
    end do 

    call thread_stats%start_time(sync_time)

    !wait until all computation is done
    !$omp barrier
    call thread_stats%stop_time(sync_time)
#else
    do i_section = i_first_local_section, i_last_local_section
#       if defined(_OPENMP_TASKS)
            !$omp task default(shared) firstprivate(i_section) mergeable
#       endif

        call traversal%sections(i_section)%stats%start_time(inner_compute_time)
        call traverse_section_wrapper(traversal%threads(i_thread), traversal%sections(i_section), grid%threads%elements(i_thread), grid%sections%elements_alloc(i_section))
        call traversal%sections(i_section)%stats%stop_time(inner_compute_time)

        !WARNING: Do not use thread_stats to track sync times here,
        !a race condition can occur otherwise.
        call traversal%sections(i_section)%stats%start_time(sync_time)
#       if !defined(_GT_NODE_MPI_TYPE) && !defined(_GT_EDGE_MPI_TYPE)
            call send_mpi_boundary(grid%sections%elements_alloc(i_section))
#       elif defined(_GT_NODE_MPI_TYPE) && !defined(_GT_EDGE_MPI_TYPE)
            call send_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_node_type_optional=traversal%mpi_node_type)
#       elif defined(_GT_EDGE_MPI_TYPE) && !defined(_GT_NODE_MPI_TYPE)
            call send_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_edge_type_optional=traversal%mpi_edge_type)
#       else
            call send_mpi_boundary(grid%sections%elements_alloc(i_section), mpi_edge_type_optional=traversal%mpi_edge_type, mpi_node_type_optional=traversal%mpi_node_type)
#       endif
        call traversal%sections(i_section)%stats%stop_time(sync_time)

#       if defined(_OPENMP_TASKS)
            !$omp end task
#       endif
    end do
    call thread_stats%start_time(sync_time)

#   if defined(_OPENMP_TASKS)
        !$omp taskwait
#   endif

    !wait until all computation is done
    !$omp barrier
    call thread_stats%stop_time(sync_time)
#endif

    !sync and call post traversal operator
    call thread_stats%start_time(sync_time)
#   if !defined(_GT_NODE_MPI_TYPE) && !defined(_GT_EDGE_MPI_TYPE)
        call collect_boundary_data(grid, edge_merge_wrapper_op, node_merge_wrapper_op)
#   elif defined(_GT_NODE_MPI_TYPE) && !defined(_GT_EDGE_MPI_TYPE)
        call collect_boundary_data(grid, edge_merge_wrapper_op, node_merge_wrapper_op, mpi_node_type_optional=traversal%mpi_node_type)
#   elif defined(_GT_EDGE_MPI_TYPE) && !defined(_GT_NODE_MPI_TYPE)
        call collect_boundary_data(grid, edge_merge_wrapper_op, node_merge_wrapper_op, mpi_edge_type_optional=traversal%mpi_edge_type)
#   else
        call collect_boundary_data(grid, edge_merge_wrapper_op, node_merge_wrapper_op, mpi_node_type_optional=traversal%mpi_node_type, mpi_edge_type_optional=traversal%mpi_edge_type)
#   endif

    do i_section = i_first_local_section, i_last_local_section
        assert_eq(i_section, grid%sections%elements_alloc(i_section)%index)
        call traversal%sections(i_section)%stats%start_time(post_compute_time)
        call post_traversal_wrapper(traversal%sections(i_section), grid%sections%elements_alloc(i_section))
        call traversal%sections(i_section)%stats%stop_time(post_compute_time)
    end do

    call grid%reverse()

    !$omp barrier
    call thread_stats%stop_time(sync_time)

    call thread_stats%start_time(barrier_time)

    select type (traversal)
        type is (_GT)
            !$omp single
            call post_traversal_grid(traversal, grid)
            !$omp end single
        class default
            assert(.false.)
    end select

    call thread_stats%stop_time(barrier_time)
    call thread_stats%stop_time(traversal_time)
    
    do i_section = i_first_local_section, i_last_local_section
        call set_stats_counters(traversal%sections(i_section)%stats, grid%sections%elements_alloc(i_section))

        grid%sections%elements_alloc(i_section)%stats = grid%sections%elements_alloc(i_section)%stats + traversal%sections(i_section)%stats

#       if defined(_ASAGI_TIMING)
            !HACK: in lack of a better method, we reduce ASAGI timing data like this for now - should be changed in the long run, so that stats belongs to the section and not the traversal
            call traversal%sections(i_section)%stats%add_time(asagi_time, grid%sections%elements_alloc(i_section)%stats%get_time(asagi_time))
#       endif

        thread_stats = thread_stats + traversal%sections(i_section)%stats
    end do

    ! used for auto-tuning LB
	! but don't consider the displacement traversal here... it is only used for one timestep!
#   ifndef _GT_SKIP_LB_TIMING
	!$omp critical
        grid%r_computation_time_since_last_LB = grid%r_computation_time_since_last_LB + thread_stats%get_time(pre_compute_time) + thread_stats%get_time(inner_compute_time) + thread_stats%get_time(post_compute_time)
    !$omp end critical
#   endif

    traversal%threads(i_thread)%stats = traversal%threads(i_thread)%stats + thread_stats
    grid%threads%elements(i_thread)%stats = grid%threads%elements(i_thread)%stats + thread_stats
end subroutine

#if defined(_GT_USE_CHAMELEON)
subroutine traverse_section_wrapper_chameleon( section_traversal,&
                                               section_metadata,&
                                               global_data,&
                                               cells,&
                                               crossed_edges_in,&
                                               crossed_edges_out,&
                                               color_edges_in,&
                                               color_edges_out,&
                                               nodes_in,&
                                               nodes_out,&
                                               boundary_edges_red,&
                                               boundary_edges_green,&
                                               boundary_nodes_red,&
                                               boundary_nodes_green) bind(C) !required? leads to problems with gnu
    use SFC_data_types
    use Grid_section
    use Cell_stream
    use Crossed_edge_stream
    use Color_edge_stream
    use Node_stream

    use Boundary_edge_stream
    use Boundary_node_stream

    type(t_section_traversal), intent(inout)        :: section_traversal
    type(t_section_metadata), intent(inout), target         :: section_metadata
    type(t_global_data), intent(inout)              :: global_data
    type(c_ptr),value,  intent(in)                     :: cells, crossed_edges_in, crossed_edges_out, color_edges_in, color_edges_out, nodes_in, nodes_out, boundary_edges_red, boundary_edges_green, boundary_nodes_red, boundary_nodes_green


    type(t_section_traversal)                       :: section_traversal_local
    type(t_thread_traversal)                        :: thread_traversal_local
    type(t_grid_thread)                             :: thread_local
    type(t_grid_section)                            :: section_local

    type(t_cell_stream_data), pointer              :: cells_ptr(:)
    type(t_crossed_edge_stream_data), pointer      :: crossed_edges_in_ptr(:)
    type(t_crossed_edge_stream_data), pointer     :: crossed_edges_out_ptr(:)
    type(t_color_edge_stream_data), pointer        :: color_edges_in_ptr(:)
    type(t_color_edge_stream_data), pointer        :: color_edges_out_ptr(:)
    type(t_node_stream_data), pointer              :: nodes_in_ptr(:)
    type(t_node_stream_data), pointer              :: nodes_out_ptr(:)
    type(t_edge_data), pointer :: boundary_edges_red_ptr(:), boundary_edges_green_ptr(:)
    type(t_node_data), pointer :: boundary_nodes_red_ptr(:), boundary_nodes_green_ptr(:)    

    logical                                        :: forward

    forward = section_metadata%forward

    section_traversal_local = section_traversal
    !section_local = section
    !thread_local = thread
    section_local%t_global_data = global_data

    call section_traversal_local%stats%start_time(inner_compute_time)
    
    call c_f_pointer(cells, cells_ptr, [section_metadata%size_cells])
    if(forward) then  
      section_local%cells%elements => cells_ptr
    else
      section_local%cells%elements => cells_ptr(section_metadata%size_cells: 1 :-1)
    end if
    section_local%cells%forward = forward
   
    call c_f_pointer(crossed_edges_in, crossed_edges_in_ptr, [section_metadata%size_crossed_edges_in])
    if(forward) then  
      section_local%crossed_edges_in%elements => crossed_edges_in_ptr
    else
      section_local%crossed_edges_in%elements => crossed_edges_in_ptr(section_metadata%size_crossed_edges_in: 1 :-1)
    end if


    call c_f_pointer(crossed_edges_out, crossed_edges_out_ptr, [section_metadata%size_crossed_edges_out])
    if(forward) then  
      section_local%crossed_edges_out%elements => crossed_edges_out_ptr
    else
      section_local%crossed_edges_out%elements => crossed_edges_out_ptr(section_metadata%size_crossed_edges_out: 1 :-1)
    end if

    call c_f_pointer(color_edges_in, color_edges_in_ptr, [section_metadata%size_color_edges_in])
    if(forward) then  
      section_local%color_edges_in%elements => color_edges_in_ptr
    else
      section_local%color_edges_in%elements => color_edges_in_ptr(section_metadata%size_color_edges_in: 1 :-1)
    end if

    call c_f_pointer(color_edges_out, color_edges_out_ptr, [section_metadata%size_color_edges_out])
    if(forward) then  
      section_local%color_edges_out%elements => color_edges_out_ptr
    else
      section_local%color_edges_out%elements => color_edges_out_ptr(section_metadata%size_color_edges_out: 1 :-1)
    end if

    call c_f_pointer(nodes_in, nodes_in_ptr, [section_metadata%size_nodes_in])
    if(forward) then  
      section_local%nodes_in%elements => nodes_in_ptr
    else
      section_local%nodes_in%elements => nodes_in_ptr(section_metadata%size_nodes_in: 1 :-1)
    end if

    call c_f_pointer(nodes_out, nodes_out_ptr, [section_metadata%size_nodes_out])
    if(forward) then  
      section_local%nodes_out%elements => nodes_out_ptr
    else
      section_local%nodes_out%elements => nodes_out_ptr(section_metadata%size_nodes_out: 1 :-1)
    end if

    call c_f_pointer(boundary_edges_red, boundary_edges_red_ptr, [section_metadata%size_boundary_edges_red])
    if(forward) then  
      section_local%boundary_edges(RED)%elements => boundary_edges_red_ptr
    else
      section_local%boundary_edges(RED)%elements => boundary_edges_red_ptr(section_metadata%size_boundary_edges_red: 1 :-1)
    end if

    call c_f_pointer(boundary_edges_green, boundary_edges_green_ptr, [section_metadata%size_boundary_edges_green])
    if(forward) then  
      section_local%boundary_edges(GREEN)%elements => boundary_edges_green_ptr
    else
      section_local%boundary_edges(GREEN)%elements => boundary_edges_green_ptr(section_metadata%size_boundary_edges_green: 1 :-1)
    end if

    call c_f_pointer(boundary_nodes_red, boundary_nodes_red_ptr, [section_metadata%size_boundary_nodes_red])
    if(forward) then  
      section_local%boundary_nodes(RED)%elements => boundary_nodes_red_ptr
    else
      section_local%boundary_nodes(RED)%elements => boundary_nodes_red_ptr(section_metadata%size_boundary_nodes_red: 1 :-1)
    end if

    call c_f_pointer(boundary_nodes_green, boundary_nodes_green_ptr, [section_metadata%size_boundary_nodes_green])
    if(forward) then  
      section_local%boundary_nodes(GREEN)%elements => boundary_nodes_green_ptr
    else
      section_local%boundary_nodes(GREEN)%elements => boundary_nodes_green_ptr(section_metadata%size_boundary_nodes_green: 1 :-1)
    end if

    call create_ringbuffer(thread_traversal_local%elements)
    call thread_local%create(section_local%max_dest_stack-section_local%min_dest_stack)
 
    call traverse_section(thread_traversal_local, section_traversal_local, thread_local, section_local)

    call thread_local%destroy()

    call section_traversal_local%stats%stop_time(inner_compute_time)
   
    section_traversal = section_traversal_local
    global_data = section_local%t_global_data

end subroutine
#endif

subroutine traverse_section_wrapper(thread_traversal, section_traversal, thread, section)
    type(t_thread_traversal), intent(inout)         :: thread_traversal
    type(t_section_traversal), intent(inout)        :: section_traversal
    type(t_grid_thread), intent(inout)              :: thread
    type(t_grid_section), intent(inout)             :: section

    type(t_section_traversal)                       :: section_traversal_local
    type(t_grid_thread)                             :: thread_local
    type(t_grid_section)                            :: section_local

    section_traversal_local = section_traversal
    thread_local = thread
    section_local = section

    call traverse_section(thread_traversal, section_traversal_local, thread_local, section_local)

    section_traversal = section_traversal_local
    thread = thread_local
    section = section_local
end subroutine

subroutine pre_traversal_wrapper(section_traversal, section)
	type(t_section_traversal), intent(inout)        :: section_traversal
	type(t_grid_section), intent(inout)			    :: section

	type(t_section_traversal)	                    :: section_traversal_local
	type(t_grid_section)		                    :: section_local

    section_traversal_local = section_traversal
    section_local = section

    call pre_traversal(section_traversal_local, section_local)

    section_traversal = section_traversal_local
    section = section_local
end subroutine

subroutine post_traversal_wrapper(section_traversal, section)
	type(t_section_traversal), intent(inout)        :: section_traversal
	type(t_grid_section), intent(inout)		        :: section

	type(t_section_traversal)	                    :: section_traversal_local
	type(t_grid_section)		                    :: section_local

    section_traversal_local = section_traversal
    section_local = section

    call post_traversal(section_traversal_local, section_local)

    section_traversal = section_traversal_local
    section = section_local
end subroutine

!> Generic iterative traversal subroutine
subroutine traverse_section(thread_traversal, section_traversal, thread, section)
	type(t_thread_traversal), target, intent(inout)	    :: thread_traversal
	type(t_section_traversal), target, intent(inout)	:: section_traversal
	type(t_grid_thread), intent(inout)					:: thread
	type(t_grid_section), intent(inout)					:: section

	! local variables
	integer (kind = GRID_SI)							:: i, i_current_element, i_next_element

#	if (_DEBUG_LEVEL > 4)
		_log_write(5, '(2X, A)') "section input state :"
		call section%print()
#	endif

#	if (_DEBUG_LEVEL > 5)
		_log_write(6, '(2X, A)') "input cells:"
		do i = lbound(section%cells%elements, 1), ubound(section%cells%elements, 1)
			_log_write(6, '(3X, I0, X, A)') i, section%cells%elements(i)%to_string()
		end do
		_log_write(6, '(A)') ""
#	endif

	select case (section%cells%get_size())
        case (1)
            !process the only element
            call init(section, thread_traversal%elements(1))
            call leaf(section_traversal, thread, section, thread_traversal%elements(1))
        case (2:)
            call init(section, thread_traversal%elements(1))

            !process first element
            call init(section, thread_traversal%elements(2))
            call leaf(section_traversal, thread, section, thread_traversal%elements(1))

            i_current_element = 2

            do
                i_next_element = mod(i_current_element, 8) + 1

                select case (thread_traversal%elements(i_current_element)%cell%geometry%i_entity_types)
                    case (INNER_OLD)
                        !init next element for the skeleton operator
                        call init(section, thread_traversal%elements(i_next_element))
                        call old_leaf(section_traversal, thread, section, thread_traversal%elements(i_current_element))
                    case (INNER_NEW)
                        !init next element for the skeleton operator
                        call init(section, thread_traversal%elements(i_next_element))
                        call new_leaf(section_traversal, thread, section, thread_traversal%elements(i_current_element))
                    case (INNER_OLD_BND)
                        !init next element for the skeleton operator
                        call init(section, thread_traversal%elements(i_next_element))
                        call old_bnd_leaf(section_traversal, thread, section, thread_traversal%elements(i_current_element))
                    case (INNER_NEW_BND)
                        !init next element for the skeleton operator
                        call init(section, thread_traversal%elements(i_next_element))
                        call new_bnd_leaf(section_traversal, thread, section, thread_traversal%elements(i_current_element))
                    case default
                        !this should happen only for the last element
                        exit
                end select

                i_current_element = i_next_element
            end do

            !process last element
            call leaf(section_traversal, thread, section, thread_traversal%elements(i_current_element))
    end select

#	if (_DEBUG_LEVEL > 4)
		_log_write(5, '(2X, A)') "section output state :"
		call section%print()
#	endif

#	if (_DEBUG_LEVEL > 5)
		_log_write(6, '(2X, A)') "output cells:"
		do i = lbound(section%cells%elements, 1), ubound(section%cells%elements, 1)
			_log_write(6, '(3X, I0, X, A)') i, section%cells%elements(i)%to_string()
		end do
		_log_write(6, '(A)') ""
#	endif
end subroutine

subroutine leaf(section_traversal, thread, section, current_element)
	type(t_section_traversal), intent(inout)	        :: section_traversal
	type(t_grid_thread), intent(inout)					:: thread
	type(t_grid_section), intent(inout)					:: section
	type(t_traversal_element), intent(inout)	        :: current_element

	call read(section_traversal, thread, section, current_element)

#	if defined(_GT_ELEMENT_OP)
		call _GT_ELEMENT_OP(section_traversal, section, current_element%t_element_base)
#	endif

	call write(section_traversal, thread, section, current_element)
end subroutine

subroutine old_leaf(section_traversal, thread, section, current_element)
	type(t_section_traversal), intent(inout)	        :: section_traversal
	type(t_grid_thread), intent(inout)					:: thread
	type(t_grid_section), intent(inout)					:: section
	type(t_traversal_element), intent(inout)	        :: current_element

	call read_oon(section_traversal, thread, section, current_element)

#	if defined(_GT_INNER_ELEMENT_OP)
		call _GT_INNER_ELEMENT_OP(section_traversal, section, current_element%t_element_base)
#	endif

	call write_oon(section_traversal, thread, section, current_element)
end subroutine

subroutine new_leaf(section_traversal, thread, section, current_element)
	type(t_section_traversal), intent(inout)	        :: section_traversal
	type(t_grid_thread), intent(inout)					:: thread
	type(t_grid_section), intent(inout)					:: section
	type(t_traversal_element), intent(inout)	        :: current_element

	call read_onn(section_traversal, thread, section, current_element)

#	if defined(_GT_INNER_ELEMENT_OP)
		call _GT_INNER_ELEMENT_OP(section_traversal, section, current_element%t_element_base)
#	endif

	call write_onn(section_traversal, thread, section, current_element)
end subroutine

subroutine old_bnd_leaf(section_traversal, thread, section, current_element)
	type(t_section_traversal), intent(inout)	        :: section_traversal
	type(t_grid_thread), intent(inout)					:: thread
	type(t_grid_section), intent(inout)					:: section
	type(t_traversal_element), intent(inout)	        :: current_element

	call read_odn(section_traversal, thread, section, current_element)

#	if defined(_GT_ELEMENT_OP)
		call _GT_ELEMENT_OP(section_traversal, section, current_element%t_element_base)
#	endif

	call write_odn(section_traversal, thread, section, current_element)
end subroutine

subroutine new_bnd_leaf(section_traversal, thread, section, current_element)
	type(t_section_traversal), intent(inout)	        :: section_traversal
	type(t_grid_thread), intent(inout)					:: thread
	type(t_grid_section), intent(inout)					:: section
	type(t_traversal_element), intent(inout)	        :: current_element

	call read_obn(section_traversal, thread, section, current_element)

#	if defined(_GT_ELEMENT_OP)
		call _GT_ELEMENT_OP(section_traversal, section, current_element%t_element_base)
#	endif

	call write_obn(section_traversal, thread, section, current_element)
end subroutine

subroutine create_ringbuffer(elements)
	type(t_traversal_element), dimension(:), target, intent(inout)		:: elements

	integer (kind = GRID_SI)											:: i

	do i = 1, size(elements)
		elements(i)%previous => elements(mod(i + size(elements) - 2, size(elements)) + 1)
		elements(i)%next => elements(mod(i, size(elements)) + 1)

		nullify(elements(i)%cell%geometry)
		nullify(elements(i)%cell%data_pers)
		elements(i)%cell%data_temp => elements(i)%cell_data_temp

		nullify(elements(i)%color_node_out%ptr)
		nullify(elements(i)%transfer_node%ptr)
		nullify(elements(i)%color_node_in%ptr)

#		if defined(_GT_EDGES)
            nullify(elements(i)%color_edge%ptr)
            elements(i)%previous_edge%ptr => elements(i)%previous%next_edge_data
            elements(i)%next_edge%ptr => elements(i)%next_edge_data
#		endif
	end do
end subroutine

#define _GT_INPUT
#define _GT_OUTPUT

#include "Tools_adaptive_traversal.f90"

#undef _GT_INPUT
#undef _GT_OUTPUT

#undef _GT_CELL_FIRST_TOUCH_OP
#undef _GT_CELL_LAST_TOUCH_OP
#undef _GT_CELL_REDUCE_OP
#undef _GT_EDGE_FIRST_TOUCH_OP
#undef _GT_EDGE_LAST_TOUCH_OP
#undef _GT_EDGE_REDUCE_OP
#undef _GT_EDGE_MPI_TYPE
#undef _GT_EDGE_MERGE_OP
#undef _GT_EDGE_WRITE_OP
#undef _GT_NODE_FIRST_TOUCH_OP
#undef _GT_NODE_LAST_TOUCH_OP
#undef _GT_NODE_REDUCE_OP
#undef _GT_NODE_MPI_TYPE
#undef _GT_NODE_MERGE_OP
#undef _GT_NODE_WRITE_OP
#undef _GT_INNER_EDGE_FIRST_TOUCH_OP
#undef _GT_INNER_EDGE_LAST_TOUCH_OP
#undef _GT_INNER_EDGE_REDUCE_OP
#undef _GT_INNER_NODE_FIRST_TOUCH_OP
#undef _GT_INNER_NODE_LAST_TOUCH_OP
#undef _GT_INNER_NODE_REDUCE_OP
#undef _GT_ELEMENT_OP
#undef _GT_INNER_ELEMENT_OP
#undef _GT_PRE_TRAVERSAL_OP
#undef _GT_POST_TRAVERSAL_OP
#undef _GT_PRE_TRAVERSAL_GRID_OP
#undef _GT_POST_TRAVERSAL_GRID_OP

#undef _GT_NODES
#undef _GT_EDGES
#undef _GT_NO_COORDS
#undef _GT_NAME

#undef t_section_traversal

#undef _GT_SKIP_LB_TIMING
