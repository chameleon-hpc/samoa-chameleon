#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

! This module defines the traversal data structure for the XDMF output traversals
module XDMF_output_base_data_types
    
    use SFC_edge_traversal
    use Tools_openmp

    use XDMF_data_types
    use XDMF_hdf5

    implicit none

    ! Per-section buffers
    type t_xdmf_section_buffer
        integer (INT32), &
            dimension(:), allocatable           :: tree
        integer (INT32), &
            dimension(:, :), allocatable        :: valsi
        real (XDMF_ISO_P), &
            dimension(:, :, :), allocatable     :: valsr
        real (XDMF_ISO_P), &
            dimension(:, :, :, :), allocatable  :: valsuv
        real (XDMF_ISO_P), &
            dimension(:, :), allocatable        :: valsg

        contains 

        procedure, pass                         :: allocate => xdmf_section_buffer_allocate
        procedure, pass                         :: deallocate => xdmf_section_buffer_deallocate
    end type

    type t_xdmf_section_buffer_ptr
        type(t_xdmf_section_buffer), &
            pointer                             :: ptr
    end type

    ! The base traversal data structure
    type t_xdmf_base_output_filter_traversal
        integer (XDMF_GRID_SI)	                :: i_sim_iteration = 0, output_iteration = 0
    end type
    type, extends(t_xdmf_base_output_filter_traversal) :: t_xdmf_base_output_traversal
        integer (XDMF_GRID_SI)	                :: num_cells = 0
        integer (XDMF_GRID_SI)	                :: grid_scale
        integer (XDMF_GRID_SI)	                :: sect_store_index = 1
        character(len = 256)					:: s_file_stamp, s_file_stamp_base
        integer                                 :: xdmf_remove_lines = 0

        type(t_xdmf_layout_descriptor)          :: root_layout_desc
        type(t_xdmf_file_descriptor)            :: root_desc
        type(t_xdmf_section_buffer_ptr)         :: sect_store
    end type

    ! This pointer data structure is needed to pass the section arrays to the core XDMF library
    type t_xdmf_base_output_filter_traversal_ptr
        type(t_xdmf_base_output_filter_traversal), pointer :: ptr
    end type
    type t_xdmf_base_output_traversal_ptr
        type(t_xdmf_base_output_traversal), pointer :: ptr
    end type

    contains

    ! This routine allocates memory for the cell attribute buffer in a section
    subroutine xdmf_section_buffer_allocate(this, sect_cells, param)
        class(t_xdmf_section_buffer), &
            intent(inout)                       :: this
        integer, intent(in)                     :: sect_cells        
        type(t_xdmf_parameter), intent(in)      :: param

        integer(XDMF_GRID_SI)                   :: sect_cells_actual
        integer                                 :: error

#       if defined (_XDMF_PATCH)
            sect_cells_actual = sect_cells * _XDMF_PATCH_ORDER_SQUARE
#       else
            sect_cells_actual = sect_cells
#       endif

        allocate(this%tree(sect_cells), stat = error); assert_eq(error, 0)
        this%tree = 0
        allocate(this%valsi(sect_cells_actual, param%hdf5_valsi_width), stat = error); assert_eq(error, 0)
        this%valsi = 0
        allocate(this%valsr(param%hdf5_attr_width, sect_cells_actual, param%hdf5_valsr_width), stat = error); assert_eq(error, 0)
        this%valsr = 0
        allocate(this%valsuv(hdf5_vector_width, param%hdf5_attr_width, sect_cells_actual, param%hdf5_valsuv_width), stat = error); assert_eq(error, 0)
        this%valsuv = 0
        allocate(this%valsg(param%hdf5_valsg_width, sect_cells_actual * param%hdf5_valst_width), stat = error); assert_eq(error, 0)
        this%valsg = 0
    end subroutine

     ! This routine frees the cell attribute buffers of a section
    subroutine xdmf_section_buffer_deallocate(this)
        class(t_xdmf_section_buffer), &
            intent(inout)                   :: this

        integer                             :: error

        deallocate(this%tree, stat = error); assert_eq(error, 0)
        deallocate(this%valsi, stat = error); assert_eq(error, 0)
        deallocate(this%valsr, stat = error); assert_eq(error, 0)
        deallocate(this%valsuv, stat = error); assert_eq(error, 0)
        deallocate(this%valsg, stat = error); assert_eq(error, 0)
    end subroutine


end module
