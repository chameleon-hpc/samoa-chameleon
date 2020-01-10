#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

! This module provides common data types used by the XDMF routines
module XDMF_data_types

    use HDF5

    use, intrinsic :: iso_fortran_env

    implicit none

    ! Needed to prevent module conflicts
#   if defined(_SINGLE_PRECISION)
        integer, PARAMETER :: XDMF_GRID_SR = kind(1.0e0)
#   elif defined(_DOUBLE_PRECISION)
        integer, PARAMETER :: XDMF_GRID_SR = kind(1.0d0)
#   elif defined(_QUAD_PRECISION)
        integer, PARAMETER :: XDMF_GRID_SR = kind(1.0q0)
#    else
#       error "No floating point precision is chosen!"
#    endif

    integer, PARAMETER :: XDMF_GRID_SI = selected_int_kind(8)
    integer, PARAMETER :: XDMF_GRID_DI = selected_int_kind(16)

    ! Sizes of the datasets
    integer, parameter			         	:: hdf5_rank = 2, hdf5_rank_v = 3
    integer(HSIZE_T), parameter				:: hdf5_tree_width = 2, hdf5_sect_tree_width = 2, hdf5_vector_width = 2
    integer(HSIZE_T), parameter, &
        dimension(hdf5_rank)			   	:: hdf5_tree_chunk_size = (/ hdf5_tree_width, 8_HSIZE_T /)
    integer(HSIZE_T), parameter             :: hdf5_vals_chunk_length = 32_HSIZE_T
    integer(HSIZE_T), parameter, &
        dimension(hdf5_rank)			   	:: hdf5_subset_stride = (/ 1, 1 /)
    integer(HSIZE_T), parameter, &
        dimension(hdf5_rank)			   	:: hdf5_subset_block = (/ 1, 1 /)
    integer(HSIZE_T), parameter, &
        dimension(hdf5_rank_v)			   	:: hdf5_subset_stride_v = (/ 1, 1, 1 /)
    integer(HSIZE_T), parameter, &
        dimension(hdf5_rank_v)			   	:: hdf5_subset_block_v = (/ 1, 1, 1 /)

    ! Names of the datasets
    character(len = 1), parameter			:: hdf5_attr_dname_nz = "a"
    character(len = 2), parameter			:: hdf5_attr_dname = hdf5_attr_dname_nz//char(0)
    character(len = 2), parameter			:: hdf5_tree_dname = "t"//char(0)
    character(len = 1), parameter			:: hdf5_valsg_dname_nz = "g"
    character(len = 1), parameter			:: hdf5_valst_dname_nz = "p"
    character(len = 2), parameter			:: hdf5_valsg_dname = hdf5_valsg_dname_nz//char(0)
    character(len = 2), parameter			:: hdf5_valst_dname = hdf5_valst_dname_nz//char(0)

    ! This parameter controls the amount of data for a cell
    type t_xdmf_parameter
        integer(HSIZE_T)                    :: hdf5_valsg_width     !< Amount of geometry data associated with each geometry data point
        integer(HSIZE_T)                    :: hdf5_valst_width     !< Amount of geometry data points per cell
        integer(HSIZE_T)                    :: hdf5_attr_width      !< Amount of attribute values per cell
        integer(HSIZE_T)                    :: hdf5_valsi_width     !< Amount of int32 values per cell
        character(len = 1), dimension(:), allocatable :: hdf5_valsi_dnames !< Names of the int32 values
        integer(HSIZE_T)                    :: hdf5_valsuv_width   !< Amount of 2D vector real32 values per attribute
        character(len = 1), dimension(:), allocatable :: hdf5_valsuv_dnames !< Names of the 2D vector real32 values
        integer(HSIZE_T)                    :: hdf5_valsr_width     !< Amount of real32 values per attribute
        character(len = 1), dimension(:), allocatable :: hdf5_valsr_dnames !< Names of the real32 values

        contains

        procedure, pass                     :: allocate => xdmf_parameter_allocate
        procedure, pass                     :: deallocate => xdmf_parameter_deallocate
    end type                                         

    ! This structure describes the location and space of a section in a HDF5 file
    type t_xdmf_section_descriptor
        integer(XDMF_GRID_SI)			    :: num_cells = 0, offset_cells = 0, offset_cells_buffer = 0
    end type

    ! This structure describes the layout of all sections of a rank in a HDF5 file
    type t_xdmf_section_descriptor_collection
        type(t_xdmf_section_descriptor), &
            dimension(:), allocatable       :: sections
        integer(XDMF_GRID_SI)               :: num_cells = 0, offset_cells = 0

        contains

        procedure, pass                     :: allocate => xdmf_file_descriptor_collection_allocate
        procedure, pass                     :: deallocate => xdmf_file_descriptor_collection_deallocate
    end type

    ! This structure describes the layout of all sections in a HDF5 file
    type t_xdmf_layout_descriptor
        type(t_xdmf_section_descriptor_collection), &
            dimension(:), allocatable       :: ranks

        contains

        procedure, pass                     :: allocate => xdmf_layout_descriptor_allocate
        procedure, pass                     :: deallocate => xdmf_layout_descriptor_deallocate
        procedure, pass                     :: scatter_to => xdmf_layout_descriptor_scatter_to
    end type

    contains

    subroutine xdmf_parameter_allocate(this)
        class(t_xdmf_parameter), &
            intent(inout)                   :: this

        integer                             :: error

        allocate(this%hdf5_valsi_dnames(this%hdf5_valsi_width), stat = error); assert_eq(error, 0)
        allocate(this%hdf5_valsr_dnames(this%hdf5_valsr_width), stat = error); assert_eq(error, 0)
        allocate(this%hdf5_valsuv_dnames(this%hdf5_valsuv_width), stat = error); assert_eq(error, 0)
    end subroutine

    subroutine xdmf_parameter_deallocate(this)
        class(t_xdmf_parameter), &
            intent(inout)                   :: this

        integer                             :: error

        deallocate(this%hdf5_valsi_dnames, stat = error); assert_eq(error, 0)
        deallocate(this%hdf5_valsr_dnames, stat = error); assert_eq(error, 0)
        deallocate(this%hdf5_valsuv_dnames, stat = error); assert_eq(error, 0)
    end subroutine

    subroutine xdmf_file_descriptor_collection_allocate(this, num_sections)
        class(t_xdmf_section_descriptor_collection), &
            intent(inout)                   :: this
        integer, intent(in)                 :: num_sections

        integer                             :: error

        allocate(this%sections(num_sections), stat = error); assert_eq(error, 0)
    end subroutine

    subroutine xdmf_file_descriptor_collection_deallocate(this)
        class(t_xdmf_section_descriptor_collection), &
            intent(inout)                    :: this

        integer                              :: error

        deallocate(this%sections, stat = error); assert_eq(error, 0)
    end subroutine

    subroutine xdmf_layout_descriptor_scatter_to(this, child)
        class(t_xdmf_layout_descriptor), &
            intent(inout)                    :: this
        class(t_xdmf_layout_descriptor), &
            intent(inout)                    :: child

        integer                              :: n, k, error

        allocate(child%ranks(size(this%ranks)), stat = error); assert_eq(error, 0)
        do n = 1, size(this%ranks)
            child%ranks(n)%num_cells = this%ranks(n)%num_cells
            child%ranks(n)%offset_cells = this%ranks(n)%offset_cells
            allocate(child%ranks(n)%sections(size(this%ranks(n)%sections)), stat = error); assert_eq(error, 0)
            do k = 1, size(this%ranks(n)%sections)
                child%ranks(n)%sections(k)%num_cells = this%ranks(n)%sections(k)%num_cells
                child%ranks(n)%sections(k)%offset_cells = this%ranks(n)%sections(k)%offset_cells
                child%ranks(n)%sections(k)%offset_cells_buffer = this%ranks(n)%sections(k)%offset_cells_buffer
            end do
        end do
    end subroutine

    subroutine xdmf_layout_descriptor_allocate(this, num_ranks)
        class(t_xdmf_layout_descriptor), &
            intent(inout)                   :: this
        integer, intent(in)                 :: num_ranks

        integer                             :: error

        allocate(this%ranks(num_ranks), stat = error); assert_eq(error, 0)
    end subroutine

    subroutine xdmf_layout_descriptor_deallocate(this)
        class(t_xdmf_layout_descriptor), &
            intent(inout)                   :: this

        integer                             :: n, error

        do n = 1, size(this%ranks)
            call this%ranks(n)%deallocate()
        end do
        deallocate(this%ranks, stat = error); assert_eq(error, 0)
    end subroutine

end module