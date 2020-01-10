#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

! This module defines the traversal data structure for the XDMF input traversals
module XDMF_initialize_dofs_base_data_types

    use SFC_edge_traversal
    use Tools_openmp

    use XDMF_data_types
    use XDMF_hdf5

    use, intrinsic :: iso_fortran_env

    implicit none

    ! The base traversal data structure
    type t_xdmf_base_initialize_dofs_traversal
        integer (XDMF_GRID_SI)	                :: grid_scale, output_iteration
        character(len = 256)					:: s_file_stamp
        logical								    :: l_load_data

        type(t_xdmf_file_descriptor)            :: root_desc

        integer (XDMF_GRID_DI)			        :: i_refinements_issued = 0
        integer (INT64)                         :: htbl_size, hash2_prime
    end type

    ! This pointer data structure is needed to pass the section arrays to the core XDMF library
    type t_xdmf_base_initialize_dofs_traversal_ptr
        type(t_xdmf_base_initialize_dofs_traversal), pointer :: ptr
    end type

    contains

end module