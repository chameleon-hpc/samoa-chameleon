#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

! This module provides the XDMF core API for the input traversal
module XDMF_initialize_dofs_base
    
    use SFC_data_types
    use Tools_openmp
    use Tools_mpi

    use XDMF_data_types
    use XDMF_initialize_dofs_base_data_types
    use XDMF_math
    use XDMF_hdf5

    use, intrinsic :: iso_fortran_env

    implicit none

    contains

#   if defined (_XDMF_PATCH)    
        subroutine xdmf_base_pre_traversal_grid_op(base, sections_ptr, grid, param_cells, param_patches)
#   else
        subroutine xdmf_base_pre_traversal_grid_op(base, sections_ptr, grid, param_cells)
#   endif
        type(t_xdmf_base_initialize_dofs_traversal), intent(inout)	:: base
        type(t_xdmf_base_initialize_dofs_traversal_ptr), dimension(:)   :: sections_ptr
        type(t_grid), intent(inout)					:: grid
        type(t_xdmf_parameter), intent(in)                              :: param_cells
#       if defined (_XDMF_PATCH)
            type(t_xdmf_parameter), intent(in)                          :: param_patches
#       endif

        integer                                                         :: result, i, hdf5_error
        integer(GRID_SI)                                                :: output_meta_iteration
        character(len = 256)					        :: file_name_h5
        logical                                                         :: hdf5_step_exists
        integer(HSIZE_T), dimension(hdf5_rank)                          :: hdf5_tree_dims, hdf5_tree_maxdims

        base%grid_scale = ishft(1_XDMF_GRID_SI, ishft(cfg%i_max_depth + 1, -1))
        output_meta_iteration = int(real(base%output_iteration, XDMF_GRID_SR) / cfg%xdmf%i_xdmfspf, XDMF_GRID_SI)
      
        ! Open file and init meta hdf5 ids
        write (file_name_h5, "(A, A, I0, A, A)") trim(base%s_file_stamp), "_", output_meta_iteration, "_xdmf.h5", char(0)
        call base%root_desc%hdf5_meta_ids%createopen(file_name_h5, base%output_iteration, .true., .false., .false., .false., hdf5_step_exists, result)
        if (result.ne.0) then
            _log_write(1, '(A, I0, A)') " XDMF: Error: Cannot open file: ", trim(base%s_file_stamp)
            assert(.false.)
            call exit()
        end if
        if (.not.hdf5_step_exists) then
            _log_write(1, '(A, I0)') " XDMF: Error: Step does not exist: ", base%output_iteration
        end if

        ! Open datasets
        call base%root_desc%hdf5_ids_cells%open(param_cells, base%root_desc%hdf5_meta_ids%step_group_cells_id)
#       if defined(_XDMF_PATCH)
            call base%root_desc%hdf5_ids_patches%open(param_patches, base%root_desc%hdf5_meta_ids%step_group_patches_id)
#       endif
        ! Read tree table size
        call h5sget_simple_extent_dims_f(base%root_desc%hdf5_ids_cells%tree%dspace_id, hdf5_tree_dims, hdf5_tree_maxdims, hdf5_error)
        base%htbl_size = hdf5_tree_dims(2)
        ! Compute hashtable secondary prime
        call close_prime(base%htbl_size - 1, .false., base%hash2_prime)

        ! Scatter computated data across all sections
        do i = 1, size(sections_ptr)
            call base%root_desc%hdf5_meta_ids%scatter_to(sections_ptr(i)%ptr%root_desc%hdf5_meta_ids)
            call sections_ptr(i)%ptr%root_desc%hdf5_ids_cells%allocate(param_cells)
            call base%root_desc%hdf5_ids_cells%scatter_to(param_cells, sections_ptr(i)%ptr%root_desc%hdf5_ids_cells)
#           if defined(_XDMF_PATCH)
                call sections_ptr(i)%ptr%root_desc%hdf5_ids_patches%allocate(param_patches)
                call base%root_desc%hdf5_ids_patches%scatter_to(param_patches, sections_ptr(i)%ptr%root_desc%hdf5_ids_patches)
#           endif
            sections_ptr(i)%ptr%s_file_stamp = base%s_file_stamp
            sections_ptr(i)%ptr%grid_scale = base%grid_scale
            sections_ptr(i)%ptr%output_iteration = base%output_iteration
            sections_ptr(i)%ptr%l_load_data = base%l_load_data
            sections_ptr(i)%ptr%htbl_size = base%htbl_size
            sections_ptr(i)%ptr%hash2_prime = base%hash2_prime
        end do
    end subroutine

#   if defined (_XDMF_PATCH)   
        subroutine xdmf_base_post_traversal_grid_op(base, sections_ptr, grid, param_cells, param_patches)
#   else
        subroutine xdmf_base_post_traversal_grid_op(base, sections_ptr, grid, param_cells)
#   endif
        type(t_xdmf_base_initialize_dofs_traversal), intent(inout)	:: base
        type(t_xdmf_base_initialize_dofs_traversal_ptr), dimension(:)   :: sections_ptr
        type(t_grid), intent(inout)					:: grid
        type(t_xdmf_parameter), intent(in)                              :: param_cells
#       if defined (_XDMF_PATCH)
            type(t_xdmf_parameter), intent(in)                          :: param_patches
#       endif

        integer                                                         :: i, error
        integer(INT32)                                                  :: refinements_issued_local

        ! Close hdf5 ids
        call base%root_desc%hdf5_ids_cells%close(param_cells)
#       if defined(_XDMF_PATCH)
            call base%root_desc%hdf5_ids_patches%close(param_patches)
#       endif
        call base%root_desc%hdf5_meta_ids%close()

        ! Collect number of refinements issued
        refinements_issued_local = 0
        do i = 1, size(sections_ptr)
            call sections_ptr(i)%ptr%root_desc%hdf5_ids_cells%deallocate()
#       if defined(_XDMF_PATCH)
            call sections_ptr(i)%ptr%root_desc%hdf5_ids_patches%deallocate()
#       endif
            refinements_issued_local = refinements_issued_local + sections_ptr(i)%ptr%i_refinements_issued
        end do
#       if defined(_MPI)
            call mpi_allreduce(refinements_issued_local, base%i_refinements_issued, 1, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, error)
#       else
            base%i_refinements_issued = refinements_issued_local
#       endif
    end subroutine

end module

