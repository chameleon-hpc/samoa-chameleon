#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

! This module defines the XDMF input traversal
#if defined(_SWE)
    module SWE_XDMF_Initialize_Dofs
        use Tools_noise

        use iso_c_binding
        use SFC_edge_traversal
        use SWE_euler_timestep
        use Samoa_swe

        use SWE_XDMF_config
        use XDMF_initialize_dofs_base

#       if defined(_SWE_PATCH)
            use SWE_PATCH
#       endif
#       if defined(_SWE_DG)
            use SWE_dg_matrices
            use SWE_data_types
#       endif   

        implicit none

        ! For traversal data, use the XDMF core data structure as a base
        type num_traversal_data
            type(t_xdmf_base_initialize_dofs_traversal) :: base
        end type

#       define _GT_EDGES
#		define _GT_NAME							t_swe_xdmf_init_dofs_traversal

#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op
#		define _GT_ELEMENT_OP					element_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

        ! This routine creates an array of wrappedsection pointers, 
        ! because fortran does not like arrays of pointers.
        subroutine ptr_wrap_sections(traversal, sections_ptr)
            type(t_swe_xdmf_init_dofs_traversal), intent(inout)			:: traversal
            type(t_xdmf_base_initialize_dofs_traversal_ptr), &
                dimension(:), allocatable, intent(out)                      :: sections_ptr

            integer                                                         :: error, i

            allocate(sections_ptr(size(traversal%sections)), stat = error); assert_eq(error, 0)
            do i = 1, size(traversal%sections)
                sections_ptr(i)%ptr => traversal%sections(i)%base
            end do
        end subroutine

        subroutine pre_traversal_grid_op(traversal, grid)
            type(t_swe_xdmf_init_dofs_traversal), intent(inout)		    :: traversal
            type(t_grid), intent(inout)							    		:: grid

            type(t_xdmf_base_initialize_dofs_traversal_ptr), &
                dimension(:), allocatable                                   :: sections_ptr
            integer                                                         :: error
            
            ! Pass this call the the core API
            call ptr_wrap_sections(traversal, sections_ptr)
            call xdmf_base_pre_traversal_grid_op(traversal%base, sections_ptr, grid, swe_xdmf_param)
            deallocate(sections_ptr, stat = error); assert_eq(error, 0)
        end subroutine

        subroutine post_traversal_grid_op(traversal, grid)
            type(t_swe_xdmf_init_dofs_traversal), intent(inout)		    :: traversal
            type(t_grid), intent(inout)							    		:: grid

            type(t_xdmf_base_initialize_dofs_traversal_ptr), &
                dimension(:), allocatable                                   :: sections_ptr
            integer                                                         :: error
            
            ! Pass this call the the core API
            call ptr_wrap_sections(traversal, sections_ptr)
            call xdmf_base_post_traversal_grid_op(traversal%base, sections_ptr, grid, swe_xdmf_param)
            deallocate(sections_ptr, stat = error); assert_eq(error, 0)
        end subroutine

        subroutine pre_traversal_op(traversal, section)
            type(t_swe_xdmf_init_dofs_traversal), intent(inout)			:: traversal
            type(t_grid_section), intent(inout)								:: section

            traversal%base%i_refinements_issued = 0
        end subroutine

        subroutine element_op(traversal, section, element)
            type(t_swe_xdmf_init_dofs_traversal), intent(inout)			:: traversal
            type(t_grid_section), intent(inout)								:: section
            type(t_element_base), intent(inout)					        	:: element

            integer															:: hdf5_error, i, point_id
            integer(INT64)                                                  :: element_hash, htbl_key, htbl_try
            integer(HID_T)													:: memspace_id
            integer (HSIZE_T), dimension(hdf5_rank)                    		:: hdf5_tree_offset, hdf5_tree_dims = (/ 2_HSIZE_T, 1_HSIZE_T /)
            integer (INT32), dimension(hdf5_tree_width)                     :: hdf5_tree_buffer

            type(t_state), dimension(_SWE_CELL_SIZE)            		    :: Q

            real (XDMF_ISO_P), dimension(swe_xdmf_param%hdf5_attr_width, 1)   :: hdf5_bh_buffer, hdf5_b_buffer
            integer(INT32), dimension(1) 									:: hdf5_l_buffer
            real (XDMF_ISO_P), dimension(hdf5_vector_width, swe_xdmf_param%hdf5_attr_width, 1) :: hdf5_f_buffer, hdf5_g_buffer

            ! Assume no refinement per default
            element%cell%geometry%refinement = 0

            ! Compute this cells hash
            call xdmf_hash_element(element, int(traversal%base%grid_scale, INT64), element_hash)

            ! Do the lookup in the tree hash table
            htbl_try = 0
            do
                ! Compute hash
                call xdmf_hashtable_hash(element_hash, traversal%base%htbl_size, traversal%base%hash2_prime, htbl_try, htbl_key)
                hdf5_tree_offset = (/ 0_HSIZE_T, htbl_key /)
                
                ! Lookup in HDF5 file
                ! This is marked critical, because HDF5 manages an internal state which is not thread-safe.
                ! Even if it was (using threadsafe-hdf5), thread-local variables like memspace_id would become inconsistent in the HDF5 state
                !$omp critical
                call h5screate_simple_f(hdf5_rank, hdf5_tree_dims, memspace_id, hdf5_error)
                call h5sselect_hyperslab_f(traversal%base%root_desc%hdf5_ids%tree%dspace_id, H5S_SELECT_SET_F, &
                    hdf5_tree_offset, hdf5_tree_dims, hdf5_error, hdf5_subset_stride, hdf5_subset_block)
                call h5dread_f(traversal%base%root_desc%hdf5_ids%tree%dset_id, H5T_NATIVE_INTEGER, hdf5_tree_buffer, hdf5_tree_dims, &
                    hdf5_error, memspace_id, traversal%base%root_desc%hdf5_ids%tree%dspace_id, &
                    xfer_prp = traversal%base%root_desc%hdf5_meta_ids%access_dset_id)
                call h5sclose_f(memspace_id, hdf5_error)
                !$omp end critical

                ! Cell found or this cell does not exist
                if ((hdf5_tree_buffer(1) - 1 .eq. element_hash) .or. &
                    (hdf5_tree_buffer(1) .eq. 0)) then
                    exit
                else
                    ! Rehash
                    htbl_try = htbl_try + 1
                end if
            end do
            
            ! Actual cell value data is only loaded in the very last traversal
            if(.not.traversal%base%l_load_data) then
                if(hdf5_tree_buffer(1) .eq. 0) then
                    element%cell%geometry%refinement = 1
                    traversal%base%i_refinements_issued = traversal%base%i_refinements_issued + 1
                end if
            else
                if(hdf5_tree_buffer(1) .eq. 0) then
                    _log_write(1, *) "XDMF: ERROR: Element data not found"
                    assert(.false.)
                end if

                hdf5_tree_buffer(2) = hdf5_tree_buffer(2) - 1
                
                ! TODO XDMF read values
                ! Load actual cell values
                ! !$omp critical
                ! call hdf5_read_chunk_real(traversal%base%root_desc%hdf5_ids%batch_valsr(swe_hdf5_valsr_b_offset)%dset_id, &
                !     traversal%base%root_desc%hdf5_ids%batch_valsr(swe_hdf5_valsr_b_offset)%dspace_id, &
                !     (/ 0_HSIZE_T, int(hdf5_tree_buffer(2), HSIZE_T) /), (/ swe_xdmf_param%hdf5_attr_width, 1_HSIZE_T /), &
                !     hdf5_rank, hdf5_subset_stride, hdf5_subset_block, &
                !     hdf5_b_buffer, traversal%base%root_desc%hdf5_meta_ids%access_dset_id)
                ! call hdf5_read_chunk_real(traversal%base%root_desc%hdf5_ids%batch_valsr(swe_hdf5_valsr_bh_offset)%dset_id, &
                !     traversal%base%root_desc%hdf5_ids%batch_valsr(swe_hdf5_valsr_bh_offset)%dspace_id, &
                !     (/ 0_HSIZE_T, int(hdf5_tree_buffer(2), HSIZE_T) /), (/ swe_xdmf_param%hdf5_attr_width, 1_HSIZE_T /), &
                !     hdf5_rank, hdf5_subset_stride, hdf5_subset_block, &
                !     hdf5_bh_buffer, traversal%base%root_desc%hdf5_meta_ids%access_dset_id)
                ! call hdf5_read_chunk_real(traversal%base%root_desc%hdf5_ids%batch_valsuv(swe_hdf5_valsuv_f_offset)%dset_id, &
                !     traversal%base%root_desc%hdf5_ids%batch_valsuv(swe_hdf5_valsuv_f_offset)%dspace_id, &
                !     (/ 0_HSIZE_T, 0_HSIZE_T, int(hdf5_tree_buffer(2), HSIZE_T) /), (/ hdf5_vector_width, &
                !     swe_xdmf_param%hdf5_attr_width, 1_HSIZE_T /), hdf5_rank_v, hdf5_subset_stride_v, hdf5_subset_block_v, &
                !     hdf5_f_buffer, traversal%base%root_desc%hdf5_meta_ids%access_dset_id)
                ! call hdf5_read_chunk_real(traversal%base%root_desc%hdf5_ids%batch_valsuv(swe_hdf5_valsuv_g_offset)%dset_id, &
                !     traversal%base%root_desc%hdf5_ids%batch_valsuv(swe_hdf5_valsuv_g_offset)%dspace_id, &
                !     (/ 0_HSIZE_T, 0_HSIZE_T, int(hdf5_tree_buffer(2), HSIZE_T) /), (/ hdf5_vector_width, &
                !     swe_xdmf_param%hdf5_attr_width, 1_HSIZE_T /), hdf5_rank_v, hdf5_subset_stride_v, hdf5_subset_block_v, &
                !     hdf5_g_buffer, traversal%base%root_desc%hdf5_meta_ids%access_dset_id)
                ! call hdf5_read_chunk_int(traversal%base%root_desc%hdf5_ids%batch_valsi(swe_hdf5_valsi_plotter_offset)%dset_id, &
                !     traversal%base%root_desc%hdf5_ids%batch_valsi(swe_hdf5_valsi_plotter_offset)%dspace_id, &
                !     (/ 0_HSIZE_T, int(hdf5_tree_buffer(2), HSIZE_T) /), (/ 1_HSIZE_T, 1_HSIZE_T /), &
                !     hdf5_rank, hdf5_subset_stride, hdf5_subset_block, &
                !     hdf5_l_buffer, traversal%base%root_desc%hdf5_meta_ids%access_dset_id)
                ! !$omp end critical

                !Write data into samoa data model
                do i = 1, swe_xdmf_param%hdf5_attr_width
                    if (hdf5_l_buffer(1) .le. 0) then
                        point_id = i
                        if (i .eq. 1) then
                            point_id = 3
                        else if (i .eq. 2) then
                            point_id = 1
                        else if (i .eq. 3) then
                            point_id = 2
                        end if
                    else
                        point_id = i
                        if (i .eq. 1) then
                            point_id = 2
                        else if (i .eq. 2) then
                            point_id = 1
                        end if
                    end if

                    ! TODO XDMF write data into cell model
                    ! element%cell%data_pers%Q(point_id)%b = hdf5_b_buffer(i, 1)
                    ! element%cell%data_pers%Q(point_id)%h = hdf5_bh_buffer(i, 1)
                    ! element%cell%data_pers%Q(point_id)%p(:) = hdf5_f_buffer(:, i, 1)
                    ! element%cell%data_pers%Q(point_id)%gradB(:) = hdf5_g_buffer(:, i, 1)
                    ! element%cell%data_pers%Q_inter(point_id)%b = element%cell%data_pers%Q(point_id)%b
                    ! element%cell%data_pers%Q_inter(point_id)%h = element%cell%data_pers%Q(point_id)%h
                    ! element%cell%data_pers%Q_inter(point_id)%p(:) = element%cell%data_pers%Q(point_id)%p(:)
                    ! element%cell%data_pers%Q_inter(point_id)%gradB(:) = element%cell%data_pers%Q(point_id)%gradB(:)
                end do
            end if
        end subroutine

    end module
#endif