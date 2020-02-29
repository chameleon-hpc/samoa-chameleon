#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

! This odule provides data types and routines to ease the usage of HDF5
module XDMF_hdf5

    use XDMF_data_types
    use Tools_mpi
    use HDF5

    use, intrinsic :: iso_fortran_env

    implicit none

    ! This structure stores the HDF5 ids of a dataset and its dataspace
    type t_xdmf_hdf5_idpair
        integer(HID_T)			         	:: dspace_id = 0, dset_id = 0
        logical                             :: is_open = .false.

        contains

        procedure, pass                     :: create => xdmf_hdf5_idpair_create
        procedure, pass                     :: open => xdmf_hdf5_idpair_open
        procedure, pass                     :: close => xdmf_hdf5_idpair_close
        procedure, pass                     :: scatter_to => xdmf_hdf5_idpair_scatter_to
    end type
                     
    ! This structure stores the HDF5 ids of a file
    type t_xdmf_hdf5_metadata 
        integer(HID_T)			         	:: file_id = 0, step_group_id = 0, step_group_cells_id = 0, access_dset_id = 0
#       if defined(_XDMF_PATCH)
            integer(HID_T)			        :: step_group_patches_id = 0
#       endif

        contains

        procedure, pass                     :: createopen => xdmf_hdf5_metadata_createopen
        procedure, pass                     :: close => xdmf_hdf5_metadata_close
        procedure, pass                     :: scatter_to => xdmf_hdf5_metadata_scatter_to
    end type

    ! This structure stores the HDF5 ids of all datasets of a step
    type t_xdmf_hdf5
        integer(HID_T)                      :: attr_group_id = 0
        type(t_xdmf_hdf5_idpair)         	:: tree, valsg, valst
        type(t_xdmf_hdf5_idpair), &
            dimension(:), allocatable       :: batch_valsi
        type(t_xdmf_hdf5_idpair), &
            dimension(:), allocatable       :: batch_valsr
        type(t_xdmf_hdf5_idpair), &
            dimension(:), allocatable       :: batch_valsuv

        contains

        procedure, pass                     :: allocate => xdmf_hdf5_allocate
        procedure, pass                     :: deallocate => xdmf_hdf5_deallocate
        procedure, pass                     :: create => xdmf_hdf5_create
        procedure, pass                     :: open => xdmf_hdf5_open
        procedure, pass                     :: close => xdmf_hdf5_close
        procedure, pass                     :: scatter_to => xdmf_hdf5_scatter_to
    end type

    ! This structure stores all HDF5 ids of a file and one step
    type t_xdmf_file_descriptor
        type(t_xdmf_hdf5_metadata)              :: hdf5_meta_ids
        type(t_xdmf_hdf5)                       :: hdf5_ids_cells
#       if defined(_XDMF_PATCH)
            type(t_xdmf_hdf5)                   :: hdf5_ids_patches
#       endif

        contains

        procedure, pass                         :: scatter_to => xdmf_file_descriptor_scatter_to
    end type

    contains

    subroutine xdmf_hdf5_idpair_create(this, loc_id, type_id, name, dims)
        class(t_xdmf_hdf5_idpair), intent(inout)                        :: this
        integer(HID_T), intent(in)                                      :: loc_id, type_id
        character(len = *), intent(in)                                  :: name 
        integer(HSIZE_T), dimension(:), intent(in)                      :: dims

        integer                                                         :: hdf5_error

        call h5screate_simple_f(size(dims), dims, this%dspace_id, hdf5_error)
        call h5dcreate_f(loc_id, name, type_id, this%dspace_id, this%dset_id, hdf5_error)
        this%is_open = .true.
    end subroutine

    subroutine xdmf_hdf5_idpair_open(this, loc_id, name)
        class(t_xdmf_hdf5_idpair), intent(inout)                        :: this
        integer(HID_T), intent(in)                                      :: loc_id
        character(len = *), intent(in)                                  :: name 

        integer                                                         :: hdf5_error

        call h5dopen_f(loc_id, name, this%dset_id, hdf5_error)
        call h5dget_space_f(this%dset_id, this%dspace_id, hdf5_error)
        this%is_open = .true.
    end subroutine

    subroutine xdmf_hdf5_idpair_close(this)
        class(t_xdmf_hdf5_idpair), intent(inout)                        :: this

        integer                                                         :: hdf5_error

        if (this%is_open) then
            call h5dclose_f(this%dset_id, hdf5_error)
            call h5sclose_f(this%dspace_id, hdf5_error)
            this%is_open = .false.
        end if
    end subroutine

    subroutine xdmf_hdf5_idpair_scatter_to(this, child)
        class(t_xdmf_hdf5_idpair), intent(inout)                        :: this
        class(t_xdmf_hdf5_idpair), intent(inout)                        :: child

        child%dspace_id = this%dspace_id
        child%dset_id = this%dset_id
        child%is_open = this%is_open
    end subroutine


    ! This routine opens or creates a HDF5 file
    subroutine xdmf_hdf5_metadata_createopen(this, file_name, step, use_cache, allow_create, use_collective, truncate, step_exists, result)
        class(t_xdmf_hdf5_metadata), intent(inout)                      :: this
        character(len = *), intent(in)                                  :: file_name
        integer(XDMF_GRID_SI), intent(in)                               :: step
        logical, intent(in)                                             :: use_cache, allow_create, use_collective, truncate
        logical, intent(out)                                            :: step_exists
        integer, intent(out)                                            :: result

        integer                                                         :: hdf5_error, i, error, num_existing_groups, hdf5_existing_group_type, hdf5_existing_step
        character(len = 256)					                        :: hdf5_step_group_name, hdf5_existing_group_name
        logical                                                         :: file_exist
        integer (HID_T)                                                 :: hdf5_plist_access_id

        ! Create dataset access policy 
        call h5pcreate_f(H5P_DATASET_XFER_F, this%access_dset_id, hdf5_error)
#       if defined(_MPI)
            ! Use parallel access for datasets
            if (use_collective) then
                call h5pset_dxpl_mpio_f(this%access_dset_id, H5FD_MPIO_COLLECTIVE_F, hdf5_error)
            else
                call h5pset_dxpl_mpio_f(this%access_dset_id, H5FD_MPIO_INDEPENDENT_F, hdf5_error)
            end if
            ! Ensure all nodes are properly initialized
            call mpi_barrier(MPI_COMM_WORLD, error); assert_eq(error, 0)
#       endif

        result = 0
        if ((.not.use_cache).or.(this%file_id.eq.0)) then
            ! Create file access policy
            call h5pcreate_f(H5P_FILE_ACCESS_F, hdf5_plist_access_id, hdf5_error)
            inquire(file=file_name, exist=file_exist)
#           if defined(_MPI)
                call mpi_barrier(MPI_COMM_WORLD, error); assert_eq(error, 0)
                ! Use parallel access for this file
                call h5pset_fapl_mpio_f(hdf5_plist_access_id, MPI_COMM_WORLD, MPI_INFO_NULL, hdf5_error)
                ! Ensure all nodes are properly initialized
                call mpi_barrier(MPI_COMM_WORLD, error); assert_eq(error, 0)
#           endif

            if (file_exist) then
                ! Open existing file
                call h5fopen_f(file_name, H5F_ACC_RDWR_F, this%file_id, hdf5_error, access_prp = hdf5_plist_access_id)
                result = hdf5_error
            else
                if (allow_create) then
                    ! Create new file
                    call h5fcreate_f(file_name, H5F_ACC_EXCL_F, this%file_id, hdf5_error, access_prp = hdf5_plist_access_id)
                    result = hdf5_error
                else
                    result = -1
                end if
            end if

            call h5pclose_f(hdf5_plist_access_id, hdf5_error)
        end if

        if (result.eq.0) then
            ! Try to find this step 
            step_exists = .false.
            ! List existing steps
            call h5gn_members_f(this%file_id, "/"//char(0), num_existing_groups, hdf5_error)
            do i = 0, num_existing_groups - 1
                call h5gget_obj_info_idx_f(this%file_id, "/"//char(0), i, hdf5_existing_group_name, hdf5_existing_group_type, hdf5_error)
                read(hdf5_existing_group_name, *) hdf5_existing_step
                if(hdf5_existing_step.eq.step) then
                    step_exists = .true.
                    exit
                end if
            end do

            write (hdf5_step_group_name, "(I0, A)") step, char(0)
            if (step_exists) then
                if (truncate) then
                    ! Remove step
                    call h5ldelete_f(this%file_id, hdf5_step_group_name, hdf5_error)
                    ! Recreate step
                    call h5gcreate_f(this%file_id, hdf5_step_group_name, this%step_group_id, hdf5_error)
                    ! Recreate substeps
                    call h5gcreate_f(this%step_group_id, hdf5_gr_cells_dname, this%step_group_cells_id, hdf5_error)
#                   if defined(_XDMF_PATCH)
                        call h5gcreate_f(this%step_group_id, hdf5_gr_patches_dname, this%step_group_patches_id, hdf5_error)
#                   endif
                else
                    ! Open step group
                    call h5gopen_f(this%file_id, hdf5_step_group_name, this%step_group_id, hdf5_error)
                    ! Open substeps
                    call h5gopen_f(this%step_group_id, hdf5_gr_cells_dname, this%step_group_cells_id, hdf5_error)
#                   if defined(_XDMF_PATCH)
                        call h5gopen_f(this%step_group_id, hdf5_gr_patches_dname, this%step_group_patches_id, hdf5_error)
#                   endif
                end if
            else
                if (allow_create) then
                    ! Create step group
                    call h5gcreate_f(this%file_id, hdf5_step_group_name, this%step_group_id, hdf5_error)
                    ! Create substeps
                    call h5gcreate_f(this%step_group_id, hdf5_gr_cells_dname, this%step_group_cells_id, hdf5_error)
#                   if defined(_XDMF_PATCH)
                        call h5gcreate_f(this%step_group_id, hdf5_gr_patches_dname, this%step_group_patches_id, hdf5_error)
#                   endif
                end if
            end if
        end if
    end subroutine

    ! This routine closes a HDF5 file
    subroutine xdmf_hdf5_metadata_close(this)
        class(t_xdmf_hdf5_metadata), intent(inout)                      :: this

        integer                                                         :: hdf5_error

        call h5gclose_f(this%step_group_id, hdf5_error)
        call h5gclose_f(this%step_group_cells_id, hdf5_error)
#       if defined(_XDMF_PATCH)
            call h5gclose_f(this%step_group_patches_id, hdf5_error)
#       endif
        call h5pclose_f(this%access_dset_id, hdf5_error)
        call h5fflush_f(this%file_id, H5F_SCOPE_LOCAL_F, hdf5_error)
        ! call h5fclose_f(this%file_id, hdf5_error)
    end subroutine

    subroutine xdmf_hdf5_metadata_scatter_to(this, child)
        class(t_xdmf_hdf5_metadata), intent(inout)                      :: this
        class(t_xdmf_hdf5_metadata), intent(inout)                      :: child

        child%file_id = this%file_id
        child%step_group_id = this%step_group_id
        child%step_group_cells_id = this%step_group_cells_id
#       if defined(_XDMF_PATCH)
            child%step_group_patches_id = this%step_group_patches_id
#       endif
    end subroutine


    subroutine xdmf_hdf5_allocate(this, param)
        class(t_xdmf_hdf5), intent(inout)                               :: this
        type(t_xdmf_parameter), intent(in)                              :: param

        integer                                                         :: error

        allocate(this%batch_valsi(param%hdf5_valsi_width), stat = error); assert_eq(error, 0)
        allocate(this%batch_valsr(param%hdf5_valsr_width), stat = error); assert_eq(error, 0)
        allocate(this%batch_valsuv(param%hdf5_valsuv_width), stat = error); assert_eq(error, 0)
    end subroutine

    subroutine xdmf_hdf5_deallocate(this)
        class(t_xdmf_hdf5), intent(inout)                               :: this

        integer                                                         :: error

        deallocate(this%batch_valsi, stat = error); assert_eq(error, 0)
        deallocate(this%batch_valsr, stat = error); assert_eq(error, 0)
        deallocate(this%batch_valsuv, stat = error); assert_eq(error, 0)
    end subroutine

    ! This routine creates the datasets for a step
    subroutine xdmf_hdf5_create(this, param, loc_id, num_cells, tree_length)
        class(t_xdmf_hdf5), intent(inout)                               :: this
        type(t_xdmf_parameter), intent(in)                              :: param
        integer(HID_T), intent(in)                                      :: loc_id
        integer(XDMF_GRID_SI), intent(in)                               :: num_cells
        integer(XDMF_GRID_DI), intent(in)                               :: tree_length

        integer                                                         :: n, hdf5_error

        call this%allocate(param)
        call h5gcreate_f(loc_id, hdf5_attr_dname, this%attr_group_id, hdf5_error)
        if (tree_length.ne.0) then
            call this%tree%create(loc_id, H5T_NATIVE_INTEGER, &
                hdf5_tree_dname, (/ hdf5_tree_width, int(tree_length, HSIZE_T) /))
        end if
        do n = 1, param%hdf5_valsi_width
            call this%batch_valsi(n)%create(this%attr_group_id, H5T_NATIVE_INTEGER, &
                param%hdf5_valsi_dnames(n), (/ 1_HSIZE_T, int(num_cells, HSIZE_T) /))
        end do
        do n = 1, param%hdf5_valsr_width
            call this%batch_valsr(n)%create(this%attr_group_id, XDMF_HDF_P, &
            param%hdf5_valsr_dnames(n), (/ param%hdf5_attr_width, int(num_cells, HSIZE_T) /))
        end do
        do n = 1, param%hdf5_valsuv_width
            call this%batch_valsuv(n)%create(this%attr_group_id, XDMF_HDF_P, &
            param%hdf5_valsuv_dnames(n), (/ hdf5_vector_width, param%hdf5_attr_width, &
            int(num_cells, HSIZE_T) /))
        end do
        call this%valsg%create(loc_id, XDMF_HDF_P, &
            hdf5_valsg_dname, (/ param%hdf5_valsg_width, int(num_cells * param%hdf5_valst_width, HSIZE_T) /))
        call this%valst%create(loc_id, H5T_NATIVE_INTEGER, &
            hdf5_valst_dname, (/ param%hdf5_valst_width, int(num_cells, HSIZE_T) /))
    end subroutine

     ! This routine opens the datasets for a step
    subroutine xdmf_hdf5_open(this, param, loc_id)
        class(t_xdmf_hdf5), intent(inout)                               :: this
        type(t_xdmf_parameter), intent(in)                              :: param
        integer(HID_T), intent(in)                                      :: loc_id

        integer                                                         :: n, hdf5_error

        call this%allocate(param)
        call h5gopen_f(loc_id, hdf5_attr_dname, this%attr_group_id, hdf5_error)
        call this%tree%open(loc_id, hdf5_tree_dname)
        do n = 1, param%hdf5_valsi_width
            call this%batch_valsi(n)%open(this%attr_group_id, param%hdf5_valsi_dnames(n))
        end do
        do n = 1, param%hdf5_valsr_width
            call this%batch_valsr(n)%open(this%attr_group_id, param%hdf5_valsr_dnames(n))
        end do
        do n = 1, param%hdf5_valsuv_width
            call this%batch_valsuv(n)%open(this%attr_group_id, param%hdf5_valsuv_dnames(n))
        end do
        call this%valsg%open(loc_id, hdf5_valsg_dname)
        call this%valst%open(loc_id, hdf5_valst_dname)
    end subroutine

    ! This routine closes the datasets for a step
    subroutine xdmf_hdf5_close(this, param)
        class(t_xdmf_hdf5), intent(inout)                               :: this
        type(t_xdmf_parameter), intent(in)                              :: param

        integer                                                         :: n, hdf5_error

        call this%tree%close()
        do n = 1, param%hdf5_valsi_width
            call this%batch_valsi(n)%close()
        end do
        do n = 1, param%hdf5_valsr_width
            call this%batch_valsr(n)%close()
        end do
        do n = 1, param%hdf5_valsuv_width
            call this%batch_valsuv(n)%close()
        end do
        call this%valsg%close()
        call this%valst%close()
        call h5gclose_f(this%attr_group_id, hdf5_error)
        call this%deallocate()
    end subroutine

    ! This routine scatters / copies a hdf5 descriptor to a target
    subroutine xdmf_hdf5_scatter_to(this, param, child)
        class(t_xdmf_hdf5), intent(inout)                           :: this
        type(t_xdmf_parameter), intent(in)                          :: param
        class(t_xdmf_hdf5), intent(inout)                           :: child

        integer                                                     :: n

        child%attr_group_id = this%attr_group_id
        call this%tree%scatter_to(child%tree)
        call this%valsg%scatter_to(child%valsg)
        call this%valst%scatter_to(child%valst)
        do n = 1, param%hdf5_valsi_width
            call this%batch_valsi(n)%scatter_to(child%batch_valsi(n))
        end do
        do n = 1, param%hdf5_valsr_width
            call this%batch_valsr(n)%scatter_to(child%batch_valsr(n))
        end do
        do n = 1, param%hdf5_valsuv_width
            call this%batch_valsuv(n)%scatter_to(child%batch_valsuv(n))
        end do
    end subroutine

    ! This routine scatters / copies a hdf5 file descriptor to a target
    subroutine xdmf_file_descriptor_scatter_to(this, param, child)
        class(t_xdmf_file_descriptor), intent(inout)    :: this
        type(t_xdmf_parameter), intent(in)              :: param
        class(t_xdmf_file_descriptor), intent(inout)    :: child

        call this%hdf5_meta_ids%scatter_to(child%hdf5_meta_ids)
        call this%hdf5_ids_cells%scatter_to(param, child%hdf5_ids_cells)
#       if defined(_XDMF_PATCH)
            call this%hdf5_ids_patches%scatter_to(param, child%hdf5_ids_patches)
#       endif
    end subroutine


    ! Writes a 2D chunk of XDMF_ISO_P float data to a hdf5 file
    subroutine hdf5_write_chunk_real(dset_id, dspace_id, offset, dims, rank, stride, block, buffer, plist_access_id)
        integer(HID_T), intent(in)                                      :: dset_id, dspace_id, plist_access_id
        integer(HSIZE_T), dimension(:), intent(in)                      :: offset, dims, stride, block
        integer, intent(in)                                             :: rank
        real(XDMF_ISO_P), dimension(*), intent(inout)                   :: buffer

        integer(HID_T)                                                  :: memspace_id
        integer                                                         :: hdf5_error

        call h5screate_simple_f(rank, dims, memspace_id, hdf5_error)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, dims, hdf5_error, stride, block)
        call h5dwrite_f(dset_id, XDMF_HDF_P, buffer, dims, hdf5_error, memspace_id, dspace_id, xfer_prp = plist_access_id)
        call h5sclose_f(memspace_id, hdf5_error)
    end subroutine

    ! Writes a 2D chunk of int32 data to a hdf5 file
    subroutine hdf5_write_chunk_int(dset_id, dspace_id, offset, dims, rank, stride, block, buffer, plist_access_id)
        integer(HID_T), intent(in)                                      :: dset_id, dspace_id, plist_access_id
        integer(HSIZE_T), dimension(:), intent(in)                      :: offset, dims, stride, block
        integer, intent(in)                                             :: rank
        integer(INT32), dimension(*), intent(inout)                     :: buffer

        integer(HID_T)                                                  :: memspace_id
        integer                                                         :: hdf5_error


        call h5screate_simple_f(rank, dims, memspace_id, hdf5_error)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, dims, hdf5_error, stride, block)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, buffer, dims, hdf5_error, memspace_id, dspace_id, xfer_prp = plist_access_id)
        call h5sclose_f(memspace_id, hdf5_error)
    end subroutine

    ! Reads a 2D chunk of XDMF_ISO_P float data from a hdf5 file
    subroutine hdf5_read_chunk_real(dset_id, dspace_id, offset, dims, rank, stride, block, buffer, plist_access_id)
        integer(HID_T), intent(in)                                      :: dset_id, dspace_id, plist_access_id
        integer(HSIZE_T), dimension(:), intent(in)                      :: offset, dims, stride, block
        integer, intent(in)                                             :: rank
        real(XDMF_ISO_P), dimension(*), intent(inout)                   :: buffer

        integer(HID_T)                                                  :: memspace_id
        integer															:: error

        call h5screate_simple_f(rank, dims, memspace_id, error)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, dims, error, stride, block)
        call h5dread_f(dset_id, XDMF_HDF_P, buffer, dims, error, memspace_id, dspace_id, xfer_prp = plist_access_id)
        call h5sclose_f(memspace_id, error)
    end subroutine

    ! Reads a 2D chunk of int32 data from a hdf5 file
    subroutine hdf5_read_chunk_int(dset_id, dspace_id, offset, dims, rank, stride, block, buffer, plist_access_id)
        integer(HID_T), intent(in)                                      :: dset_id, dspace_id, plist_access_id
        integer(HSIZE_T), dimension(hdf5_rank), intent(in)              :: offset, dims, stride, block
        integer, intent(in)                                             :: rank
        integer(INT32), dimension(*), intent(inout)                     :: buffer

        integer(HID_T)                                                  :: memspace_id
        integer															:: error

        call h5screate_simple_f(hdf5_rank, dims, memspace_id, error)
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, offset, dims, error, stride, block)
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, buffer, dims, error, memspace_id, dspace_id, xfer_prp = plist_access_id)
        call h5sclose_f(memspace_id, error)
    end subroutine

end module