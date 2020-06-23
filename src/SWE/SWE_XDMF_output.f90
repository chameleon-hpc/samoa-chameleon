#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

#if defined(_SWE)
    module SWE_XDMF_output
        
        use SFC_edge_traversal
        use Samoa_swe
        use Tools_openmp

        use SWE_XDMF_config
        use XDMF_output_base
        use XDMF_xmf  

#       if defined(_SWE_PATCH)
            use SWE_PATCH
#       endif
#       if defined(_SWE_DG)
            use SWE_dg_matrices
            use SWE_data_types
            use SWE_dg_solver
#       endif 

        implicit none

        ! For traversal data, use the XDMF core data structure as a base
        type num_traversal_data
            type(t_xdmf_base_output_traversal) :: base
        end type

#       define _GT_EDGES
#		define _GT_NAME								    t_swe_xdmf_output_traversal

#		define _GT_PRE_TRAVERSAL_OP					    pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP				    post_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP				pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP				post_traversal_grid_op

#		define _GT_ELEMENT_OP						    element_op

#       define _GT_CELL_TO_EDGE_OP				        cell_to_edge_op_dg

#		include "SFC_generic_traversal_ringbuffer.f90"

        ! This routine creates an array of wrappedsection pointers, 
        ! because fortran does not like arrays of pointers.
        subroutine ptr_wrap_sections(traversal, sections_ptr)
            type(t_swe_xdmf_output_traversal), intent(inout)				:: traversal
            type(t_xdmf_base_output_traversal_ptr), &
                dimension(:), allocatable, intent(inout)                      :: sections_ptr

            integer                                                         :: error, i

            allocate(sections_ptr(size(traversal%sections)), stat = error); assert_eq(error, 0)
            do i = 1, size(traversal%sections)
                sections_ptr(i)%ptr => traversal%sections(i)%base
            end do
        end subroutine


        subroutine pre_traversal_grid_op(traversal, grid)
            type(t_swe_xdmf_output_traversal), intent(inout)				:: traversal
            type(t_grid), intent(inout)							            :: grid

            type(t_xdmf_base_output_traversal_ptr), &
                dimension(:), allocatable                                   :: sections_ptr
            integer                                                         :: error

            ! Pass this call the the core API
            call ptr_wrap_sections(traversal, sections_ptr)
#           if defined(_SWE_PATCH)
                call xdmf_base_pre_traversal_grid_op(traversal%base, sections_ptr, grid, swe_xdmf_param_cells, swe_xdmf_param_patches)
#           else
                call xdmf_base_pre_traversal_grid_op(traversal%base, sections_ptr, grid, swe_xdmf_param_cells)
#           endif
            deallocate(sections_ptr, stat = error); assert_eq(error, 0)
        end subroutine

        subroutine post_traversal_grid_op(traversal, grid)
            type(t_swe_xdmf_output_traversal), intent(inout)				:: traversal
            type(t_grid), intent(inout)							            :: grid

            type(t_xdmf_base_output_traversal_ptr), &
                dimension(:), allocatable                                   :: sections_ptr
            integer                                                         :: error

            ! Pass this call the the core API
            call ptr_wrap_sections(traversal, sections_ptr)
#           if defined(_SWE_PATCH)
                call xdmf_base_post_traversal_grid_op(traversal%base, sections_ptr, grid, swe_xdmf_param_cells, swe_xdmf_param_patches)
#           else
                call xdmf_base_post_traversal_grid_op(traversal%base, sections_ptr, grid, swe_xdmf_param_cells)
#           endif
            deallocate(sections_ptr, stat = error); assert_eq(error, 0)

            if(rank_MPI .eq. 0) then
                ! Output the XMF file
                call write_xdmf(traversal%base, grid)
            end if
    
            traversal%base%output_iteration = traversal%base%output_iteration + 1
        end subroutine

        subroutine pre_traversal_op(traversal, section)
            type(t_swe_xdmf_output_traversal), intent(inout)				:: traversal
            type(t_grid_section), intent(inout)							    :: section

            traversal%base%sect_store_index_cells = 1
#           if defined(_SWE_PATCH)
                traversal%base%sect_store_index_patches = 1
#           endif
        end subroutine

        subroutine post_traversal_op(traversal, section)
            type(t_swe_xdmf_output_traversal), intent(inout)				:: traversal
            type(t_grid_section), intent(inout)							    :: section

        end subroutine

        subroutine element_op(traversal, section, element)
            type(t_swe_xdmf_output_traversal), intent(inout)				:: traversal
            type(t_grid_section), intent(inout)							    :: section
            type(t_element_base), intent(inout)					            :: element

            real(XDMF_ISO_P), dimension(2, 3)			                    :: position
            integer(GRID_SI)	                                            :: i, offset_cells_buffer, offset_tree_buffer
#           if defined(_SWE_DG)
                integer(GRID_SI)	                                        :: point_id
                real(GRID_SR), dimension(_SWE_DG_DOFS)                      :: norm_b_dg, norm_h_dg, norm_hu_dg, norm_hv_dg
                real(GRID_SR), dimension(_SWE_PATCH_ORDER_SQUARE)           :: norm_b_fv, norm_h_fv, norm_hu_fv, norm_hv_fv
#           endif
#           if defined(_SWE_PATCH)
                integer(GRID_SI)	                                        :: offset_patches_buffer, patch_offs
#           endif
            integer(INT64)                                                  :: element_hash
            real(XDMF_ISO_P)                                                :: new_h, new_bigh
            integer(GRID_SI)                                                :: cell_offs
#           if defined (_SWE_PATCH)
                integer(GRID_SI)                                            :: j, patch_cell_id, patch_cell_offs, row, col
#           endif
            logical                                                         :: write_cp, filter_result_cells
#           if defined(_SWE_PATCH)
                logical                                                     :: filter_result_patches
#           endif

            filter_result_cells = .true.
            filter_result_patches = .true.

            ! Compute whether to output tree (checkpoint) data
            if(cfg%xdmf%i_xdmfcpint.eq.0) then
                write_cp = .false.
            else
                write_cp = mod(int(traversal%base%output_iteration, GRID_SI), int(cfg%xdmf%i_xdmfcpint, GRID_SI)).eq.0
            end if

            ! Evaluate filter if this step is not a checkpoint
            if ((.not. write_cp) .and. ((cfg%xdmf%i_xdmffilter_index .ne. 0) .or. &
                (iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_all) .eq. xdmf_output_mode_all))) then
                filter_result_cells = .false.
                if (iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_cells) .eq. xdmf_output_mode_cells) then
                    filter_result_cells = SWE_xdmf_filter(element, cfg%xdmf%i_xdmffilter_index, cfg%xdmf%i_xmdffilter_params_count, &
                        cfg%xdmf%r_xdmffilter_params_vector, .true.)
                end if
#               if defined(_SWE_PATCH)
                    filter_result_patches = .false.
                    if (iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_patches) .eq. xdmf_output_mode_patches) then
                        filter_result_patches = SWE_xdmf_filter(element, cfg%xdmf%i_xdmffilter_index, cfg%xdmf%i_xmdffilter_params_count, &
                            cfg%xdmf%r_xdmffilter_params_vector, .false.)
                    end if
#               endif
            end if

#           if defined (_SWE_PATCH)
                if (write_cp .or. filter_result_cells .or. filter_result_patches) then
#           else
                if (write_cp .or. filter_result_cells) then
#           endif
                ! Get thread buffer offset
                offset_cells_buffer = traversal%base%root_layout_desc%ranks(rank_MPI + 1)%sections(section%index)%offset_cells_buffer
#               if defined (_SWE_PATCH)
                    offset_tree_buffer = offset_cells_buffer / _SWE_PATCH_ORDER_SQUARE
                    offset_patches_buffer = traversal%base%root_layout_desc%ranks(rank_MPI + 1)%sections(section%index)%offset_patches_buffer
#               else
                    offset_tree_buffer = offset_cells_buffer
#               endif

                ! Check buffer overflows
                if((((iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_cells) .eq. xdmf_output_mode_cells) .and. filter_result_cells) .or. write_cp) .and. &
                    ((traversal%base%sect_store_index_cells + offset_tree_buffer .gt. size(traversal%base%sect_store_cells%ptr%tree)) .or. &
                    (traversal%base%sect_store_index_cells + offset_cells_buffer .gt. size(traversal%base%sect_store_cells%ptr%valsi, 1)))) then
                        _log_write(1, '(A, I0, A, I0, A, I0, A, I0, A)') " XDMF: Warning: Writing ", &
                            traversal%base%sect_store_index_cells + offset_cells_buffer, &
                            " cells (tree: ", traversal%base%sect_store_index_cells + offset_tree_buffer, "), into buffer of size ", &
                            size(traversal%base%sect_store_cells%ptr%valsi, 1), ", (tree: ", size(traversal%base%sect_store_cells%ptr%tree), "). Skipping."
                else if((((iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_patches) .eq. xdmf_output_mode_patches) .and. filter_result_patches) .or. write_cp) .and. &
                    ((traversal%base%sect_store_index_patches + offset_patches_buffer .gt. size(traversal%base%sect_store_patches%ptr%valsi, 1)))) then
                        _log_write(1, '(A, I0, A, I0, A)') " XDMF: Warning: Writing ", &
                            traversal%base%sect_store_index_patches + offset_patches_buffer, &
                            " patches, into buffer of size ", &
                            size(traversal%base%sect_store_patches%ptr%valsi, 1), ". Skipping."
                else
                    if (write_cp) then
                        ! Compute the cells offset in the topology tree
                        call xdmf_hash_element(element, int(traversal%base%grid_scale, INT64), element_hash)
                        ! Write the existence of this cell to the tree data buffer for this section
                        if (element_hash .gt. (2_INT64**31 - 1)) then
                            _log_write(1, *) " XDMF: Error: Cell hash does not fit into 32 bit signed integer.", &
                                "This is not yet implemented. This cell will be corrupt. Reduce level of detail / depth."
                            assert(.false.)
                            element_hash = 0
                        end if
                        traversal%base%sect_store_cells%ptr%tree(traversal%base%sect_store_index_cells + offset_tree_buffer) = int(element_hash, INT32)
                    end if
#                   if defined (_SWE_PATCH)
                        if ((filter_result_cells .and. (iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_cells) .eq. xdmf_output_mode_cells)) .or. write_cp) then
                            norm_h_fv = element%cell%data_pers%H
                            norm_hu_fv = element%cell%data_pers%HU
                            norm_hv_fv = element%cell%data_pers%HV
                            norm_b_fv = element%cell%data_pers%B
                            ! Apply attribute size correction
                            ! This cell shall be output as FV
#                           if defined(_SWE_DG)    
                                if (isDG(element%cell%data_pers%troubled)) then
                                    ! DG stores H, FV stores H + B, so we normalize DG to FV
                                    call apply_phi(element%cell%data_pers%Q%H + element%cell%data_pers%Q%b, norm_h_fv)
                                    call apply_phi(element%cell%data_pers%Q%p(1), norm_hu_fv)
                                    call apply_phi(element%cell%data_pers%Q%p(2), norm_hv_fv)
                                    call apply_phi(element%cell%data_pers%Q%b, norm_b_fv)
                                end if
#                           endif
                            ! Depending on the iteration direction (sign. plotter type), the order of the cells inside
                            ! a patch may need to be flipped
                            row = 1
                            col = 1
                            do i = 1, _SWE_PATCH_ORDER_SQUARE
                                if (element%cell%geometry%i_plotter_type.gt.0) then 
                                    patch_cell_id = i
                                else
                                    patch_cell_id = (row - 1) * (row - 1) + 2 * row - col
                                end if

                                ! Compute the global offset in the datasets of this cell
                                ! Note that the triangle order inside a patch is normalized here, i.e. independent from plotter direction
                                patch_cell_offs = ((traversal%base%sect_store_index_cells - 1) * _SWE_PATCH_ORDER_SQUARE) + patch_cell_id
                                cell_offs = patch_cell_offs + offset_cells_buffer
                            
                                ! Compute and store the actual position of this subcell in the domain
                                forall(j = 1:swe_xdmf_param_cells%hdf5_valst_width) traversal%base%sect_store_cells%ptr%valsg(:, &
                                    ((cell_offs - 1) * swe_xdmf_param_cells%hdf5_valst_width) + j) = &
                                    real(cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, SWE_PATCH_geometry%coords(:, j, i)) + cfg%offset, REAL32)

                                ! Store cell values, see SWE implementation for details
                                traversal%base%sect_store_cells%ptr%valsi(cell_offs, :) = (/ &
                                    int(element%cell%geometry%i_depth, INT32), &
                                    int(rank_MPI, INT32), &
                                    int(element%cell%geometry%i_plotter_type, INT32), &
                                    int(section%index, INT32) &
#                                   if defined(_SWE_DG)
                                        , int(element%cell%data_pers%troubled, INT32) &
#                                   endif
                                    /)

                                ! Store point data in traversal buffer
                                ! h = WaterLevel = H + B
                                ! bigh = WaterHeight = WaterLevel - B = (H + B) - B
                                ! clamp WaterHeight to 0, negative value = below landmass
                                new_bigh = max(0.0_XDMF_ISO_P, real(norm_h_fv(patch_cell_id) - norm_b_fv(patch_cell_id), XDMF_ISO_P))
                                new_h = new_bigh + real(norm_b_fv(patch_cell_id), XDMF_ISO_P) ! TODO
                                do j=1, swe_xdmf_param_cells%hdf5_attr_width
                                    traversal%base%sect_store_cells%ptr%valsr(j, cell_offs, :) = (/ &
                                        real(norm_b_fv(patch_cell_id), XDMF_ISO_P), new_h, new_bigh /)
                                    traversal%base%sect_store_cells%ptr%valsuv(:, j, cell_offs, swe_hdf5_valsuv_f_offset) = &
                                        real((/ norm_hu_fv(patch_cell_id), norm_hv_fv(patch_cell_id) /), XDMF_ISO_P)
                                end do

                                ! Compute the subcells cartesian position
                                do j = 1, swe_xdmf_param_cells%hdf5_valst_width
                                    position(:, j) = real(cfg%scaling * samoa_barycentric_to_world_point(element%transform_data, SWE_PATCH_geometry%coords(:, j, i)) + cfg%offset, REAL32)
                                end do
                                forall(j = 1:swe_xdmf_param_cells%hdf5_valst_width) &
                                    traversal%base%sect_store_cells%ptr%valsg(:, &
                                    ((cell_offs - 1) * swe_xdmf_param_cells%hdf5_valst_width) + j) = position(:, j)

                                col = col + 1
                                if (col.eq.(2 * row)) then
                                    col = 1
                                    row = row + 1
                                end if
                            end do

                            traversal%base%sect_store_index_cells = traversal%base%sect_store_index_cells + 1
                        end if

                        if ((filter_result_patches .and. (iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_patches) .eq. xdmf_output_mode_patches)) .or. write_cp) then
                            norm_h_dg = element%cell%data_pers%Q%H
                            norm_hu_dg = element%cell%data_pers%Q%p(1)
                            norm_hv_dg = element%cell%data_pers%Q%p(2)
                            norm_b_dg = element%cell%data_pers%Q%b
#                           if defined(_SWE_DG)
                                ! Apply attribute size correction
                                ! This cell shall be output as DG
                                if (isFV(element%cell%data_pers%troubled)) then
                                    ! DG stores H, FV stores H + B, so we normalize FV to DG
                                    call apply_mue(element%cell%data_pers%H - element%cell%data_pers%B, norm_h_dg)
                                    call apply_mue(element%cell%data_pers%HU, norm_hu_dg)
                                    call apply_mue(element%cell%data_pers%HV, norm_hv_dg)
                                    call apply_mue(element%cell%data_pers%B, norm_b_dg)
                                end if
#                           endif
                            ! Store the patch, too
                            patch_offs = traversal%base%sect_store_index_patches + offset_patches_buffer
                            ! Store patch values, see SWE implementation for details
                            traversal%base%sect_store_patches%ptr%valsi(patch_offs, :) = (/ &
                                int(element%cell%geometry%i_depth, INT32), &
                                int(rank_MPI, INT32), &
                                int(element%cell%geometry%i_plotter_type, INT32), &
                                int(section%index, INT32) &
#                               if defined(_SWE_DG)
                                    , int(element%cell%data_pers%troubled, INT32) &
#                               endif
                                /)
                                
#                           if defined(_SWE_DG)
                                ! Depending on the iteration direction (sign. plotter type), the order of the DOFs inside
                                ! a cell may need to be flipped
                                do i = 1, _SWE_DG_DOFS
                                    if (element%cell%geometry%i_plotter_type .le. 0) then 
                                        point_id = i
                                    else
                                        point_id = mirrored_coords(i)
                                    end if
                                    
                                    ! Store point data in traversal buffer
                                    traversal%base%sect_store_patches%ptr%valsr(i, patch_offs, swe_hdf5_valsr_b_offset) = &
                                        real(norm_b_dg(point_id), XDMF_ISO_P)
                                    ! DG stores H, FV stores H + B, so we normalize DG to FV
                                    ! clamp WaterHeight to 0, negative value = below landmass
                                    new_bigh = max(0.0_XDMF_ISO_P, real(norm_h_dg(point_id), XDMF_ISO_P))
                                    new_h = new_bigh + real(norm_b_dg(point_id), XDMF_ISO_P)  ! TODO
                                    traversal%base%sect_store_patches%ptr%valsr(i, patch_offs, swe_hdf5_valsr_h_offset) = new_h ! TODO
                                    traversal%base%sect_store_patches%ptr%valsr(i, patch_offs, swe_hdf5_valsr_bh_offset) = new_bigh
                                    traversal%base%sect_store_patches%ptr%valsuv(:, i, patch_offs, swe_hdf5_valsuv_f_offset) = &
                                        real((/ norm_hu_dg(point_id), norm_hv_dg(point_id) /), XDMF_ISO_P)
                                end do
#                           endif

                            ! Compute the patchs cartesian position
                            ! Rotate the vertices by one, because thats just how it is
                            do i = 1, swe_xdmf_param_patches%hdf5_valst_width
                                position(:, i) = real((cfg%scaling * element%nodes(mod(i, swe_xdmf_param_patches%hdf5_valst_width) + 1)%ptr%position) &
                                    + cfg%offset(:), XDMF_ISO_P) 
                            end do
                            forall(i = 1:swe_xdmf_param_patches%hdf5_valst_width) &
                                traversal%base%sect_store_patches%ptr%valsg(:, &
                                ((patch_offs - 1) * swe_xdmf_param_patches%hdf5_valst_width) + i) = position(:, i)

                            traversal%base%sect_store_index_patches = traversal%base%sect_store_index_patches + 1
                        end if
#                   else
                        assert(.false.)
#                       error SWE_PATCH is missing
#                   endif
                end if
            end if
        end subroutine

        ! This routine generates the XMF file needed to index the HDF5 files
        subroutine write_xdmf(base, grid)
            type(t_xdmf_base_output_traversal), intent(inout)				:: base
            type(t_grid), intent(inout)							            :: grid
    
            character(len = 256)					                        :: file_name_h5, file_name_xmf
            integer                                                         :: xml_file_id
            integer(GRID_SI)                                                :: output_meta_iteration
            character(len = 17)                                             :: xml_time_string
            logical                                                         :: file_xmf_exist
            integer                                                         :: write_cp, i
    
            xml_file_id = 42

            write(file_name_xmf, "(A, A, A)") trim(base%s_file_stamp), "_xdmf.xmf"
    
            inquire(file=trim(file_name_xmf)//char(0), exist=file_xmf_exist)
            if(file_xmf_exist) then
                ! If file exists, remove last line and write from there
                open(unit=xml_file_id, file=trim(file_name_xmf)//char(0), status="old", position="append", access="sequential", form="formatted")
                ! If needed, remove newer timesteps
                do i = 1, base%xdmf_remove_lines
                    backspace(xml_file_id)
                end do
            else
                ! If file does not exist, write XML header
                open(unit=xml_file_id, file=trim(file_name_xmf)//char(0), status="new", action="write", access="sequential", form="formatted")
                write(xml_file_id, "(A, A, A, A, A, A)", advance="yes") '<?xml version="1.0" encoding="UTF-8"?><Xdmf Version="3.0"><Domain>', &
                    '<Information Name="CommandLine" Value="', trim(cfg%xdmf%s_krakencommand), &
                    '" /><Information Name="FileStamp" Value="', trim(base%s_file_stamp), &
                    '" /><Grid GridType="Collection" CollectionType="Temporal">'
            end if
    
            output_meta_iteration = int(real(base%output_iteration, GRID_SR) / cfg%xdmf%i_xdmfspf, GRID_SI)
            write (file_name_h5, "(A, A, I0, A, A)") trim(base%s_file_stamp), "_", output_meta_iteration, "_xdmf.h5", char(0)
    
            write(xml_file_id, "(A)", advance="no") '<Grid GridType="Collection" CollectionType="Spatial">'

            ! Compute whether to output tree (checkpoint) data
            write_cp = 0
            if((cfg%xdmf%i_xdmfcpint.ne.0) .and. &
                mod(int(base%output_iteration, GRID_SI), int(cfg%xdmf%i_xdmfcpint, GRID_SI)).eq.0) then
                write_cp = 1
            end if
            write(xml_file_id, "(A, I0, A)", advance="no") '<Information Name="CP" Value="', write_cp, '" />'
    
            ! Time
            xml_time_string(:) = " "
            if (base%i_sim_iteration.eq.0) then
                write(xml_time_string, "(F0.8)") (-1.0_GRID_SR / real(base%output_iteration + 1, GRID_SR))
            else
                write(xml_time_string, "(F0.8)") grid%r_time
            end if
            write(xml_file_id, "(A, A, A, A, I0, A, A, F0.8, A)", advance="no") '<Time Value="', trim(xml_time_string), '" />', &
                '<Attribute Name="Step" Center="Grid"><DataItem Format="XML" NumberType="Int" Dimensions="1">', &
                base%i_sim_iteration, '</DataItem></Attribute>', &
                '<Attribute Name="DeltaTime" Center="Grid"><DataItem Format="XML" NumberType="Float" Dimensions="1">', &
                grid%r_dt, '</DataItem></Attribute>'

            if ((iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_cells) .eq. xdmf_output_mode_cells) .or. write_cp .eq. 1) then
                call write_xdmf_data(base, swe_xdmf_param_cells, xml_file_id, int(base%num_cells, HSIZE_T), output_meta_iteration, hdf5_gr_cells_dname_nz, .false.)
            end if
#           if defined(_SWE_PATCH)
                if ((iand(cfg%xdmf%i_xdmfoutput_mode, xdmf_output_mode_patches) .eq. xdmf_output_mode_patches) .or. write_cp .eq. 1) then
                    call write_xdmf_data(base, swe_xdmf_param_patches, xml_file_id, int(base%num_patches, HSIZE_T), output_meta_iteration, &
                        hdf5_gr_patches_dname_nz, cfg%xdmf%l_xdmfoutput_lagrange, swe_xdmf_param_patches%hdf5_valst_width)
                end if
#           endif

            write(xml_file_id, "(A)", advance="yes") '</Grid>'
    
            write(xml_file_id, "(A)", advance="no") '</Grid></Domain></Xdmf>'
            close(unit=xml_file_id)
        end subroutine

        ! This routine generates the XMF attributes for a grid
        subroutine write_xdmf_data(base, param, xml_file_id, num_cells, output_meta_iteration, subgroup_dname_nz, lagrange, hdf5_attr_width_override)
            type(t_xdmf_base_output_traversal), intent(inout)				    :: base
            type(t_xdmf_parameter), intent(in)                                  :: param
            integer, intent(in)                                                 :: xml_file_id
            integer(HSIZE_T), intent(in)                                        :: num_cells
            integer(GRID_SI), intent(in)                                        :: output_meta_iteration
            character(*), intent(in)                                            :: subgroup_dname_nz
            logical, intent(in)                                                 :: lagrange
            integer (XDMF_GRID_DI), intent(in), optional	                    :: hdf5_attr_width_override

            character(len = 21)                                                 :: xml_dims_string
            character(len = 512)					                            :: xml_hdf5_path_string
            integer (XDMF_GRID_DI)	                                            :: hdf5_attr_width_override_p


            if ((.not. lagrange) .and. present(hdf5_attr_width_override)) then
                hdf5_attr_width_override_p = hdf5_attr_width_override
            else
                hdf5_attr_width_override_p =  param%hdf5_attr_width
            end if

            write(xml_file_id, "(A, A, A)", advance="no") '<Grid><Information Name="Layer" Value="', subgroup_dname_nz, '" />'

            ! Topology
            xml_dims_string(:) = " "
            write (xml_dims_string, "(I0, A, I0)") num_cells, " ", param%hdf5_valst_width
            xml_hdf5_path_string(:) = " "
            write (xml_hdf5_path_string, "(A, A, I0, A, I0, A, A, A, A)") trim(base%s_file_stamp_base), "_", &
                output_meta_iteration, "_xdmf.h5:/", base%output_iteration, "/", subgroup_dname_nz, "/", hdf5_valst_dname_nz
            write(xml_file_id, "(A, A, A, A, A)", advance="no") '<Topology TopologyType="Triangle"><DataItem Format="HDF" NumberType="Int" Dimensions="', &
                trim(xml_dims_string),'">', trim(xml_hdf5_path_string),'</DataItem></Topology>'
    
            ! Geometry
            xml_dims_string(:) = " "
            write (xml_dims_string, "(I0, A, I0)") (num_cells * param%hdf5_valst_width), " ", param%hdf5_valsg_width
            xml_hdf5_path_string(:) = " "
            write (xml_hdf5_path_string, "(A, A, I0, A, I0, A, A, A, A)") trim(base%s_file_stamp_base), "_", &
                output_meta_iteration, "_xdmf.h5:/", base%output_iteration, "/", subgroup_dname_nz, "/", hdf5_valsg_dname_nz
            write(xml_file_id, "(A, I0, A, A, A, A, A)", advance="no") '<Geometry GeometryType="XY"><DataItem Format="HDF" NumberType="Float" Precision="', &
                XDMF_XMF_P, '" Dimensions="', trim(xml_dims_string),'">', trim(xml_hdf5_path_string),'</DataItem></Geometry>'
    
            ! Cell attributes
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                subgroup_dname_nz, num_cells, 0_HSIZE_T, 0_HSIZE_T, swe_hdf5_attr_depth_dname_nz, "Depth", &
                .true., .true., .false., xml_file_id)
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                subgroup_dname_nz, num_cells, 0_HSIZE_T, 0_HSIZE_T, swe_hdf5_attr_rank_dname_nz, "Rank", &
                .true., .true., .false., xml_file_id)
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                subgroup_dname_nz, num_cells, 0_HSIZE_T, 0_HSIZE_T, swe_hdf5_attr_plotter_dname_nz, "Plotter", &
                .true., .true., .false., xml_file_id)
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                subgroup_dname_nz, num_cells, 0_HSIZE_T, 0_HSIZE_T, swe_hdf5_attr_section_dname_nz, "Section", &
                .true., .true., .false., xml_file_id)
#           if defined(_SWE_DG)                                             
                call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                    subgroup_dname_nz, num_cells, 0_HSIZE_T, 0_HSIZE_T, swe_hdf5_attr_troubled_dname_nz, "Troubled", &
                    .true., .true., .false., xml_file_id)
#           endif
            ! Simulation values
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                subgroup_dname_nz, num_cells, hdf5_attr_width_override_p, 0_HSIZE_T, swe_hdf5_attr_b_dname_nz, "Bathymetry", &
                .false., .false., lagrange, xml_file_id)
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                subgroup_dname_nz, num_cells, hdf5_attr_width_override_p, 0_HSIZE_T, swe_hdf5_attr_bh_dname_nz, "WaterHeight", &
                .false., .false., lagrange, xml_file_id)
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                subgroup_dname_nz, num_cells, hdf5_attr_width_override_p, 0_HSIZE_T, swe_hdf5_attr_h_dname_nz, "WaterLevel", &
                .false., .false., lagrange, xml_file_id)
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                subgroup_dname_nz, num_cells, hdf5_attr_width_override_p, 2_HSIZE_T, swe_hdf5_attr_f_dname_nz, "Momentum", &
                .false., .false., lagrange, xml_file_id)
    
            write(xml_file_id, "(A)", advance="no") '</Grid>'
        end subroutine

    end module
#endif
