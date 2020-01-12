#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

#if defined(_SWE)
    module SWE_XDMF_output
        
        use SFC_edge_traversal
        use Samoa_swe
        use SWE_euler_timestep
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

#       define _GT_CELL_TO_EDGE_OP				        cell_to_edge_op

#		include "SFC_generic_traversal_ringbuffer.f90"

        ! This routine creates an array of wrappedsection pointers, 
        ! because fortran does not like arrays of pointers.
        subroutine ptr_wrap_sections(traversal, sections_ptr)
            type(t_swe_xdmf_output_traversal), intent(inout)				:: traversal
            type(t_xdmf_base_output_traversal_ptr), &
                dimension(:), allocatable, intent(out)                      :: sections_ptr

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
            call xdmf_base_pre_traversal_grid_op(traversal%base, sections_ptr, grid, swe_xdmf_param)
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
            call xdmf_base_post_traversal_grid_op(traversal%base, sections_ptr, grid, swe_xdmf_param)
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

            traversal%base%sect_store_index = 1
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
            integer(GRID_SI)	                                            :: i, offset_cells_buffer, offset_tree_buffer, point_id
            integer(INT64)                                                  :: element_hash
            real(XDMF_ISO_P)                                                :: new_h, new_bigh
            integer(GRID_SI)                                                :: cell_offs
            logical                                                         :: write_cp, filter_result = .true.

            ! Compute whether to output tree (checkpoint) data
            if(cfg%xdmf%i_xdmfcpint.eq.0) then
                write_cp = .false.
            else
                write_cp = mod(int(traversal%base%output_iteration, GRID_SI), int(cfg%xdmf%i_xdmfcpint, GRID_SI)).eq.0
            end if

            ! Evaluate filter if this step is not a checkpoint
            if ((.not. write_cp) .and. (cfg%xdmf%i_xdmffilter_index .ne. 0)) then
                call swe_xdmf_filter(element, cfg%xdmf%i_xdmffilter_index, cfg%xdmf%i_xmdffilter_params_count, &
                    cfg%xdmf%r_xdmffilter_params_vector, filter_result)
            end if

            if (write_cp .or. filter_result) then
                ! Compute the cells cartesian position
                do i = 1, 3
                    position(:, i) = real((cfg%scaling * element%nodes(i)%ptr%position) + cfg%offset(:), XDMF_ISO_P) 
                end do

                ! Get thread buffer offset
                offset_cells_buffer = traversal%base%root_layout_desc%ranks(rank_MPI + 1)%sections(section%index)%offset_cells_buffer
                offset_tree_buffer = offset_cells_buffer

                ! Check buffer overflows
                if((traversal%base%sect_store_index + offset_tree_buffer .gt. size(traversal%base%sect_store%ptr%tree)) .or. &
                    (traversal%base%sect_store_index + offset_cells_buffer .gt. size(traversal%base%sect_store%ptr%valsi, 1))) then
                        _log_write(1, '(A, I0, A, I0, A, I0, A, I0, A)') " XDMF: Warning: Writing ", &
                            traversal%base%sect_store_index + offset_cells_buffer, &
                            " cells (tree: ", traversal%base%sect_store_index + offset_tree_buffer, "), into buffer of size ", &
                            size(traversal%base%sect_store%ptr%valsi, 1), ", (tree: ", size(traversal%base%sect_store%ptr%tree), "). Skipping."
                else
                    ! Compute the cells offset in the topology tree
                    call xdmf_hash_element(element, int(traversal%base%grid_scale, INT64), element_hash)

                    ! Write the existence of this cell to the tree data buffer for this section
                    if (element_hash .gt. (2_INT64**31 - 1)) then
                        _log_write(1, *) " XDMF: Error: Cell hash does not fit into 32 bit signed integer.", &
                            "This is not yet implemented. This cell will be corrupt. Reduce level of detail / depth."
                        assert(.false.)
                        element_hash = 0
                    end if
                    traversal%base%sect_store%ptr%tree(traversal%base%sect_store_index + offset_tree_buffer) = int(element_hash, INT32)

                    cell_offs = traversal%base%sect_store_index + offset_cells_buffer
                    ! Store cell values, see SWE implementation for details
                    traversal%base%sect_store%ptr%valsi(cell_offs, :) = (/ &
                        int(element%cell%geometry%i_depth, INT32), &
                        int(rank_MPI, INT32), &
                        int(element%cell%geometry%i_plotter_type, INT32), &
                        int(section%index, INT32) /)

                    do i = 1, swe_xdmf_param%hdf5_attr_width
                        ! Compute data storage order in cell
                        if (element%cell%geometry%i_plotter_type .le. 0) then
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
                        
                        ! Store point data in traversal buffer
                        ! TODO XDMF write data into buffer
                        ! new_h = real(element%cell%data_pers%Q(point_id)%h, XDMF_ISO_P)
                        ! if (new_h.le.cfg%dry_tolerance) then
                        !     new_h = 0
                        ! else
                        !     new_h = new_h + real(element%cell%data_pers%Q(point_id)%b, XDMF_ISO_P)
                        ! end if
                        ! new_bigh = real(element%cell%data_pers%Q(point_id)%h, XDMF_ISO_P)
                        ! traversal%base%sect_store%ptr%valsr(i, cell_offs, :) = (/ &
                        !     real(element%cell%data_pers%Q(point_id)%b, XDMF_ISO_P), new_h, new_bigh /)
                        ! traversal%base%sect_store%ptr%valsuv(:, i, cell_offs, swe_hdf5_valsuv_f_offset) = &
                        !     real(element%cell%data_pers%Q(point_id)%p(:), XDMF_ISO_P)
                        ! traversal%base%sect_store%ptr%valsuv(:, i, cell_offs, swe_hdf5_valsuv_g_offset) = &
                        !     real(element%cell%data_pers%Q(point_id)%gradB(:), XDMF_ISO_P)
                    end do
                    forall(i = 1:swe_xdmf_param%hdf5_valst_width) &
                        traversal%base%sect_store%ptr%valsg(:, &
                        ((cell_offs - 1) * swe_xdmf_param%hdf5_valst_width) + i) = position(:, i)
                end if
                traversal%base%sect_store_index = traversal%base%sect_store_index + 1
            end if
        end subroutine

        ! This routine generates the XMF file needed to index the HDF5 files
        subroutine write_xdmf(base, grid)
            type(t_xdmf_base_output_traversal), intent(inout)				:: base
            type(t_grid), intent(inout)							            :: grid
    
            character(len = 256)					                        :: file_name_h5, file_name_xmf
            integer                                                         :: xml_file_id = 42
            integer(GRID_SI)                                                :: output_meta_iteration
            integer(HSIZE_T)                                                :: num_cells
            character(len = 17)                                             :: xml_time_string
            character(len = 21)                                             :: xml_dims_string
            character(len = 512)					                        :: xml_hdf5_path_string
            logical                                                         :: file_xmf_exist
            integer                                                         :: is_cp, i
    
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
    
            write(xml_file_id, "(A)", advance="no") '<Grid>'

            ! Compute whether to output tree (checkpoint) data
            is_cp = 0
            if((cfg%xdmf%i_xdmfcpint.ne.0) .and. &
                mod(int(base%output_iteration, GRID_SI), int(cfg%xdmf%i_xdmfcpint, GRID_SI)).eq.0) then
                is_cp = 1
            end if
            write(xml_file_id, "(A, I0, A)", advance="no") '<Information Name="CP" Value="', is_cp, '" />'
    
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
    
            num_cells = int(base%num_cells, HSIZE_T)
    
            ! Topology
            xml_dims_string(:) = " "
            write (xml_dims_string, "(I0, A, I0)") num_cells, " ", swe_xdmf_param%hdf5_valst_width
            xml_hdf5_path_string(:) = " "
            write (xml_hdf5_path_string, "(A, A, I0, A, I0, A, A)") trim(base%s_file_stamp_base), "_", &
                output_meta_iteration, "_xdmf.h5:/", base%output_iteration, "/", hdf5_valst_dname_nz
            write(xml_file_id, "(A, A, A, A, A)", advance="no") '<Topology TopologyType="Triangle"><DataItem Format="HDF" NumberType="Int" Dimensions="', &
                trim(xml_dims_string),'">', trim(xml_hdf5_path_string),'</DataItem></Topology>'
    
            ! Geometry
            xml_dims_string(:) = " "
            write (xml_dims_string, "(I0, A, I0)") (num_cells * swe_xdmf_param%hdf5_valst_width), " ", swe_xdmf_param%hdf5_valsg_width
            xml_hdf5_path_string(:) = " "
            write (xml_hdf5_path_string, "(A, A, I0, A, I0, A, A)") trim(base%s_file_stamp_base), "_", &
                output_meta_iteration, "_xdmf.h5:/", base%output_iteration, "/", hdf5_valsg_dname_nz
            write(xml_file_id, "(A, I0, A, A, A, A, A)", advance="no") '<Geometry GeometryType="XY"><DataItem Format="HDF" NumberType="Float" Precision="', &
                XDMF_XMF_P, '" Dimensions="', trim(xml_dims_string),'">', trim(xml_hdf5_path_string),'</DataItem></Geometry>'
    
            ! Cell attributes
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                num_cells, 0_HSIZE_T, 0_HSIZE_T, swe_hdf5_attr_depth_dname_nz, "Depth", .true., .true., xml_file_id)
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                num_cells, 0_HSIZE_T, 0_HSIZE_T, swe_hdf5_attr_rank_dname_nz, "Rank", .true., .true., xml_file_id)
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                num_cells, 0_HSIZE_T, 0_HSIZE_T, swe_hdf5_attr_plotter_dname_nz, "Plotter", .true., .true., xml_file_id)
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                num_cells, 0_HSIZE_T, 0_HSIZE_T, swe_hdf5_attr_section_dname_nz, "Section", .true., .true., xml_file_id)
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                num_cells, swe_xdmf_param%hdf5_attr_width, 0_HSIZE_T, swe_hdf5_attr_b_dname_nz, "Bathymetry", .false., .false., xml_file_id)
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                num_cells, swe_xdmf_param%hdf5_attr_width, 0_HSIZE_T, swe_hdf5_attr_bh_dname_nz, "WaterHeight", .false., .false., xml_file_id)
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                num_cells, swe_xdmf_param%hdf5_attr_width, 0_HSIZE_T, swe_hdf5_attr_h_dname_nz, "WaterLevel", .false., .false., xml_file_id)
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                num_cells, swe_xdmf_param%hdf5_attr_width, 2_HSIZE_T, swe_hdf5_attr_f_dname_nz, "Momentum", .false., .false., xml_file_id)
            call xdmf_xmf_add_attribute(base%output_iteration, output_meta_iteration, base%s_file_stamp_base, &
                num_cells, swe_xdmf_param%hdf5_attr_width, 2_HSIZE_T, swe_hdf5_attr_g_dname_nz, "BathymetryGradient", .false., .false., xml_file_id)
    
            write(xml_file_id, "(A)", advance="yes") '</Grid>'
    
            write(xml_file_id, "(A)", advance="no") '</Grid></Domain></Xdmf>'
            close(unit=xml_file_id)
        end subroutine

    end module
#endif
