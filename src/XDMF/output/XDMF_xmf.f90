#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

! This module provides routines to ease the generation of XDMF files
module XDMF_xmf

    use XDMF_data_types

    use, intrinsic :: iso_fortran_env

    implicit none

    contains

    ! This routine generates an XML tag describing a cell attribute in XDMF
    subroutine xdmf_xmf_add_attribute(output_iteration, output_meta_iteration, &
            file_stamp_base, subgroup_dname_nz, dim_length, dim_width, dim_depth, attr_dname_nz, attr_cname, is_int, is_cell, xml_file_id)
        integer (XDMF_GRID_SI), intent(in)					            :: output_iteration, output_meta_iteration
        character (len = *), intent(in)					                :: file_stamp_base, subgroup_dname_nz, attr_dname_nz, attr_cname
        integer (XDMF_GRID_DI), intent(in)					            :: dim_length, dim_width, dim_depth
        logical, intent(in)                                             :: is_int, is_cell
        integer, intent(in)                                             :: xml_file_id

        character(len = 512)					                        :: xml_hdf5_path_string
        character(len = 32)                                             :: xml_dims_string
        character(len = 32)                                             :: xml_number_type_string
        character(len = 4)                                              :: xml_center_type_string

        xml_number_type_string = " "
        if(is_int) then
            write(xml_number_type_string, "(A)") "Int"
        else
            write(xml_number_type_string, "(A, I0)") 'Float" Precision="', XDMF_XMF_P
        end if

        xml_center_type_string = " "
        if(is_cell) then
            write(xml_center_type_string, "(A)") "Cell"
        else
            write(xml_center_type_string, "(A)") "Node"
        end if

        xml_dims_string(:) = " "
        if(dim_width.ne.0_HSIZE_T) then
            if(dim_depth.ne.0_HSIZE_T) then
                write (xml_dims_string, "(I0, A, I0, A, I0)") dim_length, " ", dim_width, " ", dim_depth
            else
                write (xml_dims_string, "(I0, A, I0)") dim_length, " ", dim_width
            end if
        else
            write (xml_dims_string, "(I0)") dim_length
        end if

        xml_hdf5_path_string(:) = " "
        write (xml_hdf5_path_string, "(A, A, I0, A, I0, A, A, A, A, A, A)") trim(file_stamp_base), &
            "_", output_meta_iteration, "_xdmf.h5:/", output_iteration, "/", subgroup_dname_nz, "/", hdf5_attr_dname_nz, "/", attr_dname_nz

        write(xml_file_id, "(A, A, A, A, A, A, A, A, A, A, A)", advance="no") '<Attribute Name="', attr_cname, &
            '" Center="', xml_center_type_string, '"><DataItem Format="HDF" NumberType="',  trim(xml_number_type_string), &
            '" Dimensions="', trim(xml_dims_string), '">', trim(xml_hdf5_path_string), '</DataItem></Attribute>'
    end subroutine

end module