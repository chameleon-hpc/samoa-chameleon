#include "Compilation_control.f90"
#include "XDMF/XDMF_compilation_control.f90"

! This module provides routines to read metadata from XDMF files
module XDMF_config

    use M_kracken

    use Tools_log
    use Tools_openmp
    use Tools_mpi
    
    use XDMF_data_types
    use XDMF_fox

    use, intrinsic :: iso_fortran_env

    implicit none

    ! Filter parameters
    integer, parameter                      :: xdmf_filter_params_width = 16
    character(len = 1), parameter           :: xdmf_filter_params_delim = " "

    ! Output mode flags
    integer, parameter                      :: xdmf_output_mode_none = 0
    integer, parameter                      :: xdmf_output_mode_cells = 1
    integer, parameter                      :: xdmf_output_mode_patches = 2
    integer, parameter                      :: xdmf_output_mode_all = 3

    ! This structure contains all XDMF-releated configuration variables
    type t_config_xdmf
        logical 				            :: l_xdmfoutput			                            !< xdmf output on/off
        integer                             :: i_xdmfoutput_mode                                !< xdmf output mode flags
        logical                             :: l_xdmfoutput_lagrange                            !< xdmf output use vtk lagrance cells for dg
        character(256) 				        :: s_xdmfinput_file			                        !< xdmf input file
        integer                             :: i_xdmfcpint                                      !< write checkpoint data every n steps
        integer                             :: i_xdmfspf                                        !< amount of steps per file
        logical                             :: l_xdmfcheckpoint                                 !< xdmf is restarting from checkpoint
        logical                             :: l_xdmfcheckpoint_loaded                          !< cached data is loaded
        character(256) 				        :: s_xdmffile_stamp                                 !< cached output file stamp
        integer(INT32)                      :: i_xdmfoutput_iteration                           !< last output iteration
        integer(INT32)                      :: i_xdmfcheckpoint_iteration                       !< last checkpoint iteration
        integer(INT32)                      :: i_xdmfsimulation_iteration                       !< cached simulation iteration
        double precision                    :: r_xdmftime                                       !< cached simulation time
        double precision                    :: r_xdmfdeltatime                                  !< cached simulation time delta
        character(1024)                     :: s_krakencommand                                  !< cached command line for checkpointing
        integer                             :: i_xdmffilter_index                               !< filter function index
        character(256)                      :: s_xdmffilter_params                              !< filter function parameters string
        integer                             :: i_xmdffilter_params_count                        !< filter funstion parameters count
        real, dimension(xdmf_filter_params_width) :: r_xdmffilter_params_vector                 !< filter function parameters vector
    end type

    contains

    ! This routine is meant to be called during initial configuration
    ! It reads the metadata of the last step found in the XDMF file
    subroutine xdmf_load_config(config, result)
        type(t_config_xdmf), intent(inout)  :: config
        integer, intent(inout)              :: result

        character (1024)                    :: fox_attribute_text
        integer                             :: i, j, k
        logical                             :: xdmf_file_exists, is_cp
        type(t_fox_Node), pointer           :: fox_doc, fox_domain, fox_child, fox_attribute_name, fox_attribute_value
        type(t_fox_NodeList), pointer       :: fox_children, fox_grandchildren
        type(t_fox_NamedNodeMap), pointer   :: fox_attributes

        if(config%l_xdmfcheckpoint.and.(.not.config%l_xdmfcheckpoint_loaded)) then
            _log_write(0, ' (" XDMF: Loading checkpoint ", A)') trim(config%s_xdmfinput_file)
            config%l_xdmfcheckpoint_loaded = .true. ! Prevent recursion
            inquire(file=trim(config%s_xdmfinput_file)//char(0), exist=xdmf_file_exists)
            if(xdmf_file_exists) then
                fox_doc => fox_parseFile(trim(config%s_xdmfinput_file)//char(0))
                fox_children => fox_getElementsByTagName(fox_doc, "Domain")
                fox_domain => fox_item(fox_children, 0)

                fox_children => fox_getChildNodes(fox_domain)
                do i = 0, fox_getLength(fox_children) - 1
                    fox_child => fox_item(fox_children, i)

                    ! Read header information
                    if (fox_getLocalName(fox_child) == "Information") then
                        fox_attributes => fox_getAttributes(fox_child)
                        fox_attribute_name => fox_getNamedItem(fox_attributes, "Name")
                        fox_attribute_value => fox_getNamedItem(fox_attributes, "Value")
                        fox_attribute_text = fox_getNodeValue(fox_attribute_value)
                        select case(fox_getNodeValue(fox_attribute_name))
                            case("CommandLine")
                                _log_write(0, ' (" XDMF: Previous command line: ", A)') trim(fox_attribute_text)
                                config%s_krakencommand = trim(fox_attribute_text)
                            case("FileStamp")
                                _log_write(0, ' (" XDMF: Previous file stamp: ", A)') trim(fox_attribute_text)
                                config%s_xdmffile_stamp = trim(fox_attribute_text)
                        end select
                    else if (fox_getLocalName(fox_child) == "Grid") then
                        fox_children => fox_getElementsByTagName(fox_child, "Grid")
                        config%i_xdmfoutput_iteration = fox_getLength(fox_children) - 1

                        ! Search for last checkpoint
                        ol: do k = fox_getLength(fox_children) - 1, 1, -1
                            fox_child => fox_item(fox_children, k)
                            fox_grandchildren => fox_getChildNodes(fox_child)
                            do j = 0, fox_getLength(fox_grandchildren) - 1
                                fox_child => fox_item(fox_grandchildren, j)
                                if (fox_getLocalName(fox_child).eq."Information") then
                                    fox_attributes => fox_getAttributes(fox_child)
                                    fox_attribute_name => fox_getNamedItem(fox_attributes, "Name")
                                    if (fox_getNodeValue(fox_attribute_name).eq."CP") then
                                        call fox_extractDataAttribute(fox_child, "Value", is_cp)
                                        if (is_cp) then
                                            exit ol
                                        end if
                                    end if
                                end if
                            end do
                        end do ol

                        ! Select the last checkpoint
                        config%i_xdmfcheckpoint_iteration = k
                        if(config%i_xdmfcheckpoint_iteration.eq.0) then
                            _log_write(0, ' (" XDMF: ERROR: No last checkpoint step found ")')
                            call exit
                        end if
                        _log_write(0, ' (" XDMF: Previous last checkpoint step: ", I0)') config%i_xdmfcheckpoint_iteration
                        fox_child => fox_item(fox_children, config%i_xdmfcheckpoint_iteration)

                        ! Read the checkpoint
                        fox_children => fox_getChildNodes(fox_child)
                        do j = 0, fox_getLength(fox_children) - 1
                            fox_child => fox_item(fox_children, j)
                            ! Read information of the last step found
                            select case(fox_getLocalName(fox_child))
                                case ("Attribute")
                                    fox_attributes => fox_getAttributes(fox_child)
                                    fox_attribute_name => fox_getNamedItem(fox_attributes, "Name")
                                    fox_child => fox_getFirstChild(fox_child)
                                    select case(fox_getNodeValue(fox_attribute_name))
                                        case("Step")
                                            call fox_extractDataContent(fox_child, config%i_xdmfsimulation_iteration)
                                            _log_write(0, ' (" XDMF: Previous last simulation step: ", I0)') config%i_xdmfsimulation_iteration
                                        case("DeltaTime")
                                            call fox_extractDataContent(fox_child, config%r_xdmfdeltatime)
                                            _log_write(0, ' (" XDMF: Previous delta time: ", F0.8)') config%r_xdmfdeltatime
                                    end select
                                case ("Time")
                                    call fox_extractDataAttribute(fox_child, "Value", config%r_xdmftime)
                                    _log_write(0, ' (" XDMF: Previous last simulation time: ", F0.8)') config%r_xdmftime
                                    if (config%r_xdmftime.le.0) then
                                        _log_write(0, ' (" XDMF: Previous last simulation time is zero or lower. Restarting simulation from beginning.")')
                                        config%l_xdmfcheckpoint = .false.
                                    end if
                            end select
                        end do

                    end if
                end do

                call fox_destroy(fox_doc)
                result = 1
                return
            else
                _log_write(0, ' (" XDMF: Error: Checkpoint file not found ")')
                result = 0
                return
            end if
        end if
        result = -1
    end subroutine

#   if defined(_MPI)
        ! This routine scatters the config data read from a XDMF file to all other ranks
        subroutine xdmf_bcast_config(config)
            type(t_config_xdmf), intent(inout)  :: config

            integer                             :: error

            call mpi_bcast(config%s_krakencommand, 1024, MPI_CHARACTER, 0, MPI_COMM_WORLD, error); assert_eq(error, 0)
            call mpi_bcast(config%s_xdmffile_stamp, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, error); assert_eq(error, 0)
            call mpi_bcast(config%i_xdmfcheckpoint_iteration, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, error); assert_eq(error, 0)
            call mpi_bcast(config%i_xdmfoutput_iteration, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, error); assert_eq(error, 0)
            call mpi_bcast(config%i_xdmfsimulation_iteration, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, error); assert_eq(error, 0)
            call mpi_bcast(config%r_xdmftime, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error); assert_eq(error, 0)
            call mpi_bcast(config%r_xdmfdeltatime, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error); assert_eq(error, 0)
            if(rank_MPI == 0) then
                _log_write(0, ' (" XDMF: Broadcasted and loading previous configuration options ")')
            end if
        end subroutine
#   endif

    ! This routine splits the filter parameters into an array
    subroutine xdmf_splitfilter_config(config)
        type(t_config_xdmf), intent(inout)  :: config

        character(len = 16), dimension(xdmf_filter_params_width) :: xdmf_filter_params_strings
        integer                             :: iicount, ilen, i, error
        integer, dimension(xdmf_filter_params_width) :: ibegin, iterm

        call delim(config%s_xdmffilter_params, xdmf_filter_params_strings, xdmf_filter_params_width, &
            iicount, ibegin, iterm, ilen, xdmf_filter_params_delim)
        config%i_xmdffilter_params_count = iicount

        config%r_xdmffilter_params_vector = 0
        do i = 1, iicount
            call string_to_real(xdmf_filter_params_strings(i), config%r_xdmffilter_params_vector(i), error)
        end do
    end subroutine

    ! This routine prints filter configuration
    subroutine xdmf_filter_write(config)
        type(t_config_xdmf), intent(inout)  :: config

        integer                             :: i
        character(len = 32)                 :: param_string
        character(len = 1024)               :: params_string

        params_string(:) = " "
        do i = 1, config%i_xmdffilter_params_count
            param_string(:) = " "
            write(param_string, "(F0.6)") config%r_xdmffilter_params_vector(i)
            if (i .eq. 1) then
                params_string = trim(param_string)
            else
                params_string = trim(params_string) // ", " // trim(param_string)
            end if
        end do
        if (config%i_xdmffilter_index .ne. 0) then
            _log_write(0, '(" XDMF: output filter index: ", I0, ", ", I0, " parameters: ", A)'), &
                config%i_xdmffilter_index, config%i_xmdffilter_params_count, trim(params_string)
        else if(iand(config%i_xdmfoutput_mode, xdmf_output_mode_all) .eq. xdmf_output_mode_all) then
            _log_write(0, '(" XDMF: output filter used only for hybrid layering")')
        else
            _log_write(0, '(" XDMF: output filter disabled")')
        end if
    end subroutine

end module
