! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE
		use SFC_edge_traversal
		use SWE_data_types

		use SWE_adapt
		use SWE_initialize_bathymetry
		use SWE_initialize_dofs
		use SWE_displace
		use SWE_output
		use SWE_xml_output
		use SWE_xml_point_output
		use SWE_point_output
#       if defined(_XDMF)
            use HDF5
            use XDMF_output_base_data_types
            use SWE_XDMF_config
            use SWE_XDMF_output
            use SWE_XDMF_output_filter
            use SWE_XDMF_initialize_dofs
            use SWE_XDMF_adapt
#       endif
		use SWE_euler_timestep
#       if defined(_SWE_PATCH)
			use SWE_PATCH
#       endif

#		if defined(_SWE_DG)
			!use SWE_dg_predictor
			use SWE_dg_solver
#		endif

		use Samoa_swe
		
		! No ASAGI -> Artificial scenario selector 
#       if !defined(_ASAGI)
			use SWE_Scenario
#       endif 

		implicit none

		PRIVATE
		PUBLIC t_swe

		type t_swe
			type(t_swe_init_b_traversal)            :: init_b
			type(t_swe_init_dofs_traversal)         :: init_dofs
			type(t_swe_displace_traversal)          :: displace
			type(t_swe_output_traversal)            :: output
			type(t_swe_xml_point_output_traversal)  :: xml_point_output
			type(t_swe_xml_output_traversal)        :: xml_output
			type(t_swe_point_output_traversal)	    :: point_output
#           if defined(_XDMF)
                type(t_swe_xdmf_output_traversal)   :: xdmf_output
                type(t_swe_xdmf_output_filter_traversal) :: xdmf_output_filter
                type(t_swe_xdmf_adaption_traversal)  :: xdmf_adaption
                type(t_swe_xdmf_init_dofs_traversal) :: xdmf_init_dofs
#           endif
			
			type(t_swe_euler_timestep_traversal)    :: euler
			type(t_swe_adaption_traversal)          :: adaption
	 
#			if defined(_SWE_DG)
	 			! type(t_swe_dg_predictor_traversal)      :: dg_predictor
	 			type(t_swe_dg_timestep_traversal)   :: dg_timestep
#			endif
	 
   			contains
	 
			procedure, pass :: create => swe_create
			procedure, pass :: run => swe_run
			procedure, pass :: destroy => swe_destroy
  		end type t_swe
  
	contains

		!> Creates all required runtime objects for the scenario
		subroutine swe_create(swe, grid, l_log, i_asagi_mode)
			class(t_swe), intent(inout)                                 :: swe
			type(t_grid), intent(inout)									:: grid
			logical, intent(in)						                    :: l_log
			integer, intent(in)											:: i_asagi_mode

			!local variables
			character(64)												:: s_log_name, s_date, s_time
			integer                                                     :: i_error
#           if defined(_XDMF)
                integer                      							:: i_xdmffile_stamp_bn_index
#           endif

#			if defined (_SWE_PATCH)
				call SWE_PATCH_geometry%init(_SWE_PATCH_ORDER)
#			endif

			call date_and_time(s_date, s_time)

#           if defined(_MPI)
				call mpi_bcast(s_date, len(s_date), MPI_CHARACTER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
				call mpi_bcast(s_time, len(s_time), MPI_CHARACTER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#           endif

			swe%output%s_file_stamp = trim(cfg%output_dir) // "/swe_" // trim(s_date) // "_" // trim(s_time)
#           if defined(_XDMF)
                call SWE_xdmf_config_load()
                if(cfg%xdmf%l_xdmfcheckpoint) then
                    swe%output%s_file_stamp = cfg%xdmf%s_xdmffile_stamp
                    i_xdmffile_stamp_bn_index = index(cfg%xdmf%s_xdmfinput_file, "_", .true.)
                    swe%xdmf_init_dofs%base%s_file_stamp = trim(cfg%xdmf%s_xdmfinput_file(:i_xdmffile_stamp_bn_index-1))
                    swe%xdmf_init_dofs%base%output_iteration = cfg%xdmf%i_xdmfcheckpoint_iteration
                    swe%xdmf_init_dofs%base%l_load_data = .false.
                    swe%xdmf_output%base%output_iteration = cfg%xdmf%i_xdmfcheckpoint_iteration + 1
                    swe%xdmf_output_filter%base%output_iteration = swe%xdmf_output%base%output_iteration
                else
                    swe%xdmf_output%base%output_iteration = 0
                    swe%xdmf_output_filter%base%output_iteration = 0
                end if
                swe%xdmf_output%base%s_file_stamp = swe%output%s_file_stamp
#           endif
			swe%xml_output%s_file_stamp = swe%output%s_file_stamp
			swe%xml_point_output%s_file_stamp = swe%output%s_file_stamp
			swe%point_output%s_file_stamp = swe%output%s_file_stamp
			s_log_name = trim(swe%output%s_file_stamp) // ".log"

			if (l_log) then
				_log_open_file(s_log_name)
			endif

			call load_scenario(grid)

			call swe%init_b%create()
			call swe%init_dofs%create()
			call swe%displace%create()
			call swe%output%create()
			call swe%xml_output%create()
			call swe%xml_point_output%create()
#           if defined(_XDMF)
                call swe%xdmf_output%create()
                call swe%xdmf_output_filter%create()
                call swe%xdmf_adaption%create()
                call swe%xdmf_init_dofs%create()
#           endif
			call swe%euler%create()
			call swe%adaption%create()
#			if defined(_SWE_DG)
				!call swe%dg_predictor%create()
				call swe%dg_timestep%create()
#			endif
		end subroutine

		subroutine load_scenario(grid)
			type(t_grid), intent(inout)             :: grid

			integer                                 :: i_error

#			if defined(_ASAGI)
				cfg%afh_bathymetry = asagi_grid_create(ASAGI_FLOAT)
				cfg%afh_displacement = asagi_grid_create(ASAGI_FLOAT)

#               if defined(_MPI)
					call asagi_grid_set_comm(cfg%afh_bathymetry, MPI_COMM_WORLD)
					call asagi_grid_set_comm(cfg%afh_displacement, MPI_COMM_WORLD)
#               endif

				call asagi_grid_set_threads(cfg%afh_bathymetry, cfg%i_threads)
				call asagi_grid_set_threads(cfg%afh_displacement, cfg%i_threads)

				!convert ASAGI mode to ASAGI parameters

				select case(cfg%i_asagi_mode)
					case (0)
						!i_asagi_hints = GRID_NO_HINT
					case (1)
						!i_asagi_hints = ieor(GRID_NOMPI, GRID_PASSTHROUGH)
						call asagi_grid_set_param(cfg%afh_bathymetry, "grid", "pass_through")
						call asagi_grid_set_param(cfg%afh_displacement, "grid", "pass_through")
					case (2)
						!i_asagi_hints = GRID_NOMPI
					case (3)
						!i_asagi_hints = ieor(GRID_NOMPI, SMALL_CACHE)
					case (4)
						!i_asagi_hints = GRID_LARGE_GRID
						call asagi_grid_set_param(cfg%afh_bathymetry, "grid", "cache")
						call asagi_grid_set_param(cfg%afh_displacement, "grid", "cache")
					case default
						try(.false., "Invalid asagi mode, must be in range 0 to 4")
				end select

				!$omp parallel private(i_error), copyin(cfg)
					i_error = asagi_grid_open(cfg%afh_bathymetry,  trim(cfg%s_bathymetry_file), 0); assert_eq(i_error, ASAGI_SUCCESS)
					i_error = asagi_grid_open(cfg%afh_displacement, trim(cfg%s_displacement_file), 0); assert_eq(i_error, ASAGI_SUCCESS)
				!$omp end parallel

				associate(afh_d => cfg%afh_displacement, afh_b => cfg%afh_bathymetry)
					cfg%scaling = max(asagi_grid_max(afh_b, 0) - asagi_grid_min(afh_b, 0), asagi_grid_max(afh_b, 1) - asagi_grid_min(afh_b, 1))

					cfg%offset = [0.5_GRID_SR * (asagi_grid_min(afh_d, 0) + asagi_grid_max(afh_d, 0)), 0.5_GRID_SR * (asagi_grid_min(afh_d, 1) + asagi_grid_max(afh_d, 1))] - 0.5_GRID_SR * cfg%scaling
					cfg%offset = min(max(cfg%offset, [asagi_grid_min(afh_b, 0), asagi_grid_min(afh_b, 1)]), [asagi_grid_max(afh_b, 0), asagi_grid_max(afh_b, 1)] - cfg%scaling)

					if (asagi_grid_dimensions(afh_d) > 2) then
						cfg%dt_eq = asagi_grid_delta(afh_d, 2)
						cfg%t_min_eq = asagi_grid_min(afh_d, 2)
						cfg%t_max_eq = asagi_grid_max(afh_d, 2)
					else
						cfg%dt_eq = 0.0_SR
						cfg%t_min_eq = 0.0_SR
						cfg%t_max_eq = 0.0_SR
					end if

					if (rank_MPI == 0) then
						_log_write(1, '(" SWE: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "]")') &
							trim(cfg%s_bathymetry_file), asagi_grid_min(afh_b, 0), asagi_grid_max(afh_b, 0),  asagi_grid_min(afh_b, 1), asagi_grid_max(afh_b, 1)
						_log_write(1, '(" SWE:  dx: ", F0.2, " dy: ", F0.2)') asagi_grid_delta(afh_b, 0), asagi_grid_delta(afh_b, 1)

						!if the data file has more than two dimensions, we assume that it contains time-dependent displacements
						if (asagi_grid_dimensions(afh_d) > 2) then
							_log_write(1, '(" SWE: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "], time: [", F0.2, ", ", F0.2, "]")') &
							trim(cfg%s_displacement_file), asagi_grid_min(afh_d, 0), asagi_grid_max(afh_d, 0),  asagi_grid_min(afh_d, 1), asagi_grid_max(afh_d, 1), asagi_grid_min(afh_d, 2), asagi_grid_max(afh_d, 2)
							_log_write(1, '(" SWE:  dx: ", F0.2, " dy: ", F0.2, " dt: ", F0.2)') asagi_grid_delta(afh_d, 0), asagi_grid_delta(afh_d, 1), asagi_grid_delta(afh_d, 2)
						else
							_log_write(1, '(" SWE: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "]")') &
							trim(cfg%s_displacement_file), asagi_grid_min(afh_d, 0), asagi_grid_max(afh_d, 0),  asagi_grid_min(afh_d, 1), asagi_grid_max(afh_d, 1)
							_log_write(1, '(" SWE:  dx: ", F0.2, " dy: ", F0.2)') asagi_grid_delta(afh_d, 0), asagi_grid_delta(afh_d, 1)
						end if

						_log_write(1, '(" SWE: computational domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "]")'), cfg%offset(1), cfg%offset(1) + cfg%scaling, cfg%offset(2), cfg%offset(2) + cfg%scaling
					end if
			   end associate
#           else
				cfg%scaling = SWE_Scenario_get_scaling()
				cfg%offset = SWE_Scenario_get_offset()
#			endif
		end subroutine

		!> Destroys all required runtime objects for the scenario
		subroutine swe_destroy(swe, grid, l_log)
			class(t_swe), intent(inout)     :: swe
			type(t_grid), intent(inout)     :: grid
			logical, intent(in)		        :: l_log

			call swe%init_b%destroy()
			call swe%init_dofs%destroy()
			call swe%displace%destroy()
			call swe%output%destroy()
			call swe%xml_output%destroy()
			call swe%xml_point_output%destroy()
			call swe%point_output%destroy()
#           if defined(_XDMF)
                call swe%xdmf_output%destroy()
                call swe%xdmf_output_filter%destroy()
                call swe%xdmf_adaption%destroy()
                call swe%xdmf_init_dofs%destroy()
				call swe_xdmf_param_cells%deallocate()
#				if defined(_SWE_PATCH)
					call swe_xdmf_param_patches%deallocate()
#				endif
#           endif
			call swe%euler%destroy()
			call swe%adaption%destroy()

#			if defined(_SWE_DG)
				call swe%dg_timestep%destroy()
#			endif

#			if defined(_ASAGI)
				call asagi_grid_close(cfg%afh_displacement)
				call asagi_grid_close(cfg%afh_bathymetry)
#			endif

			if (l_log) then
				_log_close_file()
			endif
		end subroutine

		!*********************************
		! run()-method
		!*********************************

		!> Sets the initial values of the SWE and runs the time steps
		subroutine swe_run(swe, grid)
			class(t_swe), intent(inout)                                 :: swe
			type(t_grid), intent(inout)									:: grid

			real (kind = GRID_SR)										:: r_time_next_output
			type(t_grid_info)           	                            :: grid_info, grid_info_max
			integer (kind = GRID_SI)                                    :: i_initial_step, i_time_step
			integer  (kind = GRID_SI)                                   :: i_stats_phase
#           if defined(_XDMF)
                integer                                                 :: hdf5_error
#           endif

			!init parameters
			r_time_next_output = 0.0_GRID_SR

			if (rank_MPI == 0) then
                !$omp master
#           		if defined(_XDMF)
                		if (.not. cfg%xdmf%l_xdmfcheckpoint) then
#           		endif
                    		_log_write(0, *) "SWE: setting initial values and a priori refinement.."
#           		if defined(_XDMF)
                		else
                    		_log_write(0, *) "SWE: reconstructing initial refinement from XDMF.."
                		end if
#           		endif
                	_log_write(0, *) ""
                !$omp end master
            end if
#           if defined(_XDMF)
                call h5open_f(hdf5_error)
#           endif

			call update_stats(swe, grid)
			i_stats_phase = 0
			i_initial_step = 0

#           if defined(_XDMF)
				if (.not. cfg%xdmf%l_xdmfcheckpoint) then
#           endif
				!initialize the bathymetry
				call swe%init_b%traverse(grid)
#           if defined(_XDMF)
				end if
#           endif

			!initialize dofs and set refinement conditions
			do
#               if defined(_XDMF)
                    if (.not. cfg%xdmf%l_xdmfcheckpoint) then
#               endif
                        call swe%init_dofs%traverse(grid)
#               if defined(_XDMF)
                    else
                        call swe%xdmf_init_dofs%traverse(grid)
                    end if
#               endif
				

				grid_info%i_cells = grid%get_cells(MPI_SUM, .true.)
				if (rank_MPI == 0) then
					!$omp master
#                 		if defined(_SWE_PATCH)
							_log_write(1, "(A, I0, A, I0, A, I0, A)") " SWE: ", i_initial_step, " adaptions, ", grid_info%i_cells, " patches = ", grid_info%i_cells * _SWE_PATCH_ORDER_SQUARE, " cells"
#                   	else
							_log_write(1, "(A, I0, A, I0, A)") " SWE: ", i_initial_step, " adaptions, ", grid_info%i_cells, " cells"
#                   	endif
					!$omp end master
				end if
				

#               if defined(_XDMF)
                    if (.not. cfg%xdmf%l_xdmfcheckpoint) then
#               endif
                        if (swe%init_dofs%i_refinements_issued .le. 0) then
                            exit
                        end if
#               if defined(_XDMF)
                    else
                        if (swe%xdmf_init_dofs%base%i_refinements_issued .le. 0) then
                            exit
                        end if
                    end if
#               endif
				
#               if defined(_XDMF)
                    if (.not. cfg%xdmf%l_xdmfcheckpoint) then
#               endif
                        call swe%adaption%traverse(grid)
#               if defined(_XDMF)
                    else
                        call swe%xdmf_adaption%traverse(grid)
                    end if
#               endif

#               if defined(_XDMF)
                    if (.not. cfg%xdmf%l_xdmfcheckpoint) then
#               endif
                        !output grids during initial phase if and only if t_out is 0
                        if (cfg%r_output_time_step == 0.0_GRID_SR) then
							if (cfg%l_pointoutput) then
								call swe%point_output%traverse(grid)
							end if
							
							if(cfg%l_gridoutput) then
								call swe%xml_output%traverse(grid)
							end if

							if(cfg%l_xml_pointoutput) then
								call swe%xml_point_output%traverse(grid)
							end if

#                           if defined(_XDMF)
                                if(cfg%xdmf%l_xdmfoutput) then
                                    call xdmf_output_wrapper(swe, grid, 0)
                                end if
#                           endif

                            r_time_next_output = r_time_next_output + cfg%r_output_time_step
                        end if
#               if defined(_XDMF)
                    !do not output grid when reloading checkpoint
                    end if
#               endif

				i_initial_step = i_initial_step + 1
			end do
			
#           if defined(_XDMF)
				if (cfg%xdmf%l_xdmfcheckpoint) then
					if (rank_MPI == 0) then
						!$omp master
							_log_write(0, *) "XDMF: loading values"
						!$omp end master
					end if
					swe%xdmf_init_dofs%base%l_load_data = .true.
					call swe%xdmf_init_dofs%traverse(grid)
				end if
#           endif   

			grid_info = grid%get_info(MPI_SUM, .true.)

			if (rank_MPI == 0) then
				!$omp master
				_log_write(0, *) "SWE: done."
				_log_write(0, *) ""
				call grid_info%print()
				!$omp end master
			end if

			!output initial grid
#           if defined(_XDMF)
				if (.not. cfg%xdmf%l_xdmfcheckpoint) then
#           endif
					if (cfg%i_output_time_steps > 0 .or. cfg%r_output_time_step >= 0.0_GRID_SR) then
						if (cfg%l_pointoutput) then
							call swe%point_output%traverse(grid)
						end if
						
						if(cfg%l_gridoutput) then
							call swe%xml_output%traverse(grid)
						end if

						if(cfg%l_xml_pointoutput) then
							call swe%xml_point_output%traverse(grid)
						end if

#                       if defined(_XDMF)
							if(cfg%xdmf%l_xdmfoutput) then
								call xdmf_output_wrapper(swe, grid, 0)
							end if
#                       endif

						r_time_next_output = r_time_next_output + cfg%r_output_time_step
					end if
#          	if defined(_XDMF)
				!do not output grid when reloading checkpoint
				end if
#           endif

			!print initial stats
			if (cfg%i_stats_phases >= 0) then
				call update_stats(swe, grid)

				i_stats_phase = i_stats_phase + 1
			end if

			!$omp master
#           	if defined(_XDMF)
					if (.not. cfg%xdmf%l_xdmfcheckpoint) then
#           	endif
						call swe%init_dofs%reduce_stats(MPI_SUM, .true.)
						call swe%adaption%reduce_stats(MPI_SUM, .true.)
#           	if defined(_XDMF)
					else
						call swe%xdmf_init_dofs%reduce_stats(MPI_SUM, .true.)
						call swe%xdmf_adaption%reduce_stats(MPI_SUM, .true.)
					end if
#           	endif
				call grid%reduce_stats(MPI_SUM, .true.)

				if (rank_MPI == 0) then
#           		if defined(_XDMF)
						if (.not. cfg%xdmf%l_xdmfcheckpoint) then
#           		endif
							_log_write(0, *) "SWE: running time steps.."
#           		if defined(_XDMF)
						else
							_log_write(0, *) "SWE: continuing time steps.."
						end if
#           		endif
					_log_write(0, *) ""
				end if
			!$omp end master

#           if defined(_XDMF)
                if (.not. cfg%xdmf%l_xdmfcheckpoint) then
#           endif
                    i_time_step = 0
#           if defined(_XDMF)
                else
                    i_time_step = cfg%xdmf%i_xdmfsimulation_iteration
                    grid%r_time = cfg%xdmf%r_xdmftime
                    grid%r_dt = cfg%xdmf%r_xdmfdeltatime
                end if
#           endif

#           if defined(_ASAGI)
				! during the earthquake, do small time steps that include a displacement

				do
					if ((cfg%r_max_time >= 0.0 .and. grid%r_time >= cfg%r_max_time) .or. (cfg%i_max_time_steps >= 0 .and. i_time_step >= cfg%i_max_time_steps)) then
						exit
					end if

					if (grid%r_time > cfg%t_max_eq) then
						exit
					end if

					i_time_step = i_time_step + 1



#					if defined(_SWE_DG)
						!call swe%dg_predictor%traverse(grid)
#					endif
#					if defined(_SWE_DG)
						!call swe%dg_predictor%traverse(grid)
						!print*,"predictor"
				   		call swe%adaption%traverse(grid)
#					else
						if (cfg%i_adapt_time_steps > 0 .and. mod(i_time_step, cfg%i_adapt_time_steps) == 0) then
							call swe%adaption%traverse(grid)
						end if
#					endif

#					if defined(_SWE_DG)
						!print*,"timestep"
						call swe%dg_timestep%traverse(grid)
#					elif defined(_SWE_PATCH)
						call swe%euler%traverse(grid)
#					endif
					!print*,"displace"
					!displace time-dependent bathymetry
					
					call swe%displace%traverse(grid)

					grid_info%i_cells = grid%get_cells(MPI_SUM, .true.)
					if (rank_MPI == 0) then
						!$omp master
#                       if defined (_SWE_PATCH)
							_log_write(1, '(" SWE: EQ time step: ", I0, ", sim. time:", A, ", dt:", A, ", patches : " I0, " cells: ", I0)') i_time_step, trim(time_to_hrt(DBLE(grid%r_time))), trim(time_to_hrt(DBLE(grid%r_dt))), grid_info%i_cells, grid_info%i_cells * _SWE_PATCH_ORDER_SQUARE
#                       else
							_log_write(1, '(" SWE: EQ time step: ", I0, ", sim. time:", A, ", dt:", A, ", cells: ", I0)') i_time_step, trim(time_to_hrt(DBLE(grid%r_time))), trim(time_to_hrt(DBLE(grid%r_dt))), grid_info%i_cells
#                       endif

						!$omp end master
					end if

					!output grid
					if ((cfg%i_output_time_steps > 0 .and. mod(i_time_step, cfg%i_output_time_steps) == 0) .or. &
						(cfg%r_output_time_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output)) then


						if(cfg%l_gridoutput) then
							call swe%xml_output%traverse(grid)
						end if

						if(cfg%l_xml_pointoutput) then
						   call swe%xml_point_output%traverse(grid)
						end if
						
						if (cfg%l_pointoutput) then
							call swe%point_output%traverse(grid)
						end if

#                       if defined(_XDMF)
							if(cfg%xdmf%l_xdmfoutput) then
								call xdmf_output_wrapper(swe, grid, i_time_step)
							end if
#                       endif

						r_time_next_output = r_time_next_output + cfg%r_output_time_step
					end if
				end do

				!print EQ phase stats
				if (cfg%i_stats_phases >= 0) then
					call update_stats(swe, grid)
				end if
#           endif

			!regular tsunami time steps begin after the earthquake is over

				do
					if ((cfg%r_max_time >= 0.0 .and. grid%r_time >= cfg%r_max_time) .or. (cfg%i_max_time_steps >= 0 .and. i_time_step >= cfg%i_max_time_steps))  then
						!print final solution only if not written directly before
						if (abs((r_time_next_output - cfg%r_output_time_step) - grid%r_time) >= 10e-6) then

								if(cfg%l_gridoutput) then
									call swe%xml_output%traverse(grid)
								end if

								if(cfg%l_xml_pointoutput) then
									call swe%xml_point_output%traverse(grid)
								end if

								if (cfg%l_pointoutput) then
									call swe%point_output%traverse(grid)
								end if

#                       		if defined(_XDMF)
									if(cfg%xdmf%l_xdmfoutput) then
										call xdmf_output_wrapper(swe, grid, i_time_step)
									end if
#                       		endif
						end if
						exit
					end if

					! if ((cfg%r_max_time >= 0.0 .and. grid%r_time >= cfg%r_max_time) .or. (cfg%i_max_time_steps >= 0 .and. i_time_step >= cfg%i_max_time_steps)) then
					!    !print final solution

					!      if(cfg%l_gridoutput) then
					!          call swe%xml_output%traverse(grid)
					!      end if

					!      if(cfg%l_xml_pointoutput) then
					!         call swe%xml_point_output%traverse(grid)
					!      end if

					!      if (cfg%l_pointoutput) then
					!          call swe%point_output%traverse(grid)
					!      end if

					!    exit

						
					! end if
					
					i_time_step = i_time_step + 1
				   
				   
#					if defined(_SWE_DG)
						!call swe%dg_predictor%traverse(grid)
				   		call swe%adaption%traverse(grid)
#					else
						if (cfg%i_adapt_time_steps > 0 .and. mod(i_time_step, cfg%i_adapt_time_steps) == 0) then
							call swe%adaption%traverse(grid)
						end if
#					endif

      !$omp barrier
      call swe%dg_timestep%traverse(grid)

      grid_info%i_cells = grid%get_cells(MPI_SUM, .true.)
					if (rank_MPI == 0) then
						!$omp master
#                       	if defined (_SWE_PATCH)
								_log_write(1, '(" SWE: time step: ", I0, ", sim. time:", A, ", dt:", A, ", patches : " I0, " cells: ", I0)') i_time_step, trim(time_to_hrt(DBLE(grid%r_time))), trim(time_to_hrt(DBLE(grid%r_dt))), grid_info%i_cells, grid_info%i_cells * _SWE_PATCH_ORDER_SQUARE
#                       	else
								_log_write(1, '(" SWE: time step: ", I0, ", sim. time:", A, ", dt:", A, ", cells: ", I0)') i_time_step, trim(time_to_hrt(DBLE(grid%r_time))), trim(time_to_hrt(DBLE(grid%r_dt))), grid_info%i_cells
#                       	endif
						!$omp end master
					end if
				   
					!output grid
					if ((cfg%i_output_time_steps > 0 .and. mod(i_time_step, cfg%i_output_time_steps) == 0) .or. &
							(cfg%r_output_time_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output)) then
						
						
						if (cfg%l_pointoutput) then
							call swe%point_output%traverse(grid)
						end if
						
						if(cfg%l_gridoutput) then
							call swe%xml_output%traverse(grid)
						end if

						if(cfg%l_xml_pointoutput) then
							call swe%xml_point_output%traverse(grid)
						end if
						
#                       if defined(_XDMF)
							if(cfg%xdmf%l_xdmfoutput) then
								call xdmf_output_wrapper(swe, grid, i_time_step)
							end if
#                       endif

						r_time_next_output = r_time_next_output + cfg%r_output_time_step
					end if
				   
				   
				   !print stats
				if ((cfg%r_max_time >= 0.0d0 .and. grid%r_time * cfg%i_stats_phases >= i_stats_phase * cfg%r_max_time) .or. &
					(cfg%i_max_time_steps >= 0 .and. i_time_step * cfg%i_stats_phases >= i_stats_phase * cfg%i_max_time_steps)) then
					call update_stats(swe, grid)

					i_stats_phase = i_stats_phase + 1
				end if
			 end do

			grid_info = grid%get_info(MPI_SUM, .true.)
			grid_info_max = grid%get_info(MPI_MAX, .true.)

			!$omp master
			if (rank_MPI == 0) then
				_log_write(0, '(" SWE: done.")')
				_log_write(0, '()')
				_log_write(0, '("  Cells: avg: ", I0, " max: ", I0)') grid_info%i_cells / (omp_get_max_threads() * size_MPI), grid_info_max%i_cells
				_log_write(0, '()')

				call grid_info%print()
			end if
#           if defined(_XDMF)
                if(swe%xdmf_output%base%root_desc%hdf5_meta_ids%file_id .ne. 0) then
                    call h5fclose_f(swe%xdmf_output%base%root_desc%hdf5_meta_ids%file_id, hdf5_error)
                end if
                call h5close_f(hdf5_error)
#           endif
			!$omp end master
		end subroutine

  		subroutine update_stats(swe, grid)
			class(t_swe), intent(inout)   	:: swe
			type(t_grid), intent(inout)     :: grid
			
			double precision, save          :: t_phase = huge(1.0d0)
			
			!$omp master
				!Initially, just start the timer and don't print anything
				if (t_phase < huge(1.0d0)) then
					t_phase = t_phase + get_wtime()
					
					call swe%init_dofs%reduce_stats(MPI_SUM, .true.)
					call swe%displace%reduce_stats(MPI_SUM, .true.)
#					if defined(_SWE_DG)                    
						call swe%dg_timestep%reduce_stats(MPI_SUM, .true.)       
#					else
						call swe%euler%reduce_stats(MPI_SUM, .true.)
#					endif                    
					call swe%adaption%reduce_stats(MPI_SUM, .true.)
#                   if defined(_XDMF)
                        call swe%xdmf_init_dofs%reduce_stats(MPI_SUM, .true.)
                        call swe%xdmf_output%reduce_stats(MPI_SUM, .true.)
                        call swe%xdmf_output_filter%reduce_stats(MPI_SUM, .true.)
#                   endif
					call grid%reduce_stats(MPI_SUM, .true.)
	   
					if (rank_MPI == 0) then
						_log_write(0, *) ""
						_log_write(0, *) "Phase statistics:"
						_log_write(0, *) ""
						_log_write(0, '(A, T34, A)') " Init: ", trim(swe%init_dofs%stats%to_string())
						_log_write(0, '(A, T34, A)') " Displace: ", trim(swe%displace%stats%to_string())
#						if defined(_SWE_DG)
							_log_write(0, '(A, T34, A)') " Time steps solver: ", trim(swe%dg_timestep%stats%to_string())
#						else                        
							_log_write(0, '(A, T34, A)') " Time steps: ", trim(swe%euler%stats%to_string())
#						endif                        

#						if defined(_SWE_DG)
							_log_write(0, '(A, T34, A)') "Predictor and Adaptions: ", trim(swe%adaption%stats%to_string())
#						else
							_log_write(0, '(A, T34, A)') " Adaptions: ", trim(swe%adaption%stats%to_string())
#						endif
#                   	if defined(_XDMF)
							if(cfg%xdmf%l_xdmfcheckpoint) then
								_log_write(0, '(A, T34, A)') " XDMF Input: ", trim(swe%xdmf_init_dofs%stats%to_string())
							end if
							if(cfg%xdmf%l_xdmfoutput) then
								_log_write(0, '(A, T34, A)') " XDMF Output: ", trim(swe%xdmf_output%stats%to_string())
								if (cfg%xdmf%i_xdmffilter_index .ne. 0) then
									_log_write(0, '(A, T34, A)') " XDMF Output filter: ", trim(swe%xdmf_output_filter%stats%to_string())
								end if
							end if
#                   	endif
						_log_write(0, '(A, T34, A)') " Grid: ", trim(grid%stats%to_string())
						! throughput calculations are a bit different if using patches
#                       if defined(_SWE_PATCH)                        
							_log_write(0, '(A, T34, F12.4, A)') " Element throughput: ", 1.0d-6 * dble(grid%stats%get_counter(traversed_cells)) * (_SWE_PATCH_ORDER_SQUARE)  / t_phase, " M/s"
							_log_write(0, '(A, T34, F12.4, A)') " Memory throughput: ", dble(grid%stats%get_counter(traversed_memory)) * (_SWE_PATCH_ORDER_SQUARE)  / ((1024 * 1024 * 1024) * t_phase), " GB/s"
#							if defined(_SWE_DG)
								_log_write(0, '(A, T34, F12.4, A)') " Cell update throughput solver: ", 1.0d-6 * dble(swe%dg_timestep%stats%get_counter(traversed_cells)) * (_SWE_PATCH_ORDER_SQUARE)  / t_phase, " M/s"
								_log_write(0, '(A, T34, F12.4, A)') " Total cell updates solver: ", 1.0d-6 * dble(swe%dg_timestep%stats%get_counter(traversed_cells)) * (_SWE_PATCH_ORDER_SQUARE), " millions"
								_log_write(0, '(A, T34, F12.4, A)') " Flux solver throughput solver: ", 1.0d-6 * dble(swe%dg_timestep%stats%get_counter(traversed_edges)) * (_SWE_PATCH_NUM_EDGES)   / t_phase, " M/s"
#							else                            
								_log_write(0, '(A, T34, F12.4, A)') " Cell update throughput: ", 1.0d-6 * dble(swe%euler%stats%get_counter(traversed_cells)) * (_SWE_PATCH_ORDER_SQUARE)  / t_phase, " M/s"
								_log_write(0, '(A, T34, F12.4, A)') " Total cell updates: ", 1.0d-6 * dble(swe%euler%stats%get_counter(traversed_cells)) * (_SWE_PATCH_ORDER_SQUARE), " millions"
								_log_write(0, '(A, T34, F12.4, A)') " Flux solver throughput: ", 1.0d-6 * dble(swe%euler%stats%get_counter(traversed_edges)) * (_SWE_PATCH_NUM_EDGES)   / t_phase, " M/s"
#							endif                            
#                       else
							_log_write(0, '(A, T34, F12.4, A)') " Element throughput: ", 1.0d-6 * dble(grid%stats%get_counter(traversed_cells)) / t_phase, " M/s"
							_log_write(0, '(A, T34, F12.4, A)') " Memory throughput: ", dble(grid%stats%get_counter(traversed_memory)) / ((1024 * 1024 * 1024) * t_phase), " GB/s"
							_log_write(0, '(A, T34, F12.4, A)') " Cell update throughput: ", 1.0d-6 * dble(swe%euler%stats%get_counter(traversed_cells)) / t_phase, " M/s"
							_log_write(0, '(A, T34, F12.4, A)') " Total cell updates: ", 1.0d-6 * dble(swe%euler%stats%get_counter(traversed_cells)), " millions"
							_log_write(0, '(A, T34, F12.4, A)') " Flux solver throughput: ", 1.0d-6 * dble(swe%euler%stats%get_counter(traversed_edges)) / t_phase, " M/s"
#                       endif
						_log_write(0, '(A, T34, F12.4, A)') " Asagi time:", grid%stats%get_time(asagi_time), " s"
						_log_write(0, '(A, T34, F12.4, A)') " Phase time:", t_phase, " s"
						_log_write(0, *) ""
					end if
				end if

				call swe%init_dofs%clear_stats()
				call swe%displace%clear_stats()
#				if defined(_SWE_DG)
					call swe%dg_timestep%clear_stats()
#				else
					call swe%euler%clear_stats()
#				endif
#           	if defined(_XDMF)
					call swe%xdmf_init_dofs%clear_stats()
					call swe%xdmf_output%clear_stats()
					call swe%xdmf_output_filter%clear_stats()
#           	endif               

				call swe%adaption%clear_stats()
				call grid%clear_stats()

				t_phase = -get_wtime()
			!$omp end master
		end subroutine

#       if defined(_XDMF)
            subroutine xdmf_output_wrapper(swe, grid, time_step)
                class(t_swe), intent(inout)                               :: swe
                type(t_grid), intent(inout)									:: grid
                integer (kind = GRID_SI), intent(in)						:: time_step

                integer                                                     :: i

                swe%xdmf_output%base%i_sim_iteration = time_step
                swe%xdmf_output_filter%base%i_sim_iteration = time_step
                if (cfg%xdmf%i_xdmffilter_index .ne. 0) then
                    call swe%xdmf_output_filter%traverse(grid)
                end if
                call swe%xdmf_output%traverse(grid)
            end subroutine
#       endif

	END MODULE SWE
#endif
