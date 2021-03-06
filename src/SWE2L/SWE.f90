! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE2L)
	MODULE SWE2L
		use SFC_edge_traversal
		use SWE2L_data_types

		use SWE2L_adapt
		use SWE2L_initialize_bathymetry
		use SWE2L_initialize_dofs
		use SWE2L_displace
		use SWE2L_output
		use SWE2L_xml_output
		use SWE2L_ascii_output
		use SWE2L_point_output
		use SWE2L_euler_timestep
#       if defined(_SWE_PATCH)
            use SWE_PATCH
#       endif
		use Samoa_swe2l
		
        ! No ASAGI -> Artificial scenario selector 
#       if !defined(_ASAGI)
            use SWE2L_Scenario
#       endif 

		implicit none

		PRIVATE
		PUBLIC t_swe

		type t_swe
            type(t_swe_init_b_traversal)            :: init_b
            type(t_swe_init_dofs_traversal)         :: init_dofs
            type(t_swe_displace_traversal)          :: displace
            type(t_swe_output_traversal)            :: output
            type(t_swe_xml_output_traversal)        :: xml_output
            type(t_swe_ascii_output_traversal)      :: ascii_output
	        type(t_swe_point_output_traversal)	    :: point_output

            type(t_swe_euler_timestep_traversal)    :: euler
            type(t_swe_adaption_traversal)          :: adaption

            contains

            procedure, pass :: create => swe_create
            procedure, pass :: run => swe_run
            procedure, pass :: destroy => swe_destroy
        end type

		contains

		!> Creates all required runtime objects for the scenario
		subroutine swe_create(swe, grid, l_log, i_asagi_mode)
            class(t_swe), intent(inout)                                 :: swe
			type(t_grid), intent(inout)									:: grid
			logical, intent(in)						                    :: l_log
			integer, intent(in)											:: i_asagi_mode

			!local variables
			character(64)												:: s_log_name, s_date, s_time
			character(10)						    :: s_rank
			integer                                                     :: i_error

#if defined (_SWE_PATCH)
            call SWE_PATCH_geometry%init(_SWE_PATCH_ORDER)
#endif

			call date_and_time(s_date, s_time)
			write (s_rank,'(I0)') rank_MPI

#           if defined(_MPI)
                call mpi_bcast(s_date, len(s_date), MPI_CHARACTER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
                call mpi_bcast(s_time, len(s_time), MPI_CHARACTER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#           endif

            swe%output%s_file_stamp = trim(cfg%output_dir) // "/swe_" // trim(s_date) // "_" // trim(s_time)
			swe%xml_output%s_file_stamp = trim(cfg%output_dir) // "/swe_" // trim(s_date) // "_" // trim(s_time)
            swe%point_output%s_file_stamp = trim(cfg%output_dir) // "/swe_" // trim(s_date) // "_" // trim(s_time)

			if(cfg%l_log_stats_per_process) then
				s_log_name = trim(swe%xml_output%s_file_stamp) // "_R" // trim(s_rank) // ".log"
			else
				s_log_name = trim(swe%xml_output%s_file_stamp) // ".log"
			endif

			if (l_log) then
				_log_open_file(s_log_name)
			endif

			call load_scenario(grid)

			call swe%init_b%create()
			call swe%init_dofs%create()
            call swe%displace%create()
            call swe%output%create()
            call swe%xml_output%create()
            call swe%ascii_output%create()
            call swe%euler%create()
            call swe%adaption%create()
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
                        _log_write(1, '(" SWE2L: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "]")') &
                            trim(cfg%s_bathymetry_file), asagi_grid_min(afh_b, 0), asagi_grid_max(afh_b, 0),  asagi_grid_min(afh_b, 1), asagi_grid_max(afh_b, 1)
                        _log_write(1, '(" SWE2L:  dx: ", F0.2, " dy: ", F0.2)') asagi_grid_delta(afh_b, 0), asagi_grid_delta(afh_b, 1)

                        !if the data file has more than two dimensions, we assume that it contains time-dependent displacements
                        if (asagi_grid_dimensions(afh_d) > 2) then
                            _log_write(1, '(" SWE2L: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "], time: [", F0.2, ", ", F0.2, "]")') &
                            trim(cfg%s_displacement_file), asagi_grid_min(afh_d, 0), asagi_grid_max(afh_d, 0),  asagi_grid_min(afh_d, 1), asagi_grid_max(afh_d, 1), asagi_grid_min(afh_d, 2), asagi_grid_max(afh_d, 2)
                            _log_write(1, '(" SWE2L:  dx: ", F0.2, " dy: ", F0.2, " dt: ", F0.2)') asagi_grid_delta(afh_d, 0), asagi_grid_delta(afh_d, 1), asagi_grid_delta(afh_d, 2)
                        else
                            _log_write(1, '(" SWE2L: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "]")') &
                            trim(cfg%s_displacement_file), asagi_grid_min(afh_d, 0), asagi_grid_max(afh_d, 0),  asagi_grid_min(afh_d, 1), asagi_grid_max(afh_d, 1)
                            _log_write(1, '(" SWE2L:  dx: ", F0.2, " dy: ", F0.2)') asagi_grid_delta(afh_d, 0), asagi_grid_delta(afh_d, 1)
                        end if

                        _log_write(1, '(" SWE2L: computational domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "]")'), cfg%offset(1), cfg%offset(1) + cfg%scaling, cfg%offset(2), cfg%offset(2) + cfg%scaling
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
            call swe%ascii_output%destroy()
            call swe%point_output%destroy()
            call swe%euler%destroy()
            call swe%adaption%destroy()

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

			!init parameters
			r_time_next_output = 0.0_GRID_SR

            if (rank_MPI == 0) then
                !$omp master
                _log_write(0, *) "SWE2L: setting initial values and a priori refinement.."
                _log_write(0, *) ""
                !$omp end master
            end if

            call update_stats(swe, grid)
			i_stats_phase = 0

            i_initial_step = 0

            !initialize the bathymetry
            call swe%init_b%traverse(grid)

			do
				!initialize dofs and set refinement conditions
				call swe%init_dofs%traverse(grid)

				grid_info%i_cells = grid%get_cells(MPI_SUM, .true.)
                if (rank_MPI == 0) then
                    !$omp master
#                   if defined(_SWE_PATCH)
                        _log_write(1, "(A, I0, A, I0, A, I0, A)") " SWE2L: ", i_initial_step, " adaptions, ", grid_info%i_cells, " patches = ", grid_info%i_cells * _SWE_PATCH_ORDER_SQUARE, " cells"
#                   else
                        _log_write(1, "(A, I0, A, I0, A)") " SWE2L: ", i_initial_step, " adaptions, ", grid_info%i_cells, " cells"
#                   endif
                    !$omp end master
                end if

				if (swe%init_dofs%i_refinements_issued .le. 0) then
					exit
				endif

                ! always perform LB in this phase (minimize chance of a process running out of memory)
                grid%i_steps_since_last_LB = cfg%i_lb_frequency - 1
                ! only use homogeneous LB here (same reason)
                grid%l_grid_generation = .true.

				call swe%adaption%traverse(grid)

                !output grids during initial phase if and only if t_out is 0
                if (cfg%r_output_time_step == 0.0_GRID_SR) then
                    if (cfg%l_ascii_output) then
                        call swe%ascii_output%traverse(grid)
                    end if

                    if(cfg%l_gridoutput) then
                        call swe%xml_output%traverse(grid)
                    end if

                    if (cfg%l_pointoutput) then
                        call swe%point_output%traverse(grid)
                    end if

                    r_time_next_output = r_time_next_output + cfg%r_output_time_step
                end if

				i_initial_step = i_initial_step + 1
			end do

			grid%l_grid_generation = .false.

            grid_info = grid%get_info(MPI_SUM, .true.)

            if (rank_MPI == 0) then
                !$omp master
                _log_write(0, *) "SWE2L: done."
                _log_write(0, *) ""

                call grid_info%print()
                !$omp end master
			end if

			!output initial grid
			if (cfg%i_output_time_steps > 0 .or. cfg%r_output_time_step >= 0.0_GRID_SR) then
                if (cfg%l_ascii_output) then
                    call swe%ascii_output%traverse(grid)
                end if

                if(cfg%l_gridoutput) then
                    call swe%xml_output%traverse(grid)
                end if

                if (cfg%l_pointoutput) then
                    call swe%point_output%traverse(grid)
                end if

				r_time_next_output = r_time_next_output + cfg%r_output_time_step
			end if

			!print initial stats
			if (cfg%i_stats_phases >= 0) then
                call update_stats(swe, grid)

                i_stats_phase = i_stats_phase + 1
			end if

            !$omp master
            call swe%init_dofs%reduce_stats(MPI_SUM, .true.)
            call swe%adaption%reduce_stats(MPI_SUM, .true.)
            call grid%reduce_stats(MPI_SUM, .true.)

            if (rank_MPI == 0) then
                _log_write(0, *) "SWE2L: running time steps.."
                _log_write(0, *) ""
			end if
            !$omp end master

            i_time_step = 0

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

                    if (cfg%i_adapt_time_steps > 0 .and. mod(i_time_step, cfg%i_adapt_time_steps) == 0) then
                        !refine grid
#if defined(_IPM)
	 	    call mpi_pcontrol( 1,"adaption_traversal"//char(0))	
#endif
                        call swe%adaption%traverse(grid)
#if defined(_IPM)
	 	    call mpi_pcontrol( -1,"adaption_traversal"//char(0))	
#endif
                    end if

                    !do an euler time step
#if defined(_IPM)
	 	    call mpi_pcontrol( 1,"euler_traversal"//char(0))	
#endif
                    call swe%euler%traverse(grid)
#if defined(_IPM)
		    call mpi_pcontrol( -1,"euler_traversal"//char(0))
#endif

                    !displace time-dependent bathymetry
#if defined(_IPM)
	 	    call mpi_pcontrol( 1,"displace_traversal"//char(0))	
#endif
                    call swe%displace%traverse(grid)
#if defined(_IPM)
	 	    call mpi_pcontrol( -1,"displace_traversal"//char(0))	
#endif

                    grid_info%i_cells = grid%get_cells(MPI_SUM, .true.)
                    if (rank_MPI == 0) then
                        !$omp master
#                       if defined (_SWE_PATCH)
                            _log_write(1, '(" SWE2L: EQ time step: ", I0, ", sim. time:", A, ", dt:", A, ", patches : " I0, " cells: ", I0)') i_time_step, trim(time_to_hrt(DBLE(grid%r_time))), trim(time_to_hrt(DBLE(grid%r_dt))), grid_info%i_cells, grid_info%i_cells * _SWE_PATCH_ORDER_SQUARE
#                       else
                            _log_write(1, '(" SWE2L: EQ time step: ", I0, ", sim. time:", A, ", dt:", A, ", cells: ", I0)') i_time_step, trim(time_to_hrt(DBLE(grid%r_time))), trim(time_to_hrt(DBLE(grid%r_dt))), grid_info%i_cells
#                       endif

                        !$omp end master
                    end if

                    !output grid
                    if ((cfg%i_output_time_steps > 0 .and. mod(i_time_step, cfg%i_output_time_steps) == 0) .or. &
                        (cfg%r_output_time_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output)) then

                        if (cfg%l_ascii_output) then
                            call swe%ascii_output%traverse(grid)
                        end if

                        if(cfg%l_gridoutput) then
                            call swe%xml_output%traverse(grid)
                        end if

                        if (cfg%l_pointoutput) then
                            call swe%point_output%traverse(grid)
                        end if

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
				if ((cfg%r_max_time >= 0.0 .and. grid%r_time >= cfg%r_max_time) .or. (cfg%i_max_time_steps >= 0 .and. i_time_step >= cfg%i_max_time_steps)) then
					exit
				end if

				i_time_step = i_time_step + 1

                if (cfg%i_adapt_time_steps > 0 .and. mod(i_time_step, cfg%i_adapt_time_steps) == 0) then
                    !refine grid
#if defined(_IPM)
	 	    call mpi_pcontrol( 1,"adaption_traversal"//char(0))	
#endif
                    call swe%adaption%traverse(grid)
#if defined(_IPM)
	 	    call mpi_pcontrol( -1,"adaption_traversal"//char(0))	
#endif
                end if

				!do a time step
#if defined(_IPM)
	 	    call mpi_pcontrol( 1,"euler_traversal"//char(0))	
#endif
                    call swe%euler%traverse(grid)
#if defined(_IPM)
		    call mpi_pcontrol( -1,"euler_traversal"//char(0))
#endif

                grid_info%i_cells = grid%get_cells(MPI_SUM, .true.)
                if (rank_MPI == 0) then
                    !$omp master
#                       if defined (_SWE_PATCH)
                            _log_write(1, '(" SWE2L: time step: ", I0, ", sim. time:", A, ", dt:", A, ", patches : " I0, " cells: ", I0)') i_time_step, trim(time_to_hrt(DBLE(grid%r_time))), trim(time_to_hrt(DBLE(grid%r_dt))), grid_info%i_cells, grid_info%i_cells * _SWE_PATCH_ORDER_SQUARE
#                       else
                            _log_write(1, '(" SWE2L: time step: ", I0, ", sim. time:", A, ", dt:", A, ", cells: ", I0)') i_time_step, trim(time_to_hrt(DBLE(grid%r_time))), trim(time_to_hrt(DBLE(grid%r_dt))), grid_info%i_cells
#                       endif
                    !$omp end master
                end if

				!output grid
				if ((cfg%i_output_time_steps > 0 .and. mod(i_time_step, cfg%i_output_time_steps) == 0) .or. &
				    (cfg%r_output_time_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output)) then

                    if (cfg%l_ascii_output) then
             	       call swe%ascii_output%traverse(grid)
               	    end if

                    if(cfg%l_gridoutput) then
                        call swe%xml_output%traverse(grid)
                    end if

                    if (cfg%l_pointoutput) then
                        call swe%point_output%traverse(grid)
                    end if

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
                _log_write(0, '(" SWE2L: done.")')
                _log_write(0, '()')
                _log_write(0, '("  Cells: avg: ", I0, " max: ", I0)') grid_info%i_cells / (omp_get_max_threads() * size_MPI), grid_info_max%i_cells
                _log_write(0, '()')

                call grid_info%print()
            end if
            !$omp end master
		end subroutine

        subroutine update_stats(swe, grid)
            	class(t_swe), intent(inout)   :: swe
    		type(t_grid), intent(inout)     :: grid

		double precision, save          :: t_phase = huge(1.0d0)
                character (len = 1024), save     :: str(5)

			!$omp master
                !Initially, just start the timer and don't print anything
                if (t_phase < huge(1.0d0)) then
                    t_phase = t_phase + get_wtime()
                   
		    if( cfg%l_verbose_stats) then
			    ! CAUTION: this is ugly, but ultimately necessary since threads may not be associated
			    if ( associated(swe%init_dofs%threads) ) then
			    	str(1) = swe%init_dofs%stats%to_verbose_string(swe%init_dofs%threads%stats)
			    else
				call swe%init_dofs%reduce_stats(MPI_SUM, .true.)
				str(1) = swe%init_dofs%stats%to_string()
			    endif

			    if ( associated(swe%displace%threads) ) then
			    	str(2) = swe%displace%stats%to_verbose_string(swe%displace%threads%stats)
			    else
				call swe%displace%reduce_stats(MPI_SUM, .true.)
				str(2) = swe%displace%stats%to_string()
			    endif	

			    if ( associated(swe%euler%threads) ) then
			    	str(3) = swe%euler%stats%to_verbose_string(swe%euler%threads%stats)
			    else
				call swe%euler%reduce_stats(MPI_SUM, .true.)
				str(3) = swe%euler%stats%to_string()
			    endif 			    

			    if ( associated(swe%adaption%threads) ) then
			    	str(4) = swe%adaption%stats%to_verbose_string(swe%adaption%threads%stats)
			    else
				call swe%adaption%reduce_stats(MPI_SUM, .true.)
				str(4) = swe%adaption%stats%to_string()
			    endif	

			    if (grid%threads%get_size()>0 ) then
			    	str(5) = grid%stats%to_verbose_string(grid%threads%elements%stats)
			    else
				call grid%reduce_stats(MPI_SUM, .true.)
				str(5) = grid%stats%to_string()
			    endif 

			    if (rank_MPI == 0 .or. cfg%l_log_stats_per_process) then
		                _log_write(0, *) ""
		                _log_write(0, *) "Phase statistics:"
		                _log_write(0, *) ""
				_log_write(0, *) "Format: (tsum, tmin, tmax), reduced over all threads"
		                _log_write(0, '(A, T34, A)') " Init: ", trim(str(1))
		                _log_write(0, '(A, T34, A)') " Displace: ", trim(str(2))
		                _log_write(0, '(A, T34, A)') " Time steps: ", trim(str(3))
		                _log_write(0, '(A, T34, A)') " Adaptions: ", trim(str(4))
		                _log_write(0, '(A, T34, A)') " Grid: ", trim(str(5))
		                _log_write(0, '(A, T34, F12.4, A)') " Phase time:", t_phase, " s"
		                _log_write(0, *) ""
			    endif	

		    else
			    if (.not. cfg%l_log_stats_per_process) then	
		                call swe%init_dofs%reduce_stats(MPI_SUM, .true.)
		                call swe%displace%reduce_stats(MPI_SUM, .true.)
		                call swe%euler%reduce_stats(MPI_SUM, .true.)
		                call swe%adaption%reduce_stats(MPI_SUM, .true.)
		                call grid%reduce_stats(MPI_SUM, .true.)
			    endif

		            if (rank_MPI == 0 .or. cfg%l_log_stats_per_process) then
		                _log_write(0, *) ""
		                _log_write(0, *) "Phase statistics:"
		                _log_write(0, *) ""
		                _log_write(0, '(A, T34, A)') " Init: ", trim(swe%init_dofs%stats%to_string())
		                _log_write(0, '(A, T34, A)') " Displace: ", trim(swe%displace%stats%to_string())
		                _log_write(0, '(A, T34, A)') " Time steps: ", trim(swe%euler%stats%to_string())
		                _log_write(0, '(A, T34, A)') " Adaptions: ", trim(swe%adaption%stats%to_string())
		                _log_write(0, '(A, T34, A)') " Grid: ", trim(grid%stats%to_string())
		                ! throughput calculations are a bit different if using patches
#                       if defined(_SWE_PATCH)                        
		                    _log_write(0, '(A, T34, F12.4, A)') " Element throughput: ", 1.0d-6 * dble(grid%stats%get_counter(traversed_cells)) * (_SWE_PATCH_ORDER_SQUARE)  / t_phase, " M/s"
		                    _log_write(0, '(A, T34, F12.4, A)') " Memory throughput: ", dble(grid%stats%get_counter(traversed_memory)) / ((1024 * 1024 * 1024) * t_phase), " GB/s"
		                    _log_write(0, '(A, T34, F12.4, A)') " Cell update throughput: ", 1.0d-6 * dble(swe%euler%stats%get_counter(traversed_cells)) * (_SWE_PATCH_ORDER_SQUARE)  / t_phase, " M/s"
		                    _log_write(0, '(A, T34, F12.4, A)') " Total cell updates: ", 1.0d-6 * dble(swe%euler%stats%get_counter(traversed_cells)) * (_SWE_PATCH_ORDER_SQUARE), " millions"
		                    _log_write(0, '(A, T34, F12.4, A)') " Flux solver throughput: ", 1.0d-6 * dble(swe%euler%stats%get_counter(traversed_edges)) * (_SWE_PATCH_NUM_EDGES)   / t_phase, " M/s"
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
                end if

                call swe%init_dofs%clear_stats()
                call swe%displace%clear_stats()
                call swe%euler%clear_stats()
                call swe%adaption%clear_stats()
                call grid%clear_stats()

                t_phase = -get_wtime()
            !$omp end master
        end subroutine


	END MODULE SWE2L
#endif
