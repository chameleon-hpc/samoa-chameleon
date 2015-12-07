! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE
		use SFC_edge_traversal
		use SWE_data_types

		use SWE_adapt
		use SWE_initialize
		use SWE_displace
		use SWE_output
		use SWE_xml_output
		use SWE_ascii_output
		use SWE_point_output
		use SWE_euler_timestep
#if defined(_SWE_SIMD)
		use SWE_SIMD
#endif

		use Samoa_swe

		implicit none

		PRIVATE
		PUBLIC t_swe

		type t_swe
            type(t_swe_init_traversal)              :: init
            type(t_swe_displace_traversal)          :: displace
            type(t_swe_output_traversal)            :: output
            type(t_swe_xml_output_traversal)        :: xml_output
            type(t_swe_ascii_output_traversal)      :: ascii_output                     !-------------------------
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
			type(t_grid), intent(inout)			:: grid
			logical, intent(in)		                :: l_log
			integer, intent(in)				:: i_asagi_mode

			!local variables
			character(64)					:: s_log_name, s_date, s_time
			integer                                         :: i_error

			call date_and_time(s_date, s_time)
			
#if defined (_SWE_SIMD)
			call SWE_SIMD_geometry%init(_SWE_SIMD_ORDER)
#endif

#           if defined(_MPI)
                call mpi_bcast(s_date, len(s_date), MPI_CHARACTER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
                call mpi_bcast(s_time, len(s_time), MPI_CHARACTER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
#           endif

			write (swe%output%s_file_stamp, "(A, A, A8, A, A6)") "output/swe", "_", s_date, "_", s_time
			write (swe%xml_output%s_file_stamp, "(A, A, A8, A, A6)") "output/swe", "_", s_date, "_", s_time
            write (swe%point_output%s_file_stamp, "(A, A, A8, A, A6)") "output/swe", "_", s_date, "_", s_time
			write (s_log_name, '(A, A)') TRIM(swe%xml_output%s_file_stamp), ".log"

			if (l_log) then
				_log_open_file(s_log_name)
			endif

			call load_scenario(grid, cfg%s_bathymetry_file, cfg%s_displacement_file)

			call swe%init%create()
            call swe%displace%create()
            call swe%output%create()
            call swe%xml_output%create()
            call swe%ascii_output%create()
            call swe%euler%create()
            call swe%adaption%create()
		end subroutine

		subroutine load_scenario(grid, ncd_bath, ncd_displ, scaling, offset)
			type(t_grid), intent(inout)             :: grid
            character(*), intent(in)                :: ncd_bath, ncd_displ
            double precision, optional,intent(in)   :: scaling, offset(2)

			integer						            :: i_asagi_hints
			integer                                 :: i_error

#			if defined(_ASAGI)
                select case(cfg%i_asagi_mode)
                    case (0)
                        !i_asagi_hints = GRID_NO_HINT   ! OLD VERSION OF ASAGI! REMOVE?
                    case (1)
                        !i_asagi_hints = ieor(GRID_NOMPI, GRID_PASSTHROUGH)   ! OLD VERSION OF ASAGI! REMOVE?
                    case (2)
                        !i_asagi_hints = GRID_NOMPI   ! OLD VERSION OF ASAGI! REMOVE?
                    case (3)
                        !i_asagi_hints = ieor(GRID_NOMPI, SMALL_CACHE)   ! OLD VERSION OF ASAGI! REMOVE?
                    case (4)
                        !i_asagi_hints = GRID_LARGE_GRID   ! OLD VERSION OF ASAGI! REMOVE?
                    case default
                        try(.false., "Invalid asagi mode, must be in range 0 to 4")
                end select

#               if defined(_ASAGI_NUMA)
                    cfg%afh_bathymetry = grid_create_for_numa(grid_type = ASAGI_FLOAT, hint = i_asagi_hints, levels = 1, tcount=omp_get_max_threads())
                    cfg%afh_displacement = grid_create_for_numa(grid_type = ASAGI_FLOAT, hint = ior(i_asagi_hints, GRID_HAS_TIME), levels = 1, tcount=omp_get_max_threads())

                    !$omp parallel private(i_error)
						i_error = grid_register_thread(cfg%afh_bathymetry); assert_eq(i_error, GRID_SUCCESS)
						i_error = grid_register_thread(cfg%afh_displacement); assert_eq(i_error, GRID_SUCCESS)
                        i_error = asagi_open(cfg%afh_bathymetry, trim(ncd_bath), 0); assert_eq(i_error, GRID_SUCCESS)
                        i_error = asagi_open(cfg%afh_displacement, trim(ncd_displ), 0); assert_eq(i_error, GRID_SUCCESS)
                    !$omp end parallel
#               else
                    cfg%afh_bathymetry = asagi_create(type = ASAGI_FLOAT)
                    cfg%afh_displacement = asagi_create(type = ASAGI_FLOAT)

                    i_error = asagi_open(cfg%afh_bathymetry, trim(ncd_bath), 0); assert_eq(i_error, 0) ! 0 --> GRID_SUCCESS
                    i_error = asagi_open(cfg%afh_displacement, trim(ncd_displ), 0); assert_eq(i_error, 0) ! 0 --> GRID_SUCCESS
#               endif

                associate(afh_d => cfg%afh_displacement, afh_b => cfg%afh_bathymetry)
                
                    if (present(scaling)) then
                    else
                        cfg%scaling = max(asagi_grid_max(afh_b,0) - asagi_grid_min(afh_b,0), asagi_grid_max(afh_b,1) - asagi_grid_min(afh_b,1))
                    end if

                    if (present(offset)) then
                    else
                        cfg%offset = [0.5_GRID_SR * (asagi_grid_min(afh_d,0) + asagi_grid_max(afh_d,0)), 0.5_GRID_SR * (asagi_grid_min(afh_d,1) + asagi_grid_max(afh_d,1))] - 0.5_GRID_SR * cfg%scaling
                        cfg%offset = min(max(cfg%offset, [asagi_grid_min(afh_b,0), asagi_grid_min(afh_b,1)]), [asagi_grid_max(afh_b,0), asagi_grid_max(afh_b,1)] - cfg%scaling)
                    end if

                    if (rank_MPI == 0) then
                        _log_write(1, '(" SWE: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "], time: [", F0.2, ", ", F0.2, "]")') &
                            trim(ncd_bath), asagi_grid_min(afh_b,0), asagi_grid_max(afh_b,0),  asagi_grid_min(afh_b,1), asagi_grid_max(afh_b,1),  asagi_grid_min(afh_b,2), asagi_grid_max(afh_b,2)

                        _log_write(1, '(" SWE: loaded ", A, ", domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "], time: [", F0.2, ", ", F0.2, "]")') &
                            trim(ncd_displ), asagi_grid_min(afh_d,0), asagi_grid_max(afh_d,0),  asagi_grid_min(afh_d,1), asagi_grid_max(afh_d,1),  asagi_grid_min(afh_d,2), asagi_grid_max(afh_d,2)

                        _log_write(1, '(" SWE: computational domain: [", F0.2, ", ", F0.2, "] x [", F0.2, ", ", F0.2, "]")'), cfg%offset(1), cfg%offset(1) + cfg%scaling, cfg%offset(2), cfg%offset(2) + cfg%scaling
                    end if
               end associate
#           else
                cfg%scaling = 1.0_GRID_SR
                cfg%offset = [0.0_GRID_SR, 0.0_GRID_SR]
#			endif
		end subroutine

		!> Destroys all required runtime objects for the scenario
		subroutine swe_destroy(swe, grid, l_log)
            class(t_swe), intent(inout)     :: swe
			type(t_grid), intent(inout)     :: grid
			logical, intent(in)		        :: l_log

			call swe%init%destroy()
            call swe%displace%destroy()
            call swe%output%destroy()
            call swe%xml_output%destroy()
            call swe%ascii_output%destroy()
            call swe%point_output%destroy()
            call swe%euler%destroy()
            call swe%adaption%destroy()

#			if defined(_ASAGI)
				call asagi_close(cfg%afh_displacement)
				call asagi_close(cfg%afh_bathymetry)
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
                _log_write(0, *) "SWE: setting initial values and a priori refinement.."
                _log_write(0, *) ""
                !$omp end master
            end if

            call update_stats(swe, grid)
			i_stats_phase = 0

            i_initial_step = 0

			do
				!set numerics and check for refinement
				call swe%init%traverse(grid)

                if (rank_MPI == 0) then
                    grid_info%i_cells = grid%get_cells(MPI_SUM, .false.)

                    !$omp master
                    _log_write(1, "(A, I0, A, I0, A)") " SWE: ", i_initial_step, " adaptions, ", grid_info%i_cells, " cells"
                    !$omp end master
                end if

                grid_info%i_cells = grid%get_cells(MPI_SUM, .true.)
				if (swe%init%i_refinements_issued .le. grid_info%i_cells / 100_GRID_DI) then
					exit
				endif

				call swe%adaption%traverse(grid)

				i_initial_step = i_initial_step + 1
			end do

            grid_info = grid%get_info(MPI_SUM, .true.)

            if (rank_MPI == 0) then
                !$omp master
                _log_write(0, *) "SWE: done."
                _log_write(0, *) ""

                call grid_info%print()
                !$omp end master
			end if

			!output initial grid
			if (cfg%r_output_time_step >= 0.0_GRID_SR) then
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
            call swe%init%reduce_stats(MPI_SUM, .true.)
            call swe%adaption%reduce_stats(MPI_SUM, .true.)
            call grid%reduce_stats(MPI_SUM, .true.)

            if (rank_MPI == 0) then
                _log_write(0, *) "SWE: running time steps.."
                _log_write(0, *) ""
			end if
            !$omp end master

            i_time_step = 0

#           if defined(_ASAGI)
                ! during the earthquake, do small time steps that include a displacement

                do
                    if ((cfg%r_max_time >= 0.0 .and. grid%r_time > cfg%r_max_time) .or. (cfg%i_max_time_steps >= 0 .and. i_time_step >= cfg%i_max_time_steps)) then
                        exit
                    end if

                    if (grid%r_time > asagi_grid_max(cfg%afh_displacement,3)) then
                        exit
                    end if

                    !do an euler time step
                    call swe%adaption%traverse(grid)

                    call swe%euler%traverse(grid)
                    i_time_step = i_time_step + 1

                    !displace time-dependent bathymetry
                    call swe%displace%traverse(grid)

                    if (rank_MPI == 0) then
                        grid_info%i_cells = grid%get_cells(MPI_SUM, .false.)

                        !$omp master
                        _log_write(1, '(A, I0, A, ES14.7, A, ES14.7, A, I0)') " SWE: EQ time step: ", i_time_step, ", sim. time:", grid%r_time, " s, dt:", grid%r_dt, " s, cells: ", grid_info%i_cells
                        !$omp end master
                    end if

                    !output grid
                    if (cfg%r_output_time_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output) then
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
				if ((cfg%r_max_time >= 0.0 .and. grid%r_time > cfg%r_max_time) .or. (cfg%i_max_time_steps >= 0 .and. i_time_step >= cfg%i_max_time_steps)) then
					exit
				end if
! TODO : correct adaption operators for SWE_SIMD
#ifndef _SWE_SIMD
				call swe%adaption%traverse(grid)
#endif

				!do a time step
				call swe%euler%traverse(grid)
				i_time_step = i_time_step + 1

                if (rank_MPI == 0) then
                    grid_info%i_cells = grid%get_cells(MPI_SUM, .false.)

                    !$omp master
                    _log_write(1, '(A, I0, A, ES14.7, A, ES14.7, A, I0)') " SWE: time step: ", i_time_step, ", sim. time:", grid%r_time, " s, dt:", grid%r_dt, " s, cells: ", grid_info%i_cells
                    !$omp end master
                end if

				!output grid
				if (cfg%r_output_time_step >= 0.0_GRID_SR .and. grid%r_time >= r_time_next_output) then
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
                _log_write(0, '(" SWE: done.")')
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

			!$omp master
                !Initially, just start the timer and don't print anything
                if (t_phase < huge(1.0d0)) then
                    t_phase = t_phase + get_wtime()

                    call swe%init%reduce_stats(MPI_SUM, .true.)
                    call swe%displace%reduce_stats(MPI_SUM, .true.)
                    call swe%euler%reduce_stats(MPI_SUM, .true.)
                    call swe%adaption%reduce_stats(MPI_SUM, .true.)
                    call grid%reduce_stats(MPI_SUM, .true.)

                    if (rank_MPI == 0) then
                        _log_write(0, *) ""
                        _log_write(0, *) "Phase statistics:"
                        _log_write(0, *) ""
                        _log_write(0, '(A, T34, A)') " Init: ", trim(swe%init%stats%to_string())
                        _log_write(0, '(A, T34, A)') " Displace: ", trim(swe%displace%stats%to_string())
                        _log_write(0, '(A, T34, A)') " Time steps: ", trim(swe%euler%stats%to_string())
                        _log_write(0, '(A, T34, A)') " Adaptions: ", trim(swe%adaption%stats%to_string())
                        _log_write(0, '(A, T34, A)') " Grid: ", trim(grid%stats%to_string())
                        _log_write(0, '(A, T34, F12.4, A)') " Element throughput: ", 1.0d-6 * dble(grid%stats%i_traversed_cells) / t_phase, " M/s"
                        _log_write(0, '(A, T34, F12.4, A)') " Memory throughput: ", dble(grid%stats%i_traversed_memory) / ((1024 * 1024 * 1024) * t_phase), " GB/s"
                        _log_write(0, '(A, T34, F12.4, A)') " Cell update throughput: ", 1.0d-6 * dble(swe%euler%stats%i_traversed_cells) / t_phase, " M/s"
                        _log_write(0, '(A, T34, F12.4, A)') " Flux solver throughput: ", 1.0d-6 * dble(swe%euler%stats%i_traversed_edges) / t_phase, " M/s"
                        _log_write(0, '(A, T34, F12.4, A)') " Asagi time:", grid%stats%r_asagi_time, " s"
                        _log_write(0, '(A, T34, F12.4, A)') " Phase time:", t_phase, " s"
                        _log_write(0, *) ""
                    end if
                end if

                call swe%init%clear_stats()
                call swe%displace%clear_stats()
                call swe%euler%clear_stats()
                call swe%adaption%clear_stats()
                call grid%clear_stats()

                t_phase = -get_wtime()
            !$omp end master
        end subroutine
	END MODULE SWE
#endif
