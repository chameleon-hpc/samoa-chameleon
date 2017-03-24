! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"
#include "Tools_grid.f90"

module Tools_loadbalancing
   use SFC_data_types
   use Grid
   use Tools_mpi

   public

   contains
	
#   if defined (_MPI)

	     subroutine compute_rank_load_and_throughput(grid, rank_load, rank_throughput)
        		type(t_grid), intent(inout)   :: grid
       			 integer (kind = GRID_DI), intent(inout) :: rank_load
        		double precision, intent(out) :: rank_throughput
		
			integer :: min_load
			double precision :: rank_computation_time
		
			! rank load is simply the sum of loads of all sections in the rank
			rank_load = sum(grid%sections%elements_alloc(:)%load)
		
			! if not using heterogeneous LB, set throughput as 1 and return
			if (cfg%l_lb_hh .eqv. .false.) then 
			    rank_throughput = 1.0d0
			    return 
			endif

			! if cfg%l_lb_hh_auto is not set, use Host/MIC ratio chosen by user (option -lbhhratio HOST MIC, default is 1 1)
			if (cfg%l_lb_hh_auto .eqv. .false.) then
#		               if defined (__MIC__)
				rank_throughput = cfg%r_lb_hh_ratio(2)
#		               else
				rank_throughput = cfg%r_lb_hh_ratio(1)
#		               endif
			else
			    
			    ! if at least one rank has zero load, distribute the load evenly so that it is possible to compute 
			    ! the throughput for every rank at the next iteration
			    ! do the same if still in grid generation phase
			    min_load = rank_load
			    call reduce(min_load, MPI_MIN)

			    if (min_load == 0 .or. grid%l_grid_generation) then
				rank_throughput = 1
				grid%r_computation_time_since_last_LB = 0

			    else
				! now we need to actually compute the throughput
				
				! compute throughput based on time for the last steps (from the average thread)
				rank_computation_time = grid%r_computation_time_since_last_LB / size(grid%threads%elements(:))
				
				grid%r_computation_time_since_last_LB = 0
				
				rank_throughput = rank_load / max(rank_computation_time, 1.0e-6) ! avoid division by zero
				
			    end if
			endif
    	    end subroutine

	    subroutine compute_imbalance(grid, rank_load, rank_throughput, imbalance)
            type(t_grid), intent(in)   :: grid
            integer (kind = GRID_DI), intent(in) :: rank_load
            double precision, intent(in) :: rank_throughput
            double precision, intent(out) :: imbalance
            double precision :: total_throughput
            integer :: total_load

            ! compute total_load and total_throughput (sum from all ranks)
            total_load = rank_load
            total_throughput = rank_throughput
            call reduce(total_load, MPI_SUM)
            call reduce(total_throughput, MPI_SUM)
            
            
            if (rank_load > 0) then
                ! the computation below is the same as: max_imbalance = (rank_load/rank_throughput)/(total_load/total_throughput)
                ! but with only one division.
                imbalance = (rank_load*total_throughput)/(total_load*rank_throughput)
                
                if (imbalance < 1) then ! avoid negative values
                    if (imbalance == 0) then ! avoid division by zero
                        imbalance = 1.0 
                    else
                        imbalance = 1.0 / imbalance
                    end if
                end if
                imbalance = imbalance - 1 ! 0 = perfectly balanced
            else 
                imbalance = 1 ! if this rank is empty, it is not balanced
            end if

            ! consider only the greatest value among all ranks
            call reduce(imbalance, MPI_MAX)
	    end subroutine

       !> uses a distributed algorithm to compute the new partitiion
        !> this approach scales well, but requires local methods and cannot achieve optimal load balance
        subroutine compute_partition_distributed(grid, i_rank_out_local, i_section_index_out_local, i_rank_in_local, l_early_exit)
		    type(t_grid), intent(inout)		                :: grid
		    integer, pointer, intent(inout)                 :: i_rank_out_local(:), i_section_index_out_local(:), i_rank_in_local(:)
		    logical, intent(out)                            :: l_early_exit

            integer, pointer, save                          :: i_rank_out(:) => null(), i_section_index_out(:) => null(), i_rank_in(:) => null()
            integer, allocatable, save                      :: i_sections_out(:), i_sections_in(:), i_partial_sections_in(:), i_partial_sections_out(:), i_delta_out(:), requests_out(:, :), requests_in(:)
            integer (kind = GRID_SI)						:: i_first_local_section, i_last_local_section, i_rank, i_section, i_comm
            integer									        :: i_first_rank_out, i_last_rank_out, i_first_rank_in, i_last_rank_in
            type(t_grid_section), pointer					:: section
            integer						                    :: i_error, i_sections, requests(2)
            integer	(BYTE)  		                        :: i_color
            integer (kind = GRID_DI)					    :: load, partial_load, total_load, rank_imbalance
            character (len=100)                             :: msg

	        l_early_exit = .false.
        	call grid%get_local_sections(i_first_local_section, i_last_local_section)

            !$omp single
			!switch to integer arithmetics from now on, we need exact arithmetics

			load = grid%load
            partial_load = load
			call prefix_sum(partial_load, MPI_SUM)
			total_load = load
			call reduce(total_load, MPI_SUM)

			assert_gt(total_load, 0)
			i_sections = grid%sections%get_size()

            if (i_sections > 0) then
                rank_imbalance = \
                    ((2_GRID_DI * partial_load - grid%sections%elements_alloc(i_sections)%load) * size_MPI) / (2_GRID_DI * total_load) - \
                    ((2_GRID_DI * (partial_load - load) + grid%sections%elements_alloc(1)%load) * size_MPI) / (2_GRID_DI * total_load)
            else
                rank_imbalance = 0
            end if

			call reduce(rank_imbalance, MPI_MAX)

			!allocate arrays on first call
			if (.not. associated(i_rank_out)) then
                allocate(i_rank_out(0), stat=i_error); assert_eq(i_error, 0)
                allocate(i_section_index_out(0), stat=i_error); assert_eq(i_error, 0)
            end if

            if (.not. allocated(i_sections_out)) then
                allocate(i_sections_out(0), stat=i_error); assert_eq(i_error, 0)
                allocate(i_partial_sections_out(0), stat=i_error); assert_eq(i_error, 0)
                allocate(i_delta_out(0), stat=i_error); assert_eq(i_error, 0)
                allocate(requests_out(0, 0), stat=i_error); assert_eq(i_error, 0)
            end if

            if (.not. allocated(i_sections_in)) then
                allocate(i_sections_in(0), stat=i_error); assert_eq(i_error, 0)
                allocate(i_partial_sections_in(0), stat=i_error); assert_eq(i_error, 0)
                allocate(requests_in(0), stat=i_error, errmsg=msg); if (i_error .ne. 0) then; print *, msg; end if; assert_eq(i_error, 0)
            end if

	        !$omp end single copyprivate(rank_imbalance, i_sections, load, partial_load, total_load)

	        !exit early if the imbalance cannot be improved
	        if (rank_imbalance .le. 0) then
                _log_write(2, '(4X, "load balancing: imbalance is improvable? no.")')
                l_early_exit = .true.
                return
	        end if

	        _log_write(2, '(4X, "load balancing: imbalance is improvable? yes.")')

			!$omp single
            i_first_rank_out = ((partial_load - load) * size_MPI) / total_load	!round down
            i_last_rank_out = (partial_load * size_MPI - 1) / total_load		!round up and subtract 1

            !special case: integers are rounded up below 0, so if partial_load = 0, the result is 0, but should be -1
            if (partial_load == 0) then
                i_last_rank_out = -1
            end if

			!allocate space for variables that are stored per output section
			if (size(i_rank_out) .ne. i_sections) then
                deallocate(i_rank_out, stat = i_error); assert_eq(i_error, 0)
                deallocate(i_section_index_out, stat = i_error); assert_eq(i_error, 0)

                allocate(i_rank_out(i_sections), stat=i_error); assert_eq(i_error, 0)
                allocate(i_section_index_out(i_sections), stat=i_error); assert_eq(i_error, 0)
            end if

			!allocate space for variables that are stored per output rank
            if (lbound(i_sections_out, 1) .ne. i_first_rank_out .or. ubound(i_sections_out, 1) .ne. i_last_rank_out) then
                deallocate(i_sections_out, stat=i_error); assert_eq(i_error, 0)
                deallocate(i_partial_sections_out, stat=i_error); assert_eq(i_error, 0)
                deallocate(i_delta_out, stat=i_error); assert_eq(i_error, 0)
                deallocate(requests_out, stat=i_error); assert_eq(i_error, 0)

                allocate(i_sections_out(i_first_rank_out : i_last_rank_out), stat=i_error); assert_eq(i_error, 0)
                allocate(i_partial_sections_out(i_first_rank_out : i_last_rank_out), stat=i_error); assert_eq(i_error, 0)
                allocate(i_delta_out(i_first_rank_out : i_last_rank_out), stat=i_error); assert_eq(i_error, 0)
                allocate(requests_out(i_first_rank_out : i_last_rank_out, 2), stat=i_error); assert_eq(i_error, 0)
            end if

			i_sections_out = 0
			requests_out = MPI_REQUEST_NULL
			requests = MPI_REQUEST_NULL
			!$omp end single copyprivate(i_first_rank_out, i_last_rank_out)

			_log_write(3, '(4X, "match ranks: out ranks ", I0, " to ", I0)') i_first_rank_out, i_last_rank_out

	        !compute the number of sections, that are sent to each rank
	        do i_section = i_first_local_section, i_last_local_section
	            section => grid%sections%elements_alloc(i_section)
	            assert_eq(section%index, i_section)

	            !assign this section to its respective ideal rank
	            i_rank = (((2_GRID_DI * (partial_load - load + section%partial_load) - section%load)) * size_MPI) / (2_GRID_DI * total_load)
	        	assert_ge(i_rank, i_first_rank_out)
	        	assert_le(i_rank, i_last_rank_out)

	        	i_rank_out(i_section) = i_rank

			 	!$omp atomic
                i_sections_out(i_rank) = i_sections_out(i_rank) + 1
	        end do

			!$omp barrier

	        !$omp single
	        assert_eq(sum(i_sections_out), i_sections)
			!receive 2 messages, containing the rank number of the first and the last source rank (they must be unique)

			call mpi_irecv(i_first_rank_in, 1, MPI_INTEGER, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, requests(1), i_error); assert_eq(i_error, 0)
			call mpi_irecv(i_last_rank_in, 1, MPI_INTEGER, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, requests(2), i_error); assert_eq(i_error, 0)

			!My rank must be the first input rank for all my output ranks except the first one
			!and the last input rank for all my output ranks except the last one

            do i_rank = i_first_rank_out + 1, i_last_rank_out
                call mpi_isend(rank_MPI, 1, MPI_INTEGER, i_rank, 1, MPI_COMM_WORLD, requests_out(i_rank, 1), i_error); assert_eq(i_error, 0)
                call mpi_isend(rank_MPI, 1, MPI_INTEGER, i_rank - 1, 2, MPI_COMM_WORLD, requests_out(i_rank - 1, 2), i_error); assert_eq(i_error, 0)

                _log_write(3, '("send: from: ", I0, " to: ", I0, " I am your first input rank")') rank_MPI, i_rank
                _log_write(3, '("send: from: ", I0, " to: ", I0, " I am your last input rank")') rank_MPI, i_rank - 1
            end do

			!check if i am first or last input rank for first and last output rank (meaning an exact match)
			if (load > 0 .and. (partial_load - load) * size_MPI == total_load * i_first_rank_out) then
                call mpi_isend(rank_MPI, 1, MPI_INTEGER, i_first_rank_out, 1, MPI_COMM_WORLD, requests_out(i_first_rank_out, 1), i_error); assert_eq(i_error, 0)
               _log_write(3, '("send: from: ", I0, " to: ", I0, " I am your first input rank")') rank_MPI, i_first_rank_out
			end if

			if (load > 0 .and. partial_load * size_MPI == total_load * (i_last_rank_out + 1)) then
                call mpi_isend(rank_MPI, 1, MPI_INTEGER, i_last_rank_out, 2, MPI_COMM_WORLD, requests_out(i_last_rank_out, 2), i_error); assert_eq(i_error, 0)
                _log_write(3, '("send: from: ", I0, " to: ", I0, " I am your last input rank")') rank_MPI, i_last_rank_out
			end if

			call mpi_waitall(i_last_rank_out - i_first_rank_out + 1, requests_out(i_first_rank_out:i_last_rank_out, 1), MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
			call mpi_waitall(i_last_rank_out - i_first_rank_out + 1, requests_out(i_first_rank_out:i_last_rank_out, 2), MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
			call mpi_waitall(2, requests, MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)

			_log_write(3, '(4X, "count sections: out ranks ", I0, " to ", I0, " in ranks: ", I0, " to ", I0)') i_first_rank_out, i_last_rank_out, i_first_rank_in, i_last_rank_in

			!now that I know which ranks I get data from, allocate an array that stores the number of sections I get from each source rank
            if (lbound(i_sections_in, 1) .ne. i_first_rank_in .or. ubound(i_sections_in, 1) .ne. i_last_rank_in) then
                deallocate(i_sections_in, stat=i_error); assert_eq(i_error, 0)
                deallocate(i_partial_sections_in, stat=i_error); assert_eq(i_error, 0)
                deallocate(requests_in, stat=i_error, errmsg=msg); if (i_error .ne. 0) then; print *, msg; end if; assert_eq(i_error, 0)

                allocate(i_sections_in(i_first_rank_in : i_last_rank_in), stat=i_error); assert_eq(i_error, 0)
                allocate(i_partial_sections_in(i_first_rank_in : i_last_rank_in), stat=i_error); assert_eq(i_error, 0)
                allocate(requests_in(i_first_rank_in : i_last_rank_in), stat=i_error, errmsg=msg); if (i_error .ne. 0) then; print *, msg; end if; assert_eq(i_error, 0)
            end if

			i_sections_in = 0
			requests_in = MPI_REQUEST_NULL
			!$omp end single copyprivate(i_first_rank_in, i_last_rank_in)

	        !$omp single
 		    !Communicate the number of sections sent to the destination ranks / received from the source ranks

			_log_write(3, '(4X, "send number of outgoing sections and receive number of incoming sections per rank")')

	        assert_veq(requests_in, MPI_REQUEST_NULL)
	        assert_veq(requests_out, MPI_REQUEST_NULL)

			do i_rank = i_first_rank_in, i_last_rank_in
				call mpi_irecv(i_sections_in(i_rank), 1, MPI_INTEGER, i_rank, 0, MPI_COMM_WORLD, requests_in(i_rank), i_error); assert_eq(i_error, 0)
        	end do

			do i_rank = i_first_rank_out, i_last_rank_out
                call mpi_isend(i_sections_out(i_rank), 1, MPI_INTEGER, i_rank, 0, MPI_COMM_WORLD, requests_out(i_rank, 1), i_error); assert_eq(i_error, 0)
        	end do

			call mpi_waitall(i_last_rank_out - i_first_rank_out + 1, requests_out(i_first_rank_out:i_last_rank_out, 1), MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
 	    	call mpi_waitall(i_last_rank_in - i_first_rank_in + 1, requests_in, MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)

	        assert_veq(requests_out, MPI_REQUEST_NULL)
			assert_veq(requests_in, MPI_REQUEST_NULL)

			_log_write(3, '(4X, "communicate new index of first incoming section per rank")')

			!compute prefix sum over number of input sections and output sections to find out which rank gets which section indices assigned
	        call prefix_sum(i_partial_sections_in, i_sections_in)
	        call prefix_sum(i_partial_sections_out, i_sections_out)

			!communicate this info again
			do i_rank = i_first_rank_out, i_last_rank_out
                call mpi_irecv(i_delta_out(i_rank), 1, MPI_INTEGER, i_rank, 0, MPI_COMM_WORLD, requests_out(i_rank, 1), i_error); assert_eq(i_error, 0)
        	end do

			do i_rank = i_first_rank_in, i_last_rank_in
                call mpi_isend(i_partial_sections_in(i_rank), 1, MPI_INTEGER, i_rank, 0, MPI_COMM_WORLD, requests_in(i_rank), i_error); assert_eq(i_error, 0)
        	end do

            if (.not. associated(i_rank_in) .or. size(i_rank_in) .ne. sum(i_sections_in)) then
                if (associated(i_rank_in)) then
                    deallocate(i_rank_in, stat=i_error); assert_eq(i_error, 0)
                end if

                allocate(i_rank_in(sum(i_sections_in)), stat=i_error); assert_eq(i_error, 0)
            end if

            i_section = 0
            do i_rank = i_first_rank_in, i_last_rank_in
                if (i_sections_in(i_rank) > 0) then
                    i_rank_in(i_section + 1 : i_section + i_sections_in(i_rank)) = i_rank
                    i_section = i_section + i_sections_in(i_rank)
                end if
            end do

			call mpi_waitall(i_last_rank_out - i_first_rank_out + 1, requests_out(i_first_rank_out:i_last_rank_out, 1), MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
	    	call mpi_waitall(i_last_rank_in - i_first_rank_in + 1, requests_in, MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)

			_log_write(3, '(4X, "compute new indices of all outgoing sections")')

        	i_delta_out = i_delta_out - i_partial_sections_out
	        !$omp end single

        	!compute the new section indices
	        do i_section = i_first_local_section, i_last_local_section
	            i_section_index_out(i_section) = i_delta_out(i_rank_out(i_section)) + i_section
	        end do

            !pass private copies of the array pointers back to caller function
            i_rank_out_local => i_rank_out
            i_section_index_out_local => i_section_index_out
	        i_rank_in_local => i_rank_in
        end subroutine

	!> Computes chains-on-chains assignment of n tasks with given load l to processes using the midpoint cutoff approximation.
        !> @param l Input load array
	!> @param s Output assignment of tasks to ranks 
        subroutine midpoint_cutoff(l, s)
            integer(kind = GRID_DI), intent(in)     :: l(:)
            integer(kind = GRID_SI), intent(inout)  :: s(:)

            integer                                 :: n, i, j
            integer (kind = GRID_DI)                :: L_j, L_n

            n = size(l)
            L_n = sum(l)
            L_j = 0_GRID_DI

            do j = 1, n
                i = int((size_MPI * (2_GRID_DI * L_j + l(j))) / (2_GRID_DI * L_n))
                s(j) = i

                L_j = L_j + l(j)
            end do
        end subroutine

	!> Computes chains-on-chains assignment of n tasks with given load l to processes using exact optimal binary search.
        !> @param l Input load array
	!> @param s Output assignment of tasks to ranks 
        subroutine iterative_binary(l, s)
            integer(kind = GRID_DI), intent(in)     :: l(:)
            integer(kind = GRID_SI), intent(inout)  :: s(:)

            integer(kind = GRID_SI), allocatable    :: s_test(:)
            integer                                 :: n, i_error, i, j
            integer (kind = GRID_DI)                :: test, test_max, current_min, current_max, load_i

            n = size(l)

            allocate(s_test(n), stat=i_error); assert_eq(i_error, 0)

            current_max = 0_GRID_DI
            i = 0
            load_i = 0_GRID_DI

            do j = 1, n
                if (s(j) > i) then
                    current_max = max(current_max, load_i)
                    i = s(j)
                    load_i = 0_GRID_DI
                end if

                load_i = load_i + l(j)
            end do

            current_max = max(current_max, load_i) + 1_GRID_DI
            current_min = max(sum(l) / size_MPI, maxval(l))

            test = (current_min + current_max) / 2

            do
                if (current_max <= current_min) then
                    exit
                end if

                test_max = 0_GRID_DI
                i = 0
                load_i = 0_GRID_DI

                do j = 1, n
                    if (load_i + l(j) >= test .and. i < size_MPI - 1) then
                        test_max = max(test_max, load_i)
                        i = i + 1
                        load_i = 0_GRID_DI
                    end if

                    load_i = load_i + l(j)
                    s_test(j) = i
                end do

                test_max = max(test_max, load_i)

                if (test_max < current_max) then
                    current_max = test_max
                    test = (current_min + current_max) / 2

                    s(:) = s_test(:)
                else
                    current_min = test
                    test = current_max
                end if
            end do

            deallocate(s_test, stat=i_error); assert_eq(i_error, 0)
        end subroutine

        !> uses a serial algorithm to compute the new partition
        !> this approach scales badly, but can apply global methods and can achieve optimal load balance
        subroutine compute_partition_serial(grid, i_rank_out_local, i_section_index_out_local, i_rank_in_local, l_early_exit)
		    type(t_grid), intent(inout)		                :: grid
		    integer, pointer, intent(inout)                 :: i_rank_out_local(:), i_section_index_out_local(:), i_rank_in_local(:)
		    logical, intent(out)                            :: l_early_exit

            integer, pointer, save                          :: i_rank_out(:) => null(), i_section_index_out(:) => null(), i_rank_in(:) => null()
            integer, allocatable, save                      :: all_sections(:), displacements(:), all_ranks(:), all_section_indices_out(:)
            integer, allocatable, save                      :: requests_in(:), requests_out(:)
            integer (kind = GRID_DI), allocatable, save     :: all_load(:), local_load(:)
            integer                                         :: total_sections, i_sections_out, i_sections_in, i_error, i, j, k

			l_early_exit = .false.

            !$omp single

            !gather section indices

            if (rank_MPI == 0) then
                allocate(all_sections(0 : size_MPI - 1), stat=i_error); assert_eq(i_error, 0)
                allocate(displacements(0 : size_MPI - 1), stat=i_error); assert_eq(i_error, 0)
            else
                allocate(all_sections(0), stat=i_error); assert_eq(i_error, 0)
                allocate(displacements(0), stat=i_error); assert_eq(i_error, 0)
            end if

            i_sections_out = grid%sections%get_size()

            call mpi_gather(i_sections_out, 1, MPI_INTEGER, all_sections, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)

            !gather load

            if (rank_MPI == 0) then
                call prefix_sum(displacements, all_sections)
                total_sections = displacements(size_MPI - 1)
                displacements(:) = displacements(:) - all_sections(:)

                allocate(all_load(total_sections), stat=i_error); assert_eq(i_error, 0)
                allocate(all_ranks(total_sections), stat=i_error); assert_eq(i_error, 0)
                allocate(all_section_indices_out(total_sections), stat=i_error); assert_eq(i_error, 0)
            else
                allocate(all_load(0), stat=i_error); assert_eq(i_error, 0)
                allocate(all_ranks(0), stat=i_error); assert_eq(i_error, 0)
                allocate(all_section_indices_out(0), stat=i_error); assert_eq(i_error, 0)
            end if

            allocate(local_load(i_sections_out), stat=i_error); assert_eq(i_error, 0)

            local_load(:) = grid%sections%elements_alloc(:)%load
            call mpi_gatherv(local_load, i_sections_out, MPI_INTEGER8, all_load, all_sections, displacements, MPI_INTEGER8, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)

            if (rank_MPI == 0) then
                j = 0
                do i = 0, size_MPI - 1
                    do k = 1, all_sections(i)
                        all_ranks(j + k) = i
                    end do

                    j = j + all_sections(i)
                end do

                call iterative_binary(all_load, all_ranks)

                all_sections(:) = 0

                do j = 1, total_sections
                    i = all_ranks(j)
                    all_sections(i) = all_sections(i) + 1
                    all_section_indices_out(j) = all_sections(i)
                end do
            end if

            !scatter new section count

            call mpi_scatter(all_sections, 1, MPI_INTEGER, i_sections_in, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)

            !scatter new ranks

            !allocate space for variables that are stored per input section
			if (.not. associated(i_rank_in) .or. size(i_rank_in) .ne. i_sections_in) then
                if (associated(i_rank_in)) then
                    deallocate(i_rank_in, stat = i_error); assert_eq(i_error, 0)
                    deallocate(requests_in, stat=i_error); assert_eq(i_error, 0)
                end if

                allocate(i_rank_in(i_sections_in), stat=i_error); assert_eq(i_error, 0)
                allocate(requests_in(i_sections_in), stat=i_error); assert_eq(i_error, 0)
            end if

			!allocate space for variables that are stored per output section
			if (.not. associated(i_rank_out) .or. size(i_rank_out) .ne. i_sections_out) then
                if (associated(i_rank_out)) then
                    deallocate(i_rank_out, stat = i_error); assert_eq(i_error, 0)
                    deallocate(i_section_index_out, stat = i_error); assert_eq(i_error, 0)
                    deallocate(requests_out, stat = i_error); assert_eq(i_error, 0)
                end if

                allocate(i_rank_out(i_sections_out), stat=i_error); assert_eq(i_error, 0)
                allocate(i_section_index_out(i_sections_out), stat=i_error); assert_eq(i_error, 0)
                allocate(requests_out(i_sections_out), stat=i_error); assert_eq(i_error, 0)
            end if

            if (rank_MPI == 0) then
                all_sections(0 : size_MPI - 2) = displacements(1 : size_MPI - 1) - displacements(0 : size_MPI - 2)
                all_sections(size_MPI - 1) = total_sections - displacements(size_MPI - 1)
            end if

            call mpi_scatterv(all_ranks, all_sections, displacements, MPI_INTEGER, i_rank_out, i_sections_out, MPI_INTEGER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
            call mpi_scatterv(all_section_indices_out, all_sections, displacements, MPI_INTEGER, i_section_index_out, i_sections_out, MPI_INTEGER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)

            deallocate(local_load, stat=i_error); assert_eq(i_error, 0)
            deallocate(all_load, stat=i_error); assert_eq(i_error, 0)
            deallocate(all_sections, stat=i_error); assert_eq(i_error, 0)
            deallocate(displacements, stat=i_error); assert_eq(i_error, 0)
            deallocate(all_ranks, stat=i_error); assert_eq(i_error, 0)
            deallocate(all_section_indices_out, stat=i_error); assert_eq(i_error, 0)

            do j = 1, i_sections_in
                call mpi_irecv(i_rank_in(j), 1, MPI_INTEGER, MPI_ANY_SOURCE, j, MPI_COMM_WORLD, requests_in(j), i_error); assert_eq(i_error, 0)
            end do

            do j = 1, i_sections_out
                call mpi_isend(rank_MPI, 1, MPI_INTEGER, i_rank_out(j), i_section_index_out(j), MPI_COMM_WORLD, requests_out(j), i_error); assert_eq(i_error, 0)
            end do

            call mpi_waitall(i_sections_in, requests_in, MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
            call mpi_waitall(i_sections_out, requests_out, MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
            !$omp end single

            !pass private copies of the array pointers back to caller function
            i_rank_out_local => i_rank_out
            i_section_index_out_local => i_section_index_out
            i_rank_in_local => i_rank_in
        end subroutine
        
        !> for heterogeneous hardware: uses a serial algorithm to compute the new partition based on each rank's performance
        !> this serial approach scales badly, but can apply global methods and can achieve optimal load balance
        subroutine compute_partition_serial_heterogeneous(grid, i_rank_out_local, i_section_index_out_local, i_rank_in_local, l_early_exit, my_throughput, my_load)
            type(t_grid), intent(inout)                     :: grid
            integer, pointer, intent(inout)                 :: i_rank_out_local(:), i_section_index_out_local(:), i_rank_in_local(:)
            logical, intent(out)                            :: l_early_exit
            double precision, intent(inout)                 :: my_throughput 
            integer (kind = GRID_DI), intent(inout)         :: my_load

            integer, pointer, save                          :: i_rank_out(:) => null(), i_section_index_out(:) => null(), i_rank_in(:) => null()
            integer, allocatable, save                      :: all_sections(:), displacements(:), all_ranks(:), all_section_indices_out(:)
            integer, allocatable, save                      :: requests_in(:), requests_out(:)
            integer (kind = GRID_DI), allocatable, save     :: all_load(:), local_load(:)
            integer                                         :: total_sections, i_sections_out, i_sections_in, i_error, i, j, k
            
            double precision :: total_throughput ! sum of all ranks' throughputs
            double precision, allocatable :: rank_throughput(:) ! throughputs of each rank
            double precision, allocatable :: section_position(:) ! proportional position of section
            integer (kind = GRID_DI), allocatable :: rank_load(:) ! total load of each rank
            integer (kind = GRID_DI), allocatable :: prefix_sum_load(:)
            integer (kind = GRID_DI) :: total_load

            l_early_exit = .false.

            !$omp single
            
                if (rank_MPI == 0) then
                    allocate(all_sections(0 : size_MPI - 1), stat=i_error); assert_eq(i_error, 0)
                    allocate(displacements(0 : size_MPI - 1), stat=i_error); assert_eq(i_error, 0)
                    allocate(rank_throughput(0 : size_MPI - 1), stat=i_error); assert_eq(i_error, 0)
                    allocate(rank_load(0 : size_MPI - 1), stat=i_error); assert_eq(i_error, 0)
                else
                    allocate(all_sections(0), stat=i_error); assert_eq(i_error, 0)
                    allocate(displacements(0), stat=i_error); assert_eq(i_error, 0)
                    allocate(rank_throughput(0), stat=i_error); assert_eq(i_error, 0)
                    allocate(rank_load(0), stat=i_error); assert_eq(i_error, 0)
                end if

                !gather section indices
                i_sections_out = grid%sections%get_size()
                call mpi_gather(i_sections_out, 1, MPI_INTEGER, all_sections, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
                
                ! gather load and throughput from all ranks
                call mpi_gather(my_load, 1, MPI_INTEGER8, rank_load, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
                call mpi_gather(my_throughput, 1, MPI_DOUBLE, rank_throughput, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
                
                ! apply prefix sum to rank_throughput and compute total_throughput
                if (rank_MPI == 0) then
                    call prefix_sum(rank_throughput, rank_throughput)
                    total_throughput = rank_throughput(size_MPI - 1)
                    total_load = sum(rank_load)
                endif
                
                call MPI_bcast(total_load, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
                         
                !gather load
                if (rank_MPI == 0) then
                    call prefix_sum(displacements, all_sections)
                    total_sections = displacements(size_MPI - 1)
                    displacements(:) = displacements(:) - all_sections(:)

                    allocate(all_load(total_sections), stat=i_error); assert_eq(i_error, 0)
                    allocate(prefix_sum_load(total_sections), stat=i_error); assert_eq(i_error, 0)
                    allocate(section_position(total_sections), stat=i_error); assert_eq(i_error, 0)
                    allocate(all_ranks(total_sections), stat=i_error); assert_eq(i_error, 0)
                    allocate(all_section_indices_out(total_sections), stat=i_error); assert_eq(i_error, 0)
                else
                    allocate(all_load(0), stat=i_error); assert_eq(i_error, 0)
                    allocate(prefix_sum_load(0), stat=i_error); assert_eq(i_error, 0)
                    allocate(section_position(0), stat=i_error); assert_eq(i_error, 0)
                    allocate(all_ranks(0), stat=i_error); assert_eq(i_error, 0)
                    allocate(all_section_indices_out(0), stat=i_error); assert_eq(i_error, 0)
                end if

                allocate(local_load(i_sections_out), stat=i_error); assert_eq(i_error, 0)

                local_load(:) = grid%sections%elements_alloc(:)%load
                call mpi_gatherv(local_load, i_sections_out, MPI_INTEGER8, all_load, all_sections, displacements, MPI_INTEGER8, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
                
                if (rank_MPI == 0) then
                
                    _log_write(0, '(4X, "Before LB, load: ", F6.3, " ", F6.3, " ", F6.3)'), dble(rank_load(0))/total_load, dble(rank_load(1))/total_load, dble(rank_load(2))/total_load
                    
                    ! assign each section to a rank according to the relative sections' loads and the ranks' throuhputs
                    ! first compute a relative position for each section (relative to total throughput)
                    call prefix_sum(prefix_sum_load, all_load)
                    do i = 1, total_sections
                        section_position(i) = (prefix_sum_load(i) - 0.5 * all_load(i)) * total_throughput / total_load
                    end do
                    
                    ! now assign each section to a rank
                    j = 0 ! rank
                    do i = 1, total_sections ! i = section id
                        do while (section_position(i) > rank_throughput(j)) !this is actually the prefix sum of throughputs
                            j = j + 1
                        end do
                        
                        all_ranks(i) = j
                        
                        assert_le(j, MPI_rank - 1)
                    end do

                    all_sections(:) = 0
                    rank_load = 0

                    do j = 1, total_sections
                        i = all_ranks(j)
                        all_sections(i) = all_sections(i) + 1
                        all_section_indices_out(j) = all_sections(i)
                        
                        rank_load(i) = rank_load(i) + all_load(j)
                    end do
                    _log_write(0, '(4X, "After LB, load: ", F6.3, " ", F6.3, " ", F6.3)'), dble(rank_load(0))/total_load, dble(rank_load(1))/total_load, dble(rank_load(2))/total_load
                end if
                
                !scatter new section count

                call mpi_scatter(all_sections, 1, MPI_INTEGER, i_sections_in, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
                
                !scatter new ranks

                !allocate space for variables that are stored per input section
                if (.not. associated(i_rank_in) .or. size(i_rank_in) .ne. i_sections_in) then
                    if (associated(i_rank_in)) then
                        deallocate(i_rank_in, stat = i_error); assert_eq(i_error, 0)
                        deallocate(requests_in, stat=i_error); assert_eq(i_error, 0)
                    end if

                    allocate(i_rank_in(i_sections_in), stat=i_error); assert_eq(i_error, 0)
                    allocate(requests_in(i_sections_in), stat=i_error); assert_eq(i_error, 0)
                end if

                !allocate space for variables that are stored per output section
                if (.not. associated(i_rank_out) .or. size(i_rank_out) .ne. i_sections_out) then
                    if (associated(i_rank_out)) then
                        deallocate(i_rank_out, stat = i_error); assert_eq(i_error, 0)
                        deallocate(i_section_index_out, stat = i_error); assert_eq(i_error, 0)
                        deallocate(requests_out, stat = i_error); assert_eq(i_error, 0)
                    end if

                    allocate(i_rank_out(i_sections_out), stat=i_error); assert_eq(i_error, 0)
                    allocate(i_section_index_out(i_sections_out), stat=i_error); assert_eq(i_error, 0)
                    allocate(requests_out(i_sections_out), stat=i_error); assert_eq(i_error, 0)
                end if

                if (rank_MPI == 0) then

                    all_sections(0 : size_MPI - 2) = displacements(1 : size_MPI - 1) - displacements(0 : size_MPI - 2)
                    all_sections(size_MPI - 1) = total_sections - displacements(size_MPI - 1)
                end if

                call mpi_scatterv(all_ranks, all_sections, displacements, MPI_INTEGER, i_rank_out, i_sections_out, MPI_INTEGER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)
                call mpi_scatterv(all_section_indices_out, all_sections, displacements, MPI_INTEGER, i_section_index_out, i_sections_out, MPI_INTEGER, 0, MPI_COMM_WORLD, i_error); assert_eq(i_error, 0)

                do j = 1, i_sections_in
                    call mpi_irecv(i_rank_in(j), 1, MPI_INTEGER, MPI_ANY_SOURCE, j, MPI_COMM_WORLD, requests_in(j), i_error); assert_eq(i_error, 0)
                end do

                do j = 1, i_sections_out
                    call mpi_isend(rank_MPI, 1, MPI_INTEGER, i_rank_out(j), i_section_index_out(j), MPI_COMM_WORLD, requests_out(j), i_error); assert_eq(i_error, 0)
                end do

                call mpi_waitall(i_sections_in, requests_in, MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
                call mpi_waitall(i_sections_out, requests_out, MPI_STATUSES_IGNORE, i_error); assert_eq(i_error, 0)
                
                deallocate(local_load, stat=i_error); assert_eq(i_error, 0)
                deallocate(all_load, stat=i_error); assert_eq(i_error, 0)
                deallocate(all_sections, stat=i_error); assert_eq(i_error, 0)
                deallocate(displacements, stat=i_error); assert_eq(i_error, 0)
                deallocate(all_ranks, stat=i_error); assert_eq(i_error, 0)
                deallocate(all_section_indices_out, stat=i_error); assert_eq(i_error, 0)
                deallocate(rank_throughput, stat=i_error); assert_eq(i_error, 0)
                deallocate(rank_load, stat=i_error); assert_eq(i_error, 0)
                deallocate(prefix_sum_load, stat=i_error); assert_eq(i_error, 0)
                deallocate(section_position, stat=i_error); assert_eq(i_error, 0)

            !$omp end single copyprivate(l_early_exit)
            
            !pass private copies of the array pointers back to caller function
            i_rank_out_local => i_rank_out
            i_section_index_out_local => i_section_index_out
            i_rank_in_local => i_rank_in
        end subroutine
        
#   endif

end module
