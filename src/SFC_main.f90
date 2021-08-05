! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

PROGRAM samoa
	USE SFC_traversal
#ifdef CHAMELEON
        USE chameleon_lib
#endif
	implicit none
#ifdef CHAMELEON_CALL
       integer :: ierr
#endif

    !init MPI
    call init_mpi()
#ifdef CHAMELEON_CALL
    ierr = chameleon_init()
    ierr = chameleon_determine_base_addresses(c_null_ptr)
    !$omp parallel
    ierr = chameleon_thread_init()
    !$omp end parallel
#endif

    !read config from program arguments and print it out
    call cfg%read_from_program_arguments()

    if (rank_MPI == 0) then
        call cfg%print()
    end if

    !init element transformation data
    call init_transform_data()

    !run scenario selector
    call sfc_generic()

#ifdef CHAMELEON_CALL
    ierr = chameleon_finalize()
#endif
    !finalize MPI
    call finalize_mpi()

	stop
end PROGRAM samoa
