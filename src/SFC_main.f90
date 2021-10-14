! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

PROGRAM samoa
    USE SFC_traversal
#ifdef _TRACE_EXTRAE
    USE EXTRAE_MODULE
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR, C_PTR, C_LOC
#endif
#ifdef CHAMELEON
    USE chameleon_lib
#endif
    implicit none
#ifdef CHAMELEON_CALL
    integer :: ierr
#endif
#if defined(_TRACE_EXTRAE)
    INTEGER*8, PARAMETER, DIMENSION(8)  :: values = (/ 0, 10, 11, 12, 13, 14, 15, 16 /)

    CHARACTER(KIND=C_CHAR,LEN=20)       :: evt_desc = "SamoaEvents" // C_NULL_CHAR
    CHARACTER(KIND=C_CHAR,LEN=30), DIMENSION(8), TARGET :: description_values
    TYPE(C_PTR), DIMENSION(8)           :: description_values_ptrs

    description_values(1) = "End " // C_NULL_CHAR
    description_values_ptrs(1) = C_LOC(description_values(1))

    description_values(2) = 'traverse_pre' // C_NULL_CHAR
    description_values_ptrs(2) = C_LOC(description_values(2))

    description_values(3) = 'traverse_post' // C_NULL_CHAR
    description_values_ptrs(3) = C_LOC(description_values(3))

    description_values(4) = 'traverse_collect' // C_NULL_CHAR
    description_values_ptrs(4) = C_LOC(description_values(4))

    description_values(5) = 'traverse_compute_pre' // C_NULL_CHAR
    description_values_ptrs(5) = C_LOC(description_values(5))

    description_values(6) = 'traverse_copy_boundary' // C_NULL_CHAR
    description_values_ptrs(6) = C_LOC(description_values(6))

    description_values(7) = 'traverse_compute_inner' // C_NULL_CHAR
    description_values_ptrs(7) = C_LOC(description_values(7))

    description_values(8) = 'traverse_adaptive' // C_NULL_CHAR
    description_values_ptrs(8) = C_LOC(description_values(8))

    CALL extrae_define_event_type(1000, evt_desc, 8, values, description_values_ptrs)
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
