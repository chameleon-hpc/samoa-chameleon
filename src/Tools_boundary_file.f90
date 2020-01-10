#include "Compilation_control.f90"

! This module provides routines to read a csv time series and interpolate between the records

module Tools_boundary_file

    implicit none

    ! Needed to prevent module conflicts
#   if defined(_SINGLE_PRECISION)
        integer, PARAMETER :: BND_GRID_SR = kind(1.0e0)
#   elif defined(_DOUBLE_PRECISION)
        integer, PARAMETER :: BND_GRID_SR = kind(1.0d0)
#   elif defined(_QUAD_PRECISION)
        integer, PARAMETER :: BND_GRID_SR = kind(1.0q0)
#   else
#       error "No floating point precision is chosen!"
#   endif

    integer, parameter                              :: fdid = 73

    real(BND_GRID_SR), dimension(:, :), allocatable, save  :: file_data
    integer, save                                   :: file_data_length

    logical, save                                   :: boundary_file_success
    character(len = 256), save                      :: boundary_file_last_msg

    private
    public boundary_file_parse, boundary_file_get, boundary_file_success, boundary_file_last_msg

    contains

    ! Read a csv file and allocate an array to hold its data
    subroutine boundary_file_parse(path)
        character (len = *), intent(in)		:: path

        logical                             :: exists
        integer                             :: iostat, allocstat, i
        real(BND_GRID_SR)                   :: test

        boundary_file_last_msg(:) = " "
        boundary_file_success = .false.

        inquire(file=trim(path), exist=exists)
        if(exists) then
            file_data_length = 0 
            open(fdid, file=trim(path), action='read')
            do
                read(fdid, *, iostat=iostat)
                if (iostat .ne. 0) exit
                file_data_length = file_data_length + 1
            end do
            rewind(fdid)
            if(file_data_length .gt. 0) then
                allocate(file_data(2, file_data_length), stat=allocstat)
                file_data(:, :) = 0.0_BND_GRID_SR
                if(allocstat .eq. 0) then
                    read(fdid, *) file_data(:, :)
                    boundary_file_success = .true.
                    write (boundary_file_last_msg, "(A, I0, A, F0.4, A, F0.4)") "File read successfully, ", file_data_length, &
                        " records from t=", file_data(1, 1), " to t=", file_data(1, file_data_length)
                else
                    boundary_file_last_msg = "Cannot allocate memory"
                end if
            else
                boundary_file_last_msg = "File is empty"
            end if
            close(fdid)
        else
            boundary_file_last_msg = "File does not exist"
        end if
    end subroutine

    ! Interpolate on the previously read data
    function boundary_file_get(t) result(result)
        real (BND_GRID_SR), intent(in)			:: t
        real (BND_GRID_SR)                      :: result

        integer                                 :: i

        if(t .le. file_data(1, 1)) then
            result = file_data(2, 1)
        else if(t .ge. file_data(1, file_data_length)) then
            result = file_data(2, file_data_length)
        else
            do i = 1, file_data_length - 1
                if(t .ge. file_data(1, i) .and. t .le. file_data(1, i + 1)) then
                    result = file_data(2, i) + ((t - file_data(1, i)) * &
                        ((file_data(2, i + 1) - file_data(2, i)) / (file_data(1, i + 1) - file_data(1, i))))
                    exit
                end if
            end do
        end if
    end function

end module