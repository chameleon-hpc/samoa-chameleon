program matrices
  use gen_matrices
  use elementary

  character(len=1024) :: n_arg,basis_arg
  integer :: N

!real(kind=GRID_SR),allocatable :: identity(:,:)

call get_command_argument(1,n_arg)
call get_command_argument(2,basis_arg)

write(*,*), "using order: ",trim(n_arg)
write(*,*), "basis: ",trim(basis_arg)

read(n_arg,'(I1)') N

basis_arg=trim(basis_arg)

call settings(N,basis_arg)

call generate_matrices(N)

end program matrices


