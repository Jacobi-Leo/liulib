program test
  use mpi
  use mpiConstant
  implicit none

  real(kind=WP) :: a
  integer :: i, n, nproc, id

  call mpi_init ( IERR )
  call mpi_comm_size ( mpi_comm_world, nproc, IERR )
  call mpi_comm_rank ( mpi_comm_world, id, IERR )
  call createLIUREAL()

  if ( id == 0 ) then
     a = 1.0
     do i = 1, nproc - 1
        call mpi_send ( a, 1, LIUREAL, i, 0, mpi_comm_world, IERR )
     end do
  else
     call mpi_recv ( a, 1, LIUREAL, 0, 0, mpi_comm_world, mpi_status_ignore, IERR )
  end if

  a = a / 3.0
  write(*,*) "This is processor ", id, ", with ", a
end program

