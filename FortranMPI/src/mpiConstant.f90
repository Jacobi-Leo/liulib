! Created by Jacob Zeyu LIU (liuzeyu271828@gmail.com)

module mpiConstant
  use constant
  use mpi
  implicit none
  save

  integer :: LIUREAL, IERR

contains

  subroutine createLIUREAL()
    call mpi_type_create_f90_real(precision(1.0_WP), range(1.0_WP), LIUREAL, IERR)
  end subroutine createLIUREAL

end module mpiConstant
