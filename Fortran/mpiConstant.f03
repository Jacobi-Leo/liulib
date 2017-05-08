! Created by Jacob Zeyu LIU (liuzeyu271828@gmail.com)

module mpiConstant
  use constant
  use mpi
  implicit none
  private

  integer :: LIUREAL
  public :: LIUREAL, mpiConstant_init

contains

  subroutine mpiConstant_init()
    integer :: IERROR
    mpi_type_create_f90_real(precision(1.0), range(1.0), LIUREAL, IERROR)
  end subroutine mpiConstant_init 

end module mpiConstant
