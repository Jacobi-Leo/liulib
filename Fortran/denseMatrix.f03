! Created by Jacob Zeyu LIU (liuzeyu271828@gmail.com)
! This module contains the definition of matrix and some basic
! operations, still under development

module denseMatrix
  use constant

  type :: Matrix
     integer :: nrow  ! number of rows
     integer :: ncol  ! number of columns
     real(kind=WP), dimension(nrow, ncol) :: comp  ! components of the matrix
   contains
     generic :: assignment(=) => matrixToArray, arrayToMatrix
     ! generic :: operator(+) => matrixAdd 
     ! generic :: operator(-) => matrixSubtract
     ! generic :: operator(*) => matrixTimesReal, matrixTimesInt, realTimesMatrix, intTimesMatrix, matrixTimesMatrix
     ! procedure, pass :: inv
  end type matrix

contains

  subroutine matrixToArray ( array, matrix )
    type(Matrix), intent(in) :: matrix
    real(kind=WP), intent(out), dimension(matrix%nrow, matrix%ncol) :: array
    array = matrix%coeff
  end subroutine matrixToArray

  subroutine arrayToMatrix ( matrix, array )
    type(Matrix), intent(out) :: matrix
    real(kind=WP), intent(in), dimension(:,:) :: array
    matrix%coeff = array
    matrix%nrow = size(array)
    matrix%ncol = size(array, 2)
  end subroutine arrayToMatrix

end module denseMatrix
