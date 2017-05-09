! Created by Jacob Zeyu LIU (liuzeyu271828@gmail.com)
! This module contains the definition of matrix and some basic
! operations, still under development

module denseMatrix
  use constant

  type :: Matrix
     integer :: nrow  ! number of rows
     integer :: ncol  ! number of columns
     real(kind=WP), allocatable, dimension(:, :) :: comp  ! components of the matrix
   contains
     procedure, nopass :: arrayToMatrix
     procedure, nopass :: matrixToArray
     generic :: assignment(=) => matrixToArray, arrayToMatrix
     ! generic :: operator(+) => matrixAdd
     ! generic :: operator(-) => matrixSubtract
     ! generic :: operator(*) => matrixTimesReal, matrixTimesInt, realTimesMatrix, intTimesMatrix, matrixTimesMatrix
     ! procedure, pass :: inv
     ! procedure, pass :: solve
  end type matrix

contains

  subroutine matrixToArray ( array, matrix )
    type(Matrix), intent(in) :: matrix
    real(kind=WP), intent(out), dimension(matrix%nrow, matrix%ncol) :: array
    array = matrix%comp
  end subroutine matrixToArray

  subroutine arrayToMatrix ( matrix, array )
    type(Matrix), intent(out) :: matrix
    real(kind=WP), intent(in), dimension(:,:) :: array
    integer :: m, n
    m = size(array)
    n = size(array, 2)
    allocate( matrix%comp( m, n ) )
    matrix%comp = array
    matrix%nrow = m
    matrix%ncol = n
  end subroutine arrayToMatrix

end module denseMatrix
