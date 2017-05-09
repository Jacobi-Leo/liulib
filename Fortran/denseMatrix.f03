! Created by Jacob Zeyu LIU (liuzeyu271828@gmail.com)
!
! This module contains the definition of matrix and some basic
! operations, still under development
!
! This code must be compiled with (in gfortran) -fdefault-real-8 and -std=f2003

module denseMatrix
  use constant
  implicit none
  private

  !======================================================================
  ! Error Information List
  !   info = 11    The previous operation has unmatched dimension
  !======================================================================
  type, public :: Matrix
     integer :: nrow  ! number of rows
     integer :: ncol  ! number of columns
     real(kind=WP), allocatable :: comp(:,:)  ! components of the matrix
     integer :: info  ! error information
   contains
     generic, public :: assignment(=) => arrayToMatrix
     procedure, private, pass :: arrayToMatrix
     generic, public :: operator(+) => matrixAdd
     procedure, private, pass :: matrixAdd
     ! generic :: operator(-) => matrixSubtract
     ! generic :: operator(*) => matrixTimesReal, matrixTimesInt, realTimesMatrix, intTimesMatrix, matrixTimesMatrix
     procedure, public, pass :: addMatrix
     ! procedure, pass :: inv
     ! procedure, pass :: solve
  end type matrix

contains

  subroutine arrayToMatrix ( mat, array )
    class(Matrix), intent(out) :: mat
    real(kind=WP), intent(in), dimension(:,:) :: array
    integer :: m, n
    m = size(array)
    n = size(array, 2)
    allocate( mat%comp( m, n ) )
    mat%comp = array
    mat%nrow = m
    mat%ncol = n
    mat%info = 0
  end subroutine arrayToMatrix

  type(Matrix) function matrixAdd ( self, mat )
    class(Matrix), intent(in) :: self, mat
    integer :: m1, m2, n1, n2

    m1 = self%nrow
    n1 = self%ncol
    m2 = mat%nrow
    n2 = mat%ncol

    if ( m1 /= m2 .or. n1 /= n2 ) then
       matrixAdd%info = 11
    else
       matrixAdd%nrow = m1
       matrixAdd%ncol = n1
       allocate( matrixAdd%comp(m1,n1) )
       matrixAdd%comp = self%comp + mat%comp
    end if
  end function matrixAdd

  subroutine addMatrix ( self, mat )
    class(Matrix), intent(in) :: mat
    class(Matrix), intent(inout) :: self
    integer :: m1, m2, n1, n2

    m1 = self%nrow
    n1 = self%ncol
    m2 = mat%nrow
    n2 = mat%ncol

    if ( m1 /= m2 .or. n1 /= n2 ) then
       self%info = 11
    else
       self%comp = self%comp + mat%comp
    end if
  end subroutine addMatrix

end module denseMatrix
