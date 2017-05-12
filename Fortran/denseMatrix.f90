! Created by Jacob Zeyu LIU (liuzeyu271828@gmail.com)
!
! This module contains the definition of matrix and some basic operations,
! still under development.
!
! This code must be compiled with (in gfortran) -fdefault-real-8 and -std=f2003,
! or -autodouble (in ifort).
!
! Since Fortran 2003 features are massively utilized in this code,
! gfortran version>=5.0 is required (v6.0 or higher is recommended),
! and ifort 17 is verified while lower versions not tested.


module denseMatrix
  use constant
  implicit none
  private

  !======================================================================
  ! Error Information List
  !   info = 127    The previous operation has unmatched dimension
  !======================================================================
  type, public :: Matrix
     private
     integer :: nrow = 0  ! number of rows
     integer :: ncol = 0  ! number of columns
     real(kind=WP), allocatable :: comp(:,:)  ! components of the matrix
     integer :: info = 0  ! error information
     logical :: diagonal = .false.
     logical :: bidiagonal = .false.
   contains
     !!===== Assignment and operators, as well as related. =======
     generic, public :: assignment(=) => arrayToMatrix, dArrayToMatrix
     generic, public :: operator(*) => matrixTimesReal, matrixTimesInt, realTimesMatrix, intTimesMatrix, matrixTimesMatrix
     generic, public :: operator(+) => matrixAdd, arrayAdd
     generic, public :: operator(-) => matrixSubtract
     procedure, private, pass :: arrayToMatrix
     procedure, private, pass :: dArrayToMatrix
     procedure, private, pass :: matrixAdd
     procedure, private, pass :: arrayAdd
     procedure, public, pass :: addMatrix
     procedure, private, pass :: matrixSubtract
     procedure, private, pass(self) :: realTimesMatrix
     procedure, private, pass :: matrixTimesreal
     procedure, private, pass(self) :: intTimesMatrix
     procedure, private, pass :: matrixTimesint
     procedure, private, pass :: matrixTimesMatrix
     !!===========================================================
     procedure, public, pass :: T  ! transpose
     procedure, public, pass :: writeToFile
     procedure, public, pass :: isAllocated
     procedure, public, pass :: isDiagonal
     procedure, public, pass :: isBidiagonal
     procedure, public, pass :: getNrow
     procedure, public, pass :: getNcolumn
     procedure, public, pass :: getEllement
     ! procedure, pass :: inv
     ! procedure, pass :: solve
     final :: matrixClean
  end type matrix

contains

  function isDiagonal ( self ) result ( r )
    class(Matrix), intent(in) :: self
    logical :: r
    r = self%diagonal
  end function isDiagonal

  function isBidiagonal ( self ) result ( r )
    class(Matrix), intent(in) :: self
    logical :: r
    r = self%bidiagonal
  end function isBidiagonal

  function getEllement ( self, i, j ) result ( r )
    class(Matrix), intent(in) :: self
    real(kind=WP) :: r
    integer, intent(in) :: i, j
    r = self%comp(i,j)
  end function getEllement

  function getNrow ( self ) result ( r )
    integer :: r
    class(Matrix), intent(in) :: self
    r = self%nrow
  end function getNrow

  function getNcolumn (self ) result ( r )
    integer :: r
    class(Matrix), intent(in) :: self
    r = self%ncol
  end function getNcolumn

  function isAllocated ( self ) result ( r )
    logical :: r
    class(Matrix), intent(in) :: self
    r = allocated( self%comp )
  end function isAllocated

  subroutine arrayToMatrix ( mat, array )
    class(Matrix), intent(out) :: mat
    real(kind=WP), intent(in), dimension(:,:) :: array
    integer :: m, n
    m = size(array, 1)
    n = size(array, 2)
    allocate( mat%comp( m, n ) )
    mat%comp = array
    mat%nrow = m
    mat%ncol = n
    mat%info = 0
  end subroutine arrayToMatrix

  subroutine dArrayToMatrix ( mat, array )
    class(Matrix), intent(out) :: mat
    real(kind=WP), intent(in), dimension(:) :: array
    integer :: m
    m = size(array)
    allocate( mat%comp( m, 1 ) )
    mat%comp(:,1) = array
    mat%nrow = m
    mat%ncol = 1
    mat%info = 0
  end subroutine dArrayToMatrix

  function arrayAdd ( self, array )
    type(Matrix) :: arrayAdd
    class(Matrix), intent(in) :: self
    type(Matrix) :: here
    real(kind=WP), intent(in), dimension(:,:) :: array

    here = array
    arrayAdd = self + here
  end function arrayAdd

  function matrixAdd ( self, mat )
    type(Matrix) :: matrixAdd
    class(Matrix), intent(in) :: self, mat
    integer :: m1, m2, n1, n2

    m1 = self%nrow
    n1 = self%ncol
    m2 = mat%nrow
    n2 = mat%ncol

    if ( m1 /= m2 .or. n1 /= n2 ) then
       matrixAdd%nrow = -1
       matrixAdd%ncol = -1
       matrixAdd%info = 127
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
       self%info = 127
    else
       self%comp = self%comp + mat%comp
    end if
  end subroutine addMatrix

  function matrixSubtract ( self, mat )
    type(Matrix) :: matrixSubtract
    class(Matrix), intent(in) :: self, mat
    real(kind=WP) :: neg = -1.0

    matrixSubtract = self + neg * mat
  end function matrixSubtract

  function realTimesMatrix ( a, self )
    type(Matrix) :: realTimesMatrix
    real(kind=WP), intent(in) :: a
    class(Matrix), intent(in) :: self

    realTimesMatrix%nrow = self%nrow
    realTimesMatrix%ncol = self%ncol

    if ( allocated(self%comp) ) then
       realTimesMatrix%comp = a * self%comp
    end if
  end function realTimesMatrix

  function matrixTimesReal ( self, a )
    type(Matrix) :: matrixTimesReal
    real(kind=WP), intent(in) :: a
    class(Matrix), intent(in) :: self

    matrixTimesReal%nrow = self%nrow
    matrixTimesReal%ncol = self%ncol

    if ( allocated(self%comp) ) then
       matrixTimesReal%comp = a * self%comp
    end if
  end function matrixTimesReal

  function intTimesMatrix ( i, self )
    type(Matrix) :: intTimesMatrix
    integer, intent(in) :: i
    real(kind=WP) :: a
    class(Matrix), intent(in) :: self

    a = real(i, WP)
    intTimesMatrix = a * self
  end function intTimesMatrix

  function matrixTimesInt ( self, i )
    type(Matrix) :: matrixTimesInt
    integer, intent(in) :: i
    real(kind=WP) :: a
    class(Matrix), intent(in) :: self

    a = real(i, WP)
    matrixTimesInt = a * self
  end function matrixTimesInt

  function matrixTimesMatrix ( self, mat )
    type(Matrix) :: matrixTimesMatrix
    class(Matrix), intent(in) :: self, mat

    if ( self%ncol /= mat%nrow .or. allocated( self%comp ) == .false. .or. allocated( mat%comp ) == .false. ) then
       matrixTimesMatrix%info = 127
    else
       matrixTimesMatrix%nrow = self%nrow
       matrixTimesMatrix%ncol = mat%ncol
       matrixTimesMatrix%comp = matmul( self%comp, mat%comp )
    end if
  end function matrixTimesMatrix

  subroutine T ( self )
    class(Matrix), intent(inout) :: self
    integer :: tmp

    if ( allocated(self%comp) ) then
       tmp = self%nrow
       self%nrow = self%ncol
       self%ncol = tmp
       self%comp = transpose(self%comp)
    end if

  end subroutine T

  subroutine writeToFile ( self, fileUnit )
    class(Matrix), intent(in) :: self
    integer, intent(in) :: fileUnit
    integer :: i, j

    !! This is a bad solution
    !!========================================================
    ! do i = 1, self%nrow
    !    do j = 1, self%ncol
    !       write(fileUnit, *) self%comp(i,j)
    !    end do
    !    write(fileUnit,*)
    ! end do

    ! exactly the way I want
    do i = 1, self%nrow
       write(fileUnit, *) ( self%comp(i,j), j = 1, self%ncol )
    end do

  end subroutine writeToFile

  subroutine matrixClean ( self )
    type(Matrix), intent(inout) :: self

    !! This is debug info
    !write(*,*) "In finalizer"

    if ( allocated( self%comp ) ) then
       deallocate( self%comp )
    end if
  end subroutine matrixClean

end module denseMatrix
