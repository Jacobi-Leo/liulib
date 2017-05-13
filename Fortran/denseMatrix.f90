! Created by Jacob Zeyu LIU (liuzeyu271828@gmail.com)
!
! This module contains the definition of matrix and some basic operations,
! still under development.
!
! This code must be compiled with -fdefault-real-8 (in gfortran)
! or -autodouble (in ifort).
!
! Since Fortran 2003 features are massively utilized in this code,
! gfortran version>=5.0 is required (v6.0 or higher is recommended),
! and ifort 17 is verified while lower versions not tested (according to
! official document, ifort 16 is the first edition of ifort with full
! Fortran 2003 support).
!
! LAPACK is required by the solve routine, so -llapack should be added
! to compiler flag.

module denseMatrix
  use constant
  implicit none
  private

  !======================================================================
  ! Error Information List
  !   info = 127    The previous operation has unmatched dimension
  !   info = 126    Trying to use unallocated components
  !======================================================================
  ! Error information is used for debugging and testing whether the
  ! previous operation is valid, never judge the feasibility of one
  ! operation with attribute /info/
  !======================================================================
  type, public :: Matrix
     private
     integer :: nrow = 0  ! number of rows
     integer :: ncol = 0  ! number of columns
     integer :: nupper = 0  ! number of superdiagonals
     integer :: nlower = 0  ! number of subdiagonals
     real(kind=WP), allocatable :: comp(:,:)  ! components of the matrix
     real(kind=WP), allocatable :: band(:,:)  ! form of band matrix
     integer :: info = 0  ! error information
     logical :: diagonal = .false.
     logical :: bidiagonal = .false.
     logical :: tridiagonal = .false.
     logical :: hessenberg = .false. ! upper Hessenberg considered only
     logical :: banded = .false.     ! not implemented
     logical :: symmetric = .false.
     logical :: hermitian = .false.  ! not implemented
     logical :: orthogonal = .false. ! not implemented
     logical :: positiveDefinite = .false. ! not implemented
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
     procedure, public, pass :: addMatrix     ! this is a subroutine, and other procedures are functions
     procedure, private, pass :: matrixSubtract
     procedure, private, pass(self) :: realTimesMatrix
     procedure, private, pass :: matrixTimesreal
     procedure, private, pass(self) :: intTimesMatrix
     procedure, private, pass :: matrixTimesint
     procedure, private, pass :: matrixTimesMatrix
     !!===========================================================
     ! predicates begins with /is/, and are functions
     ! procedures begins with /set/ are subroutines
     ! procedures begins with /get/ are functions
     procedure, public, pass :: isAllocated
     procedure, public, pass :: isDiagonal
     procedure, public, pass :: isBidiagonal
     procedure, public, pass :: isTridiagonal
     procedure, public, pass :: isHessenberg
     procedure, public, pass :: isSymmetric
     procedure, public, pass :: setDiagonal
     procedure, public, pass :: setBidiagonal
     procedure, public, pass :: setTridiagonal
     procedure, public, pass :: setHessenberg
     procedure, public, pass :: setSymmetric
     procedure, public, pass :: getNrow
     procedure, public, pass :: getNcolumn
     procedure, public, pass :: getEllement
     ! The following to procedures are subroutines
     procedure, public, pass :: printSpecialAttributes   ! not finished
     procedure, public, pass :: resetToGeneral           ! not finished
     !!===========================================================
     ! manipulate the shape of matrix, all subroutines
     procedure, public, pass :: pushRow
     procedure, public, pass :: popRow
     procedure, public, pass :: pushColumn
     procedure, public, pass :: popColumn
     procedure, public, pass :: T  ! transpose
     !!===========================================================
     ! subroutine to write components into file in matrix form
     procedure, public, pass :: writeToFile
     !!===========================================================
     ! linear solver
     ! procedure, pass :: inv
     ! procedure, public, pass :: solve
     !!===========================================================
     final :: matrixClean
  end type matrix

  interface gesv
     subroutine sgesv ( N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       use constant
       integer, intent(in) :: N, NRHS, LDA, LDB
       integer, intent(out) :: INFO
       real(kind=SGL), intent(inout), dimension(LDA, N) :: A
       real(kind=SGL), intent(inout), dimension(LDB, NRHS) :: B
       integer, intent(inout), dimension(N) :: IPIV
     end subroutine sgesv
     subroutine dgesv ( N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       use constant
       integer, intent(in) :: N, NRHS, LDA, LDB
       integer, intent(out) :: INFO
       real(kind=DBL), intent(inout), dimension(LDA, N) :: A
       real(kind=DBL), intent(inout), dimension(LDB, NRHS) :: B
       integer, intent(inout), dimension(N) :: IPIV
     end subroutine dgesv
     subroutine zgesv ( N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       use constant
       integer, intent(in) :: N, NRHS, LDA, LDB
       integer, intent(out) :: INFO
       complex(kind=DBL), intent(inout), dimension(LDA, N) :: A
       complex(kind=DBL), intent(inout), dimension(LDB, NRHS) :: B
       integer, intent(inout), dimension(N) :: IPIV
     end subroutine zgesv
     subroutine cgesv ( N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       use constant
       integer, intent(in) :: N, NRHS, LDA, LDB
       integer, intent(out) :: INFO
       complex(kind=SGL), intent(inout), dimension(LDA, N) :: A
       complex(kind=SGL), intent(inout), dimension(LDB, NRHS) :: B
       integer, intent(inout), dimension(N) :: IPIV
     end subroutine cgesv
  end interface gesv

contains

  subroutine solve ( self, b )  !! not finished
    class(Matrix), intent(inout) :: self, b
    integer, allocatable, dimension(:) :: ipiv

    if ( (.not. self%isAllocated()) .or. (.not. b%isAllocated()) ) then
       self%info = 126
       return
    else if ( self%ncol /= b%nrow ) then
       self%info = 127
       return
    end if

    call self%resetToGeneral()
    call b%resetToGeneral()

    associate( A => self%comp, bb => b%comp, info => self%info, &
         m => self%nrow, n => self%ncol, nrhs => b%ncol )
      if ( self%banded ) then
!!!
         return
      else
         !! The last case
         allocate( ipiv(n) )
         call gesv( n, nrhs, A, m, ipiv, bb, n, info )
      end if
    end associate
  end subroutine solve

  subroutine pushRow ( self, row )
    class(Matrix), intent(inout) :: self
    real(kind=WP), intent(in), dimension(:) :: row
    real(kind=WP), dimension(:,:), allocatable :: tmp
    if ( self%ncol /= size(row) ) then
       self%info = 127
       return
    else
       call self%resetToGeneral()
       allocate( tmp(self%nrow+1, self%ncol) )
       tmp(1:self%nrow, 1:self%ncol) = self%comp
       tmp(self%nrow+1, :) = row
       self = tmp
       deallocate( tmp )
    end if
  end subroutine pushRow

  subroutine popRow ( self )
    class(Matrix), intent(inout) :: self
    real(kind=WP), allocatable, dimension(:,:) :: tmp

    call self%resetToGeneral()
    allocate( tmp(self%nrow-1, self%ncol) )
    tmp = self%comp(1:self%nrow-1, 1:self%ncol)
    self = tmp
    deallocate( tmp )
  end subroutine popRow

  subroutine pushColumn ( self, column )
    class(Matrix), intent(inout) :: self
    real(kind=WP), intent(in), dimension(:) :: column
    call self%resetToGeneral()
    if ( self%nrow /= size(column) ) then
       self%info = 127
       return
    else
       call self%T()
       call self%pushRow(column)
       call self%T()
    end if
  end subroutine pushColumn

  subroutine popColumn ( self )
    class(Matrix), intent(inout) :: self
    call self%resetToGeneral()
    call self%T()
    call self%popRow()
    call self%T()
  end subroutine popColumn

  !================ below are not implemented ====================
  subroutine printSpecialAttributes ( self )
    class(Matrix), intent(in) :: self
  end subroutine printSpecialAttributes

  subroutine resetToGeneral ( self )
    class(Matrix), intent(inout) :: self

    if ( allocated( self%band ) ) then
       deallocate( self%band )
    end if
    !! more to do
  end subroutine resetToGeneral

  !============ The above procedures are not finished =================

  subroutine setHessenberg ( self, status )
    class(Matrix), intent(inout) :: self
    logical, intent(in) :: status
    self%hessenberg = status
  end subroutine setHessenberg

  subroutine setTridiagonal ( self, status )
    class(Matrix), intent(inout) :: self
    logical, intent(in) :: status
    self%tridiagonal = status
  end subroutine setTridiagonal

  subroutine setDiagonal ( self, status )
    class(Matrix), intent(inout) :: self
    logical, intent(in) :: status
    self%diagonal = status
  end subroutine setDiagonal

  subroutine setBidiagonal ( self, status )
    class(Matrix), intent(inout) :: self
    logical, intent(in) :: status
    self%bidiagonal = status
  end subroutine setBidiagonal

  subroutine setSymmetric ( self, status )
    class(Matrix), intent(inout) :: self
    logical, intent(in) :: status
    self%bidiagonal = status
  end subroutine setSymmetric

  function isSymmetric ( self ) result ( r )
    class(Matrix), intent(in) :: self
    logical :: r
    r = self%symmetric
  end function isSymmetric

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

  function isTridiagonal ( self ) result ( r )
    class(Matrix), intent(in) :: self
    logical :: r
    r = self%tridiagonal
  end function isTridiagonal

  function isHessenberg ( self ) result ( r )
    class(Matrix), intent(in) :: self
    logical :: r
    r = self%hessenberg
  end function isHessenberg

  function getEllement ( self, i, j ) result ( r )
    class(Matrix), intent(in) :: self
    real(kind=WP) :: r
    integer, intent(in) :: i, j
    if ( self%isAllocated() ) then
       r = self%comp(i,j)
    else
       stop '126'
    end if
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

    if ( mat%isAllocated() ) then
       deallocate( mat%comp )
    end if

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

    if ( mat%isAllocated() ) then
       deallocate( mat%comp )
    end if

    m = size(array)
    allocate( mat%comp( m, 1 ) )
    mat%comp(:,1) = array
    mat%nrow = m
    mat%ncol = 1
    mat%info = 0
  end subroutine dArrayToMatrix

  function matrixAdd ( self, mat ) result ( r )
    !==========================================
    ! If the dimension of self and mat does not
    ! match, program stopped.
    ! If the parameter self is not allocated, mat
    ! is returned with error info.
    ! If the parameter mat is not allocated, self
    ! is returned with error info.
    !
    ! The following three procedures are based on
    ! this one.
    !==========================================
    type(Matrix) :: r
    class(Matrix), intent(in) :: self, mat
    integer :: m1, m2, n1, n2

    m1 = self%nrow
    n1 = self%ncol
    m2 = mat%nrow
    n2 = mat%ncol

    if ( .not. self%isAllocated() ) then
       r = mat
       r%info = 126
    else if ( .not. mat%isAllocated() ) then
       r = self
       r%info = 126
    else if ( m1 /= m2 .or. n1 /= n2 ) then
       r%info = 127
       stop '127'
    else
       r%nrow = m1
       r%ncol = n1
       allocate( r%comp(m1,n1) )
       r%comp = self%comp + mat%comp
    end if
  end function matrixAdd

  function arrayAdd ( self, array )
    type(Matrix) :: arrayAdd
    class(Matrix), intent(in) :: self
    type(Matrix) :: here
    real(kind=WP), intent(in), dimension(:,:) :: array

    here = array
    arrayAdd = self + here
  end function arrayAdd

  subroutine addMatrix ( self, mat )
    class(Matrix), intent(in) :: mat
    class(Matrix), intent(inout) :: self
    integer :: m1, m2, n1, n2

    if ( .not. self%isAllocated() ) then
       self%info = 126
       return
    else if ( self%nrow /= mat%nrow .or. self%ncol /= mat%ncol ) then
       self%info = 127
       return
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

  function realTimesMatrix ( a, self ) result ( r )
    !=================================================
    ! This function is the basic scalar multiplication
    ! function and all other three ones just call this
    ! function.
    !=================================================
    type(Matrix) :: r
    real(kind=WP), intent(in) :: a
    class(Matrix), intent(in) :: self

    if ( .not. self%isAllocated() ) then
       r%info = 126
       return
    end if

    r%nrow = self%nrow
    r%ncol = self%ncol

    if ( self%isAllocated() ) then
       r%comp = a * self%comp
    end if
  end function realTimesMatrix

  function matrixTimesReal ( self, a )
    type(Matrix) :: matrixTimesReal
    real(kind=WP), intent(in) :: a
    class(Matrix), intent(in) :: self

    matrixTimesReal%nrow = self%nrow
    matrixTimesReal%ncol = self%ncol

    if ( self%isAllocated() ) then
       matrixTimesReal = a * self
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

    if ( self%ncol /= mat%nrow ) then
       matrixTimesMatrix%info = 127
       return
    else if ( ( .not. self%isAllocated() ) .or. ( .not. mat%isAllocated() ) ) then
       matrixTimesMatrix%info = 126
       return
    else
       matrixTimesMatrix%nrow = self%nrow
       matrixTimesMatrix%ncol = mat%ncol
       matrixTimesMatrix%comp = matmul( self%comp, mat%comp )
    end if
  end function matrixTimesMatrix

  subroutine T ( self )
    class(Matrix), intent(inout) :: self
    integer :: tmp

    if ( ( .not. self%isAllocated() ) .or. self%isSymmetric() ) then
       return
    else
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

    if ( .not. self%isAllocated() ) then
       return
    else
       !! This is a bad solution
       !!==========================================!
       ! do i = 1, self%nrow                       !
       !    do j = 1, self%ncol                    !
       !       write(fileUnit, *) self%comp(i,j)   !
       !    end do                                 !
       !    write(fileUnit,*)                      !
       ! end do                                    !
       !===========================================!
       !! exactly the way I want
       do i = 1, self%nrow
          write(fileUnit, *) ( self%comp(i,j), j = 1,self%ncol )
       end do
    end if

  end subroutine writeToFile

  subroutine matrixClean ( self )
    type(Matrix), intent(inout) :: self
    !! This is debug info
    !write(*,*) "In finalizer"

    if ( self%isAllocated() ) then
       deallocate( self%comp )
    end if
  end subroutine matrixClean

end module denseMatrix



