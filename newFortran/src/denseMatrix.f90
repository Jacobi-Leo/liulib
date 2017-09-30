! Created by Jacob Zeyu LIU (liuzeyu271828@gmail.com)
!
! This module contains the definition of matrix and some basic operations,
! still under development.
!
! This code must be compiled with -fdefault-real-8 (in gfortran)
! or -autodouble (in ifort).
!
! Since Fortran 2003 features are massively utilized in this code,
! gfortran version>=5.0 is required (v5.4 has been verified),
! and ifort 17 is also verified while lower versions not tested (according to
! official document, ifort 16 is the first edition of ifort with full
! Fortran 2003 support).
!
! LAPACK is required by the solve routine, so -llapack should be added
! to compiler flag.

module denseMatrix
  use utils
  implicit none
  private

  public :: eye, solve, inv, writeMatrix, zeros, ones

  interface gesv
     subroutine sgesv ( N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       use utils, only: SGL
       integer, intent(in) :: N, NRHS, LDA, LDB
       integer, intent(out) :: INFO
       real(kind=SGL), intent(inout), dimension(LDA, N) :: A
       real(kind=SGL), intent(inout), dimension(LDB, NRHS) :: B
       integer, intent(inout), dimension(N) :: IPIV
     end subroutine sgesv
     subroutine dgesv ( N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       use utils, only: DBL
       integer, intent(in) :: N, NRHS, LDA, LDB
       integer, intent(out) :: INFO
       real(kind=DBL), intent(inout), dimension(LDA, N) :: A
       real(kind=DBL), intent(inout), dimension(LDB, NRHS) :: B
       integer, intent(inout), dimension(N) :: IPIV
     end subroutine dgesv
     subroutine zgesv ( N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       use utils, only: DBL
       integer, intent(in) :: N, NRHS, LDA, LDB
       integer, intent(out) :: INFO
       complex(kind=DBL), intent(inout), dimension(LDA, N) :: A
       complex(kind=DBL), intent(inout), dimension(LDB, NRHS) :: B
       integer, intent(inout), dimension(N) :: IPIV
     end subroutine zgesv
     subroutine cgesv ( N, NRHS, A, LDA, IPIV, B, LDB, INFO)
       use utils, only: SGL
       integer, intent(in) :: N, NRHS, LDA, LDB
       integer, intent(out) :: INFO
       complex(kind=SGL), intent(inout), dimension(LDA, N) :: A
       complex(kind=SGL), intent(inout), dimension(LDB, NRHS) :: B
       integer, intent(inout), dimension(N) :: IPIV
     end subroutine cgesv
  end interface gesv

  interface gels
     !! There is bug here.  The array WORK is assumed-size dummy array instead of
     !! assumed-shape dummy array.  This is REALLY a bug in modern Fortran.
     !! However, if assumed-shape dummy array used, workspace query is broken,
     !! which is wired, still can't figure out why.
     subroutine sgels ( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
       use utils, only: SGL
       character, intent(in) :: TRANS
       integer, intent(in) :: M, N, NRHS, LDA, LDB, LWORK
       integer, intent(out) :: INFO
       real(kind=SGL), intent(inout), dimension(LDA, N) :: A
       real(kind=SGL), intent(inout), dimension(LDB, NRHS) :: B
       real(kind=SGL), intent(inout), dimension(*) :: WORK
     end subroutine sgels
     subroutine dgels ( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
       use utils, only: DBL
       character, intent(in) :: TRANS
       integer, intent(in) :: M, N, NRHS, LDA, LDB, LWORK
       integer, intent(out) :: INFO
       real(kind=DBL), intent(inout), dimension(LDA, N) :: A
       real(kind=DBL), intent(inout), dimension(LDB, NRHS) :: B
       real(kind=DBL), intent(inout), dimension(*) :: WORK
     end subroutine dgels
     subroutine cgels ( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
       use utils, only: SGL
       character, intent(in) :: TRANS
       integer, intent(in) :: M, N, NRHS, LDA, LDB, LWORK
       integer, intent(out) :: INFO
       complex(kind=SGL), intent(inout), dimension(LDA, N) :: A
       complex(kind=SGL), intent(inout), dimension(LDB, NRHS) :: B
       real(kind=SGL), intent(inout), dimension(*) :: WORK
     end subroutine cgels
     subroutine zgels ( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
       use utils, only: DBL
       character, intent(in) :: TRANS
       integer, intent(in) :: M, N, NRHS, LDA, LDB, LWORK
       integer, intent(out) :: INFO
       complex(kind=DBL), intent(inout), dimension(LDA, N) :: A
       complex(kind=DBL), intent(inout), dimension(LDB, NRHS) :: B
       real(kind=DBL), intent(inout), dimension(*) :: WORK
     end subroutine zgels
  end interface gels

  interface getri
     subroutine sgetri ( N, A, LDA, IPIV, WORK, LWORK, INFO )
       use utils, only: SGL
       integer, intent(in) :: N, LDA, LWORK
       integer, intent(out) ::  INFO
       real(kind=SGL), dimension( lda, n ) :: A
       integer, dimension( n ), intent(inout) ::  IPIV
       real(kind=SGL), dimension( * ), intent(inout) :: WORK
     end subroutine sgetri
     subroutine dgetri ( N, A, LDA, IPIV, WORK, LWORK, INFO )
       use utils, only: DBL
       integer, intent(in) :: N, LDA, LWORK
       integer, intent(out) ::  INFO
       real(kind=DBL), dimension( lda, n ) :: A
       integer, dimension( n ), intent(inout) ::  IPIV
       real(kind=DBL), dimension( * ), intent(inout) :: WORK
     end subroutine dgetri
     subroutine zgetri ( N, A, LDA, IPIV, WORK, LWORK, INFO )
       use utils, only: DBL
       integer, intent(in) :: N, LDA, LWORK
       integer, intent(out) ::  INFO
       complex(kind=DBL), dimension( lda, n ) ::  A
       integer, dimension( n ), intent(inout) ::  IPIV
       complex(kind=DBL), dimension( * ), intent(inout) ::  WORK
     end subroutine zgetri
     subroutine cgetri ( N, A, LDA, IPIV, WORK, LWORK, INFO )
       use utils, only: SGL
       integer, intent(in) :: N, LDA, LWORK
       integer, intent(out) ::  INFO
       complex(kind=SGL), dimension( lda, n ) ::  A
       integer, dimension( n ), intent(inout) ::  IPIV
       complex(kind=SGL), dimension( * ), intent(inout) ::  WORK
     end subroutine cgetri
  end interface getri

  interface getrf
     subroutine sgetrf ( M, N, A, LDA, IPIV, INFO )
       use utils, only: SGL
       integer, intent(in) :: M, N, LDA
       integer, intent(out) ::  INFO
       real(kind=SGL), dimension( lda, n ) :: A
       integer, dimension( n ), intent(inout) ::  IPIV
     end subroutine sgetrf
     subroutine dgetrf ( M, N, A, LDA, IPIV, INFO )
       use utils, only: DBL
       integer, intent(in) :: M, N, LDA
       integer, intent(out) ::  INFO
       real(kind=DBL), dimension( lda, n ) :: A
       integer, dimension( n ), intent(inout) ::  IPIV
     end subroutine dgetrf
     subroutine cgetrf ( M, N, A, LDA, IPIV, INFO )
       use utils, only: SGL
       integer, intent(in) :: M, N, LDA
       integer, intent(out) ::  INFO
       complex(kind=SGL), dimension( lda, n ) ::  A
       integer, dimension( n ), intent(inout) ::  IPIV
     end subroutine cgetrf
     subroutine zgetrf ( M, N, A, LDA, IPIV, INFO )
       use utils, only: DBL
       integer, intent(in) :: M, N, LDA
       integer, intent(out) ::  INFO
       complex(kind=DBL), dimension( lda, n ) ::  A
       integer, dimension( n ), intent(inout) ::  IPIV
     end subroutine zgetrf
  end interface getrf

  interface solve
     module procedure solve1
     module procedure solve2
  end interface solve

contains

  function solve1 ( A, b ) result ( x )
    real(kind=WP), dimension(:,:), intent(in) :: A
    real(kind=WP), dimension(:), intent(in) :: b
    real(kind=WP), dimension(size(A,2),1) :: x

    real(kind=WP), dimension(size(A,1), size(A,2)) :: tmpA
    ! real(kind=WP), allocatable, dimension(:) :: work
    integer, allocatable, dimension(:) :: ipiv

    integer :: info

    if ( size(A, 1) /= size(b, 1) ) then
       stop 'dimension mismatch in solving linear system!'
    end if

    if ( size(A, 1) == size(A, 2) ) then
       allocate( ipiv(size(A, 1)) )
       tmpA = A
       x(:,1) = b
       call gesv ( size(A, 1), 1, tmpA, size(A, 2), ipiv, x, size(x, 1), info )
    else
       stop 'linear fitting is not implemented...'
    end if
    ! not finished
  end function solve1

  function solve2 ( A, b ) result ( x )
    real(kind=WP), dimension(:,:), intent(in) :: A, b
    real(kind=WP), dimension(size(A,2), size(b,2)) :: x

    real(kind=WP), dimension(size(A,1), size(A,2)) :: tmpA
    ! real(kind=WP), allocatable, dimension(:) :: work
    integer, allocatable, dimension(:) :: ipiv

    integer :: info

    if ( size(A, 1) /= size(b, 1) ) then
       stop 'dimension mismatch in solving linear system!'
    end if

    if ( size(A, 1) == size(A, 2) ) then
       allocate( ipiv(size(A, 1)) )
       tmpA = A
       x = b
       call gesv ( size(A, 1), size(b,2), tmpA, size(A, 2), ipiv, x, size(x, 1), info )
    else
       stop 'linear fitting is not implemented...'
    end if
    ! not finished
  end function solve2

  function eye ( n )
    integer, intent(in) :: n
    real(kind=WP), dimension(n, n) :: eye
    integer :: i
    eye = 0.0_WP
    forall (i=1:n)
       eye(i,i) = 1.0_WP
    end forall
  end function eye

  function zeros ( n )
    integer, intent(in) :: n
    real(kind=WP), dimension(n, n) :: eye
    zeros = 0.0_WP
  end function zeros

  function ones ( n )
    integer, intent(in) :: n
    real(kind=WP), dimension(n, n) :: eye
    eye = 1.0_WP
  end function ones

  function inv(A) result(Ainv)
    real(kind=WP), dimension(:,:), intent(in) :: A
    real(kind=WP), dimension(size(A,1),size(A,2)) :: Ainv
    real(kind=WP), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call GETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call GETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end function inv

  subroutine writeMatrix ( A, unit )
    real(kind=WP), dimension(:,:), intent(in) :: A
    integer, optional, intent(in) :: unit

    integer :: i, j

    if ( present(unit) ) then
       open( unit=unit, file=outputName(unit) )
       write(unit, *) ((A(i,j), j=1, size(A,2)), i=1, size(A,1))
    else
       write(*, *) ((A(i,j), j=1, size(A,2)), new_line('*'), i=1, size(A,1))
    end if
  end subroutine writeMatrix

end module denseMatrix
