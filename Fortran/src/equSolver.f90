module equSolver
  !! buggy, not usable.
  use constant
  implicit none
  save
  private

  integer, parameter :: equN = 18  ! this number should be multiple of six
  integer, parameter :: varN = 18  !
  real(kind=DBL), parameter :: tol = epsilon(1.0_DBL) * 10._DBL
  logical :: isJacobianPresented = .false.
  real(kind=WP), dimension(equN) :: equX0
  real(kind=WP), dimension(equN) :: equX
  public :: equSolve, equX, equX0, equN, varN

contains

  subroutine equSolve ( method )
    integer, intent(in), optional :: method
    real(kind=DBL), dimension(equN) :: fvec
    real(kind=DBL), dimension(equN, varN) :: fjac
    integer :: info, lwa
    real(kind=DBL), dimension(:), allocatable :: wa
    external :: func_j, func_noj
    lwa = equN * (equN + 13) / 2 + 1
    allocate( wa(lwa) )
    equX = equX0

    if ( present(method) ) then
       select case (method)
       case (1) ! this is Newton method
       case (2) ! this Simplified Newton method
       case default
          write(*,*) "This method has not implemented."
       end select
    else        ! this is MINPACK
       if ( isJacobianPresented ) then
          call hybrj1 ( func_j, equN, equX, fvec, fjac, varN, tol, info, wa, lwa )
       else
          call hybrd1 ( func_noj, equN, equX, fvec, tol, info, wa, lwa )
       end if
    end if

  end subroutine equSolve

  subroutine func_noj ( N, X, FVEC, IFLAG )
    integer, intent(in) :: N
    integer, intent(inout) :: IFLAG
    real(kind=DBL), dimension(N), intent(in) :: X
    real(kind=DBL), dimension(N), intent(out) :: FVEC

    FVEC = f ( X )
  end subroutine func_noj

  subroutine func_j ( N, X, FVEC, FJAC, LDFJAC, IFLAG )
    integer, intent(in) :: N, LDFJAC
    integer, intent(inout) :: IFLAG
    real(kind=DBL), dimension(N), intent(in) :: X
    real(kind=DBL), dimension(N), intent(inout) :: FVEC
    real(kind=DBL), dimension(LDFJAC, N), intent(inout) :: FJAC

    if ( IFLAG == 1 ) then
       FVEC = f ( X )
    else if ( IFLAG == 2 ) then
       FJAC = jacobian ( X )
    end if
  end subroutine func_j

  function f ( x )
    real(kind=WP), dimension(equN) :: f
    real(kind=WP), intent(in), dimension(equN) :: x
    integer :: i, j, N
    real(kind=WP) :: M1 = 7.0

    N = equN

    f = 0.0

    do i = 0, N/2-1
       do j = 1, N
          if ( j - 2*i > 0 ) then
             f(i+1) = f(i+1) + x(j) * x(j-2*i)
          end if
       end do
    end do

    f(1) = f(1) - 2.0

    f(N/2+1) = sum( x ) - 2.0

    do i = 1, N/3-1
       do j = 1, N
          f(N/2+1+i) = f(N/2+1+i) + (-1)**(j+1) * x(j) * (j-1)**i
       end do
    end do

    do i = 1, N/6
       do j = 1, N
          f(5*N/6+i) = f(5*N/6+i) + (j-1)**(2*i-1) * x(j)
       end do
       f(5*N/6+i) = f(5*N/6+i) - 2.0 * M1**(2*i-1)
    end do

  end function f

  function jacobian ( x ) result ( jac )
    real(kind=WP), dimension(equN, equN) :: jac
    real(kind=WP), dimension(equN), intent(in) :: x
    jac = 0.
  end function jacobian

end module equSolver
