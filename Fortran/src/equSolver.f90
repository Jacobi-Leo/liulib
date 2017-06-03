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
    lwa = equN * (equN + 13) / 2 + 1
    allocate( wa(lwa) )
    equX = equX0

    if ( present(method) ) then
       select case (method)
       case (1) ! this is BFS method
          call BFS ()
       case (2) ! this
          write(*,*) "This method has not implemented."
       case default
          write(*,*) "This method is not valid!"
       end select
    else        ! this is MINPACK
       write(*,*) "There are known bugs, choose another method."
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

  subroutine BFS ()
    implicit none
    integer :: i, j, n=equN itmax = 200
    real(kind=WP) :: eps, reps
    x0 = equX0

    if ( WP == DBL ) then
       eps = 5d-5
       reps = 2d4
    else if ( WP == SGL ) then
       eps = 5d-2
       reps = 2d1
    else
       eps = epsilon(1.0_wp) * 10000.0_wp
       reps = 1.0_wp / eps
    end if

    do j = 1, equN
       x1 = x0
       x1(j) = x0(j) + eps
       df(:,j) = ( f(x1) - f(x0) ) * reps
    end do
    jac = df
    H0 = jac%inv()

    do i = 1, itmax
       f0 = f( x0 )
       tmp = H0 * (-1.0)
       tmp = tmp * f0
       bridge = tmp%trans()
       xtmp = bridge(1,:)
       x1 =  xtmp + x0
       dx = x1 - x0
       f1 = f( x1 )
       y = f1 - f0
       tmp = y%trans() * ( H0 * y )
       t1 = tmp%getEllement(1,1)
       tmp = dx
       tmp = tmp%trans() * y
       t2 = tmp%getEllement(1,1)
       u = 1. + t1 / t2
       tmp = dx
       m1 = tmp * tmp%trans() * u
       m2 = tmp * y%trans() * H0
       m3 = H0 * y * tmp%trans()
       tmp = tmp%trans() * y
       ttt = 1.0 / tmp%getEllement(1,1)
       H1 = ( m1 + ( (-1.0) * ( m2 + ( (-1.0) * m3 ) ) ) ) * ttt

       x0 = x1
       H0 = H1

       dx2 = sqrt( sum( dx * dx ) )

       if ( dx2 < tol ) exit

    end do

    equX = x0

  end subroutine BFS

end module equSolver
