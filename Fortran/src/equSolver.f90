module equations
  use constant
  implicit none
  save
  private

  integer :: equN = 18  ! this number should be multiple of six
  public :: equN, f

contains
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
end module equations

module equSolver
  use constant
  use equations
  use denseMatrix
  implicit none
  private
end module equSolver

