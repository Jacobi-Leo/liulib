module linearAlgebra
  use utils
  use denseMatrix
  implicit none
  public

contains
  function linspace ( a1, an, n) result (l)
    integer, intent(in) :: n
    real(kind=WP), dimension(n) :: l
    real(kind=WP), intent(in) :: a1, an

    real(kind=WP) :: step
    integer :: i

    if ( a1 < an .or. n < 2 ) then
       stop 'Invalid argument'
    end if

    step = (an - a1) / (n-1)

    do i = 1, n
       l(i) = a1 + (i-1) * step
    end do
  end function linspace

  function logspace ( a1, an, n, base ) result (l)
    integer, intent(in) :: n
    real(kind=WP), intent(in) :: a1, an
    real(kind=WP), intent(in), optional :: base
    real(kind=WP), dimension(n) :: l, ll
    real(kind=WP) :: b

    if ( .not. present( base ) ) then
       b = 10.0_WP
    else
       b = base
    end if

    ll = linspace(a1, an, n)
    l = base**ll
  end function logspace

end module linearAlgebra


