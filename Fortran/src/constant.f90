! Created by Jacob Zeyu LIU (liuzeyu271828@gmail.com)
!
! Must be compiled with -fdefault-real-8 with gfortran or
! similar parameter in ifort.

module constant
  implicit none

  integer, parameter :: SGL = selected_real_kind(p=6), DBL = selected_real_kind(p=13)
  integer, parameter :: WP = DBL
  real(kind=WP), parameter :: E = exp( 1.0_WP ), PI = 4.0_WP * atan( 1.0_WP ), EPS = epsilon( 1.0_WP )

end module constant

