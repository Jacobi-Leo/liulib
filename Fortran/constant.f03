! Created by Jacob Zeyu LIU (liuzeyu271828@gmail.com)

module constant
  implicit none

  integer, parameter :: SGL = selected_real_kind(p=6), DBL = selected_real_kind(p=13)
  integer, parameter :: WP = DBL
  real(kind=WP), parameter :: E = exp(1.0d0), PI = 4.0*atan(1.0d0)

end module constant

