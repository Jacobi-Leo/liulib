program test
  use utils
  use linearAlgebra
  implicit none

  real(kind=WP) :: a(4,4), b(4), b2(4,1)

  a = reshape([2.0, 0.0, 0.0, 0.0, 2.0 , 3.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0], [4, 4])
  b = [5.0, 4.0, 2.0, 2.0]
  b2 = reshape([5.0, 4.0, 2.0, 2.0], [4,1])

  write(*,*) solve(a, b), new_line('*')
  write(*,*) solve(a, b2), new_line('*')
  call writematrix( inv(a) )
  call writeMatrix(solve(a, eye(4)))
end program test
