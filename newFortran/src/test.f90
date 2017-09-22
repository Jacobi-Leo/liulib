program test
  use constant
  use denseMatrix
  implicit none

  real(kind=WP) :: a(3,3), b(3), b2(3,1)

  a = reshape([2.,0.,0.,2.,3.,0.,1.,1.,1.], [3,3])
  b = [5.,4.,1.]
  b2 = reshape([5.,4.,1.], [3,1])

  write(*,*) solve(a, b), new_line('*')
  write(*,*) solve(a, b2), new_line('*')
  call writematrix( inv(a) )
  call writeMatrix(solve(a, eye(3)))
end program test
