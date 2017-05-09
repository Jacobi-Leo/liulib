program test
  use constant
  use denseMatrix
  implicit none

  integer, parameter :: m = 4, n = 3
  real(kind=WP), dimension(m, n) :: A, B
  type(Matrix) :: Mat
  integer :: i, j

  forall ( i=1:m, j=1:n )
     A(i,j) = real(i*10+j, WP)
  end forall
  Mat = A



  open(unit=11, file="constantResult.txt", action="write", status="replace")
  write(11, *) SGL, DBL, WP, E, PI
  close(11)

  open(unit=11, file="denseMatrixResult.txt", action="write", status="replace")
  write(11, *) A
  write(11, *)
  write(11, *) Mat
  close(11)

end program

