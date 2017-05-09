program test
  use constant
  use denseMatrix
  implicit none

  integer, parameter :: m = 5, n = 3
  real(kind=WP), dimension(m, n) :: A, B
  type(Matrix) :: Mat, Mat2, Mat3, Mat4
  integer :: i, j

  forall ( i=1:m, j=1:n )
     A(i,j) = real(i*10+j, WP)
  end forall
  B = A + 1.0
  Mat = A
  Mat2 = A * 2.0
  Mat3 = Mat + Mat2
  Mat4 = Mat3
  call Mat4%addMatrix(Mat)


  open(unit=11, file="constantResult.txt", action="write", status="replace")
  write(11, *) SGL, DBL, WP, E, PI
  close(11)

  open(unit=11, file="denseMatrixResult.txt", action="write", status="replace")
  write(11, *) A
  write(11, *)
  write(11, *) Mat4%comp
  close(11)

end program

