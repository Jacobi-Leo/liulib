program test
  use constant
  use denseMatrix
  implicit none

  integer, parameter :: m = 5, n = 3
  real(kind=WP), dimension(m, n) :: A, B
  real(kind=WP), dimension(m) :: bb
  type(Matrix) :: Mat, Mat2, Mat3, Mat4, Mat5, Mat6, Vec
  integer :: i, j
  namelist / mylist2 / i, j

  write(*,*) Mat%getNrow()

  forall ( i=1:m, j=1:n )
     A(i,j) = real(i*10+j, WP)
  end forall
  B = A + 1.0
  bb = 1.0
  Vec = bb
  Mat = A
  ! Mat2 = A * 2.0
  Mat2 = 2 * A
  Mat3 = Mat + Mat2
  Mat4 = Mat3
  call Mat4%addMatrix(Mat)
  Mat5 = 5.0 * Mat


  open(unit=11, file="constantResult.txt", action="write", status="replace")
  write(11, *) SGL, DBL, WP, E, PI
  close(11)

  open(unit=11, file="denseMatrixResult.txt", action="write", status="replace")
  write(11, *) 'This is array A: '
  write(11, *) A
  write(11, *) 'This is matrix A: '
  call Mat%writeToFile(11)
  write(11, *) 'This is matrix Mat5: '
  write(11, *) Mat%getNrow(), Mat5%getNcolumn()
  call Mat5%writeToFile(11)
  write(11, *) 'This is vector b'
  call Vec%writeToFile(11)

  call Mat5%T()
  Mat6 = Mat5 * Mat

  write(11, *) 'This is transpose and matrix multiplication: '
  call Mat6%writeToFile(11)

  write(*,*) Mat6%isDiagonal()

  close(11)

end program
