program test
  use constant
  use denseMatrix
  implicit none

  integer, parameter :: m = 5, n = 3
  real(kind=WP), dimension(:, :), allocatable :: A, B
  real(kind=WP), dimension(:), allocatable :: bb
  type(Matrix) :: Mat, Mat2, Mat3, Mat4, Mat5, Mat6, Vec, Vec2
  integer :: i, j, ii, jj
  namelist / mylist2 / i, j

  allocate( A(m,n), B(m,n), bb(n) )

  write(*,*) Mat%getNrow()

  forall ( i=1:m, j=1:n )
     A(i,j) = real(i*10+j, WP)
  end forall
  B = A + 1.0
  do i = 1, size(bb)
     bb(i) = i
  end do
  Vec = bb
  Mat = A
  ! Mat2 = A * 2.0
  Mat2 = 2 * A
  Mat3 = Mat + Mat2
  Mat4 = Mat3
  call Mat4%addMatrix(Mat)
  Mat5 = 5 * Mat


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

  write(11, *) 'This is transpose and matrix multiplication: Mat5**T * Mat'
  call Mat6%writeToFile(11) !! this result has been verified by Matlab

  write(11, *) 'This matrix times a column vector: Mat5**T * Vec'
  Vec2 = Mat5 * Vec
  call Vec2%writeToFile(11) !! this result has been verified by Matlab

  call Mat5%T()
  call Mat5%pushColumn( bb )
  write(11, *) 'This is the extended Mat5: '
  write(11, *) Mat5%getNrow(), Mat5%getNcolumn()
  call Mat5%writeToFile(11)

  !!============== test attribute predicates ========================
  write(*,*) 'Is Mat6 diagonal? ', Mat6.isDiagonal()
  write(*,*) 'Is Mat6 allocated? ', Mat6.isAllocated()
  write(*,*) 'Is Mat6 bidiagonal? ', Mat6.isBidiagonal()
  write(*,*) 'Is Mat6 Hessenberg? ', Mat6.isHessenberg()
  write(*,*) 'Is Mat6 symmetric? ', Mat6.isSymmetric()
  write(*,*) 'Is Mat6 tridiagonal? ', Mat6.isTridiagonal()
  ii = Mat6.getNrow()
  jj = Mat6.getNcolumn()
  do i = 1, ii
     write(*,*) (Mat6.getEllement(i,j), j=1,jj)
  end do
  call Mat6.printSpecialAttributes()

  deallocate( A, B, bb )
  allocate( A(3, 3) )
  A = reshape([3,1,1,0,4,0,1,2,5], [3,3])

  write(*,*) 'array A now is: '
  write(*,*) A

  Mat = A

  write(11, *) 'Mat is now: '
  call Mat.writeToFile(11)

  allocate( bb(3) )
  bb = [6, 15, 16]
  Vec = bb
  write(11, *) 'Vec is now:'
  call Vec.writeToFile(11)

  call Mat.solve(Vec)
  write(11, *) 'Solution to equation Mat * x = Vec is  x = '
  call Vec.writeToFile(11)

  Mat = A
  Vec = bb
  deallocate( bb )
  allocate( bb(3) )
  bb = [4., 4., 5.]
  call Mat.pushRow( bb )

  deallocate( bb )
  allocate( bb(1) )
  bb = [27.]
  call Vec.pushRow( bb )

  write(11, *) 'Now   Mat = '
  call Mat.writeToFile(11)
  write(11, *) 'Now   Vec = '
  call Vec.writeToFile(11)
  write(11, *) 'Now solution is   x = '
  ! call Mat.solve(Vec)
  ! call Vec.writeToFile(11)


  close(11)

end program
