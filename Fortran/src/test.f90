module tester
  implicit none
contains

  subroutine tester_constant()
    use constant
    implicit none

    110 format ( 1X, A40, ES21.14) ! The most standard way of output double precision real
    120 format ( 1X, A40, I10.3)
    open(unit=11, file="constantResult.txt", action="write", status="replace")
    write(11, 120) "Single precision kind is ", SGL
    write(11, 120) "Double precision kind is ", DBL
    write(11, 120) "Working precision kind is ", WP
    write(11, 110) "Euler constant e is ", E
    !write(11, 110) "π is ", PI ! Unicode support is fragile.
    write(11, 110) "PI is ", PI
    close(11)
    write(*,*) "Please check file constantResult.txt..."
  end subroutine tester_constant

  subroutine tester_denseMatrix()
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

    write(*,*) "The current row number of Matrix A is:"
    write(11,*) Mat%getNrow()

    forall ( i=1:m, j=1:n )
       A(i,j) = real(i*10+j, WP)
    end forall
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

    open(unit=11, file="denseMatrixResult.txt", action="write", status="replace")
    write(11, *) 'This is array A: '
    write(11, *) A
    B = Mat
    write(11, *) 'This is array B: '
    write(11, *) B
    write(11, *) 'This is matrix A: '
    call Mat%writeToFile(11)
    write(11, *) 'This is matrix A*5: '
    write(11, *) Mat%getNrow(), Mat5%getNcolumn()
    call Mat5%writeToFile(11)
    write(11, *) 'This is vector b'
    call Vec%writeToFile(11)

    Mat5 = Mat5%trans()
    Mat6 = Mat5 * Mat

    write(11, *) 'This is transpose and matrix multiplication: Mat5**T * Mat'
    call Mat6%writeToFile(11) !! this result has been verified by Matlab

    write(11, *) 'This matrix times a column vector: Mat5**T * Vec'
    Vec2 = Mat5 * Vec
    call Vec2%writeToFile(11) !! this result has been verified by Matlab

    Mat5 = Mat5%trans()
    call Mat5%pushColumn( bb )
    write(11, *) 'This is the extended Mat5: '
    write(11, *) Mat5%getNrow(), Mat5%getNcolumn()
    call Mat5%writeToFile(11)

    !!!============== test attribute predicates ========================
    !write(*,*) 'Is Mat6 diagonal? ', Mat6%isDiagonal()
    !write(*,*) 'Is Mat6 allocated? ', Mat6%isAllocated()
    !write(*,*) 'Is Mat6 bidiagonal? ', Mat6%isBidiagonal()
    !write(*,*) 'Is Mat6 Hessenberg? ', Mat6%isHessenberg()
    !write(*,*) 'Is Mat6 symmetric? ', Mat6%isSymmetric()
    !write(*,*) 'Is Mat6 tridiagonal? ', Mat6%isTridiagonal()
    !ii = Mat6%getNrow()
    !jj = Mat6%getNcolumn()
    !do i = 1, ii
    !   write(*,*) (Mat6%getEllement(i,j), j=1,jj)
    !end do
    !call Mat6%printSpecialAttributes()
    !!!============== finished =========================================

    write(11, *) "Now, reset all variables... "
    deallocate( A, B, bb )
    allocate( A(3, 3) )
    A = reshape([3,1,1,0,4,0,1,2,5], [3,3])

    write(11,*) 'array A now is: '
    write(11,*) A

    Mat = A

    write(11, *) 'Mat A is now: '
    call Mat%writeToFile(11)

    allocate( bb(3) )
    bb = [6, 15, 16]
    Vec = bb
    write(11, *) 'Vec is now:'
    call Vec%writeToFile(11)

    Vec2 = Mat%solve(Vec)
    write(11, *) 'Solution to equation Mat * x = Vec is  x = '
    call Vec2%writeToFile(11)

    deallocate( bb )
    allocate( bb(3) )
    bb = [4., 4., 5.]
    call Mat%pushRow( bb )

    deallocate( bb )
    allocate( bb(1) )
    bb = [27.]
    call Vec%pushRow( bb )

    write(11, *) 'Now   Mat = '
    call Mat%writeToFile(11)
    write(11, *) 'Now   Vec = '
    call Vec%writeToFile(11)
    write(11, *) 'Now solution is   x = '
    Vec2 = Mat%solve(Vec)
    call Vec2%writeToFile(11)

    call Mat%popRow()
    call Vec%popRow()

    write(11,*) 'Now   Mat = '
    call Mat%writeToFile(11)

    Mat2 = Mat%inv()
    write(11,*) 'The inverse matrix of  Mat  is:'
    call Mat2%writeToFile(11)

    Mat3 = Mat * Mat2
    write(11,*) 'Mat * Mat^-1 ='
    call Mat3%writeToFile(11)

    Mat3 = Mat2 * Mat
    write(11,*) 'Mat^-1 * Mat = '
    call Mat3%writeToFile(11)

    close(11)
    write(*,*) "Please check file denseMatrixResult.txt..."
  end subroutine tester_denseMatrix

  subroutine tester_equSolver
    use equSolver
    implicit none

    integer :: i
    i = 1
    call equSolve( i )

    open(unit=11, file="equSolverResult.txt", action="write", status="replace")
    do i = 1, equN
       write(11,*) equX(i)
    end do

    close(11)
    write(*,*) 'Please check file equSolverResult.txt...'
  end subroutine tester_equSolver

end module tester

program test
  use tester
  implicit none

  integer :: comp

  write(*,*) "This is the list of my components:"
  write(*,*) "1. constant"
  write(*,*) "2. denseMatrix"
  write(*,*) "3. equSolver"
  write(*,*) "which module should I test?"
  read(*,*) comp

  select case ( comp )
  case (1)
     call tester_constant ()
  case (2)
     call tester_denseMatrix ()
  case (3)
     call tester_equSolver ()
  case default
     write(*,'(1X, A33, I2, A33)') &
          &"Wrong input! You want to test the", comp, &
          &"th module, but I don't have that!"
  end select

end program test
