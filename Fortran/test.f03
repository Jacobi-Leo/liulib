program test
  use constant
  implicit none

  open(unit=11, file="constantResult.txt", action="write", status="replace")
  write(11, *) SGL, DBL, WP, E, PI
  close(11)

end program test

