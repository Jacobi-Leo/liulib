module utils
  use constant
  implicit none
  public

contains

  function outputName ( unit ) result ( name )
    integer, intent(in) :: unit
    character(len=15) :: name
    character(len=10) :: cTemp

    if ( unit > 999999999 ) then
       stop 'unit too large!'
    end if

    write(cTemp, '(i10)') unit
    name = 'data' // trim(adjustl(cTemp)) // '.out'
  end function outputName
end module utils

