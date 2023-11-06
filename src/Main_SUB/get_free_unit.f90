! A subroutine that returns a free unit number to open a file.
subroutine get_free_unit(unit)
use iso_fortran_env, only: int8, int16
implicit none

integer(kind=int8) :: unit
integer(kind=int16) :: iostat
logical :: opened
do unit = 10, 126 ! Loop over possible unit numbers (126 is int8 max number)
  inquire (unit=unit, opened=opened, iostat=iostat) ! Check if the unit is open
  if (iostat /= 0) cycle ! Skip if there is an error
  if (.not. opened) exit ! Exit the loop if the unit is free
end do

! Return the free unit number
end subroutine