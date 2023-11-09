integer function mainsgfunc (i, ipareq)

!-------------------------------------------------------------------------------
!       This function returns the main group number corresponding to the 
!       subgroup "i" of the "ipareq" parameter table:
!           1: liquid-liquid
!           2: liquid-vapor
!           3: infinite-dilution
!       The function calculates the record number to read from by adding i 
!       to the product of ipareq and 150, which is the number of records 
!       per parameter table. It assumes that unit 14 (gruposram.mds) 
!       has been opened and positioned correctly before calling it.
!       The inputs are "i", which is the record number within a parent request, 
!       and "ipareq", which is the parameter table number.
!       The output is mainsgfunc, which is the main group number.
!-------------------------------------------------------------------------------

    use iso_fortran_env, only: int8, int32

    implicit none

    integer(kind = int8), intent(in) :: i, ipareq
    integer(kind = int8) :: main_group
    integer(kind = int32) :: read_record

	read_record = i + (ipareq - 1) * 150
    
    read (14, '(i4)', rec = read_record) main_group
    
    mainsgfunc = main_group
    return

end function mainsgfunc