module fileUnits
    use iso_fortran_env, only: int8

    implicit none

    private

    public :: &
    & flash_input_unit,&

    & intrcn32_unit, qPar150_unit, rPar150_unit, &
    & intrcn_unit, &
    & intrcnas_unit, &
    & parvolas_unit, & 
    & pareneas_unit, &
    & gruposram_unit, &
    & intrcngcalpha_unit, &
    & gruposramgc_unit, &
    & intrcngckapa_unit, &

    & output_unit, &
    & iout_unit

    ! These are module variables that can be used by any program unit 
    ! that uses this module
    integer(kind=int8) :: &
    & flash_input_unit, &

    & intrcn32_unit, &
    & qPar150_unit, &
    & rPar150_unit, &
    & intrcn_unit, &
    & intrcnas_unit, &
    & parvolas_unit, & 
    & pareneas_unit, &
    & gruposram_unit, &
    & intrcngcalpha_unit, &
    & gruposramgc_unit, &
    & intrcngckapa_unit, &

    & output_unit, &
    & iout_unit = 1


  end module fileUnits