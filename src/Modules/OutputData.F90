! MODULE: outputData
! This module store some parameters and variables needed to output results
! calc.
module outputData
    
    !use iso_fortran_env
    
    implicit none

    private 


    public :: &
        & output_console, &
        & output_lleasoccuzada, &
        & output_outputfile

    ! Allow 'write' or 'print' to be shown in console
    logical :: output_console = .true.
    logical :: output_lleasoccuzada = .true.
    logical :: output_outputfile = .true.



end module outputData