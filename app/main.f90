! Program: Main
!   * The 'implicit none' statement tells the compiler that no implicit typing 
!   is to be used in the program. This means that every variable must be 
!   explicitly declared with a type before it can be used.
!
!   * The 'call llecalas()' statement calls the subroutine `llecalas` from the 
!   main program. Some time before, 'llecalas' was the main program but now it
!   is rewritten as a subroutine in order to allow other programs to make use
!   of it.
!
!   * The 'endprogram main' statement marks the end of the main program.

!RECORDAR:
!Cada funcion/subrutina que lee de las bases de datos deben estar aisladas para 
!su posible modificacion
!Probar real128
program main
    use inputData, only: input_filename
    use outputData
    
    implicit none

    ! Allow 'write' or 'print' to be shown in console
    output_console = .true.
    input_filename = "name.dat"

    call llecalas()

    
endprogram main