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
program main
    implicit none
    
    call llecalas()

endprogram main