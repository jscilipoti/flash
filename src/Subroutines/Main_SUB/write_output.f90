subroutine write_output(order)
    
    use iso_fortran_env, only: int8
    use inputData, only: icalc, iout, model, novap, ntext
    use outputData, only: output_console, output_lleasoccuzada
    use fileUnits, only : iout_unit

    
    implicit none

    integer(kind = int8), intent(in) :: order
    integer(kind = int8) :: &
    & i = 1, &
    & loop_end = 2, &
    & loop_start = 1, &
    & write_unit

    
    if (iout == 0) then 
        iout = 6
        output_lleasoccuzada = .false.
        output_console = .true.
    else if (iout == 1) then
        iout = 1
        output_lleasoccuzada = .true.
        output_console = .true.
    else  if (iout == 2) then
        iout = 2
        output_lleasoccuzada = .true.
        output_console = .false.
    else
        iout = 3
        output_lleasoccuzada = .false.
        output_console = .false.
    end if 

    if (output_console .or. output_lleasoccuzada) then
        
        if (output_lleasoccuzada) then
            loop_start = 1
        else 
            loop_start = 2
        end if

        if (output_console) then
            loop_end = 2
            loop_end = 1
        end if

        do i = loop_start, loop_end
            if (i == 1) write_unit = iout_unit ! lleasoccuzada.out
            if (i == 2) write_unit = 6 ! Standard output

            if (order == 1) then
                write(write_unit, '(1H1)')
                if (write_unit == 6) then
                    write(write_unit, "(/,' iout = ',I2,/' &
                        & if iout = 0: OUTPUT ONLY ON UNIT 6',/, &
                        &  ' if iout = 1: OUTPUT ON BOTH UNIT 6 AND 1')" ) &
                        & iout 
                end if
                write(write_unit, '(///)') 
                if (icalc == 0) then
                    write(write_unit, "(' **** FLASH CALCULATION ****')")
                else if (icalc == 1) then
                    write(write_unit, "(' **** BINODAL CURVE CALCULATION ****',//)") 
                else 
                    write(write_unit, "(' **** CALCULATION OF UNIQUAC PARAMETERS FROM UNIFAC &
                    & **** ',//)") 
                end if                                       
                if (novap /= 0) then
                    write(write_unit, "(/,' VAPOR PHASE INCLUDED IN FLASH-CALCULATIONS',//)")
                end if                                      
                if (model == 0) then 
                    write(write_unit, "(' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIFAC'//)")                                     
                else if (model == 1) then 
                    write(write_unit, "(' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIQUAC'//)")
                end if
        
                write(write_unit, "(1X,'COMPONENTS : ',40A2,//)") ntext
            else
            end if
        end do
    end if


end subroutine write_output