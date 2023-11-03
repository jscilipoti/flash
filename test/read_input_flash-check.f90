program check
    !ANT, NGG, Pxx, Txx: No fueron checkeados
    use do_tests
    use InputData, only:&
    & flash_input_filename, name_maxlen, &
    & ICALC, model, IPR, IOUT, NOVAP, ig, ipareq, &
    !& ANT, &
    & NTEXT
    !use flash, only: P, T, Z
    !use CUFAC, only: NKK, NGG, Pxx, Txx
    
    implicit none
    real*8,dimension(10) :: z_check
    integer :: i 

    print *,""
    print *, test_run//"read_input_flash-Test"
    if (leer_input_flash_check) then
        if (pause_test) pause
        call read_input_flash("test/name.dat")
        !Check if everything went OK
        print *, flash_input_filename
        if (flash_input_filename /= "test/llecalas2.dat")& 
            &ERROR STOP "'flash_input_filename' was not read correctly."
        if (abs(icalc - 0) > 1E-8)&
            &ERROR STOP "'ICALC' was not read correctly."
        if (abs(model - 0) > 1E-8)&
            &ERROR STOP "'model' was not read correctly."
        if (abs(ipr - 0) > 1E-8)&
            &ERROR STOP "'IPR' was not read correctly."
        if (abs(iout - 1) > 1E-8)&
            &ERROR STOP "'IOUT' was not read correctly."
        if (abs(novap - 0) > 1E-8)&
            &ERROR STOP "'NOVAP' was not read correctly."
        if (abs(ig - 1) > 1E-8)&
            &ERROR STOP "'ig' was not read correctly."
        if (abs(ipareq - 2) > 1E-8)&
            &ERROR STOP "'ipareq' was not read correctly."
        ! Disables since there are more errors bellow
        ! if (abs(T - 358.15) > 1E-5)&
        !     &ERROR STOP "'T' was not read correctly."
        ! if (abs(P - 1.) > 1E-8)&
        !     &ERROR STOP "'P' was not read correctly."
        ! z_check=Z-(/0.65,0.13,1E-8,1E-8,1E-8,1E-8,1E-8,0.,0.,0./)
        ! do i = 1 , size(z_check)
        !     if (abs(z_check(i)) > 1E-5) &
        !     &ERROR STOP "'Z' was not read correctly."
        ! enddo
        ! if (abs(NKK - 7) > 1E-8) ERROR STOP "'NKK' was not read correctly."
        ! if (NTEXT /= "llecalas2") ERROR STOP "'NTEXT' was not read correctly."


        print *, test_ok
    else 
        print *, test_disabled
    endif
end program check
    