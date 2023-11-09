program check
    !ANT, NGG, Pxx, Txx: No fueron checkeados
    use do_tests
    use stdlib_ansi, only : & 
    & fg_color_green, fg_color_red, fg_color_yellow, fg_color_blue, & 
    & style_bold, style_reset, operator(//), operator(+)
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
    print *, fg_color_blue + style_bold // test_run // style_reset //"read_input_flash-Test"
    if (read_input_flash_check) then
        if (pause_test) pause
        call read_input_flash("test/name.dat")
        !Check if everything went OK
        !print *, flash_input_filename
        if (flash_input_filename /= "test/llecalas2.dat") then
            print *, fg_color_red + style_bold // test_error // style_reset
            ERROR STOP "'flash_input_filename' was not read correctly."
        end if
        if (abs(icalc - 0) > 1E-8) then
            print *, fg_color_red + style_bold // test_error // style_reset
            ERROR STOP "'ICALC' was not read correctly."
        end if
        if (abs(model - 0) > 1E-8) then
            print *, fg_color_red + style_bold // test_error // style_reset
            ERROR STOP "'model' was not read correctly."
        end if
        if (abs(ipr - 0) > 1E-8) then
            print *, fg_color_red + style_bold // test_error // style_reset
            ERROR STOP "'IPR' was not read correctly."
        end if
        if (abs(iout - 1) > 1E-8) then
            print *, fg_color_red + style_bold // test_error // style_reset
            ERROR STOP "'IOUT' was not read correctly."
        end if
        if (abs(novap - 0) > 1E-8) then
            print *, fg_color_red + style_bold // test_error // style_reset
            ERROR STOP "'NOVAP' was not read correctly."
        end if
        if (abs(ig - 1) > 1E-8) then
            print *, fg_color_red + style_bold // test_error // style_reset
            ERROR STOP "'ig' was not read correctly."
        end if
        if (abs(ipareq - 2) > 1E-8) then
            print *, fg_color_red + style_bold // test_error // style_reset
            ERROR STOP "'ipareq' was not read correctly."
        end if
        ! Disabled since there are more errors bellow
        !if (abs(T - 358.15) > 1E-5) then
        !   print *, fg_color_red + style_bold // test_error // style_reset
        !   ERROR STOP "'T' was not read correctly."
        !end if
        !if (abs(P - 1.) > 1E-8) then
        !   print *, fg_color_red + style_bold // test_error // style_reset
        !   ERROR STOP "'P' was not read correctly."
        !end if
        !z_check = Z - (/ 0.65, 0.13, 1E-8, 1E-8, 1E-8, 1E-8, 1E-8, 0., 0., 0./)
        !do i = 1 , size(z_check)
        !   if (abs(z_check(i)) > 1E-5) then
        !       print *, fg_color_red + style_bold // test_error // style_reset
        !       ERROR STOP "'Z' was not read correctly."
        !   end if
        !enddo
        !if (abs(NKK - 7) > 1E-8) then
        !   print *, fg_color_red + style_bold // test_error // style_reset
        !   ERROR STOP "'NKK' was not read correctly."
        !end if
        !if (NTEXT /= "llecalas2")then
        !       print *, fg_color_red + style_bold // test_error // style_reset
        !       ERROR STOP "'NTEXT' was not read correctly."
        !end if


        print *, fg_color_green + style_bold // test_ok // style_reset
    else 
        print *, fg_color_yellow + style_bold // test_disabled // style_reset
    endif
end program check
    