program check
    use do_tests
    use stdlib_ansi, only : & 
    & fg_color_green, fg_color_red, fg_color_yellow, fg_color_blue, & 
    & style_bold, style_reset, operator(//), operator(+)
    
    implicit none


    print *,""
    print *, fg_color_blue + style_bold // test_run // style_reset //"open_database-check"
    
    if (open_database_check) then
        if (pause_test) pause
        call open_database(0)
        close(unit=13)
        close(unit=14)
        call open_database(1)
        close(unit=13)
        close(unit=14)
        call open_database(2)
        close(unit=13)
        close(unit=14)
        close(unit=15)
        close(unit=16)
        call open_database(3)
        close(unit=13)
        close(unit=14)
        close(unit=16)
        !Check if everything went OK
        ! if ("string" /= "string")& 
        !     &ERROR STOP ""
        ! if (abs(0 - 0) > 1E-8)&
        !     &ERROR STOP ""
            
        print *, fg_color_green + style_bold // test_ok // style_reset
    else 
        print *, fg_color_yellow + style_bold // test_disabled // style_reset
    endif

!print *, "Put some tests in here!"

end program check
