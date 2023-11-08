program check
    use do_tests
    use stdlib_ansi, only : & 
    & fg_color_green, fg_color_red, fg_color_yellow, & 
    & style_bold, style_reset, operator(//), operator(+)


    implicit none


    print *,""
    print *, test_run//"check"
    
    if (.true.) then
        if (pause_test) pause
 
        !Check if everything went OK
        if ("string" /= "string")& 
            &ERROR STOP ""
        if (abs(0 - 0) > 1E-8)&
            &ERROR STOP ""
            
        print *, fg_color_green + style_bold // test_ok // style_reset

    else 
        print *, fg_color_yellow + style_bold // test_disabled // style_reset
    endif

!print *, "Put some tests in here!"

end program check
