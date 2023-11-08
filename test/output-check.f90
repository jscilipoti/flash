program check
    use do_tests
    use stdlib_ansi, only : & 
    & fg_color_green, fg_color_red, fg_color_yellow, & 
    & style_bold, style_reset, operator(//), operator(+)
    use iso_fortran_env, only: int32
    
    implicit none
    
    character(len=86),dimension(70) :: &
    &lleasoccuzada_new,lleasoccuzada_old
    character(len=70),dimension(6) :: &
    &output_new,output_old
    
    integer(kind=int32) :: i, j

    print *,""
    print *, test_run//"output-check"
    
    if (output_check) then
        if (pause_test) pause

        !call llecalas 
        
        call open_textfile("lleasoccuzada.OUT",&
            &lleasoccuzada_new,70,86)
        !print "(A)",lleasoccuzada_new
        call open_textfile("lleasoccuzada_original.OUT",&
            &lleasoccuzada_old,70,86)
        !print "(A)",lleasoccuzada_old
        call open_textfile("output.OUT",&
            &output_new,6,70)
        !print "(A)",output_new
        call open_textfile("output_original.OUT",&
            &output_old,6,70)
        !print "(A)",output_old
        
        !Check if everything went OK
        do i = 1, 70
            if (trim(lleasoccuzada_new(i)) /= trim(lleasoccuzada_old(i))) then
                print *, fg_color_red + style_bold // test_ error // style_reset
                print *,&
                    trim(lleasoccuzada_new(i)) /= trim(lleasoccuzada_old(i))
                print '(3A)',&
                    &trim(lleasoccuzada_new(i)),"//",trim(lleasoccuzada_old(i))
                ERROR STOP ""
            end if
        end do
        do i = 1, 6
            if (trim(output_new(i)) /= trim(output_old(i))) then
                print *, fg_color_red + style_bold // test_ error // style_reset
                print *,&
                    &trim(output_new(i)) /= trim(output_old(i))
                print '(3A)',&
                    &trim(output_new(i)),"//",trim(output_old(i))
                ERROR STOP ""
            end if
        end do
        
        ! if (abs(0 - 0) > 1E-8)&
        !     &ERROR STOP ""
            
        print *, fg_color_green + style_bold // test_ok // style_reset
    else 
        print *, fg_color_yellow + style_bold // test_disabled // style_reset
    endif

!print *, "Put some tests in here!"

end program check
