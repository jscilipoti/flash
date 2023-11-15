program check
    use do_tests

    use inputData, only: input_filename
    use iso_fortran_env, only: int32
    
    use stdlib_ansi, only : & 
    & fg_color_green, fg_color_red, fg_color_yellow, fg_color_blue, & 
    & style_bold, style_reset, operator(//), operator(+)
    
    implicit none
    
    integer(kind=int32), parameter :: lleasoccuzada_linelen = 86
    integer(kind=int32), parameter :: lleasoccuzada_maxlines = 93
    integer(kind=int32), parameter :: output_linelen = 70
    integer(kind=int32), parameter :: output_maxlines = 4
    
    logical :: is_OK = .true.
    integer(kind=int32) :: error_count = 0
    
    character(len=lleasoccuzada_linelen), dimension(lleasoccuzada_maxlines) :: &
    &lleasoccuzada_new,lleasoccuzada_old
    character(len=output_linelen), dimension(output_maxlines) :: &
    &output_new,output_old
    character(len=*), parameter :: PATH = "test/"
    character(len=*), parameter :: flashtype = "ALV2_flash_"
    
    integer(kind=int32) :: i, j

    input_filename = PATH//flashtype//"name.dat"
    

    print *,""
    print *, fg_color_blue + style_bold // test_run // style_reset //"ALV2_flash_check"
    
    if (.true.) then
        if (pause_test) pause

        call llecalas 
        
        call open_textfile("lleasoccuzada.OUT", &
            & lleasoccuzada_new, lleasoccuzada_maxlines, lleasoccuzada_linelen)
        !print "(A)",lleasoccuzada_new
        call open_textfile(PATH//flashtype//"lleasoccuzada.OUT",&
            &lleasoccuzada_old, lleasoccuzada_maxlines, lleasoccuzada_linelen)
        !print "(A)",lleasoccuzada_old
        call open_textfile("output.OUT",&
            &output_new, output_maxlines, output_linelen)
        !print "(A)",output_new
        call open_textfile(PATH//flashtype//"output.OUT",&
            &output_old, output_maxlines, output_linelen)
        !print "(A)",output_old
        
        !Check if everything went OK
        do i = 1, lleasoccuzada_maxlines
            if (trim(lleasoccuzada_new(i)) /= trim(lleasoccuzada_old(i))) then
                print *, fg_color_red + style_bold // test_error // style_reset
                print *, "At line: ",i, " of file: lleasoccuzada.OUT"
                print *,&
                    trim(lleasoccuzada_new(i)) /= trim(lleasoccuzada_old(i))
                print '(3A)',&
                    &trim(lleasoccuzada_new(i)),"//",trim(lleasoccuzada_old(i))
                is_OK = .false.
                error_count = error_count + 1
                !ERROR STOP ""
            end if
        end do
        do i = 1, output_maxlines
            if (trim(output_new(i)) /= trim(output_old(i))) then
                print *, fg_color_red + style_bold // test_error // style_reset
                print *, "At line: ",i, " of file: output.OUT"
                print *,&
                    &trim(output_new(i)) /= trim(output_old(i))
                print '(3A)',&
                    &trim(output_new(i)),"//",trim(output_old(i))
                is_OK = .false.
                error_count = error_count + 1
                !ERROR STOP ""
            end if
        end do
        
        ! if (abs(0 - 0) > 1E-8)&
        !     &ERROR STOP ""
                if (error_count >= 10) then
            ERROR STOP "TOO MANY ERRORS"
        end if
        if (error_count < 9) then
            ERROR STOP "LESS ERRORS!!!"
        end if
        if (is_OK) then
            print *, fg_color_green + style_bold // test_ok // style_reset
        end if    
    else 
        print *, fg_color_yellow + style_bold // test_disabled // style_reset
    endif

!print *, "Put some tests in here!"

end program check
