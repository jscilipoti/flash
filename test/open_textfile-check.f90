program check
    use do_tests
    use stdlib_ansi, only : & 
    & fg_color_green, fg_color_red, fg_color_yellow, fg_color_blue, & 
    & style_bold, style_reset, operator(//), operator(+)
    use iso_fortran_env, only: int32
    
    implicit none
    
    character(len=8) :: filename = "name.dat"
    character(len=16) :: flash_filename = "llecalas2.dat"
    character(len=36),dimension(2) :: array = (/"",""/)
    integer(kind=int32) :: doLeerBases_tmp
    
    print *,""
    print *, fg_color_blue + style_bold // test_run // style_reset //"open_file-check"
    if (pause_test) pause
    call open_textfile(filename,array,2,36)
    if (open_file_check) then
        if (filename/="name.dat") then
            print *, fg_color_red + style_bold // test_error // style_reset
            ERROR STOP "Filename has been changed."
        end if
        if (size(array)/=2) then
            print *, fg_color_red + style_bold // test_error // style_reset
            ERROR STOP "Max. lines has been changed."
        end if
        if (len(array(1))/=36) then
            print *, fg_color_red + style_bold // test_error // style_reset
            ERROR STOP "Max. char number has been changed."
        end if
        read(array(1),'(I1)') doLeerBases_tmp
        if (doLeerBases_tmp*0/=0) then
            print *, fg_color_red + style_bold // test_error // style_reset
            ERROR STOP "First line read was not a number"
        end if
        if (len_trim(array(2))/=15) then
            print *, fg_color_red + style_bold // test_error // style_reset
            ERROR STOP "Second line has been changed(15)."
        end if
        if (trim(array(2))/='"llecalas2.dat"') then
            print *, fg_color_red + style_bold // test_error // style_reset
            ERROR STOP "Second line has been changed."
        end if
    print *, fg_color_green + style_bold // test_ok // style_reset
    else
        print *, fg_color_yellow + style_bold // test_disabled // style_reset
    end if
end program check