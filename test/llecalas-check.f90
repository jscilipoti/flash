! BLoque raro de interpretar
    program check
    

    use iso_fortran_env, only: int8, real64
    use stdlib_ansi, only : & 
    & fg_color_green, fg_color_red, fg_color_yellow, & 
    & style_bold, style_reset, operator(//), operator(+)
    
    implicit none

    integer(kind=int8) ::  &
    & NF = 0, &
    & i = 0, &
    & j = 0

    integer(kind=int8), parameter :: &
    & N = 2

    real(kind=real64), dimension(N) :: &
    & YVAL = (/0.D0,1.D0/), &
    & YVAL_2 =  (/0.D0,1.D0/)

    !Codigo original
    100 NF = NF + 1                                                           
    104 do 105 i = 1, N                                                      
            if (YVAL(i) > 1.D0) then 
                !print *, "GOTO 106"  
                GOTO 106
            end if                                      
    105 CONTINUE
        !print *, "GOTO 109"                                                          
        GOTO 109                                                          
    106 do 107 i = 1, N                                                      
    107   YVAL(i) = YVAL(i) / 10.                                               
        !print *, "GOTO 104"  
        GOTO 104                                                          
    109 CONTINUE
    
    ! Print some values to compare later...
    ! print '(A5,I6)',"N " , N
    ! print '(A5,I6)',"NF " , NF
    ! print '(A5,I6)',"i " , i
    ! print '(A5,2F6.3)',"YVAL " , YVAL

    !Nuevo Código
    NF = 0
    i = 0

    1000 NF = NF + 1                                                           
    do i = 1, N                                                      
            if (YVAL_2(i) > 1.D0) then 
                do j = 1, N                                                      
                    YVAL_2(j) = YVAL_2(j) / 10.D0
                end do                                               
            end if
        end do                                      
    
    ! print '(A5,I6)',"N " , N
    ! print '(A5,I6)',"NF " , NF
    ! print '(A5,I6)',"i " , i
    ! print '(A5,2F6.3)',"YVAL " , YVAL_2
    
    ! Check if something changed
    do i = 1, N
        if (YVAL(i) /= YVAL_2(i)) then
            print *, fg_color_red + style_bold // test_ error // style_reset
            print *, "At: ", i, "| ", YVAL(i), " /= ", YVAL_2(i)
            ERROR STOP "YVAL /= YVAL_2"
        end if
    end do

        

end program check
    

