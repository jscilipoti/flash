! BLoque raro de interpretar
    program check
    

    use iso_fortran_env, only: int8, real64
    implicit none

    integer(kind=int8) ::  &
    & NF = 0, &
    & i = 0

    integer(kind=int8), parameter :: &
    & N = 2

    real(kind=real64), dimension(N) :: &
    & YVAL = 2.D0, &
    & YVAL_2 = 2.D0

    !Codigo original
    100 NF = NF + 1                                                           
    104 do 105 i = 1, N                                                      
            if (YVAL(i) > 1.D0) then 
                print *, "GOTO 106"  
                GOTO 106
            end if                                      
    105 CONTINUE
        print *, "GOTO 109"                                                          
        GOTO 109                                                          
    106 do 107 i = 1, N                                                      
    107   YVAL(i) = YVAL(i) / 10.                                               
        print *, "GOTO 104"  
        GOTO 104                                                          
    109 CONTINUE

    print '(A5,I6)',"N " , N
    print '(A5,I6)',"NF " , NF
    print '(A5,I6)',"i " , i
    print '(A5,2F6.3)',"YVAL " , YVAL

    !Nuevo Código
    NF = 0
    i = 0

    1000 NF = NF + 1                                                           
    1004 do 1005 i = 1, N                                                      
            if (YVAL_2(i) > 1.D0) then 
                print *, "GOTO 106"  
                GOTO 1006
            end if                                      
    1005 CONTINUE
        print *, "GOTO 109"                                                          
        GOTO 1009                                                          
    1006 do 1007 i = 1, N                                                      
    1007   YVAL_2(i) = YVAL_2(i) / 10.D0                                               
        print *, "GOTO 104"  
        GOTO 1004                                                          
    1009 CONTINUE
    print '(A5,I6)',"N " , N
    print '(A5,I6)',"NF " , NF
    print '(A5,I6)',"i " , i
    print '(A5,2F6.3)',"YVAL " , YVAL_2
    do i = 1, N
        if (YVAL(i) /= YVAL_2(i)) then
            print *, "ERROR"
            print *, "At: ", i, "| ", YVAL(i), " /= ", YVAL_2(i)
            ERROR STOP "YVAL /= YVAL_2"
        end if
    end do
        

end program check
    

