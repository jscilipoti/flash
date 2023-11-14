subroutine binodal_calc
    
    !implicit none
    use iso_fortran_env, only: real64

    IMPLICIT real(kind=real64) (A-H,O-Z)

    COMMON/CUFAC/N,NG,P(10,10),T
    COMMON/COUT/iout
    COMMON/CGIBBS/NF,z_max_index,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AL(10),&
    & DA(10,10),XM(10,4)
    COMMON/CY/Y13,Y21,STEP
    

    if (N /= 2 .and. N /= 3) write(*, 616)                                
    if (iout /= 6 .and. N /= 2 .and. N /= 3) write(iout, 616)               
    Y13 = Z(1)                                                          
    Y21 = Z(2)                                                          
                
    write(*,633) T                                                    
                
    if (iout /= 6) write(iout, 633) T                                   
    if (N == 3) then !if (N == 3) GOTO 12                                                
                    !12 STEP=Z(3) / 100.D0
    STEP = Z(3) / 100.D0

    if (STEP == 0.D0) STEP = 0.02D0                                         
        call BINOD
        return
    end if
    call SOLBIN
    return
    616 FORMAT(//,' * WRONG INPUT SPECIFICATION *',//) 
    633 FORMAT(///,'   TEMPERATURE =',F10.2,' DEG K')
end subroutine binodal_calc