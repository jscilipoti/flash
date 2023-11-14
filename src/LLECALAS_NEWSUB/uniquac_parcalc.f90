subroutine uniquac_parcalc
    
    !implicit none
    use iso_fortran_env, only: real64

    IMPLICIT real(kind=real64) (A-H,O-Z)

    COMMON/CUFAC/N,NG,P(10,10),T
    COMMON/COUT/iout
    COMMON/CGIBBS/NF,z_max_index,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AL(10),&
    & DA(10,10),XM(10,4)
    COMMON/CY/Y13,Y21,STEP
    COMMON/CA/XC(5),GE(5,2),GC(5,2)
    COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)
    COMMON/CIPR/ipr
    COMMON/CQT/QT(10,10),Q(10),R(10)

    DIMENSION X(2)
    DIMENSION Y(10)!,DLX(10),YVAL(30),GRAD(30),XMAT(30,30),WORK(30,5) 

    if (N /= 2) write(*, 616)                                           
    if (iout /= 6 .and. N /= 2) write(iout, 616)                          
    XC(1) = 0.D0                                                          
    XC(2) = 0.2D0                                                        
    XC(3) = 0.5D0                                                        
    XC(4) = 0.8D0                                                        
    XC(5) = 1.D0                                                        
    
    do 17 k = 1, 5                                                       
        Y(1) = XC(k)                                                        
        Y(2) = 1.D0 - XC(k)                                                   
        call unifac(1, Y, ACT1, DACT1, PACT)                                  
        GE(K, 1) = ACT1(1)                                                   
    17  GE(K, 2) = ACT1(2)

    READ(2, *) R(1), Q(1)                                               
    READ(2, *) R(2), Q(2)                                               
                                
    write(*, 627)                                                      
    
    do 14 i = 1, 2                                                       
    14   write(*, 626) i, R(i), Q(i)                                          
    
    if (iout /= 6) then !if (iout == 6) GOTO 13                                             
        write(iout, 627)                                                   
        do 11 i = 1, 2                                                       
    11  write(iout, 626) i, R(i), Q(i)                                        
    !13 CONTINUE
    end if

    X(1) = Z(1) / 300.D0                                                  
    X(2) = Z(2) / 300.D0                                                  
    do 18 i = 1, 2                                                       
        do 18 j = 1, 2                                                       
            QT(i, j) = 0.D0                                                        
    18       P(i, j) = 0.D0                                                         
    QT(1, 1) = Q(1)                                                      
    QT(2, 2) = Q(2)                                                      
    NK = 2                                                              
    NG = 2                                                              
    XLAMB = 1.D0                                                          
    call MARQ(FUNC, 2, 10, X, XLAMB, 3.D0, 1.D-7, 99)                        
    write(*, 633) T                                                    
    if (iout /= 6) write(iout, 633) T                                   
    write(*, 617) P(1, 2), P(2, 1)     
    if (IPR == 1) write(*, 618)                                         
    do 21 L = 1, 5                                                       
        do 21 i = 1, 2                                                       
            GE(L,i) = DEXP(GE(L, i))                                             
21       GC(L, i) = DEXP(GC(L, i))                                             
    if (IPR == 1) write(*, 619) ((GE(L, i), L = 1, 5), i = 1, 2)                 
    if (IPR == 1) write(*, 619) ((GC(L, i), L = 1, 5), i = 1, 2)                 
    
    if (iout /= 6) then !if (iout == 6)                                             
        write(iout, 617) P(1, 2), P(2, 1)                                     
        if (IPR == 1) &
        & write(iout, 618)                                      
        if (IPR == 1) &
        & write(iout, 619) ((GE(L, i), L = 1, 5), i = 1, 2)              
        if (IPR == 1) & 
        & write(iout, 619) ((GC(L, i), L = 1, 5), i = 1, 2)              

    end if
    return

    616 FORMAT(//,' * WRONG INPUT SPECIFICATION *',//) 
    617 FORMAT(///,' ** UNIQUAC PARAMETERS FROM UNIFAC **',//,5X,'A12/R =  ',F12.3,' K ,  A21/R = ',F12.3,' K',///)
    618 FORMAT(//,' ** COMPARISON OF ACTIVITIES CALCULATED BY UNIFAC AND UNIQUAC, RESPECTIVELY **'//)                                       
619 FORMAT(10F12.5)                          
    626 FORMAT(I5,2F15.4) 
    627 FORMAT(//,' SPECIFIED UNIQUAC R AND Q',/) 
    633 FORMAT(///,'   TEMPERATURE =',F10.2,' DEG K')

end subroutine uniquac_parcalc