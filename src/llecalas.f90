! *****************************************************************************
! *   Subroutine: llecalas(Tf, Pf, Zf)                                        *
! *****************************************************************************
!
subroutine llecalas!(Tf, Pf, Zf) 
!
! *****************************************************************************
! *                                                                           *
! *  SUBROUTINE   L L E C A L A S                                                *
! *  (asociacion incorporada para el calculo de flash y                       *
! *  curva binodal (icalc 0 y 1)                                              *
! *                                                                           *
! *                           ASOCIACION CRUZADA                              *
! *                          VERSION GENERALIZADA                             *
! *                              Junio 2023                                   *
! *   Modificada por:                                                         *
! *                        ALFONSINA ESTER ANDREATTA                          *
! *                          JOSE ANTONIO SCILIPOTI                           *
! *                           JUAN PABLO ROVEZZI                              *
! *   Revisada en Octubre del 2007 en el chequeo de estabilidad               *
! *                                                                           *
! *   Basada en las simplificaciones de los papers:                           *
! *      # Michelsen, et al. (Fluid Phase Equilibria, 180(2001)165-174 )      *
! *      # Tan, et al.  (Ind. Eng. Chem. Res, 2004, 43, 203-208).             *
! *                                                                           *
! *   Esto permitio  que todos los casos particulares de asociacion se puedan *
! *   simplificar a un unico calculo.                                         *
! *                                                                           *
! *   Valido para un maximo numero grupo asociativo de 12.                    *
! *   Con la implementacion en el calculo de la fraccion no asociada en el    *
! *   componente puro por  el metodo iterativo aqui implementado se permite   *
! *   que una molecula tenga mas de un grupo asociativo (14/07/06).           *
! *                                                                           * 
! *   El calculo se limita a que el numero maximo de sitios sea dos           * 
! *   (por razones matematicas).                                              *
! *                                                                           *                                                   
! *****************************************************************************  
! *                        DATE: 24/3/1982 /TJ                                *
! *****************************************************************************
! 
! Tf: Temperature
! Pf: Pression
! Zf: Composition for each component
!  
! Required modules:
    use InputData
    use iso_fortran_env, only: real64, int16

    ! Implicit statements need to be removed. It must use 'implicit none'
    IMPLICIT real(kind=real64) (A-H,O-Z)                                          

    ! External functions
    !EXTERNAL STABIL,GMIX ! Not used
    EXTERNAL FUNC                                         

    ! Common blocks that are meant to be removed
                            
    COMMON/CGIBBS/NF,z_max_index,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AL(10),&
    & DA(10,10),XM(10,4)                                       
    COMMON/CUFAC/N,NG,P(10,10),T                                      
    COMMON/CY/Y13,Y21,STEP                                            
    COMMON/CA/XC(5),GE(5,2),GC(5,2)                                   
                                                    
    COMMON/CQT/QT(10,10),Q(10),R(10)                                  
    COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                                             
    common/nga/nga,mass(12)




    COMMON/CMODEL/model_common 
    COMMON/CIPR/ipr_common 
    COMMON/COUT/iout_common
    COMMON/CVAP/novap_common,NDUM,IDUM(4),PRAT(10) 
    common/ig/ig_common

    !   
    DIMENSION DLX(10),YVAL(30),Y(10),GRAD(30),XMAT(30,30),WORK(30,5)  
    !DIMENSION NTEXT(36),
    DIMENSION X(2),ANT(10,3)
    dimension xmj(10),actgam(10),agam(10,4),de(10,10),pe(2,2)          

    ! These integer are meant to be removed
    integer :: &
    & ICALC_common,MODEL_common,IPR_common,IOUT_common,NOVAP_common,ig_common

    ! These integer are meant to be defined
    integer :: & 
    & z_max_index            
    
    integer(kind=int16):: NN = 0

                                                                         

    real(kind=real64) :: &

    ! The pressure of the system to be read in the input flash file
    & PP = 0.D0, &

    ! Apparentlym, this is a value to test what happends at this T
    & T1 = 0.D0, & !T1=0.

    ! The maximum value among all the molar fractions Z(i) of the components to
    ! be calculated from Z array
    & z_max = 0.D0, &

    ! The sum of all the molar fractions Z(i) to be calculated from Z array.
    & z_sum = 0.D0

!---START-------------------------------------------------------------------------

    call read_input_flash("name.dat")

    ICALC_common = icalc
    MODEL_common = model
    IPR_common = ipr
    IOUT_common = iout
    NOVAP_common = novap
    ig_common = ig


    call open_database(model)

        if (iout == 1) then
            open (unit = 1, file = 'lleasoccuzada.OUT', form = 'FORMATTED')
        end if
        
        open (unit = 3, file = 'output.OUT', form = 'FORMATTED')
        
        output = 3
        
        ! Write to screen
        write(6, 608)                                                      
        write(6, 628) iout                                                 
        write(6, 610)                                                      
        
        if(iout == 0) iout = 6                                              
        if(icalc == 0) write(6, 620)                                       
        if(icalc == 1) write(6, 621)                                       
        if(icalc == 2) write(6, 622)                                       
        if(novap /= 0) write(6, 629)                                       
        if(model == 0) write(6, 624)                                       
        if(model == 1) write(6, 625)                                       
        
        write(6,623) NTEXT                                                
        
        !if(iout == 6) GOTO 5
        if (iout /= 6) then                                              
        
        ! write to iout
        write(iout, 608)                                                   
        write(iout, 610)                                                   
        
        if(icalc == 0) write(iout,620)                                    
        if(icalc == 1) write(iout,621)                                    
        if(icalc == 2) write(iout,622)                                    
        if(novap /= 0) write(iout,629)                                    
        if(model == 0) write(iout,624)                                    
        if(model == 1) write(iout,625)                                    
        write(iout,623) NTEXT                                             
        
        end if
        !5 CONTINUE                                                          
        
        call PARIN2                                                       
        
        !-----------------------------------------------------------------------
        !Testear y luego modificar
        ! If no vapour phase is considered...
        if(novap /= 0 ) then                                             
            do 6 j=1,N                                                                                 
        6       READ(2,*) (ANT(K,j),K=1,3)                                        
            do 7 j=1,N                                                        
                ANT(1,j)=2.302585*(ANT(1,j)-2.880814)                             
        7       ANT(2,j)=2.302585*ANT(2,j)                                        
        endif                                                          
        !-----------------------------------------------------------------------
        
        T1=0.                                                             
        NN=0                                                              

        ! It checks whether the problem is a FLASH CALCULATION (0) or a 
        ! BINODAL CURVE CALCULATION OR CALCULATION OF UNIQUAC PARAMETERS 
        ! FROM UNIFAC. After that it reads the Temperature and Pressure.

        10  if (icalc == 0) then 
                READ(2,*) T, PP
            else !if (icalc .GE. 1) READ(2,*) T
                READ(2,*) T 
            end if                                   

        
        ! Exit the subroutine if the temperature has not been set 
        if(T == 0.) then 
            call close_llecalas() !GOTO 10000
            return
        end if

        ! If no vapour phase is considered...
        if(PP /= 0.D0 .and. novap /= 0) then                                 
            do i=1,N                                                        
                PRAT(i)=DLOG(PP)-ANT(1,i)+ANT(2,i)/(T-273.15+ANT(3,i))
            end do                                     
        end if
        
        ! Read the composition (Z) of each component (Zi)
        4 READ(2,*) (Z(i),i=1,N)                                            
        
        !z_sum = 0.D0                                                       
        !z_max = 0.D0                                                          
        
        ! It sums the compostion of each component to get Z_sum and the max
        ! value of composition and its index.
        do i = 1, N !15 i = 1, N                                                       
            z_sum = z_sum + z(i)                                                    
            if(Z(i) < z_max) cycle !GOTO 15                                          
            z_max = z(i)                                                         
            z_max_index = i                                                            
        !15 CONTINUE
        end do

        ! The following line has no sense since, as far as I know, T1
        ! had been initialized at the beggining as T1 = 0.
        if(T == T1) GOTO 30                                               
        
        call PARAM2                                                       
        if(icalc.NE.1) GOTO 16                                            
        if(N.NE.2.AND.N.NE.3) write(6,616)                                
        if(iout.NE.6.AND.N.NE.2.AND.N.NE.3) write(iout,616)               
        Y13=Z(1)                                                          
        Y21=Z(2)                                                          
        write(6,633) T                                                    
        if(iout.NE.6) write(iout,633) T                                   
        if(N == 3) GOTO 12                                                
        call SOLBIN                                                       
        call close_llecalas() !GOTO 10000
        return
                                                                
    12 STEP=Z(3)/100.D0                                                  
        if(STEP == 0.) STEP=.02D0                                         
        call BINOD                                                        
        call close_llecalas() !GOTO 10000
        return                                                        
    16 CONTINUE                                                          
        if(icalc.NE.2) GOTO 19                                            
        if(N.NE.2) write(6,616)                                           
        if(iout.NE.6.AND.N.NE.2) write(iout,616)                          
        XC(1)=0.                                                          
        XC(2)=.2D0                                                        
        XC(3)=.5D0                                                        
        XC(4)=.8D0                                                        
        XC(5)=1.D0                                                        
        do 17 K=1,5                                                       
            Y(1)=XC(K)                                                        
            Y(2)=1.D0-XC(K)                                                   
            call unifac(1,Y,ACT1,DACT1,PACT)                                  
            GE(K,1)=ACT1(1)                                                   
    17   GE(K,2)=ACT1(2)                                                   
        READ(2,*) R(1),Q(1)                                               
        READ(2,*) R(2),Q(2)                                               
    !C     READ(2,502) R(1),Q(1)                                             
    !C     READ(2,502) R(2),Q(2)                                             
        write(6,627)                                                      
        do 14 i=1,2                                                       
    14   write(6,626) i,R(i),Q(i)                                          
        if(iout == 6) GOTO 13                                             
        write(iout,627)                                                   
        do 11 i=1,2                                                       
    11   write(iout,626) i,R(i),Q(i)                                       
    13 CONTINUE                                                          
        X(1)=Z(1)/300.D0                                                  
        X(2)=Z(2)/300.D0                                                  
        do 18 i=1,2                                                       
            do 18 j=1,2                                                       
                QT(i,j)=0.                                                        
    18       P(i,j)=0.                                                         
        QT(1,1)=Q(1)                                                      
        QT(2,2)=Q(2)                                                      
        NK=2                                                              
        NG=2                                                              
        XLAMB=1.                                                          
        call MARQ(FUNC,2,10,X,XLAMB,3.D0,1.D-7,99)                        
        write(6,633) T                                                    
        if(iout.NE.6) write(iout,633) T                                   
        write(6,617) P(1,2),P(2,1)       !(///,' ** UNIQUAC PARAMETERS FROM UNIFAC **',//,5X,'A12/R =  ',F12.3,' K ,  A21/R = ',F12.3,' K',///)                                  
        if(IPR == 1) write(6,618)                                         
        do 21 L=1,5                                                       
            do 21 i=1,2                                                       
                GE(L,i)=DEXP(GE(L,i))                                             
    21       GC(L,i)=DEXP(GC(L,i))                                             
        if(IPR == 1) write(6,619) ((GE(L,i),L=1,5),i=1,2)                 
        if(IPR == 1) write(6,619) ((GC(L,i),L=1,5),i=1,2)                 
        if(iout == 6) GOTO 22                                             
        write(iout,617) P(1,2),P(2,1)                                     
        if(IPR == 1) write(iout,618)                                      
        if(IPR == 1) write(iout,619) ((GE(L,i),L=1,5),i=1,2)              
        if(IPR == 1) write(iout,619) ((GC(L,i),L=1,5),i=1,2)              
    22 CONTINUE                                                          
        call close_llecalas() !GOTO 10000
        return

    19 CONTINUE                                                          
        do 20 i=1,N                                                       
            do 20 j=1,N                                                       
                GAM(i,j)=0.D0                                                     
                if(j == i) GOTO 20                                                
                call GAMINF(i,j,G)                                                
                GAM(i,j)=G                                                        
    20 CONTINUE                                                          
    30 T1=T                                                              
        NN=NN+1                                                           
        write(6,602) NN                                                   
        do 35 i=1,N                                                       
    35 Z(i)=Z(i)/z_sum                                                    
        write(6,605) T,PP,z_sum,(Z(i),i=1,N)                               
        if(iout.NE.6) write(iout,602) NN                                  
        if(iout.NE.6) write(iout,605) T,PP,z_sum,(Z(i),i=1,N)              
        call unifac(1,Z,AL,DA,PACT)                                       
        SFAS(1)=1.                                                        
        GNUL=0.                                                           
        do 40 i=1,N                                                       
            XVL(i,1)=1.                                                       
            Z(i)=Z(i)+1.D-20                                                  
            DLX(i)=DLOG(Z(i))                                                 
            A(i)=AL(i)+DLX(i)                                                 
    40   GNUL=GNUL+Z(i)*AL(i)                                              
        NF=1                                                              
    50 call STIG(Y,S)                                                    
        if(S.GT.-1.D-7) GOTO 70                                           
        write(6,603)                                                      
        if(iout.NE.6) write(iout,603)                                     
        do 60 i=1,N                                                       
            YVAL(i)=1.D-5*Y(i)/Z(i)                                           
    60 CONTINUE                                                          
        GOTO 100                                                          
    70 do 75 i=1,N                                                       
    75   YVAL(i)=DLOG(Y(i))                                                
        XLAM=1.                                                           
        if(NF == 1.AND.IPR.GT.0) write(6,606)                             
        if(NF.GT.1.AND.IPR.GT.0) write(6,609) NF                          
        if(iout.NE.6.AND.NF == 1.AND.IPR.GT.0) write(iout,606)            
        if(iout.NE.6.AND.NF.GT.1.AND.IPR.GT.0) write(iout,609) NF         
        call TMSSJ(30,N,IPR,15,XLAM,1.D-12,FUN,YVAL,GRAD,XMAT,WORK,1)     
        if(FUN < -1.D-7) GOTO 80                                         
        write(6,604)         
        write(output,*) 1
            write(output,2613) (Z(j),j=1,N)
            write(output,2613) (AL(j),j=1,N)        
        write(output,*) "SYSTEM IS STABLE"                                                   
    
        ! write(7,46) T,  (xM(l,1),l=1,N)     !Alfonsina
        ! write(7,46) T,  (xM(l,2),l=1,N)     !Alfonsina
        ! write(7,*)                          !Alfonsina
    
    
    
        if(iout.NE.6) write(iout,604)                                     
        GOTO 10                                                           
    80 write(6,603)                                                      
        if(iout.NE.6) write(iout,603)                                     
        do 90 i=1,N                                                       
    90   YVAL(i)=1.D-5*DEXP(YVAL(i))/Z(i)                                  
    100 NF=NF+1                                                           
    104 do 105 i=1,N                                                      
            if(YVAL(i).GT.1.D0) GOTO 106                                      
    105 CONTINUE                                                          
        GOTO 109                                                          
    106 do 107 i=1,N                                                      
    107   YVAL(i)=YVAL(i)/10.                                               
        GOTO 104                                                          
    109 CONTINUE                                                          
        SFAS(NF)=1.                                                       
        XLAM=.2                                                           
        if(NF == 2) XLAM=.5                                               
        M=(NF-1)*N                                                        
        if(IPR.GT.0) write(6,607) NF                                      
        if(iout.NE.6.AND.IPR.GT.0) write(iout,607) NF                     
        call TMSSJ(30,M,IPR,60,XLAM,1.D-16,FUN,YVAL,GRAD,XMAT,WORK,2)     
        NT=NF*N                                                           
        NB=NT-N                                                           
        do 110 i=1,NB                                                     
    110   YVAL(NT+1-i)=YVAL(NB+1-i)                                         
        write(6,614) NF                                                   
        NVAP=0                                                            
        do 111 j=1,NF                                                     
            if(IDUM(j) == 1) NVAP=j                                           
    111 CONTINUE                                                          
        if(NVAP == 0) write(6,630)                                        
        if(NVAP.NE.0) write(6,631) NVAP                                   
        if(iout.NE.6.AND.NVAP == 0) write(iout,630)                       
        if(iout.NE.6.AND.NVAP.NE.0) write(iout,631) NVAP                  
        write(6,611) (j,SFAS(j),j=1,NF)                                   
        write(6,612) (j,j=1,NF)                                           
        if(iout.NE.6) write(iout,614) NF                                  
        if(iout.NE.6) write(iout,611)(j,SFAS(j),j=1,NF)                   
        if(iout.NE.6) write(iout,612) (j,j=1,NF)                          
        SUM=0.                                                            
        do 115 i=1,N                                                      
            DLX(i)=XVL(i,NF)*Z(i)/SFAS(NF)                                    
    115   SUM=SUM+DLX(i)                                                    
        SUM=DLOG(SUM)                                                     
        call unifac(1,DLX,A,DA,PACT)                                      
        do 120 i=1,N                                                      
            DLX(i)=DLOG(DLX(i))                                               
    120   A(i)=A(i)+DLX(i)-SUM                                              
    !c-----
        do 1130 j=1,nf
            do 1131 i=1,n
    1131       xmj(i)=xm(i,j)
            call unifac(1,xmj,actgam,de,pe)
            do 1132 i=1,n
    1132       agam(i,j)=actgam(i)
    1130 continue
        write (output,*) NF
        do i=1,NF !escribe resultados para el output que ser� le�do por excel
            write(output,2613) (XM(j,i),j=1,N)
            write(output,2613) (agam(j,i),j=1,N)  
        enddo
    !      do i=1,NF !escribe resultados para el output que ser� le�do por excel
    !        write(output,2613) (agam(j,i),j=1,N)      
    !      enddo
        
        do 130 i=1,N                                                      
            write(6,613) i,(XM(i,j),j=1,NF)     !composition        
    130   write(6,1613) i,(agam(i,j),j=1,nf) !Ln(gamma)
        
        if(iout == 6) GOTO 132                                            
        do 131 i=1,N                                                      
        write(iout,613) i,(XM(i,j),j=1,NF)    
    131 write(iout,1613) i,(agam(i,j),j=1,nf)
    46	FORMAT (2X,F12.2, 8X,F12.6, 8X,F12.6 , 8X,F12.6, 8X,F12.6, & !Alfonsina                                                 
        8X,F12.6,8X,F12.6,8X,F12.6,8X,F12.6,8X,F12.6,8X,F12.6,8X,F12.6)
    !c-----
    132 CONTINUE                                                          
        GOTO 50                                                           
    
    
    501 FORMAT(36A2)                                                      
    502 FORMAT(8F10.2)                                                    
    503 FORMAT(20I3)                                                      
    602 FORMAT(///,' * * * FLASH NUMBER',I3,' * * *',//)                  
    603 FORMAT(/,' SYSTEM IS UNSTABLE, PHASE SPLIT PERFORMED')            
    604 FORMAT(/,' * SYSTEM IS STABLE *',/)                               
    605 FORMAT(' TEMPERATURE =',F10.4,' K, PRESSURE =',F7.3,' ATM, FEED ='&
        ,F10.2,' MOLES',/,' FEED COMPOSITION (MOLE PERCENT):',/,1X,15(2PF7&
        .3))                                                              
    606 FORMAT(//,' DIFFERENTIAL STABILITY TEST FOR FEED MIXTURE:')       
    607 FORMAT(/,' PHASE SPLIT CALCULATION,',I2,' PHASES:')               
    608 FORMAT(1H1)                                                       
    609 FORMAT(//,' DIFFERENTIAL STABILITY TEST FOR',I2,'-PHASE SYSTEM')  
    610 FORMAT(///)                                                       
    611 FORMAT(/,'  PHASE FRACTIONS (PERCENT):',4(5X,I3,2PF7.3,5X))       
    612 FORMAT(/,'  COMPOSITION  ',10X,4(8X,I3,9X))                       
    613 FORMAT('   X(',I2,')            ',5(8X,F12.8))                    
    !c-----
    1613 format('  ln(G',i2,')            ',5(8x,f12.8))
    2613 format(5(2x,f12.8))
    
    !c-----
    614 FORMAT(//,' RESULT OF',I2,'-PHASE CALCULATION:')                  
    616 FORMAT(//,' * WRONG INPUT SPECIFICATION *',//)                    
    617 FORMAT(///,' ** UNIQUAC PARAMETERS FROM UNIFAC **',//,5X,'A12/R =  ',F12.3,' K ,  A21/R = ',F12.3,' K',///)                         
    618 FORMAT(//,' ** COMPARISON OF ACTIVITIES CALCULATED BY UNIFAC AND UNIQUAC, RESPECTIVELY **'//)                                       
    619 FORMAT(10F12.5)                                                   
    620 FORMAT(' **** FLASH CALCULATION ****')                            
    621 FORMAT(' **** BINODAL CURVE CALCULATION ****',//)                 
    622 FORMAT(' **** CALCULATION OF UNIQUAC PARAMETERS FROM UNIFAC **** ',//)                                                              
    623 FORMAT(1X,'COMPONENTS : ',40A2,//)                                
    624 FORMAT(' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIFAC'//)     
    625 FORMAT(' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIQUAC'//)    
    626 FORMAT(I5,2F15.4)                                                 
    627 FORMAT(//,' SPECIFIED UNIQUAC R AND Q',/)                         
    628 FORMAT(/,' iout = ',I2,/' if iout = 0: OUTPUT ONLY ON UNIT 6',/,  ' if iout = 1: OUTPUT ON BOTH UNIT 6 AND 1')                      
    629 FORMAT(/,' VAPOR PHASE INCLUDED IN FLASH-CALCULATIONS',//)        
    630 FORMAT(' NO VAPOR PHASE')                                         
    631 FORMAT(' PHASE',I2,' IS A VAPOR PHASE')                           
    633 FORMAT(///,'   TEMPERATURE =',F10.2,' DEG K')                     
    
    ! Close all opened units in this subroutine 
    call close_llecalas()
    
    ! 10000 CLOSE (UNIT=2)
    !     if (iout == 1) CLOSE (UNIT=1)
    !     close (unit=3)
    ! !c      call salida(name)
    !     !STOP

    return

end SUBROUTINE llecalas   
                      
!---END-------------------------------------------------------------------------

    ! Subroutine close_llecalas: This subroutine replaces the original 
    ! "goto 10000" statement. It closes Units: 1 (lleasoccuzada.OUT),
    ! 2 (llecalas.dat) and 3 (output.out).
    subroutine close_llecalas()
        use iso_fortran_env, only: int8
        use inputData, only: iout

        implicit none
        
        if (iout == 1) close (unit=1)
        close (unit=2)
        close (unit=3)
        return

    end subroutine close_llecalas

!-------------------------------------------------------------------------------