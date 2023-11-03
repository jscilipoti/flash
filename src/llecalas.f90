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

! Implicit statements need to be removed. It must use 'implicit none'
IMPLICIT REAL*8(A-H,O-Z)                                          

! External functions
EXTERNAL STABIL,GMIX,FUNC                                         

! Common blocks that are meant to be removed
COMMON/CVAP/NOVAP,NDUM,IDUM(4),PRAT(10)                           
COMMON/CGIBBS/NF,MAXZ,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AL(10),DA(10,10),XM(10,4)                                       
COMMON/CUFAC/N,NG,P(10,10),T                                      
COMMON/CY/Y13,Y21,STEP                                            
COMMON/CA/XC(5),GE(5,2),GC(5,2)                                   
COMMON/CIPR/IPR                                                   
COMMON/CQT/QT(10,10),Q(10),R(10)                                  
COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
COMMON/CMODEL/MODEL                                               
COMMON/COUT/IOUT                                                  
common/nga/nga,mass(12)
common/ig/ig

!   
DIMENSION DLX(10),YVAL(30),Y(10),GRAD(30),XMAT(30,30),WORK(30,5)  
DIMENSION NTEXT(36),X(2),ANT(10,3)
dimension xmj(10),actgam(10),agam(10,4),de(10,10),pe(2,2)          

integer::ICALC,MODEL,IPR,IOUT,NOVAP,ig            
character(len=36)::name, name1 
integer:: parameters 
    

    OPEN (UNIT=1,FILE ='name.dat',status='OLD',FORM='FORMATTED')
    read(1,*)parameters 
    read(1,"(A36)") name
    name = name(2:len_trim(name)-1)
    if (parameters==1)then
        call get_database_data(name)
        !stop
        return
    endif    
    CLOSE (UNIT=1)
            
    OPEN (UNIT=2,FILE=name,status='OLD',FORM='FORMATTED')
    READ(2,501) NTEXT                                                                             
    READ(2,*) ICALC,MODEL,IPR,IOUT,NOVAP,ig, ipareq 

    call open_database(model)

        IF (IOUT.EQ.1) OPEN (UNIT=1,FILE='lleasoccuzada.OUT',FORM='FORMATTED')
        OPEN (UNIT=3,FILE='output.OUT',FORM='FORMATTED')
        output=3
        WRITE(6,608)                                                      
        WRITE(6,628) IOUT                                                 
        WRITE(6,610)                                                      
        IF(IOUT.EQ.0) IOUT=6                                              
        IF(ICALC.EQ.0) WRITE(6,620)                                       
        IF(ICALC.EQ.1) WRITE(6,621)                                       
        IF(ICALC.EQ.2) WRITE(6,622)                                       
        IF(NOVAP.NE.0) WRITE(6,629)                                       
        IF(MODEL.EQ.0) WRITE(6,624)                                       
        IF(MODEL.EQ.1) WRITE(6,625)                                       
        WRITE(6,623) NTEXT                                                
        IF(IOUT.EQ.6) GOTO 5                                              
        WRITE(IOUT,608)                                                   
        WRITE(IOUT,610)                                                   
        IF(ICALC.EQ.0) WRITE(IOUT,620)                                    
        IF(ICALC.EQ.1) WRITE(IOUT,621)                                    
        IF(ICALC.EQ.2) WRITE(IOUT,622)                                    
        IF(NOVAP.NE.0) WRITE(IOUT,629)                                    
        IF(MODEL.EQ.0) WRITE(IOUT,624)                                    
        IF(MODEL.EQ.1) WRITE(IOUT,625)                                    
        WRITE(IOUT,623) NTEXT                                             
        5 CONTINUE                                                          
        CALL PARIN2                                                       
        IF(NOVAP/=0) then                                             
            DO 6 J=1,N                                                        
    !C   6 READ(2,502) (ANT(K,J),K=1,3)                                      
        6       READ(2,*) (ANT(K,J),K=1,3)                                        
            DO 7 J=1,N                                                        
                ANT(1,J)=2.302585*(ANT(1,J)-2.880814)                             
        7       ANT(2,J)=2.302585*ANT(2,J)                                        
        endif                                                          
        T1=0.                                                             
        NN=0                                                              
    !C  10 READ(2,502) T,PP                                                  
    10 IF (ICALC.EQ.0) READ(2,*) T,PP                                   
        IF (ICALC.GE.1) READ(2,*) T
        IF(T.EQ.0.) GOTO 10000 
        !IF(PP.EQ.0..OR.NOVAP.EQ.0) then !reemplazado por la siguiente estructura
        IF(PP/=0..and.NOVAP/=0) then                                 
            DO 1 I=1,N                                                        
        1       PRAT(I)=DLOG(PP)-ANT(1,I)+ANT(2,I)/(T-273.15+ANT(3,I))            
    !C   4 READ(2,502) (Z(I),I=1,N)                                          
        endif
        4 READ(2,*) (Z(I),I=1,N)                                            
        ZSUM=0.                                                           
        ZMAX=0.                                                           
        DO 15 I=1,N                                                       
            ZSUM=ZSUM+Z(I)                                                    
            IF(Z(I).LT.ZMAX) cycle !GOTO 15                                          
            ZMAX=Z(I)                                                         
            MAXZ=I                                                            
    15 CONTINUE                                                          
        IF(T.EQ.T1) GOTO 30                                               
        CALL PARAM2                                                       
        IF(ICALC.NE.1) GOTO 16                                            
        IF(N.NE.2.AND.N.NE.3) WRITE(6,616)                                
        IF(IOUT.NE.6.AND.N.NE.2.AND.N.NE.3) WRITE(IOUT,616)               
        Y13=Z(1)                                                          
        Y21=Z(2)                                                          
        WRITE(6,633) T                                                    
        IF(IOUT.NE.6) WRITE(IOUT,633) T                                   
        IF(N.EQ.3) GOTO 12                                                
        CALL SOLBIN                                                       
        GOTO 10000                                                        
    12 STEP=Z(3)/100.D0                                                  
        IF(STEP.EQ.0.) STEP=.02D0                                         
        CALL BINOD                                                        
        GOTO 10000                                                        
    16 CONTINUE                                                          
        IF(ICALC.NE.2) GOTO 19                                            
        IF(N.NE.2) WRITE(6,616)                                           
        IF(IOUT.NE.6.AND.N.NE.2) WRITE(IOUT,616)                          
        XC(1)=0.                                                          
        XC(2)=.2D0                                                        
        XC(3)=.5D0                                                        
        XC(4)=.8D0                                                        
        XC(5)=1.D0                                                        
        DO 17 K=1,5                                                       
            Y(1)=XC(K)                                                        
            Y(2)=1.D0-XC(K)                                                   
            CALL unifac(1,Y,ACT1,DACT1,PACT)                                  
            GE(K,1)=ACT1(1)                                                   
    17   GE(K,2)=ACT1(2)                                                   
        READ(2,*) R(1),Q(1)                                               
        READ(2,*) R(2),Q(2)                                               
    !C     READ(2,502) R(1),Q(1)                                             
    !C     READ(2,502) R(2),Q(2)                                             
        WRITE(6,627)                                                      
        DO 14 I=1,2                                                       
    14   WRITE(6,626) I,R(I),Q(I)                                          
        IF(IOUT.EQ.6) GOTO 13                                             
        WRITE(IOUT,627)                                                   
        DO 11 I=1,2                                                       
    11   WRITE(IOUT,626) I,R(I),Q(I)                                       
    13 CONTINUE                                                          
        X(1)=Z(1)/300.D0                                                  
        X(2)=Z(2)/300.D0                                                  
        DO 18 I=1,2                                                       
            DO 18 J=1,2                                                       
                QT(I,J)=0.                                                        
    18       P(I,J)=0.                                                         
        QT(1,1)=Q(1)                                                      
        QT(2,2)=Q(2)                                                      
        NK=2                                                              
        NG=2                                                              
        XLAMB=1.                                                          
        CALL MARQ(FUNC,2,10,X,XLAMB,3.D0,1.D-7,99)                        
        WRITE(6,633) T                                                    
        IF(IOUT.NE.6) WRITE(IOUT,633) T                                   
        WRITE(6,617) P(1,2),P(2,1)       !(///,' ** UNIQUAC PARAMETERS FROM UNIFAC **',//,5X,'A12/R =  ',F12.3,' K ,  A21/R = ',F12.3,' K',///)                                  
        IF(IPR.EQ.1) WRITE(6,618)                                         
        DO 21 L=1,5                                                       
            DO 21 I=1,2                                                       
                GE(L,I)=DEXP(GE(L,I))                                             
    21       GC(L,I)=DEXP(GC(L,I))                                             
        IF(IPR.EQ.1) WRITE(6,619) ((GE(L,I),L=1,5),I=1,2)                 
        IF(IPR.EQ.1) WRITE(6,619) ((GC(L,I),L=1,5),I=1,2)                 
        IF(IOUT.EQ.6) GOTO 22                                             
        WRITE(IOUT,617) P(1,2),P(2,1)                                     
        IF(IPR.EQ.1) WRITE(IOUT,618)                                      
        IF(IPR.EQ.1) WRITE(IOUT,619) ((GE(L,I),L=1,5),I=1,2)              
        IF(IPR.EQ.1) WRITE(IOUT,619) ((GC(L,I),L=1,5),I=1,2)              
    22 CONTINUE                                                          
        GOTO 10000                                                        
    19 CONTINUE                                                          
        DO 20 I=1,N                                                       
            DO 20 J=1,N                                                       
                GAM(I,J)=0.D0                                                     
                IF(J.EQ.I) GOTO 20                                                
                CALL GAMINF(I,J,G)                                                
                GAM(I,J)=G                                                        
    20 CONTINUE                                                          
    30 T1=T                                                              
        NN=NN+1                                                           
        WRITE(6,602) NN                                                   
        DO 35 I=1,N                                                       
    35 Z(I)=Z(I)/ZSUM                                                    
        WRITE(6,605) T,PP,ZSUM,(Z(I),I=1,N)                               
        IF(IOUT.NE.6) WRITE(IOUT,602) NN                                  
        IF(IOUT.NE.6) WRITE(IOUT,605) T,PP,ZSUM,(Z(I),I=1,N)              
        CALL unifac(1,Z,AL,DA,PACT)                                       
        SFAS(1)=1.                                                        
        GNUL=0.                                                           
        DO 40 I=1,N                                                       
            XVL(I,1)=1.                                                       
            Z(I)=Z(I)+1.D-20                                                  
            DLX(I)=DLOG(Z(I))                                                 
            A(I)=AL(I)+DLX(I)                                                 
    40   GNUL=GNUL+Z(I)*AL(I)                                              
        NF=1                                                              
    50 CALL STIG(Y,S)                                                    
        IF(S.GT.-1.D-7) GOTO 70                                           
        WRITE(6,603)                                                      
        IF(IOUT.NE.6) WRITE(IOUT,603)                                     
        DO 60 I=1,N                                                       
            YVAL(I)=1.D-5*Y(I)/Z(I)                                           
    60 CONTINUE                                                          
        GOTO 100                                                          
    70 DO 75 I=1,N                                                       
    75   YVAL(I)=DLOG(Y(I))                                                
        XLAM=1.                                                           
        IF(NF.EQ.1.AND.IPR.GT.0) WRITE(6,606)                             
        IF(NF.GT.1.AND.IPR.GT.0) WRITE(6,609) NF                          
        IF(IOUT.NE.6.AND.NF.EQ.1.AND.IPR.GT.0) WRITE(IOUT,606)            
        IF(IOUT.NE.6.AND.NF.GT.1.AND.IPR.GT.0) WRITE(IOUT,609) NF         
        CALL TMSSJ(30,N,IPR,15,XLAM,1.D-12,FUN,YVAL,GRAD,XMAT,WORK,1)     
        IF(FUN.LT.-1.D-7) GOTO 80                                         
        WRITE(6,604)         
        write(output,*) 1
            write(output,2613) (Z(j),J=1,N)
            write(output,2613) (AL(j),j=1,N)        
        write(output,*) "SYSTEM IS STABLE"                                                   
    
        write(7,46) T,  (xM(l,1),l=1,N)    !Alfonsina
        write(7,46) T,  (xM(l,2),l=1,N)!Alfonsina
        write(7,*)                    !Alfonsina
    
    
    
        IF(IOUT.NE.6) WRITE(IOUT,604)                                     
        GOTO 10                                                           
    80 WRITE(6,603)                                                      
        IF(IOUT.NE.6) WRITE(IOUT,603)                                     
        DO 90 I=1,N                                                       
    90   YVAL(I)=1.D-5*DEXP(YVAL(I))/Z(I)                                  
    100 NF=NF+1                                                           
    104 DO 105 I=1,N                                                      
            IF(YVAL(I).GT.1.D0) GOTO 106                                      
    105 CONTINUE                                                          
        GOTO 109                                                          
    106 DO 107 I=1,N                                                      
    107   YVAL(I)=YVAL(I)/10.                                               
        GOTO 104                                                          
    109 CONTINUE                                                          
        SFAS(NF)=1.                                                       
        XLAM=.2                                                           
        IF(NF.EQ.2) XLAM=.5                                               
        M=(NF-1)*N                                                        
        IF(IPR.GT.0) WRITE(6,607) NF                                      
        IF(IOUT.NE.6.AND.IPR.GT.0) WRITE(IOUT,607) NF                     
        CALL TMSSJ(30,M,IPR,60,XLAM,1.D-16,FUN,YVAL,GRAD,XMAT,WORK,2)     
        NT=NF*N                                                           
        NB=NT-N                                                           
        DO 110 I=1,NB                                                     
    110   YVAL(NT+1-I)=YVAL(NB+1-I)                                         
        WRITE(6,614) NF                                                   
        NVAP=0                                                            
        DO 111 J=1,NF                                                     
            IF(IDUM(J).EQ.1) NVAP=J                                           
    111 CONTINUE                                                          
        IF(NVAP.EQ.0) WRITE(6,630)                                        
        IF(NVAP.NE.0) WRITE(6,631) NVAP                                   
        IF(IOUT.NE.6.AND.NVAP.EQ.0) WRITE(IOUT,630)                       
        IF(IOUT.NE.6.AND.NVAP.NE.0) WRITE(IOUT,631) NVAP                  
        WRITE(6,611) (J,SFAS(J),J=1,NF)                                   
        WRITE(6,612) (J,J=1,NF)                                           
        IF(IOUT.NE.6) WRITE(IOUT,614) NF                                  
        IF(IOUT.NE.6) WRITE(IOUT,611)(J,SFAS(J),J=1,NF)                   
        IF(IOUT.NE.6) WRITE(IOUT,612) (J,J=1,NF)                          
        SUM=0.                                                            
        DO 115 I=1,N                                                      
            DLX(I)=XVL(I,NF)*Z(I)/SFAS(NF)                                    
    115   SUM=SUM+DLX(I)                                                    
        SUM=DLOG(SUM)                                                     
        CALL unifac(1,DLX,A,DA,PACT)                                      
        DO 120 I=1,N                                                      
            DLX(I)=DLOG(DLX(I))                                               
    120   A(I)=A(I)+DLX(I)-SUM                                              
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
            write(output,2613) (XM(j,i),J=1,N)
            write(output,2613) (agam(j,i),j=1,N)  
        enddo
    !      do i=1,NF !escribe resultados para el output que ser� le�do por excel
    !        write(output,2613) (agam(j,i),j=1,N)      
    !      enddo
        
        DO 130 I=1,N                                                      
            WRITE(6,613) I,(XM(I,J),J=1,NF)     !composition        
    130   write(6,1613) i,(agam(i,j),j=1,nf) !Ln(gamma)
        
        IF(IOUT.EQ.6) GOTO 132                                            
        DO 131 I=1,N                                                      
        WRITE(IOUT,613) I,(XM(I,J),J=1,NF)    
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
    628 FORMAT(/,' IOUT = ',I2,/' IF IOUT = 0: OUTPUT ONLY ON UNIT 6',/,  ' IF IOUT = 1: OUTPUT ON BOTH UNIT 6 AND 1')                      
    629 FORMAT(/,' VAPOR PHASE INCLUDED IN FLASH-CALCULATIONS',//)        
    630 FORMAT(' NO VAPOR PHASE')                                         
    631 FORMAT(' PHASE',I2,' IS A VAPOR PHASE')                           
    633 FORMAT(///,'   TEMPERATURE =',F10.2,' DEG K')                     
    10000 CLOSE (UNIT=2)
        IF (IOUT.EQ.1) CLOSE (UNIT=1)
        close (unit=3)
    !c      call salida(name)
        !STOP
        return                                                              
        end SUBROUTINE llecalas   
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                       
    !C******************************* F I N ***************************************
    !C
    