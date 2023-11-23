SUBROUTINE BINOD                                                  
    IMPLICIT REAL*8 (A-H,O-Z)                                         
    COMMON/CY/Y13,Y21,STEP                                            
    COMMON/COUT/IOUT                                                  
    DIMENSION YMAT(70,6),NITE3(70)                                    
    COMMON/CBISO/NIC1,NIC2,IC1(120),IC2(120)                          
    DIMENSION Y(4),YOLD(4),DY(4),DYOLD(4),DMAT(4,4)                   
    DIMENSION K(6)                                                    
!c-----
    dimension xandrea1(10),xandrea2(10),act1(10),act2(10),dact(10,10),pact(2,2),acttot(70,6)
    common/ig/ig
!c-----      
    DATA K/1,3,2,4,5,6/                                               
    WRITE(6,601)                                                      
    IF(IOUT.NE.6) WRITE(IOUT,601)                                     
601 FORMAT(///,10X,' ** BINODAL CURVE, CONCENTRATIONS IN MOLE PERCENT **',//,18X,'COMP. 1',13X,'COMP. 2',13X,'COMP. 3')                 
!c-----
    write(6,1601)
    if(iout.ne.6) write(iout,1601)
1601 format(18x,'gamma1',14x,'gamma2',14x,'gamma3',/)
!c-----     
    DO 1 I=1,4                                                        
    DO 1 J=1,4                                                        
1     DMAT(I,J)=0.D0                                                    
    IRUND=200                                                         
    NIC1=0                                                            
    NIC2=0                                                            
    ICOND=0                                                           
    NOLD=2                                                            
    N=0                                                               
    Y(1)=1.D0-Y13/100.D0                                              
    Y(2)=0.D0                                                         
    Y(3)=Y21/100.D0                                                   
    Y(4)=0.D0                                                         
 12 CALL SOLVE(Y,DY,NOLD,NEW,NITER,N,NT)                              
    IF(NITER.LE.NT) GOTO 16                                           
    IF(N.GT.0) GOTO 19                                                
    ICOND=-10                                                         
    WRITE(6,902)                                                      
    IF(IOUT.NE.6) WRITE(IOUT,902)                                     
902 FORMAT(/,' THE BASE LINE CALCULATION DID NOT CONVERGE IN 10 ITERATIONS'/)                                                           
    GOTO 3                                                            
 19 IF(IHALF.LT.5) GOTO 20                                            
    ICOND=-3                                                          
    GOTO 3                                                            
 20 IHALF=IHALF+1                                                     
    ST=ST/2.D0                                                        
    GOTO 17                                                           
 16 IF(DABS(Y(1)-Y(3))+DABS(Y(2)-Y(4)).GT.1.D-8) GOTO 21              
    IF(N.GT.0) GOTO 19                                                
    WRITE(6,903)                                                      
    IF(IOUT.NE.6) WRITE(IOUT,903)                                     
903 FORMAT(/,' THE CONCENTRATIONS ON THE BASE LINE ARE IDENTICAL IN THE TWO PHASES'/)                                                   
    GOTO 3                                                            
 21 N=N+1                                                             
    NITE3(N)=NITER                                                    
    IHALF=0                                                           
    DO 2 I=1,4                                                        
2     YMAT(N,I)=Y(I)                                                    
    IF(ICOND.EQ.2.AND.Y(1).LT.1.D-10) GOTO 3                          
    DYMAX=DABS(DY(NEW))                                               
    DO 4 I=1,4                                                        
4     DY(I)=DY(I)/DYMAX                                                 
    IF(N.EQ.1)GOTO 5                                                  
    STAP=DABS(Y(NEW)-YOLD(NEW))                                       
    IF(DY(NEW)*DYOLD(NEW).GT.0.D0)GOTO 6                              
    DO 7 I=1,4                                                        
7     DY(I)=-DY(I)                                                      
6     IF(NEW.EQ.NOLD)GOTO 8                                             
    RR=DY(NEW)/DYOLD(NEW)                                             
    DO 9 I=1,4                                                        
9     DYOLD(I)=DYOLD(I)*RR                                              
8     DO 10 I=1,4                                                       
    Z=(YOLD(I)-Y(I))/STAP                                             
    DMAT(I,3)=(3.D0*Z+2.D0*DY(I)+DYOLD(I))/STAP                       
10    DMAT(I,4)=(2.D0*Z+DY(I)+DYOLD(I))/STAP**2                         
5     ST=RUND(Y(NEW),DY(NEW),STEP,IRUND)                                
    DO 18 I=1,4                                                       
    DMAT(I,1)=Y(I)                                                    
    DMAT(I,2)=DY(I)                                                   
    YOLD(I)=Y(I)                                                      
18    DYOLD(I)=DY(I)                                                    
17    DO 11 I=1,4                                                       
    Y(I)=DMAT(I,4)                                                    
    DO 11 J=1,3                                                       
11    Y(I)=ST*Y(I)+DMAT(I,4-J)                                          
    IF(IHALF.GT.0)GOTO 12                                             
    CALL TERM(Y,DMAT,ICOND,NEW)                                       
    NOLD=NEW                                                          
    IF(ICOND.EQ.0.OR.ICOND.EQ.2) GOTO 12                              
3     NGIT=N                                                            
    IF(N.EQ.0) RETURN                                                 
    IF(ICOND.NE.1)GOTO 13                                             
    N=N+1                                                             
    DO 14 I=1,4                                                       
14    YMAT(N,I)=Y(I)                                                    
    NITE3(N)=0                                                        
13    IF(N.EQ.0)RETURN                                                  
    DO 15 I=1,N                                                       
    YMAT(I,5)=1.D0-YMAT(I,1)-YMAT(I,2)                                
15    YMAT(I,6)=1.D0-YMAT(I,3)-YMAT(I,4)                                
    DO 30 I=1,N                                                       

    xandrea1(1)=ymat(i,k(1))
    xandrea1(2)=ymat(i,k(3))
    xandrea1(3)=ymat(i,k(5))
    xandrea2(1)=ymat(i,k(2))
    xandrea2(2)=ymat(i,k(4))
    xandrea2(3)=ymat(i,k(6))
    call unifac(0,xandrea1,act1,dact,pact)
    call unifac(0,xandrea2,act2,dact,pact)
    acttot(i,1)=dexp(act1(1))
    acttot(i,2)=dexp(act2(1))
    acttot(i,3)=dexp(act1(2))
    acttot(i,4)=dexp(act2(2))
    acttot(i,5)=dexp(act1(3))
    acttot(i,6)=dexp(act2(3))

    WRITE(6,600) I,NITE3(I),(YMAT(I,K(J)),J=1,6)   
  WRITE(3,611) I,NITE3(I),(YMAT(I,K(J)),J=1,6)   
  
    write(8,47) YMAT (i,1),  YMAT (i,2), YMAT (i,5)    !Alfonsina
   write(8,47) YMAT (i,3),  YMAT (i,4), YMAT (i,6)  !Alfonsina
   write(8,*) 
   
 47	FORMAT (8X,F12.6 , 8X,F12.6 , 8X,F12.6) !Alfonsina                  

    if(ig.eq.1) write(6,1600) (acttot(i,j),j=1,6)
 30 continue

600 FORMAT(5X,2I3,2X,3(2P2F9.4,2X))     
611 FORMAT(5X,2I3,2X,3(2P2F9.4,2X))                                   
    IF(NIC1.GT.0) WRITE(6,900) IC1(1),IC1(NIC1)                       
1600 format(12x,6(f9.3,1x))
900 FORMAT(/,' FALSE SOLUTION IN PHASE 1 FOR CALCULATED TIE LINES',I3,' TO',I3,/)                                                       
    IF(NIC2.GT.0) WRITE(6,901) IC2(1),IC2(NIC2)                       
901 FORMAT(/,' FALSE SOLUTION IN PHASE 2 FOR CALCULATED TIE LINES',I3,' TO',I3,/)                                                       
    IF(IOUT.EQ.6) GOTO 50                                             
    DO 35 I=1,N                                                       
    WRITE(IOUT,600) I,NITE3(I),(YMAT(I,K(J)),J=1,6)                   
  


    if(ig.eq.1) write(iout,1600) (acttot(i,j),j=1,6)
 35 continue

    IF(NIC1.GT.0) WRITE(IOUT,900) IC1(1),IC1(NIC1)                    
    IF(NIC2.GT.0) WRITE(IOUT,901) IC2(1),IC2(NIC2)                    
 50 CONTINUE                                                          
    RETURN                                                            
    END