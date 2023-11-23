SUBROUTINE SOLVE(Y,DY,NOLD,NEW,NITER,N,NT)                        
    IMPLICIT REAL*8(A-H,O-Z)                                          
    COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
    COMMON/CBISO/NIC1,NIC2,IC1(120),IC2(120)                          
    DIMENSION Y(4),DY(4),AMAT(3,5)                                    
    DIMENSION INO(3)                                                  
    common/nga/nga,mass(12)
    NK=1000  !cambi� de 3 a 100 !Alfonsina                            
    NITER=0                                                           
    NT=1000  !cambi� de 3 a 100 !Alfonsina                            
    DO 1 I=1,3                                                        
  1 IF(Y(I).LT.1.D-14) NT=1000    !Cambie de 10 iteraciones a100 !Alfonsina 0
11    NITER=NITER+1                                                     
    IF(NITER.GT.NT) RETURN                                            
    DO 2 I=1,4                                                        
2     IF(Y(I).LT.0.D0)Y(I)=0.D0                                         
    DO 3 I=1,2                                                        
    Y1(I)=Y(I)                                                        
3     Y2(I)=Y(I+2)                                                      
    Y1(3)=1.D0-Y1(1)-Y1(2)                                            
    Y2(3)=1.D0-Y2(1)-Y2(2)                                            
    IF(Y1(3).LT.0.)Y1(3)=0.                                           
    IF(Y2(3).LT.0.) Y2(3)=0.                                          
    CALL unifac(3,Y1,ACT1,DACT1,PACT)                                 
    CALL unifac(3,Y2,ACT2,DACT2,PACT)                                 
    J=0                                                               
    DO 6 I=1,4                                                        
    IF(I.EQ.NOLD)GOTO 6                                               
    J=J+1                                                             
    INO(J)=I                                                          
6     CONTINUE                                                          
    DO 7 I=1,3                                                        
    DO 7 J=1,2                                                        
    AMAT(I,J)=DACT1(I,J)-DACT1(I,3)                                   
7     AMAT(I,J+2)=DACT2(I,3)-DACT2(I,J)                                 
    DO 8 I=1,3                                                        
    AMAT(I,5)=AMAT(I,NOLD)                                            
    DO 9 J=1,3                                                        
9     AMAT(I,J)=AMAT(I,INO(J))                                          
8     AMAT(I,4)=ACT1(I)-ACT2(I)                                         
    CALL GAUSL(3,5,3,2,AMAT)                                          
    RES=0.D0                                                          
    DO 10 I=1,3                                                       
    Y(INO(I))=Y(INO(I))-AMAT(I,4)                                     
    DY(INO(I))=-AMAT(I,5)                                             
10    RES=RES+AMAT(I,4)**2                                              
    IF(RES.GT.1.D-10)GOTO 11                                          
    IZ=0                                                              
    DO 14 I=1,2                                                       
    IF(Y1(I).LT.1.D-14) IZ=1                                          
 14 IF(Y2(I).LT.1.D-14) IZ=1                                          
    IF(IZ.EQ.1) GOTO 13                                               
    CALL GCON(3,Y1,ACT1,DACT1,ICVEX)                                  
    IF(ICVEX.EQ.1) GOTO 15                                            
    NIC1=NIC1+1                                                       
    IC1(NIC1)=N+1                                                     
 15 CALL GCON(3,Y2,ACT2,DACT2,ICVEX)                                  
    IF(ICVEX.EQ.1) GOTO 13                                            
    NIC2=NIC2+1                                                       
    IC2(NIC2)=N+1                                                     
13    DY(NOLD)=1.D0                                                     
    NEW=NOLD                                                          
    DYMAX=1.D0                                                        
    DO 12 I=1,4                                                       
    IF(DABS(DY(I)).LE.DYMAX)GOTO 12                                   
    NEW=I                                                             
    DYMAX=DABS(DY(I))                                                 
12    CONTINUE                                                          
    RETURN                                                            
    END