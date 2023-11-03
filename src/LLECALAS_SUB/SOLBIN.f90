SUBROUTINE SOLBIN                                                 
    IMPLICIT REAL*8(A-H,O-Z)                                          
    COMMON/CACT/X1(10),X2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
    COMMON/CY/Y13,Y21,STEP                                            
    COMMON/COUT/IOUT                                                  
    DIMENSION DMAT(2,3)                                               
    common/nga/nga,mass(12)

    NITER=0                                                           
    X1(1)=1.D0-Y13/100.D0                                             
    X2(1)=Y21/100.D0                                                  
 10 NITER=NITER+1                                                     
    IF(NITER.GT.10) GOTO 50                                           
    IF(X1(1).LT.0.D0) X1(1)=0.D0                                      
    IF(X2(1).LT.0.D0) X2(1)=0.D0                                      
    X1(2)=1.D0-X1(1)                                                  
    X2(2)=1.D0-X2(1)                                                  
    CALL unifac(3,X1,ACT1,DACT1,PACT)                                 
    CALL unifac(3,X2,ACT2,DACT2,PACT)                                 
    DO 20 I=1,2                                                       
    DMAT(I,1)=DACT1(I,1)-DACT1(I,2)                                   
    DMAT(I,2)=DACT2(I,2)-DACT2(I,1)                                   
 20 DMAT(I,3)=ACT1(I)-ACT2(I)                                         
    CALL GAUSL(2,3,2,1,DMAT)                                          
    RES=DMAT(1,3)**2+DMAT(2,3)**2                                     
    X1(1)=X1(1)-DMAT(1,3)                                             
    X2(1)=X2(1)-DMAT(2,3)                                             
    IF(RES.GT.1.D-20) GOTO 10                                         
 50 CONTINUE                                                          
    WRITE(6,603)                                                      
    IF(IOUT.NE.6) WRITE(IOUT,603)                                     
    IF(IOUT.NE.6) WRITE(IOUT,604) X1(1),X2(1),X1(2),X2(2)             
    WRITE(6,604) X1(1),X2(1),X1(2),X2(2)                              
603 FORMAT(///,5X,'** BINARY SOLUBILITIES IN MOLE FRACTIONS **',//,11X,'COMPONENT 1',15X,'COMPONENT 2',/)                               
604 FORMAT(2(2X,2P2D12.2)//)                                          
    CALL GCON(2,X1,ACT1,DACT1,ICVEX)                                  
    IF(IOUT.NE.6.AND.ICVEX.EQ.-1) WRITE(IOUT,601)                     
    IF(ICVEX.EQ.-1) WRITE(6,601)                                      
    CALL GCON(2,X2,ACT2,DACT2,ICVEX)                                  
    IF(IOUT.NE.6.AND.ICVEX.EQ.-1) WRITE(IOUT,602)                     
    IF(ICVEX.EQ.-1) WRITE(6,602)                                      
601 FORMAT(' FALSE SOLUTION IN PHASE 1')                              
602 FORMAT(' FALSE SOLUTION IN PHASE 2')                              
    RETURN                                                            
    END