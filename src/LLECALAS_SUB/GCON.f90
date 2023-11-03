SUBROUTINE GCON(NK,X,ACT,DACT,ICVEX)                              
    IMPLICIT REAL*8(A-H,O-Z)                                          
    DIMENSION X(3),DG(2),DDG(2,2),ACT(3),DACT(10,10)                  
    ICVEX=1                                                           
    DO 1 I=1,NK                                                       
  1 IF(ACT(I).LT.1.D-10) ACT(I)=1.D-10                                
    DO 5 I=1,NK                                                       
    DO 5 J=1,NK                                                       
  5 DACT(I,J)=DACT(I,J)/ACT(I)                                        
    IF(NK.EQ.3) GOTO 9                                                
    DDG(2,2)=DACT(2,2)-DACT(1,2)-DACT(2,1)+DACT(1,1)                  
    GOTO 30                                                           
9     DO 20 I=2,NK                                                      
    II=I-1                                                            
    DO 20 J=2,NK                                                      
    JJ=J-1                                                            
 20 DDG(II,JJ)=DACT(I,J)-DACT(1,J)-DACT(I,1)+DACT(1,1)                
    IF(X(1).LE.1.D-12.OR.X(2).LE.1.D-12) GOTO 30                      
    DET=DDG(1,1)*DDG(2,2)-DDG(2,1)*DDG(2,1)                           
    IF(DET.LE.0.D0.OR.DDG(1,1).LE.0.D0.OR.DDG(2,2).LE.0.D0) ICVEX=-1  
    GOTO 100                                                          
 30 CONTINUE                                                          
    IF(DDG(2,2).LE.0.D0) ICVEX=-1                                     
100 CONTINUE                                                          
    RETURN                                                            
    END