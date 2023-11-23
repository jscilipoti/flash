SUBROUTINE TERM(Y,DMAT,ICOND,NEW)                                 
    IMPLICIT REAL*8(A-H,O-Z)                                          
    DIMENSION Y(4),DMAT(4,4),A(4)                                     
    IF(Y(1).LT.1.D-14.OR.Y(3).LT.1.D-14) GOTO 1                       
    IF(Y(1)+Y(2).GT.1.D0.OR.Y(3)+Y(4).GT.1.D0) GOTO 2                 
    IF(Y(1)+Y(2)-.01D0.LT.Y(3)+Y(4).AND.Y(1)-.01D0.LT.Y(3))GOTO 3     
    RETURN                                                            
1     ICOND=2                                                           
    DS=DMAT(1,1)/(DMAT(1,1)-Y(1))                                     
    DO 5 I=1,4                                                        
5     Y(I)=DMAT(I,1)+DS*(Y(I)-DMAT(I,1))                                
    Y(1)=0.D0                                                         
    NEW=1                                                             
    RETURN                                                            
2     ICOND=-2                                                          
    RETURN                                                            
3     ICOND=1                                                           
    ND=2+NEW                                                          
    IF(ND.GT.4)ND=ND-4                                                
    DO 6 I=1,4                                                        
6     A(I)=DMAT(NEW,I)-DMAT(ND,I)                                       
    DS=0.D0                                                           
    NITER=0                                                           
7     NITER=NITER+1                                                     
    IF(NITER.LE.10)GOTO 8                                             
    ICOND=-1                                                          
    RETURN                                                            
8     F=((A(4)*DS+A(3))*DS+A(2))*DS+A(1)                                
    DF=(3.D0*A(4)*DS+2.D0*A(3))*DS+A(2)                               
    DF=-F/DF                                                          
    DS=DS+DF                                                          
    IF(DABS(DF).GT.1.D-6)GOTO 7                                       
    DO 9 I=1,4                                                        
    Y(I)=DMAT(I,4)                                                    
    DO 9 J=1,3                                                        
9     Y(I)=Y(I)*DS+DMAT(I,4-J)                                          
    RETURN                                                            
    END