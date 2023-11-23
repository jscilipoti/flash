SUBROUTINE CHOL(N,A)                                              
    IMPLICIT REAL*8(A-H,O-Z)                                          
    DIMENSION A(2,2)                                                  
    DO 50 I=1,N                                                       
    I1=I-1                                                            
    IF(I1.EQ.0) GOTO 30                                               
    DO 20 J=I,N                                                       
    DO 20 K=1,I1                                                      
 20 A(I,J)=A(I,J)-A(I,K)*A(J,K)                                       
 30 IF(A(I,I).LT.1.D-14) A(I,I)=1.D-14                                
    A(I,I)=DSQRT(A(I,I))                                              
    IF(I.EQ.N) GOTO 100                                               
    J1=I+1                                                            
    DO 50 J=J1,N                                                      
 50 A(J,I)=A(I,J)/A(I,I)                                              
100 RETURN                                                            
    END 