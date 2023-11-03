SUBROUTINE PARAM2                                                 
    IMPLICIT REAL*8(A-H,O-Z)                                          
    COMMON/CUFAC/NK,NG,P(10,10),T                                     
    COMMON/CPAR/TAU(10,10),S(10,10),F(10)                             
    COMMON/CQT/QT(10,10),Q(10),R(10)                                  
    DO 30 I=1,NG                                                      
      DO 30 J=1,NG                                                      
 30       TAU(I,J)=DEXP(-P(I,J)/T)                                          
 40 CONTINUE                                                          
    DO 50 I=1,NK                                                      
      DO 50 K=1,NG                                                      
          S(K,I)=0.D0                                                       
    DO 50 M=1,NG                                                      
 50   S(K,I)=S(K,I)+QT(M,I)*TAU(M,K)                                    
    DO 60 I=1,NK                                                      
      F(I)=1.D0                                                         
    DO 60 J=1,NG                                                      
 60   F(I)=F(I)+QT(J,I)*DLOG(S(J,I))                                    
    RETURN                                                            
    END