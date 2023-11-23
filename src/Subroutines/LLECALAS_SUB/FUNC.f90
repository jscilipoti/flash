SUBROUTINE FUNC(N,M,NDIF,X,SSQ)                                   
    IMPLICIT REAL*8(A-H,O-Z)                                          
    COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
    COMMON/CMARQ/GRAD(2),XJTJ(2,2)                                    
    COMMON/CUFAC/NK,NG,P(10,10),T                                     
    COMMON/CA/XC(5),GE(5,2),GC(5,2)                                   
    DIMENSION X(2),F(2)                                               
    JD=0                                                              
    IF(NDIF.EQ.1) JD=4                                                
    P(1,2)=X(1)*300.D0                                                
    P(2,1)=X(2)*300.D0                                                
    CALL PARAM2                                                       
    SSQ=0.                                                            
    IF(NDIF.EQ.0) GOTO 11                                             
    DO 10 I=1,2                                                       
    GRAD(I)=0.                                                        
    DO 10 J=1,2                                                       
 10 XJTJ(I,J)=0.                                                      
 11 CONTINUE                                                          
    DO 21 L=1,5                                                       
    Y1(1)=XC(L)                                                       
    Y1(2)=1.D0-XC(L)                                                  
    CALL unifac(JD,Y1,ACT1,DACT1,PACT)                                
    DO 17 I=1,2                                                       
    F(I)=ACT1(I)-GE(L,I)                                              
    GC(L,I)=ACT1(I)                                                   
 17 SSQ=SSQ+F(I)*F(I)                                                 
    IF(JD.EQ.0) GOTO 21                                               
    DO 19 I=1,2                                                       
    DO 20 J=1,2                                                       
    GRAD(J)=GRAD(J)+F(I)*PACT(I,J)                                    
    DO 20 K=1,2                                                       
 20 XJTJ(J,K)=XJTJ(J,K)+PACT(I,J)*PACT(I,K)                           
 19 CONTINUE                                                          
 21 CONTINUE                                                          
    RETURN                                                            
    END