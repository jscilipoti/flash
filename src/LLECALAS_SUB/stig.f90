SUBROUTINE STIG(Y,S)                                              
    IMPLICIT REAL*8(A-H,O-Z)                                          
    COMMON/CVAP/NOVAP,NDUM,IDUM(4),PRAT(10)                           
    COMMON/CUFAC/N,NG,P(10,10),T                                      
    COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
    COMMON/CGIBBS/NF,MAXZ,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AA(10),DA(10,10),XM(10,4)                                       
    DIMENSION Y(10),V(10),YGEM(10)                                    
    common/nga/nga,mass(12)
    JPUR=0                                                            
    AMAX=0.                                                           
    DO 10 I=1,N                                                       
    IF(A(I).LT.AMAX) GOTO 10                                          
    JPUR=I                                                            
    AMAX=A(I)                                                         
 10 CONTINUE                                                          
    RMAX=1.D5                                                         
    NGN=N                                                             
    IF(NF.GT.1) NGN=N+NF                                              
    NEG=0                                                             
    DO 100 KK=1,NGN                                                   
    JM=KK                                                             
    IF(JPUR.NE.0) JM=JPUR                                             
    IF(JM.LE.N) GOTO 30                                               
    DO 20 I=1,N                                                       
 20 Y(I)=Z(I)*(2+XVL(I,JM-N)/SFAS(JM-N))/3                            
    GOTO 40                                                           
 30 SUM=0.                                                            
    DO 35 I=1,N                                                       
    GG=A(I)-GAM(JM,I)                                                 
    IF(GG.LT.-50.D0) GG=-50.D0                                        
    Y(I)=DEXP(GG)                                                     
 35 SUM=SUM+Y(I)                                                      
 40 NA=3                                                              
    DO 43 K=1,NA                                                      
    DO 36 I=1,N                                                       
 36 Y(I)=Y(I)/SUM                                                     
    CALL unifac(1,Y,AA,DA,PACT)                                       
    IF(K.EQ.NA) GOTO 44                                               
    DO 41 I=1,N                                                       
 41 Y(I)=DEXP(A(I)-AA(I))                                             
 42 SUM=0.                                                            
    DO 43 I=1,N                                                       
 43 SUM=SUM+Y(I)                                                      
 44 CONTINUE                                                          
    YV1=0.                                                            
    DO 50 J=1,NF                                                      
 50 V(J)=0.                                                           
    DO 60 I=1,N                                                       
    GD=DLOG(Y(I))+AA(I)-A(I)                                          
    YV1=YV1+Y(I)*GD                                                   
    DO 60 J=1,NF                                                      
    K=J                                                               
    VV=XVL(I,K)*Z(I)/SFAS(K)                                          
    D=GD*(Y(I)-VV)                                                    
 60 V(J)=V(J)+D                                                       
    YV2=V(1)                                                          
    DO 70 J=1,NF                                                      
    IF(V(J).LT.YV2) YV2=V(J)                                          
 70 CONTINUE                                                          
    RT1=YV1                                                           
    IF(YV2.GT.0.) RT1=RT1-YV2/2                                       
    IF(NEG.EQ.0.AND.YV1.GT.0.) GOTO 80                                
    RT1=YV1                                                           
    IF(NEG.EQ.0) RMAX=0.                                              
    NEG=1                                                             
 80 IF(RT1.GT.RMAX) GOTO 100                                          
    S=YV1                                                             
    RMAX=RT1                                                          
    CC=DEXP(-YV1)                                                     
    DO 90 I=1,N                                                       
 90 YGEM(I)=Y(I)*CC                                                   
    IF(JPUR.NE.0) GOTO 110                                            
100 CONTINUE                                                          
110 DO 120 I=1,N                                                      
120 Y(I)=YGEM(I)                                                      
    RETURN                                                            
    END