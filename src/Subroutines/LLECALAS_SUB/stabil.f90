SUBROUTINE STABIL(N,NDIF,FUN,GRAD,XMAT,Y)                         
    IMPLICIT REAL*8(A-H,O-Z)                                          
    COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
    COMMON/CGIBBS/NF,MAXZ,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AL(10),DA(10,10),XM(10,4)                                       
    DIMENSION GRAD(30),XMAT(30,30),Y(30),YEX(10),XEX(10)              
    common/nga/nga,mass(12)

    SUM=0.                                                            
    DO 10 I=1,N                                                       
    YEX(I)=0.                                                         
    IF(Y(I).GT.-40.) YEX(I)=DEXP(Y(I))                                
 10 SUM=SUM+YEX(I)                                                    
    DO 15 I=1,N                                                       
 15 XEX(I)=YEX(I)/SUM                                                 
    JD=1                                                              
    IF(NDIF.EQ.2) JD=2                                                
    CALL unifac(JD,XEX,AL,DA,PACT)                                    
    FUN=1.                                                            
    DO 20 I=1,N                                                       
    S=Y(I)+AL(I)-A(I)                                                 
    IF(NDIF.EQ.0) GOTO 20                                             
    GRAD(I)=YEX(I)*S                                                  
 20 FUN=FUN+YEX(I)*(S-1)                                              
    IF(NDIF.LT.2) GOTO 50                                             
    DO 30 I=1,N                                                       
    S=XEX(I)                                                          
    DO 40 J=1,I                                                       
    XMAT(I,J)=S*YEX(J)*DA(I,J)                                        
 40 XMAT(J,I)=XMAT(I,J)                                               
 30 XMAT(I,I)=XMAT(I,I)+YEX(I)+GRAD(I)                                



 50 RETURN                                                            
    END