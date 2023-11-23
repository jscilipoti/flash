SUBROUTINE GMIX(NARG,NDIF,FUN,GRAD,XMAT,YVAL)                     
    IMPLICIT REAL*8(A-H,O-Z)                                          
    COMMON/CVAP/NOVAP,NDUM,IDUM(4),PRAT(10)                           
    COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                     
    COMMON/CGIBBS/NF,MAXZ,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),XX(10),DA(10,10),XM(10,4)                                       
    DIMENSION YVAL(30),GRAD(30),X(10),TM(10),FG(10),XMAT(30,30)       
    NG=NF-1                                                           
    N=NARG/NG                                                         
    JD=1                                                              
    IF(NDIF.EQ.2) JD=2                                                
    IF(NDIF.EQ.2) CALL CHECK(N,YVAL)                                  
    IF(NF.NE.NG) GOTO 20                                              
    NARG=NARG-N                                                       
    NG=NG-1                                                           
    DO 10 I=1,N                                                       
    DO 10 J=1,NG                                                      
 10 YVAL(I+(J-1)*N)=XVL(I,J)                                          
 20 FUN=-GNUL                                                         
    DO 50 I=1,N                                                       
    XVL(I,NF)=1.                                                      
    DO 30 J=1,NG                                                      
    XVL(I,J)=YVAL(I+(J-1)*N)                                          
 30 XVL(I,NF)=XVL(I,NF)-XVL(I,J)                                      
    ! write(*,*)xvl
    DO 40 J=1,NF                                                      
    IF(XVL(I,J).GT.0.) GOTO 40                                        
    FUN=0.                                                            
    GOTO 1000                                                         
 40 CONTINUE                                                          
 50 CONTINUE                                                          
    DO 200 J=1,NF                                                     
    SFAS(J)=0.                                                        
    DO 60 I=1,N                                                       
    X(I)=XVL(I,J)*Z(I)                                                
 60 SFAS(J)=SFAS(J)+X(I)                                              
    DO 65 I=1,N                                                       
    XX(I)=X(I)/SFAS(J)                                                
 65 XM(I,J)=XX(I)                                                     
    CALL unifac(JD,XX,FG,DA,PACT)                                     
    IDUM(J)=NDUM                                                      
    DO 70 I=1,N                                                       
    TM(I)=DLOG(XVL(I,J)/SFAS(J))+FG(I)                                
 70 FUN=FUN+X(I)*TM(I)                                                
    IF(NDIF.EQ.0) GOTO 200                                            
    DO 80 I=1,N                                                       
    S=Z(I)*TM(I)                                                      
    IF(J.EQ.NF) GOTO 75                                               
    GRAD(I+(J-1)*N)=S                                                 
    GOTO 80                                                           
 75 DO 76 K=1,NG                                                      
    NK=I+(K-1)*N                                                      
 76 GRAD(NK)=GRAD(NK)-S                                               
 80 CONTINUE                                                          
    IF(NDIF.EQ.1) GOTO 200                                            
    DO 100 I=1,N                                                      
    ST=Z(I)/SFAS(J)                                                   
    DO 100 L=1,N                                                      
    S=ST*(DA(I,L)-1.)*Z(L)                                            
    IF(L.EQ.I)S=S+Z(I)/XVL(I,J)                                       
    IF(J.EQ.NF) GOTO 90                                               
    XMAT(I+(J-1)*N,L+(J-1)*N)=S                                       
    GOTO 95                                                           
 90 DO 92 K=1,NG                                                      
    DO 92 M=1,K                                                       
    NK=I+(K-1)*N                                                      
    NM=L+(M-1)*N                                                      
    IF(K.NE.M) XMAT(NK,NM)=S                                          
    IF(K.EQ.M) XMAT(NK,NM)=XMAT(NK,NM)+S                              
 92 CONTINUE                                                          
 95 CONTINUE                                                          
100 CONTINUE                                                          
200 CONTINUE                                                          
1000 RETURN                                                            
    END