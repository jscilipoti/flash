SUBROUTINE SPLIT(ND,N,IDI,BETA,DEL,G,W)                           
    IMPLICIT REAL*8(A-H,O-Z)                                          
    DIMENSION G(ND,ND),W(ND,5)                                        
    IDI=0                                                             
    DO 20 I=1,N                                                       
    I1=I-1                                                            
    I2=I+1                                                            
    ID=0                                                              
    TV=G(I,I)                                                         
    SV=TV                                                             
    IF (I1.EQ.0) GO TO 25                                             
    DO 30 J=1,I1                                                      
 30 SV=SV-G(I,J)**2                                                   
 25 IF (SV.LT. DEL*DEL) ID=1                                          
    SVR=DSQRT(DABS(SV))                                               
    IF (SVR.LT. DEL) SVR=DEL                                          
    XM=0.                                                             
    IF (I.EQ.N) GO TO 35                                              
    DO 40 J=I2,N                                                      
    S=G(J,I)                                                          
    IF (I1.EQ.0) GO TO 45                                             
    DO 50 K=1,I1                                                      
 50 S=S-G(I,K)*G(J,K)                                                 
 45 S=S/SVR                                                           
    IF (DABS(S).GT. XM) XM=DABS(S)                                    
 40 G(J,I)=S                                                          
 35 IF (XM.LT. BETA) GO TO 55                                         
    ID=1                                                              
    XM=XM/BETA                                                        
    DO 60 J=I,N                                                       
 60 G(J,I)=G(J,I)/XM                                                  
    SVR=SVR*XM                                                        
 55 IF (ID.EQ.1) W(I,1)=SVR**2-SV                                     
    G(I,I)=SVR                                                        
    IDI=IDI+ID                                                        
 20 CONTINUE                                                          
    RETURN                                                            
    END                                                               
