SUBROUTINE TMSSJ(ND,N,IPR,NMAX,XLAM,GLIM,F   ,X,GD,G,W,ifunc)     
    IMPLICIT REAL*8(A-H,O-Z)                                          
    COMMON/COUT/IOUT                                                  
    integer::ifunc
    DIMENSION X(ND),GD(ND),G(ND,ND),W(ND,5)                           
    NTM=0                                                             
    DEL=1.D-7                                                         
150 NTM=NTM+1                                                         
    GNM=0.                                                            
    
    if(ifunc==1)then
      CALL STABIL(N,2,F,GD,G,X)                                       
    elseif(ifunc==2)then
      CALL GMIX(N,2,F,GD,G,X)  
    endif
    DO 5 I=1,N                                                        
  5 GNM=GNM+GD(I)**2                                                  
    IF (GNM.LT. GLIM) GO TO 200                                       
    BETA=0.                                                           
    DO 10 I=1,N                                                       
    W(I,4)=X(I)                                                       
    W(I,5)=GD(I)                                                      
    DO 10 J=1,I                                                       
    S=DABS(G(I,J))                                                    
    IF (I.EQ.J) S=S*N                                                 
    IF (S.GT. BETA) BETA=S                                            
 10 W(I,1)=0.                                                         
    BETA=DSQRT(BETA/N)                                                
    CALL SPLIT(ND,N,IDI,BETA,DEL,G,W)                                 
    XLM=0.                                                            
    NTS=0                                                             
350 NTS=NTS+1                                                         
    CALL LINE(ND,N,XLM,GD,G,W)                                        
    SMAX=0.                                                           
    GP=0.                                                             
    DP=0.                                                             
    DO 205 I=1,N                                                      
    S=W(I,3)                                                          
    IF (DABS(S).GT. SMAX) SMAX=DABS(S)                                
    GP=GP+S*GD(I)                                                     
205 DP=DP+S*S*W(I,1)                                                  
    FAC=1.                                                            
    IF (SMAX.GT. XLAM) FAC=XLAM/SMAX                                  
    DER=FAC*GP                                                        
    ALFA=1.                                                           
210 FF=ALFA*FAC                                                       
    DO 215 I=1,N                                                      
215 X(I)=W(I,4)+FF*W(I,3)                                             

     if(ifunc==1)then
      CALL STABIL(N,1,FNEW,GD,G,X)                           
    elseif(ifunc==2)then
      CALL GMIX(N,1,FNEW,GD,G,X)   
    endif



                                           
    IF(FNEW.NE.0.) GOTO 220                                           
    ALFA=ALFA/2.                                                      
    GOTO 210                                                          
220 DELS=FNEW-F                                                       
     IF (FNEW.NE. 0.  .AND. DELS.LT. 1.D-10) GO TO 230                
    IF (NTS.GT. 1) GO TO 125                                          
    DE2=(DELS-ALFA*DER)/ALFA**2/2                                     
    GT=-DER/DE2                                                       
    IF (GT.GT. ALFA) ALFA=ALFA/3.                                     
    IF (GT.LE. ALFA) ALFA=GT                                          
    GO TO 210                                                         
230 PRED=FF*GP-FF**2/2*(GP+DP)                                        
    CALPRE=DELS/PRED                                                  
    F=FNEW                                                            
    DO 345 I=1,N                                                      
345 W(I,4)=X(I)                                                       
    IF (NTS.GE.3) GO TO 125                                           
    IF (ALFA.EQ. 1.  .AND. CALPRE.GT. .8) GO TO 350                   
125 DO 126 I=1,N                                                      
126 X(I)=W(I,4)                                                       
    IF (IPR.NE.0) WRITE(6,130) NTM,F,CALPRE,GNM,ALFA                  
    IF(IOUT.NE.6.AND.IPR.NE.0) WRITE(IOUT,130) NTM,F,CALPRE,GNM,ALFA  
130 FORMAT(1X,I2,':  F = ',D16.8,'  CALPRE = ',F8.3,'  GRAD = ',D16.8,'  ALFA = ',D10.2)                                                
    IF (IPR.GT.1) WRITE(6,131) (X(I),I=1,N)                           
    IF(IOUT.NE.6.AND.IPR.GT.1) WRITE(IOUT,131) (X(I),I=1,N)           
131 FORMAT(' X-VECTOR',10F11.5,(/,9X,10F11.5))                        
    IF (IPR.GT. 1) WRITE(6,132)                                       
    IF(IOUT.NE.6.AND.IPR.GT.1) WRITE(IOUT,132)                        
132 FORMAT(' ')                                                       
    IF (NTM.LT. NMAX) GO TO 150                                       
200 IF(IPR.GT.0) WRITE(6,133) NTM,GNM,F                               
    IF(IOUT.NE.6.AND.IPR.GT.0) WRITE(IOUT,133) NTM,GNM,F              
133 FORMAT(/,'  NUMBER OF ITERATIONS = ',I2,', NORM OF GRADIENT = ',  D12.5,', F = ',D14.6)                                             
    IF(IPR.GT.0) WRITE(6,131) (X(I),I=1,N)                            
    IF(IOUT.NE.6.AND.IPR.GT.0) WRITE(IOUT,131) (X(I),I=1,N)           
    RETURN                                                            
    END