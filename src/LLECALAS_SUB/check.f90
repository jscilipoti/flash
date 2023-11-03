SUBROUTINE CHECK(N,YVAL)                                          
    IMPLICIT REAL*8(A-H,O-Z)                                          
    COMMON/CGIBBS/NF,MAXZ,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AL(10),DA(10,10),XM(10,4)                                       
    COMMON/CIPR/IPR                                                   
    COMMON/COUT/IOUT                                                  
    DIMENSION YVAL(30)                                                
    JMAX=NF                                                           
    SMAX=SFAS(NF)                                                     
    NT=NF-1                                                           
    DO 5 J=1,NT                                                       
    IF(SFAS(J).LT.SMAX) GOTO 5                                        
    SMAX=SFAS(J)                                                      
    JMAX=J                                                            
  5 CONTINUE                                                          
    IF(JMAX.EQ.NF) GOTO 20                                            
    DO 10 I=1,N                                                       
    DS=1.                                                             
    DO 15 J=1,NT                                                      
 15 DS=DS-YVAL(I+(J-1)*N)                                             
 10 YVAL(I+(JMAX-1)*N)=DS                                             
 20 IF(NF.EQ.2) GOTO 100                                              
    DO 21 I=1,N                                                       
    XVL(I,NF)=1.                                                      
    DO 21 J=1,NT                                                      
    XVL(I,J)=YVAL(I+(J-1)*N)                                          
 21 XVL(I,NF)=XVL(I,NF)-XVL(I,J)                                      
    DO 23 J=1,NF                                                      
    SFAS(J)=0.                                                        
    DO 22 I=1,N                                                       
    AL(I)=XVL(I,J)*Z(I)                                               
 22 SFAS(J)=SFAS(J)+AL(I)                                             
    DO 23 I=1,N                                                       
 23 XM(I,J)=AL(I)/SFAS(J)                                             
    DO 30 I=1,NT                                                      
    JN=I+1                                                            
    DO 30 J=JN,NF                                                     
    IF(DABS(XM(MAXZ,I)-XM(MAXZ,J)).LT.5.D-3) GOTO 40                  
 30 CONTINUE                                                          
    GOTO 100                                                          
 40 IV=I                                                              
    JV=J                                                              
    DMAX=0.                                                           
    DO 50 I=1,N                                                       
    R1=DABS(XM(I,IV)-XM(I,JV))                                        
    IF(R1.GT.DMAX) DMAX=R1                                            
 50 CONTINUE                                                          
    IF(DMAX.GT.2.5D-2) GOTO 100                                       
    WRITE(6,601)                                                      
    IF(IOUT.NE.6) WRITE(IOUT,601)                                     
    NF=NF-1                                                           
    NT=NT-1                                                           
    DO 60 I=1,N                                                       
    XVL(I,IV)=XVL(I,IV)+XVL(I,JV)                                     
 60 XVL(I,JV)=XVL(I,NF+1)                                             
100 CONTINUE                                                          
    IF(NF.LT.3) GOTO 250                                              
    MINF=0                                                            
    DO 200 I=1,NF                                                     
    IF(SFAS(I).LT.1.D-12) MINF=I                                      
200 CONTINUE                                                          
    IF(MINF.EQ.0) GOTO 250                                            
    WRITE(6,602) MINF,SFAS(MINF)                                      
    IF(IOUT.NE.6) WRITE(IOUT,602) MINF,SFAS(MINF)                     
    DO 220 I=1,NF                                                     
    IF(I.EQ.MINF) GOTO 220                                            
    DO 210 J=1,N                                                      
210 XVL(J,I)=XVL(J,I)+SFAS(I)*XVL(J,MINF)                             
220 CONTINUE                                                          
    IF(MINF.EQ.NF) GOTO 250                                           
    NNF=NF-1                                                          
    DO 230 I=1,NNF                                                    
    IF(I.LT.MINF) GOTO 230                                            
    DO 240 J=1,N                                                      
240 XVL(J,I)=XVL(J,I+1)                                               
230 CONTINUE                                                          
    NF=NF-1                                                           
250 CONTINUE                                                          
601 FORMAT(//,' * * * TWO PHASES ARE IDENTICAL * * *'//)              
602 FORMAT(/,' * THE AMOUNT OF PHASE',I2,' IS',D12.4,'. THE PHASE IS ELIMINATED *',/)                                                   
    RETURN                                                            
    END