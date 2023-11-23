SUBROUTINE MARQ(FUNC,N,M,X,XLAMB,FAC,EPSG,MAXF)                   
    IMPLICIT REAL*8(A-H,O-Z)                                          
    COMMON/COUT/IOUT                                                  
    COMMON/CMARQ/GRAD(2),XJTJ(2,2)                                    
    COMMON/CIPR/IPR                                                   
    DIMENSION X(2),Y(2),XNY(2),A(2,2),DX(2)                           
    IEVAL=0                                                           
    ISTOP=0                                                           
    IER=0                                                             
    XMAXL=XLAMB*1.D+4                                                 
    ITER=0                                                            
    IF(IOUT.NE.6.AND.IPR.EQ.1) WRITE(IOUT,603)                        
    IF(IPR.EQ.1) WRITE(6,603)                                         
    CALL FUNC(N,M,1,X,SRES)                                           
    IEVAL=IEVAL+1                                                     
    SSQ=SRES                                                          
    IF(IPR.EQ.1.AND.IOUT.NE.6) WRITE(IOUT,601) ITER,SSQ               
    IF(IPR.EQ.1) WRITE(6,601) ITER,SSQ                                
 10 CONTINUE                                                          
    IF(IEVAL.NE.1) CALL FUNC(N,M,1,X,SRES)                            
    GNORM=0.D0                                                        
    DO 140 I=1,N                                                      
140 GNORM=GNORM+GRAD(I)**2                                            
    GNORM=DSQRT(GNORM)                                                
    IF(GNORM.LT.EPSG) ISTOP=ISTOP+1                                   
    IF(ISTOP.GT.0) GOTO 1000                                          
    ITER=ITER+1                                                       
 49 CONTINUE                                                          
    IF(IEVAL.GT.MAXF) GOTO 998                                        
    DO 41 I=1,N                                                       
    DO 40 J=1,N                                                       
 40 A(I,J)=XJTJ(I,J)                                                  
 41 A(I,I)=A(I,I)+XLAMB                                               
    CALL CHOL(N,A)                                                    
    Y(1)=-GRAD(1)/A(1,1)                                              
    DO 81 I=2,N                                                       
    SUM=0.D0                                                          
    II=I-1                                                            
    DO 80 J=1,II                                                      
 80 SUM=SUM+A(I,J)*Y(J)                                               
 81 Y(I)=(-GRAD(I)-SUM)/A(I,I)                                        
    DX(N)=Y(N)/A(N,N)                                                 
    DO 85 I=2,N                                                       
    II=N-I+1                                                          
    SUM=0.D0                                                          
    III=II+1                                                          
    DO 84 J=III,N                                                     
 84 SUM=SUM+A(J,II)*DX(J)                                             
 85 DX(II)=(Y(II)-SUM)/A(II,II)                                       
    DO 90 I=1,N                                                       
 90 XNY(I)=X(I)+DX(I)                                                 
    CALL FUNC(N,M,0,XNY,SRES)                                         
    IEVAL=IEVAL+1                                                     
    SSQNY=SRES                                                        
    SQ1=0.D0                                                          
    SQ2=0.D0                                                          
    DO 110 I=1,N                                                      
    Y(I)=XNY(I)*300.D0                                                
    SQ1=SQ1+DX(I)**2                                                  
    SQ2=SQ2-DX(I)*GRAD(I)                                             
110 CONTINUE                                                          
    CCAL=SSQ-SSQNY                                                    
    CPRE=SQ2+XLAMB*SQ1                                                
    CALPRE=CCAL/(CPRE+1.D-14)                                         
    IF(IPR.EQ.1) WRITE(6,601) ITER,SSQNY,Y(1),Y(2),GNORM,XLAMB,CALPRE 
    IF(IOUT.NE.6.AND.IPR.EQ.1) WRITE(IOUT,601) ITER,SSQNY,Y(1),Y(2),  GNORM,XLAMB,CALPRE                                                
    IF(SSQ-SSQNY.GT..75*CPRE) XLAMB=XLAMB/FAC                         
    IF(SSQ-SSQNY.LT..25*CPRE) XLAMB=XLAMB*FAC                         
    IF(XLAMB.GT.XMAXL) GOTO 999                                       
    IF(SSQNY-SSQ) 120,120,119                                         
119 CONTINUE                                                          
    IF(SSQNY-SSQ.GT.DABS(SSQ)*.5D0) XLAMB=XLAMB*FAC                   
    GOTO 49                                                           
120 CONTINUE                                                          
    IF(SSQ-SSQNY.GT.DABS(SSQ)*.8D0) XLAMB=XLAMB/FAC                   
    DO 130 I=1,N                                                      
130 X(I)=XNY(I)                                                       
    SSQ=SSQNY                                                         
    GOTO 10                                                           
998 IER=1                                                             
    GOTO 1000                                                         
999 IER=2                                                             
    GOTO 1000                                                         
601 FORMAT(1X,I3,2X,D12.4,2F10.3,5X,2D12.4,3X,F10.2)                  
602 FORMAT(//,' ISTOP= ',I2,5X,'IER = ',I2,5X,'IEVAL = ',I5,//)       
603 FORMAT(///,'  ** ITERATIONCOURSE, UNIQUAC-PARAMETERS FROM UNIFAC **'/)                                                              
1000 CONTINUE                                                          
    IF(IOUT.NE.6.AND.IPR.EQ.1) WRITE(IOUT,602) ISTOP,IER,IEVAL        
    IF(IPR.EQ.1) WRITE(6,602)ISTOP, IER, IEVAL                        
    RETURN                                                            
    END