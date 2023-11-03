FUNCTION RUND(Y,DY,S,IRUND)                                       
    IMPLICIT REAL*8(A-H,O-Z)                                          
    X=Y+S*DY+1.D-8*DY**2                                              
    IX=IRUND*X                                                        
    Z=DFLOAT(IX)/IRUND-Y                                              
    RUND=DABS(Z)                                                      
    RETURN                                                            
    END