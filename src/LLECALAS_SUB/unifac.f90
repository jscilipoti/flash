SUBROUTINE unifac(NDIF,X,ACT,DACT,PACT)                           
    IMPLICIT REAL*8(A-H,O-Z)                                          
!c------
    common/asoc/nktt,igamt(20,12),nytt(20,12)   
    common/nga/nga,mass(12)
    common/grupas1/rkass(6,12,6,12),enass(6,12,6,12),deloh(6,12,6,12)!Alfonsin
!c------
    COMMON/CVAP/NOVAP,NDUM,IDUM(4),PRAT(10)                           
    COMMON/CUFAC/NK,NG,P(10,10),T                                     
    COMMON/CPAR/TAU(10,10),S(10,10),F(10)                             
    COMMON/CQT/QT(10,10),Q(10),R(10)                                  
    DIMENSION X(10),GAM(10),ACT(10),DACT(10,10),THETA(10),PHI(10),RI(10),&
    &QI(10),ETA(10),QIL(10),RIL(10),QID(10),ETAL(10),TETAR(10) 
    DIMENSION U(10,10),V(10,10),PACT(2,2),DTAU(2,2,2)                 
!c------
    dimension goh(10),xgamk(20),dxohdx(10),dxxdx(10,10),dasdx1(10,10),dasdx2(10,10),dasdx(10,10)
  common/ioh2/rngoh(12,12)

    dimension dif(12,12), dif1(10,12,12) !Alfonsina
    common/zzzas/xoh(6,12),xohi0(12,6,12),xoh_old(6,12),xohi(6,12),xohi_old(6,12), xohi0_old(12,6,12)  !Alfonsina
  dimension m_lambda(nga*2,nga*2),m_lambda1(nga*2,nga*2) !Alfonsina
  dimension psin(12) !Alfonsina
  dimension indx(20)
    double precision  m_lambda,m_lambda1,xoh,xohi0,xoh_old,xohi0_old  !Alfon
  double precision del, del1, dif, dif1, d1,psin, xgam, xnoh1 !Alfonsina
  integer order !Alfonsina
    double precision sum1, sum2, sum3, sum4, SUMA1J, sumaj !Alfonsina
    dimension xnohi0(12,12),tgt(12),dnohdx(12,12),actas(12) !Alfonsina
  dimension xnoh1(12), xnoh(12),das1(3),das3(3),dxkdni(12,6,12), dxkdnic(12,6,12), dgasdx(12)  !Alfonsina
    dimension dgasdxij (12,12), drhodx(12), drhodni(12,6,12)
  


    dk=1.381e-23
    deloh=0.0
    xnoh=0.0
    xnoh1=0.0
  xoh=0.0
    xgam=0.0
    do 7777 i=1,10
  xohi0=0
    xnohi0=0.0
    tgt(i)=0.0
7777 continue

    THETS=0.                                                          
    PHS=0.                                                            
    DO 10 I=1,NK                                                      
    THETA(I)=X(I)*Q(I)                                                
    PHI(I)=R(I)*X(I)                                                  
    THETS=THETS+THETA(I)                                              
 10 PHS=PHS+PHI(I)                                                    
    DO 20 I=1,NK                                                      
    RI(I)=R(I)/PHS                                                    
    RIL(I)=DLOG(RI(I))                                                
    QI(I)=Q(I)/THETS                                                  
 20 QIL(I)=DLOG(QI(I))                                                

    do 33 i=1,nk
    goh(i)=0.
    tgt(i)=0.0
    xnohi0=0.0
    xgam=0.0
   
!CCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    do j=1,nk
      tgt(j)=0.0
      xgam=0.0
  end do
 33 continue
    do i=1,nk
  if(nga.gt.0) then
  
      do k=1,nktt
          tgt(i)=tgt(i)+nytt(k,i)
      end do

      do j=1,nga  
          xnohi0(i,j)=rngoh(i,j)/R(i)  
      end do
      xgam=xgam+R(i)*x(i)
      end if  
  end do

    xnoh1=0d0
  do ja=1,nga
      do i=1,nk
      xnoh1(ja)=xnoh1(ja)+rngoh(i,ja)*x(i)
    end do
    end do
  
    do ja=1,nga
      xnoh(ja)=xnoh1(ja)/xgam
  end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!c------
    DO 40 I=1,NG                                                      
    ETA(I)=0.                                                         
    DO 45 J=1,NK                                                      
 45 ETA(I)=ETA(I)+S(I,J)*X(J)                                         
 40 ETAL(I)=DLOG(ETA(I))                                              
    DO 55 I=1,NG                                                      
    TETAR(I)=0.                                                       
    DO 55 J=1,NK                                                      
 55 TETAR(I)=TETAR(I)+QT(I,J)*X(J)                                    
    DO 60 I=1,NK                                                      
    QID(I)=1.-RI(I)/QI(I)                                             
    XX=F(I)+Q(I)*(1.-QIL(I))-RI(I)+RIL(I)                             
    XX=XX-5.*Q(I)*(QID(I)+RIL(I)-QIL(I))                              
    ACT(I)=XX                                                         
    DO 661 J=1,NG                                                     
    U(J,I)=S(J,I)/ETA(J)                                              
    V(J,I)=U(J,I)*TETAR(J)                                            
661 ACT(I)=ACT(I)-V(J,I)-QT(J,I)*ETAL(J)                              
!c------


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !**************************calculo de las fuerzas de asociacion*******************
    if(nga.ne.0) then
      DO J=1,NGA 
            IF(MASS(J).EQ.0) GO TO 201
              DO m=1,NGA
                IF(MASS(m).EQ.0) GO TO 101
                      DO L=1,MASS(J)
                              DO K=1,MASS(m)
                                IF(ENASS(K,m,L,J).EQ.0) THEN
                                       CONTINUE
                                  ELSE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    deloh(k,m,l,j)=(DEXP(ENASS(K,m,L,J)/T) - 1 )*RKASS(K,m,L,J)

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                                END IF
                              END DO
                      END DO

101            CONTINUE
              END DO
201          CONTINUE
      END DO


  end if  
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  !***********************************calculo de xoh**************************
!c C�lculo de la  fracci�n no asociada Paper:Ind. Eng. Chem. Res, 2004,43,203-208   
!c !Inicializaci�
    if(nga.ne.0) then
    !VALORES ORIGINALES
    !xoh=1.0d0 
    !del=1.d0
    !VALORES MODIFICADOS PORQUE SINO NO SE CUMPLIA EL DO WHILE DEBAJO Y QUEDABA
    !TRABADO ACÁ
    xoh=0.d0 
    del=0.d0
!c Iteraciones con tolerancia de 10-9
    do while (del>1.0d-10)
    xoh_old=xoh
  do m=1, nga
      do j=1,2
    sum1=0.D0
  do k=1, nga
  sum2=0.D0
  do l=1,2

  sum2=sum2+xoh_old(l,k)*deloh(l,k,j,m)
  end do
  sum1=sum1+sum2*xnoh1(k)
    end do
    xoh(j,m)=1.D0/(1.D0+sum1/xgam)          
  dif(j,m)=dabs((xoh(j,m)-xoh_old(j,m))/xoh(j,m))
  end do
  end do
  del=maxval(dif)
    end do
    end if

          write(4,*)"T=", t           
  write(4,*)"xoh (1,1)=", xoh (1,1)  
  write(4,*)"xoh (1,2)=", xoh(2,1) 

   

!cc Fin del C�lculo de la  fracci�n no asociada 
!c	!*****************************calculo de xohi0**************************************
!C	xohi0=1d0
!c	do i=1, nc
!C		do j=1, nga
!C	       do l=1,mass(j)
!C	do k=1,nga
!C	 do m=1,mass(k)
!C			If (rngoh(i,j).eq.0d0) then
!C					xohi0(i,l,j)=1.d0
!C			elseif (deloh(l,j,m,k).gt.0d0.and.rngoh(i,j).ne.0d0.and.
!C     @		mass(j).eq.2)then
!C		xohi0(i,l,j)=(-1d0+dsqrt(1d0+4d0*xnohi0(i,j)*deloh(l,j,m,k)))/
!C     @			(2d0*xnohi0(i,j)*deloh(l,j,m,k))
!C	
!C			end if
!C		 end do
!C	end do
!C	end do
!C		end do
!c	end do
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!c****************************calculo de xohi0(esta implementaci�n permite que una mol�cula
!c tenga m�s de un grupo asociativo 14/07/06**************************************
!c !Inicializaci�
   
    if(nga.ne.0) then
    xohi0=1.D0
    del1=1.D0
!c Iteraciones con tolerancia de 10-12
    dif1=0.D0 !ESTE VALOR LO INICIALIZÉ YO PORQUE ORIGINALMENTE SE INICIALIZABA
    !SOLO COMO 1E+310 Y JAMÁS SALIA DEL BUCLE SIGUIENTE
    do while (del1>1.0d-10)
    xohi0_old=xohi0
  do m=1, nga
  if	(rngoh(i,m).gt.0d0) then
      do j=1,2
    sum3=0.D0
  do k=1, nga
  sum4=0.D0
  do l=1,2 
  sum4=sum4+ xohi0_old(i,l,k)*deloh(l,k,j,m)*xnohi0(i,k)
  end do
  sum3=sum3+sum4
    end do
    xohi0(i,j,m)=1.D0/(1.D0+sum3)    
   dif1(i,j,m)=dabs((xohi0(i,j,m)-xohi0_old(i,j,m))/xohi0(i,j,m))
  end do
  else
  end if
  end do
  del1=maxval(dif1)
    end do
    end if

!c*****************************fin del calculo de xohi0**************************************
!C�lculo del gama de asociaci�n ALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C actas(M) = LOGARITMO NATURAL DEL GAMA DE ASOCIACI�N DEL COMPONENTE I  

     if(nga.ne.0) then	

    SUMAJ = 0.D0
    DO J=1,NGA 
    IF(MASS(J).NE.0) THEN      
    DO K=1,MASS(J)
  If(XOH(K,J).gt.1d-13)then
    SUMAJ = SUMAJ + RNGOH(i,j)*(dlog(XOH(K,J)/XOHI0(I,K,J))+0.5D0*(XOHi0(i,K,J)-1))+0.5D0*R(i)*xnoh(j)*(1-xoh(k,j))
  end if
    END DO
    ELSE
    CONTINUE
    END IF
    END DO
    actas(I) = SUMAJ
  end if
    act(i)=act(i)+actas(i)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 60 continue
    NDUM=0                                                            
    IF(NOVAP.EQ.0) GOTO 69                                            
    SS=0                                                              
    DO 61 I=1,NK                                                      
 61 SS=SS+X(I)*(PRAT(I)-ACT(I))                                       
    IF(SS.GT.0.) GOTO 69                                              
    NDUM=1                                                            
    DO 62 I=1,NK                                                      
    ACT(I)=PRAT(I)                                                    
    DO 62 J=1,NK                                                      
 62 DACT(I,J)=0.                                                      
    GOTO 100                                                          
 69 CONTINUE                                                          
    IF(NDIF.EQ.4) GOTO 90                                             
    IF(NDIF.LT.2) GOTO 100                                            
    DO 70 I=1,NK                                                      
    DO 70 J=I,NK                                                      
    XX=Q(I)*QI(J)*(1.-5.*QID(I)*QID(J))+(1.-RI(I))*(1.-RI(J))         
    DO 75 K=1,NG                                                      
 75 XX=XX+U(K,I)*(V(K,J)-QT(K,J))-U(K,J)*QT(K,I)                      

!********************************calculo de dxkdni Alfonsina**************************************

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!cCalcula los elementos de la matriz deltapq para el c�lculo de la derivada de la fracci�n 
!c no asociada respecto a la fracci�n molar del componente
    psin=0.0d0
  if(nga.ne.0) then
  m_lambda1=0.0d0
  m_lambda=0.0d0
    z=0; y=0
  do n=1,2
  do m=1,nga
       z=z+1
  do l=1, 2
  do k=1, nga
      y=y+1
    m_lambda(z,y)=xnoh(k)*deloh(l,k,n,m)*xoh(n,m)**2 
     if (z.eq.y)  then
    m_lambda(z,y)=m_lambda(z,y)+ 1.0d0
    end if
  end do
  end do
  y=0
  end do
  end do
    order=nga*2
  end if 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINA CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!c Calculo de  los elementos de la matriz [Yp] para el c�lculo de la derivada de la fracci�n 
!c no asociada respecto a la fracci�n molar del componente
    if(nga.ne.0) then
    do k=1,nga

  do ll=1,2
  do m=1,nga
  drhodni(j,ll,m)=(((rngoh(j,m)*xgam-xnoh1(m)*R(j)))/xgam**2)	  
  end do
  end do

    z=0     
  do ll=1,2
  do m=1,nga
  sum3=0.0d0
    do l=1,nga
  sum4=0.0d0


  do kk=1,2
  sum4=sum4+ (xoh(kk,l)*deloh(kk,l,ll,m))*drhodni(j,kk,l)
  end do
  sum3=sum3+sum4
  end do
  z=z+1
  psin(z)=-(xoh(ll,m)**2)*sum3
  end do
  end do


    N=order
  NP=order
    m_lambda1=m_lambda
    call  ludcmp(m_lambda1,N,NP,indx,d1)
     call lubksb(m_lambda1,N,NP,indx,psin)
!c colectando las derivadas en su correspondiente sub�ndice
    z=0
  do m=1,2
  do l=1, nga
  z= z+1
  dxkdni(k,m,l)=psin(z)
  end do
  end do 
  end do


  do l=1,nga
  do m=1,2
  do kk=1,nga
  if (rngoh(i,kk).ne.0) then
    dxkdnic(j,m,l)=dxkdni(kk,m,l)   
  end if
  end do
    end do
  end do

!c fin del c�lculo de la derivada de la fracci�n no asociada respecto a la 
!c fracci�n molar del componente
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C dgasdx(M) = derivada LOGARITMO NATURAL DEL GAMA DE ASOCIACI�N DEL COMPONENTE I
    if(nga.ne.0) then
 
    DO l=1,NGA 
     SUMA1J = 0.D0

    IF(MASS(l).NE.0) THEN      
    DO K=1,MASS(l)
  If(XOH(K,l).gt.1d-13)then

    SUMA1J = SUMA1J + RNGOH(i,l)*1.D0/XOH(K,l)*dxkdnic(j,k,l)+0.5D0*&
   r(i)*(drhodni(j,k,l)-xnoh(l)*dxkdnic(j,k,l)-drhodni(j,k,l)*&
   XOH(K,l))


  end if
    END DO
    ELSE
    CONTINUE
    END IF
  xx=xx+ SUMA1J
   end do   
  end if
      end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCALFONSINACCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    DACT(I,J)=XX                                                      
 70 DACT(J,I)=XX    
     IF(NDIF.LT.3) GOTO 100                                           
    DO 80 I=1,NK                                                      
    GAM(I)=DEXP(ACT(I))                                               
 80 ACT(I)=GAM(I)*X(I)                                                
    DO 85 I=1,NK                                                      
    DO 85 J=1,NK                                                      
    DACT(I,J)=ACT(I)*(DACT(I,J)-1.D0)                                 
    IF(J.EQ.I)DACT(I,J)=DACT(I,J)+GAM(I)                              
 85 CONTINUE                                                          
    GOTO 100                                                          
 90 CONTINUE                                                          
    DO 91 I=1,2                                                       
    DO 91 K=1,2                                                       
    DTAU(I,K,K)=0.                                                    
    DO 91 L=1,2                                                       
    IF(L.EQ.K) GOTO 91                                                
    H1=TETAR(L)-QT(L,I)*ETA(L)/S(L,I)                                 
    H2=QT(K,I)-S(L,I)*TETAR(K)/ETA(L)                                 
    DTAU(I,K,L)=-H1*H2/ETA(L)                                         
 91 CONTINUE                                                          
    DO 92 I=1,NK                                                      
    PACT(I,1)=-DTAU(I,1,2)*TAU(1,2)/T*300.D0                          
 92 PACT(I,2)=-DTAU(I,2,1)*TAU(2,1)/T*300.D0                          
100 RETURN                                                            
    END