SUBROUTINE GAMINF(N1,N2,GAM)                                      
    IMPLICIT REAL*8(A-H,O-Z)                                          
    COMMON/CUFAC/NK,NG,P(10,10),T                                     
    COMMON/CPAR/TAU(10,10),S(10,10),F(10)                             
    COMMON/CQT/QT(10,10),Q(10),R(10)                                  
!c------
    common/asoc/nktt,igamt(20,12),nytt(20,12)  
    common/nga/nga,mass(12)
    common/grupas1/rkass(6,12,6,12),enass(6,12,6,12), deloh(6,12,6,12) !Alfonsina
!c------
    dimension xohi0(10),xnohi0(10),tgt(10),xgamk(20),x(2)
  common/ioh2/rngoh(12,12)
  common/ioh2sis/rngoht(3,2)

!c------
    dk=1.381e-23
    deloh=0.0
    xnoh=0.0
    xnoh1=0.0
  xoh=0.0
    xgam=0.0
    xgamt=0.0
    do 7777 i=1,10
    xnohi0(i)=0.0
    tgt(i)=0.0
7777 continue
    do 8888 k=1,20
    xgamk(k)=0.0
8888 continue
!c------
!c      x(n1)=1.0
!c      x(n2)=0.0
!c------


!c      do 33 jc=1,2
!c      if(jc.eq.1) then
!c      i=n1
!c      else if(jc.eq.2) then
!c      i=n2
!c     end if
!c      goh(i)=0
!c      tgt(i)=0.0
!c      xnohi0(i)=0.0
!c      xgam=0.0
!c      do 33 j=1,nktt
!c      if(((igamt(j,i).eq.14).or.(igamt(j,i).eq.17)).and.(nytt(j,i). 
!c     *ne.0)) goh(i)=nytt(j,i)
!c   33 continue
!c      do 32 jc=1,2
!c      if(jc.eq.1) then
!c      i=n1
!c      else if(jc.eq.2) then
!c      i=n2
!c      end if
!c      if(nga.eq.1) then
!c      do 5 k=1,nktt
!c      tgt(i)=tgt(i)+nytt(k,i)
!c   5 continue
    !alterar
!c	xnoh1=xnoh1+goh(i)*x(i)
  !write (1,*) "xnoh1parte2=", xnoh1 
!c      end if
!c      xnohi0(i)=goh(i)/R(i)
  !write (1,*) "xnohi0parte2(",i,")=", xnohi0(i) 
!c   32 continue
!c      do 7 jc=1,2
!c      if(jc.eq.1) then
!c      i=n1
!c      else if(jc.eq.2) then
!c      i=n2
!c      end if
!c	!alterar
!c      xgam=xgam+R(i)*x(i)
  !write (1,*) "xgamparte2=", xgam
!c    7 continue
    !alterar
!c     xnoh=xnoh1/xgam
  !write (1,*) "xnohparte2=", xnoh 
!c------
    Q1=Q(N2)/Q(N1)                                                    
    R1=R(N2)/R(N1)                                                    
    QR=R1/Q1                                                          
    GAM=F(N2)+Q(N2)*(1.-DLOG(Q1))-R1+DLOG(R1)-5.D0*Q(N2)*(1.-QR+DLOG(R1)-DLOG(Q1))                                                      
    DO 10 I=1,NG                                                      
 10 GAM=GAM-S(I,N2)/S(I,N1)*QT(I,N1)-QT(I,N2)*DLOG(S(I,N1))           
!c------
!c      deloh=p4*(dexp(p5/t)-1.)
!c      if((xnohi0(n2).eq.0.0).or.(deloh.eq.0.0)) then
!c      xohi0(n2)=1.0
!c      else
!c      xohi0(n2)=(-1.+dsqrt(1.+4.*xnohi0(n2)*deloh))/(2.*xnohi0(n2)*
!c     *deloh)
!c      end if
!c      if((xnoh.eq.0.0).or.(deloh.eq.0.0)) then
!c      xoh=1.0
!c      else
!c      xoh=(-1.+dsqrt(1.+4.*xnoh*deloh))/(2.*xnoh*deloh)
!c      end if
!	!alterar
!c      gam=gam+goh(n2)*(2.*dlog(xoh/xohi0(n2))+xohi0(n2)-xoh)-(2.*xoh-
!c     *xoh**2.)*deloh*xnoh*(goh(n2)-R(n2)*xnoh)/(1.+2.*xnoh*xoh*deloh)
!c------
    RETURN                                                            
    END