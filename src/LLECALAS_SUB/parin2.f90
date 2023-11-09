SUBROUTINE PARIN2           
    use InputData                                      
    IMPLICIT REAL*8(A-H,O-Z)    
  !PARAMETER(NMS=2) !,NCOM=3,NSCM=10)                                      
!c------
    common/asoc/nktt,igamt(20,12),nytt(20,12)  !,p4,p5
    common/grupas1/rkass(6,12,6,12),enass(6,12,6,12), deloh(6,12,6,12)!Alfon
    common/nga/nga,mass(12) 
      common/ioh2/rngoh(12,12)

  common/ioh2sis/rngoht(3,2)
!c------
    COMMON/CUFAC/NK,NG,P(10,10),T                                     
    COMMON/CQT/QT(10,10),Q(10),R(10)                                  
    !COMMON/CMODEL/MODEL                                               
    !COMMON/COUT/IOUT                                                  
    DIMENSION RT(10,10),A(100,100),NGM(10) !,MAINSG(57)                   
    DIMENSION MS(10,10,2),NY(10,20),JH(150),IH(20)  

  DIMENSION NPUNTA(NMG),NGRUPA(NMG),NPUNT(NMG),NGRUP(NMG) !xnoh1(12),&
  !RR(NMG),
!    ,ENASS(6,12,6,12),&
!    RKASS(6,12,6,12),deloh(6,12,6,12),R(NC),xnohi0(12,12),&
!    X(NC),xnoh(12),dif(12,12),xoh_old(6,12),xoh(6,12),&
!    xohi0(12,6,12),xohi0_old(12,6,12),xxohi0(12,6,12), dif1(10,12,12),&
!    act(NC),actas(NC),NGRUPASD(NMG),NPUNTASD(NMG)
  INTEGER CS(NMG),TS(NMG,NMS),NPUNTMG(NMG),NMAINGR(NMG),J,K,CANIL(10) !,A,B,rngoh(NC,NSCM),GA,&
           !
      real*8::par             
    logical VL
    external mainsgfunc
                       
    REAL*4 RR(150),QQ(150)                                              
    REAL*4    A1(32),A2(32),A3(32),A4(32),A5(32),A6(32),A7(32),A8(32),&
   A9(32),A10(32),A11(32),A12(32),A13(32),A14(32),A15(32),A16(32),A17(32)&
   ,A18(32),A19(32),A20(32),A21(32),A22(32),A23(32),A24(32),A25(32)&
   ,A26(32),A27(32),A28(32),A29(32),A30(32),A31(32),A32(32)        
 !   DATA MAINSG/4*1,4*2,2*3,3*4,5,6,7,8,9,2*10,11,12,2*13,2*14,4*15,3*16,3*17,2*18,19,20,2*21,22,3*23,24,25,26,3*27,28,29,30,31,32/     
   DATA RR/.9011,.6744,.4469,.2195,1.3454,1.1167,.8886,1.1173,.5313,.3652,&
   1.2663,1.0396,.8121,1.,3.2499,3.2491,.92,.8952,1.6724,1.4457,&
  .998,3.168,1.3013,1.528,1.9031,1.6764,1.145,.9183,.6908,.9183,1.4654,&
  1.238,1.006,2.2564,2.0606,1.8016,2.87,2.6401,3.39,1.1562,1.870132&
  ,1.6434,1.06,2.0086,1.7818,1.5544,1.4199,2.4088,4.0013,2.9993,2.83&
  ,2.667,3.3092,2.4317,3.0856,4.0358,2.8266,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0/                    
  DATA QQ/.848,.54,.228,0.,1.176,.867,.676,.988,.4,.12,.968,.66,.348&
  ,1.2,3.128,3.124,1.4,.68,1.488,1.18,.948,2.484,1.224,1.532,1.728,1.42&
  ,1.088,.78,.468,1.1,1.264,.952,.724,1.988,1.684,1.448,2.41,2.184&
  ,2.91,.844,1.724,1.416,.816,1.868,1.56,1.248,1.104,2.248,3.568,2.113&
  ,1.833,1.553,2.86,2.192,2.736,3.2,2.472,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,&
  0.0,0.0,0.0/                      
  DATA A1/0.,292.3,156.5,104.4,328.2,-136.7,-131.9,342.4,-159.8,66.56&
  ,146.1,14.78,1744.,-320.1,1571.,73.8,27.9,21.23,89.97,-59.06,29.08&
  ,175.8,94.34,193.6,108.5,81.49,-128.8,147.3,-11.91,14.91,67.84,36.42/                                                              
   DATA A2/74.54,0.,-94.78,-269.7,470.7,-135.7,-135.7,220.6,1.,306.1,&
  517.,1.,-48.52,485.6,76.44,-24.36,-52.71,-185.1,-293.7,1.,34.78,1.&
  ,375.4,5*1.,176.7,132.1,42.73,60.82/                              
   DATA A3/-114.8,340.7,0.,-146.8,-9.21,-223.,-252.,372.8,-473.2,-78.31&
  ,-75.3,-10.44,75.49,114.8,52.13,4.68,1.,288.5,-4.7,777.8,56.41,-218.9&
  ,113.6,7.18,247.3,-50.71,-255.3,1.,-80.48,-17.78,59.16,29.77/
   DATA A4/-115.7,4102.,167.,0.,1.27,-162.6,-273.6,203.7,-470.4,-73.87&
  ,223.2,-184.9,147.3,-170.,65.69,122.9,1.,33.61,134.7,-47.13,-53.29&
  ,-15.41,-97.05,-127.1,453.4,-30.28,-124.6,3*1.,26.59,55.97      /
   DATA A5/644.6,724.4,703.9,4000.,0.,-281.1,-268.8,-122.4,-63.15,216.&
  ,-431.3,444.7,118.4,180.6,137.1,455.1,669.2,418.4,713.5,1989.,2011.&
  ,529.,483.8,332.6,-289.3,-99.56,-319.2,837.9,4*1.              /
   DATA A6/329.6,1731.,511.5,136.6,937.3,2*0.,247.,-547.,401.7,643.4,&
  -94.64,728.7,-76.64,-218.1,351.5,-186.1,-465.7,-260.3,3*1.,264.7,9*1./                                                              
   DATA A7/310.7,1731.,577.3,906.8,991.3,2*0.,104.9,-547.2,-127.6,231.4&
  ,732.3,349.1,-152.8,-218.1,351.5,-401.6,-465.7,512.2,3*1.,264.7,9*1./                                                             
   DATA A8/1300.,896.,859.4,5695.,28.73,-61.29,5.89,0.,-595.9,634.8,623.7&
  ,211.6,652.3,385.9,212.8,770.,740.4,793.2,1205.,390.7,63.48,-239.8&
  ,13.32,439.9,-424.3,1.,203.,1153.,-311.,-262.6,1.11,1.       /
   DATA A9/2255.,1.,1649.,292.6,-195.5,-153.2,-153.2,344.5,0.,-568.,3*1.&
  ,-337.3,4*1.,1616.,2*1.,-860.3,1.,-230.4,523.,1.,-222.7,5*1.  /
   DATA A10/472.6,343.7,593.7,916.7,67.07,-47.41,353.8,-171.8,-825.7,&
  0.,128.,48.93,-101.3,58.84,52.38,483.9,550.6,342.2,550.,190.5,-349.2&
  ,857.7,377.,211.6,82.77,2*1.,417.4,4*1.          /              
   DATA A11/158.1,-214.7,362.3,1218.,1409.,-344.1,-338.6,-349.9,1.,-37.36,&
   0.,-311.6,1051.,1090.,1.,-47.51,16*1./                       
   DATA A12/383.,1.,31.14,715.6,-140.3,299.3,-241.8,66.95,1.,120.3,1724.&
  ,0.,-115.7,-46.13,2*1.,808.8,203.1,70.14,5*1.,-75.23,1.,-201.9,123.2,1.,&
  -281.9,2*1./                                             
   DATA A13/139.4,1647.,461.8,339.1,-104.,244.4,-57.98,-465.7,1.,1247.&
  ,.75,1919.,0.,1417.,1402.,337.1,437.7,370.4,438.1,1349.,1.,681.4,&
  152.4,1.,-1707.,2*1.,639.7,4*1.          /                        
   DATA A14/972.4,-577.5,6.,5688.,195.6,19.57,487.1,-6.32,-898.3,258.70&
  ,-245.8,57.7,-117.6,0.,461.3,1.,-132.9,176.5,129.5,-246.3,2.41,3*1.,&
  29.86,7*1./                      
   DATA A15/662.1,289.3,32.14,213.1,262.5,1970.,1970.,64.42,1.,5.202,2*1.&
  ,-96.62,-235.7,0.,225.4,-197.7,-20.93,113.9,3*1.,-94.49,9*1. /
   DATA A16/42.14,99.61,-18.81,-114.1,62.05,-166.4,-166.4,315.9,1.,1000.&
  ,751.8,1.,19.77,1.,301.1,0.,-21.35,-157.1,11.8,13*1.       /   
   DATA A17/-243.9,337.1,2*1.,272.2,128.6,507.8,370.7,1.,-301.,1.,-347.9,&
   1670.,108.9,137.8,110.5,0.,1.,17.97,13*1.         /           
   DATA A18/7.5,4583.,-231.9,-12.14,-61.57,2*1544.,356.8,1.,12.01,1.,&
  -249.3,48.15,-209.7,-154.3,249.2,1.,0.,51.9,1.,-15.62,-216.3,4*1.,&
  -114.7,5*1. /                                                     
   DATA A19/-5.55,5831.,3000.,-141.3,-41.75,224.6,-207.,502.9,4894.,&
  -10.88,1.,61.59,43.83,54.57,47.67,62.42,56.33,-30.1,0.,-255.4,-54.86,&
  8455.,-34.68,514.6,8*1. /                                       
   DATA A20/924.8,1.,-878.1,-107.3,-597.1,2*1.,-97.27,1.,902.6,2*1.,&
  874.3,629.,4*1.,475.8,0.,-465.2,1.,794.4,1.,-241.7,1.,-906.5,5*1. /
   DATA A21/696.8,405.9,29.13,1208.,-189.3,2*1.,198.3,1.,430.6,3*1.,&
   -149.2,3*1.,70.04,492.,346.2,0.,5*1.,-169.7,5*1.        /          
   DATA A22/902.2,1.,1.64,689.6,-348.2,1.,1.,-109.8,-851.6,1010.,2*1.,&
   942.2,4*1.,-75.5,1302.,2*1.,0.,1.,175.8,164.4,1.,-944.9,5*1./    
   DATA A23/556.7,425.7,-1.77,3629.,-30.7,150.8,150.8,1538.6,1.,400.,&
   2*1.,446.3,1.,95.18,3*1.,490.9,-154.5,2*1.,0.,1.,481.3,7*1.      /
   DATA A24/575.7,1.,-11.19,-175.6,-159.,2*1.,32.92,-16.13,-328.6,8*1.,&
   534.7,2*1.,179.9,1.,0.,-246.,7*1.           /                   
   DATA A25/527.5,1.,358.9,337.7,536.6,2*1.,-269.2,-538.6,211.6,1.,&
   -278.2,572.7,343.1,5*1.,124.8,1.,125.3,139.8,963.,0.,7*1.    /      
   DATA A26/269.2,1.,363.5,1023.,53.37,20*1.,0.,6*1.          /      
   DATA A27/-300.,1.,-578.2,-390.7,183.3,2*1.,-873.6,-637.3,2*1.,-208.4,&
   5*1.,18.98,1.,-387.7,134.3,924.5,4*1.,0.,5*1.     /            
   DATA A28/-63.6,3*1.,-44.44,2*1.,1429.,1.,148.,1.,-13.91,-2.16,14*1.,0.,&
   4*1./                                                        
   DATA A29/928.3,500.7,364.2,4*1.,-364.2,20*1.,0.,3*1.    /         
   DATA A30/331.,115.4,-58.1,4*1.,-117.4,3*1.,173.8,17*1.,0.,2*1.   /
   DATA A31/561.4,784.4,21.97,238.,3*1.,18.41,22*1.,0.,1.       /    
   DATA A32/956.5,265.4,84.16,132.2,27*1.,0./                        
    DO 5 I=1,32                                                       
    A(I,1)=A1(I)                                                      
    A(I,2)=A2(I)                                                      
    A(I,3)=A3(I)                                                      
    A(I,4)=A4(I)                                                      
    A(I,5)=A5(I)                                                      
    A(I,6)=A6(I)                                                      
    A(I,7)=A7(I)                                                      
    A(I,8)=A8(I)                                                      
    A(I,9)=A9(I)                                                      
    A(I,10)=A10(I)                                                    
    A(I,11)=A11(I)                                                    
    A(I,12)=A12(I)                                                    
    A(I,13)=A13(I)                                                    
    A(I,14)=A14(I)                                                    
    A(I,15)=A15(I)                                                    
    A(I,16)=A16(I)                                                    
    A(I,17)=A17(I)                                                    
    A(I,18)=A18(I)                                                    
    A(I,19)=A19(I)                                                    
    A(I,20)=A20(I)                                                    
    A(I,21)=A21(I)                                                    
    A(I,22)=A22(I)                                                    
    A(I,23)=A23(I)                                                    
    A(I,24)=A24(I)                                                    
    A(I,25)=A25(I)                                                    
    A(I,26)=A26(I)                                                    
    A(I,27)=A27(I)                                                    
    A(I,28)=A28(I)                                                    
    A(I,29)=A29(I)                                                    
    A(I,30)=A30(I)                                                    
    A(I,31)=A31(I)                                                    
  5 A(I,32)=A32(I)                                                    
    IF(IOUT.EQ.0) IOUT=6                                              
!C     READ(2,501) NK                                                    
    READ(2,*) NK                                                      
    NC = NK                                                           
    DO 15 I=1,10                                                      
      DO 15 J=1,NK                                                      
          QT(I,J)=0.D0                                                      
 15 RT(I,J)=0.D0                                                      
    IF(MODEL.NE.1) GOTO 19                                            
    NG=NK                                                             
    DO 16 I=1,NK                                                      
!C  16 READ(2,502) RT(I,I),QT(I,I),(P(I,J),J=1,NK)                       
 16 READ(2,*) RT(I,I),QT(I,I),(P(I,J),J=1,NK)                         
 19 CONTINUE                                                          
    IF(MODEL.EQ.1) GOTO 21                                            
!C     READ(2,501) IOWNRQ,IOWNP                                          

!Lectura de par�metros R, Q, int
  read(2,*) IOWNRQ,IOWNP                                            

  if(IOWNRQ/=0)then
      do I=1,IOWNRQ                                                   
          read(2,*) K,RR(K),QQ(K)  !K= n�m de grupo                                     
      enddo
  endif            

10  continue

  if(IOWNP/=0)then
      do I=1,IOWNP                                                   
          read(2,*) J,K,par !A(J,K) !j y k= n�m de grupos
          j=mainsgfunc(j,ipareq)
          k=mainsgfunc(k,ipareq)
          a(j,k)=par
      enddo
  endif
                                                
14  CONTINUE                                                          
!c------
!ccccccccccccccccccccccccccAlfonsinaccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !rngoh=0.
  !LECTURA DE LA COMPOSICION GRUPAL ASOCIATIVA
! AGREGAR ALGO ACÁ QUE LEA 0,0 si no hay 10 pares de numero de grupo y cantidad
         do ja=1, nga
      read(2,*) (rngoh(i,ja),i=1,nc) 
  end do     
  


        !LECTURA DEL NUMERO DE SITIOS Y PARAMETROS ASOCIATIVOS

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA (basado en el Aparaest) CCCCCCCCCCCCCCCCCCCCC



  IF(NGA.GT.0) READ(2,*)(MASS(I),I=1,NGA)
  MS(:,:,:)=0                                                       
  JH(:)=0                                                           
  DO 50 I=1,NK
     ! do j=1,10
          READ(2,*) (MS(I,J,1),MS(I,J,2),j=1,size(ms(1,:,1)))    
     ! enddo
    !  READ(2,*) NGI,(MS(I,J,1),MS(I,J,2),J=1,NGI)    
      do 50 j=1,size(ms(1,:,1))
          igamt(j,i)=ms(i,j,1)
          nytt(j,i)=ms(i,j,2)
50  continue	

!****Jose S****

  NUM = 0
  NMGR = 0
  NGA = 0
  DO J=1,NC
    DO I=1,10 
        IF(MS(J,I,1).EQ.0)CYCLE
        IREC2=MS(J,I,1)+(ipareq-1)*150      
            READ(14,500,REC=IREC2)MGR,RRT,ICS,ITS1,ITS2 !500  FORMAT(2x,I2,37X,D15.8,120X,3I2)
            if(model==2) then
            IF (ICS.NE.0) THEN !GRUPOS ASOCIATIVOS
              NG = MS(J,I,1)
              CALL BUSCARAS (NG,NGRUPA,NMG,VL)
              IF(VL)THEN
                  IF(MS(J,I,1).EQ.10) THEN
                      RNGOH(J,NPUNTA(NG))=CANIL(J)
                  ELSE
                      RNGOH(J,NPUNTA(NG))=MS(J,I,2)
                  ENDIF                  
              ELSE
                  NGA=NGA+1
                  NPUNTA(NG) = NGA
                  NGRUPA(NGA) = NG
                  CS(NGA) = ICS
                  MASS(nga) = ICS
                  TS(NGA,1) = ITS1
                  TS(NGA,2) = ITS2
                  IF(MS(J,I,1).EQ.10) THEN
                      RNGOH(J,NGA)=CANIL(J)
                  ELSE
                      RNGOH(J,NGA)=MS(J,I,2)
                  ENDIF
              ENDIF
            ENDIF
          endif        	    
            IF (NPUNT(MS(J,I,1)).EQ.0) THEN !Ordena todos los grupos
              NUM = NUM + 1
              NPUNT(MS(J,I,1)) = NUM
              NGRUP(NUM) = MS(J,I,1)
             ! RR(NUM)=RRT
              CALL BUSCARAS (MGR,NMAINGR,NMG,VL)
              IF(.NOT.VL) THEN !Crea vectores NPUNTMG y NMAINGR
                  NMGR = NMGR+1
                  NPUNTMG(MGR) = NMGR
                  NMAINGR(NMGR) = MGR
              ENDIF
          END IF
      ENDDO
  ENDDO


!C....................................................................................
!C.....Genera las matrices ENASS Y RKASS con los par�metros de energ�a y volumen 
!c.....de asociaci�n, respectivamente, seg�n los grupos asociativos de los componentes
!c.....del sistema que se est� corriendo. (si el modelo elegido es A-UNIFAC)
!C....................................................................................     
  enass(:,:,:,:)=0.0
  rkass(:,:,:,:)=0.0
  if(model==2)then
    DO J=1,nga
      DO K=1,nga
          DO BB=1,CS(J)
              DO AA=1,CS(K)
                  IREC1 = MAINSGfunc(NGRUPA(K),IPAREQ)+(ipareq-1)*70
                  IF(TS(J,BB).EQ.1.OR.TS(K,AA).EQ.1)THEN !Si alguno de ambos sitios es del tipo 1
                     CALL LEEPAR (J,IREC1,IPAREQ,NGRUPA,ENASST,RKASST)
                     ENASS(AA,K,BB,J)=ENASST
                     RKASS(AA,K,BB,J)=RKASST                        
                  ELSEIF (TS(J,BB).NE.TS(K,AA)) THEN
                     CALL LEEPAR (J,IREC1,IPAREQ,NGRUPA,ENASST,RKASST)
                     ENASS(AA,K,BB,J)=ENASST
                     RKASS(AA,K,BB,J)=RKASST
                  ELSE !ELSEIF (TS(J,B).EQ.TS(K,A)) THEN
                     ENASS(AA,K,BB,J)=0.0
                     RKASS(AA,K,BB,J)=0.0     
                  ENDIF   
              ENDDO
          ENDDO
      ENDDO
    ENDDO
  endif

!****Jose S****


!!C
!!c     lectura parametros ENERG�TICOS DE asociacion
!!c
!	DO J=1,NGA
!			IF(MASS(J).EQ.0) GO TO 201
!		DO L=1,NGA
!				IF(MASS(L).EQ.0) GO TO 103
!			DO I=1,MASS(J)
!				DO K=1,MASS(L)
!						READ(2,*)ENASS(I,J,K,L) 
!				END DO
!			END DO
!  103			CONTINUE
!			END DO
!  201		   CONTINUE
!	END DO
!!C
!!c     lectura parametros VOLUM�TRICOS DE asociacion
!!C
!      DO J=1,NGA
!			IF(MASS(J).EQ.0) GO TO 3330
!		DO L=1,NGA
!				IF(MASS(L).EQ.0) GO TO 5550
!			DO I=1,MASS(J)
!				DO K=1,MASS(L)
!						READ(2,*) RKASS(I,J,K,L)
! 179						FORMAT (F16.10)
!				END DO
!			END DO
! 5550			CONTINUE 
!		END DO
! 3330			CONTINUE
!	END DO


!ccccccccccccccccccccccccccAlfonsinaccccccccccccccccccccccccccccccccccccccccccccccccccccc


    IC=1                                                              
    DO 71 I=1,NK                                                      
      DO 70 J=1,10                                                      
          IF(MS(I,J,1).EQ.0) GOTO 71                                        
          IH(IC)=MS(I,J,1)                                                  
          IF(IC.EQ.1) GOTO 69                                               
          IF(IH(IC).EQ.IH(IC-1)) GOTO 70                                    
          IF(IH(IC).GT.IH(IC-1)) GOTO 69                                    
          IF(IC.GT.2) GOTO 55                                               
          IHH=IH(1)                                                         
          IH(1)=IH(2)                                                       
          IH(2)=IHH                                                         
          GOTO 69                                                           
 55       I1=IC-1                                                           
          DO 65 I2=1,I1                                                     
              IF(IH(IC).GT.IH(I2)) GOTO 65                                      
              IF(IH(IC).EQ.IH(I2)) GOTO 70                                      
              I4=IC-I2                                                          
              DO 61 I3=1,I4                                                     
 61               IH(IC+1-I3)=IH(IC-I3)                                             
              IH(I2)=MS(I,J,1)                                                  
 65       CONTINUE                                                          
 69       IC=IC+1                                                           
          IF(IC.GT.20) WRITE(6,607)                                         
          IF(IOUT.NE.6.AND.IC.GT.20) WRITE(IOUT,607)                        
 70   CONTINUE                                                          
 71 CONTINUE                                                          
    IC=IC-1                                                           
!c------
    nktt=ic
!c------
    DO 73 I=1,IC                                                      
 73   JH(IH(I))=I                                                       
    DO 72 I=1,10                                                      
      DO 72 J=1,20                                                      
 72       NY(I,J)=0                                                         
    DO 75 I=1,NK                                                      
      DO 74 J=1,10                                                      
          IF(MS(I,J,1).EQ.0) GOTO 75                                        
          N1=MS(I,J,1)                                                      
          N2=MS(I,J,2)                                                      
          IF(N1.EQ.0) GOTO 75                                               
          N3=JH(N1)                                                         
 74   NY(I,N3)=N2                                                       
 75 CONTINUE                                                          
    I=0                                                               
    NGMGL=0                                                           
    DO 80 K=1,IC                                                      
      NSG=IH(K)                                                         
      NGMNY=MAINSGfunc(NSG,ipareq)                                                 
      IF(NGMNY.NE.NGMGL) I=I+1                                          
      NGM(I)=NGMNY                                                      
      NGMGL=NGMNY                                                       
    DO 80 J=1,NK                                                      
    RT(I,J)=RT(I,J)+NY(J,K)*RR(NSG)                                   
 80 QT(I,J)=QT(I,J)+NY(J,K)*QQ(NSG)                                   
    NG=I                                                              
    WRITE(6,608) (IH(K),K=1,IC)                                       
    WRITE(6,609) (MAINSGfunc(IH(K),ipareq),K=1,IC)                               
    WRITE(6,610)                                                      
    DO 90 I=1,NK                                                      
 90 WRITE(6,611) I,(NY(I,K),K=1,IC)                                   
    WRITE(6,699)                                                      
    IF(IOUT.EQ.6) GOTO 85                                             
    WRITE(IOUT,608) (IH(K),K=1,IC)                                    
    WRITE(IOUT,609) (MAINSGfunc(IH(K),ipareq),K=1,IC)                            
    WRITE(IOUT,610)                                                   
    DO 91 I=1,NK                                                      
 91 WRITE(IOUT,611) I,(NY(I,K),K=1,IC)                                
    WRITE(IOUT,699)                                                   
 85 CONTINUE                                                          
    DO 20 I=1,NG                                                      
    DO 20 J=1,NG                                                      
    NI=NGM(I)                                                         
    NJ=NGM(J)                                                         
 20 P(I,J)=A(NI,NJ)                                                   
    WRITE(6,612)                                                      
    DO 95 K=1,IC                                                      
    NN=IH(K)                                                          
 95 WRITE(6,613) NN,RR(NN),QQ(NN)                                     
    WRITE(6,699)                                                      
    IF(IOUT.EQ.6) GOTO 99                                             
    WRITE(IOUT,612)                                                   
    DO 96 K=1,IC                                                      
    NN=IH(K)                                                          
 96 WRITE(IOUT,613) NN,RR(NN),QQ(NN)                                  
    WRITE(IOUT,699)                                                   
 99 CONTINUE                                                          
 21 CONTINUE                                                          
    WRITE(6,604)   
                                                  
    DO 25 I=1,NG                                                      
 25 WRITE(6,603) (P(I,J),J=1,NG)                                      
    WRITE(6,699)                                                      
    IF(MODEL.EQ.0) WRITE(6,605)                                       
    IF(MODEL.EQ.1) WRITE(6,627)                                       
    IF(IOUT.EQ.6) GOTO 26                                             
    WRITE(IOUT,604)                                                   
    DO 27 I=1,NG                                                      
 27 WRITE(IOUT,603) (P(I,J),J=1,NG)                                   


!ccccccccccccccccc Escritura de los par�metros de asociaci�n ALFONSINAccccccccccccccccc
     IF (NGA.GT.0) THEN

!-------------------------------------------------------------------------------
!BLOQUE ORIGINAL
  !WRITE(1,218) (I, I=1,NC)
  !218	FORMAT(/,X,'"ASSOC GROUP COMPOSITION" ',/,23X,'COMPONENTES',/,' GRUPO 	  #  SIT ASOC  ',I5,/)
      !DO ja=1,NGA
        !WRITE(1,219) ja,MASS(ja)  , (rngoh(i,ja),i=1,nc)   
      !219	FORMAT(3X,I3,9X,I3,6X,f5.1)
        !END DO
!-------------------------------------------------------------------------------
!BLOQUE MODIFICADO PORQUE TIRA ERRORES. NO ESTA BIEN ESTO PERO SINO NO SE PUEDE SEGUIR
 WRITE(1,218) (I, I=1,NC)
218	FORMAT(/,X,'"ASSOC GROUP COMPOSITION" ',/,23X,'COMPONENTES',/,' GRUPO 	  #  SIT ASOC  ',2I5,/)  

  DO ja=1,NGA
  WRITE(1,*) ja, MASS(ja), (rngoh(i,ja),i=1,nc)    
219	FORMAT(3X,I3,9X,I3,6X,f5.1)
  END DO
!-------------------------------------------------------------------------------
  WRITE(1,220)
220	FORMAT(/,X,'PARAMETROS DE ENERGIA DE ASOCIACION (Kelvin)  ',/)  

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  DO J=1,NGA
          IF(MASS(J).EQ.0) GO TO 202
      DO L=1,NGA
              IF(MASS(L).EQ.0) GO TO 104
          DO I=1,MASS(J)
              DO K=1,MASS(L)
                      WRITE(1,221) I,J,K,L,ENASS(I,J,K,L)
221	FORMAT(X,' ENASS( ',I3,I3,I3,I3,' ) = ',F10.4)
              END DO
          END DO
104			CONTINUE
          END DO
202			CONTINUE
  END DO


  WRITE(1,222)
222	FORMAT (/,X,'PARAMETROS DE VOLUMEN DE ASOCIACI�N (cm3/mol) ',/)

  DO J=1,NGA
          IF(MASS(J).EQ.0) GO TO 301
      DO L=1,NGA
              IF(MASS(L).EQ.0) GO TO 5011
          DO I=1,MASS(J)
              DO K=1,MASS(L)
                  WRITE(1,223) I,J,K,L,RKASS(I,J,K,L)
223	FORMAT(X,' RKASS( ',I3,I3,I3,I3,' ) = ',F10.4)
              END DO
          END DO
5011			CONTINUE
          END DO
301		   CONTINUE
  END DO
  END IF
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA (basado en el Aparaest) CCCCCCCCCCCCCCCCCCCCC
!c------
    WRITE(IOUT,699)                                                   
    IF(MODEL.EQ.0) WRITE(IOUT,605)                                    
    IF(MODEL.EQ.1) WRITE(IOUT,627)                                    
 26 CONTINUE                                                          
    DO 30 I=1,NK                                                      
    Q(I)=0.D0                                                         
    R(I)=0.D0                                                         
    DO 30 K=1,NG                                                      
    Q(I)=Q(I)+QT(K,I)                                                 
 30 R(I)=R(I)+RT(K,I)                                                 
    DO 40 I=1,NK                                                      
 40 WRITE(6,606) I,R(I),Q(I)                                          
    IF(IOUT.EQ.6) GOTO 42                                             
    DO 41 I=1,NK                                                      
 41 WRITE(IOUT,606) I,R(I),Q(I)                                       
 42 CONTINUE      
500 FORMAT(2x,I2,37X,D15.8,120X,3I2)                                                     
501 FORMAT(20I3)                                                      
502 FORMAT(8F10.2)                                                    
503 FORMAT(I3,2F10.2)                                                 
504 FORMAT(2I3,F10.2)                                                 
603 FORMAT(1X,10F12.3)                                                
604 FORMAT('  INTERACTION PARAMETERS',/)                              
605 FORMAT(' UNIFAC MOLECULAR R AND Q',/)                             
606 FORMAT(I5,2F15.4)                                                 
607 FORMAT(' ** WARNING: NUMBER OF SUB GROUPS MUST NOT EXCEED 20 **') 
608 FORMAT(//,' SUB GROUPS :',20I3)                                   
609 FORMAT(' MAIN GROUPS:',20I3)                                      
610 FORMAT(' COMPONENT')                                              
611 FORMAT(6X,I2,5X,20I3)                                             
612 FORMAT(' GROUP R- AND Q-VALUES',/)                                
613 FORMAT(1X,I3,2F10.4)                                              
627 FORMAT(' SPECIFIED UNIQUAC R AND Q',/)                            
699 FORMAT(//)                                                        
!c------
1603 format('  ASSOCIATION PARAMETERS',//,10X,'K(OH)   :',F7.3,/,10X,'E(OH)/k :',F7.1,' K-1')
!c------
    RETURN                                                            
    END                   