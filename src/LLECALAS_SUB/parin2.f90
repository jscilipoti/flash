SUBROUTINE PARIN2           
    use InputData
    use iso_fortran_env, only: real64, real32, int32
    use fileUnits, only: intrcn32_unit, qPar150_unit, rPar150_unit
    
    
    IMPLICIT REAL*8(A-H, O-Z)    
  
    common/asoc/nktt, igamt(20, 12), nytt(20, 12)
    common/grupas1/rkass(6, 12, 6, 12), enass(6, 12, 6, 12), deloh(6, 12, 6, 12)!
    common/nga/nga, mass(12) 
    common/ioh2/rngoh(12, 12)

    common/ioh2sis/rngoht(3, 2)
    COMMON/CUFAC/n_comp, NG, P(10, 10), T                                     
    COMMON/CQT/QT(10, 10), Q(10), R(10)                                  

    !DIMENSION RT(10, 10), A(100, 100), NGM(10) !, MAINSG(57)
    !DIMENSION RT(10, 10), NGM(10)
    DIMENSION NGM(10)         
    DIMENSION MS(10, 10, 2), NY(10, 20), JH(150), IH(20)  

  DIMENSION NPUNTA(NMG), NGRUPA(NMG), NPUNT(NMG), NGRUP(NMG) 

  !INTEGER CS(NMG), TS(NMG, NMS), NPUNTMG(NMG), NMAINGR(NMG), j, k, CANIL(10)
  INTEGER CS(NMG), TS(NMG, NMS), NPUNTMG(NMG), NMAINGR(NMG), CANIL(10)

    real*8::intrcn_par             
    logical VL
    external mainsgfunc, get_free_unit
    integer(int32) :: get_free_unit

    !USED:
    integer(kind = int32) :: &
        & i, j, k          
    ! Here the interaction parameters are dim 100 because
    ! there is a maximum of 10 different functional groups in total
    real(kind = real64), dimension(100, 100) :: &
        & intrcnPar_matrix = 0.D0, &
        & intrcn32_data = 0.D0

    real(kind = real64), dimension(10,10) :: &
        & RT = 0.D0

    ! Arrays with the r and q parameters. Dim = 150 because there are 150 values.
    real(kind = real64), dimension(150) :: & !if real64, it prints slightly  different values. 
        & rPar_array = 0.D0, &
        & qPar_array = 0.D0, &
        & rPar150_data = 0.D0, &
        & qPar150_data = 0.D0

    !COMMON
    integer(kind = int32) :: &
    & n_comp ! Number of components in the system

    real(kind = real64) :: &
        & QT = 0.D0
   
    ! New method to fill up intrcn32_data (new intrcnPar_matrix)
    intrcn32_unit = get_free_unit()
    open(unit = intrcn32_unit, file = "src/database/intrcn32.mds", &
    &status = 'old', action = 'read', form = "formatted")
    read(intrcn32_unit, *) intrcn32_data(:32, :32)
    close(intrcn32_unit)

    ! New method to fill up rPar150_data (new rPar_array matrix)
    rPar150_unit = get_free_unit()
    open(unit = rPar150_unit, file = "src/database/rPar150.mds", &
    &status = 'old', action = 'read', form = "formatted")
    read(rPar150_unit, *) rPar150_data(:)
    close(rPar150_unit)

    ! New method to fill up qPar150_data (new qPar_array matrix)
    qPar150_unit = get_free_unit()
    open(unit = qPar150_unit, file = "src/database/qPar150.mds", &
    &status = 'old', action = 'read', form = "formatted")
    read(qPar150_unit, *) qPar150_data(:)
    close(qPar150_unit)

    intrcnPar_matrix = intrcn32_data
    rPar_array = rPar150_data
    qPar_array = qPar150_data

    if (IOUT == 0) IOUT = 6

    read(2, *) n_comp                                                      
                                                          
    ! QT and RT are filled up better above
    ! do i = 1, 10                                                      
    !     do j = 1, n_comp                                                      
    !         QT(i, j) = 0.D0                                                      
    !         RT(i, j) = 0.D0
    !   end do
    ! end do


    !if (model /= 1) !GOTO 19
    if (model == 1) then !GOTO 19 

        NG = n_comp !REVISAR PORQUE HACE ESTO                                                             
    
        do i = 1, n_comp                                                                          
            read(2, *) RT(i, i), QT(i, i), (P(i, j), j = 1, n_comp)
        end do

    end if

    !19 continue

    if (model == 1) GOTO 21                                                                                   

    ! Read the group number and parameters R and Q. 
    read(2, *) total_groups, total_int_par                                            

    if (total_groups /=  0) then
        rPar_array = 0.D0
        qPar_array = 0.D0
        do i = 1, total_groups
            k = 0                                                   
            read(2, *) k, rPar_array(k), qPar_array(k)     
        end do
    end if            

    if (total_int_par /=  0) then
        intrcnPar_matrix = 0.D0
        do i = 1, total_int_par
            j = 0
            k = 0                                                   
            read(2, *) j, k, intrcn_par ! j and k are group numbers  
            j = mainsgfunc(j, ipareq)
            k = mainsgfunc(k, ipareq)
            intrcnPar_matrix(j, k) =  intrcn_par
        end do
    end if
                                                
                                                    
  !LECTURA DE LA COMPOSICION GRUPAL ASOCIATIVA
! AGREGAR ALGO ACÁ QUE LEA 0, 0 si no hay 10 pares de numero de grupo y cantidad
         do ja = 1, nga
      read(2, *) (rngoh(i, ja), i = 1, nc) 
  end do     
  


        !LECTURA DEL NUMERO DE SITIOS Y PARAMETROS ASOCIATIVOS

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA (basado en el Aparaest) CCCCCCCCCCCCCCCCCCCCC



  if (NGA > 0) read(2, *)(MASS(i), i = 1, NGA)
  MS(:, :, :) = 0                                                       
  JH(:) = 0                                                           
  do 50 i = 1, n_comp
     ! do j = 1, 10
          read(2, *) (MS(i, j, 1), MS(i, j, 2), j = 1, size(ms(1, :, 1)))    
     ! end do
    !  read(2, *) NGI, (MS(i, j, 1), MS(i, j, 2), j = 1, NGI)    
      do 50 j = 1, size(ms(1, :, 1))
          igamt(j, i) = ms(i, j, 1)
          nytt(j, i) = ms(i, j, 2)
50  continue	

!****Jose S****

  NUM = 0
  NMGR = 0
  NGA = 0
  do j = 1, n_comp
    do i = 1, 10 
        if (MS(j, i, 1) == 0) cycle
        IREC2 = MS(j, i, 1)+(ipareq-1)*150      
            read(14, 500, REC = IREC2)MGR, RRT, ICS, ITS1, ITS2 !500  FORMAT(2x, I2, 37X, D15.8, 120X, 3I2)
            if (model == 2) then
            if (ICS /= 0) then !GRUPOS ASOCIATIVOS
              NG = MS(j, i, 1)
              call buscaras (NG, NGRUPA, NMG, VL)
              if (VL)then
                  if (MS(j, i, 1) == 10) then
                      RNGOH(j, NPUNTA(NG)) = CANIL(j)
                  else
                      RNGOH(j, NPUNTA(NG)) = MS(j, i, 2)
                  end if                  
              else
                  NGA = NGA+1
                  NPUNTA(NG) = NGA
                  NGRUPA(NGA) = NG
                  CS(NGA) = ICS
                  MASS(nga) = ICS
                  TS(NGA, 1) = ITS1
                  TS(NGA, 2) = ITS2
                  if (MS(j, i, 1) == 10) then
                      RNGOH(j, NGA) = CANIL(j)
                  else
                      RNGOH(j, NGA) = MS(j, i, 2)
                  end if
              end if
            end if
          end if        	    
            if (NPUNT(MS(j, i, 1)) == 0) then !Ordena todos los grupos
              NUM = NUM + 1
              NPUNT(MS(j, i, 1)) = NUM
              NGRUP(NUM) = MS(j, i, 1)
             ! rPar_array(NUM) = RRT
              call buscaras (MGR, NMAINGR, NMG, VL)
              if (.NOT.VL) then !Crea vectores NPUNTMG y NMAINGR
                  NMGR = NMGR+1
                  NPUNTMG(MGR) = NMGR
                  NMAINGR(NMGR) = MGR
              end if
          end if
      end do
  end do


!C....................................................................................
!C.....Genera las matrices ENASS Y RKASS con los par�metros de energ�a y volumen 
!c.....de asociaci�n, respectivamente, seg�n los grupos asociativos de los componentes
!c.....del sistema que se est� corriendo. (si el modelo elegido es A-UNIFAC)
!C....................................................................................     
  enass(:, :, :, :) = 0.0
  rkass(:, :, :, :) = 0.0
  if (model == 2) then
    do j = 1, nga
      do k = 1, nga
          do BB = 1, CS(j)
              do AA = 1, CS(k)
                  IREC1 = MAINSGfunc(NGRUPA(k), IPAREQ)+(ipareq-1)*70
                  if (TS(j, BB) == 1.OR.TS(k, AA) == 1)then !Si alguno de ambos sitios es del tipo 1
                     call LEEPAR (j, IREC1, IPAREQ, NGRUPA, ENASST, RKASST)
                     ENASS(AA, k, BB, j) = ENASST
                     RKASS(AA, k, BB, j) = RKASST                        
                  ELSEif (TS(j, BB) /= TS(k, AA)) then
                     call LEEPAR (j, IREC1, IPAREQ, NGRUPA, ENASST, RKASST)
                     ENASS(AA, k, BB, j) = ENASST
                     RKASS(AA, k, BB, j) = RKASST
                  else !ELSEif (TS(j, B) == TS(k, A)) then
                     ENASS(AA, k, BB, j) = 0.0
                     RKASS(AA, k, BB, j) = 0.0     
                  end if   
              end do
          end do
      end do
    end do
  end if

!****Jose S****


!!C
!!c     lectura parametros ENERG�TICOS DE asociacion
!!c
!	do j = 1, NGA
!			if (MASS(j) == 0) GO TO 201
!		do L = 1, NGA
!				if (MASS(L) == 0) GO TO 103
!			do i = 1, MASS(j)
!				do k = 1, MASS(L)
!						read(2, *)ENASS(i, j, k, L) 
!				end do
!			end do
!  103			continue
!			end do
!  201		   continue
!	end do
!!C
!!c     lectura parametros VOLUM�TRICOS DE asociacion
!!C
!      do j = 1, NGA
!			if (MASS(j) == 0) GO TO 3330
!		do L = 1, NGA
!				if (MASS(L) == 0) GO TO 5550
!			do i = 1, MASS(j)
!				do k = 1, MASS(L)
!						read(2, *) RKASS(i, j, k, L)
! 179						FORMAT (F16.10)
!				end do
!			end do
! 5550			continue 
!		end do
! 3330			continue
!	end do


!ccccccccccccccccccccccccccAlfonsinaccccccccccccccccccccccccccccccccccccccccccccccccccccc


    IC = 1                                                              
    do 71 i = 1, n_comp                                                      
      do 70 j = 1, 10                                                      
          if (MS(i, j, 1) == 0) GOTO 71                                        
          IH(IC) = MS(i, j, 1)                                                  
          if (IC == 1) GOTO 69                                               
          if (IH(IC) == IH(IC-1)) GOTO 70                                    
          if (IH(IC) > IH(IC-1)) GOTO 69                                    
          if (IC > 2) GOTO 55                                               
          IHH = IH(1)                                                         
          IH(1) = IH(2)                                                       
          IH(2) = IHH                                                         
          GOTO 69                                                           
 55       I1 = IC-1                                                           
          do 65 I2 = 1, I1                                                     
              if (IH(IC) > IH(I2)) GOTO 65                                      
              if (IH(IC) == IH(I2)) GOTO 70                                      
              I4 = IC-I2                                                          
              do 61 I3 = 1, I4                                                     
 61               IH(IC+1-I3) = IH(IC-I3)                                             
              IH(I2) = MS(i, j, 1)                                                  
 65       continue                                                          
 69       IC = IC+1                                                           
          if (IC > 20) write(6, 607)                                         
          if (IOUT /= 6.AND.IC > 20) write(IOUT, 607)                        
 70   continue                                                          
 71 continue                                                          
    IC = IC-1                                                           
!c------
    nktt = ic
!c------
    do 73 i = 1, IC                                                      
 73   JH(IH(i)) = i                                                       
    do 72 i = 1, 10                                                      
      do 72 j = 1, 20                                                      
 72       NY(i, j) = 0                                                         
    do 75 i = 1, n_comp                                                      
      do 74 j = 1, 10                                                      
          if (MS(i, j, 1) == 0) GOTO 75                                        
          N1 = MS(i, j, 1)                                                      
          N2 = MS(i, j, 2)                                                      
          if (N1 == 0) GOTO 75                                               
          N3 = JH(N1)                                                         
 74   NY(i, N3) = N2                                                       
 75 continue                                                          
    i = 0                                                               
    NGMGL = 0                                                           
    do 80 k = 1, IC                                                      
      NSG = IH(k)                                                         
      NGMNY = MAINSGfunc(NSG, ipareq)                                                 
      if (NGMNY /= NGMGL) i = i+1                                          
      NGM(i) = NGMNY                                                      
      NGMGL = NGMNY                                                       
    do 80 j = 1, n_comp                                                      
    RT(i, j) = RT(i, j)+NY(j, k)*rPar_array(NSG)                                   
 80 QT(i, j) = QT(i, j)+NY(j, k)*qPar_array(NSG)                                   
    NG = i                                                              
    write(6, 608) (IH(k), k = 1, IC)                                       
    write(6, 609) (MAINSGfunc(IH(k), ipareq), k = 1, IC)                               
    write(6, 610)                                                      
    do 90 i = 1, n_comp                                                      
 90 write(6, 611) i, (NY(i, k), k = 1, IC)                                   
    write(6, 699)                                                      
    if (IOUT == 6) GOTO 85                                             
    write(IOUT, 608) (IH(k), k = 1, IC)                                    
    write(IOUT, 609) (MAINSGfunc(IH(k), ipareq), k = 1, IC)                            
    write(IOUT, 610)                                                   
    do 91 i = 1, n_comp                                                      
 91 write(IOUT, 611) i, (NY(i, k), k = 1, IC)                                
    write(IOUT, 699)                                                   
 85 continue                                                          
    do 20 i = 1, NG                                                      
    do 20 j = 1, NG                                                      
    NI = NGM(i)                                                         
    NJ = NGM(j)                                                         
 20 P(i, j) = intrcnPar_matrix(NI, NJ)                                                   
    write(6, 612)                                                      
    do 95 k = 1, IC                                                      
    NN = IH(k)                                                          
 95 write(6, 613) NN, rPar_array(NN), qPar_array(NN)                                     
    write(6, 699)                                                      
    if (IOUT == 6) GOTO 99                                             
    write(IOUT, 612)                                                   
    do 96 k = 1, IC                                                      
    NN = IH(k)                                                          
 96 write(IOUT, 613) NN, rPar_array(NN), qPar_array(NN)                                  
    write(IOUT, 699)                                                   
 99 continue                                                          
 21 continue                                                          
    write(6, 604)   
                                                  
    do 25 i = 1, NG                                                      
 25 write(6, 603) (P(i, j), j = 1, NG)                                      
    write(6, 699)                                                      
    if (model == 0) write(6, 605)                                       
    if (model == 1) write(6, 627)                                       
    if (IOUT == 6) GOTO 26                                             
    write(IOUT, 604)                                                   
    do 27 i = 1, NG                                                      
 27 write(IOUT, 603) (P(i, j), j = 1, NG)                                   


!ccccccccccccccccc Escritura de los par�metros de asociaci�n ALFONSINAccccccccccccccccc
     if (NGA > 0) then

!-------------------------------------------------------------------------------
!BLOQUE ORIGINAL
  !write(1, 218) (i, i = 1, n_comp)
  !218	FORMAT(/, X, '"ASSOC GROUP COMPOSITION" ', /, 23X, 'COMPONENTES', /, ' GRUPO 	  #  SIT ASOC  ', I5, /)
      !do ja = 1, NGA
        !write(1, 219) ja, MASS(ja), (rngoh(i, ja), i = 1, nc)   
      !219	FORMAT(3X, I3, 9X, I3, 6X, f5.1)
        !end do
!-------------------------------------------------------------------------------
!BLOQUE MODifICADO PORQUE TIRA ERRORES. NO ESTA BIEN ESTO PERO SINO NO SE PUEDE SEGUIR
 write(1, 218) (i, i = 1, n_comp)
218	FORMAT(/, X, '"ASSOC GROUP COMPOSITION" ', /, 23X, 'COMPONENTES', /, ' GRUPO 	  #  SIT ASOC  ', 2I5, /)  

  do ja = 1, NGA
  write(1, *) ja, MASS(ja), (rngoh(i, ja), i = 1, nc)    
219	FORMAT(3X, I3, 9X, I3, 6X, f5.1)
  end do
!-------------------------------------------------------------------------------
  write(1, 220)
220	FORMAT(/, X, 'PARAMETROS DE ENERGIA DE ASOCIACION (Kelvin)  ', /)  

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  do j = 1, NGA
          if (MASS(j) == 0) GO TO 202
      do L = 1, NGA
              if (MASS(L) == 0) GO TO 104
          do i = 1, MASS(j)
              do k = 1, MASS(L)
                      write(1, 221) i, j, k, L, ENASS(i, j, k, L)
221	FORMAT(X, ' ENASS( ', I3, I3, I3, I3, ' ) = ', F10.4)
              end do
          end do
104			continue
          end do
202			continue
  end do


  write(1, 222)
222	FORMAT (/, X, 'PARAMETROS DE VOLUMEN DE ASOCIACI�N (cm3/mol) ', /)

  do j = 1, NGA
          if (MASS(j) == 0) GO TO 301
      do L = 1, NGA
              if (MASS(L) == 0) GO TO 5011
          do i = 1, MASS(j)
              do k = 1, MASS(L)
                  write(1, 223) i, j, k, L, RKASS(i, j, k, L)
223	FORMAT(X, ' RKASS( ', I3, I3, I3, I3, ' ) = ', F10.4)
              end do
          end do
5011			continue
          end do
301		   continue
  end do
  end if
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCC ALFONSINA (basado en el Aparaest) CCCCCCCCCCCCCCCCCCCCC
!c------
    write(IOUT, 699)                                                   
    if (model == 0) write(IOUT, 605)                                    
    if (model == 1) write(IOUT, 627)                                    
 26 continue                                                          
    do 30 i = 1, n_comp                                                      
    Q(i) = 0.D0                                                         
    R(i) = 0.D0                                                         
    do 30 k = 1, NG                                                      
    Q(i) = Q(i)+QT(k, i)                                                 
 30 R(i) = R(i)+RT(k, i)                                                 
    do 40 i = 1, n_comp                                                      
 40 write(6, 606) i, R(i), Q(i)                                          
    if (IOUT == 6) GOTO 42                                             
    do 41 i = 1, n_comp                                                      
 41 write(IOUT, 606) i, R(i), Q(i)                                       
 42 continue      
500 FORMAT(2x, I2, 37X, D15.8, 120X, 3I2)                                                     
501 FORMAT(20I3)                                                      
502 FORMAT(8F10.2)                                                    
503 FORMAT(I3, 2F10.2)                                                 
504 FORMAT(2I3, F10.2)                                                 
603 FORMAT(1X, 10F12.3)                                                
604 FORMAT('  INTERACTION PARAMETERS', /)                              
605 FORMAT(' UNIFAC MOLECULAR R AND Q', /)                             
606 FORMAT(I5, 2F15.4)                                                 
607 FORMAT(' ** WARNING: NUMBER OF SUB GROUPS MUST NOT EXCEED 20 **') 
608 FORMAT(//, ' SUB GROUPS :', 20I3)                                   
609 FORMAT(' MAIN GROUPS:', 20I3)                                      
610 FORMAT(' COMPONENT')                                              
611 FORMAT(6X, I2, 5X, 20I3)                                             
612 FORMAT(' GROUP R- AND Q-VALUES', /)                                
613 FORMAT(1X, I3, 2F10.4)                                              
627 FORMAT(' SPECifIED UNIQUAC R AND Q', /)                            
699 FORMAT(//)                                                        
!c------
1603 format('  ASSOCIATION PARAMETERS', //, 10X, 'k(OH)   :', F7.3, /, 10X, 'E(OH)/k :', F7.1, ' k-1')
!c------
    RETURN                                                            
    end
    
    ! function round(val,idigits0)
    !     use,intrinsic :: iso_fortran_env, only : dp=>real64,sp=>real32
    !     use,intrinsic :: iso_fortran_env, only : int64
    !     implicit none
        
    !     !$@(#) M_verify::round(3f): round val to specified number of significant digits
        
    !     real*8,intent(in)            :: val
    !     integer,intent(in)         :: idigits0
    !     integer(kind=int64)        :: idigits,ipow
    !     real(kind=dp)              :: aval,rnormal
    !     real*8                       :: round
    !        ! make sure a reasonable number of digits has been requested
    !        idigits=max(1,idigits0)
    !        aval=abs(val)
    !     !  select a power that will normalize the number (put it in the range 1 > abs(val) <= 0)
    !        if(aval.ge.1)then
    !           ipow=int(log10(aval)+1)
    !        else
    !           ipow=int(log10(aval))
    !        endif
    !        rnormal=val/(10.0d0**ipow)
    !        if(rnormal.eq.1)then
    !           ipow=ipow+1
    !        endif
    !        !normalize, multiply by 10*idigits to an integer, and so on
    !        round=real(anint(val*10.d0**(idigits-ipow),kind=dp),kind=dp)*10.d0**(ipow-idigits)
    !     end function round