! *****************************************************************************
! *   Subroutine: llecalas(Tf, Pf, Zf)                                        *
! *****************************************************************************
!
subroutine llecalas!(Tf, Pf, Zf) 
!
! *****************************************************************************
! *                                                                           *
! *  SUBROUTINE   L L E C A L A S                                                *
! *  (asociacion incorporada para el calculo de flash y                       *
! *  curva binodal (icalc 0 y 1)                                              *
! *                                                                           *
! *                           ASOCIACION CRUZADA                              *
! *                          VERSION GENERALIZADA                             *
! *                              Junio 2023                                   *
! *   Modificada por:                                                         *
! *                        ALFONSINA ESTER ANDREATTA                          *
! *                          JOSE ANTONIO SCILIPOTI                           *
! *                           JUAN PABLO ROVEZZI                              *
! *   Revisada en Octubre del 2007 en el chequeo de estabilidad               *
! *                                                                           *
! *   Basada en las simplificaciones de los papers:                           *
! *      # Michelsen, et al. (Fluid Phase Equilibria, 180(2001)165-174 )      *
! *      # Tan, et al.  (Ind. Eng. Chem. Res, 2004, 43, 203-208).             *
! *                                                                           *
! *   Esto permitio  que todos los casos particulares de asociacion se puedan *
! *   simplificar a un unico calculo.                                         *
! *                                                                           *
! *   Valido para un maximo numero grupo asociativo de 12.                    *
! *   Con la implementacion en el calculo de la fraccion no asociada en el    *
! *   componente puro por  el metodo iterativo aqui implementado se permite   *
! *   que una molecula tenga mas de un grupo asociativo (14/07/06).           *
! *                                                                           * 
! *   El calculo se limita a que el numero maximo de sitios sea dos           * 
! *   (por razones matematicas).                                              *
! *                                                                           *                                                   
! *****************************************************************************  
! *                        DATE: 24/3/1982 /TJ                                *
! *****************************************************************************
! 
! Tf: Temperature
! Pf: Pression
! Zf: Composition for each component
!  
! Required modules:
    use InputData
    use iso_fortran_env, only: real64, int16
    use fileUnits, only: iout_unit, output_unit

    ! Implicit statements need to be removed. It must use 'implicit none'
    IMPLICIT real(kind=real64) (A-H,O-Z)                                          

    ! External functions
    !EXTERNAL STABIL,GMIX ! Not used
    EXTERNAL FUNC                                         

    ! Common blocks that are meant to be removed
                            
    COMMON/CGIBBS/NF,z_max_index,GNUL,Z(10),A(10),XVL(10,4),SFAS(4),GAM(10,10),AL(10),&
    & DA(10,10),XM(10,4)                                       
    COMMON/CUFAC/N,NG,P(10,10),T                                      
    COMMON/CY/Y13,Y21,STEP                                            
    COMMON/CA/XC(5),GE(5,2),GC(5,2)                                   
                                                    
    COMMON/CQT/QT(10,10),Q(10),R(10)                                  
    COMMON/CACT/Y1(10),Y2(10),ACT1(10),ACT2(10),DACT1(10,10),DACT2(10,10),PACT(2,2)                                                                             
    common/nga/nga,mass(12)




    COMMON/CMODEL/model_common 
    COMMON/CIPR/ipr_common 
    COMMON/COUT/iout_common
    COMMON/CVAP/novap_common,NDUM,IDUM(4),PRAT(10) 
    common/ig/ig_common

    !   
    DIMENSION DLX(10),YVAL(30),Y(10),GRAD(30),XMAT(30,30),WORK(30,5)  
    !DIMENSION NTEXT(36),
    DIMENSION X(2),ANT(10,3)
    dimension xmj(10),actgam(10),agam(10,4),de(10,10),pe(2,2)          

    ! These integer are meant to be removed
    integer :: &
    & ICALC_common,MODEL_common,IPR_common,IOUT_common,NOVAP_common,ig_common

    ! These integer are meant to be defined
    integer :: & 
    & z_max_index            
    
    integer(kind=int16):: NN = 0

                                                                         

    real(kind=real64) :: &

    ! The pressure of the system to be read in the input flash file
    & PP = 0.D0, &

    ! Apparently, this is a value to test what happends at this T
    & T1 = 0.D0, & !T1=0.

    ! The maximum value among all the molar fractions Z(i) of the components to
    ! be calculated from Z array
    & z_max = 0.D0, &

    ! The sum of all the molar fractions Z(i) to be calculated from Z array.
    & z_sum = 0.D0

!---START------------------------------------------------------------------------

    ! Look for the input_flash file by reading the name in name.dat file.
    ! Then, read all the data in input_flash file.
    call read_input_flash(input_filename)

    ICALC_common = icalc
    MODEL_common = model
    IPR_common = ipr
    IOUT_common = iout
    NOVAP_common = novap
    ig_common = ig

    ! Now, since it has read the number of the model, it needs to open the 
    ! appropiate databases 
    call open_database(model)

    ! Check if "iout = 1" to allow "lleasoccuzada.out" output file.
    if (iout == 1) then
        open (unit = 1, file = 'lleasoccuzada.OUT', form = 'FORMATTED')
    end if
    
    ! Set the output unit for output.OUT file
    output_unit = get_free_unit()
    open (unit = output_unit, file = 'output.OUT', form = 'FORMATTED')
    
    
    
    ! Write to screen
    write(*, 608)                                                      
    write(*, 628) iout                                                 
    write(*, 610)                                                      
    
    if (iout == 0) iout = 6                                              
    if (icalc == 0) write(*, 620)                                       
    if (icalc == 1) write(*, 621)                                       
    if (icalc == 2) write(*, 622)                                       
    if (novap /= 0) write(*, 629)                                       
    if (model == 0) write(*, 624)                                       
    if (model == 1) write(*, 625)                                       
    
    write(*,623) NTEXT                                                
    
    !if (iout == 6) GOTO 5
    if (iout /= 6) then                                              
    
        ! write to iout
        write(iout, 608)                                                   
        write(iout, 610)                                                   
        
        if (icalc == 0) write(iout,620)                                    
        if (icalc == 1) write(iout,621)                                    
        if (icalc == 2) write(iout,622)                                    
        if (novap /= 0) write(iout,629)                                    
        if (model == 0) write(iout,624)                                    
        if (model == 1) write(iout,625)                                    
        write(iout,623) NTEXT                                             
    
    end if
    !5 CONTINUE                                                          
    
    call PARIN2                                                       
    
    !-----------------------------------------------------------------------
    !Testear y luego modificar
    ! If no vapour phase is considered...
    if (novap /= 0 ) then                                             
        do 6 j=1,N                                                                                 
    6       READ(2,*) (ANT(K,j),K=1,3)                                        
        do 7 j=1,N                                                        
            ANT(1,j)=2.302585*(ANT(1,j)-2.880814)                             
    7       ANT(2,j)=2.302585*ANT(2,j)                                        
    endif                                                          
    !-----------------------------------------------------------------------
    
    !T1 = 0.                                                            
    !NN = 0                                                              

    
    do while (.true.)
        ! It checks whether the problem is a FLASH CALCULATION (0) or a 
        ! BINODAL CURVE CALCULATION OR CALCULATION OF UNIQUAC PARAMETERS 
        ! FROM UNIFAC. After that it reads the Temperature and Pressure.
        10  if (icalc == 0) then 
                READ(2,*) T, PP
            else !if (icalc .GE. 1) READ(2,*) T
                READ(2,*) T 
            end if                                   

        
        ! Exit the subroutine if the temperature has not been set 
        if (T == 0.) then 
            call close_llecalas() !GOTO 10000
            return
        end if

        ! If no vapour phase is considered...
        if (PP /= 0.D0  .and.  novap /= 0) then                                 
            do i=1,N                                                        
                PRAT(i)=DLOG(PP)-ANT(1,i)+ANT(2,i)/(T-273.15+ANT(3,i))
            end do                                     
        end if
        
        ! Read the composition (Z) of each component (Zi)
        4 READ(2,*) (Z(i),i=1,N)                                            
        
        !z_sum = 0.D0                                                       
        !z_max = 0.D0                                                          
        
        ! It sums the compostion of each component to get Z_sum and the max
        ! value of composition and its index.
        do i = 1, N !15 i = 1, N                                                       
            z_sum = z_sum + z(i)                                                    
            if (Z(i) < z_max) cycle !GOTO 15                                          
            z_max = z(i)                                                         
            z_max_index = i                                                            
        !15 CONTINUE
        end do

        ! The following line has no sense since, as far as I know, T1
        ! had been initialized at the beggining as T1 = 0. Also it cuts some
        ! code related to the Flash Calculation
        !if (T == T1) GOTO 30
        if (T /= T1) then                                             
        
            call PARAM2                                                       
            
            !** BINODAL CURVE CALCULATION **********************************
            if (icalc == 1) then !if (icalc /= 1) GOTO 16                                            
            
                if (N /= 2 .and. N /= 3) write(*, 616)                                
                if (iout /= 6 .and. N /= 2 .and. N /= 3) write(iout, 616)               
                Y13 = Z(1)                                                          
                Y21 = Z(2)                                                          
                
                write(*,633) T                                                    
                
                if (iout /= 6) write(iout, 633) T                                   
                if (N == 3) then !if (N == 3) GOTO 12                                                
                    !12 STEP=Z(3) / 100.D0
                    STEP = Z(3) / 100.D0

                    if (STEP == 0.D0) STEP = 0.02D0                                         
                    call BINOD                                                        
                    
                    call close_llecalas() !GOTO 10000
                    return
                end if

                call SOLBIN                                                       
                
                call close_llecalas() !GOTO 10000
                return
                                                                    
        

            !16 CONTINUE
            end if
                                                                    
            !** CALCULATION OF UNIQUAC PARAMETERS FROM UNIFAC **************
            if (icalc == 2) then !if (icalc /= 2) GOTO 19                                            
                if (N /= 2) write(*, 616)                                           
                if (iout /= 6 .and. N /= 2) write(iout, 616)                          
                XC(1) = 0.D0                                                          
                XC(2) = 0.2D0                                                        
                XC(3) = 0.5D0                                                        
                XC(4) = 0.8D0                                                        
                XC(5) = 1.D0                                                        
                
                do 17 k = 1, 5                                                       
                    Y(1) = XC(k)                                                        
                    Y(2) = 1.D0 - XC(k)                                                   
                    call unifac(1, Y, ACT1, DACT1, PACT)                                  
                    GE(K, 1) = ACT1(1)                                                   
                17  GE(K, 2) = ACT1(2)

                READ(2, *) R(1), Q(1)                                               
                READ(2, *) R(2), Q(2)                                               
                                            
                write(*, 627)                                                      
                
                do 14 i = 1, 2                                                       
                14   write(*, 626) i, R(i), Q(i)                                          
                
                if (iout /= 6) then !if (iout == 6) GOTO 13                                             
                    write(iout, 627)                                                   
                    do 11 i = 1, 2                                                       
                11  write(iout, 626) i, R(i), Q(i)                                        
                !13 CONTINUE
                end if

                X(1) = Z(1) / 300.D0                                                  
                X(2) = Z(2) / 300.D0                                                  
                do 18 i = 1, 2                                                       
                    do 18 j = 1, 2                                                       
                        QT(i, j) = 0.D0                                                        
                18       P(i, j) = 0.D0                                                         
                QT(1, 1) = Q(1)                                                      
                QT(2, 2) = Q(2)                                                      
                NK = 2                                                              
                NG = 2                                                              
                XLAMB = 1.D0                                                          
                call MARQ(FUNC, 2, 10, X, XLAMB, 3.D0, 1.D-7, 99)                        
                write(*, 633) T                                                    
                if (iout /= 6) write(iout, 633) T                                   
                write(*, 617) P(1, 2), P(2, 1)     
                if (IPR == 1) write(*, 618)                                         
                do 21 L = 1, 5                                                       
                    do 21 i = 1, 2                                                       
                        GE(L,i) = DEXP(GE(L, i))                                             
            21       GC(L, i) = DEXP(GC(L, i))                                             
                if (IPR == 1) write(*, 619) ((GE(L, i), L = 1, 5), i = 1, 2)                 
                if (IPR == 1) write(*, 619) ((GC(L, i), L = 1, 5), i = 1, 2)                 
                
                if (iout /= 6) then !if (iout == 6) GOTO 22                                             
                    write(iout, 617) P(1, 2), P(2, 1)                                     
                    if (IPR == 1) &
                    & write(iout, 618)                                      
                    if (IPR == 1) &
                    & write(iout, 619) ((GE(L, i), L = 1, 5), i = 1, 2)              
                    if (IPR == 1) & 
                    & write(iout, 619) ((GC(L, i), L = 1, 5), i = 1, 2)              
                !22 CONTINUE
                end if

                call close_llecalas() !GOTO 10000
                return

            19 CONTINUE                                                          
            end if
            
            !** FLASH CALCULATION ******************************************
            do i = 1, N                                                       
                do j = 1, N                                                       
                    GAM(i, j) = 0.D0                                                     
                    if (j /= i) then !if (j == i) GOTO 20                                                
                        call GAMINF(i, j, G)                                                
                        GAM(i,j) = G
                    end if
                end do
            enddo                                                        
            !20  CONTINUE                                                          
        end if
    30  T1 = T                                                              
        NN = NN + 1

        write(*, 602) NN                                                   
        
        do 35 i = 1, N                                                       
    35  Z(i) = Z(i) / z_sum                                                    
        
        write(*, 605) T, PP, z_sum, (Z(i), i = 1, N)

        if (iout /= 6) write(iout, 602) NN                                  
        if (iout /= 6) write(iout, 605) T, PP, z_sum, (Z(i), i = 1, N)              
        call unifac(1, Z, AL, DA, PACT)                                       
        SFAS(1) = 1.D0                                                        
        GNUL = 0.D0                                                           
        do 40 i = 1, N                                                       
            XVL(i, 1) = 1.D0                                                       
            Z(i) = Z(i) + 1.D-20                                                  
            DLX(i) = DLOG(Z(i))                                                 
            A(i) = AL(i) + DLX(i)                                                 
    40  GNUL = GNUL + Z(i) * AL(i)                                              
        NF = 1
                                                                    
    do while (.true.) ! A loop that ends when "FUN > -1.D-7"
        50  call STIG(Y,S)                                                    
            if (S < -1.D-7) then !if (S > -1.D-7) GOTO 70                                           
                write(*, 603)                                                      
                if (iout /= 6) write(iout, 603)                                     
                do i = 1, N                                                       
                    YVAL(i) = 1.D-5 * Y(i) / Z(i)
                end do                                                                                                    
                !GOTO 100 !ELIMINAR GOTO 
            !end if
            else                                                       
        
            !70 do 75 i = 1, N
            
                do i = 1, N                                                        
                    YVAL(i) = DLOG(Y(i))
                end do
                                                        
                XLAM = 1.                                                           
                if (NF == 1 .and. IPR > 0) write(*, 606)                             
                if (NF > 1 .and. IPR > 0) write(*, 609) NF                          
                if (iout /= 6 .and. NF == 1 .and. IPR > 0) &
                    & write(iout, 606)            
                if (iout /= 6 .and. NF > 1 .and. IPR > 0) &
                    & write(iout, 609) NF         
                
                call TMSSJ(30, N, IPR, 15, XLAM, 1.D-12, FUN, YVAL, GRAD, &
                    & XMAT, WORK, 1)     
                
                !if (FUN < -1.D-7) GOTO 80                                         
                if (FUN > -1.D-7) then 
                    write(*, 604)         
                    write(output_unit, *) 1
                        write(output_unit, 2613) (Z(j), j = 1, N)
                        write(output_unit, 2613) (AL(j), j= 1, N)        
                    write(output_unit, *) "SYSTEM IS STABLE"                                                   

                    ! write(7,46) T,  (xM(l,1),l=1,N)     !Alfonsina
                    ! write(7,46) T,  (xM(l,2),l=1,N)     !Alfonsina
                    ! write(7,*)                          !Alfonsina
            
                    if (iout /= 6) write(iout, 604)                                     
                    !GOTO 10 
                    exit 
                end if

            80  write(*, 603)                                                      
                if (iout /= 6) write(iout, 603)                                     
                
                do i = 1, N                                                       
                    YVAL(i) = 1.D-5 * DEXP(YVAL(i)) / Z(i)
                end do
            end if
                                            
        100 NF = NF + 1                                                           
            do i = 1, N                                                      
                if (YVAL(i) > 1.D0) then 
                    do j = 1, N                                                      
                        YVAL(j) = YVAL(j) / 10.D0
                    end do                                               
                end if
            end do

            SFAS(NF) = 1.D0                                                       
            XLAM = 0.2D0                                                           
            
            if (NF == 2) XLAM = 0.5D0                                               
            M = (NF - 1) * N                                                        
            
            if (IPR > 0) write(*,607) NF                                      
            
            if (iout /= 6 .and. IPR > 0) write(iout,607) NF                     
            
            call TMSSJ(30, M, IPR, 60, XLAM, 1.D-16, FUN, YVAL, GRAD, &
            & XMAT, WORK, 2)     
            
            NT = NF * N                                                           
            NB = NT - N                                                           
            
            do i = 1, NB                                                     
                YVAL(NT + 1 - i) = YVAL(NB + 1 - i)
            end do

            write(*, 614) NF                                                   
            
            NVAP = 0                                                            
            do j = 1, NF                                                     
                if (IDUM(j) == 1) NVAP = j
            end do                                           
                                                                
            if (NVAP == 0) write(*, 630)                                        
            if (NVAP /= 0) write(*, 631) NVAP                                   
            if (iout /= 6 .and. NVAP == 0) write(iout, 630)                       
            if (iout /= 6 .and. NVAP /= 0) write(iout, 631) NVAP                  
            write(*, 611) (j, SFAS(j), j = 1, NF)                                   
            write(*, 612) (j, j = 1, NF)                                           
            if (iout /= 6) write(iout, 614) NF                                  
            if (iout /= 6) write(iout, 611)(j, SFAS(j), j = 1, NF)                   
            if (iout /= 6) write(iout, 612) (j, j = 1, NF)                          
            
            SUM = 0.D0                                                            
            
            do i = 1, N                                                      
                DLX(i) = XVL(i, NF) * Z(i) / SFAS(NF)                                    
                SUM = SUM + DLX(i)
            end do                                                    
            
            SUM = DLOG(SUM)                                                     
            call unifac(1, DLX, A, DA, PACT)                                      
            do 120 i = 1, N                                                      
                DLX(i) = DLOG(DLX(i))                                               
        120   A(i) = A(i) + DLX(i) - SUM                                              
        !c-----
            do 1130 j = 1, nf
                do 1131 i = 1, n
        1131       xmj(i) = xm(i, j)
                call unifac(1, xmj, actgam, de, pe)
                do 1132 i = 1, n
        1132       agam(i, j) = actgam(i)
        1130 continue
            write(output_unit,*) NF
            ! Print the output_unit to be read by an Excel Sheet
            do i = 1, NF
                write(output_unit, 2613) (XM(j, i),j = 1, N)
                write(output_unit, 2613) (agam(j, i),j = 1, N)  
            end do
            
            do i = 1, N                                                      
                write(*, 613) i, (XM(i, j), j = 1, NF)    ! Composition        
                write(*, 1613) i, (agam(i, j), j = 1, nf) ! Ln(gamma)
            end do
            !if (iout == 6) GOTO 132
            if (iout /= 6) then                                          
                do i=1,N                                                      
                    write(iout, 613) i, (XM(i, j), j = 1, NF)    
                    write(iout, 1613) i, (agam(i, j), j = 1, nf)
                end do
                                                                    
            ! 46 FORMAT( &
            !   & 2X, F12.2, 8X,F12.6, &
            !   & 8X, F12.6, 8X, F12.6, &
            !   & 8X, F12.6, 8X, F12.6, &
            !   & 8X, F12.6, 8X, F12.6, &
            !   & 8X, F12.6, 8X, F12.6, &
            !   & 8X, F12.6, 8X, F12.6) ! Alfonsina
            end if
            !132 CONTINUE
            !GOTO 50
    end do
end do 

!-------------------------------------------------------------------------------
!
!     FORMATS USED IN THIS PROGRAM:
!
!-------------------------------------------------------------------------------
501 FORMAT(36A2)                                                      
502 FORMAT(8F10.2)                                                    
503 FORMAT(20I3)                                                      
602 FORMAT(///,' * * * FLASH NUMBER',I3,' * * *',//)                  
603 FORMAT(/,' SYSTEM IS UNSTABLE, PHASE SPLIT PERFORMED')            
604 FORMAT(/,' * SYSTEM IS STABLE *',/)                               
605 FORMAT(' TEMPERATURE =',F10.4,' K, PRESSURE =',F7.3,' ATM, FEED ='&
    ,F10.2,' MOLES',/,' FEED COMPOSITION (MOLE PERCENT):',/,1X,15(2PF7&
    .3))                                                              
606 FORMAT(//,' DIFFERENTIAL STABILITY TEST FOR FEED MIXTURE:')       
607 FORMAT(/,' PHASE SPLIT CALCULATION,',I2,' PHASES:')               
608 FORMAT(1H1)                                                       
609 FORMAT(//,' DIFFERENTIAL STABILITY TEST FOR',I2,'-PHASE SYSTEM')  
610 FORMAT(///)                                                       
611 FORMAT(/,'  PHASE FRACTIONS (PERCENT):',4(5X,I3,2PF7.3,5X))       
612 FORMAT(/,'  COMPOSITION  ',10X,4(8X,I3,9X))                       
613 FORMAT('   X(',I2,')            ',5(8X,F12.8))                    
!c-----
1613 format('  ln(G',i2,')            ',5(8x,f12.8))
2613 format(5(2x,f12.8))

!c-----
614 FORMAT(//,' RESULT OF',I2,'-PHASE CALCULATION:')                  
616 FORMAT(//,' * WRONG INPUT SPECIFICATION *',//)                    
617 FORMAT(///,' ** UNIQUAC PARAMETERS FROM UNIFAC **',//,5X,'A12/R =  ',F12.3,' K ,  A21/R = ',F12.3,' K',///)                         
618 FORMAT(//,' ** COMPARISON OF ACTIVITIES CALCULATED BY UNIFAC AND UNIQUAC, RESPECTIVELY **'//)                                       
619 FORMAT(10F12.5)                                                   
620 FORMAT(' **** FLASH CALCULATION ****')                            
621 FORMAT(' **** BINODAL CURVE CALCULATION ****',//)                 
622 FORMAT(' **** CALCULATION OF UNIQUAC PARAMETERS FROM UNIFAC **** ',//)                                                              
623 FORMAT(1X,'COMPONENTS : ',40A2,//)                                
624 FORMAT(' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIFAC'//)     
625 FORMAT(' MODEL USED FOR LIQUID PHASE NON-IDEALITY: UNIQUAC'//)    
626 FORMAT(I5,2F15.4)                                                 
627 FORMAT(//,' SPECIFIED UNIQUAC R AND Q',/)                         
628 FORMAT(/,' iout = ',I2,/' if iout = 0: OUTPUT ONLY ON UNIT 6',/,  ' if iout = 1: OUTPUT ON BOTH UNIT 6 AND 1')                      
629 FORMAT(/,' VAPOR PHASE INCLUDED IN FLASH-CALCULATIONS',//)        
630 FORMAT(' NO VAPOR PHASE')                                         
631 FORMAT(' PHASE',I2,' IS A VAPOR PHASE')                           
633 FORMAT(///,'   TEMPERATURE =',F10.2,' DEG K')                     

! Close all opened units in this subroutine 
call close_llecalas()

! 10000 CLOSE (UNIT=2)
!     if (iout == 1) CLOSE (UNIT = 1)
!     close (unit=3)
! !c      call salida(name)
!     !STOP

return

end SUBROUTINE llecalas   
                    
!---END-------------------------------------------------------------------------

! Subroutine close_llecalas: This subroutine replaces the original 
! "goto 10000" statement. It closes Units: 1 (lleasoccuzada.OUT),
! 2 (llecalas.dat) and output_unit (output.out).
subroutine close_llecalas()
    use iso_fortran_env, only: int8
    use fileUnits, only: output_unit
    use inputData, only: iout

    implicit none
    
    if (iout == 1) close (unit = 1)
    close (unit = 2)
    close (unit = output_unit)
    return

end subroutine close_llecalas

!-------------------------------------------------------------------------------