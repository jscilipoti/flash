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
    use iso_fortran_env, only: real64, int16, int32
    use fileUnits, only: iout_unit, output_unit
    use outputData, only: output_console, output_lleasoccuzada

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
    DIMENSION X(2)!,ANT(10,3)
    dimension xmj(10),actgam(10),agam(10,4),de(10,10),pe(2,2)          

    ! These integer are meant to be removed
    integer :: &
    & ICALC_common,MODEL_common,IPR_common,IOUT_common,NOVAP_common,ig_common

    ! These integer are meant to be defined
    integer :: & 
    & z_max_index            
    
    ! The number of flash Calculations done
    integer(kind=int32):: flashcalc_loop = 0

                                                                         

    real(kind=real64) :: &

    ! The pressure of the system to be read in the input flash file
    & PP = 0.D0, &

    ! The maximum value among all the molar fractions Z(i) of the components to
    ! be calculated from Z array
    & z_max = 0.D0, &

    ! The sum of all the molar fractions Z(i) to be calculated from Z array.
    & z_sum = 0.D0

!---START------------------------------------------------------------------------

    ! Look for the input_flash file by reading the name in name.dat file.
    ! Then, read all the data in input_flash file.
    call read_input_flash(input_filename)

    IOUT_common = iout
    !ICALC_common = icalc
    !MODEL_common = model
    !IPR_common = ipr
    !NOVAP_common = novap
    !ig_common = ig

    ! Now, since it has read the number of the model, it needs to open the 
    ! appropiate databases 
    !call open_database(model)

    ! Check if "iout = 1" to allow "lleasoccuzada.out" output file.
    if (iout == 1) then
        !iout_unit = get_free_unit()
        open (unit = 1, file = 'lleasoccuzada.OUT', form = 'FORMATTED')
        !open (unit = iout_unit, file = 'lleasoccuzada.OUT', form = 'FORMATTED')
    end if
    
    ! Set the output unit for output.OUT file
    output_unit = get_free_unit()
    open (unit = output_unit, file = 'output.OUT', form = 'FORMATTED')
    
    
    call write_output(1)

                                                      
    
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

    
    mainloop : do while (.true.)
        ! It reads the Temperature and Pressure if the problem is a 
        ! FLASH CALCULATION (0) or it reads only the Temperature if the problem
        ! is a BINODAL CURVE CALCULATION (1) or UNIQUAC PARAMETERS CALCULATION
        ! FROM UNIFAC (2).
        select case (icalc)
            case (0)
                READ(2,*) T, PP
            case (1)
                READ(2,*) T 
            case (2)
                READ(2,*) T 
            case default
                return
        end select                                  

        ! Exit the subroutine if the temperature has not been set or there is no
        ! more inputs for T in the input file: "0,0"
        if (T == 0.) then 
            call close_llecalas()
            return
        end if

        ! If no vapour phase is considered...
        if (PP /= 0.D0  .and.  novap /= 0) then                                 
            do i = 1, N                                                        
                PRAT(i) = DLOG(PP) - ANT(1,i) + ANT(2,i) / (T-273.15 + ANT(3,i))
            end do                                     
        end if
        
        ! Read the composition (Z) of each component (Zi)
        4 READ(2,*) (Z(i), i = 1, N)                                                                                                   
        
        ! It sums the compostion of each component to get Z_sum and the max
        ! value of composition and its index.
        do i = 1, N                                                       
            z_sum = z_sum + z(i)                                                    
            if (Z(i) < z_max) cycle                                        
            z_max = z(i)                                                         
            z_max_index = i                                                            
        end do                                         
        
        call PARAM2                                                       
        
        problem : select case (icalc)
            case (0) ! *** FLASH CALCULATION ***********************************
                do i = 1, N                                                       
                    do j = 1, N                                                       
                        GAM(i, j) = 0.D0                                                     
                        if (j /= i) then                                               
                            call GAMINF(i, j, G)                                                
                            GAM(i,j) = G
                        end if
                    end do
                end do                                                                                                                 
                                                                    
                flashcalc_loop = flashcalc_loop + 1

                do i = 1, N                                                       
                    Z(i) = Z(i) / z_sum
                end do   
                
                ! Write
                write(*, "(///,' * * * FLASH NUMBER',I3,' * * *',//)") flashcalc_loop  
                write(*, "(' TEMPERATURE =',F10.4,' K, PRESSURE =',F7.3,' ATM, FEED =' &
                    & ,F10.2,' MOLES',/,' FEED COMPOSITION (MOLE PERCENT):',/,1X,15(2PF7&
                    & .3))") T, PP, z_sum, (Z(i), i = 1, N)
                if (iout /= 6) write(iout, "(///,' * * * FLASH NUMBER',I3,' * * *',//)") flashcalc_loop                                  
                if (iout /= 6) write(iout, "(' TEMPERATURE =',F10.4,' K, PRESSURE =',F7.3,' ATM, FEED =' &
                    & ,F10.2,' MOLES',/,' FEED COMPOSITION (MOLE PERCENT):',/,1X,15(2PF7 &
                    & .3))") T, PP, z_sum, (Z(i), i = 1, N)
                ! ----

                call unifac(1, Z, AL, DA, PACT)

                SFAS(1) = 1.D0                                                        
                GNUL = 0.D0

                do i = 1, N                                                       
                    XVL(i, 1) = 1.D0                                                       
                    Z(i) = Z(i) + 1.D-20                                                  
                    DLX(i) = DLOG(Z(i))                                                 
                    A(i) = AL(i) + DLX(i)                                                 
                    GNUL = GNUL + Z(i) * AL(i)
                end do                                              
                
                NF = 1
                                                                        
                flash_exit : do while (.true.) 
                ! A loop that ends when "FUN >-1.D-7"
                    call STIG(Y,S)                                                    
                    if (S < -1.D-7) then !if (S > -1.D-7) GOTO 70                                           
                        
                        ! Write
                        write(*, "(/,' SYSTEM IS UNSTABLE, PHASE SPLIT PERFORMED')")                                                      
                        if (iout /= 6) write(iout, "(/,' SYSTEM IS UNSTABLE, PHASE SPLIT PERFORMED')")
                        ! ----                                     
                        
                        do i = 1, N                                                       
                            YVAL(i) = 1.D-5 * Y(i) / Z(i)
                        end do                                                 
                    else                                                       
                        do i = 1, N                                                        
                            YVAL(i) = DLOG(Y(i))
                        end do
                                                                
                        XLAM = 1.D0

                        ! Write
                        if (NF == 1 .and. IPR > 0) write(*, "(//,' DIFFERENTIAL STABILITY TEST FOR FEED MIXTURE:')")                             
                        if (NF > 1 .and. IPR > 0) write(*, "(//,' DIFFERENTIAL STABILITY TEST FOR',I2,'-PHASE SYSTEM')") NF                          
                        if (iout /= 6 .and. NF == 1 .and. IPR > 0) &
                            & write(iout, "(//,' DIFFERENTIAL STABILITY TEST FOR FEED MIXTURE:')")            
                        if (iout /= 6 .and. NF > 1 .and. IPR > 0) &
                            & write(iout, "(//,' DIFFERENTIAL STABILITY TEST FOR',I2,'-PHASE SYSTEM')") NF         
                        ! ----

                        call TMSSJ(30, N, IPR, 15, XLAM, 1.D-12, FUN, YVAL, GRAD, &
                            & XMAT, WORK, 1)     
                        
                        !if (FUN < -1.D-7) GOTO 80                                         
                        if (FUN > -1.D-7) then
                            ! Write 
                            write(*, "(/,' * SYSTEM IS STABLE *',/)")         
                            write(output_unit, *) 1
                                write(output_unit, "(5(2x,f12.8))") (Z(j), j = 1, N)
                                write(output_unit, "(5(2x,f12.8))") (AL(j), j= 1, N)        
                            write(output_unit, *) "SYSTEM IS STABLE"                                                   
                            if (iout /= 6) write(iout, "(/,' * SYSTEM IS STABLE *',/)")
                            ! ----

                            exit 
                        end if

                        ! Write
                        write(*, "(/,' SYSTEM IS UNSTABLE, PHASE SPLIT PERFORMED')")                                                      
                        if (iout /= 6) write(iout, "(/,' SYSTEM IS UNSTABLE, PHASE SPLIT PERFORMED')")
                        ! ----                                     
                        
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
                    
                    ! Write
                    if (IPR > 0) write(*,"(/,' PHASE SPLIT CALCULATION,',I2,' PHASES:')") NF                                      
                    if (iout /= 6 .and. IPR > 0) write(iout,"(/,' PHASE SPLIT CALCULATION,',I2,' PHASES:')") NF
                    ! -----                     
                    
                    call TMSSJ(30, M, IPR, 60, XLAM, 1.D-16, FUN, YVAL, GRAD, &
                    & XMAT, WORK, 2)     
                    
                    NT = NF * N                                                           
                    NB = NT - N                                                           
                    
                    do i = 1, NB                                                     
                        YVAL(NT + 1 - i) = YVAL(NB + 1 - i)
                    end do

                                                                       
                    
                    NVAP = 0                                                            
                    do j = 1, NF                                                     
                        if (IDUM(j) == 1) NVAP = j
                    end do                                           
                    
                    ! Write
                    write(*, "(//,' RESULT OF',I2,'-PHASE CALCULATION:')") NF
                    if (NVAP == 0) write(*, "(' NO VAPOR PHASE')")                                        
                    if (NVAP /= 0) write(*, "(' PHASE',I2,' IS A VAPOR PHASE')") NVAP                                   
                    if (iout /= 6 .and. NVAP == 0) write(iout, "(' NO VAPOR PHASE')")                       
                    if (iout /= 6 .and. NVAP /= 0) write(iout, "(' PHASE',I2,' IS A VAPOR PHASE')") NVAP                  
                    write(*, "(/,'  PHASE FRACTIONS (PERCENT):',4(5X,I3,2PF7.3,5X))") (j, SFAS(j), j = 1, NF)                                   
                    write(*, "(/,'  COMPOSITION  ',10X,4(8X,I3,9X))") (j, j = 1, NF)                                           
                    if (iout /= 6) write(iout, "(//,' RESULT OF',I2,'-PHASE CALCULATION:')") NF                                  
                    if (iout /= 6) write(iout, "(/,'  PHASE FRACTIONS (PERCENT):',4(5X,I3,2PF7.3,5X))")(j, SFAS(j), j = 1, NF)                   
                    if (iout /= 6) write(iout, "(/,'  COMPOSITION  ',10X,4(8X,I3,9X))") (j, j = 1, NF)                          
                    ! -----

                    SUM = 0.D0                                                            
                    
                    do i = 1, N                                                      
                        DLX(i) = XVL(i, NF) * Z(i) / SFAS(NF)                                    
                        SUM = SUM + DLX(i)
                    end do                                                    
                    
                    SUM = DLOG(SUM)                                                     
                    call unifac(1, DLX, A, DA, PACT)                                      
                    do i = 1, N                                                      
                        DLX(i) = DLOG(DLX(i))                                               
                        A(i) = A(i) + DLX(i) - SUM
                    end do                                              
  
                    do j = 1, nf
                        do i = 1, n
                            xmj(i) = xm(i, j)
                        end do

                        call unifac(1, xmj, actgam, de, pe)
                        
                        do i = 1, n
                            agam(i, j) = actgam(i)
                        end do
                    end do
                    
                    write(output_unit,*) NF
                    ! Print the output_unit to be read by an Excel Sheet
                    do i = 1, NF
                        write(output_unit, "(5(2x,f12.8))") (XM(j, i),j = 1, N)
                        write(output_unit, "(5(2x,f12.8))") (agam(j, i),j = 1, N)  
                    end do
                    
                    do i = 1, N                                                      
                        write(*, "('   X(',I2,')            ',5(8X,F12.8))") i, (XM(i, j), j = 1, NF)    ! Composition        
                        write(*, "('  ln(G',i2,')            ',5(8x,f12.8))") i, (agam(i, j), j = 1, nf) ! Ln(gamma)
                    end do
                    if (iout /= 6) then                                          
                        do i = 1, N                                                      
                            write(iout, "('   X(',I2,')            ',5(8X,F12.8))") i, (XM(i, j), j = 1, NF)    
                            write(iout, "('  ln(G',i2,')            ',5(8x,f12.8))") i, (agam(i, j), j = 1, nf)
                        end do
                    end if
                end do flash_exit
            case (1) ! *** BINODAL CURVE CALCULATION ***************************
                call binodal_calc
                call close_llecalas()
                return
            case (2) ! *** CALCULATION OF UNIQUAC PARAMETERS FROM UNIFAC *******                                           
                call uniquac_parcalc
                call close_llecalas()
                return
            case default
                return
        end select problem
    end do mainloop                   

! Close all opened units in this subroutine 
call close_llecalas()

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