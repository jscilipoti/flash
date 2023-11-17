subroutine read_input_unifac
    ! This Fortran subroutine opens a file and reads its contents. 
    ! It first reads the number of parameters and then reads the name of the 
    ! file which contains the data for the flash calculation. 
    ! The subroutine then removes the leading and trailing quotes from the name.
    ! If the number of parameters is 1, the subroutine reads the bases from the 
    ! file using the `get_database_data()` subroutine. Finally, 
    ! the file is closed.

    use iso_fortran_env, only: int16, int32

    use inputData, only: &
        & flash_input_filename, &
        & NTEXT,&
        & icalc, model, IPR, IOUT, NOVAP, ig, ipareq
        

    use unifac_input_module, only: &
        & n_comp, &
        & total_groups, total_int_par, &
        & rPar_array, qPar_array



    use fileUnits, only: flash_input_unit
    
    
    implicit none    
    integer(kind=int16) :: &
    & i, j, k, stat

    ! An external function to get a free unit to open a file.
    integer(kind=int32) :: get_free_unit
    external :: get_free_unit
    
    call read_input_flash("name.dat")

    close(flash_input_unit)
    close(2)
    flash_input_unit = get_free_unit()
    open(unit = flash_input_unit, file = flash_input_filename, status = 'old', &
        & form = 'formatted', action = 'read', iostat = stat)

    ! Read the first line of the flash calc input file where the name is.
    read(flash_input_unit,'(36A2)') NTEXT

    ! Read some info related to the calculation (see InputData module...)   
    read(flash_input_unit,*) ICALC, model, IPR, IOUT, NOVAP, ig, ipareq     

    ! Read the number of components in system
    read(flash_input_unit,*) n_comp

    ! Read the group number and interaction parameters
    read(flash_input_unit,*) total_groups, total_int_par

    ! Read the R and Q parameters
    if (total_groups /=  0) then
        do i = 1, total_groups
            k = 0                                                   
            read(flash_input_unit,*) k, rPar_array(k), qPar_array(k)     
        end do
    end if            

    ! if (total_int_par /=  0) then
    !     intrcnPar_matrix = 0.D0
    !     do i = 1, total_int_par
    !         j = 0
    !         k = 0                                                   
    !         read(2,*) j, k, intrcn_par ! j and k are group numbers  
    !         j = mainsgfunc(j, ipareq)
    !         k = mainsgfunc(k, ipareq)
    !         intrcnPar_matrix(j, k) =  intrcn_par
    !     end do
    ! end if
    
    
end subroutine read_input_unifac
