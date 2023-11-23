subroutine read_input_unifac
    ! This Fortran subroutine gets the info of the kind of input file it needs
    ! to read and its parameters required by the UNIFAC model.


    use iso_fortran_env, only: int16, int32

    use inputData, only: &
        & flash_input_filename, &
        & NTEXT,&
        & icalc, model, IPR, IOUT, NOVAP, ig, ipareq
        

    use unifac_input_module, only: &
        & n_comp, & ! The number of components in system
        & total_groups, & ! Total number of groups
        & total_int_par, & ! Total number of interaction parameters
        & rPar_array, & ! Array with the R parameters.
        & qPar_array, & ! Array with the r and q parameters. 
        & intrcnPar_matrix, & ! The interaction parameters matrix
        & intrcn_par ! A variable with a intereaction parameter



    use fileUnits, only: flash_input_unit
    
    
    implicit none    
    integer(kind=int32) :: &
    & i, j, k, stat

    ! An external function to get a free unit to open a file.
    integer(kind=int32) :: get_free_unit
    external :: get_free_unit
    ! An external function to get the main subgroup id.
    integer(kind=int32) :: mainsgfunc
    external :: mainsgfunc

    ! --- START ----------------------------------------------------------------

    ! It calls "read_input_flash" subroutine in order to get the name of the 
    ! input data file and some necessary information such as the model applied.
    call read_input_flash("name.dat")
    call open_database(model)

    ! Close all the files opened by "read_input_flash" 
    close(flash_input_unit)
    close(2)

    ! Reopen the flash_input file.
    flash_input_unit = get_free_unit()
    open(unit = flash_input_unit, file = flash_input_filename, status = 'old', &
        & form = 'formatted', action = 'read', iostat = stat)
    if (stat /= 0) then ! check for errors
        print *, 'Error opening file ', flash_input_filename
        !ERROR STOP 'Error opening file.'
    end if

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

    if (total_int_par /=  0) then
        do i = 1, total_int_par
            j = 0
            k = 0                                                   
            read(flash_input_unit,*) j, k, intrcn_par ! j & k are group numbers
            j = mainsgfunc(j, ipareq)
            k = mainsgfunc(k, ipareq)
            intrcnPar_matrix(j, k) =  intrcn_par
        end do
    end if
    
    
end subroutine read_input_unifac
