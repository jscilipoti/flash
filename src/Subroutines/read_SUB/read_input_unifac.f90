subroutine read_input_unifac(input_filename)
    ! This Fortran subroutine gets the info of the kind of input file it needs
    ! to read and its parameters required by the UNIFAC model.


    use iso_fortran_env, only: int32, real64

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
        & sk_idGroups_matrix, sk_nGroups_matrix, &
        & skeletal_matrix



    use fileUnits, only: flash_input_unit
    
    
    implicit none    
    common/asoc/nktt, sk_idGroups_matrix_common, sk_nGroups_matrix_common
    COMMON/CUFAC/n_comp_common, NG, P(10, 10), T   

    integer(kind=int32) :: &
    & i, j, k, stat

    ! An external function to get a free unit to open a file.
    integer(kind=int32) :: get_free_unit
    external :: get_free_unit
    ! An external function to get the main subgroup id.
    integer(kind=int32) :: mainsgfunc
    external :: mainsgfunc

    real(kind = real64) :: &
        & intrcn_par ! A variable with a interaction parameter

    ! The name of the input_file to be read. It must contain the name of the 
    ! flash input file and the value of search_parameters_flag.
    character(len=*), intent(in) :: input_filename

    ! COMMON USED (REMOVE)

    integer(kind=int32) :: &
        & n_comp_common = 0
    integer(kind = int32), dimension(20,12) :: &
        & sk_idGroups_matrix_common, &
        & sk_nGroups_matrix_common

    ! COMMON NOT USED (REMOVE)
        
    integer(kind = int32) :: nktt, ng
    real(kind = real64) :: P, T




    ! --- START ----------------------------------------------------------------

    ! It calls "read_input_flash" subroutine in order to get the name of the 
    ! input data file and some necessary information such as the model applied.
    close(flash_input_unit)
    close(2)
    call read_input_flash(input_filename)
    call open_database(model)

    ! Close all the files opened by "read_input_flash" 
    close(flash_input_unit)
    !close(2)

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
    ! Read the interaction parameters
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
    ! Read the skeletal matrix. It shows the information related to how many
    ! functional groups and which ones are present in each compound.
    do i = 1, n_comp
        read(flash_input_unit,*) &
            & (skeletal_matrix(i, j, 1), skeletal_matrix(i, j, 2), &
            & j = 1, size(skeletal_matrix(1, :, 1))) 
        
        do j = 1, size(skeletal_matrix(1, :, 1))
            ! Save the IDs of each group of one compound.
            sk_idGroups_matrix(j, i) = skeletal_matrix(i, j, 1)
            ! Save the numbers of each group in one compound.
            sk_nGroups_matrix(j, i) = skeletal_matrix(i, j, 2)
        end do
    end do
    
    ! Save Common variables (REMOVE)
    n_comp_common = n_comp
    sk_idGroups_matrix_common = sk_idGroups_matrix
    sk_nGroups_matrix_common = sk_nGroups_matrix


end subroutine read_input_unifac
