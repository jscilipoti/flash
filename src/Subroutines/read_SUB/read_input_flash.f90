subroutine read_input_flash(input_filename)
    ! This Fortran subroutine opens a file and reads its contents. 
    ! It first reads the number of parameters and then reads the name of the 
    ! file which contains the data for the flash calculation. 
    ! The subroutine then removes the leading and trailing quotes from the name.
    ! If the number of parameters is 1, the subroutine reads the bases from the 
    ! file using the `get_database_data()` subroutine. Finally, 
    ! the file is closed.

    !use CUFAC, only: NKK, NGG, Pxx, Txx
    !use flash, only: P, T, Z
    use iso_fortran_env, only: int16, int8, int32

    use inputData, only:&
        & flash_input_filename, name_maxlen,&
        & icalc, model, IPR, IOUT, NOVAP, ig, ipareq,&
        & NTEXT


    use fileUnits, only: flash_input_unit
    
    
    implicit none    
    
    integer(kind=int16) :: i, j, k, stat
    ! The number of variables to get from the input_file
    integer(kind=int16), parameter :: file_vars = 2
    ! A variable to be read from the input_file in order to allow or not further
    ! parameter seach.
    integer(kind=int8) :: search_parameters_flag = 0
    ! An external function to get a free unit to open a file.
    integer(kind=int32) :: get_free_unit
    external :: get_free_unit
    
    ! The name of the input_file to be read. It must contain the name of the 
    ! flash input file and the value of search_parameters_flag.
    character(len=*), intent(in) :: input_filename
    ! The raw data from the input_file
    character(len=name_maxlen), dimension(file_vars) :: file_data
    
    ! Open the file for reading
    call open_textfile(input_filename, file_data, file_vars, name_maxlen)

     ! Read the number of parameters and the name from the file
    read(file_data(1),'(I1)') search_parameters_flag
    flash_input_filename = file_data(2)

    ! Remove the leading and trailing quotes from the name
    flash_input_filename = &
        & flash_input_filename(2:len_trim(flash_input_filename)-1)

    ! If the number of search_parameters_flag is 1, read the bases from the file
    if (search_parameters_flag == 1) then
        call get_database_data()
        stop
    end if
    
    !Old way to open flashInput file
    open(unit=2, file=flash_input_filename, status='old', form='formatted',&
        &action='read', iostat=stat)
    if (stat /= 0) then ! check for errors
        print *, 'Error opening file ', flash_input_filename
        ERROR STOP 'Error opening file.'
    end if
    !New way to open flashInput file
    
    flash_input_unit = get_free_unit()
    open(unit=flash_input_unit, file=flash_input_filename, status='old',&
        &form='formatted', action='read', iostat=stat)
    if (stat /= 0) then ! check for errors
        print *, 'Error opening file ', flash_input_filename
        !ERROR STOP 'Error opening file.'
    end if

    ! Read the first line of the flash calc input file where the name is.
    read(2,'(36A2)') NTEXT

    ! Read some info related to the calculation (see InputData module...)   
    read(2,*) ICALC, model, IPR, IOUT, NOVAP, ig, ipareq     
    
    !Disabled since there are more errors below
    call open_database(model)


end subroutine read_input_flash
