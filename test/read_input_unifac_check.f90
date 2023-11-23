program check
    use do_tests

    use iso_fortran_env, only: int32

    use inputData, only: &
        & flash_input_filename, &
        & icalc, model, IPR, IOUT, NOVAP, ig, ipareq, &
        & NTEXT

        use unifac_input_module, only: &
        & n_comp, & ! The number of components in system
        & total_groups, & ! Total number of groups
        & total_int_par, & ! Total number of interaction parameters
        & rPar_array, & ! Array with the R parameters.
        & qPar_array, & ! Array with the r and q parameters. 
        & intrcnPar_matrix, & ! The interaction parameters matrix
        & sk_idGroups_matrix, sk_nGroups_matrix, &
        & skeletal_matrix
    
    use stdlib_ansi, only : & 
    & fg_color_green, fg_color_red, fg_color_yellow, fg_color_blue, & 
    & style_bold, style_reset, operator(//), operator(+)
    
    implicit none
    
    logical :: is_OK = .true.
    logical :: print_to_screen = .true.
    integer(kind=int32) :: error_count = 0
    
    character(len=*), parameter :: PATH = ""
    character(len=*), parameter :: flashtype = ""
    
    integer(kind=int32) :: i, j, k
    integer, dimension(8)  :: group_indexs 

    ! An external function to get the main subgroup id.
    integer(kind=int32) :: mainsgfunc
    external :: mainsgfunc

    

    print *,""
    print *, fg_color_blue + style_bold // test_run // style_reset //"read_input_unifac_check"
    
    if (.true.) then
        if (pause_test) pause
        
        call read_input_unifac("name.dat")

        if (print_to_screen) then
            print '(36A2)', NTEXT
            print '(7(I1,","))', ICALC, model, IPR, IOUT, NOVAP, ig, ipareq
            print '(I1)', n_comp
            print '(I1,",",I2)', total_groups, total_int_par
            ! Print RQ parameters
            k = 0 
            do i = 1, 150
                k = k + 1
                if (rPar_array(k) /= 0) then
                    print '(I2,2(",",F5.3))', k, rPar_array(k), qPar_array(k)
                end if     
            end do
            ! Print interaction parameters
            group_indexs = (/1,2,3,15,22,17,38,43/)
            do i = 1, total_groups
                do j = 1, total_groups
                    print '(2(I2,","),F10.3)', group_indexs(i), &
                        & group_indexs(j), &
                        & intrcnPar_matrix(mainsgfunc(group_indexs(i), ipareq), mainsgfunc(group_indexs(j), ipareq)) 
                end do
            end do
            ! Print Skeletal matrix
            do i = 1, n_comp
                !do j = 1, size(skeletal_matrix(1, :, 1))
                    ! Save the IDs of each group of one compound.
                    j=1
                    print "(20(I3,','))", &
                    & (sk_idGroups_matrix(j, i), sk_nGroups_matrix(j, i), &
                    j = 1, size(skeletal_matrix(1, :, 1))) 
                    ! Save the numbers of each group in one compound.
                !end do
            end do

        end if 
        
        !Check if everything went OK
        ! do i = 1, lleasoccuzada_maxlines
        !     if (trim(lleasoccuzada_new(i)) /= trim(lleasoccuzada_old(i))) then
        !         print *, fg_color_red + style_bold // test_error // style_reset
        !         print *, "At line: ",i, " of file: lleasoccuzada.OUT"
        !         print *,&
        !             trim(lleasoccuzada_new(i)) /= trim(lleasoccuzada_old(i))
        !         print '(3A)',&
        !             &trim(lleasoccuzada_new(i)),"//",trim(lleasoccuzada_old(i))
        !         is_OK = .false.
        !         error_count = error_count + 1
        !         !ERROR STOP ""
        !     end if
        ! end do


        ! Print the result of the test:
        ! Errors:
        if (error_count >= 1) then
            ERROR STOP "TOO MANY ERRORS"
        end if
        if (error_count < 0) then
            ERROR STOP "LESS ERRORS!!!"
        end if
        ! Test passed:
        if (is_OK) then
            print *, fg_color_green + style_bold // test_ok // style_reset
        end if    
    else
        ! Test was previusly disabled: 
        print *, fg_color_yellow + style_bold // test_disabled // style_reset
    end if

end program check
