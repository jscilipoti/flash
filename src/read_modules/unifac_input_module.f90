module unifac_input_module
    use iso_fortran_env, only: real64, int32
    
    implicit none
    
    private

    public :: &

        & n_comp, & ! The number of components in system
        & total_groups, & ! Total number of groups
        & total_int_par, & ! Total number of interaction parameters
        & rPar_array, & ! Array with the R parameters.
        & qPar_array, & ! Array with the r and q parameters. 
        & intrcnPar_matrix, & ! The interaction parameters matrix
        & intrcn_par ! A variable with a intereaction parameter

  
 
    ! These are module variables that can be used by any program unit that uses 
    ! this module
    integer(kind=int32) :: &
        & n_comp = 0, &
        & total_groups = 0, & 
        & total_int_par = 0 

    real(kind = real64) :: &
        & intrcn_par ! A variable with a intereaction parameter
    
    real(kind = real64), dimension(150) :: &  
        ! Dim = 150 because there is a max number of 150 
        ! R and Q parameters.
        & rPar_array = 0.D0, &
        & qPar_array = 0.D0

    ! Here the interaction parameters matrix has "dim=(100,100)" because
    ! there is a maximum of 10 different functional groups in total
        real(kind = real64), dimension(100, 100) :: &
        & intrcnPar_matrix = 0.D0

    
  
  end module unifac_input_module