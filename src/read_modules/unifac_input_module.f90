module unifac_input_module
    use iso_fortran_env, only: real64, int32
    
    implicit none
    
    private

    public :: &

        & n_comp, & ! The number of components in system
        & total_groups, & ! Total number of groups
        & total_int_par, & ! Total number of interaction parameters
        & rPar_array, & ! Array with the R parameters.
        & qPar_array ! Array with the r and q parameters. 

  
 
    ! These are module variables that can be used by any program unit that uses 
    ! this module
    integer(kind=int32) :: &
        & n_comp = 0, &
        & total_groups = 0, & 
        & total_int_par = 0 
    
    real(kind = real64), dimension(150) :: &  
        ! Dim = 150 because there is a max number of 150 
        ! R and Q parameters.
        & rPar_array = 0.D0, &
        & qPar_array = 0.D0

    
  
  end module unifac_input_module