module CUFAC
    use iso_fortran_env, only: real64, int32
    
    implicit none
    
    private

    public :: NKK,NGG,Pxx,Txx
  
 
    ! These are module variables that can be used by any program unit that uses 
    ! this module
    integer(kind=int32) :: NKK = 0.0, NGG = 0.0
    real(kind=real64) :: Txx = 0.D0
    real(kind=real64), dimension(10,10) :: Pxx = 0.D0
    
  
  end module CUFAC