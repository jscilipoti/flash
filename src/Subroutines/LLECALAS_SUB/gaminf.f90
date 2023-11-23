SUBROUTINE gaminf(n1,n2,gam)
    use iso_fortran_env, only: real64, int32                                     
    
    implicit none                                          
    
    integer(kind = int32) :: i, NG
    real(kind=real64), dimension(10,10) :: S, QT
    real(kind=real64), dimension(10) :: F, Q, R
    
    integer(kind = int32), intent(in) :: n1, n2
    real(kind = real64), intent(out) :: gam

    ! ELIMINAR:
    integer(kind = int32) :: NK
    real(kind=real64), dimension(10,10) :: P, TAU
    real(kind=real64) :: T
    COMMON/CUFAC/NK,NG,P,T                                     
    COMMON/CPAR/TAU,S,F                             
    COMMON/CQT/QT,Q,R                                  
    !common/grupas1/rkass(6,12,6,12),enass(6,12,6,12), deloh(6,12,6,12) !Alfonsina

    !deloh =     0.D0
    ! ----------------------

                                                              
    GAM = & 
        & F(n2) + &
        & Q(n2) * (1.D0 - DLOG(Q(n2) / Q(n1))) - &
        & (R(n2) / R(n1)) + DLOG((R(n2) / R(n1))) - &
        & 5.D0 * Q(n2) * (1.D0 - ((R(n2) / R(n1)) / (Q(n2) / Q(n1))) + &
        & DLOG((R(n2) / R(n1))) - DLOG(Q(n2) / Q(n1)))

    do i = 1, NG                                                      
        GAM = GAM - S(i, n2) / S(i, n1) * QT(i, n1) - QT(i,n2) * DLOG(S(i, n1))           
    end do

    return
end subroutine gaminf                                                        
