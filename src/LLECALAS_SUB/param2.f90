subroutine param2
    use iso_fortran_env, only: int32, real64                                                
    !IMPLICIT REAL * 8(A-H, O-Z)
    implicit none

    COMMON/CUFAC/nk, ng, P, T                                     
    COMMON/CPAR/TAU, S, F                            
    COMMON/CQT/QT, Q, R

    integer(kind = int32) :: i, j, k
    integer(kind = int32) :: ng, nk
    real(kind = real64), dimension(10,10) :: P, QT, TAU, S
    real(kind = real64), dimension(10) :: F
    real(kind = real64) :: T

    ! NOT USED
    real(kind = real64), dimension(10) :: Q, R
    ! ----

    do i = 1, ng                                                      
        do j = 1, ng                                                      
            TAU(i, j) = DEXP(-P(i, j) / T)
        end do
    end do

    do i = 1, nk                                                      
        do j = 1, ng                                                      
            S(j, i) = 0.D0                                                       
            do k = 1, ng                                                      
                S(j, i) = S(j, i) + QT(k, i) * TAU(k, j)
            end do
        end do
    end do                                    
    do i = 1, nk                                                      
        F(i) = 1.D0                                                         
        do j = 1, ng                                                      
            F(i) = F(i) + QT(j, i) * DLOG(S(j, i))
        end do
    end do
    
    return
end subroutine param2

