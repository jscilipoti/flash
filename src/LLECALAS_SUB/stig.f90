subroutine STIG(y, S)                                              
    use iso_fortran_env, only: real64, int32
    
    !IMPLICIT REAL * 8(A-H, O-Z)   
    implicit none                                       
                               
    COMMON/CUFAC/n, NG, P, T                                      
    COMMON/CACT/Y1, Y2, ACT1, ACT2, DACT1, DACT2, pact                                                     
    COMMON/CGIBBS/nf, MAXZ, GNUL, Z, A, XVL, SFAS, GAM, aa, da, XM                                       
    
    !real(kind = real64), dimension(10) :: y, V, YGEM
                                        
    !USED:
    integer(kind = int32) :: i, j, jm, k, kk, na, neg, ngn
    integer(kind = int32) :: jpur = 0
    real(kind = real64) :: GD, GG, D, CC, amax
    real(kind = real64) :: & 
        !& amax = 1.D0, &
        !& CC = 1.D0, &
        !& D = 1.D0, &
        !& GD = 1.D0, &
        !& GG = 1.D0, &
        & rmax = 1.D0, &
        & RT1 = 1.D0, &
        & SUM = 1.D0, &
        & VV = 1.D0, &
        & YV1 = 1.D0, &
        & YV2 = 1.D0
    real(kind = real64), dimension(10) :: &
        &   V = 0.D0, &
        &   YGEM = 0.D0

    real(kind = real64), dimension(10) :: Y
    real(kind = real64) :: S
    
    !COMMON
    integer(kind = int32) :: n, nf
    real(kind = real64), dimension(4) :: SFAS
    real(kind = real64), dimension(10) :: A, AA, Z
    real(kind = real64), dimension(2, 2) :: pact
    real(kind = real64), dimension(10, 4) :: XVL
    real(kind = real64), dimension(10, 10) :: GAM, da

    
    
    !NOT USED
    integer(kind = int32) :: ng, maxz, gnul
    real(kind = real64) :: T
    real(kind = real64), dimension(10) :: Y1, Y2, ACT1, ACT2
    real(kind = real64), dimension(10, 4) :: XM
    real(kind = real64), dimension(10, 10) :: P, DACT1, DACT2

    do i = 1, n                                                       
        if(A(i) > amax) then                                         
            jpur = i                                                            
            amax = A(i)
        end if
    end do                                                         
                                                          
    rmax = 1.0D5                                                         
    ngn = n                                                             
    if(nf > 1) then
        ngn = n + nf
    end if                                              
    neg = 0

    mainloop : do kk = 1, ngn                                                   
        jm = kk                                                             
        if(jpur /= 0) jm = jpur                                             
        
        if(.not.(jm <= n)) then                                                
        
            do i = 1, n                                                       
                y(i) = Z(i) * (2.D0 + XVL(i, jm-n) / SFAS(jm - n)) / 3.D0
            end do                            
            
        else                                                           
    
            SUM = 0.D0                                                            
        
            do i = 1, n                                                       
                GG = A(i)-GAM(jm, i)                                                 
                if(GG < -50.D0) GG = -50.D0                                        
                y(i) = DEXP(GG)                                                     
                SUM = SUM + y(i)
            end do    

        end if
        
        na = 3                                                              
        
        inner: do k = 1, na                                                      
            do i = 1, n                                                       
                y(i) = y(i) / SUM
            end do

            call unifac(1, y, aa, da, pact)

            if(k == na) exit inner

            do i = 1, n
                y(i) = DEXP(A(i)-aa(i))
            end do                                             
            
            SUM = 0.D0                                                            
            
            do i = 1, n                                                       
                SUM = SUM + y(i)
            end do
        end do inner
                                                         
        
        YV1 = 0.D0                                                            
        
        do j = 1, nf                                                      
            V(j) = 0.D0
        end do

        do i = 1, n                                                       
            GD = DLOG(y(i)) + aa(i)-A(i)                                          
            YV1 = YV1 + y(i) * GD                                                   
            do j = 1, nf                                                      
                k = j                                                               
                VV = XVL(i, k) * Z(i)/SFAS(k)                                          
                D = GD * (y(i)-VV)                                                    
                V(j) = V(j) + D
            end do
        end do                                                       
        
        YV2 = V(1)                                                          
        
        do j = 1, nf                                                      
            if(V(j) < YV2) then 
                YV2 = V(j)
            end if                                          
        end do

        RT1 = YV1

        if(YV2 > 0.) then 
            RT1 = RT1 - YV2 / 2
        end if

        if(.not.(neg == 0 .and. YV1 > 0.)) then                                
            RT1 = YV1                                                           
            if(neg == 0) rmax = 0.D0                                              
            neg = 1
        end if

        if(RT1 > rmax) cycle                                          
        S = YV1                                                             
        rmax = RT1                                                          
        CC = DEXP(-YV1)                                                     
        
        do i = 1, n                                                       
            YGEM(i) = y(i) * CC
        end do                                                   
        
        if(jpur /= 0) then 
            exit mainloop 
        end if

    end do mainloop                                                          


    do i = 1, n                                                      
        y(i) = YGEM(i)
    end do                                                      
    
    return

end subroutine stig