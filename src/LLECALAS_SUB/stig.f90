subroutine STIG(Y, S)                                              
    use iso_fortran_env, only: real64, int32
    
    implicit none                                       
                               
    COMMON/CUFAC/n, NG, P, T                                      
    COMMON/CACT/Y1, Y2, ACT1, ACT2, DACT1, DACT2, pact                                                     
    COMMON/CGIBBS/nf, MAXZ, GNUL, Z, A, XVL, SFAS, GAM, aa, da, XM                                       
                                           
    !USED:
    integer(kind = int32) :: i, j, k, l
    integer(kind = int32) :: &
        & j_max = 0, &
        & A_max_index = 0, &
        & k_max = 3, &
        & neg = 0, &
        & ngn = 0

    real(kind = real64) :: & 
        A_max = 0.D0, &
        & CC = 0.D0, &
        & D = 0.D0, &
        & GD = 0.D0, &
        & GG = 0.D0, &
        & GG_minval = -50.D0, &
        & rmax = 1.0D5, &
        & RT1 = 0.D0, &
        & sum = 0.D0, &
        & VV = 0.D0, &
        & YV1 = 0.D0, &
        & YV2 = 0.D0

    real(kind = real64), dimension(10) :: &
        &   V = 0.D0, &
        &   YGEM = 0.D0

    real(kind = real64), dimension(10), intent(out) :: Y
    real(kind = real64), intent(out) :: S
    
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
    
    

    ! Search for the max value of A(i) and its index
    do i = 1, n                                                       
        if(A(i) > A_max) then                                         
            A_max_index = i                                                            
            A_max = A(i)
        end if
    end do                                                         
                                                          
    Y = 0.D0
    S = 0.D0
    rmax = 1.0D5
    neg = 0                                                         
    ngn = n

    if(nf > 1) then
        ngn = n + nf
    end if                                              
    

    mainloop : do l = 1, ngn                                                   
        j_max = l                                                             
        if(A_max_index /= 0) j_max = A_max_index                                             
        
        if(.not.(j_max <= n)) then                                                
        
            do i = 1, n                                                       
                Y(i) = Z(i) * (2.D0 + XVL(i, j_max-n) / SFAS(j_max - n)) / 3.D0
            end do                            
            
        else                                                           
            ! Sum of all Y(i)
            sum = 0.D0                                                            
            do i = 1, n                                                       
                GG = A(i) - GAM(j_max, i)                                                 
                if(GG < GG_minval) GG = GG_minval                                       
                Y(i) = DEXP(GG)                                                     
                sum = sum + Y(i)
            end do    

        end if
        
        k_max = 3                                                              
        
        inner: do k = 1, k_max                                                      
            ! Get the fractions of each Y(i)
            do i = 1, n                                                       
                Y(i) = Y(i) / sum
            end do

            call unifac(1, Y, aa, da, pact)

            if(k == k_max) exit inner

            do i = 1, n
                Y(i) = DEXP(A(i)-aa(i))
            end do                                             
            
            ! Sum of all Y(i)
            sum = 0.D0                                                            
            do i = 1, n                                                       
                sum = sum + Y(i)
            end do
        end do inner
                                                         
        
        YV1 = 0.D0                                                            
        
        do j = 1, nf                                                      
            V(j) = 0.D0
        end do

        do i = 1, n                                                       
            GD = DLOG(Y(i)) + aa(i)-A(i)                                          
            YV1 = YV1 + Y(i) * GD                                                   
            do j = 1, nf                                                      
                k = j                                                               
                VV = XVL(i, k) * Z(i)/SFAS(k)                                          
                D = GD * (Y(i)-VV)                                                    
                V(j) = V(j) + D
            end do
        end do                                                       
        
        YV2 = V(1)                                                          
        
        ! Get the min value from all V(j)
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
            YGEM(i) = Y(i) * CC
        end do                                                   
        
        if(A_max_index /= 0) then 
            exit mainloop 
        end if

    end do mainloop                                                          


    do i = 1, n                                                      
        Y(i) = YGEM(i)
    end do                                                      
    
    return

end subroutine stig