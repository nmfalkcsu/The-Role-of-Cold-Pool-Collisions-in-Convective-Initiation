module base_state_nocturnal

    implicit none

    contains

    subroutine base_state_nocturnal()

        use constants
        use constants_b
        use model_arrays
        !use sounding
        implicit none

        integer k

        real log_u

        real, dimension(NZ) :: &
            THETAB_z, &
            THETAB_w_z, &
            RVB_z, &
            RVB_w_z

        real N_squared, theta_0, P_sfc

        N_squared = 0.02**2
        theta_0 = 300
        P_sfc = 100000.

        !make z grids
        do k = 1, NZ
            w_grid_z(k) = delta_z*(k-2)
            s_grid_z(k) = delta_z*(k-1.5)
        enddo

        !make x grids
        do k = 1, NX
            u_grid_x(k) = delta_x*(k-2)
            s_grid_x(k) = delta_x*(k-1.5)
        enddo

        !create base state theta
        do k = 2, NZ
            if (s_grid_z(k)<=1000) then
                THETAB(:,k) = theta_0*exp(N_squared*s_grid_z(k)/g)
            else 
                THETAB(:,k) = theta_0*exp(1000*N_squared/g)       
            endif
        enddo
        THETAB(:,1) = THETAB(:,2)

        !create base state r
        RVB(:,:) = 0
        RVB(:,1) = RVB(:,2)
        
        !create base state thetav
                                                    !RVB(:,:) = 0.00
        THETAVB = THETAB*(1+0.61*RVB)
        
                                                !THETAB(:,:) = 300
                                                !THETAVB(:,:) = 300

        !create base state exener function
        PIB_w(:,2) = (P_sfc/P_0)**(Rd/c_p)
        PIB(:,1) = PIB_w(:,2) + g/(c_p*(THETAVB(:,1)))*(w_grid_z(2)-s_grid_z(1))
        do k = 2, (NZ-1)
            PIB(:,k) = PIB(:,k-1) - g/(c_p*(THETAVB(:,k)+THETAVB(:,k-1))/2)*(w_grid_z(k+1)-w_grid_z(k))
        enddo
        PIB(:,NZ) = PIB(:,k-1) - g/(c_p*(THETAVB(:,k)+THETAVB(:,k-1))/2)*(delta_z)
        !create base state temperature
        TB = THETAB*PIB

        !create base state pressure
        PB = P_0*(PIB)**(c_p/Rd)
        
        !create base state density at scalar levels
        RHOU = P_0*(PIB)**(c_v/Rd)/(Rd*THETAVB)

        !create base state relative humitidy
        RVSB = (380/PB)*exp((17.27*(TB-273.0))/(TB-36.0))
        RHB = 100*RVB/RVSB
        


        !create base state density at w levels
        do k = 2, NZ
            if (w_grid_z(k)<=1000) then
                THETAB_w(:,k) = theta_0*exp(N_squared*w_grid_z(k)/g)
            else 
                THETAB_w(:,k) = theta_0*exp(1000*N_squared/g)       
            endif
        enddo
        THETAB_w(:,1) = THETAB_w(:,2)

        !create base state r
        RVB_w(:,:) = 0
        RVB_w(:,1) = RVB_w(:,2)
                                                                            !RVB_w(:,:) = 0
        THETAVB_w = THETAB_w*(1+0.61*RVB_w)
                                                                            !THETAVB_w(:,:) = 300
        PIB_w(:,1) = PIB_W(:,2) + g/(c_p*(THETAVB_w(:,1)))*(w_grid_z(2)-w_grid_z(1))
        do k = 3, NZ
            PIB_w(:,k) = PIB_w(:,k-1) - g/(c_p*(THETAVB_w(:,k)+THETAVB_w(:,k-1))/2)*(s_grid_z(k)-s_grid_z(k-1))
        enddo
        RHOW = P_0*(PIB_w)**(c_v/Rd)/(Rd*THETAVB_w)
    

        UB(:,:) = 0

        !UB(:,:) = 10

        
        !do k=2,NZ-1
        !    if ((k-1)<10) then
        !        UB(:,k) = k-1
        !    else
        !        UB(:,k) = 10
        !    endif
        !enddo
        !UB(:,1) = UB(:,2)
        !UB(:,NZ) = UB(:,NZ-1)

        !do k = 2, NZ
        !    log_u = (u_star/vk_constant)*log((s_grid_z(k)-zero_plane_displacement)/roughness_length)
        !    if(log_u<0) then
        !        UB(:,k) = 0
        !    else if(log_u>max_u) then
        !        UB(:,k) = max_u
        !    else
        !        UB(:,k) = log_u
        !    endif
        !enddo


    end subroutine base_state_nocturnal
    
end module base_state_nocturnal