module base_state_dry_isentropic

    implicit none

    contains

    subroutine base_state_dry_isentropic()

        use constants
        use constants_b
        use model_arrays
        implicit none

        integer k

        real P_sfc

        !Code to make dry isentropic (theta = 300 K) base state

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

        !define surface pressure
        P_sfc = 100000.

        !create base state theta
        THETAB(:,:) = 300
        THETAB(:,1) = THETAB(:,2)
        THETAB_w(:,:) = 300
        THETAB_w(:,1) = THETAB_w(:,2)

        !create base state rv
        RVB(:,:) = 0
        RVB(:,1) = RVB(:,2)
        RVB_w(:,:) = 0
        rvb_w(:,1) = rvb_w(:,2)
        
        !create base state thetav
        THETAVB = THETAB*(1+0.61*RVB)

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
        THETAVB_w = THETAB_w*(1+0.61*RVB_w)

        PIB_w(:,1) = PIB_W(:,2) + g/(c_p*(THETAVB_w(:,1)))*(w_grid_z(2)-w_grid_z(1))
        do k = 3, NZ
            PIB_w(:,k) = PIB_w(:,k-1) - g/(c_p*(THETAVB_w(:,k)+THETAVB_w(:,k-1))/2)*(s_grid_z(k)-s_grid_z(k-1))
        enddo
        RHOW = P_0*(PIB_w)**(c_v/Rd)/(Rd*THETAVB_w)


    end subroutine base_state_dry_isentropic
    
end module base_state_dry_isentropic