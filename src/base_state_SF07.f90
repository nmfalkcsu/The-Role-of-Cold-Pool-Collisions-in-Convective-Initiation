module base_state_SF07

    implicit none

    contains

    subroutine base_state_SF07()

        use constants
        use constants_b
        use model_arrays
        use sounding_SF07
        implicit none

        integer k

        real log_u, P_sfc

        !real, dimension(NZ) :: &
        !    THETAB_z, &
        !    THETAB_w_z, &
        !    RVB_z, &
        !    RVB_w_z

        !Code to make soudning from SF07 S07
        P_sfc = 1.0055*100*1000;
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
            thetab(:,k) = thetab_z(k)
        enddo
        THETAB(:,1) = THETAB(:,2)

        !create base state rv
        do k = 2, NZ
            RVB(:,k) = RVB_z(k)
        enddo
        RVB(:,1) = RVB(:,2)
        
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
        do k = 2, NZ
            thetab_w(:,k) = thetab_w_z(k)
        enddo
        THETAB_w(:,1) = THETAB_w(:,2)


        do k = 2, NZ
            RVB_w(:,k) = RVB_w_z(k)
        enddo
        RVB_w(:,1) = RVB_w(:,2)

        THETAVB_w = THETAB_w*(1+0.61*RVB_w)

        PIB_w(:,1) = PIB_W(:,2) + g/(c_p*(THETAVB_w(:,1)))*(w_grid_z(2)-w_grid_z(1))
        do k = 3, NZ
            PIB_w(:,k) = PIB_w(:,k-1) - g/(c_p*(THETAVB_w(:,k)+THETAVB_w(:,k-1))/2)*(s_grid_z(k)-s_grid_z(k-1))
        enddo
        RHOW = P_0*(PIB_w)**(c_v/Rd)/(Rd*THETAVB_w)
    

        UB(:,:) = 0



    end subroutine base_state_SF07
    
end module base_state_SF07