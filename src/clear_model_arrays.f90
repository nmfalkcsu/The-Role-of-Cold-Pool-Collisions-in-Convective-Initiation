module clear_model_arrays

    implicit none

    contains

    subroutine clear_model_arrays()

        use model_arrays

        implicit none

            THETAB(:,:) = 0
            PIB(:,:) = 0
            RVB(:,:) = 0
            THETAVB(:,:) = 0
            RHOU(:,:) = 0
            TB(:,:) = 0
            PB(:,:) = 0
            RVSB(:,:) = 0
            RHB(:,:) = 0

            PIB_w(:,:) = 0
            THETAB_w(:,:) = 0
            THETAVB_w(:,:) = 0
            RVB_w(:,:) = 0
            RHOW(:,:) = 0

            TH_p_p(:,:) = 0
            TH_p_c(:,:) = 0
            TH_p_f(:,:) = 0
            PI_p_p(:,:) = 0
            PI_p_c(:,:) = 0
            PI_p_f(:,:) = 0
            RV_p_p(:,:) = 0
            RV_p_c(:,:) = 0
            RV_p_f(:,:) = 0
            u_p(:,:) = 0
            u_c(:,:) = 0
            u_f(:,:) = 0
            w_p(:,:) = 0
            w_c(:,:) = 0
            w_f(:,:) = 0
            P_p(:,:) = 0

            Kmx(:,:) = 0
            Kmz(:,:) = 0
            Khx(:,:) = 0
            Khz(:,:) = 0
            rho(:,:) = 0

            rv_b(:,:) = 0
            rc_p(:,:) = 0
            rc_c(:,:) = 0
            rc_f(:,:) = 0


            w_grid_z(:) = 0
            s_grid_z(:) = 0
      
            u_grid_x(:) = 0
            s_grid_x(:) = 0



    end subroutine clear_model_arrays
    
end module clear_model_arrays