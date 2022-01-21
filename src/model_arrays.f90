module model_arrays
      use constants_b
      use constants
      
   implicit none



      real, dimension(NX,NZ) :: &
            THETAB, &
            PIB, &
            RVB, &
            THETAVB, &
            RHOU, &
            TB, &
            PB, &
            RVSB, &
            RHB, &

            PIB_w, &
            THETAB_w, &
            THETAVB_w, &
            RVB_w, &
            RHOW, &

            TH_p_p, &
            TH_p_c, &
            TH_p_f, &
            PI_p_p, &
            PI_p_c, &
            PI_p_f, &
            RV_p_p, &
            RV_p_c, &
            RV_p_f, &
            u_p, &
            u_c, &
            u_f, &
            w_p, &
            w_c, &
            w_f, &
            P_p, &

            Kmx, &
            Kmz, &
            Khx, &
            Khz, &
            rho, &

            rv_b, &
            rc_p, &
            rc_c, &
            rc_f, &

            Vt, &
            rr_p, &
            rr_c, &
            rr_f, &
            pre_lambda, &
            ub, &
            lambda

      real :: &
            n0, gamma45, kvt



       real, dimension(NZ) :: &
            w_grid_z, &
            s_grid_z
      
      real, dimension(NX) :: &
            u_grid_x, &
            s_grid_x


end module model_arrays