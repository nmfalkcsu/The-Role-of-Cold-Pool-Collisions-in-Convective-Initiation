module model_integration

    implicit none

    contains

    subroutine model_integration()

        use constants
        use constants_b
        use model_arrays
        use save_vars_no_mphys
        use solve_prog
        use get_k
        use store_saved_timeseries
        use save_saved_timeseries
        use save_base_state
        use damping

        implicit none

        rd2x = 1/(2*delta_x)
        rd2z = 1/(2*delta_z)
        rdx = 1/(delta_x)
        rdz = 1/(delta_z)
    
        l = 0

        call store_saved_timeseries
        call save_vars_no_mphys
        call save_base_state
        call get_k


        do l = 1, num_timesteps
            call solve_prog

            if (l>1) call damping

            if(IARF == 1) then
                if(l>1) then
                    u_c = u_c+ARF_gamma*(u_f-2*u_c+u_p)
                    w_c = w_c+ARF_gamma*(w_f-2*w_c+w_p)
                    th_p_c = th_p_c+ARF_gamma*(th_p_f-2*th_p_c+th_p_p)
                    pi_p_c = pi_p_c+ARF_gamma*(pi_p_f-2*pi_p_c+pi_p_p)
                endif
            endif

            u_p = u_c 
            u_c = u_f
            pi_p_p = pi_p_c 
            pi_p_c = pi_p_f
            w_p = w_c 
            w_c = w_f
            th_p_p = th_p_c 
            th_p_c = th_p_f
            rv_p_p = rv_p_c
            rv_p_c = rv_p_f
            rc_p = rc_c
            rc_c = rc_f

            if(mod((l*dt),save_freq)==0) then
                call save_vars_no_mphys
            endif

            call store_saved_timeseries
        enddo

        call save_saved_timeseries

    end subroutine model_integration
    
end module model_integration