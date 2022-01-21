module store_saved_timeseries

    implicit none

    contains

    subroutine store_saved_timeseries()

        use model_arrays
        use constants_b
        use saved_timeseries
        
        w_max_ts(l+1) = maxval(w_c)
        w_max_loc_ts(:,l+1) = maxloc(w_c)
        th_at_w_max_loc_ts(l+1) = th_p_c(w_max_loc_ts(1,l+1),w_max_loc_ts(2,l+1))
        th_min_ts(l+1) = minval(th_p_c)
        u_abs_max_ts(l+1) = maxval(abs(u_c))

    end subroutine store_saved_timeseries

end module store_saved_timeseries