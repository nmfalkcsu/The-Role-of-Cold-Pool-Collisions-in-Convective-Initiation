module saved_timeseries
      use constants_b
      
   implicit none

      real, dimension(num_timesteps+1) :: &
            w_max_ts, &
            th_at_w_max_loc_ts, &
            th_min_ts, &
            u_abs_max_ts

      integer, dimension(2,num_timesteps+1) :: &
            w_max_loc_ts
            
end module saved_timeseries