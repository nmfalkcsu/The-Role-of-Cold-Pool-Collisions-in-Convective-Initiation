module save_saved_timeseries

    implicit none

    contains

    subroutine save_saved_timeseries()

        use looping_stuff
        use saved_timeseries
        use constants_b

        integer k

        character(len=1024) :: filename
        character(len=7) :: specific

        write (specific, "(I2.2,a1,I3.3,a1)") theta_p_loop,'_',dfc_loop,'_'
        !print*,specific
        write (filename, "(a33,a7,a5,a4)") file_path,specific,'final','.dat'


        open(1,file=filename)
        write(1,*), 'w_max_ts,w_max_x_ts,w_max_y_ts,th_at_w_max_loc_ts,th_min_ts,u_abs_max_ts'
        do k = 1, num_timesteps+1
                write(1,'(ES15.8E3,a1,I5,a1,I5,a1,ES15.8E3,a1,ES15.8E3,a1,ES15.8E3)'), &
                w_max_ts(k),',',w_max_loc_ts(1,k),',',w_max_loc_ts(2,k),',',th_at_w_max_loc_ts(k),',',th_min_ts(k), &
                ',',u_abs_max_ts(k)
        enddo
        close(1)



    end subroutine save_saved_timeseries
end module save_saved_timeseries