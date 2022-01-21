module save_vars_no_mphys

    implicit none

    contains

    subroutine save_vars_no_mphys()

        use constants
        use constants_b
        use model_arrays
        use looping_stuff

        integer i,k,miliseconds

        character(len=1024) :: filename
        character(len=7) :: specific

        write (6,*) theta_p_loop*-1,'K',(dfc_loop-24)*200,'m',nint(l*dt/60),'min'
        call flush(6)

        miliseconds = nint(l*dt*1000)
        write (specific, "(I2.2,a1,I3.3,a1)") theta_p_loop,'_',dfc_loop,'_'
        !print*,specific
        write (filename, "(a33,a7,I13.13,a4)") file_path,specific,miliseconds,'.dat'


        open(1,file=filename)
        write(1,'(I4,a1,I4)') NX,',',NZ
        write(1,*), 'th_p_c,u_c,w_c,pi_p_c,rv_p_c,rc_c'
        do k = 1, NZ
            do i = 1, NX
                write(1,'(ES15.8E3,a1,ES15.8E3,a1,ES15.8E3,a1,ES15.8E3,a1,ES15.8E3,a1,ES15.8E3)'), &
                Th_p_c(i,k),',',u_c(i,k), &
                ',',w_c(i,k),',',pi_p_c(i,k),',',rv_p_c(i,k),',',rc_c(i,k)
            enddo
        enddo
        close(1)


    end subroutine save_vars_no_mphys
    
end module save_vars_no_mphys