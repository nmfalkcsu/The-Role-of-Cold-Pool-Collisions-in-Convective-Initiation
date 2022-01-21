module save_base_state

    implicit none

    contains

    subroutine save_base_state()

        use constants
        use constants_b
        use model_arrays
        use looping_stuff

        integer i,k

        character(len=1024) :: filename
        character(len=7) :: specific

        write (specific, "(I2.2,a1,I3.3,a1)") theta_p_loop,'_',dfc_loop,'_'
        !print*,specific
        write (filename, "(a33,a7,a4,a4)") file_path,specific,'base','.dat'


        open(1,file=filename)
        write(1,'(I4,a1,I4)') NX,',',NZ
        write(1,*), 'rvb,thetab,thetab_w,rvb_w,pib_w,pib,pb'
        do k = 1, NZ
            do i = 1, NX
                write(1,'(ES15.8E3,a1,ES15.8E3,a1,ES15.8E3,a1,ES15.8E3,a1,ES15.8E3,a1,ES15.8E3,a1,ES15.8E3)'), &
                rvb(i,k),',',thetab(i,k),',',thetab_w(i,k),',',rvb_w(i,k),',',&
                pib_w(i,k),',',pib(i,k),',',PB(i,k)
            enddo
        enddo
        close(1)


    end subroutine save_base_state
    
end module save_base_state