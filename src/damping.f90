module damping

    implicit none

    contains

    subroutine damping()

        use constants
        use constants_b
        use model_arrays

        real :: coef, coef2

        real, dimension(NX,NZ) :: th_p_tend, u_tend, w_tend

        integer :: i, j, k, zstart


        zstart = nz-damp_levels
        coef = 1. / ( damp_tau*(s_grid_z(NZ)-s_grid_z(zstart)))

        th_p_tend = (th_p_f-th_p_p)/(2*dt)
        u_tend = (u_f-u_p)/(2*dt)
        w_tend = (w_f-w_p)/(2*dt)

        do k=zstart+1,nz
            do i=1,nx
                coef2 = coef * ( s_grid_z(k)-s_grid_z(zstart) )
                th_p_tend(i,k) = th_p_tend(i,k) - coef2*th_p_c(i,k)
                u_tend(i,k) = u_tend(i,k) - coef2*u_c(i,k)
                w_tend(i,k) = w_tend(i,k) - coef2*w_c(i,k)
            enddo
        enddo

        th_p_f = th_p_p + 2*dt*th_p_tend
        u_f = u_p + 2*dt*u_tend    
        w_f = w_p + 2*dt*w_tend

    end subroutine damping
    
end module damping