module solve_prog

    implicit none

    contains

    subroutine solve_prog()

        use constants
        use constants_b
        use model_arrays

        implicit none

        integer i,k
        
        if(l==1) then
            do i = 2,NX-1
                do k = 2,NZ-1
                    u_f(i,k) = u_c(i,k) + dt*(-1*rdx*(((u_c(i+1,k)+u_c(i,k))*.5)**2-((u_c(i,k)+u_c(i-1,k))*.5)**2.) &
                    - rdz/rhou(i,k)*((.25*rhow(i,k+1)*(w_c(i-1,k+1)+w_c(i,k+1))*(u_c(i,k)+u_c(i,k+1)))-(.25*rhow(i,k)*(w_c(i-1,k)+w_c(i,k))*(u_c(i,k-1)+u_c(i,k)))) &
                    - c_p*thetavb(i,k)*rdx*(pi_p_c(i,k)-pi_p_c(i-1,k)) &
                    +Kmx(i,k)*(u_c(i+1,k)-2*u_c(i,k)+u_c(i-1,k))*rdx**2 &
                    +Kmz(i,k)*(u_c(i,k+1)-2*u_c(i,k)+u_c(i,k-1))*rdz**2)
                enddo
            enddo
        else
            do i = 2,NX-1
                do k = 2,NZ-1
                    u_f(i,k) = u_p(i,k) + 2*dt*(-1*rdx*(((u_c(i+1,k)+u_c(i,k))*.5)**2-((u_c(i,k)+u_c(i-1,k))*.5)**2.) &
                    - rdz/rhou(i,k)*((.25*rhow(i,k+1)*(w_c(i-1,k+1)+w_c(i,k+1))*(u_c(i,k)+u_c(i,k+1)))-(.25*rhow(i,k)*(w_c(i-1,k)+w_c(i,k))*(u_c(i,k-1)+u_c(i,k)))) &
                    - c_p*thetavb(i,k)*rdx*(pi_p_c(i,k)-pi_p_c(i-1,k)) &
                    +Kmx(i,k)*(u_p(i+1,k)-2*u_p(i,k)+u_p(i-1,k))*rdx**2 &
                    +Kmz(i,k)*(u_p(i,k+1)-2*u_p(i,k)+u_p(i,k-1))*rdz**2)
                enddo
            enddo
        endif
        u_f(2:nx-1,1) = u_f(2:nx-1,2)
        u_f(2:nx-1,nz) = u_f(2:nx-1,nz-1)
        u_f(1,:) = u_f(nx-1,:)
        u_f(nx,:) = u_f(2,:)
        

        if(l==1) then
            do i = 2,NX-1
                do k = 2,NZ-1
                    pi_p_f(i,k) = pi_p_c(i,k) + dt*(-(c_s**2)/(rhou(i,k)*c_p*thetavb(i,k)**2)*(rhou(i,k)*thetavb(i,k)*rdx*(u_c(i+1,k)-u_c(i,k)) &
                    + rdz*.5*(w_c(i,k+1)*((thetavb(i,k+1)+thetavb(i,k))*(rhow(i,k+1)))-w_c(i,k)*(thetavb(i,k)+thetavb(i,k-1))*(rhow(i,k)))) &
                    +Khx(i,k)*(pi_p_c(i+1,k)-2*pi_p_c(i,k)+pi_p_c(i-1,k))*rdx**2 &
                    +Khz(i,k)*(pi_p_c(i,k+1)-2*pi_p_c(i,k)+pi_p_c(i,k-1))*rdz**2) 
                enddo
            enddo
        else
            do i = 2,NX-1
                do k = 2,NZ-1
                    pi_p_f(i,k) = pi_p_p(i,k) + 2*dt*(-(c_s**2)/(rhou(i,k)*c_p*thetavb(i,k)**2)*(rhou(i,k)*thetavb(i,k)*rdx*(u_c(i+1,k)-u_c(i,k)) &
                    + rdz*.5*(w_c(i,k+1)*(thetavb(i,k+1)+thetavb(i,k))*(rhow(i,k+1))-w_c(i,k)*(thetavb(i,k)+thetavb(i,k-1))*(rhow(i,k)))) &
                    +Khx(i,k)*(pi_p_p(i+1,k)-2*pi_p_p(i,k)+pi_p_p(i-1,k))*rdx**2 &
                    +Khz(i,k)*(pi_p_p(i,k+1)-2*pi_p_p(i,k)+pi_p_p(i,k-1))*rdz**2)
                enddo
            enddo
        endif
        pi_p_f(2:nx-1,1) = pi_p_f(2:nx-1,2)
        pi_p_f(2:nx-1,nz) = pi_p_f(2:nx-1,nz-1)
        pi_p_f(1,:) = pi_p_f(nx-1,:)
        pi_p_f(nx,:) = pi_p_f(2,:)


        if(l==1) then
            do i = 2,NX-1
                do k = 3,NZ-1
                    w_f(i,k) = w_c(i,k) +  dt*(-1*rdx*.25*((u_c(i+1,k)+u_c(i+1,k-1))*(w_c(i+1,k)+w_c(i,k)) - (u_c(i,k)+u_c(i,k-1))*(w_c(i,k)+w_c(i-1,k))) &
                    -(1/rhow(i,k))*rdz*((rhou(i,k)*((w_c(i,k+1)+w_c(i,k))*.5)**2)-(rhou(i,k-1)*((w_c(i,k)+w_c(i,k-1))*.5)**2)) &
                    -c_p*.5*(thetavb(i,k)+thetavb(i,k-1))*rdz*(pi_p_c(i,k)-pi_p_c(i,k-1)) &
                    +g*((th_p_c(i,k)+th_p_c(i,k-1))/(thetab(i,k)+thetab(i,k-1)) + 0.61*.5*(rv_p_c(i,k)+rv_p_c(i,k-1)))  &
                    +Kmx(i,k)*(w_c(i+1,k)-2*w_c(i,k)+w_c(i-1,k))*rdx**2 &
                    +Kmz(i,k)*(w_c(i,k+1)-2*w_c(i,k)+w_c(i,k-1))*rdz**2)
                enddo
            enddo
        else
            do i = 2,NX-1
                do k = 3,NZ-1
                    w_f(i,k) = w_p(i,k) + 2*dt*(-1*rdx*.25*((u_c(i+1,k)+u_c(i+1,k-1))*(w_c(i+1,k)+w_c(i,k)) - (u_c(i,k)+u_c(i,k-1))*(w_c(i,k)+w_c(i-1,k))) &
                    -(1/rhow(i,k))*rdz*((rhou(i,k)*((w_c(i,k+1)+w_c(i,k))*.5)**2)-(rhou(i,k-1)*((w_c(i,k)+w_c(i,k-1))*.5)**2)) &
                    -c_p*.5*(thetavb(i,k)+thetavb(i,k-1))*rdz*(pi_p_c(i,k)-pi_p_c(i,k-1)) &
                    +g*((th_p_c(i,k)+th_p_c(i,k-1))/(thetab(i,k)+thetab(i,k-1)) + 0.61*.5*(rv_p_c(i,k)+rv_p_c(i,k-1)))  &
                    +Kmx(i,k)*(w_p(i+1,k)-2*w_p(i,k)+w_p(i-1,k))*rdx**2 &
                    +Kmz(i,k)*(w_p(i,k+1)-2*w_p(i,k)+w_p(i,k-1))*rdz**2)
                enddo
            enddo
        endif    

        w_f(:,1) = 0
        w_f(:,2) = 0
        w_f(:,nz) = 0
        w_f(1,:) = w_f(nx-1,:)
        w_f(nx,:) = w_f(2,:)
        
        
        if(l==1) then
            do i = 2,NX-1
                do k = 2,NZ-1
                    th_p_f(i,k) = th_p_c(i,k) + dt*(-1*rdx*.5*(u_c(i+1,k)*(th_p_c(i+1,k)+th_p_c(i,k))-u_c(i,k)*(th_p_c(i,k)+th_p_c(i-1,k))) & 
                    -(1/rhou(i,k))*rdz*.5*((rhow(i,k+1)*w_c(i,k+1)*(th_p_c(i,k+1)+th_p_c(i,k)))-(rhow(i,k)*w_c(i,k)*(th_p_c(i,k)+th_p_c(i,k-1)))) &
                    -.5*(rdz/rhou(i,k))*(rhow(i,k)*w_c(i,k)*(thetab(i,k)-thetab(i,k-1))+rhow(i,k+1)*w_c(i,k+1)*(thetab(i,k+1)-thetab(i,k))) &
                    +Khx(i,k)*(th_p_c(i+1,k)-2*th_p_c(i,k)+th_p_c(i-1,k))*rdx**2 &
                    +Khz(i,k)*(th_p_c(i,k+1)-2*th_p_c(i,k)+th_p_c(i,k-1))*rdz**2) 
                enddo
            enddo
        else
            do i = 2,NX-1
                do k = 2,NZ-1
                    th_p_f(i,k) = th_p_p(i,k) + 2*dt*(-1*rdx*.5*(u_c(i+1,k)*(th_p_c(i+1,k)+th_p_c(i,k))-u_c(i,k)*(th_p_c(i,k)+th_p_c(i-1,k))) & 
                    -(1/rhou(i,k))*rdz*.5*((rhow(i,k+1)*w_c(i,k+1)*(th_p_c(i,k+1)+th_p_c(i,k)))-(rhow(i,k)*w_c(i,k)*(th_p_c(i,k)+th_p_c(i,k-1)))) &
                    -.5*(rdz/rhou(i,k))*(rhow(i,k)*w_c(i,k)*(thetab(i,k)-thetab(i,k-1))+rhow(i,k+1)*w_c(i,k+1)*(thetab(i,k+1)-thetab(i,k))) &
                    +Khx(i,k)*(th_p_p(i+1,k)-2*th_p_p(i,k)+th_p_p(i-1,k))*rdx**2 &
                    +Khz(i,k)*(th_p_p(i,k+1)-2*th_p_p(i,k)+th_p_p(i,k-1))*rdz**2) 
                enddo
            enddo
        endif
        th_p_f(2:nx-1,1) = th_p_f(2:nx-1,2)
        th_p_f(2:nx-1,nz) = th_p_f(2:nx-1,nz-1)
        th_p_f(1,:) = th_p_f(nx-1,:)
        th_p_f(nx,:) = th_p_f(2,:)
        
        if(l==1) then
            do i = 2,NX-1
                do k = 2,NZ-1
                    rv_p_f(i,k) = rv_p_c(i,k) + dt*(-1*rdx*.5*(u_c(i+1,k)*(rv_p_c(i+1,k)+rv_p_c(i,k))-u_c(i,k)*(rv_p_c(i,k)+rv_p_c(i-1,k))) & 
                    -(1/rhou(i,k))*rdz*.5*((rhow(i,k+1)*w_c(i,k+1)*(rv_p_c(i,k+1)+rv_p_c(i,k)))-(rhow(i,k)*w_c(i,k)*(rv_p_c(i,k)+rv_p_c(i,k-1)))) &
                    -.5*(rdz/rhou(i,k))*(rhow(i,k)*w_c(i,k)*(rvb(i,k)-rvb(i,k-1))+rhow(i,k+1)*w_c(i,k+1)*(rvb(i,k+1)-rvb(i,k))) &
                    +Khx(i,k)*(rv_p_c(i+1,k)-2*rv_p_c(i,k)+rv_p_c(i-1,k))*rdx**2 &
                    +Khz(i,k)*(rv_p_c(i,k+1)-2*rv_p_c(i,k)+rv_p_c(i,k-1))*rdz**2) 
                enddo
            enddo
        else
            do i = 2,NX-1
                do k = 2,NZ-1
                    rv_p_f(i,k) = rv_p_p(i,k) + 2*dt*(-1*rdx*.5*(u_c(i+1,k)*(rv_p_c(i+1,k)+rv_p_c(i,k))-u_c(i,k)*(rv_p_c(i,k)+rv_p_c(i-1,k))) & 
                    -(1/rhou(i,k))*rdz*.5*((rhow(i,k+1)*w_c(i,k+1)*(rv_p_c(i,k+1)+rv_p_c(i,k)))-(rhow(i,k)*w_c(i,k)*(rv_p_c(i,k)+rv_p_c(i,k-1)))) &
                    -.5*(rdz/rhou(i,k))*(rhow(i,k)*w_c(i,k)*(rvb(i,k)-rvb(i,k-1))+rhow(i,k+1)*w_c(i,k+1)*(rvb(i,k+1)-rvb(i,k))) &
                    +Khx(i,k)*(rv_p_p(i+1,k)-2*rv_p_p(i,k)+rv_p_p(i-1,k))*rdx**2 &
                    +Khz(i,k)*(rv_p_p(i,k+1)-2*rv_p_p(i,k)+rv_p_p(i,k-1))*rdz**2) 
                enddo
            enddo
        endif
        rv_p_f(2:nx-1,1) = rv_p_f(2:nx-1,2)
        rv_p_f(2:nx-1,nz) = rv_p_f(2:nx-1,nz-1)
        rv_p_f(1,:) = rv_p_f(nx-1,:)
        rv_p_f(nx,:) = rv_p_f(2,:)

        if(l==1) then
            do i = 2,NX-1
                do k = 2,NZ-1
                    rc_f(i,k) = rc_c(i,k) + dt*(-1*rdx*.5*(u_c(i+1,k)*(rc_c(i+1,k)+rc_c(i,k))-u_c(i,k)*(rc_c(i,k)+rc_c(i-1,k))) & 
                    -(1/rhou(i,k))*rdz*.5*((rhow(i,k+1)*w_c(i,k+1)*(rc_c(i,k+1)+rc_c(i,k)))-(rhow(i,k)*w_c(i,k)*(rc_c(i,k)+rc_c(i,k-1)))) &
                    +Khx(i,k)*(rc_c(i+1,k)-2*rc_c(i,k)+rc_c(i-1,k))*rdx**2 &
                    +Khz(i,k)*(rc_c(i,k+1)-2*rc_c(i,k)+rc_c(i,k-1))*rdz**2) 
                enddo
            enddo
        else
            do i = 2,NX-1
                do k = 2,NZ-1
                    rc_f(i,k) = rc_p(i,k) + 2*dt*(-1*rdx*.5*(u_c(i+1,k)*(rc_c(i+1,k)+rc_c(i,k))-u_c(i,k)*(rc_c(i,k)+rc_c(i-1,k))) & 
                    -(1/rhou(i,k))*rdz*.5*((rhow(i,k+1)*w_c(i,k+1)*(rc_c(i,k+1)+rc_c(i,k)))-(rhow(i,k)*w_c(i,k)*(rc_c(i,k)+rc_c(i,k-1)))) &
                    +Khx(i,k)*(rc_p(i+1,k)-2*rc_p(i,k)+rc_p(i-1,k))*rdx**2 &
                    +Khz(i,k)*(rc_p(i,k+1)-2*rc_p(i,k)+rc_p(i,k-1))*rdz**2) 
                enddo
            enddo
        endif
        rc_f(2:nx-1,1) = rc_f(2:nx-1,2)
        rc_f(2:nx-1,nz) = rc_f(2:nx-1,nz-1)
        rc_f(1,:) = rc_f(nx-1,:)
        rc_f(nx,:) = rc_f(2,:)
        

    end subroutine solve_prog
    
end module solve_prog