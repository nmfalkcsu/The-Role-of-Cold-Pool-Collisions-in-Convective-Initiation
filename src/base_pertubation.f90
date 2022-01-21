module base_pertubation

    implicit none

    contains

    subroutine base_pertubation()

        use constants
        use constants_b
        use model_arrays
        use looping_stuff

        implicit none

        integer i,k, dfc, AMP

        real imid, rad

        !set distance from center for bubbles and temperature perturbation
        dfc = dfc_loop
        AMP = -1*theta_p_loop

        if(dfc>=0) then
            !create first bubble
            imid = (NX/2.0 - dfc)+.5
            do i = 1, NX
                do k = 1, NZ
                    rad = sqrt( ((s_grid_z(k)-zcnt)/radz)**2 + ((delta_x*(i-imid))/radx)**2)
                    if (rad<1.0) then
                        TH_p_p(i,k) = 0.5*AMP*(cos(rad*trigpi)+1.0)
                        rc_p(i,k) = 1
                        rc_c(i,k) = 1
                    else
                        TH_p_p(i,k) = 0
                        rc_p(i,k) = 0
                        rc_c(i,k) = 0
                    endif
                enddo 
            enddo
            !create second bubble
            imid = (NX/2.0 + dfc)+.5
            do i = 1, NX
                do k = 1, NZ
                    rad = sqrt( ((s_grid_z(k)-zcnt)/radz)**2 + ((delta_x*(i-imid))/radx)**2)
                    if (rad<1.0) then
                        TH_p_p(i,k) = 0.5*AMP*(cos(rad*trigpi)+1.0)
                        rc_p(i,k) = 1
                        rc_c(i,k) = 1
                    endif
                enddo 
            enddo

        else
            !create only one bubble
            imid = (NX/2.0)+.5
            do i = 1, NX
                do k = 1, NZ
                    rad = sqrt( ((s_grid_z(k)-zcnt)/radz)**2 + ((delta_x*(i-imid))/radx)**2)
                    if (rad<1.0) then
                        TH_p_p(i,k) = 0.5*AMP*(cos(rad*trigpi)+1.0)
                        rc_p(i,k) = 1
                        rc_c(i,k) = 1
                    else
                        TH_p_p(i,k) = 0
                        rc_p(i,k) = 0
                        rc_c(i,k) = 0
                    endif
                enddo 
            enddo

        endif

    
        !make bottom zero gradient
        Th_p_p(:,1) = Th_p_p(:,2)

        !fix pi and pressure
        do i = 1, NX
            do k = NZ,1,-1
                if(k==NZ) then
                    pi_p_p(i,k) = 0
                else
                    pi_p_p(i,k) = pi_p_p(i,k+1) - g*TH_p_p(i,k)/(c_p*thetab(i,k)**2)*delta_z
                endif
                p_p(i,k) = pi_p_p(i,k)*c_p*rhou(i,k)*thetavb(i,k)
            enddo
        enddo
        
        pi_p_p(:,1) = pi_p_p(:,2)
        p_p(:,1) = p_p(:,2)

        Th_p_c = Th_p_p
        pi_p_c = pi_p_p

    end subroutine base_pertubation
    
end module base_pertubation