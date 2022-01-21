module constants_b

    implicit none

    integer, parameter :: &
        NZ = 160,&
        NX = 803

    integer :: &
        l = 0

    real, parameter :: &
        !grid specification
        delta_x = 100,&
        delta_z = 50,&

        !time information
        dt = .25, &
        save_freq = 60,&

        ! for thermal pertubation
        RADX = 2500.0, &
        RADZ = 1000.0,  &
        ZCNT = 0000.0, &
        
        !speed of sound
        c_s = 50.0, &

        !A-R filter info
        IARF = 1,&
        ARF_gamma = .2, &

        !damping info
        damp_tau = 60, & !seconds

        !turbulence info
        k_factor = 1.0/2.0 !1.0/2.0 - standard | 1.0/4.0 - half | 1.0/1.0 - double (diffusion)

    character(len=33), parameter :: &
            file_path = '/squall/nfalk/HFCACPG/out_prod21/'

    integer,parameter :: &
        num_timesteps = 4*60*60,&
        damp_levels = 6
    

    real :: &
        rd2x = 0, &
        rd2z = 0, &
        rdx = 0, &
        rdz = 0, &
        rho_bar = 0.0

end module constants_b