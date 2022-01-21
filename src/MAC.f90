program MAC
    use base_state_dry_isentropic
    use base_pertubation
    use model_integration
    use looping_stuff
    use clear_model_arrays

    implicit none
    character(len=3) :: tp_min_str, tp_max_str, dfc_min_str, dfc_max_str
    integer :: tp_min, tp_max, dfc_min, dfc_max
    call getarg(1,tp_min_str)
    call getarg(2,tp_max_str)
    call getarg(3,dfc_min_str)
    call getarg(4,dfc_max_str)
    read(tp_min_str,*) tp_min
    read(tp_max_str,*) tp_max
    read(dfc_min_str,*) dfc_min
    read(dfc_max_str,*) dfc_max

    do theta_p_loop=tp_min,tp_max!1,20
        do dfc_loop=dfc_min,dfc_max!24,160,4
            call base_state_dry_isentropic
            call base_pertubation
            call model_integration
            call clear_model_arrays
        enddo
    enddo

end program MAC