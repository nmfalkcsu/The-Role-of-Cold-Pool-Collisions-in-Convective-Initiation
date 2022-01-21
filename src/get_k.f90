module get_k
        implicit none

        contains

        subroutine get_k()
                use constants
                use constants_b
                use model_arrays

                implicit none

                Kmx(:,:) = delta_x*k_factor
                Kmz(:,:) = delta_z*k_factor
                Khx(:,:) = Kmx*3
                Khz(:,:) = Kmz*3

        end subroutine

end module