module stochastic_type
    use types, only : dp

    implicit none
    
    private 

    public stochastic_struct
    
    type stochastic_struct
        real(dp) :: nelec_calc
        complex(dp), allocatable :: Mn(:), Mn_Dos(:)
        real(dp), allocatable :: Sqrt_Density_Cn(:), Density_Cn(:), OMDensity_Cn(:), Entropy_Cn(:), DOS_Cn(:)
        real(dp) :: Emin, Emax, Ebar, deltaE
        logical :: do_stoc
        integer :: k_point
    end type

end module