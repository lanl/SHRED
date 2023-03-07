module odp_type
    use types, only : op, dp
    use, intrinsic :: iso_c_binding

    implicit none
    
    private 
    
    public :: orbital_struct, spin_DenPot_struct, DenPot_struct, field_struct, deterministic, stochastic, full, &
              lowest, find_min, eigen_vector, time_dependent

    integer, parameter :: deterministic=1, time_dependent=3, eigen_vector=1
    integer, parameter :: stochastic=2
    integer, parameter :: full=1, lowest=1, find_min=0


    
    !Field Struct is the Fourier transformed function in 3D space
    !The data types depend on FFT package, here FFTW3
    type field_struct
        complex(C_DOUBLE_COMPLEX), pointer :: R(:,:,:)
        complex(C_DOUBLE_COMPLEX), pointer :: G(:,:,:)
        real(dp), pointer :: R_gamma(:,:,:)
        type(C_ptr) :: p_R, p_G, p_R_gamma
        integer :: grid
    end type

    !Type of Kohn Sham Orbitals/Bands 
    type orbital_struct
        integer  :: band, k_point, n_spinor, spin, degeneracy
        real(dp), allocatable :: occ(:), eig(:), filter(:), weight(:)
        type(field_struct), allocatable :: of(:)
        integer, allocatable :: seed_array(:)
        complex(dp), allocatable :: of_G_mix(:,:,:,:,:)
        integer :: type, stoc_type, det_type, td_type, color
    endtype

    type DenPot_struct
        type(field_struct) :: of
    endtype

    type spin_DenPot_struct
        integer  :: n_s !number of spin components (functional derivatives may have more than 2)
        type(field_struct), allocatable :: of(:)
    endtype

    

    

end module