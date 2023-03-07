module system_type
    use types, only : dp

    implicit none
    
    private 

    public system_struct

type system_struct
    integer :: n_kpoints, n_kpoints_local
    integer :: n_spin, n_spin_local !trivial collinear spin
    integer :: n_spinor !spinor components of WF
    integer :: n_elements, n_proj
    integer :: n_atoms
    integer :: n_orbitals, n_orbitals_local !includes buffers
    integer :: n_orbitals_deterministic_lowest
    integer :: n_deterministic, n_find_min, n_deterministic_local
    integer :: n_orbitals_stochastic_full
    integer :: n_stochastic, n_stochastic_local
    integer :: n_orbitals_buffer_lowest, n_orbitals_smoother_lowest

    integer :: deterministic_lowest_start, deterministic_lowest_end
    integer :: deterministic_start, deterministic_end
    integer :: stochastic_full_start, stochastic_full_end
    integer :: stochastic_start, stochastic_end
    integer :: find_min_start, find_min_end
    integer :: stoch_gen_type

    real(dp) :: Temperature, Temperature_cheap
    real(dp) :: density
    real(dp) :: mass
    integer :: spin_type
    real(dp), allocatable :: nelec(:), chemical_potential(:), k_point(:,:)
    real(dp) :: epsatm_all
    logical :: orbital_free=.false.
    integer :: cells(3)
endtype

end module