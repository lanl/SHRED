module tasks_type
    use types, only : dp

    implicit none
    
    private 

    public tasks_struct

type tasks_struct
    logical :: density_from_file
    logical :: initial_orbital_free
    logical :: read_orbitals_from_file
    logical :: run_Kubo_Greenwood
    logical :: calc_DOS
    logical :: run_Time_Dependent
    logical :: run_Stopping_Power
    logical :: test_Hamiltonian, test_current, test_S
    logical :: calculate_current
    real(dp) :: wall_time
    logical :: restart
    logical :: Cheby_Filter, RMM_DIIS
endtype

end module