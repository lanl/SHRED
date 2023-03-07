program main
    use simulation_type, only: simulation_struct
    use initialization_finalization, only: initialize_sim, finalize_sim, init_den_pot_orb
    use read_in
    use Simple_Tests, only: test_apply_Hamiltonian, test_current
    use SCF, only : Self_consistent_field
    use Time_Dependent, only: Initialize_Time_Dependent, Time_Dependent_Propagation
    use Kubo_Greenwood, only: Mixed_Kubo_Greenwood
    use parallel_mod
    use Density_of_States, only : Calculate_Density_of_States

implicit none

type(simulation_struct) :: sim
real(dp) :: time1, time2

!Initialize Simulation
call initialize_sim(sim)
if(sim%parallel%myid.eq.0) print *, 'Initialization complete'
!call test_grid_transfer(sim)
if(.false.) call test_new_FFT(sim)
if(.false.) stop
!Start performing tasks
if(sim%parallel%myid.eq.0) print *, 'Start tasks'
call init_den_pot_orb(sim)
!Simple Tests can only perform one, they all end in stops
if(sim%tasks%test_Hamiltonian) call test_apply_Hamiltonian(sim)
if(sim%tasks%test_current) call test_current(sim)

if( .not. sim%tasks%read_orbitals_from_file  &
   .and. .not.sim%system%orbital_free) then
    call cpu_time(time1)
    if(sim%parallel%myid.eq.0) print *, 'Initial SCF'
    call Self_consistent_field(sim, sim%grids)
    call cpu_time(time2)
    if(sim%parallel%myid.eq.0) print *, 'SCF Time: ', time2-time1
endif

if(sim%tasks%calc_DOS) then
    call Calculate_Density_of_States(sim)
endif

if(sim%tasks%run_Kubo_Greenwood) then
    call Mixed_Kubo_Greenwood(sim)
endif

if(sim%tasks%run_Time_Dependent) then
    call Initialize_Time_Dependent(sim)
    call Time_Dependent_Propagation(sim)
endif

!Finalize the Simulation
if(sim%parallel%myid.eq.0) print *, 'Begin Finalization: Goodbye'
call parallel_wait(sim%parallel)
call finalize_sim(sim)

end program
