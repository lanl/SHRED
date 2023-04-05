module Time_Dependent
    use types, only : dp
    use constants, only : i_, pi

    implicit none

    public :: Initialize_Time_Dependent, Time_Dependent_Propagation, Current_stuff, &
     print_energies_pressures, print_positions_and_velocities

    contains

    subroutine Initialize_Time_Dependent(sim)
        use simulation_type, only : simulation_struct
        use Ewald_mod, only : Ewald_Calc, init_ewald, free_ewald
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use External_Field, only : Initialize_External_Field_Potential, Initialize_External_Field_Potential_restart
        use Update_Potentials_and_Energies, only : Calculate_all_static_Energies,Update_electronic_from_orbtials_or_density_TD
        use initialization_finalization,only : initialize_atoms_temperature, initialize_time_dependent_0, &
        finalize_all_paw, initialize_paw, initialize_elements, initialize_atoms, finalize_all_HGH
        use Stopping_Power, only : Stopping_Power_Initialization
        use fft,only : real_to_recip
        use td_type, only : TD_BOMD, TD_RT
        use Non_Local_ion, only : pseudopotential_update
        use parallel_mod, only : parallel_task, parallel_wait
        use Stress_Pressure, only : Stress_and_Pressure
        use allocate_memory, only : allocate_required_memory_atoms_and_elements
        use Density, only : calc_core_Charge_Density
        use inputs, only : read_pseudopotentials
        use Non_Local_ion_EnFoSt, only : NL_PAW_Density_matrix

        type(simulation_struct), intent(inout) :: sim
        integer :: s, i

        if(sim%parallel%myid.eq.0) print *, 'Initializing Time Dependent'; call parallel_wait(sim%parallel); flush(6)
        call initialize_time_dependent_0(sim)

        sim%td%it=0
        if(sim%parallel%myid.eq.0) then
            call initialize_atoms_temperature(sim%atoms, sim%elements, sim%td, sim%system)
        endif
            call parallel_task('bcast', sim%atoms(:)%P(1), sim%parallel, 'all', root=0)
            call parallel_task('bcast', sim%atoms(:)%P(2), sim%parallel, 'all', root=0)
            call parallel_task('bcast', sim%atoms(:)%P(3), sim%parallel, 'all', root=0)

        if(sim%system%orbital_free.and.((sim%td%type.eq.TD_RT).or.(abs(sim%tf%lambda).gt.tiny(1.0_dp)))) then
            if(size(sim%orbitals(1,1,1)%of).gt.1) then
                print *, 'Need to figure out spinors for orbital free time dependent'
                stop
            endif
            do s=1,size(sim%density%of)
                sim%orbitals(1,1,s)%of(1)%R=sqrt(real(sim%density%of(s)%R))/sqrt(sim%system%nelec(s))
                sim%orbitals(1,1,s)%occ=sim%system%nelec(s)
                sim%orbitals(1,1,s)%weight=1.0
                call real_to_recip(sim%orbitals(1,1,s)%of(1), sim%grids)
            enddo
        endif
    
        call Ewald_Calc(sim%atoms, sim%elements, sim%grids, sim%ewald, sim%energies, .true., sim%stress) 
        call Calculate_all_static_Energies(sim%orbitals, sim%grids, sim%system, &
        sim%parallel, sim%density, sim%potentials, sim%energies, sim%atoms,  &
        sim%elements, sim%xc, sim%tf, sim%all_PAW, sim%main_output, 0)
        call Stress_and_Pressure(sim%orbitals, sim%grids, sim%system, sim%parallel, sim%density, &
        sim%potentials, sim%energies, sim%stress, sim%pressures, sim%atoms, sim%elements, &
        sim%xc, sim%tf, sim%all_PAW, sim%main_output, 0)

        !initialize Stopping Power Calculation
        if(sim%stopping%do_stopping) then
            if(sim%parallel%myid.eq.0) print *, 'Initializing Stopping Power'; call parallel_wait(sim%parallel); flush(6)

            call finalize_all_HGH(sim%elements,sim%atoms)
            do i=1,size(sim%elements)
                if(allocated(sim%elements(i)%Ven0G)) deallocate(sim%elements(i)%Ven0G)
                if(allocated(sim%elements(i)%dVen0GdG2)) deallocate(sim%elements(i)%dVen0GdG2)
                if(allocated(sim%elements(i)%counts)) deallocate(sim%elements(i)%counts)
                if(allocated(sim%elements(i)%displs)) deallocate(sim%elements(i)%displs)
            enddo
            call free_ewald(sim%ewald)
            call Stopping_Power_Initialization(sim%orbitals, sim%stopping, sim%density,  &
            sim%system, sim%all_PAW, sim%atoms, sim%elements, sim%parallel, sim %grids, sim%density_t0, sim%orbitals_t0)
            call read_pseudopotentials(sim%elements, sim%grids, sim%system, sim%xc, sim%parallel, sim%all_PAW)
            call initialize_elements(sim%elements(:), &
                                    sim%parallel, sim%grids, sim%all_PAW)
            call initialize_atoms(sim%atoms(sim%stopping%stop_at:sim%stopping%stop_at), &
                                    sim%elements(:), sim%parallel)
            call allocate_required_memory_atoms_and_elements(sim)
            call init_ewald(sim%ewald, sim%grids(1)%box_length, sim%atoms)
            call Ewald_Calc(sim%atoms, sim%elements, sim%grids, sim%ewald, sim%energies, .true., sim%stress)
            call pseudopotential_update(sim%atoms, sim%elements, sim%grids(2), sim%grids(1), sim%parallel, sim%all_PAW)   
            if(sim%parallel%myid.eq.0) print *, 'Stopping Power Initialized'; call parallel_wait(sim%parallel)  
        else
            sim%stopping%stop_at=sim%system%n_atoms+1
        endif

        call update_Forces_after_Ewald(sim, print_force=.true.)

        sim%td%last_time = sim%td%time 

        call Initialize_External_Field_Potential(sim%td, sim%potentials, sim%grids)
        if(sim%parallel%myid.eq.0) print *, 'Initial A:', sim%td%A(:)
        if(sim%parallel%myid.eq.0) print *, 'Initial Ex:', sim%td%Ex(:)

        if(sim%parallel%myid.eq.0) then
            if(sim%tasks%calculate_current) then
                open(newunit=sim%current_output_unit,  file="Current.data", status="unknown", recl=1024)
                write(sim%current_output_unit, *) "1:Time", ' 2-4: Current (a.u.)', ' 5-7:Dipole(a.u.)'
            endif

            
            open(newunit=sim%energies%output,  file="Energies.data", status="unknown", recl=1024)
            open(newunit=sim%pressures%output,  file="Pressures.data", status="unknown", recl=1024)

            write(sim%energies%output, *) "1:Time 2:Total Free(eV) 3: Total Pot 4:Hartree 5:-kTEntropy  &
            & 6:xc 7:Ion Local 8:Ion_Non_local 9:Core 10:External &
            & 11:kinetic 12:kinetic_local 13:kinetic_dynamic 14:nuclear_Pot 15:nuclear_kinetic" 

            write(sim%pressures%output, *) "1:Time 2:Total Free(GPa) 3: Total Pot 4:Hartree 5:-kTEntropy  &
            & 6:xc 7:Ion Local 8:Ion_Non_local 9:Core 10:External &
            & 11:kinetic 12:kinetic_local 13:kinetic_dynamic 14:nuclear_Pot 15:nuclear_kinetic" 

        
            open(newunit=sim%positions_output_1,  file="Positions.data", status="unknown", recl=1024)
            open(newunit=sim%positions_output_2,  file="Relative_Positions.data", status="unknown", recl=1024)
            open(newunit=sim%velocities_output,  file="Velocities.data", status="unknown", recl=1024)
            open(newunit=sim%forces_output,  file="Forces.data", status="unknown", recl=1024)  
            open(newunit=sim%trajectory_output,  file="trajectory.data", status="unknown", recl=1024)  

        endif


        if((sim%tasks%calculate_current.or.sim%tf%dynamic_kedp)) then
            call Current_stuff(sim)
        endif

        if(sim%td%type.eq.TD_RT) then
            call  SIL_Initialize(sim%td%Q, sim%orbitals, sim%system, sim%parallel, sim%grids, sim%td%SIL_Rank)
        endif

    end subroutine

    subroutine Time_Dependent_Propagation(sim)
        use simulation_type, only : simulation_struct
        use Non_Local_ion, only : pseudopotential_update
        use Non_Local_ion_EnFoSt, only : NL_PAW_Density_matrix
        use Update_Potentials_and_Energies, only : Update_ion_and_TD_Potentials, Update_electronic_from_orbtials_or_density_TD
        use Stopping_Power, only : Stopping_Power_Calculation
        use td_type, only : TD_RT, TD_BOMD
        use SCF, only : Self_consistent_field
        use Ewald_mod, only : Ewald_Calc
        use Orbital_Free_Min, only : orbital_free_energy_min
        use parallel_mod, only : parallel_task
        use odp_type, only : orbital_struct
        use Density, only : calc_core_Charge_Density
        use initialization_finalization, only : init_den_pot_orb
        use fft, only : real_to_recip
        use constants, only : FMT_int, FMT_real
        use initialization_finalization, only : deterministic_vector_initializer
        type(simulation_struct), intent(inout), target :: sim

        real(dp) :: last_cycle_time=0.0_dp, time_start, delta, time_save
        integer :: it, ix, iy, iz, s, u, i, j
        character(6) :: id1, id2
        real(dp) ::  dt_scale


        dt_scale=1.0_dp
        if(sim%td%type.eq.TD_RT) dt_scale=0.5_dp
    
        do it=sim%td%it, sim%td%nt
            call cpu_time(time_start) 
            if(sim%parallel%myid.eq.0) print *, "Starting Timestep:", it
            sim%atoms(:)%update_me=.not.sim%atoms(:)%frozen_R
            sim%atoms(:)%update_me_force=.not.sim%atoms(:)%frozen_V

            !Velocity Verlet Momentum  step up to 1/2 step
            if(any(.not.sim%atoms(:)%frozen_V)) &
            call VVerlet_step_P(sim%atoms, sim%elements, sim%td, sim%parallel, (sim%td%time-sim%td%last_time))

            !Velocity Verlet Position  step up to 1/2 step
            if(any(.not.sim%atoms(:)%frozen_R)) &
            call VVerlet_step_R(sim%atoms,sim%elements, sim%grids, sim%parallel, &
                    (sim%td%time-sim%td%last_time)*dt_scale)

                    !Velocity Verlet Force adjustment new atom positions, original electrons
            if(.not.all(sim%atoms(:)%frozen_R))   then
                !Calculate new PP, electric field and local ion-potentials 
                    if(.not.all(sim%atoms(:(sim%stopping%stop_at-1))%frozen_R)) &
                        call Ewald_Calc(sim%atoms, sim%elements, sim%grids, sim%ewald, sim%energies, .true.)
                    call pseudopotential_update(sim%atoms, sim%elements,  sim%grids(2), sim%grids(1), sim%parallel, sim%all_PAW)
                    call calc_core_Charge_Density(sim%all_PAW, sim%grids, sim%atoms, sim%elements, sim%parallel)
                    call  NL_PAW_Density_matrix(sim%orbitals(:,:,:), sim%atoms, sim%elements, sim%grids, &
                                                sim%parallel, sim%all_PAW) 
            endif
            time_save=sim%td%time
            sim%td%time=sim%td%last_time + 0.5_dp*(sim%td%time-sim%td%last_time) !Electric field needs to update up to td%time, but may only want half step
            call Update_ion_and_TD_Potentials(sim%grids, sim%system, sim%parallel, sim%potentials, &
                                    sim%td, sim%atoms, sim%elements, sim%all_PAW, sim%fine_to_coarse)
            sim%td%time=time_save

               
            !step the orbitals / electrons
            if(abs(sim%td%time-sim%td%last_time).gt.tiny(1.0_dp)) then
                if(sim%td%type.eq.TD_RT) then
                    call Real_Time_Step(sim,sim%td%time-sim%td%last_time)
                    if(sim%tasks%calculate_current.or.sim%tf%dynamic_kedp) then
                        call Current_stuff(sim)
                    endif
                    call Update_electronic_from_orbtials_or_density_TD( sim%grids, sim%system, sim%parallel, &
                    sim%coarse_density, sim%density, sim%current_density, sim%potentials, sim%energies, sim%atoms, sim%elements,  &
                    sim%xc, sim%tf, sim%all_PAW, sim%fine_to_coarse, sim%coarse_to_fine, sim%orbitals)

                else if(sim%td%type.eq.TD_BOMD) then
                    if(sim%system%orbital_free) then
                        call  orbital_free_energy_min(sim%system, sim%density, sim%all_PAW, sim%grids,  &
                            sim%parallel, sim%potentials, sim%xc, sim%energies, sim%tf)
                            if(abs(sim%tf%lambda).gt.tiny(1.0_dp)) then
                                do s=1,size(sim%density%of)
                                    sim%orbitals(1,1,s)%of(1)%R=sqrt(real(sim%density%of(s)%R))/sqrt(sim%system%nelec(s))
                                    call real_to_recip(sim%orbitals(1,1,s)%of(1), sim%grids)
                                enddo
                            endif
                    else
                        if(mod(it,20).eq.0.and.(it.gt.0)) then
                             call init_den_pot_orb(sim) !start from fresh density to prevent infinite growth of convergence error
                             do i=1,size(sim%orbitals,3); do j=1, size(sim%orbitals,2)
                                call deterministic_vector_initializer(sim%orbitals(:,j,i),sim%grids, sim%parallel)
                             enddo;enddo
                        endif
                        call Self_consistent_field(sim, sim%grids)
                    endif
                endif
            endif

            sim%atoms(:)%update_me=.not.sim%atoms(:)%frozen_R
            sim%atoms(:)%update_me_force=.not.sim%atoms(:)%frozen_V

            !Velocity Verlet Force adjustment new atom positions
            if(sim%td%type.eq.TD_RT) then
                if(any(.not.sim%atoms(:)%frozen_R)) &
                call VVerlet_step_R(sim%atoms,sim%elements, sim%grids, sim%parallel, &
                (sim%td%time-sim%td%last_time)*0.5_dp)
                if(.not.all(sim%atoms(:)%frozen_R))   then
                    !Calculate new PP, electric field and local ion-potentials
                    if(.not.all(sim%atoms(:(sim%stopping%stop_at-1))%frozen_R)) & 
                        call Ewald_Calc(sim%atoms, sim%elements, sim%grids, sim%ewald, sim%energies, .true.)
                    call pseudopotential_update(sim%atoms, sim%elements,  sim%grids(2), sim%grids(1), sim%parallel, sim%all_PAW)
                    call calc_core_Charge_Density(sim%all_PAW, sim%grids, sim%atoms, sim%elements, sim%parallel)
                    call  NL_PAW_Density_matrix(sim%orbitals(:,:,:), sim%atoms, sim%elements, sim%grids, sim%parallel, sim%all_PAW) 
                endif
                call Update_ion_and_TD_Potentials(sim%grids, sim%system, sim%parallel, sim%potentials, &
                sim%td, sim%atoms, sim%elements, sim%all_PAW, sim%fine_to_coarse)
            endif

            call update_Forces_after_Ewald(sim, print_force=.false.)

             !Velocity Verlet Momentum  step up to full step
            if(any(.not.sim%atoms(:)%frozen_V)) &
            call VVerlet_step_P(sim%atoms, sim%elements, sim%td, sim%parallel, (sim%td%time-sim%td%last_time))

            !To-Do: Output Calculations & plots
            if(sim%stopping%do_stopping.and.sim%stopping%remove_t0) then
                if(sim%elements(sim%stopping%stop_el)%n_proj.gt.0) then
                    call Stopping_Power_Calculation(sim%orbitals, sim%stopping, sim%density, sim%potentials, &
                    sim%system, sim%td, sim%atoms, sim%elements, sim%parallel, sim%grids, sim%all_PAW, &
                    sim%density_t0, sim%orbitals_t0)
                else
                    call Stopping_Power_Calculation(sim%orbitals, sim%stopping, sim%density,sim%potentials,  &
                    sim%system, sim%td, sim%atoms, sim%elements, sim%parallel, sim%grids, sim%all_PAW, sim%density_t0)
                endif
            else if(sim%stopping%do_stopping) then
                call Stopping_Power_Calculation(sim%orbitals, sim%stopping, sim%density,sim%potentials,  &
                sim%system, sim%td, sim%atoms, sim%elements, sim%parallel, sim%grids, sim%all_PAW)
            endif

            call print_energies_pressures(sim, it)            
            call print_positions_and_velocities(sim)
            
            sim%td%last_time=sim%td%time
            sim%td%time=min(sim%td%time+sim%td%dt,sim%td%total_time)

            if(sim%parallel%myid.eq.0) then
                call cpu_time(last_cycle_time) 
                last_cycle_time=last_cycle_time-time_start
                print *, 'TD cycle time: ', last_cycle_time
            endif
            
            if(.false.) then
            if((mod(floor(sim%td%time/sim%td%dt),1000).eq.0)) then
                if(sim%parallel%myid_space.eq.0) then
                    write(id1,FMT_int) sim%parallel%my_space_group
                    write(id2,FMT_int) floor(sim%td%time/sim%td%dt)
                    open(newunit=u, file= "density_yz"//trim(id1)//"_"//trim(id2)//".data", status="unknown", recl=1024)
                     if(sim%stopping%do_stopping)then
                        ix=floor(sim%atoms(sim%stopping%stop_at)%R(1)/sim%grids(1)%dR(1)+0.5_dp)
                     else
                        ix=1
                     endif
                    delta=0.0_dp
                    if(sim%stopping%do_stopping.and.sim%stopping%remove_t0) then
                        do s=1,size(sim%density%of(:))
                            delta=delta+real(sim%density%of(s)%R(ix,iy,iz))-real(sim%density_t0%of(s)%R(ix,iy,iz))
                        enddo
                    else
                        do s=1,size(sim%density%of(:))
                            delta=delta+real(sim%density%of(s)%R(ix,iy,iz))
                        enddo
                    endif
                    do iy=1, sim%grids(1)%Nr_local(2); do iz=1, sim%grids(1)%Nr_local(3)
                        write(u, FMT_real) sim%grids(1)%R(ix,iy,iz,2), sim%grids(1)%R(ix,iy,iz,3), delta
                    enddo; enddo
                    close(u)
                endif
            endif
            endif

        enddo

        if(sim%td%type.eq.TD_RT) call SIL_Finalize(sim%td%Q, sim%grids)

    end subroutine


    subroutine VVerlet_step_R(atoms, elements, grids, parallel, dt)
        use atom_type, only : atom_struct 
        use element_type, only : element_struct  
        use parallel_type, only : parallel_struct  
        use parallel_mod, only : parallel_task
        use Thermostat_mod, only : IsoKinetic_Thermostat
        use grids_type, only : grid_struct


        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout), target :: elements(:)
        type(element_struct), pointer :: element
        type(grid_struct), intent(inout) :: grids(:)
        type(parallel_struct), intent(in):: parallel

        real(dp), intent(in) :: dt
        integer :: i        

        !this is a combination and equivalent to the steps:
        !2) R_1-> R_0+P_0/M*dt (+ box cycling)
        if(parallel%myid.eq.0) then
            do i=1,size(atoms)
                element=>elements(atoms(i)%element)
                if(.not.atoms(i)%frozen_R) atoms(i)%R(:)=atoms(i)%R(:)+(atoms(i)%P(:)/element%M)*dt
                atoms(i)%R(:)=atoms(i)%R(:)-grids(1)%box_length(:)*floor(atoms(i)%R(:)/grids(1)%box_length(:))
            enddo
        endif   

        do i=1,size(atoms)
            if(.not.atoms(i)%frozen_R) call parallel_task('bcast', atoms(i)%R(:), parallel, 'all', root=0)
        enddo

    end subroutine

    subroutine VVerlet_step_P(atoms, elements, td, parallel, dt)
        use atom_type, only : atom_struct 
        use element_type, only : element_struct  
        use td_type, only : td_struct  , IsoKinetic
        use parallel_type, only : parallel_struct  
        use parallel_mod, only : parallel_task
        use Thermostat_mod, only : IsoKinetic_Thermostat

        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout), target :: elements(:)
        type(td_struct), intent(in):: td
        type(parallel_struct), intent(in):: parallel
        real(dp), allocatable :: atoms_P(:,:)
        real(dp), intent(in) :: dt

        integer :: i

        allocate(atoms_P(3,size(atoms)))    
        if(parallel%myid.eq.0) then
            do i=1,size(atoms)
                atoms_P(:,i)=atoms(i)%P(:)
            enddo
            if(td%thermostat.eq.IsoKinetic) then
                    call IsoKinetic_Thermostat(atoms, atoms_P, elements, 0.5*dt)
            else
                    do i=1,size(atoms)
                        if(atoms(i)%frozen_V) cycle
                        atoms_P(:,i)=atoms(i)%P(:)+0.5_dp*atoms(i)%F(:)*dt
                    enddo
            endif
        endif
 
        !Should test if sycronizing these is really important or if we can do the atom propagation on each processor and avoid 
        call parallel_task('bcast', atoms_P, parallel, 'all', root=0)
        do i=1,size(atoms)
            atoms(i)%P(:)=atoms_P(:,i)
        enddo
        deallocate(atoms_P)    

        
    end subroutine


    subroutine Real_Time_Step(sim, dt)
        use simulation_type, only : simulation_struct, potentials_struct
        use odp_type, only : field_struct, orbital_struct
        use grids_type, only : grid_struct
        use odp, only : allocate_field, deallocate_field, allocate_orbital, deallocate_orbital
        use element_type, only : PAW
        use Non_Local_Ion, only : Calculate_Projector_overlaps
        use Non_Local_ion_EnFoSt, only :  NL_PAW_Density_matrix
        use Current, only : Current_Density_Calculation
        use Update_Potentials_and_Energies, only : Update_electronic_from_orbtials_or_density_TD
        use grids_mod, only : Gspace_grid_to_grid_transfer
        use fft, only: real_to_recip, recip_to_real
        use parallel_mod, only: parallel_task
        use Density, only : calc_Compensation_Charge_Density
        use current, only : orbital_current

        type(simulation_struct), intent(inout),target :: sim
        real(dp) :: dt, factor
        integer :: i,j,k,s, spin, at, e, ip, jp, step, dir
        type(field_struct), allocatable :: phi_save(:)
        type(grid_struct), pointer :: grid
        type(orbital_struct), allocatable :: J_orbital(:)

        complex(dp), allocatable :: PYpsi(:,:,:)

        if(sim%td%ETRS_steps.gt.0) then  
            do at=1, size(sim%atoms)
                e=sim%atoms(at)%element
                if(sim%elements(sim%atoms(at)%element)%PP_type.ne.PAW) cycle
                if(.not.allocated(sim%atoms(at)%PAW%rho_ij)) cycle 
                sim%atoms(at)%PAW%rho_ij_save=sim%atoms(at)%PAW%rho_ij
            enddo
        endif

    do step=1, sim%td%ETRS_steps + 1
        if(step.le.sim%td%ETRS_steps) then
            do s=1, size(sim%coarse_density_pred%of(:))
                sim%coarse_density_pred%of(s)%R=0.0_dp
            enddo
            do at=1, size(sim%atoms)
                e=sim%atoms(at)%element
                if(sim%elements(sim%atoms(at)%element)%PP_type.ne.PAW) cycle
                if(.not.allocated(sim%atoms(at)%PAW%rho_ij)) cycle 
                sim%atoms(at)%PAW%rho_ij(:,:,:)=0.0_dp
            enddo
            if(sim%tf%dynamic_kedp) then
                allocate(J_orbital(3))
                do dir=1,3;do s=1,sim%current_density(dir)%n_s
                    sim%coarse_current_density_pred(dir)%of(s)%R=0.0_dp
                enddo; enddo
            endif
        endif
    
        do i = 1, size(sim%orbitals,3);do j = 1, size(sim%orbitals,2);do k = 1, size(sim%orbitals,1)
            if(sum(abs(sim%orbitals(k,j,i)%occ(:))).lt.tiny(1.0_dp)) cycle  

            if(sim%tf%dynamic_kedp) then
                do dir=1,3
                    J_orbital(dir)%band=sim%orbitals(k,j,i)%band
                    J_orbital(dir)%k_point=sim%orbitals(k,j,i)%k_point
                    J_orbital(dir)%spin=sim%orbitals(k,j,i)%spin
                    J_orbital(dir)%degeneracy=sim%orbitals(k,j,i)%degeneracy
                    J_orbital(dir)%type=sim%orbitals(k,j,i)%type
                    J_orbital(dir)%n_spinor=sim%orbitals(k,j,i)%n_spinor
                    call allocate_orbital(J_orbital(dir), sim%grids,  &
                                            sim%orbitals(k,j,i)%k_point+2, sim%parallel)
                    J_orbital(dir)%occ=sim%orbitals(k,j,i)%occ
                    J_orbital(dir)%weight=sim%orbitals(k,j,i)%weight
                    J_orbital(dir)%filter=sim%orbitals(k,j,i)%filter
                enddo  
            endif

            if(step.le.sim%td%ETRS_steps) then
                allocate(phi_save(size(sim%orbitals(k,j,i)%of(:))))
                phi_save(:)%grid=sim%orbitals(k,j,i)%of(:)%grid
                grid=>sim%grids(phi_save(1)%grid)
                do s=1, size(sim%orbitals(k,j,i)%of(:))
                    call allocate_field(phi_save(s), grid, sim%parallel)
                    phi_save(s)%R=sim%orbitals(k,j,i)%of(s)%R
                    phi_save(s)%G=sim%orbitals(k,j,i)%of(s)%G
                enddo
                
            endif

            !Main Forward Step
            call SIL_Propagate(sim%orbitals(k,j,i), sim%td%Q(:,:,j,i), sim%td%SIL_rank, sim%grids,  &
            sim%potentials, sim%atoms, sim%elements, &
            sim%parallel, dt, sim%td%err_allowed)   

            
            
            if(step.le.sim%td%ETRS_steps) then !Should this be le?
                !Add contribution to predicted density
                do s=1, size(sim%orbitals(k,j,i)%of(:))
                    spin=s*i
                    sim%coarse_density_pred%of(spin)%R= sim%coarse_density_pred%of(s)%R + &
                        sim%orbitals(k,j,i)%weight(s)*sim%orbitals(k,j,i)%occ(s)*&
                            abs(sim%orbitals(k,j,i)%of(s)%R)**2
                enddo
                !Add contribution to predicted current density
                if(sim%tf%dynamic_kedp) then
                    call orbital_current(sim%orbitals(k,j,i), J_orbital, sim%parallel, sim%potentials, sim%atoms, sim%elements, &
                                         sim%all_PAW, sim%grids)
                    do s=1,size(sim%orbitals(k,j,i)%of)
                        !i and s are both spin variables, only one should be > 1
                        sim%coarse_current_density_pred(dir)%of(s*i)%R =  sim%coarse_current_density_pred(dir)%of(s*i)%R + &
                                                    real(conjg(sim%orbitals(k,j,i)%of(s)%R)*J_orbital(dir)%of(s)%R)*&
                                                    sim%orbitals(k,j,i)%occ(s)*sim%orbitals(k,j,i)%weight(s)
                    enddo
                endif
                !Restore orbital
                do s=1, size(sim%orbitals(k,j,i)%of(:))
                    sim%orbitals(k,j,i)%of(s)%R=phi_save(s)%R
                    sim%orbitals(k,j,i)%of(s)%G=phi_save(s)%G
                    call deallocate_field(phi_save(s),sim%grids(phi_save(s)%grid))
                enddo
                deallocate(phi_save)
            endif
            if(sim%tf%dynamic_kedp) then
                do dir=1,3
                    call deallocate_orbital(J_orbital(dir),sim%grids)
                enddo
            endif
        enddo;enddo;enddo;
        if(sim%tf%dynamic_kedp) deallocate(J_orbital)


         !compute Predition of Potentials if there is another ERTS step to do
         if(step.le.sim%td%ETRS_steps) then
            !Need to average predicted densities w/ original densities 
            !We could save densities from previous ERTS steps also and mix them to help convergence,
            !But for now I see it as just setting the number of steps as a fixed input

            do s=1, size(sim%coarse_density_pred%of(:))
                call parallel_task('sum', sim%coarse_density_pred%of(s)%R(:,:,:), sim%parallel, 'space') 
                sim%coarse_density_pred%of(s)%R(:,:,:)=(sim%coarse_density_pred%of(s)%R(:,:,:)+ &
                                                        sim%coarse_density     %of(s)%R(:,:,:))*0.5_dp
                call real_to_recip(sim%coarse_density_pred%of(s), sim%grids)
                sim%coarse_density_pred%of(s)%G=sim%coarse_density_pred%of(s)%G*sim%grids(2)%cutden
                call recip_to_real(sim%coarse_density_pred%of(s), sim%grids)
                sim%coarse_density_pred%of(s)%R=abs(sim%coarse_density_pred%of(s)%R) ! G cutting may cause negative if near zero
                call Gspace_grid_to_grid_transfer(sim%grids(2), sim%grids(1), sim%parallel, sim%coarse_density_pred%of(s)%G, &
                    sim%density%of(s)%G, sim%coarse_to_fine)
                call recip_to_real(sim%density%of(s), sim%grids)
                 ! add the compensation charge density if needed
                sim%density%of(s)%R=abs(sim%density%of(s)%R) ! G cutting may cause negative if near zero
                if(sim%parallel%myid.eq.0) print *, 'ERTS_step', step, 'Nelectrons of spin ', s, ' by density integeration', &
                    sim%density%of(s)%G(1,1,1)*product(sim%grids(1)%box_length)
            enddo

            if(sim%tf%dynamic_kedp) then
                do dir=1,3; do s=1,sim%coarse_current_density_pred(dir)%n_s
                    call parallel_task('sum', sim%coarse_current_density_pred(dir)%of(s)%R(:,:,:), sim%parallel, 'space') 
                    sim%coarse_current_density_pred(dir)%of(s)%R=(sim%coarse_current_density_pred(dir)%of(s)%R+ &
                                                                  sim%coarse_current_density(dir)%of(s)%R)*0.5_dp
                    call real_to_recip(sim%coarse_current_density_pred(dir)%of(s), sim%grids)
                    sim%coarse_current_density_pred(dir)%of(s)%G=sim%coarse_current_density_pred(dir)%of(s)%G*sim%grids(2)%cutden
                    call recip_to_real(sim%coarse_current_density_pred(dir)%of(s), sim%grids)
                    call Gspace_grid_to_grid_transfer(sim%grids(2), sim%grids(1), sim%parallel, &
                    sim%coarse_current_density_pred(dir)%of(s)%G, sim%current_density(dir)%of(s)%G, sim%coarse_to_fine)
                    call recip_to_real(sim%current_density(dir)%of(s), sim%grids)
                enddo; enddo
            endif


            !Calculate new potentials & PAW Dij' from average densities (not from orbitals)
            call Update_electronic_from_orbtials_or_density_TD(sim%grids, sim%system, sim%parallel, &
            sim%coarse_density_pred, sim%density, sim%current_density, sim%potentials, sim%energies, &
            sim%atoms, sim%elements, sim%xc, sim%tf, sim%all_PAW, sim%fine_to_coarse, sim%coarse_to_fine)
         endif
         
    enddo
        !


    end subroutine

    subroutine SIL_Propagate(orbital, Q, SIL_rank, grids, potentials, atoms, elements, parallel, dt, err_allowed)
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use odp_type, only : orbital_struct, field_struct
        use odp, only : allocate_field, deallocate_field
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct

        type(parallel_struct), intent(in) :: parallel
        type(orbital_struct), intent(inout) :: orbital
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout),target :: grids(:)
        type(field_struct), intent(inout) :: Q(:,:)
        integer, intent(inout) :: SIL_rank

        real(dp),intent(in) :: dt, err_allowed
    
        call SIL_NC_Propagate(orbital, Q, SIL_rank, grids, potentials, atoms, elements, parallel, dt, err_allowed)

    end subroutine

    subroutine SIL_NC_Propagate(orbital, Q, SIL_rank, grids, potentials, atoms, elements, parallel, dt, err_allowed)
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use odp_type, only : orbital_struct, field_struct
        use odp, only : allocate_field, deallocate_field

        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only:  recip_to_real
        use linalg, only: deig_tri
        use Apply_Hamiltonian, only: Apply_H
        use operations_3D, only:  integrate_3D_G

        type(parallel_struct), intent(in) :: parallel
        type(orbital_struct), intent(inout) :: orbital
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout),target :: grids(:)
        type(grid_struct), pointer :: grid
        type(field_struct), intent(inout) :: Q(:,:)

        real(dp),intent(in) :: dt, err_allowed

        real(dp), allocatable ::  a(:),U(:,:), b(:), a_tmp(:), b_tmp(:)
        complex(dp), allocatable :: c_m(:)
        integer, intent(inout) :: SIL_rank
        integer :: s,n,m, n2, K_max
        real(dp) :: dt_remain, dt_check
        logical :: first_pass, converged


            !Short Iteratiive Lancoz method
            !Only Apply H needed, like Chebychev, but lowest Eiegenvector isn't required.

            allocate(a(SIL_rank))
            allocate(b(size(a)))
            allocate(U(size(a),size(a)))
            allocate(c_m(size(a)))

            allocate(a_tmp(size(a)))
            allocate(b_tmp(size(a)))

            dt_remain=dt
            dt_check=dt
            first_pass=.true.
            do while(dt_remain.gt.tiny(1.0_dp))
                converged=.false.
                a(:)=0.0_dp
                b(:)=0.0_dp
                K_max=size(a)
                do n=1, size(a)
                    if(n.eq.1) then
                        do s=1, size(orbital%of(:))
                            Q(s,n)%G=orbital%of(s)%G(:,:,:)
                            Q(s,n)%R=orbital%of(s)%R(:,:,:)
                        enddo
                        call Apply_H(Q(:,n), Q(:,n+1), grids, potentials, atoms, elements, parallel, calc_R=.false.)
                    else
                        call Apply_H(Q(:,n), Q(:,n+1), grids, potentials, atoms, elements, parallel, calc_R=.true.)
                    endif
                    do s=1, size(orbital%of(:))
                        grid=>grids(Q(s,n)%grid)
                        a(n)= a(n)+ real(integrate_3D_G(Q(s,n+1)%G*conjg(Q(s,n)%G), grid, parallel))
                    enddo
                    if(n.lt.size(a)) then
                        if(size(a).gt.1) then
                            do s=1, size(orbital%of(:))
                                grid=>grids(Q(s,n)%grid)
                                Q(s,n+1)%G=Q(s,n+1)%G-a(n)*Q(s,n)%G
                                if(n.gt.1) Q(s,n+1)%G=Q(s,n+1)%G-b(n-1)*Q(s,n-1)%G
                                b(n)=b(n)+real(integrate_3D_G(Q(s,n+1)%G*conjg(Q(s,n+1)%G), grid, parallel))
                            enddo
                            b(n)=sqrt(b(n))
                            do s=1, size(orbital%of(:))
                                Q(s,n+1)%G=Q(s,n+1)%G/b(n)
                            enddo
                        endif
                    endif
                    if(n.eq.1) cycle
                    a_tmp=a(:)
                    b_tmp=b(:)
                    call deig_tri(a_tmp(:n),b_tmp(:(n-1)),U(:n,:n))!a becomes eigenvalues, U holds eigenvectors in the Q subpace
                    do m=1,n
                        c_m(m)=0.0_dp
                        do n2=1,n
                            c_m(m)= c_m(m) + U(1,n2)*exp(-i_*dt_check*a_tmp(n2))*U(m,n2)
                        enddo
                    enddo
                    if(abs(c_m(n)).lt.err_allowed) then
                        converged=.true.
                        K_max=n
                        exit
                    endif
                enddo
    
                do while(first_pass.and..not.converged) !only reduce the dt_size based on first application of the Hamiltonian
                    if(abs(c_m(size(a))).gt.err_allowed) then
                        if(parallel%myid_k.eq.0) print *, c_m(:)
                        if(parallel%myid_k.eq.0) print *, '|c_max|=', abs(c_m(size(a))), err_allowed, dt_check
                        dt_check=dt_check*0.5
                    else
                        first_pass=.false.
                        exit
                    endif
                    if(dt_check.lt.0.0624*dt) then
                        print *, 'Warning: SIL step down to 1/16th dt'
                        first_pass=.false.
                        exit
                    endif
                    do m=1,size(a)
                        c_m(m)=0.0_dp
                        do n=1,size(a)
                            c_m(m)= c_m(m) + U(1,n)*exp(-i_*dt_check*a_tmp(n))*U(m,n)
                        enddo
                    enddo
                enddo

                do s=1, orbital%n_spinor
                    orbital%of(s)%G=0.0_dp
                    do m=1,K_max
                        orbital%of(s)%G=orbital%of(s)%G + c_m(m)*Q(s,m)%G
                    enddo
                    orbital%of(s)%G=orbital%of(s)%G*grids(orbital%of(s)%grid)%cutwf
                    call recip_to_real(orbital%of(s), grids)
                enddo

                dt_remain=dt_remain-dt_check

            enddo
                
            

            deallocate(a)
            deallocate(b)
            deallocate(U)
            deallocate(c_m)

    end subroutine

    subroutine SIL_Initialize(Q, orbitals, system, parallel, grids, SIL_Rank)
        use odp, only : allocate_field
        use odp_type, only : orbital_struct, field_struct, eigen_vector, time_dependent
        use system_type, only : system_struct
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct


        type(system_struct), intent(in) :: system
        type(parallel_struct), intent(in) :: parallel
        type(field_struct), intent(inout), allocatable :: Q(:,:,:,:)
        type(orbital_struct), intent(inout) :: orbitals(:,:,:)
        type(grid_struct), intent(inout),target :: grids(:)
        type(grid_struct), pointer :: grid

        integer, intent(in) :: SIL_Rank
        integer i,j,n,s

        
        if(.not.allocated(Q)) then
            allocate(Q(system%n_spinor, SIL_Rank+1, system%n_kpoints / parallel%n_k_groups, &
                       system%n_spin/ parallel%n_spin_groups))

            do i=1, size(Q,4); do j=1, size(Q,3); do n=1, size(Q,2) ;do s=1, size(Q,1)
                    Q(s,n,j,i)%grid=orbitals(1,j,i)%of(s)%grid
                     grid=>grids(Q(s,n,j,i)%grid)
                    call allocate_field(Q(s,n,j,i), grid, parallel)
            enddo; enddo; enddo; enddo
        endif 

        do i=1, size(orbitals,3); do j=1, size(orbitals,2); do n=1, size(orbitals,1)
            if(orbitals(n,j,i)%td_type.eq.eigen_vector) orbitals(n,j,i)%td_type=time_dependent
        enddo; enddo; enddo
    end subroutine

    subroutine SIL_Finalize(Q,grids)
        use odp_type, only : field_struct
        use odp, only : deallocate_field
        use grids_type, only : grid_struct
        type(grid_struct), intent(inout) :: grids(:)

        type(field_struct), intent(inout), allocatable :: Q(:,:,:,:)
        integer i,j,n,s

        if(allocated(Q)) then
            do i=1, size(Q,4); do j=1, size(Q,3); do n=1, size(Q,2) ;do s=1, size(Q,1)
                    call deallocate_field(Q(s,n,j,i),grids(Q(s,n,j,i)%grid))
            enddo; enddo; enddo; enddo
            deallocate(Q)
        endif   
    end subroutine

    subroutine update_Forces_after_Ewald(sim, print_force)
        use simulation_type, only : simulation_struct
        use Local_Ion, only : calculate_local_ion_force
        use Non_Local_ion_EnFoSt, only : NL_PP_Force_RL
        use Density, only : Force_or_Stress_Comp_charge, calculate_core_xc_force
        use parallel_mod, only: parallel_wait
        type(simulation_struct), intent(inout) :: sim
        logical, intent(in) :: print_force
        real(dp) :: FCM(3), stress(3,3)
        integer :: n_not_frozen, i

        if(any(.not.sim%atoms(:)%frozen_V)) then
            n_not_frozen=0
            sim%elements(:)%all_frozen_V=.true.
            sim%elements(:)%all_frozen_R=.true.
            do i=1,size(sim%atoms)
                if(.not.sim%atoms(i)%frozen_V) then
                     n_not_frozen= n_not_frozen + 1
                     sim%elements(sim%atoms(i)%element)%all_frozen_V=.false.
                endif
                if(.not.sim%atoms(i)%frozen_R) then
                    n_not_frozen= n_not_frozen + 1
                    sim%elements(sim%atoms(i)%element)%all_frozen_R=.false.
               endif
            enddo
            call  calculate_local_ion_force(sim%density, sim%system, sim%atoms, sim%elements, sim%parallel, sim%grids)
            sim%atoms(:)%F_nonloc(1)=0.0_dp
            sim%atoms(:)%F_nonloc(2)=0.0_dp
            sim%atoms(:)%F_nonloc(3)=0.0_dp
            call NL_PP_Force_RL(sim%orbitals, sim%atoms, sim%elements, sim%potentials, sim%grids, sim%parallel, sim%all_PAW)
            sim%atoms(:)%F(1)=sim%atoms(:)%F_nonloc(1) + sim%atoms(:)%F_loc(1) + sim%atoms(:)%F_nn(1) 
            sim%atoms(:)%F(2)=sim%atoms(:)%F_nonloc(2) + sim%atoms(:)%F_loc(2) + sim%atoms(:)%F_nn(2) 
            sim%atoms(:)%F(3)=sim%atoms(:)%F_nonloc(3) + sim%atoms(:)%F_loc(3) + sim%atoms(:)%F_nn(3)

            if(sim%all_PAW%N_PAW_atoms.gt.0) then
                call Force_or_Stress_Comp_charge(sim%all_PAW, sim%potentials, sim%grids,sim%atoms, sim%parallel, &
                    calc_F=.true.,calc_S=.false.,stress=stress)
                call calculate_core_xc_force(sim%potentials, sim%system, sim%atoms, sim%elements, sim%parallel, sim%grids)
                sim%atoms(:)%F(1)=sim%atoms(:)%F(1) + sim%atoms(:)%F_comp(1) + sim%atoms(:)%F_xc_core(1)
                sim%atoms(:)%F(2)=sim%atoms(:)%F(2) + sim%atoms(:)%F_comp(2) + sim%atoms(:)%F_xc_core(2)
                sim%atoms(:)%F(3)=sim%atoms(:)%F(3) + sim%atoms(:)%F_comp(3) + sim%atoms(:)%F_xc_core(3)
            endif

            FCM(1)=sum(sim%atoms(:)%F(1)) 
            FCM(2)=sum(sim%atoms(:)%F(2)) 
            FCM(3)=sum(sim%atoms(:)%F(3)) 
          !  if(sim%parallel%myid.eq.0) print *, 'FCM: ', FCM(:)

        if(print_force) then
           do i=1,size(sim%atoms)
                if(sim%parallel%myid.eq.0) print *, 'Floc reduced: ', i, sim%atoms(i)%F_loc(:)*sim%grids(1)%Box_Length(:)
                if(sim%parallel%myid.eq.0) print *, 'Fnonloc: ', i, sim%atoms(i)%F_nonloc(:)*sim%grids(1)%Box_Length(:)
                if(sim%parallel%myid.eq.0) print *, 'Fnn: ', i, sim%atoms(i)%F_nn(:)*sim%grids(1)%Box_Length(:)
            if(sim%all_PAW%N_PAW_atoms.gt.0) then
                if(sim%parallel%myid.eq.0) print *, 'F_comp: ', i, sim%atoms(i)%F_comp(:)*sim%grids(1)%Box_Length(:)
                if(sim%parallel%myid.eq.0) print *, 'F_xc_core: ', i, sim%atoms(i)%F_xc_core(:)*sim%grids(1)%Box_Length(:)
            endif
            if(sim%parallel%myid.eq.0) print *, 'Floc (eV/A): ', i, sim%atoms(i)%F_loc(:)*51.42208619083232_dp
            if(sim%parallel%myid.eq.0) print *, 'Fnonloc (eV/A): ', i, sim%atoms(i)%F_nonloc(:)*51.42208619083232_dp
            if(sim%parallel%myid.eq.0) print *, 'Fnn (eV/A): ', i, sim%atoms(i)%F_nn(:)*51.42208619083232_dp
                if(sim%all_PAW%N_PAW_atoms.gt.0) then
            if(sim%parallel%myid.eq.0) print *, 'F_comp (eV/A): ', i, sim%atoms(i)%F_comp(:)*51.42208619083232_dp
            if(sim%parallel%myid.eq.0) print *, 'F_xc_core (eV/A): ', i, sim%atoms(i)%F_xc_core(:)*51.42208619083232_dp
                endif
                if(sim%parallel%myid.eq.0) print *, 'F_total(eV/A): ', i, sim%atoms(i)%F(:)*51.42208619083232_dp

            enddo
        endif
        call parallel_wait(sim%parallel)
        do i=1,size(sim%atoms)
            if(.not.sim%atoms(i)%frozen_V) then
                sim%atoms(i)%F(:)=sim%atoms(i)%F(:)-FCM(:)/n_not_frozen
            endif
        enddo
        if(print_force) then
            do i=1,size(sim%atoms)
        if(sim%parallel%myid.eq.0) print *, 'F_total_CM_removed(eV/A): ', i, sim%atoms(i)%F(:)*51.42208619083232_dp
            enddo
        endif
    endif
    end subroutine

    subroutine Current_stuff(sim)
        use simulation_type, only : simulation_struct
        use Current, only : Current_Density_Calculation
        use Thomas_Fermi, only : calc_Dynamic_KEDP_and_energy
        use operations_3D, only : integrate_3D_R

        type(simulation_struct), intent(inout) :: sim
        real(dp) , allocatable :: total_current(:,:), total_dipole(:,:)
        integer :: i, j

        if(sim%tasks%calculate_current.or.sim%tf%dynamic_kedp) then
            call Current_Density_Calculation(sim%orbitals(:,:,:), sim%parallel, sim%potentials, &
                 sim%atoms, sim%elements, sim%all_PAW, sim%grids, &
                sim%coarse_current_density, sim%current_density, sim%coarse_to_fine)
            allocate(total_current(3, sim%current_density(1)%n_s))
            allocate(total_dipole(3, sim%current_density(1)%n_s))
            do i=1,3;do j=1,sim%current_density(i)%n_s
                total_current(i,j)=real(integrate_3D_R(sim%current_density(i)%of(j)%R, &
                    sim%grids(1), sim%parallel))/product(sim%grids(1)%box_length)
                 total_dipole(i,j)=real(integrate_3D_R(sim%grids(1)%R(:,:,:,i)*sim%density%of(j)%R(:,:,:), &
                  sim%grids(1), sim%parallel))/sim%system%nelec(j)
            enddo;enddo
            if(sim%parallel%myid.eq.0) write(sim%current_output_unit, *) sim%td%time, sum(total_current(1,:)), &
            sum(total_current(2,:)), sum(total_current(3,:)), &
            sum(total_dipole(1,:)),sum(total_dipole(2,:)),sum(total_dipole(3,:))
            if(sim%parallel%myid.eq.0) flush(sim%current_output_unit)
            if(sim%tf%dynamic_kedp) then
                call calc_Dynamic_KEDP_and_energy(sim%potentials%kinetic_dynamic, sim%energies%kinetic_dynamic, &
                    sim%current_density, sim%density, sim%system, sim%grids, sim%parallel)
                !No need to put on coarse because it is only for orbital Free which only uses dense grid
                if(sim%parallel%myid.eq.0) print *, 'Dynamic KEDF Energy (a.u.): ' , sim%energies%kinetic_dynamic
            endif
        endif

    end subroutine

    subroutine print_energies_pressures(sim, it)
        use constants, only : Ha2eV, au2GPa
        use simulation_type, only : simulation_struct
        use Update_Potentials_and_Energies, only : Calculate_all_static_Energies
        use Stress_Pressure, only : Stress_and_Pressure
        use Entropy, only : Calc_Entropy
        use Ewald_mod, only : Ewald_Calc
        use constants, only : FMT_real

        type(simulation_struct), intent(inout) :: sim 
        integer :: it


        if(mod(it,sim%td%p_print_mod).eq.0.or.mod(it,sim%td%e_print_mod).eq.0) then
        call Ewald_Calc(sim%atoms, sim%elements, sim%grids, sim%ewald, sim%energies, .true., sim%stress)
        endif
        if(mod(it,sim%td%e_print_mod).eq.0) then
        call Calculate_all_static_Energies(sim%orbitals, sim%grids, sim%system, sim%parallel, sim%density, &
        sim%potentials, sim%energies, sim%atoms, sim%elements, sim%xc, sim%tf, sim%all_PAW, &
        sim%main_output, it)
        call Calc_Entropy(sim)
        endif

        if(mod(it,sim%td%p_print_mod).eq.0) then
        call Stress_and_Pressure(sim%orbitals, sim%grids, sim%system, sim%parallel, &
                sim%density, sim%potentials, sim%energies, sim%stress, sim%pressures, sim%atoms,  &
                sim%elements, sim%xc, sim%tf, sim%all_PAW, sim%main_output, it)
        endif
            
        if(sim%parallel%myid.eq.0) then
            if(sim%all_PAW%N_PAW_atoms.gt.0) then
                if(mod(it,sim%td%e_print_mod).eq.0) &
                write(sim%energies%output, FMT_real) sim%td%time, sim%energies%total*Ha2eV, &
                sim%energies%total*Ha2eV-sim%energies%nuclear_kinetic*Ha2eV, sim%energies%hartree*Ha2eV, &
                sim%entropy*Ha2eV, sim%energies%xc*Ha2eV, sim%energies%ion_local*Ha2eV, &
                (sim%energies%ion_non_local + sim%all_PAW%epaw)*Ha2eV, sim%energies%core*Ha2eV, sim%energies%external_field*Ha2eV, &
                sim%energies%kinetic*Ha2eV, sim%energies%kinetic_local*Ha2eV, sim%energies%kinetic_dynamic*Ha2eV, &
                sim%energies%nuclear_pot*Ha2eV, sim%energies%nuclear_kinetic*Ha2eV
                flush(sim%energies%output)
                if(mod(it,sim%td%p_print_mod).eq.0) &
                write(sim%pressures%output, FMT_real) sim%td%time, sim%pressures%total*au2GPa, sim%pressures%hartree*au2GPa, &
                0.0_dp, sim%pressures%xc*au2GPa, sim%pressures%ion_local*au2GPa, &
                sim%pressures%ion_non_local*au2GPa, sim%pressures%core*au2GPa, sim%pressures%external_field*au2GPa, &
                sim%pressures%kinetic*au2GPa, sim%pressures%kinetic_local*au2GPa, sim%pressures%kinetic_dynamic*au2GPa, &
                sim%pressures%nuclear_pot*au2GPa, sim%pressures%nuclear_kinetic*au2GPa, &
                sim%pressures%xc_core*au2GPa, sim%pressures%comp*au2GPa
                flush(sim%pressures%output)
            else
                if(mod(it,sim%td%e_print_mod).eq.0) &
                write(sim%energies%output, FMT_real) sim%td%time, sim%energies%total*Ha2eV,&
                sim%energies%total*Ha2eV-sim%energies%nuclear_kinetic*Ha2eV, sim%energies%hartree*Ha2eV, &
                sim%entropy*Ha2eV, sim%energies%xc*Ha2eV, sim%energies%ion_local*Ha2eV, &
                sim%energies%ion_non_local*Ha2eV, sim%energies%core*Ha2eV, sim%energies%external_field*Ha2eV, &
                sim%energies%kinetic*Ha2eV, sim%energies%kinetic_local*Ha2eV, sim%energies%kinetic_dynamic*Ha2eV, &
                sim%energies%nuclear_pot*Ha2eV, sim%energies%nuclear_kinetic*Ha2eV
                flush(sim%energies%output)
                if(mod(it,sim%td%p_print_mod).eq.0) &
                write(sim%pressures%output, FMT_real) sim%td%time, sim%pressures%total*au2GPa, sim%pressures%hartree*au2GPa, &
                0.0_dp, sim%pressures%xc*au2GPa, sim%pressures%ion_local*au2GPa, &
                sim%pressures%ion_non_local*au2GPa, sim%pressures%core*au2GPa, sim%pressures%external_field*au2GPa, &
                sim%pressures%kinetic*au2GPa, sim%pressures%kinetic_local*au2GPa, sim%pressures%kinetic_dynamic*au2GPa, &
                sim%pressures%nuclear_pot*au2GPa, sim%pressures%nuclear_kinetic*au2GPa
                flush(sim%pressures%output)
            endif
        endif

        
    end subroutine

    subroutine print_positions_and_velocities(sim)
        use simulation_type, only : simulation_struct
        use constants, only : FMT_real
        type(simulation_struct), intent(inout) :: sim
        integer :: i

        if(sim%parallel%myid.eq.0) then

            write(sim%positions_output_1, *) 'Time (a.u): ', sim%td%time
            write(sim%positions_output_2, *) 'Time (a.u): ', sim%td%time
            write(sim%velocities_output, *) 'Time (a.u): ', sim%td%time
            write(sim%forces_output, *) 'Time (a.u): ', sim%td%time

            do i=1, size(sim%atoms)
                write(sim%positions_output_1, FMT_real) sim%atoms(i)%R(:)
                write(sim%positions_output_2, FMT_real) sim%atoms(i)%R(:)/sim%grids(1)%box_length(:)
                write(sim%velocities_output, FMT_real) sim%atoms(i)%P(:)/sim%elements(sim%atoms(i)%element)%M
                write(sim%forces_output, FMT_real) sim%atoms(i)%F(:)
                write(sim%trajectory_output, FMT_real) sim%td%time, sim%atoms(i)%R(:), &
                sim%atoms(i)%P(:)/sim%elements(sim%atoms(i)%element)%M, sim%atoms(i)%F(:)
            enddo

            write(sim%positions_output_1, *) 
            write(sim%positions_output_2, *) 
            write(sim%velocities_output, *) 
            write(sim%forces_output, *) 

            flush(sim%positions_output_1)
            flush(sim%positions_output_2)
            flush(sim%velocities_output)
            flush(sim%trajectory_output)

        endif
        
    end subroutine

end module
