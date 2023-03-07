module inputs
    use types, only : dp
    use simulation_type, only : simulation_struct
    use parallel_mod

    implicit none

    public 

    contains

    subroutine read_input_numerics(numerics_in, tasks, parallel, elements, system)
        use constants, only : Ha2eV,density2gcm3,bohr2ang
        use numerics_type, only : numerics_struct
        use parallel_type, only : parallel_struct
        use tasks_type, only : tasks_struct
        use element_type, only : element_struct, PAW
        use system_type, only : system_struct

        type(parallel_struct), intent(in) ::  parallel
        type(numerics_struct), intent(inout) :: numerics_in
        type(tasks_struct), intent(in) :: tasks
        type(element_struct), intent(in) :: elements(:)
        type(system_struct), intent(in) :: system

        integer :: u
        integer :: lobpcg_inner_steps, diis_n_steps, lobpcg_n_init, lobpcg_PC_type,diis_PC_type
        real(dp) :: lobpcg_soft_lock_thresh
        integer :: pulay_n, pulay_k, pulay_max_iter, cheby_inner_steps, n_scf_delay
        real(dp) ::  pulay_alpha, pulay_L2_eps, pulay_eig_eps, pulay_L2, diis_min_resnorm, kerker_A, kerker_B, kerker_C
        real(dp) :: kerker_Amin
        logical :: precondition_density, kerker_plasma
        namelist/numerics/ lobpcg_inner_steps, lobpcg_soft_lock_thresh, pulay_n, pulay_k, pulay_max_iter, &
        pulay_alpha, pulay_L2_eps, pulay_eig_eps, pulay_L2, cheby_inner_steps, precondition_density, &
        diis_n_steps, lobpcg_n_init, lobpcg_PC_type, diis_PC_type, diis_min_resnorm, kerker_A, kerker_B, kerker_C, &
        kerker_plasma, kerker_Amin, n_scf_delay
     
        lobpcg_soft_lock_thresh=1.0E-8_dp
        diis_min_resnorm=1.0E-8_dp
        lobpcg_inner_steps=8
        pulay_n=5
        pulay_k=3
        pulay_max_iter=100
        pulay_alpha=-1.0_dp
        pulay_L2_eps=1.0E-5_dp
        pulay_eig_eps=1.0E-6_dp
        precondition_density=.true.
        cheby_inner_steps=8
        diis_n_steps=3
        lobpcg_PC_type=1
        diis_PC_type=1
        kerker_plasma=.false.
        kerker_A=0.6 !~(VASP Amix)/0.7
        if(any(elements(:)%PP_type.eq.PAW)) kerker_A=0.4 !VASP A
        kerker_B=kerker_A/0.1 !VASP Amin, but smoothly cutoff
        kerker_C=-1.0_dp !will be shift from VASP Bmix to compensate for smooth cutoff
        kerker_Amin=0.0_dp
        n_scf_delay=0

        if(tasks%Cheby_Filter) then
            lobpcg_n_init = 1
        else if(tasks%RMM_DIIS) then
            lobpcg_n_init = 4
        else
            lobpcg_n_init = -1
        endif

        if(parallel%myid.eq.0) then
            open(newunit=u, file="input", status="old")
            read(u,nml=numerics)
            close(u)
            numerics_in%lobpcg%inner_steps=lobpcg_inner_steps
            numerics_in%lobpcg%soft_lock_thresh=lobpcg_soft_lock_thresh
            numerics_in%pulay%n=pulay_n
            numerics_in%pulay%k=pulay_k
            numerics_in%pulay%max_iter=pulay_max_iter
            if(pulay_alpha.ge.0) then
                print *, 'Warning pulay_alpha is depricated, just use kerker_A'
            endif
            numerics_in%pulay%n_skip=n_scf_delay
            numerics_in%pulay%L2_eps=pulay_L2_eps
            numerics_in%pulay%eig_eps=pulay_eig_eps
            numerics_in%precondition_density=precondition_density
            numerics_in%cheby%inner_steps=cheby_inner_steps
            numerics_in%diis%n_steps=diis_n_steps
            if(lobpcg_n_init.lt.0.or.((.not.tasks%Cheby_Filter).and.(.not.tasks%RMM_DIIS))) lobpcg_n_init=pulay_max_iter
            numerics_in%lobpcg%n_init=lobpcg_n_init
            numerics_in%lobpcg%PC_type=lobpcg_PC_type
            numerics_in%diis%PC_type=diis_PC_type
            numerics_in%diis%min_resnorm=diis_min_resnorm
            numerics_in%kerker%plasma=kerker_plasma
            numerics_in%kerker%A=kerker_A
            numerics_in%kerker%B=kerker_B
            if(kerker_C.lt.tiny(1.0_dp)) kerker_C=(2.0_dp+kerker_A)/3.0_dp
            numerics_in%kerker%C=kerker_C
            numerics_in%kerker%Amin=kerker_Amin
            if(numerics_in%kerker%plasma) then
                print *, 'Using Plasma Preconditioner (setting kerker_C=(Fermi_velocity)^-1'
            endif
        endif
        call parallel_task('bcast',numerics_in%diis%min_resnorm, parallel, 'all', root=0)
        call parallel_task('bcast',numerics_in%lobpcg%PC_type, parallel, 'all', root=0)
        call parallel_task('bcast',numerics_in%diis%PC_type, parallel, 'all', root=0)
        call parallel_task('bcast',numerics_in%lobpcg%n_init, parallel, 'all', root=0)
        call parallel_task('bcast',numerics_in%cheby%inner_steps, parallel, 'all', root=0)
        call parallel_task('bcast',numerics_in%lobpcg%inner_steps, parallel, 'all', root=0)
        call parallel_task('bcast',numerics_in%lobpcg%soft_lock_thresh, parallel, 'all', root=0)
        call parallel_task('bcast',numerics_in%pulay%n, parallel, 'all', root=0)
        call parallel_task('bcast',numerics_in%pulay%k, parallel, 'all', root=0)
        call parallel_task('bcast',numerics_in%pulay%n_skip, parallel, 'all', root=0)
        call parallel_task('bcast',numerics_in%pulay%max_iter, parallel, 'all', root=0)
        call parallel_task('bcast',numerics_in%pulay%L2_eps, parallel, 'all', root=0)
        call parallel_task('bcast',numerics_in%pulay%eig_eps, parallel, 'all', root=0)
        call parallel_task('bcast',numerics_in%precondition_density, parallel, 'all', root=0)
        call parallel_task('bcast',numerics_in%diis%n_steps, parallel, 'all', root=0)

        call parallel_task('bcast', numerics_in%kerker%plasma, parallel, 'all', root=0)
        call parallel_task('bcast', numerics_in%kerker%A, parallel, 'all', root=0)
        call parallel_task('bcast', numerics_in%kerker%B, parallel, 'all', root=0)
        call parallel_task('bcast', numerics_in%kerker%C, parallel, 'all', root=0)
        call parallel_task('bcast', numerics_in%kerker%Amin, parallel, 'all', root=0)

        if(parallel%myid.eq.0) then
            print *, 'Numerical Methods Values: '
            if(tasks%Cheby_Filter) then
                print *, 'Deterministic Eigenvector Chebychev Filtering===='
                print *, 'Expansion Terms: ',numerics_in%cheby%inner_steps
                print *, 'Initial LOBPCG===='
                print *, 'Number of initial LOBPCG iterations before Cheby:', numerics_in%lobpcg%n_init
            else if(tasks%RMM_DIIS) then
                print *, 'Deterministic Eigenvector RMM DIIS===='
                print *, '# DIIS steps (max Rank of Risiduals overlap matrix): ', numerics_in%diis%n_steps
                if(numerics_in%diis%n_steps.lt.2) then 
                    numerics_in%diis%n_steps=2
                    print *, '# DIIS steps Set to at least 2: ', numerics_in%diis%n_steps
                endif
                print *, 'Initial LOBPCG===='
                print *, 'Number of initial LOBPCG iterations before RMM DIIS:', numerics_in%lobpcg%n_init
            else
                print *, 'LOBPCG===='
            endif
            print *, 'Inner Iterations: ', numerics_in%lobpcg%inner_steps
            print *, 'Soft Locking Threshold: ', numerics_in%lobpcg%soft_lock_thresh
            print *, 'SCF Pulay Mixing===='
            print *, 'Pulay n : ', numerics_in%pulay%n
            print *, 'Pulay k : ', numerics_in%pulay%k
            print *, 'Max SCF Iterations: ', numerics_in%pulay%max_iter
            !print *, 'Pulay alpha : ', numerics_in%pulay%alpha
            print *, 'Pulay Skip mixing for first # steps: ',numerics_in%pulay%n_skip
            print *, 'Pulay potential L2_norm convegence : ',numerics_in%pulay%L2_eps
            print *, 'Pulay eigenvalues convegence : ',numerics_in%pulay%eig_eps
            print *, 'Using Preconditioner : K(|G|)=MAX(A_MIN,A*(B^-1 + C^2*G^2)/(1+C^2*G^2))'
            if(numerics_in%kerker%plasma) then
                print *, 'A:', numerics_in%kerker%A, 'B:', numerics_in%kerker%B, 'C: inverse Fermi momentum'
            else
                print *, 'A:', numerics_in%kerker%A,'B:', numerics_in%kerker%B,'C:', numerics_in%kerker%C
                print *, 'ABinit Equivalence-> ', ' diemac: ', numerics_in%kerker%B, &
                                                  ' dielng: ', numerics_in%kerker%C
                print *, 'VASP Equivalence-> ', 'Bmix: ', 1.0_dp/numerics_in%kerker%C**2
                print *, 'Smooth minimum (A/B): ', numerics_in%kerker%A/numerics_in%kerker%B, &
                         'Hard minimum (Amin):', numerics_in%kerker%Amin
            endif
        endif

        allocate(numerics_in%pulay%R_i(numerics_in%pulay%n+1,system%n_spin))
        allocate(numerics_in%pulay%F_i(numerics_in%pulay%n+1,system%n_spin))
    end subroutine

    subroutine read_input_tasks(tasks, parallel)
        use constants, only : Ha2eV,density2gcm3,bohr2ang
        use tasks_type, only : tasks_struct
        use parallel_type, only : parallel_struct
        type(parallel_struct), intent(in) ::  parallel
        type(tasks_struct), intent(inout) :: tasks
        integer :: u
        logical :: read_density_from_file, initial_orbital_free, read_orbitals_from_file, &
                   run_time_dependent, run_Kubo_Greenwood, run_Stopping_Power, test_Hamiltonian, &
                   calculate_current, test_current, Cheby_Filter, RMM_DIIS, test_S, calc_DOS
        real(dp) :: wall_time_hours

        namelist/task_list/ read_density_from_file, initial_orbital_free, read_orbitals_from_file, &
        run_Time_Dependent, run_Kubo_Greenwood, run_Stopping_Power, test_Hamiltonian, calculate_current, &
        test_current, wall_time_hours, Cheby_Filter,RMM_DIIS, test_S, calc_DOS
     
        read_density_from_file=.false.
        read_orbitals_from_file=.false.
        run_Time_Dependent=.false.
        run_Kubo_Greenwood=.false.
        calc_DOS=.false.

        run_Stopping_Power=.false.
        test_Hamiltonian=.false.
        test_S=.false.
        calculate_current=.false.
        initial_orbital_free=.false.
        test_current=.false.
        Cheby_Filter = .false.
        RMM_DIIS=.false.
        wall_time_hours=1000.0_dp

        if(parallel%myid.eq.0) then
            open(newunit=u, file="input", status="old")
            read(u,nml=task_list)
            close(u)
            tasks%density_from_file=read_density_from_file
            tasks%initial_orbital_free=initial_orbital_free
            tasks%read_orbitals_from_file=read_orbitals_from_file
            tasks%run_Time_Dependent=run_Time_Dependent
            tasks%run_Kubo_Greenwood=run_Kubo_Greenwood
            tasks%run_Stopping_Power=run_Stopping_Power
            tasks%test_Hamiltonian=test_Hamiltonian
            tasks%calculate_current=calculate_current
            tasks%test_current=test_current
            tasks%test_S=test_S
            tasks%wall_time=wall_time_hours*3600 ! to seconds
            tasks%Cheby_Filter=Cheby_Filter
            tasks%RMM_DIIS=RMM_DIIS
            tasks%calc_DOS=calc_DOS
        endif
        
        call parallel_task('bcast',tasks%RMM_DIIS, parallel, 'all', root=0)
        call parallel_task('bcast',tasks%Cheby_Filter, parallel, 'all', root=0)
        call parallel_task('bcast',tasks%density_from_file, parallel, 'all', root=0)
        call parallel_task('bcast',tasks%initial_orbital_free, parallel, 'all', root=0)
        call parallel_task('bcast',tasks%read_orbitals_from_file, parallel, 'all', root=0)
        call parallel_task('bcast',tasks%run_Time_Dependent, parallel, 'all', root=0)
        call parallel_task('bcast',tasks%run_Kubo_Greenwood, parallel, 'all', root=0)
        call parallel_task('bcast',tasks%run_Stopping_Power, parallel, 'all', root=0)
        call parallel_task('bcast',tasks%test_Hamiltonian, parallel, 'all', root=0)
        call parallel_task('bcast',tasks%test_current, parallel, 'all', root=0)
        call parallel_task('bcast',tasks%test_S, parallel, 'all', root=0)
        call parallel_task('bcast',tasks%calculate_current, parallel, 'all', root=0)
        call parallel_task('bcast',tasks%wall_time, parallel, 'all', root=0)
        call parallel_task('bcast',tasks%calc_DOS, parallel, 'all', root=0)


        if(tasks%run_Stopping_Power.and.(.not.tasks%run_Time_Dependent)) then
            print *, 'Stopping Power Requires run_Time_Dependent=.true.'
            tasks%run_Time_Dependent=.true.
        endif


        if(parallel%myid.eq.0) then
            print *, 'Tasks '
            print *, 'Read Density from File:', tasks%density_from_file
            print *, 'Read Orbitals from File:', tasks%read_orbitals_from_file
            print *, 'Initially calculate Orbital Free Density:', tasks%initial_orbital_free

            if(tasks%RMM_DIIS.and.tasks%Cheby_Filter) then
                print *, 'Choose either Cheby_Filter, RMM_DIIS, or neither, not both'
                stop
            endif
            if(tasks%RMM_DIIS) then
                print *, 'Using RMM_DIIS for eigensolver'
            else if (tasks%Cheby_Filter) then
                print *, 'Using Chebychev Filtering for eigensolver'
            else
                print *, 'Using LOBPCG for eigensolver'
            endif

            print *, 'Post (initial) -SCF calculations: '
            print *, 'Running a Time Dependent Calculation:', tasks%run_Time_Dependent
            print *, 'Running a Kubo-Greenwood Calculation:', tasks%run_Kubo_Greenwood
            print *, 'You have allowed ', wall_time_hours, ' real time hours for simulation'  
            print *
        endif
    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads the input file to set the physical system parameters for the 
!simulaitons including:
! kT (Temperature), density 
!Input: None
!Output: None
!Module Variables updated- 
!        system: density, kT, num_atom_types, num_atoms_of_type, PP_type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_input_system(system, parallel)
        use constants, only : Ha2eV,density2gcm3,bohr2ang
        use system_type, only : system_struct
        use parallel_type, only : parallel_struct
        type(parallel_struct), intent(inout) ::  parallel
        type(system_struct), intent(inout) :: system
        real(dp) :: Temperature_eV, Density_gcc, k_shift(3)
        real(dp), allocatable :: kx(:), ky(:), kz(:), k_weight(:)
        logical :: Gamma, Time_Reversal_Symmetry

        integer :: u, Number_of_elements, Number_of_stoch_vect, Number_of_KS_orbitals, &
                      Number_of_buffer_orbitals, Number_of_smoother_orbitals, &
                      N_kpoints, k ,k2,i, x,y,z
        integer :: Spin_type
        integer :: Replica_x, Replica_y, Replica_z,  k_divisions(3) , reason, stoch_gen_type
        character(len=8)   :: header


        namelist/system_parameters/ Temperature_eV, Density_gcc, Number_of_elements, &
        Number_of_KS_orbitals, Number_of_stoch_vect, &
        Spin_type, Number_of_buffer_orbitals, Number_of_smoother_orbitals, &
         N_kpoints, Replica_x, Replica_y, Replica_z, k_divisions, &
        k_shift, Gamma, Time_Reversal_Symmetry, stoch_gen_type
     
        Temperature_eV=-1
        Density_gcc=-1
        Number_of_elements=-1
        Spin_type=1
        Number_of_stoch_vect = 0
        Number_of_KS_orbitals = 0
        Number_of_buffer_orbitals = 0
        Number_of_smoother_orbitals = 0
        stoch_gen_type=1

        N_kpoints = -1
        Gamma=.true.
        k_divisions(:)=-1
        k_shift=0.0_dp
        Time_Reversal_Symmetry=.false.
        system%n_kpoints=0

        Replica_x=0
        Replica_y=0
        Replica_z=0

        if(parallel%myid.eq.0) then
            open(newunit=u, file="input", status="old")
            read(u,nml=system_parameters)
            close(u)   

            !Setting up k-point grid
            if(N_kpoints.gt.0) then
                allocate(kx(N_kpoints))
                allocate(ky(N_kpoints))
                allocate(kz(N_kpoints))
                allocate(k_weight(N_kpoints))

                if(all(k_divisions.gt.0)) then
                    print *, 'You are trying to specify manual K_points (N_k_points>0), and &
                    & a MP grid (k_divisions>0), which you cannot do, stop'
                    stop
                endif
                open(newunit=u, file="input", status="old")
                do
                    read(u,'(a)', IOSTAT=reason) header
                    if(reason.gt.0) then
                        print *, 'Something went wrong looking for k-point list'
                        stop
                    else if (reason.lt.0) then
                        print *, 'No k-point list found, but ',N_kpoints,' manually input k-points specified'
                        stop
                    else
                        if ( &
                            header(1:8).eq.'$kpoints'.or.&
                            header(1:8).eq.'&kpoints' &
                            ) exit ! breaking infinite loop
                    endif
                enddo
                k=0
                do i=1, N_kpoints
                        read(u,*) kx(i), ky(i), kz(i), k_weight(i)
                        print *, kx(i), ky(i), kz(i), k_weight(i)
                enddo
                system%n_kpoints=N_kpoints
                allocate(system%k_point(4,system%n_kpoints))
                do k=1, system%n_kpoints
                        system%k_point(1,k)=kx(k)
                        system%k_point(2,k)=ky(k)
                        system%k_point(3,k)=kz(k)
                        system%k_point(4,k)=k_weight(k)
                enddo
                close(u)
            else if(all(k_divisions.gt.0)) then
                N_kpoints=product(k_divisions)
                allocate(kx(N_kpoints))
                allocate(ky(N_kpoints))
                allocate(kz(N_kpoints))
                allocate(k_weight(N_kpoints))
                k_weight=1.0_dp/N_kpoints
                k=1
                kx=0.0_dp
                ky=0.0_dp
                kz=0.0_dp
                do x=0, k_divisions(1)-1
                    do y=0, k_divisions(2)-1
                        do z=0, k_divisions(3)-1
                            if((mod(k_divisions(1),2).eq.0).and.(.not.Gamma)) then
                                kx(k)=(x-0.5_dp+k_shift(1))/k_divisions(1) 
                            else
                                kx(k)=(x+k_shift(1)-1.0_dp)/k_divisions(1) 
                            endif
                            if((mod(k_divisions(2),2).eq.0).and.(.not.Gamma)) then
                                ky(k)=(y-0.5_dp+k_shift(2))/k_divisions(2) 
                            else
                                ky(k)=(y+k_shift(2)-1.0_dp)/k_divisions(2) 
                            endif
                            if((mod(k_divisions(3),2).eq.0).and.(.not.Gamma)) then
                                kz(k)=(z-0.5_dp+k_shift(3))/k_divisions(3) 
                            else
                                kz(k)=(z+k_shift(3)-1.0_dp)/k_divisions(3)
                            endif
                    k=k+1
                enddo;enddo; enddo

                if(Time_Reversal_Symmetry) then
                    do k=1, N_kpoints
                        do k2=k+1,N_kpoints
                            if((abs(kx(k)+kx(k2)).lt.1E-10).and. &
                               (abs(ky(k)+ky(k2)).lt.1E-10).and. &
                               (abs(kz(k)+kz(k2)).lt.1E-10)) then
                                k_weight(k2)=k_weight(k2)+k_weight(k)
                                k_weight(k)=0.0_dp
                                exit
                            endif
                        enddo
                        if(k_weight(k).gt.tiny(1.0_dp)) then
                            system%n_kpoints=system%n_kpoints+1
                        endif
                    enddo
                    allocate(system%k_point(4,system%n_kpoints))
                    k2=0
                    do k=1, N_kpoints
                        if(k_weight(k).gt.tiny(1.0_dp)) then
                            k2=k2+1
                            system%k_point(1,k2)=kx(k)
                            system%k_point(2,k2)=ky(k)
                            system%k_point(3,k2)=kz(k)
                            system%k_point(4,k2)=k_weight(k)
                            print *, k2, system%k_point(4,k2), system%n_kpoints
                        endif
                    enddo
                    N_kpoints=system%n_kpoints
                else
                    system%n_kpoints=N_kpoints
                    allocate(system%k_point(4,system%n_kpoints))
                    do k=1, system%n_kpoints
                            system%k_point(1,k)=kx(k)
                            system%k_point(2,k)=ky(k)
                            system%k_point(3,k)=kz(k)
                            system%k_point(4,k)=k_weight(k)
                    enddo
                endif
            else
                print *, 'No K-points requested running Gamma point Calculation'
                system%n_kpoints=1
                allocate(system%k_point(4,system%n_kpoints))
                system%k_point(1:3,1)=0.0_dp
                system%k_point(4,1)=1.0_dp
                N_kpoints=1
            endif

            if(abs(sum(system%k_point(4,:))-1.0_dp).gt.1E-10_dp) then
                print *, sum(system%k_point(4,:))
                print *, 'k-point weights do not add to 1'
                stop
            endif                
            system%cells(1)=Replica_x+1
            system%cells(2)=Replica_y+1
            system%cells(3)=Replica_z+1
            print *, 'Replicating atom positions & provided velocities ', Replica_x, 'times in x direction'
            print *, 'Replicating atom positions & provided velocities ', Replica_y, 'times in y direction'
            print *, 'Replicating atom positions & provided velocities ', Replica_z, 'times in z direction'

            if(any(system%cells.lt.1)) then 
                print *, 'Cannot have negative replicas'; stop
            endif

            Number_of_KS_orbitals=Number_of_KS_orbitals*product(system%cells)
            Number_of_buffer_orbitals=Number_of_buffer_orbitals*product(system%cells)
            Number_of_smoother_orbitals=Number_of_smoother_orbitals*product(system%cells)
            !Stochastic orbitals essentially are scaling with system size so no need to increase number

            system%spin_type = Spin_type
            if(system%spin_type.eq.1) print *, 'Restricted Spin'
            if(system%spin_type.eq.2) print *, 'Collinnear Fully UnRestricted Spin'
            if(system%spin_type.eq.3) print *, 'Colinear Restricted-Open Spin (spinor form)'
            if(system%spin_type.eq.4) print *, 'Noncolinear Spin (spinor form)'
            if(system%spin_type.eq.5) then
                print *, 'Dirac-Kohn-Sham, you wish we had that'
                stop
            endif

            if(parallel%n_spin_groups.gt.2) then
                print *, 'the # of parallel spin groups cannot be larger than 2'
                stop
            endif
            if((system%spin_type.eq.1).and.(parallel%n_spin_groups.gt.1)) then
                print *, 'the # of spin groups is larger than 1, but spin type is restricted'
                stop

            endif
            if((system%spin_type.gt.2).and.(parallel%n_spin_groups.gt.1)) then
                print *, 'you have selected a spinor representation'
                print *, 'can only parallelize over collinear, fully unrestricted, spins'
                stop
            endif

            if(system%n_kpoints.lt.parallel%n_k_groups) then
                print *, '# of k-points must be at least as large as the number of k-point groups'
                stop
            endif

            if(Number_of_elements.le.0) then
                print *, 'Number_of_elements needs to be provided, >0'
            endif

            if((Number_of_KS_orbitals.eq.0).and. (Number_of_stoch_vect.eq.0)) then
                print *, 'No Kohn Sham orbitals requested, performing Orbital-Free DFT'
                system%orbital_free=.true.
                if(parallel%n_band_groups.gt.1) then
                    print *, 'Cannot parallelize over bands in an orbital-free calcualtion'
                    stop
                endif
                if(system%n_kpoints.gt.1) then
                    print *, 'K-points do not apply to orbital free calcuations'
                    stop
                endif
                if(Number_of_buffer_orbitals.gt.0) then
                    print *, 'buffer orbitals do not apply to orbital free calcuations'
                    Number_of_buffer_orbitals=0
                endif
                if(Number_of_smoother_orbitals.gt.0) then
                    print *, 'smoother orbitals do not apply to orbital free calcuations'
                    Number_of_smoother_orbitals=0
                endif
            endif

            if((Number_of_KS_orbitals.lt.0).or. (Number_of_stoch_vect.lt.0) &
           .or.(Number_of_buffer_orbitals.lt.0) .or. (Number_of_smoother_orbitals.lt.0)) then
                print *, '# of orbitals cannot be negative, Number_of_KS_orbitals: ', Number_of_KS_orbitals
                print *, '# of orbitals cannot be negative, Number_of_stoch_vect: ', Number_of_stoch_vect
                print *, '# of orbitals cannot be negative, Number_of_buffer_orbitals: ',Number_of_buffer_orbitals
                print *, '# of orbitals cannot be negative, Number_of_smoother_orbitals: ', Number_of_smoother_orbitals
                stop
            endif

            if(Number_of_smoother_orbitals.gt.Number_of_KS_orbitals) then
                print *, 'cannot smooth over more orbitals than total, Number_of_KS_orbitals: ', Number_of_KS_orbitals, &
                'Number_of_smoother_orbitals: ', Number_of_smoother_orbitals
            endif

            if(mod(Number_of_stoch_vect,parallel%n_band_groups).ne.0) then
                Number_of_stoch_vect=(Number_of_stoch_vect/parallel%n_band_groups)*parallel%n_band_groups &
                                     + parallel%n_band_groups
            endif

            system%n_orbitals_deterministic_lowest=Number_of_KS_orbitals+Number_of_buffer_orbitals 

            if(Number_of_stoch_vect.gt.0.and.system%n_orbitals_deterministic_lowest.eq.0) then
                !add a buffer for min eigenvalue calculation
                system%n_find_min=parallel%n_band_groups_sub
            else
                system%n_find_min=0
            endif
            
            if(mod(system%n_orbitals_deterministic_lowest,parallel%n_band_groups).ne.0) then
                print *, 'Increasing buffer from ', Number_of_buffer_orbitals
                system%n_orbitals_deterministic_lowest= &
                    (system%n_orbitals_deterministic_lowest/parallel%n_band_groups)*parallel%n_band_groups &
                    + parallel%n_band_groups
                Number_of_buffer_orbitals=system%n_orbitals_deterministic_lowest-Number_of_KS_orbitals
                print *, 'to ', Number_of_buffer_orbitals, ' due to # of groups requested'
            endif


           

            system%n_deterministic= system%n_orbitals_deterministic_lowest
            system%n_orbitals_buffer_lowest=Number_of_buffer_orbitals
            system%n_orbitals_smoother_lowest=Number_of_smoother_orbitals
            system%n_orbitals_stochastic_full = Number_of_stoch_vect
            system%n_stochastic = system%n_orbitals_stochastic_full
            system%n_orbitals=system%n_deterministic+system%n_stochastic

           !if(parallel%my_sub.eq.0) then 
           !     system%n_orbitals=system%n_orbitals+system%n_find_min
            !endif
            system%Temperature=Temperature_eV/Ha2eV
            system%density=Density_gcc/density2gcm3
            system%n_elements=Number_of_elements

            system%stoch_gen_type=stoch_gen_type

        endif

        call parallel_task('bcast',system%n_kpoints, parallel, 'all', root=0)
        call parallel_task('bcast',system%Temperature, parallel, 'all', root=0)
        call parallel_task('bcast',system%density, parallel, 'all', root=0)
        call parallel_task('bcast',system%n_elements, parallel, 'all', root=0)

        call parallel_task('bcast',system%n_orbitals_deterministic_lowest, parallel, 'all', root=0)
        call parallel_task('bcast',system%n_orbitals_stochastic_full, parallel, 'all', root=0)
        call parallel_task('bcast',system%n_orbitals_buffer_lowest, parallel, 'all', root=0)
        call parallel_task('bcast',system%n_orbitals_smoother_lowest, parallel, 'all', root=0)
        
        call parallel_task('bcast',system%n_deterministic, parallel, 'all', root=0)
        call parallel_task('bcast',system%n_stochastic, parallel, 'all', root=0)

        call parallel_task('bcast',system%n_orbitals, parallel, 'all', root=0)
        call parallel_task('bcast',system%n_find_min, parallel, 'all', root=0)


        call parallel_task('bcast',system%spin_type , parallel, 'all', root=0)
        call parallel_task('bcast',system%orbital_free , parallel, 'all', root=0)
        call parallel_task('bcast',system%cells, parallel, 'all', root=0)

        call parallel_task('bcast',system%stoch_gen_type, parallel, 'all', root=0)


        if(parallel%myid.ne.0) then
            allocate(system%k_point(4,system%n_kpoints))
        endif
        
        call parallel_task('bcast',system%k_point, parallel, 'all', root=0)

        if(system%spin_type.eq.1) then
            system%n_spinor=1
            system%n_spin=1
        else if(system%spin_type.eq.2) then
            system%n_spin=2
            system%n_spinor=1
        else
            print *, 'Spin type not available'
            system%n_spin=1
            !sim%system%n_spinor >= 2
            !Not implemented yet
            stop 
        endif

        allocate(system%chemical_potential(system%n_spin))
        system%chemical_potential=0.0_dp
        
        system%n_spin_local=system%n_spin/parallel%n_spin_groups
        system%n_kpoints_local=system%n_kpoints/parallel%n_k_groups
        system%n_orbitals_local=system%n_orbitals/parallel%n_band_groups
        system%n_deterministic_local=system%n_deterministic/parallel%n_band_groups
        system%n_stochastic_local=system%n_stochastic/parallel%n_band_groups


        if(parallel%myid.eq.0) then
            print *, 'System Options'
            print *, 'Temperature (eV):', system%Temperature*Ha2eV
            print *, '# of element types:', system%n_elements
            print *, '#lowest Eigenspectrum'
            print *, '# of deterministic KS orbitals:',system%n_orbitals_deterministic_lowest
            print *, '# of buffer:',system%n_orbitals_buffer_lowest
            print *, '# of smoothing orbitals:',system%n_orbitals_smoother_lowest

            print *, '# of stochastic vectors (full):', system%n_orbitals_stochastic_full

            print *, '# orbitals total:', system%n_orbitals
            print *, '# orbitals (no buffer):', system%n_orbitals - system%n_orbitals_buffer_lowest


            print *, 'Spin Type:', system%spin_type
            print *, 'number of k-points:', system%n_kpoints
            do k=1, system%n_kpoints
                print *, 'k: ', k, system%k_point(:,k)
            enddo
            print *
         endif

         system%n_proj=0
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads the input file to set the atom types (pseudopotentials), bumber of such atoms, and the
! pseudopotential application method for the simulaitons including:
! Number_of_atoms_of_element, PP_file, PP_type, Projector_Application
!Input: element, n_elements
!Output: None
!Module Variables updated- 
!        system: num_atoms_of_type(:), PP_type(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine read_input_elements(elements, n_elements, system, parallel, tasks)
        use constants, only : u2au, density2gcm3
        use parallel_mod, only: parallel_wait
        use element_type, only : element_struct, Real_local
        use parallel_type, only : parallel_struct
        use system_type, only : system_struct
        use grids_type, only : grid_struct
        use tasks_type, only : tasks_struct

        type(parallel_struct), intent(in) ::  parallel
        type(element_struct), intent(inout) :: elements(:)
        type(system_struct), intent(inout) :: system
        type(tasks_struct), intent(in) :: tasks
        
        integer, intent(in) :: n_elements
        integer, allocatable :: Number_of_atoms_of_element(:)
        character(len=256), allocatable :: PP_file(:), PP_type(:), Projector_Application(:) 
        real(dp), allocatable :: Mass_amu(:), radius(:) 
        integer, allocatable :: OH_factor(:)
        logical, allocatable :: frozen_positions(:), frozen_velocities(:)
        logical, allocatable :: thawed_positions(:), thawed_velocities(:)
        logical, allocatable:: recip_apply(:)

        integer :: i, u
        namelist/element_list/ Number_of_atoms_of_element, PP_file, PP_type, Projector_Application, Mass_amu, &
                               radius, OH_factor, frozen_positions, frozen_velocities, thawed_positions, thawed_velocities, &
                               recip_apply

     

        allocate(Number_of_atoms_of_element(n_elements))
        allocate(PP_file(n_elements))
        allocate(PP_type(n_elements))
        allocate(Projector_Application(n_elements))
        allocate(Mass_amu(n_elements))
        allocate(radius(n_elements))
        allocate(OH_factor(n_elements))
        allocate(frozen_positions(n_elements))
        allocate(frozen_velocities(n_elements))
        allocate(thawed_positions(n_elements))
        allocate(thawed_velocities(n_elements))
        allocate(recip_apply(n_elements))

        radius=3.5_dp
        OH_factor=2!5
        Projector_Application(:)='None'

        frozen_positions=.false.
        frozen_velocities=.false.
        thawed_positions=.false.
        thawed_velocities=.false.
        recip_apply=.false.

        if(.not.tasks%run_Time_Dependent) then
            frozen_positions=.true. 
            frozen_velocities=.true.
        endif

        if(parallel%myid.eq.0) then
            open(newunit=u, file="input", status="old")
            read(u,nml=element_list)
            if(any(thawed_positions.and.frozen_positions)) then
                print *, "all positions cannot be frozen and thawed for any element, &
                & set either thawed_positions or frozen_positions or both to false"
                stop
            endif
            if(any(thawed_velocities.and.frozen_velocities)) then
                print *, "all velocities cannot be frozen and thawed for any element, &
                & set either thawed_velocities or thawed_positions or both to false"
                stop
            endif
            elements(:)%n_atoms_of_element=Number_of_atoms_of_element(:)*product(system%cells)
            elements(:)%M=Mass_amu(:)*u2au
            elements(:)%all_frozen_R=frozen_positions(:)
            elements(:)%all_frozen_V=frozen_velocities(:)
            elements(:)%all_thawed_R=thawed_positions(:)
            elements(:)%all_thawed_V=thawed_velocities(:)
            elements(:)%recip_apply=recip_apply

            do i=1, n_elements
                call read_single_element(PP_type(i), Projector_Application(i), PP_file(i), elements(i))
            enddo
            close(u)
        endif


        do i=1, n_elements
            call parallel_task('bcast',elements(i)%n_atoms_of_element, parallel, 'all', root=0)
            call parallel_task('bcast',elements(i)%PP_type, parallel, 'all', root=0)
            call parallel_task('bcast',elements(i)%PP_file, parallel, 'all', root=0)
            call parallel_task('bcast',elements(i)%PA_type, parallel, 'all', root=0)
            call parallel_task('bcast',elements(i)%M, parallel, 'all', root=0)
            call parallel_task('bcast',elements(i)%all_frozen_R, parallel, 'all', root=0)
            call parallel_task('bcast',elements(i)%all_frozen_V, parallel, 'all', root=0)
            call parallel_task('bcast',elements(i)%all_thawed_R, parallel, 'all', root=0)
            call parallel_task('bcast',elements(i)%all_thawed_V, parallel, 'all', root=0)
            call parallel_task('bcast',elements(i)%recip_apply, parallel, 'all', root=0)
        enddo
            
            system%n_atoms=sum(elements(:)%n_atoms_of_element)
            system%mass=sum(elements(:)%M*elements(:)%n_atoms_of_element)

        do i=1, n_elements
            if(elements(i)%PA_type.eq.Real_local) then
                allocate(elements(i)%RL)
                if(parallel%myid.eq.0) then
                    elements(i)%RL%radius=radius(i)
                    elements(i)%RL%OH_factor=OH_factor(i)
                endif
                call parallel_task('bcast',elements(i)%RL%radius, parallel, 'all', root=0)
                call parallel_task('bcast',elements(i)%RL%OH_factor, parallel, 'all', root=0)
            endif
            !PP_files are only seen by reading thread  
        enddo

        if(parallel%myid.eq.0) then
            print *, 'Element Options'
            do i=1, n_elements
                print *, 'Pseudopotential file of element: ', trim(elements(i)%PP_file)
            enddo
            print *, 'Total # of atoms', system%n_atoms
            print *, '# of atoms of each element:', elements(:)%n_atoms_of_element
            print *, 'Mass of each element:', elements(:)%M/u2au
            print *, 'Pseudopotential type of element: ', elements(:)%PP_type
            print *, 'Pseudopotential Application of element: ', elements(:)%PA_type
            print *
            do i=1, n_elements
                if(elements(i)%PA_type.eq.Real_local) then
                    print *, 'Atomic Sphere Radius of element ', i, ': ', elements(i)%RL%radius
                endif
            enddo
        endif

    end subroutine

    subroutine read_single_element(PP_type, Projector_Application, PP_file, element)
        use element_type, only : element_struct, Real_local, Euler_Spline, Reciprocal, local_only, &
             None, HGH, PAW
        type(element_struct), intent(inout) :: element
        character(*) :: PP_file, PP_type, Projector_Application 
        
        if(PP_type.eq.'None') then
            element%PP_type=None
        else  if(PP_type.eq.'HGH') then
            element%PP_type=HGH
        else  if(PP_type.eq.'PAW') then
            element%PP_type=PAW
        else
            print *, 'Error Unrecognized PP_type: ', PP_type
            print *, 'Options are: None, HGH, PAW'
            stop
        endif

        !The Pseduopotential Projector Application Scheme
        if(Projector_Application.eq.'Real_local'.or. &
           Projector_Application.eq.'Real_Local'.or. &
           Projector_Application.eq.'REAL_LOCAL') then
            element%PA_type=Real_local
        else if (Projector_Application.eq.'Euler_Spline') then
            element%PA_type=Euler_Spline
            print *, 'Euler_Spline not implemented'
            stop
        else if (Projector_Application.eq.'Reciprocal') then
            element%PA_type=Reciprocal
        else if (Projector_Application.eq.'None') then
            element%PA_type=local_only
        else
            print *, 'Projector_Application= ', trim(Projector_Application), ' not supported'
            print *, 'Options are: Real_local, Euler_Spline, Reciprocal'
            stop
        endif
        element%PP_file=trim(PP_file)
    end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads the input file to set the atom positions and momenta parameters for the 
!simulaitons including:
! kT (Temperature), Box_Length(:)
!Input: element, n_elements
!Output: None
!Module Variables updated- 
!        system: num_atoms_of_type(:), PP_type(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_input_atoms(atoms, elements, system, grid, parallel)
        use constants, only : bohr2ang
        use atom_type, only : atom_struct
        use parallel_type, only : parallel_struct
        use system_type, only : system_struct
        use element_type, only : element_struct
        use grids_type, only : grid_struct

        type(parallel_struct), intent(in) ::  parallel
        type(system_struct), intent(inout) :: system
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout) :: grid(:)
        logical :: not_frozen


        character(len=11)   :: header
        integer :: u, reason, i,j,k, rx,ry,rz, atoms_cell
        real(dp) :: Rmin, dR, shift, Rj_near(3), cell_length(3)
        logical :: check
        real(dp), allocatable :: R_base(:,:), P_base(:,:)


        if(parallel%myid.eq.0) then
            !Make the element/atom list
            k=0
            do i=1, system%n_elements
                do j=1, elements(i)%n_atoms_of_element
                    k=k+1
                    atoms(k)%element=i
                    atoms(k)%atom_of_element=j
                enddo
            enddo

            atoms_cell= size(atoms)/product(system%cells)
            cell_length(:) = grid(1)%Box_Length(:)/system%cells(:)
            allocate(R_base(atoms_cell,3))
            allocate(P_base(atoms_cell,3))

            open(newunit=u, file="input", status="old")
            do
                read(u,'(a)', IOSTAT=reason) header
                if(reason.gt.0) then
                    print *, 'Something went wrong looking for positions of atoms'
                    stop
                else if (reason.lt.0) then
                    print *, 'No atom positions found, will generate random positions'
                    Rmin=(0.25*(product(cell_length(:))/atoms_cell)**(1.0_dp/3))
                    do i=1, atoms_cell
                        call random_number(R_base(i,:))
                        R_base(i,:)=R_base(i,:)*cell_length(:)
                    enddo
                    check=.true.
                    do while(check)
                        check=.false.
                        do i=1, atoms_cell
                            do j=i+1, atoms_cell
                                Rj_near(:)=R_base(j,:)
                                !compare to the nearest periodic image
                                Rj_near(:)=R_base(j,:) - &
                                    idnint((R_base(j,:)-R_base(i,:))/cell_length(:)) &
                                        *cell_length(:)
                                dR=sqrt(sum((R_base(i,:)-Rj_near(:))**2))
                                if(dR.lt.Rmin) then
                                    check=.true.
                                    call random_number(shift) 
                                    shift=shift*2.0_dp - 1.0_dp !random between -1 and 1
                                    !shift atom j
                                    R_base(j,:)=R_base(j,:)+(1.0_dp+shift)*(R_base(j,:)-R_base(i,:))*Rmin
                                    !if moved outside box, adjust
                                    R_base(j,:)=R_base(j,:)-floor(R_base(j,:)/cell_length(:)) &
                                        *cell_length(:)
                                endif
                            enddo
                        enddo
                    enddo
                    do i=1, atoms_cell
                        R_base(i,:)=R_base(i,:)/cell_length(:)
                    enddo
                    exit
                else
                    if ( &
                        header(1:10).eq.'$positions'.or.&
                        header(1:10).eq.'&positions' &
                        ) exit ! breaking infinite loop
                endif
            enddo
            k=0

            do i=1, atoms_cell
                if(reason.eq.0) then
                    print *, elements(atoms(k+1)%element)%all_frozen_R,  &
                             elements(atoms(k+1)%element)%all_thawed_R, &
                             .not.elements(atoms(k+1)%element)%all_frozen_R.and. &
                        .not.elements(atoms(k+1)%element)%all_thawed_R


                    if(.not.elements(atoms(k+1)%element)%all_frozen_R.and. &
                        .not.elements(atoms(k+1)%element)%all_thawed_R) then
                        read(u,*) R_base(i,1),R_base(i,2),R_base(i,3), not_frozen
                    else
                        read(u,*) R_base(i,1),R_base(i,2),R_base(i,3)
                        if(elements(atoms(k+1)%element)%all_thawed_R) not_frozen=.true.
                        if(elements(atoms(k+1)%element)%all_frozen_R) not_frozen=.false.
                    endif
                endif
                do rx=1,system%cells(1);do ry=1,system%cells(2); do rz=1,system%cells(3)
                    k=k+1
                    atoms(k)%R(:)=R_base(i,:)+[rx-1,ry-1,rz-1]*1.0_dp
                    atoms(k)%R(:)=atoms(k)%R(:)/system%cells(:)
                    if(reason.eq.0) then
                        atoms(k)%frozen_R=.not.not_frozen
                    else
                        atoms(k)%frozen_R=elements(atoms(k)%element)%all_frozen_R
                    endif
                    atoms(k)%R=atoms(k)%R(:)*grid(1)%Box_Length(:)
                    if(elements(atoms(k)%element)%all_frozen_R.and.not_frozen) then
                        atoms(k)%frozen_R=.true.
                        print *, 'You have selected to freeze all atoms of element: ', atoms(i)%element, &
                        'in the &elements namelits, ignoring any thawed selection for atom: ', i, &
                        ' in the &positions'
                    endif
                enddo; enddo; enddo
            enddo
            close(u)
            open(newunit=u, file="input", status="old")
            do
                read(u,'(a)', IOSTAT=reason) header
                if(reason.gt.0) then
                    print *, 'Something went wrong looking for velocities of atoms'
                    stop
                else if (reason.lt.0) then
                    print *, 'No atom velocites found, will set to 0, or randomize if Isokinetic ensemble used (at start of TD)'
                    atoms(:)%frozen_V=elements(atoms(:)%element)%all_frozen_V
                    atoms(:)%P(1)=0.0_dp
                    atoms(:)%P(2)=0.0_dp
                    atoms(:)%P(3)=0.0_dp
                    atoms(:)%F(1)=0.0_dp
                    atoms(:)%F(2)=0.0_dp
                    atoms(:)%F(3)=0.0_dp
                    exit
                else
                    if ( &
                        header(1:11).eq.'$velocities'.or.&
                        header(1:11).eq.'&velocities' &
                        ) exit ! breaking infinite loop
                endif
            enddo
            if(reason.eq.0) then
                k=0
                atoms_cell= size(atoms)/product(system%cells)
                do i=1, atoms_cell
                    if(.not.elements(atoms(k+1)%element)%all_frozen_V.and. &
                       .not.elements(atoms(k+1)%element)%all_thawed_V) then
                        read(u,*) P_base(i,1),P_base(i,2),P_base(i,3), not_frozen
                    else
                        read(u,*) P_base(i,1),P_base(i,2),P_base(i,3)
                        if(elements(atoms(k+1)%element)%all_thawed_V) not_frozen=.true.
                        if(elements(atoms(k+1)%element)%all_frozen_V) not_frozen=.false.
                    endif
                    do rx=1,system%cells(1);do ry=1,system%cells(2); do rz=1,system%cells(3)
                            k=k+1
                            atoms(k)%P(:)=P_base(i,:)
                            atoms(k)%frozen_V=.not.not_frozen
                    enddo; enddo; enddo
                enddo
                k=0
                do i=1, system%n_elements
                    do j=1, elements(i)%n_atoms_of_element
                        k=k+1
                        atoms(k)%P(:)=atoms(k)%P(:)*elements(i)%M
                        atoms(k)%F(:)=0.0_dp
                    enddo
                enddo
            endif
            close(u)
        endif

        do i=1, system%n_atoms
            call parallel_task('bcast',atoms(i)%P(:),parallel,'all', root=0)
            call parallel_task('bcast',atoms(i)%R(:),parallel,'all', root=0)
            call parallel_task('bcast',atoms(i)%element,parallel,'all', root=0)
            call parallel_task('bcast',atoms(i)%atom_of_element,parallel,'all', root=0)
            call parallel_task('bcast',atoms(i)%frozen_R,parallel,'all', root=0)
            call parallel_task('bcast',atoms(i)%frozen_V,parallel,'all', root=0)
        enddo

        atoms(:)%F(1)=0.0_dp; atoms(:)%F(2)=0.0_dp; atoms(:)%F(3)=0.0_dp

        
        if(parallel%myid.eq.0) then
            print *, "Atom Positions (absolute, Angstrom)"
            k=0
            do i=1, system%n_elements
                do j=1, elements(i)%n_atoms_of_element
                    k=k+1
                    print *, atoms(k)%R(1)*bohr2ang, atoms(k)%R(2)*bohr2ang,atoms(k)%R(3)*bohr2ang
                enddo
            enddo
            print *, "Atom Momentum"
            k=0
            do i=1, system%n_elements
                do j=1, elements(i)%n_atoms_of_element
                    k=k+1
                    print *, atoms(k)%P(1), atoms(k)%P(2),atoms(k)%P(3)
                enddo
            enddo
            print *
        endif

    end subroutine

    subroutine read_input_grids(grids, system, parallel)
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use system_type, only : system_struct
        use constants, only : Ha2eV,pi, bohr2ang, density2gcm3

        type(grid_struct), intent(inout) :: grids(:)
        type(parallel_struct), intent(inout) ::  parallel
        type(system_struct), intent(inout) :: system
        real(dp) :: Ecut, Length_angst(3), Ecut_fine, ecutsm
        integer :: Nx, Ny, Nz, u, Nx_fine, Ny_fine, Nz_fine, g
        logical :: pruned_WF_FFT, internal_FFT, use_time_reversal_sym, debug_TRS
        namelist /simulation_grid/ Nx, Ny, Nz, Ecut, Length_angst, Nx_fine, Ny_fine, Nz_fine, Ecut_fine, ecutsm, & 
            pruned_WF_FFT, internal_FFT, use_time_reversal_sym, debug_TRS

        Nx=-1
        Ny=-1
        Nz=-1
        Length_angst(:)=-1
        Nx_fine=-1
        Ny_fine=-1
        Nz_fine=-1
        Ecut_fine=-1.0
        ecutsm=-1.0_dp
        pruned_WF_FFT=.true.
        internal_FFT=.true.
    
        use_time_reversal_sym=.false.
        debug_TRS=.false.

        Ecut=-1.0
        if(parallel%myid.eq.0) then
            open(newunit=u, file="input", status="old")
                read(u,nml=simulation_grid)
            close(u)

            Ecut=Ecut/Ha2eV
            Ecut_fine=Ecut_fine/Ha2eV
            ecutsm=ecutsm/Ha2eV
            if(pruned_WF_FFT) internal_FFT=.true.

            if(all(Length_angst(:).gt.0.0_dp)) then
                if(system%density.gt.0.0_dp) then
                    print *, 'Specifying both Density and Box dimensions, ignoring input density'
                endif
                system%density=system%mass/product(Length_angst(:)/bohr2ang*system%cells(:))
            else
                if(system%density.lt.0.0_dp) then
                print *, 'Neither Density nor Box dimensions specified in input, stop'
                    stop
                endif
                Length_angst(:)=bohr2ang*(system%density/system%mass*product(system%cells))**(-1.0_dp/3.0_dp)
            endif

            grids(:)%ecutsm=ecutsm
            grids(:)%gamma=.false.

            do g=3, system%n_kpoints+2
                grids(g)%Box_Length(:)=Length_angst(:)/bohr2ang*system%cells(:)
                grids(g)%k(1)=system%k_point(1,g-2)
                grids(g)%k(2)=system%k_point(2,g-2)
                grids(g)%k(3)=system%k_point(3,g-2)
                if(internal_FFT .and. all(abs(grids(g)%k).lt.tiny(1.0_dp))) grids(g)%gamma=.true.
            enddo

            grids(2)%Box_Length=grids(3)%Box_Length
            grids(2)%k(:)=0.0_dp
            if(internal_FFT) grids(2)%gamma=.true.

            grids(1)%Box_Length=grids(2)%Box_Length
            grids(1)%k(:)=0.0_dp
            if(internal_FFT) grids(1)%gamma=.true.

            if(.not.use_time_reversal_sym) grids(3:)%gamma=.false.
            if(debug_TRS) grids(1:2)%gamma=.false.

            grids(:)%FFT_type=1
            do g=1, system%n_kpoints+2
                if(pruned_WF_FFT.and.g.gt.2) then
                    grids(g)%FFT_type=21
                else if (internal_FFT) then
                    grids(g)%FFT_type=11
                endif
                if(parallel%nproc_fft_xy.gt.1) then
                    grids(g)%FFT_type=max(grids(g)%FFT_type, 11)
                    grids(g)%FFT_type=grids(g)%FFT_type+1
                endif
            enddo

            do g=2, system%n_kpoints+2
                if(Nx.gt.0) then
                    if(Ny.lt.0) Ny=Nx
                    if(Nz.lt.0) Nz=Nx
                    grids(g)%Ng(1) = Nx*system%cells(1)
                    grids(g)%Ng(2) = Ny*system%cells(2)
                    grids(g)%Ng(3) = Nz*system%cells(3)
                    if(Ecut.lt.0) then
                        grids(g)%Ecut=(minval(grids(g)%Ng)/grids(g)%Box_Length(minloc(grids(g)%Ng,1))*pi)**2/8.0_dp
                    else 
                        grids(g)%Ecut=Ecut
                    endif
                else if(Ecut.gt.0) then 
                    grids(g)%Ecut=Ecut
                    if(Nx.lt.0) then
                        grids(g)%Ng(:)=floor(sqrt(Ecut*2.0_dp)*grids(g)%Box_Length(:)*2.0_dp/pi+0.5)
                    endif
                else   
                    print *, 'Need to specity Nx or Ecut to establish grid size'
                endif
            enddo
            
            if(Nx_fine.lt.0 .and. Ecut_fine.lt.0.0) then
                grids(1)%Ecut  = grids(2)%Ecut
                grids(1)%Ng(1) = grids(2)%Ng(1)*system%cells(1)
                grids(1)%Ng(2) = grids(2)%Ng(2)*system%cells(2)
                grids(1)%Ng(3) = grids(2)%Ng(3)*system%cells(3)
                
            else if(Nx_fine.gt.0) then
                if(Ny_fine.lt.0) Ny_fine=Nx_fine
                if(Nz_fine.lt.0) Nz_fine=Nx_fine
                grids(1)%Ng(1) = Nx_fine*system%cells(1)
                grids(1)%Ng(2) = Ny_fine*system%cells(2)
                grids(1)%Ng(3) = Nz_fine*system%cells(3)
                if(Ecut_fine.lt.0) then
                    grids(1)%Ecut=(minval(grids(1)%Ng)/grids(1)%Box_Length(minloc(grids(1)%Ng(:),1))*pi)**2/8.0_dp
                else 
                    grids(1)%Ecut=Ecut_fine
                endif
            else 
                grids(1)%Ecut=Ecut_fine
                if(Nx_fine.lt.0) then
                    grids(1)%Ng(:)=floor(sqrt(Ecut_fine*2.0_dp)*grids(1)%Box_Length(:)*2.0_dp/pi+0.5)
                endif
            endif
        endif
            call parallel_task('bcast', system%density, parallel, 'all', root=0)
        do g=1, system%n_kpoints+2
            call parallel_task('bcast', grids(g)%Ng(:), parallel, 'all', root=0)
            call parallel_task('bcast', grids(g)%Ecut,  parallel, 'all', root=0)
            call parallel_task('bcast', grids(g)%Box_Length, parallel, 'all', root=0)
            call parallel_task('bcast', grids(g)%k(:), parallel, 'all', root=0)
            call parallel_task('bcast', grids(g)%ecutsm, parallel, 'all', root=0)
            call parallel_task('bcast', grids(g)%FFT_type, parallel, 'all', root=0)
            call parallel_task('bcast', grids(g)%gamma, parallel, 'all', root=0)
        enddo

        if(parallel%myid.eq.0) then
            if(system%density.gt.0.0_dp) print *, 'Density (g/cc):', system%density*density2gcm3
            print *, 'Box Lengths (Angstrom):', grids(1)%Box_Length(:)*bohr2ang
            print *, "Fine Grid: ", grids(1)%Ng(1), grids(1)%Ng(2), grids(1)%Ng(3)
            print *, "Fine Ecut (eV): ", grids(1)%Ecut*Ha2eV
            print *, "k-point:", grids(1)%k(1), grids(1)%k(2), grids(1)%k(3)
            print *, "Den Coarse Grid: ", grids(2)%Ng(1), grids(2)%Ng(2), grids(2)%Ng(3)
            print *, "Ecut (eV): ", grids(2)%Ecut*Ha2eV
            print *, "k-point:", grids(2)%k(1), grids(2)%k(2), grids(2)%k(3)

            do g=3, system%n_kpoints+2
                if(grids(g)%FFT_type.gt.20)  print *, "WF grid will be pruned in recipricoal space"
                print *, "WF Coarse Grids: ", g, grids(g)%Ng(1), grids(g)%Ng(2), grids(g)%Ng(3)
                print *, "Ecut (eV): ", grids(g)%Ecut*Ha2eV
                print *, "k-point:", grids(g)%k(1), grids(g)%k(2), grids(g)%k(3)
            enddo
           
            print *
        endif


    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine reads the input file to set the physical system parameters for the 
!simulaitons including:
! kT (Temperature), density 
!Input: None
!Output: None
!Module Variables updated- 
!        system: density, kT, num_atom_types, num_atoms_of_type, PP_type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_input_xc(xc, system, parallel)
        use xc_type, only : xc_struct
        use system_type, only : system_struct
        use parallel_type, only : parallel_struct
        use xc_interface, only : init_xc
        
        type(parallel_struct), intent(in) ::  parallel
        type(xc_struct), allocatable, intent(inout) :: xc(:)
        type(system_struct), intent(inout) :: system
        character(len=256) :: exchange, correlation, exchange_correlation
        integer :: u
        
        namelist/xc_functional/ exchange_correlation, correlation, exchange
     
        exchange_correlation= 'NONE'
        exchange= 'LDA_X'
        correlation='LDA_C_PZ'

        if(parallel%myid.eq.0) then
            open(newunit=u, file="input", status="old")
            read(u,nml=xc_functional)
            close(u)
            if(exchange_correlation.ne.'NONE') then
                exchange='NONE'
                correlation='NONE'
            endif
        endif
        call parallel_task('bcast',exchange_correlation, parallel, 'all', root=0)
        call parallel_task('bcast',exchange, parallel, 'all', root=0)
        call parallel_task('bcast',correlation, parallel, 'all', root=0)

        
        if(exchange_correlation.ne.'NONE') then
            allocate(xc(1))
            if(exchange_correlation.eq.'LDA_XC_KSDT') then
                call init_xc(xc(1), trim(exchange_correlation), system%Temperature)
            else
                call init_xc(xc(1), trim(exchange_correlation))
            endif
        else
            allocate(xc(2))
            call init_xc(xc(1), trim(correlation))
            call init_xc(xc(2), trim(exchange))
        endif

        if(parallel%myid.eq.0) then
            print *, 'Exchange Correlation Functional:'
            if(exchange_correlation.eq.'NONE') then
                print *, trim(exchange), " ", trim(correlation)
            else
                print *, trim(exchange_correlation)
            endif
            print *
         endif

    end subroutine

    subroutine read_input_tf(tf, tasks, parallel)
        use tf_type, only : tf_struct
        use parallel_type, only : parallel_struct
        use tasks_type, only : tasks_struct

        type(tasks_struct), intent(inout) ::  tasks
        type(parallel_struct), intent(in) ::  parallel
        type(tf_struct), intent(inout) :: tf
        real(dp) :: lambda, gamma, energy_eps, brent_eps
        integer :: max_iter, update_type
        integer :: u, IOstatus
        logical :: dynamic_kedp
        
        namelist/thomas_fermi/ lambda, gamma, energy_eps, brent_eps, max_iter, update_type, dynamic_kedp
     
        lambda=1.0
        gamma=1.0
        energy_eps=1e-12_dp
        brent_eps=1e-3_dp
        max_iter=2000
        update_type=2
        dynamic_kedp=.false.
        
        if(parallel%myid.eq.0) then
            open(newunit=u, file="input", status="old")
            read(u,nml=thomas_fermi, IOSTAT=IOstatus)
                if(IOstatus < 0) then
                    close(u)
                else
                    tf%lambda=lambda
                    tf%gamma= gamma
                    tf%energy_eps=energy_eps
                    tf%brent_eps=brent_eps
                    tf%max_iter=max_iter
                    tf%update_type=update_type
                    tf%dynamic_kedp=dynamic_kedp
                    close(u)
                endif
            if(tf%dynamic_kedp.and.(.not.tasks%calculate_current)) then
                print *, 'Dynamic KEDP requires current calculation, turning "calculate_current" task on'
                tasks%calculate_current=.true.
            endif
        endif
        call parallel_task('bcast',tf%lambda, parallel, 'all', root=0)
        call parallel_task('bcast',tf%gamma, parallel, 'all', root=0)
        call parallel_task('bcast',tf%energy_eps, parallel, 'all', root=0)
        call parallel_task('bcast',tf%brent_eps, parallel, 'all', root=0)
        call parallel_task('bcast',tf%max_iter, parallel, 'all', root=0)
        call parallel_task('bcast',tf%update_type, parallel, 'all', root=0)
        call parallel_task('bcast',tf%dynamic_kedp, parallel, 'all', root=0)
        call parallel_task('bcast',tasks%calculate_current, parallel, 'all', root=0)

        if(parallel%myid.eq.0) then
            print *, 'orbital free parameters'
            print *, 'lambda: ', tf%lambda
            print *, 'gamma: ', tf%gamma
            print *, 'energy_eps: ', tf%energy_eps
            print *, 'brent_eps: ', tf%brent_eps
            print *, 'max_iter: ', tf%max_iter
            print *, 'update_type: ', tf%update_type
            print *, 'dynamic_kedp: ', tf%dynamic_kedp
            print *
        endif

    end subroutine

    subroutine read_input_KG(KG, parallel)
        use KG_type, only : KG_struct, Cheb, SIL
        use parallel_type, only : parallel_struct
        use constants, only : Ha2eV, pi
        
        type(parallel_struct), intent(in) ::  parallel
        type(KG_struct), intent(inout) :: KG
        real(dp) :: delta_eV, omega_max_eV, dw_eV, min_energy_gap_au, err_allowed
        integer :: u, IOstatus, SIL_rank, Nc_Cheby
        logical :: filter_stoc,project_out_J, SIL_or_C

        namelist/Kubo_Greenwood/ delta_eV, omega_max_eV, dw_eV, min_energy_gap_au, &
            SIL_rank, err_allowed, filter_stoc, project_out_J, Nc_Cheby
     
        delta_eV=0.5_dp
        omega_max_eV=100.0_dp
        dw_eV=0.5_dp
        min_energy_gap_au=1.0E-7_dp
        SIL_rank=5
        err_allowed=1.0E-4_dp
        filter_stoc=.false.
        project_out_J=.false.
        if(parallel%myid.eq.0) then
            open(newunit=u, file="input", status="old")
            print *, delta_eV

            read(u,nml=Kubo_Greenwood, IOSTAT=IOstatus)
                if(IOstatus < 0) then
                    print *, 'Warning: No KG namespace found, all default parameters'
                    close(u)
                else
                    close(u)
                endif
                if(filter_stoc) project_out_J=.false.
                KG%delta=delta_eV/Ha2eV
                KG%dw=dw_eV/Ha2eV
                KG%nw=ceiling(omega_max_eV/dw_eV)
                KG%energy_gap_min=min_energy_gap_au
                KG%SIL_rank=SIL_rank
                KG%filter_stoc=filter_stoc
                KG%project_out_J=project_out_J
                KG%err_allowed=err_allowed
        endif
        call parallel_task('bcast',KG%SIL_rank, parallel, 'all', root=0)
        call parallel_task('bcast',KG%err_allowed, parallel, 'all', root=0)
        call parallel_task('bcast',KG%delta, parallel, 'all', root=0)
        call parallel_task('bcast',KG%nw, parallel, 'all', root=0)
        call parallel_task('bcast',KG%dw, parallel, 'all', root=0)
        call parallel_task('bcast',KG%energy_gap_min, parallel, 'all', root=0)
        call parallel_task('bcast',KG%filter_stoc, parallel, 'all', root=0)
        call parallel_task('bcast',KG%project_out_J, parallel, 'all', root=0)

        if(KG%SIL_rank.lt.2) then
             KG%SIL_rank=2
             if(parallel%myid.eq.0) print *, 'Propagator Krylov subspace must be at least 2 (X and AX)'
        endif
        if(parallel%myid.eq.0) then
            print *, 'Kubo Greenwood parameters'
            print *, 'delta: ', KG%delta*Ha2eV
            print *, 'omega_max: ', KG%dw*KG%nw*Ha2eV
            print *, 'energy_gap_min: ', KG%energy_gap_min
            print *
        endif
    end subroutine

    subroutine read_input_DOS(DOS, parallel)
        use Density_of_states_type, only : Density_of_states_struct
        use parallel_type, only : parallel_struct
        use constants, only : Ha2eV, pi
        
        type(parallel_struct), intent(in) ::  parallel
        type(Density_of_states_struct), intent(inout) :: DOS
        real(dp) :: delta_eV, omega_max_eV, dw_eV
        integer :: u, IOstatus
        
        namelist/Density_of_states/ delta_eV, omega_max_eV, dw_eV
     
        delta_eV=0.5_dp
        omega_max_eV=100.0_dp
        dw_eV=0.5_dp
       

        if(parallel%myid.eq.0) then
            open(newunit=u, file="input", status="old")
            print *, delta_eV

            read(u,nml=Density_of_states, IOSTAT=IOstatus)
                if(IOstatus < 0) then
                    print *, 'Warning: No DOS namespace found, all default parameters'
                    close(u)
                else
                    close(u)
                endif
                DOS%a=(delta_eV/Ha2eV)**2
                DOS%dw=dw_eV/Ha2eV
                DOS%nw=ceiling(omega_max_eV/dw_eV)
               
        endif

        call parallel_task('bcast',DOS%a, parallel, 'all', root=0)
        call parallel_task('bcast',DOS%nw, parallel, 'all', root=0)
        call parallel_task('bcast',DOS%dw, parallel, 'all', root=0)

        if(parallel%myid.eq.0) then
            print *, 'Density of States (Gaussian Approximation) parameters'
            print *, 'a: ', DOS%a*Ha2eV
            print *, 'omega_max: ', DOS%dw*DOS%nw*Ha2eV
            print *
        endif
    end subroutine

    subroutine read_input_td(td, system, grids, elements, parallel, tasks)
        use td_type, only : td_struct, TD_BOMD, TD_RT, TD_None, IsoKinetic, IsoEnergy
        use parallel_type, only : parallel_struct
        use system_type, only: system_struct
        use grids_type, only : grid_struct
        use constants, only : pi
        use element_type, only: element_struct
        use tasks_type, only : tasks_struct

        type(parallel_struct), intent(in) ::  parallel
        type(tasks_struct), intent(in) ::  tasks
        type(system_struct), intent(in) ::  system
        type(grid_struct), intent(in) ::  grids(:)
        type(td_struct), intent(inout) :: td
        type(element_struct), intent(in) :: elements(:)
        integer :: u, IOstatus
        character(len=256), allocatable :: BOMD_or_RealTime, Thermostat


        real(dp) :: dt, total_time, time, err_allowed
        real(dp) :: pulse_E0(3)
        real(dp) :: pulse_t_max_field, pulse_w_field, pulse_tw, pulse_t0, pulse_phase_field
        integer  :: nt, SIL_rank, ETRS_steps, pressure_print_mod, energy_print_mod
        namelist/Time_Dependent/ dt, total_time, pulse_E0, pulse_t_max_field, pulse_w_field, &
        pulse_t0, pulse_phase_field, nt, BOMD_or_RealTime, Thermostat, SIL_rank,pulse_tw, err_allowed, ETRS_steps, &
        pressure_print_mod, energy_print_mod

        time=0.0_dp
        dt=-1
        total_time=-1
        nt=-1
        pulse_E0=0.0_dp
        pulse_t_max_field=-1.0_dp
        pulse_w_field=0.0_dp
        pulse_tw=1.0_dp
        pulse_t0=-10*pulse_tw
        pulse_phase_field=0.0_dp
        BOMD_or_RealTime='None'
        Thermostat='IsoEnergy'
        err_allowed=1.0E-4_dp
        SIL_rank=5
        ETRS_steps=0 
        pressure_print_mod=1
        energy_print_mod=1

        if(parallel%myid.eq.0) then
            open(newunit=u, file="input", status="old")
            read(u,nml=Time_Dependent, IOSTAT=IOstatus)
            if(IOstatus < 0) then
                if(tasks%run_Time_Dependent) then
                        print *, 'Need to specify the &Time_Dependent namelist to run time dependent calculation requested in &
                                & task namelist'
                        stop
                endif
                close(u)
            else
                if(BOMD_or_RealTime.eq.'RealTime'.or. &
                    BOMD_or_RealTime.eq.'Real Time'.or. &
                    BOMD_or_RealTime.eq.'realtime'.or. &
                    BOMD_or_RealTime.eq.'real time'.or. &
                    BOMD_or_RealTime.eq.'Real time'.or. &
                    BOMD_or_RealTime.eq.'Real'.or. &
                    BOMD_or_RealTime.eq.'real') then
                        td%type=TD_RT
                else if(BOMD_or_RealTime.eq.'BOMD'.or. &
                    BOMD_or_RealTime.eq.'bomd'.or. &
                    BOMD_or_RealTime.eq.'Bomd') then
                        td%type=TD_BOMD
                else if(BOMD_or_RealTime.eq.'None'.or. &
                        BOMD_or_RealTime.eq.'NONE'.or. &
                        BOMD_or_RealTime.eq.'none') then
                        td%type=TD_None
                else
                    print *, 'Unrecognized TD_type: ', trim(BOMD_or_RealTime)
                    print *, 'Types are "Real Time" for electron dynamics or Ehrenfest'
                    print *, 'Types are "BOMD" for Born-Oppenheimer Molecular Dynamics'
                    stop
                endif
                
            endif


            if((dt.lt.0.and.total_time.lt.0).or. &
                (dt.lt.0.and.nt.lt.0)) then
                print *, 'dt not specified (nor it total time and number of steps'
                if(td%type.eq.TD_BOMD) then
                    print *, 'dt estimated from density, average element mass, and minimum element mass and temperature'
                    dt=(0.75_dp/pi*product(grids(1)%Box_Length(:))/system%n_atoms)**(1.0_dp/3) / &
                                sqrt(3.0_dp*system%temperature/minval(elements(:)%M))
                    dt=dt/25.0_dp
                else if(td%type.eq.TD_RT) then
                    dt=0.25_dp*pi/grids(2)%Ecut
                    print *, 'dt estimated as pi/Ecut (for wf)'
                else
                    dt=0.0_dp
                endif
                print *, 'dt (Atomic units): ', dt
            endif

            td%dt=dt
            td%total_time=total_time
            td%time=time
            td%nt=nt
            td%E0=pulse_E0(:)
            td%t_max_field=pulse_t_max_field
            td%w_field=pulse_w_field
            td%t0=pulse_t0
            td%tw=pulse_tw
            td%SIL_rank=SIL_rank
            td%ETRS_steps=ETRS_steps
            td%phase_field=pulse_phase_field
            td%err_allowed=err_allowed
            td%p_print_mod=pressure_print_mod
            td%e_print_mod=energy_print_mod

            if(td%SIL_rank.lt.2) then
                td%SIL_rank=2
                if(parallel%myid.eq.0) print *, 'Propagator Krylov subspace must be at least 2 (X and AX)'
            endif
            if(Thermostat.eq.'IsoKinetic'.or. &
                Thermostat.eq.'Isokinetic'.or. &
                Thermostat.eq.'isokinetic'.or. &
                Thermostat.eq.'IsoKineti'.or. &
                Thermostat.eq.'Isokineti'.or. &
                Thermostat.eq.'isokineti') then
                td%thermostat=IsoKinetic
            else if(Thermostat.eq.'IsoEnergy'.or. &
                    Thermostat.eq.'Isoenergy'.or. &
                    Thermostat.eq.'isoenergy') then
                td%thermostat=IsoEnergy
            else 
                print *, 'Unrecognized Thermostat, options are IsoKinetic and IsoEnergy'
                stop
            endif
            
            close(u)
            if(pulse_tw.le.tiny(1.0_dp)) then
                print *, "Pulse width needs to be a positive finite value"
                stop
            endif 
            if((dt.ge.0).and.(total_time.ge.0).and.(nt.ge.0)) then
                print *, 'The time step, number of time steps, and total time cannot be simultaneously defined'
                print *, 'recalculating the number of time steps'
                td%nt=ceiling((td%total_time-td%time)/td%dt)
            else if((dt.gt.0).and.(total_time.ge.0)) then
                td%nt=ceiling((td%total_time-td%time)/td%dt)
            else if((nt.ge.0).and.(total_time.ge.0)) then
                td%dt=(td%total_time-td%time)/td%nt
            else if((nt.ge.0).and.(dt.ge.0)) then
                td%total_time=td%time + td%nt*td%dt
            endif
        endif
        call parallel_task('bcast',td%dt, parallel, 'all', root=0)
        call parallel_task('bcast',td%total_time, parallel, 'all', root=0)
        call parallel_task('bcast',td%time, parallel, 'all', root=0)
        call parallel_task('bcast',td%nt, parallel, 'all', root=0)
        call parallel_task('bcast',td%E0, parallel, 'all', root=0)
        call parallel_task('bcast',td%t_max_field, parallel, 'all', root=0)
        call parallel_task('bcast',td%w_field, parallel, 'all', root=0)
        call parallel_task('bcast',td%t0, parallel, 'all', root=0)
        call parallel_task('bcast',td%tw, parallel, 'all', root=0)
        call parallel_task('bcast',td%phase_field, parallel, 'all', root=0)
        call parallel_task('bcast',td%type, parallel, 'all', root=0)
        call parallel_task('bcast',td%SIL_rank, parallel, 'all', root=0)
        call parallel_task('bcast',td%thermostat, parallel, 'all', root=0)
        call parallel_task('bcast',td%err_allowed, parallel, 'all', root=0)
        call parallel_task('bcast',td%ETRS_steps, parallel, 'all', root=0)
        call parallel_task('bcast',td%Nc, parallel, 'all', root=0)
        call parallel_task('bcast',td%dt, parallel, 'all', root=0)
        call parallel_task('bcast',td%p_print_mod, parallel, 'all', root=0)
        call parallel_task('bcast',td%e_print_mod, parallel, 'all', root=0)

        if(parallel%myid.eq.0) then
            print *, 'Time Dependent: '
            print *, 'Number of time steps: ', td%nt
            print *, 'Start and End time: ', td%time, td%total_time
            print *, 'Time Step Size: ', td%dt
            if(td%type.eq.TD_RT)   print *, 'Dynamics are Real Time Electron or Nonadiabatic'
            if(td%type.eq.TD_BOMD) print *, 'Dynamics are Born Oppenheimer MD'
            if(td%type.eq.TD_None) print *, 'There is no Time dependent'
            print *, 'Print energies every', td%e_print_mod, 'th step'
            print *, 'Print pressure every', td%p_print_mod, 'th step'

            if(any(abs(td%E0).gt.0.0_dp)) then
                print *, 'Time Dependent Pulse: '
                print *, 'E0: ', td%E0
                print *, 'Max Pulse Time: ', td%t_max_field
                print *, 'Field frequency (omega): ', td%w_field
                print *, 'Pulse center time: ', td%t0
                print *, 'Pulse width time: ', td%tw
                print *, 'Field Phase shift: ', td%phase_field
                print *
            endif
        endif

    end subroutine

    subroutine read_input_stopping(stopping, system, grids, parallel)
        use constants, only : u2au, bohr2ang
        use stopping_type, only : stopping_struct
        use system_type, only : system_struct

        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use element_type, only : Real_Local

        type(parallel_struct), intent(in) ::  parallel
        type(stopping_struct), intent(inout) :: stopping
        type(system_struct), intent(inout) :: system
        type(grid_struct), intent(inout), target ::  grids(:)

        real(dp) :: Projectile_R(3), Projectile_V(3), Projectile_Mass_amu, Projectile_radius 
        real(dp) :: t0
        character(len=256) :: Projectile_PP_file, Projectile_PP_type, Projectile_Projector_Application 

        logical :: Constant_Velocity , remove_BO
        integer :: u, IOstatus, Projectile_OH_factor
        
        namelist/stopping_power/ Projectile_PP_file, Projectile_PP_type, Projectile_Projector_Application, &
        Projectile_Mass_amu, Projectile_radius, Projectile_R, Projectile_V, &
        Constant_Velocity, remove_BO, t0, Projectile_OH_factor
     
        Projectile_radius=3.5_dp
        Projectile_Projector_Application(:)='None'
        Constant_Velocity=.true.
        remove_BO=.false.
        t0=0.0_dp
        Projectile_OH_factor=2

        if(parallel%myid.eq.0) then
            open(newunit=u, file="input", status="old")
            read(u,nml=stopping_power, IOSTAT=IOstatus)
                if(IOstatus < 0) then
                    close(u)
                    print *, 'Stopping Power Task listed, but no stopping_power inputs found'
                    stop
                else
                    stopping%element(1)%n_atoms_of_element=1
                    stopping%element(1)%M=Projectile_Mass_amu*u2au
                    call read_single_element(Projectile_PP_type, Projectile_Projector_Application, &
                        Projectile_PP_file, stopping%element(1))
                    close(u)
                endif
            stopping%atom(1)%R(:)=Projectile_R(:)*grids(1)%Box_Length(:)/system%cells(:)
            stopping%atom(1)%P(:)=Projectile_V(:)*stopping%element(1)%M
            stopping%atom(1)%F(:)=0.0_dp
            stopping%atom(1)%frozen_R=.false.
            stopping%atom(1)%frozen_V=Constant_Velocity
            stopping%atom(1)%element=system%n_elements+1
            stopping%atom(1)%atom_of_element=1
            stopping%remove_t0=remove_BO
            stopping%t0=t0
        endif

        call parallel_task('bcast',stopping%element(1)%n_atoms_of_element, parallel, 'all', root=0)
        call parallel_task('bcast',stopping%element(1)%PP_type, parallel, 'all', root=0)
        call parallel_task('bcast',stopping%element(1)%PP_file, parallel, 'all', root=0)
        call parallel_task('bcast',stopping%element(1)%PA_type, parallel, 'all', root=0)
        call parallel_task('bcast',stopping%element(1)%M, parallel, 'all', root=0)
        call parallel_task('bcast',stopping%atom(1)%P(:),parallel,'all', root=0)
        call parallel_task('bcast',stopping%atom(1)%R(:),parallel,'all', root=0)
        call parallel_task('bcast',stopping%atom(1)%F(:),parallel,'all', root=0)
        call parallel_task('bcast',stopping%atom(1)%element,parallel,'all', root=0)
        call parallel_task('bcast',stopping%atom(1)%atom_of_element,parallel,'all', root=0)
        call parallel_task('bcast',stopping%atom(1)%frozen_R,parallel,'all', root=0)
        call parallel_task('bcast',stopping%atom(1)%frozen_V,parallel,'all', root=0)
        call parallel_task('bcast',stopping%remove_t0,parallel,'all', root=0)
        call parallel_task('bcast',stopping%t0,parallel,'all', root=0)

        if(stopping%element(1)%PA_type.eq.Real_Local) then
            allocate(stopping%element(1)%RL)
            if(parallel%myid.eq.0) then
                stopping%element(1)%RL%radius=Projectile_radius
                stopping%element(1)%RL%OH_factor=Projectile_OH_factor
            endif
            call parallel_task('bcast',stopping%element(1)%RL%radius, parallel, 'all', root=0)
            call parallel_task('bcast',stopping%element(1)%RL%OH_factor, parallel, 'all', root=0)
        endif

        if(parallel%myid.eq.0) then
            print *, 'Projectile Element Options'
            print *, 'Pseudopotential file of Projectile: ', trim(stopping%element(1)%PP_file)
            print *, 'Mass:', stopping%element(1)%M/u2au
            print *, 'Pseudopotential type of Projectile: ', stopping%element(1)%PP_type
            print *, 'Pseudopotential Application of Projectile: ', stopping%element(1)%PA_type
            print *
            if(stopping%element(1)%PA_type.eq.Real_local) then
                print *, 'Atomic Sphere Radius of Projectile: ', stopping%element(1)%RL%radius
            endif
            print *, "Atom Positions (absolute, Angstrom)"
            print *, stopping%atom(1)%R(1)*bohr2ang, stopping%atom(1)%R(2)*bohr2ang, stopping%atom(1)%R(3)*bohr2ang
            print *, 'Veloctiy is constant?', stopping%atom(1)%frozen_V
            print *, "Atom Momentum"
                    print *, stopping%atom(1)%P(1), stopping%atom(1)%P(2),stopping%atom(1)%P(3)
            print *
        endif

    end subroutine

    subroutine read_pseudopotentials(elements, grid, system, xc, parallel, all_PAW)
        use constants, only: pi
        use constants, only : u2au
        use element_type, only : element_struct, None, HGH, PAW, Real_Local
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use system_type, only : system_struct
        use simulation_type, only : all_PAW_struct

        use linalg, only: eigh
        use xc_type, only : xc_struct
        use m_pawxmlps, only: rdpawpsxml, paw_setup_free, paw_setuploc
        use m_libpaw_libxc_funcs, only: libpaw_libxc_init
        use xc_f03_lib_m
        use m_pawpsp, only : pawpsp_main
        use m_pawtab, only: pawtab_nullify


        type(xc_struct), intent(inout), allocatable :: xc(:)
        type(parallel_struct), intent(in) ::  parallel
        type(element_struct), intent(inout), target :: elements(:)
        type(element_struct),  pointer :: element

        type(grid_struct), intent(inout), target ::  grid(:)
        type(system_struct), intent(inout), target ::  system
        type(all_PAW_struct), intent(inout), target :: all_PAw

        integer :: u
        integer ::  i, j, l, m, k, nnonloc, e, n_elements_total

        n_elements_total=size(elements)
        !Count PAW elements
        all_PAW%N_PAW_elements=0
        do e = 1, n_elements_total
            element=>elements(e)
            if(element%PP_type.eq.PAW) all_PAW%N_PAW_elements=all_PAW%N_PAW_elements+1
        enddo

        !Setup PAW targets & pointers, basic PAW variables
        if(all_PAW%N_PAW_elements.gt.0) then
            
            allocate(all_PAW%list(all_PAW%N_PAW_elements))
            allocate(all_PAW%rad(all_PAW%N_PAW_elements))
            allocate(all_PAW%tab(all_PAW%N_PAW_elements))
            allocate(all_PAW%nattyp(all_PAW%N_PAW_elements))
            allocate(all_PAW%znucl(all_PAW%N_PAW_elements))

            all_PAW%pawxcdev=1
            all_PAW%ixc=0
            if(allocated(xc)) then
                all_PAW%ixc=xc(1)%id
                if(xc(1)%family.eq.XC_FAMILY_LDA) then 
                    all_PAW%xclevel=1
                else
                    all_PAW%xclevel=2
                endif
            else
                all_PAW%xclevel=0
            endif
            if(allocated(xc).and.size(xc).gt.1) then
                all_PAW%ixc=all_PAW%ixc+xc(2)%id*1000
                if(all_PAW%xclevel.eq.1) then
                    if(xc(2)%family.ne.XC_FAMILY_LDA) then
                        all_PAW%xclevel=2
                    endif 
                endif
            endif
            all_PAW%ixc=-all_PAW%ixc
            if(any(xc(:)%family.eq.XC_FAMILY_LDA)) then 
                all_PAW%xclevel=1
            endif
            all_PAW%nspden=system%n_spinor*system%n_spin
            all_PAW%spnorb=0
            if(system%n_spinor.ne.1) all_PAW%spnorb=1

            i=0
            do e = 1, n_elements_total
                element=>elements(e)
                if(element%PP_type.eq.PAW) then
                    allocate(element%PAW)
                    i=i+1
                    all_PAW%list(i)=e
                    element%PAW%index=i
                    call pawtab_nullify(all_PAW%tab(i))
                    element%PAW%tab=>all_PAW%tab(i)
                    element%PAW%rad=>all_PAW%rad(i)
                endif
            enddo
        endif

        !Start reading pseudopotential files
        do e=1, n_elements_total
            element=>elements(e)
            element%n_proj=0
            if(element%PP_type.eq.HGH) then
                if(parallel%myid==0) open(newunit=u, file=trim(element%PP_file), status="old")
                if(.not.allocated(element%HGH)) allocate(element%HGH)
                if(parallel%myid==0) then
                    element%HGH%C(:)=0.0_dp
                    read(u,*) ! Header
                    read(u,*) element%Znuc, element%Zion
                    read(u,*) !only need for NL
                    read(u,*) element%HGH%rloc, element%HGH%nloc
                    backspace u; read(u,*) element%HGH%rloc, element%HGH%nloc, element%HGH%C(1:element%HGH%nloc)
                    element%HGH%epsatm=2.d0*pi*element%HGH%rloc**2*element%Zion &
                    +(2.d0*pi)**(1.5d0)*element%HGH%rloc**3*&
                    & (element%HGH%C(1)+3.d0*element%HGH%C(2)+15.d0*element%HGH%C(3)+105.d0*element%HGH%C(4))
                endif

                call parallel_task('bcast',element%HGH%C, parallel, 'all', root=0)
                call parallel_task('bcast',element%HGH%nloc, parallel, 'all', root=0)
                call parallel_task('bcast',element%HGH%rloc, parallel, 'all', root=0)
                call parallel_task('bcast',element%Zion, parallel, 'all', root=0)
                call parallel_task('bcast',element%Znuc, parallel, 'all', root=0)
                call parallel_task('bcast',element%HGH%epsatm, parallel, 'all', root=0)
                
                if(parallel%myid==0) then
                    print *, 'HGH Local Pseudopotential Parameters'
                    print *, 'Z nuclei:', element%Znuc, ' Z pseudo: ', element%Zion
                    print *, 'rloc: ', element%HGH%rloc
                    print *, 'C: ', element%HGH%C(1:element%HGH%nloc)
                    print *, 'epsatm: ', element%HGH%epsatm
                endif

                element%HGH%n_proj=0
                element%HGH%H=0.0_dp
                element%HGH%K=0.0_dp
                element%HGH%EKB=0.0_dp
                element%HGH%KKB=0.0_dp
                element%HGH%nlh=0.0
                element%HGH%nlk=0.0
                element%HGH%rl=0.0
                if(parallel%myid.eq.0) close(u)
                if(parallel%myid.eq.0) open(newunit=u, file=trim(element%PP_file), status="old")

                if(parallel%myid.eq.0) then
                    print *, 'Read NL'
                    read(u,*) ! Header
                    read(u,*) !Already Read this
                    read(u,*) !Already Read this
                    read(u,*) !Already Read this
                    read(u,*) nnonloc
                    print *, '# NonLoc in ', trim(element%PP_file), ': ', nnonloc
                    do i=1,nnonloc
                            read(u,*) element%HGH%rl(i), element%HGH%nlh(i), element%HGH%H(1,1:element%HGH%nlh(i),i)
                            element%HGH%nlk(i)=element%HGH%nlh(i)
                            print *, '# N_terms in proj ', i ,': ', element%HGH%nlh(i), element%HGH%nlk(i)
                            do j=2,element%HGH%nlh(i)
                                    read(u,*) element%HGH%H(j,j:element%HGH%nlh(i),i) 
                            enddo
                            if(i.gt.1) then
                                    do j=1,element%HGH%nlk(i)
                                            read(u,*) element%HGH%K(j,j:element%HGH%nlk(i),i) 
                                    enddo
                            endif
                    enddo
                    close(u)

                    do l=1,nnonloc;do i=1,3;do j=1,3
                        if(abs(element%HGH%H(j,i,l)).lt.1E-8) element%HGH%H(j,i,l)=0.0_dp 
                        if(abs(element%HGH%K(j,i,l)).lt.1E-8) element%HGH%K(j,i,l)=0.0_dp
                    enddo;enddo;enddo
                    do l=1,nnonloc
                        element%HGH%H(2,1,l)=element%HGH%H(1,2,l)
                        element%HGH%H(3,1,l)=element%HGH%H(1,3,l)
                        element%HGH%H(3,2,l)=element%HGH%H(2,3,l)
                        element%HGH%K(2,1,l)=element%HGH%K(1,2,l)
                        element%HGH%K(3,1,l)=element%HGH%K(1,3,l)
                        element%HGH%K(3,2,l)=element%HGH%K(2,3,l)
                    enddo          
                endif
                call parallel_task('bcast',nnonloc, parallel, 'all', root=0)
                call parallel_task('bcast',element%HGH%rl, parallel, 'all', root=0)
                call parallel_task('bcast',element%HGH%nlh, parallel, 'all', root=0)
                call parallel_task('bcast',element%HGH%nlk, parallel, 'all', root=0)
                call parallel_task('bcast',element%HGH%H, parallel, 'all', root=0)
                call parallel_task('bcast',element%HGH%K, parallel, 'all', root=0)

                do l=1,nnonloc
                    do j=1,element%HGH%nlh(l)
                            if(abs(element%HGH%H(j,j,l)).lt.1E-8) then
                                    element%HGH%nlh(l)=j-1
                                    exit
                            endif
                    enddo
                    do j=1,element%HGH%nlk(l)
                            if(abs(element%HGH%K(j,j,l)).lt.1E-8) then
                                    element%HGH%nlk(l)=j-1
                                    exit
                            endif
                    enddo
                    do j=1,element%HGH%nlh(l);do i=1,element%HGH%nlh(l)
                            if((abs(element%HGH%H(i,j,l)).gt.1E-8).or.(abs(element%HGH%H(i,j,l)).gt.1E-8)) then
                            if(parallel%myid.eq.0)        print *, "l: ", l-1, " , i: ",i, ", j: ",j,&
                                                ", H:", element%HGH%H(i,j,l)
                            else
                                            
                            endif
                    enddo;enddo
                    do j=1,element%HGH%nlk(l);do i=1,element%HGH%nlk(l)
                        if((abs(element%HGH%K(i,j,l)).gt.1E-8).or.(abs(element%HGH%K(i,j,l)).gt.1E-8)) then
                            if(parallel%myid.eq.0) print *, "l: ", l-1, " , i: ",i, ", j: ",j,&
                                                ", K:", element%HGH%K(i,j,l)
                        else
                                        
                        endif
                    enddo;enddo

                    if(system%spin_type.le.1) element%HGH%nlk(l)=0

                    if(element%HGH%nlh(l).gt.0) call eigh(element%HGH%H(1:element%HGH%nlh(l),1:element%HGH%nlh(l),l),&
                        element%HGH%EKB(1:element%HGH%nlh(l),l), element%HGH%U_HE(1:element%HGH%nlh(l),1:element%HGH%nlh(l),l))
                    if(element%HGH%nlk(l).gt.0) call eigh(element%HGH%K(1:element%HGH%nlk(l),1:element%HGH%nlk(l),l),&
                        element%HGH%KKB(1:element%HGH%nlk(l),l), element%HGH%U_HK(1:element%HGH%nlk(l),1:element%HGH%nlk(l),l))

                    element%HGH%n_proj=element%HGH%n_proj + (element%HGH%nlh(l) + element%HGH%nlk(l))*(2*(l-1)+1) 
                    !if(parallel%myid.eq.0) print *, l-1, element%HGH%n_proj
                enddo

                if(element%HGH%n_proj.gt.0) then
                    allocate(element%HGH%l(element%HGH%n_proj))
                    allocate(element%HGH%m(element%HGH%n_proj))
                    allocate(element%HGH%i(element%HGH%n_proj))
                    allocate(element%HGH%h_or_k(element%HGH%n_proj))

                    k=0
                    do l=1,nnonloc; do i= 1, element%HGH%nlh(l); do m=-(l-1),(l-1)
                                k=k+1                    
                                element%HGH%l(k)=l-1
                                element%HGH%m(k)=m
                                element%HGH%i(k)=i
                                element%HGH%h_or_k(k)='h'
                    enddo; enddo; enddo
                    if(system%spin_type.gt.1) then
                        do l=1,nnonloc;do i= 1, element%HGH%nlk(l); do m=-(l-1),(l-1)
                                k=k+1                    
                                element%HGH%l(k)=l-1
                                element%HGH%m(k)=m
                                element%HGH%i(k)=i
                                element%HGH%h_or_k(k)='k'
                        enddo;enddo; enddo
                    endif
                endif
                element%n_proj= element%n_proj + element%HGH%n_proj   

            else if(element%PP_type.eq.PAW) then
                element%PAW%filepsp=trim(element%PP_file)
                element%PAW%ixc=all_PAW%ixc
                element%PAW%xclevel=all_PAW%xclevel
                element%PAW%nspden=all_PAW%nspden
                element%PAW%spnorb=all_PAW%spnorb
                element%PAW%pawxcdev=all_PAW%pawxcdev
                element%PAW%usewvl=all_PAW%usewvl
                element%PAW%usexcnhat=all_PAW%usexcnhat
                element%PAW%xc_denpos=all_PAW%xc_denpos

                if (element%PAW%xml_file) then
                    !LibPAW internally calls paw_setuploc so it needs to be called here
                    element%PAW%setup=>paw_setuploc
                    call rdpawpsxml(element%PAW%filepsp, element%PAW%setup)
                ! call rdpawpsxml(element%PAW%filepsp, paw_setuploc)
                    element%PAW%lnmax=element%PAW%setup%valence_states%nval
                    element%Zion=element%PAW%setup%atom%zval
                    element%Znuc=element%PAW%setup%atom%znucl
                else
                    print *, 'Not setup to read non XML PAW files'
                    stop
                end if
                call libpaw_libxc_init(element%PAW%ixc,element%PAW%nspden)

                element%PAW%mqgrid_ff = 3001
                element%PAW%mqgrid_vl = 3001
                allocate(element%PAW%qgrid_ff(element%PAW%mqgrid_ff))
                allocate(element%PAW%qgrid_vl(element%PAW%mqgrid_vl))
                do i=1,element%PAW%mqgrid_vl
                    element%PAW%qgrid_vl(i)=(i-1)*(1.2_dp/pi*sqrt(grid(1)%Ecut*2)/(element%PAW%mqgrid_vl-1))
                enddo
                do i=1,element%PAW%mqgrid_ff
                    element%PAW%qgrid_ff(i)=(i-1)*(1.2_dp/pi*sqrt(grid(1)%Ecut*2)/(element%PAW%mqgrid_ff-1))
                enddo
                allocate(element%PAW%ffspl(element%PAW%mqgrid_ff,2,element%PAW%lnmax))
                allocate(element%PAW%vlspl(element%PAW%mqgrid_vl,2))
                !compute and store real space projectors (has_tproj_in=2)

                !pawpsp_main has to read paw_setuploc as the psxml, this seems like a bug of the library
                !but it's ok to just use paw_setuploc (a global variable! )
                call pawpsp_main(pawrad=element%PAW%rad,pawtab=element%PAW%tab,&
                & filpsp=element%PAW%filepsp,usewvl=element%PAW%usewvl,icoulomb=element%PAW%icoulomb,ixc=element%PAW%ixc, &
                xclevel=element%PAW%xclevel,pawxcdev=element%PAW%pawxcdev,usexcnhat=element%PAW%usexcnhat,&
                & qgrid_ff=element%PAW%qgrid_ff,qgrid_vl=element%PAW%qgrid_vl,ffspl=element%PAW%ffspl,vlspl=element%PAW%vlspl, &
                epsatm=element%PAW%epsatm,xcccrc=element%PAW%xcccrc,zionpsp=element%Zion,znuclpsp=element%Znuc,&
                psxml=element%PAW%setup,comm_mpi=parallel%comm_all,xc_denpos=element%PAW%xc_denpos, has_tproj_in=2)

                element%PAW%n_proj=element%PAW%tab%lmn_size
                element%n_proj=element%n_proj+element%PAW%n_proj
                element%PAW%radius=element%PAW%tab%rpaw
                if(element%PA_type.eq.Real_local) then
                    element%RL%radius=element%PAW%tab%rpaw
                endif

                      
                !call paw_setup_free(element%PAW%setup)
                element%PAW%setup=>Null()
                call paw_setup_free(paw_setuploc)
            else if(element%PP_type.eq.None) then
                    if(parallel%myid==0) then
                        if(parallel%myid.eq.0) open(newunit=u, file=trim(element%PP_file), status="old")
                        read(u,*) ! Header
                        read(u,*) element%Znuc
                        element%Zion=element%Znuc
                    endif
                    
                    call parallel_task('bcast',element%Zion, parallel, 'all', root=0)
                    call parallel_task('bcast',element%Znuc, parallel, 'all', root=0)

                call parallel_task('bcast',element%Zion, parallel, 'all', root=0)
                call parallel_task('bcast',element%Znuc, parallel, 'all', root=0)

                if(parallel%myid==0) then
                    print *, 'Coulombic Potential Used'
                    print *, 'Z:', element%Zion
                endif           
            endif
        enddo
    end subroutine

end module
