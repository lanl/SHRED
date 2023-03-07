module initialization_finalization
    use types, only : dp

    implicit none

    public initialize_sim, finalize_sim, initialize_atoms, initialize_elements, init_den_pot_orb, &
     initialize_atoms_temperature, initialize_time_dependent_0, finalize_all_paw

contains

    subroutine initialize_sim(sim)
        use, intrinsic :: iso_c_binding

        use simulation_type, only: simulation_struct
        use constants, only : Ha2eV,density2gcm3,bohr2ang, pi
        use parallel_mod, only: initialize_parallel
        use allocate_memory, only: allocate_required_memory_atoms_and_elements, allocate_required_memory_density_and_orbitals
        use grids_mod, only: setup_grids, map_between_distributed_grids, setup_grids_new
        use fft, only: real_to_recip
        use inputs
        use Non_Local_ion, only: pseudopotential_update
        use Ewald_mod, only: Ewald_Calc, init_ewald, free_ewald
        use odp_type, only : deterministic, stochastic
        use Density, only : calc_core_Charge_Density

        type(simulation_struct), intent(inout), target :: sim
        integer :: i
        !Random number generator Seed                                 
        integer,allocatable :: seed_array(:)                   
        integer :: seedsize
        logical :: restart_td
        
        call cpu_time(sim%wall_time_start) !start the clock
        !==================================================================
        call initialize_parallel(sim%parallel)
        if(sim%parallel%myid.eq.0) &
            open(newunit=sim%main_output,  file="main.log", status="unknown", recl=1024)  

        call random_init(.false.,.false.)
        !call random_seed()
        !Get all processes their own RNG
        call random_seed(size=seedsize)
        allocate(seed_array(seedsize))
        call random_seed(get=seed_array)
        !seed_array=1 !set for debug
        if(sim%parallel%myid.ne.83) then
            seed_array(:)=abs( mod((seed_array(:)*181)*((sim%parallel%myid-83)*359), 104729) )
         else
            seed_array(:)=abs( mod((seed_array(:)*193)*((sim%parallel%myid-97)*389), 104729) )
         endif
         call random_seed(put=seed_array)
        !==================================================================
        !Setup the Calculations
        call read_input_tasks(sim%tasks, sim%parallel)
        sim%stopping%do_stopping=sim%tasks%run_Stopping_Power

        !==================================================================

        !Read the rest of input file
        call read_input_system(sim%system, sim%parallel)
        sim%coarse_density%n_s=sim%system%n_spinor*sim%system%n_spin
        sim%density%n_s=sim%system%n_spinor*sim%system%n_spin

        if(sim%system%orbital_free) then
            if(.not.sim%tasks%density_from_file .and. &
               .not.sim%tasks%read_orbitals_from_file) sim%tasks%initial_orbital_free=.true.
                if(sim%tasks%run_Kubo_Greenwood.and.sim%parallel%myid.eq.0) &
                 print *, 'No Kubo Greenwood for Orbital Free'
                sim%tasks%run_Kubo_Greenwood=.false.
        endif
        !==================================================================
        allocate(sim%elements(sim%system%n_elements))
        call read_input_elements(sim%elements,sim%system%n_elements,sim%system, sim%parallel, sim%tasks)            
        !==================================================================
        call read_input_tf(sim%tf, sim%tasks, sim%parallel)
        !==================================================================
        call read_input_xc(sim%xc, sim%system, sim%parallel)
        !==================================================================
        allocate(sim%grids(sim%system%n_kpoints+2))
        call read_input_grids(sim%grids, sim%system, sim%parallel) 
        if(sim%tasks%run_Kubo_Greenwood.or.sim%tasks%run_Stopping_Power) sim%grids(3:)%gamma=.false.        
        do i=1,size(sim%grids)
            sim%grids(i)%work%grid=i
            if(sim%grids(i)%FFT_type.eq.1)  call setup_grids(sim%grids(i), sim%parallel)
            if(sim%grids(i)%FFT_type.eq.11) call setup_grids_new(sim%grids(i), sim%parallel, .false.)
            if(sim%grids(i)%FFT_type.eq.21) call setup_grids_new(sim%grids(i), sim%parallel, .true.)
        enddo
        call map_between_distributed_grids(sim%grids(2), sim%grids(1), sim%parallel, sim%coarse_to_fine)
        call map_between_distributed_grids(sim%grids(1), sim%grids(2), sim%parallel, sim%fine_to_coarse)
        !==================================================================

        !All the ks_orbitals
        sim%system%deterministic_start=1
        sim%system%deterministic_end=sim%system%deterministic_start -1 &
                                            + sim%system%n_deterministic/sim%parallel%n_band_groups
        !The lowest eigenvector KS orbitals
        sim%system%deterministic_lowest_start=1
        sim%system%deterministic_lowest_end=sim%system%n_orbitals_deterministic_lowest/sim%parallel%n_band_groups

        !All the stochastic orbitals
        sim%system%stochastic_start=sim%system%deterministic_end + 1
        sim%system%stochastic_end=sim%system%stochastic_start-1 &
                                            + sim%system%n_stochastic/sim%parallel%n_band_groups
        !The full filter stochastic orbitals
        sim%system%stochastic_full_start=sim%system%deterministic_end + 1
        sim%system%stochastic_full_end=sim%system%stochastic_full_start-1 &
                                            + sim%system%n_orbitals_stochastic_full/sim%parallel%n_band_groups
        
        !orbitals just for throwing into LOBPCG to get the lowest eigenvalue
        sim%system%find_min_start=sim%system%n_orbitals_local+1
        if(sim%parallel%my_sub.eq.0) then
            sim%system%find_min_end=sim%system%n_orbitals_local+1
        else
            sim%system%find_min_end=sim%system%n_orbitals_local
        endif   
    

        if(sim%parallel%myid.eq.0) print *, 'Allocate Memory for tasks: orbitals and densities'
        call allocate_required_memory_density_and_orbitals(sim)

        allocate(sim%atoms(sim%system%n_atoms))
        call read_input_atoms(sim%atoms, sim%elements, sim%system, sim%grids, sim%parallel)

        allocate(sim%system%nelec(sim%system%n_spin))
        sim%system%nelec=0.0_dp
        call read_pseudopotentials(sim%elements,sim%grids,sim%system, sim%xc, sim%parallel, sim%all_PAW)
        do i=1, sim%system%n_elements
                sim%system%nelec(:)=sim%system%nelec(:)+sim%elements(i)%n_atoms_of_element*sim%elements(i)%Zion
        enddo
        if(sim%system%spin_type.gt.1) sim%system%nelec = sim%system%nelec*0.5_dp
        if(sim%parallel%myid.eq.0) print *, 'Total Number of Valence Electrons: ', sim%system%nelec
        !==================================================================
        call initialize_elements(sim%elements,sim%parallel, sim%grids, sim%all_PAW)
        !==================================================================
        call initialize_atoms(sim%atoms, sim%elements, sim%parallel)
        call init_ewald(sim%ewald, sim%grids(1)%box_length, sim%atoms)
        !==================================================================
        call read_input_numerics(sim%numerics, sim%tasks, sim%parallel, sim%elements, sim%system)

        if(sim%parallel%myid.eq.0) print *, 'initialize all_PAW'
        call initialize_PAW(sim%atoms, sim%elements, sim%parallel, sim%system, sim%grids, sim%all_PAW, sim%numerics)
        sim%atoms(:)%update_me=.true.
        sim%atoms(:)%update_me_force=.true.
        if(sim%parallel%myid.eq.0) print *, 'allocate element & atomic functions'
        call allocate_required_memory_atoms_and_elements(sim)

        if(sim%parallel%myid.eq.0) print *, 'Initial Pseudopotenital Setup'
        call pseudopotential_update(sim%atoms, sim%elements, sim%grids(2), sim%grids(1), sim%parallel, sim%all_PAW)
        if(sim%parallel%myid.eq.0) print *, 'Pseudopotenital Setup Complete'
        call calc_core_Charge_Density(sim%all_PAW, sim%grids, sim%atoms, sim%elements, sim%parallel)

        if(sim%parallel%myid.eq.0) print *, 'Ground State initialization Complete'

    end subroutine

    subroutine initialize_time_dependent_0(sim)
        use simulation_type, only: simulation_struct
        use allocate_memory, only: allocate_required_memory_TD
        use inputs
        use td_type, only : TD_RT
        use constants, only : pi
        type(simulation_struct), intent(inout), target :: sim

        !==================================================================
        if(sim%tasks%run_Time_Dependent) call read_input_td(sim%td,sim%system, sim%grids, sim%elements, sim%parallel, sim%tasks)
        if(sim%td%type.eq.TD_RT.and.sim%all_PAW%N_PAW_atoms.gt.0) then
            print *, 'Real Time TD-DFT not compatible with PAW pseudopotentials, use HGH only'
        endif
        !==================================================================
        if(sim%tasks%run_Stopping_Power) then
            call read_input_stopping(sim%stopping, sim%system, sim%grids, sim%parallel)

            if(abs(sim%td%dt-pi/sim%grids(2)%Ecut).lt.tiny(1.0_dp).and. &
                (minval(0.1*sim%grids(1)%dR(:))/sqrt(sum(sim%stopping%atom(1)%P(:)**2))*sim%stopping%element(1)%M)  &
                    .lt. sim%td%dt) then  
                    if(sim%parallel%myid.eq.0) &
                         print *, 'Default electron dynamics time-step detected ', sim%td%dt
                        sim%td%dt=minval(0.1*sim%grids(1)%dR(:))/sqrt(sum(sim%stopping%atom(1)%P(:)**2))*sim%stopping%element(1)%M
                    if(sim%parallel%myid.eq.0) &
                         print *, 'Replacing with Default Stopping power time step &
                         & (0.1 * minimum (grid dR)/projectile speed): ' &
                        , sim%td%dt
            endif
        endif
        !==================================================================

        sim%atoms(:)%update_me=.true.
        sim%atoms(:)%update_me_force=.not.sim%atoms(:)%frozen_V
        
        if(sim%system%orbital_free) then
            if((sim%td%type.eq.TD_RT).or.(abs(sim%tf%lambda).gt.tiny(1.0_dp))) then
                sim%system%n_orbitals_deterministic_lowest=1
                sim%system%n_deterministic=1
                sim%system%n_orbitals=1
                sim%system%n_orbitals_local=1

                  !All the ks_orbitals
                    sim%system%deterministic_start=1
                    sim%system%deterministic_end=1
        
            endif
            if(.not.sim%tasks%density_from_file .and. &
               .not.sim%tasks%read_orbitals_from_file) sim%tasks%initial_orbital_free=.true.
        else if(sim%tf%dynamic_kedp) then
            sim%tf%dynamic_kedp=.false.
            if(sim%parallel%myid.eq.0) print *, 'KEDP is only for Orbital_Free'
        endif

        if(sim%parallel%myid.eq.0) print *, 'Allocate Memory for  time dependent tasks'
        call allocate_required_memory_TD(sim)

    end subroutine

    subroutine initialize_Kubo_Greenwood(sim)
        use simulation_type, only: simulation_struct
        use inputs

        type(simulation_struct), intent(inout), target :: sim

        if(sim%tasks%run_Kubo_Greenwood) call read_input_KG(sim%KG, sim%parallel)

        if(sim%system%orbital_free) then
            if(sim%tasks%run_Kubo_Greenwood.and.sim%parallel%myid.eq.0) &
             print *, 'No Kubo Greenwood for Orbital Free'
            sim%tasks%run_Kubo_Greenwood=.false.
        endif

    end subroutine

    subroutine initialize_DOS(sim)
        use simulation_type, only: simulation_struct
        use inputs

        type(simulation_struct), intent(inout), target :: sim

        if(sim%tasks%calc_DOS) call read_input_DOS(sim%DOS, sim%parallel)

        if(sim%system%orbital_free) then
            if(sim%tasks%calc_DOS.and.sim%parallel%myid.eq.0) &
             print *, 'No DOS for Orbital Free'
            sim%tasks%calc_DOS=.false.
        endif

    end subroutine

    subroutine finalize_sim(sim)
        use simulation_type, only: simulation_struct
        use parallel_mod, only: finalize_parallel
        type(simulation_struct), intent(inout), target :: sim
        
        call finalize_parallel(sim%parallel)
    end subroutine

    subroutine initialize_elements(elements, parallel, grids, all_PAW)
        use types, only : dp
        use constants, only : pi
        use element_type, only : element_struct, Real_Local, HGH, PAW, None
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use simulation_type, only : all_PAW_struct
        use m_paw_init, only : pawinit
        use m_pawtab, only: pawtab_print
        use splines, only: spline3ders_unsorted, spline3_unsorted
        use m_pawrad, only : simp_gen, nderiv_gen
        use linalg, only: eigh
        use parallel_mod, only: parallel_wait
        use m_paw_onsite, only : pawnabla_init, paw_x_nabla_init,paw_dipole_init
        use m_pawrad, only: pawrad_ifromr
        use constants, only : FMT_real, FMT_int, FMT_signed_int
        use grids_mod, only : allocate_local_fields_G

        type(parallel_struct), intent(in) :: parallel
        type(element_struct), intent(inout), target :: elements(:)
        type(element_struct), pointer :: element
        type(all_PAW_struct), intent(inout) :: all_PAw

        type(grid_struct), intent(inout), target :: grids(:)
        real(dp), pointer :: w2(:,:,:), C(:), rloc, q2vq(:), dq2vq(:), d2q2vq(:), &
             g2vec(:), tvale_local(:)
        integer :: i, j, k, e, n_elements_total, mpsang, n_lm, lm, imax
        real(dp), target, allocatable :: d2Ven0G(:,:,:), d2Den0G(:,:,:)
        real(dp), allocatable :: L_ij(:,:) , U_ij(:,:), lam_i(:), Q_ij(:,:) ,  sij_lm(:,:), tzeta(:,:), tproj_lm(:,:), &
                                 ffspl_lm(:,:,:), ffzeta(:,:,:)
        integer, allocatable ::  lmn_of_lm(:)
        real(dp),allocatable :: dphi(:,:),dtphi(:,:)
        integer :: u
        real(dp) :: tmp

        n_elements_total=size(elements)
    !Initalize all PAW elements PAW objects
        if(all_PAW%N_PAW_elements.gt.0) then
            mpsang=0
            do i=1,size(all_PAW%tab(:))
                mpsang=max(mpsang, (maxval(all_PAW%tab(i)%orbitals)+1))
            enddo
            !Options just taken as ABINIT defaults
            call pawinit(1.0_dp, 1, 4.4_dp/pi**2, 0.0_dp,10,10,mpsang,&
            &                  13,0,12,all_PAW%ang,all_PAW%rad(:), &
                                all_PAW%spnorb,all_PAW%tab(:),all_PAW%pawxcdev,all_PAW%ixc,all_PAW%usepotzero)
            if(parallel%myid.eq.0) call pawtab_print(all_PAW%tab(:), unit=6)
            call paw_dipole_init(mpsang,all_PAW%N_PAW_elements,all_PAW%rad(:),all_PAW%tab(:))
            call pawnabla_init(mpsang,all_PAW%N_PAW_elements,all_PAW%rad(:),all_PAW%tab(:))
            call paw_x_nabla_init(mpsang,all_PAW%N_PAW_elements,all_PAW%rad(:),all_PAW%tab(:))

            if(.true..and.parallel%myid.eq.0) then
                imax=pawrad_ifromr(all_PAW%rad(1),all_PAW%tab(1)%rpaw)
                print *, all_PAW%rad(1)%mesh_size, all_PAW%tab(1)%partialwave_mesh_size,  &
                        size(all_PAW%tab(1)%phi(:,1)), size(all_PAW%tab(1)%tphi(:,1)), &
                         size(all_PAW%tab(1)%tproj(:,1)), imax
                         call simp_gen(tmp, all_PAW%tab(1)%tproj(:imax,1)*all_PAW%tab(1)%tphi(:imax,1),all_PAW%rad(1))
                print *, 'simpgen norm: ', tmp 
                do e = 1, size(all_PAW%tab)
                    allocate( dphi(all_PAW%rad(e)%mesh_size,size(all_PAW%tab(e)%phi,2)))
                    allocate(dtphi(all_PAW%rad(e)%mesh_size,size(all_PAW%tab(e)%tphi,2)))
                    do i=1,size(all_PAW%tab(e)%phi,2)
                        call nderiv_gen(dphi(:,i),all_PAW%tab(e)%phi(:,i),all_PAW%rad(e))
                    enddo
                    do i=1,size(all_PAW%tab(e)%phi,2)
                        call nderiv_gen(dtphi(:,i),all_PAW%tab(e)%tphi(:,i),all_PAW%rad(e))
                    enddo
                    deallocate(dphi)
                    deallocate(dtphi)
                    close(u)
                enddo
                !do i=1, all_PAW%tab(1)%mqgrid
                !    print *, elements(1)%PAW%qgrid_vl(i), all_PAW%tab(1)%tcorespl(i,1)
                !enddo

               ! print *, all_PAW%tab(1)%core_mesh_size, size(all_PAW%tab(1)%coredens(:)), size(all_PAW%tab(1)%tcoredens(:,1)), &
               !          all_PAW%tab(1)%rcore
               ! do ir=1, all_PAW%tab(1)%core_mesh_size
               !     print *, all_PAW%rad(1)%rad(ir),  all_PAW%tab(1)%coredens(ir), all_PAW%tab(1)%tcoredens(ir,1)
               ! enddo
            endif


            do e = 1, n_elements_total
                element=>elements(e)
                if(element%PP_type.ne.PAW) cycle
                    imax=pawrad_ifromr(element%PAW%rad,element%PAW%tab%rpaw)

                    if(parallel%myid.eq.0) then
                        print *, 'projector index'
                print *, '        lmn', '           l           ','m           ','n           ','lm          ','ln          ','spin'
                        do i=1, element%PAW%tab%lmn_size
                            write(6,FMT_signed_int) i, element%PAW%tab%indlmn(:,i)
                        enddo
                    endif
                    if(parallel%myid.eq.0) then
                        print *, 'projector^2 index'
                print *, '        klmn', '         klm         ','kln         ','|i_l-j_l|  ', 'i_l+j_l    ', &
                        'i_lm        ','j_ln        ', 'i_lmn      ', 'j_lmn        '
                        do i=1, element%PAW%tab%lmn2_size
                            write(6,FMT_signed_int) i, element%PAW%tab%indklmn(:,i)
                        enddo
                    endif
                    allocate(element%PAW%sij(element%PAW%tab%lmn_size,element%PAW%tab%lmn_size))
                    allocate(element%PAW%eijkl(element%PAW%tab%lmn_size,element%PAW%tab%lmn_size, &
                                            element%PAW%tab%lmn_size,element%PAW%tab%lmn_size))
                    allocate(element%PAW%dij0(element%PAW%tab%lmn_size,element%PAW%tab%lmn_size))
                    allocate(element%PAW%tproj_RO(imax,element%PAW%tab%lmn_size))
                    allocate(element%PAW%ffspl_Ort(size(element%PAW%ffspl,1),2,element%PAW%tab%lmn_size))

                    allocate(element%PAW%o_i(element%PAW%tab%lmn_size))

                    if(parallel%myid.eq.0) print *, 'sij, i/j=> lmn, j<=i'
                    k=0
                    do j=1, element%PAW%tab%lmn_size
                        do i=1, j
                            k=k+1
                            if(parallel%myid.eq.0) write(6,'(E15.5)',advance="no") element%PAW%tab%sij(k)
                            element%PAW%sij(i,j)=element%PAW%tab%sij(k)
                            element%PAW%sij(j,i)=element%PAW%tab%sij(k)
                        enddo
                        if(parallel%myid.eq.0) write(6,*)
                    enddo
                    if(parallel%myid.eq.0) print *, 'nabla_ij.nabla_ij, i/j=> lmn, j<=i'
                    k=0
                    do j=1, element%PAW%tab%lmn_size
                        do i=1, element%PAW%tab%lmn_size
                            if(parallel%myid.eq.0) write(6,'(E15.5)',advance="no") sum(element%PAW%tab%nabla_ij(:,i,j)**2)
                        enddo
                        if(parallel%myid.eq.0) write(6,*)
                    enddo
                    element%PAW%eijkl=-666
                   ! if(parallel%myid.eq.0) print *, 'eij,kl, ij/kl=> lmn2, kl<=ij, j<=i, l<=k'
                    do j=1, element%PAW%tab%lmn2_size
                        do i=1, j
                    !        if(parallel%myid.eq.0) write(6,'(E15.5)',advance="no") element%PAW%tab%eijkl(i,j)
                            element%PAW%eijkl(element%PAW%tab%indklmn(7,i), &
                                            element%PAW%tab%indklmn(8,i), &
                                            element%PAW%tab%indklmn(7,j), &
                                            element%PAW%tab%indklmn(8,j))  = element%PAW%tab%eijkl(i,j)
                            element%PAW%eijkl(element%PAW%tab%indklmn(8,i), &
                                            element%PAW%tab%indklmn(7,i), &
                                            element%PAW%tab%indklmn(7,j), &
                                            element%PAW%tab%indklmn(8,j))  = element%PAW%tab%eijkl(i,j)
                            element%PAW%eijkl(element%PAW%tab%indklmn(7,i), &
                                            element%PAW%tab%indklmn(8,i), &
                                            element%PAW%tab%indklmn(8,j), &
                                            element%PAW%tab%indklmn(7,j))  = element%PAW%tab%eijkl(i,j)
                            element%PAW%eijkl(element%PAW%tab%indklmn(8,i), &
                                            element%PAW%tab%indklmn(7,i), &
                                            element%PAW%tab%indklmn(8,j), &
                                            element%PAW%tab%indklmn(7,j))  = element%PAW%tab%eijkl(i,j)

                            element%PAW%eijkl(element%PAW%tab%indklmn(7,j), &
                                            element%PAW%tab%indklmn(8,j), &
                                            element%PAW%tab%indklmn(7,i), &
                                            element%PAW%tab%indklmn(8,i))  = element%PAW%tab%eijkl(i,j)
                            element%PAW%eijkl(element%PAW%tab%indklmn(8,j), &
                                            element%PAW%tab%indklmn(7,j), &
                                            element%PAW%tab%indklmn(7,i), &
                                            element%PAW%tab%indklmn(8,i))  = element%PAW%tab%eijkl(i,j)
                            element%PAW%eijkl(element%PAW%tab%indklmn(7,j), &
                                            element%PAW%tab%indklmn(8,j), &
                                            element%PAW%tab%indklmn(8,i), &
                                            element%PAW%tab%indklmn(7,i))  = element%PAW%tab%eijkl(i,j)
                            element%PAW%eijkl(element%PAW%tab%indklmn(8,j), &
                                            element%PAW%tab%indklmn(7,j), &
                                            element%PAW%tab%indklmn(8,i), &
                                            element%PAW%tab%indklmn(7,i))  = element%PAW%tab%eijkl(i,j)
                        enddo
                    !    if(parallel%myid.eq.0) write(6,*)
                    enddo

                    if(parallel%myid.eq.0) print *, 'dij0, i/j=> lmn, j<=i'
                    k=0
                    do j=1, element%PAW%tab%lmn_size
                        do i=1, j
                            k=k+1
                            if(parallel%myid.eq.0) write(6,'(E15.5)',advance="no") element%PAW%tab%dij0(k)
                            element%PAW%dij0(i,j)=element%PAW%tab%dij0(k)
                            element%PAW%dij0(j,i)=element%PAW%tab%dij0(k)
                        enddo
                        if(parallel%myid.eq.0) write(6,*)
                    enddo

                    !Rotate Blocks of lm
                    do lm=1, maxval(element%PAW%tab%indlmn(4,:))
                        !print *, 'Rotate tproj to RO tproj: lm block', lm

                        n_lm=0
                        do i=1, element%PAW%tab%lmn_size
                            if(element%PAW%tab%indlmn(4,i).ne.lm) cycle
                                n_lm=n_lm+1
                        enddo

                        allocate(L_ij(n_lm,n_lm))
                        allocate(U_ij(n_lm,n_lm))
                        allocate(lam_i(n_lm))
                        allocate(Q_ij(n_lm,n_lm))
                        allocate(sij_lm(n_lm,n_lm))
                        allocate(tzeta(imax,n_lm))
                        allocate(tproj_lm(imax,n_lm))
                        allocate(ffspl_lm(size(element%PAW%ffspl,1),2,n_lm))
                        allocate(ffzeta(size(element%PAW%ffspl,1),2,n_lm))

                        allocate(lmn_of_lm(n_lm))

                        L_ij=0.0_dp
                        !From t_proj to t_proj_RO
                        n_lm=0
                        do i=1, element%PAW%tab%lmn_size
                            if(element%PAW%tab%indlmn(4,i).ne.lm) cycle
                            n_lm=n_lm+1
                            tproj_lm(:,n_lm)=element%PAW%tab%tproj(:imax,element%PAW%tab%indlmn(5,i))
                            ffspl_lm(:,:,n_lm)=element%PAW%ffspl(:,:,element%PAW%tab%indlmn(5,i))
                            lmn_of_lm(n_lm)=i
                        enddo

                        do i=1,n_lm
                            do j=1, n_lm
                                tzeta(:imax,1)=tproj_lm(:,i)*tproj_lm(:,j)
                                call simp_gen(L_ij(i,j),tzeta(:,1),element%PAW%rad)
                                sij_lm(i,j)=element%PAW%sij(lmn_of_lm(i),lmn_of_lm(j))
                            enddo
                        enddo
                        
                        call eigh(L_ij,lam_i, U_ij)
                        !L_ij(i,j)-sum(U_ij(:,i)*lam_i(:)*U_ij(:,j))=0
                        ! lam_i(i)*delta_ij = sum(U_ij(:,i)*Q_ij(:,j)) [Q_ij(i,j)=sum(L_ij(i,:)*U_ij(:,j))]

                        ! |zeta>=(|P>.dot.U)/sqrt[lam]
                        call dgemm('N','N', imax, n_lm, n_lm, 1.0_dp, &
                        tproj_lm, imax, U_ij, n_lm, &
                        0.0_dp,  tzeta, imax)
                        do i=1, n_lm
                            tzeta(:imax,i)=tzeta(:imax,i)/sqrt(lam_i(i))
                        enddo

                        call dgemm('N','N', size(ffspl_lm,1), n_lm, n_lm, 1.0_dp, &
                        ffspl_lm(:,1,:), size(ffspl_lm,1), U_ij, n_lm, &
                        0.0_dp,  ffzeta(:,1,:), size(ffspl_lm,1))
                        call dgemm('N','N', size(ffspl_lm,1), n_lm, n_lm, 1.0_dp, &
                        ffspl_lm(:,2,:), size(ffspl_lm,1), U_ij, n_lm, &
                        0.0_dp,  ffzeta(:,2,:), size(ffspl_lm,1))
                        do i=1, n_lm
                            ffzeta(:,:,i)=ffzeta(:,:,i)/sqrt(lam_i(i))
                        enddo

                        !Q_ij = s.U
                        call dgemm('N','N', n_lm, n_lm, n_lm, 1.0_dp, &
                                sij_lm, n_lm, U_ij, n_lm, &
                                0.0_dp,  Q_ij, n_lm)
                        !L_ij = Ut.s.U 
                        
                        call dgemm('T','N', n_lm, n_lm, n_lm, 1.0_dp, &
                                U_ij, n_lm, Q_ij, n_lm, &
                                0.0_dp,  L_ij, n_lm)
                        !L_ij => lam^1/2 L_ij lam^1/2
                        do i=1, n_lm
                            L_ij(:,i)=L_ij(:,i)*sqrt(lam_i(:))*sqrt(lam_i(i))
                        enddo

                        call eigh(L_ij,lam_i, U_ij)

                        ! |nu>=(|zeta>.dot.U)
                        call dgemm('N','N', imax, n_lm, n_lm, 1.0_dp, &
                        tzeta, imax, U_ij, n_lm, &
                        0.0_dp,  tproj_lm, imax)

                        call dgemm('N','N', size(ffspl_lm,1), n_lm, n_lm, 1.0_dp, &
                        ffspl_lm(:,1,:), size(ffspl_lm,1), U_ij, n_lm, &
                        0.0_dp,  ffspl_lm(:,1,:), size(ffspl_lm,1))
                        call dgemm('N','N', size(ffspl_lm,1), n_lm, n_lm, 1.0_dp, &
                        ffspl_lm(:,2,:), size(ffspl_lm,1), U_ij, n_lm, &
                        0.0_dp,  ffspl_lm(:,2,:), size(ffspl_lm,1))

                        do i=1, n_lm
                            j=lmn_of_lm(i)
                            element%PAW%tproj_RO(:,j)=tproj_lm(:,i)
                            element%PAW%ffspl_Ort(:,:,j)=ffspl_lm(:,:,i)
                            element%PAW%o_i(j)=max(lam_i(i),-0.997_dp)
                        enddo

                        deallocate(L_ij)
                        deallocate(U_ij)
                        deallocate(Q_ij)
                        deallocate(sij_lm)
                        deallocate(lmn_of_lm)
                        deallocate(lam_i)
                        deallocate(tzeta)
                        deallocate(tproj_lm)
                        deallocate(ffzeta)
                        deallocate(ffspl_lm)
                    enddo

                    if(parallel%myid.eq.0) print *, 's_RO_i, i=> lmn'
                    do j=1, element%PAW%tab%lmn_size
                            if(parallel%myid.eq.0) write(6,'(E15.5)') element%PAW%o_i(j)

                    enddo
            enddo
        endif

        !Set the local pseudopotentials
        w2=>grids(1)%G2

        do e = 1, n_elements_total
            element=>elements(e)
            call allocate_local_fields_G(element%Ven0G, grids(1))
            call allocate_local_fields_G(element%dVen0GdG2,grids(1))

            element%Ven0G= 0.0_dp
            element%dVen0GdG2 =0.0_dp

            if(element%PP_type.eq.HGH) then

                C=>element%HGH%C
                rloc=>element%HGH%rloc
                where(w2>tiny(1.0_dp))
                    element%Ven0G = -4.0_dp*pi*element%Zion/w2*exp(-w2*rloc**2/2.0_dp) &
                        + sqrt(8.0_dp*pi**3)*rloc**3*exp(-w2*rloc**2/2.0_dp) &
                        *(C(1) +C(2)*(3.0_dp-w2*rloc**2) &
                        +C(3)*(15.0_dp-10.0_dp*w2*rloc**2 + w2**2*rloc**4) &
                        +C(4)*(105.0_dp-105.0_dp*w2*rloc**2 + 21.0_dp*w2**2*rloc**4 -w2**3*rloc**6) &
                        )
                    element%dVen0GdG2 = -rloc**2/2.0_dp*element%Ven0G &
                        +4.0_dp*pi*element%Zion/w2**2*exp(-w2*rloc**2/2.0_dp) &
                        + sqrt(8.0_dp*pi**3)*rloc**3*exp(-w2*rloc**2/2.0_dp) & 
                        *(C(2)*(-rloc**2) &
                        +C(3)*(-10.0_dp*rloc**2 + 2.0*w2*rloc**4) &
                        +C(4)*(-105.0_dp*rloc**2 + 42.0_dp*w2*rloc**4 -3.0_dp*w2**2*rloc**6) &
                        )
                elsewhere
                end where

            else if(element%PP_type.eq.PAW) then
                call allocate_local_fields_G(d2Ven0G,grids(1))
                d2Ven0G=0.0_dp

                q2vq(1:product(grids(1)%Ng_local))=>element%Ven0G(:,:,:)
                g2vec(1:product(grids(1)%Ng_local))=>grids(1)%G2(:,:,:)
                dq2vq(1:product(grids(1)%Ng_local))=>element%dVen0GdG2(:,:,:)
                d2q2vq(1:product(grids(1)%Ng_local))=>d2Ven0G(:,:,:)

                call spline3ders_unsorted(element%PAW%qgrid_vl(:),element%PAW%vlspl(:,1),sqrt(g2vec)*0.5_dp/pi, &
                q2vq, dq2vq, d2q2vq)

                deallocate(d2Ven0G)

                flush(6)
                call parallel_wait(parallel)
                where(grids(1)%G2.gt.tiny(1.0_dp))            
                    element%Ven0G=element%Ven0G/grids(1)%G2*4.0_dp*pi*pi
                    element%dVen0GdG2=element%dVen0GdG2-sqrt(grids(1)%G2)*element%Ven0G/pi
                    element%dVen0GdG2=element%dVen0GdG2*pi/sqrt(grids(1)%G2)/grids(1)%G2
                endwhere

                call allocate_local_fields_G(element%Den0G,grids(1))
                call allocate_local_fields_G(element%dDen0GdG2,grids(1))
                if(element%PAW%tab%usetcore.eq.1) call allocate_local_fields_G(d2Den0G,grids(1))
                call allocate_local_fields_G(element%Den1G,grids(1))

                element%Den1G= 0.0_dp
                element%Den0G= 0.0_dp
                element%dDen0GdG2 =0.0_dp
                g2vec(1:product(grids(1)%Ng_local))=>grids(1)%G2(:,:,:)

                if(element%PAW%tab%usetcore.eq.1) then
                    q2vq(1:product(grids(1)%Ng_local))=>element%Den0G(:,:,:)
                    dq2vq(1:product(grids(1)%Ng_local))=>element%dDen0GdG2(:,:,:)
                    d2q2vq(1:product(grids(1)%Ng_local))=>d2Den0G(:,:,:)

                    call spline3ders_unsorted(element%PAW%qgrid_vl(:),element%PAW%tab%tcorespl(:,1),sqrt(g2vec)*0.5_dp/pi, &
                    q2vq, dq2vq, d2q2vq)
                    where(grids(1)%G2.gt.tiny(1.0_dp))            
                    element%dDen0GdG2(:,:,:)=element%dDen0GdG2(:,:,:)*0.25_dp/pi/sqrt(grids(1)%G2)
                    endwhere
                    deallocate(d2Den0G)
                endif

                if(element%PAW%tab%has_tvale.eq.1) then
                    tvale_local(1:product(grids(1)%Ng_local))=>element%Den1G(:,:,:)
                    tvale_local = spline3_unsorted(element%PAW%qgrid_vl(:),element%PAW%tab%tvalespl(:,1),sqrt(g2vec)*0.5_dp/pi)
                endif

            else if(element%PP_type.eq.None) then
                where(w2>tiny(1.0_dp))
                element%Ven0G = -4.0_dp*pi*element%Zion/w2
                element%dVen0GdG2= 4.0_dp*pi*element%Zion/w2**2
                elsewhere
                end where
            endif 
            element%Ven0G = element%Ven0G*grids(1)%cutden
            element%dVen0GdG2 = element%dVen0GdG2*grids(1)%cutden 
            !This needs explicitely doing since Paw splines will give value but 
            !convention is to be 0 (and integration with compenesation density will give non-zero energy/stress shift) 
            if(grids(1)%G2(1,1,1).lt.tiny(1.0_dp)) element%Ven0G(1,1,1)=0.0_dp
        enddo

        !Distribute atoms of elements over same node processors consecutively

        do e = 1, n_elements_total
            element=>elements(e)
            allocate(element%counts(parallel%nproc_node))
            allocate(element%displs(parallel%nproc_node))
            element%displs=0
            element%counts(:)=element%n_atoms_of_element/parallel%nproc_node
            element%counts(:mod(element%n_atoms_of_element,parallel%nproc_node))= &
            element%counts(:mod(element%n_atoms_of_element,parallel%nproc_node)) +1
            element%my_natom=element%counts(parallel%myid_node+1)
            element%my_atom_min=0
            element%my_atom_max=-1
            element%displs(1)=0
    
            do j=2, parallel%nproc_node
                element%displs(j)=element%displs(j-1) + element%counts(j-1)
            enddo
            element%my_atom_min=element%displs(parallel%myid_node+1) + 1
            element%my_atom_max=element%my_atom_min + (element%my_natom -1)
            print *, e , element%my_atom_min, element%my_atom_max
        enddo

        do e = 1, n_elements_total
            element=>elements(e)
            if(element%PA_type.eq.Real_local) then
                element%RL%max_sphere(:) = &
                    2*(ceiling(element%RL%radius/grids(2)%dR(:)))+1
        
                if(parallel%myid.eq.0) print *, 'Element: ', e, ' Sphere size, rough_grid: ', &
                                                    element%RL%max_sphere(:)
                element%RL%max_sphere_dense(:) = &
                    1+(element%RL%max_sphere(:)-1)*element%RL%OH_factor
        
                if(parallel%myid.eq.0) print *, 'Element: ', e, ' Sphere size, dense_grid (interpolation): ', &
                                                    element%RL%max_sphere_dense(:)                                    
            endif
            if(element%PP_type.eq.PAW) then
                element%PAW%max_sphere_fine(:) = &
                2*(ceiling(element%PAW%radius/grids(1)%dR(:)))+1

                if(parallel%myid.eq.0)  print *, 'Element: ', e, ' Sphere size, fine_grid: ', &
                                element%PAW%max_sphere_fine(:)
            endif
        enddo

        flush(6)

    end subroutine

    subroutine initialize_atoms(atoms, elements, parallel)
        use atom_type, only : atom_struct
        use element_type, only : element_struct
        use simulation_type, only : all_PAW_struct
        use system_type, only : system_struct

        use parallel_type, only : parallel_struct
        use parallel_mod, only : parallel_wait
        use m_paral_atom, only : get_my_natom, get_my_atmtab
        use m_paw_ij, only : paw_ij_init
        use m_pawcprj, only : pawcprj_alloc, pawcprj_getdim
        use m_pawfgrtab, only : pawfgrtab_init
        use m_paw_occupancies, only : initrhoij
        use numerics_type, only : numerics_struct

        type(parallel_struct), intent(in) :: parallel
        type(element_struct), intent(in) :: elements(:)
        type(atom_struct), intent(inout), target :: atoms(:)
        type(atom_struct), pointer :: atom

        integer :: i , n_atoms_total
            if(parallel%myid.eq.0) print *, 'Initialize Atoms'
            atoms(:)%mine=.false.
            n_atoms_total=size(atoms)

            do i = 1, n_atoms_total
                    atom=>atoms(i)
                call parallel_wait(parallel)
                if( (atom%atom_of_element.le.elements(atom%element)%my_atom_max) .and. &
                    (atom%atom_of_element.ge.elements(atom%element)%my_atom_min) ) &
                    atom%mine=.true.
            enddo

    end subroutine

    subroutine initialize_atoms_temperature(atoms, elements, td, system)
        use types, only : dp
        use constants, only: pi
        use atom_type, only : atom_struct
        use td_type, only : td_struct, IsoKinetic
        use element_type, only : element_struct
        use parallel_type, only : parallel_struct
        use system_type, only : system_struct
        type(system_struct), intent(in) :: system

        type(element_struct), intent(inout) :: elements(:)
        type(atom_struct), intent(inout) :: atoms(:)
        type(td_struct), intent(in) :: td

        integer :: n_not_frozen(size(elements)), i, e
        real(dp) :: PCM(3), theta, phi, norm

        PCM(:)=0.0_dp
        if(any(.not.atoms(:)%frozen_V)) then
            n_not_frozen(:)=0
            do i=1,size(atoms)
                if(.not.atoms(i)%frozen_V) then
                    e=atoms(i)%element
                     n_not_frozen(e)= n_not_frozen(e) + 1
                     PCM(:)=PCM(:)+atoms(i)%P(:)
                endif
            enddo
        endif

        elements(:)%isokin_ec=0.0_dp
        do i = 1, size(atoms)
            if(.not.atoms(i)%frozen_V) atoms(i)%P(:)=atoms(i)%P(:)-PCM(:)/sum(n_not_frozen)
            if(td%thermostat.eq.IsoKinetic) then
                if(.not.atoms(i)%frozen_V) then
                    e=atoms(i)%element
                    elements(e)%isokin_ec=elements(e)%isokin_ec+sum(atoms(i)%P(:)**2)/elements(e)%M
                endif
            endif
        enddo

        if(td%thermostat.eq.IsoKinetic) then
            do e = 1, size(elements)
                !if P was left undefined Randomize it and set it to temperature
                if(elements(e)%isokin_ec.lt.tiny(1.0_dp).and..not.all(atoms(:)%frozen_V)) then
                    elements(e)%isokin_ec=0.0_dp
                    do i = 1, size(atoms)
                        if(e.eq.atoms(i)%element.and..not.atoms(i)%frozen_V) then
                            theta=0.0_dp;phi=0.0_dp;norm=0.0_dp
                            call random_number(theta)
                            theta=theta*2.0_dp*pi
                            call random_number(phi)
                            phi=phi*pi-0.5_dp*pi
                            call random_number(norm)
                            atoms(i)%P(1)=norm*cos(phi)*cos(theta)
                            atoms(i)%P(2)=norm*cos(phi)*sin(theta)
                            atoms(i)%P(3)=norm*sin(phi)
                            print *, 'Random P for atom' ,i,  atoms(i)%P(:), n_not_frozen(e), elements(e)%M
                        endif
                    enddo
                    PCM(:)=0.0_dp
                    do i = 1, size(atoms)
                        if(e.eq.atoms(i)%element.and..not.atoms(i)%frozen_V) then
                            PCM(:)=PCM(:)+atoms(i)%P(:)
                        endif
                    enddo           
                    do i = 1, size(atoms)
                        if(e.eq.atoms(i)%element.and..not.atoms(i)%frozen_V) then
                            atoms(i)%P(:)=atoms(i)%P(:)-PCM(:)/n_not_frozen(e)
                            elements(e)%isokin_ec=elements(e)%isokin_ec+sum(atoms(i)%P(:)**2)/elements(e)%M
                        endif
                    enddo
                endif
                do i = 1, size(atoms)
                    if(e.eq.atoms(i)%element.and..not.atoms(i)%frozen_V.and.(n_not_frozen(e).gt.1)) then
                        atoms(i)%P(:)=atoms(i)%P(:)/sqrt(elements(e)%isokin_ec/n_not_frozen(e)/3.0_dp) &
                        *sqrt(system%temperature)
                    endif
                enddo
                elements(e)%isokin_ec=0.0_dp
                do i = 1, size(atoms)
                    if(e.eq.atoms(i)%element.and..not.atoms(i)%frozen_V) then
                        elements(e)%isokin_ec=elements(e)%isokin_ec+sum(atoms(i)%P(:)**2)/elements(e)%M
                    endif
                enddo
            enddo
        endif
            
       end subroutine

       subroutine init_den_pot_orb(sim)
        use simulation_type, only: simulation_struct
        use Orbital_Free_Min, only : orbital_free_energy_min, orbital_free_energy
        use Local_ion, only : Calc_Local_PP_Energy, calculate_local_ion_potential, Calc_Core_Energy
        use Update_Potentials_and_Energies, only : Update_BOMD_Potentials
        use grids_mod, only : Gspace_grid_to_grid_transfer
        use read_in, only : read_in_density
        use fft, only: real_to_recip, recip_to_real
        use Stochastic_Mod, only : stochastic_vector_initializer
        use parallel_mod, only : parallel_task
        use odp_type, only : orbital_struct,field_struct
        use density, only : calc_Compensation_Charge_Density, calc_valance_Charge_Density_for_guess
        use Non_Local_ion, only : Apply_S_power
        use odp, only : allocate_field, deallocate_field

        type(simulation_struct), intent(inout), target :: sim
        type(orbital_struct), pointer :: ks_orbitals(:), stoc_orbitals(:)

        integer :: i, s_loc, k_loc
        real(dp) :: nval
        nval=0.0_dp
        !Initialization of Density
        if(sim%tasks%density_from_file) then
            if(sim%parallel%myid.eq.0) print *, 'Read Density from file'
            call read_in_density(sim%coarse_density,sim%grids, 2, sim%parallel)
                do i = 1, sim%density%n_s
                    call Gspace_grid_to_grid_transfer(sim%grids(2), sim%grids(1), sim%parallel, sim%coarse_density%of(i)%G, &
                                                    sim%density%of(i)%G, sim%coarse_to_fine)
                    call recip_to_real(sim%density%of(i), sim%grids)
                enddo
        else if(sim%tasks%read_orbitals_from_file) then
            if(sim%parallel%myid.eq.0) print *, 'Reading in Orbitals not supported yet'
            stop
        else
            if(sim%all_PAW%N_PAW_atoms.gt.0) then
                call calc_valance_Charge_Density_for_guess(sim%all_PAW, sim%grids, sim%atoms, sim%elements, sim%parallel)
                do i=1, sim%density%n_s 
                    sim%density%of(i)%R=sim%all_PAW%tnvale%of(i)%R
                    sim%density%of(i)%G=sim%all_PAW%tnvale%of(i)%G
                    if(sim%parallel%myid.eq.0) then
                        nval=real(sim%density%of(i)%G(1,1,1)*product(sim%grids(sim%density%of(i)%grid)%box_length))
                        print *, 'ne_paw_guess, spin', nval, i, 'nelec:', sim%system%nelec(i), 'diff:', sim%system%nelec(i)-nval
                    endif
                    call parallel_task('bcast', nval, sim%parallel, 'all', root=0)
                    sim%density%of(i)%R= sim%density%of(i)%R + &
                         (sim%system%nelec(i)-nval)/product(sim%grids(sim%density%of(i)%grid)%box_length)
                    call real_to_recip(sim%density%of(i), sim%grids)
                    call Gspace_grid_to_grid_transfer(sim%grids(1), sim%grids(2), sim%parallel, sim%density%of(i)%G, &
                                                    sim%coarse_density%of(i)%G, sim%fine_to_coarse)
                    call recip_to_real(sim%coarse_density%of(i), sim%grids)
                enddo
            else
                do i=1, sim%density%n_s 
                    sim%density%of(i)%R=sim%system%nelec(i)/product(sim%grids(sim%density%of(i)%grid)%box_length)
                    call real_to_recip(sim%density%of(i), sim%grids)
                    sim%coarse_density%of(i)%R=sim%system%nelec(i)/product(sim%grids(sim%coarse_density%of(i)%grid)%box_length)
                    call real_to_recip(sim%coarse_density%of(i), sim%grids)
                enddo
            endif
        endif    
        
        !Start with an Orbital-Free Density calcualtion?
        if(sim%tasks%initial_orbital_free) then
            !Update the Potentials
            call Update_BOMD_Potentials(sim%grids, sim%system, sim%parallel, sim%density, &
            sim%potentials, sim%energies, sim%atoms, sim%elements, sim%xc, sim%tf, sim%all_PAW, sim%fine_to_coarse)
            if(sim%parallel%myid.eq.0) print *, 'Orbital Free Calculation'
            if(sim%parallel%myid_space.eq.0) then
                call  orbital_free_energy_min(sim%system, sim%density, sim%all_PAW, sim%grids,  sim%parallel, &
                 sim%potentials, sim%xc, sim%energies, sim%tf)
            endif
            do i = 1, sim%density%n_s
                call parallel_task('bcast', sim%density%of(i)%G, sim%parallel, 'space', root=0)
                call parallel_task('bcast', sim%density%of(i)%R, sim%parallel, 'space', root=0)

                call Gspace_grid_to_grid_transfer(sim%grids(1), sim%grids(2), sim%parallel, sim%density%of(i)%G, &
                                                    sim%coarse_density%of(i)%G, sim%fine_to_coarse)
                call recip_to_real(sim%coarse_density%of(i), sim%grids)
            enddo
        endif

        if(sim%all_PAW%N_PAW_atoms.gt.0) then
            call calc_Compensation_Charge_Density(sim%all_PAW, sim%grids, sim%atoms, sim%parallel)
            do i = 1, sim%density%n_s
                sim%density%of(i)%G=sim%density%of(i)%G + sim%all_PAW%rho_comp%of(i)%G
                sim%density%of(i)%R=sim%density%of(i)%R + sim%all_PAW%rho_comp%of(i)%R
            enddo
        endif

        !Update the Potentials
        if(sim%parallel%myid.eq.0) print *, 'Initial Potential Calcualtion'
        call Update_BOMD_Potentials(sim%grids, sim%system, sim%parallel, sim%density, &
        sim%potentials, sim%energies, sim%atoms, sim%elements, sim%xc, sim%tf, sim%all_PAW, sim%fine_to_coarse)

        if(.not.sim%tasks%read_orbitals_from_file) then
            if(sim%parallel%myid.eq.0) print *, 'Initialize Orbitals'
            do s_loc=1, size(sim%orbitals(:,:,:),3); do k_loc=1, size(sim%orbitals(:,:,:),2)

                if((sim%parallel%my_sub.eq.0).and.(sim%system%n_find_min.gt.0)) then
                    ks_orbitals(1:sim%system%n_find_min/sim%parallel%n_band_groups_sub)=>&
                        sim%orbitals(sim%system%find_min_start:sim%system%find_min_end, k_loc, s_loc)
                    call deterministic_vector_initializer(ks_orbitals(:),sim%grids, sim%parallel)
                endif
                
                if(sim%system%n_deterministic.gt.0) then
                    ks_orbitals(1:sim%system%n_deterministic/sim%parallel%n_band_groups)=>&
                        sim%orbitals(sim%system%deterministic_start:sim%system%deterministic_end, k_loc, s_loc)
                    call deterministic_vector_initializer(ks_orbitals(:),sim%grids, sim%parallel)
                endif

                if(sim%system%n_stochastic.gt.0) then
                    stoc_orbitals(1:sim%system%n_stochastic/sim%parallel%n_band_groups)=> &
                        sim%orbitals(sim%system%stochastic_start:sim%system%stochastic_end, k_loc, s_loc)
                        sim%stoc_norm=0.0_dp
                        sim%stoc_occ=0.0_dp
                    call stochastic_vector_initializer(stoc_orbitals(:), sim%system, &
                                                        sim%grids, sim%parallel, sim%stoc_norm, sim%stoc_occ, &
                                                        sim%all_PAW, sim%atoms, sim%elements)
                endif                    
            enddo; enddo
            if(sim%system%n_stochastic.gt.0) then 
                call parallel_task('sum', sim%stoc_norm(:,:,:), sim%parallel, 'space')
                call parallel_task('sum', sim%stoc_occ(:,:,:), sim%parallel, 'space')
            endif
        endif
        if(sim%parallel%myid.eq.0) print *, 'Initialization of Density, Potential & orbitals done'
            
    end subroutine

    subroutine deterministic_vector_initializer(orbitals,grids, parallel)
        use types, only : dp
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct, inter_grid_struct
        use odp_type, only : orbital_struct
        use operations_3D, only : integrate_3D_G
        use fft, only: real_to_recip
        use grids_mod, only : allocate_local_fields_R


        type(parallel_struct), intent(in) :: parallel
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid

        real(dp), allocatable :: random_number_field(:,:,:)
        integer:: i, s
        call allocate_local_fields_R(random_number_field, grids(orbitals(1)%of(1)%grid))
        do i=1,size(orbitals)
            do s=1,size(orbitals(i)%of)
                grid=>grids(orbitals(i)%of(s)%grid)
                call random_number(random_number_field)
                orbitals(i)%of(s)%R=random_number_field
                call real_to_recip(orbitals(i)%of(s), grids)
                orbitals(i)%of(s)%G=orbitals(i)%of(s)%G*grid%cutwf
                orbitals(i)%of(s)%G=orbitals(i)%of(s)%G/ &
                    sqrt(real(integrate_3D_G(conjg(orbitals(i)%of(s)%G)*orbitals(i)%of(s)%G, grid, parallel)))
            enddo
        enddo

        deallocate(random_number_field)
    end subroutine

    subroutine initialize_PAW(atoms, elements, parallel, system, grids, all_PAW, numerics)
        use atom_type, only : atom_struct
        use element_type, only : element_struct, PAW
        use simulation_type, only : all_PAW_struct
        use system_type, only : system_struct
        use grids_type, only :grid_struct
        use parallel_type, only : parallel_struct
        use parallel_mod, only : parallel_wait
        use m_paral_atom, only : get_my_natom, get_my_atmtab
        use m_paw_ij, only : paw_ij_init
        use m_pawcprj, only : pawcprj_alloc, pawcprj_getdim
        use m_pawfgrtab, only : pawfgrtab_init
        use m_paw_an, only : paw_an_init
        use m_paw_occupancies, only : initrhoij
        use numerics_type, only : numerics_struct
        use m_pawtab, only : pawtab_get_lsize
        use m_paw_denpot, only : pawdenpot
        use m_pawrhoij, only : pawrhoij_filter
        use m_pawrhoij, only: pawrhoij_unpack

        type(parallel_struct), intent(in) :: parallel
        type(element_struct), intent(in) :: elements(:)
        type(atom_struct), intent(inout), target :: atoms(:)
        type(all_PAW_struct), intent(inout), target :: all_PAW
        type(system_struct), intent(in) ::  system
        type(numerics_struct), intent(in) :: numerics
        type(grid_struct), intent(in) :: grids(:)

        integer, allocatable :: lpawu(:), lexexch(:)
        real(dp), allocatable :: spinat(:,:)
        integer, allocatable  :: l_size_atm(:)
        real(dp) :: nucdipmom(3,size(atoms))

        real(dp) :: compch_sph,epaw,epawdc
        integer :: i , j, k, e, n_atoms_total, n_paw_ij_total, klmn,ip, jp, at,all_proj

        n_atoms_total=size(atoms)
        nucdipmom=0.0_dp
        all_PAW%epaw=0.0_dp
        all_PAW%vpotzero=0.0_dp

        !!!MORE PAW INITIALIZATION

        all_PAW%N_PAW_atoms=0
        do i = 1, n_atoms_total
            if(elements(atoms(i)%element)%PP_type.ne.PAW) cycle
            all_PAW%N_PAW_atoms=all_PAW%N_PAW_atoms+1
            all_PAW%nattyp(elements(atoms(i)%element)%PAW%index)= &
            all_PAW%nattyp(elements(atoms(i)%element)%PAW%index) + 1
        enddo

        do e=1, size(elements)
            if(elements(e)%PP_type.ne.PAW) cycle
            all_PAW%znucl(elements(e)%PAW%index)=elements(e)%Znuc
        enddo

        if(all_PAW%N_PAW_atoms.gt.0) then
            allocate(all_PAW%ij(all_PAW%N_PAW_atoms))
            allocate(all_PAW%typat(all_PAW%N_PAW_atoms))
            allocate(all_PAW%all_at(all_PAW%N_PAW_atoms))
        else 
            if(parallel%myid.eq.0) print *,  'No PAW atoms'
            return
        endif
        j=0
        n_paw_ij_total=0

        all_PAW%tncore%n_s=all_PAW%nspden
        all_PAW%tnvale%n_s=all_PAW%nspden
        all_PAW%rho_comp%n_s=all_PAW%nspden
        all_PAW%gr_rho_comp(1:3)%n_s=all_PAW%nspden

        all_proj=0
        do i = 1, n_atoms_total
            e=atoms(i)%element
            if(elements(atoms(i)%element)%PP_type.ne.PAW) cycle
            j=j+1
            all_PAW%typat(j)=elements(atoms(i)%element)%PAW%index
            all_PAW%all_at(j)=i

            all_proj=all_proj+elements(e)%n_proj

            allocate(atoms(i)%PAW)
            allocate(atoms(i)%PAW%obar_i(elements(e)%PAW%tab%lmn_size))
            allocate(atoms(i)%PAW%Dij(elements(e)%PAW%tab%lmn_size,elements(e)%PAW%tab%lmn_size, elements(e)%PAW%nspden))
            allocate(atoms(i)%PAW%rho_ij(elements(e)%PAW%tab%lmn_size,elements(e)%PAW%tab%lmn_size, elements(e)%PAW%nspden))
            allocate(atoms(i)%PAW%Ei(elements(e)%PAW%tab%lmn_size))
            allocate(atoms(i)%PAW%Uij(elements(e)%PAW%tab%lmn_size,elements(e)%PAW%tab%lmn_size))
            allocate(atoms(i)%PAW%Sij(elements(e)%PAW%tab%lmn_size,elements(e)%PAW%tab%lmn_size))
            allocate(atoms(i)%PAW%Sinv_block(elements(e)%PAW%tab%lmn_size,elements(e)%PAW%tab%lmn_size))

            atoms(i)%PAW%Dij=0.0_dp
            atoms(i)%PAW%ij=>all_PAW%ij(j)
            atoms(i)%PAW%rad=>elements(atoms(i)%element)%PAW%rad
            atoms(i)%PAW%tab=>elements(atoms(i)%element)%PAW%tab
            atoms(i)%PAW%typat=>all_PAW%typat(j)

            n_paw_ij_total=n_paw_ij_total+elements(atoms(i)%element)%PAW%tab%lmn2_size
            atoms(i)%PAW%index=j
            atoms(i)%PAW%an=>null()
            atoms(i)%PAW%fgrtab=>null()
            atoms(i)%PAW%rhoij=>null()
            atoms(i)%PAW%cprj=>null()
            atoms(i)%PAW%dimcprj=>null()
        enddo

        allocate(all_PAW%proj_overlap(all_proj,all_proj))
        allocate(all_PAW%S_inv_full(all_proj,all_proj))

        allocate(all_PAW%Rhoij_mix_R_i(numerics%pulay%n+1, all_PAW%nspden))
        allocate(all_PAW%Rhoij_mix_F_i(numerics%pulay%n+1, all_PAW%nspden))
        allocate(all_PAW%Rhoij_mix_Ri(n_paw_ij_total, numerics%pulay%n, all_PAW%nspden))
        allocate(all_PAW%Rhoij_mix_Fi(n_paw_ij_total, numerics%pulay%n, all_PAW%nspden))
        allocate(all_PAW%Rhoij_mix_0(n_paw_ij_total,all_PAW%nspden))
        allocate(all_PAW%Rhoij_mix_m1(n_paw_ij_total,all_PAW%nspden))
        allocate(all_PAW%Rhoij_mix_m2(n_paw_ij_total,all_PAW%nspden))
        allocate(all_PAW%Rhoij_mix_fm1(n_paw_ij_total,all_PAW%nspden))


        all_PAW%paral_atom=.true.
        all_PAW%my_atmtab=>null()
        call get_my_natom(parallel%comm_space,all_PAW%my_natom, all_PAW%N_PAW_atoms)
        call get_my_atmtab(parallel%comm_space,all_PAW%my_atmtab,all_PAW%my_atmtab_allocated, all_PAW%paral_atom, &
                            all_PAW%N_PAW_atoms)
        call pawtab_get_lsize(Pawtab=all_PAW%tab,l_size_atm=l_size_atm,natom=all_PAW%N_PAW_atoms,typat=all_PAW%typat, &
                            mpi_atmtab=all_PAW%my_atmtab)
        

        allocate(all_PAW%an(all_PAW%my_natom))
        allocate(all_PAW%my_ij(all_PAW%my_natom))
        allocate(all_PAW%fgrtab(all_PAW%my_natom))
        allocate(all_PAW%rhoij(all_PAW%my_natom))
        allocate(all_PAW%cprj(all_PAW%my_natom,system%n_orbitals))
        allocate(all_PAW%dimcprj(all_PAW%my_natom))

        if(all_PAW%my_natom.gt.0) then
            j=0
            k=0
            n_paw_ij_total=0
            do i = 1, n_atoms_total
                if(elements(atoms(i)%element)%PP_type.ne.PAW) cycle
                j=j+1
                if(all_PAW%paral_atom) then
                    if(all_PAW%my_atmtab_allocated) then; if(j.eq.all_PAW%my_atmtab(k+1)) then
                        k=k+1
                        atoms(i)%PAW%an=>all_PAW%an(k)
                        atoms(i)%PAW%ij=>all_PAW%my_ij(k)
                        atoms(i)%PAW%fgrtab=>all_PAW%fgrtab(k)
                        atoms(i)%PAW%rhoij=>all_PAW%rhoij(k)
                        atoms(i)%PAW%cprj(1:system%n_orbitals)=>all_PAW%cprj(k,1:system%n_orbitals)
                        atoms(i)%PAW%dimcprj=>all_PAW%dimcprj(k)
                    endif; endif
                else
                        atoms(i)%PAW%an=>all_PAW%an(j)
                        atoms(i)%PAW%ij=>all_PAW%my_ij(j)
                        atoms(i)%PAW%fgrtab=>all_PAW%fgrtab(j)
                        atoms(i)%PAW%rhoij=>all_PAW%rhoij(j)
                        atoms(i)%PAW%cprj(1:system%n_orbitals)=>all_PAW%cprj(j,1:system%n_orbitals)
                        atoms(i)%PAW%dimcprj=>all_PAW%dimcprj(j)
                endif
                if(k.eq.all_PAW%my_natom) exit
            enddo
        endif
        if(parallel%myid.eq.0) Print *, 'LibPAW an init'
        call paw_an_init(all_PAW%an,natom=all_PAW%N_PAW_atoms, ntypat=all_PAW%N_PAW_elements, &
                nkxc1=0,nk3xc1=0,nspden=all_PAW%nspden,cplex=1,pawxcdev=all_PAW%pawxcdev,typat=all_PAW%typat, &
                pawang=all_PAW%ang,pawtab=all_PAW%tab,&
                has_vhartree=1,has_vxc=1,has_vxctau=0,has_vxcval=0,has_kxc=0,has_k3xc=0,has_vxc_ex=0, & ! optional arguments
                mpi_atmtab=all_PAW%my_atmtab,comm_atom=parallel%comm_space) ! optional arguments (parallelism)
        if(parallel%myid.eq.0) Print *, 'LibPAW IJ init'
        call paw_ij_init(all_PAW%my_ij,cplex=1,nspinor=system%n_spinor,nsppol=system%n_spin,nspden=all_PAW%nspden, &
            pawspnorb=all_PAW%spnorb,natom=all_PAW%N_PAW_atoms, ntypat=all_PAW%N_PAW_elements, &
            typat=all_PAW%typat,pawtab=all_PAW%tab,&
            has_dij=1,has_dij0=1,has_dijfock=0,has_dijfr=1,has_dijhartree=1,has_dijhat=1,& ! Optional
            has_dijxc=1,has_dijxc_hat=1,has_dijxc_val=0,has_dijnd=0,has_dijso=0,has_dijU=0, has_dijexxc=0, &  ! Optional
            has_exexch_pot=0,has_pawu_occ=0,&!nucdipmom, & ! Optional
            mpi_atmtab=all_PAW%my_atmtab,comm_atom=parallel%comm_space) ! optional arguments (parallelism)
        if(parallel%myid.eq.0) Print *, 'LibPAW cprj init'
        call pawcprj_getdim(dimcprj=all_PAW%dimcprj,natom=all_PAW%my_natom,nattyp=all_PAW%nattyp, &
            ntypat=all_PAW%N_PAW_elements,typat=all_PAW%typat,Pawtab=all_PAW%tab,sort_mode='R')
        call pawcprj_alloc(cprj=all_PAW%cprj,ncpgr=0,nlmn=all_PAW%dimcprj)
        if(parallel%myid.eq.0) Print *, 'LibPAW Fine Grid table init'
        call pawfgrtab_init(all_PAW%fgrtab,cplex=1,l_size_atm=l_size_atm, &
            nspden=all_PAW%nspden,typat=all_PAW%typat, &
            mpi_atmtab=all_PAW%my_atmtab,comm_atom=parallel%comm_space) ! optional arguments (parallelism)
        allocate(lpawu(all_PAW%N_PAW_elements))
        allocate(lexexch(all_PAW%N_PAW_elements))
        allocate(spinat(3,all_PAW%N_PAW_atoms))
        lpawu=-1
        lexexch=-1
        spinat=-1
        if(parallel%myid.eq.0) Print *, 'LibPAW rho_ij init'
        call initrhoij(cpxocc=1,lexexch=lexexch,lpawu=lpawu,my_natom=all_PAW%my_natom,natom=all_PAW%N_PAW_atoms, &
            nspden=all_PAW%nspden,nspinor=system%n_spinor,nsppol=system%n_spin,&
            ntypat=all_PAW%N_PAW_elements,pawrhoij=all_PAW%rhoij,pawspnorb=all_PAW%spnorb, &
            pawtab=all_PAW%tab,qphase=1,spinat=spinat,typat=all_PAW%typat, &
            ngrhoij=0,nlmnmix=0,use_rhoij_=1,use_rhoijres=1,& ! optional arguments
            mpi_atmtab=all_PAW%my_atmtab,comm_atom=parallel%comm_space) ! optional arguments (parallelism)

        deallocate(lpawu)
        deallocate(lexexch)
        deallocate(spinat)

        if(parallel%myid.eq.0) Print *, 'LibPAW rho_ij unpack'
        call pawrhoij_unpack(all_PAW%rhoij)
        do at=1, size(atoms)
            if(elements(atoms(at)%element)%PP_type.ne.PAW) cycle
            if(associated(atoms(at)%PAW%rhoij)) then !The Libpaw struct
                e=atoms(at)%element
                do klmn=1,atoms(at)%PAW%rhoij%lmn2_size
                        ip=elements(e)%PAW%tab%indklmn(7,klmn)
                        jp=elements(e)%PAW%tab%indklmn(8,klmn)
                        atoms(at)%PAW%rho_ij(ip,jp,:)=atoms(at)%PAW%rhoij%rhoij_(klmn,:)
                enddo
            endif
        enddo
        if(parallel%myid.eq.0) Print *, 'LibPAW denpot init'

        call pawdenpot(compch_sph,epaw,epawdc,ipert=0,ixc=all_PAW%ixc, &
        my_natom=all_PAW%my_natom,natom=all_PAW%N_PAW_atoms,nspden=all_PAW%nspden,ntypat=all_PAW%N_PAW_elements, &
        nucdipmom=nucdipmom(:,:all_PAW%N_PAW_atoms), nzlmopt=-1,option=0,paw_an=all_PAW%an,paw_an0=all_PAW%an, &
        paw_ij=all_PAW%my_ij,pawang=all_PAW%ang,pawprtvol=0,pawrad=all_PAW%rad,pawrhoij=all_PAW%rhoij, &
        pawspnorb=all_PAW%spnorb,pawtab=all_PAW%tab,pawxcdev=all_PAW%pawxcdev,spnorbscl=0.0_dp, &
        xclevel=all_PAW%xclevel,xc_denpos=all_PAW%xc_denpos,ucvol=product(grids(1)%Box_Length),znucl=all_PAW%znucl,&
        mpi_atmtab=all_PAW%my_atmtab,comm_atom=parallel%comm_space, & ! optional arguments (parallelism)
        vpotzero=all_PAW%vpotzero)!,hyb_mixing,hyb_mixing_sr) ! optional arguments

        if(parallel%myid.eq.0) Print *, 'LibPAW init done'

    end subroutine

    subroutine finalize_all_paw(all_PAW, elements, atoms, grids)
        use simulation_type, only : all_PAW_struct
        use atom_type, only : atom_struct
        use element_type, only : element_struct
        use odp, only : deallocate_DenPot
        use grids_type, only :grid_struct

        type(all_PAW_struct), intent(inout), target :: all_PAW
        type(element_struct), intent(inout) :: elements(:)
        type(atom_struct), intent(inout) :: atoms(:)
        type(grid_struct), intent(in) :: grids(:)

        integer :: i 
            if(allocated(all_PAW%an)) deallocate(all_PAW%an)
            if(allocated(all_PAW%my_ij)) deallocate(all_PAW%my_ij)
            if(allocated(all_PAW%fgrtab)) deallocate(all_PAW%fgrtab)
            if(allocated(all_PAW%rhoij)) deallocate(all_PAW%rhoij)
            if(allocated(all_PAW%cprj)) deallocate(all_PAW%cprj)
            if(allocated(all_PAW%dimcprj)) deallocate(all_PAW%dimcprj)

            if(allocated(all_PAW%rad)) deallocate(all_PAW%rad) 
            if(allocated(all_PAW%tab)) deallocate(all_PAW%tab) 
            if(allocated(all_PAW%list)) deallocate(all_PAW%list)
            if(allocated(all_PAW%nattyp))deallocate(all_PAW%nattyp)
            if(allocated(all_PAW%znucl))deallocate(all_PAW%znucl)

            if(allocated(all_PAW%proj_overlap))deallocate(all_PAW%proj_overlap)
            if(allocated(all_PAW%S_inv_full))deallocate(all_PAW%S_inv_full)


            if(allocated(all_PAW%Rhoij_mix_R_i))deallocate(all_PAW%Rhoij_mix_R_i)
            if(allocated(all_PAW%Rhoij_mix_F_i)) deallocate(all_PAW%Rhoij_mix_F_i)
            if(allocated(all_PAW%Rhoij_mix_Ri))deallocate(all_PAW%Rhoij_mix_Ri)
            if(allocated(all_PAW%Rhoij_mix_Fi))deallocate(all_PAW%Rhoij_mix_Fi)
            if(allocated(all_PAW%Rhoij_mix_0))deallocate(all_PAW%Rhoij_mix_0)
            if(allocated(all_PAW%Rhoij_mix_m1))deallocate(all_PAW%Rhoij_mix_m1)
            if(allocated(all_PAW%Rhoij_mix_m2))deallocate(all_PAW%Rhoij_mix_m2)
            if(allocated(all_PAW%Rhoij_mix_fm1)) deallocate(all_PAW%Rhoij_mix_fm1)

            if(all_PAW%N_PAW_atoms.gt.0) then
                deallocate(all_PAW%ij)
                deallocate(all_PAW%typat)
                deallocate(all_PAW%all_at)
            endif

            do i=1,size(elements)
                if(allocated(elements(i)%PAW)) deallocate(elements(i)%PAW)
                if(allocated(elements(i)%Den0G)) deallocate(elements(i)%Den0G)
                if(allocated(elements(i)%dDen0GdG2)) deallocate(elements(i)%dDen0GdG2)
                if(allocated(elements(i)%Den1G)) deallocate(elements(i)%Den1G)
            enddo
            do i=1,size(atoms)
                if(allocated(atoms(i)%PAW)) deallocate(atoms(i)%PAW)
            enddo

            if(allocated(all_PAW%tncore%of)) call deallocate_DenPot(all_PAW%tncore,grids)
            if(allocated(all_PAW%tncore%of)) call deallocate_DenPot(all_PAW%rho_comp,grids)
            do i=1,3
                if(allocated(all_PAW%gr_rho_comp(i)%of)) call deallocate_DenPot(all_PAW%gr_rho_comp(i),grids)
            enddo


    end subroutine

    subroutine finalize_all_HGH(elements, atoms)
        use atom_type, only : atom_struct
        use element_type, only : element_struct

        type(element_struct), intent(inout) :: elements(:)
        type(atom_struct), intent(inout) :: atoms(:)
        integer :: i

        do i=1,size(elements)
            if(allocated(elements(i)%HGH)) deallocate(elements(i)%HGH)
        enddo
     

    end subroutine

    

end module