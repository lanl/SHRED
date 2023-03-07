module Kubo_Greenwood
    use types, only : dp

    implicit none

    public :: Mixed_Kubo_Greenwood


    contains

    subroutine Mixed_Kubo_Greenwood(sim)
        use simulation_type, only: simulation_struct
        use initialization_finalization, only : initialize_Kubo_Greenwood
        use parallel_mod,only : parallel_task, parallel_wait
        use Time_Dependent, only : SIL_Propagate, SIL_Initialize
        use odp_type, only : orbital_struct
        use constants, only : pi, i_, FMT_real
        use simulation_type, only: simulation_struct
        use parallel_mod,only : parallel_task, parallel_wait
        use odp_type, only : orbital_struct, field_struct
        use constants, only : pi, i_
        use odp, only: allocate_field, deallocate_field, allocate_orbital, deallocate_orbital 
        use grids_type, only : grid_struct
        use operations_3D, only: integrate_3D_G
        use Non_Local_ion, only : Apply_projectors_RL, Calculate_Projector_overlaps, Apply_S_inverse
        use Apply_Hamiltonian, only: Apply_H, Apply_SinvH
        use fft, only: recip_to_real, real_to_recip
        use Eigen_LOBPCG, only : Orbitals_Project_out
        use odp_type, only : deterministic, stochastic, lowest, full, find_min, eigen_vector
        use Current,only: Dipole_Matrix_Elements, J_to_orbitals, J_Matrix_Elements, orbitals_both_current
        use KG_type, only : SIL

        type(simulation_struct), intent(inout), target :: sim
        type(orbital_struct), pointer :: stoc_orbitals(:), ks_orbitals(:)
        type(orbital_struct), allocatable :: J_orbitals(:,:,:,:,:), unfiltered_stoc_orbitals(:,:,:), J_orbitals_t(:,:,:,:,:)
        type(grid_struct), pointer :: grid

        real(dp) :: Cmn_t_det(4), Cmn_t_mix(4), Lmn_t(4), Cmn_t_stoc(4), Cmn_t(4)
        real(dp) :: dt, time

        complex(dp), allocatable  :: DS_0(:,:,:,:,:), SD_t(:,:,:,:,:), DS(:,:,:,:), DD_prefactor(:,:,:,:)
        complex(dp), allocatable  :: DFS_0(:,:,:,:,:), FSD_t(:,:,:,:,:)
        integer, allocatable :: my_dd_b(:), my_dd_a(:), my_ds_b(:), my_ds_a(:)
        integer :: i,j,dir, k, s, kp, it, my_count_det2, my_count_detstoc, a,b,l, m
        complex(dp) :: tmp1, tmp2

        if(sim%all_PAW%N_PAW_atoms.gt.0.and.sim%system%n_stochastic.gt.0) then
            print *, 'Stochastic/Mixed DFT KG not compatible with PAW pseudopotentials, use HGH only'
        endif

        call initialize_Kubo_Greenwood(sim)

        if(sim%parallel%myid.eq.0) print *, 'Start Kubo Greenwood: ', 'N_deterministic: ', sim%system%n_deterministic, &
            'N_stochastic: ', sim%system%n_stochastic
        
        my_count_det2=0
        j=0
        i=sim%parallel%myid_k+1
        allocate(my_dd_b(sim%system%n_deterministic**2/sim%parallel%nproc_k+1))
        allocate(my_dd_a(sim%system%n_deterministic**2/sim%parallel%nproc_k+1))
        do a=1, sim%system%n_deterministic; do b=1,sim%system%n_deterministic
            j=j+1
            if(j.eq.i) then
                my_count_det2=my_count_det2+1
                i=i+sim%parallel%nproc_k
                my_dd_a(my_count_det2)=a
                my_dd_b(my_count_det2)=b
            endif
        enddo;enddo

        my_count_detstoc=0
        j=0
        i=sim%parallel%myid_k+1
        allocate(my_ds_b(sim%system%n_deterministic*sim%system%n_stochastic/sim%parallel%nproc_k+1))
        allocate(my_ds_a(sim%system%n_deterministic*sim%system%n_stochastic/sim%parallel%nproc_k+1))
        do a=1, sim%system%n_deterministic; do b=1,sim%system%n_stochastic
            j=j+1
            if(j.eq.i) then
                my_count_detstoc=my_count_detstoc+1
                i=i+sim%parallel%nproc_k
                my_ds_a(my_count_detstoc)=a
                my_ds_b(my_count_detstoc)=b
            endif
        enddo;enddo

        !Prefactor & initialization  Phase=========================================================================
        !Deterministic-Deterministic block
        if(sim%system%n_deterministic.gt.0) then
            allocate(sim%DPM(sim%system%n_deterministic, sim%system%n_deterministic, &
                    3, sim%system%n_kpoints/sim%parallel%n_k_groups, &
                    sim%system%n_spin/sim%parallel%n_spin_groups))  
            allocate(DD_prefactor(my_count_det2, 4, sim%system%n_kpoints/sim%parallel%n_k_groups, &
                    sim%system%n_spin/sim%parallel%n_spin_groups))
            DD_prefactor=0.0_dp    
            sim%DPM=0.0_dp
            do i=1, size(sim%orbitals,3);  do j=1, size(sim%orbitals,2) 
                kp=sim%orbitals(1,j,i)%k_point
                s=sim%orbitals(1,j,i)%spin

                ks_orbitals(1:sim%system%n_deterministic/sim%parallel%n_band_groups) => &
                sim%orbitals(sim%system%deterministic_start:sim%system%deterministic_end,j,i)
                call Dipole_Matrix_Elements(ks_orbitals, ks_orbitals, &
                sim%grids, sim%parallel, sim%potentials, sim%atoms, sim%elements, sim%all_PAW, sim%DPM(:,:,:,j,i))
                do k=1, my_count_det2 
                    if(my_dd_a(k).eq.my_dd_b(k)) cycle
                    tmp1= 2.0_dp*sim%all_occ(my_dd_a(k),kp,s)*sim%all_weight(kp,s) &
                    *sum(sim%DPM(my_dd_a(k),my_dd_b(k),:,j,i)*sim%DPM(my_dd_b(k),my_dd_a(k),:,j,i))
                    tmp2=0.5_dp*(sim%all_eigs(my_dd_a(k),kp,s)+sim%all_eigs(my_dd_b(k),kp,s))- &
                    (sim%system%chemical_potential(i)-sim%entropy/sum(sim%system%nelec))

                    DD_prefactor(k,1,j,i)=DD_prefactor(k,1,j,i)+ tmp1
                    DD_prefactor(k,2,j,i)=DD_prefactor(k,2,j,i)- tmp1*tmp2
                    DD_prefactor(k,3,j,i)=DD_prefactor(k,3,j,i)- tmp1*tmp2
                    DD_prefactor(k,4,j,i)=DD_prefactor(k,4,j,i)+ tmp1*tmp2**2
                enddo
            enddo; enddo
        endif

        if(sim%system%n_stochastic.gt.0) then
            if(sim%parallel%myid.eq.0.and.sim%KG%filter_stoc) &
            print *, "Initialize stochastic Kubo Greenwood with (1-FD)"
            call pure_Stochastic_Kubo_Greenwood_initialize(sim%orbitals, sim%stochastic, sim%KG,  &
            sim%system, sim%potentials, sim%atoms, sim%elements, sim%parallel, sim%all_PAW, sim%grids,&
            J_orbitals, J_orbitals_t, sim%system%chemical_potential(:)-sim%entropy/sum(sim%system%nelec))

            if(sim%system%n_deterministic.gt.0) then
                if(sim%parallel%myid.eq.0) print *, "Initialize Mixed stochastic Kubo Greenwood"
                call Mixed_Kubo_Greenwood_initialize(sim%orbitals, sim%stochastic, sim%KG, sim%system, &
                sim%potentials, sim%atoms, sim%elements, sim%all_PAW, &
                sim%parallel, sim%grids, unfiltered_stoc_orbitals, sim%all_filter)  

                allocate(DS(sim%system%n_deterministic, sim%system%n_stochastic,3,2))
                allocate(DS_0(my_count_detstoc,3,2,sim%system%n_kpoints_local,sim%system%n_spin_local))
                allocate(SD_t(my_count_detstoc,3,2,sim%system%n_kpoints_local,sim%system%n_spin_local))

                if(sim%KG%project_out_J) then
                    allocate(DFS_0(my_count_detstoc,3,2,sim%system%n_kpoints_local,sim%system%n_spin_local))
                    allocate(FSD_t(my_count_detstoc,3,2,sim%system%n_kpoints_local,sim%system%n_spin_local))
                endif

                do i=1, size(sim%orbitals,3);  do j=1, size(sim%orbitals,2) 
                    kp=sim%orbitals(1,j,i)%k_point
                    s=sim%orbitals(1,j,i)%spin

                    ks_orbitals(1:sim%system%n_deterministic/sim%parallel%n_band_groups) => &
                    sim%orbitals(sim%system%deterministic_start:sim%system%deterministic_end,j,i)

                    call Dipole_Matrix_Elements(ks_orbitals, unfiltered_stoc_orbitals(:,j,i), &
                    sim%grids, sim%parallel, sim%potentials, sim%atoms, sim%elements, sim%all_PAW, DS(:,:,:,1), &
                     DS(:,:,:,2), sim%system%chemical_potential(i)-sim%entropy/sum(sim%system%nelec))

                    do dir=1,3; do k=1, my_count_detstoc
                        DS_0(k,dir,:,j,i)=  DS(my_ds_a(k),my_ds_b(k),dir,:)*sim%stoc_norm(my_ds_b(k),kp,s) &
                                *2.0_dp*sim%all_occ(my_ds_a(k),kp,s)*sim%all_weight(kp,s)
                    enddo; enddo

                    if(sim%KG%project_out_J) then
                        do l=1,2; do dir=1,3
                            call Orbitals_Project_out(ks_orbitals(:), J_orbitals(:,dir,l,j,i), &
                            sim%grids, sim%parallel,sim%atoms, sim%elements, sim%all_PAW, sim%all_filter(:,kp,s))
                        enddo; enddo 

                        stoc_orbitals(1:sim%system%n_stochastic/sim%parallel%n_band_groups)=> &
                            sim%orbitals(sim%system%stochastic_start:sim%system%stochastic_end, j, i)

                            call Dipole_Matrix_Elements(ks_orbitals, stoc_orbitals(:), &
                                sim%grids, sim%parallel, sim%potentials, sim%atoms, sim%elements, sim%all_PAW, DS(:,:,:,1), &
                                DS(:,:,:,2), sim%system%chemical_potential(i)-sim%entropy/sum(sim%system%nelec))

                        do dir=1,3; do k=1, my_count_detstoc
                            DFS_0(k,dir,:,j,i)=  DS(my_ds_a(k),my_ds_b(k),dir,:)*sim%stoc_occ(my_ds_b(k),kp,s) &
                                    *2.0_dp*sim%all_weight(kp,s)
                        enddo; enddo
                    endif
                enddo;enddo
            endif
        endif

        !Start time Dependent phase
        call  SIL_Initialize(sim%KG%Q, sim%orbitals, sim%system, sim%parallel, sim%grids, sim%KG%SIL_Rank)

        dt=min(&
            pi/(sim%KG%nw*sim%KG%dw), &
            pi/(1.10_dp*sim%grids(2)%Ecut-minval(sim%Minimum_eigenvalue(:,:))-sim%grids(2)%Ecut*0.1_dp))

        if(sim%parallel%myid.eq.0) print *, &
            'dt:', dt, 'w_max:', sim%KG%nw*sim%KG%dw
        if(sim%parallel%myid.eq.0)   print *, 'dt(w_max):' , pi/(sim%KG%nw*sim%KG%dw), &
            'dt(Ecut):', pi/(1.10_dp*sim%grids(2)%Ecut-minval(sim%Minimum_eigenvalue(:,:))-sim%grids(2)%Ecut*0.1_dp), 'DeltaE:', &
             (1.10_dp*sim%grids(2)%Ecut-minval(sim%Minimum_eigenvalue(:,:))-sim%grids(2)%Ecut*0.1_dp)
        if(sim%parallel%myid.eq.0) print *, 'nt:', int(5.25/sim%KG%delta/dt)
        if(sim%parallel%myid.eq.0) print *, 'SIL_Rank: ', sim%KG%SIL_Rank, sim%KG%err_allowed
        if(sim%parallel%myid.eq.0) open(newunit=sim%KG%output_unit,  file="Kubo_Greenwood_Cmn_t.data", status="unknown", recl=1024)
        if(sim%parallel%myid.eq.0) write(sim%KG%output_unit,*) 'time ', 'Re_C11(t) ', 'Re_C12(t) ' , 'Re_C21(t) ', 'Re_C22(t) ', &
                                                                       'Im_C11(t) ', 'Im_C12(t) ', 'Im_C21(t) ', 'Im_C22(t) '

        time=0.0_dp
        Lmn_t=0.0_dp
        Cmn_t_det=0.0_dp
        Cmn_t_stoc=0.0_dp
        Cmn_t_mix=0.0_dp
        Cmn_t=0.0_dp

        do it=1, int(5.25_dp/sim%KG%delta/dt)
            if(sim%parallel%myid.eq.0) print *, it, "/", int(5.25_dp/sim%KG%delta/dt),  &
            (1.0_dp*it)/int(5.25_dp/sim%KG%delta/dt)*100, '%'

            Lmn_t=Lmn_t - dt*Cmn_t*0.5_dp*exp(-0.25_dp*sim%KG%delta**2*time**2)
            time=dt*(it-1)
            Cmn_t=0.0_dp
            Cmn_t_det=0.0_dp
            Cmn_t_stoc=0.0_dp
            Cmn_t_mix=0.0_dp
            do i=1, size(sim%orbitals,3);  do j=1, size(sim%orbitals,2) 
                kp=sim%orbitals(1,j,i)%k_point
                s=sim%orbitals(1,j,i)%spin

                if(sim%system%n_deterministic.gt.0) then
                    ks_orbitals(1:sim%system%n_deterministic/sim%parallel%n_band_groups) => &
                    sim%orbitals(sim%system%deterministic_start:sim%system%deterministic_end,j,i)
                endif
                if(sim%system%n_stochastic.gt.0) then
                    stoc_orbitals(1:sim%system%n_stochastic/sim%parallel%n_band_groups)=> &
                        sim%orbitals(sim%system%stochastic_start:sim%system%stochastic_end, j, i)
                endif

                !Update stochastic orbitals
                if(sim%system%n_stochastic.gt.0) then
                    if(it.gt.1) then
                        do k = 1, size(stoc_orbitals);
                            !Main Forward Step
                            call SIL_Propagate(stoc_orbitals(k), sim%KG%Q(:,:,j,i), sim%KG%SIL_rank, sim%grids,  &
                                sim%potentials, sim%atoms, sim%elements, &
                                sim%parallel, dt, sim%KG%err_allowed)
                            if(sim%system%n_deterministic.gt.0) then
                                call SIL_Propagate(unfiltered_stoc_orbitals(k,j,i), sim%KG%Q(:,:,j,i), sim%KG%SIL_rank, &
                                    sim%grids,sim%potentials, sim%atoms, sim%elements, &
                                    sim%parallel, dt, sim%KG%err_allowed)
                            endif
                        enddo
                        do l=1,2; do dir=1,3; do k = 1, size(stoc_orbitals)
                            call SIL_Propagate(J_orbitals(k,dir,l,j,i), sim%KG%Q(:,:,j,i), sim%KG%SIL_rank, sim%grids,  &
                                sim%potentials, sim%atoms, sim%elements, &
                                sim%parallel, dt, sim%KG%err_allowed)
                        enddo; enddo; enddo
                    endif
                endif

                !Deterministic TD-phases
                if(sim%system%n_deterministic.gt.0) then
                    do l=1,4
                        Cmn_t_det(l)=Cmn_t_det(l)+aimag(sum(DD_prefactor(:,l,j,i)* &
                             exp(-i_*(sim%all_eigs(my_dd_a(:my_count_det2),kp,s)&
                                    -sim%all_eigs(my_dd_b(:my_count_det2),kp,s))*time)))
                    enddo
                    !Deterministic-Stochastic Mixed
                    if(sim%system%n_stochastic.gt.0) then

                        call Dipole_Matrix_Elements(ks_orbitals, unfiltered_stoc_orbitals(:,j,i), &
                        sim%grids, sim%parallel, sim%potentials, sim%atoms, sim%elements, sim%all_PAW, DS(:,:,:,1), &
                         DS(:,:,:,2), sim%system%chemical_potential(i)-sim%entropy/sum(sim%system%nelec))

                        do dir=1,3; do k=1, my_count_detstoc
                            SD_t(k,dir,:,j,i)= conjg(DS(my_ds_a(k),my_ds_b(k),dir,:))*exp(-i_*sim%all_eigs(my_ds_a(k),kp,s)*time)
                        enddo; enddo
                        do dir=1,3
                            Cmn_t_mix(1)=Cmn_t_mix(1)+aimag(sum(DS_0(:,dir,1,j,i)*SD_t(:,dir,1,j,i))) 
                            Cmn_t_mix(2)=Cmn_t_mix(2)-aimag(sum(DS_0(:,dir,1,j,i)*SD_t(:,dir,2,j,i)))  
                            Cmn_t_mix(3)=Cmn_t_mix(3)-aimag(sum(DS_0(:,dir,2,j,i)*SD_t(:,dir,1,j,i)))  
                            Cmn_t_mix(4)=Cmn_t_mix(4)+aimag(sum(DS_0(:,dir,2,j,i)*SD_t(:,dir,2,j,i))) 
                        enddo

                        if(sim%KG%project_out_J) then
                            call J_Matrix_Elements(ks_orbitals, stoc_orbitals(:), &
                            sim%grids, sim%parallel, sim%potentials, sim%atoms, sim%elements, sim%all_PAW, DS(:,:,:,1), &
                            DS(:,:,:,2), sim%system%chemical_potential(i)-sim%entropy/sum(sim%system%nelec), .false.)

                            do dir=1,3; do k=1, my_count_detstoc
                                FSD_t(k,dir,:,j,i)= -conjg(DS(my_ds_a(k),my_ds_b(k),dir,:))&
                                                        *exp(-i_*sim%all_eigs(my_ds_a(k),kp,s)*time)
                            enddo; enddo
                            do dir=1,3
                                Cmn_t_mix(1)=Cmn_t_mix(1)+aimag(sum(DFS_0(:,dir,1,j,i)*FSD_t(:,dir,1,j,i))) 
                                Cmn_t_mix(2)=Cmn_t_mix(2)-aimag(sum(DFS_0(:,dir,1,j,i)*FSD_t(:,dir,2,j,i)))  
                                Cmn_t_mix(3)=Cmn_t_mix(3)-aimag(sum(DFS_0(:,dir,2,j,i)*FSD_t(:,dir,1,j,i)))  
                                Cmn_t_mix(4)=Cmn_t_mix(4)+aimag(sum(DFS_0(:,dir,2,j,i)*FSD_t(:,dir,2,j,i))) 
                            enddo
                        endif

                    endif

                endif

                if(sim%system%n_stochastic.gt.0) then
                    stoc_orbitals(1:sim%system%n_stochastic/sim%parallel%n_band_groups)=> &
                    sim%orbitals(sim%system%stochastic_start:sim%system%stochastic_end, j, i)
                    grid=>sim%grids(stoc_orbitals(1)%of(1)%grid)

                    
                    call orbitals_both_current(stoc_orbitals, sim%potentials, sim%atoms, sim%elements, &
                                               sim%parallel,sim% all_PAW, sim%grids, &
                    J_orbitals_t(:,:,:,j,i), sim%system%chemical_potential(i)-sim%entropy/sum(sim%system%nelec))
                     
                     do k = 1, size(stoc_orbitals)
                     do dir=1,3; do s=1, size(stoc_orbitals(k)%of)                    
                        Cmn_t_stoc(1) = Cmn_t_stoc(1) +&
                            2.0_dp*aimag(integrate_3D_G(conjg(J_orbitals(k,dir,1,j,i)%of(s)%G)*J_orbitals_t(k,dir,1,j,i)%of(s)%G,&
                            grid, sim%parallel))*stoc_orbitals(k)%occ(s)*stoc_orbitals(k)%weight(s)
                        Cmn_t_stoc(2) = Cmn_t_stoc(2) -&
                            2.0_dp*aimag(integrate_3D_G(conjg(J_orbitals(k,dir,1,j,i)%of(s)%G)*J_orbitals_t(k,dir,2,j,i)%of(s)%G,&
                            grid, sim%parallel))*stoc_orbitals(k)%occ(s)*stoc_orbitals(k)%weight(s)
                        Cmn_t_stoc(3) = Cmn_t_stoc(3) -&
                            2.0_dp*aimag(integrate_3D_G(conjg(J_orbitals(k,dir,2,j,i)%of(s)%G)*J_orbitals_t(k,dir,1,j,i)%of(s)%G,&
                            grid, sim%parallel))*stoc_orbitals(k)%occ(s)*stoc_orbitals(k)%weight(s)
                        Cmn_t_stoc(4) = Cmn_t_stoc(4) + &
                            2.0_dp*aimag(integrate_3D_G(conjg(J_orbitals(k,dir,2,j,i)%of(s)%G)*J_orbitals_t(k,dir,2,j,i)%of(s)%G,&
                            grid, sim%parallel))*stoc_orbitals(k)%occ(s)*stoc_orbitals(k)%weight(s)
                        enddo; enddo
                    enddo
                endif

            enddo; enddo

            Cmn_t=Cmn_t_det+Cmn_t_mix
            call parallel_task('sum', Cmn_t, sim%parallel, 'all')
            call parallel_task('sum', Cmn_t_stoc, sim%parallel, 'space')
            Cmn_t=Cmn_t+Cmn_t_stoc

            Cmn_t(:)=Cmn_t(:)*1.0_dp/3.0_dp/product(sim%grids(1)%Box_Length(:))
            if(sim%parallel%myid.eq.0) write(sim%KG%output_unit, FMT_real) time, Cmn_t(:), Lmn_t(:)
            if(sim%parallel%myid.eq.0) flush(sim%KG%output_unit)                                                           
            Lmn_t=Lmn_t(:) -dt*Cmn_t*0.5_dp*exp(-0.25_dp*sim%KG%delta**2*time**2)                       
        enddo

        if(allocated(unfiltered_stoc_orbitals)) then
            do i = 1, size(unfiltered_stoc_orbitals,3)
            do j = 1, size(unfiltered_stoc_orbitals,2);do k = 1, size(unfiltered_stoc_orbitals,1)
                call deallocate_orbital(unfiltered_stoc_orbitals(k,j,i),sim%grids)
            enddo;enddo;enddo
            deallocate(unfiltered_stoc_orbitals)
        endif

        if(allocated(J_orbitals)) then
            do i = 1, size(J_orbitals,5);do j = 1, size(J_orbitals,4);do k = 1, size(J_orbitals,3)
                do l = 1, size(J_orbitals,2);do m = 1, size(J_orbitals,1)
                    call deallocate_orbital(J_orbitals(m,l,k,j,i),sim%grids)
                enddo;enddo
            enddo;enddo;enddo
            deallocate(J_orbitals)
        endif

        if(allocated(J_orbitals_t)) then
            do i = 1, size(J_orbitals_t,5);do j = 1, size(J_orbitals_t,4);do k = 1, size(J_orbitals_t,3)
                do l = 1, size(J_orbitals_t,2);do m = 1, size(J_orbitals_t,1)
                    call deallocate_orbital(J_orbitals_t(m,l,k,j,i),sim%grids)
                enddo;enddo
            enddo;enddo;enddo
            deallocate(J_orbitals_t)
        endif

    end subroutine

    subroutine Mixed_Kubo_Greenwood_initialize(orbitals,stochastic, KG,   &
        system, potentials, atoms, elements, all_PAW, parallel, &
        grids, unfiltered_stoc, all_filter)
    use KG_type, only : KG_struct
    use system_type, only : system_struct
    use parallel_type, only : parallel_struct
    use odp_type, only : field_struct, orbital_struct
    use operations_3D, only : integrate_3D_G
    use Local_Ion, only : calculate_local_ion_force
    use odp_type, only: orbital_struct
    use grids_type, only: grid_struct
    use Stochastic_Mod, only : stochastic_FD_filter, stochastic_vector_restore, stochastic_root_OMFD_filter
    use stochastic_type, only : stochastic_struct
    use odp, only: allocate_orbital
    use simulation_type, only : potentials_struct, all_PAW_struct
    use element_type, only : element_struct
    use atom_type, only : atom_struct
    use Eigen_LOBPCG, only : Orbitals_Project_out

    type(parallel_struct), intent(in) :: parallel
    type(grid_struct), intent(inout) :: grids(:)
    type(system_struct), intent(inout) :: system
    type(atom_struct), intent(inout) :: atoms(:)
    type(element_struct), intent(inout) :: elements(:)
    type(all_PAW_struct), intent(inout) :: all_PAW
    type(potentials_struct), intent(in) :: potentials
    type(KG_struct),intent(inout) :: KG
    type(stochastic_struct), intent(inout) :: stochastic(:,:)

    type(orbital_struct), intent(in), target :: orbitals(:,:,:)
    type(orbital_struct), intent(inout), allocatable :: unfiltered_stoc(:,:,:)
    type(orbital_struct), pointer ::  stoc_orbitals(:), ks_orbitals(:)

    real(dp) , intent(in) :: all_filter(:,:,:)

    real(dp), allocatable :: dummy_1(:,:,:), dummy_2(:,:,:)

    integer :: s, i, j, k, kp, dir   

    allocate(unfiltered_stoc( system%n_stochastic/ parallel%n_band_groups, &
                          system%n_kpoints / parallel%n_k_groups, &
                          system%n_spin/ parallel%n_spin_groups))
                          
    do i = 1, size(unfiltered_stoc,3);do j = 1, size(unfiltered_stoc,2)
    
        stoc_orbitals(1:system%n_stochastic/parallel%n_band_groups)=> &
        orbitals(system%stochastic_start:system%stochastic_end, j, i)

        do k = 1, size(unfiltered_stoc,1)
            unfiltered_stoc(k,j,i)%band=stoc_orbitals(k)%band
            unfiltered_stoc(k,j,i)%k_point=stoc_orbitals(k)%k_point
            unfiltered_stoc(k,j,i)%spin=stoc_orbitals(k)%spin
            unfiltered_stoc(k,j,i)%degeneracy=stoc_orbitals(k)%degeneracy
            unfiltered_stoc(k,j,i)%n_spinor=stoc_orbitals(k)%n_spinor
            call allocate_orbital(unfiltered_stoc(k,j,i), grids, unfiltered_stoc(k,j,i)%k_point+2, parallel)
            unfiltered_stoc(k,j,i)%occ=stoc_orbitals(k)%occ
            unfiltered_stoc(k,j,i)%weight=stoc_orbitals(k)%weight
            unfiltered_stoc(k,j,i)%filter=stoc_orbitals(k)%filter
            unfiltered_stoc(k,j,i)%type=stoc_orbitals(k)%type
            if(allocated(unfiltered_stoc(k,j,i)%seed_array)) deallocate(unfiltered_stoc(k,j,i)%seed_array)
            allocate(unfiltered_stoc(k,j,i)%seed_array(size(stoc_orbitals(k)%seed_array)))
            unfiltered_stoc(k,j,i)%seed_array=stoc_orbitals(k)%seed_array
        enddo

        kp=stoc_orbitals(1)%k_point
        s=stoc_orbitals(1)%spin

        !unfiltered_stoc restored to unfiltered state, and orthogonalized to KS

        call stochastic_vector_restore(unfiltered_stoc(:,j,i), system, grids, parallel, &
            dummy_1, dummy_2, all_PAW, atoms, elements)
            ks_orbitals(1:system%n_deterministic/parallel%n_band_groups) => &
            orbitals(system%deterministic_start:system%deterministic_end,j,i)
        call Orbitals_Project_out(ks_orbitals(:), unfiltered_stoc(:,j,i), &
                grids, parallel, atoms, elements, all_PAW, all_filter(:,kp,s))  
                
        if(KG%filter_stoc) then
            do dir=1,3
                call stochastic_root_OMFD_filter(unfiltered_stoc(:,j,i), stochastic(j,i), &
                system, grids, potentials, atoms, elements, parallel, all_PAW)
            enddo
        endif
    enddo; enddo

end subroutine


subroutine pure_Stochastic_Kubo_Greenwood_initialize(orbitals, stochastic, KG,  &
        system, potentials, atoms, elements, parallel, all_PAW, grids, J_orbitals, J_orbitals_t, he)
    use KG_type, only : KG_struct
    use system_type, only : system_struct
    use parallel_type, only : parallel_struct
    use odp_type, only : field_struct, orbital_struct
    use operations_3D, only : integrate_3D_G
    use Local_Ion, only : calculate_local_ion_force
    use odp_type, only : orbital_struct
    use grids_type, only : grid_struct
    use Stochastic_Mod, only :  stochastic_vector_restore, stochastic_OMFD_filter
    use stochastic_type, only : stochastic_struct
    use simulation_type, only : potentials_struct, all_PAW_struct
    use element_type, only : element_struct, PAW
    use atom_type, only : atom_struct
    use Eigen_LOBPCG, only : Orbitals_Project_out
    use constants, only : pi, i_
    use Non_Local_ion, only : Apply_projectors_RL, Calculate_Projector_overlaps, Apply_S_inverse
    use fft, only: recip_to_real, real_to_recip
    use Apply_Hamiltonian, only: Apply_H, Apply_SinvH
    use odp, only: allocate_field, deallocate_field, allocate_orbital, deallocate_orbital
    use Current,only: J_to_orbitals, orbitals_both_current


    type(parallel_struct), intent(in) :: parallel
    type(grid_struct), intent(inout) :: grids(:)
    type(system_struct), intent(inout) :: system
    type(orbital_struct), intent(in), target :: orbitals(:,:,:)
    type(orbital_struct), intent(inout), allocatable :: J_orbitals(:,:,:,:,:), J_orbitals_t(:,:,:,:,:)
    type(orbital_struct), pointer ::  stoc_orbitals(:)
    type(KG_struct),intent(inout) :: KG
    type(potentials_struct), intent(in) :: potentials
    type(atom_struct), intent(inout) :: atoms(:)
    type(element_struct), intent(inout) :: elements(:)
    type(stochastic_struct), intent(inout) :: stochastic(:,:)
    type(all_PAW_struct), intent(inout) :: all_PAW

    type(field_struct), allocatable :: Sinv_J(:)

    real(dp) , intent(in) ::  he(:)

    integer :: i, j, k, dir,l,s
    
    if(parallel%myid.eq.0) then
        open(newunit=KG%output_unit_time,  file="Kubo_Lmn_time.data", status="unknown", recl=1024)
        write(KG%output_unit_time, *) "1:Time", '2-4:R(a.u.)', '5-7:P(a.u.)[3]', '8:SP(t)(a.u.)', '9:<SP> (a.u)', &
                                        '10:rho(a.u.)', '11:rho/<rho>', '12:F(a.u.)' 
    endif    

    allocate(J_orbitals( system%n_stochastic/ parallel%n_band_groups,3,2, &
                          system%n_kpoints / parallel%n_k_groups, &
                          system%n_spin/ parallel%n_spin_groups))

    allocate(J_orbitals_t( system%n_stochastic/ parallel%n_band_groups,3,2, &
                           system%n_kpoints / parallel%n_k_groups, &
                           system%n_spin/ parallel%n_spin_groups))
                     
    do i = 1, size(J_orbitals,5);do j = 1, size(J_orbitals,4)
        stoc_orbitals(1:system%n_stochastic/parallel%n_band_groups)=> &
        orbitals(system%stochastic_start:system%stochastic_end, j, i)
        do l=1,2; do dir=1,3
            do k = 1, size(J_orbitals,1)
                J_orbitals(k,dir,l,j,i)%band=stoc_orbitals(k)%band
                J_orbitals(k,dir,l,j,i)%k_point=stoc_orbitals(k)%k_point
                J_orbitals(k,dir,l,j,i)%spin=stoc_orbitals(k)%spin
                J_orbitals(k,dir,l,j,i)%degeneracy=stoc_orbitals(k)%degeneracy
                J_orbitals(k,dir,l,j,i)%n_spinor=stoc_orbitals(k)%n_spinor
                call allocate_orbital(J_orbitals(k,dir,l,j,i), grids, J_orbitals(k,dir,l,j,i)%k_point+2, parallel)
                J_orbitals(k,dir,l,j,i)%occ=stoc_orbitals(k)%occ
                J_orbitals(k,dir,l,j,i)%weight=stoc_orbitals(k)%weight
                J_orbitals(k,dir,l,j,i)%filter=stoc_orbitals(k)%filter
                J_orbitals(k,dir,l,j,i)%type=stoc_orbitals(k)%type
                if(allocated(J_orbitals(k,dir,l,j,i)%seed_array)) deallocate(J_orbitals(k,dir,l,j,i)%seed_array)
                allocate(J_orbitals(k,dir,l,j,i)%seed_array(size(stoc_orbitals(k)%seed_array)))
                J_orbitals(k,dir,l,j,i)%seed_array=stoc_orbitals(k)%seed_array      
            enddo
        enddo; enddo

        do l=1,2; do dir=1,3
            do k = 1, size(J_orbitals_t,1)
                J_orbitals_t(k,dir,l,j,i)%band=stoc_orbitals(k)%band
                J_orbitals_t(k,dir,l,j,i)%k_point=stoc_orbitals(k)%k_point
                J_orbitals_t(k,dir,l,j,i)%spin=stoc_orbitals(k)%spin
                J_orbitals_t(k,dir,l,j,i)%degeneracy=stoc_orbitals(k)%degeneracy
                J_orbitals_t(k,dir,l,j,i)%n_spinor=stoc_orbitals(k)%n_spinor
                call allocate_orbital(J_orbitals_t(k,dir,l,j,i), grids, J_orbitals_t(k,dir,l,j,i)%k_point+2, parallel)
                J_orbitals_t(k,dir,l,j,i)%occ=stoc_orbitals(k)%occ
                J_orbitals_t(k,dir,l,j,i)%weight=stoc_orbitals(k)%weight
                J_orbitals_t(k,dir,l,j,i)%filter=stoc_orbitals(k)%filter
                J_orbitals_t(k,dir,l,j,i)%type=stoc_orbitals(k)%type
                if(allocated(J_orbitals_t(k,dir,l,j,i)%seed_array)) deallocate(J_orbitals_t(k,dir,l,j,i)%seed_array)
                allocate(J_orbitals_t(k,dir,l,j,i)%seed_array(size(stoc_orbitals(k)%seed_array)))
                J_orbitals_t(k,dir,l,j,i)%seed_array=stoc_orbitals(k)%seed_array      
            enddo
        enddo; enddo

        call orbitals_both_current(stoc_orbitals, potentials, atoms, elements, parallel, all_PAW, grids, &
                                    J_orbitals(:,:,:,j,i), he(i))
                                           
        if(any(elements%PP_type.eq.PAW)) then
            allocate(Sinv_J(size(J_orbitals(1,1,1,j,i)%of(:))))
            Sinv_J(:)%grid=J_orbitals(1,1,1,j,i)%of(:)%grid

            do s=1, size(Sinv_J)
                call allocate_field(Sinv_J(s), grids(Sinv_J(s)%grid), parallel)
            enddo

            do l=1,2; do dir=1,3
                do k = 1, size(J_orbitals,1)
                    call Apply_S_inverse(J_orbitals(k,dir,l,j,i)%of, atoms, elements, grids, &
                    parallel, Sinv_J, all_PAW, with_CG=.true.)
                    do s=1, size(Sinv_J)
                        J_orbitals(k,dir,l,j,i)%of(s)%R=Sinv_J(s)%R
                        J_orbitals(k,dir,l,j,i)%of(s)%G=Sinv_J(s)%G
                    enddo
                enddo
            enddo;enddo

            do s=1, size(Sinv_J)
                call deallocate_field(Sinv_J(s),grids(Sinv_J(s)%grid))
            enddo
            deallocate(Sinv_J)
        endif

        if(KG%filter_stoc) then
            do dir=1,3
              call stochastic_OMFD_filter(J_orbitals(:,dir,1,j,i), stochastic(j,i), &
              system, grids, potentials, atoms, elements, parallel, all_PAW)
              call stochastic_OMFD_filter(J_orbitals(:,dir,2,j,i), stochastic(j,i), &
              system, grids, potentials, atoms, elements, parallel, all_PAW)
            enddo
        endif
          
    enddo; enddo

end subroutine


end module