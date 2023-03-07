module SCF  

use types, only : dp
    use constants, only : pi
    implicit none
    
    
    public :: Self_consistent_field

    contains

    subroutine Self_consistent_field(sim, grids)

        use simulation_type, only : simulation_struct
        use parallel_mod, only : parallel_wait, parallel_task
        use grids_type, only : grid_struct

        use grids_mod, only : Gspace_grid_to_grid_transfer
        use fft, only: real_to_recip, recip_to_real
        use Eigen_LOBPCG, only : LOBPCG, Orbitals_Project_out, Orbitals_Orthogonalize
        use Eigen_ChebFilter, only : Eig_ChebyFilter
        use Eigen_RMM_DIIS, only: Eig_RMM_DIIS
        use FFT2Matrix, only : orbitals_FFT2Matrix, orbitals_Matrix2FFT
        use Chemical_Potential, only : Find_Chemical_Potential
        use Local_ion, only : Calc_Local_PP_Energy,Calc_Core_Energy,core_pot_0
        use xc_mod, only : calc_XC
        use Hartree, only : calc_Hartree_potential, Hartree_Energy
        use Density, only : calc_density, calc_Compensation_Charge_Density
        use Kinetic, only : Kinetic_Energy
        use Stress_Pressure, only : Stress_and_Pressure
        use Entropy, only: Calc_Entropy
        use Stochastic_Mod, only : stochastic_rootFD_filter, stochastic_vector_restore
        use Ewald_mod, only : Ewald_Calc
        use Update_Potentials_and_Energies, only : Calculate_all_static_Energies
        use odp_type, only : orbital_struct, field_struct
        use odp, only : allocate_field, deallocate_field
        use element_type, only : PAW
        use Non_Local_ion, only : Apply_S_power
        use Non_Local_ion_EnFoSt, only : NL_PAW_Density_matrix
        use PAW_Pseudopotentials, only : Calc_Dij
        use initialization_finalization, only : deterministic_vector_initializer, init_den_pot_orb
        use operations_3D, only : integrate_3D_R
        use grids_mod, only : allocate_local_fields_G
        logical :: test_S=.false.


        type(simulation_struct), intent(inout), target :: sim
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid_density, grid_orbitals
        type(orbital_struct), pointer :: ks_orbitals_lowest(:), ks_orbitals(:), &
                                         stoc_orbitals_full(:), stoc_orbitals(:)
        complex(dp), allocatable, target :: xlast(:,:,:,:), xbeforelast(:,:,:,:), flast(:,:,:,:), &
                                            Ri(:,:,:,:,:), Fi(:,:,:,:,:), x0(:,:,:,:)
        real(dp) :: kf_2(sim%density%n_s), q21, q2max, q2min, ratio

        integer :: i, s, k, k_loc, s_loc, iter, comm_save, myid_save, nproc_save, n, stoc_iter, check_start, check_end, &
                    comm_save_2, myid_save_2, nproc_save_2,s2, ij, ij_all, at, ip, jp,k2

        real(dp) :: delta_eigs(sim%system%n_deterministic,sim%system%n_kpoints,sim%system%n_spin), delta_eigs_max
        real(dp) ::  last_eigs(sim%system%n_deterministic,sim%system%n_kpoints,sim%system%n_spin), mu, ne_check
        real(dp), allocatable :: last_mu(:), before_last_mu(:), flast_mu(:), min_eigs(:), inv_diel(:,:,:), mix_metric(:,:,:), &
                                 F_mu(:,:), R_mu(:,:)

        logical, allocatable :: success(:,:)
        
        call core_pot_0(sim%system, sim%elements, sim%grids, sim%all_PAW)

        !grids passed into SCF subroutined, should be sim$grids, passed in expllicitly for target label
        grid_density=> grids(sim%density%of(1)%grid)
        grid_orbitals=> grids(sim%coarse_density%of(1)%grid)

        allocate(success(sim%system%n_kpoints,sim%system%n_spin))

        call allocate_local_fields_G(x0, grid_density, sim%density%n_s)
        call allocate_local_fields_G(xlast, grid_density,sim%density%n_s)
        call allocate_local_fields_G(xbeforelast, grid_density,sim%density%n_s)
        call allocate_local_fields_G(flast ,grid_density, sim%density%n_s)
        call allocate_local_fields_G(Ri, grid_density, sim%numerics%pulay%n, sim%density%n_s)
        call allocate_local_fields_G(Fi, grid_density, sim%numerics%pulay%n, sim%density%n_s)
        call allocate_local_fields_G(inv_diel, grid_density)
        call allocate_local_fields_G(mix_metric, grid_density)

        allocate(last_mu(sim%density%n_s))
        allocate(before_last_mu(sim%density%n_s))
        allocate(flast_mu(sim%density%n_s))
        allocate(F_mu(sim%numerics%pulay%n+1, sim%density%n_s))
        allocate(R_mu(sim%numerics%pulay%n+1, sim%density%n_s))

        do s=1,sim%density%n_s
            xlast(:,:,:,s)=sim%density%of(s)%G !Input to SCF cycle
            last_mu(s)=sim%system%chemical_potential(s)
            if(sim%all_PAW%N_PAW_atoms.gt.0) then
                ij_all=0
                do at=1, size(sim%atoms)
                    if(sim%elements(sim%atoms(at)%element)%PP_type.ne.PAW) cycle
                    do ij=1, sim%elements(sim%atoms(at)%element)%PAW%tab%lmn2_size
                        ip=sim%elements(sim%atoms(at)%element)%PAW%tab%indklmn(7,ij)
                        jp=sim%elements(sim%atoms(at)%element)%PAW%tab%indklmn(8,ij)
                        ij_all=ij_all+1
                        sim%all_PAW%Rhoij_mix_m1(ij_all,s)=real(sim%atoms(at)%PAW%rho_ij(ip,jp,s))
                    enddo
                enddo
            endif
        enddo

        xbeforelast=0.0
        if(sim%numerics%kerker%plasma) then
            kf_2(:)=(4.0_dp/pi)*(3.0_dp*(pi**2.0_dp)*sim%density%n_s*sim%system%nelec(:)&
                /product(grid_density%box_length))**(1.0_dp/3.0_dp)
            sim%numerics%kerker%C=(sum(kf_2(:))/sim%density%n_s)**(-0.5_dp)
        endif
        !Calculate the Approximate inverse dielectric for density mixing 
        inv_diel =  grid_density%cutden * &
                    max(sim%numerics%kerker%A*(1.0_dp/sim%numerics%kerker%B + &
                    sim%numerics%kerker%C**2*grid_density%G2)/(1.0_dp+sim%numerics%kerker%C**2*grid_density%G2), &
                    sim%numerics%kerker%Amin) 

        !Calculate metric for pulay density mixing 
        q2min=minval(grid_density%dG)**2
        q2max=grid_density%Ecut*8.0_dp
        ratio=min(20.0_dp,q2max/q2min)
        q21=q2max*q2min*(ratio-1.0_dp)/(q2max-q2min)          
        where(grid_density%G2.gt.q2min)
            mix_metric=(grid_density%G2+q21)/grid_density%G2
        elsewhere
            mix_metric=(q2min+q21)/q2min
        endwhere

        !Set density and residual differences pointers
        Ri=0.0_dp
        Fi=0.0_dp
        do s=1,sim%density%n_s
            do i=1, sim%numerics%pulay%n
                sim%numerics%pulay%R_i(i,s)%ppointer(1:product(grid_density%Ng_local))=>Ri(:,:,:,i,s)
                sim%numerics%pulay%F_i(i,s)%ppointer(1:product(grid_density%Ng_local))=>Fi(:,:,:,i,s)
            enddo
        enddo
        if(sim%all_PAW%N_PAW_atoms.gt.0) then
            sim%all_PAW%Rhoij_mix_Ri=0.0_dp
            sim%all_PAW%Rhoij_mix_Fi=0.0_dp
            do s=1,sim%density%n_s
                do i=1, sim%numerics%pulay%n
                    sim%all_PAW%Rhoij_mix_R_i(i,s)%ppointer(1:size(sim%all_PAW%Rhoij_mix_Ri,1))=>sim%all_PAW%Rhoij_mix_Ri(:,i,s)
                    sim%all_PAW%Rhoij_mix_F_i(i,s)%ppointer(1:size(sim%all_PAW%Rhoij_mix_Fi,1))=>sim%all_PAW%Rhoij_mix_Fi(:,i,s)
                enddo
            enddo
        endif

        sim%numerics%pulay%L2=10.0_dp*sim%numerics%pulay%L2_eps
        delta_eigs_max=0.0

        if(sim%system%n_deterministic.gt.0) then
            delta_eigs=10.0_dp*sim%numerics%pulay%eig_eps
            last_eigs=0.0_dp
            sim%all_eigs(:,:,:)=0.0_dp
        endif

        iter=0
        stoc_iter=0
        last_mu=-10*sim%numerics%pulay%eig_eps

        if(sim%system%n_deterministic.gt.0) then
            do s_loc=1, size(sim%orbitals(:,:,:),3)
                s=sim%orbitals(1,1,s_loc)%spin
            do k_loc=1, size(sim%orbitals(:,:,:),2)
                k=sim%orbitals(1,k_loc, s_loc)%k_point
                ks_orbitals(1:sim%system%n_deterministic/sim%parallel%n_band_groups) => &
                                    sim%orbitals(sim%system%deterministic_start:sim%system%deterministic_end,k_loc,s_loc)
                !call deterministic_vector_initializer(ks_orbitals(:),sim%grids, sim%parallel)
                !call Orbitals_Orthogonalize(ks_orbitals, sim%grids, sim%parallel, &
                !sim%atoms, sim%elements, sim%all_PAW)
            enddo;enddo
        endif
        do while(sim%numerics%pulay%L2.gt.sim%numerics%pulay%L2_eps .or. delta_eigs_max.gt.sim%numerics%pulay%eig_eps &
                    .or.any((sim%system%chemical_potential(:)-last_mu(:)).gt.sim%numerics%pulay%eig_eps))

            if(sim%parallel%myid.eq.0) print *, 'SCF iteration: ', iter+1     
            last_mu=sim%system%chemical_potential
            if(sim%system%n_deterministic.gt.0) then
                sim%all_occ(:,:,:)=0.0_dp
                sim%all_filter(:,:,:)=0.0_dp
                sim%all_eigs(:,:,:)=0.0_dp
            endif
           
            sim%Minimum_eigenvalue=0.0_dp

            do s_loc=1, size(sim%orbitals(:,:,:),3)
                s=sim%orbitals(1,1,s_loc)%spin
            do k_loc=1, size(sim%orbitals(:,:,:),2)
                k=sim%orbitals(1,k_loc, s_loc)%k_point
                !Calculate New Deterministic orbitals
                if(sim%system%n_find_min.gt.0) then
                    !No deterministic orbitals, just looking for lowest eigenvector, so use sub group
                    if(sim%parallel%my_sub.eq.0) then
                        allocate(min_eigs(sim%parallel%n_band_groups_sub))
                        ks_orbitals_lowest(1:1)=>&
                        sim%orbitals((sim%system%n_orbitals_local+1):(sim%system%n_orbitals_local+1), k_loc,s_loc)

                        comm_save=sim%parallel%comm_k; myid_save=sim%parallel%myid_k; nproc_save=sim%parallel%nproc_k
                        comm_save_2=sim%parallel%comm_diff_b; myid_save_2=sim%parallel%myid_diff_b
                        nproc_save_2=sim%parallel%nproc_diff_b
                        
                        sim%parallel%comm_k=sim%parallel%comm_sub
                        sim%parallel%myid_k=sim%parallel%myid_sub
                        sim%parallel%nproc_k=sim%parallel%nproc_sub
                        sim%parallel%comm_diff_b=sim%parallel%comm_diff_b_sub
                        sim%parallel%myid_diff_b=sim%parallel%myid_diff_b_sub
                        sim%parallel%nproc_diff_b=sim%parallel%nproc_diff_b_sub

                        call LOBPCG(ks_orbitals_lowest(:1), sim%grids, sim%parallel, &
                            sim%numerics, sim%potentials, sim%atoms, sim%elements, &
                            min_eigs(:), sim%all_PAW, &
                            sim%parallel%n_band_groups_sub)

                        sim%parallel%comm_k=comm_save; sim%parallel%myid_k=myid_save; sim%parallel%nproc_k=nproc_save
                        sim%parallel%comm_diff_b=comm_save_2; sim%parallel%myid_diff_b=myid_save_2
                        sim%parallel%nproc_diff_b=nproc_save_2
                        
                        if(sim%parallel%myid_sub.eq.0) print *, 'Lowest Eigenvalues:', min_eigs(:)
                        sim%Minimum_eigenvalue(k,s)=min_eigs(1)
                        ks_orbitals_lowest(1)%eig= min_eigs(sim%parallel%my_band_group+1)
                        ks_orbitals_lowest(1)%occ=0.0_dp
                        ks_orbitals_lowest(1)%filter=0.0_dp
                        ks_orbitals_lowest(1)%weight=0.0_dp
                        deallocate(min_eigs)
                    endif
                    call parallel_task('bcast',  sim%Minimum_eigenvalue, sim%parallel, 'k', root=0)

                else if(iter.lt.sim%numerics%lobpcg%n_init) then
                    ks_orbitals_lowest(1:(sim%system%n_orbitals_deterministic_lowest/sim%parallel%n_band_groups))=>&
                        sim%orbitals(sim%system%deterministic_lowest_start:sim%system%deterministic_lowest_end, k_loc,s_loc)
                    call LOBPCG(ks_orbitals_lowest(:), sim%grids, sim%parallel,  &
                                sim%numerics, sim%potentials, sim%atoms, sim%elements, &
                                sim%all_eigs(:sim%system%n_orbitals_deterministic_lowest,k,s), &
                                sim%all_PAW, sim%parallel%n_band_groups)
                    sim%Minimum_eigenvalue(k,s)=sim%all_eigs(1,k,s)
                else if(sim%tasks%Cheby_Filter) then
                    ks_orbitals_lowest(1:(sim%system%n_orbitals_deterministic_lowest/sim%parallel%n_band_groups))=>&
                        sim%orbitals(sim%system%deterministic_lowest_start:sim%system%deterministic_lowest_end, k_loc,s_loc)
                    call Eig_ChebyFilter(ks_orbitals_lowest(:), sim%grids, sim%parallel,&
                                sim%numerics, sim%potentials, sim%atoms, sim%elements,&
                                sim%all_eigs(:sim%system%n_orbitals_deterministic_lowest,k,s), sim%all_PAW, &
                                sim%parallel%n_band_groups)
                    sim%Minimum_eigenvalue(k,s)=sim%all_eigs(1,k,s)
                else if(sim%tasks%RMM_DIIS) then
                    ks_orbitals_lowest(1:(sim%system%n_orbitals_deterministic_lowest/sim%parallel%n_band_groups))=>&
                        sim%orbitals(sim%system%deterministic_lowest_start:sim%system%deterministic_lowest_end, k_loc,s_loc)
                    call Eig_RMM_DIIS(ks_orbitals_lowest(:), sim%grids, sim%parallel,&
                                sim%numerics%diis, sim%potentials, sim%atoms, sim%elements,&
                                sim%all_eigs(:sim%system%n_orbitals_deterministic_lowest,k,s),&
                                sim%all_PAW, sim%parallel%n_band_groups)
                    sim%Minimum_eigenvalue(k,s)=minval(sim%all_eigs(:,k,s))
                endif

                if(sim%system%n_deterministic.gt.0) then
                    !Apply Smoothing to KS orbitals truncation
                    check_start=sim%system%deterministic_lowest_start -1 + &
                                sim%system%n_orbitals_deterministic_lowest &
                                - sim%system%n_orbitals_buffer_lowest &
                                - sim%system%n_orbitals_smoother_lowest + 1
                    sim%all_filter(:(check_start-1),k,s) = 1.0_dp
                    check_end=check_start + sim%system%n_orbitals_smoother_lowest
                    if(check_end.gt.check_start) then
                        sim%all_filter(check_start:check_end,k,s) = &
                            (exp(sim%all_eigs(check_start,k,s)-sim%all_eigs(check_start:check_end,k,s)) &
                            -exp(sim%all_eigs(check_start,k,s)-sim%all_eigs(check_end,k,s))) &
                            /(1.0_dp-exp(-(sim%all_eigs(check_end,k,s)-sim%all_eigs(check_start,k,s))))
                    endif
                    sim%all_filter((check_end+1):,k,s) = 0.0_dp

                endif
            enddo; enddo
            
            if(test_S.and.iter.eq.100) stop
            if(test_S) iter=iter+1; if(test_S) cycle

            call parallel_task('sum', sim%Minimum_eigenvalue(:,:), sim%parallel, 'diff_k')
            call parallel_task('sum', sim%Minimum_eigenvalue(:,:), sim%parallel, 'diff_s')
            if(sim%system%n_deterministic.gt.0) then
                call parallel_task('sum', sim%all_eigs(:,:,:), sim%parallel, 'diff_k')
                call parallel_task('sum', sim%all_eigs(:,:,:), sim%parallel, 'diff_s')
                call parallel_task('sum', sim%all_filter(:,:,:), sim%parallel, 'diff_k')
                call parallel_task('sum', sim%all_filter(:,:,:), sim%parallel, 'diff_s')
            endif
        

            !Prepare Stochastic orbitals
            do s_loc=1, size(sim%orbitals(:,:,:),3)
                s=sim%orbitals(1,1,s_loc)%spin
                do k_loc=1, size(sim%orbitals(:,:,:),2)
                    k=sim%orbitals(1,k_loc,s_loc)%k_point
                    success(k,s)=.false.
                    if(sim%system%n_stochastic.gt.0) then
                        sim%stochastic(k_loc,s_loc)%Emin=sim%Minimum_eigenvalue(k,s)-sim%grids(2)%Ecut*0.1_dp
                        sim%stochastic(k_loc,s_loc)%Emax=1.10_dp*sim%grids(2)%Ecut                
                        sim%stochastic(k_loc,s_loc)%deltaE= &
                                (sim%stochastic(k_loc,s_loc)%Emax-sim%stochastic(k_loc,s_loc)%Emin)*0.5_dp
                        sim%stochastic(k_loc,s_loc)%Ebar= &
                                (sim%stochastic(k_loc,s_loc)%Emax+sim%stochastic(k_loc,s_loc)%Emin)*0.5_dp
                        if(sim%parallel%myid_k.eq.0) print *, 'ChebyFilter: ', &
                            'Emin: ', sim%stochastic(k_loc,s_loc)%Emin, &
                            'Emax: ', sim%stochastic(k_loc,s_loc)%Emax, &
                            'DeltaE: ', sim%stochastic(k_loc,s_loc)%deltaE, &  
                            'Ebar: ', sim%stochastic(k_loc,s_loc)%Ebar
                    endif
                enddo!end parallel loop on k-points
            enddo !end parallel loop on colinear spin


            sim%all_occ(:,:,:)=0.0_dp

            do s_loc=1, size(sim%orbitals(:,:,:),3); 
                s=sim%orbitals(1,1, s_loc)%spin
                !Calculate Stochastic Orbitals
                success(:,s)=.false.
                do while(.not.all(success(:,s)))
                    if(sim%system%n_stochastic.gt.0) then
                        stoc_iter=stoc_iter+1
                        do k_loc=1, size(sim%orbitals(:,:,:),2)
                            k=sim%orbitals(1,k_loc, s_loc)%k_point

                            ks_orbitals(1:sim%system%n_deterministic/sim%parallel%n_band_groups) => &
                                sim%orbitals(sim%system%deterministic_start:sim%system%deterministic_end,k_loc,s_loc)

                            stoc_orbitals(1:sim%system%n_stochastic/sim%parallel%n_band_groups) => &
                                sim%orbitals(sim%system%stochastic_start:sim%system%stochastic_end,k_loc,s_loc)
                            
                            stoc_orbitals_full(1:sim%system%n_orbitals_stochastic_full/sim%parallel%n_band_groups) => &
                                sim%orbitals(sim%system%stochastic_full_start:sim%system%stochastic_full_end,k_loc,s_loc)

                                
                            if(.not.success(k,s)) then
                                call stochastic_vector_restore(stoc_orbitals_full(:), sim%system, sim%grids, sim%parallel,  &
                                    sim%stoc_norm, sim%stoc_occ, sim%all_PAW, sim%atoms, sim%elements)
                                !Orthogonalize Stochasic to Determinisic
                                if(sim%system%n_deterministic.gt.0) then
                                    call Orbitals_Project_out(ks_orbitals(:), stoc_orbitals(:), &
                                        sim%grids, sim%parallel, sim%atoms, sim%elements, sim%all_PAW, &
                                        sim%all_filter(:,k,s))
                                endif
                                call stochastic_rootFD_filter(stoc_orbitals_full(:), sim%stochastic(k_loc,s_loc), &
                                    sim%system, sim%grids,sim%potentials, sim%atoms, sim%elements, sim%parallel, sim%all_PAW)

                                call stochastic_ne_calc_scf(stoc_orbitals_full(:), sim%stochastic(k_loc,s_loc), &
                                    sim%grids, sim%atoms, sim%elements, sim%parallel, sim%all_PAW)

                                    
                                if(sim%parallel%myid_k.eq.0) print *, 'Ne_Stochastic:', sim%stochastic(k_loc,s_loc)%nelec_calc 
                                if(((sim%stochastic(k_loc,s_loc)%nelec_calc/sim%all_weight(k,s)) &
                                        .gt.(1000*sim%system%nelec(s))) &
                                        .or.(sim%stochastic(k_loc,s_loc)%nelec_calc.lt.0.0_dp)) then

                                        sim%stochastic(k_loc,s_loc)%Emin=sim%stochastic(k_loc,s_loc)%Emin-1.0_dp
                                        sim%stochastic(k_loc,s_loc)%Emax=sim%stochastic(k_loc,s_loc)%Emax+1.0_dp
                                        sim%stochastic(k_loc,s_loc)%deltaE= &
                                        (sim%stochastic(k_loc,s_loc)%Emax-sim%stochastic(k_loc,s_loc)%Emin)*0.5_dp
                                        sim%stochastic(k_loc,s_loc)%Ebar= &
                                        (sim%stochastic(k_loc,s_loc)%Emax+sim%stochastic(k_loc,s_loc)%Emin)*0.5_dp
                                        if(sim%parallel%myid_k.eq.0) print *, ' New ChebyFilter: ', &
                                            'ne_stoc: ',sim%stochastic(k_loc,s_loc)%nelec_calc, &
                                            'mu: ',sim%system%chemical_potential(:), &
                                            'Emin: ', sim%stochastic(k_loc,s_loc)%Emin, &
                                            'Emax: ', sim%stochastic(k_loc,s_loc)%Emax, &
                                            'DeltaE: ', sim%stochastic(k_loc,s_loc)%deltaE, &  
                                            'Ebar: ', sim%stochastic(k_loc,s_loc)%Ebar
                                    success(k,s)=.false.
                                else
                                    success(k,s)=.true.
                                endif
                            endif
                        enddo
                    else
                        success(:,s)=.true.
                    endif
                    if(any(.not.success(:,s))) cycle !back to start of the while

                    call Find_Chemical_Potential( sim%Minimum_eigenvalue(:,s), &
                        sim%all_eigs(:,:,s), sim%all_occ(:,:,s), sim%all_filter(:,:,s), sim%all_degeneracy(:,:,s), & 
                        sim%all_weight(:,s), sim%stochastic(:,s_loc), sim%system%temperature, &
                        sim%system%chemical_potential(s), sim%system, &
                        sim%grids(:), sim%system%nelec(s), sim%parallel, success(:,s))
                    if(any(.not.success(:,s))) cycle !back to start of the while

                    do k_loc=1, size(sim%orbitals(:,:,:),2)
                        if(sim%system%n_deterministic.gt.0) then
                            ks_orbitals(1:sim%system%n_deterministic/sim%parallel%n_band_groups)=>&
                                sim%orbitals(sim%system%deterministic_start:sim%system%deterministic_end, k_loc,s_loc)
                            do n= 1, size(ks_orbitals)
                                k=ks_orbitals(n)%k_point
                                ks_orbitals(n)%eig= sim%all_eigs(ks_orbitals(n)%band,k,s)
                                ks_orbitals(n)%occ= sim%all_occ(ks_orbitals(n)%band,k,s)
                                ks_orbitals(n)%filter= sim%all_filter(ks_orbitals(n)%band,k,s)
                                ks_orbitals(n)%occ=ks_orbitals(n)%occ*ks_orbitals(n)%filter
                            enddo
                        endif
                        if(sim%system%n_stochastic.gt.0) then
                            stoc_orbitals(1:sim%system%n_stochastic/sim%parallel%n_band_groups) => &
                                sim%orbitals(sim%system%stochastic_start:sim%system%stochastic_end,k_loc,s_loc)
                            do n= 1, size(stoc_orbitals)
                                k=stoc_orbitals(n)%k_point
                                stoc_orbitals(n)%eig= 0.0_dp !This kills the eigenvalue to NL force in PAW (which needs to be calculated differently)
                            enddo
                        endif
                    enddo

                    if(sim%parallel%myid_spin.eq.0) print *, 'mu: ', sim%system%chemical_potential(s), &
                                                            ' it worked?=', .not.any(.not.success(:, s)), 's=', s
                    
                enddo !end while loop
                mu=sim%system%chemical_potential(s)
                sim%system%chemical_potential(:)=0.0_dp
                sim%system%chemical_potential(s)=mu
            enddo !end parallel loop on colinear spin

            call parallel_task('sum', sim%all_occ(:,:,:), sim%parallel, 'diff_s')
            call parallel_task('sum', sim%system%chemical_potential(:), sim%parallel, 'diff_s') 

            if(sim%parallel%myid.eq.0) then
                print *, 'Eigs, Occs, filter, k, s'
                do s=1, sim%system%n_spin
                    print *, 'Chemical Potential : ', s, sim%system%chemical_potential(s)!- sim%system%epsatm_all
                    do k=1, sim%system%n_kpoints
                        do i=1, sim%system%n_deterministic
                            print *, sim%all_eigs(i,k,s), &! - sim%system%epsatm_all, &
                            sim%all_occ(i,k,s), &
                            sim%all_filter(i,k,s), 'k=', k , 's=', s
                        enddo
                    enddo
                enddo
            endif
            
            if(sim%system%n_deterministic.gt.0) then
                delta_eigs=sim%all_filter*(sim%all_eigs-last_eigs)*sim%all_occ/sum(sim%all_occ)*size(sim%all_occ)
                last_eigs=sim%all_eigs
                delta_eigs_max=maxval(abs(delta_eigs))
            endif

            if(iter.ge.sim%numerics%pulay%n_skip) then
                if(any(sim%elements(:)%PP_type.eq.PAW)) then
                    call  NL_PAW_Density_matrix(sim%orbitals(:,:,:), sim%atoms, sim%elements, sim%grids, sim%parallel, sim%all_PAW) 
                    call calc_Compensation_Charge_Density(sim%all_PAW, sim%grids, sim%atoms, sim%parallel)               
                endif   
                !Calc Density & new potentials
                call calc_density(sim%coarse_density, sim%density,  &
                                sim%orbitals(:,:,:), sim%grids, sim%parallel, sim%coarse_to_fine, sim%all_PAW)
            endif
            !Density Mixing, w/ PAW rho_ij mixing and chemical potential mixing
            do s=1,sim%density%n_s
                if(iter.ge.sim%numerics%pulay%n_skip) then
                    if(sim%parallel%myid.eq.0) print *, 'Pulay Mix Density & Rhoij if PAW'
                    if(sim%all_PAW%N_PAW_atoms.gt.0) then
                        call pulay_mixing(s, sim%density%of(s)%G, &
                                        xlast(:,:,:,s), xbeforelast(:,:,:,s), flast(:,:,:,s), &
                                        inv_diel, max(sim%numerics%kerker%A, sim%numerics%kerker%Amin), &
                                        mix_metric, sim%system%chemical_potential(s), last_mu(s), before_last_mu(s), flast_mu(s), &
                                        F_mu(:,s), R_mu(:,s), &
                                        sim%numerics%pulay, iter-sim%numerics%pulay%n_skip,  &
                                        sim%parallel, sim%grids, sim%all_PAW, sim%atoms, sim%elements)
                    else
                        call pulay_mixing(s, sim%density%of(s)%G, &
                                        xlast(:,:,:,s), xbeforelast(:,:,:,s), flast(:,:,:,s), &
                                        inv_diel,max(sim%numerics%kerker%A, sim%numerics%kerker%Amin), &
                                        mix_metric, sim%system%chemical_potential(s), last_mu(s), before_last_mu(s), flast_mu(s), &
                                        F_mu(:,s), R_mu(:,s), &
                                        sim%numerics%pulay, iter-sim%numerics%pulay%n_skip, sim%parallel, sim%grids)
                    endif  
                else
                    xlast(:,:,:,s)=sim%density%of(s)%G
                    last_mu(s)=sim%system%chemical_potential(s)
                    if(sim%all_PAW%N_PAW_atoms.gt.0) then
                        ij_all=0
                        do at=1, size(sim%atoms)
                            if(sim%elements(sim%atoms(at)%element)%PP_type.ne.PAW) cycle
                            do ij=1, sim%elements(sim%atoms(at)%element)%PAW%tab%lmn2_size
                                ip=sim%elements(sim%atoms(at)%element)%PAW%tab%indklmn(7,ij)
                                jp=sim%elements(sim%atoms(at)%element)%PAW%tab%indklmn(8,ij)
                                ij_all=ij_all+1
                                sim%all_PAW%Rhoij_mix_m1(ij_all,s)=real(sim%atoms(at)%PAW%rho_ij(ip,jp,s))
                            enddo
                        enddo
                    endif
                endif !skip first steps
            enddo               
        
            if(iter.ge.sim%numerics%pulay%n_skip) then
                !Add Compensation Chareg density for any PAW atoms
                call calc_Hartree_potential(sim%potentials%hartree, sim%density, sim%grids)
                call calc_XC(sim%density, sim%xc, sim%all_PAW, &
                            sim%parallel, sim%grids, sim%energies%xc, sim%potentials%xc)
                    
                do s=1,sim%density%n_s
                    sim%potentials%total_local_fine%of(s)%R =real(sim%potentials%ion_local%of%R) + &
                                                            real(sim%potentials%hartree%of%R) + &
                                                            real(sim%potentials%xc%of(s)%R)

                    if(sim%all_PAW%usepotzero.gt.0)  sim%potentials%total_local_fine%of(s)%R = &
                                                    sim%potentials%total_local_fine%of(s)%R + sum(sim%all_PAW%vpotzero)

                    call real_to_recip(sim%potentials%total_local_fine%of(s), sim%grids)
                    sim%potentials%total_local_fine%of(s)%G=sim%potentials%total_local_fine%of(s)%G*grid_density%cutden
                    call recip_to_real(sim%potentials%total_local_fine%of(s), sim%grids)
                    sim%potentials%total_local_fine%of(s)%R=real(sim%potentials%total_local_fine%of(s)%R)

                    !Transfer the Fine Grid potential to the coarse WF grid
                    call Gspace_grid_to_grid_transfer(grid_density, grid_orbitals, sim%parallel, &
                            sim%potentials%total_local_fine%of(s)%G, &
                            sim%potentials%total_local%of(s)%G, sim%fine_to_coarse)
                    sim%potentials%total_local%of(s)%G=sim%potentials%total_local%of(s)%G*grid_orbitals%cutden
                    sim%potentials%total_local_fine%of(s)%R=real(sim%potentials%total_local_fine%of(s)%R)
                    call recip_to_real(sim%potentials%total_local%of(s), sim%grids)
                    sim%potentials%total_local%of(s)%R=real(sim%potentials%total_local%of(s)%R)
                enddo

                call Calc_Dij(sim%potentials, sim%atoms, sim%elements, sim%all_PAW, sim%grids, sim%parallel)
            endif
            iter=iter+1
            if(iter.eq.sim%numerics%pulay%max_iter) then
                print *, 'Maximum number of SCF iterations reached'
                stop
            endif
        enddo

        if(test_S) stop
    
        do s=1,sim%density%n_s
            do i=1, sim%numerics%pulay%n
                ! nullify(sim%numerics%pulay%R_i(i,s)%ppointer(1:product(grid_density%Ng_local)))
                ! nullify(sim%numerics%pulay%F_i(i,s)%ppointer(1:product(grid_density%Ng_local)))
            enddo
        enddo

        deallocate(xlast)
        deallocate(xbeforelast)
        deallocate(flast)
        deallocate(Ri)
        deallocate(Fi)
        deallocate(flast_mu)
        deallocate(before_last_mu)
        deallocate(last_mu)
        deallocate(R_mu)
        deallocate(F_mu)

        if(sim%parallel%myid.eq.0) then
            print *, 'Final : Eigs, Occs, filter, k, s'
            write(sim%main_output, *) 'Final : Eigs, Occs, filter, k, s'

            do s=1, sim%system%n_spin
                print *, 'Chemical Potential : ', s, sim%system%chemical_potential(s) !- sim%system%epsatm_all
                write(sim%main_output, *)  'Chemical Potential : ', s, sim%system%chemical_potential(s) !- sim%system%epsatm_all

                do k=1, sim%system%n_kpoints
                    do i=1, sim%system%n_deterministic
                    print *, sim%all_eigs(i,k,s),&! - sim%system%epsatm_all, &
                    sim%all_occ(i,k,s), &
                    sim%all_filter(i,k,s), 'k=', k , 's=', s
                    write(sim%main_output, *) sim%all_eigs(i,k,s),&! - sim%system%epsatm_all, &
                    sim%all_occ(i,k,s), &
                    sim%all_filter(i,k,s), 'k=', k , 's=', s
                    enddo
                enddo
            enddo
        endif

        if(.not.sim%tasks%run_Time_Dependent) then
            call Ewald_Calc(sim%atoms, sim%elements, sim%grids, sim%ewald, sim%energies, .false., sim%stress)
            call Calculate_all_static_Energies(sim%orbitals, sim%grids, sim%system, &
                sim%parallel, sim%density, sim%potentials, sim%energies, sim%atoms,  &
                sim%elements, sim%xc, sim%tf, sim%all_PAW, sim%main_output, 0) 
            call Stress_and_Pressure(sim%orbitals, sim%grids, sim%system, sim%parallel, &
            sim%density, sim%potentials, sim%energies, sim%stress, sim%pressures, sim%atoms,  &
            sim%elements, sim%xc, sim%tf, sim%all_PAW, sim%main_output, 0)
            call Calc_Entropy(sim)
        endif

        

    end subroutine
    
    subroutine pulay_mixing(s, x0, xlast, xbeforelast, flast, inv_diel, A, mix_metric, &
                            mu_new, mu_last, mu_beforelast, mu_flast, F_mu, R_mu, &
                            pulay, iter, parallel, grids, all_PAW, atoms, elements)
        use numerics_type, only : pulay_numerics
        use parallel_mod,only : parallel_task
        use parallel_type, only : parallel_struct
        use linalg, only : solve
        use grids_type, only : grid_struct
        use atom_type, only : atom_struct
        use element_type, only : element_struct, PAW
        use simulation_type, only: all_PAW_struct
        use m_paral_atom, only : get_my_natom, get_my_atmtab
        use m_pawrhoij, only : pawrhoij_filter
        use m_paw_denpot, only : pawdenpot
        use Density, only : calc_Compensation_Charge_Density

        ! Finds "x" so that x = g(x)
        ! Implements Algorithm 1. from [1], we use the same notation, except we call
        ! the independent variable "x" instead of "rho":
        ! do i = 1, 2, ...
        !   f_i = g(x_i) - x_i
        !   if modulo(i+1, k) == 0
        !     x_i = x_i + alpha*f_i - (R_i+alpha*F_i)*(F_i^T F_i)^-1 F_i^T f_i
        !   else
        !     x_i = x_i + alpha*f_i
        ! until |f_i| < tol
        !
        ! alpha-> kerker type preconditioner
        ! [1] Banerjee, A. S., Suryanarayana, P., Pask, J. E. (2016). Periodic Pulay
        ! method for robust and efficient convergence acceleration of self-consistent
        ! field iterations. Chemical Physics Letters, 647, 31–35.
        ! https://doi.org/10.1016/j.cplett.2016.01.033
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        type(pulay_numerics), intent(inout) :: pulay
        type(parallel_struct), intent(in) :: parallel
        type(all_PAW_struct), intent(inout), optional :: all_PAW
        type(atom_struct), intent(inout), optional :: atoms(:)
        type(element_struct), intent(inout), optional :: elements(:)


        complex(dp), intent(inout), contiguous, target :: x0(:,:,:), xlast(:,:,:), xbeforelast(:,:,:), &
                                                          flast(:,:,:)
        real(dp), intent(inout) :: mu_new, mu_last, mu_beforelast, mu_flast, F_mu(:), R_mu(:)
        real(dp), intent(in), contiguous, target :: inv_diel(:,:,:), mix_metric(:,:,:)
        real(dp), intent(in) :: A
        complex(dp), pointer :: x_new(:), x_last(:), x_before_last(:), f_last(:)
        real(dp), pointer :: inv_diel_p(:), mix_metric_p(:)
        complex(dp) :: fi(size(x0))
        complex(dp) :: FTf(pulay%n), FF(pulay%n,pulay%n), W(pulay%n)
        real(dp), allocatable :: rhoij_fi(:)
        real(dp) :: mu_fi
        real(dp) :: fnorm, xnorm
        integer, intent(in) :: iter, s
        integer :: i, j, n_rhoij, ij,ij_all, at, e, m, usexcnhat, ip, jp
        real(dp) :: compch_sph,epawdc
        real(dp), allocatable :: nucdipmom(:,:)

        grid=>grids(1) ! Mix the fine density 
        !First Rotate the pointers
        pulay%F_i(pulay%n+1,s)%ppointer(1:size(x0))=>pulay%F_i(1,s)%ppointer(1:size(x0))
        pulay%R_i(pulay%n+1,s)%ppointer(1:size(x0))=>pulay%R_i(1,s)%ppointer(1:size(x0))
        F_mu(pulay%n+1) = F_mu(1)
        R_mu(pulay%n+1) = R_mu(1)
        do i=1, pulay%n 
            pulay%F_i(i,s)%ppointer(1:size(x0))=>pulay%F_i(i+1,s)%ppointer(1:size(x0))
            pulay%R_i(i,s)%ppointer(1:size(x0))=>pulay%R_i(i+1,s)%ppointer(1:size(x0))
            F_mu(i) = F_mu(i+1)
            R_mu(i) = R_mu(i+1)
        enddo
        if(present(all_PAW)) then
            n_rhoij=size(all_PAW%Rhoij_mix_0,1)
            allocate(rhoij_fi(n_rhoij))
            allocate(nucdipmom(3,all_PAW%N_PAW_atoms))
            all_PAW%Rhoij_mix_F_i(pulay%n+1,s)%ppointer(1:n_rhoij)=>all_PAW%Rhoij_mix_F_i(1,s)%ppointer(1:n_rhoij)
            all_PAW%Rhoij_mix_R_i(pulay%n+1,s)%ppointer(1:n_rhoij)=>all_PAW%Rhoij_mix_R_i(1,s)%ppointer(1:n_rhoij)
            do i=1, pulay%n 
                all_PAW%Rhoij_mix_F_i(i,s)%ppointer(1:n_rhoij)=>all_PAW%Rhoij_mix_F_i(i+1,s)%ppointer(1:n_rhoij)
                all_PAW%Rhoij_mix_R_i(i,s)%ppointer(1:n_rhoij)=>all_PAW%Rhoij_mix_R_i(i+1,s)%ppointer(1:n_rhoij)
            enddo
            ij_all=0
            do at=1, size(all_PAW%ij)
                if(elements(atoms(at)%element)%PP_type.ne.PAW) cycle
                do ij=1, elements(atoms(at)%element)%PAW%tab%lmn2_size
                    ip=elements(atoms(at)%element)%PAW%tab%indklmn(7,ij)
                    jp=elements(atoms(at)%element)%PAW%tab%indklmn(8,ij)
                    ij_all=ij_all+1
                    all_PAW%Rhoij_mix_0(ij_all,s)=real(atoms(at)%PAW%rho_ij(ip,jp,s))
                enddo
            enddo
        endif

        x_new(1:size(x0))=>x0(:,:,:)
        x_last(1:size(x0))=>xlast(:,:,:)
        x_before_last(1:size(x0))=>xbeforelast(:,:,:)
        f_last(1:size(x0))=>flast(:,:,:)
        inv_diel_p(1:size(x0))=>inv_diel(:,:,:)
        mix_metric_p(1:size(x0))=>mix_metric(:,:,:)

        if(iter.eq.0) then
             !x_last set at start of SCF cycle
             x_before_last(:)=x_last(:)!0.0_dp
             f_last(:)=0.0_dp
             mu_beforelast=mu_last
             mu_flast=0.0_dp
        endif

        fi(:)=x_new(:)-x_last(:)
        fnorm=sum(real(conjg(fi(:))*fi(:)))
        xnorm=sum(real(conjg(x_last(:))*x_last(:)))
        mu_fi=mu_new-mu_last
        !print *, fnorm
        call parallel_task('sum', fnorm, parallel, 'band')
        call parallel_task('sum', xnorm, parallel, 'band')

        fnorm=sqrt(fnorm*product(grid%box_length(:)))
        xnorm=sqrt(xnorm*product(grid%box_length(:)))
        if(xnorm.lt.1E-12_dp) xnorm=1E-12_dp
        if(iter.gt.0) pulay%L2=fnorm/xnorm
        if(parallel%myid_k.eq.0) print *, 'L2 Conv: ', pulay%L2, pulay%L2_eps
        pulay%F_i(pulay%n,s)%ppointer(:)=fi(:)-f_last(:) 
        pulay%R_i(pulay%n,s)%ppointer(:)=x_last(:)-x_before_last(:)

        F_mu(pulay%n)=mu_fi-mu_flast
        R_mu(pulay%n)=mu_last-mu_beforelast

        if(present(all_PAW)) then
            if(iter.eq.0) then
                all_PAW%Rhoij_mix_m2(:,s)=all_PAW%Rhoij_mix_m1(:,s)!0.0_dp
                all_PAW%Rhoij_mix_fm1(:,s)=0.0_dp
            endif
            rhoij_fi(:)=all_PAW%Rhoij_mix_0(:,s)-all_PAW%Rhoij_mix_m1(:,s)

            all_PAW%Rhoij_mix_F_i(pulay%n,s)%ppointer(:)=rhoij_fi(:)-all_PAW%Rhoij_mix_fm1(:,s)
            all_PAW%Rhoij_mix_R_i(pulay%n,s)%ppointer(:)=all_PAW%Rhoij_mix_m1(:,s)-all_PAW%Rhoij_mix_m2(:,s)
        endif

        !linear mix / G1 term
        x_new(:)=x_last(:) + inv_diel_p(:)*fi(:)
        if(present(all_PAW)) then
            all_PAW%Rhoij_mix_0(:,s)= all_PAW%Rhoij_mix_m1(:,s) + A*rhoij_fi(:)
        endif
        mu_new=mu_last+ A*mu_fi

        m=max(pulay%n-iter+1,1)
        !Linear mix only on the first step to generate meaningful f_last for pulay%F_i(n,s)  
        if(iter.gt.1.and.(mod((iter+1),pulay%k).eq.0)) then
            FF=0.0_dp
            do i=m, pulay%n
            do j=m, pulay%n
                FF(i,j)=sum(conjg(pulay%F_i(i,s)%ppointer(:))*mix_metric_p(:)*pulay%F_i(j,s)%ppointer(:))
            enddo
            enddo
            call parallel_task('sum', FF, parallel, 'band')   
            FF=FF*product(grid%box_length(:)) 

            FTf=0.0_dp
            do i=m, pulay%n
                FTf(i)=sum(conjg(pulay%F_i(i,s)%ppointer(:))*mix_metric_p(:)*fi(:))
            enddo
            call parallel_task('sum', FTf, parallel, 'band')
            FTf=FTf*product(grid%box_length(:))

            W(m:pulay%n)=solve(FF(m:pulay%n,m:pulay%n),FTf(m:pulay%n)) !W=(FT.F)^(-1).(FT.fi)

            !Use same weights for density, rho_ij and chemical potential
            do i=m, pulay%n
                x_new(:)=x_new(:)-(pulay%R_i(i,s)%ppointer(:) +  &
                        inv_diel_p(:)*pulay%F_i(i,s)%ppointer(:))*W(i)
            enddo
            do i=m, pulay%n
                mu_new=mu_new-(R_mu(i) + A*F_mu(i))*real(W(i))
            enddo
            if(present(all_PAW)) then
                do i=m, pulay%n
                    all_PAW%Rhoij_mix_0(:,s)= all_PAW%Rhoij_mix_0(:,s) &
                    - (all_PAW%Rhoij_mix_R_i(i,s)%ppointer(:) + A*all_PAW%Rhoij_mix_F_i(i,s)%ppointer(:))*real(W(i))
                enddo
            endif

            if(parallel%myid.eq.0) print *, 'mu_new, mu_old, W(:):', mu_new, mu_last, real(W(m:pulay%n))

        endif

        x_before_last=x_last
        x_last=x_new
        f_last=fi

        mu_beforelast=mu_last
        mu_last=mu_new
        mu_flast=mu_fi

        if(present(all_PAW)) then
            if(.not.present(elements).or..not.present(atoms)) then
                print *, 'need elements and atoms into pulay mix with all_PAW'
                stop
            endif
            all_PAW%Rhoij_mix_m2(:,s)=all_PAW%Rhoij_mix_m1(:,s)
            all_PAW%Rhoij_mix_m1(:,s)=all_PAW%Rhoij_mix_0(:,s)
            all_PAW%Rhoij_mix_fm1(:,s)=rhoij_fi(:)
            deallocate(rhoij_fi)
            ij_all=0
            do at=1, size(atoms)
                if(elements(atoms(at)%element)%PP_type.ne.PAW) cycle
                    do ij=1, elements(atoms(at)%element)%PAW%tab%lmn2_size
                        ip=elements(atoms(at)%element)%PAW%tab%indklmn(7,ij)
                        jp=elements(atoms(at)%element)%PAW%tab%indklmn(8,ij)
                        ij_all=ij_all+1
                        atoms(at)%PAW%rho_ij(ip,jp,s)=all_PAW%Rhoij_mix_0(ij_all,s)
                    enddo
            enddo

            nucdipmom=0.0_dp
            do at=1, size(atoms)
                e=atoms(at)%element
                if(associated(atoms(at)%PAW%rhoij)) then !The Libpaw struct
                    do ij=1,atoms(at)%PAW%rhoij%lmn2_size
                        ip=elements(atoms(at)%element)%PAW%tab%indklmn(7,ij)
                        jp=elements(atoms(at)%element)%PAW%tab%indklmn(8,ij)
                        atoms(at)%PAW%rhoij%rhoij_(ij,:)=real(atoms(at)%PAW%rho_ij(ip,jp,s))
                    enddo
                    if(at.eq.1.and.parallel%myid.eq.0) then
                        print *, 'Atom 1 PAW Rho_ij (After mix)'
                        do ip=1, atoms(at)%PAW%rhoij%lmn_size
                                print *, real(atoms(at)%PAW%rho_ij(ip,:ip,:))
                        enddo
                     endif
                    call pawrhoij_filter(rhoij=atoms(at)%PAW%rhoij%rhoijp, &
                    rhoijselect=atoms(at)%PAW%rhoij%rhoijselect,nselect=atoms(at)%PAW%rhoij%nrhoijsel, &
                    cplex=atoms(at)%PAW%rhoij%cplex_rhoij,qphase=atoms(at)%PAW%rhoij%qphase, &
                    lmn2_size=atoms(at)%PAW%rhoij%lmn2_size, nspden=atoms(at)%PAW%rhoij%nspden,&
                    rhoij_input=atoms(at)%PAW%rhoij%rhoij_)
                endif
            enddo

            !Recalculate PAW onsite quanteties
            all_PAW%paral_atom=.true.
            all_PAW%my_atmtab=>null()
            call get_my_natom(parallel%comm_space,all_PAW%my_natom, all_PAW%N_PAW_atoms)
            call get_my_atmtab(parallel%comm_space,all_PAW%my_atmtab,all_PAW%my_atmtab_allocated, all_PAW%paral_atom, &
                                all_PAW%N_PAW_atoms)
            nucdipmom=0.0_dp
            call pawdenpot(compch_sph, all_PAW%epaw,epawdc,ipert=0,ixc=all_PAW%ixc, &
             my_natom=all_PAW%my_natom,natom=all_PAW%N_PAW_atoms,nspden=all_PAW%nspden,ntypat=all_PAW%N_PAW_elements, &
             nucdipmom=nucdipmom, nzlmopt=-1,option=0,paw_an=all_PAW%an,paw_an0=all_PAW%an, &
             paw_ij=all_PAW%my_ij,pawang=all_PAW%ang,pawprtvol=0,pawrad=all_PAW%rad,pawrhoij=all_PAW%rhoij, &
             pawspnorb=all_PAW%spnorb,pawtab=all_PAW%tab,pawxcdev=all_PAW%pawxcdev,spnorbscl=0.0_dp, &
             xclevel=all_PAW%xclevel,xc_denpos=all_PAW%xc_denpos,ucvol=product(grid%Box_Length),znucl=all_PAW%znucl,&
             mpi_atmtab=all_PAW%my_atmtab,comm_atom=parallel%comm_space)!,hyb_mixing,hyb_mixing_sr) ! optional arguments
    
            if(parallel%myid.eq.0) print *, 'compch_sph: from Abi spherical (After Mix)', compch_sph

            !recalculated Compensation_Charge_Density if it must be removed for XC potential
            usexcnhat=maxval(all_PAW%tab(1:all_PAW%N_PAW_elements)%usexcnhat)
            if(all_PAW%usexcnhat.eq.1) usexcnhat=1
            if(usexcnhat==0)  call calc_Compensation_Charge_Density(all_PAW, grids, atoms, parallel)
           
        endif
    end subroutine

    subroutine stochastic_ne_calc_scf(orbitals, stochastic, grids, atoms, elements, parallel, all_PAW)
        use odp_type, only :  orbital_struct, field_struct
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use system_type, only : system_struct
        use stochastic_type, only : stochastic_struct
        use simulation_type, only : potentials_struct, all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use operations_3D, only: integrate_3D_R
        use parallel_mod, only: parallel_task
        use Apply_Hamiltonian, only: Apply_S
        use odp, only : allocate_field, deallocate_field
        use Eigen_LOBPCG, only : Orbitals_Project_out

        type(orbital_struct), intent(inout) :: orbitals(:)
        type(stochastic_struct), intent(inout) :: stochastic
        type(parallel_struct), intent(in) :: parallel
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(all_PAW_struct), intent(inout) :: all_PAW
        type(field_struct), allocatable :: Sphi(:)
        integer  :: i, s

        stochastic%nelec_calc=0.0_dp
            do i= 1, size(orbitals)
                if(all_PAW%N_PAW_atoms.gt.0) then
                    allocate(Sphi(size(orbitals(i)%of)))
                        do s= 1, size(orbitals(i)%of)
                            Sphi%grid=orbitals(i)%of(s)%grid
                            grid=>grids(orbitals(i)%of(s)%grid)
                            call allocate_field(Sphi(s), grid, parallel)
                        end do
                    call Apply_S(orbitals(i)%of, Sphi, grids, atoms, elements, parallel, calc_R=.false.)
                    do s= 1, size(orbitals(i)%of)
                        grid=>grids(orbitals(i)%of(s)%grid)
                        stochastic%nelec_calc = stochastic%nelec_calc &
                                                + orbitals(i)%occ(s)*orbitals(i)%weight(s)* &
                                                    integrate_3D_R(abs(orbitals(i)%of(s)%R*conjg(Sphi(s)%R)), grid,  parallel)
                    enddo
                    do s = 1, size(Sphi)
                        call deallocate_field(Sphi(s),grids(Sphi(s)%grid))
                    enddo
                    deallocate(Sphi)
                else
                    do s= 1, size(orbitals(i)%of)
                        grid=>grids(orbitals(i)%of(s)%grid)
                        stochastic%nelec_calc = stochastic%nelec_calc &
                                                 + orbitals(i)%occ(s)*orbitals(i)%weight(s)* &
                                                    integrate_3D_R(abs(orbitals(i)%of(s)%R)**2, grid,  parallel)
                    enddo
                endif
            enddo
            call parallel_task('sum', stochastic%nelec_calc, parallel, 'diff_b')
        end subroutine

end module