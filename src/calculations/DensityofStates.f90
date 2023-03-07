module Density_of_States
    use types, only : dp

    implicit none

    public :: Calculate_Density_of_States


    contains

    subroutine Calculate_Density_of_States(sim)
        use simulation_type, only: simulation_struct
        use constants, only : pi
        use Stochastic_Mod, only : stochastic_delta, stochastic_vector_restore
        use odp_type, only : orbital_struct, field_struct
        use Chebychev, only : calculate_Chebychev_coefficients_FD_Gaussian, calculate_Chebychev_coefficients_Gaussian
        use odp, only : allocate_orbital, deallocate_orbital
        use Eigen_LOBPCG, only : Orbitals_Project_out
        use initialization_finalization, only : initialize_DOS
        use parallel_mod, only : parallel_task
        type(simulation_struct), intent(inout), target :: sim
        integer :: i,j,k,jj,k_loc,s_loc,s,n
        real(dp) :: w
        type(orbital_struct), pointer :: stoc_orbitals(:), ks_orbitals(:)
        type(orbital_struct), allocatable :: delta_orbitals(:)

        if(sim%parallel%myid.eq.0)  print *, sim%parallel%myid, 'Start Density of States'
        call initialize_DOS(sim)
        if(allocated(sim%DOS%dos)) deallocate(sim%DOS%dos)
        allocate(sim%DOS%dos(sim%DOS%nw))
        if(allocated(sim%DOS%fdos)) deallocate(sim%DOS%fdos)
        allocate(sim%DOS%fdos(sim%DOS%nw))
        sim%DOS%dos(:)=0.0_dp
        sim%DOS%fdos(:)=0.0_dp

        if(sim%system%n_deterministic.gt.0) then                
            if(allocated(sim%DOS%det_dos)) deallocate(sim%DOS%det_dos)
            allocate(sim%DOS%det_dos(sim%DOS%nw))
            if(allocated(sim%DOS%det_fdos)) deallocate(sim%DOS%det_fdos)
            allocate(sim%DOS%det_fdos(sim%DOS%nw))
            sim%DOS%det_dos(:)=0.0_dp
            sim%DOS%det_fdos(:)=0.0_dp
            do k=1,sim%system%n_spin; do j=1,sim%system%n_kpoints
                do i=1, sim%system%n_deterministic
                    do jj=1, sim%DOS%nw

                        w=jj*sim%DOS%dw+minval(sim%Minimum_eigenvalue)-10*sqrt(sim%DOS%a)
                        sim%DOS%det_dos(jj)= sim%DOS%det_dos(jj)  &
                        + sim%all_weight(j,k)*sim%all_degeneracy(i,j,k)*sim%all_filter(i,j,k) &
                        *exp(-(sim%all_eigs(i,j,k) - w)**2/sim%DOS%a)/(sqrt(pi*sim%DOS%a))
                    
                        if((sim%all_eigs(i,j,k)-sim%system%chemical_potential(k))/sim%system%Temperature &
                            .gt.50.0_dp) cycle

                        sim%DOS%det_fdos(jj)= sim%DOS%det_fdos(jj)  &
                        + sim%all_weight(j,k)*sim%all_degeneracy(i,j,k)*sim%all_filter(i,j,k) &
                        *exp(-(sim%all_eigs(i,j,k) - w)**2/sim%DOS%a)/(sqrt(pi*sim%DOS%a)) &
                        /(exp((sim%all_eigs(i,j,k)-sim%system%chemical_potential(k))/sim%system%Temperature)+1.0_dp)

                    enddo
                enddo
            enddo; enddo
            sim%DOS%dos=sim%DOS%det_dos
            sim%DOS%fdos=sim%DOS%det_fdos

        endif

        if(sim%system%n_stochastic.gt.0) then 
            if(allocated(sim%DOS%stoc_dos)) deallocate(sim%DOS%stoc_dos)
            allocate(sim%DOS%stoc_dos(sim%DOS%nw))
            if(allocated(sim%DOS%stoc_fdos)) deallocate(sim%DOS%stoc_fdos)
            allocate(sim%DOS%stoc_fdos(sim%DOS%nw))
            sim%DOS%stoc_dos(:)=0.0_dp
            sim%DOS%stoc_fdos(:)=0.0_dp
            do s_loc=1, size(sim%orbitals(:,:,:),3)
                s=sim%orbitals(1,1,s_loc)%spin
                do k_loc=1, size(sim%orbitals(:,:,:),2)
                    k=sim%orbitals(1,k_loc,s_loc)%k_point

                    stoc_orbitals(1:sim%system%n_stochastic/sim%parallel%n_band_groups) => &
                                sim%orbitals(sim%system%stochastic_start:sim%system%stochastic_end,k_loc,s_loc)
                    
                    allocate(delta_orbitals( sim%system%n_stochastic/ sim%parallel%n_band_groups))
                    do i = 1, size(delta_orbitals)
                        delta_orbitals(i)%band=stoc_orbitals(k)%band
                        delta_orbitals(i)%k_point=stoc_orbitals(k)%k_point
                        delta_orbitals(i)%spin=stoc_orbitals(k)%spin
                        delta_orbitals(i)%degeneracy=stoc_orbitals(k)%degeneracy
                        delta_orbitals(i)%n_spinor=stoc_orbitals(k)%n_spinor
                        call allocate_orbital(delta_orbitals(i), sim%grids, delta_orbitals(i)%k_point+2, sim%parallel)
                        delta_orbitals(i)%occ=stoc_orbitals(k)%occ
                        delta_orbitals(i)%weight=stoc_orbitals(k)%weight
                        delta_orbitals(i)%filter=stoc_orbitals(k)%filter
                        delta_orbitals(i)%type=stoc_orbitals(k)%type
                        if(allocated(delta_orbitals(i)%seed_array)) deallocate(delta_orbitals(i)%seed_array)
                        allocate(delta_orbitals(i)%seed_array(size(stoc_orbitals(k)%seed_array)))
                        delta_orbitals(i)%seed_array=stoc_orbitals(k)%seed_array
                    enddo

                    call stochastic_vector_restore(delta_orbitals(:), sim%system, &
                                sim%grids, sim%parallel, sim%stoc_norm, sim%stoc_occ, &
                                sim%all_PAW, sim%atoms, sim%elements)

                    if(sim%system%n_deterministic.gt.0) then
                        ks_orbitals(1:sim%system%n_deterministic/sim%parallel%n_band_groups) => &
                        sim%orbitals(sim%system%deterministic_start:sim%system%deterministic_end,k_loc,s_loc)

                        call Orbitals_Project_out(ks_orbitals(:), delta_orbitals(:), &
                            sim%grids, sim%parallel, sim%atoms, sim%elements, sim%all_PAW, sim%all_filter(:,k,s))
                    endif

                    call stochastic_delta(delta_orbitals(:), sim%stochastic(k_loc,s_loc), sim%DOS, sim%grids, &
                             sim%potentials, sim%atoms, sim%elements, sim%parallel, sim%all_PAW)
                    
                    do jj=1, sim%DOS%nw
                        w=jj*sim%DOS%dw+minval(sim%Minimum_eigenvalue)-10*sqrt(sim%DOS%a)
                        call calculate_Chebychev_coefficients_Gaussian(sim%stochastic(k_loc,s_loc)%Ebar, &
                             sim%stochastic(k_loc,s_loc)%deltaE, sim%stochastic(k_loc,s_loc)%DOS_Cn, &
                             size(sim%stochastic(k_loc,s_loc)%Mn_Dos(:)), w, sim%DOS%a)

                        do n=1, size(sim%stochastic(k_loc,s_loc)%Mn_Dos(:))
                                sim%DOS%stoc_dos(jj)= sim%DOS%stoc_dos(jj)&
                                    +real(sim%stochastic(k_loc,s_loc)%Mn_Dos(n))*sim%stochastic(k_loc,s_loc)%DOS_Cn(n)
                        enddo
                        
                        call calculate_Chebychev_coefficients_FD_Gaussian(sim%stochastic(k_loc,s_loc)%Ebar, &
                            sim%stochastic(k_loc,s_loc)%deltaE, sim%stochastic(k_loc,s_loc)%DOS_Cn, &
                            size(sim%stochastic(k_loc,s_loc)%Mn_Dos(:)), &
                            w, sim%DOS%a, sim%system%chemical_potential(s), sim%system%Temperature)

                        do n=1, size(sim%stochastic(k_loc,s_loc)%Mn_Dos(:))
                                sim%DOS%stoc_fdos(jj)= sim%DOS%stoc_fdos(jj)&
                                    +real(sim%stochastic(k_loc,s_loc)%Mn_Dos(n))*sim%stochastic(k_loc,s_loc)%DOS_Cn(n)
                        enddo
                    enddo

                    do i = 1, size(delta_orbitals)
                        call deallocate_orbital(delta_orbitals(i),sim%grids)
                    enddo
                    deallocate(delta_orbitals)
                enddo
            enddo
            call parallel_task('sum',  sim%DOS%stoc_dos, sim%parallel, 'diff_k')
            call parallel_task('sum',  sim%DOS%stoc_dos, sim%parallel, 'diff_s')
            call parallel_task('sum',  sim%DOS%stoc_fdos, sim%parallel, 'diff_k')
            call parallel_task('sum',  sim%DOS%stoc_fdos, sim%parallel, 'diff_s')

            sim%DOS%dos=sim%DOS%dos+sim%DOS%stoc_dos
            sim%DOS%fdos=sim%DOS%fdos+sim%DOS%stoc_fdos

        endif
            

    if(sim%parallel%myid.eq.0) call DOS_Output(sim%DOS, minval(sim%Minimum_eigenvalue))

    end subroutine



    subroutine DOS_Output(DOS, minimum)
        use Density_of_states_type, only : Density_of_states_struct
        use constants, only : FMT_real

        type(Density_of_states_struct),intent(inout) ::DOS
        integer :: jj
        real(dp) :: w, minimum

        open(newunit=DOS%output_unit,  file="DOS.data", status="unknown", recl=1024)
        write(DOS%output_unit, *) "1:Freq(a.u.)", '2:DOS(1/a.u.)', '3:FDOS(1/a.u.)'
        do jj=1, DOS%nw
            w=jj*DOS%dw+minimum-10*sqrt(DOS%a)
            if(allocated(DOS%det_dos).and.allocated(DOS%stoc_dos)) then
                write(DOS%output_unit, FMT_real) w, DOS%dos(jj), DOS%fdos(jj), DOS%det_dos(jj),  &
                DOS%det_fdos(jj), DOS%stoc_dos(jj), DOS%stoc_fdos(jj)
            else
                write(DOS%output_unit, FMT_real) w, DOS%dos(jj), DOS%fdos(jj)
            endif
        enddo
        close(DOS%output_unit)
    end subroutine

    
end module