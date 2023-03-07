module Stochastic_Mod
    use types, only : dp
    use constants, only: i_, pi

    implicit none

    public :: stochastic_vector_initializer, stochastic_vector_restore, stochastic_FD_filter, Stoc_Entropy,&
              stochastic_rootFD_filter, stochastic_OMFD_filter, stochastic_root_OMFD_filter


    contains

    subroutine stochastic_vector_initializer(orbitals, system, grids, parallel, stoc_norm, stoc_occ, all_PAW, atoms, elements)
        use odp_type, only : field_struct, orbital_struct
        use grids_type, only : grid_struct
        use system_type, only : system_struct
        use parallel_type, only : parallel_struct
        use operations_3D, only:  integrate_3D_G
        use fft, only: real_to_recip, recip_to_real
        use simulation_type, only : all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use parallel_mod, only : parallel_task, parallel_wait
        use Non_Local_ion, only : Apply_S_power
        use odp, only : allocate_field, deallocate_field
        use grids_mod,only : allocate_local_fields_R, allocate_local_fields_G

        type(parallel_struct), intent(in) :: parallel
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(system_struct), intent(in) :: system
        type(all_PAW_struct), intent(inout) :: all_PAW
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        
        real(dp), intent(inout) :: stoc_norm(:,:,:),stoc_occ(:,:,:)
        integer :: rand_size
        integer, allocatable ::  seed_array(:), seed_array_save(:)

        if(size(orbitals).lt.1) return

        !Save the incoming seed array on each rank
        call random_seed(size=rand_size)
        allocate(seed_array_save(rand_size))
        call random_seed(get=seed_array_save)

        !Fresh single stream RNG for all orbitals on all proc       
        call random_init(.false.,.false.)
        !call random_seed()
        call random_seed(size=rand_size)
        call parallel_task('bcast', rand_size, parallel, 'all', root=0)
        allocate(seed_array(rand_size))
        call random_seed(get=seed_array)
        call parallel_task('bcast', seed_array, parallel, 'all', root=0)
        call random_seed(put=seed_array)

        call stochastc_vector_fill(orbitals, system, grids, parallel, stoc_norm,stoc_occ, &
                                     all_PAW, atoms, elements, seed_array, .true.)

        !Restore the incoming seed array on each rank
        call random_seed(put=seed_array_save)
        deallocate(seed_array_save)
        deallocate(seed_array)

    end subroutine

    subroutine stochastic_vector_restore(orbitals, system, grids, parallel, stoc_norm, stoc_occ, all_PAW, atoms, elements)
        use odp_type, only : field_struct, orbital_struct
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use operations_3D, only:  integrate_3D_G
        use fft, only: real_to_recip, recip_to_real
        use simulation_type, only : all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use Non_Local_ion, only : Apply_S_power
        use odp, only : allocate_field, deallocate_field
        use grids_mod,only : allocate_local_fields_R, allocate_local_fields_G
        use system_type, only : system_struct

        type(parallel_struct), intent(in) :: parallel
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(all_PAW_struct), intent(inout) :: all_PAW
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(system_struct), intent(in) :: system
        real(dp), intent(inout) :: stoc_norm(:,:,:),stoc_occ(:,:,:)

        integer, allocatable :: seed_array_save(:)

        integer ::  rand_size

        !Save the incoming seed array on each rank
        call random_seed(size=rand_size)
        allocate(seed_array_save(rand_size))
        call random_seed(get=seed_array_save)

        !stored single stream RNG for all orbitals on all proc
        call random_seed(put=orbitals(1)%seed_array)

        call stochastc_vector_fill(orbitals, system, grids, parallel, stoc_norm, stoc_occ, &
                                    all_PAW, atoms, elements, orbitals(1)%seed_array, .false.)

        !Restore the incoming seed array on each rank
        call random_seed(put=seed_array_save)
        deallocate(seed_array_save)
    end subroutine

    subroutine stochastc_vector_fill(orbitals, system, grids, parallel, &
            stoc_norm,stoc_occ, all_PAW, atoms, elements, seed_array, init)
        use odp_type, only : field_struct, orbital_struct
        use grids_type, only : grid_struct
        use system_type, only : system_struct
        use parallel_type, only : parallel_struct
        use operations_3D, only:  integrate_3D_G
        use fft, only: real_to_recip, recip_to_real
        use simulation_type, only : all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use parallel_mod, only : parallel_task, parallel_wait
        use Non_Local_ion, only : Apply_S_power
        use odp, only : allocate_field, deallocate_field
        use grids_mod,only : allocate_local_fields_R, allocate_local_fields_G

        type(parallel_struct), intent(in) :: parallel
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        type(system_struct), intent(in) :: system
        type(all_PAW_struct), intent(inout) :: all_PAW
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        
        logical, intent(in) :: init
        real(dp), intent(inout) :: stoc_norm(:,:,:),stoc_occ(:,:,:)
        integer :: n_orbitals
        integer :: i,s ,k, s2, ix,iy,iz, ix_local,iy_local,iz_local,my_shift(4)
        real(dp), allocatable :: random_number_field(:,:,:,:)
        integer, intent(in) ::  seed_array(:)
        real(dp) :: norm
        type(field_struct), allocatable :: Sinv(:)
        grid=>grids(orbitals(1)%of(1)%grid)

        if(system%stoch_gen_type.eq.2) then
            call allocate_local_fields_G(random_number_field, grid, size(orbitals(1)%of))
        else
            call allocate_local_fields_R(random_number_field, grid, size(orbitals(1)%of))
        endif
        !call random_seed()
        
        i=1

        do while(i.le.size(orbitals))
        !do i=1, size(orbitals)
            if(init) orbitals(i)%seed_array=seed_array !All proc and bands same seed array, staggered use, wont actually use %seed_array 


            if(i.le.system%n_orbitals_stochastic_full) then 
                n_orbitals=system%n_orbitals_stochastic_full
            else
                print *, 'Orbital arangement messed up'
                stop
            endif

            grid=>grids(orbitals(i)%of(1)%grid)
            
            if(init) orbitals(i)%occ(:)=2.0_dp/(1.0_dp*system%n_spinor*system%n_spin)
            if(init) orbitals(i)%filter(:)=1.0_dp
            !Assign the random number field and store the seed_array used to generate it
            !call random_seed(size=rand_size)
            !if(allocated(orbitals(i)%seed_array)) deallocate(orbitals(i)%seed_array)
            !allocate(orbitals(i)%seed_array(rand_size))
            !call random_seed(get=orbitals(i)%seed_array)
            !seed_array=1 !set for debug
            
            !if(parallel%myid.ne.83) then
            !    orbitals(i)%seed_array(:)= &
            !     abs( mod((orbitals(i)%seed_array(:)*181)*((parallel%myid-83)*359), 104729) )
            !else
            !    orbitals(i)%seed_array(:)= &
            !     abs( mod((orbitals(i)%seed_array(:)*193)*((parallel%myid-97)*389), 104729) )
            !endif
            !call random_seed(put=orbitals(i)%seed_array)

            !call random_number(random_number_field)
                !orbitals(i)%of(s)%R(:,:,:)=sqrt(product(grid%Nr/grid%box_length)) &
                !                            *exp(i_*pi*2.0_dp*random_number_field(:,:,:))       
                
                do k=1, parallel%nproc
                    call random_number(random_number_field)
                    if(k.eq.(parallel%myid+1)) then
                        do s=1, size(orbitals(i)%of)
                            if(.not.grid%gamma) then
                                orbitals(i)%of(s)%R(:,:,:)=sqrt(product(grid%Nr/grid%box_length)) &
                                                *exp(i_*pi*2.0_dp*random_number_field(:,:,:,s))
                            else
                                orbitals(i)%of(s)%R(:,:,:)=(random_number_field(:,:,:,s)-0.5_dp)
                                orbitals(i)%of(s)%R(:,:,:)=orbitals(i)%of(s)%R(:,:,:)/abs(orbitals(i)%of(s)%R(:,:,:))
                                orbitals(i)%of(s)%R(:,:,:)=orbitals(i)%of(s)%R(:,:,:)*sqrt(product(grid%Nr/grid%box_length))
                            endif
                        enddo
                    endif
                enddo

                if(system%stoch_gen_type.eq.1) then
                    do s=1, size(orbitals(i)%of)
                        call real_to_recip(orbitals(i)%of(s), grids)
                    enddo
                endif

                do s=1, size(orbitals(i)%of)
                    orbitals(i)%of(s)%G=orbitals(i)%of(s)%G*grid%cutwf
                    norm=real(integrate_3D_G(conjg(orbitals(i)%of(s)%G)*orbitals(i)%of(s)%G, grid, parallel))
                    if(init)    stoc_norm(orbitals(i)%band,orbitals(i)%k_point,orbitals(i)%spin*s)=norm
                    if(init)  orbitals(i)%occ(s)=orbitals(i)%occ(s)*norm 
                    if(norm.gt.tiny(1.0_dp)) then
                        orbitals(i)%of(s)%G=orbitals(i)%of(s)%G/sqrt(norm)
                        call recip_to_real(orbitals(i)%of(s), grids)
                    endif
                enddo
                if(init) then
                    if(sum(stoc_norm(orbitals(i)%band,orbitals(i)%k_point,orbitals(i)%spin:)).le.tiny(1.0_dp)) cycle
                endif
                if(all_PAW%N_PAW_atoms.gt.0) then
                    allocate(Sinv(size(orbitals(i)%of)))
                    do s2=1, size(Sinv)
                        Sinv(s2)%grid=orbitals(i)%of(s2)%grid
                        call allocate_field(Sinv(s2), grids(Sinv(s2)%grid), parallel)
                    enddo
                    call Apply_S_power(orbitals(i)%of, atoms, elements,  &
                            grids, parallel, Sinv, -0.5_dp)
                    do s2 = 1, size(Sinv)
                        orbitals(i)%of(s2)%R=Sinv(s2)%R
                        call real_to_recip(orbitals(i)%of(s2), grids)
                        orbitals(i)%of(s2)%G=orbitals(i)%of(s2)%G &
                            *grids(Sinv(s2)%grid)%cutwf
                        call recip_to_real(orbitals(i)%of(s2), grids)
                        call deallocate_field(Sinv(s2),grids(Sinv(s2)%grid))
                    enddo
                    deallocate(Sinv)
                endif
                i=i+1 !advance orbital
        end do 
        if(init) then
        do i=1, size(orbitals)
            do s=1, size(orbitals(i)%of)
                stoc_occ(orbitals(i)%band,orbitals(i)%k_point,orbitals(i)%spin*s)=orbitals(i)%occ(s)
            enddo
        enddo
        endif

        deallocate(random_number_field)

    end subroutine

    subroutine stochastic_rootFD_filter(orbitals, stochastic, system, grids, potentials, atoms, elements, parallel, all_PAW)
        use odp_type, only :  orbital_struct, field_struct
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use system_type, only : system_struct
        use stochastic_type, only : stochastic_struct
        use simulation_type, only : potentials_struct, all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use Chebychev, only : cheby_exp_and_app, calculate_Chebychev_coefficients_sqrt_fd
        use operations_3D, only: integrate_3D_R
        use parallel_mod, only: parallel_task
        use Apply_Hamiltonian, only: Apply_S
        use odp, only : allocate_field, deallocate_field

        type(orbital_struct), intent(inout) :: orbitals(:)
        type(stochastic_struct), intent(inout) :: stochastic
        type(system_struct), intent(in) :: system
        type(parallel_struct), intent(in) :: parallel
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(all_PAW_struct), intent(inout) :: all_PAW
        integer  :: Nc, Nc_max, i, power


            Nc_max=min(4096,5*(floor(stochastic%deltaE/system%temperature)+1)+10)
            Nc_max=max(Nc_max,128) !lower bound on Nc 
            !Nc_max needs to be power of 2 
            power=1
            do while(power<Nc_max)
                power=power*2
            enddo
            Nc_max=power

            if(allocated(stochastic%Sqrt_Density_Cn)) deallocate(stochastic%Sqrt_Density_Cn)
            if(allocated(stochastic%Mn)) deallocate(stochastic%Mn)

            allocate(stochastic%Sqrt_Density_Cn(Nc_max))
            call calculate_Chebychev_coefficients_sqrt_fd(stochastic%Ebar, stochastic%deltaE, stochastic%Sqrt_Density_Cn, Nc_max,&
                                    system%Temperature, system%chemical_potential(orbitals(1)%spin))
            Nc=Nc_max
            do i=2, Nc_max
            if(all(abs(stochastic%Sqrt_Density_Cn(i:)).lt.(1.0_dp/sqrt(2.0_dp)*1E-9))) then
                    Nc=i-1
                    exit
            endif
            enddo
            allocate(stochastic%Mn(Nc))
            call cheby_exp_and_app(orbitals, &
                    stochastic%Ebar, stochastic%deltaE, stochastic%Mn, stochastic%Sqrt_Density_Cn(:Nc), &
                                                     grids, potentials, atoms, elements, parallel, all_PAW)
            
    end subroutine

    subroutine stochastic_FD_filter(orbitals, stochastic, system, grids, potentials, atoms, elements, parallel, all_PAW)
        use odp_type, only :  orbital_struct, field_struct
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use system_type, only : system_struct
        use stochastic_type, only : stochastic_struct
        use simulation_type, only : potentials_struct, all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use Chebychev, only : cheby_exp_and_app, calculate_Chebychev_coefficients_fd
        use operations_3D, only: integrate_3D_R
        use parallel_mod, only: parallel_task

        type(orbital_struct), intent(inout) :: orbitals(:)
        type(stochastic_struct), intent(inout) :: stochastic
        type(system_struct), intent(in) :: system
        type(parallel_struct), intent(in) :: parallel
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(all_PAW_struct), intent(inout) :: all_PAW

        integer  :: Nc, Nc_max, i, power

            Nc_max=min(4096,5*(floor(stochastic%deltaE/system%temperature)+1)+10)
            Nc_max=max(Nc_max,128) !lower bound on Nc 
            !Nc_max needs to be power of 2 
            power=1
            do while(power<Nc_max)
                power=power*2
            enddo
            Nc_max=power

            if(allocated(stochastic%Density_Cn)) deallocate(stochastic%Density_Cn)
            if(allocated(stochastic%Mn)) deallocate(stochastic%Mn)

            allocate(stochastic%Density_Cn(Nc_max))
            call calculate_Chebychev_coefficients_fd(stochastic%Ebar, stochastic%deltaE, stochastic%Density_Cn, Nc_max,&
                                    system%Temperature, system%chemical_potential(orbitals(1)%spin))
            Nc=Nc_max
            do i=2, Nc_max
            if(all(abs(stochastic%Density_Cn(i:)).lt.(1E-9))) then
                    Nc=i-1
                    exit
            endif
            enddo
            allocate(stochastic%Mn(Nc))

            call cheby_exp_and_app(orbitals, &
                    stochastic%Ebar, stochastic%deltaE, stochastic%Mn, stochastic%Density_Cn(:Nc), &
                                                     grids, potentials, atoms, elements, parallel, all_PAW)
    end subroutine

    subroutine stochastic_OMFD_filter(orbitals, stochastic, system, grids, potentials, atoms, elements, parallel, all_PAW)
        use odp_type, only :  orbital_struct, field_struct
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use system_type, only : system_struct
        use stochastic_type, only : stochastic_struct
        use simulation_type, only : potentials_struct, all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use Chebychev, only : cheby_exp_and_app, calculate_Chebychev_coefficients_omfd
        use operations_3D, only: integrate_3D_R
        use parallel_mod, only: parallel_task

        type(orbital_struct), intent(inout) :: orbitals(:)
        type(stochastic_struct), intent(inout) :: stochastic
        type(system_struct), intent(in) :: system
        type(parallel_struct), intent(in) :: parallel
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(all_PAW_struct), intent(inout) :: all_PAW

        integer  :: Nc, Nc_max, i, power

            Nc_max=min(4096,5*(floor(stochastic%deltaE/system%temperature)+1)+10)
            Nc_max=max(Nc_max,128) !lower bound on Nc 
            !Nc_max needs to be power of 2 
            power=1
            do while(power<Nc_max)
                power=power*2
            enddo
            Nc_max=power

            if(allocated(stochastic%OMDensity_Cn)) deallocate(stochastic%OMDensity_Cn)
            if(allocated(stochastic%Mn)) deallocate(stochastic%Mn)

            allocate(stochastic%OMDensity_Cn(Nc_max))
            call calculate_Chebychev_coefficients_omfd(stochastic%Ebar, stochastic%deltaE, stochastic%OMDensity_Cn, Nc_max,&
                                    system%Temperature, system%chemical_potential(orbitals(1)%spin))
            Nc=Nc_max
            do i=2, Nc_max
            if(all(abs(stochastic%OMDensity_Cn(i:)).lt.(1E-9))) then
                    Nc=i-1
                    exit
            endif
            enddo
            allocate(stochastic%Mn(Nc))

            call cheby_exp_and_app(orbitals, &
                    stochastic%Ebar, stochastic%deltaE, stochastic%Mn, stochastic%OMDensity_Cn(:Nc), &
                                                     grids, potentials, atoms, elements, parallel, all_PAW)
    end subroutine

    subroutine stochastic_root_OMFD_filter(orbitals, stochastic, system, grids, potentials, atoms, elements, parallel, all_PAW)
        use odp_type, only :  orbital_struct, field_struct
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use system_type, only : system_struct
        use stochastic_type, only : stochastic_struct
        use simulation_type, only : potentials_struct, all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use Chebychev, only : cheby_exp_and_app, calculate_Chebychev_coefficients_root_omfd
        use operations_3D, only: integrate_3D_R
        use parallel_mod, only: parallel_task

        type(orbital_struct), intent(inout) :: orbitals(:)
        type(stochastic_struct), intent(inout) :: stochastic
        type(system_struct), intent(in) :: system
        type(parallel_struct), intent(in) :: parallel
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(all_PAW_struct), intent(inout) :: all_PAW

        integer  :: Nc, Nc_max, i, power

            Nc_max=min(4096,5*(floor(stochastic%deltaE/system%temperature)+1)+10)
            Nc_max=max(Nc_max,128) !lower bound on Nc 
            !Nc_max needs to be power of 2 
            power=1
            do while(power<Nc_max)
                power=power*2
            enddo
            Nc_max=power

            if(allocated(stochastic%OMDensity_Cn)) deallocate(stochastic%OMDensity_Cn)
            if(allocated(stochastic%Mn)) deallocate(stochastic%Mn)

            allocate(stochastic%OMDensity_Cn(Nc_max))
            call calculate_Chebychev_coefficients_root_omfd(stochastic%Ebar, stochastic%deltaE, stochastic%OMDensity_Cn, Nc_max,&
                                    system%Temperature, system%chemical_potential(orbitals(1)%spin))
            Nc=Nc_max
            do i=2, Nc_max
            if(all(abs(stochastic%OMDensity_Cn(i:)).lt.(1E-9))) then
                    Nc=i-1
                    exit
            endif
            enddo
            allocate(stochastic%Mn(Nc))

            call cheby_exp_and_app(orbitals, &
                    stochastic%Ebar, stochastic%deltaE, stochastic%Mn, stochastic%OMDensity_Cn(:Nc), &
                                                     grids, potentials, atoms, elements, parallel, all_PAW)
    end subroutine
    
    subroutine stochastic_delta(orbitals, stochastic, dos, grids, potentials, atoms, elements, parallel, all_PAW)
        use odp_type, only :  orbital_struct, field_struct
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use stochastic_type, only : stochastic_struct
        use simulation_type, only : potentials_struct, all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use Chebychev, only : cheby_exp_and_app, calculate_Chebychev_coefficients_Gaussian
        use operations_3D, only: integrate_3D_R
        use parallel_mod, only: parallel_task
        use density_of_states_type, only : density_of_states_struct

        type(orbital_struct), intent(inout) :: orbitals(:)
        type(stochastic_struct), intent(inout) :: stochastic
        type(parallel_struct), intent(in) :: parallel
        type(density_of_states_struct), intent(in) :: dos
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(all_PAW_struct), intent(inout) :: all_PAW

        integer  :: Nc, Nc_max, i, power

            Nc_max=min(32768,10*(floor(stochastic%deltaE/sqrt(dos%a))+1)+10)
            Nc_max=max(Nc_max,128) !lower bound on Nc 
            !Nc_max needs to be power of 2 
            power=1
            do while(power<Nc_max)
                power=power*2
            enddo
            Nc_max=power

            if(allocated(stochastic%Mn_Dos)) deallocate(stochastic%Mn_Dos)
            if(allocated(stochastic%DOS_Cn)) deallocate(stochastic%DOS_Cn)

            allocate(stochastic%DOS_Cn(Nc_max))
            call calculate_Chebychev_coefficients_Gaussian(stochastic%Ebar, stochastic%deltaE, stochastic%DOS_Cn, Nc_max,&
                                    0.0_dp, dos%a)
            Nc=Nc_max
            do i=2, Nc_max
            if(all(abs(stochastic%DOS_Cn(i:)).lt.(1/sqrt(2.0_dp)*1E-9))) then
                    Nc=i-1
                    exit
            endif
            enddo
            allocate(stochastic%Mn_Dos(Nc))

            call cheby_exp_and_app(orbitals, &
                    stochastic%Ebar, stochastic%deltaE, stochastic%Mn_Dos, stochastic%DOS_Cn(:Nc), &
                                                     grids, potentials, atoms, elements, parallel, all_PAW)
            
    end subroutine

    function Stochastic_Entropy_ks(stochastic, T_au, mu) result(entrop)
        use stochastic_type, only : stochastic_struct
        use Chebychev, only: calculate_Chebychev_coefficients_entropy
        type(stochastic_struct), intent(inout) :: stochastic
        real(dp), intent(in) :: mu, T_au
        real(dp):: entrop
        integer  :: Nc_max, n, power

        if(allocated(stochastic%Mn)) then
            Nc_max=min(4096,5*(floor(stochastic%deltaE/T_au)+1)+10)
            Nc_max=max(Nc_max,128) !lower bound on Nc 
            !Nc_max needs to be power of 2 
            power=1
            do while(power<Nc_max)
                power=power*2
            enddo
            Nc_max=power

            if(allocated(stochastic%Entropy_Cn)) deallocate(stochastic%Entropy_Cn)
            
            allocate(stochastic%Entropy_Cn(Nc_max))

            call calculate_Chebychev_coefficients_entropy(stochastic%Ebar, stochastic%deltaE, stochastic%Entropy_Cn, Nc_max, &
                                    T_au, mu)
            entrop=0.0_dp
            do n=1, min(size(stochastic%Mn(:)),Nc_max)
                entrop=entrop+real(stochastic%Mn(n))*stochastic%Entropy_Cn(n)
            enddo
        endif
    end function

    function Stoc_Entropy(sim) result (Entropy)
        use simulation_type, only : simulation_struct
        type(simulation_struct), intent(inout) :: sim
        real(dp) :: Entropy
        integer:: i,j
        Entropy=0.0_dp
        if(sim%system%n_orbitals_stochastic_full.gt.0) then
        do j=1,size(sim%stochastic,2);do i=1,size(sim%stochastic,1)
            Entropy=Entropy+Stochastic_Entropy_ks(sim%stochastic(i,j), &
                sim%system%temperature, sim%system%chemical_potential(j))
        enddo; enddo
        endif
        Entropy=sim%system%temperature*Entropy
    end function
        
end module