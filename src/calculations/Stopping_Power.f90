module Stopping_Power
    use types, only : dp
    use constants , only : i_

    implicit none

    public :: Stopping_Power_Calculation, Stopping_Power_Initialization


    contains

    subroutine Stopping_Power_Calculation(orbitals, stopping, density, potentials, &
            system, td, atoms, elements, parallel, grids, all_PAW, density_0, orbitals_0)
        use stopping_type, only : stopping_struct
        use td_type, only : td_struct
        use atom_type, only : atom_struct
        use element_type, only : element_struct
        use grids_type, only : grid_struct
        use system_type, only : system_struct
        use parallel_type, only : parallel_struct
        use odp_type, only : field_struct, orbital_struct, spin_DenPot_struct
        use operations_3D, only : integrate_3D_G
        use Local_Ion, only : calculate_local_ion_force
        use Non_Local_ion_EnFoSt, only : NL_PP_Force_RL
        use simulation_type, only : potentials_struct
        use simulation_type, only : all_PAW_struct
        use constants, only : FMT_real

        type(parallel_struct), intent(in) :: parallel
        type(spin_DenPot_struct), intent(in) :: density
        type(spin_DenPot_struct), intent(in), optional ::  density_0
        type(system_struct), intent(in) :: system
        type(orbital_struct), intent(in) :: orbitals(:,:,:)
        type(orbital_struct), intent(in), optional :: orbitals_0(:,:,:)
        type(element_struct), intent(inout) :: elements(:)
        type(atom_struct), intent(inout) :: atoms(:)
        type(stopping_struct), intent(inout) :: stopping
        type(grid_struct), intent(inout) :: grids(:)
        type(td_struct), intent(in) :: td
        type(potentials_struct), intent(inout) :: potentials 
        type(all_PAW_struct), intent(inout) :: all_PAW

        logical :: updates(size(atoms))
        real(dp) :: F_in(3), F_0(3)
        integer :: s
        
        updates=atoms(:)%update_me_force
        F_in(:)=atoms(stopping%stop_at)%F(:)
        !Remove T=0 forces if t=0 density/orbitals present and all background atoms frozen 
        F_0(:)=0.0_dp
        if(all(atoms(:(stopping%stop_at-1))%frozen_R)) then
            atoms(:)%update_me_force=.false.
            atoms(stopping%stop_at)%update_me_force=.true.
            atoms(stopping%stop_at)%F_loc(:)=0.0_dp
            if(present(density_0)) call  calculate_local_ion_force(density_0, system, atoms(:), elements(:), parallel, grids)
            F_0(:)=F_0(:)+atoms(stopping%stop_at)%F_loc(:)
            atoms(stopping%stop_at)%F_nonloc(:)=0.0_dp
            if(present(orbitals_0)) call NL_PP_Force_RL(orbitals_0, atoms, elements, potentials, grids, parallel, all_PAW)
            F_0(:)=F_0(:)+atoms(stopping%stop_at)%F_nonloc(:) 
            print *, 'F0:', F_0(:), atoms(stopping%stop_at)%F_loc(:), atoms(stopping%stop_at)%F_nonloc(:), &
                             atoms(stopping%stop_at)%F_nn(:)
            if(.not.atoms(stopping%stop_at)%frozen_V.and.present(density_0)) then
                F_0(:)=F_0(:)+atoms(stopping%stop_at)%F_nn(:)
            endif
        endif
        if(atoms(stopping%stop_at)%frozen_V) then
            atoms(:)%update_me_force=.false.
            atoms(stopping%stop_at)%update_me_force=.true.
            call  calculate_local_ion_force(density, system, atoms(:), elements(:), parallel, grids)
            atoms(stopping%stop_at)%F_nonloc(1)=0.0_dp
            atoms(stopping%stop_at)%F_nonloc(2)=0.0_dp
            atoms(stopping%stop_at)%F_nonloc(3)=0.0_dp
            call NL_PP_Force_RL(orbitals, atoms, elements, potentials, grids, parallel, all_PAW)

            atoms(stopping%stop_at)%F(:)= atoms(stopping%stop_at)%F_nonloc(:)  &
                + atoms(stopping%stop_at)%F_loc(:)
    
            if(.not.all(atoms(:(stopping%stop_at-1))%frozen_R).or..not.present(density_0)) &
                atoms(stopping%stop_at)%F(:)=atoms(stopping%stop_at)%F(:) + atoms(stopping%stop_at)%F_nn(:)
            print *, "F_full:", atoms(stopping%stop_at)%F(:), &
                                atoms(stopping%stop_at)%F_loc(:), &
                                atoms(stopping%stop_at)%F_nonloc(:), atoms(stopping%stop_at)%F_nn(:)
        endif
        atoms(stopping%stop_at)%F(:)=atoms(stopping%stop_at)%F(:)-F_0(:)
        stopping%power=-sum(atoms(stopping%stop_at)%F(:)*atoms(stopping%stop_at)%P(:))/sqrt(sum(atoms(stopping%stop_at)%P(:)**2))
        do s=1, size(density%of)
            stopping%density_t(s)=real(integrate_3D_G(density%of(s)%G* &
            exp(i_*(grids(1)%G(:,:,:,1)*atoms(stopping%stop_at)%R(1) + grids(1)%G(:,:,:,2)*atoms(stopping%stop_at)%R(2) + &
                     grids(1)%G(:,:,:,3)*atoms(stopping%stop_at)%R(3))),grids(1), parallel))/product(grids(1)%box_length)
        enddo
        !stopping%density_t=stopping%density_t/product(grids(1)%Nr)
        if((td%time-stopping%t0).gt.tiny(1.0_dp)) then
            if(stopping%t0.gt.td%last_time) then
                stopping%power_average=stopping%power
            else
                stopping%power_average=stopping%power_average*(td%last_time-stopping%t0)
                stopping%power_average=stopping%power_average+stopping%power*(td%time-td%last_time)
                stopping%power_average=stopping%power_average/(td%time-stopping%t0)
            endif
        else
            stopping%power_average=0.0_dp
        endif

        if(parallel%myid.eq.0) then
            print *, 'Projectile data:'
            print *, 'Time: ', td%time, 'R: ', atoms(stopping%stop_at)%R
            print *, 'P: ', atoms(stopping%stop_at)%P
            print *, 'F: ', atoms(stopping%stop_at)%F

            write(stopping%output_unit, FMT_real) &
                                            (td%time-stopping%t0), atoms(stopping%stop_at)%R(:), atoms(stopping%stop_at)%P(:),  &
                                            stopping%power, stopping%power_average, sum(stopping%density_t), &
                                            sum(stopping%density_t)/sum(system%nelec)*product(grids(1)%box_length),  &
                                            atoms(stopping%stop_at)%F(:)
            flush(stopping%output_unit)
        endif

        atoms(:)%update_me_force=updates(:)
        atoms(stopping%stop_at)%F(:)=F_in(:)


    end subroutine  

    subroutine Stopping_Power_Initialization(orbitals, stopping, density,  &
        system, all_PAW, atoms, elements, parallel, grids, density_t0, orbitals_t0)
    use stopping_type, only : stopping_struct
    use td_type, only : td_struct
    use atom_type, only : atom_struct
    use element_type, only : element_struct
    use system_type, only : system_struct
    use parallel_type, only : parallel_struct
    use odp_type, only : field_struct, orbital_struct, spin_DenPot_struct, stochastic
    use Local_Ion, only : calculate_local_ion_force
    use Non_Local_ion_EnFoSt, only : NL_PP_Force_RL
    use odp_type, only : orbital_struct
    use grids_type, only : grid_struct
    use Stochastic_Mod, only : stochastic_vector_restore
    use initialization_finalization, only : deterministic_vector_initializer, initialize_elements, initialize_atoms
    use odp, only: allocate_DenPot, allocate_orbital
    use simulation_type, only : all_PAW_struct

    type(parallel_struct), intent(in) :: parallel
    type(grid_struct), intent(inout) :: grids(:)
    type(spin_DenPot_struct), intent(in) :: density
    type(spin_DenPot_struct), intent(inout) ::  density_t0
    type(system_struct), intent(inout) :: system
    type(orbital_struct), intent(in) :: orbitals(:,:,:)
    type(orbital_struct), intent(inout), allocatable,target :: orbitals_t0(:,:,:)
    type(element_struct), intent(inout), allocatable :: elements(:)
    type(atom_struct), intent(inout), allocatable :: atoms(:)
    type(stopping_struct), intent(inout) :: stopping
    type(orbital_struct), pointer :: ks_orbitals_t0(:),  stoc_orbitals_t0(:)
    type(all_PAW_struct), intent(inout) :: all_PAW

    type(atom_struct), allocatable :: new_atoms(:)
    type(element_struct), allocatable :: new_elements(:)
    integer :: s, i, j, k
    real(dp), allocatable :: dummy_1(:,:,:), dummy_2(:,:,:)
    
    allocate(stopping%density_t(size(density%of)))
    if(parallel%myid.eq.0) then
        open(newunit=stopping%output_unit,  file="Stopping_Power.data", status="unknown")
        write(stopping%output_unit, *) "1:Time", '2-4:R(a.u.)', '5-7:P(a.u.)[3]', '8:SP(t)(a.u.)', '9:<SP> (a.u)', &
                                        '10:rho(a.u.)', '11:rho/<rho>', '12:F(a.u.)' 
    endif

    stopping%power_average=0.0_dp

    system%n_atoms=system%n_atoms+1
    system%n_elements=system%n_elements+1

    stopping%stop_at=system%n_atoms
    stopping%stop_el=system%n_elements

    !Concatenate the projectile onto the atoms and elements list 
    !Heres where you wish this was really python
    allocate(new_atoms(system%n_atoms))
    allocate(new_elements(system%n_elements))

    new_atoms(:(system%n_atoms-1))=atoms(:)
    new_elements(:(system%n_elements-1))=elements(:)

    new_atoms(stopping%stop_at)=stopping%atom(1)
    new_elements(stopping%stop_el)=stopping%element(1)

    deallocate(atoms)
    deallocate(elements)

    allocate(atoms(system%n_atoms))
    allocate(elements(system%n_elements))

    atoms(:)=new_atoms(:)
    elements(:)=new_elements(:)

    if(stopping%remove_t0) then

        density_t0%n_s=density%n_s
        call allocate_DenPot(density_t0, grids, 1, parallel)


        do s=1 , size(density_t0%of)
            density_t0%of(s)%R=density%of(s)%R
            density_t0%of(s)%G=density%of(s)%G
        enddo

        if(elements(stopping%stop_el)%n_proj.gt.0) then
            allocate(orbitals_t0( system%n_orbitals_local, &
                                  system%n_kpoints_local, &
                                  system%n_spin_local))
                          
            do i = 1, size(orbitals_t0,3);do j = 1, size(orbitals_t0,2);do k = 1, size(orbitals_t0,1);
                orbitals_t0(k,j,i)%band=orbitals(k,j,i)%band
                orbitals_t0(k,j,i)%k_point=orbitals(k,j,i)%k_point
                orbitals_t0(k,j,i)%spin=orbitals(k,j,i)%spin
                orbitals_t0(k,j,i)%degeneracy=orbitals(k,j,i)%degeneracy
                orbitals_t0(k,j,i)%n_spinor=orbitals(k,j,i)%n_spinor
                call allocate_orbital(orbitals_t0(k,j,i), grids, orbitals(k,j,i)%k_point+2, parallel)
                orbitals_t0(k,j,i)%occ=orbitals(k,j,i)%occ
                orbitals_t0(k,j,i)%weight=orbitals(k,j,i)%weight
                orbitals_t0(k,j,i)%filter=orbitals(k,j,i)%filter

                orbitals_t0(k,j,i)%type=orbitals(k,j,i)%type
                if(orbitals_t0(k,j,i)%type.eq.stochastic) then
                    if(allocated(orbitals_t0(i,k,s)%seed_array)) deallocate(orbitals_t0(i,k,s)%seed_array)
                    allocate(orbitals_t0(i,k,s)%seed_array(size(orbitals(i,k,s)%seed_array)))
                    orbitals_t0(i,k,s)%seed_array=orbitals(i,k,s)%seed_array
                endif

            enddo; enddo;enddo

            do i = 1, size(orbitals_t0,3);do j = 1, size(orbitals_t0,2);
                if(system%n_deterministic_local.gt.0) then
                    ks_orbitals_t0(1:system%n_deterministic_local)=>&
                    orbitals_t0(system%deterministic_start:system%deterministic_end, j, i)

                    call deterministic_vector_initializer(ks_orbitals_t0(:),grids, parallel)
                endif
                
                if(system%n_stochastic_local.gt.0) then
                    stoc_orbitals_t0(1:system%n_stochastic_local)=> &
                        orbitals_t0(system%stochastic_start:system%stochastic_end, j, i)

                    call stochastic_vector_restore(stoc_orbitals_t0(:), system, grids, parallel, &
                                dummy_1(:,:,:), dummy_2(:,:,:), all_PAW, atoms, elements)
                endif                      
            enddo; enddo

            do i = 1, size(orbitals,3);do j = 1, size(orbitals,2);do k = 1, system%n_orbitals_local
                do s = 1, size(orbitals(k,j,i)%of)
                    orbitals_t0(k,j,i)%of(s)%R=orbitals(k,j,i)%of(s)%R
                    orbitals_t0(k,j,i)%of(s)%G=orbitals(k,j,i)%of(s)%G
                enddo
            enddo;enddo;enddo
        endif
    endif

end subroutine 

end module
