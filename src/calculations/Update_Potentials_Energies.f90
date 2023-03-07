module Update_Potentials_and_Energies  

use types, only : dp
    use constants, only : pi
    implicit none
    
    
    public :: Update_ion_and_TD_Potentials, Calculate_all_static_Energies, &
    Update_electronic_from_orbtials_or_density_TD, fine_to_coarse_all_potentials, Update_BOMD_Potentials

    contains

    subroutine Update_ion_and_TD_Potentials(grids, system, parallel, potentials, td, atoms, elements, all_PAW, fine_to_coarse)
        use grids_type, only : grid_struct, inter_grid_struct
        use system_type, only : system_struct
        use td_type, only : td_struct
        use xc_type, only : xc_struct
        use odp_type, only : field_struct, orbital_struct, spin_DenPot_struct
        use simulation_type, only : all_PAW_struct, potentials_struct, energies_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use parallel_type, only : parallel_struct

        use Local_ion, only : calculate_local_ion_potential
        use External_Field, only : External_Field_Potential
        use PAW_Pseudopotentials, only: calc_Dij


        type(parallel_struct), intent(in) :: parallel
        type(system_struct), intent(in) :: system
        type(potentials_struct), intent(inout) :: potentials 
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(td_struct), intent(inout) :: td
        type(grid_struct), intent(inout), target :: grids(:)
        type(inter_grid_struct) ::  fine_to_coarse 
        type(all_PAW_struct), intent(inout) :: all_PAW
        integer :: s

        if(.not.all(atoms(:)%frozen_R)) call calculate_local_ion_potential(potentials%ion_local, system, atoms, elements, grids)
        call External_Field_Potential(td, potentials, grids)

        do s=1,size(potentials%total_local_fine%of)
            potentials%total_local_fine%of(s)%R=real(potentials%ion_local%of%R) + &
                                            potentials%hartree%of%R + &
                                            potentials%external_field%of%R + &
                                            potentials%xc%of(s)%R
            if(all_PAW%usepotzero.gt.0)  potentials%total_local_fine%of(s)%R = &
                potentials%total_local_fine%of(s)%R + sum(all_PAW%vpotzero)                                             

            if(system%orbital_free) potentials%total_local_fine%of(s)%R=potentials%total_local_fine%of(s)%R + &
                                potentials%kinetic_local%of(s)%R + &
                                potentials%kinetic_dynamic%of(s)%R
        enddo

        call fine_to_coarse_all_potentials(potentials, grids, parallel, fine_to_coarse)
        call Calc_Dij(potentials, atoms, elements, all_PAW, grids, parallel)

    end subroutine

    subroutine Update_BOMD_Potentials(grids, system, parallel, fine_density, &
        potentials, energies, atoms, elements, xc, tf, all_PAW, fine_to_coarse)
        use grids_type, only : grid_struct, inter_grid_struct
        use system_type, only : system_struct
        use xc_type, only : xc_struct
        use odp_type, only : field_struct, spin_DenPot_struct
        use simulation_type, only : all_PAW_struct, potentials_struct, energies_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use parallel_type, only : parallel_struct
        use tf_type, only: tf_struct

        use Local_ion, only : calculate_local_ion_potential
        use xc_mod, only : calc_XC
        use Hartree, only : calc_Hartree_potential
        use External_Field, only : External_Field_Potential
        use Density, only : calc_density
        use Thomas_Fermi, only: calc_Thomas_Fermi_Potential
        use PAW_Pseudopotentials, only: calc_Dij

        type(parallel_struct), intent(in) :: parallel
        type(spin_DenPot_struct), intent(inout) :: fine_density
        type(system_struct), intent(in) :: system
        type(xc_struct), intent(inout) :: xc(:)
        type(potentials_struct), intent(inout) :: potentials 
        type(energies_struct), intent(inout) :: energies
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(all_PAW_struct), intent(inout) :: all_PAW
        type(tf_struct), intent(in) :: tf

        type(inter_grid_struct) ::  fine_to_coarse
        integer :: s

        if(system%orbital_free) then
             call calc_Thomas_Fermi_Potential(fine_density, potentials, system, tf)
        else
            do s=1,size(potentials%kinetic_local%of)
                potentials%kinetic_local%of(s)%G=0.0_dp
                potentials%kinetic_local%of(s)%R=0.0_dp
            enddo
        endif
        call calc_Hartree_potential(potentials%hartree, fine_density, grids)
        call calc_XC(fine_density, xc, all_PAW, parallel, grids, energies%xc, potentials%xc)
        call calculate_local_ion_potential(potentials%ion_local, system, atoms, elements, grids)
        do s=1,size(potentials%total_local_fine%of)
            potentials%total_local_fine%of(s)%R=real(potentials%ion_local%of%R) + &
                                            potentials%hartree%of%R + &
                                            potentials%xc%of(s)%R
            if(all_PAW%usepotzero.gt.0)  potentials%total_local_fine%of(s)%R = &
                                            potentials%total_local_fine%of(s)%R + sum(all_PAW%vpotzero)  
            if(system%orbital_free) potentials%total_local_fine%of(s)%R=potentials%total_local_fine%of(s)%R + &
                                                                    potentials%kinetic_local%of(s)%R
        enddo
        call fine_to_coarse_all_potentials(potentials, grids, parallel, fine_to_coarse)
        call Calc_Dij(potentials, atoms, elements, all_PAW, grids, parallel)

    end subroutine

    subroutine Calculate_all_static_Energies(orbitals, grids, system, parallel, fine_density, &
        potentials, energies, atoms, elements, xc, tf, all_PAW, out_unit, step)
        use grids_type, only : grid_struct, inter_grid_struct
        use system_type, only : system_struct
        use xc_type, only : xc_struct
        use odp_type, only : field_struct, orbital_struct, spin_DenPot_struct
        use simulation_type, only : potentials_struct
        use simulation_type, only : energies_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use parallel_type, only : parallel_struct
        use tf_type, only: tf_struct

        use Local_ion, only : Calc_Local_PP_Energy,Calc_Core_Energy
        use Kinetic, only : Kinetic_Energy
        use xc_mod, only : calc_XC
        use Hartree, only : Hartree_Energy
        use External_Field, only : External_Field_Energy
        use Non_Local_ion_EnFoSt, only: NL_PP_Energy_RL
        use Thomas_Fermi, only: Thomas_Fermi_Energy
        use simulation_type, only : all_PAW_struct

        type(parallel_struct), intent(inout) :: parallel
        type(spin_DenPot_struct), intent(inout) :: fine_density
        type(system_struct), intent(in) :: system
        type(xc_struct), intent(inout) :: xc(:)
        type(tf_struct), intent(inout) :: tf
        type(orbital_struct), intent(inout) :: orbitals(:,:,:)
        type(potentials_struct), intent(inout) :: potentials 
        type(energies_struct), intent(inout) :: energies
        type(element_struct), intent(inout) :: elements(:)
        type(atom_struct), intent(inout) :: atoms(:)
        type(all_PAW_struct), intent(inout) :: all_PAw

        type(grid_struct), intent(inout), target :: grids(:)
        integer, intent(in) :: out_unit, step

        energies%hartree=Hartree_Energy(potentials%hartree, fine_density, parallel, grids)
        if(parallel%myid.eq.0) print *, 'Hartree Energy (a.u.): ', energies%hartree
        if(parallel%myid.eq.0) write(out_unit,*) 'Hartree Energy (a.u.): ', energies%hartree

        call calc_XC(fine_density, xc, all_PAW, parallel, grids, energies%xc)
        if(parallel%myid.eq.0) print *, 'XC (a.u.): ' ,  energies%xc
        if(parallel%myid.eq.0) write(out_unit,*) 'XC (a.u.): ' ,  energies%xc

        energies%ion_local=Calc_Local_PP_Energy(potentials%ion_local, fine_density, grids, parallel)
        if(parallel%myid.eq.0) print *, 'Local Ion (a.u.): ', energies%ion_local
        if(parallel%myid.eq.0) write(out_unit,*) 'Local Ion (a.u.): ', energies%ion_local

        energies%ion_non_local=NL_PP_Energy_RL(orbitals, atoms, elements, grids, parallel)
        if(parallel%myid.eq.0) print *, 'Non Local Ion (non-PAW) (a.u.): ', energies%ion_non_local
        if(parallel%myid.eq.0) write(out_unit,*) 'Non Local Ion (non-PAW) (a.u.): ', energies%ion_non_local

        energies%core=Calc_Core_Energy(system, elements, grids, all_PAW)
        if(parallel%myid.eq.0) print *, 'Core (a.u.): ', energies%core
        if(parallel%myid.eq.0) write(out_unit,*) 'Core (a.u.): ', energies%core

        energies%kinetic=Kinetic_Energy(orbitals, parallel, grids)
        energies%kinetic_local=0.0_dp
        if(parallel%myid.eq.0) print *, 'Kinetic Energy (a.u.): ', energies%kinetic
        if(parallel%myid.eq.0) write(out_unit,*) 'Kinetic Energy (a.u.): ', energies%kinetic

        if(system%orbital_free) then
            energies%kinetic_local=Thomas_Fermi_Energy(fine_density, &
             grids, parallel, system, tf)
             if(parallel%myid.eq.0) print *, 'Local Kinetic Energy (a.u.): ', energies%kinetic_local
             if(parallel%myid.eq.0) write(out_unit,*) 'Local Kinetic Energy (a.u.): ', energies%kinetic_local
        endif
        energies%external_field=External_Field_Energy(potentials%external_field, fine_density, grids, parallel)
        if(parallel%myid.eq.0) print *, 'External_Field (a.u.): ' ,  energies%external_field
        if(parallel%myid.eq.0) write(out_unit,*) 'External_Field (a.u.): ' ,  energies%external_field

        if(parallel%myid.eq.0) print *, 'Nuclear Potential (a.u.): ' ,  energies%nuclear_pot
        if(parallel%myid.eq.0) write(out_unit,*) 'Nuclear Potential (a.u.): ' ,  energies%nuclear_pot

        energies%nuclear_kinetic= &
            0.5_dp*sum((atoms(:)%P(1)**2+atoms(:)%P(2)**2+atoms(:)%P(3)**2)/elements(atoms(:)%element)%M)

        if(parallel%myid.eq.0) print *, 'Nuclear Kinetic(a.u.): ' , energies%nuclear_kinetic
        if(parallel%myid.eq.0) write(out_unit,*) 'Nuclear Kinetic(a.u.): ' , energies%nuclear_kinetic

        energies%total=energies%external_field + energies%kinetic_local + energies%kinetic + &
        energies%core + energies%ion_non_local + energies%ion_local +energies%xc  + energies%hartree + &
        energies%nuclear_kinetic + energies%nuclear_pot
        
        if(all_PAW%N_PAW_atoms.gt.0) then
            if(parallel%myid.eq.0) print *, 'PAW on-site (a.u.): ' , all_PAW%epaw
            if(parallel%myid.eq.0) write(out_unit,*) 'PAW on-site (a.u.): ' , all_PAW%epaw
            energies%total=energies%total+all_PAW%epaw

        endif

        if(parallel%myid.eq.0) print *, 'Total Energy (a.u.): ' , energies%total
        if(parallel%myid.eq.0) write(out_unit,*)  'Total Energy (a.u.): ' , energies%total
        if(parallel%myid.eq.0) print *, 'Total Ion Potential Energy : ' , energies%total-energies%nuclear_kinetic
        if(parallel%myid.eq.0) write(out_unit,*) 'Total Ion Potential Energy : ' , energies%total-energies%nuclear_kinetic
        ! ENTROPY

    end subroutine

    subroutine Update_electronic_from_orbtials_or_density_TD(grids, system, parallel, coarse_density, fine_density, &
        current_density, potentials, energies, atoms, elements, xc, tf, all_PAW, fine_to_coarse, coarse_to_fine, orbitals)
        use grids_type, only : grid_struct, inter_grid_struct
        use system_type, only : system_struct
        use xc_type, only : xc_struct
        use tf_type, only : tf_struct
        use odp_type, only : field_struct, orbital_struct, spin_DenPot_struct
        use simulation_type, only : all_PAW_struct, potentials_struct, energies_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use parallel_type, only : parallel_struct

        use Local_ion, only : Calc_Local_PP_Energy, Calc_Core_Energy
        use Kinetic, only : Kinetic_Energy
        use xc_mod, only : calc_XC
        use Hartree, only : calc_Hartree_potential, Hartree_Energy
        use External_Field, only : External_Field_Energy
        use Density, only : calc_density, calc_Compensation_Charge_Density
        use Thomas_Fermi, only: Thomas_Fermi_Energy, calc_Thomas_Fermi_Potential, calc_Dynamic_KEDP_and_energy
        use PAW_Pseudopotentials, only: calc_Dij
        use simulation_type, only : all_PAW_struct
        use Non_Local_ion_EnFoSt, only : NL_PAW_Density_matrix, NL_PP_Energy_RL

        type(parallel_struct), intent(inout) :: parallel
        type(spin_DenPot_struct), intent(inout) :: fine_density, coarse_density, current_density(3)
        type(system_struct), intent(in) :: system
        type(xc_struct), intent(inout) :: xc(:)
        type(tf_struct), intent(inout) :: tf
        type(orbital_struct), intent(inout), optional :: orbitals(:,:,:)
        type(potentials_struct), intent(inout) :: potentials 
        type(energies_struct), intent(inout) :: energies
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(all_PAW_struct), intent(inout) :: all_PAW
        type(inter_grid_struct) ::  fine_to_coarse, coarse_to_fine
        integer :: s

        if(all_PAW%N_PAW_atoms.gt.0) then
            if(present(orbitals)) call  NL_PAW_Density_matrix(orbitals(:,:,:), atoms, elements, grids, parallel, all_PAW) 
            if(present(orbitals)) call calc_Compensation_Charge_Density(all_PAW, grids, atoms, parallel)               
        endif   

        !Calc Density & new potentials just in case they are not already updated
        if(present(orbitals)) call calc_density(coarse_density, fine_density,  &
            orbitals, grids, parallel, coarse_to_fine, all_PAW)
      
        energies%kinetic_dynamic=0.0_dp
        if(system%orbital_free) then
                call calc_Thomas_Fermi_Potential(fine_density, potentials, system, tf)
                if(tf%dynamic_kedp) then
                    call calc_Dynamic_KEDP_and_energy(potentials%kinetic_dynamic, energies%kinetic_dynamic, &
                        current_density, fine_density, system, grids, parallel)
                    !No need to put on coarse because it is only for orbital Free which only uses dense grid
                endif
        else
            do s=1,size(potentials%kinetic_local%of)
                potentials%kinetic_local%of(s)%G=0.0_dp
                potentials%kinetic_dynamic%of(s)%G=0.0_dp
                potentials%kinetic_local%of(s)%R=0.0_dp
                potentials%kinetic_dynamic%of(s)%R=0.0_dp
            enddo
        endif
        call calc_Hartree_potential(potentials%hartree, fine_density, grids)
        call calc_XC(fine_density, xc, all_PAW, parallel, grids, energies%xc, potentials%xc)

        if(present(orbitals)) then !Don't calculate energies on potential iterations

            energies%hartree=Hartree_Energy(potentials%hartree, fine_density, parallel, grids)
            energies%ion_local=Calc_Local_PP_Energy(potentials%ion_local, fine_density, grids, parallel)
            energies%core=Calc_Core_Energy(system, elements, grids, all_PAW)
            energies%kinetic=Kinetic_Energy(orbitals, parallel, grids)
            energies%external_field=External_Field_Energy(potentials%external_field, fine_density, grids, parallel)
            energies%ion_non_local=NL_PP_Energy_RL(orbitals, atoms, elements, grids, parallel)
            if(system%orbital_free) energies%kinetic_local=Thomas_Fermi_Energy(fine_density, &
                grids, parallel, system, tf)

            if(parallel%myid.eq.0.and.system%orbital_free) print *, 'Local Kinetic Energy (a.u.): ', energies%kinetic_local
            if(parallel%myid.eq.0) print *, 'Kinetic Energy (a.u.): ', energies%kinetic
            if(parallel%myid.eq.0) print *, 'Hartree Energy (a.u.): ', energies%hartree
            if(parallel%myid.eq.0) print *, 'Local Ion (a.u.): ', energies%ion_local
            if(parallel%myid.eq.0) print *, 'Non Local (non-PAW) Ion (a.u.): ', energies%ion_non_local
            if(parallel%myid.eq.0) print *, 'Core (a.u.): ', energies%core
            if(parallel%myid.eq.0) print *, 'XC (a.u.): ' ,  energies%xc
            if(parallel%myid.eq.0) print *, 'External_Field (a.u.): ' ,  energies%external_field
            if(system%orbital_free) then
                if(parallel%myid.eq.0) print *, 'Kinetic Local (a.u.): ', energies%kinetic_local
                if(parallel%myid.eq.0) print *, 'Dynamic KEDF Energy (a.u.): ' , energies%kinetic_dynamic
            endif

            energies%total=energies%external_field + energies%kinetic_local + energies%kinetic + &
            energies%core + energies%ion_non_local + energies%ion_local +energies%xc  + energies%hartree + &
            energies%kinetic_dynamic

            if(all_PAW%N_PAW_atoms.gt.0) then
                if(parallel%myid.eq.0) print *, 'PAW on-site / NL (a.u.): ' , all_PAW%epaw
                energies%total=energies%total+all_PAW%epaw
            endif
        endif
        do s=1,size(potentials%total_local_fine%of)
            potentials%total_local_fine%of(s)%R=real(potentials%ion_local%of%R) + &
                                            potentials%hartree%of%R + &
                                            potentials%external_field%of%R + &
                                            potentials%xc%of(s)%R
            if(all_PAW%usepotzero.gt.0)  potentials%total_local_fine%of(s)%R = &
                                            potentials%total_local_fine%of(s)%R + sum(all_PAW%vpotzero)  
            if(system%orbital_free) potentials%total_local_fine%of(s)%R=potentials%total_local_fine%of(s)%R + &
                                potentials%kinetic_local%of(s)%R + &
                                potentials%kinetic_dynamic%of(s)%R
        enddo

        call fine_to_coarse_all_potentials(potentials, grids, parallel, fine_to_coarse)
        call Calc_Dij(potentials, atoms, elements, all_PAW, grids, parallel)


    end subroutine
    
    subroutine fine_to_coarse_all_potentials(potentials, grids, parallel, fine_to_coarse)
        use grids_type, only : grid_struct, inter_grid_struct
        use grids_mod, only : Gspace_grid_to_grid_transfer
        use simulation_type, only : potentials_struct
        use parallel_type, only : parallel_struct
        use system_type, only : system_struct
        use fft, only: real_to_recip, recip_to_real

        type(potentials_struct), intent(inout) :: potentials 
        type(grid_struct), intent(inout) :: grids(:)
        type(inter_grid_struct), intent(inout) ::  fine_to_coarse
        type(parallel_struct), intent(in) :: parallel

        integer :: s

        do s=1, size(potentials%total_local_fine%of)
            call real_to_recip(potentials%total_local_fine%of(s), grids)
            potentials%total_local_fine%of(s)%G=potentials%total_local_fine%of(s)%G*grids(1)%cutden
            call recip_to_real(potentials%total_local_fine%of(s), grids)
            call Gspace_grid_to_grid_transfer(grids(1), grids(2), parallel, potentials%total_local_fine%of(s)%G, &
                            potentials%total_local%of(s)%G, fine_to_coarse)
            call recip_to_real(potentials%total_local%of(s), grids)
            potentials%total_local%of(s)%R=real(potentials%total_local%of(s)%R)
        enddo

    end subroutine

end module