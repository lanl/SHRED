module Stress_Pressure
    use types, only : dp

    implicit none

    public :: Stress_and_Pressure


    contains

    subroutine Stress_and_Pressure(orbitals, grids, system, parallel, density_in, &
        potentials, energies, stress, pressures, atoms, elements, xc, tf, all_PAW, output, step)
        use constants ,only : au2GPa
        use grids_type, only : grid_struct
        use system_type, only : system_struct
        use xc_type, only : xc_struct
        use odp_type, only : field_struct, orbital_struct, spin_DenPot_struct
        use simulation_type, only : potentials_struct, energies_struct, stress_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use parallel_type, only : parallel_struct
        use tf_type, only: tf_struct

        use Local_ion, only : Local_PP_Stress
        use Non_Local_ion_EnFoSt, only : NL_PP_Stress_RL
        use Kinetic, only : Kinetic_Stress
        use xc_mod, only : calc_XC
        use Hartree, only : Hartree_Stress
        use External_Field, only : External_Field_Energy
        use simulation_type, only : all_PAW_struct
        use Density, only: core_xc_Stress, Force_or_Stress_Comp_charge
        use Thomas_Fermi, only: Thomas_Fermi_Energy, calc_Thomas_Fermi_Potential
        use operations_3D, only : integrate_3D_R

        type(parallel_struct), intent(inout) :: parallel
        type(system_struct), intent(in) :: system
        type(spin_DenPot_struct), intent(inout) :: density_in
        type(xc_struct), intent(inout) :: xc(:)
        type(orbital_struct), intent(inout) :: orbitals(:,:,:)
        type(potentials_struct), intent(inout) :: potentials 
        type(energies_struct), intent(inout) :: energies, pressures
        type(stress_struct), intent(inout) :: stress
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(all_PAW_struct), intent(inout) :: all_PAw
        type(tf_struct), intent(inout) :: tf
        integer, intent(in) :: output, step
        integer :: dir,i

        if(parallel%myid.eq.0) print *,  'Stress and Pressure'
        if(parallel%myid.eq.0) write(output,*)  'Stress and Pressure'

        stress%kinetic_local=0.0_dp
        if(system%orbital_free) then
            grid=>grids(density_in%of(1)%grid)
            call calc_Thomas_Fermi_Potential(density_in, potentials, system, tf)
            stress%kinetic_local(1,1)=Thomas_Fermi_Energy(density_in, grids, parallel, system, tf)
            do i=1,density_in%n_s
                stress%kinetic_local(1,1)=stress%kinetic_local(1,1) - &
                real(integrate_3D_R(potentials%kinetic_local%of(i)%R*density_in%of(i)%R, grid, parallel))
            enddo     
            stress%kinetic_local(1,1)=stress%kinetic_local(1,1)
            stress%kinetic_local(2,2)=stress%kinetic_local(1,1)
            stress%kinetic_local(3,3)=stress%kinetic_local(1,1)
        endif

        stress%external_field=0.0_dp ! Probably dont stress/pressure from external field, but could add at some point
        stress%hartree=Hartree_Stress(potentials%hartree, density_in, parallel, grids)
        stress%ion_local=Local_PP_Stress(density_in, grids, atoms, elements, parallel)
        stress%ion_non_local=0.0_dp
        stress%xc_core=0.0_dp
        stress%comp=0.0_dp
        call NL_PP_Stress_RL(orbitals, atoms, elements, potentials, grids, parallel, all_PAW,stress%ion_non_local)

        stress%kinetic= Kinetic_Stress(orbitals, parallel, grids)
        if(system%orbital_free) then
            stress%kinetic=stress%kinetic*tf%lambda
        endif
        call calc_XC(density_in, xc, all_PAW, parallel, grids, energies%xc, potentials%xc, stress%xc)
        stress%core=0.0_dp

        stress%core(1,1)=-energies%core
        stress%core(2,2)=-energies%core
        stress%core(3,3)=-energies%core

        if(all_PAW%N_PAW_atoms.gt.0) then
            call Force_or_Stress_Comp_charge(all_PAW,  potentials, grids, atoms, parallel, &
                 calc_F=.false., calc_S=.true., stress=stress%comp)
            stress%xc_core=core_xc_Stress(potentials, grids, atoms, elements, parallel)
           ! stress%ion_non_local=stress%ion_non_local+stress%comp
        endif

        stress%nuclear_kinetic=0.0_dp
        do i=1, size(atoms)
            if(all(atoms(:)%frozen_R)) then
                do dir=1,3
                    stress%nuclear_kinetic(dir,dir)=stress%nuclear_kinetic(dir,dir)-system%temperature
                enddo
            else
                if(.not.atoms(i)%frozen_V) then
                do dir=1,3
                stress%nuclear_kinetic(dir,:)=stress%nuclear_kinetic(dir,:)-atoms(i)%P(dir)*atoms(i)%P(:) &
                    /elements(atoms(i)%element)%M
                enddo
                endif
            endif
        enddo

        stress%core = stress%core/product(grids(1)%box_length(:))
        stress%hartree = stress%hartree/product(grids(1)%box_length(:))
        stress%ion_local = stress%ion_local/product(grids(1)%box_length(:))
        stress%kinetic = stress%kinetic/product(grids(1)%box_length(:))
        stress%kinetic_local = stress%kinetic_local/product(grids(1)%box_length(:))
        stress%xc= stress%xc/product(grids(1)%box_length(:))
        stress%nuclear_pot= stress%nuclear_pot/product(grids(1)%box_length(:))
        stress%nuclear_kinetic= stress%nuclear_kinetic/product(grids(1)%box_length(:))
        stress%ion_non_local = stress%ion_non_local/product(grids(1)%box_length(:))

        stress%total=stress%core + stress%hartree + stress%ion_local + stress%ion_non_local &
                     + stress%kinetic + stress%kinetic_local +stress%xc
        !Dont forget to add the stress of the nuclei previously calculated in Ewald subroutine
        stress%total=stress%total+stress%nuclear_pot + stress%nuclear_kinetic
        if(all_PAW%N_PAW_atoms.gt.0) then
            stress%xc_core=stress%xc_core/product(grids(1)%box_length(:))
            stress%comp=stress%comp/product(grids(1)%box_length(:))
            stress%total= stress%total+stress%xc_core
            stress%total= stress%total+stress%comp
        endif

        if(parallel%myid.eq.0) then
            write(output,*) 'Step:', step

            print *, 'Total Stress: ',  stress%total(1,1),stress%total(2,2),stress%total(3,3), &
                                        stress%total(1,2), stress%total(1,3), stress%total(2,3)
            write(output,*) 'Total Stress: ',  stress%total(1,1),stress%total(2,2),stress%total(3,3), &
                                        stress%total(1,2), stress%total(1,3), stress%total(2,3)
            stress%total= stress%total-stress%nuclear_kinetic
            print *, 'Potential Stress', stress%total(1,1),stress%total(2,2),stress%total(3,3), &
                                        stress%total(1,2), stress%total(1,3), stress%total(2,3)
            write(output,*) 'Potential Stress: ',  stress%total(1,1),stress%total(2,2),stress%total(3,3), &
                                        stress%total(1,2), stress%total(1,3), stress%total(2,3)
            stress%total= stress%total+stress%nuclear_kinetic

           

            print *, '-Hartree Stress: ', stress%hartree(1,1),stress%hartree(2,2),stress%hartree(3,3), &
                                        stress%hartree(1,2), stress%hartree(1,3), stress%hartree(2,3)
            write(output,*) '-Hartree Stress: ', stress%hartree(1,1),stress%hartree(2,2),stress%hartree(3,3), &
                                        stress%hartree(1,2), stress%hartree(1,3), stress%hartree(2,3)
        
            print *, '-Local Ion Stress: ', stress%ion_local(1,1),stress%ion_local(2,2),stress%ion_local(3,3), &
                                        stress%ion_local(1,2), stress%ion_local(1,3), stress%ion_local(2,3)
            write(output,*) '-Local Ion Stress: ', stress%ion_local(1,1),stress%ion_local(2,2),stress%ion_local(3,3), &
                                        stress%ion_local(1,2), stress%ion_local(1,3), stress%ion_local(2,3)
                
            print *, '-Nonlocal Ion Stress: ', stress%ion_non_local(1,1),stress%ion_non_local(2,2),stress%ion_non_local(3,3), &
                                               stress%ion_non_local(1,2), stress%ion_non_local(1,3), stress%ion_non_local(2,3), &
                                               stress%ion_non_local(2,1), stress%ion_non_local(3,1), stress%ion_non_local(3,2)
            
            write(output,*) '-Nonlocal Ion Stress: ', &
                stress%ion_non_local(1,1),stress%ion_non_local(2,2),stress%ion_non_local(3,3), &
                stress%ion_non_local(1,2), stress%ion_non_local(1,3), stress%ion_non_local(2,3), &
                stress%ion_non_local(2,1), stress%ion_non_local(3,1), stress%ion_non_local(3,2)

            print *, '-Kinetic electron Stress: ', stress%kinetic(1,1),stress%kinetic(2,2),stress%kinetic(3,3), &
                                                    stress%kinetic(1,2), stress%kinetic(1,3), stress%kinetic(2,3)
            
            write(output,*) '-Kinetic electron Stress: ', stress%kinetic(1,1),stress%kinetic(2,2),stress%kinetic(3,3), &
                                                    stress%kinetic(1,2), stress%kinetic(1,3), stress%kinetic(2,3)

            print *, '-Local Kinetic electron Stress: ',  &
                                     stress%kinetic_local(1,1),stress%kinetic_local(2,2),stress%kinetic_local(3,3)

            write(output,*)  '-Local Kinetic electron Stress: ',  &
                                     stress%kinetic_local(1,1),stress%kinetic_local(2,2),stress%kinetic_local(3,3)
            print *, '-XC Stress: ', stress%xc(1,1),stress%xc(2,2),stress%xc(3,3), &
                                     stress%xc(1,2), stress%xc(1,3), stress%xc(2,3)
            write(output,*)  '-XC Stress: ', stress%xc(1,1),stress%xc(2,2),stress%xc(3,3), &
                                     stress%xc(1,2), stress%xc(1,3), stress%xc(2,3)

            print *, '-Ion core Stress: ', stress%core(1,1),stress%core(2,2),stress%core(3,3)
            write(output,*)  '-Ion core Stress: ', stress%core(1,1),stress%core(2,2),stress%core(3,3)

            if(all_PAW%N_PAW_atoms.gt.0) then
                print *, '-Core XC Stress: ', stress%xc_core(1,1),stress%xc_core(2,2),stress%xc_core(3,3), &
                                              stress%xc_core(1,2), stress%xc_core(1,3), stress%xc_core(2,3)
                write(output,*) '-Core XC Stress: ', stress%xc_core(1,1),stress%xc_core(2,2),stress%xc_core(3,3), &
                                              stress%xc_core(1,2), stress%xc_core(1,3), stress%xc_core(2,3)
                print *, '-Comp Charge Stress: ', stress%comp(1,1),stress%comp(2,2),stress%comp(3,3), &
                                              stress%comp(1,2), stress%comp(1,3), stress%comp(2,3), &
                                              stress%comp(2,1), stress%comp(3,1), stress%comp(3,2)
                write(output,*) '-Comp Charge Stress: ', stress%comp(1,1),stress%comp(2,2),stress%comp(3,3), &
                                              stress%comp(1,2), stress%comp(1,3), stress%comp(2,3), &
                                              stress%comp(2,1), stress%comp(3,1), stress%comp(3,2)
                                              
                print *, '-Comp + NL Stress: ', stress%comp(1,1) + stress%ion_non_local(1,1), &
                                                stress%comp(2,2) + stress%ion_non_local(2,2), &
                                                stress%comp(3,3) + stress%ion_non_local(3,3), &
                                                stress%comp(1,2) + stress%ion_non_local(1,2), &
                                                stress%comp(1,3) + stress%ion_non_local(1,3), &
                                                stress%comp(2,3) + stress%ion_non_local(2,3), &
                                                stress%comp(2,1) + stress%ion_non_local(2,1), &
                                                stress%comp(3,1) + stress%ion_non_local(3,1), &
                                                stress%comp(3,2) + stress%ion_non_local(3,2) 

                write(output,*) '-Comp + NL Stress: ', stress%comp(1,1) + stress%ion_non_local(1,1), &
                                                stress%comp(2,2) + stress%ion_non_local(2,2), &
                                                stress%comp(3,3) + stress%ion_non_local(3,3), &
                                                stress%comp(1,2) + stress%ion_non_local(1,2), &
                                                stress%comp(1,3) + stress%ion_non_local(1,3), &
                                                stress%comp(2,3) + stress%ion_non_local(2,3), &
                                                stress%comp(2,1) + stress%ion_non_local(2,1), &
                                                stress%comp(3,1) + stress%ion_non_local(3,1), &
                                                stress%comp(3,2) + stress%ion_non_local(3,2)                              
            endif
            
            print *, '-Nuclear Kinetic: ',  &
                stress%nuclear_kinetic(1,1),stress%nuclear_kinetic(2,2),stress%nuclear_kinetic(3,3), &
                stress%nuclear_kinetic(1,2), stress%nuclear_kinetic(1,3), stress%nuclear_kinetic(2,3)
            write(output,*) '-Nuclear Kinetic: ', &
                stress%nuclear_kinetic(1,1),stress%nuclear_kinetic(2,2),stress%nuclear_kinetic(3,3), &
                stress%nuclear_kinetic(1,2), stress%nuclear_kinetic(1,3), stress%nuclear_kinetic(2,3)
            print *, '-Nuclear Potential: ', stress%nuclear_pot(1,1),stress%nuclear_pot(2,2),stress%nuclear_pot(3,3), &
                                             stress%nuclear_pot(1,2), stress%nuclear_pot(1,3), stress%nuclear_pot(2,3)
            write(output,*) '-Nuclear Potential: ', &
                stress%nuclear_pot(1,1),stress%nuclear_pot(2,2),stress%nuclear_pot(3,3), &
                stress%nuclear_pot(1,2), stress%nuclear_pot(1,3), stress%nuclear_pot(2,3)
        endif


        pressures%kinetic_local=0.0_dp
        pressures%hartree=0.0_dp
        pressures%ion_local=0.0_dp
        pressures%ion_non_local=0.0_dp
        pressures%kinetic=0.0_dp
        pressures%xc=0.0_dp
        pressures%core=0.0_dp
        pressures%nuclear_pot=0.0_dp
        pressures%nuclear_kinetic=0.0_dp
        pressures%external_field=0.0_dp
        pressures%xc_core=0.0_dp
        pressures%total=0.0_dp

        do dir=1,3
            pressures%kinetic_local=pressures%kinetic_local - stress%kinetic_local(dir,dir)/3.0_dp
            pressures%hartree=pressures%hartree - stress%hartree(dir,dir)/3.0_dp
            pressures%ion_local=pressures%ion_local - stress%ion_local(dir,dir)/3.0_dp
            pressures%ion_non_local=pressures%ion_non_local - stress%ion_non_local(dir,dir)/3.0_dp
            pressures%kinetic=pressures%kinetic - stress%kinetic(dir,dir)/3.0_dp
            pressures%xc=pressures%xc - stress%xc(dir,dir)/3.0_dp
            pressures%core=pressures%core - stress%core(dir,dir)/3.0_dp
            pressures%nuclear_pot=pressures%nuclear_pot - stress%nuclear_pot(dir,dir)/3.0_dp
            pressures%nuclear_kinetic=pressures%nuclear_kinetic - &
                stress%nuclear_kinetic(dir,dir)/3.0_dp
            pressures%external_field=pressures%external_field - &
                stress%external_field(dir,dir)/3.0_dp
            pressures%total=pressures%total - stress%total(dir,dir)/3.0_dp
            if(all_PAW%N_PAW_atoms.gt.0) then
                pressures%xc_core=pressures%xc_core - stress%xc_core(dir,dir)/3.0_dp
                pressures%comp=pressures%comp - stress%comp(dir,dir)/3.0_dp
            endif
        enddo

        if(parallel%myid.eq.0) then
            print *, 'Total Pressure (GPA, a.u.): ', pressures%total*au2GPa , pressures%total
            pressures%total=pressures%total - pressures%nuclear_kinetic
            print *, 'Electronic & Nuc Pot Pressure: ', pressures%total*au2GPa , pressures%total
            pressures%total=pressures%total - pressures%nuclear_pot
            print *, 'Electronic Pressure: ', pressures%total*au2GPa , pressures%total
            pressures%total=pressures%total + pressures%nuclear_pot
            pressures%total=pressures%total + pressures%nuclear_kinetic

            print *, 'Nuclear Pressure: ', (pressures%nuclear_pot+pressures%nuclear_kinetic)*au2GPa , &
                                           (pressures%nuclear_pot+pressures%nuclear_kinetic)
            print *, '-Nuclear Potential: ', pressures%nuclear_pot*au2GPa, pressures%nuclear_pot 
            print *, '-Nuclear Kinetic: ', pressures%nuclear_kinetic*au2GPa, pressures%nuclear_kinetic

            print *, '-Hartree Pressure: ', pressures%hartree*au2GPa , pressures%hartree
            print *, '-Local Ion Pressure: ', pressures%ion_local*au2GPa , pressures%ion_local
            print *, '-Nonlocal Ion Pressure: ', pressures%ion_non_local*au2GPa , pressures%ion_non_local
            print *, '-Kinetic electron Pressure: ', pressures%kinetic*au2GPa , pressures%kinetic
            print *, '-Local Kinetic electron Pressure: ', pressures%kinetic_local*au2GPa , pressures%kinetic_local
            print *, '-XC Pressure: ', pressures%xc*au2GPa , pressures%xc
            print *, '-Ion core Pressure: ', pressures%core*au2GPa , pressures%core
            if(all_PAW%N_PAW_atoms.gt.0) then
                print *, '-Core XC Pressure: ', pressures%xc_core*au2GPa , pressures%xc_core
                print *, '-Comp Charge Pressure: ', pressures%comp*au2GPa , pressures%comp
                print *, '-Comp Charge + NL Pressure: ', (pressures%comp+pressures%ion_non_local)*au2GPa , &
                     pressures%comp+pressures%ion_non_local

            endif

            write(output,*)  'Total Pressure (GPA, a.u.): ', pressures%total*au2GPa , pressures%total
            pressures%total=pressures%total - pressures%nuclear_kinetic
            write(output,*)  'Electronic & Nuc Pot Pressure: ', pressures%total*au2GPa , pressures%total
            pressures%total=pressures%total - pressures%nuclear_pot
            write(output,*)  'Electronic Pressure: ', pressures%total*au2GPa , pressures%total
            pressures%total=pressures%total + pressures%nuclear_pot
            pressures%total=pressures%total + pressures%nuclear_kinetic

            write(output,*) 'Nuclear Pressure: ', (pressures%nuclear_pot+pressures%nuclear_kinetic)*au2GPa , &
                                           (pressures%nuclear_pot+pressures%nuclear_kinetic)
            write(output,*) '-Nuclear Potential: ', pressures%nuclear_pot*au2GPa, pressures%nuclear_pot
            write(output,*) '-Nuclear Kinetic: ', pressures%nuclear_kinetic*au2GPa, pressures%nuclear_kinetic

            write(output,*) '-Hartree Pressure: ', pressures%hartree*au2GPa , pressures%hartree
            write(output,*) '-Local Ion Pressure: ', pressures%ion_local*au2GPa , pressures%ion_local
            write(output,*) '-Nonlocal Ion Pressure: ', pressures%ion_non_local*au2GPa , pressures%ion_non_local
            write(output,*) '-Kinetic electron Pressure: ', pressures%kinetic*au2GPa , pressures%kinetic
            write(output,*) '-Local Kinetic electron Pressure: ', pressures%kinetic_local*au2GPa , pressures%kinetic_local
            write(output,*) '-XC Pressure: ', pressures%xc*au2GPa , pressures%xc
            write(output,*) '-Ion core Pressure: ', pressures%core*au2GPa , pressures%core
            if(all_PAW%N_PAW_atoms.gt.0) then
                write(output,*) '-Core XC Pressure: ', pressures%xc_core*au2GPa , pressures%xc_core
                write(output,*) '-Comp Charge Pressure: ', pressures%comp*au2GPa , pressures%comp
                write(output,*) '-Comp Charge + NL Pressure: ', (pressures%comp+pressures%ion_non_local)*au2GPa , &
                     pressures%comp+pressures%ion_non_local

            endif
            write(output,*) 'Nuclear Pressure: ', (pressures%nuclear_pot+pressures%nuclear_kinetic)*au2GPa , &
                                           (pressures%nuclear_pot+pressures%nuclear_kinetic)
            write(output,*) '-Nuclear Potential: ', pressures%nuclear_pot*au2GPa, pressures%nuclear_pot
            write(output,*) '-Nuclear Kinetic: ', pressures%nuclear_kinetic*au2GPa, pressures%nuclear_kinetic

            !write(output,*) 'External Field Pressure: ', pressures%external_field*au2GPa 
        endif
        
        flush(6)
        if(parallel%myid.eq.0) flush(output)
    end subroutine

end module
