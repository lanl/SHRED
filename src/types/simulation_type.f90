module simulation_type
    use types, only : dp
    use element_type, only : element_struct
    use system_type, only : system_struct
    use atom_type, only : atom_struct
    use odp_type, only : orbital_struct, spin_DenPot_struct, DenPot_struct
    use parallel_type, only : parallel_struct
    use stochastic_type, only : stochastic_struct
    use stopping_type, only : stopping_struct
    use numerics_type, only : numerics_struct, pulay_real_pointers
    use grids_type, only : grid_struct, inter_grid_struct
    use xc_type, only: xc_struct
    use tf_type, only: tf_struct
    use tasks_type, only: tasks_struct
    use shared_type, only: shared_struct
    use td_type, only: td_struct
    use KG_type, only : KG_struct
    use Density_of_states_type, only : Density_of_states_struct
    use ewald_type, only : ewald_struct
    use m_pawrad, only: pawrad_type
    use m_pawtab, only: pawtab_type
    use m_pawang, only: pawang_type
    use m_paw_an, only : paw_an_type
    use m_paw_ij, only : paw_ij_type
    use m_pawfgrtab,only : pawfgrtab_type
    use m_pawrhoij, only :pawrhoij_type
    use m_pawcprj, only : pawcprj_type


    
    implicit none
    
    private

    public simulation_struct, potentials_struct, energies_struct, stress_struct, all_PAW_struct

type all_PAW_struct
integer :: N_PAW_elements=0, N_PAW_atoms, my_atoms
!By element
type(pawrad_type), allocatable :: rad(:)
type(pawtab_type), allocatable :: tab(:)
!By atom
type(paw_an_type), allocatable :: an(:)
type(paw_ij_type), allocatable  :: my_ij(:), ij(:)
type(pawfgrtab_type), allocatable  :: fgrtab(:)
type(pawrhoij_type), allocatable  :: rhoij(:)
type(pawcprj_type), allocatable :: cprj(:,:)
real(dp), allocatable :: Rhoij_mix_Ri(:,:,:), Rhoij_mix_Fi(:,:,:), Rhoij_mix_0(:,:), &
                         Rhoij_mix_m1(:,:), Rhoij_mix_m2(:,:), Rhoij_mix_fm1(:,:) 
type(pulay_real_pointers), allocatable :: Rhoij_mix_R_i(:,:), Rhoij_mix_F_i(:,:) 

integer, allocatable :: typat(:), dimcprj(:), nattyp(:), all_at(:)
logical :: my_atmtab_allocated, paral_atom
!arrays
integer, pointer :: my_atmtab(:)
integer :: my_natom

integer, allocatable :: list(:)
type(pawang_type) :: ang
integer  :: icoulomb=0
integer  :: ixc, nspden, lnmax, spnorb
integer  :: xclevel
integer  :: pawxcdev=1,usewvl=0,usexcnhat=-1, usepotzero=0
real(dp) :: xc_denpos=1.0E-14
type(spin_DenPot_struct) :: tncore, tnvale,  rho_comp, gr_rho_comp(3)
real(dp), allocatable :: znucl(:), proj_overlap(:,:), S_inv_full(:,:)
real(dp) :: epaw, vpotzero(2)

endtype

type potentials_struct
    type(DenPot_struct) :: hartree
    type(spin_DenPot_struct) :: xc
    type(DenPot_struct) :: ion_local
    type(DenPot_struct) :: external_field
    type(spin_DenPot_struct) :: kinetic_local
    type(spin_DenPot_struct) :: kinetic_dynamic
    type(spin_DenPot_struct) :: total_local
    type(spin_DenPot_struct) :: total_local_fine
end type

type energies_struct
    real(dp) :: hartree
    real(dp) :: nuclear_pot, nuclear_kinetic
    real(dp) :: kinetic
    real(dp) :: kinetic_local, kinetic_dynamic
    real(dp) :: xc
    real(dp) :: ion_local
    real(dp) :: ion_non_local
    real(dp) :: core
    real(dp) :: external_field
    real(dp) :: total
    real(dp) :: comp(3,3)
    real(dp) :: xc_core(3,3)
    integer :: output
end type

type stress_struct
    real(dp) :: comp(3,3)
    real(dp) :: hartree(3,3)
    real(dp) :: nuclear_pot(3,3), nuclear_kinetic(3,3)
    real(dp) :: kinetic(3,3)
    real(dp) :: kinetic_local(3,3), kinetic_dynamic(3,3)
    real(dp) :: xc(3,3), xc_core(3,3)
    real(dp) :: ion_local(3,3)
    real(dp) :: ion_non_local(3,3)
    real(dp) :: core(3,3)
    real(dp) :: external_field(3,3)
    real(dp) :: total(3,3)
    integer :: output
end type


type simulation_struct
    type(system_struct) :: system
    type(parallel_struct) :: parallel
    type(numerics_struct) :: numerics
    type(grid_struct), allocatable :: grids(:)
    type(atom_struct), allocatable :: atoms(:)
    type(element_struct), allocatable :: elements(:)

    type(orbital_struct), allocatable :: orbitals(:,:,:), orbitals_t0(:,:,:), orbitals_find_min(:,:,:)

    type(stochastic_struct), allocatable :: stochastic(:,:)
    type(xc_struct), allocatable :: xc(:)
    type(tf_struct):: tf
    type(td_struct):: td
    type(ewald_struct) :: ewald
    type(stopping_struct):: stopping

    type(tasks_struct) :: tasks
    type(spin_DenPot_struct) :: density, coarse_density, density_t0, coarse_density_pred
    type(spin_DenPot_struct) :: current_density(3), coarse_current_density(3), &
                                 coarse_current_density_pred(3)
    type(potentials_struct) :: potentials
    type(energies_struct) :: energies, pressures
    type(stress_struct) :: stress
    type(KG_struct) :: KG
    type(all_PAW_struct) :: all_PAW
    type(Density_of_states_struct) :: dos

    type(inter_grid_struct) :: coarse_to_fine, fine_to_coarse
    real(dp), allocatable :: stoc_norm(:,:,:), stoc_occ(:,:,:)
    real(dp), allocatable :: all_eigs(:,:,:), all_occ(:,:,:), all_filter(:,:,:), all_weight(:,:), Minimum_eigenvalue(:,:)
    real(dp) :: entropy
    integer, allocatable :: all_degeneracy(:,:,:)
    complex(dp), allocatable :: DPM(:,:,:,:,:)

    integer :: current_output_unit, positions_output_1, positions_output_2, velocities_output, forces_output
    integer :: trajectory_output, main_output

    real(dp) :: wall_time_start, wall_time


endtype

end module