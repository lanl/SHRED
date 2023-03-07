module allocate_memory
use types, only : dp
use xc_type

implicit none

private 

public :: allocate_required_memory_density_and_orbitals, allocate_required_memory_atoms_and_elements, &
          allocate_atomic_projectors, allocate_required_memory_TD

contains

subroutine allocate_required_memory_density_and_orbitals(sim)
    use simulation_type, only: simulation_struct
    use odp, only: allocate_DenPot, allocate_orbital
    use parallel_mod, only : parallel_task
    use shared_type, only : allocate_shared
    use Non_local_ion, only : allocate_and_fill_Ono_Hirose_matrix
    use odp_type, only : deterministic, stochastic, lowest, full, find_min, eigen_vector
    type(simulation_struct), intent(inout) :: sim
    integer :: i, j,k, n_min

    !Allocate the Electron density
    sim%coarse_density%n_s=sim%density%n_s
    sim%potentials%xc%n_s=sim%density%n_s
    sim%potentials%kinetic_local%n_s=sim%density%n_s
    sim%potentials%total_local%n_s=sim%density%n_s
    sim%potentials%total_local_fine%n_s=sim%density%n_s


    allocate(sim%Minimum_eigenvalue(sim%system%n_kpoints,sim%system%n_spin))

    call allocate_DenPot(sim%density, sim%grids, 1, sim%parallel)
    call allocate_DenPot(sim%coarse_density, sim%grids, 2, sim%parallel)

    !Allocate the potentials
    !Real-space potentials defined by density (on fine grid)
    call allocate_DenPot(sim%potentials%hartree, sim%grids, 1, sim%parallel)
    call allocate_DenPot(sim%potentials%xc, sim%grids, 1, sim%parallel)
    call allocate_DenPot(sim%potentials%ion_local, sim%grids, 1, sim%parallel)
    call allocate_DenPot(sim%potentials%kinetic_local, sim%grids, 1, sim%parallel)
    call allocate_DenPot(sim%potentials%external_field, sim%grids, 1, sim%parallel)
    call allocate_DenPot(sim%potentials%total_local_fine, sim%grids, 1, sim%parallel)
    !Collected total local Real-space potential which multiplies WF's
    call allocate_DenPot(sim%potentials%total_local, sim%grids, 2, sim%parallel)
    

    do i=1,sim%density%n_s
        sim%potentials%xc%of(i)%R=0.0_dp
        sim%potentials%kinetic_local%of(i)%R=0.0_dp
        sim%potentials%total_local%of(i)%R=0.0_dp
        sim%potentials%total_local_fine%of(i)%R=0.0_dp
    enddo

    sim%potentials%hartree%of%R=0.0_dp
    sim%potentials%ion_local%of%R=0.0_dp
    sim%potentials%external_field%of%R=0.0_dp

    do i=1,sim%density%n_s
        sim%potentials%xc%of(i)%G=0.0_dp
        sim%potentials%kinetic_local%of(i)%G=0.0_dp
        sim%potentials%total_local%of(i)%G=0.0_dp
        sim%potentials%total_local_fine%of(i)%G=0.0_dp
    enddo
    sim%potentials%hartree%of%G=0.0_dp
    sim%potentials%ion_local%of%G=0.0_dp
    sim%potentials%external_field%of%G=0.0_dp

    allocate(sim%all_weight(sim%system%n_kpoints, sim%system%n_spin ))
    do k=1, size( sim%all_weight,1)
        sim%all_weight(k,:)=sim%system%k_point(4,k)
    enddo

    if(sim%system%n_deterministic.gt.0) then
        allocate(sim%all_eigs  (sim%system%n_deterministic, sim%system%n_kpoints, sim%system%n_spin ))
        allocate(sim%all_occ   (sim%system%n_deterministic, sim%system%n_kpoints, sim%system%n_spin ))
        allocate(sim%all_filter(sim%system%n_deterministic, sim%system%n_kpoints, sim%system%n_spin ))
        allocate(sim%all_degeneracy(sim%system%n_deterministic, sim%system%n_kpoints, sim%system%n_spin ))
        
        sim%all_filter=0.0_dp
        sim%all_occ=0.0_dp
        sim%all_eigs=0.0_dp
        
        if(sim%system%spin_type.eq.1) then
            sim%all_degeneracy=2
        else if(sim%system%spin_type.eq.2) then
            sim%all_degeneracy=1
        else
            !not supported
            stop
        endif
    else
        allocate(sim%all_eigs  (0, sim%system%n_kpoints, sim%system%n_spin ))
        allocate(sim%all_occ   (0, sim%system%n_kpoints, sim%system%n_spin ))
        allocate(sim%all_filter(0, sim%system%n_kpoints, sim%system%n_spin ))
        allocate(sim%all_degeneracy(0, sim%system%n_kpoints, sim%system%n_spin ))
    endif

    if(sim%system%n_stochastic.gt.0) then 
        allocate(sim%stoc_norm  (sim%system%n_stochastic, sim%system%n_kpoints, sim%system%n_spin*sim%system%n_spinor ))
        allocate(sim%stoc_occ   (sim%system%n_stochastic, sim%system%n_kpoints, sim%system%n_spin*sim%system%n_spinor ))
    endif

    if(sim%parallel%my_sub.eq.0) then
        n_min=sim%system%n_find_min/sim%parallel%n_band_groups_sub
    else
        n_min=0
    endif

    allocate(sim%orbitals(sim%system%n_orbitals_local+n_min, &
                          sim%system%n_kpoints_local, &
                          sim%system%n_spin_local))

    do i=1, sim%system%n_spin_local  
        do j=1, sim%system%n_kpoints_local 
            do k=1, sim%system%n_orbitals_local+n_min
                !not all orbitals should need to have same  # spin components just lazy for now
                !When they are eventually set to be potentially different, 
                !then all orbitals in the same band group do need to have the same spin for distribution. 
                
                sim%orbitals(k,j,i)%n_spinor=sim%system%n_spinor
                !call parallel_task('max', sim%orbitals(k,j,i)%n_spinor, sim%parallel, 'band')   
            enddo
        enddo
    enddo

    allocate(sim%stochastic(sim%system%n_kpoints_local,sim%system%n_spin_local ))

    if(sim%system%n_orbitals_stochastic_full.gt.0) then
        sim%stochastic(:,:)%do_stoc=.true.
    else
        sim%stochastic(:,:)%do_stoc=.false.
    endif

    do i=1, sim%system%n_spin_local 
        do j=1, sim%system%n_kpoints_local
            do k=1, sim%system%n_orbitals_local+n_min
                sim%orbitals(k,j,i)%band=sim%parallel%my_band_group+(k-1)*sim%parallel%n_band_groups+1
                sim%orbitals(k,j,i)%k_point=sim%parallel%my_k_group+(j-1)*sim%parallel%n_k_groups+1
                sim%orbitals(k,j,i)%spin=sim%parallel%my_spin_group+(i-1)*sim%parallel%n_spin_groups+1
                if(sim%system%spin_type.eq.1) then
                    sim%orbitals(k,j,i)%degeneracy=2
                else if(sim%system%spin_type.eq.2) then
                    sim%orbitals(k,j,i)%degeneracy=1
                else
                    !not supported
                    stop
                endif
                if(sim%system%n_find_min.gt.0.and.sim%system%n_orbitals_deterministic_lowest.gt.0) then
                    print *, 'find_min orbitals allocated but already have deterministic lowest, this is a bug'
                    stop
                end if

                if(k.gt. sim%system%n_orbitals_local) then
                    sim%orbitals(k,j,i)%type=deterministic
                    sim%orbitals(k,j,i)%det_type=find_min
                    sim%orbitals(k,j,i)%td_type=eigen_vector
                    sim%orbitals(k,j,i)%band=sim%orbitals(k,j,i)%band-sim%system%n_orbitals
                else if(k.le.(sim%system%n_deterministic/sim%parallel%n_band_groups)) then
                    sim%orbitals(k,j,i)%type=deterministic
                    sim%orbitals(k,j,i)%td_type=eigen_vector
                    if(k.le.(sim%system%n_orbitals_deterministic_lowest/sim%parallel%n_band_groups)) then
                        sim%orbitals(k,j,i)%det_type=lowest
                    endif
                    sim%orbitals(k,j,i)%stoc_type=-1
                else
                    !Stochastic bands are independent of deterministic & find_min 
                    sim%orbitals(k,j,i)%band=sim%orbitals(k,j,i)%band-sim%system%n_deterministic
                    sim%orbitals(k,j,i)%type=stochastic 
                    sim%orbitals(k,j,i)%td_type=stochastic
                    if(k.le.(sim%system%n_orbitals_stochastic_full+sim%system%n_deterministic)/sim%parallel%n_band_groups) then
                        sim%orbitals(k,j,i)%stoc_type=full
                    endif
                    sim%orbitals(k,j,i)%det_type=-1
                endif
                call allocate_orbital(sim%orbitals(k,j,i), sim%grids, sim%orbitals(k,j,i)%k_point+2, sim%parallel)
                 sim%orbitals(k,j,i)%weight(:)=sim%all_weight(sim%orbitals(k,j,i)%k_point,sim%orbitals(k,j,i)%spin)
            enddo
        enddo
    enddo

end subroutine

subroutine allocate_required_memory_atoms_and_elements(sim)
    use simulation_type, only: simulation_struct
    use Non_local_ion, only : allocate_and_fill_Ono_Hirose_matrix, Calculate_Projector_GNL, &
                              Calculate_Projector_GNL_ORT,  allocate_and_fill_maps
    use odp, only : allocate_DenPot
    use grids_mod, only : allocate_local_fields_G
    use element_type, only : Reciprocal
    use Non_Local_ion, only : Calculate_Projector_GNL_ORT, Calculate_Projector_GNL
    type(simulation_struct), intent(inout) :: sim
    integer :: i, k

    call allocate_atomic_projectors(sim%atoms, sim%elements, sim%grids, sim%parallel, derivatives=.true.)
    if(sim%all_PAW%N_PAW_atoms.gt.0) then
        if(.not.allocated(sim%all_PAW%tncore%of)) then
            call allocate_DenPot(sim%all_PAW%tncore, sim%grids, 1, sim%parallel)
        endif
        if(.not.allocated(sim%all_PAW%tnvale%of)) then
            call allocate_DenPot(sim%all_PAW%tnvale, sim%grids, 1, sim%parallel)
        endif
        if(.not.allocated(sim%all_PAW%rho_comp%of)) then
            call allocate_DenPot(sim%all_PAW%rho_comp, sim%grids, 1, sim%parallel)
        endif
        do i=1,3
            if(.not.allocated(sim%all_PAW%gr_rho_comp(i)%of)) then
                call allocate_DenPot(sim%all_PAW%gr_rho_comp(i), sim%grids, 1, sim%parallel)
            endif
        enddo
        
    endif

    do i=1, size(sim%elements)
            if(allocated(sim%elements(i)%GNL)) cycle
            allocate(sim%elements(i)%GNL(size(sim%grids)))
            do k=1, size(sim%grids)
                call allocate_local_fields_G(sim%elements(i)%GNL(k)%projector_G, sim%grids(k), sim%elements(i)%n_proj)
                call Calculate_Projector_GNL(sim%elements(i), sim%grids(k), k)
                call allocate_local_fields_G(sim%elements(i)%GNL(k)%projector_G_ort, sim%grids(k), sim%elements(i)%n_proj)
                call Calculate_Projector_GNL_ORT(sim%elements(i), sim%grids(k), k)
            enddo
    enddo
    call allocate_and_fill_maps(sim%grids, sim%elements)
    call allocate_and_fill_Ono_Hirose_matrix(sim%grids, sim%elements, sim%parallel)

end subroutine

subroutine allocate_required_memory_TD(sim)
    use simulation_type, only: simulation_struct
    use odp, only: allocate_DenPot, allocate_orbital
    use odp_type, only : deterministic, lowest
    use element_type, only : PAW
    use td_type, only : TD_RT

    type(simulation_struct), intent(inout) :: sim
    integer :: i, j, k,e

    sim%potentials%kinetic_dynamic%n_s=sim%density%n_s
    sim%current_density%n_s=sim%density%n_s
    sim%coarse_current_density%n_s=sim%density%n_s
    call allocate_DenPot(sim%potentials%kinetic_dynamic, sim%grids, 1, sim%parallel)
    if(sim%tasks%calculate_current) then
        do i=1,3
            call allocate_DenPot(sim%current_density(i), sim%grids, 1, sim%parallel)
            call allocate_DenPot(sim%coarse_current_density(i), sim%grids, 2, sim%parallel)
        enddo
    endif

    if(sim%td%ETRS_steps.gt.0.and.sim%td%type.eq.TD_RT) then
        sim%coarse_density_pred%n_s=sim%density%n_s
        call allocate_DenPot(sim%coarse_density_pred, sim%grids, 2, sim%parallel)

        do i = 1, size(sim%atoms)
            e=sim%atoms(i)%element
            if(sim%elements(sim%atoms(i)%element)%PP_type.ne.PAW) cycle
            allocate(sim%atoms(i)%PAW%rho_ij_save(sim%elements(e)%PAW%tab%lmn_size,&
                sim%elements(e)%PAW%tab%lmn_size, sim%elements(e)%PAW%nspden))
        enddo

        if(sim%tf%dynamic_kedp) then 
            do i=1,3
                call allocate_DenPot(sim%coarse_current_density_pred(i), sim%grids, 2, sim%parallel)
            enddo
        endif
    endif

    if(sim%system%orbital_free) then; if(((sim%td%type.eq.TD_RT).or.(abs(sim%tf%lambda).gt.tiny(1.0_dp)))) then
        deallocate(sim%orbitals)
        allocate(sim%orbitals(1, &
        sim%system%n_kpoints_local, &
        sim%system%n_spin_local))


        do i=1, sim%system%n_spin_local 
            do j=1, sim%system%n_kpoints_local
                    sim%orbitals(1,j,i)%n_spinor=sim%system%n_spinor
                    sim%orbitals(1,j,i)%band=1
                    sim%orbitals(1,j,i)%type=deterministic
                    sim%orbitals(1,j,i)%det_type=lowest
                    sim%orbitals(1,j,i)%stoc_type=-1

                    sim%orbitals(1,j,i)%k_point=sim%parallel%my_k_group+(j-1)*sim%parallel%n_k_groups+1
                    sim%orbitals(1,j,i)%spin=sim%parallel%my_spin_group+(i-1)*sim%parallel%n_spin_groups+1
                    if(sim%system%spin_type.eq.1) then
                        sim%orbitals(1,j,i)%degeneracy=2
                    else if(sim%system%spin_type.eq.2) then
                        sim%orbitals(1,j,i)%degeneracy=1
                    else
                        !not supported
                        stop
                    endif
                    call allocate_orbital(sim%orbitals(1,j,i), sim%grids, sim%orbitals(1,j,i)%k_point+2, sim%parallel)
                    sim%orbitals(1,j,i)%weight(:)=sim%all_weight(sim%orbitals(1,j,i)%k_point,sim%orbitals(1,j,i)%spin)
            enddo
        enddo

    endif;endif

end subroutine

subroutine allocate_atomic_projectors(atoms, elements, grids, parallel, derivatives)
    use system_type, only : system_struct
    use atom_type, only : atom_struct
    use element_type, only : element_struct, Real_Local, PAW
    use parallel_type, only : parallel_struct
    use shared_type, only : allocate_shared
    use grids_type, only : grid_struct

    type(atom_struct), intent(inout) :: atoms(:)
    type(element_struct), intent(inout) :: elements(:)
    type(parallel_struct), intent(in) :: parallel
    type(grid_struct), intent(in) :: grids(:)
    logical, intent(in), optional :: derivatives
    integer, allocatable :: arrayshape(:)
    integer i,j
    logical :: derivatives_
    
    derivatives_=.false.
    if(present(derivatives)) derivatives_=derivatives

    do j=1, size(atoms)
        i=atoms(j)%element

        if(elements(i)%PP_type.eq.PAW) then
            allocate(atoms(j)%PAW%map_fine(5, &
                min(product(grids(1)%Nr_local(:)),product(elements(i)%PAW%max_sphere_fine(:)))))
            allocate(atoms(j)%PAW%Rs_fine(3, &
                min(product(grids(1)%Nr_local(:)),product(elements(i)%PAW%max_sphere_fine(:)))))
        endif

        if(elements(i)%PA_type.eq.Real_local) then
            if(allocated(atoms(j)%RL)) deallocate(atoms(j)%RL)
            allocate(atoms(j)%RL)
            atoms(j)%RL%ns=0
            if(elements(i)%n_proj.gt.0) then
                allocate(atoms(j)%RL%map(5, &
                    min(product(grids(2)%Nr_local(:)),product(elements(i)%RL%max_sphere(:)))))
                allocate(atoms(j)%RL%Rs(3, &
                    min(product(grids(2)%Nr_local(:)),product(elements(i)%RL%max_sphere(:)))))

                allocate(arrayshape(2))
                arrayshape=(/product(elements(i)%RL%max_sphere),elements(i)%n_proj/)
                call allocate_shared(atoms(j)%RL%projector, arrayshape, &
                                    atoms(j)%RL%proj_shared, parallel)
                deallocate(arrayshape)

                if(derivatives) then
                    allocate(arrayshape(3))
                    arrayshape=(/product(elements(i)%RL%max_sphere),3,elements(i)%n_proj/)
                    call allocate_shared(atoms(j)%RL%deriv_projector, arrayshape, &
                                        atoms(j)%RL%deriv_proj_shared, parallel)
                    deallocate(arrayshape)
                endif

                if(elements(i)%PP_type.eq.PAW) then
                    allocate(arrayshape(2))
                    arrayshape=(/product(elements(i)%RL%max_sphere),elements(i)%n_proj/)
                    call allocate_shared(atoms(j)%PAW%projector_ort, arrayshape, &
                                        atoms(j)%PAW%proj_shared_ort, parallel)
                    deallocate(arrayshape)
                    if(derivatives) then
                        allocate(arrayshape(3))
                        arrayshape=(/product(elements(i)%RL%max_sphere),3,elements(i)%n_proj/)
                        call allocate_shared(atoms(j)%PAW%deriv_projector_ort, arrayshape, &
                                            atoms(j)%PAW%deriv_proj_shared_ort, parallel)
                        deallocate(arrayshape)
                    endif
                endif
            endif
        endif
    enddo

    
end subroutine

end module