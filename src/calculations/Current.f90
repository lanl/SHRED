module Current
    use types, only : dp
    use constants , only : i_

    implicit none

    public ::  Current_Density_Calculation, orbital_current, orbitals_both_current, &
                 Dipole_Matrix_Elements, J_to_orbitals, J_Matrix_Elements
    private :: Add_to_local_Current_Density

    contains

    subroutine Current_Density_Calculation(orbitals, parallel, potentials,  atoms, elements, all_PAW, grids, &
                                           current_density, fine_current_density, coarse_to_fine)
        use td_type, only : td_struct
        use grids_type, only : grid_struct, inter_grid_struct
        use parallel_type, only : parallel_struct
        use odp_type, only : field_struct, orbital_struct, spin_DenPot_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only : recip_to_real, real_to_recip
        use parallel_mod, only : parallel_task
        use grids_mod, only: Gspace_grid_to_grid_transfer
        use simulation_type, only : potentials_struct, all_PAW_struct

        type(parallel_struct), intent(in) :: parallel
        type(orbital_struct), intent(inout) :: orbitals(:,:,:)
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(spin_DenPot_struct), intent(inout) :: current_density(:), fine_current_density(:)
        type(inter_grid_struct), intent(inout) ::  coarse_to_fine
        type(potentials_struct), intent(in) :: potentials 
        type(all_PAW_struct), intent(inout) :: all_PAW

        integer :: dir, i
        do dir=1,3;do i=1,current_density(dir)%n_s
                current_density(dir)%of(i)%R=0.0_dp
        enddo; enddo

        call Add_to_local_Current_Density(orbitals,  parallel, potentials, atoms, elements, all_PAW, grids, current_density)

        do dir=1,3; do i=1,current_density(dir)%n_s
            call parallel_task('sum', current_density(dir)%of(i)%R(:,:,:), parallel, 'space') 
            call real_to_recip(current_density(dir)%of(i), grids)
            current_density(dir)%of(i)%G=current_density(dir)%of(i)%G*grids(2)%cutden
            call recip_to_real(current_density(dir)%of(i), grids)
            call Gspace_grid_to_grid_transfer(grids(2), grids(1), parallel, current_density(dir)%of(i)%G, &
                  fine_current_density(dir)%of(i)%G, coarse_to_fine)
            call recip_to_real(fine_current_density(dir)%of(i), grids)
        enddo; enddo

    end subroutine


    subroutine Add_to_local_Current_Density(orbitals, parallel, potentials,  atoms, elements, all_PAW, grids, current_density)
        use constants, only : i_
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use odp_type, only : field_struct, orbital_struct, spin_DenPot_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use Non_Local_ion, only : Apply_projectors_RL, Calculate_Projector_overlaps
        use odp, only: allocate_orbital, deallocate_orbital
        use simulation_type, only : potentials_struct, all_PAW_struct

        type(parallel_struct), intent(in) :: parallel
        type(orbital_struct), intent(inout) :: orbitals(:,:,:)
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(spin_DenPot_struct), intent(inout) :: current_density(:)
        type(potentials_struct), intent(in) :: potentials 
        type(all_PAW_struct), intent(inout) :: all_PAW

        type(orbital_struct), allocatable :: J_orbital(:)

        integer :: i, j, k, dir, s

        allocate(J_orbital(3))

        
        do i=1,size(orbitals,3); do j=1,size(orbitals,2);do k=1,size(orbitals,1)
            if(sum(abs(orbitals(k,j,i)%occ(:))).lt.tiny(1.0_dp)) cycle

            do dir=1,3
                J_orbital(dir)%band=orbitals(k,j,i)%band
                J_orbital(dir)%k_point=orbitals(k,j,i)%k_point
                J_orbital(dir)%spin=orbitals(k,j,i)%spin
                J_orbital(dir)%degeneracy=orbitals(k,j,i)%degeneracy
                J_orbital(dir)%type=orbitals(k,j,i)%type
                J_orbital(dir)%n_spinor=orbitals(k,j,i)%n_spinor
                call allocate_orbital(J_orbital(dir), grids, orbitals(k,j,i)%k_point+2, parallel)
                J_orbital(dir)%occ=orbitals(k,j,i)%occ
                J_orbital(dir)%weight=orbitals(k,j,i)%weight
                J_orbital(dir)%filter=orbitals(k,j,i)%filter
            enddo  
    

            call orbital_current(orbitals(k,j,i), J_orbital, parallel, potentials, atoms, elements, all_PAW, &
                                 grids)

            do s=1,size(orbitals(k,j,i)%of)
                !i and s are both spin variables, only one should be > 1
                do dir=1,3
                    current_density(dir)%of(s*i)%R =  current_density(dir)%of(s*i)%R + &
                                             real(conjg(orbitals(k,j,i)%of(s)%R)*J_orbital(dir)%of(s)%R)*&
                                             orbitals(k,j,i)%occ(s)*orbitals(k,j,i)%weight(s)
                enddo
            enddo
            do dir=1,3
                call deallocate_orbital(J_orbital(dir),grids)
            enddo
        enddo;enddo;enddo

        deallocate(J_orbital)

            
    end subroutine

    subroutine orbital_current(orbital, P_orbital, parallel, potentials, atoms, elements, &
                                 all_PAW, grids, H_field_in, H_field_out)
        use constants, only : i_
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use odp_type, only : field_struct, orbital_struct, spin_DenPot_struct
        use element_type, only : element_struct, HGH, PAW
        use atom_type, only : atom_struct
        use Non_Local_ion, only : Apply_projectors_RL, Calculate_Projector_overlaps
        use fft, only : recip_to_real, real_to_recip
        use odp, only : allocate_field, deallocate_field
        use simulation_type, only : potentials_struct, all_PAW_struct
        use Apply_Hamiltonian, only : Apply_SinvH, Apply_HSinv

        type(parallel_struct), intent(in) :: parallel
        type(orbital_struct), intent(inout) :: orbital, P_orbital(:)
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(potentials_struct), intent(in) :: potentials 
        type(all_PAW_struct), intent(inout) :: all_PAW

        type(grid_struct), pointer :: grid
        type(field_struct), allocatable :: tmp_field(:,:), tmp_field_3(:,:)
        type(field_struct), allocatable, target :: H_field_tmp(:)
        type(field_struct), intent(inout), optional, target :: H_field_in(:), H_field_out(:)
        type(field_struct), pointer :: H_field(:)
        complex(dp),  allocatable :: Pypsi(:,:,:),dPypsi(:,:,:,:)
        integer :: dir, s
        

        if(present(H_field_in)) then
            H_field(1:size(H_field_in))=>H_field_in(1:size(H_field_in))
        else if(present(H_field_out)) then
            H_field(1:size(H_field_out))=>H_field_out(1:size(H_field_out))
            !Dont need H_psi
        endif

        allocate(tmp_field(3,size(orbital%of)))
        do s=1,size(tmp_field,2)
            do dir=1,3
                grid=>grids(orbital%of(s)%grid)
                call allocate_field(tmp_field(dir,s), grid, parallel)
                tmp_field(dir,s)%grid=orbital%of(s)%grid
                tmp_field(dir,s)%R=0.0_dp
                tmp_field(dir,s)%G=0.0_dp
                P_orbital(dir)%of(s)%R=0.0_dp
                P_orbital(dir)%of(s)%G=0.0_dp
            enddo
        enddo

        allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(orbital%of)))
        do s=1,size(tmp_field,2)
            do dir=1,3
                P_orbital(dir)%of(s)%R=0.0_dp
                P_orbital(dir)%of(s)%G=0.0_dp
            enddo
        enddo
        if(any(elements%PP_type.eq.HGH)) then
            !Phat= (G+k) i_*(-[r,V_NL]
            PYpsi=0.0_dp
            allocate(dPYpsi(3,maxval(elements(:)%n_proj), size(atoms), size(orbital%of)))
            dPYpsi=0.0_dp
            do s=1, size(orbital%of)
                call Calculate_Projector_overlaps(orbital%of(s),  &
                PYpsi(:,:,s), atoms, elements, grid, parallel, dPYpsi(:,:,:,s), skip_PAW=.true.) 
            enddo
            call Apply_projectors_RL(orbital%of, atoms, elements, &
            grid, PYpsi, dPYpsi, midrVpsi=tmp_field, skip_PAW=.true.)
            do s=1,size(orbital%of)
                do dir=1,3
                    P_orbital(dir)%of(s)%R=tmp_field(dir,s)%R
                enddo 
            enddo
            deallocate(dPYpsi)
        endif

        if(any(elements%PP_type.eq.PAW)) then
            PYpsi=0.0_dp               
                !Basic PAW P
                !Phat= (G+k) -i_*|Pi>(nabla_ij +i_*k*sij(no k, use phases))<Pj|
                do s=1, size(orbital%of)
                    call Calculate_Projector_overlaps(orbital%of(s),  &
                    PYpsi(:,:,s), atoms, elements, grid, parallel, skip_HGH=.true.) 
                enddo

                call Apply_projectors_RL(orbital%of, atoms, elements, &
                    grid, PYpsi, nabla_psi=tmp_field, skip_HGH=.true.) 
                do s=1,size(orbital%of)
                    do dir=1,3
                        P_orbital(dir)%of(s)%R=P_orbital(dir)%of(s)%R-tmp_field(dir,s)%R
                    enddo 
                enddo
        endif 

        deallocate(PYpsi)
        do s=1,size(tmp_field,2)
            do dir=1,3
                call deallocate_field(tmp_field(dir,s),grids(tmp_field(dir,s)%grid))
                if(allocated(tmp_field_3)) call deallocate_field(tmp_field_3(dir,s),grids(tmp_field_3(dir,s)%grid))
            enddo
            if(allocated(H_field_tmp))  call deallocate_field(H_field_tmp(s),grids(H_field_tmp(s)%grid))
        enddo
        if(associated(H_field)) nullify(H_field)
        if(allocated(H_field_tmp)) deallocate(H_field_tmp)

        deallocate(tmp_field)
        if(allocated(tmp_field_3)) deallocate(tmp_field_3)


        do s=1,size(orbital%of)
            grid=>grids(orbital%of(s)%grid)
            do dir=1,3
                if(any(elements%PP_type.eq.PAW).or.any(elements%PP_type.eq.HGH)) then
                    P_orbital(dir)%of(s)%R=i_*P_orbital(dir)%of(s)%R
                    call real_to_recip(P_orbital(dir)%of(s), grids)
                else
                    P_orbital(dir)%of(s)%G=0.0_dp
                endif

                P_orbital(dir)%of(s)%G=P_orbital(dir)%of(s)%G+grid%G(:,:,:,dir)*orbital%of(s)%G
                P_orbital(dir)%of(s)%G=P_orbital(dir)%of(s)%G*grid%cutwf
                call recip_to_real(P_orbital(dir)%of(s), grids)
            enddo 
        enddo

    end subroutine

    subroutine orbital_energy_current(orbital, P_orbital, EP_orbital, enthalpy, parallel, &
                                         potentials, atoms, elements, all_PAW, grids, H_orbital_)
        use constants, only : i_
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use odp_type, only : field_struct, orbital_struct, spin_DenPot_struct
        use element_type, only : element_struct, PAW
        use atom_type, only : atom_struct
        use Non_Local_ion, only : Apply_projectors_RL, Calculate_Projector_overlaps
        use fft, only : recip_to_real, real_to_recip
        use odp, only: allocate_orbital, deallocate_orbital
        use simulation_type, only : potentials_struct, all_PAW_struct
        use Apply_Hamiltonian, only : Apply_SinvH, Apply_HSinv, Apply_H

        type(parallel_struct), intent(in) :: parallel
        type(orbital_struct), intent(inout) :: orbital, P_orbital(:), EP_orbital(:)
        type(orbital_struct), intent(inout), optional, target :: H_orbital_
        type(orbital_struct), allocatable, target :: H_orbital_tmp
        type(orbital_struct), pointer :: H_orbital
        real(dp), intent(in) :: enthalpy

        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(potentials_struct), intent(in) :: potentials 
        type(all_PAW_struct), intent(inout) :: all_PAW

        type(orbital_struct), allocatable :: orbital_tmp(:)
        integer :: dir, s

        if(present(H_orbital_)) then
            H_orbital=>H_orbital_
        else  
            allocate(H_orbital_tmp)
            H_orbital_tmp%band=orbital%band
            H_orbital_tmp%k_point=orbital%k_point
            H_orbital_tmp%spin=orbital%spin
            H_orbital_tmp%degeneracy=orbital%degeneracy
            H_orbital_tmp%type=orbital%type
            H_orbital_tmp%n_spinor=orbital%n_spinor
            call allocate_orbital(H_orbital_tmp, grids, orbital%k_point+2, parallel)
            H_orbital_tmp%occ=orbital%occ
            H_orbital_tmp%weight=orbital%weight
            H_orbital_tmp%filter=orbital%filter
            if(any(elements%PP_type.eq.PAW)) then
                call Apply_SinvH(orbital%of, H_orbital_tmp%of, grids, potentials, atoms, elements, &
                parallel, all_PAW, calc_G=.false., with_CG=.false.)
            else
                call Apply_H(orbital%of, H_orbital_tmp%of, grids, potentials,&
                atoms, elements, parallel, calc_R=.false.)
                do s=1,size(H_orbital_tmp%of)
                    call recip_to_real(H_orbital_tmp%of(s), grids)
                enddo
            endif
            H_orbital=>H_orbital_tmp
        endif

        allocate(orbital_tmp(3))
        do dir=1,3
            orbital_tmp(dir)%band=orbital%band
            orbital_tmp(dir)%k_point=orbital%k_point
            orbital_tmp(dir)%spin=orbital%spin
            orbital_tmp(dir)%degeneracy=orbital%degeneracy
            orbital_tmp(dir)%type=orbital%type
            orbital_tmp(dir)%n_spinor=orbital%n_spinor
            call allocate_orbital(orbital_tmp(dir), grids, orbital%k_point+2, parallel)
            orbital_tmp(dir)%occ=orbital%occ
            orbital_tmp(dir)%weight=orbital%weight
            orbital_tmp(dir)%filter=orbital%filter
        enddo  

        do dir=1,3
            do s=1,size(EP_orbital(dir)%of)
                EP_orbital(dir)%of(s)%R=0.0_dp
                EP_orbital(dir)%of(s)%G=0.0_dp
            enddo
        enddo

         !P(S^-1)H
        call orbital_current(H_orbital, EP_orbital, parallel, potentials, atoms, elements, all_PAW, grids)

        !H(S^-1)P
        do dir=1,3
            if(any(elements%PP_type.eq.PAW)) then
                call Apply_HSinv(P_orbital(dir)%of(:), orbital_tmp(dir)%of, grids, potentials, atoms, elements, &
                parallel, all_PAW, calc_R=.false., with_CG=.false.)
            else
                call Apply_H(P_orbital(dir)%of(:), orbital_tmp(dir)%of, grids, potentials,&
                atoms, elements, parallel, calc_R=.false.)
                do s=1,size(orbital_tmp(dir)%of)
                    call recip_to_real(orbital_tmp(dir)%of(s), grids)
                enddo
            endif
            do s=1,size(orbital_tmp(dir)%of)
                EP_orbital(dir)%of(s)%R=EP_orbital(dir)%of(s)%R+orbital_tmp(dir)%of(s)%R
                EP_orbital(dir)%of(s)%G=EP_orbital(dir)%of(s)%G+orbital_tmp(dir)%of(s)%G
            enddo
        enddo

        do dir=1,3
            do s=1,size(EP_orbital(dir)%of)
                EP_orbital(dir)%of(s)%R=0.5_dp*EP_orbital(dir)%of(s)%R-P_orbital(dir)%of(s)%R*enthalpy
                EP_orbital(dir)%of(s)%G=0.5_dp*EP_orbital(dir)%of(s)%G-P_orbital(dir)%of(s)%G*enthalpy
            enddo
        enddo

        do dir=1,3
            call deallocate_orbital(orbital_tmp(dir),grids)
        enddo
        deallocate(orbital_tmp)
   
        if(allocated(H_orbital_tmp)) call deallocate_orbital(H_orbital_tmp,grids)
        if(allocated(H_orbital_tmp)) deallocate(H_orbital_tmp)
        if(associated(H_orbital)) nullify(H_orbital)

    end subroutine

    subroutine orbitals_both_current(orbitals, potentials, atoms, elements, parallel, all_PAW, grids, J_orbitals, enthalpy)

        use parallel_type, only : parallel_struct
        use odp_type, only : field_struct, orbital_struct
        use grids_type, only : grid_struct
        use simulation_type, only : potentials_struct, all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use odp, only: allocate_field, deallocate_field, allocate_orbital, deallocate_orbital

        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(orbital_struct), intent(inout) :: J_orbitals(:,:,:)
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(field_struct), allocatable ::  tmp_field_1(:)
        type(orbital_struct), allocatable :: tmp_orbital_1
        type(all_PAW_struct), intent(inout) :: all_PAW
        real(dp) , intent(in) ::  enthalpy

        integer :: s, k

        grid=>grids(orbitals(1)%of(1)%grid)
        
        !Orbital = | \chi >
        do k=1, size(orbitals)
            call  orbital_current(orbitals(k), J_orbitals(k,:,1), parallel, potentials, atoms, elements, all_PAW, grids)
            call orbital_energy_current(orbitals(k), J_orbitals(k,:,1), J_orbitals(k,:,2), &
                        enthalpy, parallel, potentials, atoms, elements, all_PAW, grids)
        enddo

    end subroutine


    subroutine Dipole_Matrix_Elements(orbitals_l, orbitals_r, &
        grids, parallel, potentials, atoms, elements, all_PAW, DPM, HPM, enthalpy)
        use constants, only : pi, i_
        use grids_type, only : grid_struct
        use odp_type, only : field_struct, orbital_struct, deterministic, stochastic
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: recip_to_real, real_to_recip
        use lapack, only: zgemm

        use parallel_mod,only : parallel_task, parallel_wait
        use parallel_type, only : parallel_struct
        use numerics_type, only : numerics_struct
        use Apply_Hamiltonian, only: Apply_H, Apply_SinvH, Apply_HSinv
        use FFT2Matrix, only : orbitals_FFT2Matrix, orbitals_Matrix2FFT, grids_G2_FFT2Matrix, &
                            orbital_FFT2Matrix
        use odp, only: allocate_field, deallocate_field, allocate_orbital, deallocate_orbital 
        use Non_Local_ion, only : Apply_projectors_RL, Calculate_Projector_overlaps 
        use simulation_type, only : all_PAW_struct


        type(parallel_struct), intent(in) :: parallel
        type(orbital_struct), intent(inout) :: orbitals_l(:), orbitals_r(:)
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(all_PAW_struct), intent(inout) :: all_PAW
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        complex(dp),intent(inout) :: DPM(:,:,:)
        complex(dp),intent(inout), optional ::  HPM(:,:,:)
        real(dp), intent(in), optional :: enthalpy
        logical :: symmetric
        complex(dp), allocatable :: Y(:,:), AY(:,:)
        integer :: nYl, nYr, nX, nband_block_l, nband_block_r, nG
        integer:: i
        integer :: orb_grid
        type(orbital_struct), allocatable:: J_orbitals_r(:,:,:)
        integer :: dir

        if(parallel%myid.eq.0) print *, 'Calculating Momentum (DPM) and/or Energy Momentum (HPM) Matrix Elements'

        grid=>grids(orbitals_l(1)%of(1)%grid)

        nband_block_l=size(orbitals_l)
        nX=parallel%n_band_groups
        nYl=nX*nband_block_l
        orb_grid=orbitals_l(1)%of(1)%grid
        if(orb_grid.ne.orbitals_r(1)%of(1)%grid) then
            print *, 'Trying to calculate Matrix element between different grids'
        endif

        nband_block_r=size(orbitals_r)
        nX=parallel%n_band_groups
        nYr=nX*nband_block_r

        do i=1,nband_block_l
            if(orbitals_l(i)%of(1)%grid.ne.orb_grid) then
                print *, 'All orbitals passed into Egensolver need to be on the same grid (same k-point)'
            endif
        enddo
        
        do i=1,nband_block_r
            if(orbitals_l(i)%of(1)%grid.ne.orb_grid) then
                print *, 'All orbitals passed into Egensolver need to be on the same grid (same k-point)'
            endif
        enddo


        grid=>grids(orb_grid)
        if(grid%reduction_option.eq.2) then
            nG=product(grid%Ng_small)
        else
            nG=product(grid%Ng_local)
        endif
        if(mod(nG,nX).ne.0) nG=(nG/nX+1)*nX

        if(present(HPM)) then
            allocate(J_orbitals_r(size(orbitals_r),3,2))
        else
            allocate(J_orbitals_r(size(orbitals_r),3,1))
        endif
        do dir=1,3
        do i = 1, size(orbitals_r)
            J_orbitals_r(i,dir,:)%band=orbitals_r(i)%band
            J_orbitals_r(i,dir,:)%k_point=orbitals_r(i)%k_point
            J_orbitals_r(i,dir,:)%spin=orbitals_r(i)%spin
            J_orbitals_r(i,dir,:)%degeneracy=orbitals_r(i)%degeneracy
            J_orbitals_r(i,dir,:)%type=orbitals_r(i)%type
            J_orbitals_r(i,dir,:)%n_spinor=orbitals_r(i)%n_spinor

            call allocate_orbital(J_orbitals_r(i,dir,1), grids, orbitals_r(i)%k_point+2, parallel)
            J_orbitals_r(i,dir,1)%occ=orbitals_r(i)%occ
            J_orbitals_r(i,dir,1)%weight=orbitals_r(i)%weight
            J_orbitals_r(i,dir,1)%filter=orbitals_r(i)%filter
            if(present(HPM)) then
                call allocate_orbital(J_orbitals_r(i,dir,2), grids, orbitals_r(i)%k_point+2, parallel)
                J_orbitals_r(i,dir,2)%occ=orbitals_r(i)%occ
                J_orbitals_r(i,dir,2)%weight=orbitals_r(i)%weight
                J_orbitals_r(i,dir,2)%filter=orbitals_r(i)%filter
            endif
        enddo
        enddo  

        !Move orbitals into Matrix Form for Overlaps

        nYl=nX*nband_block_l
        if(maxval(orbitals_l(:)%n_spinor).ne.maxval(orbitals_r(:)%n_spinor)) then
            print *, 'Spinor and nonspinor cannot overlap in Kubo Greenwood'
        endif

        allocate(Y(maxval(orbitals_l(:)%n_spinor)*nG/nX,nYl))
        allocate(AY(size(Y,1),nYr))

        !<psi_l|
        call orbitals_FFT2Matrix(orbitals_l(:), Y(:,:), grids, parallel, nX, &
                                grids(orb_grid)%reduction_option)

        !P|psi_r>
        do i=1, size(orbitals_r)
            call orbital_current(orbitals_r(i), J_orbitals_r(i,:,1), parallel, potentials, atoms, elements, all_PAW, &
                                 grids)
        enddo 
        !DPM=<psi_l|P|psi_r>               
        do dir=1,3
            call orbitals_FFT2Matrix(J_orbitals_r(:,dir,1), AY(:,:), grids, parallel, nX, &
            grids(orb_grid)%reduction_option)
            call zgemm('C','N', nYl, nYr, size(Y,1), (1.0_dp,0.0_dp), Y(:,:),  &
                size(Y,1), AY(:,:), size(AY,1), (0.0_dp,0.0_dp),  DPM(:,:,dir), nYl)
        enddo
        call parallel_task('sum', DPM(:,:,:), parallel, 'k')

        if(present(HPM)) then
            do i=1, size(orbitals_r)
                !1/2 (P H + H P) - enthalpy P |psi_r>
                call orbital_energy_current(orbitals_r(i), J_orbitals_r(i,:,1), J_orbitals_r(i,:,2), &
                         enthalpy, parallel, potentials, atoms, elements, all_PAW, grids)
            enddo
            !HPM = <psi_l| 1/2 (P H + H P) - enthalpy P  |psi_r>
            do dir=1,3
                call orbitals_FFT2Matrix(J_orbitals_r(:,dir,2), AY(:,:), grids, parallel, nX, &
                grids(orb_grid)%reduction_option)
                call zgemm('C','N', nYl, nYr, size(Y,1), (1.0_dp,0.0_dp), Y(:,:),  &
                    size(Y,1), AY(:,:), size(AY,1), (0.0_dp,0.0_dp),  HPM(:,:,dir), nYl)
            enddo
            call parallel_task('sum', HPM(:,:,:), parallel, 'k')
        endif

        do dir=1,3
            do i = 1, size(orbitals_r)
                call deallocate_orbital(J_orbitals_r(i,dir,1),grids)
                if(present(HPM)) call deallocate_orbital(J_orbitals_r(i,dir,2),grids)
            enddo
        enddo 
        
        deallocate(J_orbitals_r)
        deallocate(Y)
        deallocate(AY)
        if(parallel%myid.eq.0) print *, 'Done'
    end subroutine


    subroutine J_to_orbitals(orbitals, potentials, atoms, elements, parallel, all_PAW, grids, J_orbitals, he)

        use KG_type, only : KG_struct
        use system_type, only : system_struct
        use parallel_type, only : parallel_struct
        use odp_type, only : field_struct, orbital_struct
        use Local_Ion, only : calculate_local_ion_force
        use odp_type, only : orbital_struct
        use grids_type, only : grid_struct
        use Stochastic_Mod, only :  stochastic_vector_restore, stochastic_OMFD_filter
        use stochastic_type, only : stochastic_struct
        use odp, only: allocate_orbital
        use simulation_type, only : potentials_struct, all_PAW_struct
        use element_type, only : element_struct, PAW, HGH
        use atom_type, only : atom_struct
        use Eigen_LOBPCG, only : Orbitals_Project_out
        use constants, only : pi, i_
        use Non_Local_ion, only : Apply_projectors_RL, Calculate_Projector_overlaps, Apply_S_inverse
        use fft, only: recip_to_real, real_to_recip
        use Apply_Hamiltonian, only: Apply_H, Apply_SinvH
        use odp, only: allocate_field, deallocate_field, allocate_orbital, deallocate_orbital

        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(orbital_struct), intent(inout) :: J_orbitals(:,:,:)
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(field_struct), allocatable ::  tmp_field_1(:,:), tmp_field_2(:)
        complex(dp),  allocatable :: Pypsi(:,:,:),dPypsi(:,:,:,:)
        type(all_PAW_struct), intent(inout) :: all_PAW
        real(dp) , intent(in) ::  he

        integer :: s, k, dir

        grid=>grids(orbitals(1)%of(1)%grid)
        
        allocate(tmp_field_1( 3, maxval(orbitals(:)%n_spinor) ))
        do s=1,size(tmp_field_1,2); do dir=1,3
            call allocate_field(tmp_field_1(dir,s), grid, parallel)
            tmp_field_1(dir,s)%grid=orbitals(1)%of(1)%grid
        enddo; enddo
        
        allocate(tmp_field_2( maxval(orbitals(:)%n_spinor) ))
        do s=1,size(tmp_field_2)
            call allocate_field(tmp_field_2(s), grid, parallel)
            tmp_field_2(s)%grid=orbitals(1)%of(1)%grid
        enddo

        !Orbital = | \chi >
        do k=1, size(orbitals)
            allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(orbitals(k)%of)))
            PYpsi=0.0_dp
            !P|chi>
            do s=1,size(orbitals(k)%of)
                do dir=1,3
                    J_orbitals(k,dir,1)%of(s)%G= 0.0_dp
                enddo 
            enddo
            !HGH
            if(any(elements%PP_type.eq.HGH)) then
                PYpsi=0.0_dp
                allocate(dPYpsi(3,maxval(elements(:)%n_proj), size(atoms), size(orbitals(k)%of)))
                dPYpsi=0.0_dp
                do s=1, size(orbitals(k)%of)
                    call Calculate_Projector_overlaps(orbitals(k)%of(s),  &
                    PYpsi(:,:,s), atoms, elements, grid, parallel, dPYpsi(:,:,:,s)) 
                enddo
                call Apply_projectors_RL(orbitals(k)%of, atoms, elements, &
                grid, PYpsi, dPYpsi, midrVpsi=tmp_field_1, skip_PAW=.true.)
                do s=1, size(orbitals(k)%of)
                    do dir=1,3
                        call real_to_recip(tmp_field_1(dir,s), grids)
                        J_orbitals(k,dir,1)%of(s)%G= J_orbitals(k,dir,1)%of(s)%G + &
                            i_*tmp_field_1(dir,s)%G*grid%cutwf
                    enddo 
                enddo
            endif

             !PAW
            if(any(elements%PP_type.eq.PAW)) then
                do s=1, size(orbitals(k)%of)
                    call Calculate_Projector_overlaps(orbitals(k)%of(s),  &
                    PYpsi(:,:,s),  atoms, elements, grid, parallel, skip_HGH=.true.) 
                enddo
                call Apply_projectors_RL(orbitals(k)%of, atoms, elements, &
                    grid, PYpsi, nabla_psi=tmp_field_1, skip_HGH=.true.)
                    do s=1, size(orbitals(k)%of)
                        do dir=1,3
                        call real_to_recip(tmp_field_1(dir,s), grids)
                        J_orbitals(k,dir,1)%of(s)%G= J_orbitals(k,dir,1)%of(s)%G - &
                            i_*tmp_field_1(dir,s)%G*grid%cutwf
                    enddo 
                enddo

                PYpsi=0.0_dp
                do s=1, size(orbitals(k)%of)
                    call Calculate_Projector_overlaps(orbitals(k)%of(s),  &
                    PYpsi(:,:,s), atoms, elements, grid, parallel, skip_HGH=.true.) 
                enddo
                call Apply_projectors_RL(orbitals(k)%of, atoms, elements, &
                    grid, PYpsi, nabla_psi=tmp_field_1, skip_HGH=.true.) 
                do s=1,size(orbitals(k)%of)
                    do dir=1,3
                        call real_to_recip(tmp_field_1(dir,s), grids)
                        J_orbitals(k,dir,1)%of(s)%G= J_orbitals(k,dir,1)%of(s)%G - &
                        i_*tmp_field_1(dir,s)%G*grid%cutwf
                    enddo 
                enddo
                call Apply_projectors_RL(orbitals(k)%of, atoms, elements, &
                grid, PYpsi, midrVpsi=tmp_field_1, skip_HGH=.true.) 
                do s=1,size(orbitals(k)%of)
                    do dir=1,3
                        call real_to_recip(tmp_field_1(dir,s), grids) !fix to do 1 FFT, use 2 fields
                        J_orbitals(k,dir,1)%of(s)%G= J_orbitals(k,dir,1)%of(s)%G + &
                        i_*tmp_field_1(dir,s)%G*grid%cutwf
                    enddo 
                enddo
                call Apply_SinvH(orbitals(k)%of, tmp_field_2, grids, potentials, atoms, elements, &
                                 parallel, all_PAW, calc_G=.false., PYpsi_in=PYpsi, with_CG=.false.)
                PYpsi=0.0_dp
                do s=1, size(orbitals(k)%of)
                    call Calculate_Projector_overlaps(tmp_field_2(s),  &
                    PYpsi(:,:,s), atoms, elements, grid, parallel, skip_HGH=.true.) 
                enddo
                call Apply_projectors_RL(tmp_field_2, atoms, elements, &
                                 grid, PYpsi, midrSpsi=tmp_field_1, skip_HGH=.true.) 
                do s=1,size(orbitals(k)%of)
                    do dir=1,3
                        call real_to_recip(tmp_field_1(dir,s), grids) !fix to do 1 FFT, use 2 fields
                        J_orbitals(k,dir,1)%of(s)%G= J_orbitals(k,dir,1)%of(s)%G - &
                        i_*tmp_field_1(dir,s)%G*grid%cutwf
                    enddo 
                enddo
                call Apply_projectors_RL(tmp_field_2, atoms, elements, &
                                 grid, PYpsi, dipole_psi=tmp_field_1, skip_HGH=.true.) 
                do s=1,size(orbitals(k)%of)
                    do dir=1,3
                        call real_to_recip(tmp_field_1(dir,s), grids) !fix to do 1 FFT, use 2 fields
                        J_orbitals(k,dir,1)%of(s)%G= J_orbitals(k,dir,1)%of(s)%G - &
                        i_*tmp_field_1(dir,s)%G*grid%cutwf
                    enddo 
                enddo

            endif

            !Plane Wave
            do s=1, size(orbitals(k)%of)
                do dir=1,3
                    J_orbitals(k,dir,1)%of(s)%G= J_orbitals(k,dir,1)%of(s)%G + &
                        grid%G(:,:,:,dir)*orbitals(k)%of(s)%G*grid%cutwf
                        call recip_to_real(J_orbitals(k,dir,1)%of(s), grids)
                enddo 
            enddo              

            !PH|sqrt(fd) chi>
            if(any(elements%PP_type.eq.PAW)) then
                call Apply_SinvH(orbitals(k)%of(:), tmp_field_2(:), grids, potentials, &
                      atoms, elements, parallel, all_PAW, calc_G=.false., PYpsi_in=PYpsi)
            else
                call Apply_H(orbitals(k)%of(:), tmp_field_2(:), grids, potentials, &
                      atoms, elements, parallel,  calc_R=.false., PYpsi_in=PYpsi)
            endif
            
            do s=1,size(tmp_field_2)
                call recip_to_real(tmp_field_2(s), grids)
            enddo
            
            do s=1,size(orbitals(k)%of)
                do dir=1,3
                    J_orbitals(k,dir,2)%of(s)%G= 0.0_dp
                enddo 
            enddo
            PYpsi=0.0_dp
            if(any(elements%PP_type.eq.HGH)) then
                dPYpsi=0.0_dp
                do s=1, size(tmp_field_2)
                    call Calculate_Projector_overlaps(tmp_field_2(s),  &
                    PYpsi(:,:,s), atoms, elements, grid, parallel, dPYpsi(:,:,:,s)) 
                enddo

                call Apply_projectors_RL(tmp_field_2, atoms, elements, &
                grid, PYpsi, dPYpsi, midrVpsi=tmp_field_1, skip_PAW=.true.)
                do s=1, size(tmp_field_2)
                    do dir=1,3
                        call real_to_recip(tmp_field_1(dir,s), grids)
                        J_orbitals(k,dir,2)%of(s)%G= J_orbitals(k,dir,2)%of(s)%G + &
                            i_*tmp_field_1(dir,s)%G*grid%cutwf
                    enddo 
                enddo
            endif

            !PAW
            if(any(elements%PP_type.eq.PAW)) then
                do s=1, size(tmp_field_2)
                    call Calculate_Projector_overlaps(tmp_field_2(s),  &
                    PYpsi(:,:,s),  atoms, elements, grid, parallel, skip_HGH=.true.) 
                enddo
                call Apply_projectors_RL(tmp_field_2, atoms, elements, &
                    grid, PYpsi, nabla_psi=tmp_field_1, skip_HGH=.true.)
                    do s=1, size(tmp_field_2)
                        do dir=1,3
                        call real_to_recip(tmp_field_1(dir,s), grids)
                        J_orbitals(k,dir,2)%of(s)%G= J_orbitals(k,dir,2)%of(s)%G - &
                            i_*tmp_field_1(dir,s)%G*grid%cutwf
                    enddo 
                enddo
            endif

            !Plane Wave
            do s=1, size(tmp_field_2)
                do dir=1,3
                    J_orbitals(k,dir,2)%of(s)%G= J_orbitals(k,dir,2)%of(s)%G + &
                        grid%G(:,:,:,dir)*tmp_field_2(s)%G*grid%cutwf
                    call recip_to_real(J_orbitals(k,dir,2)%of(s), grids)
                enddo 
            enddo   

            !HP|sqrt(fd) chi>
            do dir=1,3
                if(any(elements%PP_type.eq.PAW)) then
                    call Apply_S_inverse(J_orbitals(k,dir,1)%of, atoms, elements, grids, &
                    parallel, tmp_field_1(1,:), all_PAW, with_CG=.false.)
                    call Apply_H(tmp_field_1(1,:), tmp_field_2(:), grids, potentials, &
                          atoms, elements, parallel,  calc_R=.false.)
                else
                    call Apply_H(J_orbitals(k,dir,1)%of, tmp_field_2(:), grids, potentials, &
                          atoms, elements, parallel,  calc_R=.false.)
                endif
                do s=1,size(orbitals(k)%of)
                    J_orbitals(k,dir,2)%of(s)%G = J_orbitals(k,dir,2)%of(s)%G + &
                                                    tmp_field_2(s)%G
                    !(PH/2+HP/2-heP)|sqrt(fd) chi>
                    J_orbitals(k,dir,2)%of(s)%G= 0.5_dp*J_orbitals(k,dir,2)%of(s)%G &
                        - he*J_orbitals(k,dir,1)%of(s)%G
                    J_orbitals(k,dir,2)%of(s)%G = J_orbitals(k,dir,2)%of(s)%G*grid%cutwf
                    call recip_to_real(J_orbitals(k,dir,2)%of(s), grids)
                enddo 
            enddo

            deallocate(PYpsi)
            if(any(elements%PP_type.eq.HGH)) deallocate(dPYpsi)
        enddo

        do s=1,size(tmp_field_1,2); do dir=1,3
            call deallocate_field(tmp_field_1(dir,s),grids(tmp_field_1(dir,s)%grid))
        enddo; enddo
        deallocate(tmp_field_1)

        do s=1,size(tmp_field_2)
            call deallocate_field(tmp_field_2(s),grids(tmp_field_2(s)%grid))
        enddo
        deallocate(tmp_field_2)

    end subroutine

    subroutine J_Matrix_Elements(orbitals_l, orbitals_r, &
        grids, parallel, potentials, atoms, elements, all_PAW, DPM, HPM, he, symmetric_)
    use constants, only : pi, i_
    use grids_type, only : grid_struct
    use odp_type, only : field_struct, orbital_struct, deterministic, stochastic
    use simulation_type, only : potentials_struct
    use element_type, only : element_struct, HGH, PAW
    use atom_type, only : atom_struct
    use fft, only: recip_to_real, real_to_recip
    use lapack, only: zgemm

    use parallel_mod,only : parallel_task, parallel_wait
    use parallel_type, only : parallel_struct
    use numerics_type, only : numerics_struct
    use Apply_Hamiltonian, only: Apply_H, Apply_SinvH, Apply_HSinv
    use FFT2Matrix, only : orbitals_FFT2Matrix, orbitals_Matrix2FFT, grids_G2_FFT2Matrix, &
                           orbital_FFT2Matrix
    use odp, only: allocate_field, deallocate_field, allocate_orbital, deallocate_orbital 
    use Non_Local_ion, only : Apply_projectors_RL, Calculate_Projector_overlaps 
    use simulation_type, only : all_PAW_struct


    type(parallel_struct), intent(in) :: parallel
    type(orbital_struct), intent(inout) :: orbitals_l(:), orbitals_r(:)
    type(potentials_struct), intent(in) :: potentials
    type(atom_struct), intent(inout) :: atoms(:)
    type(element_struct), intent(inout) :: elements(:)
    type(all_PAW_struct), intent(inout) :: all_PAW
    type(grid_struct), intent(inout), target :: grids(:)
    type(grid_struct), pointer:: grid
    complex(dp),intent(inout) :: DPM(:,:,:)
    complex(dp),intent(inout), optional ::  HPM(:,:,:)
    real(dp), intent(in), optional :: he
    logical, intent(in), optional :: symmetric_
    logical :: symmetric
    complex(dp), allocatable :: Y(:,:), AY(:,:)
    integer :: nYl, nYr, nX, nband_block_l, nband_block_r, nG
    integer:: i
    integer :: orb_grid
    type(field_struct), allocatable ::  tmp_field_1(:,:), tmp_field_2(:), tmp_field_3(:,:), tmp_field_4(:)
    type(orbital_struct), allocatable:: J_orbitals_r(:,:)
    integer :: dir, s
    complex(dp),  allocatable :: Pypsi(:,:,:),dPypsi(:,:,:,:)

    if(present(symmetric_)) then
        symmetric= symmetric_
    else 
        symmetric=.false.
    endif

    if(parallel%myid.eq.0) print *, 'Starting Kubo Greenwood'

    grid=>grids(orbitals_l(1)%of(1)%grid)
    allocate(tmp_field_1( 3, maxval(orbitals_r(:)%n_spinor) ))
    do s=1,size(tmp_field_1,2); do dir=1,3
        call allocate_field(tmp_field_1(dir,s), grid, parallel)
        tmp_field_1(dir,s)%grid=orbitals_r(1)%of(1)%grid
    enddo; enddo
    allocate(tmp_field_3( 3, maxval(orbitals_r(:)%n_spinor) ))
    do s=1,size(tmp_field_3,2); do dir=1,3
        call allocate_field(tmp_field_3(dir,s), grid, parallel)
        tmp_field_3(dir,s)%grid=orbitals_r(1)%of(1)%grid
    enddo; enddo
    allocate(tmp_field_2( maxval(orbitals_r(:)%n_spinor) ))
    do s=1,size(tmp_field_2)
        call allocate_field(tmp_field_2(s), grid, parallel)
        tmp_field_2(s)%grid=orbitals_r(1)%of(1)%grid
    enddo
    if(present(HPM)) then
        allocate(tmp_field_4( maxval(orbitals_r(:)%n_spinor) ))
        do s=1,size(tmp_field_4)
            call allocate_field(tmp_field_4(s), grid, parallel)
            tmp_field_4(s)%grid=orbitals_r(1)%of(1)%grid
        enddo
    endif

    nband_block_l=size(orbitals_l)
    nX=parallel%n_band_groups
    nYl=nX*nband_block_l
    orb_grid=orbitals_l(1)%of(1)%grid
    if(orb_grid.ne.orbitals_r(1)%of(1)%grid) then
        print *, 'Trying to calculate Matrix element between different grids'
    endif

    nband_block_r=size(orbitals_r)
    nX=parallel%n_band_groups
    nYr=nX*nband_block_r

    if(symmetric.and.(nYr.ne.nYl)) then
        print *, 'You say DPM is symmetric, but dont pass same orbitals to Bra and Ket'
    endif

    do i=1,nband_block_l
        if(orbitals_l(i)%of(1)%grid.ne.orb_grid) then
            print *, 'All orbitals passed into Egensolver need to be on the same grid (same k-point)'
        endif
    enddo
    
    do i=1,nband_block_r
        if(orbitals_l(i)%of(1)%grid.ne.orb_grid) then
            print *, 'All orbitals passed into Egensolver need to be on the same grid (same k-point)'
        endif
    enddo


    grid=>grids(orb_grid)
    if(grid%reduction_option.eq.2) then
        nG=product(grid%Ng_small)
    else
        nG=product(grid%Ng_local)
    endif
    if(mod(nG,nX).ne.0) nG=(nG/nX+1)*nX


    allocate(J_orbitals_r(size(orbitals_r),3))

    do dir=1,3
    do i = 1, size(orbitals_r)
        J_orbitals_r(i,dir)%band=orbitals_r(i)%band
        J_orbitals_r(i,dir)%k_point=orbitals_r(i)%k_point
        J_orbitals_r(i,dir)%spin=orbitals_r(i)%spin
        J_orbitals_r(i,dir)%degeneracy=orbitals_r(i)%degeneracy
        J_orbitals_r(i,dir)%type=orbitals_r(i)%type
        J_orbitals_r(i,dir)%n_spinor=orbitals_r(i)%n_spinor
        call allocate_orbital(J_orbitals_r(i,dir), grids, orbitals_r(i)%k_point+2, parallel)
        J_orbitals_r(i,dir)%occ=orbitals_r(i)%occ
        J_orbitals_r(i,dir)%weight=orbitals_r(i)%weight
        J_orbitals_r(i,dir)%filter=orbitals_r(i)%filter
    enddo
    enddo  

    !Move orbitals into Matrix Form for Overlaps

    nYl=nX*nband_block_l
    if(maxval(orbitals_l(:)%n_spinor).ne.maxval(orbitals_r(:)%n_spinor)) then
        print *, 'Spinor and nonspinor cannot overlap in Kubo Greenwood'
    endif

    allocate(Y(maxval(orbitals_l(:)%n_spinor)*nG/nX,nYl))
    allocate(AY(size(Y,1),nYr))

    !<psi_l|
    call orbitals_FFT2Matrix(orbitals_l(:), Y(:,:), grids, parallel, nX, &
                            grids(orb_grid)%reduction_option)

    !P|psi_r>
    do i=1, size(orbitals_r)
        allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(orbitals_r(i)%of)))

        do s=1,size(orbitals_r(i)%of)
            do dir=1,3
                J_orbitals_r(i,dir)%of(s)%G= 0.0_dp
                J_orbitals_r(i,dir)%of(s)%R= 0.0_dp
            enddo 
        enddo

        !HGH
        if(any(elements%PP_type.eq.HGH)) then
            PYpsi=0.0_dp
            allocate(dPYpsi(3,maxval(elements(:)%n_proj), size(atoms), size(orbitals_r(i)%of)))
            dPYpsi=0.0_dp
            do s=1, size(orbitals_r(i)%of)
                call Calculate_Projector_overlaps(orbitals_r(i)%of(s),  &
                PYpsi(:,:,s), atoms, elements, grid, parallel, dPYpsi(:,:,:,s), skip_PAW=.true.) 
            enddo

            call Apply_projectors_RL(orbitals_r(i)%of, atoms, elements, &
                grid, PYpsi, dPYpsi, midrVpsi=tmp_field_1, skip_PAW=.true.)
            do s=1,size(orbitals_r(i)%of)
                do dir=1,3
                    J_orbitals_r(i,dir)%of(s)%R=tmp_field_1(dir,s)%R
                enddo 
            enddo
            deallocate(dPYpsi)
        endif
         !PAW
        if(any(elements%PP_type.eq.PAW)) then
            PYpsi=0.0_dp
            allocate(dPYpsi(3,maxval(elements(:)%n_proj), size(atoms), size(orbitals_r(i)%of)))
            do s=1, size(orbitals_r(i)%of)
                call Calculate_Projector_overlaps(orbitals_r(i)%of(s),  &
                PYpsi(:,:,s), atoms, elements, grid, parallel, dPYpsi(:,:,:,s), skip_HGH=.true.) 
            enddo

          !Basic PAW P
                 !Phat= (G+k) -i_*|Pi>(nabla_ij +i_*k*sij)<Pj|

                !call Apply_projectors_RL(orbitals_r(i)%of, atoms, elements, &
                !    grid, PYpsi, nabla_psi=current_field, skip_HGH=.true.) 

                !Better PAW P=i[H,r] (less sensitive to PAW projector set)
                !Phat= (G+k) i_*(-[r,V_NL] +( r|Pi>sij<Pj| - |Pi>rij<Pj| ) * S^-1 H +
                !                       + H S^-1 * ( |Pi>rij<Pj| - |Pi>sij<Pj|r  ) )

                !-[r,V_NL]|Psi>
            call Apply_projectors_RL(orbitals_r(i)%of, atoms, elements, &
            grid, PYpsi, dPYpsi, midrVpsi=tmp_field_1, skip_HGH=.true.) 
            do s=1,size(orbitals_r(i)%of)
                do dir=1,3
                    J_orbitals_r(i,dir)%of(s)%R=J_orbitals_r(i,dir)%of(s)%R+tmp_field_1(dir,s)%R
                enddo 
            enddo

            ! |Pi>sij<Pj|r|Psi>
            call Apply_projectors_RL(orbitals_r(i)%of, atoms, elements, &
            grid, PYpsi, dPYpsi, SRpsi=tmp_field_3, skip_HGH=.true.) 

            ! |Pi>rij<Pj|Psi>
            call Apply_projectors_RL(orbitals_r(i)%of, atoms, elements, &
            grid, PYpsi, dipole_psi=tmp_field_1, skip_HGH=.true.) 

            do s=1,size(orbitals_r(i)%of)
                do dir=1,3
                    tmp_field_1(dir,s)%R=tmp_field_1(dir,s)%R-tmp_field_3(dir,s)%R
                enddo 
            enddo

            !H S^-1 * ( |Pi>rij<Pj| - |Pi>sij<Pj|r  ) |Psi>
            do dir=1,3
                call Apply_HSinv(tmp_field_1(dir,:), tmp_field_3(dir,:), grids, potentials, atoms, elements, &
                parallel, all_PAW, calc_R=.false., with_CG=.false.)
            enddo

            do s=1,size(orbitals_r(i)%of)
                do dir=1,3
                    J_orbitals_r(i,dir)%of(s)%R=J_orbitals_r(i,dir)%of(s)%R+tmp_field_3(dir,s)%R
                enddo
            enddo

            !S^-1 H |Psi>
            call Apply_SinvH(orbitals_r(i)%of, tmp_field_2, grids, potentials, atoms, elements, &
            parallel, all_PAW, calc_G=.false., PYpsi_in=PYpsi, with_CG=.false.)

            PYpsi=0.0_dp
            do s=1, size(orbitals_r(i)%of)
                call Calculate_Projector_overlaps(tmp_field_2(s),  &
                PYpsi(:,:,s), atoms, elements, grid, parallel, skip_HGH=.true.) 
            enddo
            !r|Pi>sij<Pj|S^-1 H |Psi>
            call Apply_projectors_RL(tmp_field_2, atoms, elements, &
                             grid, PYpsi, RSpsi=tmp_field_3, skip_HGH=.true.) 
            ! |Pi>rij<Pj|S^-1 H |Psi>
            call Apply_projectors_RL(tmp_field_2, atoms, elements, &
                             grid, PYpsi, dipole_psi=tmp_field_1, skip_HGH=.true.) 

            do s=1,size(orbitals_r(i)%of)
                do dir=1,3
                    J_orbitals_r(i,dir)%of(s)%R=J_orbitals_r(i,dir)%of(s)%R + &
                                            tmp_field_3(dir,s)%R - tmp_field_1(dir,s)%R
                enddo
            enddo

            deallocate(dPYpsi)
        endif

        !Plane Wave
        do s=1,size(orbitals_r(i)%of)
            do dir=1,3
                if(any(elements%PP_type.eq.PAW).or.any(elements%PP_type.eq.HGH)) then
                    J_orbitals_r(i,dir)%of(s)%R=i_*J_orbitals_r(i,dir)%of(s)%R
                    call real_to_recip(J_orbitals_r(i,dir)%of(s), grids)
                else
                    J_orbitals_r(i,dir)%of(s)%G=0.0_dp
                endif
                
                J_orbitals_r(i,dir)%of(s)%G= J_orbitals_r(i,dir)%of(s)%G + &
                     grid%G(:,:,:,dir)*orbitals_r(i)%of(s)%G
                J_orbitals_r(i,dir)%of(s)%G=J_orbitals_r(i,dir)%of(s)%G*grid%cutwf
            enddo 
        enddo
        deallocate(PYpsi)
    enddo 
    !DPM=<psi_l|P|psi_r>               
    do dir=1,3
        call orbitals_FFT2Matrix(J_orbitals_r(:,dir), AY(:,:), grids, parallel, nX, &
        grids(orb_grid)%reduction_option)
        call zgemm('C','N', nYl, nYr, size(Y,1), (1.0_dp,0.0_dp), Y(:,:),  &
            size(Y,1), AY(:,:), size(AY,1), (0.0_dp,0.0_dp),  DPM(:,:,dir), nYl)
    enddo

    if(present(HPM)) then
        !<psi_l| H
        do i=1, size(orbitals_l)
            if(orbitals_l(i)%type.eq.deterministic) then
                do s=1,size(orbitals_l(i)%of)
                    orbitals_l(i)%of(s)%G=orbitals_l(i)%of(s)%G*orbitals_l(i)%eig(1)
                enddo
            else
                if(any(elements%PP_type.eq.PAW)) then
                    call Apply_SinvH(orbitals_l(i)%of(:), tmp_field_4(:), grids, potentials, &
                        atoms, elements, parallel, all_PAW, calc_G=.false.)
                else
                    call Apply_H(orbitals_l(i)%of(:), tmp_field_4(:), grids, potentials,&
                         atoms, elements, parallel, calc_R=.false.)
                endif
                do s=1,size(orbitals_l(i)%of)
                        if(any(elements%PP_type.eq.PAW)) then
                                call real_to_recip(tmp_field_4(s), grids)
                        endif
                    orbitals_l(i)%of(s)%G= tmp_field_4(s)%G*grid%cutwf
                    !don't do inverse FFT or overwrite %R field (save for restoration)                                
                enddo
            endif
        enddo
        
        call orbitals_FFT2Matrix(orbitals_l(:), Y(:,:), grids, parallel, nX, &
                            grids(orb_grid)%reduction_option)

        !<psi_l| H P |psi_r>
        do dir=1,3
            call orbitals_FFT2Matrix(J_orbitals_r(:,dir), AY(:,:), grids, parallel, nX, &
                     grids(orb_grid)%reduction_option)
            call zgemm('C','N', nYl, nYr, size(Y,1), (1.0_dp,0.0_dp), Y(:,:),  &
                size(Y,1), AY(:,:), size(AY,1), (0.0_dp,0.0_dp),  HPM(:,:,dir), nYl)
        enddo
        
        !<psi_l|
        do i=1, size(orbitals_l)
            if(orbitals_l(i)%type.eq.deterministic) then
                do s=1,size(orbitals_l(i)%of)
                    orbitals_l(i)%of(s)%G=orbitals_l(i)%of(s)%G/orbitals_l(i)%eig(1)
                enddo
            else
                do s=1,size(orbitals_l(i)%of)
                     ! FFT of saved %R field to restore %G
                    call real_to_recip(orbitals_l(i)%of(s), grids)
                enddo
            endif
        enddo

        if(symmetric) then
            !<psi_l| P H |psi_r>
            do dir=1,3
                call orbitals_FFT2Matrix(J_orbitals_r(:,dir), AY(:,:), grids, parallel, nX, &
                     grids(orb_grid)%reduction_option)

                call zgemm('C','N', nYr, nYl, size(AY,1), (0.5_dp,0.0_dp), AY(:,:),  &
                    size(AY,1), Y(:,:), size(Y,1), (0.5_dp,0.0_dp),  HPM(:,:,dir), nYr)
            enddo
        else
            !<psi_l|
            call orbitals_FFT2Matrix(orbitals_l(:), Y(:,:), grids, parallel, nX, &
                                grids(orb_grid)%reduction_option)
            do i=1, size(orbitals_r)
                allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(orbitals_r(i)%of)))
    
                ! H |psi_r>
                if(orbitals_r(i)%type.eq.deterministic) then
                    do s=1,size(orbitals_r(i)%of)
                        orbitals_r(i)%of(s)%G=orbitals_r(i)%of(s)%G*orbitals_r(i)%eig(1)
                        orbitals_r(i)%of(s)%R=orbitals_r(i)%of(s)%R*orbitals_r(i)%eig(1)
                    enddo
                else
                    PYpsi=0.0_dp
                    do s=1, size(orbitals_r(i)%of)
                        call Calculate_Projector_overlaps(orbitals_r(i)%of(s),  &
                        PYpsi(:,:,s), atoms, elements, grid, parallel) 
                    enddo
                    if(any(elements%PP_type.eq.PAW)) then
                        call Apply_SinvH(orbitals_r(i)%of(:), tmp_field_2(:), grids, potentials, &
                            atoms, elements, parallel, all_PAW, calc_G=.false., Pypsi_in=PYpsi)
                    else
                        call Apply_H(orbitals_r(i)%of(:), tmp_field_2(:),  &
                            grids, potentials, atoms, elements, parallel, calc_R=.false., Pypsi_in=PYpsi)
                    endif
                    do s=1,size(orbitals_r(i)%of)
                            if(any(elements%PP_type.eq.PAW)) then
                                call real_to_recip(tmp_field_2(s), grids)
                            endif
                            tmp_field_2(s)%R=orbitals_r(i)%of(s)%R
                            !Store %R for restore
                            orbitals_r(i)%of(s)%G= tmp_field_2(s)%G*grid%cutwf
                            call recip_to_real(orbitals_r(i)%of(s), grids)
                    enddo
                endif 
                ! P H |psi_r>
                PYpsi=0.0_dp
                do s=1,size(orbitals_r(i)%of)
                    do dir=1,3
                        J_orbitals_r(i,dir)%of(s)%G= 0.0_dp
                    enddo 
                enddo
                if(any(elements%PP_type.eq.HGH)) then
                    allocate(dPYpsi(3,maxval(elements(:)%n_proj), size(atoms), size(orbitals_r(i)%of)))
                    dPYpsi=0.0_dp

                    do s=1, size(orbitals_r(i)%of)
                        call Calculate_Projector_overlaps(orbitals_r(i)%of(s),  &
                        PYpsi(:,:,s), atoms, elements, grid, parallel, dPYpsi(:,:,:,s), skip_PAW=.true.) 
                    enddo
    
                    call Apply_projectors_RL(orbitals_r(i)%of, atoms, elements, &
                        grid, PYpsi, dPYpsi, midrVpsi=tmp_field_1, skip_PAW=.true.)
                    do s=1,size(orbitals_r(i)%of)
                        do dir=1,3
                            call real_to_recip(tmp_field_1(dir,s), grids)
                            J_orbitals_r(i,dir)%of(s)%G= J_orbitals_r(i,dir)%of(s)%G + i_*tmp_field_1(dir,s)%G
                        enddo 
                    enddo
                    deallocate(dPYpsi)
                endif
                !PAW
                if(any(elements%PP_type.eq.PAW)) then
                    PYpsi=0.0_dp
                    do s=1, size(orbitals_r(i)%of)
                        call Calculate_Projector_overlaps(orbitals_r(i)%of(s),  &
                        PYpsi(:,:,s), atoms, elements, grid, parallel, skip_HGH=.true.) 
                    enddo
                    call Apply_projectors_RL(orbitals_r(i)%of, atoms, elements, &
                        grid, PYpsi, nabla_psi=tmp_field_1, skip_HGH=.true.)
                    do s=1,size(orbitals_r(i)%of)
                        do dir=1,3
                            call real_to_recip(tmp_field_1(dir,s), grids)
                            J_orbitals_r(i,dir)%of(s)%G= J_orbitals_r(i,dir)%of(s)%G - i_*tmp_field_1(dir,s)%G
                        enddo 
                    enddo
                endif

                !Plane Wave
                do s=1,size(orbitals_r(i)%of)
                    do dir=1,3
                        J_orbitals_r(i,dir)%of(s)%G= J_orbitals_r(i,dir)%of(s)%G + &
                            grid%G(:,:,:,dir)*orbitals_r(i)%of(s)%G
                        J_orbitals_r(i,dir)%of(s)%G=J_orbitals_r(i,dir)%of(s)%G*grid%cutwf
                    enddo 
                enddo

                do s=1,size(orbitals_r(i)%of)
                    if(orbitals_r(i)%type.eq.deterministic) then
                        orbitals_r(i)%of(s)%R = orbitals_r(i)%of(s)%R/orbitals_r(i)%eig(1)
                    else
                        orbitals_r(i)%of(s)%R = tmp_field_2(s)%R
                    endif
                    call real_to_recip(orbitals_r(i)%of(s), grids) 
                enddo
                                   
                deallocate(PYpsi)
            enddo 

            ! <psi_l| P H |psi_r>
            !HPM = <psi_l| 1/2 (P H + H P) |psi_r> - he*<psi_l| P |psi_r>
            do dir=1,3
                call orbitals_FFT2Matrix(J_orbitals_r(:,dir), AY(:,:), grids, parallel, nX, &
                grids(orb_grid)%reduction_option)
                call zgemm('C','N', nYl, nYr, size(Y,1), (0.5_dp,0.0_dp), Y(:,:),  &
                    size(Y,1), AY(:,:), size(AY,1), (0.5_dp,0.0_dp),  HPM(:,:,dir), nYl)
            enddo

        endif
        do dir=1,3
            HPM(:,:,dir)=HPM(:,:,dir)-he*DPM(:,:,dir)
        enddo
        call parallel_task('sum', HPM(:,:,:), parallel, 'k')
    endif

    call parallel_task('sum', DPM(:,:,:), parallel, 'k')

    do dir=1,3
        do i = 1, size(orbitals_r)
            call deallocate_orbital(J_orbitals_r(i,dir),grids)
        enddo
    enddo 

    deallocate(J_orbitals_r)
    do s=1,size(tmp_field_1,2); do dir=1,3
        call deallocate_field(tmp_field_1(dir,s),grids(tmp_field_1(dir,s)%grid))
        call deallocate_field(tmp_field_3(dir,s),grids(tmp_field_3(dir,s)%grid))
    enddo; enddo
    do s=1,size(tmp_field_2)
        call deallocate_field(tmp_field_2(s),grids(tmp_field_2(s)%grid))
    enddo
    deallocate(tmp_field_1)
    deallocate(tmp_field_2)
    deallocate(tmp_field_3)
    deallocate(Y)
    deallocate(AY)
    
end subroutine


end module