module Density
    use types, only : op, dp
    use constants, only : pi
    implicit none
    public :: calc_density, calc_Compensation_Charge_Density, calc_core_Charge_Density, calc_valance_Charge_Density_for_guess, &
    Force_or_Stress_Comp_charge, calculate_core_xc_force, core_xc_Stress, calc_coarse_density

    contains

    subroutine calc_density(coarse_density, fine_density, psi, grids, parallel, coarse_to_fine, all_PAW)
        use odp_type, only: spin_DenPot_struct,orbital_struct
        use grids_type, only: grid_struct, inter_grid_struct
        use parallel_type, only: parallel_struct
        use parallel_mod, only: parallel_task
        use fft, only: real_to_recip, recip_to_real
        use grids_mod, only : Gspace_grid_to_grid_transfer
        use operations_3D, only: integrate_3D_R
        use simulation_type, only : all_PAW_struct


        type(orbital_struct), intent(in) :: psi(:,:,:)
        type(spin_DenPot_struct), intent(inout) :: coarse_density, fine_density
        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(inout) :: grids(:)
        type(inter_grid_struct), intent(inout) :: coarse_to_fine
        type(all_PAW_struct), intent(inout), target :: all_PAW

        integer :: i

        call calc_coarse_density(coarse_density, psi, grids, parallel)

        do i=1,coarse_density%n_s
            call Gspace_grid_to_grid_transfer(grids(2), grids(1), parallel, coarse_density%of(i)%G, &
                  fine_density%of(i)%G, coarse_to_fine)
            call recip_to_real(fine_density%of(i), grids)
            fine_density%of(i)%R=abs(fine_density%of(i)%R) ! G cutting may cause negative if near zero
            ! add the compensation charge density if needed
            if(all_PAW%N_PAW_atoms.gt.0) then
                fine_density%of(i)%G=fine_density%of(i)%G + all_PAW%rho_comp%of(i)%G
                fine_density%of(i)%R=fine_density%of(i)%R + all_PAW%rho_comp%of(i)%R
            endif
            if(parallel%myid.eq.0) print *, 'Nelectrons of spin ', i, ' by density integeration', &
                fine_density%of(i)%G(1,1,1)*product(grids(1)%box_length)
        enddo

    end subroutine

    subroutine calc_coarse_density(coarse_density, psi, grids, parallel)
        use odp_type, only: spin_DenPot_struct,orbital_struct
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use parallel_mod, only: parallel_task
        use fft, only: real_to_recip, recip_to_real


        type(orbital_struct), intent(in) :: psi(:,:,:)
        type(spin_DenPot_struct), intent(inout) :: coarse_density
        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(inout) :: grids(:)

        integer :: i,j,k,s, spin

        do i=1,coarse_density%n_s
            coarse_density%of(i)%R=0.0_dp
        enddo

        !Need to generalize this for spinor which gives density matrix with off-diagonal components

        !Note that i and s both relate account for spin,
        !but only one of the two should be able to be > 1 at a time
        do i=1,size(psi,3);do j=1,size(psi,2);do k=1,size(psi,1); do s=1, psi(k,j,i)%n_spinor  
            if(abs(psi(k,j,i)%occ(s)).lt.tiny(1.0_dp)) cycle   
            spin=s*i
            coarse_density%of(spin)%R= coarse_density%of(s)%R + &
                    psi(k,j,i)%weight(s)*psi(k,j,i)%occ(s)*&
                        abs(psi(k,j,i)%of(s)%R)**2
        enddo; enddo; enddo; enddo

        do i=1,coarse_density%n_s
            call parallel_task('sum', coarse_density%of(i)%R(:,:,:), parallel, 'space') 
            call real_to_recip(coarse_density%of(i), grids)
            coarse_density%of(i)%G=coarse_density%of(i)%G*grids(2)%cutden
            call recip_to_real(coarse_density%of(i), grids)
            coarse_density%of(i)%R=abs(coarse_density%of(i)%R) ! G cutting may cause negative if near zero
        enddo

    end subroutine

    subroutine G_space_Fine_density_interp(coarse_density, fine_density, grids, parallel, coarse_to_fine, all_PAW)
        use odp_type, only: spin_DenPot_struct
        use grids_type, only: grid_struct, inter_grid_struct
        use parallel_type, only: parallel_struct
        use parallel_mod, only: parallel_task
        use fft, only: real_to_recip, recip_to_real
        use grids_mod, only : Gspace_grid_to_grid_transfer
        use operations_3D, only: integrate_3D_R
        use simulation_type, only : all_PAW_struct

        type(spin_DenPot_struct), intent(in) :: coarse_density
        type(spin_DenPot_struct), intent(inout) :: fine_density
        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(inout) :: grids(:)
        type(inter_grid_struct), intent(inout) :: coarse_to_fine
        type(all_PAW_struct), intent(inout), target :: all_PAW

        integer :: i

        do i=1,coarse_density%n_s
            call Gspace_grid_to_grid_transfer(grids(2), grids(1), parallel, coarse_density%of(i)%G, &
                  fine_density%of(i)%G, coarse_to_fine)
            call recip_to_real(fine_density%of(i), grids)
            fine_density%of(i)%R=abs(fine_density%of(i)%R) ! G cutting may cause negative if near zero
            ! add the compensation charge density if needed
            if(all_PAW%N_PAW_atoms.gt.0) then
                fine_density%of(i)%G=fine_density%of(i)%G + all_PAW%rho_comp%of(i)%G
                fine_density%of(i)%R=fine_density%of(i)%R + all_PAW%rho_comp%of(i)%R
            endif
            if(parallel%myid.eq.0) print *, 'Nelectrons of spin ', i, ' by density integeration', &
                fine_density%of(i)%G(1,1,1)*product(grids(1)%box_length)
        enddo
    end subroutine

    subroutine calc_Compensation_Charge_Density(all_PAW, grids, atoms, parallel)
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use parallel_mod, only: parallel_task
        use atom_type, only : atom_struct
        use simulation_type, only : all_PAW_struct
        use m_paral_atom, only : get_my_natom, get_my_atmtab
        use m_paw_nhat, only : pawmknhat
        use fft, only: real_to_recip, recip_to_real

        type(parallel_struct), intent(in) :: parallel
        type(atom_struct), intent(inout) :: atoms(:)
        type(grid_struct), intent(inout) :: grids(:)

        type(all_PAW_struct), intent(inout) :: all_PAw
        integer ::   i, j, nhatgrdim
        integer :: ix_local,iy_local,iz_local, ic
        real(dp), allocatable :: rho(:,:), gr_rho(:,:,:)
 
        real(dp) :: gprimd(3,3),qphon(3),rprimd(3,3),xred(3,all_PAW%N_PAW_atoms)


        !Divide the atoms up by same space (spin, k-point and band redundant here)

        allocate(rho(product(grids(1)%Nr_local),all_paw%rho_comp%n_s))
        rho = 0.0_dp

        nhatgrdim=0
        allocate(gr_rho(product(grids(1)%Nr_local),all_paw%gr_rho_comp(1)%n_s,3*nhatgrdim))

        all_PAW%paral_atom=.true.
        all_PAW%my_atmtab=>null()
        call get_my_natom(parallel%comm_space,all_PAW%my_natom, all_PAW%N_PAW_atoms)
        call get_my_atmtab(parallel%comm_space,all_PAW%my_atmtab,all_PAW%my_atmtab_allocated, all_PAW%paral_atom, &
                            all_PAW%N_PAW_atoms)

        gprimd=0.0_dp
        gprimd(1,1)=1.0_dp/grids(1)%Box_length(1)
        gprimd(2,2)=1.0_dp/grids(1)%Box_length(2)
        gprimd(2,2)=1.0_dp/grids(1)%Box_length(3)
        qphon(:)=0
        rprimd=0.0_dp
        rprimd(1,1)=grids(1)%Box_length(1)
        rprimd(2,2)=grids(1)%Box_length(2)
        rprimd(2,2)=grids(1)%Box_length(3)

        do i=1,all_PAW%N_PAW_atoms
            j=all_PAW%all_at(i)
            xred(:,i)=atoms(j)%R(:)/grids(1)%Box_length(:)
        enddo

        call pawmknhat(cplex=1,ider=0,idir=0,ipert=0,gprimd=gprimd,&
          my_natom=all_PAW%my_natom, natom=all_PAW%N_PAW_atoms, nfft=product(grids(1)%Nr_local), &
          nhatgrdim=nhatgrdim,nspden=all_PAW%nspden,ntypat=all_PAW%N_PAW_elements,pawang=all_PAW%ang,&
          pawfgrtab=all_PAW%fgrtab,pawgrnhat=gr_rho,pawnhat=rho,pawrhoij=all_PAW%rhoij, &
          pawrhoij0=all_PAW%rhoij,  pawtab=all_PAW%tab,qphon=qphon,rprimd=rprimd,xred=xred,&
          mpi_atmtab=all_PAW%my_atmtab,comm_atom=parallel%comm_space)

        do i=1,all_paw%rho_comp%n_s
            do iz_local=1, grids(1)%Nr_local(3)
                do iy_local=1, grids(1)%Nr_local(2)
                    do ix_local=1, grids(1)%Nr_local(1)
                        ic=ix_local + grids(1)%Nr_local(1)*(iy_local-1+grids(1)%Nr_local(2)*(iz_local-1))
                        all_PAW%rho_comp%of(i)%R(ix_local,iy_local,iz_local)=rho(ic,i)
                    enddo
                enddo
            enddo
            call real_to_recip(all_PAW%rho_comp%of(i), grids)
            all_PAW%rho_comp%of(i)%G=all_PAW%rho_comp%of(i)%G*grids(1)%cutden
            if(parallel%myid.eq.0) print *, 'Compensation Charge FFT intgeral:', &
            all_PAW%rho_comp%of(i)%G(1,1,1)*product(grids(1)%box_length)
            call recip_to_real(all_PAW%rho_comp%of(i), grids)
        enddo

        deallocate(rho)
        deallocate(gr_rho)
    end subroutine

    subroutine Force_or_Stress_Comp_charge(all_PAW,  potentials, grids, atoms, parallel, calc_F, calc_S, stress)
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use parallel_mod, only: parallel_task, parallel_wait
        use atom_type, only : atom_struct
        use simulation_type, only : all_PAW_struct, potentials_struct
        use m_paral_atom, only : get_my_natom, get_my_atmtab
        use m_paw_nhat, only : pawmknhat
        use fft, only: real_to_recip, recip_to_real
        use operations_3D, only: integrate_3D_R

        type(parallel_struct), intent(in) :: parallel
        type(atom_struct), intent(inout) :: atoms(:)
        type(grid_struct), intent(inout) :: grids(:)
        type(potentials_struct), intent(inout) :: potentials

        type(all_PAW_struct), intent(inout) :: all_PAw
        logical, intent(in) :: calc_F, calc_S
        integer ::   i, j, s
        real(dp), intent(out) :: stress(3,3)
        real(dp), allocatable ::  vtrial(:,:), intLM(:,:,:), intRLM(:,:,:,:)

        real(dp) ::  tmp(3), rho_s, rq, diag_s
        integer :: ix_local,iy_local,iz_local, ic, dir, ifft
        integer :: irhoij, klmn, klm, lmin, lmax, ll, ilm, isel, i_p, j_p,usexcnhat

        !Divide the atoms up by same space (spin, k-point and band redundant here)
        usexcnhat=maxval(all_PAW%tab(1:all_PAW%N_PAW_elements)%usexcnhat)
        if(all_PAW%usexcnhat.eq.1) usexcnhat=1

        allocate(vtrial(product(grids(1)%Nr_local),all_paw%rho_comp%n_s))
        vtrial=0.0_dp
        do s=1, all_paw%rho_comp%n_s
            do iz_local=1, grids(1)%Nr_local(3)
                do iy_local=1, grids(1)%Nr_local(2)
                    do ix_local=1, grids(1)%Nr_local(1)
                        ic=ix_local + grids(1)%Nr_local(1)*(iy_local-1+grids(1)%Nr_local(2)*(iz_local-1))
                        vtrial(ic,s)=real(potentials%total_local_fine%of(s)%R(ix_local,iy_local,iz_local))
                        if(usexcnhat==0) vtrial(ic,s)=vtrial(ic,s)- &
                            real(potentials%xc%of(s)%R(ix_local,iy_local,iz_local))
                    enddo
                enddo
            enddo
        enddo

        if(calc_F) then
            do i=1,size(atoms)
                if(.not.atoms(i)%update_me_force) cycle
                atoms(i)%F_comp(:)=0.0_dp
            enddo 
        endif
        if(calc_S) stress(:,:)=0.0_dp
        all_PAW%paral_atom=.true.
        all_PAW%my_atmtab=>null()
        call get_my_natom(parallel%comm_space,all_PAW%my_natom, all_PAW%N_PAW_atoms)
        call get_my_atmtab(parallel%comm_space,all_PAW%my_atmtab,all_PAW%my_atmtab_allocated, all_PAW%paral_atom, &
                            all_PAW%N_PAW_atoms)
        do i=1, all_PAW%my_natom
            if(all_PAW%paral_atom) then
                j=all_PAW%all_at(all_PAW%my_atmtab(i))
            else
                j=all_PAW%all_at(i)
            endif
            if(.not.atoms(j)%update_me_force.and.(.not.calc_S)) cycle
            if(calc_F) allocate(intLM(3,atoms(j)%PAW%an%lm_size,all_PAW%nspden))
            if(calc_S) allocate(intRLM(3,3,atoms(j)%PAW%an%lm_size,all_PAW%nspden))
            if(calc_F) intLM=0.0_dp
            if(calc_S) intRLM=0.0_dp
            !Calcualte for each LM integral V(r)dGSLM(r-R_i)/dRi
            !integration here only over local space (sum later)
            do s=1, all_PAW%nspden; do ilm=1, atoms(j)%PAW%an%lm_size
                do ic=1,atoms(j)%PAW%fgrtab%nfgd
                    ifft=atoms(j)%PAW%fgrtab%ifftsph(ic)
                    tmp(:)=vtrial(ifft,s)*atoms(j)%PAW%fgrtab%gylmgr(:,ic,ilm)
                    if(calc_F) intLM(:,ilm,s)=intLM(:,ilm,s) + tmp(:)
                    if(calc_S) then
                        do dir=1,3
                            intRLM(:,dir,ilm,s)=intRLM(:,dir,ilm,s) &
                                + 0.5_dp*tmp(dir)*atoms(j)%PAW%fgrtab%rfgd(:,ic) &
                                + 0.5_dp*tmp(:)*atoms(j)%PAW%fgrtab%rfgd(dir,ic)
                        enddo
                    endif
                enddo
                if(calc_F) intLM(:,ilm,s)=intLM(:,ilm,s)*product(grids(1)%Box_length/grids(1)%Nr) 
                if(calc_S) intRLM(:,:,ilm,s)=intRLM(:,:,ilm,s)*product(grids(1)%Box_length/grids(1)%Nr)
            enddo; enddo;

            !Calcualte for rho_ij, rho_ij*qij^LM*int
            do irhoij=1,atoms(j)%PAW%rhoij%nrhoijsel !looping through packed (nonzero) rhoij
                klmn=atoms(j)%PAW%rhoij%rhoijselect(irhoij) 
                klm =atoms(j)%PAW%tab%indklmn(1,klmn)
                lmin=atoms(j)%PAW%tab%indklmn(3,klmn)
                lmax=atoms(j)%PAW%tab%indklmn(4,klmn)
                i_p=atoms(j)%PAW%tab%indklmn(7,klmn)
                j_p=atoms(j)%PAW%tab%indklmn(8,klmn)
                do s=1, all_PAW%nspden
                    rho_s = real(atoms(j)%PAW%rho_ij(i_p,j_p,s))
                    if(i_p.ne.j_p) rho_s = 2*rho_s
                    do ll=lmin,lmax,2; do ilm=ll**2+1,(ll+1)**2
                        isel=all_PAW%ang%gntselect(ilm,klm)
                        if (isel.le.0) cycle
                        rq=rho_s*atoms(j)%PAW%tab%qijl(ilm,klmn)
                        if(calc_F) atoms(j)%F_comp(:)=atoms(j)%F_comp(:) + intLM(:,ilm,s)*rq
                        if(calc_S) stress(:,:)= stress(:,:) + intRLM(:,:,ilm,s)*rq
                    end do; end do
                enddo;
              end do 
              if(calc_F) deallocate(intLM)
              if(calc_S) deallocate(intRLM)
        enddo
        !Sum over atoms, or fill 0s
        if(calc_S) then
            call parallel_task('sum', stress(:,:) , parallel, 'all')
        endif
        if(calc_F) then
            do i=1,size(atoms)
                if(.not.atoms(i)%update_me_force) cycle
                call parallel_task('sum',  atoms(i)%F_comp(:), parallel, 'all')
            enddo
        endif

        if(calc_S) then
            diag_s=0.0_dp
            do s=1, all_PAW%nspden
                if(usexcnhat==0) then
                    diag_s=diag_s+real(integrate_3D_R( &
                    (potentials%total_local_fine%of(s)%R - &
                    potentials%xc%of(s)%R)*all_PAW%rho_comp%of(s)%R, grids(1), parallel))
                else
                    diag_s=diag_s+real(integrate_3D_R( &
                    potentials%total_local_fine%of(s)%R*all_PAW%rho_comp%of(s)%R, grids(1), parallel))
                endif
            enddo
            stress(1,1)=stress(1,1)+diag_s
            stress(2,2)=stress(2,2)+diag_s
            stress(3,3)=stress(3,3)+diag_s

        endif



    end subroutine 

    subroutine calc_core_Charge_Density(all_PAW, grids, atoms, elements, parallel)
        use odp_type, only: spin_DenPot_struct,orbital_struct
        use grids_type, only: grid_struct, inter_grid_struct
        use parallel_type, only: parallel_struct
        use parallel_mod, only: parallel_task
        use atom_type, only : atom_struct
        use element_type, only : element_struct, PAW
        use simulation_type, only : all_PAW_struct
        use m_paral_atom, only : get_my_natom, get_my_atmtab
        use m_pawrad, only: pawrad_ifromr
        use splines, only: spline3
        use grids_type, only: grid_struct, inter_grid_struct
        use constants, only : i_
        use fft, only: recip_to_real


        type(parallel_struct), intent(in) :: parallel
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(all_PAW_struct), intent(inout) :: all_PAw

        integer :: e, i,k

        if(all_PAW%N_PAW_atoms.lt.1) return
        if(.not.any(atoms(:)%update_me)) return
        
        do i=1, size(all_PAW%tncore%of)
            all_PAW%tncore%of(i)%G=0.0_dp
            grid=>grids(all_PAW%tncore%of(i)%grid)
            do k=1,size(atoms)
                e=atoms(k)%element
                if(elements(e)%PP_type.ne.PAW) cycle
                if(elements(e)%PAW%tab%usetcore.eq.1) then
                    all_PAW%tncore%of(i)%G(:,:,:) = all_PAW%tncore%of(i)%G(:,:,:) + elements(e)%Den0G(:,:,:) * &
                        exp(-i_*(grid%G(:,:,:,1)*atoms(k)%R(1) + grid%G(:,:,:,2)*atoms(k)%R(2) &
                            + grid%G(:,:,:,3)*atoms(k)%R(3)))
                endif
            enddo
            if(grid%gamma) then
                where(abs(grid%G(:,:,:,1)).gt.tiny(1.0_dp))
                    all_PAW%tncore%of(i)%G(:,:,:) = all_PAW%tncore%of(i)%G(:,:,:)*sqrt(2.0_dp)
                endwhere
            endif
            all_PAW%tncore%of(i)%G=all_PAW%tncore%of(i)%G/product(grid%Box_Length(:))/size(all_PAW%tncore%of)*grid%cutden
        enddo
            

        do i=1, size(all_PAW%tncore%of)
            call recip_to_real(all_PAW%tncore%of(i), grids)
            if(parallel%myid.eq.0) print *,  '<tncore>:', all_PAW%tncore%of(i)%G(1,1,1)*product(grid%Box_Length(:))
        enddo
        
    end subroutine

    subroutine calculate_core_xc_force(potentials, system, atoms, elements, parallel, grids)
        use constants, only: i_
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use parallel_type, only : parallel_struct
        use system_type, only : system_struct
        use grids_type, only: grid_struct
        use operations_3D, only: integrate_3D_G
        use simulation_type, only: potentials_struct

        type(parallel_struct), intent(in) :: parallel
        type(system_struct), intent(in) :: system
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(in) :: elements(:)
        type(potentials_struct), intent(in) :: potentials
        integer :: i, k, j, dir
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer:: grid
        integer :: s
        complex(dp) :: factor(size(grids(potentials%xc%of(1)%grid)%G2,1), &
                              size(grids(potentials%xc%of(1)%grid)%G2,2), &
                              size(grids(potentials%xc%of(1)%grid)%G2,3))
        grid=>grids(potentials%xc%of(1)%grid)

        k=0
        do i = 1, system%n_elements
        do j = 1, elements(i)%n_atoms_of_element
            k=k+1
            atoms(k)%F_xc_core=0.0_dp
            if(.not.allocated(atoms(k)%PAW)) cycle
            if(elements(i)%PAW%tab%usetcore.ne.1) cycle
            if(atoms(k)%update_me_force) then
                factor=exp(-i_*(grid%G(:,:,:,1)*atoms(k)%R(1) + grid%G(:,:,:,2)*atoms(k)%R(2) + grid%G(:,:,:,3)*atoms(k)%R(3)))
                factor=factor*grid%cutden
                if(grid%gamma) then
                    where(abs(grid%G(:,:,:,1)).gt.tiny(1.0_dp)) 
                        factor= factor*sqrt(2.0_dp)
                    endwhere
                endif
            do dir=1,3
                atoms(k)%F_xc_core(dir)=0.0_dp
                do s=1,size(potentials%xc%of)
                    atoms(k)%F_xc_core(dir)= atoms(k)%F_xc_core(dir) +&
                    real(integrate_3D_G( (i_ * elements(i)%Den0G(:,:,:)  &
                    * conjg(potentials%xc%of(s)%G(:,:,:)) * grid%G(:,:,:,dir) * &
                    exp(-i_*(grid%G(:,:,:,1)*atoms(k)%R(1) + grid%G(:,:,:,2)*atoms(k)%R(2) + grid%G(:,:,:,3)*atoms(k)%R(3))) ), &
                    grid, parallel))/product(grid%Box_Length(:))
                enddo
            enddo
            endif
        enddo
        enddo
            
    end subroutine

    function core_xc_Stress(potentials, grids, atoms, elements, parallel) result (S_core_xc)
        use constants, only : i_
        use simulation_type, only: potentials_struct
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use operations_3D, only: integrate_3D_G
        use element_type, only : element_struct
        use atom_type, only : atom_struct

        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(atom_struct), intent(in) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(potentials_struct), intent(in) :: potentials
        
        real(dp) :: S_core_xc(3,3)
        integer :: i, dir, dir2, e, a
        grid=>grids(potentials%xc%of(1)%grid)

        S_core_xc= 0.0_dp
        do dir=1,3
            do dir2=dir,3
                grid%work%G(:,:,:)=0.0_dp
                do a=1,size(atoms)
                    if(.not.allocated(atoms(a)%PAW)) cycle
                    e=atoms(a)%element
                    if(elements(e)%PAW%tab%usetcore.ne.1) cycle
                    grid%work%G= grid%work%G + exp(-i_*(grid%G(:,:,:,1)*atoms(a)%R(1)   + &
                                                        grid%G(:,:,:,2)*atoms(a)%R(2)   + &
                                                        grid%G(:,:,:,3)*atoms(a)%R(3))) * &
                        (elements(e)%dDen0GdG2(:,:,:)*grid%G(:,:,:,dir)*grid%G(:,:,:,dir2)*2.0_dp)
                enddo
                if(grid%gamma) then
                    where(abs(grid%G(:,:,:,1)).gt.tiny(1.0_dp)) 
                    grid%work%G = grid%work%G*sqrt(2.0_dp)
                    endwhere
                endif
                do i= 1, size(potentials%xc%of)
                    S_core_xc(dir,dir2) =  S_core_xc(dir,dir2) - real(integrate_3D_G( &
                            grid%work%G*conjg(potentials%xc%of(i)%G)*grid%cutden  &
                        , grid, parallel))/product(grid%Box_Length(:))
                enddo
            enddo        
        enddo
        do dir=1,3
            do dir2=dir,3
                S_core_xc(dir2,dir)=S_core_xc(dir,dir2) 
            enddo
        enddo
        
    end function

    subroutine calc_valance_Charge_Density_for_guess(all_PAW, grids, atoms, elements, parallel)
        use odp_type, only: spin_DenPot_struct,orbital_struct
        use grids_type, only: grid_struct, inter_grid_struct
        use parallel_type, only: parallel_struct
        use parallel_mod, only: parallel_task
        use atom_type, only : atom_struct
        use element_type, only : element_struct
        use simulation_type, only : all_PAW_struct
        use m_paral_atom, only : get_my_natom, get_my_atmtab
        use m_pawrad, only: pawrad_ifromr
        use splines, only: spline3
        use grids_type, only: grid_struct, inter_grid_struct
        use constants, only : i_
        use fft, only: recip_to_real


        type(parallel_struct), intent(in) :: parallel
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(all_PAW_struct), intent(inout) :: all_PAw

        integer :: e, i,k

        if(all_PAW%N_PAW_atoms.lt.1) return
        if(.not.any(atoms(:)%update_me)) return
        
        do i=1, size(all_PAW%tnvale%of)
            all_PAW%tnvale%of(i)%G=0.0_dp
            grid=>grids(all_PAW%tnvale%of(i)%grid)
            all_PAW%tnvale%of(i)%G(:,:,:) = 0.0_dp
            do k=1,size(atoms)
                e=atoms(k)%element
                all_PAW%tnvale%of(i)%G(:,:,:) = all_PAW%tnvale%of(i)%G(:,:,:) + elements(e)%Den1G(:,:,:) * &
                    exp(-i_*(grid%G(:,:,:,1)*atoms(k)%R(1) + grid%G(:,:,:,2)*atoms(k)%R(2) &
                        + grid%G(:,:,:,3)*atoms(k)%R(3)))
                
            enddo
            
            if(grid%gamma) then
                where(abs(grid%G(:,:,:,1)).gt.tiny(1.0_dp)) 
                    all_PAW%tnvale%of(i)%G = all_PAW%tnvale%of(i)%G*sqrt(2.0_dp)
                endwhere
            endif

            all_PAW%tnvale%of(i)%G=all_PAW%tnvale%of(i)%G*grid%cutden/product(grid%Box_Length(:))/size(all_PAW%tnvale%of)
        enddo

        do i=1, size(all_PAW%tnvale%of)
            call recip_to_real(all_PAW%tnvale%of(i), grids)
            if(parallel%myid.eq.0) print *,  '<tnvale>:', all_PAW%tnvale%of(i)%G(1,1,1)*product(grid%Box_Length(:))
        enddo
        
    end subroutine

end module