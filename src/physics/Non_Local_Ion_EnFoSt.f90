module Non_Local_ion_EnFoSt
    use types, only : op, dp
    use constants, only : pi, i_

    implicit none
    public :: NL_PP_Energy_RL, NL_PP_Force_RL, NL_PP_Stress_RL, NL_PAW_Density_matrix

    contains

subroutine NL_PAW_Density_matrix(psis, atoms, elements, grids, parallel, all_PAW, rho_ij_done)                     
        use fft, only: recip_to_real, real_to_recip
        use odp_type, only: orbital_struct
        use  parallel_mod, only: parallel_task
        use odp_type, only : field_struct
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct, PAW
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        use m_pawrhoij, only : pawrhoij_filter
        use m_paral_atom, only : get_my_natom, get_my_atmtab
        use simulation_type, only : all_PAW_struct
        use m_paw_denpot, only : pawdenpot
        use Non_Local_Ion, only : Calculate_Projector_overlaps

        type(orbital_struct), intent(in), target :: psis(:,:,:)
        type(field_struct), pointer :: psi(:)
        type(parallel_struct), intent(in) :: parallel
        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer:: grid
        type(element_struct), intent(inout), target :: elements(:)
        type(all_PAW_struct), intent(inout), target :: all_PAW

        complex(dp), allocatable :: PYpsi(:,:,:)
        real(dp) :: factor, nucdipmom(3,all_PAW%N_PAW_atoms)
        integer :: ii, j,k, s, at, e, ip, jp, klmn, jl
        real(dp) :: compch_sph,epawdc
        logical, intent(in), optional :: rho_ij_done
        logical :: rho_ij_done_

        if(all_PAW%N_PAW_atoms.lt.1) return

        rho_ij_done_=.false.
        nucdipmom=0.0_dp
        if(present(rho_ij_done)) rho_ij_done_=rho_ij_done
        if(.not.rho_ij_done_) then
                do at=1, size(atoms)
                        if(allocated(atoms(at)%PAW)) atoms(at)%PAW%rho_ij(:,:,:)=0.0_dp
                enddo
                do ii=1, size(psis,3); do j=1, size(psis,2); do k=1, size(psis,1)
                        if(sum(abs(psis(k,j,ii)%occ(:))).lt.tiny(1.0_dp)) cycle
                        psi=>psis(k,j,ii)%of
                        grid=>grids(psi(1)%grid)
                        allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(psi)))
                        PYpsi=0.0_dp
                        if(any(elements(:)%n_proj.ne.0))  then 
                                do s=1, size(psi)
                                        call Calculate_Projector_overlaps(psi(s), PYpsi(:,:,s), atoms, elements, grid, parallel)
                                enddo
                        endif
                        do s=1,size(psi)
                                factor=psis(k,j,ii)%occ(s)*psis(k,j,ii)%weight(s)
                                do at=1, size(atoms)
                                        e=atoms(at)%element
                                        if(elements(e)%PP_type.ne.PAW) cycle
                                        do jp=1, elements(e)%n_proj
                                                if(elements(e)%PAW%tab%indlmn(6,jp).ne.s) cycle
                                                jl=elements(e)%PAW%tab%indlmn(1,jp)
                                                do ip=1, elements(e)%n_proj
                                                        if(elements(e)%PAW%tab%indlmn(6,ip).ne.s) cycle
                                                        
                                                        atoms(at)%PAW%rho_ij(ip,jp,s)=atoms(at)%PAW%rho_ij(ip,jp,s) + &
                                                                PYpsi(jp,at,s)*conjg(PYpsi(ip,at,s))*factor
                                                enddo
                                        enddo
                                enddo
                        enddo
                        deallocate(PYpsi)
                enddo; enddo; enddo

                do at=1, size(atoms)
                        e=atoms(at)%element
                        if(elements(e)%PP_type.eq.PAW) then
                                call parallel_task('sum', atoms(at)%PAW%rho_ij, parallel, 'space')
                        endif
                enddo
        endif
        do at=1, size(atoms)
                e=atoms(at)%element
                if(.not.allocated(elements(e)%PAW)) cycle
                if(associated(atoms(at)%PAW%rhoij)) then !The Libpaw struct
                        do klmn=1,atoms(at)%PAW%rhoij%lmn2_size
                                ip=elements(e)%PAW%tab%indklmn(7,klmn)
                                jp=elements(e)%PAW%tab%indklmn(8,klmn)
                                !Rho_ij real
                                ! As long as hamiltonian is real only real part of  rho_ij contribute since only Sums  
                                ! with rho_ij*D_ij or rho_ij*rho_kl*Eijkl over all indicies contribute and D & E are real matricies 
                                ! Complex Hamiltonian (Spin Orbit coupling) will require using complex potentials, 
                                ! complex D, and complex Rho_ij
                                atoms(at)%PAW%rhoij%rhoij_(klmn,:)=real(atoms(at)%PAW%rho_ij(ip,jp,:))
                        enddo
                         if(at.eq.1.and.parallel%myid.eq.0) then
                                print *, 'Atom 1 PAW Rho_ij'
                                do s=1,size(atoms(at)%PAW%rho_ij,3); do ip=1, atoms(at)%PAW%rhoij%lmn_size
                                        print *, real(atoms(at)%PAW%rho_ij(ip,:ip,:))
                                enddo; enddo
                        endif
                        
                        call pawrhoij_filter(rhoij=atoms(at)%PAW%rhoij%rhoijp, &
                        rhoijselect=atoms(at)%PAW%rhoij%rhoijselect,nselect=atoms(at)%PAW%rhoij%nrhoijsel, &
                        cplex=atoms(at)%PAW%rhoij%cplex_rhoij,qphase=atoms(at)%PAW%rhoij%qphase, &
                        lmn2_size=atoms(at)%PAW%rhoij%lmn2_size, nspden=atoms(at)%PAW%rhoij%nspden,&
                        rhoij_input=atoms(at)%PAW%rhoij%rhoij_)
                endif
        enddo

        all_PAW%paral_atom=.true.
        all_PAW%my_atmtab=>null()
        call get_my_natom(parallel%comm_space,all_PAW%my_natom, all_PAW%N_PAW_atoms)
        call get_my_atmtab(parallel%comm_space,all_PAW%my_atmtab,all_PAW%my_atmtab_allocated, all_PAW%paral_atom, &
                            all_PAW%N_PAW_atoms)
        call pawdenpot(compch_sph, all_PAW%epaw,epawdc,ipert=0,ixc=all_PAW%ixc, &
         my_natom=all_PAW%my_natom,natom=all_PAW%N_PAW_atoms,nspden=all_PAW%nspden,ntypat=all_PAW%N_PAW_elements, &
         nucdipmom=nucdipmom, nzlmopt=-1,option=0,paw_an=all_PAW%an,paw_an0=all_PAW%an, &
         paw_ij=all_PAW%my_ij,pawang=all_PAW%ang,pawprtvol=0,pawrad=all_PAW%rad,pawrhoij=all_PAW%rhoij, &
         pawspnorb=all_PAW%spnorb,pawtab=all_PAW%tab,pawxcdev=all_PAW%pawxcdev,spnorbscl=0.0_dp, &
         xclevel=all_PAW%xclevel,xc_denpos=all_PAW%xc_denpos,ucvol=product(grids(1)%Box_Length),znucl=all_PAW%znucl,&
         mpi_atmtab=all_PAW%my_atmtab,comm_atom=parallel%comm_space)!,hyb_mixing,hyb_mixing_sr) ! optional arguments

         if(parallel%myid.eq.0) print *, 'compch_sph: from Abi spherical', compch_sph

        return 

end subroutine

function NL_PP_Energy_RL(psis, atoms, elements, grids, parallel) result (Een_NL)                     
        use fft, only: recip_to_real, real_to_recip
        use odp_type, only: orbital_struct, field_struct
        use  parallel_mod, only: parallel_task
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct, HGH, Real_local, PAW
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        use Non_Local_Ion, only : Calculate_Projector_overlaps 
        type(orbital_struct), intent(in), target :: psis(:,:,:)

        type(field_struct), pointer :: psi(:)

        type(parallel_struct), intent(in) :: parallel
        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer:: grid

        type(element_struct), intent(inout), target :: elements(:)
        real(dp) :: KB
        type(atom_struct), pointer :: atom
        type(element_struct), pointer :: element
        complex(dp), allocatable :: PYpsi(:,:,:)
        integer :: ii, j,k, l, at, i_p, m, s, i
        real(dp) ::  Een_NL
        real(dp), allocatable :: Een_NL_psi(:)
        complex(dp), allocatable :: Mij(:,:,:), MijPYpsi(:,:)

        Een_NL=0.0_dp
        do ii=1, size(psis,3); do j=1, size(psis,2); do k=1, size(psis,1)
                if(sum(abs(psis(k,j,ii)%occ(:))).lt.tiny(1.0_dp)) cycle
                psi=>psis(k,j,ii)%of
                grid=>grids(psi(1)%grid)
                if(any(elements(:)%n_proj.ne.0))  then 
                        allocate(Een_NL_psi(size(psi)))
                        Een_NL_psi=0.0_dp
                        allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(psi)))
                        PYpsi=0.0_dp
                        do s=1, size(psi)
                                call Calculate_Projector_overlaps(psi(s), PYpsi(:,:,s), atoms, elements, grid, parallel) 
                        enddo
                        do at=1, size(atoms)
                        atom=>atoms(at)
                        element=>elements(atom%element)
                        if(element%PA_type.ne.Real_local) cycle
                        if(element%PP_type.eq.HGH) then
                                do i_p=1, element%n_proj
                                i=element%HGH%i(i_p)
                                l=element%HGH%l(i_p)
                                m=element%HGH%m(i_p)
                                do s=1, size(psi)
                                        if(element%HGH%h_or_k(i_p).eq.'h') then
                                                        KB=element%HGH%EKB(i,l+1)
                                        else
                                                        KB=element%HGH%KKB(i,l+1)*(s*1.0_dp-1.5_dp)*m
                                        endif
                                        Een_NL_psi(s)= Een_NL_psi(s)+ abs(PYpsi(i_p,at,s))**2*KB
                                enddo
                                enddo               
                        else if(element%PP_type.eq.PAW) then
                                !if(all(abs(grid%k).lt.tiny(1.0_dp)).and.all(abs(grid%A).lt.tiny(1.0_dp))) cycle
                                cycle !This is double counted
                                allocate(Mij(element%n_proj,element%n_proj, size(psi)))
                                do s=1, size(psi); do i=1, element%n_proj  
                                        Mij(:,i,s)=-i_*((grid%A(1)+grid%k(1))*element%PAW%tab%nabla_ij(1,:,i) + &
                                                        (grid%A(2)+grid%k(2))*element%PAW%tab%nabla_ij(2,:,i) + &
                                                        (grid%A(3)+grid%k(3))*element%PAW%tab%nabla_ij(3,:,i))+ &
                                                sum((grid%A(:)+grid%k(:))**2)*element%PAW%sij(:,i)*0.5_dp
                                        Mij(:,i,s)=0.0_dp
                                enddo; enddo
                                allocate(MijPYpsi(element%n_proj, size(psi)))
                                do s=1, size(psi)
                                        do i_p=1, element%n_proj
                                                MijPYpsi(i_p,s)=sum(Mij(:,i_p,s)*PYpsi(:element%n_proj,at,s))
                                        enddo
                                        Een_NL_psi(s)= Een_NL_psi(s)+ real(sum(conjg(MijPYpsi(:,s))*PYpsi(:element%n_proj,at,s)))
                                enddo                                
                                deallocate(Mij)
                                deallocate(MijPYpsi)
                        else
                                print *, 'Unknown PP type: ', element%PP_type; stop
                        endif       
                        enddo
                        deallocate(PYpsi)
                        do s=1, size(psi);
                                Een_NL=Een_NL+Een_NL_psi(s)*psis(k,j,ii)%occ(s)*psis(k,j,ii)%weight(s)
                        enddo
                        deallocate(Een_NL_psi)
                endif

        enddo; enddo; enddo
        call parallel_task('sum', Een_NL, parallel, 'space')

        return 

end function

subroutine NL_PP_Force_RL(psis, atoms, elements, potentials, grids, parallel, all_PAW)                     
        use fft, only: recip_to_real, real_to_recip
        use odp_type, only: orbital_struct, eigen_vector,field_struct
        use  parallel_mod, only: parallel_task
        use odp, only : allocate_field, deallocate_field
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct, HGH, Real_local, PAW
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        use simulation_type, only : potentials_struct
        use Apply_Hamiltonian, only: Apply_SinvH
        use Non_local_ion, only : Calculate_Projector_overlaps, Calculate_deriv_Projector_overlaps
        use simulation_type, only : all_PAW_struct
        use operations_3D, only : integrate_3D_R

        type(orbital_struct), intent(in), target :: psis(:,:,:)

        type(field_struct), pointer :: psi(:), Hpsi(:)

        type(parallel_struct), intent(in) :: parallel
        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        type(potentials_struct), intent(in) :: potentials
        type(all_PAW_struct), intent(inout) :: all_PAW

        type(element_struct), intent(inout), target :: elements(:)
        real(dp) :: KB
        type(atom_struct), pointer :: atom
        type(element_struct), pointer :: element
        complex(dp), allocatable :: PYpsi(:,:,:), dPYpsi(:,:,:,:), PYHpsi(:,:,:)
        integer :: ii, j,k, l, at, i_p, m, s, s2, i, j_p, dir
        real(dp) :: Fen_NL(size(atoms),3)
        real(dp), allocatable :: Fen_NL_psi(:,:,:)
        complex(dp), allocatable :: MijPYpsi(:,:), TijPYpsi(:,:,:), SijPYHpsi(:,:)
        complex(dp), allocatable :: Mij(:,:,:)
        complex(dp) :: check
        Fen_NL=0.0_dp
        do ii=1, size(psis,3); do j=1, size(psis,2); do k=1, size(psis,1)
                if(sum(abs(psis(k,j,ii)%occ(:))).lt.tiny(1.0_dp)) cycle
                psi=>psis(k,j,ii)%of
                grid=>grids(psi(1)%grid)
                if(any(elements(:)%n_proj.ne.0))  then 
                        if(psis(k,j,ii)%td_type.ne.eigen_vector.and.all_PAW%N_PAW_atoms.gt.0) then
                                allocate(Hpsi(size(psi)))
                                do i=1, size(psi)
                                call allocate_field(Hpsi(i), grid, parallel)
                                enddo
                                Hpsi%grid=psi%grid
                                call Apply_SinvH(psi(:), Hpsi, grids, potentials, atoms, elements, &
                                                                parallel, all_PAW, calc_G=.false., with_CG=.true.)
                                allocate(TijPYpsi(3,maxval(elements(:)%n_proj), size(psi)))
                                allocate(SijPYHpsi(maxval(elements(:)%n_proj), size(psi)))
                                allocate(PYHpsi(maxval(elements(:)%n_proj), size(atoms), size(psi)))
                                PYHpsi=0.0_dp
                        endif
                        allocate(Fen_NL_psi(size(atoms),3,size(psi)))
                        Fen_NL_psi=0.0_dp
                        allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(psi)))
                        allocate(dPYpsi(3,maxval(elements(:)%n_proj), size(atoms), size(psi)))
                        PYpsi=0.0_dp
                        dPYpsi=0.0_dp
                        do s=1, size(psi)
                                call Calculate_Projector_overlaps(psi(s), PYpsi(:,:,s), atoms, elements, grid, parallel)
                                call Calculate_deriv_Projector_overlaps(psi(s), dPYpsi(:,:,:,s), atoms, elements, grid, parallel) 
                                if(psis(k,j,ii)%td_type.ne.eigen_vector.and.all_PAW%N_PAW_atoms.gt.0) then
                                        call Calculate_Projector_overlaps(Hpsi(s), PYHpsi(:,:,s), atoms, elements, grid, parallel)
                                endif
                        enddo
                        if(psis(k,j,ii)%td_type.ne.eigen_vector.and.all_PAW%N_PAW_atoms.gt.0) then
                                do i=1, size(psi)
                                        call deallocate_field(Hpsi(i),grids(Hpsi(i)%grid))
                                enddo
                                deallocate(Hpsi)
                        endif

                        do at=1, size(atoms)
                                atom=>atoms(at)
                                element=>elements(atom%element)
                                if(element%PA_type.ne.Real_local) cycle
                                if(element%PP_type.eq.PAW) then
                                        allocate(Mij(element%n_proj,element%n_proj, size(psi)))

                                       ! if(any(abs(grid%k).gt.tiny(1.0_dp)).or.any(abs(grid%A).gt.tiny(1.0_dp))) then
                                       !         do s=1, size(psi); do i=1, element%n_proj  
                                       !                 Mij(:,i,s)=atom%PAW%Dij(:,i,s)  &
                                       !                          -i_*((grid%A(1)+grid%k(1))*element%PAW%tab%nabla_ij(1,:,i) + &
                                       !                              (grid%A(2)+grid%k(2))*element%PAW%tab%nabla_ij(2,:,i) + &
                                       !                              (grid%A(3)+grid%k(3))*element%PAW%tab%nabla_ij(3,:,i))+ &
                                       !                         sum((grid%A(:)+grid%k(:))**2)*element%PAW%sij(:,i)*0.5_dp
                                       !         enddo; enddo
                                       ! else
                                                Mij=atom%PAW%Dij +i_*0.0_dp
                                        !endif
                                        
                                        if(psis(k,j,ii)%td_type.eq.eigen_vector) then
                                                do s=1, size(psi); do j_p=1, element%n_proj  
                                                        Mij(:,j_p,s)=Mij(:,j_p,s)-psis(k,j,ii)%eig(s)*element%PAW%sij(:,j_p)
                                                enddo; enddo
                                        endif

                                        allocate(MijPYpsi(element%n_proj, size(psi)))
                                        do s=1, size(psi)
                                                do i_p=1, element%n_proj
                                                        MijPYpsi(i_p,s)=sum(Mij(i_p,:,s)*PYpsi(:element%n_proj,at,s))
                                                enddo
                                        enddo
                                        if(psis(k,j,ii)%td_type.ne.eigen_vector) then
                                                do s=1, size(psi); do i_p=1, element%n_proj
                                                        do dir=1,3               
                                                        TijPYpsi(dir,i_p,s)=-sum(element%PAW%tab%nabla_ij(dir,i_p,:) &
                                                                                        *PYpsi(:element%n_proj,at,s))
                                                        enddo 
                                                        SijPYHpsi(i_p,s)=-sum(element%PAW%sij(i_p,:) &
                                                                *PYHpsi(:element%n_proj,at,s))
                                                enddo; enddo
                                        endif
                                endif
                                do i_p=1, element%n_proj
                                        if(element%PP_type.eq.HGH) then
                                                i=element%HGH%i(i_p)
                                                l=element%HGH%l(i_p)
                                                m=element%HGH%m(i_p)
                                                do s=1, size(psi); do s2=1, size(psi)
                                                        if(element%HGH%h_or_k(i_p).eq.'h') then
                                                                if(s.eq.s2) then
                                                                        KB=element%HGH%EKB(i,l+1)
                                                                else
                                                                        cycle
                                                                endif
                                                        else
                                                                if(s.eq.s2) then
                                                                        KB=element%HGH%KKB(i,l+1)*(s*1.0_dp-1.5_dp)*m
                                                                else
                                                                        KB=element%HGH%KKB(i,l+1) &
                                                                        *0.5_dp*sqrt(1.0_dp*(l+m+1))*sqrt(1.0_dp*(l-m))
                                                                endif
                                                        endif
                                                        Fen_NL_psi(at,:,s)= Fen_NL_psi(at,:,s)  &
                                                        -2.0_dp*real(conjg(dPYpsi(:,i_p,at,s))*PYpsi(i_p,at,s))*KB
                                                enddo; enddo
                                        else if(element%PP_type.eq.PAW) then
                                                do s=1, size(psi)
                                                        
                                                        Fen_NL_psi(at,:,s)= Fen_NL_psi(at,:,s) &
                                                        -2.0_dp*real(conjg(dPYpsi(:,i_p,at,s))*MijPYpsi(i_p,s))
                                                        if(psis(k,j,ii)%td_type.ne.eigen_vector) then     
                                                        Fen_NL_psi(at,:,s)= Fen_NL_psi(at,:,s) &
                                                                -2.0_dp*real(conjg(dPYpsi(:,i_p,at,s))*SijPYHpsi(i_p,s)) &
                                                                -2.0_dp*real(conjg(TijPYpsi(:,i_p,s))*PYHpsi(i_p,at,s))
                                                        endif
                                                enddo
                                        else
                                                print *, 'Unknown PP type: ', element%PP_type; stop
                                        endif       
                                enddo
                                if(element%PP_type.eq.PAW ) deallocate(Mij)
                                if(element%PP_type.eq.PAW )deallocate(MijPYpsi)
                                
                        enddo
                        deallocate(PYpsi)
                        deallocate(dPYpsi)
                        if(psis(k,j,ii)%td_type.ne.eigen_vector.and.all_PAW%N_PAW_atoms.gt.0) then
                                deallocate(PYHpsi)
                                deallocate(TijPYpsi)
                                deallocate(SijPYHpsi)
                        endif
                        do s=1, size(psi);
                                Fen_NL(:,:)= Fen_NL(:,:)+ Fen_NL_psi(:,:,s)*psis(k,j,ii)%occ(s)*psis(k,j,ii)%weight(s)
                        enddo
                        deallocate(Fen_NL_psi)
                endif
        enddo; enddo; enddo

        call parallel_task('sum', Fen_NL, parallel, 'space')
        atoms(:)%F_nonloc(1)=atoms(:)%F_nonloc(1)+Fen_NL(:,1)
        atoms(:)%F_nonloc(2)=atoms(:)%F_nonloc(2)+Fen_NL(:,2)
        atoms(:)%F_nonloc(3)=atoms(:)%F_nonloc(3)+Fen_NL(:,3)
        return 

end subroutine

subroutine NL_PP_Stress_RL(psis, atoms, elements, potentials, grids, parallel, all_PAW, SNL)                     
        use fft, only: recip_to_real, real_to_recip
        use odp_type, only: orbital_struct, eigen_vector,field_struct
        use  parallel_mod, only: parallel_task, parallel_wait
        use odp, only : allocate_field, deallocate_field
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct, HGH, Real_local,PAW
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        use simulation_type, only : potentials_struct
        use Apply_Hamiltonian, only: Apply_SinvH, Apply_H
        use Non_local_ion, only : Calculate_Projector_overlaps, Calculate_deriv_Projector_R_overlaps
        use operations_3D, only : Integrate_3D_R
        use simulation_type, only : all_PAW_struct

        type(orbital_struct), intent(in), target :: psis(:,:,:)

        type(field_struct), pointer :: psi(:)
        type(field_struct), allocatable :: Hpsi(:)
        type(parallel_struct), intent(in) :: parallel
        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        type(potentials_struct), intent(in) :: potentials
        type(element_struct), intent(inout), target :: elements(:)
        real(dp) ,intent(inout) :: SNL(3,3)
        type(all_PAW_struct), intent(inout) :: all_PAW

        real(dp) :: KB
        type(atom_struct), pointer :: atom
        type(element_struct), pointer :: element
        complex(dp), allocatable :: PYpsi(:,:,:), dPYpsi(:,:,:,:,:), PYHpsi(:,:,:)
        complex(dp), allocatable :: MijPYpsi(:,:), TijPYpsi(:,:,:,:), SijPYHpsi(:,:)
        complex(dp), allocatable :: Mij(:,:,:)

        integer :: ii, j,k, l, at, i_p, m, s, i, dir,dir2, j_p
        real(dp), allocatable :: Sen_NL_psi(:,:,:,:)
        real(dp) :: SNL_here(3,3), energy

        SNL_here=0.0_dp
        do ii=1, size(psis,3); do j=1, size(psis,2); do k=1, size(psis,1)
                if(sum(abs(psis(k,j,ii)%occ(:))).lt.tiny(1.0_dp)) cycle
                psi=>psis(k,j,ii)%of
                grid=>grids(psi(1)%grid)
                if(any(elements(:)%n_proj.ne.0))  then 
                        allocate(Sen_NL_psi(size(atoms),3,3,size(psi)))
                        Sen_NL_psi=0.0_dp
                        allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(psi)))
                        allocate(dPYpsi(3,3,maxval(elements(:)%n_proj), size(atoms), size(psi)))
                        PYpsi=0.0_dp
                        dPYpsi=0.0_dp
                        if(psis(k,j,ii)%td_type.ne.eigen_vector.and. &
                                    any(elements(:)%PP_type.eq.PAW)) then
                                allocate(Hpsi(size(psi)))
                                do i=1, size(psi)
                                call allocate_field(Hpsi(i), grid, parallel)
                                enddo
                                Hpsi%grid=psi%grid
                                
                                call Apply_SinvH(psi(:), Hpsi, grids, potentials, atoms, elements, &
                                                                parallel, all_PAW, calc_G=.false., with_CG=.true.)
                                allocate(TijPYpsi(3, 3, maxval(elements(:)%n_proj), size(psi)))
                                allocate(SijPYHpsi(maxval(elements(:)%n_proj), size(psi)))
                                allocate(PYHpsi(maxval(elements(:)%n_proj), size(atoms), size(psi)))
                                PYHpsi=0.0_dp
                        else
                        endif
                        do s=1, size(psi)
                                call Calculate_Projector_overlaps(psi(s), PYpsi(:,:,s), atoms, elements, grid, parallel)
                                call Calculate_deriv_Projector_R_overlaps(psi(s), dPYpsi(:,:,:,:,s),  &
                                        atoms, elements, grid, parallel)
                                if(allocated(Hpsi)) then
                                        call Calculate_Projector_overlaps(Hpsi(s), PYHpsi(:,:,s), atoms, elements, grid, parallel)
                                        call deallocate_field(Hpsi(s),grids(Hpsi(s)%grid))
                                endif
                        enddo
                        if(allocated(HPsi))  deallocate(Hpsi)

                        do at=1, size(atoms)
                        atom=>atoms(at)
                        element=>elements(atom%element)

                                if(element%PA_type.ne.Real_local) cycle

                                if(element%PP_type.eq.PAW) then
                                        allocate(Mij(element%n_proj,element%n_proj, size(psi)))
                                        !if(any(abs(grid%k).gt.tiny(1.0_dp)).or.any(abs(grid%A).gt.tiny(1.0_dp))) then
                                       !         do s=1, size(psi); do i=1, element%n_proj  
                                       !                 Mij(:,i,s)=atom%PAW%Dij(:,i,s) &
                                       !                 -i_*((grid%A(1)+grid%k(1))*element%PAW%tab%nabla_ij(1,:,i) + &
                                       !                 (grid%A(2)+grid%k(2))*element%PAW%tab%nabla_ij(2,:,i) + &
                                       !                 (grid%A(3)+grid%k(3))*element%PAW%tab%nabla_ij(3,:,i))+ &
                                       !                         sum((grid%A(:)+grid%k(:))**2)*element%PAW%sij(:,i)*0.5_dp
                                       !         enddo; enddo
                                       ! else
                                                Mij=atom%PAW%Dij +i_*0.0_dp
                                       ! endif
                                        
                                        if(psis(k,j,ii)%td_type.eq.eigen_vector) then
                                                do s=1, size(psi); do j_p=1, element%n_proj  
                                                        Mij(:,j_p,s)=Mij(:,j_p,s)-psis(k,j,ii)%eig(s)*element%PAW%sij(:,j_p)
                                                enddo; enddo
                                        endif
                                        allocate(MijPYpsi(element%n_proj, size(psi)))
                                        do s=1, size(psi)
                                                do i_p=1, element%n_proj
                                                        MijPYpsi(i_p,s)=sum(Mij(:,i_p,s)*PYpsi(:element%n_proj,at,s))
                                                enddo
                                        enddo
                                        if(allocated(TijPYpsi)) then
                                                do s=1, size(psi); do i_p=1, element%n_proj
                                                        !Tij+Tij^t= - delta_dir,dir2 Sij
                                                        do dir=1,3; do dir2=1,3             
                                                        TijPYpsi(dir2,dir,i_p,s)=sum( &
                                                        element%PAW%tab%x_nabla_ij(dir2,dir,i_p,:) &
                                                        *PYpsi(:element%n_proj,at,s))
                                                      ! if(dir.eq.dir2) print *, element%PAW%tab%x_nabla_ij(dir2,dir,i_p,1) + &
                                                      !          element%PAW%tab%x_nabla_ij(dir2,dir,1,i_p) &
                                                      !          +element%PAW%sij(i_p,1)
                                                        enddo; enddo
                                                        
                                                        SijPYHpsi(i_p,s)=-sum(element%PAW%sij(:,i_p) &
                                                                *PYHpsi(:element%n_proj,at,s))
                                                enddo; enddo
                                        endif
                                endif

                                do i_p=1, element%n_proj
                                !Terms with dPYpsi have opposite sign from paper formulas because it is derivative wrt atom position 
                                !not r-Ra
                                        if(element%PP_type.eq.HGH) then
                                                i=element%HGH%i(i_p)
                                                l=element%HGH%l(i_p)
                                                m=element%HGH%m(i_p)
                                                do s=1, size(psi)
                                                        if(element%HGH%h_or_k(i_p).eq.'h') then
                                                                        KB=element%HGH%EKB(i,l+1)
                                                        
                                                        else
                                                                        KB=element%HGH%KKB(i,l+1)*(s*1.0_dp-1.5_dp)*m
                                                        endif

                                        Sen_NL_psi(at,:,:,s)= Sen_NL_psi(at,:,:,s) -2.0_dp*real(conjg(dPYpsi(:,:,i_p,at,s)) &
                                                                                        *PYpsi(i_p,at,s))*KB
                                                        do dir=1,3 
                                                Sen_NL_psi(at,dir,dir,s)=Sen_NL_psi(at,dir,dir,s)-KB*abs(PYpsi(i_p,at,s))**2
                                                        enddo                                        
                                                enddo
                                        else if(element%PP_type.eq.PAW) then
                                                
                                                do s=1, size(psi)
                                                        Sen_NL_psi(at,:,:,s) =  Sen_NL_psi(at,:,:,s) &
                                                                -2.0_dp*real(conjg(dPYpsi(:,:,i_p,at,s))*MijPYpsi(i_p,s))
                                                        do dir=1,3 
                                                                Sen_NL_psi(at,dir,dir,s)=Sen_NL_psi(at,dir,dir,s) &
                                                                        +real(conjg(PYpsi(i_p,at,s))*MijPYpsi(i_p,s))
                                                        enddo
                                                        if(psis(k,j,ii)%td_type.ne.eigen_vector) then
                                                                Sen_NL_psi(at,:,:,s) =  Sen_NL_psi(at,:,:,s) &
                                                                -2.0_dp*real(conjg(dPYpsi(:,:,i_p,at,s))*SijPYHpsi(i_p,s)) &
                                                                +2.0_dp*real(conjg(TijPYpsi(:,:,i_p,s))*PYHpsi(i_p,at,s))
                                                                do dir=1,3 
                                                                 !       Sen_NL_psi(at,dir,dir,s)=Sen_NL_psi(at,dir,dir,s)  &
                                                                 !       +2.0_dp*real(conjg(PYpsi(i_p,at,s))*SijPYHpsi(i_p,s))
                                                                enddo
                                                        endif
                                                enddo
                                        else
                                                print *, 'Unknown PP type: ', element%PP_type; stop
                                        endif       
                                enddo   
                                if(element%PP_type.eq.PAW) then
                                        deallocate(Mij)  
                                        deallocate(MijPYpsi)            
                                endif     
                        enddo
                        deallocate(PYpsi)
                        deallocate(dPYpsi)
                        if(allocated(TijPYpsi))         deallocate(TijPYpsi)
                        if(allocated(SijPYHpsi))        deallocate(SijPYHpsi)
                        if(allocated(PYHpsi))           deallocate(PYHpsi)
                        do s=1, size(psi); 
                                do at=1, size(atoms); 
                                SNL_here(:,:)= SNL_here(:,:)+ Sen_NL_psi(at,:,:,s)*psis(k,j,ii)%occ(s)*psis(k,j,ii)%weight(s)
                                enddo;
                         enddo
                        deallocate(Sen_NL_psi)
                endif
        enddo; enddo; enddo
        call parallel_task('sum', SNL_here, parallel, 'space')
        SNL=SNL+SNL_here
        return 

end subroutine

end module
