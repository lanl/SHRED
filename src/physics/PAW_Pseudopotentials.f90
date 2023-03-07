module PAW_Pseudopotentials
    use types, only : op, dp
    use constants, only : pi
    implicit none
    public ::  Calculate_Projector_PAW_RL, Calculate_deriv_Projector_PAW_RL, Calculate_Projector_PAW_RO_RL, Calc_Dij, &
    Calculate_Projector_PAW_GNL, Calculate_Projector_PAW_GNL_ORT

    contains

    subroutine Calculate_Projector_PAW_RL(Rvec, element, projector)
        use  element_type, only : element_struct
        use Spherical_Harmonic_mod, only : real_spherical_harmonic2
        use splines, only: spline3
        use  atom_type, only : atom_struct
        use linalg, only : eigh
        use m_pawrad, only: pawrad_ifromr

        type(element_struct), intent(in) :: element

        real(dp), intent(in) :: Rvec(:,:)
        real(dp),intent(inout) :: projector(:,:)
        integer  :: k, i
        integer ::  ll, mm, ln, imax 
        real(dp), allocatable :: R(:), radial_spline(:)
        
        imax=pawrad_ifromr(element%PAW%rad,element%PAW%tab%rpaw)

        allocate(R(size(Rvec,2)))
        allocate(radial_spline(size(Rvec,2)))

        R(:)=sqrt(Rvec(1,:)**2+Rvec(2,:)**2+Rvec(3,:)**2)
        if(size(projector,1).ne.size(R)) then
            print *, 'Projector size not equal to size of atom sphere'
            stop
        endif

        if(R(1).lt.(tiny(1.0_dp))) R(1)=1.0E-12_dp

        do k=1, element%n_proj
            ln=element%PAW%tab%indlmn(5,k) 
            ll=element%PAW%tab%indlmn(1,k)
            mm=element%PAW%tab%indlmn(2,k)
            radial_spline(:)=spline3(element%PAW%rad%rad(:imax), &
                                element%PAW%tab%tproj(:imax,ln),R(:))
            do i=1, size(R)
                projector(i,k) = radial_spline(i)/R(i)*real_spherical_harmonic2(ll, mm, Rvec(:,i))
           enddo
        enddo

        deallocate(R)
        deallocate(radial_spline)
        
    end subroutine

    subroutine Calculate_Projector_PAW_GNL(Gvec, element, projector)
        use  element_type, only : element_struct
        use Spherical_Harmonic_mod, only : real_spherical_harmonic2
        use splines, only: spline3_unsorted
        use  atom_type, only : atom_struct
        use linalg, only : eigh
        use m_pawrad, only: pawrad_ifromr

        type(element_struct), intent(in) :: element

        real(dp), intent(in) :: Gvec(:,:)
        real(dp),intent(inout) :: projector(:,:)
        integer  :: k, i
        integer ::  ll, mm, ln, imax 
        real(dp), allocatable :: G(:), radial_spline(:)
        
        imax=element%PAW%mqgrid_ff

        allocate(G(size(Gvec,2)))
        allocate(radial_spline(size(Gvec,2)))

        G(:)=sqrt(Gvec(1,:)**2+Gvec(2,:)**2+Gvec(3,:)**2)
        if(size(projector,1).ne.size(G)) then
            print *, 'Projector size not equal to size of reciprocal space grid'
            stop
        endif

        if(G(1).lt.(tiny(1.0_dp))) G(1)=1.0E-12_dp

        do k=1, element%n_proj
            ln=element%PAW%tab%indlmn(5,k) 
            ll=element%PAW%tab%indlmn(1,k)
            mm=element%PAW%tab%indlmn(2,k)
            radial_spline(:)=spline3_unsorted(element%PAW%qgrid_ff(:imax), &
                                         element%PAW%ffspl(:imax,1,ln), G(:)/(2.0_dp*pi))
            do i=1, size(G)
                projector(i,k) = radial_spline(i)*real_spherical_harmonic2(ll, mm, Gvec(:,i)/(2.0_dp*pi))*(4.0_dp*pi)
           enddo
        enddo

        deallocate(G)
        deallocate(radial_spline)
        
    end subroutine

    subroutine Calculate_Projector_PAW_GNL_ORT(Gvec, element, projector)
        use  element_type, only : element_struct
        use Spherical_Harmonic_mod, only : real_spherical_harmonic2
        use splines, only: spline3_unsorted
        use  atom_type, only : atom_struct
        use linalg, only : eigh
        use m_pawrad, only: pawrad_ifromr

        type(element_struct), intent(in) :: element

        real(dp), intent(in) :: Gvec(:,:)
        real(dp),intent(inout) :: projector(:,:)
        integer  :: k, i
        integer ::  ll, mm, imax 
        real(dp), allocatable :: G(:), radial_spline(:)
        
        imax=element%PAW%mqgrid_ff

        allocate(G(size(Gvec,2)))
        allocate(radial_spline(size(Gvec,2)))

        G(:)=sqrt(Gvec(1,:)**2+Gvec(2,:)**2+Gvec(3,:)**2)
        if(size(projector,1).ne.size(G)) then
            print *, 'Projector size not equal to size of reciprocal space grid'
            stop
        endif

        if(G(1).lt.(tiny(1.0_dp))) G(1)=1.0E-12_dp

        do k=1, element%n_proj
            ll=element%PAW%tab%indlmn(1,k)
            mm=element%PAW%tab%indlmn(2,k)
            radial_spline(:)=spline3_unsorted(element%PAW%qgrid_ff(:imax), &
                                         element%PAW%ffspl_Ort(:imax,1,k), G(:)/(2.0_dp*pi))
            do i=1, size(G)
                projector(i,k) = radial_spline(i)*real_spherical_harmonic2(ll, mm, Gvec(:,i)/(2.0_dp*pi))*(4.0_dp*pi)
           enddo
        enddo

        deallocate(G)
        deallocate(radial_spline)
        
    end subroutine

    subroutine Calculate_Projector_PAW_RO_RL(Rvec, element, projector)
        use  element_type, only : element_struct
        use Spherical_Harmonic_mod, only : real_spherical_harmonic2
        use splines, only: spline3
        use  atom_type, only : atom_struct
        use linalg, only : eigh
        use m_pawrad, only: pawrad_ifromr

        type(element_struct), intent(in) :: element

        real(dp), intent(in) :: Rvec(:,:)
        real(dp),intent(inout) :: projector(:,:)
        integer  :: k, i
        integer ::  ll, mm, imax 
        real(dp), allocatable :: R(:), radial_spline(:)
        
        imax=pawrad_ifromr(element%PAW%rad,element%PAW%tab%rpaw)

        allocate(R(size(Rvec,2)))
        allocate(radial_spline(size(Rvec,2)))

        R(:)=sqrt(Rvec(1,:)**2+Rvec(2,:)**2+Rvec(3,:)**2)
        if(size(projector,1).ne.size(R)) then
            print *, 'Projector size not equal to size of atom sphere'
            stop
        endif
        if(R(1).lt.(tiny(1.0_dp))) R(1)=1.0E-12_dp

        do k=1, element%n_proj
            ll=element%PAW%tab%indlmn(1,k)
            mm=element%PAW%tab%indlmn(2,k)
            radial_spline(:)=spline3(element%PAW%rad%rad(:imax), &
                                element%PAW%tproj_RO(:imax,k),R(:))
            do i=1, size(R)
                projector(i,k) = radial_spline(i)/R(i)*real_spherical_harmonic2(ll, mm, Rvec(:,i))
           enddo
        enddo

        deallocate(R)
        deallocate(radial_spline)
        
    end subroutine

    subroutine Calculate_deriv_Projector_PAW_RL(Rvec, element, projector)
        use  element_type, only : element_struct
        use  atom_type, only : atom_struct
        use splines, only: spline3ders
        use Spherical_Harmonic_mod, only : real_spherical_harmonic_and_derivatives
        use m_pawrad, only: pawrad_ifromr


        type(element_struct), intent(in) :: element
        real(dp), intent(in) :: Rvec(:,:)
        real(dp),intent(inout) :: projector(:,:,:)

        integer :: i, k, imax
        
        integer ::  ll, mm, ln
        real(dp) :: R(size(Rvec,2)),R2(size(Rvec,2))
        real(dp) :: pr_of_r(size(Rvec,2)),dpr_of_r_dr(size(Rvec,2)), d2pr_of_r_dr(size(Rvec,2))
        real(dp) :: Ylm, dYlm(3)

        R2=Rvec(1,:)**2+Rvec(2,:)**2+Rvec(3,:)**2
        projector=0.0_dp
        R=sqrt(R2)

        imax=pawrad_ifromr(element%PAW%rad,element%PAW%tab%rpaw)
        if(R(1).lt.(tiny(1.0_dp))) R(1)=1.0E-12_dp

        do k=1, element%n_proj
            ln=element%PAW%tab%indlmn(5,k) 
            ll=element%PAW%tab%indlmn(1,k)
            mm=element%PAW%tab%indlmn(2,k)
            !PAW projectors are stored as proj(r)*r on radial grid
            call spline3ders(element%PAW%rad%rad(:imax), &
                    element%PAW%tab%tproj(:imax,ln), R(:), &
                    pr_of_r(:), dpr_of_r_dr(:), d2pr_of_r_dr(:))
            do i=1, size(Rvec,2)
                call real_spherical_harmonic_and_derivatives(ll,mm, Rvec(:,i),Ylm, dYlm)
                !The minus here is because it is derivative wrt R_atom and Rvec = R-R_atom
                projector(i,k,:)=-dYlm(:)*pr_of_r(i)/R(i) &
                                 +Rvec(:,i)/R(i)**3*pr_of_r(i)*Ylm &
                                 -Rvec(:,i)/R(i)**2*dpr_of_r_dr(i)*Ylm
            enddo
        enddo  


    end subroutine    
    
    subroutine Calc_Dij(potentials, atoms, elements, all_PAW, grids, parallel)
        use parallel_type, only: parallel_struct
        use parallel_mod, only: parallel_task, parallel_wait
        use atom_type, only : atom_struct
        use element_type, only : element_struct, PAW
        use simulation_type, only : all_PAW_struct, potentials_struct
        use grids_type,only: grid_struct
        use m_pawdij, only : pawdij
        use m_paw_ij, only : paw_ij_gather
        use m_paral_atom, only : get_my_natom, get_my_atmtab

        type(parallel_struct), intent(in) :: parallel
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(all_PAW_struct), intent(inout) :: all_PAW
        type(potentials_struct), intent(in) :: potentials
        type(grid_struct), intent(in) :: grids(:)
        real(dp), allocatable :: vtrial(:,:), vxc(:,:)
        real(dp) :: gprimd(3,3), qphon(3), xred(3,all_PAW%N_PAW_atoms)
        integer :: ix_local,iy_local,iz_local, ic, at, s, klmn, i,j, e

        if(all_PAW%N_PAW_atoms.le.0) return

        gprimd=0.0_dp
        gprimd(1,1)=1.0_dp/grids(1)%Box_length(1)
        gprimd(2,2)=1.0_dp/grids(1)%Box_length(2)
        gprimd(2,2)=1.0_dp/grids(1)%Box_length(3)

        allocate(vtrial(product(grids(1)%Nr_local(:)),size(potentials%total_local_fine%of(:))))
        allocate(vxc(product(grids(1)%Nr_local(:)),size(potentials%total_local_fine%of(:))))

        do s=1, size(potentials%total_local_fine%of(:))
            do iz_local=1, grids(1)%Nr_local(3)
                do iy_local=1, grids(1)%Nr_local(2)
                    do ix_local=1, grids(1)%Nr_local(1)
                        ic=ix_local + grids(1)%Nr_local(1)*(iy_local-1+grids(1)%Nr_local(2)*(iz_local-1))
                        vtrial(ic,s)=real(potentials%total_local_fine%of(s)%R(ix_local,iy_local,iz_local))
                        vxc(ic,s)=real(potentials%xc%of(s)%R(ix_local,iy_local,iz_local))
                    enddo
                enddo
            enddo
        enddo

        i=0
        do at=1,size(atoms)
            if(allocated(atoms(at)%PAW)) then
                i=i+1
                xred(:,i)=atoms(at)%R(:)/grids(1)%Box_length(:)
            endif
        enddo


        qphon(:)=0.0_dp

        if(any(elements(:)%PP_type.eq.PAW)) then
            all_PAW%paral_atom=.true.
            all_PAW%my_atmtab=>null()
            call get_my_natom(parallel%comm_space,all_PAW%my_natom, all_PAW%N_PAW_atoms)
            call get_my_atmtab(parallel%comm_space,all_PAW%my_atmtab,all_PAW%my_atmtab_allocated, all_PAW%paral_atom, &
                                all_PAW%N_PAW_atoms)
            all_PAW%my_ij(:)%has_dij=1
            all_PAW%my_ij(:)%has_dijxc=1
            all_PAW%my_ij(:)%has_dijxc_hat=0
            if(any(all_PAW%tab(:)%usexcnhat==1).and.all_PAW%usexcnhat.ne.0) all_PAW%my_ij(:)%has_dijxc_hat=1
            all_PAW%my_ij(:)%has_dijhat=1
            all_PAW%my_ij(:)%has_dijhartree=1

            call pawdij(cplex=1, enunit=0,gprimd=gprimd,ipert=0,my_natom=all_PAW%my_natom,natom = all_PAW%N_PAW_atoms, &
                      nfft=product(grids(1)%Nr_local),nfftot=product(grids(1)%Nr(:)), nspden=all_PAW%nspden, &
                      ntypat=all_PAW%N_PAW_elements, paw_an=all_PAW%an,paw_ij=all_PAW%my_ij,pawang=all_PAW%ang, &
                      pawfgrtab=all_PAW%fgrtab,pawprtvol=0,pawrad=all_PAW%rad,pawrhoij=all_PAW%rhoij,pawspnorb=all_PAW%spnorb, &
                      pawtab=all_PAW%tab, pawxcdev=all_PAW%pawxcdev,qphon=qphon,spnorbscl=1.0_dp, &
                      ucvol=product(grids(1)%Box_length), charge=0.0_dp,vtrial=vtrial,vxc=vxc,xred=xred, &
                      mpi_atmtab=all_PAW%my_atmtab,comm_atom=parallel%comm_space,mpi_comm_grid=parallel%comm_bands)
           
            call paw_ij_gather(paw_ij_in=all_PAW%my_ij,paw_ij_gathered=all_PAW%ij,master=-1,comm_atom=parallel%comm_space)

            do at=1, size(atoms)
                e=atoms(at)%element
                if(elements(e)%PP_type.ne.PAW) cycle
                do klmn= 1, size(atoms(at)%PAW%ij%dij)
                    i=elements(e)%PAW%tab%indklmn(7,klmn)
                    j=elements(e)%PAW%tab%indklmn(8,klmn)
                    !FIX-ME:NEED TO GO THROUGH AND FIGURE OUT SPIN with PAW
                    atoms(at)%PAW%Dij(j,i,:)=atoms(at)%PAW%ij%dij(klmn,:)
                    atoms(at)%PAW%Dij(i,j,:)=atoms(at)%PAW%Dij(j,i,:)
                enddo
            enddo
        endif
        deallocate(vtrial)
        deallocate(vxc)
    end subroutine

end module