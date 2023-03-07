module Non_Local_ion
    use types, only : op, dp
    use constants, only : pi, i_

    implicit none
    public :: pseudopotential_update, allocate_and_fill_Ono_Hirose_matrix, Apply_projectors_RL, &
              Calculate_Projector_overlaps, Apply_S_power, Apply_S_inverse, Apply_S_rotation, Apply_S_RR, &
              Calculate_Projector_GNL, Calculate_Projector_GNL_ORT

    contains

    subroutine pseudopotential_update(atoms, elements, grid, fine_grid, parallel, all_PAW)
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct, Real_local, local_only, PAW, Reciprocal
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        use shared_type, only : shared_fence
        use  parallel_mod, only : parallel_wait, parallel_task
        use m_paw_finegrid, only : pawgylm
        use simulation_type, only : all_PAW_struct

        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(in) :: grid, fine_grid
        type(parallel_struct), intent(in) :: parallel
        type(all_PAW_struct), intent(inout) :: all_PAW


        type(element_struct), intent(inout) :: elements(:)
        integer  :: e, j,  k
        integer :: i_coarse, i_dense

        real(dp), allocatable, target :: dense_projector(:,:,:), coarse_projector(:,:,:)
        real(dp), allocatable :: projector_dense(:,:,:) , projector_dense_rot(:,:)
        real(dp), allocatable :: projector_coarse(:,:,:)

        real(dp), allocatable :: Rs_dense(:,:)
        real(dp), allocatable :: Oij(:,:)
        integer, allocatable :: map_dense(:,:), map_coarse(:,:)
        integer :: n_P, dir, ns_dense, ns_coarse

        do e=1, size(elements)
            if(elements(e)%PA_type.eq.local_only.or.elements(e)%n_proj.lt.1) cycle
            if(elements(e)%PA_type.eq.Real_local) then
                !=====================================================================================================================
                !First, We calculate the projectors on a full sphere around each atom, result is stores once on each node in shared memory
                !=====================================================================================================================

                do j=1, size(atoms)
                        if(e.eq.atoms(j)%element) then
                                call shared_fence(atoms(j)%RL%proj_shared, parallel)
                                if(elements(e)%PP_type.eq.PAW) &
                                        call shared_fence(atoms(j)%PAW%proj_shared_ort, parallel)
                        endif
                enddo

                if(any(atoms(:)%update_me_force)) then
                        do j=1, size(atoms)
                                if(e.eq.atoms(j)%element) then
                                        call shared_fence(atoms(j)%RL%deriv_proj_shared, parallel)
                                        if(elements(e)%PP_type.eq.PAW) &
                                                call shared_fence(atoms(j)%PAW%deriv_proj_shared_ort, parallel)
                                endif
                        enddo
                endif
                allocate(dense_projector(elements(e)%RL%max_sphere_dense(1), &
                                         elements(e)%RL%max_sphere_dense(2), &
                                         elements(e)%RL%max_sphere_dense(3)))
                allocate(coarse_projector(elements(e)%RL%max_sphere(1), &
                                         elements(e)%RL%max_sphere(2), &
                                         elements(e)%RL%max_sphere(3)))

                allocate(Rs_dense(3, product(elements(e)%RL%max_sphere_dense)))
                allocate(map_dense(4, product(elements(e)%RL%max_sphere_dense)))
                allocate(map_coarse(4, product(elements(e)%RL%max_sphere)))

                do j=1, size(atoms) 
                        if(atoms(j)%element.ne.e) cycle
                        if(.not.atoms(j)%update_me) cycle
                        if(elements(e)%PP_type.eq.PAW) atoms(j)%PAW%obar_i=0.0_dp
                        if(.not.atoms(j)%mine) cycle
                        !Set 1 for projector add 3 for derivatives
                        n_P=1
                        if(atoms(j)%update_me_force) n_P=4

                        !map 3D cube around atom on fine grid to 1D array of the sphere around atom 
                        !Careful that the full grid here is indeed the coarse grid
                        call map_atom_sphere(ns_dense, map_dense, elements(e), atoms(j), &
                                 grid, elements(e)%RL%dense_grid, elements(e)%RL%max_sphere_dense, Rs_dense, elements(e)%RL%radius)

                        !map 3D cube around atom on fine grid to 1D array of the sphere around atom 
                        !Careful not to sort the coarse elements, i.e. dont give an Rs, since you wont sort them 
                        ! when you make the local map, otherwise we need another map for sorted
                        call map_atom_sphere(ns_coarse, map_coarse, elements(e), atoms(j), &
                                 grid, elements(e)%RL%coarse_grid, elements(e)%RL%max_sphere, radius=elements(e)%RL%radius)

                        !Calculate projector (with Diagonal D or EKB) on fine sphere (1D array) 
                        allocate(projector_dense(ns_dense,elements(e)%n_proj, n_P))
                        call Calculate_Projector_RL(Rs_dense(:,:ns_dense), &
                                        elements(e), projector_dense(:,:,1))

                        if(atoms(j)%update_me_force) call Calculate_deriv_Projector_RL( &
                                Rs_dense(:,:ns_dense), elements(e), projector_dense(:,:,2:n_P))


                        !Smooth the projector from fine to the coarse grid
                        !Have to project back onto the cube first, then back to the coarse 1D sphere
                        do k=1,elements(e)%n_proj
                                dense_projector=0.0_dp
                                do i_dense=1, ns_dense
                                        dense_projector(map_dense(1,i_dense),map_dense(2,i_dense),map_dense(3,i_dense)) &
                                                = projector_dense(i_dense,k,1)
                                enddo
                                
                                call Ono_Hirose(elements(e)%RL%Beta_x, elements(e)%RL%Beta_y, elements(e)%RL%Beta_z, &
                                        elements(e)%RL%max_sphere, elements(e)%RL%max_sphere_dense,&
                                        dense_projector, coarse_projector)

                                do i_coarse=1, ns_coarse
                                        atoms(j)%RL%projector(i_coarse,k)= &
                                        coarse_projector(map_coarse(1,i_coarse),map_coarse(2,i_coarse),map_coarse(3,i_coarse)) 
                                enddo

                                if(atoms(j)%update_me_force) then
                                        do dir=1,3
                                                dense_projector=0.0_dp
                                                do i_dense=1, ns_dense
                                                dense_projector(map_dense(1,i_dense),map_dense(2,i_dense),map_dense(3,i_dense))= &
                                                                projector_dense(i_dense,k,dir+1)
                                                enddo
                                                call Ono_Hirose(elements(e)%RL%Beta_x, elements(e)%RL%Beta_y, &
                                                        elements(e)%RL%Beta_z, &
                                                        elements(e)%RL%max_sphere, elements(e)%RL%max_sphere_dense,&
                                                        dense_projector, coarse_projector)
                                                do i_coarse=1, ns_coarse
                                                        atoms(j)%RL%deriv_projector(i_coarse, dir,k)= &
                                coarse_projector(map_coarse(1,i_coarse),map_coarse(2,i_coarse),map_coarse(3,i_coarse)) 
                                                enddo
                                        enddo
                                endif
                        enddo

                        deallocate(projector_dense)

                        !For PAW also need the projectors in the S-Diagonal Rotation
                        if(elements(e)%PP_type.eq.PAW) then
                                !Rotate Fine first
                                allocate(projector_dense_rot(ns_dense,elements(e)%n_proj))

                                call Calculate_Projector_RL(Rs_dense(:,:ns_dense), &
                                        elements(e), projector_dense_rot(:,:), RO=.true.)

                                !Smooth the rotated dense projector to the coarse grid
                                do k=1,elements(e)%n_proj
                                        dense_projector=0.0_dp
                                        do i_dense=1, ns_dense
                                                dense_projector(map_dense(1,i_dense),map_dense(2,i_dense),map_dense(3,i_dense))=&
                                                        projector_dense_rot(i_dense,k)
                                        enddo

                                        call Ono_Hirose(elements(e)%RL%Beta_x, elements(e)%RL%Beta_y, elements(e)%RL%Beta_z, &
                                                elements(e)%RL%max_sphere, elements(e)%RL%max_sphere_dense,&
                                                dense_projector, coarse_projector)

                                        do i_coarse=1, ns_coarse
                                                atoms(j)%PAW%projector_ort(i_coarse,k)= &
                                        coarse_projector(map_coarse(1,i_coarse),map_coarse(2,i_coarse),map_coarse(3,i_coarse)) 
                                        enddo
                                enddo
                                deallocate(projector_dense_rot)

                                !atoms(j)%PAW%obar_i(:)=elements(e)%PAW%o_i(:)
                                if(.true.) then
                                        !Redo S rotation on smooth grid (as described in the paper)
                                        allocate(Oij(elements(e)%n_proj,elements(e)%n_proj))
                                        Oij=0.0_dp
                                        do k=1,elements(e)%n_proj
                                                Oij(k,k)=elements(e)%PAW%o_i(k)
                                        enddo

                                        allocate(projector_coarse(ns_coarse,elements(e)%n_proj,1))
                                        projector_coarse(:ns_coarse,:,1)= &
                                                atoms(j)%PAW%projector_ort(:ns_coarse,:)
                                        call Projector_S_transformation(grid%dR, ns_coarse, elements(e)%n_proj, &
                                                projector_coarse(:ns_coarse,:,1), atoms(j)%PAW%projector_ort(:ns_coarse,:), &
                                                Oij, atoms(j)%PAW%obar_i)
                                        !Avoid Singularity in S^-1/2 application
                                        atoms(j)%PAW%obar_i(:)=max(atoms(j)%PAW%obar_i(:),-0.997_dp)
                                        deallocate(projector_coarse)
                                        deallocate(Oij)
                                endif
                        endif  
                enddo
                do j=1, size(atoms) 
                        if(atoms(j)%element.ne.e) cycle
                        if(.not.atoms(j)%update_me) cycle
                        if(elements(e)%PP_type.eq.PAW) then
                                call parallel_task('sum', atoms(j)%PAW%obar_i , parallel, 'node')
                        endif
                enddo 

                deallocate(dense_projector)
                deallocate(coarse_projector)
                deallocate(Rs_dense)
                deallocate(map_dense)
                deallocate(map_coarse)

                do j=1, size(atoms)
                        if(e.eq.atoms(j)%element) then
                                call shared_fence(atoms(j)%RL%proj_shared, parallel)
                                if(elements(e)%PP_type.eq.PAW) &
                                        call shared_fence(atoms(j)%PAW%proj_shared_ort, parallel)
                        endif
                enddo

                if(any(atoms(:)%update_me_force)) then
                        do j=1, size(atoms)
                                if(e.eq.atoms(j)%element) then
                                        call shared_fence(atoms(j)%RL%deriv_proj_shared, parallel)
                                        if(elements(e)%PP_type.eq.PAW) &
                                                call shared_fence(atoms(j)%PAW%deriv_proj_shared_ort, parallel)

                                endif
                        enddo
                endif

                !=====================================================================================================================
                !Second, We create the map for each processors local space to the point on the atomic sphere
                !Possibly we should sacrifice the memory and have each processor have it's own copy of it's local portion of projector
                !Note that this section is done by every processor while the prior section was parallelized over atoms (by processors on same !node)
                !=====================================================================================================================

                !Calculate the local map, everyone needs the coarse map (since WF are distributed)
                do j=1, size(atoms);if(e.eq.atoms(j)%element) then
                        call map_atom_to_local(atoms(j)%RL%ns,atoms(j)%RL%map,atoms(j)%RL%Rs, &
                                elements(e), atoms(j), grid, elements(e)%RL%coarse_grid, &
                                elements(e)%RL%max_sphere, radius=elements(e)%RL%radius)
                endif;enddo
            else if(elements(e)%PA_type.eq.Reciprocal) then
                do j=1, size(atoms) 
                        if(atoms(j)%element.ne.e) cycle
                        if(.not.atoms(j)%update_me) cycle
                        if(elements(e)%PP_type.eq.PAW) atoms(j)%PAW%obar_i(:)=elements(e)%PAW%o_i(:)
                enddo
            endif

            !Fine-Grid Map needed for both recip and real projectors since it defines fgrtab
            if(elements(e)%PP_type.eq.PAW) then
                do j=1, size(atoms);if(e.eq.atoms(j)%element) then
                        !Only one proce per space needs the fine grid, distributed otherwise over atoms
                        if(.not.associated(atoms(j)%PAW%fgrtab)) cycle

                        call map_atom_to_local(atoms(j)%PAW%ns_fine, atoms(j)%PAW%map_fine,atoms(j)%PAW%Rs_fine, &
                        elements(e), atoms(j), fine_grid, elements(e)%PAW%fine_grid, &
                        elements(e)%PAW%max_sphere_fine, radius=elements(e)%PAW%radius)

                        if(allocated(atoms(j)%PAW%fgrtab%gylm)) deallocate(atoms(j)%PAW%fgrtab%gylm)
                        if(allocated(atoms(j)%PAW%fgrtab%gylmgr)) deallocate(atoms(j)%PAW%fgrtab%gylmgr)
                        if(allocated(atoms(j)%PAW%fgrtab%gylmgr2)) deallocate(atoms(j)%PAW%fgrtab%gylmgr2)
                        if(allocated(atoms(j)%PAW%fgrtab%ifftsph)) deallocate(atoms(j)%PAW%fgrtab%ifftsph)
                        if(allocated(atoms(j)%PAW%fgrtab%rfgd)) deallocate(atoms(j)%PAW%fgrtab%rfgd)


                        allocate(atoms(j)%PAW%fgrtab%gylm(atoms(j)%PAW%ns_fine,atoms(j)%PAW%an%lm_size))
                        allocate(atoms(j)%PAW%fgrtab%gylmgr(3,atoms(j)%PAW%ns_fine,atoms(j)%PAW%an%lm_size))
                        allocate(atoms(j)%PAW%fgrtab%gylmgr2(6,atoms(j)%PAW%ns_fine,atoms(j)%PAW%an%lm_size))
                        allocate(atoms(j)%PAW%fgrtab%ifftsph(atoms(j)%PAW%ns_fine))
                        allocate(atoms(j)%PAW%fgrtab%rfgd(3,atoms(j)%PAW%ns_fine))

                        atoms(j)%PAW%fgrtab%nfgd=atoms(j)%PAW%ns_fine
                        atoms(j)%PAW%fgrtab%ifftsph(:)=atoms(j)%PAW%map_fine(5,:atoms(j)%PAW%ns_fine)
                        atoms(j)%PAW%fgrtab%rfgd=atoms(j)%PAW%Rs_fine(:,:atoms(j)%PAW%ns_fine)


                        call pawgylm(gylm=atoms(j)%PAW%fgrtab%gylm, &
                        gylmgr=atoms(j)%PAW%fgrtab%gylmgr, &
                        gylmgr2=atoms(j)%PAW%fgrtab%gylmgr2, lm_size=atoms(j)%PAW%an%lm_size, &
                        nfgd=atoms(j)%PAW%fgrtab%nfgd, optgr0=1,optgr1=1,optgr2=1,pawtab=elements(e)%PAW%tab, &
                        rfgd=atoms(j)%PAW%fgrtab%rfgd)

                        atoms(j)%PAW%fgrtab%gylm_allocated=1
                        atoms(j)%PAW%fgrtab%gylmgr_allocated=1
                        atoms(j)%PAW%fgrtab%gylmgr2_allocated=1
                        atoms(j)%PAW%fgrtab%rfgd_allocated=1


                endif;enddo
        endif

        enddo    
        call PAW_projector_projector_overlaps(elements, atoms, grid, parallel, all_PAW)

    end subroutine

    subroutine Calculate_Projector_RL(Rvec, element, projector, RO)
        use  element_type, only : element_struct, HGH, PAW
        use HGH_Pseudopotentials, only : Calculate_Projector_HGH_RL
        use PAW_Pseudopotentials, only: Calculate_Projector_PAW_RL, Calculate_Projector_PAW_RO_RL

        use  atom_type, only : atom_struct

        type(element_struct), intent(inout) :: element

        real(dp), intent(in) :: Rvec(:,:)
        real(dp),intent(inout) :: projector(:,:)
        logical, optional :: RO
        logical :: RO_
        integer :: i,k
        RO_=.false.
        if(present(RO)) RO_=RO
        if(element%PP_type.eq.HGH) then
                do k=1, size(projector,2); do i=1,size(projector,1)
                call Calculate_Projector_HGH_RL(Rvec(:,i), element, element%HGH%l(k), element%HGH%m(k), &
                                                element%HGH%i(k), element%HGH%h_or_k(k), projector(i,k))
                enddo; enddo
        else if(element%PP_type.eq.PAW) then
                if(RO_) then
                        call Calculate_Projector_PAW_RO_RL(Rvec, element, projector)
                else
                        call Calculate_Projector_PAW_RL(Rvec, element, projector)
                endif
        else
                print *, ' Unrecognized pseudopotential type:',   element%PP_type
                stop
        endif
    end subroutine

    subroutine Calculate_Projector_GNL(element, grid, kpt)
        use  element_type, only : element_struct, HGH, PAW
        use HGH_Pseudopotentials, only : Calculate_Projector_HGH_GNL
        use PAW_Pseudopotentials, only : Calculate_Projector_PAW_GNL

        use  grids_type, only : grid_struct

        type(element_struct), intent(inout) :: element
        type(grid_struct), intent(in) ::  grid
        integer, intent(in) :: kpt
        integer :: k
        real(dp), allocatable :: Gvec(:,:), projector(:,:)

        if(element%PP_type.eq.HGH) then
                do k=1, element%n_proj
                call Calculate_Projector_HGH_GNL(element, element%HGH%l(k), element%HGH%m(k), &
                                                element%HGH%i(k), element%HGH%h_or_k(k), &
                                                element%GNL(kpt)%projector_G(:,:,:,k), grid%G2, grid%G)
                enddo
        else if(element%PP_type.eq.PAW) then
                allocate(Gvec(3,product(grid%Ng_local)))
                allocate(projector(product(grid%Ng_local),element%n_proj))
                do k=1,3
                        Gvec(k,:)=reshape(grid%G(:,:,:,k),[size(grid%G(:,:,:,1))])
                enddo
                call Calculate_Projector_PAW_GNL(Gvec, element, projector)
                do k=1,element%n_proj
                        element%GNL(kpt)%projector_G(:,:,:,k)= &
                                reshape(projector(:,k),[size(grid%G,1),size(grid%G,2),size(grid%G,3)])
                enddo
                deallocate(Gvec)
                deallocate(projector)
        else
                print *, ' Unrecognized pseudopotential type:',   element%PP_type
                stop
        endif

        do k=1,element%n_proj
                element%GNL(kpt)%projector_G(:,:,:,k)= &
                element%GNL(kpt)%projector_G(:,:,:,k)*grid%cutwf/(product(grid%Box_Length))
                if(grid%gamma) then
                        where(abs(grid%G(:,:,:,1)).gt.tiny(1.0_dp)) 
                                element%GNL(kpt)%projector_G(:,:,:,k)= &
                                element%GNL(kpt)%projector_G(:,:,:,k)*sqrt(2.0_dp) 
                        end where
                endif
        enddo
    end subroutine

    subroutine Calculate_Projector_GNL_ORT(element, grid, kpt)
        use  element_type, only : element_struct, HGH, PAW
        use PAW_Pseudopotentials, only : Calculate_Projector_PAW_GNL, Calculate_Projector_PAW_GNL_ORT

        use  grids_type, only : grid_struct

        type(element_struct), intent(inout) :: element
        type(grid_struct), intent(in) ::  grid
        integer, intent(in) :: kpt
        integer :: k
        real(dp), allocatable :: Gvec(:,:), projector(:,:)

        if(element%PP_type.eq.HGH) then
        else if(element%PP_type.eq.PAW) then
                allocate(Gvec(3,product(grid%Ng_local)))
                allocate(projector(product(grid%Ng_local),element%n_proj))
                do k=1,3
                        Gvec(k,:)=reshape(grid%G(:,:,:,k),[size(grid%G(:,:,:,1))])
                enddo
                call Calculate_Projector_PAW_GNL_ORT(Gvec, element, projector)
                do k=1,element%n_proj
                        element%GNL(kpt)%projector_G_ort(:,:,:,k)= &
                                reshape(projector(:,k),[size(grid%G,1),size(grid%G,2),size(grid%G,3)])
                enddo
                deallocate(Gvec)
                deallocate(projector)
                
                do k=1,element%n_proj
                        element%GNL(kpt)%projector_G_ort(:,:,:,k)= &
                        element%GNL(kpt)%projector_G_ort(:,:,:,k)*grid%cutwf/(product(grid%Box_Length))
                        if(grid%gamma) then
                                where(abs(grid%G(:,:,:,1)).gt.tiny(1.0_dp)) 
                                        element%GNL(kpt)%projector_G_ort(:,:,:,k)= &
                                        element%GNL(kpt)%projector_G_ort(:,:,:,k)*sqrt(2.0_dp) 
                                end where
                        endif
                enddo
                
        else
                print *, ' Unrecognized pseudopotential type:',   element%PP_type
                stop
        endif

       

    end subroutine

    subroutine Calculate_deriv_Projector_RL(Rvec, element, projector)
        use  element_type, only : element_struct, HGH, PAW
        use HGH_Pseudopotentials, only : Calculate_deriv_Projector_HGH_RL
        use PAW_Pseudopotentials, only: Calculate_deriv_Projector_PAW_RL

        use  atom_type, only : atom_struct

        type(element_struct), intent(inout) :: element

        real(dp), intent(in) :: Rvec(:,:)
        real(dp),intent(inout) :: projector(:,:,:)
        integer :: i,k

        if(element%PP_type.eq.HGH) then
                do k=1, size(projector,2); do i=1,size(projector,1)
                call Calculate_deriv_Projector_HGH_RL(Rvec(:,i), element, element%HGH%l(k), element%HGH%m(k), &
                                                element%HGH%i(k), element%HGH%h_or_k(k), projector(i,k,:))
                enddo; enddo
        else if(element%PP_type.eq.PAW) then
                call Calculate_deriv_Projector_PAW_RL(Rvec, element, projector)
        else
                print *, ' Unrecognized pseudopotential type:',   element%PP_type
                stop
        endif
    end subroutine

    subroutine Apply_projectors_RL(psi, atoms, elements, grid, PYPsi, dPYPsi, Vpsi, midrVpsi, Spsi, midrSpsi, &
                                 nabla_psi, dipole_psi, RSpsi, SRpsi, skip_PAW, skip_HGH)                     
        use fft, only: recip_to_real, real_to_recip
        use odp_type, only : field_struct
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct, HGH, PAW, Real_local
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        type(field_struct), intent(in) :: psi(:)
        type(field_struct), intent(inout), optional :: Vpsi(:), Spsi(:)
        type(field_struct), intent(inout), optional :: midrVpsi(:,:), midrSpsi(:,:), nabla_psi(:,:), dipole_psi(:,:)
        type(field_struct), intent(inout), optional :: RSpsi(:,:), SRpsi(:,:)

        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(in) :: grid
        type(element_struct), intent(inout), target :: elements(:)
        real(dp), allocatable :: KB(:)
        type(atom_struct), pointer :: atom
        type(element_struct), pointer :: element
        complex(dp), intent(in) :: PYPsi(:,:,:)
        complex(dp), intent(in), optional :: dPYpsi(:,:,:,:)
        complex(dp), allocatable :: SPYpsi(:,:), DijPYpsi(:,:), DijdPYpsi(:,:,:),SdPYpsi(:,:,:), &
                                        vpsi_vec(:),  dvpsi_vec(:,:), factor(:),  &
                                        TijPYpsi(:,:,:),XijPYpsi(:,:,:)!, Mij(:,:,:)
        integer :: ii, ix, iy, iz, i_sphere, i, l, at, i_p, m, s, s2, dir
        logical, optional, intent(in) :: skip_PAW, skip_HGH
        logical :: skip_PAW_, skip_HGH_
        real(dp), allocatable :: proj_packed(:,:)

        skip_PAW_=.false.; skip_HGH_=.false.; 
        if(present(skip_PAW)) skip_PAW_=skip_PAW
        if(present(skip_HGH)) skip_HGH_=skip_HGH
        do s=1, size(psi)
                if(present(Spsi)) Spsi(s)%R=psi(s)%R
        enddo
        if(all(elements(:)%n_proj.eq.0))  return
        if((.not.present(Vpsi).and..not.present(midrVpsi).and..not.present(Spsi) &
         .and..not.present(midrSpsi).and..not.present(nabla_psi).and..not.present(dipole_psi) &
         .and..not.present(RSpsi).and..not.present(SRpsi))) then
                print *, 'Ask for NL apply but dont ask for output, what are you trying to do?'
                stop
        endif

        do s=1, size(psi)
                if(present(Vpsi)) Vpsi(s)%R=0.0_dp
                if(present(midrVpsi)) then
                        midrVpsi(1,s)%R=0.0_dp
                        midrVpsi(2,s)%R=0.0_dp
                        midrVpsi(3,s)%R=0.0_dp
                endif
                if(present(midrSpsi)) then
                        midrSpsi(1,s)%R=0.0_dp
                        midrSpsi(2,s)%R=0.0_dp
                        midrSpsi(3,s)%R=0.0_dp
                endif
                if(present(nabla_psi)) then
                        nabla_psi(1,s)%R=0.0_dp
                        nabla_psi(2,s)%R=0.0_dp
                        nabla_psi(3,s)%R=0.0_dp
                endif
                if(present(dipole_psi)) then
                        dipole_psi(1,s)%R=0.0_dp
                        dipole_psi(2,s)%R=0.0_dp
                        dipole_psi(3,s)%R=0.0_dp
                endif
                if(present(RSpsi)) then
                        RSpsi(1,s)%R=0.0_dp
                        RSpsi(2,s)%R=0.0_dp
                        RSpsi(3,s)%R=0.0_dp
                endif
                if(present(SRpsi)) then
                        SRpsi(1,s)%R=0.0_dp
                        SRpsi(2,s)%R=0.0_dp
                        SRpsi(3,s)%R=0.0_dp
                endif
        enddo

        do at=1, size(atoms)
          atom=>atoms(at)
          element=>elements(atom%element)
          if(element%PA_type.ne.Real_local) cycle
          if(atom%RL%ns.lt.1) cycle
          if(skip_HGH_) then
                if(element%PP_type.eq.HGH) cycle
          endif
          if(skip_PAW_) then
                if(element%PP_type.eq.PAW) cycle
          endif
          allocate(vpsi_vec(atom%RL%ns))
          if(present(midrVpsi).or.present(midrSpsi).or.present(nabla_psi).or.present(dipole_psi) &
                .or.present(SRpsi).or.present(RSpsi)) allocate(dvpsi_vec(3,atom%RL%ns))

          if(present(Spsi).or.present(midrSpsi).or.present(RSpsi)) then
                allocate(SPYpsi(element%n_proj, size(psi)))
                do s=1, size(psi); do i=1, element%n_proj
                        SPYpsi(i,s)=sum(element%PAW%sij(:,i)*PYpsi(:element%n_proj,at,s))
                enddo;enddo
          endif

          if(present(midrSpsi).or.present(SRpsi)) then
                allocate(SdPYpsi(3, element%n_proj, size(psi)))
                do s=1, size(psi); do i=1, element%n_proj; do dir=1,3               
                        SdPYpsi(dir,i,s)=sum(element%PAW%sij(:,i)*dPYpsi(dir,:element%n_proj,at,s))
                enddo; enddo; enddo
          endif

          if(present(nabla_psi)) then
                allocate(TijPYpsi(3, element%n_proj, size(psi)))
                do s=1, size(psi); do i=1, element%n_proj; do dir=1,3               
                        TijPYpsi(dir,i,s)=sum((element%PAW%tab%nabla_ij(dir,i,:) &
                !+i_*(grid%A(dir)+grid%k(dir))*element%PAW%sij(:,i) )*PYpsi(:element%n_proj,at,s))
                 )*PYpsi(:element%n_proj,at,s))
                enddo; enddo; enddo
          endif

          if(present(dipole_psi)) then
                allocate(XijPYpsi(3, element%n_proj, size(psi)))
                do s=1, size(psi); do i=1, element%n_proj; do dir=1,3               
                        XijPYpsi(dir,i,s)=sum((element%PAW%tab%dipole_ij(dir,i,:))*PYpsi(:element%n_proj,at,s))
                enddo; enddo; enddo
          endif

          if(element%PP_type.eq.PAW) then
                if(present(Vpsi).or.present(midrVpsi)) then
                        allocate(DijPYpsi(element%n_proj, size(psi)))
                       ! allocate(Mij(element%n_proj,element%n_proj, size(psi)))

                        !if(any(abs(grid%k).gt.tiny(1.0_dp)).or.any(abs(grid%A).gt.tiny(1.0_dp))) then
                         !       do s=1, size(psi); do i=1, element%n_proj  
                          !              Mij(:,i,s)=atom%PAW%Dij(:,i,s)   &
                           !                      -i_*((grid%A(1)+grid%k(1))*element%PAW%tab%nabla_ij(1,:,i) + &
                            !                         (grid%A(2)+grid%k(2))*element%PAW%tab%nabla_ij(2,:,i) + &
                             !                        (grid%A(3)+grid%k(3))*element%PAW%tab%nabla_ij(3,:,i))+ &
                              !                   sum((grid%A(:)**2+grid%k(:))**2)*element%PAW%sij(:,i)*0.5_dp
                       !         enddo; enddo
                       ! else
                        !        Mij=atom%PAW%Dij +i_*0.0_dp
                       ! endif                        
                        do s=1, size(psi)
                                do i=1, element%n_proj
                                       ! DijPYpsi(i,s)=sum(Mij(i,:,s)*PYpsi(:element%n_proj,at,s))
                                        DijPYpsi(i,s)=sum(atom%PAW%Dij(i,:,s)*PYpsi(:element%n_proj,at,s))
                                enddo
                        enddo
                endif
                if(present(midrVpsi)) then
                        allocate(DijdPYpsi(3, element%n_proj, size(psi)))
                        do s=1, size(psi); do i=1, element%n_proj; do dir=1,3               
                              !  DijdPYpsi(dir,i,s)=sum(Mij(i,:,s)*dPYpsi(dir,:element%n_proj,at,s))
                                DijdPYpsi(dir,i,s)=sum(atom%PAW%Dij(i,:,s)*dPYpsi(dir,:element%n_proj,at,s))
                        enddo; enddo; enddo
                endif
                !if(present(Vpsi).or.present(midrVpsi))  deallocate(Mij)

          endif

          allocate(proj_packed(element%n_proj,atom%RL%ns))
          do ii= 1, atom%RL%ns
                i_sphere=atom%RL%map(4,ii)
                proj_packed(:,ii)=atoms(at)%RL%projector(i_sphere,:)
          enddo

          if(element%PP_type.eq.HGH) then
                if((any(abs(grid%k).gt.tiny(1.0_dp)).or.any(abs(grid%A).gt.tiny(1.0_dp)))) then
                        allocate(factor(atom%RL%ns))
                        factor(:)=exp(-i_*((grid%k(1)+grid%A(1))*atom%RL%Rs(1,:atom%RL%ns) + &
                                        (grid%k(2)+grid%A(2))*atom%RL%Rs(2,:atom%RL%ns) + &
                                        (grid%k(3)+grid%A(3))*atom%RL%Rs(3,:atom%RL%ns)))
                endif
                allocate(KB(element%n_proj))
                KB(:)=0.0_dp
                do s=1, size(psi); do s2=1, size(psi)
                        do i_p=1,element%n_proj
                                i=element%HGH%i(i_p)
                                l=element%HGH%l(i_p)
                                m=element%HGH%m(i_p)
                                if(element%HGH%h_or_k(i_p).eq.'h') then
                                        if(s.eq.s2) then
                                                KB(i_p)=element%HGH%EKB(i,l+1)
                                        else
                                                cycle
                                        endif
                                else
                                        if(s.eq.s2) then
                                                KB(i_p)=element%HGH%KKB(i,l+1)*(s*1.0_dp-1.5_dp)*m
                                        else
                                                KB(i_p)=element%HGH%KKB(i,l+1)*0.5_dp* &
                                                        sqrt(1.0_dp*(l+m+1))*sqrt(1.0_dp*(l-m))
                                        endif
                                endif
                        enddo
                if(present(Vpsi)) then
                        if(allocated(factor)) then
                                do ii= 1, atom%RL%ns
                                        vpsi_vec(ii)=sum(proj_packed(:,ii)*PYpsi(:,at,s)*KB(:))*factor(ii)
                                enddo
                        else
                                do ii= 1, atom%RL%ns
                                        vpsi_vec(ii)=sum(proj_packed(:,ii)*PYpsi(:,at,s)*KB(:))
                                enddo
                        endif
                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                                Vpsi(s2)%R(ix,iy,iz)=Vpsi(s2)%R(ix,iy,iz)+vpsi_vec(ii)
                        enddo
                endif
                if(present(midrVpsi))  then
                        if(allocated(factor)) then
                                do ii= 1, atom%RL%ns; do dir=1,3
                                        dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)*factor(ii)* &
                                                (dPYpsi(dir,:,at,s) - atom%RL%Rs(dir,ii)*PYpsi(:,at,s))*KB(:))
                                enddo; enddo
                        else
                                do ii= 1, atom%RL%ns; do dir=1,3
                                        dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)* &
                                                (dPYpsi(dir,:,at,s)- atom%RL%Rs(dir,ii)*PYpsi(:,at,s))*KB(:))
                                enddo; enddo
                        endif
                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                                do dir=1,3
                                        midrVpsi(dir,s2)%R(ix,iy,iz)=midrVpsi(dir,s2)%R(ix,iy,iz)+&
                                        dvpsi_vec(dir,ii)
                                enddo
                        enddo
                endif
                enddo; enddo
                if(allocated(factor)) deallocate(factor)
                deallocate(KB)
          else if(element%PP_type.eq.PAW) then
                if((any(abs(grid%k).gt.tiny(1.0_dp)).or.any(abs(grid%A).gt.tiny(1.0_dp)))) then
                        allocate(factor(atom%RL%ns))
                        factor(:)=exp(-i_*((grid%k(1)+grid%A(1))*atom%RL%Rs(1,:atom%RL%ns) + &
                                        (grid%k(2)+grid%A(2))*atom%RL%Rs(2,:atom%RL%ns) + &
                                        (grid%k(3)+grid%A(3))*atom%RL%Rs(3,:atom%RL%ns)))
                endif
                do s=1, size(psi)
                if(present(Vpsi)) then
                        if(allocated(factor)) then
                                do ii= 1, atom%RL%ns
                                        vpsi_vec(ii)=sum(proj_packed(:,ii)*DijPYpsi(:,s))*factor(ii)
                                enddo
                        else
                                do ii= 1, atom%RL%ns
                                        vpsi_vec(ii)=sum(proj_packed(:,ii)*DijPYpsi(:,s))
                                enddo
                        endif

                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                                Vpsi(s)%R(ix,iy,iz)=Vpsi(s)%R(ix,iy,iz)+vpsi_vec(ii)
                        enddo
                endif
                if(present(midrVpsi))  then
                        if(allocated(factor)) then
                                do ii= 1, atom%RL%ns; do dir=1,3
                                        dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)* &
                                        (DijdPYpsi(dir,:,s) - atom%RL%Rs(dir,ii)*DijPYpsi(:,s)))*factor(ii)
                                enddo; enddo
                        else
                                do ii= 1, atom%RL%ns; do dir=1,3
                                        dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)* &
                                                (DijdPYpsi(dir,:,s) - atom%RL%Rs(dir,ii)*DijPYpsi(:,s)))
                                enddo; enddo
                        endif

                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                                do dir=1,3
                                        midrVpsi(dir,s)%R(ix,iy,iz)= midrVpsi(dir,s)%R(ix,iy,iz)+&
                                        dvpsi_vec(dir,ii)
                                enddo
                        enddo
                endif
                if(present(Spsi))  then
                        if(allocated(factor)) then
                                do ii= 1, atom%RL%ns
                                        vpsi_vec(ii)=sum(proj_packed(:,ii)*SPYpsi(:,s))*factor(ii)
                                enddo
                        else
                                do ii= 1, atom%RL%ns
                                        vpsi_vec(ii)=sum(proj_packed(:,ii)*SPYpsi(:,s))
                                enddo
                        endif

                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                                Spsi(s)%R(ix,iy,iz)=Spsi(s)%R(ix,iy,iz)+vpsi_vec(ii)
                        enddo
                endif
                if(present(midrSpsi))  then
                        if(allocated(factor)) then
                                do ii= 1, atom%RL%ns; do dir=1,3
                                        dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)* (&
                                                SdPYpsi(dir,:,s) &
                                                 - atom%RL%Rs(dir,ii)*SPYpsi(:,s)))*factor(ii)
                                enddo; enddo
                        else      
                                do ii= 1, atom%RL%ns; do dir=1,3
                                        dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)* (&
                                                SdPYpsi(dir,:,s) &
                                                - atom%RL%Rs(dir,ii)*SPYpsi(:,s)))
                                enddo; enddo
                        endif
                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                                do dir=1,3
                                        midrSpsi(dir,s)%R(ix,iy,iz)=midrSpsi(dir,s)%R(ix,iy,iz) + &
                                        dvpsi_vec(dir,ii)
                                enddo
                        enddo
                endif
                if(present(RSpsi))  then
                        if(allocated(factor)) then
                                do ii= 1, atom%RL%ns; do dir=1,3
                                        dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)* (&
                                                atom%RL%Rs(dir,ii)*SPYpsi(:,s)))*factor(ii)
                                enddo; enddo
                        else      
                                do ii= 1, atom%RL%ns; do dir=1,3
                                        dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)* (&
                                                atom%RL%Rs(dir,ii)*SPYpsi(:,s)))
                                enddo; enddo
                        endif
                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                                do dir=1,3
                                        RSpsi(dir,s)%R(ix,iy,iz)=RSpsi(dir,s)%R(ix,iy,iz) + &
                                        dvpsi_vec(dir,ii)
                                enddo
                        enddo
                endif
                if(present(SRpsi))  then
                        if(allocated(factor)) then
                                do ii= 1, atom%RL%ns; do dir=1,3
                                        dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)* (&
                                                SdPYpsi(dir,:,s)))*factor(ii)
                                enddo; enddo
                        else      
                                do ii= 1, atom%RL%ns; do dir=1,3
                                        dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)* (&
                                                SdPYpsi(dir,:,s)))
                                enddo; enddo
                        endif
                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                                do dir=1,3
                                        SRpsi(dir,s)%R(ix,iy,iz)=SRpsi(dir,s)%R(ix,iy,iz) + &
                                        dvpsi_vec(dir,ii)
                                enddo
                        enddo
                endif
                if(present(nabla_psi))  then
                        if(allocated(factor)) then
                                do ii= 1, atom%RL%ns; do dir=1,3
                                        dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)*TijPYpsi(dir,:,s))*factor(ii)
                                enddo; enddo
                        else
                                do ii= 1, atom%RL%ns; do dir=1,3
                                        dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)*TijPYpsi(dir,:,s))  
                                enddo; enddo
                        endif

                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                                do dir=1,3
                                        nabla_psi(dir,s)%R(ix,iy,iz)=nabla_psi(dir,s)%R(ix,iy,iz) + &
                                        dvpsi_vec(dir,ii)
                                enddo
                        enddo
                endif
                if(present(dipole_psi))  then
                        if(allocated(factor)) then
                                do ii= 1, atom%RL%ns; do dir=1,3
                                        dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)*XijPYpsi(dir,:,s))*factor(ii)
                                enddo; enddo
                        else
                                do ii= 1, atom%RL%ns; do dir=1,3
                                        dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)*XijPYpsi(dir,:,s))  
                                enddo; enddo
                        endif

                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                                do dir=1,3
                                        dipole_psi(dir,s)%R(ix,iy,iz)=dipole_psi(dir,s)%R(ix,iy,iz) + &
                                        dvpsi_vec(dir,ii)
                                enddo
                        enddo
                endif
                enddo
                if(allocated(factor)) deallocate(factor)
          else
                print *, 'Unrecognized PP type', element%PP_type
          endif
          if(element%PP_type.eq.PAW) then
                if(present(Vpsi).or.present(midrVpsi)) deallocate(DijPYpsi)     
                if(present(midrVpsi)) deallocate(DijdPYpsi) 
                if(present(nabla_psi))   deallocate(TijPYpsi) 
                if(present(dipole_psi))   deallocate(XijPYpsi) 
                if(present(Spsi).or.present(midrSpsi).or.present(RSpsi))    deallocate(SPYpsi) 
                if(present(midrSpsi).or.present(SRpsi)) deallocate(SdPYpsi)
          endif
          deallocate(vpsi_vec)
          if(present(midrVpsi).or.present(midrSpsi) &
         .or.present(nabla_psi).or.present(dipole_psi).or.present(RSpsi).or.present(SRpsi)) deallocate(dvpsi_vec)
          deallocate(proj_packed)        
        enddo
        
end subroutine

subroutine Apply_projectors_Recip(psi, atoms, elements, grid, PYPsi, dPYPsi, Vpsi, midrVpsi, Spsi, midrSpsi, &
        nabla_psi, dipole_psi, RSpsi, SRpsi, skip_PAW, skip_HGH)                     
        use fft, only: recip_to_real, real_to_recip
        use odp_type, only : field_struct
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct, HGH, PAW, Reciprocal
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        use grids_mod, only : allocate_local_fields_G
        type(field_struct), intent(in) :: psi(:)
        type(field_struct), intent(inout), optional :: Vpsi(:), Spsi(:)
        type(field_struct), intent(inout), optional :: midrVpsi(:,:), midrSpsi(:,:), nabla_psi(:,:), dipole_psi(:,:)
        type(field_struct), intent(inout), optional :: RSpsi(:,:), SRpsi(:,:)

        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(in) :: grid
        type(element_struct), intent(inout), target :: elements(:)
        real(dp), allocatable :: KB(:)
        type(atom_struct), pointer :: atom
        type(element_struct), pointer :: element
        complex(dp), intent(in) :: PYPsi(:,:,:)
        complex(dp), intent(in), optional :: dPYpsi(:,:,:,:)
        complex(dp), allocatable :: SPYpsi(:,:), DijPYpsi(:,:), DijdPYpsi(:,:,:),SdPYpsi(:,:,:), &
                        factor(:,:,:), TijPYpsi(:,:,:),XijPYpsi(:,:,:)!, Mij(:,:,:)
        integer :: i, l, at, i_p, m, s, s2, dir, j, il, jl
        logical, optional, intent(in) :: skip_PAW, skip_HGH
        logical :: skip_PAW_, skip_HGH_

        skip_PAW_=.false.; skip_HGH_=.false.; 
        if(present(skip_PAW)) skip_PAW_=skip_PAW
        if(present(skip_HGH)) skip_HGH_=skip_HGH

        if(all(elements(:)%n_proj.eq.0))  return
        if((.not.present(Vpsi).and..not.present(midrVpsi).and..not.present(Spsi) &
        .and..not.present(midrSpsi).and..not.present(nabla_psi).and..not.present(dipole_psi) &
        .and..not.present(RSpsi).and..not.present(SRpsi))) then
        print *, 'Ask for NL apply but dont ask for output, what are you trying to do?'
        stop
        endif

        do s=1, size(psi)
        if(present(Vpsi)) Vpsi(s)%G=0.0_dp
        if(present(Spsi)) Spsi(s)%G=0.0_dp
        if(present(midrVpsi)) then
        midrVpsi(1,s)%G=0.0_dp
        midrVpsi(2,s)%G=0.0_dp
        midrVpsi(3,s)%G=0.0_dp
        endif
        if(present(midrSpsi)) then
        midrSpsi(1,s)%G=0.0_dp
        midrSpsi(2,s)%G=0.0_dp
        midrSpsi(3,s)%G=0.0_dp
        endif
        if(present(nabla_psi)) then
        nabla_psi(1,s)%G=0.0_dp
        nabla_psi(2,s)%G=0.0_dp
        nabla_psi(3,s)%G=0.0_dp
        endif
        if(present(dipole_psi)) then
        dipole_psi(1,s)%G=0.0_dp
        dipole_psi(2,s)%G=0.0_dp
        dipole_psi(3,s)%G=0.0_dp
        endif
        if(present(RSpsi)) then
        RSpsi(1,s)%G=0.0_dp
        RSpsi(2,s)%G=0.0_dp
        RSpsi(3,s)%G=0.0_dp
        endif
        if(present(SRpsi)) then
        SRpsi(1,s)%G=0.0_dp
        SRpsi(2,s)%G=0.0_dp
        SRpsi(3,s)%G=0.0_dp
        endif
        enddo

        do at=1, size(atoms)
        atom=>atoms(at)
        element=>elements(atom%element)
        if(element%PA_type.ne.Reciprocal) cycle
        if(skip_HGH_) then
        if(element%PP_type.eq.HGH) cycle
        endif
        if(skip_PAW_) then
        if(element%PP_type.eq.PAW) cycle
        endif
        if(present(Spsi).or.present(midrSpsi).or.present(RSpsi)) then
                allocate(SPYpsi(element%n_proj, size(psi)))
                do s=1, size(psi); do i=1, element%n_proj
                        SPYpsi(i,s)=sum(element%PAW%sij(:,i)*PYpsi(:element%n_proj,at,s))
                enddo;enddo
          endif

          if(present(midrSpsi).or.present(SRpsi)) then
                allocate(SdPYpsi(3, element%n_proj, size(psi)))
                do s=1, size(psi); do i=1, element%n_proj; do dir=1,3               
                        SdPYpsi(dir,i,s)=sum(element%PAW%sij(:,i)*dPYpsi(dir,:element%n_proj,at,s))
                enddo; enddo; enddo
          endif

          if(present(nabla_psi)) then
                allocate(TijPYpsi(3, element%n_proj, size(psi)))
                do s=1, size(psi); do i=1, element%n_proj; do dir=1,3               
                        TijPYpsi(dir,i,s)=sum((element%PAW%tab%nabla_ij(dir,i,:) &
                !+i_*(grid%A(dir)+grid%k(dir))*element%PAW%sij(:,i) )*PYpsi(:element%n_proj,at,s))
                 )*PYpsi(:element%n_proj,at,s))
                enddo; enddo; enddo
          endif

          if(present(dipole_psi)) then
                allocate(XijPYpsi(3, element%n_proj, size(psi)))
                do s=1, size(psi); do i=1, element%n_proj; do dir=1,3               
                        XijPYpsi(dir,i,s)=sum((element%PAW%tab%dipole_ij(dir,i,:))*PYpsi(:element%n_proj,at,s))
                enddo; enddo; enddo
          endif

          if(element%PP_type.eq.PAW) then
                if(present(Vpsi).or.present(midrVpsi)) then
                        allocate(DijPYpsi(element%n_proj, size(psi)))
                        do s=1, size(psi)
                                do i=1, element%n_proj
                                        DijPYpsi(i,s)=sum(atom%PAW%Dij(i,:,s)*PYpsi(:element%n_proj,at,s))
                                enddo
                        enddo
                endif
                if(present(midrVpsi)) then
                        allocate(DijdPYpsi(3, element%n_proj, size(psi)))
                        do s=1, size(psi); do i=1, element%n_proj; do dir=1,3               
                                DijdPYpsi(dir,i,s)=sum(atom%PAW%Dij(i,:,s)*dPYpsi(dir,:element%n_proj,at,s))
                        enddo; enddo; enddo
                endif
          endif
        if(element%PP_type.eq.HGH) then
                call allocate_local_fields_G(factor,grid)
                factor(:,:,:)=exp(-i_*((grid%G(:,:,:,1)-grid%A(1))*atom%R(1) +&
                                       (grid%G(:,:,:,2)-grid%A(2))*atom%R(2) +&
                                       (grid%G(:,:,:,3)-grid%A(3))*atom%R(3)  ))
                allocate(KB(element%n_proj))
                KB(:)=0.0_dp
                do s=1, size(psi); do s2=1, size(psi)
                do i_p=1,element%n_proj
                i=element%HGH%i(i_p)
                l=element%HGH%l(i_p)
                m=element%HGH%m(i_p)
                if(element%HGH%h_or_k(i_p).eq.'h') then
                        if(s.eq.s2) then
                                KB(i_p)=element%HGH%EKB(i,l+1)
                        else
                                cycle
                        endif
                else
                        if(s.eq.s2) then
                                KB(i_p)=element%HGH%KKB(i,l+1)*(s*1.0_dp-1.5_dp)*m
                        else
                                KB(i_p)=element%HGH%KKB(i,l+1)*0.5_dp* &
                                        sqrt(1.0_dp*(l+m+1))*sqrt(1.0_dp*(l-m))
                        endif
                endif
                enddo
                if(present(Vpsi)) then
                        do i_p=1,element%n_proj
                                Vpsi(s2)%G(:,:,:)=Vpsi(s2)%G(:,:,:)+ &
                                factor(:,:,:)*element%GNL(psi(s)%grid)%projector_G(:,:,:,i_p)*KB(i_p)*PYPsi(i_p,at,s)
                        enddo
                endif
                if(present(midrVpsi))  then
                ! midrVpsi(dir,s2)%R(ix,iy,iz)=midrVpsi(dir,s2)%R(ix,iy,iz)+&
                !        dvpsi_vec(dir,ii)
                        print *, 'Need midrVpsi for recip'; flush(6) 
                        stop
                endif
                enddo; enddo
                deallocate(factor)
                deallocate(KB)
        else if(element%PP_type.eq.PAW) then
                call allocate_local_fields_G(factor,grid)
                factor(:,:,:)=exp(-i_*((grid%G(:,:,:,1)-grid%A(1))*atom%R(1) +&
                                       (grid%G(:,:,:,2)-grid%A(2))*atom%R(2) +&
                                       (grid%G(:,:,:,3)-grid%A(3))*atom%R(3)  ))

                do s=1, size(psi)
                if(present(Vpsi)) then
                        do i_p=1,element%n_proj
                                il=element%PAW%tab%indlmn(1,i_p)
                                Vpsi(s)%G(:,:,:)=Vpsi(s)%G(:,:,:)+ (-i_)**il*&
                                factor(:,:,:)*element%GNL(psi(s)%grid)%projector_G(:,:,:,i_p)*DijPYpsi(i_p,s)
                        enddo
                endif
                if(present(midrVpsi))  then
                        ! midrVpsi(dir,s)%R(ix,iy,iz)=midrVpsi(dir,s)%R(ix,iy,iz)+&
                !        dvpsi_vec(dir,ii)
                        print *, 'Need midrVpsi for recip'; flush(6) 
                        stop
                endif
                if(present(Spsi))  then
                        do i_p=1,element%n_proj
                                il=element%PAW%tab%indlmn(1,i_p)
                                Spsi(s)%G(:,:,:)=Spsi(s)%G(:,:,:)+(-i_)**il*&
                                factor(:,:,:)*element%GNL(psi(s)%grid)%projector_G(:,:,:,i_p)*SPYpsi(i_p,s)
                        enddo
                endif
                if(present(midrSpsi))  then
                        print *, 'Need midrSpsi for recip'; flush(6) 
                        stop
                endif
                if(present(RSpsi))  then
                        print *, 'Need RSpsi for recip'; flush(6) 
                        stop
                endif
                if(present(SRpsi))  then
                        print *, 'Need SRpsi for recip'; flush(6) 
                        stop
                endif
                if(present(nabla_psi))  then
                        do i_p=1,element%n_proj;do dir=1,3
                                il=element%PAW%tab%indlmn(1,i_p)
                                nabla_psi(dir,s)%G(:,:,:)=nabla_psi(dir,s)%G(:,:,:)+ (-i_)**il*&
                                factor(:,:,:)*element%GNL(psi(s)%grid)%projector_G(:,:,:,i_p)*TijPYpsi(dir,i_p,s)
                        enddo; enddo
                endif

                if(present(dipole_psi))  then
                        do i_p=1,element%n_proj;do dir=1,3
                                il=element%PAW%tab%indlmn(1,i_p)
                                dipole_psi(dir,s)%G(:,:,:)=dipole_psi(dir,s)%G(:,:,:)+ (-i_)**il*&
                                factor(:,:,:)*element%GNL(psi(s)%grid)%projector_G(:,:,:,i_p)*XijPYpsi(dir,i_p,s)
                        enddo; enddo
                endif
                enddo
                deallocate(factor)
        else
        print *, 'Unrecognized PP type', element%PP_type
        endif
        if(element%PP_type.eq.PAW) then
        if(present(Vpsi).or.present(midrVpsi)) deallocate(DijPYpsi)     
        if(present(midrVpsi)) deallocate(DijdPYpsi) 
        if(present(nabla_psi))   deallocate(TijPYpsi) 
        if(present(dipole_psi))   deallocate(XijPYpsi) 
        if(present(Spsi).or.present(midrSpsi).or.present(RSpsi))    deallocate(SPYpsi) 
        if(present(midrSpsi).or.present(SRpsi)) deallocate(SdPYpsi)
        endif
        enddo

end subroutine

subroutine Apply_projectors_RL_gamma(psi, atoms, elements, grid, PYPsi_, dPYPsi_, Vpsi, midrVpsi, Spsi, midrSpsi, &
        nabla_psi, dipole_psi, RSpsi, SRpsi, skip_PAW, skip_HGH)                     
    use fft, only: recip_to_real, real_to_recip
    use odp_type, only : field_struct
    use  atom_type, only : atom_struct
    use  element_type, only : element_struct, HGH, PAW, Real_local
    use  grids_type, only : grid_struct
    use  parallel_type, only : parallel_struct
    type(field_struct), intent(in) :: psi(:)
    type(field_struct), intent(inout), optional :: Vpsi(:), Spsi(:)
    type(field_struct), intent(inout), optional :: midrVpsi(:,:), midrSpsi(:,:), nabla_psi(:,:), dipole_psi(:,:)
    type(field_struct), intent(inout), optional :: RSpsi(:,:), SRpsi(:,:)
    
    type(atom_struct), intent(inout), target :: atoms(:)
    type(grid_struct), intent(in) :: grid
    type(element_struct), intent(inout), target :: elements(:)
    real(dp), allocatable :: KB(:)
    type(atom_struct), pointer :: atom
    type(element_struct), pointer :: element
    complex(dp), intent(in) :: PYPsi_(:,:,:)
    complex(dp), intent(in), optional :: dPYpsi_(:,:,:,:)
    real(dp), allocatable :: PYPsi(:,:,:)
    real(dp), allocatable :: dPYpsi(:,:,:,:)
    
    real(dp), allocatable :: SPYpsi(:,:), DijPYpsi(:,:), DijdPYpsi(:,:,:),SdPYpsi(:,:,:), &
               vpsi_vec(:),  dvpsi_vec(:,:), factor(:),  &
               TijPYpsi(:,:,:),XijPYpsi(:,:,:), Mij(:,:,:)
    integer :: ii, ix, iy, iz, i_sphere, i, l, at, i_p, m, s, s2, dir
    logical, optional, intent(in) :: skip_PAW, skip_HGH
    logical :: skip_PAW_, skip_HGH_
    real(dp), allocatable :: proj_packed(:,:)
    
    !Cast to reals
    allocate(PYPsi(size(PYPsi_,1),size(PYPsi_,2),size(PYPsi_,3)))
    if(present(dPYpsi_)) allocate(dPYPsi(size(dPYPsi_,1),size(dPYPsi_,2),size(dPYPsi_,3),size(dPYPsi_,4)))
    PYPsi=real(PYPsi_)
    if(present(dPYpsi_)) dPYpsi=real(dPYpsi_)
    
    skip_PAW_=.false.; skip_HGH_=.false.; 
    if(present(skip_PAW)) skip_PAW_=skip_PAW
    if(present(skip_HGH)) skip_HGH_=skip_HGH
   
    do s=1, size(psi)
        if(present(Spsi)) Spsi(s)%R=psi(s)%R
    enddo

    if(all(elements(:)%n_proj.eq.0))  return
    if((.not.present(Vpsi).and..not.present(midrVpsi).and..not.present(Spsi) &
    .and..not.present(midrSpsi).and..not.present(nabla_psi).and..not.present(dipole_psi) &
    .and..not.present(RSpsi).and..not.present(SRpsi))) then
    print *, 'Ask for NL apply but dont ask for output, what are you trying to do?'
    stop
    endif
    
    do s=1, size(psi)
    if(present(Vpsi)) Vpsi(s)%R=0.0_dp
    if(present(midrVpsi)) then
    midrVpsi(1,s)%R=0.0_dp
    midrVpsi(2,s)%R=0.0_dp
    midrVpsi(3,s)%R=0.0_dp
    endif
    if(present(midrSpsi)) then
    midrSpsi(1,s)%R=0.0_dp
    midrSpsi(2,s)%R=0.0_dp
    midrSpsi(3,s)%R=0.0_dp
    endif
    if(present(nabla_psi)) then
    nabla_psi(1,s)%R=0.0_dp
    nabla_psi(2,s)%R=0.0_dp
    nabla_psi(3,s)%R=0.0_dp
    endif
    if(present(dipole_psi)) then
    dipole_psi(1,s)%R=0.0_dp
    dipole_psi(2,s)%R=0.0_dp
    dipole_psi(3,s)%R=0.0_dp
    endif
    if(present(RSpsi)) then
    RSpsi(1,s)%R=0.0_dp
    RSpsi(2,s)%R=0.0_dp
    RSpsi(3,s)%R=0.0_dp
    endif
    if(present(SRpsi)) then
    SRpsi(1,s)%R=0.0_dp
    SRpsi(2,s)%R=0.0_dp
    SRpsi(3,s)%R=0.0_dp
    endif
    enddo
    
    do at=1, size(atoms)
    atom=>atoms(at)
    element=>elements(atom%element)
    if(element%PA_type.ne.Real_local) cycle
    if(atom%RL%ns.lt.1) cycle
    if(skip_HGH_) then
    if(element%PP_type.eq.HGH) cycle
    endif
    if(skip_PAW_) then
    if(element%PP_type.eq.PAW) cycle
    endif
    allocate(vpsi_vec(atom%RL%ns))
    if(present(midrVpsi).or.present(midrSpsi).or.present(nabla_psi).or.present(dipole_psi) &
    .or.present(SRpsi).or.present(RSpsi)) allocate(dvpsi_vec(3,atom%RL%ns))
    
    if(present(Spsi).or.present(midrSpsi).or.present(RSpsi)) then
    allocate(SPYpsi(element%n_proj, size(psi)))
    do s=1, size(psi); do i=1, element%n_proj
    SPYpsi(i,s)=sum(element%PAW%sij(:,i)*PYpsi(:element%n_proj,at,s))
    enddo;enddo
    endif
    
    if(present(midrSpsi).or.present(SRpsi)) then
    allocate(SdPYpsi(3, element%n_proj, size(psi)))
    do s=1, size(psi); do i=1, element%n_proj; do dir=1,3               
    SdPYpsi(dir,i,s)=sum(element%PAW%sij(:,i)*dPYpsi(dir,:element%n_proj,at,s))
    enddo; enddo; enddo
    endif
    
    if(present(nabla_psi)) then
    allocate(TijPYpsi(3, element%n_proj, size(psi)))
    do s=1, size(psi); do i=1, element%n_proj; do dir=1,3               
    TijPYpsi(dir,i,s)=sum(element%PAW%tab%nabla_ij(dir,i,:)*PYpsi(:element%n_proj,at,s))
    enddo; enddo; enddo
    endif
    
    if(present(dipole_psi)) then
    allocate(XijPYpsi(3, element%n_proj, size(psi)))
    do s=1, size(psi); do i=1, element%n_proj; do dir=1,3               
    XijPYpsi(dir,i,s)=sum((element%PAW%tab%dipole_ij(dir,i,:))*PYpsi(:element%n_proj,at,s))
    enddo; enddo; enddo
    endif
    
    if(element%PP_type.eq.PAW) then
    if(present(Vpsi).or.present(midrVpsi)) then
    allocate(DijPYpsi(element%n_proj, size(psi)))
    allocate(Mij(element%n_proj,element%n_proj, size(psi)))
    Mij=atom%PAW%Dij
    do s=1, size(psi)
       do i=1, element%n_proj
               DijPYpsi(i,s)=sum(Mij(i,:,s)*PYpsi(:element%n_proj,at,s))
       enddo
    enddo
    endif
    if(present(midrVpsi)) then
    allocate(DijdPYpsi(3, element%n_proj, size(psi)))
    do s=1, size(psi); do i=1, element%n_proj; do dir=1,3               
       DijdPYpsi(dir,i,s)=sum(Mij(i,:,s)*dPYpsi(dir,:element%n_proj,at,s))
    enddo; enddo; enddo
    endif
    if(present(Vpsi).or.present(midrVpsi))  deallocate(Mij)
    
    endif
    
    allocate(proj_packed(element%n_proj,atom%RL%ns))
    do ii= 1, atom%RL%ns
    i_sphere=atom%RL%map(4,ii)
    proj_packed(:,ii)=atoms(at)%RL%projector(i_sphere,:)
    enddo
    
    if(element%PP_type.eq.HGH) then
    if((any(abs(grid%k).gt.tiny(1.0_dp)).or.any(abs(grid%A).gt.tiny(1.0_dp)))) then
        print *, "Error, into a Apply_Projectors for real WF (Gamma) but have finite k or A"
        stop
    endif
    allocate(KB(element%n_proj))
    KB(:)=0.0_dp
    do s=1, size(psi); do s2=1, size(psi)
    do i_p=1,element%n_proj
       i=element%HGH%i(i_p)
       l=element%HGH%l(i_p)
       m=element%HGH%m(i_p)
       if(element%HGH%h_or_k(i_p).eq.'h') then
               if(s.eq.s2) then
                       KB(i_p)=element%HGH%EKB(i,l+1)
               else
                       cycle
               endif
       else
               if(s.eq.s2) then
                       KB(i_p)=element%HGH%KKB(i,l+1)*(s*1.0_dp-1.5_dp)*m
               else
                       KB(i_p)=element%HGH%KKB(i,l+1)*0.5_dp* &
                               sqrt(1.0_dp*(l+m+1))*sqrt(1.0_dp*(l-m))
               endif
       endif
    enddo
    if(present(Vpsi)) then
       do ii= 1, atom%RL%ns
               vpsi_vec(ii)=sum(proj_packed(:,ii)*PYpsi(:,at,s)*KB(:))
       enddo
        do ii= 1, atom%RL%ns
        ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
        Vpsi(s2)%R(ix,iy,iz)=Vpsi(s2)%R(ix,iy,iz)+vpsi_vec(ii)
        enddo
    endif
    if(present(midrVpsi))  then
       do ii= 1, atom%RL%ns; do dir=1,3
               dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)* &
                       (dPYpsi(dir,:,at,s)- atom%RL%Rs(dir,ii)*PYpsi(:,at,s))*KB(:))
       enddo; enddo
        do ii= 1, atom%RL%ns
        ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
        do dir=1,3
                midrVpsi(dir,s2)%R(ix,iy,iz)=midrVpsi(dir,s2)%R(ix,iy,iz)+&
                dvpsi_vec(dir,ii)
        enddo
        enddo
    endif
    enddo; enddo
    deallocate(KB)
    else if(element%PP_type.eq.PAW) then
    if((any(abs(grid%k).gt.tiny(1.0_dp)).or.any(abs(grid%A).gt.tiny(1.0_dp)))) then
        print *, "Error, into a Apply_Projectors for real WF (Gamma) but have finite k or A"
        stop
    endif
    do s=1, size(psi)
    if(present(Vpsi)) then
       do ii= 1, atom%RL%ns
               vpsi_vec(ii)=sum(proj_packed(:,ii)*DijPYpsi(:,s))
       enddo
        do ii= 1, atom%RL%ns
        ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
        Vpsi(s)%R(ix,iy,iz)=Vpsi(s)%R(ix,iy,iz)+vpsi_vec(ii)
        enddo
    endif
    if(present(midrVpsi))  then
       do ii= 1, atom%RL%ns; do dir=1,3
               dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)* &
                       (DijdPYpsi(dir,:,s) - atom%RL%Rs(dir,ii)*DijPYpsi(:,s)))
       enddo; enddo
        do ii= 1, atom%RL%ns
        ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
        do dir=1,3
                midrVpsi(dir,s)%R(ix,iy,iz)= midrVpsi(dir,s)%R(ix,iy,iz)+&
                dvpsi_vec(dir,ii)
        enddo
        enddo
    endif
    if(present(Spsi))  then
       do ii= 1, atom%RL%ns
               vpsi_vec(ii)=sum(proj_packed(:,ii)*SPYpsi(:,s))
       enddo
        do ii= 1, atom%RL%ns
        ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
        Spsi(s)%R(ix,iy,iz)=Spsi(s)%R(ix,iy,iz)+vpsi_vec(ii)
        enddo
    endif
    if(present(midrSpsi))  then
        do ii= 1, atom%RL%ns; do dir=1,3
               dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)* (&
                       SdPYpsi(dir,:,s) &
                       - atom%RL%Rs(dir,ii)*SPYpsi(:,s)))
        enddo; enddo
        do ii= 1, atom%RL%ns
                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                do dir=1,3
                        midrSpsi(dir,s)%R(ix,iy,iz)=midrSpsi(dir,s)%R(ix,iy,iz) + &
                        dvpsi_vec(dir,ii)
                enddo
        enddo
    endif
    if(present(RSpsi))  then    
        do ii= 1, atom%RL%ns; do dir=1,3
               dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)* (&
                       atom%RL%Rs(dir,ii)*SPYpsi(:,s)))
        enddo; enddo
        do ii= 1, atom%RL%ns
                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                do dir=1,3
                        RSpsi(dir,s)%R(ix,iy,iz)=RSpsi(dir,s)%R(ix,iy,iz) + &
                        dvpsi_vec(dir,ii)
                enddo
        enddo
    endif
    if(present(SRpsi))  then    
        do ii= 1, atom%RL%ns; do dir=1,3
               dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)* (&
                       SdPYpsi(dir,:,s)))
        enddo; enddo
        do ii= 1, atom%RL%ns
                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                do dir=1,3
                        SRpsi(dir,s)%R(ix,iy,iz)=SRpsi(dir,s)%R(ix,iy,iz) + &
                        dvpsi_vec(dir,ii)
                enddo
        enddo
    endif
    if(present(nabla_psi))  then
        do ii= 1, atom%RL%ns; do dir=1,3
               dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)*TijPYpsi(dir,:,s))  
        enddo; enddo
    
        do ii= 1, atom%RL%ns
                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                do dir=1,3
                        nabla_psi(dir,s)%R(ix,iy,iz)=nabla_psi(dir,s)%R(ix,iy,iz) + &
                        dvpsi_vec(dir,ii)
                enddo
        enddo
    endif
    if(present(dipole_psi))  then
        do ii= 1, atom%RL%ns; do dir=1,3
               dvpsi_vec(dir,ii)=sum(proj_packed(:,ii)*XijPYpsi(dir,:,s))  
        enddo; enddo
    
        do ii= 1, atom%RL%ns
                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                do dir=1,3
                        dipole_psi(dir,s)%R(ix,iy,iz)=dipole_psi(dir,s)%R(ix,iy,iz) + &
                        dvpsi_vec(dir,ii)
                enddo
        enddo
    endif
    enddo
    if(allocated(factor)) deallocate(factor)
    else
    print *, 'Unrecognized PP type', element%PP_type
    endif
    if(element%PP_type.eq.PAW) then
    if(present(Vpsi).or.present(midrVpsi)) deallocate(DijPYpsi)     
    if(present(midrVpsi)) deallocate(DijdPYpsi) 
    if(present(nabla_psi))   deallocate(TijPYpsi) 
    if(present(dipole_psi))   deallocate(XijPYpsi) 
    if(present(Spsi).or.present(midrSpsi).or.present(RSpsi))    deallocate(SPYpsi) 
    if(present(midrSpsi).or.present(SRpsi)) deallocate(SdPYpsi)
    endif
    deallocate(vpsi_vec)
    if(present(midrVpsi).or.present(midrSpsi) &
    .or.present(nabla_psi).or.present(dipole_psi).or.present(RSpsi).or.present(SRpsi)) deallocate(dvpsi_vec)
    deallocate(proj_packed)        
    enddo
    
    end subroutine


subroutine Apply_S_power(psi, atoms, elements, grids, parallel, S2psi, power)                     
        use fft, only: recip_to_real, real_to_recip
        use  parallel_mod, only: parallel_task, parallel_wait
        use odp_type, only : field_struct
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct, PAW, Real_local, Reciprocal
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        use grids_mod, only : allocate_local_fields_G
        type(field_struct), intent(in) :: psi(:)
        type(field_struct), intent(inout) ::  S2psi(:)

        type(parallel_struct), intent(in) :: parallel
        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer :: grid

        type(element_struct), intent(inout), target :: elements(:)
        type(atom_struct), pointer :: atom
        type(element_struct), pointer :: element
        real(dp) :: power
        complex(dp), allocatable :: PYpsi(:,:,:), SPYpsi(:,:), vpsi_vec(:), factor_G(:,:,:), factor(:)
        integer :: ii, ix, iy, iz, i_sphere, at, s, i_p, il
        real(dp), allocatable :: proj_packed(:,:)
 
        grid=> grids(psi(1)%grid)

        do s=1, size(psi)
                S2psi(s)%R=0.0_dp
                S2psi(s)%G=psi(s)%G
        enddo
        if(all(elements(:)%n_proj.eq.0))  return

        allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(psi)))
        PYpsi=0.0_dp

        do s=1, size(psi)
                call Calculate_Projector_overlaps(psi(s), PYpsi(:,:,s), atoms, elements, grid, parallel, &
                         RO=.true.,  skip_HGH=.true.) 
        enddo

        do at=1, size(atoms)
          atom=>atoms(at)
          element=>elements(atom%element)
          if(element%PP_type.ne.PAW) cycle
          allocate(SPYpsi(element%n_proj,size(psi)))
          do s=1, size(psi)
                SPYpsi(:,s)=((1.0_dp+atom%PAW%obar_i(:element%n_proj))**(power)-1.0_dp)*PYpsi(:element%n_proj,at,s)
          enddo

          if(element%PA_type.eq.Real_local) then
               
                allocate(proj_packed(element%n_proj,atom%RL%ns))
                allocate(vpsi_vec(atom%RL%ns))

                do ii= 1, atom%RL%ns
                        i_sphere=atom%RL%map(4,ii)
                        proj_packed(:,ii)=atom%PAW%projector_ort(i_sphere,:)
                enddo          
                
                do s=1, size(psi);
                        if((any(abs(grid%k).gt.tiny(1.0_dp)).or.any(abs(grid%A).gt.tiny(1.0_dp)))) then
                                allocate(factor(atom%RL%ns))
                                factor(:)=exp(-i_*((grid%k(1)+grid%A(1))*atom%RL%Rs(1,:atom%RL%ns) + &
                                                (grid%k(2)+grid%A(2))*atom%RL%Rs(2,:atom%RL%ns) + &
                                                (grid%k(3)+grid%A(3))*atom%RL%Rs(3,:atom%RL%ns)))

                                do ii= 1, atom%RL%ns
                                        vpsi_vec(ii)=sum(proj_packed(:,ii)*SPYpsi(:,s)) &
                                                *factor(ii)
                                enddo
                                deallocate(factor)
                        else
                                do ii= 1, atom%RL%ns
                                        vpsi_vec(ii)=sum(proj_packed(:,ii)*SPYpsi(:,s))
                                enddo
                        endif

                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii)
                                iy=atom%RL%map(2,ii)
                                iz=atom%RL%map(3,ii)
                                S2psi(s)%R(ix,iy,iz)=S2psi(s)%R(ix,iy,iz)+vpsi_vec(ii)
                        enddo
                enddo
                deallocate(vpsi_vec)
                deallocate(proj_packed)
          else if(element%PA_type.eq.Reciprocal) then
                call allocate_local_fields_G(factor_G,grid)
                factor_G(:,:,:)=exp(-i_*((grid%G(:,:,:,1)-grid%A(1))*atom%R(1) +&
                                         (grid%G(:,:,:,2)-grid%A(2))*atom%R(2) +&
                                         (grid%G(:,:,:,3)-grid%A(3))*atom%R(3)  ))
                do s=1, size(psi); do i_p=1,element%n_proj
                        il=element%PAW%tab%indlmn(1,i_p)
                        S2psi(s)%G(:,:,:)=S2psi(s)%G(:,:,:)+ (-i_)**il*&
                        factor_G(:,:,:)*element%GNL(psi(s)%grid)%projector_G_ort(:,:,:,i_p)*SPYpsi(i_p,s)
                enddo; enddo
                deallocate(factor_G)
          endif
          deallocate(SPYpsi)

        enddo
        deallocate(PYpsi)

        do s=1, size(psi)
                if(any(elements(:)%PA_type.eq.Real_Local)) then
                        grid%work%R=S2psi(s)%R
                        call real_to_recip(grid%work, grids)
                        S2psi(s)%G=(S2psi(s)%G+grid%work%G)*grid%cutwf
                else
                        S2psi(s)%G=S2psi(s)%G*grid%cutwf
                endif
                call recip_to_real(S2psi(s), grids)
        enddo

end subroutine

subroutine Apply_projectors_Recip2(psi, atoms, elements, grid, PYPsi,  Vpsi, skip_PAW, skip_HGH)                     
use fft, only: recip_to_real, real_to_recip
use odp_type, only : field_struct
use  atom_type, only : atom_struct
use  element_type, only : element_struct, HGH, PAW
use  grids_type, only : grid_struct
use  parallel_type, only : parallel_struct
use grids_mod, only : allocate_local_fields_G

type(field_struct), intent(in) :: psi(:)
type(field_struct), intent(inout), optional :: Vpsi(:)

type(atom_struct), intent(inout), target :: atoms(:)
type(grid_struct), intent(in) :: grid
type(element_struct), intent(inout), target :: elements(:)
real(dp), allocatable :: KB(:)
type(atom_struct), pointer :: atom
type(element_struct), pointer :: element
complex(dp), intent(in) :: PYPsi(:,:,:)
complex(dp), allocatable :: factor(:,:,:)
integer :: i, l, at, i_p, m, s, s2
logical, optional, intent(in) :: skip_PAW, skip_HGH
logical :: skip_PAW_, skip_HGH_

skip_PAW_=.false.; skip_HGH_=.false.; 
if(present(skip_PAW)) skip_PAW_=skip_PAW
if(present(skip_HGH)) skip_HGH_=skip_HGH

if(all(elements(:)%n_proj.eq.0))  return
if(.not.present(Vpsi)) then
print *, 'Ask for NL apply but dont ask for output, what are you trying to do?'
stop
endif

do s=1, size(psi)
if(present(Vpsi)) Vpsi(s)%G=0.0_dp
enddo

do at=1, size(atoms)
atom=>atoms(at)
element=>elements(atom%element)
if(skip_HGH_) then
if(element%PP_type.eq.HGH) cycle
endif
if(skip_PAW_) then
if(element%PP_type.eq.PAW) cycle
endif

if(element%PP_type.eq.HGH) then
        call allocate_local_fields_G(factor,grid)
        factor(:,:,:)=exp(-i_*((grid%G(:,:,:,1)-grid%A(1))*atom%R(1) +&
                               (grid%G(:,:,:,2)-grid%A(2))*atom%R(2) +&
                               (grid%G(:,:,:,3)-grid%A(3))*atom%R(3)  ))
        allocate(KB(element%n_proj))
        KB(:)=0.0_dp
        do s=1, size(psi); do s2=1, size(psi)
        do i_p=1,element%n_proj
        i=element%HGH%i(i_p)
        l=element%HGH%l(i_p)
        m=element%HGH%m(i_p)
        if(element%HGH%h_or_k(i_p).eq.'h') then
                if(s.eq.s2) then
                        KB(i_p)=element%HGH%EKB(i,l+1)
                else
                        cycle
                endif
        else
                if(s.eq.s2) then
                        KB(i_p)=element%HGH%KKB(i,l+1)*(s*1.0_dp-1.5_dp)*m
                else
                        KB(i_p)=element%HGH%KKB(i,l+1)*0.5_dp* &
                                sqrt(1.0_dp*(l+m+1))*sqrt(1.0_dp*(l-m))
                endif
        endif
        Vpsi(s2)%G(:,:,:)=Vpsi(s2)%G(:,:,:)+factor(:,:,:)*element%GNL(psi(s)%grid)%projector_G(:,:,:,i_p)*KB(i_p)*PYPsi(i_p,at,s)
        enddo
        enddo; enddo

        deallocate(factor)
        deallocate(KB)
else if(element%PP_type.eq.PAW) then
        call allocate_local_fields_G(factor,grid)
        factor(:,:,:)=exp(-i_*((grid%G(:,:,:,1)-grid%A(1)-grid%k(1))*atom%R(1) +&
                               (grid%G(:,:,:,2)-grid%A(2)-grid%k(2))*atom%R(2) +&
                               (grid%G(:,:,:,3)-grid%A(3)-grid%k(3))*atom%R(3)  ))

        Vpsi(s2)%G(:,:,:)=Vpsi(s2)%G(:,:,:)+factor(:,:,:)*element%GNL(psi(s)%grid)%projector_G(:,:,:,i_p)*KB(i_p)*PYPsi(i_p,at,s)

else
print *, 'Unrecognized PP type', element%PP_type
        stop
endif
enddo

end subroutine

subroutine Apply_S_inverse(psi, atoms, elements, grids, parallel, SIpsi, all_PAW, with_CG)                     
        use fft, only: recip_to_real, real_to_recip
        use  parallel_mod, only: parallel_task, parallel_wait
        use odp_type, only : field_struct
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct, Real_local, Reciprocal
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        use simulation_type, only : all_PAW_struct
        use operations_3D, only : integrate_3D_R
        use odp, only : allocate_field, deallocate_field
        use fft, only: real_to_recip, recip_to_real
        use grids_mod, only : allocate_local_fields_G

        type(field_struct), intent(in) :: psi(:)
        type(field_struct), intent(inout) ::  SIpsi(:)
        type(field_struct), allocatable ::  p(:), res(:), AP(:)

        type(parallel_struct), intent(in) :: parallel
        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(element_struct), intent(inout), target :: elements(:)
        type(atom_struct), pointer :: atom
        type(element_struct), pointer :: element
        type(all_PAW_struct), intent(inout) :: all_PAW
        logical, intent(in), optional :: with_CG

        complex(dp), allocatable :: PYpsi(:,:,:), SPYpsi(:,:), vpsi_vec(:)
        complex(dp), allocatable :: PYpsi_all(:,:), factor(:), factor_G(:,:,:)
        integer :: ii, ix, iy, iz, i_sphere, at, s, si, ei, all_proj, i_p, i, j, i_p_local, il
        real(dp), allocatable :: proj_packed(:,:)
        logical :: with_CG_
        real(dp) :: resnorm, beta
        complex(dp) :: alpha

        grid=> grids(psi(1)%grid)
        with_CG_=.false.
        if(present(with_CG)) with_CG_=with_CG

        do s=1, size(psi)
                SIpsi(s)%R=0.0_dp
                SIpsi(s)%G=psi(s)%G
        enddo
        if(all_PAW%N_PAW_atoms.lt.1) return
        if(all(elements(:)%n_proj.eq.0))  return

        allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(psi)))
        PYpsi=0.0_dp
        do s=1, size(psi)
                call Calculate_Projector_overlaps(psi(s), PYpsi(:,:,s), atoms, elements, grid, parallel, skip_HGH=.true.) 
        enddo

        all_proj=0
        do at=1, all_PAW%N_PAW_atoms
                atom=>atoms(all_PAW%all_at(at))
                element=>elements(atom%element)
                all_proj=all_proj+element%n_proj
        enddo
        allocate(PYpsi_all(all_proj, size(psi)))
         ei=0
         do at=1, all_PAW%N_PAW_atoms
                atom=>atoms(all_PAW%all_at(at))
                element=>elements(atom%element)
                si=ei+1
                ei=ei+element%n_proj
                do s=1, size(psi)
                        PYpsi_all(si:ei,s)=PYpsi(1:element%n_proj,all_PAW%all_at(at),s)
                enddo
        enddo
        deallocate(PYpsi)
        allocate(SPYpsi(all_proj, size(psi)))
        do s=1, size(psi); do i_p=1, all_proj
                SPYpsi(i_p,s)=sum(all_PAW%S_inv_full(:,i_p)*PYpsi_all(:,s))
        enddo;enddo
        deallocate(PYpsi_all)

        ei=0
        do at=1, all_PAW%N_PAW_atoms
          atom=>atoms(all_PAW%all_at(at))
          element=>elements(atom%element)
          si=ei+1
          ei=ei+element%n_proj

          if(element%PA_type.eq.Real_Local) then
                allocate(proj_packed(element%n_proj,atom%RL%ns))
                allocate(vpsi_vec(atom%RL%ns))

                do ii= 1, atom%RL%ns
                        i_sphere=atom%RL%map(4,ii)
                        proj_packed(:,ii)=atom%RL%projector(i_sphere,:)
                enddo          
                
                do s=1, size(psi);
                        if((any(abs(grid%k).gt.tiny(1.0_dp)).or.any(abs(grid%A).gt.tiny(1.0_dp)))) then
                                allocate(factor(atom%RL%ns))
                                factor(:)=exp(-i_*((grid%k(1)+grid%A(1))*atom%RL%Rs(1,:atom%RL%ns) + &
                                                (grid%k(2)+grid%A(2))*atom%RL%Rs(2,:atom%RL%ns) + &
                                                (grid%k(3)+grid%A(3))*atom%RL%Rs(3,:atom%RL%ns)))

                                do ii= 1, atom%RL%ns
                                        vpsi_vec(ii)=sum(proj_packed(:,ii)*SPYpsi(si:ei,s)) &
                                                *factor(ii)
                                enddo
                                deallocate(factor)
                        else
                                do ii= 1, atom%RL%ns
                                        vpsi_vec(ii)=sum(proj_packed(:,ii)*SPYpsi(si:ei,s))
                                enddo
                        endif
                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii)
                                iy=atom%RL%map(2,ii)
                                iz=atom%RL%map(3,ii)
                                SIpsi(s)%R(ix,iy,iz)=SIpsi(s)%R(ix,iy,iz)+vpsi_vec(ii)
                        enddo
                enddo
                deallocate(vpsi_vec)
                deallocate(proj_packed)
          else if(element%PA_type.eq.Reciprocal) then
                call allocate_local_fields_G(factor_G,grid)
                factor_G(:,:,:)=exp(-i_*((grid%G(:,:,:,1)-grid%A(1))*atom%R(1) +&
                                         (grid%G(:,:,:,2)-grid%A(2))*atom%R(2) +&
                                         (grid%G(:,:,:,3)-grid%A(3))*atom%R(3)  ))
                do s=1, size(psi);do i_p=si,ei
                        i_p_local=i_p-si+1
                        il=element%PAW%tab%indlmn(1,i_p_local)
                        SIpsi(s)%G(:,:,:)=SIpsi(s)%G(:,:,:)+ (-i_)**il*&
                        factor_G(:,:,:)*element%GNL(psi(s)%grid)%projector_G(:,:,:,i_p_local)*SPYpsi(i_p,s)
                        
                enddo; enddo
                deallocate(factor_G)
          endif
        enddo
        deallocate(SPYpsi)

        do i=1, size(psi)
                if(any(elements(:)%PA_type.eq.Real_Local)) then
                        grid%work%R=SIpsi(i)%R
                        call real_to_recip(grid%work, grids)
                        SIpsi(i)%G=(SIpsi(i)%G+grid%work%G)*grid%cutwf
                else
                        SIpsi(i)%G=SIpsi(i)%G*grid%cutwf
                endif
                call recip_to_real(SIpsi(i), grids)
        enddo

        if(with_CG_) then
                !!Conjugent Gradient Refinement====================================================
                !==================================================================================
                !(S^-1 psi) = psi' solve S(psi')=psi, use RO inv as guess and preconditioner
                allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(psi)))
                allocate(p(size(psi)))
                allocate(res(size(psi)))
                allocate(Ap(size(psi)))
                do i=1, size(psi)
                        call allocate_field(p(i), grid, parallel)
                        call allocate_field(Ap(i), grid, parallel)
                        call allocate_field(res(i), grid, parallel)
                enddo

                p%grid=psi%grid
                res%grid=psi%grid
                Ap%grid=psi%grid

                !r_0= psi-S(S^-1psi)

                call Apply_S_RR(SIpsi, res, grids, atoms, elements, parallel)

                do i=1, size(psi)
                        call real_to_recip(res(i), grids)
                        res(i)%G=res(i)%G*grid%cutwf
                        call recip_to_real(res(i), grids)
                enddo

                resnorm=0.0_dp
                do i=1, size(psi)
                        !r=b-Ax=psi-S(S^-1psi)
                        res(i)%R=res(i)%R-psi(i)%R
                        !p_0=R_0
                        p(i)%R=res(i)%R
                        resnorm=resnorm+abs(integrate_3D_R(conjg(res(i)%R)*res(i)%R, grid, parallel))
                enddo


                do j=1, 100
                        !if(parallel%myid.eq.0)  print *, 'step:', j, resnorm
                        if(sqrt(abs(resnorm)).lt.(1E-12_dp)) exit

                        !AP=Sp
                        call Apply_S_RR(p, Ap, grids, atoms, elements, parallel)
                      !  do i=1, size(psi)
                      !          call real_to_recip(Ap(i), grids)
                      !          Ap(i)%G=Ap(i)%G*grid%cutwf
                      !          call recip_to_real(Ap(i), grids)
                      !  enddo
                        alpha=0.0_dp
                        do i=1, size(psi)
                                alpha=alpha+integrate_3D_R(Ap(i)%R*conjg(p(i)%R),grid, parallel)
                        enddo
                        alpha=resnorm/alpha
                        beta=resnorm                        
                        do i=1, size(psi)
                                !x=x+alpha*p
                                SIpsi(i)%R=SIpsi(i)%R-alpha*p(i)%R
                        enddo
                        resnorm=0.0_dp
                        do i=1, size(psi)
                                !r=r-alpha*Ap
                                res(i)%R=res(i)%R-alpha*Ap(i)%R
                                resnorm=resnorm+abs(integrate_3D_R(conjg(res(i)%R)*res(i)%R, grid, parallel))
                        enddo
                        beta=resnorm/beta
                        do i=1, size(psi)
                                !p=r+beta*p
                                p(i)%R=beta*p(i)%R
                                p(i)%R=res(i)%R+p(i)%R
                        enddo
                enddo
                do i=1, size(psi)
                        call deallocate_field(p(i),grids(p(i)%grid))
                        call deallocate_field(Ap(i),grids(Ap(i)%grid))
                        call deallocate_field(res(i),grids(res(i)%grid))
                enddo
                deallocate(p)
                deallocate(res)
                deallocate(Ap)
                deallocate(PYpsi)

                do i=1, size(psi)
                        call real_to_recip(SIpsi(i), grids)
                        SIpsi(i)%G=SIpsi(i)%G*grid%cutwf
                        call recip_to_real(SIpsi(i), grids)
                enddo
        endif


end subroutine

subroutine Apply_S_RR(field, Sfield, grids, atoms, elements, parallel, PYpsi_in)
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct

        use odp_type, only : field_struct
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct, Reciprocal
        use atom_type, only : atom_struct
        use fft, only: real_to_recip, recip_to_real

        type(parallel_struct), intent(in) :: parallel
        type(field_struct), intent(inout) :: field(:)
        type(field_struct), intent(inout) :: Sfield(:)
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        complex(dp), intent(in), optional, target :: PYPsi_in(:,:,:)
        complex(dp), pointer :: PYPsi(:,:,:)

        integer :: s

        grid=>grids(field(1)%grid)

        if(any(elements(:)%PA_type.eq.Reciprocal)) then
                do s=1,size(field)
                        call real_to_recip(field(s), grids)
                enddo
        endif

        if(present(PYpsi_in)) then
            PYpsi=>PYpsi_in
        else
            allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(field)))
            PYpsi=0.0_dp
            do s=1, size(field)
                call Calculate_Projector_overlaps(field(s),  &
                PYpsi(:,:,s), atoms, elements, grid, parallel) 
             enddo
        endif
        
        if(grid%gamma) then
                call Apply_projectors_RL_gamma(field, atoms, elements, grid, PYpsi, Spsi=Sfield) 
            else
                call Apply_projectors_RL(field, atoms, elements, grid, PYpsi, Spsi=Sfield) 
            endif
        if(any(elements(:)%PA_type.eq.Reciprocal)) then
                call Apply_projectors_Recip(field, atoms, elements, grid, PYpsi, Spsi=Sfield)
                do s=1, size(field)
                        grid%work%G=Sfield(s)%G*grid%cutwf
                        call recip_to_real(grid%work, grids)
                        Sfield(s)%R=Sfield(s)%R+grid%work%R
                enddo
        endif
        if(.not.present(PYpsi_in)) then
            deallocate(PYpsi)
        else
            nullify(PYpsi)
        endif
    end subroutine

    subroutine Apply_S_m_one_half(field, Sfield, grids, atoms, elements, parallel, calc_R)
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct

        use odp_type, only : field_struct
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: real_to_recip, recip_to_real
        use odp, only : allocate_field, deallocate_field

        type(parallel_struct), intent(in) :: parallel
        type(field_struct), intent(inout) :: field(:)
        type(field_struct), intent(inout) :: Sfield(:)
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        logical, intent(in) :: calc_R

        type(grid_struct), pointer:: grid
        complex(dp), pointer :: PYPsi(:,:,:)
        type(field_struct), allocatable:: X_nm1_field(:), X_n_field(:)
        integer :: s, n

        real(dp) :: coefficients(100)

        coefficients=(/-0.5_dp,0.375_dp,-0.3125_dp,0.273438_dp,-0.246094_dp,0.225586_dp,-0.209473_dp,&
        0.196381_dp,-0.185471_dp,0.176197_dp,-0.168188_dp,0.16118_dp,-0.154981_dp,&
        0.149446_dp,-0.144464_dp,0.13995_dp,-0.135834_dp,0.132061_dp,-0.128585_dp,&
        0.125371_dp,-0.122386_dp,0.119604_dp,-0.117004_dp,0.114567_dp,-0.112275_dp,&
        0.110116_dp,-0.108077_dp,0.106147_dp,-0.104317_dp,0.102578_dp,-0.100924_dp,&
        0.0993468_dp,-0.0978415_dp,0.0964027_dp,-0.0950255_dp,0.0937057_dp,-0.0924394_dp,&
        0.0912231_dp,-0.0900535_dp,0.0889279_dp,-0.0878434_dp,0.0867976_dp,-0.0857884_dp,&
        0.0848135_dp,-0.0838711_dp,0.0829595_dp,-0.0820769_dp,0.081222_dp,-0.0803932_dp,&
        0.0795892_dp,-0.078809_dp,0.0780512_dp,-0.0773148_dp,0.076599_dp,-0.0759026_dp,&
        0.0752249_dp,-0.074565_dp,0.0739222_dp,-0.0732958_dp,0.072685_dp,-0.0720892_dp,&
        0.0715078_dp,-0.0709403_dp,0.0703861_dp,-0.0698447_dp,0.0693155_dp,-0.0687983_dp,&
        0.0682924_dp,-0.0677975_dp,0.0673132_dp,-0.0668392_dp,0.066375_dp,-0.0659204_dp,&
        0.065475_dp,-0.0650385_dp,0.0646106_dp,-0.0641911_dp,0.0637796_dp,-0.0633759_dp,&
        0.0629798_dp,-0.0625911_dp,0.0622094_dp,-0.0618347_dp,0.0614666_dp,-0.061105_dp,&
        0.0607498_dp,-0.0604006_dp,0.0600574_dp,-0.05972_dp,0.0593883_dp,-0.059062_dp,&
        0.058741_dp,-0.0584252_dp,0.0581144_dp,-0.0578085_dp,0.0575074_dp,-0.057211_dp,&
        0.0569191_dp,-0.0566316_dp, 0.0563485_dp /)

        grid=>grids(field(1)%grid)

        !S^-1/2=(1+X)^-1/2
        !X=S-1

        allocate(X_nm1_field(size(field)))
        allocate(X_n_field(size(field)))

        do s=1,size(X_nm1_field)
            grid=>grids(field(s)%grid)
            call allocate_field(X_nm1_field(s), grid, parallel)
            call allocate_field(X_n_field(s), grid, parallel)
            X_nm1_field(s)%grid=field(s)%grid
            X_n_field(s)%grid=field(s)%grid
        enddo
        allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(field)))

        if(calc_R) then
                do s=1,size(field)
                    call recip_to_real(field(s), grids)
                enddo
        endif

        do s=1,size(field)
                X_nm1_field(s)%R=field(s)%R
                Sfield(s)%R=field(s)%R
        enddo

        do n=1,10
                PYpsi=0.0_dp
                !S(X^(n-1))
                do s=1, size(field)
                        call Calculate_Projector_overlaps(X_nm1_field(s),  &
                        PYpsi(:,:,s), atoms, elements, grid, parallel) 
                enddo
                call Apply_projectors_RL(X_nm1_field, atoms, elements, grid, PYpsi, Spsi=X_n_field) 
                
                do s=1,size(field)
                        !X^n=(S-1)(X^n-1)
                        X_n_field(s)%R=X_n_field(s)%R-X_nm1_field(s)%R
                        Sfield(s)%R=Sfield(s)%R+coefficients(n)*X_n_field(s)%R
                        X_nm1_field(s)%R=X_n_field(s)%R
                enddo
        enddo
        
        do s=1,size(field)
                call real_to_recip(Sfield(s), grids)
                Sfield(s)%G=Sfield(s)%G*grids(Sfield(s)%grid)%cutwf
                call recip_to_real(Sfield(s), grids)
        enddo

        deallocate(PYpsi)
        do s=1,size(field)
                call deallocate_field(X_nm1_field(s),grids(X_nm1_field(s)%grid))
                call deallocate_field(X_n_field(s),grids(X_n_field(s)%grid))
        enddo
        deallocate(X_nm1_field)
        deallocate(X_n_field)

    end subroutine

subroutine Apply_S_rotation(psi, atoms, elements, grid, parallel, dR, S2psi, all_PAW)                     
        use odp_type, only : field_struct
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct, Real_local, PAW
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        use simulation_type, only : all_PAW_struct

        type(field_struct), intent(in) :: psi(:)
        type(field_struct), intent(inout) ::  S2psi(:)

        type(parallel_struct), intent(in) :: parallel
        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(in) :: grid
        type(element_struct), intent(inout), target :: elements(:)
        type(atom_struct), pointer :: atom
        type(element_struct), pointer :: element
        type(all_PAW_struct), intent(inout) :: all_PAW
        real(dp), intent(in) :: dR(:,:)
        complex(dp), allocatable :: PYpsi(:,:,:), vpsi_vec(:), dPYpsi(:,:,:,:), &
        MijPYpsi(:,:)
        integer :: ii, ix, iy, iz, i_sphere, at, s, i_p, dir, i
        real(dp), allocatable :: proj_packed(:,:), tmp(:)

        do s=1, size(psi)
                S2psi(s)%R=0.0_dp
        enddo
        if(all(elements(:)%n_proj.eq.0))  return
        if(all_PAW%N_PAW_atoms.lt.1) return

        allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(psi)))
        allocate(dPYpsi(3,maxval(elements(:)%n_proj), size(atoms), size(psi)))
        dPYpsi=0.0_dp

        PYpsi=0.0_dp
        do s=1, size(psi)
                call Calculate_Projector_overlaps(psi(s), PYpsi(:,:,s), atoms, elements, grid, parallel, skip_HGH=.true.)
                call Calculate_deriv_Projector_overlaps(psi(s), dPYpsi(:,:,:,s), atoms, elements, grid, parallel) 
        enddo
        do at=1, size(atoms)
                if(all(abs(dR(:,at)).lt.tiny(1.0_dp))) cycle
                atom=>atoms(at)
                element=>elements(atom%element)
                if(element%PA_type.ne.Real_local) cycle
                if(element%PP_type.ne.PAW) cycle
                if(atom%RL%ns.lt.1) cycle
                allocate(MijPYpsi(element%n_proj, size(psi)))
                allocate(tmp(element%n_proj))       
                do s=1, size(psi); do i_p=1, element%n_proj
                        tmp(:)=0.0_dp!element%PAW%sij(:,i_p)
                        do dir=1,3    
                                tmp(:)=tmp(:)-element%PAW%tab%nabla_ij(dir,i_p,:)*dR(dir,at)
                        enddo
                        MijPYpsi(i_p,s)=sum(tmp(:)*PYpsi(:element%n_proj,at,s))
                enddo; enddo
                deallocate(tmp)
                do s=1, size(psi); do i=i_p, element%n_proj; do dir=1,3               
                        MijPYpsi(i_p,s)=MijPYpsi(i_p,s)+sum(element%PAW%sij(:,i_p)*&
                                                dPYpsi(dir,:element%n_proj,at,s))*dR(dir,at)
                enddo; enddo; enddo

                allocate(vpsi_vec(atom%RL%ns))
                allocate(proj_packed(element%n_proj,atom%RL%ns))
                do ii= 1, atom%RL%ns
                      i_sphere=atom%RL%map(4,ii)
                      proj_packed(:,ii)=atoms(at)%RL%projector(i_sphere,:)
                enddo
                do s=1, size(psi); 
                        do ii= 1, atom%RL%ns
                                vpsi_vec(ii)=sum(proj_packed(:,ii)*MijPYpsi(:,s))
                        enddo
                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii);iy=atom%RL%map(2,ii);iz=atom%RL%map(3,ii)
                                S2psi(s)%R(ix,iy,iz)=S2psi(s)%R(ix,iy,iz)+vpsi_vec(ii)
                        enddo
                enddo
                deallocate(vpsi_vec)
                deallocate(proj_packed)
                deallocate(MijPYpsi)
        enddo

        deallocate(PYpsi)
        deallocate(dPYpsi)
end subroutine

subroutine Calculate_Projector_overlaps(psi, PYpsi, atoms, elements, grid, parallel, PYirpsi, skip_PAW, skip_HGH, RO)
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        use odp_type, only : field_struct

        type(field_struct),intent(in) :: psi
        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(in) :: grid
        type(element_struct), intent(inout), target :: elements(:)
        type(parallel_struct), intent(in) :: parallel
        complex(dp),intent(inout) :: PYpsi(:,:)
        complex(dp),intent(inout) , optional ::  PYirpsi(:,:,:)

        real(dp), allocatable :: psi_R(:,:,:)
        real(dp),allocatable :: PYpsi_R(:,:)
        real(dp),allocatable ::  PYirpsi_R(:,:,:)

        logical, optional, intent(in) :: skip_PAW, skip_HGH, RO
        logical :: skip_PAW_, skip_HGH_, RO_

        skip_PAW_=.false.; skip_HGH_=.false.; RO_=.false.
        if(present(skip_PAW)) skip_PAW_=skip_PAW
        if(present(skip_HGH)) skip_HGH_=skip_HGH
        if(present(RO)) RO_=RO


        !Real Space PP's
        if(grid%gamma) then
                allocate(psi_R(size(psi%R,1),size(psi%R,2),size(psi%R,3)))
                allocate(PYpsi_R(size(PYpsi,1),size(PYpsi,2)))

                psi_R=real(psi%R)
                if(present(PYirpsi)) then
                        allocate(PYirpsi_R(size(PYirpsi,1),size(PYirpsi,2),size(PYirpsi,3)))
                        call Calculate_Projector_overlaps_gamma(&
                                psi_R, PYpsi_R, atoms, elements, grid, parallel, PYirpsi_R, skip_PAW, skip_HGH, RO)
                        PYirpsi=PYirpsi_R
                        deallocate(PYirpsi_R)
                else
                        call Calculate_Projector_overlaps_gamma(&
                                psi_R, PYpsi_R, atoms, elements, grid, parallel, skip_PAW=skip_PAW_, skip_HGH=skip_HGH_, RO=RO_)
                endif
                PYpsi=PYpsi_R
                deallocate(psi_R)
                deallocate(PYpsi_R)
        else
                if(present(PYirpsi)) then
                        call Calculate_Projector_overlaps_complex(&
                                psi%R, PYpsi, atoms, elements, grid, parallel, PYirpsi, skip_PAW, skip_HGH, RO)
                else
                        call Calculate_Projector_overlaps_complex(&
                                psi%R, PYpsi, atoms, elements, grid, parallel, skip_PAW=skip_PAW_, skip_HGH=skip_HGH_, RO=RO_)
                endif               
        endif
        !Reciprocal Space PP's (must be after Real-Space)
        call Calculate_Projector_overlaps_Recip( &
                psi%G, psi%grid, PYpsi, atoms, elements, grid, parallel, skip_PAW=skip_PAW_, skip_HGH=skip_HGH_, RO=RO_)

end subroutine

subroutine Calculate_Projector_overlaps_complex(psi, PYpsi, atoms, elements, grid, parallel, PYirpsi, skip_PAW, skip_HGH, RO) 
        use lapack, only : dgemv
        use fft, only: recip_to_real, real_to_recip
        use  parallel_mod, only: parallel_task, parallel_wait
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct, Real_local, HGH, PAW
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct


        complex(dp),intent(in) :: psi(:,:,:)
        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(in) :: grid
        type(element_struct), intent(inout), target :: elements(:)
        type(parallel_struct), intent(in) :: parallel
        complex(dp),intent(inout) :: PYpsi(:,:)
        complex(dp),intent(inout) , optional ::  PYirpsi(:,:,:)
        complex(dp), allocatable :: factor(:), Rfactor(:,:), PYpsi_tmp(:),PYpsi_save(:,:),PYirpsi_save(:,:,:)
        real(dp), allocatable :: proj_packed(:,:)

        type(atom_struct), pointer :: atom
        type(element_struct), pointer :: element
        integer :: ii, ix, iy, iz, i_sphere, at, i_p, dir

        logical, optional, intent(in) :: skip_PAW, skip_HGH, RO
        logical :: skip_PAW_, skip_HGH_, RO_
        real(dp) :: time_1, time_2

        call cpu_time(time_1)
        skip_PAW_=.false.; skip_HGH_=.false.; RO_=.false.
        if(present(skip_PAW)) skip_PAW_=skip_PAW
        if(present(skip_HGH)) skip_HGH_=skip_HGH
        if(present(RO)) RO_=RO

        if((skip_HGH_.and.any(elements(:)%PP_type.eq.HGH)) .or. &
           (skip_PAW_.and.any(elements(:)%PP_type.eq.PAW))) then
           allocate(PYpsi_save(size(PYpsi,1),size(PYpsi,2)))
           PYpsi_save=PYpsi
           if(present(PYirpsi)) then
                allocate(PYirpsi_save(size(PYirpsi,1),size(PYirpsi,2),size(PYirpsi,3)))
                PYirpsi_save=PYirpsi
           endif
        endif
        do at=1, size(atoms)
                atom=>atoms(at)
                element=>elements(atom%element)
                if(skip_HGH_) then
                        if(element%PP_type.eq.HGH) cycle
                endif
                if(skip_PAW_) then
                        if(element%PP_type.eq.PAW) cycle
                endif
                if(element%PA_type.ne.Real_local) cycle
                PYpsi(:,at)=0.0_dp
                if(present(PYirpsi)) PYirpsi(:,:,at)=0.0_dp

                if(atom%RL%ns.lt.1) cycle
                allocate(factor(atom%RL%ns))
                allocate(proj_packed(atom%RL%ns,element%n_proj))
                if((any(abs(grid%k).gt.tiny(1.0_dp)).or.any(abs(grid%A).gt.tiny(1.0_dp)))) then
                        factor(:)=exp(i_*((grid%k(1)+grid%A(1))*atom%RL%Rs(1,:atom%RL%ns) + &
                                        (grid%k(2)+grid%A(2))*atom%RL%Rs(2,:atom%RL%ns) + &
                                        (grid%k(3)+grid%A(3))*atom%RL%Rs(3,:atom%RL%ns)))
                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii)
                                iy=atom%RL%map(2,ii)
                                iz=atom%RL%map(3,ii)
                                factor(ii)=factor(ii)*psi(ix,iy,iz)
                        enddo

                else
                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii)
                                iy=atom%RL%map(2,ii)
                                iz=atom%RL%map(3,ii)
                                factor(ii)=psi(ix,iy,iz)
                        enddo
                endif
                if(RO_) then 
                        do ii= 1, atom%RL%ns
                                i_sphere=atom%RL%map(4,ii)
                                proj_packed(ii,:)=atoms(at)%PAW%projector_ort(i_sphere,:)
                        enddo
                else
                        do ii= 1, atom%RL%ns
                                i_sphere=atom%RL%map(4,ii)
                                proj_packed(ii,:)=atoms(at)%RL%projector(i_sphere,:)
                        enddo
                endif
                do i_p=1, element%n_proj
                        PYpsi(i_p,at)= sum(proj_packed(:,i_p)*factor(:))
                enddo


                if(present(PYirpsi)) then
                        allocate(Rfactor(atom%RL%ns,3))
                        allocate(PYpsi_tmp(element%n_proj))
                        do dir= 1, 3
                                Rfactor(:,dir)=atom%RL%Rs(dir,:atom%RL%ns)*factor(:)
                        enddo
                        do dir=1,3
                                do i_p=1, element%n_proj
                                        PYpsi_tmp(i_p) = sum(proj_packed(:,i_p)*Rfactor(:,dir))
                                enddo
                                PYirpsi(dir,:element%n_proj,at)=PYpsi_tmp(:)
                        enddo
                        deallocate(Rfactor)
                        deallocate(PYpsi_tmp)
                endif
        
                deallocate(factor) 
                deallocate(proj_packed)  
        enddo !atoms
        PYpsi(:,:)=PYpsi(:,:)*product(grid%box_length(:)/grid%Nr(:))
        if(present(PYirpsi)) PYirpsi(:,:,:)=PYirpsi(:,:,:)*product(grid%box_length(:)/grid%Nr(:))

        !since threads only saw their own local contributions we need to sum (we
        !wait untill here so that the threads can work through all the atoms
        !separatly and hopefully take ~the same time to finish (else for any one
        !atom most threads don't contibute
        call parallel_task('sum', PYpsi, parallel, 'band')  
        if(present(PYirpsi)) call parallel_task('sum', PYirpsi, parallel, 'band')

        if((skip_HGH_.and.any(elements(:)%PP_type.eq.HGH)) .or. &
           (skip_PAW_.and.any(elements(:)%PP_type.eq.PAW))) then
                do at=1, size(atoms)
                        atom=>atoms(at)
                        element=>elements(atom%element)
                        if(skip_HGH_) then
                                if(element%PP_type.ne.HGH) cycle
                        endif
                        if(skip_PAW_) then
                                if(element%PP_type.ne.PAW) cycle
                        endif
                        PYpsi(:,at)=PYpsi_save(:,at)
                        if(present(PYirpsi)) PYirpsi(:,:,at)=PYirpsi_save(:,:,at)
                enddo
                deallocate(PYpsi_save)
                if(present(PYirpsi)) deallocate(PYirpsi_save)
        endif
        call cpu_time(time_2)
        !print *, parallel%myid, 'PYpsi time: ',time_2-time_1

end subroutine 

subroutine Calculate_Projector_overlaps_gamma(psi, PYpsi, atoms, elements, grid, parallel, PYirpsi, skip_PAW, skip_HGH, RO) 
        use lapack, only : dgemv
        use fft, only: recip_to_real, real_to_recip
        use  parallel_mod, only: parallel_task, parallel_wait
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct, Real_local, HGH, PAW
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct


        real(dp),intent(in) :: psi(:,:,:)
        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(in) :: grid
        type(element_struct), intent(inout), target :: elements(:)
        type(parallel_struct), intent(in) :: parallel
        real(dp),intent(inout) :: PYpsi(:,:)
        real(dp),intent(inout) , optional ::  PYirpsi(:,:,:)
        real(dp), allocatable :: factor(:), Rfactor(:,:), PYpsi_tmp(:),PYpsi_save(:,:),PYirpsi_save(:,:,:)
        real(dp), allocatable :: proj_packed(:,:)

        type(atom_struct), pointer :: atom
        type(element_struct), pointer :: element
        integer :: ii, ix, iy, iz, i_sphere, at, i_p, dir

        logical, optional, intent(in) :: skip_PAW, skip_HGH, RO
        logical :: skip_PAW_, skip_HGH_, RO_
        real(dp) :: time_1, time_2

        call cpu_time(time_1)
        skip_PAW_=.false.; skip_HGH_=.false.; RO_=.false.
        if(present(skip_PAW)) skip_PAW_=skip_PAW
        if(present(skip_HGH)) skip_HGH_=skip_HGH
        if(present(RO)) RO_=RO

        if((skip_HGH_.and.any(elements(:)%PP_type.eq.HGH)) .or. &
           (skip_PAW_.and.any(elements(:)%PP_type.eq.PAW))) then
           allocate(PYpsi_save(size(PYpsi,1),size(PYpsi,2)))
           PYpsi_save=PYpsi
           if(present(PYirpsi)) then
                allocate(PYirpsi_save(size(PYirpsi,1),size(PYirpsi,2),size(PYirpsi,3)))
                PYirpsi_save=PYirpsi
           endif
        endif
        do at=1, size(atoms)
                atom=>atoms(at)
                element=>elements(atom%element)
                if(skip_HGH_) then
                        if(element%PP_type.eq.HGH) cycle
                endif
                if(skip_PAW_) then
                        if(element%PP_type.eq.PAW) cycle
                endif
                if(element%PA_type.ne.Real_local) cycle
                PYpsi(:,at)=0.0_dp
                if(present(PYirpsi)) PYirpsi(:,:,at)=0.0_dp

                if(atom%RL%ns.lt.1) cycle
                allocate(factor(atom%RL%ns))
                allocate(proj_packed(atom%RL%ns,element%n_proj))
                if((any(abs(grid%k).gt.tiny(1.0_dp)).or.any(abs(grid%A).gt.tiny(1.0_dp)))) then
                        print *, "ERRROR: you have finite k or A on a Gamma point grid"; flush(6)
                        stop
                else
                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii)
                                iy=atom%RL%map(2,ii)
                                iz=atom%RL%map(3,ii)
                                factor(ii)=psi(ix,iy,iz)
                        enddo
                endif
                if(RO_) then 
                        do ii= 1, atom%RL%ns
                                i_sphere=atom%RL%map(4,ii)
                                proj_packed(ii,:)=atoms(at)%PAW%projector_ort(i_sphere,:)
                        enddo
                else
                        do ii= 1, atom%RL%ns
                                i_sphere=atom%RL%map(4,ii)
                                proj_packed(ii,:)=atoms(at)%RL%projector(i_sphere,:)
                        enddo
                endif
                do i_p=1, element%n_proj
                        PYpsi(i_p,at)= sum(proj_packed(:,i_p)*factor(:))
                enddo

                if(present(PYirpsi)) then
                        allocate(Rfactor(atom%RL%ns,3))
                        allocate(PYpsi_tmp(element%n_proj))
                        do dir= 1, 3
                                Rfactor(:,dir)=atom%RL%Rs(dir,:atom%RL%ns)*factor(:)
                        enddo
                        do dir=1,3
                                do i_p=1, element%n_proj
                                        PYpsi_tmp(i_p) = sum(proj_packed(:,i_p)*Rfactor(:,dir))
                                enddo
                                PYirpsi(dir,:element%n_proj,at)=PYpsi_tmp(:)
                        enddo
                        deallocate(Rfactor)
                        deallocate(PYpsi_tmp)
                endif
        
                deallocate(factor) 
                deallocate(proj_packed)  
        enddo !atoms
        PYpsi(:,:)=PYpsi(:,:)*product(grid%box_length(:)/grid%Nr(:))
        if(present(PYirpsi)) PYirpsi(:,:,:)=PYirpsi(:,:,:)*product(grid%box_length(:)/grid%Nr(:))

        !since threads only saw their own local contributions we need to sum (we
        !wait untill here so that the threads can work through all the atoms
        !separatly and hopefully take ~the same time to finish (else for any one
        !atom most threads don't contibute
        call parallel_task('sum', PYpsi, parallel, 'band')  
        if(present(PYirpsi)) call parallel_task('sum', PYirpsi, parallel, 'band')

        if((skip_HGH_.and.any(elements(:)%PP_type.eq.HGH)) .or. &
           (skip_PAW_.and.any(elements(:)%PP_type.eq.PAW))) then
                do at=1, size(atoms)
                        atom=>atoms(at)
                        element=>elements(atom%element)
                        if(skip_HGH_) then
                                if(element%PP_type.ne.HGH) cycle
                        endif
                        if(skip_PAW_) then
                                if(element%PP_type.ne.PAW) cycle
                        endif
                        PYpsi(:,at)=PYpsi_save(:,at)
                        if(present(PYirpsi)) PYirpsi(:,:,at)=PYirpsi_save(:,:,at)
                enddo
                deallocate(PYpsi_save)
                if(present(PYirpsi)) deallocate(PYirpsi_save)
        endif
        call cpu_time(time_2)
        !print *, parallel%myid, 'PYpsi time: ',time_2-time_1

end subroutine

subroutine Calculate_Projector_overlaps_Recip(psi, k, PYpsi, atoms, elements, grid, parallel, skip_PAW, skip_HGH, RO) 
        use lapack, only : dgemv
        use fft, only: recip_to_real, real_to_recip
        use  parallel_mod, only: parallel_task, parallel_wait
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct, HGH, PAW, Reciprocal, Real_local
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        use operations_3D,only : integrate_3D_G
        use grids_mod, only : allocate_local_fields_G

        complex(dp),intent(in) :: psi(:,:,:)
        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(in) :: grid
        type(element_struct), intent(inout), target :: elements(:)
        type(parallel_struct), intent(in) :: parallel
        complex(dp),intent(inout) :: PYpsi(:,:)
        complex(dp), allocatable :: psi_shift(:,:,:), PYpsi_save(:,:)
        integer, intent(in) :: k

        type(atom_struct), pointer :: atom
        type(element_struct), pointer :: element
        integer :: at, i_p, il

        logical, optional, intent(in) :: skip_PAW, skip_HGH, RO
        logical :: skip_PAW_, skip_HGH_, RO_
        real(dp) :: time_1, time_2
        complex(dp) :: i_factor

        call allocate_local_fields_G(psi_shift,grid)

        call cpu_time(time_1)
        skip_PAW_=.false.; skip_HGH_=.false.; RO_=.false.
        if(present(skip_PAW)) skip_PAW_=skip_PAW
        if(present(skip_HGH)) skip_HGH_=skip_HGH
        if(present(RO)) RO_=RO

        if((skip_HGH_.and.any(elements(:)%PP_type.eq.HGH)) .or. &
           (skip_PAW_.and.any(elements(:)%PP_type.eq.PAW)) .or. &
           any(elements(:)%PA_type.eq.Real_local)) then
           allocate(PYpsi_save(size(PYpsi,1),size(PYpsi,2)))
           PYpsi_save=PYpsi
        endif
        do at=1, size(atoms)
                atom=>atoms(at)
                element=>elements(atom%element)
                if(skip_HGH_) then
                        if(element%PP_type.eq.HGH) cycle
                endif
                if(skip_PAW_) then
                        if(element%PP_type.eq.PAW) cycle
                endif
                if(element%PA_type.ne.Reciprocal) cycle

                PYpsi(:,at)=0.0_dp
                psi_shift(:,:,:)=psi(:,:,:)*exp(i_*((grid%G(:,:,:,1)-grid%A(1))*atom%R(1) +&
                                                    (grid%G(:,:,:,2)-grid%A(2))*atom%R(2) +&
                                                    (grid%G(:,:,:,3)-grid%A(3))*atom%R(3)  ))
                do i_p=1, element%n_proj
                        i_factor=1.0_dp
                        if(element%PP_type.eq.PAW) then
                                il=element%PAW%tab%indlmn(1,i_p)
                                i_factor=i_**il 
                        endif
                        if(RO_) then 
                                PYpsi(i_p,at)=integrate_3D_G( i_factor * &
                                           element%GNL(k)%projector_G_ort(:,:,:,i_p)* &
                                           psi_shift(:,:,:), grid, parallel)
                        else
                                PYpsi(i_p,at)=integrate_3D_G( i_factor * &
                                        element%GNL(k)%projector_G(:,:,:,i_p)* &
                                        psi_shift(:,:,:), grid, parallel)
                        endif
                enddo


        enddo !atoms
        deallocate(psi_shift)


        if((skip_HGH_.and.any(elements(:)%PP_type.eq.HGH)) .or. &
           (skip_PAW_.and.any(elements(:)%PP_type.eq.PAW)) .or. &
           any(elements(:)%PA_type.eq.Real_local)) then
                do at=1, size(atoms)
                        atom=>atoms(at)
                        element=>elements(atom%element)
                        if(element%PA_type.eq.Real_local) then
                                PYpsi(:,at)=PYpsi_save(:,at)
                        endif
                        if(skip_HGH_) then
                                if(element%PP_type.eq.HGH) PYpsi(:,at)=PYpsi_save(:,at)
                        endif
                        if(skip_PAW_) then
                                if(element%PP_type.eq.PAW) PYpsi(:,at)=PYpsi_save(:,at)
                        endif
                enddo
                deallocate(PYpsi_save)
        endif
        call cpu_time(time_2)
        !print *, parallel%myid, 'PYpsi time: ',time_2-time_1

end subroutine 

subroutine Calculate_deriv_Projector_overlaps(psi, PYpsi, atoms, elements, grid, parallel) 
        use lapack, only : dgemv
        use fft, only: recip_to_real, real_to_recip
        use  parallel_mod, only: parallel_task, parallel_wait
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct, Real_local
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        use odp_type, only : field_struct

        type(field_struct),intent(in) :: psi
        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(in) :: grid
        type(element_struct), intent(inout), target :: elements(:)
        type(parallel_struct), intent(in) :: parallel
        complex(dp),intent(inout) :: PYpsi(:,:,:)
        complex(dp), allocatable :: factor(:)
        type(atom_struct), pointer :: atom
        type(element_struct), pointer :: element
        real(dp), allocatable :: proj_packed(:,:,:)

        integer :: ii, ix, iy, iz, i_sphere, at, i_p, dir
        do at=1, size(atoms)
                atom=>atoms(at)
                element=>elements(atom%element)
                if(element%PA_type.ne.Real_local) cycle
                PYpsi(:,:,at)=0.0_dp
                if(atom%RL%ns.lt.1) cycle
                allocate(factor(atom%RL%ns))
                allocate(proj_packed(atom%RL%ns,3,element%n_proj))

                if((any(abs(grid%k).gt.tiny(1.0_dp)).or.any(abs(grid%A).gt.tiny(1.0_dp)))) then
                     factor(:)=exp(i_*((grid%k(1)+grid%A(1))*atom%RL%Rs(1,:atom%RL%ns) + &
                                        (grid%k(2)+grid%A(2))*atom%RL%Rs(2,:atom%RL%ns) + &
                                        (grid%k(3)+grid%A(3))*atom%RL%Rs(3,:atom%RL%ns)))
                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii)
                                iy=atom%RL%map(2,ii)
                                iz=atom%RL%map(3,ii)
                                factor(ii)=factor(ii)*psi%R(ix,iy,iz)
                        enddo

                else
                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii)
                                iy=atom%RL%map(2,ii)
                                iz=atom%RL%map(3,ii)
                                factor(ii)=psi%R(ix,iy,iz)
                        enddo
                endif                        
                do ii= 1, atom%RL%ns
                        i_sphere=atom%RL%map(4,ii)
                        proj_packed(ii,:,:)=atom%RL%deriv_projector(i_sphere,:,:)
                enddo

                do i_p=1, element%n_proj; do dir=1,3
                                PYpsi(dir,i_p,at)=sum(proj_packed(:,dir,i_p)*factor(:))                                               
                enddo; enddo
                deallocate(factor)
                deallocate(proj_packed)       
                ! print *, parallel%myid, 'PYpsi time: ',time_2-time_1
        enddo
        PYpsi(:,:,:)=PYpsi(:,:,:)*product(grid%box_length(:)/grid%Nr(:))


        !since threads only saw their own local contributions we need to sum (we
        !wait untill here so that the threads can work through all the atoms
        !separatly and hopefully take ~the same time to finish (else for any one
        !atom most threads don't contibute
        call parallel_task('sum', PYpsi, parallel, 'band')   

end subroutine 

subroutine Calculate_deriv_Projector_R_overlaps(psi, PYpsi, atoms, elements, grid, parallel) 
        use lapack, only : dgemv
        use fft, only: recip_to_real, real_to_recip
        use  parallel_mod, only: parallel_task, parallel_wait
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct, Real_local
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        use odp_type, only : field_struct

        type(field_struct),intent(in) :: psi
        type(atom_struct), intent(inout), target :: atoms(:)
        type(grid_struct), intent(in) :: grid
        type(element_struct), intent(inout), target :: elements(:)
        type(parallel_struct), intent(in) :: parallel
        complex(dp),intent(inout) :: PYpsi(:,:,:,:)
        complex(dp), allocatable :: factor(:), Rfactor(:,:)
        real(dp), allocatable :: proj_packed(:,:,:)
        type(atom_struct), pointer :: atom
        type(element_struct), pointer :: element
        integer :: ii, ix, iy, iz, i_sphere, at, i_p, dir, dir2
        do at=1, size(atoms)
                atom=>atoms(at)
                element=>elements(atom%element)
                if(element%PA_type.ne.Real_local) cycle
                PYpsi(:,:,:,at)=0.0_dp
                if(atom%RL%ns.lt.1) cycle
                allocate(factor(atom%RL%ns))
                allocate(Rfactor(atom%RL%ns,3))
                allocate(proj_packed(atom%RL%ns,3,element%n_proj))
                if((any(abs(grid%k).gt.tiny(1.0_dp)).or.any(abs(grid%A).gt.tiny(1.0_dp)))) then
                        factor(:)=exp(i_*((grid%k(1)+grid%A(1))*atom%RL%Rs(1,:atom%RL%ns) + &
                                        (grid%k(2)+grid%A(2))*atom%RL%Rs(2,:atom%RL%ns) + &
                                        (grid%k(3)+grid%A(3))*atom%RL%Rs(3,:atom%RL%ns)))
                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii)
                                iy=atom%RL%map(2,ii)
                                iz=atom%RL%map(3,ii)
                                factor(ii)=factor(ii)*psi%R(ix,iy,iz)
                        enddo

                else
                        do ii= 1, atom%RL%ns
                                ix=atom%RL%map(1,ii)
                                iy=atom%RL%map(2,ii)
                                iz=atom%RL%map(3,ii)
                                factor(ii)=psi%R(ix,iy,iz)
                        enddo
                endif
                do dir= 1, 3
                        Rfactor(:,dir)=atom%RL%Rs(dir,:atom%RL%ns)*factor(:)
                enddo
                do ii= 1, atom%RL%ns
                        i_sphere=atom%RL%map(4,ii)
                        proj_packed(ii,:,:)=atom%RL%deriv_projector(i_sphere,:,:)
                enddo        

                do i_p=1, element%n_proj; do dir=1,3; do dir2=1,3
                        PYpsi(dir,dir2,i_p,at)=sum(proj_packed(:,dir,i_p)*Rfactor(:,dir2))
                enddo; enddo; enddo

                deallocate(factor)
                deallocate(Rfactor)        
                deallocate(proj_packed)       
        enddo
        PYpsi(:,:,:,:)=PYpsi(:,:,:,:)*product(grid%box_length(:)/grid%Nr(:))

        !since threads only saw their own local contributions we need to sum (we
        !wait untill here so that the threads can work through all the atoms
        !separatly and hopefully take ~the same time to finish (else for any one
        !atom most threads don't contibute
        call parallel_task('sum', PYpsi, parallel, 'band')   

end subroutine 

subroutine allocate_and_fill_maps(grids, elements)
        use  element_type, only : element_struct, Real_local, PAW
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        use  parallel_mod, only : parallel_wait

        use splines, only: spline3
        use shared_type, only : allocate_shared, shared_fence

        type(grid_struct), intent(in) :: grids(:)
        type(element_struct), intent(inout), target :: elements(:)
        real(dp) :: X_low(3), X_high(3)
        integer i,j, dir

        do i=1, size(elements)
                if(elements(i)%PA_type.ne.Real_local) cycle
                if(allocated(elements(i)%RL%coarse_grid)) cycle

                X_low(:)=floor(-elements(i)%RL%radius/grids(2)%dR(:))*grids(2)%dR(:)
                X_high(:)=X_low(:) + (elements(i)%RL%max_sphere(:)-1)*grids(2)%dR(:)
                allocate(elements(i)%RL%coarse_grid(maxval(elements(i)%RL%max_sphere(:)),3))
                do dir=1,3
                        do j = 1, elements(i)%RL%max_sphere(dir)
                                elements(i)%RL%coarse_grid(j,dir)=(j-1)*grids(2)%dR(dir) + X_low(dir)
                                 !print *, 'coarse:', dir, j, elements(i)%RL%coarse_grid(j,dir)
                        enddo
                enddo

                elements(i)%RL%dense_dR(:)=(X_high(:)-X_low(:))/(elements(i)%RL%max_sphere_dense(:)-1)
                allocate(elements(i)%RL%dense_grid(maxval(elements(i)%RL%max_sphere_dense(:)),3))
                do dir=1,3
                        do j = 1, elements(i)%RL%max_sphere_dense(dir)
                                elements(i)%RL%dense_grid(j,dir)=(j-1)*elements(i)%RL%dense_dR(dir) + X_low(dir)
                                ! print *, 'dense:', dir, j, elements(i)%RL%dense_grid(j,dir)
                        enddo
                enddo
        enddo
        do i=1, size(elements)
                if(elements(i)%PP_type.ne.PAW) cycle

                if(allocated(elements(i)%PAW%fine_grid)) cycle
                X_low(:)=floor(-elements(i)%PAW%radius/grids(1)%dR(:))*grids(1)%dR(:)
                X_high(:)=X_low(:) + (elements(i)%PAW%max_sphere_fine(:)-1)*grids(1)%dR(:)
                allocate(elements(i)%PAW%fine_grid(maxval(elements(i)%PAW%max_sphere_fine(:)),3))
                do dir=1,3
                        do j = 1, elements(i)%PAW%max_sphere_fine(dir)
                                elements(i)%PAW%fine_grid(j,dir)=(j-1)*grids(1)%dR(dir) + X_low(dir)
                           !    print *, 'fine:', dir, j, elements(i)%PAW%fine_grid(j,dir)
                        enddo
                enddo
        enddo


end subroutine

subroutine allocate_and_fill_Ono_Hirose_matrix(grids, elements, parallel)
        use  element_type, only : element_struct, Real_local
        use  grids_type, only : grid_struct
        use  parallel_type, only : parallel_struct
        use  parallel_mod, only : parallel_wait

        use splines, only: spline3
        use shared_type, only : allocate_shared, shared_fence

        type(grid_struct), intent(in) :: grids(:)
        type(element_struct), intent(inout), target :: elements(:)
        type(parallel_struct), intent(in) :: parallel
        integer, allocatable :: arrayshape(:)
        real(dp), allocatable :: deltas(:), Beta_dir(:,:)
        integer i,j, dir

        do i=1, size(elements)
                if(elements(i)%PA_type.ne.Real_local) cycle

                !Beta matrix for dense temporary to coarse grid
                allocate(arrayshape(2))

                arrayshape=(/elements(i)%RL%max_sphere_dense(1),elements(i)%RL%max_sphere(1)/)
                call allocate_shared(elements(i)%RL%Beta_x, &
                arrayshape, elements(i)%RL%beta_x_shared, parallel)
                
                arrayshape=(/elements(i)%RL%max_sphere_dense(2),elements(i)%RL%max_sphere(2)/)
                call allocate_shared(elements(i)%RL%Beta_y, &
                arrayshape, elements(i)%RL%beta_y_shared, parallel)

                arrayshape=(/elements(i)%RL%max_sphere_dense(3),elements(i)%RL%max_sphere(3)/)
                call allocate_shared(elements(i)%RL%Beta_z, &
                arrayshape, elements(i)%RL%beta_z_shared, parallel)
                deallocate(arrayshape)

                call shared_fence(elements(i)%RL%beta_x_shared, parallel)
                call shared_fence(elements(i)%RL%beta_y_shared, parallel)
                call shared_fence(elements(i)%RL%beta_z_shared, parallel)
                if(parallel%myid_node.eq.0) then
                do dir=1,3
                        allocate(Beta_dir(elements(i)%RL%max_sphere_dense(dir), elements(i)%RL%max_sphere(dir)))
                        allocate(deltas(elements(i)%RL%max_sphere(dir)))
                        do j = 1, elements(i)%RL%max_sphere(dir)
                                deltas=0.0_dp
                                deltas(j)=1.0_dp
                                Beta_dir(:,j)=spline3(elements(i)%RL%coarse_grid(:elements(i)%RL%max_sphere(dir),dir),&
                                                deltas, elements(i)%RL%dense_grid(:elements(i)%RL%max_sphere_dense(dir),dir)) &
                                                *elements(i)%RL%dense_dR(dir)/grids(2)%dR(dir)
                        enddo


                        if(dir.eq.1) elements(i)%RL%Beta_x(:,:)=Beta_dir(:,:)
                        if(dir.eq.2) elements(i)%RL%Beta_y(:,:)=Beta_dir(:,:)
                        if(dir.eq.3) elements(i)%RL%Beta_z(:,:)=Beta_dir(:,:)

                        deallocate(Beta_dir)
                        deallocate(deltas)
                enddo
                endif
                call shared_fence(elements(i)%RL%beta_x_shared, parallel)
                call shared_fence(elements(i)%RL%beta_y_shared, parallel)
                call shared_fence(elements(i)%RL%beta_z_shared, parallel)

        enddo

end subroutine

subroutine Ono_Hirose(Beta_x, Beta_y, Beta_z, max_coarse, max_fine, dense_projector, coarse_projector)
        real(dp), intent(in) :: Beta_x(:,:), Beta_y(:,:), Beta_z(:,:)
        integer, intent(in) :: max_coarse(3), max_fine(3)
        real(dp), intent(in) :: dense_projector(:,:,:)
        real(dp), intent(inout) :: coarse_projector(:,:,:)
        integer :: iz_coarse, iy_coarse, ix_coarse
        integer :: iz_fine, iy_fine
        real(dp) :: smooth_projector_x(max_fine(2),max_fine(3)), smooth_projector_xy(max_fine(3)) 

        
        do ix_coarse=1, max_coarse(1)
                do iy_fine=1, max_fine(2); do iz_fine=1, max_fine(3)
                        smooth_projector_x(iy_fine,iz_fine)=sum(dense_projector(:,iy_fine,iz_fine)*Beta_x(:,ix_coarse))
                enddo;enddo
                do iy_coarse=1, max_coarse(2)
                        do iz_fine=1, max_fine(3)
                                smooth_projector_xy(iz_fine)=sum(smooth_projector_x(:,iz_fine)*Beta_y(:,iy_coarse))
                        enddo
                        do iz_coarse=1, max_coarse(3)
                                coarse_projector(ix_coarse,iy_coarse,iz_coarse)= &
                                        sum(smooth_projector_xy(:)*Beta_z(:,iz_coarse))
                        enddo
                enddo
        enddo

end subroutine

subroutine Projector_S_transformation(dR, ns, np, projector_in, projector_out, sij, &
                lam_i, der_projector_in, der_projector_out)
        use lapack, only : dgemm, dsyrk
        use linalg, only : eigh

        integer, intent(in) :: ns, np
        real(dp), intent(in) ::  sij(np,np), dR(3)
        real(dp), intent(in) :: projector_in(ns,np)
        real(dp), intent(inout) :: projector_out(ns,np)
        real(dp), intent(inout)  :: lam_i(np)
        real(dp) :: tmp_ij(np,np), tmp2_ij(np,np), Uij(np,np), projector_tmp(ns,np) 
        real(dp), allocatable  :: der_projector_tmp(:,:,:)

        real(dp), optional, intent(in) :: der_projector_in(ns,np,3)
        real(dp), optional, intent(inout) :: der_projector_out(ns,np,3)
        integer :: i, dir

        if(present(der_projector_in)) allocate(der_projector_tmp(ns,np,3))
        tmp_ij=0.0_dp
        !tmp_ij=Lij=< p_i | p_j >
        call dsyrk('U','C', np, ns, 1.0_dp,&
                        projector_in,  ns, 0.0_dp,  tmp_ij, np)  
        tmp_ij=tmp_ij*product(dR)

        call eigh(tmp_ij,lam_i, Uij)

        ! |eta>=(|P>.dot.U)/sqrt[lam]
        call dgemm('N','N', ns, np, np, 1.0_dp, &
                projector_in, ns, Uij, np, &
                0.0_dp,  projector_tmp, ns)

        if(present(der_projector_in)) then
                do dir=1,3
                        call dgemm('N','N', ns, np, np, 1.0_dp, &
                        der_projector_in(:,:,dir),  &
                        ns, Uij, np, &
                        0.0_dp,  der_projector_tmp(:,:,dir), ns)
                enddo
        endif
        do i=1, np
                projector_tmp(:,i)=projector_tmp(:,i)/sqrt(lam_i(i))
        enddo

        if(present(der_projector_in)) then
                do dir=1,3; do i=1, np
                        der_projector_tmp(:,i,dir) = &
                                der_projector_tmp(:,i,dir)/sqrt(lam_i(i))
                enddo; enddo
        endif

        !tmp2_ij = s.dot.U !Non-diag

        call dgemm('N','N', np, np, np, 1.0_dp, &
                sij, np, Uij, np, &
                0.0_dp,  tmp2_ij, np)
                
        !tmp_ij = U^dagger .dot. tmp2_ij
        call dgemm('T','N', np, np, np, 1.0_dp, &
                Uij, np, tmp2_ij, np, &
                0.0_dp,  tmp_ij, np)
        
        do i=1, np
                tmp_ij(:,i)=tmp_ij(:,i)*sqrt(lam_i(:))*sqrt(lam_i(i))
        enddo
        call eigh(tmp_ij,lam_i, Uij)

        do i=1, np
                lam_i(i) =max(lam_i(i),-0.997_dp)
        enddo
        ! |nu>=(|eta>.dot.U)
        call dgemm('N','N', ns, np, np, 1.0_dp, &
                projector_tmp, ns, Uij, np, &
                0.0_dp,  projector_out, ns)

        if(present(der_projector_in)) then
                do dir=1,3
                        call dgemm('N','N', ns, np, np, 1.0_dp, &
                        der_projector_tmp(:,:,dir), ns, Uij, np, &
                        0.0_dp,  der_projector_out(:,:,dir), ns)
                enddo
        endif

end subroutine

subroutine map_atom_sphere(ns, map, element, atom, grid, sphere, max_sphere, Rs, radius)
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct
        use  grids_type, only : grid_struct
        type(atom_struct), intent(inout)  :: atom
        type(element_struct), intent(inout)  :: element
        type(grid_struct), intent(in)  :: grid

        integer, intent(inout) :: ns, map(:,:)
        integer,intent(in) ::  max_sphere(3)
        real(dp), intent(in) :: sphere(:,:), radius
        real(dp), intent(inout), optional :: Rs(:,:)
        integer :: i, ix, iy, iz
        real(dp)  :: X_dense(3), X_close(3), sphere_center(3)
        real(dp), allocatable :: R2(:)

        if(present(Rs)) allocate(R2(size(Rs,2)))

        X_close(:)=floor(atom%R(:)/grid%dR(:)+0.5_dp)*grid%dR(:)
        sphere_center(1)=sphere(max_sphere(1)/2+1,1)
        sphere_center(2)=sphere(max_sphere(2)/2+1,2)
        sphere_center(3)=sphere(max_sphere(3)/2+1,3)

        ns=0
        i=0
        do iz=1, max_sphere(3)
                X_dense(3)= sphere(iz,3) + X_close(3) - sphere_center(3)
        do iy=1, max_sphere(2)
                X_dense(2)=sphere(iy,2) + X_close(2) - sphere_center(2)
        do ix=1, max_sphere(1)
                X_dense(1)=sphere(ix,1) + X_close(1) - sphere_center(1)
                i=i+1
                if(sqrt(sum((X_dense(:)-atom%R(:))**2)).le.radius) then
                        ns=ns+1 
                        map(1,ns)=ix
                        map(2,ns)=iy
                        map(3,ns)=iz
                        map(4,ns)=i
                        if(present(Rs)) Rs(:,ns)=X_dense(:)-atom%R(:)
                        if(present(Rs)) R2(ns)=sum((X_dense(:)-atom%R(:))**2)
                endif
        enddo;enddo;enddo
        !Quicksort the map (needed for spline)
        if(present(Rs)) call Sort(R2(:),Rs(:,:),map(:,:),1,ns)
end subroutine

subroutine map_atom_to_local(ns, map, Rs, element, atom, grid, sphere, max_sphere, radius)
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct
        use  grids_type, only : grid_struct
        type(atom_struct), intent(inout)  :: atom
        type(element_struct), intent(inout)  :: element
        type(grid_struct), intent(in)  :: grid
        integer, intent(in) ::  max_sphere(3)
        real(dp), intent(in) :: sphere(:,:), radius
        integer, intent(inout) :: ns, map(:,:)
        real(dp), intent(inout) :: Rs(:,:)

        integer :: ns_global
        integer, allocatable :: map_global(:,:)
        integer :: i, i_sphere(3), i_local(3), i_global(3), dir
        real(dp)  ::  atom_R_tmp(3)
        real(dp), allocatable :: R2(:)

        allocate(map_global(4, product(max_sphere)))

        call map_atom_sphere(ns_global, map_global, element, atom, grid, sphere, max_sphere, radius=radius)
        ns=0
        do i=1, ns_global
                i_sphere(1:3)=map_global(1:3,i)
                !Global point = sphere point - sphere_center + atom_nearist point
                i_global(:)=i_sphere(:) + 1 + floor(atom%R(:)/grid%dR(:)+0.5_dp) - (max_sphere(:)/2 +1)
                atom_R_tmp(:)=atom%R(:)
                do dir=1,3
                        if(i_global(dir).le.0) then 
                                i_global(dir)= i_global(dir) + grid%Nr(dir)
                                atom_R_tmp(dir)=atom_R_tmp(dir)+grid%Box_Length(dir)
                        else if(i_global(dir).gt.grid%Nr(dir)) then
                                i_global(dir)=i_global(dir)-grid%Nr(dir)
                                atom_R_tmp(dir)=atom_R_tmp(dir)-grid%Box_Length(dir)
                        endif
                enddo
                !local to global and check if in proc's range
                i_local(:)=i_global(:)-grid%myxyz_r(:)*grid%Nr_local(:)
                if(all(i_local(:).gt.0).and.all(i_local(:).le.grid%Nr_local(:))) then
                        ns=ns+1
                        map(1:3,ns)=i_local(1:3)
                        map(4,ns)=i
                        map(5,ns)= &
                        i_local(1) + grid%Nr_local(1)*(i_local(2)-1+grid%Nr_local(2)*(i_local(3)-1))
                        Rs(:,ns)=grid%R(i_local(1),i_local(2),i_local(3),:)-atom_R_tmp(:)
                endif
        enddo
        deallocate(map_global)
        if(ns.eq.0) return

        allocate(R2(ns))
        R2(1:ns)=(Rs(1,1:ns)**2+Rs(2,1:ns)**2+Rs(3,1:ns)**2)

        call Sort(R2(:), Rs(:,:),map(:,:),1,ns)
end subroutine

subroutine PAW_projector_projector_overlaps(elements, atoms, grid, parallel, all_PAW)
        use  atom_type, only : atom_struct
        use element_type, only : element_struct, Real_local, Reciprocal
        use  grids_type, only : grid_struct
        use simulation_type, only : all_PAW_struct
        use parallel_type, only : parallel_struct
        use m_paral_atom, only : get_my_natom, get_my_atmtab
        use  parallel_mod, only: parallel_task, parallel_wait
        use linalg, only: inv
        use lapack, only : dgemm
        use grids_mod, only : allocate_local_fields_R, allocate_local_fields_G
        use operations_3D, only : integrate_3D_G

        type(parallel_struct), intent(in) :: parallel
        type(atom_struct), intent(inout),target  :: atoms(:)
        type(element_struct), intent(inout),target  :: elements(:)
        type(grid_struct), intent(in)  :: grid
        type(all_PAW_struct), intent(inout), target :: all_PAW
        type(atom_struct), pointer :: atom_i, atom_j
        type(element_struct), pointer :: element_i, element_j
        integer :: i,j
        integer :: ii, ix, iy, iz, i_sphere, si,ei, sj,ej,j_p, i_p, all_proj, il,jl
        real(dp) :: delR, tmpR(3)
        real(dp), allocatable ::  proj_grid(:,:,:,:), factor(:,:), proj_packed(:,:), PTP_block(:,:), &
                                  AV(:,:), tmp(:,:)
        complex(dp), allocatable :: factor_G(:,:,:)
        integer, allocatable :: atom_blocks(:,:)
        real(dp) :: norm
        if(all_PAW%N_PAW_atoms.lt.1) return
        all_PAW%proj_overlap=0.0_dp

        all_PAW%paral_atom=.true.
        all_PAW%my_atmtab=>null()
        call get_my_natom(parallel%comm_space,all_PAW%my_natom, all_PAW%N_PAW_atoms)
        call get_my_atmtab(parallel%comm_space,all_PAW%my_atmtab,all_PAW%my_atmtab_allocated, all_PAW%paral_atom, &
                            all_PAW%N_PAW_atoms)
        
        allocate(atom_blocks(all_PAW%N_PAW_atoms,all_PAW%N_PAW_atoms))
        atom_blocks=0

        call allocate_local_fields_G(factor_G,grid)

        ei=0
        do i=1, all_PAW%N_PAW_atoms
                atom_i=>atoms(all_PAW%all_at(i))
                element_i=>elements(atom_i%element)
                si=ei+1
                ei=ei+element_i%n_proj
                if(all_PAW%paral_atom) then
                        if(.not.all_PAW%my_atmtab_allocated) cycle
                        if(all(all_PAW%my_atmtab(:).ne.i)) cycle
                endif
                atom_blocks(i,i)=1
                if(element_i%PA_type.ne.Real_local) cycle
                if(atom_i%RL%ns.lt.1) cycle
                if(allocated(proj_grid)) deallocate(proj_grid)
                call allocate_local_fields_R(proj_grid,grid,element_i%n_proj)

                !Cast atom i to the grid
                proj_grid=0.0_dp
                do ii= 1, atom_i%RL%ns
                        ix=atom_i%RL%map(1,ii)
                        iy=atom_i%RL%map(2,ii)
                        iz=atom_i%RL%map(3,ii)
                        i_sphere=atom_i%RL%map(4,ii)
                        proj_grid(ix,iy,iz,:)=atom_i%RL%projector(i_sphere,:) 
                        do i_p=1, element_i%n_proj
                                do j_p=1, i_p
                                        all_PAW%proj_overlap((si+j_p-1),(si+i_p-1))= &
                                        all_PAW%proj_overlap((si+j_p-1),(si+i_p-1)) + &
                                                atom_i%RL%projector(i_sphere,i_p)*atom_i%RL%projector(i_sphere,j_p) &
                                                *product(grid%box_length(:)/grid%Nr(:))
                        enddo; enddo
                enddo
                do i_p=1, element_i%n_proj
                        do j_p=1, i_p
                                all_PAW%proj_overlap((si+i_p-1),(si+j_p-1))=all_PAW%proj_overlap((si+j_p-1),(si+i_p-1))

                enddo; enddo
                
                !Look for overlapping j's
                ej=0
                do j=1, i-1
                        atom_j=>atoms(all_PAW%all_at(j))
                        element_j=>elements(atom_j%element)
                        sj=ej+1
                        ej=ej+element_j%n_proj
                        if(element_j%PA_type.ne.Real_local) cycle
                        if(atom_j%RL%ns.lt.1) cycle
                        if(allocated(factor)) deallocate(factor)
                        allocate(factor(atom_j%RL%ns,element_i%n_proj))
                        if(allocated(proj_packed)) deallocate(proj_packed)
                        allocate(proj_packed(atom_j%RL%ns,element_j%n_proj))

                        tmpR=atom_i%R(:)-atom_j%R(:)
                        where(tmpR(:).gt.0.5*grid%Box_Length(:))
                                tmpR(:)=tmpR(:)-grid%Box_Length(:)
                        endwhere
                        where(tmpR(:).lt.-0.5*grid%Box_Length(:))
                                tmpR(:)=tmpR(:)+grid%Box_Length(:)
                        endwhere
                        delR=sqrt(sum(tmpR(:)**2))
                        if(delR.gt.(element_i%PAW%radius+element_j%PAW%radius)) cycle
                        ! print *,  'Atom :', all_PAW%all_at(i), 'And',  all_PAW%all_at(j), 'Overlap with distance', delR
                        atom_blocks(i,j)=1
                        atom_blocks(j,i)=1
                        factor=0.0_dp
                        do ii= 1, atom_j%RL%ns
                                ix=atom_j%RL%map(1,ii)
                                iy=atom_j%RL%map(2,ii)
                                iz=atom_j%RL%map(3,ii)
                                factor(ii,:)=proj_grid(ix,iy,iz,:)
                        enddo
                        do ii= 1, atom_j%RL%ns
                                i_sphere=atom_j%RL%map(4,ii)
                                proj_packed(ii,:)=atom_j%RL%projector(i_sphere,:)
                        enddo
                        do j_p=1, element_j%n_proj; do i_p=1, element_i%n_proj
                                all_PAW%proj_overlap(si -1+ i_p, &
                                                sj -1+ j_p)= &
                                sum(factor(:,i_p)*proj_packed(:,j_p))*product(grid%box_length(:)/grid%Nr(:))
                                all_PAW%proj_overlap(sj -1+ j_p, &
                                                si -1+ i_p) = &
                                all_PAW%proj_overlap(si -1+ i_p, &
                                                sj -1+ j_p)
                        enddo;enddo
                enddo
        enddo

        call parallel_task('sum', all_PAW%proj_overlap, parallel, 'band') 

        ei=0
        do i=1, all_PAW%N_PAW_atoms
                atom_i=>atoms(all_PAW%all_at(i))
                element_i=>elements(atom_i%element)
                si=ei+1
                ei=ei+element_i%n_proj
                if(all_PAW%paral_atom) then
                        if(.not.all_PAW%my_atmtab_allocated) cycle
                        if(all(all_PAW%my_atmtab(:).ne.i)) cycle
                endif
                atom_blocks(i,i)=1
                if(element_i%PA_type.eq.Real_local) then
                        if(atom_i%RL%ns.lt.1) cycle
                        ej=0
                        do j=1, i-1
                                atom_j=>atoms(all_PAW%all_at(j))
                                element_j=>elements(atom_j%element)
                                sj=ej+1
                                ej=ej+element_j%n_proj

                                if(element_j%PA_type.ne.Reciprocal) cycle
                                tmpR=atom_i%R(:)-atom_j%R(:)
                                where(tmpR(:).gt.0.5*grid%Box_Length(:))
                                        tmpR(:)=tmpR(:)-grid%Box_Length(:)
                                endwhere
                                where(tmpR(:).lt.-0.5*grid%Box_Length(:))
                                        tmpR(:)=tmpR(:)+grid%Box_Length(:)
                                endwhere
                                delR=sqrt(sum(tmpR(:)**2))

                                if(delR.gt.(element_i%PAW%radius+element_j%PAW%radius)) cycle
                        ! print *,  'Atom :', all_PAW%all_at(i), 'And',  all_PAW%all_at(j), 'Overlap with distance', delR
                                atom_blocks(i,j)=1
                                atom_blocks(j,i)=1
                                !<P_i|P_j>
                                factor_G(:,:,:)=exp(i_*(grid%G(:,:,:,1)*tmpR(1) +&
                                                         grid%G(:,:,:,2)*tmpR(2) +&
                                                         grid%G(:,:,:,3)*tmpR(3)  ))                            
                                do j_p=1, element_j%n_proj; do i_p=1, element_i%n_proj
                                        il=element_i%PAW%tab%indlmn(1,i_p)
                                        jl=element_j%PAW%tab%indlmn(1,j_p)
                                all_PAW%proj_overlap(si -1+ i_p, &
                                                     sj -1+ j_p)= &
                                         real(integrate_3D_G(i_**(il-jl)*factor_G* &
                                          element_i%GNL(2)%projector_G(:,:,:,i_p)* &
                                          element_j%GNL(2)%projector_G(:,:,:,j_p), grid, parallel))
                                all_PAW%proj_overlap(sj -1+ j_p, &
                                                     si -1+ i_p) = &
                                all_PAW%proj_overlap(si -1+ i_p, &
                                                     sj -1+ j_p)
                                enddo;enddo
                        enddo

                else if(element_i%PA_type.eq.Reciprocal) then

                        ej=0
                        do j=1, i
                                atom_j=>atoms(all_PAW%all_at(j))
                                element_j=>elements(atom_j%element)
                                sj=ej+1
                                ej=ej+element_j%n_proj


                                tmpR=atom_i%R(:)-atom_j%R(:)
                                where(tmpR(:).gt.0.5*grid%Box_Length(:))
                                        tmpR(:)=tmpR(:)-grid%Box_Length(:)
                                endwhere
                                where(tmpR(:).lt.-0.5*grid%Box_Length(:))
                                        tmpR(:)=tmpR(:)+grid%Box_Length(:)
                                endwhere
                                delR=sqrt(sum(tmpR(:)**2))

                                if(delR.gt.(element_i%PAW%radius+element_j%PAW%radius)) cycle
                        ! print *,  'Atom :', all_PAW%all_at(i), 'And',  all_PAW%all_at(j), 'Overlap with distance', delR
                                atom_blocks(i,j)=1
                                atom_blocks(j,i)=1
                                !<P_i|P_j>
                                factor_G(:,:,:)=exp(i_*(grid%G(:,:,:,1)*tmpR(1) +&
                                                         grid%G(:,:,:,2)*tmpR(2) +&
                                                         grid%G(:,:,:,3)*tmpR(3)  ))                               
                                do j_p=1, element_j%n_proj; do i_p=1, element_i%n_proj
                                        il=element_i%PAW%tab%indlmn(1,i_p)
                                        jl=element_j%PAW%tab%indlmn(1,j_p)
                                all_PAW%proj_overlap(si -1+ i_p, &
                                                     sj -1+ j_p)= &
                                         real(integrate_3D_G(i_**(il-jl)*factor_G* &
                                          element_i%GNL(2)%projector_G(:,:,:,i_p)* &
                                          element_j%GNL(2)%projector_G(:,:,:,j_p), grid, parallel)) !&
                                        !/product(grid%box_length(:))**2
                                all_PAW%proj_overlap(sj -1+ j_p, &
                                                     si -1+ i_p) = &
                                all_PAW%proj_overlap(si -1+ i_p, &
                                                     sj -1+ j_p)
                                enddo;enddo
                        enddo
                endif
        enddo
        if(allocated(proj_grid)) deallocate(proj_grid)
        if(allocated(factor_G)) deallocate(factor_G)

        !Sum over different band/k groups which worked on different atoms
        call parallel_task('sum', atom_blocks, parallel, 'space') 
        call parallel_task('sum', all_PAW%proj_overlap, parallel, 'space') 

        if(parallel%myid.eq.0) then
                !print *, 'Atom Connections'
                !do i=1, all_PAW%N_PAW_atoms
                !        print *, atom_blocks(i,:)
                !enddo
        endif
        all_PAW%S_inv_full=0.0_dp
        ei=0
        do i=1, all_PAW%N_PAW_atoms
                atom_i=>atoms(all_PAW%all_at(i))
                element_i=>elements(atom_i%element)
                si=ei+1
                ei=ei+element_i%n_proj 
                if(allocated(PTP_block)) deallocate(PTP_block)
                allocate(PTP_block(element_i%n_proj,element_i%n_proj))
                PTP_block=0.0_dp
                atom_i%PAW%Sinv_block=0.0_dp

                if(all_PAW%paral_atom) then
                        if(.not.all_PAW%my_atmtab_allocated) cycle
                        if(all(all_PAW%my_atmtab(:).ne.i)) cycle
                endif
                do i_p=1, element_i%n_proj
                        do j_p=1, i_p
                                PTP_block(i_p,j_p)=all_PAW%proj_overlap((si+j_p-1),(si+i_p-1))
                                PTP_block(j_p,i_p)=PTP_block(i_p,j_p)
                enddo; enddo
                atom_i%PAW%Sinv_block=inv(element_i%PAW%sij(:,:))
                PTP_block=PTP_block+atom_i%PAW%Sinv_block
                atom_i%PAW%Sinv_block=inv(PTP_block)
        enddo
        ei=0
        do i=1, all_PAW%N_PAW_atoms
                atom_i=>atoms(all_PAW%all_at(i))
                element_i=>elements(atom_i%element)
                call parallel_task('sum', atom_i%PAW%Sinv_block, parallel, 'space') 
                !Initial guess at Matrix inverse
                si=ei+1
                ei=ei+element_i%n_proj
                do i_p=1, element_i%n_proj
                        do j_p=1, element_i%n_proj
                                all_PAW%S_inv_full((si+j_p-1),(si+i_p-1))=atom_i%PAW%Sinv_block(j_p,i_p)
                enddo; enddo
        enddo

        if(parallel%myid.eq.0) then
                all_proj=size(all_PAW%proj_overlap,1)
              !  print *, 'Full Projector overlap matrix'
              !  do i_p=1,all_proj
              !          print *, all_PAW%proj_overlap(i_p,:i_p)
              !  enddo
              
                !Form Sij^-1 + PTP
                ei=0
                do i=1, all_PAW%N_PAW_atoms
                        atom_i=>atoms(all_PAW%all_at(i))
                        element_i=>elements(atom_i%element)
                        si=ei+1
                        ei=ei+element_i%n_proj 
                        !PTP block used to store Sij^-1
                        if(allocated(PTP_block)) deallocate(PTP_block)
                        allocate(PTP_block(element_i%n_proj,element_i%n_proj))
                        PTP_block=0.0_dp
                        PTP_block=inv(element_i%PAW%sij(:,:))
                        do i_p=1, element_i%n_proj
                                do j_p=1, element_i%n_proj
                                        all_PAW%proj_overlap((si+j_p-1),(si+i_p-1)) = &
                                        all_PAW%proj_overlap((si+j_p-1),(si+i_p-1)) + &
                                        PTP_block(j_p,i_p)
                        enddo; enddo
                enddo

                if(.false.) then !Just use Lapack
                        all_PAW%S_inv_full=inv(all_PAW%proj_overlap)
                else !Schulz Algorithm for iterative inversion
                        allocate(AV(all_proj,all_proj))
                        allocate(tmp(all_proj,all_proj))
                        call dgemm('N','N', all_proj, all_proj, all_proj, 1.0_dp, &
                        all_PAW%proj_overlap, all_proj, all_PAW%S_inv_full, all_proj, &
                        0.0_dp,  AV, all_proj)
                        norm=0.0_dp
                        do i_p=1,all_proj
                                AV(i_p,i_p)=AV(i_p,i_p)-1.0_dp
                                norm=sum(AV(:,i_p)**2)
                                AV(i_p,i_p)=AV(i_p,i_p)+1.0_dp
                        enddo
                        do while(norm.gt.1E-16_dp)
                                tmp=0.0_dp
                                do i=1, all_proj
                                        tmp(i,i)=2.0_dp
                                enddo
                                tmp=tmp-AV
                                       
                                call dgemm('N','N', all_proj, all_proj, all_proj, 1.0_dp, &
                                all_PAW%S_inv_full, all_proj, tmp, all_proj, &
                                0.0_dp,  AV, all_proj)
                                all_PAW%S_inv_full=AV
      
                                call dgemm('N','N', all_proj, all_proj, all_proj, 1.0_dp, &
                                all_PAW%proj_overlap, all_proj, all_PAW%S_inv_full, all_proj, &
                                0.0_dp,  AV, all_proj)
                                norm=0.0_dp
                                do i_p=1,all_proj
                                        AV(i_p,i_p)=AV(i_p,i_p)-1.0_dp
                                        norm=sum(AV(:,i_p)**2)
                                        AV(i_p,i_p)=AV(i_p,i_p)+1.0_dp
                                enddo
                        enddo
        
                        deallocate(AV)
                        deallocate(tmp)
                endif 
                all_PAW%S_inv_full=-all_PAW%S_inv_full
        endif
        call parallel_task('bcast', all_PAW%S_inv_full, parallel, 'all', root=0)

end subroutine

subroutine map_atom_sphere_real_grid(ns, map, grid, max_sphere, R, radius, Rs)
        use  atom_type, only : atom_struct
        use  element_type, only : element_struct
        use  grids_type, only : grid_struct
        type(grid_struct), intent(in)  :: grid

        integer, intent(inout) :: ns, map(:,:)
        integer, intent(in) :: max_sphere(3)
        real(dp), intent(in) :: R(3), radius
        real(dp), intent(inout) :: Rs(:,:)
        integer :: i, ix, iy, iz, ix_local, iy_local, iz_local, ix_global, iy_global, iz_global, &
                   i_atom(3), i_low(3), i_high(3), local_min(3), local_max(3)
        real(dp)  ::  atom_R_tmp(3)

        ns=0
        map=-1
        Rs=0.0_dp

        local_min(:)=grid%myxyz_r(:)*grid%Nr_local(:)+1
        local_max(:)=(grid%myxyz_r(:)+1)*grid%Nr_local(:)

        i_atom(:)= floor(R(:)/grid%dR(:)+0.5_dp)  
        i_low(:)=   i_atom(:) - max_sphere(:)/2 
        i_high(:)=  i_atom(:) + max_sphere(:)/2
                

        !Atoms cube is too high in some direction (wrap arround also doesn't fall in local space)
        if(Any(i_low(:)>local_max(:).and.(i_high(:)-grid%Nr(:))<local_min(:))) return
        !Atoms cube is too  low in some direction (wrap arround also doesn't fall in local space)
        if(Any(i_high(:)<local_min(:).and.(i_low(:)+grid%Nr(:))>local_max(:))) return

        i=0
        do iz=1, max_sphere(3)
                atom_R_tmp(3)=R(3)
                !Global point = sphere point - sphere_center + atom_nearist point
                iz_global= iz + 1 + i_atom(3) - (max_sphere(3)/2 + 1)
                if(iz_global.le.0) then 
                        iz_global=iz_global+grid%Nr(3)
                        atom_R_tmp(3)=atom_R_tmp(3)+grid%Box_Length(3)
                else if(iz_global.gt.grid%Nr(3)) then
                        iz_global=iz_global-grid%Nr(3)
                        atom_R_tmp(3)=atom_R_tmp(3)-grid%Box_Length(3)
                endif
                iz_local=iz_global-grid%myxyz_r(3)*grid%Nr_local(3)
                if( (iz_local.gt.grid%Nr_local(3)).or.(iz_local.le.0) ) then
                        i= i + max_sphere(2)*max_sphere(1)
                        cycle
                endif
                do iy=1, max_sphere(2)
                        iy_global=iy + 1 + i_atom(2) - (max_sphere(2)/2 + 1)
                        atom_R_tmp(2)=R(2)
                        if(iy_global.le.0) then 
                                iy_global=iy_global+grid%Nr(2)
                                atom_R_tmp(2)=atom_R_tmp(2)+grid%Box_Length(2)
                        else if(iy_global.gt.grid%Nr(2)) then
                                iy_global=iy_global-grid%Nr(2)
                                atom_R_tmp(2)=atom_R_tmp(2)-grid%Box_Length(2)
                        endif
                        iy_local=iy_global-grid%myxyz_r(2)*grid%Nr_local(2)
                        if( (iy_local.gt.grid%Nr_local(2)).or.(iy_local.le.0) ) then
                                i=i + max_sphere(1)
                                cycle
                        endif
                        do ix=1, max_sphere(1)
                                ix_global=ix + 1 + i_atom(1) - (max_sphere(1)/2 + 1)
                                atom_R_tmp(1)=R(1)
                                if(ix_global.le.0) then 
                                        ix_global=ix_global+grid%Nr(1)
                                        atom_R_tmp(1)=atom_R_tmp(1)+grid%Box_Length(1)
                                else if(ix_global.gt.grid%Nr(1)) then
                                        ix_global=ix_global-grid%Nr(1)
                                        atom_R_tmp(1)=atom_R_tmp(1)-grid%Box_Length(1)
                                endif 
                                ix_local=ix_global-grid%myxyz_r(1)*grid%Nr_local(1)
                                i=i+1
                                if( (ix_local.gt.grid%Nr_local(1)).or.(ix_local.le.0) ) then
                                        cycle
                                else
                                        if(sqrt(sum( (grid%R(ix_local,iy_local,iz_local,:)-atom_R_tmp(:))**2)) &
                                                                .lt.radius) then
                                        ns=ns+1 
                                        map(1,ns)=ix_local            
                                        map(2,ns)=iy_local
                                        map(3,ns)=iz_local
                                        !This i needs to be the same for this n_s as given by
                                        !map_atom_sphere for this atom on the coarse sphere 
                                        map(4,ns)=i
                                        map(5,ns)=ix_local + grid%Nr_local(1)*(iy_local-1+grid%Nr_local(2)*(iz_local-1))
                                        Rs(:,ns)= &
                                                grid%R(ix_local,iy_local,iz_local,:)-atom_R_tmp(:)
                                        endif
                                endif
                        enddo
                enddo
        enddo
end subroutine

recursive subroutine Sort(x, y1, y2, left, right)
implicit none
real(dp),intent(inout)  :: x(:), y1(:,:)
integer,intent(inout) :: y2(:,:)
real(dp)  :: pivot, temp_x, temp_y1(size(y1,1))
integer :: temp_y2(size(y2,1))
integer :: left , right
integer :: i, j

pivot = x( (left+right) / 2 )
i = left
j = right
do
        do while (x(i) < pivot)
                i=i+1
        end do
        do while (pivot < x(j))
                j=j-1
        end do
        if (i >= j) exit
        temp_x = x(i);  x(i) = x(j);  x(j) = temp_x
        temp_y1(:) = y1(:,i);  y1(:,i) = y1(:,j);  y1(:,j) = temp_y1(:)
        temp_y2(:) = y2(:,i);  y2(:,i) = y2(:,j);  y2(:,j) = temp_y2(:)

        i=i+1; j=j-1
end do
if (j+1 < right)  call Sort(x, y1, y2, j+1, right)
if (left < i-1) call Sort(x, y1, y2, left, i-1)
end subroutine Sort


end module
