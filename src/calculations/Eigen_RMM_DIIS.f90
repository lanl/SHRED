module Eigen_RMM_DIIS
    use types, only : dp

    implicit none

    public :: Eig_RMM_DIIS


    contains

    subroutine Eig_RMM_DIIS(orbitals, grids, parallel, diis, potentials, atoms, elements, &
                      all_eigs, all_PAW, n_band_groups)
        use grids_type, only : grid_struct
        use odp_type, only : field_struct, orbital_struct
        use simulation_type, only : potentials_struct, all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: recip_to_real
        use Apply_Hamiltonian, only: Apply_H
        use FFT2Matrix, only : field_FFT2Matrix

        use parallel_mod,only : parallel_task, parallel_wait
        use parallel_type, only : parallel_struct
        use numerics_type, only : diis_numerics

        use FFT2Matrix, only : orbitals_FFT2Matrix, orbitals_Matrix2FFT, field_Matrix2FFT
        use odp, only: allocate_field, deallocate_field, allocate_orbital, deallocate_orbital 
        use linalg, only : eigh
        use Eigen_LOBPCG, only : orthogonalize, compute_Ritz_orthobasis
        use Eigen_ChebFilter, only : compute_Ritz_nonorthobasis


        type(parallel_struct), intent(in) :: parallel
        type(diis_numerics), intent(in) :: diis
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        real(dp) , intent(inout) :: all_eigs(:)
        type(all_PAW_struct), intent(in) :: all_PAW

        integer, intent(in), optional :: n_band_groups
        integer :: nX

        if(present(n_band_groups)) then
            nX=n_band_groups
        else
            nX=parallel%n_band_groups
        endif

        if(all_PAW%N_PAW_atoms.gt.0) then 
            call Eig_RMM_DIIS_Gen(orbitals, grids, parallel, diis, potentials, atoms, elements, &
            all_eigs, nX)
        else
            call Eig_RMM_DIIS_NC(orbitals, grids, parallel, diis, potentials, atoms, elements, &
            all_eigs, nX)
        endif

    end subroutine 

    subroutine Eig_RMM_DIIS_NC(orbitals, grids, parallel, diis, potentials, atoms, elements, all_eigs, n_band_groups)
        use grids_type, only : grid_struct
        use odp_type, only : field_struct, orbital_struct
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: recip_to_real
        use Apply_Hamiltonian, only: Apply_H
        use FFT2Matrix, only : field_FFT2Matrix

        use parallel_mod,only : parallel_task, parallel_wait
        use parallel_type, only : parallel_struct
        use numerics_type, only : diis_numerics

        use FFT2Matrix, only : orbitals_FFT2Matrix, orbitals_Matrix2FFT, field_Matrix2FFT
        use odp, only: allocate_field, deallocate_field, allocate_orbital, deallocate_orbital 
        use operations_3D, only : integrate_3D_G
        use linalg, only : eigh, solve
        use Eigen_LOBPCG, only : compute_Ritz_orthobasis, orthogonalize
        use Eigen_ChebFilter, only : compute_Ritz_nonorthobasis
        use grids_mod, only : allocate_local_fields_G


        type(parallel_struct), intent(in) :: parallel
        type(diis_numerics), intent(in) :: diis
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        real(dp) , intent(inout) :: all_eigs(:)

        integer, intent(in) :: n_band_groups

        type(field_struct), allocatable :: phi(:,:), Res(:,:), Hphi(:), HRes(:), KRes(:)
        real(dp), allocatable :: Kpc(:,:,:), lam(:), lam_band, GP(:,:,:)
        complex(dp), allocatable :: Y(:,:),  AY(:,:), Ywork(:,:)
        complex(dp) :: OR(diis%n_steps+1,diis%n_steps+1), v_diis(diis%n_steps+1), &
                       alpha(diis%n_steps+1), olap
        complex(dp) :: a,b,c,d,e,f
        complex(dp) :: alpha_cg0, A1(2,2), M1(2,2), U1(2,2)
        real(dp) :: norm, l1(2)

        integer, allocatable :: block_n_spinor(:)
        integer :: band_block, nY, nX, nband_block, nG
        integer:: i, j, k
        integer :: orb_grid
        LOGICAL :: final_Ritz

        if(parallel%myid.eq.0) print *, 'Starting RMM DIIS Eigensolver'
        nband_block=size(orbitals)
        nX=n_band_groups

        orb_grid=orbitals(1)%of(1)%grid
        do i=1,nband_block
            if(orbitals(i)%of(1)%grid.ne.orb_grid) then
                print *, 'All orbitals passed into EgenXolver need to be on the same grid (same k-point)'
            endif
        enddo

        grid=>grids(orbitals(1)%of(1)%grid)
        if(grid%reduction_option.eq.2) then
            !if(parallel%myid.eq.0) print *, 'Removing some of the zeros in x and y directions'
            nG=product(grid%Ng_small)
        else
            nG=product(grid%Ng_local)
        endif

        if(mod(nG,nX).ne.0) nG=(nG/nX+1)*nX
        call allocate_local_fields_G(Kpc,grid)
        if(diis%PC_type.eq.1) then
            call allocate_local_fields_G(GP,grid)
            if(grid%ecutsm.gt.0.0_dp) then
                GP=grid%G2*grid%p*grid%cutwf
                Kpc=(27.0_dp+9.0_dp*GP+3.0_dp*GP**2+GP**3)
                Kpc=Kpc/(Kpc+GP**4)
            else
                Kpc=(27.0_dp+9.0_dp*grid%G2+3.0_dp*grid%G2**2+grid%G2**3)
                Kpc=Kpc/(Kpc+grid%G2**4)
            endif
            if(grid%ecutsm.gt.0.0_dp) then
               deallocate(GP)
            endif
        endif
        all_eigs(:)=0.0_dp
        !Move orbitals into Matrix Form for LAPACK
        !Keep them stored in the unused portion of YXP vector
        !deallocate the orbitals to save space

        nY=nX*nband_block
        allocate(Y(maxval(orbitals(:)%n_spinor)*nG/nX,nY))
        allocate(AY(size(Y,1),nY))
        allocate(lam(nY))
        allocate(block_n_spinor(nband_block))
        block_n_spinor(:)=orbitals(:)%n_spinor
        !Initial Rayleigh-Ritz (subspace diagonalization)==================================================
        call orbitals_FFT2Matrix(orbitals(:), Y(:,:), grids, parallel, n_band_groups, &
        grids(orb_grid)%reduction_option)
        nY=0
        do band_block=1, nband_block
            allocate(Hphi(block_n_spinor(band_block)))

            do i=1, block_n_spinor(band_block)
                call allocate_field(Hphi(i), grid, parallel)
            enddo
            call Apply_H(orbitals(band_block)%of(:), Hphi, grids, potentials, atoms, elements, parallel, calc_R=.true.)
            call field_FFT2Matrix(Hphi, AY(:,(nY+1):(nY+nX)), grid, parallel, nX, grid%reduction_option)

            deallocate(Hphi)
            nY=nY+nX
        enddo
        allocate(Ywork(size(Y,1),nY))
        call compute_Ritz_orthobasis(Y(:,1:nY), AY(:,1:nY),nY, Ywork(:,1:nY), lam,  parallel, rot_A=.true., gamma=grid%gamma)
        deallocate(Ywork)
        all_eigs(1:nY)=lam(1:nY)  
        deallocate(lam)
        do band_block=1, nband_block
            orbitals(band_block)%eig=all_eigs(orbitals(band_block)%band)
        enddo

        call orbitals_Matrix2FFT(orbitals(:nband_block), Y(:,1:nY), grids, parallel, n_band_groups, &
        grid%reduction_option)

        deallocate(Y)
        !Initial Rayleigh-Ritz Done=============================================================================

        nY=0
        do band_block=1, nband_block
            allocate(Hphi(block_n_spinor(band_block)))
            allocate(HRes(block_n_spinor(band_block)))
            allocate(KRes(block_n_spinor(band_block)))
            allocate(Res(block_n_spinor(band_block),  diis%n_steps+1))
            allocate(phi(block_n_spinor(band_block),  diis%n_steps+1))
            do i=1, block_n_spinor(band_block)
                call allocate_field(Hphi(i), grid, parallel)
                call allocate_field(HRes(i), grid, parallel)
                call allocate_field(KRes(i), grid, parallel)
            enddo
            do j=1,   diis%n_steps+1
                do i=1, block_n_spinor(band_block)
                call allocate_field(Res(i,j), grid, parallel)
                call allocate_field(phi(i,j), grid, parallel)
                enddo
            enddo

            phi%grid=orb_grid
            Hphi%grid=orb_grid
 
            Res%grid=orb_grid
            HRes%grid=orb_grid
            KRes%grid=orb_grid

            do i=1, block_n_spinor(band_block)
                phi(i,1)%G=orbitals(band_block)%of(i)%G
                phi(i,1)%R=orbitals(band_block)%of(i)%R
            enddo

            OR(:,:)=0.0_dp
            
            !Zeroth Step Find Alpha (Lambda) 
            lam_band=orbitals(band_block)%eig(1)
            !Already calculated Hphi during initial Rayleigh-Ritz
            call field_Matrix2FFT(Hphi, AY(:,(nY+1):(nY+nX)), grids(orb_grid), parallel, nX, grids(orb_grid)%reduction_option)
            do i=1, block_n_spinor(band_block)
                Res(i,1)%G=Hphi(i)%G-lam_band*phi(i,1)%G
                KRes(i)%G=Res(i,1)%G
            enddo

            call precondition_Res(diis%pc_type, Kpc, grids(orb_grid), parallel, KRes(:))
            call Apply_H(KRes(:), HRes, grids, potentials, atoms, elements, parallel, calc_R=.true.)

        
            a=0.0_dp;b=0.0_dp;c=0.0_dp;d=0.0_dp; OR(1,1)=0.0_dp
            do i=1, block_n_spinor(band_block)
                a=a+integrate_3D_G(conjg(KRes(i)%G)*Hphi(i)%G, grids(orb_grid), parallel)
                b=b+integrate_3D_G(conjg(KRes(i)%G)*HRes(i)%G, grids(orb_grid), parallel)
                c=c+integrate_3D_G(conjg(KRes(i)%G)*phi(i,1)%G, grids(orb_grid), parallel)
                d=d+integrate_3D_G(conjg(KRes(i)%G)*Res(i,1)%G, grids(orb_grid), parallel)
                OR(1,1)=OR(1,1)+real(integrate_3D_G(conjg(Res(i,1)%G)*Res(i,1)%G, grids(orb_grid), parallel))
            enddo

            if(.true.) then !via Rayleigh Quotient minimization (KRESSE WAY)
                A1(1,1)=lam_band
                A1(1,2)=a
                A1(2,1)=a
                A1(2,2)=b
                M1(1,1)=1.0_dp*block_n_spinor(band_block)
                M1(1,2)=c
                M1(2,1)=c
                M1(2,2)=d
                call eigh(A1, M1, l1, U1) 
              !  if(parallel%myid.eq.0) print *, U1(:,1)
              !  if(parallel%myid.eq.0) print *, U1(:,2)

                alpha_cg0=U1(2,1)/U1(1,1)
                alpha_cg0=max(abs(alpha_cg0),0.1_dp)*alpha_cg0/abs(alpha_cg0)
                alpha_cg0=min(abs(alpha_cg0),1.0_dp)*alpha_cg0/abs(alpha_cg0)
            else ! Abinit Way - Re{<R_0|(H - e_0S)} |K R_0>} / |(H - e_0S) |K R_0>|**2
                !Tmporary use of Res(2) space
                do i=1, block_n_spinor(band_block)
                    Res(i,2)%G=HRes(i)%G-lam_band*Res(i,1)%G
                enddo
                e=0.0_dp;f=0.0_dp
                do i=1, block_n_spinor(band_block)
                    e=e + integrate_3D_G(conjg(Res(i,1)%G)*Res(i,2)%G, grids(orb_grid), parallel)
                    f=f + integrate_3D_G(conjg(Res(i,2)%G)*Res(i,2)%G, grids(orb_grid), parallel)
                enddo
                alpha_cg0=-(real(e)/abs(f))
            endif 

           ! if(parallel%myid.eq.0)  print *, 'alpha_cg0:', real(alpha_cg0)
            do i=1, block_n_spinor(band_block)
                phi(i,2)%G= phi(i,1)%G +  KRes(i)%G*alpha_cg0
            enddo
            !S-Orthogonormalize
            olap=1.0_dp+alpha_cg0*conjg(c) !<phi_1|S|phi_2>
            norm=real(olap + c*conjg(alpha_cg0)) + abs(alpha_cg0)**2*abs(d)

            do i=1, block_n_spinor(band_block)
                phi(i,2)%G=phi(i,2)%G/sqrt(norm)
            enddo

            call Apply_H(phi(:,2), Hphi, grids, potentials, atoms, elements, parallel, calc_R=.true.)

            do i=1, block_n_spinor(band_block)
                !Hphi(i)%G=(1.0_dp-olap)*Hphi(i)%G+alpha_cg0*HRes(i)%G
                !Sphi(i)%G=(1.0_dp-olap)*Sphi(i)%G+alpha_cg0*SRes(i)%G
                !Hphi(i)%G=Hphi(i)%G/sqrt(norm)
                !Sphi(i)%G=Sphi(i)%G/sqrt(norm)
                call deallocate_field(HRes(i),grids(HRes(i)%grid))
            enddo
            deallocate(HRes)

            lam_band=0.0_dp
            do i=1, block_n_spinor(band_block)
                lam_band=lam_band + real(integrate_3D_G(conjg(phi(i,2)%G)*Hphi(i)%G, grids(orb_grid), parallel))
            enddo
            
           ! lam_band=(1.0_dp-e)*lam_band+(2-e)*alpha_cg0*a+b*alpha_cg0**2
           ! lam_band=lam_band/norm
           ! if(parallel%myid.eq.0) print *, 'Initial Redisual Norm:', OR(1,1)
            do j=2,  diis%n_steps
              !  if(parallel%myid.eq.0)print *, j, 'lam_band', lam_band
                do i=1, block_n_spinor(band_block)
                    Res(i,j)%G=Hphi(i)%G-lam_band*phi(i,j)%G
                enddo
                do k=1,j
                    OR(k,j)=0.0_dp
                    do i=1, block_n_spinor(band_block)
                        OR(k,j)=OR(k,j)+integrate_3D_G(conjg(Res(i,k)%G)*Res(i,j)%G, grids(orb_grid), parallel)
                    enddo
                    OR(j,k)=conjg(OR(k,j))
                enddo
                              
                if(abs(OR(j,j)).lt.max(diis%min_resnorm,0.3_dp*real(OR(1,1)))) then 
                     !End on Trial per suggestion of: https://doi.org/10.1103/PhysRevB.54.11169
                    call precondition_Res(diis%pc_type, Kpc, grids(orb_grid), parallel, Res(:,j))
                    do i=1, block_n_spinor(band_block)                        
                        orbitals(band_block)%of(i)%G=phi(i,j)%G+alpha_cg0*Res(i,j)%G
                    enddo
                    exit
                endif

                OR(:j,j+1)=1.0_dp
                OR(j+1,:j)=1.0_dp
                OR(j+1,j+1)=0.0_dp
                v_diis(:)=0.0_dp
                v_diis(j+1)=1.0_dp
                alpha(:(j+1))=solve(OR(:(j+1),:(j+1)), v_diis(:(j+1)))
                if(parallel%myid.eq.0) print *, 'mixing co:',  alpha(:j), 'sum', sum(alpha(:j)) 
                alpha(:j)= alpha(:j) + (1.0_dp-sum(alpha(:j)))/j
                if(parallel%myid.eq.0) print *, 'lag. mult:',  alpha(j+1)

                do i=1, block_n_spinor(band_block)
                    phi(i,j+1)%G=0.0_dp
                    Res(i,j+1)%G=0.0_dp
                    do k=1,j
                        !Phi_bar(j) & Res_bar(j)
                        phi(i,j+1)%G=phi(i,j+1)%G + alpha(k)*phi(i,k)%G
                        Res(i,j+1)%G=Res(i,j+1)%G + alpha(k)*Res(i,k)%G
                    enddo
                    KRes(i)%G=Res(i,j+1)%G
                enddo
                call precondition_Res(diis%pc_type, Kpc, grids(orb_grid), parallel, KRes(:))

                do i=1, block_n_spinor(band_block)
                    phi(i,j+1)%G=phi(i,j+1)%G+alpha_cg0*KRes(i)%G
                enddo
                norm=0.0
                do i=1, block_n_spinor(band_block)
                    norm=norm+real(integrate_3D_G(conjg(phi(i,j+1)%G)*phi(i,j+1)%G, grids(orb_grid), parallel))
                enddo
                do i=1, block_n_spinor(band_block)
                    phi(i,j+1)%G=phi(i,j+1)%G/sqrt(norm)
                enddo

                call Apply_H(phi(:,j+1), Hphi, grids, potentials, atoms, elements, parallel, calc_R=.true.)
                lam_band=0.0_dp
                do i=1, block_n_spinor(band_block)
                    lam_band=lam_band + real(integrate_3D_G(conjg(phi(i,j+1)%G)*Hphi(i)%G, grids(orb_grid), parallel))
                enddo

                if(j.eq.diis%n_steps) then
                    !End on Trial per suggestion of: https://doi.org/10.1103/PhysRevB.54.11169
                    do i=1, block_n_spinor(band_block)
                        KRes(i)%G=Hphi(i)%G-lam_band*phi(i,j+1)%G
                    enddo
                    call precondition_Res(diis%pc_type, Kpc, grids(orb_grid), parallel, KRes(:))
                    do i=1, block_n_spinor(band_block)
                        orbitals(band_block)%of(i)%G=phi(i,j+1)%G + alpha_cg0*KRes(i)%G
                    enddo
                    exit
                endif
            enddo
            do i = 1, size(Hphi)
                call deallocate_field(Hphi(i),grids(Hphi(i)%grid))
                call deallocate_field(KRes(i),grids(KRes(i)%grid))
            enddo

            do j=1,diis%n_steps+1
                do i = 1, size(phi,1)
                    call deallocate_field(phi(i,j),grids(phi(i,j)%grid))
                    call deallocate_field(Res(i,j),grids(Res(i,j)%grid))
                enddo
            enddo

            deallocate(Hphi)
            deallocate(phi)
            deallocate(Res)
            deallocate(KRes)

            nY=nY+nX
        enddo        
        deallocate(Kpc)

        do band_block=1, nband_block
            do i=1, size(orbitals(band_block)%of(:))
            orbitals(band_block)%of(i)%G=orbitals(band_block)%of(i)%G*grid%cutwf
            enddo
        enddo

        allocate(Y(maxval(orbitals(:)%n_spinor)*nG/nX,nY))
        call orbitals_FFT2Matrix(orbitals(:), Y(:,:), grids, parallel, n_band_groups, &
        grids(orb_grid)%reduction_option)
        call orthogonalize(Y(:,1:nY), nY, parallel, AY(:,1:nY), gamma=grids(orb_grid)%gamma)
        call orbitals_Matrix2FFT(orbitals(:nband_block), Y(:,1:nY), grids, parallel, n_band_groups, &
        grid%reduction_option)

        final_Ritz=.false.

        if(final_Ritz) then
            allocate(lam(nY))
             ! Rayleigh-Ritz (subspace diagonalization)==================================================
            nY=0
            do band_block=1, nband_block
                allocate(Hphi(block_n_spinor(band_block)))
                do i=1, block_n_spinor(band_block)
                    call allocate_field(Hphi(i), grid, parallel)
                    Hphi(i)%grid=orbitals(band_block)%of(i)%grid
                enddo
                call Apply_H(orbitals(band_block)%of(:), Hphi, grids, potentials, atoms, elements, parallel, calc_R=.true.)
                call field_FFT2Matrix(Hphi, AY(:,(nY+1):(nY+nX)), grid, parallel, nX, grid%reduction_option)
                do i=1, block_n_spinor(band_block)
                    call deallocate_field(Hphi(i),grids(Hphi(i)%grid))
                enddo
                deallocate(Hphi)
                nY=nY+nX
            enddo

            allocate(Ywork(size(Y,1),nY))
            call compute_Ritz_orthobasis(Y(:,1:nY), AY(:,1:nY),nY, Ywork(:,1:nY), lam,  parallel, rot_A=.false., gamma=grid%gamma)
            deallocate(Ywork)
            all_eigs(1:nY)=lam(1:nY)  
            deallocate(lam)
            do band_block=1, nband_block
                orbitals(band_block)%eig=all_eigs(orbitals(band_block)%band)
            enddo
            call orbitals_Matrix2FFT(orbitals(:nband_block), Y(:,1:nY), grids, parallel, n_band_groups, &
            grid%reduction_option)
            ! Rayleigh-Ritz Done=============================================================================
        else
            nY=0
            all_eigs(:)=0.0_dp
            do band_block=1, nband_block
                allocate(Hphi(block_n_spinor(band_block)))
                do i=1, block_n_spinor(band_block)
                    call allocate_field(Hphi(i), grid, parallel)
                enddo
                Hphi%grid=orb_grid
                call Apply_H(orbitals(band_block)%of(:), Hphi, grids, potentials, atoms, elements, parallel, calc_R=.true.)
                lam_band=0.0_dp
                do i=1, block_n_spinor(band_block)
                    lam_band=lam_band  &
                    + real(integrate_3D_G(conjg(orbitals(band_block)%of(i)%G)*Hphi(i)%G, grids(orb_grid), parallel))
                enddo 
                all_eigs(orbitals(band_block)%band)=lam_band
                orbitals(band_block)%eig(:)=all_eigs(orbitals(band_block)%band) 
                do i = 1, size(Hphi)
                    call deallocate_field(Hphi(i),grids(Hphi(i)%grid))
                enddo
                deallocate(Hphi)     
                nY=nY+nX
            enddo
            call parallel_task('sum', all_eigs(:), parallel, 'diff_b')
        endif
        deallocate(AY)
        deallocate(Y)
    end subroutine

    subroutine Eig_RMM_DIIS_Gen(orbitals, grids, parallel, diis, potentials, atoms, elements, all_eigs, n_band_groups)
        use grids_type, only : grid_struct
        use odp_type, only : field_struct, orbital_struct
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: recip_to_real
        use Apply_Hamiltonian, only: Apply_H, Apply_S
        use FFT2Matrix, only : field_FFT2Matrix

        use parallel_mod,only : parallel_task, parallel_wait
        use parallel_type, only : parallel_struct
        use numerics_type, only : diis_numerics

        use FFT2Matrix, only : orbitals_FFT2Matrix, orbitals_Matrix2FFT, field_Matrix2FFT
        use odp, only: allocate_field, deallocate_field, allocate_orbital, deallocate_orbital 
        use operations_3D, only : integrate_3D_G
        use linalg, only : eigh, solve
        use Eigen_LOBPCG, only : compute_Ritz_orthobasis, B_orthogonalize
        use Eigen_ChebFilter, only : compute_Ritz_nonorthobasis
        use Non_Local_ion, only: Calculate_Projector_overlaps
        use grids_mod, only : allocate_local_fields_G

        type(parallel_struct), intent(in) :: parallel
        type(diis_numerics), intent(in) :: diis
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        real(dp) , intent(inout) :: all_eigs(:)

        integer, intent(in) :: n_band_groups

        type(field_struct), allocatable :: phi(:,:), Res(:,:), Hphi(:), HRes(:), Sphi(:), SRes(:), KRes(:)
        real(dp), allocatable :: Kpc(:,:,:), lam(:), lam_band, GP(:,:,:)
        complex(dp), allocatable :: Y(:,:),  AY(:,:), BY(:,:), Ywork(:,:), PYpsi(:,:,:)
        complex(dp) :: OR(diis%n_steps+1,diis%n_steps+1), v_diis(diis%n_steps+1), &
                       alpha(diis%n_steps+1), olap
        complex(dp) :: a,b,c,d,e,f
        complex(dp) :: alpha_cg0, A1(2,2), M1(2,2), U1(2,2)
        real(dp) :: norm, l1(2)

        integer, allocatable :: block_n_spinor(:)
        integer :: band_block, nY, nX, nband_block, nG
        integer:: i, j, k,s
        integer :: orb_grid
        logical :: final_Ritz
        if(parallel%myid.eq.0) print *, 'Starting RMM DIIS Eigensolver'
        nband_block=size(orbitals)
        nX=n_band_groups

        orb_grid=orbitals(1)%of(1)%grid
        do i=1,nband_block
            if(orbitals(i)%of(1)%grid.ne.orb_grid) then
                print *, 'All orbitals passed into EgenXolver need to be on the same grid (same k-point)'
            endif
        enddo

        grid=>grids(orbitals(1)%of(1)%grid)
        if(grid%reduction_option.eq.2) then
            !if(parallel%myid.eq.0) print *, 'Removing some of the zeros in x and y directions'
            nG=product(grid%Ng_small)
        else
            nG=product(grid%Ng_local)
        endif

        if(mod(nG,nX).ne.0) nG=(nG/nX+1)*nX
        call allocate_local_fields_G(Kpc,grid)
        if(diis%PC_type.eq.1) then
            call allocate_local_fields_G(GP,grid)
            if(grid%ecutsm.gt.0.0_dp) then
                GP=grid%G2*grid%p*grid%cutwf
                Kpc=(27.0_dp+9.0_dp*GP+3.0_dp*GP**2+GP**3)
                Kpc=Kpc/(Kpc+GP**4)
            else
                Kpc=(27.0_dp+9.0_dp*grid%G2+3.0_dp*grid%G2**2+grid%G2**3)
                Kpc=Kpc/(Kpc+grid%G2**4)
            endif
            if(grid%ecutsm.gt.0.0_dp) then
               deallocate(GP)
            endif
        endif
        all_eigs(:)=0.0_dp
        !Move orbitals into Matrix Form for LAPACK
        !Keep them stored in the unused portion of YXP vector
        !deallocate the orbitals to save space

        nY=nX*nband_block
        allocate(Y(maxval(orbitals(:)%n_spinor)*nG/nX,nY))
        allocate(AY(size(Y,1),nY))
        allocate(BY(size(Y,1),nY))
        allocate(lam(nY))
        allocate(block_n_spinor(nband_block))
        block_n_spinor(:)=orbitals(:)%n_spinor
        !Initial Rayleigh-Ritz (subspace diagonalization)==================================================

        call orbitals_FFT2Matrix(orbitals(:), Y(:,:), grids, parallel, n_band_groups, &
        grids(orb_grid)%reduction_option)
        nY=0
        do band_block=1, nband_block
            allocate(Sphi(block_n_spinor(band_block)))
            allocate(Hphi(block_n_spinor(band_block)))

            do i=1, block_n_spinor(band_block)
                call allocate_field(Sphi(i), grid, parallel)
                call allocate_field(Hphi(i), grid, parallel)
            enddo
            Sphi%grid=orb_grid
            Hphi%grid=orb_grid

            allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(orbitals(band_block)%of)))
            PYpsi=0.0_dp
            do s=1, size(orbitals(band_block)%of)
                call recip_to_real(orbitals(band_block)%of(s), grids)
                call Calculate_Projector_overlaps(orbitals(band_block)%of(s),  &
                PYpsi(:,:,s), atoms, elements, grid, parallel) 
             enddo
            call Apply_H(orbitals(band_block)%of(:), Hphi, grids, potentials, atoms, elements, parallel, &
                            calc_R=.false.,PYpsi_in=PYpsi)
            call Apply_S(orbitals(band_block)%of(:), Sphi, grids, atoms, elements, parallel, calc_R=.false.,PYpsi_in=PYpsi)
            call field_FFT2Matrix(Hphi, AY(:,(nY+1):(nY+nX)), grid, parallel, nX, grid%reduction_option)
            call field_FFT2Matrix(Sphi, BY(:,(nY+1):(nY+nX)), grid, parallel, nX, grid%reduction_option)

            deallocate(PYpsi)
            do i = 1, size(Sphi)
                call deallocate_field(Sphi(i),grids(Sphi(i)%grid))
                call deallocate_field(Hphi(i),grids(Hphi(i)%grid))
            enddo
            deallocate(Sphi)
            deallocate(Hphi)
            nY=nY+nX
        enddo
        !Make X (new) B_orthogonal to themselves, Rotate AY
        call B_orthogonalize(Y(:,1:nY), BY(:,1:nY), nY, parallel, AY(:,1:nY), gamma=grids(orb_grid)%gamma)
        allocate(Ywork(size(Y,1),nY))
        call compute_Ritz_orthobasis(Y(:,1:nY), AY(:,1:nY),nY, Ywork(:,1:nY), &
             lam,  parallel, rot_A=.true., BX=BY(:,1:nY), gamma=grids(orb_grid)%gamma)
        deallocate(Ywork)
        all_eigs(1:nY)=lam(1:nY)  
        deallocate(lam)
        do band_block=1, nband_block
            orbitals(band_block)%eig=all_eigs(orbitals(band_block)%band)
        enddo

        call orbitals_Matrix2FFT(orbitals(:nband_block), Y(:,1:nY), grids, parallel, n_band_groups, &
        grid%reduction_option)
        deallocate(Y)
        !Initial Rayleigh-Ritz Done=============================================================================

        nY=0
        do band_block=1, nband_block
            allocate(Hphi(block_n_spinor(band_block)))
            allocate(Sphi(block_n_spinor(band_block)))
            allocate(HRes(block_n_spinor(band_block)))
            allocate(KRes(block_n_spinor(band_block)))
            allocate(SRes(block_n_spinor(band_block)))
            allocate(Res(block_n_spinor(band_block),  diis%n_steps+1))
            allocate(phi(block_n_spinor(band_block),  diis%n_steps+1))
            allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(orbitals(band_block)%of)))

            do i=1, block_n_spinor(band_block)
                call allocate_field(Hphi(i), grid, parallel)
                call allocate_field(Sphi(i), grid, parallel)
                call allocate_field(HRes(i), grid, parallel)
                call allocate_field(SRes(i), grid, parallel)
                call allocate_field(KRes(i), grid, parallel)
            enddo
            do j=1,   diis%n_steps+1
                do i=1, block_n_spinor(band_block)
                call allocate_field(Res(i,j), grid, parallel)
                call allocate_field(phi(i,j), grid, parallel)
                enddo
            enddo

            phi%grid=orb_grid
            Hphi%grid=orb_grid
            Sphi%grid=orb_grid
 
            Res%grid=orb_grid
            HRes%grid=orb_grid
            SRes%grid=orb_grid
            KRes%grid=orb_grid

            do i=1, block_n_spinor(band_block)
                phi(i,1)%G=orbitals(band_block)%of(i)%G
                phi(i,1)%R=orbitals(band_block)%of(i)%R
            enddo

            OR(:,:)=0.0_dp
            
            !Zeroth Step Find Alpha (Lambda) 
            lam_band=orbitals(band_block)%eig(1)
            !Already calculated Hphi during initial Rayleigh-Ritz
            call field_Matrix2FFT(Hphi, AY(:,(nY+1):(nY+nX)), grids(orb_grid), parallel, nX, grids(orb_grid)%reduction_option)
            call field_Matrix2FFT(Sphi, BY(:,(nY+1):(nY+nX)), grids(orb_grid), parallel, nX, grids(orb_grid)%reduction_option)
            do i=1, block_n_spinor(band_block)
                Res(i,1)%G=Hphi(i)%G-lam_band*Sphi(i)%G
                KRes(i)%G=Res(i,1)%G
            enddo

            call precondition_Res(diis%pc_type, Kpc, grids(orb_grid), parallel, KRes(:))

            PYpsi=0.0_dp
            do s=1, size(KRes)
                call recip_to_real(KRes(s), grids)
                call Calculate_Projector_overlaps(KRes(s),  &
                PYpsi(:,:,s), atoms, elements, grid, parallel) 
            enddo
            call Apply_H(KRes(:), HRes, grids, potentials, atoms, elements, parallel, &
                        calc_R=.false., PYpsi_in=PYpsi)
            call Apply_S(KRes(:), SRes, grids, atoms, elements, parallel, calc_R=.false.,PYpsi_in=PYpsi)

        
            a=0.0_dp;b=0.0_dp;c=0.0_dp;d=0.0_dp; OR(1,1)=0.0_dp
            do i=1, block_n_spinor(band_block)
                a=a+integrate_3D_G(conjg(KRes(i)%G)*Hphi(i)%G, grids(orb_grid), parallel)
                b=b+integrate_3D_G(conjg(KRes(i)%G)*HRes(i)%G, grids(orb_grid), parallel)
                c=c+integrate_3D_G(conjg(KRes(i)%G)*Sphi(i)%G, grids(orb_grid), parallel)
                d=d+integrate_3D_G(conjg(KRes(i)%G)*SRes(i)%G, grids(orb_grid), parallel)
                OR(1,1)=OR(1,1)+real(integrate_3D_G(conjg(Res(i,1)%G)*Res(i,1)%G, grids(orb_grid), parallel))
            enddo

            if(.true.) then !via Rayleigh Quotient minimization (KRESSE WAY)
                A1(1,1)=lam_band
                A1(1,2)=a
                A1(2,1)=a
                A1(2,2)=b
                M1(1,1)=1.0_dp*block_n_spinor(band_block)
                M1(1,2)=c
                M1(2,1)=c
                M1(2,2)=d
                call eigh(A1, M1, l1, U1) 
              !  if(parallel%myid.eq.0) print *, U1(:,1)
              !  if(parallel%myid.eq.0) print *, U1(:,2)

                alpha_cg0=U1(2,1)/U1(1,1)
                alpha_cg0=max(abs(alpha_cg0),0.1_dp)*alpha_cg0/abs(alpha_cg0)
                alpha_cg0=min(abs(alpha_cg0),1.0_dp)*alpha_cg0/abs(alpha_cg0)
            else ! Abinit Way - Re{<R_0|(H - e_0S)} |K R_0>} / |(H - e_0S) |K R_0>|**2
                !Tmporary use of Res(2) space
                do i=1, block_n_spinor(band_block)
                    Res(i,2)%G=HRes(i)%G-lam_band*SRes(i)%G
                enddo
                e=0.0_dp;f=0.0_dp
                do i=1, block_n_spinor(band_block)
                    e=e + integrate_3D_G(conjg(Res(i,1)%G)*Res(i,2)%G, grids(orb_grid), parallel)
                    f=f + integrate_3D_G(conjg(Res(i,2)%G)*Res(i,2)%G, grids(orb_grid), parallel)
                enddo
                alpha_cg0=-(real(e)/abs(f))
            endif 

           ! if(parallel%myid.eq.0)  print *, 'alpha_cg0:', real(alpha_cg0)
            do i=1, block_n_spinor(band_block)
                phi(i,2)%G= phi(i,1)%G +  KRes(i)%G*alpha_cg0
            enddo
            !S-Orthogonormalize
            olap=1.0_dp+alpha_cg0*conjg(c) !<phi_1|S|phi_2>
            norm=real(olap + c*conjg(alpha_cg0)) + abs(alpha_cg0)**2*abs(d)

            do i=1, block_n_spinor(band_block)
                phi(i,2)%G=phi(i,2)%G/sqrt(norm)
            enddo

            PYpsi=0.0_dp
            do s=1, size(phi(:,2))
                call recip_to_real(phi(s,2), grids)
                call Calculate_Projector_overlaps(phi(s,2),  &
                PYpsi(:,:,s), atoms, elements, grid, parallel) 
            enddo
            call Apply_H(phi(:,2), Hphi, grids, potentials, atoms, elements, parallel, &
                calc_R=.false., PYpsi_in=PYpsi)
            call Apply_S(phi(:,2), Sphi, grids, atoms, elements, parallel, calc_R=.false., PYpsi_in=PYpsi)

            do i=1, block_n_spinor(band_block)
                !Hphi(i)%G=(1.0_dp-olap)*Hphi(i)%G+alpha_cg0*HRes(i)%G
                !Sphi(i)%G=(1.0_dp-olap)*Sphi(i)%G+alpha_cg0*SRes(i)%G
                !Hphi(i)%G=Hphi(i)%G/sqrt(norm)
                !Sphi(i)%G=Sphi(i)%G/sqrt(norm)
                call deallocate_field(HRes(i),grids(HRes(i)%grid))
                call deallocate_field(SRes(i),grids(SRes(i)%grid))
            enddo
            deallocate(HRes)
            deallocate(SRes)

            lam_band=0.0_dp
            do i=1, block_n_spinor(band_block)
                lam_band=lam_band + real(integrate_3D_G(conjg(phi(i,2)%G)*Hphi(i)%G, grids(orb_grid), parallel))
            enddo
            
           ! lam_band=(1.0_dp-e)*lam_band+(2-e)*alpha_cg0*a+b*alpha_cg0**2
           ! lam_band=lam_band/norm
           ! if(parallel%myid.eq.0) print *, 'Initial Redisual Norm:', OR(1,1)
            do j=2,  diis%n_steps
                do i=1, block_n_spinor(band_block)
                    Res(i,j)%G=Hphi(i)%G-lam_band*Sphi(i)%G
                enddo
                do k=1,j
                    OR(k,j)=0.0_dp
                    do i=1, block_n_spinor(band_block)
                        OR(k,j)=OR(k,j)+integrate_3D_G(conjg(Res(i,k)%G)*Res(i,j)%G, grids(orb_grid), parallel)
                    enddo
                    OR(j,k)=conjg(OR(k,j))
                enddo

                              
                if(abs(OR(j,j)).lt.max(diis%min_resnorm,0.3_dp*real(OR(1,1)))) then 
                     !End on Trial per suggestion of: https://doi.org/10.1103/PhysRevB.54.11169
                    call precondition_Res(diis%pc_type, Kpc, grids(orb_grid), parallel, Res(:,j))
                    do i=1, block_n_spinor(band_block)                        
                        orbitals(band_block)%of(i)%G=phi(i,j)%G+alpha_cg0*Res(i,j)%G
                    enddo
                    PYpsi=0.0_dp
                    do s=1, size(orbitals(band_block)%of)
                        call recip_to_real(orbitals(band_block)%of(s), grids)
                        call Calculate_Projector_overlaps(orbitals(band_block)%of(s),  &
                        PYpsi(:,:,s), atoms, elements, grid, parallel) 
                    enddo
                    call Apply_H(orbitals(band_block)%of(:), Hphi, grids, potentials, atoms, elements, parallel, &
                                    calc_R=.false.,PYpsi_in=PYpsi)
                    call Apply_S(orbitals(band_block)%of(:), Sphi, grids, atoms, elements, parallel, calc_R=.false.,PYpsi_in=PYpsi)                    
                    call field_FFT2Matrix(Sphi, BY(:,(nY+1):(nY+nX)), grid, parallel, nX, grid%reduction_option)
                    call field_FFT2Matrix(Hphi, AY(:,(nY+1):(nY+nX)), grid, parallel, nX, grid%reduction_option)
                    exit
                endif

                OR(:j,j+1)=1.0_dp
                OR(j+1,:j)=1.0_dp
                OR(j+1,j+1)=0.0_dp
                v_diis(:)=0.0_dp
                v_diis(j+1)=1.0_dp
                alpha(:(j+1))=solve(OR(:(j+1),:(j+1)), v_diis(:(j+1)))
                !if(parallel%myid.eq.0) print *, 'mixing co:',  alpha(:j), 'sum', sum(alpha(:j)) 
                alpha(:j)= alpha(:j) + (1.0_dp-sum(alpha(:j)))/j
                !if(parallel%myid.eq.0) print *, 'lag. mult:',  alpha(j+1)

                do i=1, block_n_spinor(band_block)
                    phi(i,j+1)%G=0.0_dp
                    Res(i,j+1)%G=0.0_dp
                    do k=1,j
                        !Phi_bar(j) & Res_bar(j)
                        phi(i,j+1)%G=phi(i,j+1)%G + alpha(k)*phi(i,k)%G
                        Res(i,j+1)%G=Res(i,j+1)%G + alpha(k)*Res(i,k)%G
                    enddo
                    KRes(i)%G=Res(i,j+1)%G
                enddo
                call precondition_Res(diis%pc_type, Kpc, grids(orb_grid), parallel, KRes(:))
                do i=1, block_n_spinor(band_block)
                    phi(i,j+1)%G=phi(i,j+1)%G+alpha_cg0*KRes(i)%G
                enddo

                PYpsi=0.0_dp
                do s=1, size(phi(:,j+1))
                    call recip_to_real(phi(s,j+1), grids)
                    call Calculate_Projector_overlaps(phi(s,j+1),  &
                    PYpsi(:,:,s), atoms, elements, grid, parallel) 
                enddo
                call Apply_S(phi(:,j+1), Sphi, grids, atoms, elements, parallel, calc_R=.false., PYpsi_in=PYpsi)
                norm=0.0
                do i=1, block_n_spinor(band_block)
                    norm=norm+real(integrate_3D_G(conjg(phi(i,j+1)%G)*Sphi(i)%G, grids(orb_grid), parallel))
                enddo
                do i=1, block_n_spinor(band_block)
                    phi(i,j+1)%G=phi(i,j+1)%G/sqrt(norm)
                    phi(i,j+1)%R=phi(i,j+1)%R/sqrt(norm)
                    Sphi(i)%G=Sphi(i)%G/sqrt(norm)
                    Sphi(i)%R=Sphi(i)%R/sqrt(norm)
                    PYpsi(:,:,i)=PYpsi(:,:,i)/sqrt(norm)
                enddo
                call Apply_H(phi(:,j+1), Hphi, grids, potentials, atoms, elements, parallel, calc_R=.false., PYpsi_in=PYpsi)
                lam_band=0.0_dp
                do i=1, block_n_spinor(band_block)
                    lam_band=lam_band + real(integrate_3D_G(conjg(phi(i,j+1)%G)*Hphi(i)%G, grids(orb_grid), parallel))
                enddo

                if(j.eq.diis%n_steps) then
                    !End on Trial per suggestion of: https://doi.org/10.1103/PhysRevB.54.11169
                    do i=1, block_n_spinor(band_block)
                        KRes(i)%G=Hphi(i)%G-lam_band*Sphi(i)%G
                    enddo
                    call precondition_Res(diis%pc_type, Kpc, grids(orb_grid), parallel, KRes(:))
                    do i=1, block_n_spinor(band_block)
                        orbitals(band_block)%of(i)%G=phi(i,j+1)%G + alpha_cg0*KRes(i)%G
                    enddo
                    PYpsi=0.0_dp
                    do s=1, size(orbitals(band_block)%of)
                        call recip_to_real(orbitals(band_block)%of(s), grids)
                        call Calculate_Projector_overlaps(orbitals(band_block)%of(s),  &
                        PYpsi(:,:,s), atoms, elements, grid, parallel) 
                    enddo
                    call Apply_H(orbitals(band_block)%of(:), Hphi, grids, potentials, atoms, elements, parallel, &
                                    calc_R=.false.,PYpsi_in=PYpsi)
                    call Apply_S(orbitals(band_block)%of(:), Sphi, grids, atoms, elements, parallel, calc_R=.false.,PYpsi_in=PYpsi)                    
                    call field_FFT2Matrix(Sphi, BY(:,(nY+1):(nY+nX)), grid, parallel, nX, grid%reduction_option)
                    call field_FFT2Matrix(Hphi, AY(:,(nY+1):(nY+nX)), grid, parallel, nX, grid%reduction_option)
                    exit
                endif
            enddo
    
            do i = 1, size(Hphi)
                call deallocate_field(Hphi(i),grids(Hphi(i)%grid))
                call deallocate_field(Sphi(i),grids(Sphi(i)%grid))
                call deallocate_field(KRes(i),grids(KRes(i)%grid))
            enddo
            do j=1,diis%n_steps+1
                do i = 1, size(phi,1)
                    call deallocate_field(phi(i,j),grids(phi(i,j)%grid))
                    call deallocate_field(Res(i,j),grids(Res(i,j)%grid))
                enddo
            enddo
            deallocate(Hphi)
            deallocate(phi)
            deallocate(Sphi)
            deallocate(Res)
            deallocate(KRes)
            deallocate(PYpsi)

            nY=nY+nX
        enddo        
        deallocate(Kpc)

        do band_block=1, nband_block
            do i=1, size(orbitals(band_block)%of(:))
            orbitals(band_block)%of(i)%G=orbitals(band_block)%of(i)%G*grid%cutwf
            enddo
        enddo

        allocate(Y(maxval(orbitals(:)%n_spinor)*nG/nX,nY))
        call orbitals_FFT2Matrix(orbitals(:nband_block), Y(:,:), grids, parallel, n_band_groups, &
                                 grids(orb_grid)%reduction_option)
        call B_orthogonalize(Y(:,1:nY), BY(:,1:nY), nY, parallel, AY(:,1:nY), gamma=grids(orb_grid)%gamma)
        deallocate(BY)

        final_Ritz=.false.
        if(final_Ritz) then
            allocate(lam(nY))
             ! Rayleigh-Ritz (subspace diagonalization)==================================================
            allocate(Ywork(size(Y,1),nY))
            call compute_Ritz_orthobasis(Y(:,1:nY), AY(:,1:nY),nY, Ywork(:,1:nY), &
                lam,  parallel, rot_A=.false., gamma=grids(orb_grid)%gamma)
            deallocate(Ywork)
            all_eigs(1:nY)=lam(1:nY)  
            deallocate(lam)
            do band_block=1, nband_block
                orbitals(band_block)%eig=all_eigs(orbitals(band_block)%band)
            enddo
            call orbitals_Matrix2FFT(orbitals(:nband_block), Y(:,1:nY), grids, parallel, n_band_groups, &
            grid%reduction_option)
            ! Rayleigh-Ritz Done=============================================================================
        else
            do i=1, nY
                all_eigs(i)= real(sum(conjg(Y(:,i))*AY(:,i)))
            enddo
            call parallel_task('sum', all_eigs(:), parallel, 'k')

            do band_block=1, nband_block
                orbitals(band_block)%eig=all_eigs(orbitals(band_block)%band)
            enddo
        endif
        deallocate(AY)
        deallocate(Y)
    
    end subroutine

    subroutine precondition_Res(PC_type, Kpc, grid, parallel, Res)
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use odp_type, only : field_struct
        use operations_3D, only : integrate_3D_G

        type(parallel_struct), intent(in) :: parallel
        integer, intent(in) ::PC_type
        real(dp),intent(inout) :: Kpc(:,:,:)
        type(field_struct), intent(inout) :: Res(:)
        type(grid_struct), intent(inout) :: grid
        integer :: i
        real(dp):: Ekinx3
        real(dp), allocatable :: GP(:,:,:)

        if(PC_type.eq.2) then
            Ekinx3=0.0_dp
            if(grid%ecutsm.gt.0.0_dp) then
                allocate(GP(size(grid%G2,1),size(grid%G2,2),size(grid%G2,3)))
                GP=grid%G2*grid%p*grid%cutwf
                do i=1, size(Res)
                    Ekinx3 = Ekinx3 + &
                    1.5*real(integrate_3D_G(abs(Res(i)%G)**2*GP, grid, parallel))
                enddo
            else
                do i=1, size(Res)
                    Ekinx3 = Ekinx3 + &
                    1.5*real(integrate_3D_G(abs(Res(i)%G)**2*grid%G2*grid%cutwf, grid, parallel))
                enddo
            endif
            grid%G2=grid%G2/Ekinx3
            if(grid%ecutsm.gt.0.0_dp) then
                Kpc=27.0_dp+GP*(18.0_dp+GP*(12.0_dp+8*GP))
                Kpc=Kpc/(Kpc+16.0_dp*GP**4)
            else
                Kpc=27.0_dp+grid%G2*(18.0_dp+grid%G2*(12.0_dp+8*grid%G2))
                Kpc=Kpc/(Kpc+16.0_dp*grid%G2**4)
            endif
            if(grid%ecutsm.gt.0.0_dp) then
                deallocate(GP)
            endif
            grid%G2=grid%G2*Ekinx3
        endif
        do i=1, size(Res)
            Res(i)%G=Res(i)%G*Kpc
        enddo

    end subroutine

end module
