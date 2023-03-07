module Eigen_ChebFilter
    use types, only : dp

    implicit none

    public :: Eig_ChebyFilter_NC


    contains

    subroutine Eig_ChebyFilter(orbitals, grids, parallel, numerics, potentials, atoms, elements, all_eigs, all_PAW, n_band_groups)
        use grids_type, only : grid_struct
        use odp_type, only : field_struct, orbital_struct
        use simulation_type, only : potentials_struct, all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: recip_to_real

        use parallel_mod,only : parallel_task, parallel_wait
        use parallel_type, only : parallel_struct
        use numerics_type, only : numerics_struct

        use FFT2Matrix, only : orbitals_FFT2Matrix, orbitals_Matrix2FFT, grids_G2_FFT2Matrix, &
                                orbital_FFT2Matrix
        use odp, only: allocate_field, deallocate_field, allocate_orbital, deallocate_orbital 


        type(parallel_struct), intent(in) :: parallel
        type(numerics_struct), intent(in) :: numerics
        type(all_PAW_struct), intent(inout) :: all_PAW
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout) :: grids(:)
        real(dp) , intent(inout) :: all_eigs(:)
        integer, intent(in), optional :: n_band_groups
        integer :: nX

        if(present(n_band_groups)) then
            nX=n_band_groups
        else
            nX=parallel%n_band_groups
        endif

        if(all_PAW%N_PAW_atoms.gt.0) then 
            call Eig_ChebyFilter_Gen(orbitals, grids, parallel, numerics, potentials, atoms, elements, all_eigs, nX, all_PAW)
        else
            call Eig_ChebyFilter_NC(orbitals, grids, parallel, numerics, potentials, atoms, elements, all_eigs, nX)
        endif

    end subroutine 

    subroutine Eig_ChebyFilter_NC(orbitals, grids, parallel, numerics, potentials, atoms, elements, all_eigs, n_band_groups)
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
        use numerics_type, only : numerics_struct

        use FFT2Matrix, only : orbitals_FFT2Matrix, orbitals_Matrix2FFT, grids_G2_FFT2Matrix, &
                               orbital_FFT2Matrix
        use odp, only: allocate_field, deallocate_field, allocate_orbital, deallocate_orbital 
        use operations_3D, only : integrate_3D_G


        type(parallel_struct), intent(in) :: parallel
        type(numerics_struct), intent(in) :: numerics
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        real(dp) , intent(inout) :: all_eigs(:)

        integer, intent(in) :: n_band_groups

        type(field_struct), allocatable :: C0(:), C1(:), Hpsi(:)
        real(dp), allocatable :: lam(:)
        complex(dp), allocatable :: Y(:,:),  AY(:,:), Ywork(:,:)

        integer, allocatable :: block_n_spinor(:)
        integer :: band_block, nY, nX, nband_block, nG
        integer:: i, n
        integer :: orb_grid

        real(dp) :: c, r, lam_plus, lam_minus, sig, tau, a_l, sig_2

        if(parallel%myid.eq.0) print *, 'Starting Chebychev Filter EigenXolver'
        nband_block=size(orbitals)
        nX=n_band_groups

        orb_grid=orbitals(1)%of(1)%grid
        do i=1,nband_block
            if(orbitals(i)%of(1)%grid.ne.orb_grid) then
                print *, 'All orbitals passed into EgenXolver need to be on the same grid (same k-point)'
            endif
        enddo

        !Setup Preconditioner
        grid=>grids(orbitals(1)%of(1)%grid)
        if(grid%reduction_option.eq.2) then
            !if(parallel%myid.eq.0) print *, 'Removing some of the zeros in x and y directionX'
            nG=product(grid%Ng_small)
        else
            nG=product(grid%Ng_local)
        endif

        if(mod(nG,nX).ne.0) nG=(nG/nX+1)*nX

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
        
        lam=0.0_dp
        do band_block=1, nband_block
            allocate(Hpsi(block_n_spinor(band_block)))
            do i=1, block_n_spinor(band_block)
                call allocate_field(Hpsi(i), grid, parallel)
            enddo
            Hpsi%grid=orb_grid

            call Apply_H(orbitals(band_block)%of(:), Hpsi, grids, potentials, atoms, elements, parallel, calc_R=.true.)
            do i=1, block_n_spinor(band_block)
                lam(orbitals(band_block)%band)=lam(orbitals(band_block)%band)+ &
                    real(integrate_3D_G(conjg(orbitals(band_block)%of(i)%G)*Hpsi(i)%G, grids(orb_grid), parallel))
            enddo

            do i = 1, size(Hpsi)
                call deallocate_field(Hpsi(i),grids(Hpsi(i)%grid))
            enddo

            deallocate(Hpsi)
        enddo
        call parallel_task('sum', lam(:), parallel, 'diff_b')

        lam_plus=grid%Ecut
        lam_minus=maxval(lam)
        a_l=minval(lam)

        r=lam_plus-lam_minus
        r=r*0.5_dp
        c=lam_plus+lam_minus
        c=c*0.5_dp
        sig=r/(c-a_l)
        tau=2.0_dp/sig

        if(parallel%myid_k.eq.0) print *, 'Filter max:', lam_plus, 'Filter min', lam_minus
        nY=0
        do band_block=1, nband_block
            allocate(Hpsi(block_n_spinor(band_block)))
            allocate(C0(block_n_spinor(band_block)))
            allocate(C1(block_n_spinor(band_block)))

            do i=1, block_n_spinor(band_block)
                call allocate_field(Hpsi(i), grid, parallel)
                call allocate_field(C0(i), grid, parallel)
                call allocate_field(C1(i), grid, parallel)
            enddo
            Hpsi%grid=orb_grid
            C0%grid=orb_grid
            C1%grid=orb_grid

            

            do i=1, block_n_spinor(band_block)
                C0(i)%G=orbitals(band_block)%of(i)%G
            enddo
            call Apply_H(orbitals(band_block)%of(:), Hpsi, grids, potentials, atoms, elements, parallel, calc_R=.true.)

            do i=1, block_n_spinor(band_block)
                C1(i)%G=(Hpsi(i)%G-c*C0(i)%G)*sig/r
            enddo
            do n=2, numerics%cheby%inner_steps
                sig_2=1.0_dp/(tau-sig)

                call Apply_H(C1(:), Hpsi, grids, potentials, atoms, elements, parallel, calc_R=.true.)
                do i=1, block_n_spinor(band_block)
                    orbitals(band_block)%of(i)%G=2.0_dp*sig_2/r*(Hpsi(i)%G-c*C1(i)%G)-sig_2*sig*C0(i)%G
                    C0(i)%G=C1(i)%G
                    C1(i)%G=orbitals(band_block)%of(i)%G
                    sig=sig_2
                enddo
            enddo

            call Apply_H(orbitals(band_block)%of(:), Hpsi, grids, potentials, atoms, elements, parallel, calc_R=.true.)

            call field_FFT2Matrix(Hpsi, AY(:,(nY+1):(nY+nX)), grid, parallel, nX, grid%reduction_option)

            nY=nY+nX
            do i = 1, size(Hpsi)
                call deallocate_field(Hpsi(i),grids(Hpsi(i)%grid))
                call deallocate_field(C0(i),grids(C0(i)%grid))
                call deallocate_field(C1(i),grids(C1(i)%grid))
            enddo

            deallocate(C0)
            deallocate(Hpsi)
            deallocate(C1)
        enddo

        call orbitals_FFT2Matrix(orbitals(:), Y(:,:), grids, parallel, nX, &
                                 grids(orb_grid)%reduction_option)
    
        !The Big Ritz
        allocate(Ywork(size(Y,1),nY))
        call compute_Ritz_nonorthobasis(nY, lam, Y(:,1:nY), AY(:,1:nY), Ywork(:,1:nY), parallel, gamma=grids(orb_grid)%gamma)
        deallocate(Ywork)
        deallocate(AY)

        all_eigs(1:nY)=lam(1:nY)  

        deallocate(lam)

        call orbitals_Matrix2FFT(orbitals(:nband_block), Y(:,1:nY), grids, parallel, nX, &
                                 grid%reduction_option)
        deallocate(Y)

        do band_block=1, nband_block
            orbitals(band_block)%eig(:)=all_eigs(orbitals(band_block)%band)
            do i=1, size(orbitals(band_block)%of(:))
                orbitals(band_block)%of(i)%G=orbitals(band_block)%of(i)%G*grid%cutwf
                call recip_to_real(orbitals(band_block)%of(i), grids)
            enddo
        enddo

    end subroutine

    subroutine Eig_ChebyFilter_Gen(orbitals, grids, parallel, numerics, potentials, atoms, elements, &
                                     all_eigs, n_band_groups, all_PAW)
        use grids_type, only : grid_struct
        use odp_type, only : field_struct, orbital_struct
        use simulation_type, only : potentials_struct, all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: recip_to_real, real_to_recip
        use Apply_Hamiltonian, only: Apply_SinvH, Apply_H, Apply_S
        use Non_Local_ion, only : Apply_S_power,Calculate_Projector_overlaps
        use FFT2Matrix, only : field_FFT2Matrix
        use Eigen_LOBPCG, only: compute_Ritz_orthobasis, B_orthogonalize

        use parallel_mod,only : parallel_task, parallel_wait
        use parallel_type, only : parallel_struct
        use numerics_type, only : numerics_struct

        use FFT2Matrix, only : orbitals_FFT2Matrix, orbitals_Matrix2FFT, grids_G2_FFT2Matrix, &
                               orbital_FFT2Matrix
        use odp, only: allocate_field, deallocate_field, allocate_orbital, deallocate_orbital 
        use operations_3D, only : integrate_3D_R


        type(parallel_struct), intent(in) :: parallel
        type(numerics_struct), intent(in) :: numerics
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        real(dp) , intent(inout) :: all_eigs(:)
        type(all_PAW_struct), intent(inout) :: all_PAW

        integer, intent(in) :: n_band_groups

        type(field_struct), allocatable :: C0(:), C1(:), Hpsi(:)
        real(dp), allocatable :: lam(:)
        complex(dp), allocatable :: Y(:,:),  AY(:,:), BY(:,:),Ywork(:,:), PYpsi(:,:,:)

        integer, allocatable :: block_n_spinor(:)
        integer :: band_block, nY, nX, nband_block, nG
        integer:: i, n
        integer :: orb_grid

        real(dp) :: c, r, lam_plus, lam_minus, sig, tau, a_l, sig_2

        if(parallel%myid.eq.0) print *, 'Starting Chebychev Filter EigenXolver'
        nband_block=size(orbitals)
        nX=n_band_groups

        orb_grid=orbitals(1)%of(1)%grid
        do i=1,nband_block
            if(orbitals(i)%of(1)%grid.ne.orb_grid) then
                print *, 'All orbitals passed into EgenXolver need to be on the same grid (same k-point)'
            endif
        enddo

        !Setup Preconditioner
        grid=>grids(orbitals(1)%of(1)%grid)
        if(grid%reduction_option.eq.2) then
            !if(parallel%myid.eq.0) print *, 'Removing some of the zeros in x and y directionX'
            nG=product(grid%Ng_small)
        else
            nG=product(grid%Ng_local)
        endif

        if(mod(nG,nX).ne.0) nG=(nG/nX+1)*nX

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
        
        lam=0.0_dp
        do band_block=1, nband_block
            allocate(Hpsi(block_n_spinor(band_block)))
            do i=1, block_n_spinor(band_block)
                call allocate_field(Hpsi(i), grid, parallel)
            enddo
            Hpsi%grid=orb_grid

            call Apply_SinvH(orbitals(band_block)%of(:), Hpsi, grids, potentials, atoms, elements, parallel, all_PAW, calc_G=.true.)

            do i=1, block_n_spinor(band_block)
                lam(orbitals(band_block)%band)=lam(orbitals(band_block)%band)+ &
                    real(integrate_3D_R(conjg(orbitals(band_block)%of(i)%R)*Hpsi(i)%R, grids(orb_grid), parallel))
            enddo

            do i = 1, size(Hpsi)
                call deallocate_field(Hpsi(i),grids(Hpsi(i)%grid))
            enddo

            deallocate(Hpsi)
        enddo
        call parallel_task('sum', lam(:), parallel, 'diff_b')

        lam_plus=grid%Ecut
        lam_minus=maxval(lam)
        a_l=minval(lam)

        r=lam_plus-lam_minus
        r=r*0.5_dp
        c=lam_plus+lam_minus
        c=c*0.5_dp
        sig=r/(c-a_l)
        tau=2.0_dp/sig

        if(parallel%myid_k.eq.0) print *, 'Filter max:', lam_plus, 'Filter min', lam_minus
        nY=0
        do band_block=1, nband_block
            allocate(Hpsi(block_n_spinor(band_block)))
            allocate(C0(block_n_spinor(band_block)))
            allocate(C1(block_n_spinor(band_block)))
            allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(orbitals(band_block)%of)))

            do i=1, block_n_spinor(band_block)
                call allocate_field(Hpsi(i), grid, parallel)
                call allocate_field(C0(i), grid, parallel)
                call allocate_field(C1(i), grid, parallel)
            enddo
            Hpsi%grid=orb_grid
            C0%grid=orb_grid
            C1%grid=orb_grid

            do i=1, block_n_spinor(band_block)
                C0(i)%R=orbitals(band_block)%of(i)%R
            enddo
            call Apply_SinvH(orbitals(band_block)%of(:), Hpsi, grids, potentials, atoms, elements, parallel, all_PAW, calc_G=.true.)

            do i=1, block_n_spinor(band_block)
                C1(i)%R=(Hpsi(i)%R-c*C0(i)%R)*sig/r
            enddo
            do n=2, numerics%cheby%inner_steps
                sig_2=1.0_dp/(tau-sig)

                call Apply_SinvH(C1(:), Hpsi, grids, potentials, atoms, elements, parallel, all_PAW, calc_G=.true.)
                do i=1, block_n_spinor(band_block)
                    orbitals(band_block)%of(i)%R=2.0_dp*sig_2/r*(Hpsi(i)%R-c*C1(i)%R)-sig_2*sig*C0(i)%R
                    C0(i)%R=C1(i)%R
                    C1(i)%R=orbitals(band_block)%of(i)%R
                    sig=sig_2
                enddo
            enddo

            do i = 1, size(Hpsi)
                call real_to_recip(orbitals(band_block)%of(i), grids)
                orbitals(band_block)%of(i)%G=orbitals(band_block)%of(i)%G*grid%cutwf
                call recip_to_real(orbitals(band_block)%of(i), grids)
                call Calculate_Projector_overlaps(orbitals(band_block)%of(i),  &
                PYpsi(:,:,i), atoms, elements, grid, parallel) 
            enddo

            call Apply_S(orbitals(band_block)%of(:), Hpsi, grids, atoms, elements, parallel,  calc_R=.false., PYpsi_in=PYpsi)          
            call field_FFT2Matrix(Hpsi, BY(:,(nY+1):(nY+nX)), grid, parallel, nX, grid%reduction_option)
            call Apply_H(orbitals(band_block)%of(:), Hpsi, grids, potentials, atoms, elements, parallel, &
                                  calc_R=.false., PYpsi_in=PYpsi)          
            call field_FFT2Matrix(Hpsi, AY(:,(nY+1):(nY+nX)), grid, parallel, nX, grid%reduction_option)

            nY=nY+nX
            do i = 1, size(Hpsi)
                call deallocate_field(Hpsi(i),grids(Hpsi(i)%grid))
                call deallocate_field(C0(i),grids(C0(i)%grid))
                call deallocate_field(C1(i),grids(C1(i)%grid))
            enddo

            deallocate(C0)
            deallocate(Hpsi)
            deallocate(C1)
            deallocate(PYpsi)
        enddo

        call orbitals_FFT2Matrix(orbitals(:), Y(:,:), grids, parallel, n_band_groups, &
                                 grids(orb_grid)%reduction_option)
        call B_orthogonalize(Y(:,:), BY(:,:), nY, parallel, AY(:,:), gamma=grids(orb_grid)%gamma)
        deallocate(BY)

        !print *, 'The Big Ritz'
        allocate(Ywork(size(Y,1),nY))
        call compute_Ritz_orthobasis(Y(:,:), AY(:,:), nY, Ywork(:,:), lam(:), &
            parallel, rot_A=.false., gamma=grids(orb_grid)%gamma)
        !print *, 'The Big Ritz-Done'
        deallocate(Ywork)
        deallocate(AY)

        all_eigs(1:nY)=lam(1:nY)  
        deallocate(lam)
        call orbitals_Matrix2FFT(orbitals(:nband_block), Y(:,1:nY), grids, parallel, nX, &
                                 grid%reduction_option)
        deallocate(Y)

        do band_block=1, nband_block
            orbitals(band_block)%eig(:)=all_eigs(orbitals(band_block)%band)
            do i=1, size(orbitals(band_block)%of(:))
                orbitals(band_block)%of(i)%G=orbitals(band_block)%of(i)%G*grid%cutwf
                call recip_to_real(orbitals(band_block)%of(i), grids)
            enddo
        enddo

    end subroutine

    subroutine compute_Ritz_nonorthobasis(nX, lam, X, AX, Xwork, parallel, gamma)
        use lapack, only : zgemm, ztrttp, ztpttr, zherk
        use parallel_type, only : parallel_struct
        use parallel_mod,only : parallel_task, parallel_wait

        use linalg, only : eigh

        type(parallel_struct), intent(in) :: parallel

        integer, intent(in) :: nX
        real(dp), intent(inout) :: lam(:)
        complex(dp), intent(inout) :: X(:,:), AX(:,:), Xwork(:,:)
        logical, intent(in) :: gamma

        complex(dp), allocatable :: gramA(:,:), U(:,:), gramB(:,:), Oxx(:,:) 
        integer :: ierr, nG, nXt, i, j
        complex(dp), allocatable :: send(:)
       
        nG=size(X,1)
        nXt=nX*(nX+1)/2

        allocate(gramA(nX,nX), gramB(nX,nX), U(nX,nX))
        allocate(send(nXt),Oxx(nX,nX))

        Oxx=0.0_dp
        call zgemm('C','N', nX, nX, nG, (1.0_dp,0.0_dp), X(:,1:nX),  &
                     nG, AX(:,1:nX), nG, (0.0_dp,0.0_dp),  Oxx, nX)

        call parallel_task('sum', Oxx, parallel, 'k')

        do i=1, nX
            do j=i,nX
                Oxx(j,i)=(Oxx(j,i)+conjg(Oxx(i,j)))*0.5_dp
                Oxx(i,j)=conjg(Oxx(j,i))
            enddo
        enddo
       
        if(gamma) then
            gramA=real(Oxx)
        else
            gramA=Oxx
        endif
        ! if(parallel%myid.eq.0) print *, parallel%myid, 'gramA: '
        ! do i=1, size(gramA,1)
        !     if(parallel%myid.eq.0) print *, real(gramA(i,:))
        ! enddo

       ! call zherk('U','C', nX, nG, (1.0_dp,0.0_dp), X(:,1:nX),  &
       ! nG, (0.0_dp,0.0_dp),  Oxx, nX)
       ! call ztrttp('U',nX, Oxx, nX, send, ierr) 
       ! call parallel_task('sum', send, parallel, 'k')        
       ! call ztpttr('U',nX, send, Oxx, nX, ierr) 

        Oxx=0.0_dp
        call zgemm('C','N', nX, nX, nG, (1.0_dp,0.0_dp), X(:,1:nX),  &
                     nG, X(:,1:nX), nG, (0.0_dp,0.0_dp),  Oxx, nX)

        call parallel_task('sum', Oxx, parallel, 'k')

        do i=1, nX
            do j=i,nX
                Oxx(j,i)=(Oxx(j,i)+conjg(Oxx(i,j)))*0.5_dp
                Oxx(i,j)=conjg(Oxx(j,i))
            enddo
        enddo
       
        if(gamma) then
            gramB=real(Oxx)
        else
            gramB=Oxx
        endif

        ! if(parallel%myid.eq.0) print *, parallel%myid, 'gramB: ', gramB

        deallocate(send,Oxx)
        call eigh(gramA, gramB, lam(1:nX), U)

        !U is now the Ritz vectors lam is the ritz
        !|X> = |Xin>*U
        call zgemm('N','N', size(X,1), nX, nX, (1.0_dp,0.0_dp), X,  &
        size(X,1), U, size(U,1), (0.0_dp,0.0_dp),  Xwork, size(Xwork,1))
        X=Xwork

       ! if(parallel%myid.eq.0) print *, parallel%myid, 'U: ', U

       ! if(parallel%myid.eq.0) print *, parallel%myid, 'lam: ', lam

        deallocate(gramA,gramB,U)
end subroutine

    

end module