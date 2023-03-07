module Eigen_LOBPCG
    use types, only : dp

    implicit none

    public :: LOBPCG, Orbitals_Project_out, Orbitals_Orthogonalize


    contains

    subroutine LOBPCG(orbitals, grids, parallel, numerics, potentials, atoms, elements, &
        all_eigs, all_PAW, n_band_groups)
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
        type(all_PAW_struct), intent(in) :: all_PAW
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout) :: grids(:)
        real(dp) , intent(inout) :: all_eigs(:)
        integer, intent(in), optional :: n_band_groups
        integer :: nX
        logical :: test_S=.false.

        if(present(n_band_groups)) then
            nX=n_band_groups
        else
            nX=parallel%n_band_groups
        endif

        if(all_PAW%N_PAW_atoms.gt.0) then 
            call LOBPCG_Gen(orbitals, grids, parallel, numerics, potentials, atoms, elements, &
                all_eigs, nX)
        else
            call LOBPCG_NC(orbitals, grids, parallel, numerics, potentials, atoms, elements, &
            all_eigs, nX)
                
        endif

    end subroutine 

    subroutine LOBPCG_NC(orbitals, grids, parallel, numerics, potentials, atoms, elements, &
             all_eigs, n_band_groups)
        use grids_type, only : grid_struct
        use odp_type, only : field_struct, orbital_struct
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: recip_to_real
        use lapack, only : zgemm

        use parallel_mod,only : parallel_task, parallel_wait
        use parallel_type, only : parallel_struct
        use numerics_type, only : numerics_struct

        use FFT2Matrix, only : orbitals_FFT2Matrix, orbitals_Matrix2FFT, grids_G2_FFT2Matrix, &
                               orbital_FFT2Matrix
        use odp, only: allocate_field, deallocate_field, allocate_orbital, deallocate_orbital 
        logical :: test_S=.false.


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


        type(field_struct), allocatable :: psi(:), Hpsi(:), work(:), Hs2psi(:)
        real(dp), allocatable :: G2_2(:,:), Kpc(:,:), lam(:), resnorm(:), Ekin(:)
        complex(dp), allocatable :: XPW(:,:), AXPW(:,:),  Xwork(:,:), Xwork2(:,:)
        complex(dp), allocatable :: YXP(:,:),  AY(:,:), Ywork(:,:)
        complex(dp), allocatable :: CX(:,:), CP(:,:)
        integer, allocatable :: block_n_spinor(:)
        integer :: band_block, nY, nX, nband_block, nG, nGblock
        integer:: i,sW,eW,sP,eP, nc, n, j, skip
        integer :: orb_grid
        logical :: do_p=.false.

        real(dp), allocatable:: occ(:,:), eigs(:,:), weight(:,:), filter(:,:)

        if(parallel%myid.eq.0) print *, 'Starting LOBPCG'
        nband_block=size(orbitals)
        nX=n_band_groups
        orb_grid=orbitals(1)%of(1)%grid

        do i=1,nband_block
            if(orbitals(i)%of(1)%grid.ne.orb_grid) then
                print *, 'All orbitals passed into Egensolver need to be on the same grid (same k-point)'
            endif
        enddo

        !Setup Preconditioner
        grid=>grids(orbitals(1)%of(1)%grid)
        if(grid%reduction_option.eq.2) then
            !if(parallel%myid.eq.0) print *, 'Removing some of the zeros in x and y directions'
            nG=product(grid%Ng_small)
        else
            nG=product(grid%Ng_local)
        endif

        if(mod(nG,nX).ne.0) nG=(nG/nX+1)*nX
        allocate(G2_2(nG/nX,nX))
        !Set Preconditioner
        call grids_G2_FFT2Matrix(grid, G2_2(:,1:nX), parallel, n_band_groups, grid%reduction_option)
        G2_2=G2_2*0.5_dp

        if(numerics%lobpcg%PC_type.eq.1) then
            allocate(Kpc(size(G2_2,1),nX))
            Kpc=27.0_dp+G2_2*(18.0_dp + G2_2*(12.0_dp + 8*G2_2))
            Kpc=2.0_dp*Kpc/(Kpc+16.0_dp*G2_2**4)
            deallocate(G2_2)
        else
            allocate(Ekin(nX))
            allocate(Kpc(size(G2_2,1),nX))
        endif

        all_eigs(:)=0.0_dp
        !Move orbitals into Matrix Form for LAPACK
        !Keep them stored in the unused portion of YXP vector
        !deallocate the orbitals to save space

        nY=nX*nband_block
        allocate(YXP(maxval(orbitals(:)%n_spinor)*nG/nX,nY+nX))
        allocate(AY(size(YXP,1),nY))
        allocate(lam(nY+2*nX))

        allocate(block_n_spinor(nband_block))
        block_n_spinor(:)=orbitals(:)%n_spinor

        call orbitals_FFT2Matrix(orbitals(:nband_block), YXP(:,(nX+1):), grids, parallel, nX, &
                                 grids(orb_grid)%reduction_option)

        allocate(occ(maxval(orbitals(:)%n_spinor),nband_block))
        allocate(weight(maxval(orbitals(:)%n_spinor),nband_block))
        allocate(eigs(maxval(orbitals(:)%n_spinor),nband_block))
        allocate(filter(maxval(orbitals(:)%n_spinor),nband_block))
        do band_block=1, nband_block
            occ(:block_n_spinor(band_block),band_block)=orbitals(band_block)%occ(:)
            weight(:block_n_spinor(band_block),band_block)=orbitals(band_block)%weight(:)
            filter(:block_n_spinor(band_block),band_block)=orbitals(band_block)%filter(:)
            eigs(:block_n_spinor(band_block),band_block)=orbitals(band_block)%eig(:)
            !Dont deallocate the orbitals array itself, retain the index values
            call deallocate_orbital(orbitals(band_block),grids)
        enddo
        allocate(resnorm(nX))

        allocate(CX(nX*3,nX))
        allocate(CP(nX*2,nX))
        nY=0
      
        !Big Loop over Locally Optimized Blocks
        do band_block=1, nband_block
            nGblock=nG*block_n_spinor(band_block)/nX

            allocate(XPW(nGblock,nX*3))
            allocate(Xwork(nGblock,nX))
            allocate(Xwork2(nGblock,nX))
            allocate(AXPW(nGblock,nX*3))
            XPW=0
            AXPW=0
            Xwork=0
            AXPW=0

            allocate(Hpsi(block_n_spinor(band_block)))
            allocate( psi(block_n_spinor(band_block)))
            allocate(work(block_n_spinor(band_block)))
            
            do i=1, block_n_spinor(band_block)
                call allocate_field(Hpsi(i), grid, parallel)
                call allocate_field(psi(i), grid, parallel)
                call allocate_field(work(i), grid, parallel)
            enddo
            psi%grid=orb_grid
            Hpsi%grid=orb_grid
            work%grid=orb_grid

            XPW(:,1:nX)=YXP(:nGblock,(nY+nX+1):(nY+2*nX))

            !Make X (new) orthogonal to Y (found)
            if(nY.gt.0) then
                call projection_removal_orthobasis(XPW(:,1:nX), nX, YXP(:nGblock,1:nY), nY, parallel, gamma=grid%gamma)
            endif
            !Make X (new) orthogonal to themselves
            call orthogonalize(XPW(:,1:nX), nX, parallel, gamma=grid%gamma)

            !|AX>= H|X>
            call av(XPW(:,1:nX), nX, AXPW(:,1:nX), parallel, grids, &
                orb_grid, potentials, atoms, elements, psi, Hpsi)

            !Initial Ritz
            call compute_Ritz_orthobasis(XPW(:,1:nX), AXPW(:,1:nX), nX, Xwork, lam(1:nX), parallel, rot_A=.true., gamma=grid%gamma)

            !YXP(:nGblock,(nY+1):(nY+nX))=XPW(:,1:nX) 
            do_p=.false.
            nc=0
            sP=0
            eP=-1
            sW=nX+1
            eW=2*nX
            !Calculate Residual
            !W=A.X-lam*X
            do i=1, nX
                XPW(:,nX+i)=AXPW(:,i)-XPW(:,i)*lam(i)
            enddo

            do i=1, nX
                resnorm(i)=real(sum(conjg(XPW(:,nX+i))*XPW(:,nX+i)))
            enddo
            call parallel_task('sum', resnorm, parallel, 'k')   
         !   if(parallel%myid.eq.0) print *, band_block, 'resnorm ', resnorm(:nX), 'lam ', lam(:nX); flush(6)
            
            do i=1, nX
                if(resnorm(i).lt.numerics%lobpcg%soft_lock_thresh) then
                     nc=nc+1
                else
                    exit
                endif
            enddo
            if(nc.eq.nX) then
                skip=numerics%lobpcg%inner_steps
            else 
                skip=0
            endif


            !Inner Loop for number of refinement steps
            do n=1, numerics%lobpcg%inner_steps-skip
                    if(numerics%lobpcg%PC_type.eq.2) then
                        do i=sW, eW
                            j=i-sW+1
                            Ekin(j)=real(sum(conjg(XPW(:,i))*G2_2(:,j)*XPW(:,i)))
                        enddo
                        call parallel_task('sum', Ekin, parallel, 'k')   
                        Ekin=Ekin*1.5_dp
                        do i=sW, eW
                            j=i-sW+1
                            G2_2(:,j)=G2_2(:,j)/Ekin(j)
                        enddo
                        Kpc=27.0_dp+G2_2*(18.0_dp+G2_2*(12.0_dp+8.0_dp*G2_2))
                        do i=sW, eW
                            j=i-sW+1
                            Kpc(:,j)=Kpc(:,j)/(Kpc(:,j)+16.0_dp*G2_2(:,j)**4)
                            G2_2(:,j)=G2_2(:,j)*Ekin(j)
                        enddo
                    endif
                    if(.not.test_S) XPW(:,sW:eW)=XPW(:,sW:eW)*Kpc(:nGblock,1:(eW-sW+1))
                    if(nY.gt.0) then
                        !call orthogonalize(XPW(:,sW:eW), nX-nc, parallel)
                        call projection_removal_orthobasis(XPW(:,sW:eW), eW-sW+1, YXP(:nGblock,1:nY), nY, parallel,gamma=grid%gamma)
                    endif
                    call orthogonalize(XPW(:,sW:eW), nX-nc, parallel, gamma=grid%gamma)
                    if(do_p) call orthogonalize(XPW(:,sP:eP), nX-nc, parallel, AXPW(:,sP:eP), gamma=grid%gamma)

                    !|AW>= H|W>
                    !Need to apply to all nX (even unused) to prevent rotation errors
                   
                    call av(XPW(:,sW:(sW-1+nX)), nX, AXPW(:,sW:(sW-1+nX)), parallel, grids, &
                        orb_grid, potentials, atoms, elements, psi, Hpsi)

                    call Rayleigh_Ritz_traditional(nX, nc, eW,&
                                                   lam, XPW, AXPW, CX, CP, parallel, do_p, grid%gamma)

                    call zgemm('N','N', size(XPW,1), nX, eW, (1.0_dp,0.0_dp), XPW(:,1:eW),  &
                    size(XPW,1), CX(1:eW,1:nX), size(CX(1:eW,1:nX),1), (0.0_dp,0.0_dp), Xwork, &
                    size(Xwork,1))

                    XPW(:,1:nX)=Xwork(:,1:nX) !XPW(:,(nX+1):eW) (P & W) remains the same

                    call zgemm('N','N', size(AXPW,1), nX, eW, (1.0_dp,0.0_dp), AXPW(:,1:eW),  &
                    size(AXPW,1), CX(1:eW,1:nX), size(CX(1:eW,1:nX),1), (0.0_dp,0.0_dp),  Xwork, size(Xwork,1))
                    AXPW(:,1:nX)=Xwork(:,1:nX)

                    !new W
                    !Calculate Residual & find lockable vectors
                    do i=1, nX
                        Xwork(:,i)=AXPW(:,i)-XPW(:,i)*lam(i)
                        resnorm(i)=real(sum(conjg(Xwork(:,i))*Xwork(:,i)))
                    enddo
                    call parallel_task('sum', resnorm, parallel, 'k')   
               !    if(parallel%myid.eq.0) print *, band_block, 'resnorm ', resnorm(:nX), 'lam ', lam(:nX)
                    do i=nc+1, nX
                        if(resnorm(i).lt.numerics%lobpcg%soft_lock_thresh) then
                             nc=nc+1
                        else
                            exit
                        endif
                    enddo
                    if(nc.eq.nX) then
                        exit
                    else 
                    endif

                    do_p=.true.
                    sP=nX+1
                    eP=sP-1+nX-nc
                    
                    !only keep those corresponding to unconverged
                    call zgemm('N','N', size(XPW,1), nX-nc, eW-nX, (1.0_dp,0.0_dp), XPW(:,sP:eW),  &
                    size(XPW,1), CP(1:(eW-nX),(nc+1):nX), size(CP(1:(eW-nX),(nc+1):nX),1), &
                    (0.0_dp,0.0_dp), Xwork2(:,(nc+1):nX), size(Xwork2(:,(nc+1):nX),1))

                    !Now we can update P & W to the new values
                    XPW(:,sP:eP)=Xwork2(:,(nc+1):nX)

                    call zgemm('N','N', size(AXPW,1), nX-nc, eW-nX, (1.0_dp,0.0_dp), AXPW(:,sP:eW),  &
                    size(AXPW,1), CP(1:(eW-nX),(nc+1):nX), size(CP(1:(eW-nX),(nc+1):nX),1), &
                    (0.0_dp,0.0_dp), Xwork2(:,(nc+1):nX), size(Xwork2(:,(nc+1):nX),1))

                    AXPW(:,sP:eP)=Xwork2(:,(nc+1):nX)

                    sW=eP+1
                    eW=sW-1+nX-nc

                    XPW(:,sW:eW)=Xwork(:,(nc+1):nX)

            enddo

            YXP(:nGblock,(nY+1):(nY+nX))= XPW(:,1:nX)
             AY(:nGblock,(nY+1):(nY+nX))=AXPW(:,1:nX)

            if( maxval(block_n_spinor(:)).eq.2 .and. block_n_spinor(band_block).eq.1 ) then
                YXP((nGblock+1):(2*nGblock),(nY+1):(nY+nX)) = XPW(:,1:nX)
                AY(:nGblock,(nY+1):(nY+nX)) = AXPW(:,1:nX) 
            else if( maxval(block_n_spinor(:)).eq.4 .and. block_n_spinor(band_block).eq.1 ) then
                YXP((nGblock+1):(2*nGblock),(nY+1):(nY+nX)) = XPW(:,1:nX)
                AY(:nGblock,(nY+1):(nY+nX)) = AXPW(:,1:nX) 
                YXP((nGblock*2+1):(4*nGblock),(nY+1):(nY+nX)) = YXP((nGblock+1):(2*nGblock),(nY+1):(nY+nX))
                AY((nGblock*2+1):(4*nGblock),(nY+1):(nY+nX)) =  AY((nGblock+1):(2*nGblock),(nY+1):(nY+nX))
            endif

            nY=nY+nX !permanently add the converged vectors to the hard-lock list
            !Deallocate the work spaces before the big Ritz
            do i = 1, size(psi)
                call deallocate_field(work(i),grids(work(i)%grid))
                call deallocate_field(Hpsi(i),grids(Hpsi(i)%grid))
                call deallocate_field(psi(i),grids(psi(i)%grid))
            enddo
            deallocate(work)
            deallocate(Hpsi)
            deallocate(psi)

            deallocate(XPW)
            deallocate(Xwork)
            deallocate(Xwork2)
            deallocate(AXPW)

        enddo ! End loop on blocks

        deallocate(resnorm)
        deallocate(CX)
        deallocate(CP)
        deallocate(Kpc)
        if(numerics%lobpcg%PC_type.eq.2) then
            deallocate(Ekin)
            deallocate(G2_2)
        endif
        !The Big Ritz
        allocate(Ywork(size(YXP,1),nY))
        call compute_Ritz_orthobasis(YXP(:,1:nY), AY(:,1:nY), nY, Ywork(:,1:nY), lam(1:nY), &
            parallel, rot_A=.true., gamma=grid%gamma)

        deallocate(Ywork)
        deallocate(AY)

        all_eigs(1:nY)=lam(1:nY)
  
        deallocate(lam)
        do band_block=1, nband_block
            call allocate_orbital(orbitals(band_block), grids, orb_grid, parallel)
            orbitals(band_block)%eig(:)=all_eigs(orbitals(band_block)%band)
            !Need to adjust eig for spinor components if definable like in restricted open-shell collinear
            orbitals(band_block)%filter(:)=filter(:block_n_spinor(band_block),band_block)
            orbitals(band_block)%occ(:)=occ(:block_n_spinor(band_block),band_block)
            orbitals(band_block)%weight(:)=weight(:block_n_spinor(band_block),band_block)
        enddo
        deallocate(filter)
        deallocate(occ)
        deallocate(weight)
        deallocate(eigs)


        call orbitals_Matrix2FFT(orbitals(:nband_block), YXP(:,1:nY), grids, parallel, nX, &
                                 grid%reduction_option)

        do band_block=1, nband_block
            do i=1, size(orbitals(band_block)%of(:))
                orbitals(band_block)%of(i)%G=orbitals(band_block)%of(i)%G*grid%cutwf
                call recip_to_real(orbitals(band_block)%of(i), grids)
            enddo
        enddo

        if(test_S) then
            if(parallel%myid.eq.0) then
                do i=   1,nY
                    print *, i, all_eigs(i)
                enddo
            endif
        endif
        deallocate(YXP)
    end subroutine

    subroutine LOBPCG_Gen(orbitals, grids, parallel, numerics, potentials, atoms, elements, &
        all_eigs, n_band_groups)
   use grids_type, only : grid_struct
   use odp_type, only : field_struct, orbital_struct
   use simulation_type, only : potentials_struct
   use element_type, only : element_struct
   use atom_type, only : atom_struct
   use fft, only: recip_to_real

   use parallel_mod,only : parallel_task, parallel_wait
   use parallel_type, only : parallel_struct
   use numerics_type, only : numerics_struct

   use FFT2Matrix, only : orbitals_FFT2Matrix, orbitals_Matrix2FFT, grids_G2_FFT2Matrix, &
                          orbital_FFT2Matrix
   use odp, only: allocate_field, deallocate_field, allocate_orbital, deallocate_orbital 
   use lapack, only : zgemm

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

   type(field_struct), allocatable :: psi(:), Hpsi(:), Spsi(:), work(:), Hs2psi(:)
   real(dp), allocatable :: G2_2(:,:), Kpc(:,:), lam(:), resnorm(:), Ekin(:)
   complex(dp), allocatable :: XPW(:,:), AXPW(:,:), BXPW(:,:), Xwork(:,:), Xwork2(:,:)
   complex(dp), allocatable :: YXP(:,:),  AY(:,:), BY(:,:), Ywork(:,:)
   complex(dp), allocatable :: CX(:,:), CP(:,:)
   integer, allocatable :: block_n_spinor(:)
   integer :: band_block, nY, nX, nband_block, nG, nGblock
   integer:: i,j,sW,eW,sP,eP, nc, n, skip
   integer :: orb_grid
   logical :: do_p=.false.

   real(dp), allocatable:: occ(:,:), eigs(:,:), weight(:,:), filter(:,:)

   if(parallel%myid.eq.0) print *, 'Starting Generalized LOBPCG'
   nband_block=size(orbitals)
   nX=n_band_groups

   orb_grid=orbitals(1)%of(1)%grid

   do i=1,nband_block
       if(orbitals(i)%of(1)%grid.ne.orb_grid) then
           print *, 'All orbitals passed into Egensolver need to be on the same grid (same k-point)'
       endif
   enddo

   !Setup Preconditioner
   grid=>grids(orbitals(1)%of(1)%grid)
   if(grid%reduction_option.eq.2) then
       !if(parallel%myid.eq.0) print *, 'Removing some of the zeros in x and y directions'
       nG=product(grid%Ng_small)
   else
       nG=product(grid%Ng_local)
   endif

   if(mod(nG,nX).ne.0) nG=(nG/nX+1)*nX
   allocate(G2_2(nG/nX,nX))
   !Set Preconditioner
   call grids_G2_FFT2Matrix(grid, G2_2(:,1:nX), parallel, n_band_groups, grid%reduction_option)
   G2_2=G2_2*0.5_dp

   if(numerics%lobpcg%PC_type.eq.1) then
       allocate(Kpc(size(G2_2,1),nX))
       Kpc=27.0_dp+G2_2*(18.0_dp + G2_2*(12.0_dp + 8*G2_2))
       Kpc=Kpc/(Kpc+16.0_dp*G2_2**4)
       deallocate(G2_2)
   else
       allocate(Ekin(nX))
       allocate(Kpc(size(G2_2,1),nX))
   endif

   all_eigs(:)=0.0_dp
   !Move orbitals into Matrix Form for LAPACK
   !Keep them stored in the unused portion of YXP vector
   !deallocate the orbitals to save space

   nY=nX*nband_block
   allocate(YXP(maxval(orbitals(:)%n_spinor)*nG/nX,nY+nX))
   allocate(AY(size(YXP,1),nY))
   allocate(BY(size(YXP,1),nY))
   !allocate(Oyy(nY,nY))
   !Oyy=0.0_dp

   allocate(lam(nY+2*nX))

   allocate(block_n_spinor(nband_block))
   block_n_spinor(:)=orbitals(:)%n_spinor

   call orbitals_FFT2Matrix(orbitals(:nband_block), YXP(:,(nX+1):), grids, parallel, nX, &
                            grids(orb_grid)%reduction_option)

   allocate(occ(maxval(orbitals(:)%n_spinor),nband_block))
   allocate(weight(maxval(orbitals(:)%n_spinor),nband_block))
   allocate(eigs(maxval(orbitals(:)%n_spinor),nband_block))
   allocate(filter(maxval(orbitals(:)%n_spinor),nband_block))
   do band_block=1, nband_block
       occ(:block_n_spinor(band_block),band_block)=orbitals(band_block)%occ(:)
       weight(:block_n_spinor(band_block),band_block)=orbitals(band_block)%weight(:)
       filter(:block_n_spinor(band_block),band_block)=orbitals(band_block)%filter(:)
       eigs(:block_n_spinor(band_block),band_block)=orbitals(band_block)%eig(:)
       !Dont deallocate the orbitals array itself, retain the index values
       call deallocate_orbital(orbitals(band_block),grids)
   enddo
   allocate(resnorm(nX))
   allocate(CX(nX*3,nX))
   allocate(CP(nX*2,nX))
   nY=0
 
   !Big Loop over Locally Optimized Blocks
   do band_block=1, nband_block
       nGblock=nG*block_n_spinor(band_block)/nX

       allocate(XPW(nGblock,nX*3))
       allocate(Xwork(nGblock,nX))
       allocate(Xwork2(nGblock,nX))
       allocate(AXPW(nGblock,nX*3))
       allocate(BXPW(nGblock,nX*3))
       XPW=0
       AXPW=0
       Xwork=0
       AXPW=0
       BXPW=0
       allocate(Hpsi(block_n_spinor(band_block)))
       allocate(Spsi(block_n_spinor(band_block)))
       allocate( psi(block_n_spinor(band_block)))
       allocate(work(block_n_spinor(band_block)))
       
       do i=1, block_n_spinor(band_block)
           call allocate_field(Hpsi(i), grid, parallel)
           call allocate_field(Spsi(i), grid, parallel)
           call allocate_field(psi(i), grid, parallel)
           call allocate_field(work(i), grid, parallel)
       enddo
       psi%grid=orb_grid
       Hpsi%grid=orb_grid
       Spsi%grid=orb_grid

       work%grid=orb_grid

       XPW(:,1:nX)=YXP(:nGblock,(nY+nX+1):(nY+2*nX))

       !Make X (new) B_orthogonal to Y (found)
       if(nY.gt.0) then
           call projection_removal_nonorthobasis(XPW(:,1:nX), nX, YXP(:nGblock,1:nY), BY(:nGblock,1:nY), &
                             nY, parallel, gamma=grid%gamma)
       endif

       if(.true.) then
            call av_bv(XPW(:,1:nX), nX, AXPW(:,1:nX), BXPW(:,1:nX), parallel, grids, &
            orb_grid, potentials, atoms, elements, psi, Hpsi, Spsi)
             !Make X (new) B_orthogonal to themselves and rotate AX-> AX'
            call B_orthogonalize(XPW(:,1:nX), BXPW(:,1:nX), nX, parallel,  AXPW(:,1:nX), grid%gamma)
       else
            !|BX>= S|X>
            call bv(XPW(:,1:nX), nX, BXPW(:,1:nX), parallel, grids, &
            orb_grid, atoms, elements, psi, Spsi)

            !Make X (new) B_orthogonal to themselves and rotate AX-> AX'
            call B_orthogonalize(XPW(:,1:nX), BXPW(:,1:nX), nX, parallel, gamma=grid%gamma)

            !|AX>= H|X>
            call av(XPW(:,1:nX), nX, AXPW(:,1:nX), parallel, grids, &
                orb_grid, potentials, atoms, elements, psi, Hpsi)
       endif
       !Initial Ritz
       call compute_Ritz_orthobasis(XPW(:,1:nX), AXPW(:,1:nX), nX, Xwork, lam(1:nX), parallel, &
                         rot_A=.true., BX=BXPW(:,1:nX),gamma=grid%gamma)

       !YXP(:nGblock,(nY+1):(nY+nX))=XPW(:,1:nX) 
       do_p=.false.
       nc=0
       sP=0
       eP=-1
       sW=nX+1
       eW=2*nX
       !Calculate Residual
       !W=A.X-lam*X
       do i=1, nX
           XPW(:,nX+i)=AXPW(:,i)-BXPW(:,i)*lam(i)
           resnorm(i)=real(sum(conjg(XPW(:,nX+i))*XPW(:,nX+i)))
        enddo
        call parallel_task('sum', resnorm, parallel, 'k') 
     !  if(parallel%myid.eq.0) print *, band_block, ' Initial resnorm ', resnorm(:nX), 'lam ', lam(:nX); flush(6)
        do i=1, nX
            if(resnorm(i).lt.numerics%lobpcg%soft_lock_thresh) then
                 nc=nc+1
            else
                exit
            endif
        enddo
        if(nc.eq.nX) then
            skip=numerics%lobpcg%inner_steps
        else 
            skip=0
        endif
       !Inner Loop for number of refinement steps
       do n=1, numerics%lobpcg%inner_steps-skip
            ! print *,  'inner loop:', n,  nc, nX, sW, eW, sP, eP
             if(numerics%lobpcg%PC_type.eq.2) then
                do i=sW, eW
                    j=i-sW+1
                    Ekin(j)=real(sum(conjg(XPW(:,i))*G2_2(:,j)*XPW(:,i)))
                enddo
                call parallel_task('sum', Ekin, parallel, 'k')   
                Ekin=Ekin*1.5_dp
                do i=sW, eW
                    j=i-sW+1
                    G2_2(:,j)=G2_2(:,j)/Ekin(j)
                enddo
                Kpc=27.0_dp+G2_2*(18.0_dp+G2_2*(12.0_dp+8.0_dp*G2_2))
                do i=sW, eW
                    j=i-sW+1
                    Kpc(:,j)=Kpc(:,j)/(Kpc(:,j)+16.0_dp*G2_2(:,j)**4)
                    G2_2(:,j)=G2_2(:,j)*Ekin(j)
                enddo
            endif
            XPW(:,sW:eW)=XPW(:,sW:eW)*Kpc(:nGblock,1:(eW-sW+1))
            if(nY.gt.0) then
                !call orthogonalize(XPW(:,sW:eW), nX-nc, parallel)
                call projection_removal_nonorthobasis(XPW(:,sW:eW), eW-sW+1, YXP(:nGblock,1:nY), BY(:nGblock,1:nY), &
                     nY, parallel, gamma=grid%gamma)
            endif

            if(.true.) then
                call av_bv(XPW(:,sW:(sW-1+nX)), nX, AXPW(:,sW:(sW-1+nX)), BXPW(:,sW:(sW-1+nX)), parallel, grids, &
                orb_grid, potentials, atoms, elements, psi, Hpsi, Spsi)
                 !Make X (new) B_orthogonal to themselves and rotate AX-> AX'
                 call B_orthogonalize(XPW(:,sW:eW), BXPW(:,sW:eW), nX-nc, parallel,  AXPW(:,sW:(sW-1+nX)), grid%gamma)

                if(do_p) then
                    if(nc.gt.0) then
                        Xwork=1.0_dp
                        Xwork(:,1:(eP-sP+1))=XPW(:,sP:eP)
                        !upper P's here will be the W's if nc >0, based on ordering so need to only take the result from real P's 
                        call bv(Xwork(:,1:nX), nX, Xwork2(:,1:nX), parallel, grids, &
                        orb_grid, atoms, elements, psi, Spsi)
                        BXPW(:,sP:eP)=Xwork2(:,1:(eP-sP+1))
                    else
                        call bv(XPW(:,sP:eP), nX, BXPW(:,sP:eP), parallel, grids, &
                        orb_grid, atoms, elements, psi, Spsi)
                    endif
                    call B_orthogonalize(XPW(:,sP:eP), BXPW(:,sP:eP), nX-nc, parallel, AXPW(:,sP:eP), grid%gamma)
                endif
            else
                !|BX>= S|X>, store projector overlaps
                call bv(XPW(:,sW:(sW-1+nX)), nX, BXPW(:,sW:(sW-1+nX)), parallel, grids, &
                        orb_grid, atoms, elements, psi, Spsi)
                call B_orthogonalize(XPW(:,sW:eW), BXPW(:,sW:eW), nX-nc, parallel, gamma=grid%gamma)
                if(do_p) then
                    if(nc.gt.0) then
                        Xwork=0.0_dp
                        Xwork(:,1:(eP-sP+1))=XPW(:,sP:eP)
                        !upper P's here will be the W's if nc >0, based on ordering so need to only take the result from real P's 
                        call bv(Xwork(:,1:nX), nX, Xwork2(:,1:nX), parallel, grids, &
                        orb_grid, atoms, elements, psi, Spsi)
                        BXPW(:,sP:eP)=Xwork2(:,1:(eP-sP+1))
                    else
                        call bv(XPW(:,sP:eP), nX, BXPW(:,sP:eP), parallel, grids, &
                        orb_grid, atoms, elements, psi, Spsi)
                    endif
                    call B_orthogonalize(XPW(:,sP:eP), BXPW(:,sP:eP), nX-nc, parallel, AXPW(:,sP:eP), grid%gamma)
                endif
                !|AW>= H|W>
                !Need to apply to all nX (even unused) to prevent rotation errors
                call av(XPW(:,sW:(sW-1+nX)), nX, AXPW(:,sW:(sW-1+nX)), &
                    parallel, grids, orb_grid, potentials, atoms, elements, psi, Hpsi)
            endif
            call Gen_Rayleigh_Ritz_traditional(nX, nc, eW,&
                                            lam, XPW, AXPW, BXPW, CX, CP, parallel, do_p, grid%gamma)

            call zgemm('N','N', size(XPW,1), nX, eW, (1.0_dp,0.0_dp), XPW(:,1:eW),  &
            size(XPW,1), CX(1:eW,1:nX), size(CX(1:eW,1:nX),1), (0.0_dp,0.0_dp), Xwork, &
            size(Xwork,1))

            XPW(:,1:nX)=Xwork(:,1:nX) !XPW(:,(nX+1):eW) (P & W) remains the same

            call zgemm('N','N', size(AXPW,1), nX, eW, (1.0_dp,0.0_dp), AXPW(:,1:eW),  &
            size(AXPW,1), CX(1:eW,1:nX), size(CX(1:eW,1:nX),1), (0.0_dp,0.0_dp),  Xwork, size(Xwork,1))
            AXPW(:,1:nX)=Xwork(:,1:nX)

            call zgemm('N','N', size(BXPW,1), nX, eW, (1.0_dp,0.0_dp), BXPW(:,1:eW),  &
            size(BXPW,1), CX(1:eW,1:nX), size(CX(1:eW,1:nX),1), (0.0_dp,0.0_dp),  Xwork, size(Xwork,1))
            BXPW(:,1:nX)=Xwork(:,1:nX)

            !new W
            !Calculate Residual & find lockable vectors
            do i=1, nX
                Xwork(:,i)=AXPW(:,i)-BXPW(:,i)*lam(i)
                resnorm(i)=real(sum(conjg(Xwork(:,i))*Xwork(:,i)))
            enddo
            call parallel_task('sum', resnorm, parallel, 'k')
     !       if(parallel%myid.eq.0) print *, band_block, 'resnorm ', resnorm(:nX), 'lam ', lam(:nX); flush(6)
            do i=nc+1, nX
                if(resnorm(i).lt.numerics%lobpcg%soft_lock_thresh) then
                    nc=nc+1
                else
                    exit
                endif
            enddo
            if(nc.eq.nX) then
                exit
            else 
            endif
            do_p=.true.
            sP=nX+1
            eP=sP-1+nX-nc
            
            !only keep those corresponding to unconverged
            call zgemm('N','N', size(XPW,1), nX-nc, eW-nX, (1.0_dp,0.0_dp), XPW(:,sP:eW),  &
            size(XPW,1), CP(1:(eW-nX),(nc+1):nX), size(CP(1:(eW-nX),(nc+1):nX),1), &
            (0.0_dp,0.0_dp), Xwork2(:,(nc+1):nX), size(Xwork2(:,(nc+1):nX),1))

            !Now we can update P & W to the new values
            XPW(:,sP:eP)=Xwork2(:,(nc+1):nX)

            call zgemm('N','N', size(AXPW,1), nX-nc, eW-nX, (1.0_dp,0.0_dp), AXPW(:,sP:eW),  &
            size(AXPW,1), CP(1:(eW-nX),(nc+1):nX), size(CP(1:(eW-nX),(nc+1):nX),1), &
            (0.0_dp,0.0_dp), Xwork2(:,(nc+1):nX), size(Xwork2(:,(nc+1):nX),1))

            AXPW(:,sP:eP)=Xwork2(:,(nc+1):nX)
            
            call zgemm('N','N', size(BXPW,1), nX-nc, eW-nX, (1.0_dp,0.0_dp), BXPW(:,sP:eW),  &
            size(BXPW,1), CP(1:(eW-nX),(nc+1):nX), size(CP(1:(eW-nX),(nc+1):nX),1), &
            (0.0_dp,0.0_dp), Xwork2(:,(nc+1):nX), size(Xwork2(:,(nc+1):nX),1))

            BXPW(:,sP:eP)=Xwork2(:,(nc+1):nX)

            sW=eP+1
            eW=sW-1+nX-nc

            XPW(:,sW:eW)=Xwork(:,(nc+1):nX)

       enddo
       YXP(:nGblock,(nY+1):(nY+nX))= XPW(:,1:nX)
        AY(:nGblock,(nY+1):(nY+nX))=AXPW(:,1:nX)
        BY(:nGblock,(nY+1):(nY+nX))=BXPW(:,1:nX)

       if( maxval(block_n_spinor(:)).eq.2 .and. block_n_spinor(band_block).eq.1 ) then
           YXP((nGblock+1):(2*nGblock),(nY+1):(nY+nX)) = XPW(:,1:nX)
           AY(:nGblock,(nY+1):(nY+nX)) = AXPW(:,1:nX) 
           BY(:nGblock,(nY+1):(nY+nX)) = BXPW(:,1:nX) 

       else if( maxval(block_n_spinor(:)).eq.4 .and. block_n_spinor(band_block).eq.1 ) then
           YXP((nGblock+1):(2*nGblock),(nY+1):(nY+nX)) = XPW(:,1:nX)
           AY(:nGblock,(nY+1):(nY+nX)) = AXPW(:,1:nX) 
           BY(:nGblock,(nY+1):(nY+nX)) = BXPW(:,1:nX) 
           YXP((nGblock*2+1):(4*nGblock),(nY+1):(nY+nX)) = YXP((nGblock+1):(2*nGblock),(nY+1):(nY+nX))
           AY((nGblock*2+1):(4*nGblock),(nY+1):(nY+nX)) =  AY((nGblock+1):(2*nGblock),(nY+1):(nY+nX))
           BY((nGblock*2+1):(4*nGblock),(nY+1):(nY+nX)) =  BY((nGblock+1):(2*nGblock),(nY+1):(nY+nX))
       endif

       nY=nY+nX !permanently add the converged vectors to the hard-lock list
       !Deallocate the work spaces before the big Ritz
       do i = 1, size(psi)
           call deallocate_field(work(i),grids(work(i)%grid))
           call deallocate_field(Hpsi(i),grids(Hpsi(i)%grid))
           call deallocate_field(Spsi(i),grids(Spsi(i)%grid))
           call deallocate_field(psi(i),grids(psi(i)%grid))
       enddo
       deallocate(work)
       deallocate(Hpsi)
       deallocate(Spsi)
       deallocate(psi)

       deallocate(XPW)
       deallocate(Xwork)
       deallocate(Xwork2)
       deallocate(AXPW)
       deallocate(BXPW)
   enddo ! End loop on blocks

   deallocate(resnorm)
   deallocate(CX)
   deallocate(CP)
   deallocate(Kpc)
   if(numerics%lobpcg%PC_type.eq.2) then
       deallocate(Ekin)
       deallocate(G2_2)
   endif
   !The Big Ritz
   !print *, 'The Big Ritz-Pre-Ortho'
   call B_orthogonalize(YXP(:,1:nY), BY(:,1:nY), nY, parallel,  AY(:,1:nY), grid%gamma)
   allocate(Ywork(size(YXP,1),nY))
   !print *, 'The Big Ritz'
   call compute_Ritz_orthobasis(YXP(:,1:nY), AY(:,1:nY), nY,  &
    Ywork(:,1:nY), lam(1:nY), parallel, rot_A=.true., gamma=grid%gamma )
   !print *, 'The Big Ritz-Done'

   deallocate(Ywork)
   deallocate(AY)
   deallocate(BY)

   all_eigs(1:nY)=lam(1:nY)  
   deallocate(lam)
   do band_block=1, nband_block
       call allocate_orbital(orbitals(band_block), grids, orb_grid, parallel)
       orbitals(band_block)%eig(:)=all_eigs(orbitals(band_block)%band)
       !Need to adjust eig for spinor components if definable like in restricted open-shell collinear
       orbitals(band_block)%filter(:)=filter(:block_n_spinor(band_block),band_block)
       orbitals(band_block)%occ(:)=occ(:block_n_spinor(band_block),band_block)
       orbitals(band_block)%weight(:)=weight(:block_n_spinor(band_block),band_block)
   enddo
   deallocate(filter)
   deallocate(occ)
   deallocate(weight)
   deallocate(eigs)

   !print *, 'The Final Rotation'

   call orbitals_Matrix2FFT(orbitals(:nband_block), YXP(:,1:nY), grids, parallel, nX, &
                            grid%reduction_option)
   !print *, 'The Final Rotation-done'

   do band_block=1, nband_block
       do i=1, size(orbitals(band_block)%of(:))
           orbitals(band_block)%of(i)%G=orbitals(band_block)%of(i)%G*grid%cutwf
           call recip_to_real(orbitals(band_block)%of(i), grids)
       enddo
   enddo

   deallocate(YXP)

   !print *, 'LOBPCG cycle Done'
end subroutine

subroutine Rayleigh_Ritz_traditional(nX, nc, eW, &
                                    lam, XPW, AXPW, CX, CP, parallel, do_p, gamma)
    use lapack, only : zgemm, ztrttp, ztpttr, zherk
    use parallel_type, only : parallel_struct
    use parallel_mod,only : parallel_task, parallel_wait
    use constants, only : i_
    use linalg, only : eigh

    type(parallel_struct), intent(in) :: parallel

    integer, intent(in) :: nX, nc
    integer, intent(inout) :: eW
    real(dp), intent(inout) :: lam(:)
    complex(dp), intent(inout) :: XPW(:,:), AXPW(:,:)
    complex(dp), intent(out) :: CX(:,:),CP(:,:)
    logical, intent(in) :: do_p

    complex(dp), allocatable :: gramA(:,:), U(:,:), gramB(:,:), Oss(:,:) 
    integer :: ierr, nS, np, nG, nSt, i, j
    complex(dp), allocatable :: send(:)
    logical, intent(in) :: gamma
    
    nS=eW 
    nSt=nS*(nS+1)/2
    np=0
    if(do_p) np=nX-nc
    nG=size(XPW,1)

    allocate(gramA(nS,nS), gramB(nS,nS), U(nS,nS))
    allocate(send(nSt),Oss(nS,nS))

    Oss=0.0_dp
    call zgemm('C','N', nS, nS, nG, (1.0_dp,0.0_dp), XPW(:,1:nS),  &
    nG, AXPW(:,1:nS), nG, (0.0_dp,0.0_dp),  Oss, nS)

    call parallel_task('sum', Oss, parallel, 'k')

   ! if(parallel%myid.eq.0) print *, parallel%myid, 'Oss: '
    !do i=1, size(Oss,1)
     !   if(parallel%myid.eq.0) print *, Oss(i,:)
    !enddo

    do i=1, nS
        do j=i,nS
            Oss(j,i)=(Oss(j,i)+conjg(Oss(i,j)))*0.5_dp
            Oss(i,j)=conjg(Oss(j,i))
        enddo
    enddo

    !call ztrttp('U',nS, Oss, nS, send, ierr)
    !Oss=0.0_dp 
    !call parallel_task('sum', send, parallel, 'k')     
    !call ztpttr('U',nS, send, Oss, nS, ierr) 

    if(gamma) then
        gramA=real(Oss)
    else
        gramA=Oss
    endif
   ! if(parallel%myid.eq.0) print *, parallel%myid, 'gramA: '
   ! do i=1, size(gramA,1)
   !     if(parallel%myid.eq.0) print *, gramA(i,:)
    !enddo

   !  Oss=0.0_dp
   ! call zherk('U','C', nS, nG, (1.0_dp,0.0_dp), XPW(:,1:nS),  &
   ! nG, (0.0_dp,0.0_dp),  Oss, nS)
   ! call ztrttp('U',nS, Oss, nS, send, ierr) 
   ! Oss=0.0_dp
   ! call parallel_task('sum', send, parallel, 'k')        
   ! call ztpttr('U',nS, send, Oss, nS, ierr) 

    Oss=0.0_dp
    call zgemm('C','N', nS, nS, nG, (1.0_dp,0.0_dp), XPW(:,1:nS),  &
    nG, XPW(:,1:nS), nG, (0.0_dp,0.0_dp),  Oss, nS)

    call parallel_task('sum', Oss, parallel, 'k')

    do i=1, nS
        do j=i,nS
            Oss(j,i)=(Oss(j,i)+conjg(Oss(i,j)))*0.5_dp
            Oss(i,j)=conjg(Oss(j,i))
        enddo
    enddo

    if(gamma) then
        gramB=real(Oss)
    else
        gramB=Oss
    endif   
     !if(parallel%myid.eq.0) print *, parallel%myid, 'gramB: ', gramB

    deallocate(send,Oss)

    call eigh(gramA, gramB, lam(1:nS), U)

    ! if(parallel%myid.eq.0) print *, parallel%myid, 'U: ', U

   !if(parallel%myid.eq.0) print *, parallel%myid, 'lam: ', lam(1:nS)

    !if(parallel%myid.eq.0) stop
    CX(1:nS,1:nX)=U(1:nS, 1:nX)
    CP=0.0_dp
    CP(1:(nS-nX),1:nX)=U((nX+1):nS, 1:nX)
    
    deallocate(gramA,gramB,U)
end subroutine

subroutine Gen_Rayleigh_Ritz_traditional(nX, nc, eW, &
    lam, XPW, AXPW, BXPW, CX, CP, parallel, do_p, gamma)
use lapack, only : zgemm, ztrttp, ztpttr
use parallel_type, only : parallel_struct
use parallel_mod,only : parallel_task, parallel_wait

use linalg, only : eigh

type(parallel_struct), intent(in) :: parallel

integer, intent(in) :: nX, nc
integer, intent(inout) :: eW
real(dp), intent(inout) :: lam(:)
complex(dp), intent(inout) :: XPW(:,:), AXPW(:,:), BXPW(:,:)
complex(dp), intent(out) :: CX(:,:),CP(:,:)
logical, intent(in) :: do_p, gamma

complex(dp), allocatable :: gramA(:,:), U(:,:), gramB(:,:), Oss(:,:) 
integer :: ierr, nS, np, nG, nSt, i, j
complex(dp), allocatable :: send(:)

nS=eW 
nSt=nS*(nS+1)/2
np=0
if(do_p) np=nX-nc
nG=size(XPW,1)

allocate(gramA(nS,nS), gramB(nS,nS), U(nS,nS))
allocate(send(nSt),Oss(nS,nS))

Oss=0.0_dp
call zgemm('C','N', nS, nS, nG, (1.0_dp,0.0_dp), XPW(:,1:nS),  &
nG, AXPW(:,1:nS), nG, (0.0_dp,0.0_dp),  Oss, nS)

!call ztrttp('U',nS, Oss, nS, send, ierr) 
!Oss=0.0_dp     
!call parallel_task('sum', send, parallel, 'k')
!call ztpttr('U',nS, send, Oss, nS, ierr) 
call parallel_task('sum', Oss, parallel, 'k')

do i=1, nS
    do j=i,nS
        Oss(j,i)=(Oss(j,i)+conjg(Oss(i,j)))*0.5_dp
        Oss(i,j)=conjg(Oss(j,i))
    enddo
enddo

if(gamma) then
    gramA=real(Oss)
else
    gramA=Oss
endif

!if(parallel%myid.eq.0) print *, parallel%myid, 'gramA: '
!do i=1, size(gramA,1)
!    if(parallel%myid.eq.0) print *, real(gramA(i,:))
!enddo

Oss=0.0_dp
call zgemm('C','N', nS, nS, nG, (1.0_dp,0.0_dp), XPW(:,1:nS),  &
nG, BXPW(:,1:nS), nG, (0.0_dp,0.0_dp),  Oss, nS)

!call ztrttp('U',nS, Oss, nS, send, ierr) 
!Oss=0.0_dp       
!call parallel_task('sum', send, parallel, 'k') 
!call ztpttr('U',nS, send, Oss, nS, ierr) 
call parallel_task('sum', Oss, parallel, 'k')

do i=1, nS
    do j=i,nS
        Oss(j,i)=(Oss(j,i)+conjg(Oss(i,j)))*0.5_dp
        Oss(i,j)=conjg(Oss(j,i))
    enddo
enddo

if(gamma) then
    gramB=real(Oss)
else
    gramB=Oss
endif  

deallocate(send,Oss)

call eigh(gramA, gramB, lam(1:nS), U)

!if(parallel%myid.eq.0) print *, parallel%myid, 'gramB: '
!do i=1, size(gramB,1)
!    if(parallel%myid.eq.0) print *, real(gramB(i,:))
!enddo


!if(parallel%myid.eq.0) print *, parallel%myid, 'U: ', U(1:nS,1:nS)

!if(parallel%myid.eq.0) print *, parallel%myid, 'lam: ', lam(1:nS)

CX(1:nS,1:nX)=U(1:nS, 1:nX)
CP=0.0_dp
CP(1:(nS-nX),1:nX)=U((nX+1):nS, 1:nX)

deallocate(gramA,gramB,U)
end subroutine

    subroutine compute_Ritz_orthobasis(Xin, AX, nX, Xwork, lam, parallel, rot_A, BX, gamma)
        use parallel_type, only : parallel_struct
        use parallel_mod,only : parallel_task

        use linalg, only : eigh
        use lapack, only : zgemm

        type(parallel_struct), intent(in) :: parallel
        complex(dp), intent(inout) :: Xin(:,:) 
        complex(dp), intent(out) :: Xwork(:,:) 
        complex(dp), intent(inout) :: AX(:,:)
        complex(dp), intent(inout), optional :: BX(:,:)
        logical, intent(in) :: rot_A
        real(dp), intent(out) :: lam(:) 
        integer, intent(in) :: nX
        integer :: i, j
        complex(dp) :: Oxx(nX,nX), U(nX,nX)
        logical, intent(in) :: gamma
        !Form Oxx= <X|AX> = <X|H|X>
        call zgemm('C','N', nX, nX, size(Xin,1), (1.0_dp,0.0_dp), Xin,  &
         size(Xin,1), AX, size(AX,1), (0.0_dp,0.0_dp),  Oxx, size(Oxx,1))

        call parallel_task('sum', Oxx, parallel, 'k')

        do i=1, nX
            do j=i,nX
                Oxx(j,i)=(Oxx(j,i)+conjg(Oxx(i,j)))*0.5_dp
                Oxx(i,j)=conjg(Oxx(j,i))
            enddo
        enddo
        if(gamma) Oxx=real(Oxx)
        
        call eigh(Oxx,lam,U)

        !U is now the Ritz vectors lam is the ritz
        !|X> = |Xin>*Qxx
        call zgemm('N','N', size(Xin,1), nX, nX, (1.0_dp,0.0_dp), Xin,  &
        size(Xin,1), U, size(U,1), (0.0_dp,0.0_dp),  Xwork, size(Xwork,1))
        Xin=Xwork

        if(rot_A) then
            !|AX> = |AX>*Qxx
            call zgemm('N','N', size(AX,1), nX, nX, (1.0_dp,0.0_dp), AX,  &
                        size(AX,1), U, size(U,1), (0.0_dp,0.0_dp),  Xwork, size(Xwork,1))
            AX=Xwork
        endif
        if(present(BX)) then
            !|BX> = |BX>*Qxx
            call zgemm('N','N', size(BX,1), nX, nX, (1.0_dp,0.0_dp), BX,  &
                        size(BX,1), U, size(U,1), (0.0_dp,0.0_dp),  Xwork, size(Xwork,1))
            BX=Xwork
        endif
    end subroutine

    subroutine av(X, nX, AX, parallel, grids, g, potentials, atoms, elements, psi, Hpsi)
        use parallel_type, only : parallel_struct
        use Apply_Hamiltonian, only: Apply_H, Apply_S
        use odp_type, only : field_struct
        use FFT2Matrix, only : field_FFT2Matrix, field_Matrix2FFT
        use grids_type, only : grid_struct
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: real_to_recip, recip_to_real
        use parallel_mod, only : parallel_wait, parallel_task


        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(field_struct), intent(inout) :: psi(:), Hpsi(:)

        complex(dp), intent(in), contiguous :: X(:,:) 
        complex(dp), intent(inout), contiguous :: AX(:,:) 
        integer, intent(in) :: g, nX
        integer :: s
        logical :: test_S=.false.
        grid=>grids(g)

        call field_Matrix2FFT(psi, X(:,1:nX), grid, parallel, nX, grid%reduction_option)

        if(test_S) then
            call Apply_S(psi, Hpsi, grids, atoms, elements, parallel, calc_R=.true.)
            do s=1,size(psi)
                Hpsi(s)%G=Hpsi(s)%G-psi(s)%G
            enddo 
        else
            call Apply_H(psi, Hpsi, grids, potentials, atoms, elements, parallel, calc_R=.true.)
        endif
        
        call field_FFT2Matrix(Hpsi, AX(:,1:nX), grid, parallel, nX, grid%reduction_option)

    end subroutine

    subroutine bv(X, nX, BX, parallel, grids, g, atoms, elements, psi, Spsi)
        use parallel_type, only : parallel_struct
        use Apply_Hamiltonian, only: Apply_S
        use odp_type, only : field_struct
        use FFT2Matrix, only : field_FFT2Matrix, field_Matrix2FFT
        use grids_type, only : grid_struct
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: real_to_recip, recip_to_real
        use parallel_mod, only : parallel_wait, parallel_task
        use Non_Local_ion, only: Calculate_Projector_overlaps

        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(field_struct), intent(inout) :: psi(:), Spsi(:)

        complex(dp), intent(in), contiguous :: X(:,:) 
        complex(dp), intent(inout), contiguous :: BX(:,:)

        integer, intent(in) :: g, nX
        grid=>grids(g)

        call field_Matrix2FFT(psi, X(:,1:nX), grid, parallel, nX, grid%reduction_option)

        call Apply_S(psi, Spsi, grids, atoms, elements, parallel, calc_R=.true.)
       
        call field_FFT2Matrix(Spsi, BX(:,1:nX), grid, parallel, nX, grid%reduction_option)

    end subroutine

    subroutine av_bv(X, nX, AX, BX, parallel, grids, g, potentials, atoms, elements, psi, Hpsi, Spsi)
        use parallel_type, only : parallel_struct
        use Apply_Hamiltonian, only: Apply_S, Apply_H
        use odp_type, only : field_struct
        use FFT2Matrix, only : field_FFT2Matrix, field_Matrix2FFT
        use grids_type, only : grid_struct
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: real_to_recip, recip_to_real
        use parallel_mod, only : parallel_wait, parallel_task
        use Non_Local_ion, only: Calculate_Projector_overlaps

        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(field_struct), intent(inout) :: psi(:), Spsi(:), Hpsi(:)
        type(potentials_struct), intent(in) :: potentials

        complex(dp), intent(in), contiguous :: X(:,:) 
        complex(dp), intent(inout), contiguous :: AX(:,:), BX(:,:)
        complex(dp), allocatable :: PYPsi(:,:,:)

        integer, intent(in) :: g, nX
        integer :: s
        grid=>grids(g)

        call field_Matrix2FFT(psi, X(:,1:nX), grid, parallel, nX, grid%reduction_option)

        do s=1,size(psi)
            call recip_to_real(psi(s), grids)
        enddo

        allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(psi)))
        PYpsi=0.0_dp
        do s=1, size(psi)
            call Calculate_Projector_overlaps(psi(s),  &
            PYpsi(:,:,s), atoms, elements, grid, parallel) 
        enddo
            

        call Apply_S(psi, Spsi, grids, atoms, elements, parallel, calc_R=.false.,PYpsi_in=PYpsi)
        call Apply_H(psi, Hpsi, grids, potentials, atoms, elements, parallel, calc_R=.false.,PYpsi_in=PYpsi)
       
        call field_FFT2Matrix(Spsi, BX(:,1:nX), grid, parallel, nX, grid%reduction_option)
        call field_FFT2Matrix(Hpsi, AX(:,1:nX), grid, parallel, nX, grid%reduction_option)

        deallocate(PYpsi)

    end subroutine

    subroutine Orbitals_Project_out(orbitals_known, orbitals_ort, grids, parallel, &
                                     atoms, elements, all_PAW, filter)
        use lapack, only : zgemm
        use parallel_mod,only : parallel_task
        use grids_type, only : grid_struct
        use odp_type, only : field_struct, orbital_struct
        use FFT2Matrix, only : orbitals_FFT2Matrix, orbitals_Matrix2FFT
        use parallel_type, only : parallel_struct
        use FFT, only : recip_to_real
        use Apply_Hamiltonian, only: Apply_S
        use simulation_type, only : all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use odp, only: allocate_field, deallocate_field
        use FFT2Matrix, only : field_FFT2Matrix
        use operations_3D, only : integrate_3D_R
        type(parallel_struct), intent(in) :: parallel
        type(orbital_struct), intent(inout) :: orbitals_known(:), orbitals_ort(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        real(dp), intent(in), optional :: filter(:)
        type(field_struct), allocatable :: Sphi(:)
        type(all_PAW_struct), intent(in) :: all_PAW
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)

        complex(dp), allocatable :: Y1(:,:), Y2(:,:), Y3(:,:), O12(:,:)
        integer, allocatable :: block_n_spinor(:)

        integer :: orb_grid, nX, nY1, nY2, band_block, i, nG, k
        real(dp) :: olap

        nX=parallel%n_band_groups

        orb_grid=orbitals_known(1)%of(1)%grid
        do i=1, size(orbitals_known)
            if(orbitals_known(i)%of(1)%grid.ne.orb_grid) then
                print *, 'All orbitals passed into Orthogonalization need to be on the same grid (same k-point)'
            endif
        enddo

        do i=1,size(orbitals_ort)
            if(orbitals_ort(i)%of(1)%grid.ne.orb_grid) then
                print *, 'All orbitals passed into Orthogonalization need to be on the same grid (same k-point)'
            endif
        enddo
        
        grid=>grids(orb_grid)
        if(grid%reduction_option.eq.2) then
            !if(parallel%myid.eq.0) print *, 'Removing some of the zeros in x and y directions'
            nG=product(grid%Ng_small)
        else
            nG=product(grid%Ng_local)
        endif
        if(mod(nG,nX).ne.0) nG=(nG/nX+1)*nX

        nY1=nX*size(orbitals_known)
        nY2=nX*size(orbitals_ort)

        allocate(Y1(maxval(orbitals_known(:)%n_spinor)*nG/nX,nY1))

        if((all_PAW%N_PAW_atoms.gt.0)) then 
            allocate(Y3(maxval(orbitals_known(:)%n_spinor)*nG/nX,nY1))

            allocate(block_n_spinor(size(orbitals_known)))
            block_n_spinor(:)=orbitals_known(:)%n_spinor
            nY1=0
            do band_block=1, size(orbitals_known)
                allocate(Sphi(block_n_spinor(band_block)))
    
                do i=1, block_n_spinor(band_block)
                    call allocate_field(Sphi(i), grid, parallel)
                enddo
                Sphi%grid=orb_grid
                call Apply_S(orbitals_known(band_block)%of(:), Sphi, grids, atoms, elements, parallel, calc_R=.false.)
                call field_FFT2Matrix(Sphi, Y1(:,(nY1+1):(nY1+nX)), grid, parallel, nX, grid%reduction_option)
    
                do i = 1, size(Sphi)
                    call deallocate_field(Sphi(i),grids(Sphi(i)%grid))
                enddo
                deallocate(Sphi)
                nY1=nY1+nX
            enddo
            call orbitals_FFT2Matrix(orbitals_known(:), Y3(:,:), grids, parallel, nX, &
                grids(orb_grid)%reduction_option)
        else
            call orbitals_FFT2Matrix(orbitals_known(:), Y1(:,:), grids, parallel, nX, &
                grids(orb_grid)%reduction_option)
        endif

        allocate(Y2(maxval(orbitals_ort(:)%n_spinor)*nG/nX,nY2))

        call orbitals_FFT2Matrix(orbitals_ort(:), Y2(:,:), grids, parallel, nX, &
                                 grids(orb_grid)%reduction_option)

        if((all_PAW%N_PAW_atoms.gt.0)) then 
            if(present(filter)) then
                call projection_removal_nonorthobasis(Y2, nY2, Y3, Y1, nY1, parallel, filter, gamma=grid%gamma)
            else
                call projection_removal_nonorthobasis(Y2, nY2, Y3, Y1, nY1, parallel, gamma=grid%gamma)
            endif
        else
            if(present(filter)) then
                call projection_removal_orthobasis(Y2, nY2, Y1, nY1, parallel, filter, gamma=grid%gamma)
            else
                call projection_removal_orthobasis(Y2, nY2, Y1, nY1, parallel, gamma=grid%gamma)
            endif
        endif

        call orbitals_Matrix2FFT(orbitals_ort(:), Y2(:,:), grids, parallel, nX, &
            grid%reduction_option)
        
        do band_block=1, size(orbitals_ort)
            do i=1, size(orbitals_ort(band_block)%of(:))
                orbitals_ort(band_block)%of(i)%G=orbitals_ort(band_block)%of(i)%G*grid%cutwf
                call recip_to_real(orbitals_ort(band_block)%of(i), grids)
            enddo
        enddo

    end subroutine

    subroutine projection_removal_orthobasis(X, nX, Y, nY, parallel, filter, gamma)
        use lapack, only : zgemm
        use parallel_mod,only : parallel_task
        use parallel_type, only : parallel_struct
        type(parallel_struct), intent(in) :: parallel
        complex(dp), intent(inout) :: X(:,:) 
        complex(dp), intent(in) :: Y(:,:) 
        integer, intent(in) :: nY
        integer, intent(in) :: nX
        complex(dp), allocatable :: Oyx(:,:)
        real(dp), intent(in), optional :: filter(:)
        logical, intent(in):: gamma
        integer :: i

                allocate(Oyx(nY,nX)) 

                !Project out the Y
                !overlap with  Y vectors <Y|X>
                call zgemm('C','N', nY, nX, size(Y,1), (1.0_dp,0.0_dp), Y,  &
                            size(Y,1), X, size(X,1), (0.0_dp,0.0_dp),  Oyx, size(Oyx,1))

                call parallel_task('sum', Oyx, parallel, 'k')
                if(gamma) Oyx=real(Oyx)

                if(present(filter)) then
                    do i=1, size(Oyx,1)
                        Oyx(i,:)=filter(i)*Oyx(i,:)
                    enddo
                endif
                ! |X> = |X> - |Y><Y|X>
                call zgemm('N','N', size(X,1), size(Oyx,2), nY, -(1.0_dp,0.0_dp), Y,  &
                            size(Y,1), Oyx, size(Oyx,1), (1.0_dp,0.0_dp),  X, size(X,1))
                
                deallocate(Oyx) 
           
    end subroutine

    subroutine projection_removal_nonorthobasis(X, nX, Y, BY, nY, parallel, filter, gamma)
        use lapack, only : zgemm
        use parallel_mod,only : parallel_task
        use parallel_type, only : parallel_struct
        use linalg, only: solve
        type(parallel_struct), intent(in) :: parallel
        complex(dp), intent(inout) :: X(:,:) 
        complex(dp), intent(in) :: Y(:,:), BY(:,:)
        integer, intent(in) :: nY
        integer, intent(in) :: nX
        complex(dp), allocatable :: Oyx(:,:), M2(:,:)
        real(dp), intent(in), optional :: filter(:)
        logical, intent(in):: gamma
        integer :: i

                allocate(Oyx(nY,nX)) 
                !Project out the Y
                !overlap with  Y vectors <Y|B|X>
                call zgemm('C','N', nY, nX, size(BY,1), (1.0_dp,0.0_dp), BY,  &
                            size(BY,1), X, size(X,1), (0.0_dp,0.0_dp),  Oyx, size(Oyx,1))
                call parallel_task('sum', Oyx, parallel, 'k')

                if(gamma) Oyx=real(Oyx)

            !    allocate(Oyy(nY,nY)) 
                !overlap with  Y vectors <Y|B|Y>
             !   call zgemm('C','N', nY, nY, size(BY,1), (1.0_dp,0.0_dp), BY,  &
             !               size(BY,1), Y, size(Y,1), (0.0_dp,0.0_dp),  Oyy, size(Oyy,1))
             !   call parallel_task('sum', Oyy, parallel, 'k')
                
                !M2=[<Y|B|Y>]^-1<Y|B|X>
                !do i=1,nX
                !    M2(:,i)=solve(Oyy(:nY,:nY),Oyx(:,i))
                !enddo
                !Oyx=M2
                !deallocate(Oyy)
                !deallocate(M2)
                
                if(present(filter)) then
                    do i=1, size(Oyx,1)
                        Oyx(i,:)=filter(i)*Oyx(i,:)
                    enddo
                endif

                ! |X> = |X> - |Y>[<Y|B|Y>]^-1<Y|B|X>
                call zgemm('N','N', size(X,1), size(Oyx,2), nY, -(1.0_dp,0.0_dp), Y,  &
                            size(Y,1), Oyx, size(Oyx,1), (1.0_dp,0.0_dp),  X, size(X,1))
                
                deallocate(Oyx) 
           
    end subroutine

    subroutine Orbitals_Orthogonalize(orbitals, grids, parallel, &
        atoms, elements, all_PAW)
        use grids_type, only : grid_struct
        use odp_type, only : field_struct, orbital_struct
        use FFT2Matrix, only : orbitals_FFT2Matrix, orbitals_Matrix2FFT
        use parallel_type, only : parallel_struct
        use FFT, only : recip_to_real
        use Apply_Hamiltonian, only: Apply_S
        use simulation_type, only : all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use odp, only: allocate_field, deallocate_field
        use FFT2Matrix, only : field_FFT2Matrix

        type(parallel_struct), intent(in) :: parallel
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        type(field_struct), allocatable :: Sphi(:)
        type(all_PAW_struct), intent(in) :: all_PAW
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)

        complex(dp), allocatable :: Y(:,:), BY(:,:)
        integer, allocatable :: block_n_spinor(:)

        integer :: orb_grid, nX, nY, band_block, i, nG

        nX=parallel%n_band_groups

        orb_grid=orbitals(1)%of(1)%grid
        do i=1, size(orbitals)
        if(orbitals(i)%of(1)%grid.ne.orb_grid) then
        print *, 'All orbitals passed into Orthogonalization need to be on the same grid (same k-point)'
        endif
        enddo

        grid=>grids(orb_grid)
        if(grid%reduction_option.eq.2) then
        !if(parallel%myid.eq.0) print *, 'Removing some of the zeros in x and y directions'
        nG=product(grid%Ng_small)
        else
        nG=product(grid%Ng_local)
        endif
        if(mod(nG,nX).ne.0) nG=(nG/nX+1)*nX

        nY=nX*size(orbitals)

        allocate(Y(maxval(orbitals(:)%n_spinor)*nG/nX,nY))
        if(all_PAW%N_PAW_atoms.gt.0) then 
            allocate(BY(maxval(orbitals(:)%n_spinor)*nG/nX,nY))
            allocate(block_n_spinor(size(orbitals)))
            block_n_spinor(:)=orbitals(:)%n_spinor
            nY=0
            do band_block=1, size(orbitals)
                allocate(Sphi(block_n_spinor(band_block)))
                do i=1, block_n_spinor(band_block)
                    call allocate_field(Sphi(i), grid, parallel)
                enddo
                Sphi%grid=orb_grid
                call Apply_S(orbitals(band_block)%of(:), Sphi, grids, atoms, elements, parallel, calc_R=.false.)
                call field_FFT2Matrix(Sphi, BY(:,(nY+1):(nY+nX)), grid, parallel, nX, grid%reduction_option)

                do i = 1, size(Sphi)
                    call deallocate_field(Sphi(i),grids(Sphi(i)%grid))
                enddo
                deallocate(Sphi)
                nY=nY+nX
            enddo
        endif

        call orbitals_FFT2Matrix(orbitals(:), Y(:,:), grids, parallel, nX, &
            grids(orb_grid)%reduction_option)

        if(all_PAW%N_PAW_atoms.gt.0) then 
            call B_orthogonalize(Y, BY, nY, parallel, gamma=grid%gamma)
            deallocate(BY)
        else
            call orthogonalize(Y, nY, parallel, gamma=grid%gamma)
        endif

        call orbitals_Matrix2FFT(orbitals(:), Y(:,:), grids, parallel, nX, &
            grid%reduction_option)
end subroutine

    subroutine orthogonalize(Xin, nX, parallel, vec, gamma)
        use lapack, only : zpotrf, zgemm, ztrsm
        use parallel_mod,only : parallel_task
        use parallel_type, only : parallel_struct

        type(parallel_struct), intent(in) :: parallel
        integer, intent(in) :: nX
        complex(dp), intent(inout) :: Xin(:,:) 
        complex(dp), intent(inout), optional :: vec(:,:) 
        complex(dp) :: Oxx(nX,nX)
        logical, intent(in) :: gamma
        integer :: ierr, i, j
       
        !overlap of vectors <X|X>
        call zgemm('C','N', nX, nX, size(Xin,1), (1.0_dp,0.0_dp), Xin,  &
                    size(Xin,1), Xin, size(Xin,1), (0.0_dp,0.0_dp),  Oxx, size(Oxx,1))
        call parallel_task('sum', Oxx, parallel, 'k')

        do i=1,nX; do j=i,nX
            Oxx(i,j)=(Oxx(i,j) + conjg(Oxx(j,i)))*0.5_dp
            Oxx(j,i)=conjg(Oxx(i,j))
        enddo; enddo
        if(gamma) Oxx=real(Oxx)

        !O=chol(O), X=X*O^-1(psi2->orthonormal) (This is repeated on every rank, can it be parallelized with scalapck?)
        call zpotrf('U', nX, Oxx, nX, ierr)
        call ztrsm('R', 'U', 'N', 'N', size(Xin,1), nX, (1.0_dp,0.0_dp), Oxx, nX, Xin, size(Xin,1)) 

        !V=V*O^-1
        if(present(vec)) then
            call ztrsm('R', 'U', 'N', 'N', size(vec,1), nX, (1.0_dp,0.0_dp), Oxx, nX, vec, size(vec,1)) 
        endif
        
    end subroutine 

    subroutine B_orthogonalize(Xin, BXin, nX, parallel, vec, gamma)
        use lapack, only : zpotrf, zgemm, ztrsm
        use parallel_mod,only : parallel_task
        use parallel_type, only : parallel_struct

        type(parallel_struct), intent(in) :: parallel
        integer, intent(in) :: nX
        complex(dp), intent(inout) :: Xin(:,:), BXin(:,:)  
        complex(dp), intent(inout), optional :: vec(:,:) 
        complex(dp) :: Oxx(nX,nX)
        integer :: ierr, i, j
        logical, intent(in) :: gamma
       
        !overlap of vectors <X|B|X>
        call zgemm('C','N', nX, nX, size(Xin,1), (1.0_dp,0.0_dp), Xin,  &
                    size(Xin,1), BXin, size(BXin,1), (0.0_dp,0.0_dp),  Oxx, size(Oxx,1))
        call parallel_task('sum', Oxx, parallel, 'k')

        do i=1,nX; do j=i,nX
            Oxx(i,j)=(Oxx(i,j) + conjg(Oxx(j,i)))*0.5_dp
            Oxx(j,i)=conjg(Oxx(i,j))
        enddo; enddo
        
        if(gamma) Oxx=real(Oxx)

      !  do i=1,nX
      !      if(parallel%myid_band.eq.0) print *, i, Oxx(i,:)
      !  enddo

        !O=chol(O), X=X*O^-1(psi2->orthonormal) (This is repeated on every rank, can it be parallelized with scalapck?)
        call zpotrf('U', nX, Oxx, nX, ierr)
               ! if(myid.eq.0) print *, 'R: '
               ! do i=1, nX
               ! if(myid.eq.0) print *, Oxx(i,:)
               ! enddo
        call ztrsm('R', 'U', 'N', 'N', size(Xin,1), nX, (1.0_dp,0.0_dp), Oxx, nX, Xin, size(Xin,1))        
        !BX=BX*O^-1
        call ztrsm('R', 'U', 'N', 'N', size(BXin,1), nX, (1.0_dp,0.0_dp), Oxx, nX, BXin, size(BXin,1))     
        !V=V*O^-1
        if(present(vec)) then
            call ztrsm('R', 'U', 'N', 'N', size(vec,1), nX, (1.0_dp,0.0_dp), Oxx, nX, vec, size(vec,1)) 
        endif
                
    end subroutine


end module
