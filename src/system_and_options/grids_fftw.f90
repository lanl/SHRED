module grids_mod
    use types, only : dp

    implicit none
    
    public :: setup_grids, map_between_distributed_grids, Gspace_grid_to_grid_transfer, allocate_local_fields_G, &
              allocate_full_fields_G, allocate_local_fields_R, allocate_full_fields_R


    interface allocate_local_fields_G
        module procedure allocate_local_field0_G
        module procedure allocate_local_field1_G
        module procedure allocate_local_field2_G
        module procedure allocate_local_field0_G_real
        module procedure allocate_local_field1_G_real
    end interface

    interface allocate_full_fields_G
        module procedure allocate_full_field0_G
    end interface

    interface allocate_local_fields_R
        module procedure allocate_local_field0_R_real
        module procedure allocate_local_field1_R_real
        module procedure allocate_local_field2_R_real
        module procedure allocate_local_field0_R_complex
        module procedure allocate_local_field1_R_complex
        module procedure allocate_local_field2_R_complex
    end interface

    interface allocate_full_fields_R
        module procedure allocate_full_field0_R
    end interface

    contains

    subroutine setup_grids(grid, mpi)
        use, intrinsic :: iso_c_binding
        use constants, only : pi
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use parallel_mod, only : parallel_wait 
        include 'fftw3-mpi.f03'

        type(grid_struct), intent(inout) :: grid
        type(parallel_struct), intent(in) :: mpi
        real(dp), allocatable :: X(:,:,:)

        !For FFTW
        integer(C_INTPTR_T) :: size_x, size_y, size_z, local_z, local_j_offset
        
        integer :: ijk_global(3), i,j,k

        !FFTW using full grids, same order in real and recip space
        grid%Nr=grid%Ng

        grid%FFT_comm=mpi%comm_bands
        grid%nproc_FFT=mpi%nproc_band

        !Setup FFTW w/ MPI initialization 
            size_x=grid%Ng(1)
            size_y=grid%Ng(2)
            size_z=grid%Ng(3)
            local_j_offset=0
            local_z=size_z
            call fftw_mpi_init() !figure out what this does
           
        !Adjust the grid to conform to FFTW requirements, number of z-direction slices must be a product of the
        ! processors available for distribution (nproc_bands), on the communicator comm_bands 
            if(mod(int(size_z),mpi%nproc_band).ne.0) then
                if(mpi%myid.eq.0) print *, 'Adjusting # of z direction grid points to divide &
                &    by available processors, Ecut unchanged)'
                size_z=(size_z/mpi%nproc_band)*mpi%nproc_band + mpi%nproc_band
                grid%Ng(3)=int(size_z)
                grid%Nr(3)=int(size_z)
                if(mpi%myid.eq.0) print *, "Ng_z:,",  grid%Ng(3)
            endif

        
        !FFTW setup distributed FFT's
            grid%alloc_local = fftw_mpi_local_size_3d(size_z, size_y, size_x,  mpi%comm_bands, &
            local_z, local_j_offset)
            grid%alloc_local_R=grid%alloc_local
            grid%alloc_local_G=grid%alloc_local
    
            if(mpi%myid.eq.0) then
                if((size_y.ne.size_z).and.(size_z.eq.mpi%nproc_band)) then
                        print *, 'warning: If you are running on intel Cluster FFT w/ FFTW3 interface'
                        print *, 'you are probably about to crash due to intel bug when Ny.ne.Nz and Nz=nprocs_fft'
                        print *, 'consider adjusting grid parameters/ number cores, '
                        print *, 'rotating system (if Nx.ne.Ny), or running on gcc mkl openmpi'
                endif
            endif
            call parallel_wait(mpi)

            grid%Ng_local(1)=int(size_x)
            grid%Ng_local(2)=int(size_y)
            grid%Ng_local(3)=int(local_z)

            grid%myxyz(1)=0
            grid%myxyz(2)=0
            grid%myxyz(3)=int(local_j_offset/grid%Ng_local(3))
            
            !Full FFT, grid not pruned in G space
            grid%Nr_local=grid%Ng_local
            grid%myxyz_r=grid%myxyz

            !print *, 'local:', grid%Ng_local, grid%Nr_local


        !Allocating and setting up work space sets plan for all FFT's on this grid
            grid%work%p_R=fftw_alloc_complex(grid%alloc_local)
            call c_f_pointer(grid%work%p_R, grid%work%R, [size_x,size_y,local_z])
            grid%work%p_G=fftw_alloc_complex(grid%alloc_local)
            call c_f_pointer(grid%work%p_G, grid%work%G, [size_x,size_y,local_z])
           
        !Create the FFTW3 plans, any functions that are defined on this grid can use these same plans
            grid%plan_forward = &
                fftw_mpi_plan_dft_3d(size_z, size_y, size_x, grid%work%R, grid%work%G, mpi%comm_bands, &
                                FFTW_FORWARD, FFTW_MEASURE)
            grid%plan_reverse = &
                fftw_mpi_plan_dft_3d(size_z, size_y, size_x, grid%work%G, grid%work%R, mpi%comm_bands, &
                                FFTW_BACKWARD, FFTW_MEASURE)

        !Start defining the physical grid (real & reciprocal space values, both locally (processor dependent) &
        ! and Globaly (full grid)
        allocate(grid%G(grid%Ng_local(1),grid%Ng_local(2),grid%Ng_local(3),3))
        allocate(grid%R(grid%Ng_local(1),grid%Ng_local(2),grid%Ng_local(3),3))
        allocate(grid%G2(grid%Ng_local(1),grid%Ng_local(2),grid%Ng_local(3)))
        if(grid%ecutsm.gt.0.0_dp) then
            allocate(grid%p(grid%Ng_local(1),grid%Ng_local(2),grid%Ng_local(3)))
            allocate(grid%pprime(grid%Ng_local(1),grid%Ng_local(2),grid%Ng_local(3)))
            allocate(grid%pdoubleprime(grid%Ng_local(1),grid%Ng_local(2),grid%Ng_local(3)))
            allocate(X(grid%Ng_local(1),grid%Ng_local(2),grid%Ng_local(3)))
        endif

        grid%dG(:)=2.0_dp*pi/grid%Box_Length(:)
        grid%dR(:)=grid%Box_Length(:)/grid%Ng(:)
        grid%k(:)=grid%k(:)*grid%dG(:)
        
        if(mpi%myid.eq.(mpi%nproc-1))print *, 'Real Spacing (bohr): ', grid%dR(:) 
        if(mpi%myid.eq.(mpi%nproc-1))print *, 'Reciprocal Spacing (1/bohr): ', grid%dG(:) 
        do k = 1, grid%Ng_local(3); do j = 1, grid%Ng_local(2); do i = 1, grid%Ng_local(1)
            ijk_global = [i, j, k] + grid%myxyz(:)*grid%Ng_local(:)
            grid%R(i,j,k,:) = (ijk_global-1) * grid%dR(:)
            grid%G(i,j,k,:) = grid%dG(:) * (ijk_global - 1 - grid%Ng(:)*nint((ijk_global-1.5_dp)/grid%Ng)) + grid%k(:)
            grid%G2(i,j,k)  = sum(grid%G(i,j,k,:)**2)
        end do; end do; end do
        if(grid%ecutsm.gt.0.0_dp) then
            X=(grid%Ecut-0.5_dp*grid%G2)/grid%ecutsm
            where(X.gt.1.0_dp)
                grid%p=1.0_dp
                grid%pprime=0.0_dp
                grid%pdoubleprime=0.0_dp
            elsewhere(X.lt.tiny(1.0_dp)) !p should be infinite at X, but 
                grid%p=0.0_dp
                grid%pprime=0.0_dp
                grid%pdoubleprime=0.0_dp
            elsewhere
                grid%p=1.0_dp/(X**2*(3.0_dp+X*(1.0_dp+X*(-6.0_dp+3.0_dp*X))))
                grid%pprime=-3.0*(X-1.0_dp)**2*X*(2.0_dp+5.0_dp*X)*grid%p**2
                grid%pdoubleprime=6.0_dp*X**2* &
                    (9.0_dp+X*(8.0_dp+X*(-52.0_dp+X*(-3.0_dp+X*(137.0_dp+X*(-144.0+45.0_dp*X))))))*grid%p**3
            endwhere
            deallocate(X)
        endif

        !Set up the Energy cutoff of the plane wave orbital & density/potential coefficients
        allocate(grid%cutwf(grid%Ng_local(1),grid%Ng_local(2),grid%Ng_local(3)))
        allocate(grid%cutden(grid%Ng_local(1),grid%Ng_local(2),grid%Ng_local(3)))

        where(grid%G2.le.(grid%Ecut*2.0_dp))
            grid%cutwf = 1.0_dp
        elsewhere
            grid%cutwf = 0.0_dp
        end where
        where(grid%G2.le.(grid%Ecut*8.0_dp))
            grid%cutden = 1.0_dp
        elsewhere
            grid%cutden = 0.0_dp
        end where

        !Vector potential 
        grid%A =0.0_dp

        !set up the reduced grid removing 0's for the wave-function (those compoenent already set to zero by cutwf)
        grid%Ng_small(1)=grid%Ng_local(1)/2 + mod(grid%Ng_local(1),2)
        grid%Ng_small(2)=grid%Ng_local(2)/2 + mod(grid%Ng_local(2),2)
        grid%Ng_small(3)=grid%Ng_local(3)

        grid%reduction_option=2
    end subroutine

    

    subroutine setup_grids_new(grid, mpi, pruned)
        use, intrinsic :: iso_c_binding
        use constants, only : pi
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use parallel_mod, only : parallel_wait 
        include 'fftw3-mpi.f03'

        type(grid_struct), intent(inout) :: grid
        type(parallel_struct), intent(in) :: mpi
        real(dp), allocatable :: X(:,:,:)
        logical, intent(in) :: pruned

        !For FFTW
        integer(C_INTPTR_T) :: size_x, size_y, size_z, local_x, local_z, size_tmp
        
        integer :: ijk_global(3), i,j,k, rank, cut

        !Setup FFTW w/ MPI initialization 
           
        grid%Nr=grid%Ng
        grid%FFT_comm=mpi%comm_bands
        grid%nproc_FFT=mpi%nproc_band
        

        cut=1
        if(pruned) cut=2

        call fftw_mpi_init() !figure out what this does
           
        !Adjust the real grid to conform to FFT requirements, number of z-direction slices must be a product of the
        ! processors available for distribution (comm_band)
            size_z=grid%Nr(3)
            if(mod(int(size_z),grid%nproc_FFT).ne.0) then
                if(mpi%myid.eq.0) print *, 'Adjusting # of z direction grid points to divide &
                &    by available processors in real space decomp, Ecut unchanged)'
                size_z=(size_z/grid%nproc_FFT)*grid%nproc_FFT + grid%nproc_FFT
                grid%Nr(3)=int(size_z)
                if(mpi%myid.eq.0) print *, "Nr_z:,",  grid%Nr(3)
            endif
            grid%Nr_local(3)=grid%Nr(3)/grid%nproc_FFT
            grid%Ng(3)=grid%Nr(3)/cut
            grid%Ng_local(3)=grid%Ng(3)
            local_z=grid%Ng_local(3)

        !Adjust the grid to conform to FFT requirements, 1/cut number of x-direction slices must be a product of the
        ! processors available for distribution (comm_band)

            if(grid%gamma) then
                size_x=(grid%Nr(1)/2)/cut!+1)/cut
            else
                size_x=grid%Nr(1)/cut
            endif

            if(mod(int(size_x),grid%nproc_FFT).ne.0) then
                if(mpi%myid.eq.0) print *, 'Adjusting # of x direction grid points to divide &
                &    by available processors in reduced recip space decomp for only real, Ecut unchanged)'
                size_x=(size_x/grid%nproc_FFT)*grid%nproc_FFT+ grid%nproc_FFT
                grid%Ng(1)=int(size_x)
                if(mpi%myid.eq.0) print *, "Ng_x:,",  grid%Ng(1)
            else
                grid%Ng(1)=int(size_x)
            endif

            grid%Nr_local(1)=grid%Nr(1)
            grid%Ng_local(1)=grid%Ng(1)/grid%nproc_FFT
            local_x=grid%Ng_local(1)

            grid%Nr_local(2)=grid%Nr(2)
            grid%Ng(2)=grid%Nr(2)/cut
            grid%Ng_local(2)=grid%Ng(2)

            grid%myxyz(1)=mpi%myid_band
            grid%myxyz(2)=0
            grid%myxyz(3)=0

            grid%myxyz_r(1)=0
            grid%myxyz_r(2)=0
            grid%myxyz_r(3)=mpi%myid_band
        
            !FFTW setup distributed FFT's
            !Allocating and setting up work space sets plan for all FFT's on this grid

            size_z=grid%Nr_local(3);size_y=grid%Nr_local(2);size_x=grid%Nr_local(1);
            grid%alloc_local_R=size_z*size_y*size_x
            grid%work%p_R=fftw_alloc_complex(grid%alloc_local_R)
            call c_f_pointer(grid%work%p_R, grid%work%R, [size_x,size_y,size_z])

            if(grid%gamma) then
                grid%alloc_local_R_gamma=size_z*size_y*size_x
                grid%work%p_R_gamma=fftw_alloc_real(grid%alloc_local_R_gamma)
                call c_f_pointer(grid%work%p_R_gamma, grid%work%R_gamma, [size_x,size_y,size_z])
            endif

            size_z=grid%Ng_local(3);size_y=grid%Ng_local(2);size_x=grid%Ng_local(1);
            grid%alloc_local_G=size_z*size_y*size_x
            grid%work%p_G=fftw_alloc_complex(grid%alloc_local_G)
            call c_f_pointer(grid%work%p_G, grid%work%G, [size_z,size_y,size_x])
           
            !Working Space arrays
            grid%ncache=max(grid%Nr(1),grid%Nr(2),grid%Nr(3),1024) !use abinit ncache size
            grid%lotx=grid%ncache/grid%Nr(1)
            grid%loty=grid%ncache/grid%Nr(2)
            grid%lotz=grid%ncache/grid%Nr(3)

            size_tmp=grid%ncache
            grid%p_zw=fftw_alloc_complex(size_tmp)
            call c_f_pointer(grid%p_zw, grid%zw, [grid%ncache])

            !1D parallel
            size_tmp=product([grid%Nr(2),grid%Ng(1),grid%Nr_local(3)])
            grid%p_yXz2_out=fftw_alloc_complex(size_tmp)

            size_tmp=product([grid%Nr(3),grid%Ng(2)*grid%Ng_local(1)])
            grid%p_zYX_out=fftw_alloc_complex(size_tmp)

            size_tmp=product([grid%Ng(2)*grid%Ng(1), grid%Nr_local(3)])
            grid%p_YXz_out=fftw_alloc_complex(size_tmp)
            grid%p_buff=fftw_alloc_complex(size_tmp)

            call c_f_pointer(grid%p_yXz2_out, grid%outyXz2,[grid%Nr(2),grid%Ng(1),grid%Nr_local(3)])
            call c_f_pointer(grid%p_zYX_out,  grid%outzYX, [grid%Nr(3),grid%Ng(2)*grid%Ng_local(1)])
            call c_f_pointer(grid%p_YXz_out,  grid%outYXz, [grid%Ng(2)*grid%Ng(1), grid%Nr_local(3)])
            call c_f_pointer(grid%p_buff,     grid%buff,   [grid%Ng(2)*grid%Ng(1), grid%Nr_local(3)])

        !Create the FFTW3 plans, any functions that are defined on this grid can use these same plans
        !==================================================================

            rank=1
            !Normal FFT of fast variable (x of F(x,y))

            if(grid%gamma) then
                grid%plan_fft_xyz_x_lot = fftw_plan_many_dft_r2c(rank, [grid%Nr(1)], int(grid%lotx), &
                grid%work%R_gamma, [grid%Nr(1),grid%Nr(2),grid%Nr_local(3)], 1, grid%Nr(1),                          &
                grid%zw, [grid%ncache], 1, grid%Nr(1)/2+1, &
                FFTW_MEASURE)

                if(mod(grid%Nr_local(3)*grid%Nr(2),int(grid%lotx)).gt.0) &
                grid%plan_fft_xyz_x_left = fftw_plan_many_dft_r2c(rank, [grid%Nr(1)], &
                mod( grid%Nr_local(3)*grid%Nr(2),int(grid%lotx)), &
                grid%work%R_gamma, [grid%Nr(1),grid%Nr(2), grid%Nr_local(3)], 1, grid%Nr(1),                          &
                grid%zw, [grid%ncache], 1, grid%Nr(1)/2+1, &
                FFTW_MEASURE)
            else
                grid%plan_fft_xyz_x_lot = fftw_plan_many_dft(rank, [grid%Nr(1)], int(grid%lotx), &
                grid%work%R, [grid%Nr(1),grid%Nr(2),grid%Nr_local(3)], 1, grid%Nr(1),                          &
                grid%zw, [grid%ncache], 1, grid%Nr(1), &
                FFTW_FORWARD, FFTW_MEASURE)

                if(mod(grid%Nr_local(3)*grid%Nr(2),int(grid%lotx)).gt.0) &
                grid%plan_fft_xyz_x_left = fftw_plan_many_dft(rank, [grid%Nr(1)], &
                mod( grid%Nr_local(3)*grid%Nr(2),int(grid%lotx)), &
                grid%work%R, [grid%Nr(1),grid%Nr(2), grid%Nr_local(3)], 1, grid%Nr(1),                          &
                grid%zw, [grid%ncache], 1, grid%Nr(1), &
                FFTW_FORWARD, FFTW_MEASURE)
            endif
    
            !Normal FFT of fast variable (y of F(y,Gx))
            grid%plan_fft_yXz_y_lot = fftw_plan_many_dft(rank, [grid%Nr(2)], int(grid%loty), &
            grid%outyXz, [grid%Nr(2),grid%Ng(1),grid%Nr_local(3)], 1, grid%Nr(2),                          &
            grid%zw, [grid%ncache], 1, grid%Nr(2), &
            FFTW_FORWARD, FFTW_MEASURE)
    
            if(mod(grid%Ng(1),int(grid%loty)).gt.0) &
            grid%plan_fft_yXz_y_left = fftw_plan_many_dft(rank, [grid%Nr(2)], mod(grid%Ng(1),int(grid%loty)), &
            grid%outyXz, [grid%Nr(2),grid%Ng(1),grid%Nr_local(3)], 1, grid%Nr(2),                          &
            grid%zw, [grid%ncache], 1, grid%Nr(2), &
            FFTW_FORWARD, FFTW_MEASURE)
    
            !Normal FFT of fast variable (z of F(z,GyGx))
            grid%plan_fft_zYX_z_lot = fftw_plan_many_dft(rank, [grid%Nr(3)], int(grid%lotz), &
            grid%outzYX, [grid%Nr(3),grid%Ng(2)*grid%Ng_local(1)], 1, grid%Nr(3),                          &
            grid%zw, [grid%ncache], 1, grid%Nr(3), &
            FFTW_FORWARD, FFTW_MEASURE)
    
            if(mod(grid%Ng_local(1)*grid%Ng(2),int(grid%lotz)).gt.0) &
            grid%plan_fft_zYX_z_left = fftw_plan_many_dft(rank, [grid%Nr(3)], mod(grid%Ng_local(1)*grid%Ng(2),int(grid%lotz)), &
            grid%outzYX, [grid%Nr(3),grid%Ng(2)*grid%Ng_local(1)], 1, grid%Nr(3),                          &
            grid%zw, [grid%ncache], 1, grid%Nr(3), &
            FFTW_FORWARD, FFTW_MEASURE)     
            
            !==================================================================

            !Inverse FFT of fast variable (x of F(x,y))
            if(grid%gamma) then
                grid%plan_Ifft_xyz_x_lot = fftw_plan_many_dft_c2r(rank, [grid%Nr(1)], int(grid%lotx), &
                grid%zw, [grid%ncache], 1, grid%Nr(1)/2+1, &
                grid%work%R_gamma, [grid%Nr(1),grid%Nr(2),grid%Nr_local(3)], 1, grid%Nr(1),                          &
                FFTW_MEASURE)

                if(mod(grid%Nr_local(3)*grid%Nr(2),int(grid%lotx)).gt.0) &
                grid%plan_Ifft_xyz_x_left = fftw_plan_many_dft_c2r(rank, [grid%Nr(1)], &
                mod( grid%Nr_local(3)*grid%Nr(2),int(grid%lotx)), &
                grid%zw, [grid%ncache], 1, grid%Nr(1)/2+1, &
                grid%work%R_gamma, [grid%Nr(1),grid%Nr(2), grid%Nr_local(3)], 1, grid%Nr(1),                          &
                FFTW_MEASURE)
            else
                grid%plan_Ifft_xyz_x_lot = fftw_plan_many_dft(rank, [grid%Nr(1)], int(grid%lotx), &
                grid%zw, [grid%ncache], 1, grid%Nr(1), &
                grid%work%R, [grid%Nr(1),grid%Nr(2),grid%Nr_local(3)], 1, grid%Nr(1),                          &
                FFTW_BACKWARD, FFTW_MEASURE)

                if(mod(grid%Nr_local(3)*grid%Nr(2),int(grid%lotx)).gt.0) &
                grid%plan_Ifft_xyz_x_left = fftw_plan_many_dft(rank, [grid%Nr(1)],&
                mod( grid%Nr_local(3)*grid%Nr(2),int(grid%lotx)), &
                grid%zw, [grid%ncache], 1, grid%Nr(1), &
                grid%work%R, [grid%Nr(1),grid%Nr(2), grid%Nr_local(3)], 1, grid%Nr(1),                          &
                FFTW_BACKWARD, FFTW_MEASURE)
            endif
    
            !Inverse FFT of fast variable (y of F(y,Gx))
            grid%plan_Ifft_yXz_y_lot = fftw_plan_many_dft(rank, [grid%Nr(2)], int(grid%loty), &
            grid%zw, [grid%ncache], 1, grid%Nr(2), &
            grid%outyXz, [grid%Nr(2),grid%Ng(1),grid%Nr_local(3)], 1, grid%Nr(2),                          &
            FFTW_BACKWARD, FFTW_MEASURE)
    
            if(mod(grid%Ng(1),int(grid%loty)).gt.0) &
            grid%plan_Ifft_yXz_y_left = fftw_plan_many_dft(rank, [grid%Nr(2)], mod(grid%Ng(1),int(grid%loty)), &
            grid%zw, [grid%ncache], 1, grid%Nr(2), &
            grid%outyXz, [grid%Nr(2),grid%Ng(1),grid%Nr_local(3)], 1, grid%Nr(2),                          &
            FFTW_BACKWARD, FFTW_MEASURE)
    
            !Inverse FFT of fast variable (z of F(z,GyGx))
            grid%plan_Ifft_zYX_z_lot = fftw_plan_many_dft(rank, [grid%Nr(3)], int(grid%lotz), &
            grid%zw, [grid%ncache], 1, grid%Nr(3), &
            grid%outzYX, [grid%Nr(3),grid%Ng(2)*grid%Ng_local(1)], 1, grid%Nr(3),                          &
            FFTW_BACKWARD, FFTW_MEASURE)
    
            if(mod(grid%Ng_local(1)*grid%Ng(2),int(grid%lotz)).gt.0) &
            grid%plan_Ifft_zYX_z_left = fftw_plan_many_dft(rank, [grid%Nr(3)], mod(grid%Ng_local(1)*grid%Ng(2),int(grid%lotz)), &
            grid%zw, [grid%ncache], 1, grid%Nr(3), &
            grid%outzYX, [grid%Nr(3),grid%Ng(2)*grid%Ng_local(1)], 1, grid%Nr(3),                          &
            FFTW_BACKWARD, FFTW_MEASURE)         

        !Start defining the physical grid (real & reciprocal space values, both locally (processor dependent) &
        ! and Globaly (full grid)
        allocate(grid%R(grid%Nr_local(1),grid%Nr_local(2),grid%Nr_local(3),3))
        allocate(grid%G(grid%Ng_local(3),grid%Ng_local(2),grid%Ng_local(1),3))
        allocate(grid%G2(grid%Ng_local(3),grid%Ng_local(2),grid%Ng_local(1)))
        if(grid%ecutsm.gt.0.0_dp) then
            allocate(grid%p(grid%Ng_local(3),grid%Ng_local(2),grid%Ng_local(1)))
            allocate(grid%pprime(grid%Ng_local(3),grid%Ng_local(2),grid%Ng_local(1)))
            allocate(grid%pdoubleprime(grid%Ng_local(3),grid%Ng_local(2),grid%Ng_local(1)))
            allocate(X(grid%Ng_local(3),grid%Ng_local(2),grid%Ng_local(1)))
        endif

        grid%dG(:)=2.0_dp*pi/grid%Box_Length(:)
        grid%dR(:)=grid%Box_Length(:)/grid%Nr(:)
        grid%k(:)=grid%k(:)*grid%dG(:)
        
        if(mpi%myid.eq.(mpi%nproc-1))print *, 'Real Spacing (bohr): ', grid%dR(:) 
        do k = 1, grid%Nr_local(3); do j = 1, grid%Nr_local(2); do i = 1, grid%Nr_local(1)
            ijk_global = [i, j, k] + grid%myxyz_r(:)*grid%Nr_local(:)
            grid%R(i,j,k,:) = (ijk_global-1) * grid%dR(:)
          
        end do; end do; end do

        if(mpi%myid.eq.(mpi%nproc-1))print *, 'Reciprocal Spacing (1/bohr): ', grid%dG(:) 

        if(grid%gamma) then
            do i = 1, grid%Ng_local(1); do j = 1, grid%Ng_local(2); do k = 1, grid%Ng_local(3);
                ijk_global = [i, j, k] + grid%myxyz(:)*grid%Ng_local(:)
                grid%G(k, j, i, 1) = grid%dG(1) * (ijk_global(1) - 1)
                grid%G(k, j, i, 2:3) = grid%dG(2:3) * &
                                        (ijk_global(2:3) - 1 - grid%Ng(2:3)*nint((ijk_global(2:3)-1.5_dp)/grid%Ng(2:3)))
                grid%G2(k, j, i)  = sum(grid%G(k,j,i,:)**2)
            end do; end do; end do
        else
            do i = 1, grid%Ng_local(1); do j = 1, grid%Ng_local(2); do k = 1, grid%Ng_local(3);
                ijk_global = [i, j, k] + grid%myxyz(:)*grid%Ng_local(:)
                grid%G(k, j, i,:) = grid%dG(:) * (ijk_global - 1 - grid%Ng(:)*nint((ijk_global-1.5_dp)/grid%Ng)) + grid%k(:)
                grid%G2(k, j, i)  = sum(grid%G(k,j,i,:)**2)
            end do; end do; end do
        endif

        if(grid%ecutsm.gt.0.0_dp) then
            X=(grid%Ecut-0.5_dp*grid%G2)/grid%ecutsm
            where(X.gt.1.0_dp)
                grid%p=1.0_dp
                grid%pprime=0.0_dp
                grid%pdoubleprime=0.0_dp
            elsewhere(X.lt.tiny(1.0_dp)) !p should be infinite at X, but 
                grid%p=0.0_dp
                grid%pprime=0.0_dp
                grid%pdoubleprime=0.0_dp
            elsewhere
                grid%p=1.0_dp/(X**2*(3.0_dp+X*(1.0_dp+X*(-6.0_dp+3.0_dp*X))))
                grid%pprime=-3.0*(X-1.0_dp)**2*X*(2.0_dp+5.0_dp*X)*grid%p**2
                grid%pdoubleprime=6.0_dp*X**2* &
                    (9.0_dp+X*(8.0_dp+X*(-52.0_dp+X*(-3.0_dp+X*(137.0_dp+X*(-144.0+45.0_dp*X))))))*grid%p**3
            endwhere
            deallocate(X)
        endif

        !Set up the Energy cutoff of the plane wave orbital & density/potential coefficients
        allocate(grid%cutwf(grid%Ng_local(3),grid%Ng_local(2),grid%Ng_local(1)))
        allocate(grid%cutden(grid%Ng_local(3),grid%Ng_local(2),grid%Ng_local(1)))

        where(grid%G2.le.(grid%Ecut*2.0_dp))
            grid%cutwf = 1.0_dp
        elsewhere
            grid%cutwf = 0.0_dp
        end where
        where(grid%G2.le.(grid%Ecut*8.0_dp))
            grid%cutden = 1.0_dp
        elsewhere
            grid%cutden = 0.0_dp
        end where

        !Vector potential 
        grid%A =0.0_dp

        !set up the reduced grid removing 0's for the wave-function (those compoenent already set to zero by cutwf)
        grid%Ng_small(1)=grid%Ng_local(1)
        grid%Ng_small(2)=grid%Ng_local(2)/(3-cut) + mod(grid%Ng_local(2),3-cut)
        grid%Ng_small(3)=grid%Ng_local(3)/(3-cut) + mod(grid%Ng_local(3),3-cut)

        grid%reduction_option=3-cut
    end subroutine

    subroutine map_between_distributed_grids(grid_from, grid_to, mpi, inter_grid)
        use, intrinsic :: iso_c_binding
        use system_type, only : system_struct
        use grids_type, only : grid_struct, inter_grid_struct
        use parallel_type, only : parallel_struct
        use parallel_mod, only: parallel_task, parallel_wait

        type(parallel_struct), intent(in) :: mpi
        type(grid_struct), intent(in) :: grid_from, grid_to
        type(inter_grid_struct), intent(inout) :: inter_grid

        integer, allocatable :: who_has_local(:,:,:), shift_by_proc(:,:), map(:,:,:,:), last_index(:), packed_map(:,:), &
        recv_countv(:), send_countv(:) 
        integer:: i,j,k,ijk_global_from(3),ijk_global_to(3), dir
        !map=ix_from_local, iy_from_local, iz_from_local, 1-4 (1-3=ix-iz_to_local, 4= proc(band))

        !Check that the grids have the same grid spacing in reciprocal space (same real space box size) 
        if(any(abs(grid_from%dG(:)-grid_to%dG(:)).gt.tiny(1.0_dp))) then
            print *, 'Trying to map grids from 2 different size boxes in real pace'
            stop
        endif

        if(grid_from%FFT_type.gt.1) then
            allocate(map(grid_from%Ng_local(3),grid_from%Ng_local(2),grid_from%Ng_local(1),4))
        else 
            allocate(map(grid_from%Ng_local(1),grid_from%Ng_local(2),grid_from%Ng_local(3),4))
        endif

        allocate(inter_grid%count_(mpi%nproc_band,mpi%nproc_band))
        allocate(inter_grid%sdisp(mpi%nproc_band))
        allocate(inter_grid%rdisp(mpi%nproc_band))
        allocate(send_countv(mpi%nproc_band))
        allocate(recv_countv(mpi%nproc_band))

        inter_grid%count_=0
        inter_grid%sdisp=0
        inter_grid%rdisp=0

        !First figure out who locally owns different space on the recieving grid
        if(grid_to%FFT_type.gt.1) then 
           ! print *, grid_to%Ng(:), grid_to%Ng_local(:), grid_to%myxyz(:)*grid_to%Ng_local(:)
            allocate(who_has_local(grid_to%Ng(3),grid_to%Ng(2),grid_to%Ng(1)))
        else
            allocate(who_has_local(grid_to%Ng(1),grid_to%Ng(2),grid_to%Ng(3)))
        endif
        who_has_local=0
        do k = 1, grid_to%Ng_local(3); do j = 1, grid_to%Ng_local(2); do i = 1, grid_to%Ng_local(1)
            ijk_global_to = [i,j,k] + grid_to%myxyz(:)*grid_to%Ng_local(:)
            if(grid_to%FFT_type.gt.1) then
                who_has_local(ijk_global_to(3),ijk_global_to(2),ijk_global_to(1))=mpi%myid_band
           else
                who_has_local(ijk_global_to(1),ijk_global_to(2),ijk_global_to(3))=mpi%myid_band
           endif
        enddo; enddo; enddo
        call parallel_task('sum', who_has_local, mpi, 'band')

        !Second collect the grid shifts fro each processor
        allocate(shift_by_proc(mpi%nproc_band,3))
        shift_by_proc=0
        shift_by_proc(mpi%myid_band+1,1)=grid_to%myxyz(1)
        shift_by_proc(mpi%myid_band+1,2)=grid_to%myxyz(2)
        shift_by_proc(mpi%myid_band+1,3)=grid_to%myxyz(3)
        call parallel_task('sum', shift_by_proc, mpi, 'band')
        

        do k = 1, grid_from%Ng_local(3); do j = 1, grid_from%Ng_local(2); do i = 1, grid_from%Ng_local(1)
            !real spacing index
            ijk_global_from = [i,j,k] + grid_from%myxyz(:)*grid_from%Ng_local(:)

            !from regular order real to flipped recip space index 
            if(grid_from%gamma) then
                ijk_global_from(2:3) = ijk_global_from(2:3) &
                - 1 - grid_from%Ng(2:3)*nint((ijk_global_from(2:3)-1.5_dp)/grid_from%Ng(2:3))

                ijk_global_from(1) = ijk_global_from(1) -1
            else
                ijk_global_from = ijk_global_from - 1 - grid_from%Ng(:)*nint((ijk_global_from-1.5_dp)/grid_from%Ng)
            endif

            ijk_global_to(:)=-1
            do dir=1,3
                !only points on _from grid that fit onto _to grid
                if(grid_to%gamma.and.dir.eq.1) then
                    if((ijk_global_from(dir).ge.0) .and. (ijk_global_from(dir).lt.grid_to%Ng(dir))) then
                        ijk_global_to(dir)=ijk_global_from(dir) + 1
                    endif
                else                    
                    if(mod(grid_to%Ng(dir),2).eq.0) then
                        if(abs(ijk_global_from(dir)).lt.((grid_to%Ng(dir)/2))) then !if it fits
                            if((ijk_global_from(dir).lt.0)) then !if it is negative
                                ijk_global_to(dir)=ijk_global_from(dir) + grid_to%Ng(dir) + 1
                            else 
                                ijk_global_to(dir)=ijk_global_from(dir) + 1
                            endif
                        endif
                    else
                        if(abs(ijk_global_from(dir)).le.((grid_to%Ng(dir)/2))) then !if it fits
                            if((ijk_global_from(dir).lt.0)) then !if it is negative
                                ijk_global_to(dir)=ijk_global_from(dir) + grid_to%Ng(dir) + 1
                            else 
                                ijk_global_to(dir)=ijk_global_from(dir) + 1
                            endif
                        endif
                    endif
                endif
            enddo

            if(grid_from%FFT_type.gt.1) then 
                if(grid_to%FFT_type.gt.1) then
                    if(all(ijk_global_to.ge.0)) then
                        map(k,j,i,4)=who_has_local(ijk_global_to(3),ijk_global_to(2),ijk_global_to(1))
                        inter_grid%count_(mpi%myid_band+1, map(k,j,i,4)+1)=inter_grid%count_(mpi%myid_band+1,map(k,j,i,4)+1)+1
                        map(k,j,i,1)=ijk_global_to(1) - shift_by_proc(map(k,j,i,4)+1,1)*grid_to%Ng_local(1)  
                        map(k,j,i,2)=ijk_global_to(2) - shift_by_proc(map(k,j,i,4)+1,2)*grid_to%Ng_local(2)
                        map(k,j,i,3)=ijk_global_to(3) - shift_by_proc(map(k,j,i,4)+1,3)*grid_to%Ng_local(3)
                    else
                        map(k,j,i,4)=-1 !The _to grid is smaller than the _from grid, 
                                        !so there is nowhere to put this point, it will be discarded on the _to grid
                    endif
                else
                    if(all(ijk_global_to.ge.0)) then
                        map(k,j,i,4)=who_has_local(ijk_global_to(1),ijk_global_to(2),ijk_global_to(3))
                        inter_grid%count_(mpi%myid_band+1, map(k,j,i,4)+1)=inter_grid%count_(mpi%myid_band+1,map(k,j,i,4)+1)+1
                        map(k,j,i,1)=ijk_global_to(1) - shift_by_proc(map(k,j,i,4)+1,1)*grid_to%Ng_local(1)  
                        map(k,j,i,2)=ijk_global_to(2) - shift_by_proc(map(k,j,i,4)+1,2)*grid_to%Ng_local(2)
                        map(k,j,i,3)=ijk_global_to(3) - shift_by_proc(map(k,j,i,4)+1,3)*grid_to%Ng_local(3)
                    else
                        map(k,j,i,4)=-1 !The _to grid is smaller than the _from grid, 
                                        !so there is nowhere to put this point, it will be discarded on the _to grid
                    endif
                endif
            else
                if(grid_to%FFT_type.gt.1) then
                    if(all(ijk_global_to.ge.0)) then
                        map(i,j,k,4)=who_has_local(ijk_global_to(3),ijk_global_to(2),ijk_global_to(1))
                        inter_grid%count_(mpi%myid_band+1, map(i,j,k,4)+1)=inter_grid%count_(mpi%myid_band+1,map(i,j,k,4)+1)+1
                        map(i,j,k,1)=ijk_global_to(1) - shift_by_proc(map(i,j,k,4)+1,1)*grid_to%Ng_local(1)  
                        map(i,j,k,2)=ijk_global_to(2) - shift_by_proc(map(i,j,k,4)+1,2)*grid_to%Ng_local(2)
                        map(i,j,k,3)=ijk_global_to(3) - shift_by_proc(map(i,j,k,4)+1,3)*grid_to%Ng_local(3)
                    else
                        map(i,j,k,4)=-1 !The _to grid is smaller than the _from grid, 
                                        !so there is nowhere to put this point, it will be discarded on the _to grid
                    endif
                else
                    if(all(ijk_global_to.ge.0)) then
                        map(i,j,k,4)=who_has_local(ijk_global_to(1),ijk_global_to(2),ijk_global_to(3))
                        inter_grid%count_(mpi%myid_band+1, map(i,j,k,4)+1)=inter_grid%count_(mpi%myid_band+1,map(i,j,k,4)+1)+1
                        map(i,j,k,1)=ijk_global_to(1) - shift_by_proc(map(i,j,k,4)+1,1)*grid_to%Ng_local(1)  
                        map(i,j,k,2)=ijk_global_to(2) - shift_by_proc(map(i,j,k,4)+1,2)*grid_to%Ng_local(2)
                        map(i,j,k,3)=ijk_global_to(3) - shift_by_proc(map(i,j,k,4)+1,3)*grid_to%Ng_local(3)
                    else
                        map(i,j,k,4)=-1 !The _to grid is smaller than the _from grid, 
                                        !so there is nowhere to put this point, it will be discarded on the _to grid
                    endif
                endif
            endif
        end do; end do; end do
        call parallel_task('sum', inter_grid%count_, mpi, 'band')
        do i = 1, mpi%nproc_band
            inter_grid%sdisp(i)=sum(inter_grid%count_(mpi%myid_band+1,:(i-1)))
            inter_grid%rdisp(i)=sum(inter_grid%count_(:(i-1),mpi%myid_band+1))
        enddo

        !inter_grid%pack_index, where to put local data in the send array for alltoallv distribution
        if(grid_from%FFT_type.gt.1) then 
            allocate(inter_grid%pack_index(grid_from%Ng_local(3),grid_from%Ng_local(2),grid_from%Ng_local(1)))
        else
            allocate(inter_grid%pack_index(grid_from%Ng_local(1),grid_from%Ng_local(2),grid_from%Ng_local(3)))
        endif
        allocate(last_index(mpi%nproc_band))
        allocate(packed_map(sum(inter_grid%count_(mpi%myid_band+1,:)),3))
        allocate(inter_grid%unpack_index(sum(inter_grid%count_(:,mpi%myid_band+1)),3))

        inter_grid%pack_index=-1
        last_index(:)=inter_grid%sdisp(:)
        do k = 1, grid_from%Ng_local(3); do j = 1, grid_from%Ng_local(2); do i = 1, grid_from%Ng_local(1)
            if(grid_from%FFT_type.gt.1) then 
                if(map(k,j,i,4).ge.0) then
                    last_index(map(k,j,i,4)+1)=last_index(map(k,j,i,4)+1)+1
                    inter_grid%pack_index(k,j,i)=last_index(map(k,j,i,4)+1)
                endif
            else
                if(map(i,j,k,4).ge.0) then
                    last_index(map(i,j,k,4)+1)=last_index(map(i,j,k,4)+1)+1
                    inter_grid%pack_index(i,j,k)=last_index(map(i,j,k,4)+1)
                endif
            endif
        end do; end do; end do
          !pack the map so it can be sent to the reciever for unpacking
        do k = 1, grid_from%Ng_local(3); do j = 1, grid_from%Ng_local(2); do i = 1, grid_from%Ng_local(1)
            if(grid_from%FFT_type.gt.1) then 
                if(inter_grid%pack_index(k,j,i).gt.0) then
                    packed_map(inter_grid%pack_index(k,j,i),1)=map(k,j,i,1)
                    packed_map(inter_grid%pack_index(k,j,i),2)=map(k,j,i,2)
                    packed_map(inter_grid%pack_index(k,j,i),3)=map(k,j,i,3)
                endif
            else
                if(inter_grid%pack_index(i,j,k).gt.0) then
                    packed_map(inter_grid%pack_index(i,j,k),1)=map(i,j,k,1)
                    packed_map(inter_grid%pack_index(i,j,k),2)=map(i,j,k,2)
                    packed_map(inter_grid%pack_index(i,j,k),3)=map(i,j,k,3)
                endif
            endif
        enddo; end do; end do
        deallocate(shift_by_proc)
        deallocate(who_has_local)
        deallocate(map)

        if(mpi%nproc_band.gt.1) then
            send_countv(:)=inter_grid%count_(mpi%myid_band+1,:)
            recv_countv(:)=inter_grid%count_(:,mpi%myid_band+1)
            do dir=1,3
                call parallel_task('alltoall', packed_map(:,dir), mpi, 'band', inter_grid%unpack_index(:,dir), &
                send_countv=send_countv, recv_countv=recv_countv, &
                sdispls=inter_grid%sdisp(:), rdispls=inter_grid%rdisp(:))
            enddo
        else
            inter_grid%unpack_index(:,:)=packed_map(:,:)
        endif
         
        deallocate(packed_map)
        
    end subroutine

    subroutine Gspace_grid_to_grid_transfer(grid_from, grid_to, mpi, data_in, data_out, inter_grid)
        use, intrinsic :: iso_c_binding
        use system_type, only : system_struct
        use grids_type, only : grid_struct, inter_grid_struct
        use parallel_type, only : parallel_struct
        use parallel_mod, only: parallel_task

        type(parallel_struct), intent(in) :: mpi
        type(grid_struct), intent(in) :: grid_from, grid_to
        type(inter_grid_struct), intent(in) :: inter_grid
        complex(dp), intent(in) :: data_in(:,:,:)
        complex(dp), intent(inout) :: data_out(:,:,:)
        complex(dp), allocatable :: packed_data_send(:), packed_data_recieve(:)
        integer, allocatable ::recv_countv(:), send_countv(:) 

        integer:: i,j,k
        allocate(send_countv(mpi%nproc_band))
        allocate(recv_countv(mpi%nproc_band))
        allocate(packed_data_send(sum(inter_grid%count_(mpi%myid_band+1,:))))
        allocate(packed_data_recieve(size(inter_grid%unpack_index(:,1))))

        if(grid_from%FFT_type.gt.1) then 
            do i = 1, grid_from%Ng_local(1); do j = 1, grid_from%Ng_local(2); do k = 1, grid_from%Ng_local(3)
                if(inter_grid%pack_index(k,j,i).gt.0) then
                    packed_data_send(inter_grid%pack_index(k,j,i))=data_in(k,j,i)
                endif
            enddo; end do; end do
        else
            do k = 1, grid_from%Ng_local(3); do j = 1, grid_from%Ng_local(2); do i = 1, grid_from%Ng_local(1)
                if(inter_grid%pack_index(i,j,k).gt.0) then
                    packed_data_send(inter_grid%pack_index(i,j,k))=data_in(i,j,k)
                endif
            enddo; end do; end do
        endif

        if(mpi%nproc_band.gt.1) then
            send_countv(:)=inter_grid%count_(mpi%myid_band+1,:)
            recv_countv(:)=inter_grid%count_(:,mpi%myid_band+1)
            call parallel_task('alltoall', packed_data_send(:), mpi, 'band', packed_data_recieve(:), &
            send_countv=send_countv, recv_countv=recv_countv, &
            sdispls=inter_grid%sdisp(:), rdispls=inter_grid%rdisp(:))
        else
            packed_data_recieve=packed_data_send
        endif
        data_out=0.0_dp

        do i = 1, size(packed_data_recieve)
            if(all(inter_grid%unpack_index(i,:).gt.0))  then
                if(grid_to%FFT_type.gt.1) then 
                    data_out(inter_grid%unpack_index(i,3), &
                        inter_grid%unpack_index(i,2), &
                        inter_grid%unpack_index(i,1)) = packed_data_recieve(i)
                else
                    data_out(inter_grid%unpack_index(i,1), &
                        inter_grid%unpack_index(i,2), &
                        inter_grid%unpack_index(i,3)) = packed_data_recieve(i)
                endif
                    
            else
                print *, 'Somehow Data is sent off grid ', i, inter_grid%unpack_index(i,:)
                flush(6)
                stop
            endif
        enddo
        data_out=data_out*grid_to%cutden
        
    end subroutine
    
   

    subroutine allocate_local_field0_G(psi, grid)
        use grids_type, only : grid_struct

        type(grid_struct), intent(in) :: grid
        complex(dp), allocatable, intent(inout) :: psi(:,:,:)
        if(grid%FFT_type.gt.1) then
            allocate(psi(grid%Ng_local(3),grid%Ng_local(2),grid%Ng_local(1)))
        else
            allocate(psi(grid%Ng_local(1),grid%Ng_local(2),grid%Ng_local(3)))
        endif
    end subroutine

    subroutine allocate_local_field0_G_real(psi, grid)
        use grids_type, only : grid_struct

        type(grid_struct), intent(in) :: grid
        real(dp), allocatable, intent(inout) :: psi(:,:,:)
        if(grid%FFT_type.gt.1) then
            allocate(psi(grid%Ng_local(3),grid%Ng_local(2),grid%Ng_local(1)))
        else
            allocate(psi(grid%Ng_local(1),grid%Ng_local(2),grid%Ng_local(3)))       
        endif

    end subroutine

    subroutine allocate_local_field1_G(psi, grid, N)
        use grids_type, only : grid_struct

        type(grid_struct), intent(in) :: grid
        complex(dp), allocatable, intent(inout) :: psi(:,:,:,:)
        integer , intent(in) :: N

        if(grid%FFT_type.gt.1) then
            allocate(psi(grid%Ng_local(3),grid%Ng_local(2),grid%Ng_local(1),N))
        else
            allocate(psi(grid%Ng_local(1),grid%Ng_local(2),grid%Ng_local(3),N))  
        endif


    end subroutine

    subroutine allocate_local_field1_G_real(psi, grid, N)
        use grids_type, only : grid_struct

        type(grid_struct), intent(in) :: grid
        real(dp), allocatable, intent(inout) :: psi(:,:,:,:)
        integer , intent(in) :: N

        if(grid%FFT_type.gt.1) then
            allocate(psi(grid%Ng_local(3),grid%Ng_local(2),grid%Ng_local(1),N))
        else
            allocate(psi(grid%Ng_local(1),grid%Ng_local(2),grid%Ng_local(3),N))  
        endif


    end subroutine

    subroutine allocate_local_field2_G(psi, grid, N, M)
        use grids_type, only : grid_struct

        type(grid_struct), intent(in) :: grid
        complex(dp), allocatable, intent(inout) :: psi(:,:,:,:,:)
        integer , intent(in) :: N, M

        if(grid%FFT_type.gt.1) then
            allocate(psi(grid%Ng_local(3),grid%Ng_local(2),grid%Ng_local(1),N, M))
        else
            allocate(psi(grid%Ng_local(1),grid%Ng_local(2),grid%Ng_local(3),N, M))
        endif


    end subroutine

    subroutine allocate_local_field0_R_real(psi, grid)
        use grids_type, only : grid_struct

        type(grid_struct), intent(in) :: grid
        real(dp), allocatable, intent(inout) :: psi(:,:,:)

        allocate(psi(grid%Nr_local(1),grid%Nr_local(2),grid%Nr_local(3)))
    end subroutine

    subroutine allocate_local_field1_R_real(psi, grid,N)
        use grids_type, only : grid_struct

        type(grid_struct), intent(in) :: grid
        real(dp), allocatable, intent(inout) :: psi(:,:,:,:)
        integer , intent(in) :: N

        allocate(psi(grid%Nr_local(1),grid%Nr_local(2),grid%Nr_local(3),N))
    end subroutine

    subroutine allocate_local_field2_R_real(psi, grid,N, M)
        use grids_type, only : grid_struct

        type(grid_struct), intent(in) :: grid
        real(dp), allocatable, intent(inout) :: psi(:,:,:,:,:)
        integer , intent(in) :: N,M

        allocate(psi(grid%Nr_local(1),grid%Nr_local(2),grid%Nr_local(3),N,M))
    end subroutine

    subroutine allocate_local_field0_R_complex(psi, grid)
        use grids_type, only : grid_struct

        type(grid_struct), intent(in) :: grid
        complex(dp), allocatable, intent(inout) :: psi(:,:,:)

        allocate(psi(grid%Nr_local(1),grid%Nr_local(2),grid%Nr_local(3)))
    end subroutine

    subroutine allocate_local_field1_R_complex(psi, grid,N)
        use grids_type, only : grid_struct

        type(grid_struct), intent(in) :: grid
        complex(dp), allocatable, intent(inout) :: psi(:,:,:,:)
        integer , intent(in) :: N

        allocate(psi(grid%Nr_local(1),grid%Nr_local(2),grid%Nr_local(3),N))
    end subroutine

    subroutine allocate_local_field2_R_complex(psi, grid,N, M)
        use grids_type, only : grid_struct

        type(grid_struct), intent(in) :: grid
        complex(dp), allocatable, intent(inout) :: psi(:,:,:,:,:)
        integer , intent(in) :: N,M

        allocate(psi(grid%Nr_local(1),grid%Nr_local(2),grid%Nr_local(3),N,M))
    end subroutine

    subroutine allocate_full_field0_G(psi, grid)
        use grids_type, only : grid_struct

        type(grid_struct), intent(in) :: grid
        complex(dp), allocatable, intent(inout) :: psi(:,:,:)
        
        if(grid%FFT_type.gt.1) then
            allocate(psi(grid%Ng(3),grid%Ng(2),grid%Ng(1)))
        else
            allocate(psi(grid%Ng(1),grid%Ng(2),grid%Ng(3)))          
        endif

    end subroutine

    subroutine allocate_full_field0_R(psi, grid)
        use grids_type, only : grid_struct

        type(grid_struct), intent(in) :: grid
        real(dp), allocatable, intent(inout) :: psi(:,:,:)

        allocate(psi(grid%Nr(1),grid%Nr(2),grid%Nr(3)))
    end subroutine

end module
