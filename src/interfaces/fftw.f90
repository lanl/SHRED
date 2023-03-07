module fft
use types, only : dp

implicit none

private 

public :: real_to_recip, recip_to_real, FFT3D_pruned_1D_decomp

interface real_to_recip
    module procedure real_to_recip_fftw_in_place
    module procedure real_to_recip_fftw_save
end interface

interface recip_to_real
    module procedure recip_to_real_fftw_in_place
    module procedure recip_to_real_fftw_save
end interface

contains

    subroutine real_to_recip_fftw_in_place(in_fftw, grids)
        use, intrinsic :: iso_c_binding
        use odp_type, only: field_struct
        use grids_type, only : grid_struct
        include 'fftw3-mpi.f03'
        type(field_struct), intent(inout) :: in_fftw
        type(grid_struct), intent(in), target:: grids(:)
        type(grid_struct), pointer:: grid
        grid=>grids(in_fftw%grid)
  
        if(grid%gamma) then
            in_fftw%R_gamma=real(in_fftw%R)
            if (grid%FFT_type.eq.11) then
                call FFT3D_1D_decomp_gamma(grid, in_fftw%R_gamma, in_fftw%G) 
            !Small Recip space FFT, transpose
            else if (grid%FFT_type.eq.21) then
                call FFT3D_pruned_1D_decomp_gamma(grid, in_fftw%R_gamma, in_fftw%G) 
            else 
                print *, 'Error FFT grid not initialized real_to_recip'
                stop
            endif
            where(abs(grid%G(:,:,:,1)).gt.tiny(1.0_dp)) 
                in_fftw%G=in_fftw%G*sqrt(2.0_dp) 
            end where
        else
            !FULL GRID FFT, no transpose
            if(grid%FFT_type.eq.1) then 
                call fftw_mpi_execute_dft(grid%plan_forward, in_fftw%R, in_fftw%G)
            else if (grid%FFT_type.eq.11) then
                call FFT3D_1D_decomp(grid, in_fftw%R, in_fftw%G) 
            !Small Recip space FFT, transpose
            else if (grid%FFT_type.eq.21) then
                call FFT3D_pruned_1D_decomp(grid, in_fftw%R, in_fftw%G) 
            else 
                print *, 'Error FFT grid not initialized real_to_recip'
                stop
            endif
        endif
        in_fftw%G=in_fftw%G/product(grid%Nr(:))

    end subroutine
    
    subroutine recip_to_real_fftw_in_place(in_fftw, grids)
        use, intrinsic :: iso_c_binding
        use odp_type, only:         
        use odp_type, only: field_struct
        use grids_type, only : grid_struct
        include 'fftw3-mpi.f03'
        type(field_struct), intent(inout) :: in_fftw
        type(grid_struct), intent(in), target:: grids(:)
        type(grid_struct), pointer:: grid
        grid=>grids(in_fftw%grid)        
        if(grid%gamma) then

            where(abs(grid%G(:,:,:,1)).gt.tiny(1.0_dp))  
                in_fftw%G=in_fftw%G/sqrt(2.0_dp) 
            end where

            if (grid%FFT_type.eq.11) then
                call IFFT3D_1D_decomp_gamma(grid, in_fftw%G, in_fftw%R_gamma)
            else if (grid%FFT_type.eq.21) then
                call IFFT3D_pruned_1D_decomp_gamma(grid, in_fftw%G, in_fftw%R_gamma)
            else 
                print *, 'Error FFT grid not initialized recip_to_real'
                stop
            endif
            in_fftw%R=in_fftw%R_gamma

            where(abs(grid%G(:,:,:,1)).gt.tiny(1.0_dp)) 
                in_fftw%G=in_fftw%G*sqrt(2.0_dp) 
            end where
        else
            if(grid%FFT_type.eq.1) then 
                call fftw_mpi_execute_dft(grid%plan_reverse, in_fftw%G, in_fftw%R)
            else if (grid%FFT_type.eq.11) then
                call IFFT3D_1D_decomp(grid, in_fftw%G, in_fftw%R)
            else if (grid%FFT_type.eq.21) then
                call IFFT3D_pruned_1D_decomp(grid, in_fftw%G, in_fftw%R)
            else 
                print *, 'Error FFT grid not initialized recip_to_real'
                stop
            endif
        endif
    end subroutine

    subroutine real_to_recip_fftw_save(in_fftw, result_fftw, grids)
        use, intrinsic :: iso_c_binding
        use odp_type, only: field_struct
        use grids_type, only : grid_struct
        include 'fftw3-mpi.f03'
        type(field_struct), intent(in) :: in_fftw
        type(field_struct), intent(inout) :: result_fftw
        type(grid_struct), intent(in), target:: grids(:)
        type(grid_struct), pointer:: grid
        grid=>grids(in_fftw%grid)
        if(in_fftw%grid.ne.result_fftw%grid) then
            print *, 'Attempting to transfer between different grids in FFT, stopping'
            stop
        endif            
        result_fftw%R=in_fftw%R

        if(grid%gamma) then
            in_fftw%R_gamma=real(in_fftw%R)
            if (grid%FFT_type.eq.11) then
                call FFT3D_1D_decomp_gamma(grid, in_fftw%R_gamma, result_fftw%G)
                !Small Recip space FFT, transpose
            else if (grid%FFT_type.eq.21) then
                call FFT3D_pruned_1D_decomp_gamma(grid, in_fftw%R_gamma, result_fftw%G)
            else 
                print *, 'Error FFT grid not initialized real_to_recip'
                stop
            endif
            where(abs(grid%G(:,:,:,1)).gt.tiny(1.0_dp))  
                result_fftw%G=result_fftw%G*sqrt(2.0_dp) 
            end where
        else
            if(grid%FFT_type.eq.1) then 
                call fftw_mpi_execute_dft(grid%plan_forward, in_fftw%R, result_fftw%G)
            else if (grid%FFT_type.eq.11) then
                call FFT3D_1D_decomp(grid, in_fftw%R, result_fftw%G)
            else if (grid%FFT_type.eq.21) then
                call FFT3D_pruned_1D_decomp(grid, in_fftw%R, result_fftw%G)
            else 
                print *, 'Error FFT grid not initialized real_to_recip'
                stop
            endif
        endif

        result_fftw%G=result_fftw%G/product(grid%Nr(:))
    end subroutine
    
    subroutine recip_to_real_fftw_save(in_fftw, result_fftw, grids)
        use, intrinsic :: iso_c_binding
        use odp_type, only: field_struct
        use grids_type, only : grid_struct

        include 'fftw3-mpi.f03'
        type(field_struct), intent(in) :: in_fftw
        type(field_struct), intent(inout) :: result_fftw
        type(grid_struct), intent(in), target:: grids(:)
        type(grid_struct), pointer:: grid
        grid=>grids(in_fftw%grid)
        if(in_fftw%grid.ne.result_fftw%grid) then
            print *, 'Attempting to transfer between different grids in FFT, stopping'
            stop
        endif   

        result_fftw%G=in_fftw%G
        if(grid%gamma) then

            where(abs(grid%G(:,:,:,1)).gt.tiny(1.0_dp)) 
                result_fftw%G=in_fftw%G/sqrt(2.0_dp) 
            end where

            if (grid%FFT_type.eq.11) then
                call IFFT3D_1D_decomp_gamma(grid, result_fftw%G, result_fftw%R_gamma)
            else if (grid%FFT_type.eq.21) then
                call IFFT3D_pruned_1D_decomp_gamma(grid, result_fftw%G, result_fftw%R_gamma)
            else 
                print *, 'Error FFT grid not initialized recip_to_real'
                stop
            endif
            result_fftw%R=result_fftw%R_gamma
            
        else
            if(grid%FFT_type.eq.1) then 
                call fftw_mpi_execute_dft(grid%plan_reverse, in_fftw%G, result_fftw%R)
            else if (grid%FFT_type.eq.11) then
                call IFFT3D_1D_decomp(grid, in_fftw%G, result_fftw%R)
            else if (grid%FFT_type.eq.21) then
                call IFFT3D_pruned_1D_decomp(grid, in_fftw%G, result_fftw%R)
            else 
                print *, 'Error FFT grid not initialized recip_to_real'
                stop
            endif
        endif
    end subroutine

    
    subroutine FFT3D_pruned_1D_decomp(grid, R_in, G_out)
        use, intrinsic :: iso_c_binding
        use grids_type,only: grid_struct
        use constants, only : i_
        use mpi
        include 'fftw3-mpi.f03'
        complex(C_DOUBLE_COMPLEX), intent(inout) :: R_in(:,:,:)
        complex(C_DOUBLE_COMPLEX), intent(inout) :: G_out(:,:,:)
        type(grid_struct) :: grid

        integer :: i, j,k
        integer :: nx_lot, nxy_lot, nyz_lot
        
        integer :: ix, iy, iz, ixy, iyz
        integer :: ix1,ix2,iz1, iz2, iy1, iy2
   

        integer, allocatable :: requests(:)
        integer ::  ierr

    

        allocate(requests(grid%Nr_local(3)))

        grid%outYXz=0.0_dp
        grid%outzYX=0.0_dp
        grid%outyXz2=0.0_dp

        iz1=1; iz2=1; iy1=0; iy2=0
        do iyz=1, grid%Nr_local(3)*grid%Nr(2), int(grid%lotx)
            nyz_lot=min(int(grid%lotx),grid%Nr_local(3)*grid%Nr(2)-iyz+1)
            iz=(iyz-1)/grid%Nr(2) + 1
            iy=iyz-(iz-1)*grid%Nr(2)

            if(nyz_lot.lt.grid%lotx) then
                call fftw_execute_dft(grid%plan_fft_xyz_x_left, R_in(1:,iy,iz), grid%zw)
            else
                call fftw_execute_dft(grid%plan_fft_xyz_x_lot, R_in(1:,iy,iz), grid%zw)
            endif 

            k=0    
            do i=0,nyz_lot-1
                iy2=iy2+1
                if(iy2.gt.grid%Nr(2)) then 
                    iy2=1; iz2=iz2+1
                endif
                !posive G's
                do j=1,(grid%Ng(1)+1)/2 !nx/4+1
                    k=k+1
                    grid%outyXz2(iy2,j,iz2) = grid%zw(k)!/grid%Ng(1)
                enddo
                k=k+grid%Nr(1)-grid%Ng(1)
                do j=(grid%Ng(1)+1)/2 + 1, grid%Ng(1)!nx/4+2,grid%Ng(1)
                    k=k+1
                    grid%outyXz2(iy2,j,iz2) = grid%zw(k)!/grid%Ng(1)
                enddo
            enddo
        enddo

     !   do ix=1,grid%Ng(1)
       !     print *, ix, grid%outyXz2(1,ix,1)/grid%Nr(1); flush(6)
       ! enddo
        
        do iz=1, grid%Nr_local(3)
            ix1=0; ix2=0 
            do ix=1, grid%Ng(1), int(grid%loty) !ix is first x index of lot for y transforms
                nx_lot=min(int(grid%loty), grid%Ng(1)-ix+1)

                if(nx_lot.lt.grid%loty) then
                    call fftw_execute_dft(grid%plan_fft_yXz_y_left, grid%outyXz2(:,ix,iz), grid%zw)

                else
                    call fftw_execute_dft(grid%plan_fft_yXz_y_lot, grid%outyXz2(:,ix,iz), grid%zw)
                endif 
                !Do lot2 y FFT's in a group
                !output to zw in order:
                k=0
                do i=0, nx_lot-1
                    ix2=ix2+1
                    if(ix2.gt.grid%Ng(1)) ix2=1 
                    do j=1,(grid%Ng(2)+1)/2!grid%Nr(2)/4
                        k=k+1
                        grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)=grid%zw(k)
                    enddo
                    k=k+grid%Nr(2)-grid%Ng(2)
                    do j=(grid%Ng(2)+1)/2+1, grid%Ng(2)!grid%Nr(2)/4+1, grid%Ng(2)
                        k=k+1
                        grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)=grid%zw(k)
                    enddo
                enddo
            enddo
            call mpi_Ialltoall( grid%outYXz(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX,  &
            grid%buff(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX, grid%FFT_comm, requests(iz),ierr)
        enddo

        call MPI_WAITALL(grid%Nr_local(3), requests, MPI_STATUSES_IGNORE, ierr)
        
        do iz=1, grid%Nr_local(3)
            do i=1,grid%nproc_FFT
                do ixy=1, int(grid%Ng_local(1)*grid%Ng(2))
                    grid%outzYX(iz+(i-1)*grid%Nr_local(3),ixy)= grid%buff((i-1)*int(grid%Ng_local(1)*grid%Ng(2))+ixy,iz)
                enddo
            enddo
        enddo


        !do iy=1,grid%Ng(2)
        !    print *, iy, grid%outzYX(1,iy+(3-1)*grid%Ng(2))/grid%Nr(1)/grid%Nr(2); flush(6)
        !enddo

        iy=0;ix=1
        do ixy=1, grid%Ng_local(1)*grid%Ng(2), int(grid%lotz)
            nxy_lot=min(int(grid%lotz),grid%Ng_local(1)*grid%Ng(2)-ixy+1)
    
            if(nxy_lot.lt.grid%lotz) then
                call fftw_execute_dft(grid%plan_fft_zYX_z_left,  grid%outzYX(:,ixy), grid%zw)
            else
                call fftw_execute_dft(grid%plan_fft_zYX_z_lot,  grid%outzYX(:,ixy), grid%zw)
            endif 
            !Do lotx x FFT's in a group
            !output to zw in order:
            k=0    
            do i=0,nxy_lot-1
                iy=iy+1
                if(iy.gt.grid%Ng(2)) then 
                    iy=1; ix=ix+1
                endif
                !posive G's
                do j=1,(grid%Ng(3)+1)/2!grid%Nr(3)/4+1
                    k=k+1
                    G_out(j,iy,ix) = grid%zw(k)!/grid%Ng(1)
                enddo
                k=k+grid%Nr(3)-grid%Ng(3)
                do j=(grid%Ng(3)+1)/2 + 1,grid%Ng(3)!grid%Nr(3)/4+2,grid%Ng(3)
                    k=k+1
                    G_out(j,iy,ix) = grid%zw(k)!/grid%Ng(1)
                enddo
            enddo
        enddo

        deallocate(requests)
    end subroutine

    subroutine IFFT3D_pruned_1D_decomp(grid, G_in, R_out)
        use, intrinsic :: iso_c_binding
        use grids_type,only: grid_struct
        use constants, only : i_
        use mpi
        include 'fftw3-mpi.f03'
        complex(C_DOUBLE_COMPLEX) :: G_in(:,:,:)
        complex(C_DOUBLE_COMPLEX) :: R_out(:,:,:)
        type(grid_struct) :: grid

        integer :: i, j,k
        integer :: nx_lot, nxy_lot, nyz_lot
        
        integer :: ix, iy, iz, ixy, iyz
        integer :: ix1,ix2,iz1, iz2, iy1, iy2


        integer, allocatable :: requests(:)
        integer :: STATUS(MPI_STATUS_SIZE)
        integer ::  ierr
        logical, allocatable :: mpi_flag(:)

        allocate(requests(grid%Nr_local(3)))
        allocate(mpi_flag(grid%Nr_local(3)))

        mpi_flag(:)=.false.
        grid%outYXz=0.0_dp
        grid%outzYX=0.0_dp
        grid%outyXz2=0.0_dp

        iy=0;ix=1
        do ixy=1, grid%Ng_local(1)*grid%Ng(2), int(grid%lotz)
            nxy_lot=min(int(grid%lotz),grid%Ng_local(1)*grid%Ng(2)-ixy+1)
    
            !Do lotx x FFT's in a group
            !output to zw in order:
            k=0
            grid%zw=0.0_dp    
            do i=0,nxy_lot-1
                iy=iy+1
                if(iy.gt.grid%Ng(2)) then 
                    iy=1; ix=ix+1
                endif
                !posive G's
                do j=1,(grid%Ng(3)+1)/2!grid%Nr(3)/4+1
                    k=k+1
                    grid%zw(k)=G_in(j,iy,ix)
                enddo
                k=k+grid%Nr(3)-grid%Ng(3)
                do j=(grid%Ng(3)+1)/2 + 1,grid%Ng(3)!grid%Nr(3)/4+2,grid%Ng(3)
                    k=k+1
                    grid%zw(k)=G_in(j,iy,ix)
                enddo
            enddo
            if(nxy_lot.lt.grid%lotz) then
                call fftw_execute_dft(grid%plan_Ifft_zYX_z_left,  grid%zw, grid%outzYX(:,ixy))
            else
                call fftw_execute_dft(grid%plan_Ifft_zYX_z_lot,  grid%zw, grid%outzYX(:,ixy))
            endif 
        enddo

    !    if(sim%parallel%myid.eq.0.and.l.eq.1) then 
    !        do iy=1,grid%Ng(2)
    !           print *, iy, grid%outzYX(1,iy+(3-1)*grid%Ng(2)); flush(6)
    !        enddo
    !    endif
    !    call parallel_wait(sim%parallel)

        do iz=1, grid%Nr_local(3)
            do i=1,grid%nproc_FFT
                do ixy=1, int(grid%Ng_local(1)*grid%Ng(2))
                    grid%buff((i-1)*int(grid%Ng_local(1)*grid%Ng(2))+ixy,iz)=grid%outzYX(iz+(i-1)*grid%Nr_local(3),ixy)
                enddo
            enddo
            call mpi_Ialltoall(grid%buff(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX,  &
            grid%outYXz(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX, grid%FFT_comm, requests(iz),ierr)
        enddo

        mpi_flag=.false.
        do while(.not.all(mpi_flag(:)))
        do iz=1, grid%Nr_local(3)
            !call MPI_WAIT(requests(iz), STATUS, ierr)
            if(mpi_flag(iz)) cycle !already done
            call MPI_TEST(requests(iz), mpi_flag(iz),STATUS,ierr)
            if(.not.mpi_flag(iz)) cycle !not ready yet
            ix1=0; ix2=0 
            do ix=1, grid%Ng(1), int(grid%loty) !ix is first x index of lot for y transforms
                nx_lot=min(int(grid%loty), grid%Ng(1)-ix+1)
                k=0
                grid%zw=0.0_dp
                do i=0, nx_lot-1
                    ix2=ix2+1
                    if(ix2.gt.grid%Ng(1)) ix2=1 
                    do j=1,(grid%Ng(2)+1)/2!grid%Nr(2)/4
                        k=k+1
                        grid%zw(k)=grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)
                    enddo
                    k=k+grid%Nr(2)-grid%Ng(2)
                    do j=(grid%Ng(2)+1)/2+1, grid%Ng(2)!grid%Nr(2)/4+1, grid%Ng(2)
                        k=k+1
                        grid%zw(k)=grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)
                    enddo
                enddo
                iy=(ix-1)
                if(nx_lot.lt.grid%loty) then
                    call fftw_execute_dft(grid%plan_Ifft_yXz_y_left, grid%zw,  grid%outyXz2(:,ix,iz))

                else
                    call fftw_execute_dft(grid%plan_Ifft_yXz_y_lot, grid%zw,  grid%outyXz2(:,ix,iz))
                endif
            enddo
        enddo
        enddo

    !   if(sim%parallel%myid.eq.0.and.l.eq.1) then 
    !        do ix=1,grid%Ng(1)
    !            print *, ix, grid%outyXz2(1,ix,1); flush(6)
    !        enddo
    !    endif
    !    call parallel_wait(sim%parallel)

        iz1=1; iz2=1; iy1=0; iy2=0
        do iyz=1, grid%Nr_local(3)*grid%Nr(2), int(grid%lotx)
            nyz_lot=min(int(grid%lotx),grid%Nr_local(3)*grid%Nr(2)-iyz+1)
            k=0  
            grid%zw=0.0_dp  
            do i=0,nyz_lot-1
                iy2=iy2+1
                if(iy2.gt.grid%Nr(2)) then 
                    iy2=1; iz2=iz2+1
                endif
                !posive G's
                do j=1,(grid%Ng(1)+1)/2 !nx/4+1
                    k=k+1
                    grid%zw(k)=grid%outyXz2(iy2,j,iz2)
                enddo
                k=k+grid%Nr(1)-grid%Ng(1)
                do j=(grid%Ng(1)+1)/2 + 1, grid%Ng(1)!nx/4+2,grid%Ng(1)
                    k=k+1
                    grid%zw(k)=grid%outyXz2(iy2,j,iz2)
                enddo
            enddo

            iz=(iyz-1)/grid%Nr(2) + 1
            iy=iyz-(iz-1)*grid%Nr(2)
            if(nyz_lot.lt.grid%lotx) then
                call fftw_execute_dft(grid%plan_Ifft_xyz_x_left, grid%zw, R_out(1:,iy,iz))
            else
                call fftw_execute_dft(grid%plan_Ifft_xyz_x_lot, grid%zw, R_out(1:,iy,iz))
            endif 
        enddo

        deallocate(requests)
        deallocate(mpi_flag)

    end subroutine

    subroutine FFT3D_1D_decomp(grid, R_in, G_out)
        use, intrinsic :: iso_c_binding
        use grids_type,only: grid_struct
        use constants, only : i_
        use mpi
        include 'fftw3-mpi.f03'
        complex(C_DOUBLE_COMPLEX), intent(inout) :: R_in(:,:,:)
        complex(C_DOUBLE_COMPLEX), intent(inout) :: G_out(:,:,:)
        type(grid_struct) :: grid

        integer :: i, j,k
        integer :: nx_lot, nxy_lot, nyz_lot
        
        integer :: ix, iy, iz, ixy, iyz
        integer :: ix1,ix2,iz1, iz2, iy1, iy2
   

        integer, allocatable :: requests(:)
        integer ::  ierr

    

        allocate(requests(grid%Nr_local(3)))

        grid%outYXz=0.0_dp
        grid%outzYX=0.0_dp
        grid%outyXz2=0.0_dp
        G_out=0.0_dp

        iz1=1; iz2=1; iy1=0; iy2=0
        do iyz=1, grid%Nr_local(3)*grid%Nr(2), int(grid%lotx)
            nyz_lot=min(int(grid%lotx),grid%Nr_local(3)*grid%Nr(2)-iyz+1)
            iz=(iyz-1)/grid%Nr(2) + 1
            iy=iyz-(iz-1)*grid%Nr(2)

            if(nyz_lot.lt.grid%lotx) then
                call fftw_execute_dft(grid%plan_fft_xyz_x_left, R_in(1:,iy,iz), grid%zw)
            else
                call fftw_execute_dft(grid%plan_fft_xyz_x_lot, R_in(1:,iy,iz), grid%zw)
            endif 

            k=0    
            do i=0,nyz_lot-1
                iy2=iy2+1
                if(iy2.gt.grid%Nr(2)) then 
                    iy2=1; iz2=iz2+1
                endif
                !positive G's
                do j=1, (grid%Nr(1)+1)/2
                    k=k+1
                    grid%outyXz2(iy2,j,iz2)=grid%zw(k)
                enddo
                !negative G's
                do j=(grid%Nr(1)+1)/2 + 1 + grid%Ng(1) -  grid%Nr(1), grid%Ng(1)
                    k=k+1
                    grid%outyXz2(iy2,j,iz2)=grid%zw(k)
                enddo
            enddo
        enddo
        
        do iz=1, grid%Nr_local(3)
            ix1=0; ix2=0 
            do ix=1, grid%Ng(1), int(grid%loty) !ix is first x index of lot for y transforms
                nx_lot=min(int(grid%loty), grid%Ng(1)-ix+1)

                if(nx_lot.lt.grid%loty) then
                    call fftw_execute_dft(grid%plan_fft_yXz_y_left, grid%outyXz2(:,ix,iz), grid%zw)

                else
                    call fftw_execute_dft(grid%plan_fft_yXz_y_lot, grid%outyXz2(:,ix,iz), grid%zw)
                endif 
                !Do lot2 y FFT's in a group
                !output to zw in order:
                k=0
                do i=0, nx_lot-1
                    ix2=ix2+1
                    if(ix2.gt.grid%Ng(1)) ix2=1 
                    do j=1, (grid%Nr(2)+1)/2
                        k=k+1
                        grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)=grid%zw(k)
                    enddo
                    do j=(grid%Nr(2)+1)/2 + 1 + grid%Ng(2) -  grid%Nr(2), grid%Ng(2)
                        k=k+1
                        grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)=grid%zw(k)
                    enddo
                enddo
            enddo
            call mpi_Ialltoall( grid%outYXz(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX,  &
            grid%buff(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX, grid%FFT_comm, requests(iz),ierr)
        enddo

        call MPI_WAITALL(grid%Nr_local(3), requests, MPI_STATUSES_IGNORE, ierr)
        
        do iz=1, grid%Nr_local(3)
            do i=1,grid%nproc_FFT
                do ixy=1, int(grid%Ng_local(1)*grid%Ng(2))
                    grid%outzYX(iz+(i-1)*grid%Nr_local(3),ixy)= grid%buff((i-1)*int(grid%Ng_local(1)*grid%Ng(2))+ixy,iz)
                enddo
            enddo
        enddo

        !do iy=1,grid%Ng(2)
        !    print *, iy, grid%outzYX(1,iy+(3-1)*grid%Ng(2))/grid%Nr(1)/grid%Nr(2); flush(6)
        !enddo

        iy=0;ix=1
        do ixy=1, grid%Ng_local(1)*grid%Ng(2), int(grid%lotz)
            nxy_lot=min(int(grid%lotz),grid%Ng_local(1)*grid%Ng(2)-ixy+1)
    
            if(nxy_lot.lt.grid%lotz) then
                call fftw_execute_dft(grid%plan_fft_zYX_z_left,  grid%outzYX(:,ixy), grid%zw)
            else
                call fftw_execute_dft(grid%plan_fft_zYX_z_lot,  grid%outzYX(:,ixy), grid%zw)
            endif 
            !Do lotx x FFT's in a group
            !output to zw in order:
            k=0    
            do i=0,nxy_lot-1
                iy=iy+1
                if(iy.gt.grid%Ng(2)) then 
                    iy=1; ix=ix+1
                endif
                !positive G's
                do j=1, (grid%Nr(3)+1)/2
                    k=k+1
                    G_out(j,iy,ix)=grid%zw(k)
                enddo
                !negative G's
                do j=(grid%Nr(3)+1)/2 + 1 + grid%Ng(3) -  grid%Nr(3), grid%Ng(3)
                    k=k+1
                    G_out(j,iy,ix)=grid%zw(k)
                enddo
            enddo
        enddo

        deallocate(requests)
    end subroutine

    subroutine IFFT3D_1D_decomp(grid, G_in, R_out)
        use, intrinsic :: iso_c_binding
        use grids_type,only: grid_struct
        use constants, only : i_
        use mpi
        include 'fftw3-mpi.f03'
        complex(C_DOUBLE_COMPLEX) :: G_in(:,:,:)
        complex(C_DOUBLE_COMPLEX) :: R_out(:,:,:)
        type(grid_struct) :: grid

        integer :: i, j,k
        integer :: nx_lot, nxy_lot, nyz_lot
        
        integer :: ix, iy, iz, ixy, iyz
        integer :: ix1,ix2,iz1, iz2, iy1, iy2


        integer, allocatable :: requests(:)
        integer :: STATUS(MPI_STATUS_SIZE)
        integer ::  ierr
        logical, allocatable :: mpi_flag(:)

        allocate(requests(grid%Nr_local(3)))
        allocate(mpi_flag(grid%Nr_local(3)))

        mpi_flag(:)=.false.
        grid%outYXz=0.0_dp
        grid%outzYX=0.0_dp
        grid%outyXz2=0.0_dp

        iy=0;ix=1
        do ixy=1, grid%Ng_local(1)*grid%Ng(2), int(grid%lotz)
            nxy_lot=min(int(grid%lotz),grid%Ng_local(1)*grid%Ng(2)-ixy+1)
    
            !Do lotx x FFT's in a group
            !output to zw in order:
            k=0
            grid%zw=0.0_dp    
            do i=0,nxy_lot-1
                iy=iy+1
                if(iy.gt.grid%Ng(2)) then 
                    iy=1; ix=ix+1
                endif
                !positive G's
                do j=1, (grid%Nr(3)+1)/2
                    k=k+1
                    grid%zw(k)=G_in(j,iy,ix)
                enddo
                !negative G's
                do j=(grid%Nr(3)+1)/2 + 1 + grid%Ng(3) -  grid%Nr(3), grid%Ng(3)
                    k=k+1
                    grid%zw(k)=G_in(j,iy,ix)
                enddo
            enddo
            if(nxy_lot.lt.grid%lotz) then
                call fftw_execute_dft(grid%plan_Ifft_zYX_z_left,  grid%zw, grid%outzYX(:,ixy))
            else
                call fftw_execute_dft(grid%plan_Ifft_zYX_z_lot,  grid%zw, grid%outzYX(:,ixy))
            endif 
        enddo

    !    if(sim%parallel%myid.eq.0.and.l.eq.1) then 
    !        do iy=1,grid%Ng(2)
    !           print *, iy, grid%outzYX(1,iy+(3-1)*grid%Ng(2)); flush(6)
    !        enddo
    !    endif
    !    call parallel_wait(sim%parallel)

        do iz=1, grid%Nr_local(3)
            do i=1,grid%nproc_FFT
                do ixy=1, int(grid%Ng_local(1)*grid%Ng(2))
                    grid%buff((i-1)*int(grid%Ng_local(1)*grid%Ng(2))+ixy,iz)=grid%outzYX(iz+(i-1)*grid%Nr_local(3),ixy)
                enddo
            enddo
            call mpi_Ialltoall(grid%buff(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX,  &
            grid%outYXz(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX, grid%FFT_comm, requests(iz),ierr)
        enddo

        mpi_flag=.false.
        do while(.not.all(mpi_flag(:)))
        do iz=1, grid%Nr_local(3)
            !call MPI_WAIT(requests(iz), STATUS, ierr)
            if(mpi_flag(iz)) cycle !already done
            call MPI_TEST(requests(iz), mpi_flag(iz),STATUS,ierr)
            if(.not.mpi_flag(iz)) cycle !not ready yet
            ix1=0; ix2=0 
            do ix=1, grid%Ng(1), int(grid%loty) !ix is first x index of lot for y transforms
                nx_lot=min(int(grid%loty), grid%Ng(1)-ix+1)
                k=0
                grid%zw=0.0_dp
                do i=0, nx_lot-1
                    ix2=ix2+1
                    if(ix2.gt.grid%Ng(1)) ix2=1 
                    do j=1, (grid%Nr(2)+1)/2
                        k=k+1
                        grid%zw(k)=grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)
                    enddo
                    do j=(grid%Nr(2)+1)/2 + 1 + grid%Ng(2) -  grid%Nr(2), grid%Ng(2)
                        k=k+1
                        grid%zw(k)=grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)
                    enddo
                enddo
                iy=(ix-1)
                if(nx_lot.lt.grid%loty) then
                    call fftw_execute_dft(grid%plan_Ifft_yXz_y_left, grid%zw,  grid%outyXz2(:,ix,iz))

                else
                    call fftw_execute_dft(grid%plan_Ifft_yXz_y_lot, grid%zw,  grid%outyXz2(:,ix,iz))
                endif
            enddo
        enddo
        enddo

    !   if(sim%parallel%myid.eq.0.and.l.eq.1) then 
    !        do ix=1,grid%Ng(1)
    !            print *, ix, grid%outyXz2(1,ix,1); flush(6)
    !        enddo
    !    endif
    !    call parallel_wait(sim%parallel)

        iz1=1; iz2=1; iy1=0; iy2=0
        do iyz=1, grid%Nr_local(3)*grid%Nr(2), int(grid%lotx)
            nyz_lot=min(int(grid%lotx),grid%Nr_local(3)*grid%Nr(2)-iyz+1)
            k=0  
            grid%zw=0.0_dp  
            do i=0,nyz_lot-1
                iy2=iy2+1
                if(iy2.gt.grid%Nr(2)) then 
                    iy2=1; iz2=iz2+1
                endif
                !positive G's
                do j=1, (grid%Nr(1)+1)/2
                    k=k+1
                    grid%zw(k)=grid%outyXz2(iy2,j,iz2)
                enddo
                !negative G's
                do j=(grid%Nr(1)+1)/2 + 1 + grid%Ng(1) -  grid%Nr(1), grid%Ng(1)
                    k=k+1
                    grid%zw(k)=grid%outyXz2(iy2,j,iz2)
                enddo
            enddo

            iz=(iyz-1)/grid%Nr(2) + 1
            iy=iyz-(iz-1)*grid%Nr(2)
            if(nyz_lot.lt.grid%lotx) then
                call fftw_execute_dft(grid%plan_Ifft_xyz_x_left, grid%zw, R_out(1:,iy,iz))
            else
                call fftw_execute_dft(grid%plan_Ifft_xyz_x_lot, grid%zw, R_out(1:,iy,iz))
            endif 
        enddo

        deallocate(requests)
        deallocate(mpi_flag)

    end subroutine

    subroutine FFT3D_pruned_1D_decomp_gamma(grid, R_in, G_out)
        use, intrinsic :: iso_c_binding
        use grids_type,only: grid_struct
        use constants, only : i_
        use mpi
        include 'fftw3-mpi.f03'
        real(dp), intent(inout) :: R_in(:,:,:)
        complex(C_DOUBLE_COMPLEX), intent(inout) :: G_out(:,:,:)
        type(grid_struct) :: grid

        integer :: i, j,k
        integer :: nx_lot, nxy_lot, nyz_lot
        
        integer :: ix, iy, iz, ixy, iyz
        integer :: ix1,ix2,iz1, iz2, iy1, iy2
   

        integer, allocatable :: requests(:)
        integer ::  ierr

    

        allocate(requests(grid%Nr_local(3)))

        grid%outYXz=0.0_dp
        grid%outzYX=0.0_dp
        grid%outyXz2=0.0_dp

        iz1=1; iz2=1; iy1=0; iy2=0
        do iyz=1, grid%Nr_local(3)*grid%Nr(2), int(grid%lotx)
            nyz_lot=min(int(grid%lotx),grid%Nr_local(3)*grid%Nr(2)-iyz+1)
            iz=(iyz-1)/grid%Nr(2) + 1
            iy=iyz-(iz-1)*grid%Nr(2)

            if(nyz_lot.lt.grid%lotx) then
                call fftw_execute_dft_r2c(grid%plan_fft_xyz_x_left, R_in(1:,iy,iz), grid%zw)
            else
                call fftw_execute_dft_r2c(grid%plan_fft_xyz_x_lot, R_in(1:,iy,iz), grid%zw)
            endif 

            k=0    
            do i=0,nyz_lot-1
                iy2=iy2+1
                if(iy2.gt.grid%Nr(2)) then 
                    iy2=1; iz2=iz2+1
                endif
                !posive G's
                do j=1, grid%Ng(1) !nx/4+1
                    k=k+1
                    grid%outyXz2(iy2,j,iz2) = grid%zw(k)!/grid%Ng(1)
                enddo
                k=k+(grid%Nr(1)/2+1)-grid%Ng(1)
            enddo
        enddo

     !   do ix=1,grid%Ng(1)
       !     print *, ix, grid%outyXz2(1,ix,1)/grid%Nr(1); flush(6)
       ! enddo
        
        do iz=1, grid%Nr_local(3)
            ix1=0; ix2=0 
            do ix=1, grid%Ng(1), int(grid%loty) !ix is first x index of lot for y transforms
                nx_lot=min(int(grid%loty), grid%Ng(1)-ix+1)

                if(nx_lot.lt.grid%loty) then
                    call fftw_execute_dft(grid%plan_fft_yXz_y_left, grid%outyXz2(:,ix,iz), grid%zw)

                else
                    call fftw_execute_dft(grid%plan_fft_yXz_y_lot, grid%outyXz2(:,ix,iz), grid%zw)
                endif 
                !Do lot2 y FFT's in a group
                !output to zw in order:
                k=0
                do i=0, nx_lot-1
                    ix2=ix2+1
                    if(ix2.gt.grid%Ng(1)) ix2=1 
                    do j=1,(grid%Ng(2)+1)/2!grid%Nr(2)/4
                        k=k+1
                        grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)=grid%zw(k)
                    enddo
                    k=k+grid%Nr(2)-grid%Ng(2)
                    do j=(grid%Ng(2)+1)/2+1, grid%Ng(2)!grid%Nr(2)/4+1, grid%Ng(2)
                        k=k+1
                        grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)=grid%zw(k)
                    enddo
                enddo
            enddo
            call mpi_Ialltoall( grid%outYXz(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX,  &
            grid%buff(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX, grid%FFT_comm, requests(iz),ierr)
        enddo

        call MPI_WAITALL(grid%Nr_local(3), requests, MPI_STATUSES_IGNORE, ierr)
        
        do iz=1, grid%Nr_local(3)
            do i=1,grid%nproc_FFT
                do ixy=1, int(grid%Ng_local(1)*grid%Ng(2))
                    grid%outzYX(iz+(i-1)*grid%Nr_local(3),ixy)= grid%buff((i-1)*int(grid%Ng_local(1)*grid%Ng(2))+ixy,iz)
                enddo
            enddo
        enddo


        !do iy=1,grid%Ng(2)
        !    print *, iy, grid%outzYX(1,iy+(3-1)*grid%Ng(2))/grid%Nr(1)/grid%Nr(2); flush(6)
        !enddo

        iy=0;ix=1
        do ixy=1, grid%Ng_local(1)*grid%Ng(2), int(grid%lotz)
            nxy_lot=min(int(grid%lotz),grid%Ng_local(1)*grid%Ng(2)-ixy+1)
    
            if(nxy_lot.lt.grid%lotz) then
                call fftw_execute_dft(grid%plan_fft_zYX_z_left,  grid%outzYX(:,ixy), grid%zw)
            else
                call fftw_execute_dft(grid%plan_fft_zYX_z_lot,  grid%outzYX(:,ixy), grid%zw)
            endif 
            !Do lotx x FFT's in a group
            !output to zw in order:
            k=0    
            do i=0,nxy_lot-1
                iy=iy+1
                if(iy.gt.grid%Ng(2)) then 
                    iy=1; ix=ix+1
                endif
                !posive G's
                do j=1,(grid%Ng(3)+1)/2!grid%Nr(3)/4+1
                    k=k+1
                    G_out(j,iy,ix) = grid%zw(k)!/grid%Ng(1)
                enddo
                k=k+grid%Nr(3)-grid%Ng(3)
                do j=(grid%Ng(3)+1)/2 + 1,grid%Ng(3)!grid%Nr(3)/4+2,grid%Ng(3)
                    k=k+1
                    G_out(j,iy,ix) = grid%zw(k)!/grid%Ng(1)
                enddo
            enddo
        enddo

        deallocate(requests)
    end subroutine

    subroutine IFFT3D_pruned_1D_decomp_gamma(grid, G_in, R_out)
        use, intrinsic :: iso_c_binding
        use grids_type,only: grid_struct
        use constants, only : i_
        use mpi
        include 'fftw3-mpi.f03'
        complex(C_DOUBLE_COMPLEX) :: G_in(:,:,:)
        real(dp) :: R_out(:,:,:)
        type(grid_struct) :: grid

        integer :: i, j,k
        integer :: nx_lot, nxy_lot, nyz_lot
        
        integer :: ix, iy, iz, ixy, iyz
        integer :: ix1,ix2,iz1, iz2, iy1, iy2


        integer, allocatable :: requests(:)
        integer :: STATUS(MPI_STATUS_SIZE)
        integer ::  ierr
        logical, allocatable :: mpi_flag(:)

        allocate(requests(grid%Nr_local(3)))
        allocate(mpi_flag(grid%Nr_local(3)))

        mpi_flag(:)=.false.
        grid%outYXz=0.0_dp
        grid%outzYX=0.0_dp
        grid%outyXz2=0.0_dp

        iy=0;ix=1
        do ixy=1, grid%Ng_local(1)*grid%Ng(2), int(grid%lotz)
            nxy_lot=min(int(grid%lotz),grid%Ng_local(1)*grid%Ng(2)-ixy+1)
    
            !Do lotx x FFT's in a group
            !output to zw in order:
            k=0
            grid%zw=0.0_dp    
            do i=0,nxy_lot-1
                iy=iy+1
                if(iy.gt.grid%Ng(2)) then 
                    iy=1; ix=ix+1
                endif
                !posive G's
                do j=1,(grid%Ng(3)+1)/2!grid%Nr(3)/4+1
                    k=k+1
                    grid%zw(k)=G_in(j,iy,ix)
                enddo
                k=k+grid%Nr(3)-grid%Ng(3)
                do j=(grid%Ng(3)+1)/2 + 1,grid%Ng(3)!grid%Nr(3)/4+2,grid%Ng(3)
                    k=k+1
                    grid%zw(k)=G_in(j,iy,ix)
                enddo
            enddo
            if(nxy_lot.lt.grid%lotz) then
                call fftw_execute_dft(grid%plan_Ifft_zYX_z_left,  grid%zw, grid%outzYX(:,ixy))
            else
                call fftw_execute_dft(grid%plan_Ifft_zYX_z_lot,  grid%zw, grid%outzYX(:,ixy))
            endif 
        enddo

    !    if(sim%parallel%myid.eq.0.and.l.eq.1) then 
    !        do iy=1,grid%Ng(2)
    !           print *, iy, grid%outzYX(1,iy+(3-1)*grid%Ng(2)); flush(6)
    !        enddo
    !    endif
    !    call parallel_wait(sim%parallel)

        do iz=1, grid%Nr_local(3)
            do i=1,grid%nproc_FFT
                do ixy=1, int(grid%Ng_local(1)*grid%Ng(2))
                    grid%buff((i-1)*int(grid%Ng_local(1)*grid%Ng(2))+ixy,iz)=grid%outzYX(iz+(i-1)*grid%Nr_local(3),ixy)
                enddo
            enddo
            call mpi_Ialltoall(grid%buff(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX,  &
            grid%outYXz(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX, grid%FFT_comm, requests(iz),ierr)
        enddo

        mpi_flag=.false.
        do while(.not.all(mpi_flag(:)))
        do iz=1, grid%Nr_local(3)
            !call MPI_WAIT(requests(iz), STATUS, ierr)
            if(mpi_flag(iz)) cycle !already done
            call MPI_TEST(requests(iz), mpi_flag(iz),STATUS,ierr)
            if(.not.mpi_flag(iz)) cycle !not ready yet
            ix1=0; ix2=0 
            do ix=1, grid%Ng(1), int(grid%loty) !ix is first x index of lot for y transforms
                nx_lot=min(int(grid%loty), grid%Ng(1)-ix+1)
                k=0
                grid%zw=0.0_dp
                do i=0, nx_lot-1
                    ix2=ix2+1
                    if(ix2.gt.grid%Ng(1)) ix2=1 
                    do j=1,(grid%Ng(2)+1)/2!grid%Nr(2)/4
                        k=k+1
                        grid%zw(k)=grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)
                    enddo
                    k=k+grid%Nr(2)-grid%Ng(2)
                    do j=(grid%Ng(2)+1)/2+1, grid%Ng(2)!grid%Nr(2)/4+1, grid%Ng(2)
                        k=k+1
                        grid%zw(k)=grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)
                    enddo
                enddo
                iy=(ix-1)
                if(nx_lot.lt.grid%loty) then
                    call fftw_execute_dft(grid%plan_Ifft_yXz_y_left, grid%zw,  grid%outyXz2(:,ix,iz))

                else
                    call fftw_execute_dft(grid%plan_Ifft_yXz_y_lot, grid%zw,  grid%outyXz2(:,ix,iz))
                endif
            enddo
        enddo
        enddo

    !   if(sim%parallel%myid.eq.0.and.l.eq.1) then 
    !        do ix=1,grid%Ng(1)
    !            print *, ix, grid%outyXz2(1,ix,1); flush(6)
    !        enddo
    !    endif
    !    call parallel_wait(sim%parallel)

        iz1=1; iz2=1; iy1=0; iy2=0
        do iyz=1, grid%Nr_local(3)*grid%Nr(2), int(grid%lotx)
            nyz_lot=min(int(grid%lotx),grid%Nr_local(3)*grid%Nr(2)-iyz+1)
            k=0  
            grid%zw=0.0_dp  
            do i=0,nyz_lot-1
                iy2=iy2+1
                if(iy2.gt.grid%Nr(2)) then 
                    iy2=1; iz2=iz2+1
                endif
                !posive G's
                do j=1, grid%Ng(1) !nx/4+1
                    k=k+1
                    grid%zw(k)=grid%outyXz2(iy2,j,iz2)
                enddo
                k=k+grid%Nr(1)/2+1-grid%Ng(1)
            enddo

            iz=(iyz-1)/grid%Nr(2) + 1
            iy=iyz-(iz-1)*grid%Nr(2)
            if(nyz_lot.lt.grid%lotx) then
                call fftw_execute_dft_c2r(grid%plan_Ifft_xyz_x_left, grid%zw, R_out(1:,iy,iz))
            else
                call fftw_execute_dft_c2r(grid%plan_Ifft_xyz_x_lot, grid%zw, R_out(1:,iy,iz))
            endif 
        enddo

        deallocate(requests)
        deallocate(mpi_flag)

    end subroutine

    subroutine FFT3D_1D_decomp_gamma(grid, R_in, G_out)
        use, intrinsic :: iso_c_binding
        use grids_type,only: grid_struct
        use constants, only : i_
        use mpi
        include 'fftw3-mpi.f03'
        real(dp), intent(inout) :: R_in(:,:,:)
        complex(C_DOUBLE_COMPLEX), intent(inout) :: G_out(:,:,:)
        type(grid_struct) :: grid

        integer :: i, j,k
        integer :: nx_lot, nxy_lot, nyz_lot
        
        integer :: ix, iy, iz, ixy, iyz
        integer :: ix1,ix2,iz1, iz2, iy1, iy2
   

        integer, allocatable :: requests(:)
        integer ::  ierr

    

        allocate(requests(grid%Nr_local(3)))

        grid%outYXz=0.0_dp
        grid%outzYX=0.0_dp
        grid%outyXz2=0.0_dp
        G_out=0.0_dp

        iz1=1; iz2=1; iy1=0; iy2=0
        do iyz=1, grid%Nr_local(3)*grid%Nr(2), int(grid%lotx)
            nyz_lot=min(int(grid%lotx),grid%Nr_local(3)*grid%Nr(2)-iyz+1)
            iz=(iyz-1)/grid%Nr(2) + 1
            iy=iyz-(iz-1)*grid%Nr(2)

            if(nyz_lot.lt.grid%lotx) then
                call fftw_execute_dft_r2c(grid%plan_fft_xyz_x_left, R_in(1:,iy,iz), grid%zw)
            else
                call fftw_execute_dft_r2c(grid%plan_fft_xyz_x_lot, R_in(1:,iy,iz), grid%zw)
            endif 

            k=0    
            do i=0,nyz_lot-1
                iy2=iy2+1
                if(iy2.gt.grid%Nr(2)) then 
                    iy2=1; iz2=iz2+1
                endif
                !positive G's
                do j=1, (grid%Nr(1))/2!+1
                    k=k+1
                    grid%outyXz2(iy2,j,iz2)=grid%zw(k)
                enddo
                k=k+1
            enddo
        enddo

        !do ix=1,grid%Ng(1)
        !    print *, 'step 1:', ix, grid%outyXz2(1,ix,:)/grid%Nr(1); flush(6)
        !enddo
        !do iz=1,grid%Nr_local(3)
        !    print *, 'step 1b:', ix, grid%outyXz2(1,2,iz)/grid%Nr(1); flush(6)
       ! enddo



        
        do iz=1, grid%Nr_local(3)
            ix1=0; ix2=0 
            do ix=1, grid%Ng(1), int(grid%loty) !ix is first x index of lot for y transforms
                nx_lot=min(int(grid%loty), grid%Ng(1)-ix+1)

                if(nx_lot.lt.grid%loty) then
                    call fftw_execute_dft(grid%plan_fft_yXz_y_left, grid%outyXz2(:,ix,iz), grid%zw)

                else
                    call fftw_execute_dft(grid%plan_fft_yXz_y_lot, grid%outyXz2(:,ix,iz), grid%zw)
                endif 
                !Do lot2 y FFT's in a group
                !output to zw in order:
                k=0
                do i=0, nx_lot-1
                    ix2=ix2+1
                    if(ix2.gt.grid%Ng(1)) ix2=1 
                    do j=1, (grid%Nr(2)+1)/2
                        k=k+1
                        grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)=grid%zw(k)
                    enddo
                    do j=(grid%Nr(2)+1)/2 + 1 + grid%Ng(2) -  grid%Nr(2), grid%Ng(2)
                        k=k+1
                        grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)=grid%zw(k)
                    enddo
                enddo
            enddo
            call mpi_Ialltoall( grid%outYXz(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX,  &
            grid%buff(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX, grid%FFT_comm, requests(iz),ierr)
        enddo

        call MPI_WAITALL(grid%Nr_local(3), requests, MPI_STATUSES_IGNORE, ierr)
        
        do iz=1, grid%Nr_local(3)
            do i=1,grid%nproc_FFT
                do ixy=1, int(grid%Ng_local(1)*grid%Ng(2))
                    grid%outzYX(iz+(i-1)*grid%Nr_local(3),ixy)= grid%buff((i-1)*int(grid%Ng_local(1)*grid%Ng(2))+ixy,iz)
                enddo
            enddo
        enddo

        !do iy=1,grid%Ng(2)
        !    print *, iy, grid%outzYX(1,iy+(2-1)*grid%Ng(2))/grid%Nr(1)/grid%Nr(2); flush(6)
        !enddo

        iy=0;ix=1
        do ixy=1, grid%Ng_local(1)*grid%Ng(2), int(grid%lotz)
            nxy_lot=min(int(grid%lotz),grid%Ng_local(1)*grid%Ng(2)-ixy+1)
    
            if(nxy_lot.lt.grid%lotz) then
                call fftw_execute_dft(grid%plan_fft_zYX_z_left,  grid%outzYX(:,ixy), grid%zw)
            else
                call fftw_execute_dft(grid%plan_fft_zYX_z_lot,  grid%outzYX(:,ixy), grid%zw)
            endif 
            !Do lotx x FFT's in a group
            !output to zw in order:
            k=0    
            do i=0,nxy_lot-1
                iy=iy+1
                if(iy.gt.grid%Ng(2)) then 
                    iy=1; ix=ix+1
                endif
                !positive G's
                do j=1, (grid%Nr(3)+1)/2
                    k=k+1
                    G_out(j,iy,ix)=grid%zw(k)
                enddo
                !negative G's
                do j=(grid%Nr(3)+1)/2 + 1 + grid%Ng(3) -  grid%Nr(3), grid%Ng(3)
                    k=k+1
                    G_out(j,iy,ix)=grid%zw(k)
                enddo
            enddo
        enddo

        deallocate(requests)
    end subroutine

    subroutine IFFT3D_1D_decomp_gamma(grid, G_in, R_out)
        use, intrinsic :: iso_c_binding
        use grids_type,only: grid_struct
        use constants, only : i_
        use mpi
        include 'fftw3-mpi.f03'
        complex(C_DOUBLE_COMPLEX) :: G_in(:,:,:)
        real(dp) :: R_out(:,:,:)
        type(grid_struct) :: grid

        integer :: i, j,k
        integer :: nx_lot, nxy_lot, nyz_lot
        
        integer :: ix, iy, iz, ixy, iyz
        integer :: ix1,ix2,iz1, iz2, iy1, iy2


        integer, allocatable :: requests(:)
        integer :: STATUS(MPI_STATUS_SIZE)
        integer ::  ierr
        logical, allocatable :: mpi_flag(:)

        allocate(requests(grid%Nr_local(3)))
        allocate(mpi_flag(grid%Nr_local(3)))

        mpi_flag(:)=.false.
        grid%outYXz=0.0_dp
        grid%outzYX=0.0_dp
        grid%outyXz2=0.0_dp

        iy=0;ix=1
        do ixy=1, grid%Ng_local(1)*grid%Ng(2), int(grid%lotz)
            nxy_lot=min(int(grid%lotz),grid%Ng_local(1)*grid%Ng(2)-ixy+1)
    
            !Do lotx x FFT's in a group
            !output to zw in order:
            k=0
            grid%zw=0.0_dp    
            do i=0,nxy_lot-1
                iy=iy+1
                if(iy.gt.grid%Ng(2)) then 
                    iy=1; ix=ix+1
                endif
                !positive G's
                do j=1, (grid%Nr(3)+1)/2
                    k=k+1
                    grid%zw(k)=G_in(j,iy,ix)
                enddo
                !negative G's
                do j=(grid%Nr(3)+1)/2 + 1 + grid%Ng(3) -  grid%Nr(3), grid%Ng(3)
                    k=k+1
                    grid%zw(k)=G_in(j,iy,ix)
                enddo
            enddo
            if(nxy_lot.lt.grid%lotz) then
                call fftw_execute_dft(grid%plan_Ifft_zYX_z_left,  grid%zw, grid%outzYX(:,ixy))
            else
                call fftw_execute_dft(grid%plan_Ifft_zYX_z_lot,  grid%zw, grid%outzYX(:,ixy))
            endif 
        enddo

       ! do iy=1,grid%Ng(2)
        !    print *, iy, grid%outzYX(1,iy+(2-1)*grid%Ng(2)); flush(6)
       ! enddo

        do iz=1, grid%Nr_local(3)
            do i=1,grid%nproc_FFT
                do ixy=1, int(grid%Ng_local(1)*grid%Ng(2))
                    grid%buff((i-1)*int(grid%Ng_local(1)*grid%Ng(2))+ixy,iz)=grid%outzYX(iz+(i-1)*grid%Nr_local(3),ixy)
                enddo
            enddo
            call mpi_Ialltoall(grid%buff(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX,  &
            grid%outYXz(:,iz), int(grid%Ng_local(1)*grid%Ng(2)), MPI_DOUBLE_COMPLEX, grid%FFT_comm, requests(iz),ierr)
        enddo

        mpi_flag=.false.
        do while(.not.all(mpi_flag(:)))
        do iz=1, grid%Nr_local(3)
            !call MPI_WAIT(requests(iz), STATUS, ierr)
            if(mpi_flag(iz)) cycle !already done
            call MPI_TEST(requests(iz), mpi_flag(iz),STATUS,ierr)
            if(.not.mpi_flag(iz)) cycle !not ready yet
            ix1=0; ix2=0 
            do ix=1, grid%Ng(1), int(grid%loty) !ix is first x index of lot for y transforms
                nx_lot=min(int(grid%loty), grid%Ng(1)-ix+1)
                k=0
                grid%zw=0.0_dp
                do i=0, nx_lot-1
                    ix2=ix2+1
                    if(ix2.gt.grid%Ng(1)) ix2=1 
                    do j=1, (grid%Nr(2)+1)/2
                        k=k+1
                        grid%zw(k)=grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)
                    enddo
                    do j=(grid%Nr(2)+1)/2 + 1 + grid%Ng(2) -  grid%Nr(2), grid%Ng(2)
                        k=k+1
                        grid%zw(k)=grid%outYXz(j+(ix2-1)*grid%Ng(2),iz)
                    enddo
                enddo
                iy=(ix-1)
                if(nx_lot.lt.grid%loty) then
                    call fftw_execute_dft(grid%plan_Ifft_yXz_y_left, grid%zw,  grid%outyXz2(:,ix,iz))

                else
                    call fftw_execute_dft(grid%plan_Ifft_yXz_y_lot, grid%zw,  grid%outyXz2(:,ix,iz))
                endif
            enddo
        enddo
        enddo

        !do ix=1,grid%Ng(1)
        !    print *, ix, grid%outyXz2(1,ix,1); flush(6)
       ! enddo

        iz1=1; iz2=1; iy1=0; iy2=0
        do iyz=1, grid%Nr_local(3)*grid%Nr(2), int(grid%lotx)
            nyz_lot=min(int(grid%lotx),grid%Nr_local(3)*grid%Nr(2)-iyz+1)
            k=0  
            grid%zw=0.0_dp  
            do i=0,nyz_lot-1
                iy2=iy2+1
                if(iy2.gt.grid%Nr(2)) then 
                    iy2=1; iz2=iz2+1
                endif
                !positive G's
                do j=1, (grid%Nr(1))/2!+1
                    k=k+1
                    grid%zw(k)=grid%outyXz2(iy2,j,iz2)
                enddo
                k=k+1
            enddo

            iz=(iyz-1)/grid%Nr(2) + 1
            iy=iyz-(iz-1)*grid%Nr(2)

            !print *, 'HERE', iy, iz, size(R_out,1), size(R_out,2), size(R_out,3)  
            if(nyz_lot.lt.grid%lotx) then
                call fftw_execute_dft_c2r(grid%plan_Ifft_xyz_x_left, grid%zw, R_out(1:,iy,iz))
            else
                call fftw_execute_dft_c2r(grid%plan_Ifft_xyz_x_lot, grid%zw, R_out(1:,iy,iz))
            endif 

        enddo

        deallocate(requests)
        deallocate(mpi_flag)

    end subroutine

end module