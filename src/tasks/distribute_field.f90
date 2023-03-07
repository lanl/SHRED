module distribute_field
use types, only : dp
implicit none

interface distribute_R
        module procedure distribute_real
        module procedure distribute_complex
end interface

interface collate_R
        module procedure collate_real
        module procedure collate_complex
end interface

    contains
    subroutine distribute_real(root, dsource, dtarget, comm, grid, parallel)
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use parallel_mod, only: parallel_task

        ! Distribute data from global dsource on rank 'root' to parallel dtarget
        type(grid_struct), intent(in) ::  grid
        type(parallel_struct), intent(in) ::  parallel
        character(len=*), intent(in) :: comm
        integer, intent(in) :: root
        real(dp), intent(in) :: dsource(:,:,:)
        real(dp), intent(out) :: dtarget(:,:,:)
        real(dp) :: tmp(size(dtarget,1),size(dtarget,2),size(dtarget,3))
        integer :: target_id, max_id, xyz(3)

        max_id=parallel%nproc_band
        do target_id = 0, max_id-1
            xyz=grid%myxyz_r
            !Send position requested from id to root
            if(target_id.ne.root) then
                if (parallel%myid_band == root) then
                    call parallel_task('recv', xyz, parallel, comm, root=target_id, tag=0)
                else if (parallel%myid_band == target_id) then
                    call parallel_task('send', xyz, parallel, comm, to=root, tag=0)
                endif
            endif
            !Send data at positions from root to id
            if(target_id.ne.root) then
                if (parallel%myid_band == root) then
                    tmp = dsource(xyz(1)*grid%Nr_local(1)+1:(xyz(1)+1)*grid%Nr_local(1), &
                                  xyz(2)*grid%Nr_local(2)+1:(xyz(2)+1)*grid%Nr_local(2), &
                                  xyz(3)*grid%Nr_local(3)+1:(xyz(3)+1)*grid%Nr_local(3))
                    call parallel_task('send', tmp, parallel, comm, to=target_id, tag=0)
                else if (parallel%myid_band == target_id) then
                    call parallel_task('recv', dtarget, parallel, comm, root=root, tag=0)
                end if
            else
                if (parallel%myid_band == root) then
                    dtarget = dsource(xyz(1)*grid%Nr_local(1)+1:(xyz(1)+1)*grid%Nr_local(1), &
                                      xyz(2)*grid%Nr_local(2)+1:(xyz(2)+1)*grid%Nr_local(2), &
                                      xyz(3)*grid%Nr_local(3)+1:(xyz(3)+1)*grid%Nr_local(3))
                endif
            endif
        enddo
    end subroutine

    subroutine distribute_complex(root, dsource, dtarget, comm, grid, parallel)
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use parallel_mod, only: parallel_task

        ! Distribute data from global dsource on rank 'root' to parallel dtarget
        type(grid_struct), intent(in) ::  grid
        type(parallel_struct), intent(in) ::  parallel
        character(len=*), intent(in) :: comm
        integer, intent(in) :: root
        complex(dp), intent(in) :: dsource(:,:,:)
        complex(dp), intent(out) :: dtarget(:,:,:)
        complex(dp) :: tmp(size(dtarget,1),size(dtarget,2),size(dtarget,3))
        integer :: target_id, max_id, xyz(3)

        max_id=parallel%nproc_band
        do target_id = 0, max_id-1
            xyz=grid%myxyz_r
            !Send position requested from id to root
            if(target_id.ne.root) then
                if (parallel%myid_band == root) then
                    call parallel_task('recv', xyz, parallel, comm, root=target_id, tag=0)
                else if (parallel%myid_band == target_id) then
                    call parallel_task('send', xyz, parallel, comm, to=root, tag=0)
                endif
            endif
            !Send data at positions from root to id
            if(target_id.ne.root) then
                if (parallel%myid_band == root) then
                    tmp = dsource(xyz(1)*grid%Nr_local(1)+1:(xyz(1)+1)*grid%Nr_local(1), &
                                  xyz(2)*grid%Nr_local(2)+1:(xyz(2)+1)*grid%Nr_local(2), &
                                  xyz(3)*grid%Nr_local(3)+1:(xyz(3)+1)*grid%Nr_local(3))
                    call parallel_task('send', tmp, parallel, comm, to=target_id, tag=0)
                else if (parallel%myid_band == target_id) then
                    call parallel_task('recv', dtarget, parallel, comm, root=root, tag=0)
                end if
            else
                if (parallel%myid_band == root) then
                    dtarget = dsource(xyz(1)*grid%Nr_local(1)+1:(xyz(1)+1)*grid%Nr_local(1), &
                                      xyz(2)*grid%Nr_local(2)+1:(xyz(2)+1)*grid%Nr_local(2), &
                                      xyz(3)*grid%Nr_local(3)+1:(xyz(3)+1)*grid%Nr_local(3))
                endif
            endif
        enddo
    end subroutine

    subroutine distribute_G(root, dsource, dtarget, comm, grid, parallel)
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use parallel_mod, only: parallel_task

        ! Distribute data from global dsource on rank 'root' to parallel dtarget
        type(grid_struct), intent(in) ::  grid
        type(parallel_struct), intent(in) ::  parallel
        character(len=*), intent(in) :: comm
        integer, intent(in) :: root
        complex(dp), intent(in) :: dsource(:,:,:)
        complex(dp), intent(out) :: dtarget(:,:,:)
        complex(dp) :: tmp(size(dtarget,1),size(dtarget,2),size(dtarget,3))
        integer :: target_id, max_id, xyz(3), i,j,k, ijk_global(3), Ng_local(3), Ng(3)

        Ng_local = [size(dtarget,1), size(dtarget,2), size(dtarget,3)]
        Ng  = [size(dsource,1), size(dsource,2), size(dsource,3)]
        
        max_id=parallel%nproc_band
        do target_id = 0, max_id-1
            xyz=grid%myxyz
            !Send position requested from id to root
            if(target_id.ne.root) then
                if (parallel%myid_band == root) then
                    call parallel_task('recv', xyz, parallel, comm, root=target_id, tag=0)
                else if (parallel%myid_band == target_id) then
                    call parallel_task('send', xyz, parallel, comm, to=root, tag=0)
                endif
            endif
            !Send data at positions from root to id
            if(target_id.ne.root) then
                if (parallel%myid_band == root) then

                    if(grid%FFT_type.gt.1) then
                        do i = 1, Ng_local(1); do j = 1, Ng_local(2); do k = 1, Ng_local(3)
                            ijk_global = [i, j, k] + xyz(:)*Ng_local(:) &
                                - Ng(:)*nint((ijk_global-1.5_dp)/Ng)
                            tmp(k,j,i)=dsource(ijk_global(3),ijk_global(2),ijk_global(1)) 
                        end do; end do; end do
                    else
                        do k = 1, Ng_local(3); do j = 1, Ng_local(2); do i = 1, Ng_local(1)
                            ijk_global = [i, j, k] + xyz(:)*Ng_local(:) &
                                - Ng(:)*nint((ijk_global-1.5_dp)/Ng)
                            tmp(i,j,k)=dsource(ijk_global(1),ijk_global(2),ijk_global(3)) 
                        end do; end do; end do
                    endif
                    

                    call parallel_task('send', tmp, parallel, comm, to=target_id, tag=0)
                else if (parallel%myid_band == target_id) then
                    call parallel_task('recv', dtarget, parallel, comm, root=root, tag=0)
                end if
            else
                if (parallel%myid_band == root) then
                    if(grid%FFT_type.gt.1) then
                        do i = 1, Ng_local(1); do j = 1, Ng_local(2); do k = 1, Ng_local(3)
                            ijk_global = [i, j, k] + xyz(:)*Ng_local(:) &
                                - Ng(:)*nint((ijk_global-1.5_dp)/Ng)
                            dtarget(k,j,i)=dsource(ijk_global(3),ijk_global(2),ijk_global(1)) 
                        end do; end do; end do
                    else
                        do k = 1, Ng_local(3); do j = 1, Ng_local(2); do i = 1, Ng_local(1)
                            ijk_global = [i, j, k] + xyz(:)*Ng_local(:) &
                                - Ng(:)*nint((ijk_global-1.5_dp)/Ng)
                            dtarget(i,j,k)=dsource(ijk_global(1),ijk_global(2),ijk_global(3)) 
                        end do; end do; end do
                    endif
                endif
            endif
        enddo
    end subroutine

    subroutine collate_real(root, dsource, dtarget, comm, grid, parallel)
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use parallel_mod, only: parallel_task

        ! Collect data from local dsources on ranks to global dtarget on 'root'
        type(grid_struct), intent(in) ::  grid
        type(parallel_struct), intent(in) ::  parallel
        character(len=*), intent(in) :: comm
        integer, intent(in) :: root
        real(dp), intent(in) :: dsource(:,:,:)
        real(dp), intent(out) :: dtarget(:,:,:)
        real(dp) :: tmp(size(dsource,1),size(dsource,2),size(dsource,3))
        integer :: source_id, max_id, xyz(3)

        max_id=parallel%nproc_band
        do source_id = 0, max_id-1
            xyz=grid%myxyz_r
            !Send position requested from id to root
            if(source_id.ne.root) then
                if (parallel%myid_band == root) then
                    call parallel_task('recv', xyz, parallel, comm, root=source_id, tag=0)
                else if (parallel%myid_band == source_id) then
                    call parallel_task('send', xyz, parallel, comm, to=root, tag=0)
                endif
            endif
            !Send data at positions from root to id
            if(source_id.ne.root) then
                if (parallel%myid_band == source_id) then
                    tmp = dsource
                    call parallel_task('send', tmp, parallel, comm, to=root, tag=0)
                else if (parallel%myid_band == root) then
                    call parallel_task('recv', tmp, parallel, comm, root=source_id, tag=0)
                    dtarget(xyz(1)*grid%Nr_local(1)+1:(xyz(1)+1)*grid%Nr_local(1), &
                            xyz(2)*grid%Nr_local(2)+1:(xyz(2)+1)*grid%Nr_local(2), &
                            xyz(3)*grid%Nr_local(3)+1:(xyz(3)+1)*grid%Nr_local(3))=tmp
                end if
            else
                if (parallel%myid_band == root) then
                    dtarget(xyz(1)*grid%Nr_local(1)+1:(xyz(1)+1)*grid%Nr_local(1), &
                            xyz(2)*grid%Nr_local(2)+1:(xyz(2)+1)*grid%Nr_local(2), &
                            xyz(3)*grid%Nr_local(3)+1:(xyz(3)+1)*grid%Nr_local(3))= dsource
                endif
            endif
        enddo
    end subroutine

    subroutine collate_complex(root, dsource, dtarget, comm, grid, parallel)
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use parallel_mod, only: parallel_task

        ! Collect data from local dsources on ranks to global dtarget on 'root'
        type(grid_struct), intent(in) ::  grid
        type(parallel_struct), intent(in) ::  parallel
        character(len=*), intent(in) :: comm
        integer, intent(in) :: root
        complex(dp), intent(in) :: dsource(:,:,:)
        complex(dp), intent(out) :: dtarget(:,:,:)
        complex(dp) :: tmp(size(dsource,1),size(dsource,2),size(dsource,3))
        integer :: source_id, max_id, xyz(3)

        max_id=parallel%nproc_band
        do source_id = 0, max_id-1
            xyz=grid%myxyz_r
            !Send position requested from id to root
            if(source_id.ne.root) then
                if (parallel%myid_band == root) then
                    call parallel_task('recv', xyz, parallel, comm, root=source_id, tag=0)
                else if (parallel%myid_band == source_id) then
                    call parallel_task('send', xyz, parallel, comm, to=root, tag=0)
                endif
            endif
            !Send data at positions from root to id
            if(source_id.ne.root) then
                if (parallel%myid_band == source_id) then
                    tmp = dsource
                    call parallel_task('send', tmp, parallel, comm, to=root, tag=0)
                else if (parallel%myid_band == root) then
                    call parallel_task('recv', tmp, parallel, comm, root=source_id, tag=0)
                    dtarget(xyz(1)*grid%Nr_local(1)+1:(xyz(1)+1)*grid%Nr_local(1), &
                            xyz(2)*grid%Nr_local(2)+1:(xyz(2)+1)*grid%Nr_local(2), &
                            xyz(3)*grid%Nr_local(3)+1:(xyz(3)+1)*grid%Nr_local(3))=tmp
                end if
            else
                if (parallel%myid_band == root) then
                    dtarget(xyz(1)*grid%Nr_local(1)+1:(xyz(1)+1)*grid%Nr_local(1), &
                            xyz(2)*grid%Nr_local(2)+1:(xyz(2)+1)*grid%Nr_local(2), &
                            xyz(3)*grid%Nr_local(3)+1:(xyz(3)+1)*grid%Nr_local(3))= dsource
                endif
            endif
        enddo
    end subroutine

    subroutine collate_G(root, dsource, dtarget, comm, grid, parallel)
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use parallel_mod, only: parallel_task

        ! Collect data from local dsources on ranks to global dtarget on 'root'
        type(grid_struct), intent(in) ::  grid
        type(parallel_struct), intent(in) ::  parallel
        character(len=*), intent(in) :: comm
        integer, intent(in) :: root
        complex(dp), intent(in) :: dsource(:,:,:)
        complex(dp), intent(out) :: dtarget(:,:,:)
        complex(dp) :: tmp(size(dsource,1),size(dsource,2),size(dsource,3))
        integer :: source_id, max_id, xyz(3), i,j,k, ijk_global(3), Ng_local(3), Ng(3)

        Ng_local = [size(dsource,1), size(dsource,2), size(dsource,3)]
        Ng  = [size(dtarget,1), size(dtarget,2), size(dtarget,3)]

        max_id=parallel%nproc_band
        do source_id = 0, max_id-1
            xyz=grid%myxyz
            !Send position requested from id to root
            if(source_id.ne.root) then
                if (parallel%myid_band == root) then
                    call parallel_task('recv', xyz, parallel, comm, root=source_id, tag=0)
                else if (parallel%myid_band == source_id) then
                    call parallel_task('send', xyz, parallel, comm, to=root, tag=0)
                endif
            endif
            !Send data at positions from root to id
            if(source_id.ne.root) then
                if (parallel%myid_band == source_id) then
                    tmp = dsource
                    call parallel_task('send', tmp, parallel, comm, to=root, tag=0)

                else if (parallel%myid_band == root) then
                    call parallel_task('recv', tmp, parallel, comm, root=source_id, tag=0)
                    if(grid%FFT_type.gt.1) then
                        do i = 1, Ng_local(1); do j = 1, Ng_local(2); do k = 1, Ng_local(3)
                            ijk_global = [i, j, k] + xyz(:)*Ng_local(:)
                            ijk_global = ijk_global - Ng(:)*nint((ijk_global-1.5_dp)/Ng)
                            dtarget(ijk_global(3),ijk_global(2),ijk_global(1)) = tmp(k,j,i)
                        end do; end do; end do
                    else
                        do k = 1, Ng_local(3); do j = 1, Ng_local(2); do i = 1, Ng_local(1)
                            ijk_global = [i, j, k] + xyz(:)*Ng_local(:)
                            ijk_global = ijk_global - Ng(:)*nint((ijk_global-1.5_dp)/Ng)
                            dtarget(ijk_global(1),ijk_global(2),ijk_global(3)) = tmp(i,j,k)
                        end do; end do; end do
                    endif
                end if
            else
                if (parallel%myid_band == root) then
                    if(grid%FFT_type.gt.1) then
                        do i= 1, Ng_local(1); do j = 1, Ng_local(2); do k = 1, Ng_local(3)
                            ijk_global = [i, j, k] + xyz(:)*Ng_local(:)
                            ijk_global = ijk_global - Ng(:)*nint((ijk_global-1.5_dp)/Ng)
                            dtarget(ijk_global(3),ijk_global(2),ijk_global(1)) = dsource(k,j,i)
                        end do; end do; end do 
                    else
                        do k = 1, Ng_local(3); do j = 1, Ng_local(2); do i = 1, Ng_local(1)
                            ijk_global = [i, j, k] + xyz(:)*Ng_local(:) 
                            ijk_global = ijk_global - Ng(:)*nint((ijk_global-1.5_dp)/Ng)
                            dtarget(ijk_global(1),ijk_global(2),ijk_global(3)) = dsource(i,j,k)
                        end do; end do; end do 
                    endif
                    
                endif
            endif
        enddo
    end subroutine
        
end module