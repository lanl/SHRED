module read_in
use types, only : dp
implicit none

public :: read_in_density

contains

subroutine read_in_density(density,grids,i_grid,parallel)
    use odp_type, only: spin_DenPot_struct
    use grids_type, only: grid_struct
    use fft, only: real_to_recip 
    use distribute_field, only: distribute_R
    use parallel_type, only : parallel_struct
    use parallel_mod, only: parallel_task
    use fft, only : real_to_recip, recip_to_real
    use grids_mod, only :  allocate_full_fields_R, allocate_local_fields_R


    type(spin_DenPot_struct), intent(inout) :: density
    type(grid_struct), intent(in), target :: grids(:)
    type(grid_struct), pointer :: grid

    type(parallel_struct), intent(in) :: parallel
    integer, intent(in) :: i_grid
    

    real(dp), allocatable :: Den_in(:,:,:)
    real(dp), allocatable :: Den_local(:,:,:)

    integer :: u, i

    grid=>grids(i_grid)

    if(parallel%myid.eq.0) then
        open(newunit=u, file='density.dat', status="old")
        call allocate_full_fields_R(Den_in, grid)
    endif
    call allocate_local_fields_R(Den_local, grid)

    do i=1,density%n_s
        if(parallel%myid.eq.0) then
            read(u,*) Den_in
        endif
        if((parallel%myid_diff_k.eq.0).and.(parallel%myid_diff_s.eq.0).and.(parallel%my_band_group.eq.0)) then
            call distribute_R(0, Den_in, Den_local, 'band', grid, parallel)
            density%of(i)%R=Den_local
        endif

        call parallel_task('bcast',density%of(i)%R, parallel, 'diff_b', root=0)
        call parallel_task('bcast',density%of(i)%R, parallel, 'diff_k', root=0)
        call parallel_task('bcast',density%of(i)%R, parallel, 'diff_s', root=0)

        call real_to_recip(density%of(i), grids)
        density%of(i)%G=density%of(i)%G*grid%cutden
        call recip_to_real(density%of(i), grids)


    enddo

end subroutine

end module