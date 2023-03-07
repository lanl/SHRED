module odp
    use types, only : op, dp
implicit none

private 

public :: allocate_field, deallocate_field, allocate_DenPot, deallocate_DenPot, allocate_orbital, deallocate_orbital

interface allocate_DenPot
    module procedure allocate_no_spin_DenPot
    module procedure allocate_spin_DenPot
end interface

interface deallocate_DenPot
    module procedure deallocate_no_spin_DenPot
    module procedure deallocate_spin_DenPot
end interface



contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine allocates a field and FFTW plan for FFT/inverse FFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine allocate_field(psi, grid, parallel)

    use, intrinsic :: iso_c_binding
    use odp_type, only: field_struct
    use grids_type, only: grid_struct
    use parallel_type, only: parallel_struct
   
    include 'fftw3-mpi.f03'
    type(parallel_struct), intent(in) :: parallel
    type(grid_struct), intent(in) :: grid
    type(field_struct), intent(inout) :: psi

    integer(C_INTPTR_T) :: local_z, local_y, local_x

    local_x=grid%Nr_local(1)
    local_y=grid%Nr_local(2)
    local_z=grid%Nr_local(3)

    psi%p_R=fftw_alloc_complex(grid%alloc_local_R)
    call c_f_pointer(psi%p_R, psi%R, [local_x,local_y,local_z])

    if(grid%gamma) then
        psi%p_R_gamma=fftw_alloc_real(grid%alloc_local_R_gamma)
        call c_f_pointer(psi%p_R_gamma, psi%R_gamma, [local_x,local_y,local_z])
    endif

    local_x=grid%Ng_local(1)
    local_y=grid%Ng_local(2)
    local_z=grid%Ng_local(3)

    psi%p_G=fftw_alloc_complex(grid%alloc_local_G)
    if(grid%FFT_type.gt.1) then
        call c_f_pointer(psi%p_G, psi%G, [local_z,local_y,local_x])
    else
        call c_f_pointer(psi%p_G, psi%G, [local_x,local_y,local_z])
    endif


end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocates a field and destroys FFTW plan for FFT/inverse FFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine deallocate_field(psi, grid)
    use, intrinsic :: iso_c_binding
    use odp_type, only: field_struct
    use grids_type, only: grid_struct

    include 'fftw3-mpi.f03'
    type(field_struct), intent(inout) :: psi
    type(grid_struct), intent(in) :: grid


    !deallocate(psi%R)
    !deallocate(psi%G)
    call fftw_free(psi%p_R)
    call fftw_free(psi%p_G)
    if(grid%gamma) then
        call fftw_free(psi%p_R_gamma)
    endif
    nullify(psi%R)
    nullify(psi%G)
    nullify(psi%R_gamma)
end subroutine

subroutine allocate_no_spin_DenPot(den, grids, i_grid, mpi)

    use, intrinsic :: iso_c_binding
    use odp_type, only: DenPot_struct
    use grids_type, only: grid_struct
    use parallel_type, only: parallel_struct
   
    include 'fftw3-mpi.f03'
    type(parallel_struct), intent(in) :: mpi
    type(grid_struct), intent(in) :: grids(:)
    integer, intent(in) :: i_grid
    type(DenPot_struct), intent(inout) :: den
        
        den%of%grid=i_grid
        call allocate_field(den%of, grids(den%of%grid), mpi)

end subroutine

subroutine allocate_spin_DenPot(den, grids, i_grid, mpi)

    use, intrinsic :: iso_c_binding
    use odp_type, only: spin_DenPot_struct
    use grids_type, only: grid_struct
    use parallel_type, only: parallel_struct
    type(parallel_struct), intent(in) :: mpi
    type(grid_struct), intent(in) :: grids(:)
    integer, intent(in) :: i_grid
    type(spin_DenPot_struct), intent(inout) :: den
    integer :: i

    allocate(den%of(den%n_s))
    do i=1, den%n_s
        den%of(i)%grid=i_grid
        call allocate_field(den%of(i), grids(den%of(i)%grid), mpi)
    enddo
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine deallocates a field and destroys FFTW plan for FFT/inverse FFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine deallocate_spin_DenPot(den, grids)
    use, intrinsic :: iso_c_binding
    use odp_type, only: spin_DenPot_struct
    use grids_type, only: grid_struct
    type(spin_DenPot_struct), intent(inout) :: den
    type(grid_struct), intent(in) :: grids(:)
    integer :: i

    do i=1, den%n_s
        call deallocate_field(den%of(i), grids(den%of(i)%grid))
    enddo
    deallocate(den%of)
end subroutine

subroutine deallocate_no_spin_DenPot(den, grids)
    use, intrinsic :: iso_c_binding
    use odp_type, only: DenPot_struct
    use grids_type, only: grid_struct
    type(DenPot_struct), intent(inout) :: den
    type(grid_struct), intent(in) :: grids(:)

    call deallocate_field(den%of,grids(den%of%grid)) 
    
end subroutine

subroutine allocate_orbital(orb, grids, i_grid, parallel)
    use, intrinsic :: iso_c_binding
    use odp_type, only: orbital_struct
    use grids_type, only: grid_struct
    use parallel_type, only: parallel_struct
   
    include 'fftw3-mpi.f03'
    type(parallel_struct), intent(in) :: parallel
    type(grid_struct), intent(in) :: grids(:)
    type(orbital_struct), intent(inout) :: orb
    integer, intent(in) :: i_grid
    integer :: i
   
    allocate(orb%of(orb%n_spinor))
    allocate(orb%occ(orb%n_spinor))
    allocate(orb%weight(orb%n_spinor))
    allocate(orb%eig(orb%n_spinor))
    allocate(orb%filter(orb%n_spinor))

    do i=1, orb%n_spinor
        orb%of(i)%grid=i_grid
        call allocate_field(orb%of(i), grids(orb%of(i)%grid), parallel)
    enddo
end subroutine

subroutine deallocate_orbital(orb, grids)
    use, intrinsic :: iso_c_binding
    use odp_type, only: orbital_struct
    use grids_type, only: grid_struct

    type(orbital_struct), intent(inout) :: orb
    type(grid_struct), intent(in) :: grids(:)

    integer :: i

    do i=1, orb%n_spinor
        call deallocate_field(orb%of(i),grids(orb%of(i)%grid))
    enddo
    deallocate(orb%of)
    deallocate(orb%occ)
    deallocate(orb%weight)
    deallocate(orb%eig)
    deallocate(orb%filter)
end subroutine

end module