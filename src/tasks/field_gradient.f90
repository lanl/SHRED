module gradient_mod
    use types, only: dp
    use constants, only : i_
    implicit none

    public :: gradient

    interface gradient
        module procedure gradient_fields
        module procedure gradient_density
        module procedure gradient_density_nospin
        module procedure gradient_orbitals
        module procedure gradient_orbital
    end interface


contains

    subroutine gradient_fields(field, grad_den, grids)
    use odp_type, only: field_struct
    use grids_type, only: grid_struct
    use fft, only : recip_to_real
    type(field_struct), intent(in) :: field(:)
    real(dp), intent(inout) :: grad_den(:,:,:,:,:)
    type(grid_struct), intent(inout), target :: grids(:)
    type(grid_struct),  pointer :: grid
    integer :: dir, i
    do i=1, size(field)
        grid=>grids(field(i)%grid)
        do dir=1,3
            grid%work%G=-i_*grid%G(:,:,:,dir)*field(i)%G
            call recip_to_real(grid%work, grids)
            grad_den(:,:,:,dir,i)=real(grid%work%R)
        enddo
    enddo

end subroutine
    subroutine gradient_density(den, grad_den, grids)
        use odp_type, only: spin_DenPot_struct
        use grids_type, only: grid_struct
        use fft, only : recip_to_real
        type(spin_DenPot_struct), intent(in) :: den
        real(dp), intent(inout) :: grad_den(:,:,:,:,:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct),  pointer :: grid
        integer :: dir, i
        do i=1,den%n_s
            grid=>grids(den%of(i)%grid)
            do dir=1,3
                grid%work%G=-i_*grid%G(:,:,:,dir)*den%of(i)%G
                call recip_to_real(grid%work, grids)
                grad_den(:,:,:,dir,i)=real(grid%work%R)
            enddo
        enddo

    end subroutine

    subroutine gradient_density_nospin(den, grad_den, grids)
        use odp_type, only: DenPot_struct
        use grids_type, only: grid_struct
        use fft, only : recip_to_real
        type(DenPot_struct), intent(in) :: den
        real(dp), intent(inout) :: grad_den(:,:,:,:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct),  pointer :: grid
        integer :: dir
        
        grid=>grids(den%of%grid)

        do dir=1,3
            grid%work%G=-i_*grid%G(:,:,:,dir)*den%of%G
            call recip_to_real(grid%work, grids)
            grad_den(:,:,:,dir)=real(grid%work%R)
        enddo

    end subroutine

    subroutine gradient_orbitals(psi, grad_psi, grids)
        use odp_type, only: orbital_struct
        use grids_type, only: grid_struct
        use fft, only : recip_to_real
        type(orbital_struct), intent(inout) :: psi(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct),  pointer :: grid
        complex(dp), intent(inout) :: grad_psi(:,:,:,:,:,:)
        integer :: dir, i, s
        
        do i=1,size(psi(:))
            do s=1, psi(i)%n_spinor
                grid=>grids(psi(i)%of(s)%grid)
                do dir=1,3
                    grid%work%G=-i_*(grid%G(:,:,:,dir))*psi(i)%of(s)%G
                    call recip_to_real(grid%work, grids)
                    grad_psi(:,:,:,dir,s,i)=grid%work%R
                enddo
            enddo
        enddo

    end subroutine

    subroutine gradient_orbital(psi, grad_psi, grids)
        use odp_type, only: orbital_struct
        use grids_type, only: grid_struct
        use fft, only : recip_to_real
        type(orbital_struct), intent(in) :: psi
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct),  pointer :: grid
        complex(dp), intent(inout) :: grad_psi(:,:,:,:,:)
        integer :: dir, s
       

        do s=1, psi%n_spinor
            grid=>grids(psi%of(s)%grid)
            do dir=1,3
                grid%work%G=-i_*(grid%G(:,:,:,dir))*psi%of(s)%G
                call recip_to_real(grid%work, grids)
                grad_psi(:,:,:,dir,s)=grid%work%R
            enddo
        enddo

    end subroutine
end module