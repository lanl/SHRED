module Kinetic
    use types, only : op, dp
    use constants, only : pi
    implicit none
    public :: Kinetic_Energy, Kinetic_Stress

    contains

    function Kinetic_Energy(psi, parallel, grids) result (Ekin)
        use odp_type, only: orbital_struct
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use operations_3D, only: integrate_3D_G
        use parallel_mod, only: parallel_task
        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(orbital_struct), intent(in) :: psi(:,:,:)

        real(dp):: Ekin
        integer :: i,j,k,s
        
        Ekin=0.0_dp
        !Note that i and s both relate account for spin,
        !but only one of the two should be able to be > 1 at a time
        do i=1, size(psi,3); do j=1, size(psi,2); do k=1, size(psi,1); do s=1,psi(k,j,i)%n_spinor
                if(abs(psi(k,j,i)%occ(s)).lt.tiny(1.0_dp)) cycle   
                grid=>grids(psi(k,j,i)%of(s)%grid)
                if(grid%ecutsm.gt.0.0_dp) then
                    Ekin = Ekin + psi(k,j,i)%occ(s)*psi(k,j,i)%weight(s)* &
                    0.5_dp*integrate_3D_G(abs(psi(k,j,i)%of(s)%G)**2*grid%G2*grid%p, grid, parallel)
                else
                    Ekin = Ekin + psi(k,j,i)%occ(s)*psi(k,j,i)%weight(s)* &
                    0.5_dp*integrate_3D_G(abs(psi(k,j,i)%of(s)%G)**2*grid%G2, grid, parallel)
                endif
        end do;end do;end do;end do

        call parallel_task('sum', Ekin, parallel, 'space')

        
    end function

    function Kinetic_Stress(psi, parallel, grids) result (Skin)
        use odp_type, only: orbital_struct
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use operations_3D, only: integrate_3D_G
        use parallel_mod, only: parallel_task
        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(orbital_struct), intent(in) :: psi(:,:,:)

        real(dp):: Skin(3,3)
        integer :: i,j,k,s, dir, dir2
        
        Skin=0.0_dp
        !Note that i and s both relate account for spin,
        !but only one of the two should be able to be > 1 at a time
        do i=1, size(psi,3); do j=1, size(psi,2); do k=1, size(psi,1); do s=1,psi(k,j,i)%n_spinor
                if(abs(psi(k,j,i)%occ(s)).lt.tiny(1.0_dp)) cycle   
                grid=>grids(psi(k,j,i)%of(s)%grid)
                if(grid%ecutsm.gt.0.0_dp) then
                    do dir = 1,3; do dir2=dir,3
                        Skin(dir,dir2) = Skin(dir,dir2) - psi(k,j,i)%occ(s)*psi(k,j,i)%weight(s)* &
                        integrate_3D_G(abs(psi(k,j,i)%of(s)%G)**2*grid%G(:,:,:,dir)*grid%G(:,:,:,dir2)*&
                        (grid%p - 0.5_dp*grid%G2(:,:,:)*grid%pprime/grid%ecutsm), grid, parallel)
                    enddo; enddo
                else
                    do dir = 1,3; do dir2=dir,3
                        Skin(dir,dir2) = Skin(dir,dir2) - psi(k,j,i)%occ(s)*psi(k,j,i)%weight(s)* &
                            integrate_3D_G(abs(psi(k,j,i)%of(s)%G)**2*grid%G(:,:,:,dir)*grid%G(:,:,:,dir2), grid, parallel)
                    enddo; enddo
                endif
        end do;end do;end do;end do
        
        do dir = 1,3; do dir2=dir,3
            Skin(dir2,dir) = Skin(dir,dir2)
        enddo; enddo

        call parallel_task('sum', Skin, parallel, 'space')

        
    end function

end module