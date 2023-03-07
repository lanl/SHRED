module Hartree
    use types, only : op, dp
    use constants, only : pi
    implicit none
    public :: calc_Hartree_potential, Hartree_Energy, Hartree_Stress

    contains

    subroutine calc_Hartree_potential(hartree, density, grids)
        ! Calculates VeeG = 4*pi*neG / G2, but skips G2=0 point
        use odp_type, only: spin_DenPot_struct, DenPot_struct
        use grids_type, only: grid_struct
        use fft, only: recip_to_real 

        implicit none

        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(spin_DenPot_struct), intent(in) :: density
        type(DenPot_struct), intent(inout) :: hartree

        integer :: i 
        
        if(any(density%of(:)%grid.ne.hartree%of%grid)) then
            print *, 'Attempting to calculate Hartree from density on different grid, stopping'
            stop
        endif
    
        grid=>grids(hartree%of%grid)
        hartree%of%G= 0.0_dp
        do i=1, density%n_s
            where(grid%G2>tiny(1.0_dp))
                hartree%of%G = hartree%of%G + 4.0_dp*pi*density%of(i)%G/grid%G2
            elsewhere
                hartree%of%G = 0.0_dp
            end where 
        enddo   
        hartree%of%G = hartree%of%G * grid%cutden    
        call recip_to_real(hartree%of, grids)
            
       
    end subroutine

    function Hartree_Energy(hartree, density, parallel, grids) result (Eee)
        use odp_type, only: DenPot_struct, spin_DenPot_struct
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use operations_3D, only: integrate_3D_R
        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(spin_DenPot_struct), intent(in) :: density
        type(DenPot_struct), intent(in) :: hartree
        real(dp):: Eee
        integer :: i 

        if(any(density%of(:)%grid.ne.hartree%of%grid)) then
            print *, 'Attempting to calculate Hartree from density on different grid, stopping'
            stop
        endif
        
        Eee=0.0_dp
        do i=1, density%n_s
            grid=>grids(density%of(i)%grid)
            Eee = Eee + real(integrate_3D_R(real(density%of(i)%R)*real(hartree%of%R), grid, parallel))/2.0_dp
        enddo
        
    end function

    function Hartree_Stress(hartree, density, parallel, grids) result (See)
        use odp_type, only: DenPot_struct, spin_DenPot_struct
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use operations_3D, only: integrate_3D_G
        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(spin_DenPot_struct), intent(in) :: density
        type(DenPot_struct), intent(in) :: hartree
        real(dp):: See(3,3), delta(3,3)
        logical :: check=.false.
        integer :: i, dir, dir2

        if(any(density%of(:)%grid.ne.hartree%of%grid)) then
            print *, 'Attempting to calculate Hartree from density on different grid, stopping'
            stop
        endif
        
        See=0.0_dp
        do i=1, density%n_s
            grid=>grids(density%of(i)%grid)
            if(grid%G2(1,1,1).lt.tiny(1.0_dp)) then
                check=.true.
                grid%G2(1,1,1)=grid%G2(1,1,1)+1.0_dp
            endif
            delta=0.0_dp
            do dir=1,3
            delta(dir,dir)=1.0_dp
                do dir2=dir,3
                    See(dir,dir2)=See(dir,dir2)+real(integrate_3D_G( &
                conjg(density%of(i)%G)*hartree%of%G*(2.0_dp*grid%G(:,:,:,dir)*grid%G(:,:,:,dir2)/grid%G2(:,:,:)-delta(dir,dir2)) &
                    , grid, parallel))/2.0_dp
                enddo
            enddo
            if(check) grid%G2(1,1,1)=grid%G2(1,1,1)-1.0_dp
            check=.false.
        enddo
        do dir=1,3
            do dir2=dir,3
            See(dir2,dir)=See(dir,dir2)
            !print *,  'Hartree Stress:', dir, dir2, See(dir,dir2)
            enddo
        enddo
        
    end function

end module