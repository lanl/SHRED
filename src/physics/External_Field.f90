module External_Field
    use types, only : op, dp
    use constants, only : pi
    implicit none
    public :: External_Field_Potential, External_Field_Energy, &
    Initialize_External_Field_Potential,  Initialize_External_Field_Potential_restart

    contains

    !Note that this potential and Energy only accounts for the Length Gauge portion of an electric field, Vector potential
    !contribution will be combined with the Kinetic Energy
    subroutine External_Field_Potential(td, potentials, grids)
        use td_type, only : td_struct
        use simulation_type, only : potentials_struct
        use grids_type, only : grid_struct, inter_grid_struct

        type(td_struct), intent(inout) :: td
        type(potentials_struct), intent(inout) :: potentials 
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid

        real(dp) :: A_last(3)

        integer :: g,i

        !Update the Field
        if(td%time.lt.td%t_max_field) then
            A_last(:)=td%A(:)
            if(td%w_field.gt.0) then
                td%A(:)= td%A(:)-td%Ex(:)*(td%time-td%last_time)*0.5_dp
                td%Ex(:) = -td%E0(:) * sin(pi*(td%time-td%t0)/td%tw)**2*cos(td%w_field*td%time+td%phase_field)
                td%A(:) = td%A(:)-td%Ex(:)*(td%time-td%last_time)*0.5_dp
            else
                td%Ex(:) = -td%E0(:) * exp(-(td%time-td%t0)**2/(2*td%tw**2)) / (sqrt(2*pi)*td%tw)
                td%A(:) = td%E0(:) * (erf(sqrt(2._dp)*(td%time - td%t0)/(2*td%tw))/2 + 1._dp/2)
            endif
            do g=3,size(grids) !grids 1 and 2 belong to density which doesnt know about Vector Potential
                grid=>grids(g)
                grid%G2(:,:,:) = 0.0_dp
                do i=1,3
                    grid%G(:,:,:,i) = grid%G(:,:,:,i) -A_last(i) + td%A(i) !(G+A->G+A')
                    grid%G2(:,:,:) = grid%G2(:,:,:) + abs(grid%G(:,:,:,i))**2
                enddo
                grid%A(:)=td%A(:)
            enddo
            do i=1,3
                !Pure Velocity Gauge right now
                potentials%external_field%of%R=0.0_dp
                potentials%external_field%of%G=0.0_dp
            enddo
        endif
    end subroutine

    subroutine Initialize_External_Field_Potential(td, potentials, grids)
        use td_type, only : td_struct
        use simulation_type, only : potentials_struct
        use grids_type, only : grid_struct, inter_grid_struct

        type(td_struct), intent(inout) :: td
        type(potentials_struct), intent(inout) :: potentials 
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid

        integer :: g,i, it, nt_int=10000
        real(dp) :: time_int, dt_int
        

        !Initialize the Field and Vector Potential
        if(td%w_field.gt.0.0_dp) then
            td%A(:)=0.0
            dt_int=(td%time-100.0_dp*td%t0)/nt_int
            time_int=-100.0_dp*td%t0
            do it=1,nt_int
                time_int=time_int+dt_int
                td%Ex(:) = -td%E0(:) * sin(pi*(time_int-td%t0)/td%tw)**2*cos(td%w_field*time_int+td%phase_field)
                td%A(:)=td%A(:)-td%Ex(:)*dt_int
            enddo            
        else
            if((td%time-td%t0)**2/(2*td%tw**2).lt.50.0_dp) then
                td%Ex(:) = -td%E0(:) * exp(-(td%time-td%t0)**2/(2*td%tw**2)) / (sqrt(2*pi)*td%tw)
            else
                td%Ex(:) = 0.0_dp
            endif
            td%A(:) = td%E0(:) * (erf(sqrt(2._dp)*(td%time - td%t0)/(2*td%tw))/2 + 1._dp/2)
        endif
                
        do g=3,size(grids) !grids 1 and 2 belong to density which doesnt know about Vector Potential
            grid=>grids(g)
            grid%G2(:,:,:) = 0.0_dp
            do i=1,3
                grid%G(:,:,:,i) = grid%G(:,:,:,i) + td%A(i) !(G->G+A(t))
                grid%G2(:,:,:) = grid%G2(:,:,:) + abs(grid%G(:,:,:,i))**2
            enddo
            grid%A(:)=td%A(:)
        enddo
        
        do i=1,3
            !Pure Velocity Gauge right now
            potentials%external_field%of%R=0.0_dp
            potentials%external_field%of%G=0.0_dp
        enddo
    end subroutine

    subroutine Initialize_External_Field_Potential_restart(td, potentials, grids)
        use td_type, only : td_struct
        use simulation_type, only : potentials_struct
        use grids_type, only : grid_struct, inter_grid_struct

        type(td_struct), intent(inout) :: td
        type(potentials_struct), intent(inout) :: potentials 
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid

        integer :: g,i

        !Initialize the Field and Vector Potential
        do g=3,size(grids) !grids 1 and 2 belong to density which doesnt know about Vector Potential
            grid=>grids(g)
            grid%G2(:,:,:) = 0.0_dp
            do i=1,3
                grid%G(:,:,:,i) = grid%G(:,:,:,i) + td%A(i) !(G->G+A(t))
                grid%G2(:,:,:) = grid%G2(:,:,:) + abs(grid%G(:,:,:,i))**2
            enddo
            grid%A(:)=td%A(:)
        enddo
        
        do i=1,3
            !Pure Velocity Gauge right now
            potentials%external_field%of%R=0.0_dp
            potentials%external_field%of%G=0.0_dp
        enddo
    end subroutine

    function External_Field_Energy(external_potential, density, grids, parallel) result (Een_local)
        use odp_type, only: spin_DenPot_struct, DenPot_struct
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use operations_3D, only: integrate_3D_R
        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(in) :: grids(:)
        type(spin_DenPot_struct), intent(in) :: density
        type(DenPot_struct), intent(in) :: external_potential
        real(dp) :: Een_local
        integer :: i

        if(any(density%of(:)%grid.ne.external_potential%of%grid)) then
            print *, 'Attempting to calculate Local Ion energy from density and potential on different grid, stopping'
            stop
        endif

        Een_local= 0.0_dp
        do i=1, density%n_s
            Een_local= Een_local + real(integrate_3D_R( &
                real(density%of(i)%R)*real(external_potential%of%R), grids(external_potential%of%grid), parallel))
        enddo
        
    end function

end module