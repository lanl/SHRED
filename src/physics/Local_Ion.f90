module Local_Ion
    use types, only : op, dp
    use constants, only : pi
    implicit none
    public :: calculate_local_ion_potential, Calc_Local_PP_Energy, &
    Calc_Core_Energy, calculate_local_ion_force, Local_PP_Stress, core_pot_0

    contains

    subroutine calculate_local_ion_potential(local_PP, system, atoms, elements, grids)
        use constants, only: i_
        use odp_type, only: DenPot_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use system_type, only : system_struct
        use grids_type, only: grid_struct
        use fft, only: recip_to_real 

        type(system_struct), intent(in) :: system
        type(atom_struct), intent(in) :: atoms(:)
        type(element_struct), intent(in) :: elements(:)
        type(DenPot_struct), intent(inout) :: local_PP
        integer :: e, k
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer:: grid
        grid=>grids(local_PP%of%grid)
        local_PP%of%G = 0
        do k=1,size(atoms)
            e=atoms(k)%element
            local_PP%of%G(:,:,:) = local_PP%of%G(:,:,:) + elements(e)%Ven0G(:,:,:) * &
                    exp(-i_*(grid%G(:,:,:,1)*atoms(k)%R(1) + grid%G(:,:,:,2)*atoms(k)%R(2) &
                        + grid%G(:,:,:,3)*atoms(k)%R(3)))
        enddo
        local_PP%of%G=local_PP%of%G*grid%cutden/product(grid%Box_Length(:))
        if(grid%gamma) then
            where(abs(grid%G(:,:,:,1)).gt.tiny(1.0_dp)) 
            local_PP%of%G = local_PP%of%G*sqrt(2.0_dp)
            endwhere
        endif
        call recip_to_real(local_PP%of, grids)
    end subroutine

    subroutine calculate_local_ion_force(density, system, atoms, elements, parallel, grids)
        use constants, only: i_
        use odp_type, only: spin_DenPot_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use parallel_type, only : parallel_struct
        use system_type, only : system_struct
        use grids_type, only: grid_struct
        use operations_3D, only: integrate_3D_G, integrate_3D_R
        use fft, only: recip_to_real

        type(parallel_struct), intent(in) :: parallel
        type(system_struct), intent(in) :: system
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(in) :: elements(:)
        type(spin_DenPot_struct), intent(in) :: density
        integer :: i, k, j, dir
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer:: grid
        integer :: s
        complex(dp) :: factor(size(grids(density%of(1)%grid)%G2,1), &
                              size(grids(density%of(1)%grid)%G2,2), &
                              size(grids(density%of(1)%grid)%G2,3))
        grid=>grids(density%of(1)%grid)

        k=0
        do i = 1, system%n_elements
        do j = 1, elements(i)%n_atoms_of_element
            k=k+1
            if(atoms(k)%update_me_force) then
                factor=exp(-i_*(grid%G(:,:,:,1)*atoms(k)%R(1) + grid%G(:,:,:,2)*atoms(k)%R(2) + grid%G(:,:,:,3)*atoms(k)%R(3)))
                factor=factor*grid%cutden
                if(grid%gamma) then
                    where(abs(grid%G(:,:,:,1)).gt.tiny(1.0_dp)) 
                        factor= factor*sqrt(2.0_dp)
                    endwhere
                endif
                do dir=1,3
                    atoms(k)%F_loc(dir)=0.0_dp
                    do s=1,size(density%of)
                        atoms(k)%F_loc(dir)= atoms(k)%F_loc(dir) + &
                        real(integrate_3D_G( (i_ * elements(i)%Ven0G(:,:,:) * conjg(density%of(s)%G(:,:,:)) &
                         * grid%G(:,:,:,dir) * factor(:,:,:) ), &
                        grid, parallel))/product(grid%Box_Length(:))
                    enddo
                enddo
            endif
        enddo
        enddo
            
    end subroutine

    function Calc_Core_Energy(system, elements, grids, all_PAW) result (Een_core)
        use element_type, only : element_struct, HGH, PAW
        use atom_type, only : atom_struct
        use system_type, only : system_struct    
        use grids_type, only : grid_struct    
        use simulation_type, only : all_PAW_struct    

        real(dp) :: Een_core
        type(system_struct), intent(in) :: system
        type(grid_struct), intent(in) :: grids(:)
        type(element_struct), intent(in) :: elements(:)
        type(all_PAW_struct), intent(in) :: all_PAW

        integer i,j

        Een_core=0.0_dp
        do i = 1, system%n_elements
        do j = 1, elements(i)%n_atoms_of_element
                if(elements(i)%PP_type.eq.HGH) then
                    Een_core = Een_core + elements(i)%HGH%epsatm
                else if(elements(i)%PP_type.eq.PAW.and.all_PAW%usepotzero.lt.2) then
                    Een_core = Een_core + elements(i)%PAW%epsatm 
                endif
        enddo
        enddo
        Een_core=Een_core*sum(system%nelec(:))/product(grids(1)%Box_Length(:))
        
    end function

    subroutine core_pot_0(system, elements, grids, all_PAW)
        use element_type, only : element_struct, HGH,PAW
        use atom_type, only : atom_struct
        use system_type, only : system_struct    
        use grids_type, only : grid_struct    
        use simulation_type, only : all_PAW_struct    

        type(system_struct), intent(inout) :: system
        type(grid_struct), intent(in) :: grids(:)
        type(element_struct), intent(in) :: elements(:)
        type(all_PAW_struct), intent(in) :: all_PAW

        integer i,j

        system%epsatm_all=0.0_dp
        do i = 1, system%n_elements
        do j = 1, elements(i)%n_atoms_of_element
                if(elements(i)%PP_type.eq.HGH) then
                    system%epsatm_all = system%epsatm_all + elements(i)%HGH%epsatm
                endif
        enddo
        enddo
        system%epsatm_all=system%epsatm_all
        
    end subroutine

    function Calc_Local_PP_Energy(local_PP, density, grids, parallel) result (Een_local)
        use odp_type, only: spin_DenPot_struct, DenPot_struct
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use operations_3D, only: integrate_3D_R
        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(in) :: grids(:)
        type(spin_DenPot_struct), intent(in) :: density
        type(DenPot_struct), intent(in) :: local_PP
        real(dp) :: Een_local
        integer :: i

        if(any(density%of(:)%grid.ne.local_PP%of%grid)) then
            print *, 'Attempting to calculate Local Ion energy from density and potential on different grid, stopping'
            stop
        endif

        Een_local= 0.0_dp
        do i=1, density%n_s
            Een_local= Een_local + real(integrate_3D_R( &
                real(density%of(i)%R)*real(local_PP%of%R), grids(local_PP%of%grid), parallel))
        enddo
        
    end function

    function Local_PP_Stress(density, grids, atoms, elements, parallel) result (Sen_local)
        use constants, only : i_
        use odp_type, only: spin_DenPot_struct, DenPot_struct
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use operations_3D, only: integrate_3D_G
        use element_type, only : element_struct
        use atom_type, only : atom_struct

        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(atom_struct), intent(in) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(spin_DenPot_struct), intent(in) :: density
        real(dp) :: Sen_local(3,3)
        real(dp) ::  delta(3,3)
        integer :: i, dir, dir2, e, a

        grid=>grids(density%of(1)%grid)

        Sen_local= 0.0_dp
        delta=0.0_dp
        do dir=1,3
            delta(dir,dir)=1.0_dp
            do dir2=dir,3
                    grid%work%G(:,:,:)=0.0_dp
                    do a=1,size(atoms)
                        e=atoms(a)%element
                        grid%work%G= grid%work%G + exp(-i_*(grid%G(:,:,:,1)*atoms(a)%R(1)   + &
                                                            grid%G(:,:,:,2)*atoms(a)%R(2)   + &
                                                            grid%G(:,:,:,3)*atoms(a)%R(3))) * &
        (elements(e)%dVen0GdG2(:,:,:)*grid%G(:,:,:,dir)*grid%G(:,:,:,dir2)*2.0_dp+delta(dir,dir2)*elements(e)%Ven0G)
                    enddo
                    if(grid%gamma) then
                        where(abs(grid%G(:,:,:,1)).gt.tiny(1.0_dp)) 
                            grid%work%G=grid%work%G*sqrt(2.0_dp)
                        end where
                    endif
                    do i= 1, size(density%of)
                        Sen_local(dir,dir2) =  Sen_local(dir,dir2) - real(integrate_3D_G( &
                                grid%work%G*conjg(density%of(i)%G)*grid%cutden  &
                            , grid, parallel))/product(grid%Box_Length(:))
                    enddo
            enddo        
        enddo
        do dir=1,3
            do dir2=dir,3
                Sen_local(dir2,dir)=Sen_local(dir,dir2) 
            enddo
        enddo
        
    end function

end module