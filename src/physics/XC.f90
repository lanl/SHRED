module xc_mod
    use xc_type
    use types, only: dp
    implicit none

    public

    contains

    subroutine calc_XC(density_in, xc, all_PAW, parallel, grids, xc_energy, exchange_correlation, xc_stress)
        use constants, only: i_
        use odp_type, only: spin_DenPot_struct, field_struct
        use odp, only:  allocate_field, deallocate_field
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use xc_type, only: xc_struct
        use fft, only: recip_to_real, real_to_recip 
        use operations_3D, only: integrate_3D_R
        use xc_interface, only: update_xc
        use gradient_mod, only: gradient
        use parallel_mod, only: parallel_wait
        use atom_type, only : atom_struct
        use element_type, only : element_struct
        use simulation_type, only : all_PAW_struct
        use grids_mod,only : allocate_local_fields_R

        type(parallel_struct), intent(in) :: parallel
        type(spin_DenPot_struct), intent(inout), target :: density_in
        type(all_PAW_struct), intent(inout) :: all_PAw
        type(xc_struct), intent(in) :: xc(:) !1-correlation/exchnage_correlation 2-exchange
        type(spin_DenPot_struct), optional, intent(inout) :: exchange_correlation
        real(dp), optional, intent(out) :: xc_energy, xc_stress(3,3)
        real(dp), allocatable:: exc(:,:,:), vxc(:,:,:,:), sxc(:,:,:,:), rho_total(:,:,:)! XC energy density, potential, stress tensor
        real(dp), allocatable, dimension(:,:,:,:) :: rho, sigma
        real(dp), allocatable :: grad_den(:,:,:,:,:)
        real(dp) :: tmp
        integer :: i, j, k, dir, dir2, index,usexcnhat
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        type(field_struct), allocatable, target :: density_xc_store(:)
        type(field_struct), pointer :: density_xc(:)


        grid=>grids(density_in%of(1)%grid)

        if(all_PAW%N_PAW_atoms.gt.0) then
            allocate(density_xc_store(density_in%n_s))
            do i=1,density_in%n_s
                density_xc_store(i)%grid=1
                call allocate_field(density_xc_store(i), grid, parallel)
            enddo
            density_xc(1:density_in%n_s)=>density_xc_store(1:density_in%n_s)

            usexcnhat=maxval(all_PAW%tab(1:all_PAW%N_PAW_elements)%usexcnhat)
            if(all_PAW%usexcnhat.eq.1) usexcnhat=1
            if(usexcnhat==0) then
                do i=1,density_in%n_s
                    density_xc(i)%R(:,:,:)=density_in%of(i)%R(:,:,:) &
                        -all_PAW%rho_comp%of(i)%R +all_PAW%tncore%of(i)%R(:,:,:)
                    density_xc(i)%G(:,:,:)=density_in%of(i)%G(:,:,:) &
                        -all_PAW%rho_comp%of(i)%G +all_PAW%tncore%of(i)%G(:,:,:)
                enddo
            else
                do i=1,density_in%n_s
                    density_xc(i)%R(:,:,:)=density_in%of(i)%R(:,:,:) &
                        +all_PAW%tncore%of(i)%R(:,:,:)
                    density_xc(i)%G(:,:,:)=density_in%of(i)%G(:,:,:) &
                        +all_PAW%tncore%of(i)%G(:,:,:)
                enddo
            endif
        else
            density_xc(1:density_in%n_s)=>density_in%of(1:density_in%n_s)
        endif

        call allocate_local_fields_R(rho,grid,density_in%n_s)

        do i=1,density_in%n_s
            rho(:,:,:,i)=real(density_xc(i)%R(:,:,:))
        enddo

        if((.not.present(exchange_correlation)).and.present(xc_stress))  then
            print *, 'Need to calculate exchange_correlation potential with stress'
            stop
        endif

        if(present(exchange_correlation)) call allocate_local_fields_R(vxc,grid,density_in%n_s)
        if(present(xc_energy)) call allocate_local_fields_R(exc,grid)
        if( any(xc(:)%family.eq.XC_FAMILY_GGA) ) then
            call allocate_local_fields_R(grad_den,grid, 3, density_in%n_s)
            call allocate_local_fields_R(sigma,grid,2*(density_in%n_s-1)+1)
            if(present(exchange_correlation)) then
                call allocate_local_fields_R(sxc,grid,2*(density_in%n_s-1)+1)
            endif
            call gradient(density_xc, grad_den, grids)
            sigma=0.0
            do dir=1,3
                index=0
                do j=1,density_in%n_s
                    do k=j,density_in%n_s
                        index=index+1
                        sigma(:,:,:,index)=sigma(:,:,:,index)+grad_den(:,:,:,dir,j)*grad_den(:,:,:,dir,k)
                    enddo
                enddo
            enddo


            if(allocated(exc) .and. allocated(vxc) .and. allocated(sxc)) then
                call update_xc(rho, xc, sigma, exc=exc, vxc=vxc, sxc=sxc)
            elseif(allocated(exc) .and. allocated(sxc)) then
                call update_xc(rho, xc, sigma, exc=exc, sxc=sxc)
            elseif(allocated(vxc) .and. allocated(sxc)) then
                call update_xc(rho, xc, sigma, vxc=vxc, sxc=sxc)
            elseif(allocated(exc) .and. allocated(vxc) ) then
                call update_xc(rho, xc, sigma, exc=exc,  vxc=vxc)
            elseif(allocated(exc)) then
                call update_xc(rho, xc, sigma, exc=exc)
            elseif(allocated(vxc)) then
                call update_xc(rho, xc, sigma, vxc=vxc)
            else
            endif
        else
        
            if(allocated(exc) .and. allocated(vxc) .and. allocated(sxc)) then
            call update_xc(rho, xc, exc=exc, vxc=vxc, sxc=sxc)
            elseif(allocated(exc) .and. allocated(sxc)) then
            call update_xc(rho, xc, exc=exc, sxc=sxc)
            elseif(allocated(vxc) .and. allocated(sxc)) then
            call update_xc(rho, xc, vxc=vxc, sxc=sxc)
            elseif(allocated(exc) .and. allocated(vxc) ) then
            call update_xc(rho, xc, exc=exc,  vxc=vxc)
            elseif(allocated(exc)) then
            call update_xc(rho, xc, exc=exc)
            elseif(allocated(vxc)) then
            call update_xc(rho, xc, vxc=vxc)
            else
            endif

        endif

        if(present(xc_stress)) then
            xc_stress=0.0_dp
        endif

        if( any(xc(:)%family.eq.XC_FAMILY_GGA) ) then
            if(present(xc_stress)) then
                call allocate_local_fields_R(rho_total,grid)
                rho_total=0.0_dp
                do j=1,density_in%n_s
                    rho_total=rho_total+rho(:,:,:,j)
                enddo

                do dir=1,3; do dir2=dir,3
                    index=0
                    do j=1,density_in%n_s
                        do k=j,density_in%n_s
                            index=index+1
                          ! Abinit does not seem to cutoff here, based on comparison of pressure, so I'm following
                            ! but it seems like there could be Wrap around error
                          !  grid%work%R=sxc(:,:,:,index)
                          !  call real_to_recip(grid%work, grids)
                          !  grid%work%G=grid%work%G*grid%cutden
                          !  call recip_to_real(grid%work, grids)
                          !  sxc(:,:,:,index)=real(grid%work%R) 
                            xc_stress(dir,dir2)=xc_stress(dir,dir2)-&
                                real(integrate_3D_R( &
                                2.0_dp*sxc(:,:,:,index)*grad_den(:,:,:,dir2,j)*grad_den(:,:,:,dir,k)&
                                ,grid, parallel))
                        enddo
                    enddo
                    xc_stress(dir2,dir)= xc_stress(dir,dir2)
                enddo; enddo
            endif
            if(allocated(vxc)) then
                do dir=1,3
                    index=0
                    do j=1,density_in%n_s
                        do k=j,density_in%n_s
                            index=index+1
                            if(k.eq.j) then
                                grid%work%R=2.0_dp*grad_den(:,:,:,dir,j)*sxc(:,:,:,index)
                                call real_to_recip(grid%work, grids)
                                grid%work%G=i_*grid%G(:,:,:,dir)*grid%work%G*grid%cutden
                                call recip_to_real(grid%work, grids)
                                vxc(:,:,:,j)=vxc(:,:,:,j)+real(grid%work%R)
                            endif
                        enddo
                    enddo
                enddo
            endif
            deallocate(grad_den)
        endif

        if(allocated(vxc)) then
            do j=1,density_in%n_s                   
                exchange_correlation%of(j)%R=vxc(:,:,:,j)
                call real_to_recip(exchange_correlation%of(j), grids)
            enddo
        endif

        if(present(xc_stress))then
            do j=1,density_in%n_s
                tmp = real(integrate_3D_R((exc(:,:,:)-vxc(:,:,:,j))*rho(:,:,:,j), grid, parallel))
                xc_stress(1,1)=xc_stress(1,1) + tmp
                xc_stress(2,2)=xc_stress(2,2) + tmp
                xc_stress(3,3)=xc_stress(3,3) +tmp
            enddo
        endif   

        if(present(xc_energy)) then
            if(.not.allocated(rho_total)) call allocate_local_fields_R(rho_total,grid)
            rho_total=0.0
            do j=1,density_in%n_s
                rho_total=rho_total+rho(:,:,:,j)
            enddo
            xc_energy=0.0_dp
            xc_energy=xc_energy + real(integrate_3D_R(exc(:,:,:)*rho_total, grid, parallel))
        endif

        if(all_PAW%N_PAW_atoms.gt.0) then
            do i=1,density_in%n_s
                call deallocate_field(density_xc_store(i),grids(density_xc_store(i)%grid))
            enddo
            deallocate(density_xc_store)
        endif
        if(allocated(vxc)) deallocate(vxc)
        if(allocated(exc)) deallocate(exc)
        if(allocated(sxc)) deallocate(sxc)

        return
    end subroutine

end module