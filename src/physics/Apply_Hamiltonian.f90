module Apply_Hamiltonian
    use types, only : dp

    implicit none

    public :: Apply_H, Apply_S,  Apply_H_and_S, Apply_SinvH, Apply_HSinv


    contains

    subroutine Apply_H(field, Hfield, grids, potentials, atoms, elements, parallel, calc_R, PYpsi_in)
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct

        use odp_type, only : field_struct
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: real_to_recip, recip_to_real
        use Non_Local_ion, only : Apply_projectors_RL, Apply_projectors_RL_gamma, calculate_projector_overlaps, &
                                  Apply_projectors_Recip
        type(parallel_struct), intent(in) :: parallel
        type(field_struct), intent(inout) :: field(:)
        type(field_struct), intent(inout) :: Hfield(:)
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        complex(dp), intent(in), optional, target :: PYPsi_in(:,:,:)
        complex(dp), pointer :: PYPsi(:,:,:)

        logical, intent(in) :: calc_R
        integer :: s, i
        real(dp) :: time1, time2

        !if(any(elements(:)%recip_apply)) then
        !    call Apply_H_RecipNL(field, Hfield, grids, potentials, atoms, elements, parallel, calc_R=calc_R)
        !    return
        !endif

        if(calc_R) then
            do s=1,size(field)
                call recip_to_real(field(s), grids)
            enddo
        endif

        grid=>grids(field(1)%grid)
        do s=1, size(field)
            Hfield(s)%G=0.0_dp
            Hfield(s)%R=0.0_dp
        enddo

    !
     !  do i=1, size(grid%work%G,3)
    !        if(grid%gamma.and.abs(grid%G2(1,1,i)).gt.tiny(1.0_dp)) &
    !        print *, 'inG3:', i, field(1)%G(1,1,i)/sqrt(2.0), grid%G2(1,1,i)
    !        if(.not.grid%gamma.or.abs(grid%G2(1,1,i)).lt.tiny(1.0_dp))&
    !          print *, 'inG3:', i, field(1)%G(1,1,i), grid%G2(1,1,i)
    !     enddo
    !   do i=1, size(grid%work%G,2)
    !         print *, 'inG2:', i, field(1)%G(1,i,1), grid%G2(1,i,1)
    !    enddo
    !    do i=1, size(grid%work%G,1)
    !         print *, 'inG1:', i, field(1)%G(i,1,1), grid%G2(i,1,1)
    !    enddo
    !   do i=1, size(grid%work%G,3)
    !    if(grid%gamma.and.abs(grid%G2(i,i,i)).gt.tiny(1.0_dp))&
    !      print *, 'inGD:', i, field(1)%G(i,i,i)/sqrt(2.0), grid%G2(i,i,i) 
    !    if(.not.grid%gamma.or.abs(grid%G2(i,i,i)).lt.tiny(1.0_dp))&
     !     print *, 'inGD:', i, field(1)%G(1,1,i), grid%G2(1,1,i)
     ! enddo


        if(present(PYpsi_in)) then
            PYpsi=>PYpsi_in
        else
            allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(field)))
            PYpsi=0.0_dp
            do s=1, size(field)
                call Calculate_Projector_overlaps(field(s),  &
                PYpsi(:,:,s), atoms, elements, grid, parallel) 
             enddo
        endif
        call cpu_time(time1)
        !Apply any real-space projectors
        if(grid%gamma) then
            call Apply_projectors_RL_gamma(field, atoms, elements, grid, PYpsi, Vpsi=Hfield) 
        else
            call Apply_projectors_RL(field, atoms, elements, grid, PYpsi, Vpsi=Hfield) 
        endif
        !Apply any recip-space projectors
            call Apply_projectors_Recip(field, atoms, elements, grid, PYpsi, Vpsi=Hfield) 
        call cpu_time(time2)
        !print *, 'Rl time:', time2-time1
        call cpu_time(time1)
        do s=1, size(field)

            Hfield(s)%R=  Hfield(s)%R + field(s)%R*real(potentials%total_local%of(s)%R)
            grid%work%R =  Hfield(s)%R
            
            call real_to_recip(grid%work, grids)
           ! do i=1, size(grid%work%G,3)
           !     print *, 'midGD:', i,  Hfield(s)%G(1,1,i), &
          !           grid%work%G(1,1,i), grid%G2(1,1,i)*0.5_dp*field(s)%G(1,1,i), grid%G2(1,1,i)
           ! enddo
            Hfield(s)%G = Hfield(s)%G + grid%work%G
            
            if(grid%ecutsm.gt.0.0_dp) then
                Hfield(s)%G = Hfield(s)%G + grid%G2*0.5_dp*field(s)%G*grid%p
            else
                Hfield(s)%G = Hfield(s)%G + grid%G2*0.5_dp*field(s)%G

            endif
           
            Hfield(s)%G = Hfield(s)%G*grid%cutwf
      !      do i=1, size(grid%work%G,3)
     !           if(grid%gamma.and.abs(grid%G2(1,1,i)).gt.tiny(1.0_dp)) &
     !           print *, 'G3:', i, Hfield(s)%G(1,1,i)/sqrt(2.0), grid%G2(1,1,i)
     !           if(.not.grid%gamma.or.abs(grid%G2(1,1,i)).lt.tiny(1.0_dp)) &
     !           print *, 'G3:', i, Hfield(s)%G(1,1,i), grid%G2(1,1,i)
     !       enddo
     !       do i=1, size(grid%work%G,2)
     !            print *, 'G2:', i, Hfield(s)%G(1,i,1), grid%G2(1,i,1)
     !       enddo
     !       do i=1, size(grid%work%G,1)
     !            print *, 'G1:', i, Hfield(s)%G(i,1,1), grid%G2(i,1,1)
     !       enddo
     !       do i=1, size(grid%work%G,3)
     !           if(grid%gamma.and.abs(grid%G2(i,i,i)).gt.tiny(1.0_dp))  &
     !           print *, 'GD:', i, Hfield(s)%G(i,i,i)/sqrt(2.0_dp), grid%G2(i,i,i)
     !           if(.not.grid%gamma.or.abs(grid%G2(i,i,i)).lt.tiny(1.0_dp)) &
      !           print *, 'GD:', i, Hfield(s)%G(1,1,i), grid%G2(1,1,i)
      !     enddo
        enddo
        
        if(.not.present(PYpsi_in)) then
            deallocate(PYpsi)
        else
            nullify(PYpsi)
        endif

        call cpu_time(time2)
      !  print *, 'Loc time:', time2-time1
    end subroutine

    subroutine Apply_H_RecipNL(field, Hfield, grids, potentials, atoms, elements, parallel, calc_R)
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct

        use odp_type, only : field_struct
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: real_to_recip, recip_to_real
        use Non_Local_ion, only : Apply_projectors_Recip, Calculate_Projector_overlaps_Recip
        type(parallel_struct), intent(in) :: parallel
        type(field_struct), intent(inout) :: field(:)
        type(field_struct), intent(inout) :: Hfield(:)
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        complex(dp), pointer :: PYPsi(:,:,:)

        logical, intent(in) :: calc_R
        integer :: s
        real(dp) :: time1, time2

        if(calc_R) then
            do s=1,size(field)
                call recip_to_real(field(s), grids)
            enddo
        endif

        grid=>grids(field(1)%grid)
        do s=1, size(field)
            Hfield(s)%G=0.0_dp
            Hfield(s)%R=0.0_dp
        enddo

        allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(field)))
        PYpsi=0.0_dp

        do s=1, size(field)
            call Calculate_Projector_overlaps_Recip(field(s)%G, field(s)%grid,  &
            PYpsi(:,:,s), atoms, elements, grid, parallel) 
        enddo

        call cpu_time(time1)
        call Apply_projectors_Recip(field, atoms, elements, grid, PYpsi, Vpsi=Hfield) 

        call cpu_time(time2)
        !print *, 'Rl time:', time2-time1
        call cpu_time(time1)
        do s=1, size(field)

            grid%work%R = field(s)%R*potentials%total_local%of(s)%R
            call real_to_recip(grid%work, grids)
            Hfield(s)%G = grid%work%G + Hfield(s)%G 
            if(grid%ecutsm.gt.0.0_dp) then
                Hfield(s)%G = Hfield(s)%G + grid%G2*0.5_dp*field(s)%G*grid%p
            else
                Hfield(s)%G = Hfield(s)%G + grid%G2*0.5_dp*field(s)%G

            endif
            Hfield(s)%G = Hfield(s)%G*grid%cutwf
        enddo
        
        deallocate(PYpsi)

        call cpu_time(time2)
      !  print *, 'Loc time:', time2-time1
    end subroutine

    subroutine Apply_SinvH(field, SIHfield, grids, potentials, atoms, elements, parallel, all_PAW, calc_G, PYpsi_in, with_CG)
        !Unlike other Apply_H this ends in real space no reciprocal, so it asks if you need G at input
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct

        use odp_type, only : field_struct
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: real_to_recip, recip_to_real
        use Non_Local_ion, only : Apply_projectors_RL, Apply_projectors_RL_gamma, &
                                  calculate_projector_overlaps, Apply_S_power,Apply_S_inverse, Apply_projectors_Recip
        use simulation_type, only : all_PAW_struct
        use odp, only : allocate_field, deallocate_field

        type(parallel_struct), intent(in) :: parallel
        type(field_struct), intent(inout) :: field(:)
        type(field_struct), intent(inout) :: SIHfield(:)
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        type(all_PAW_struct), intent(inout) :: all_PAW
        logical, intent(in), optional :: with_CG
        complex(dp), intent(in), optional, target :: PYPsi_in(:,:,:)
        logical, intent(in) :: calc_G
        type(field_struct), allocatable:: Hfield(:)

        complex(dp), pointer :: PYPsi(:,:,:)
        logical :: with_CG_

        integer :: s
        real(dp) :: time1, time2

        with_CG_=.false.
        if(present(with_CG)) with_CG_=with_CG

        if(calc_G) then
            do s=1,size(field)
                call real_to_recip(field(s), grids)
            enddo
        endif
        
        grid=>grids(field(1)%grid)

        allocate(Hfield(size(field)))
        do s=1, size(field)
            call allocate_field(Hfield(s), grid, parallel)
        enddo
        Hfield%grid=field%grid

        do s=1, size(field)
            Hfield(s)%G=0.0_dp
            Hfield(s)%R=0.0_dp
        enddo

        if(present(PYpsi_in)) then
            PYpsi=>PYpsi_in
        else
            allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(field)))
            PYpsi=0.0_dp
            do s=1, size(field)
                call Calculate_Projector_overlaps(field(s),  &
                PYpsi(:,:,s), atoms, elements, grid, parallel) 
             enddo
        endif

        call cpu_time(time1)
        if(grid%gamma) then
            call Apply_projectors_RL_gamma(field, atoms, elements, grid, PYpsi, Vpsi=Hfield) 
        else
            call Apply_projectors_RL(field, atoms, elements, grid, PYpsi, Vpsi=Hfield) 
        endif
        call Apply_projectors_Recip(field, atoms, elements, grid, PYpsi, Vpsi=Hfield) 

        call cpu_time(time2)
        !print *, 'Rl time:', time2-time1
        call cpu_time(time1)
        do s=1, size(field)
            Hfield(s)%R= Hfield(s)%R + field(s)%R*potentials%total_local%of(s)%R
            grid%work%R =  Hfield(s)%R

            call real_to_recip(grid%work, grids)
            Hfield(s)%G = Hfield(s)%G + grid%work%G

            if(grid%ecutsm.gt.0.0_dp) then
                Hfield(s)%G= grid%G2*0.5_dp*field(s)%G*grid%p + Hfield(s)%G
            else
                Hfield(s)%G = grid%G2*0.5_dp*field(s)%G + Hfield(s)%G
            endif

            Hfield(s)%G =Hfield(s)%G*grid%cutwf
            call recip_to_real(Hfield(s), grids)
        enddo

        !call Apply_S_power(Hfield(:), atoms, elements, grid, parallel, SIHfield(:), -1.0_dp) 
        call Apply_S_inverse(Hfield(:), atoms, elements, grids, parallel, SIHfield(:), all_PAW, with_CG) 
       
        do s=1, size(field)
            call deallocate_field(Hfield(s),grids(Hfield(s)%grid))
        enddo
        deallocate(Hfield)
        call cpu_time(time2)

        if(.not.present(PYpsi_in)) then
            deallocate(PYpsi)
        else
            nullify(PYpsi)
        endif

        !print *, 'Loc time:', time2-time1
    end subroutine

    subroutine Apply_HSinv(field, HSIfield, grids, potentials, atoms, elements, parallel, all_PAW, calc_R, with_CG)
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct

        use odp_type, only : field_struct
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: real_to_recip, recip_to_real
        use Non_Local_ion, only : Apply_projectors_RL,Apply_projectors_RL_gamma, Apply_projectors_Recip, &
                                  calculate_projector_overlaps, Apply_S_power,Apply_S_inverse
        use simulation_type, only : all_PAW_struct
        use odp, only : allocate_field, deallocate_field

        type(parallel_struct), intent(in) :: parallel
        type(field_struct), intent(inout) :: field(:)
        type(field_struct), intent(inout) :: HSIfield(:)
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        type(all_PAW_struct), intent(inout) :: all_PAW
        logical, intent(in), optional :: with_CG
        logical, intent(in) :: calc_R
        type(field_struct), allocatable:: Sfield(:)

        complex(dp), pointer :: PYPsi(:,:,:)
        logical :: with_CG_

        integer :: s
        real(dp) :: time1, time2

        with_CG_=.false.
        if(present(with_CG)) with_CG_=with_CG

        if(calc_R) then
            do s=1,size(field)
                call recip_to_real(field(s), grids)
            enddo
        endif
        
        grid=>grids(field(1)%grid)

        allocate(Sfield(size(field)))
        do s=1, size(field)
            call allocate_field(Sfield(s), grid, parallel)
        enddo
        Sfield%grid=field%grid

        do s=1, size(field)
            HSIfield(s)%G=0.0_dp
            HSIfield(s)%R=0.0_dp
        enddo
        call Apply_S_inverse(field(:), atoms, elements, grids, parallel, Sfield(:), all_PAW, with_CG) 

        allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(field)))
        PYpsi=0.0_dp
        do s=1, size(field)
            call Calculate_Projector_overlaps(Sfield(s),  &
            PYpsi(:,:,s), atoms, elements, grid, parallel) 
        enddo

        call cpu_time(time1)
        if(grid%gamma) then
            call Apply_projectors_RL_gamma(Sfield, atoms, elements, grid, PYpsi, Vpsi=HSIfield) 
        else
            call Apply_projectors_RL(Sfield, atoms, elements, grid, PYpsi, Vpsi=HSIfield) 
        endif
            call Apply_projectors_Recip(Sfield, atoms, elements, grid, PYpsi, Vpsi=HSIfield) 

        call cpu_time(time2)
        !print *, 'Rl time:', time2-time1
        call cpu_time(time1)
        do s=1, size(field)
            HSIfield(s)%R= HSIfield(s)%R + Sfield(s)%R*potentials%total_local%of(s)%R
            grid%work%R=HSIfield(s)%R
            call real_to_recip(grid%work, grids)

            HSIfield(s)%G = HSIfield(s)%G + grid%work%G

            if(grid%ecutsm.gt.0.0_dp) then
                HSIfield(s)%G= grid%G2*0.5_dp*Sfield(s)%G*grid%p + HSIfield(s)%G
            else
                HSIfield(s)%G = grid%G2*0.5_dp*Sfield(s)%G + HSIfield(s)%G
            endif

            HSIfield(s)%G =HSIfield(s)%G*grid%cutwf
            call recip_to_real(HSIfield(s), grids)
        enddo

        do s=1, size(field)
            call deallocate_field(Sfield(s),grids(Sfield(s)%grid))
        enddo
        deallocate(Sfield)
        call cpu_time(time2)

        deallocate(PYpsi)
    
        !print *, 'Loc time:', time2-time1
    end subroutine

    subroutine Apply_S(field, Sfield, grids, atoms, elements, parallel, calc_R, PYpsi_in)
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct

        use odp_type, only : field_struct
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: real_to_recip, recip_to_real
        use Non_Local_ion, only : Apply_projectors_RL, Apply_projectors_RL_gamma, calculate_projector_overlaps,&
                                  Apply_projectors_Recip
        type(parallel_struct), intent(in) :: parallel
        type(field_struct), intent(inout) :: field(:)
        type(field_struct), intent(inout) :: Sfield(:)
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        logical, intent(in) :: calc_R
        complex(dp), intent(in), optional, target :: PYPsi_in(:,:,:)
        complex(dp), pointer :: PYPsi(:,:,:)

        integer :: s
        real(dp) :: time1, time2

        grid=>grids(field(1)%grid)

        if(calc_R) then
            do s=1,size(field)
                call recip_to_real(field(s), grids)
            enddo
        endif

        call cpu_time(time1)

        if(present(PYpsi_in)) then
            PYpsi=>PYpsi_in
        else
            allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(field)))
            PYpsi=0.0_dp
            do s=1, size(field)
                call Calculate_Projector_overlaps(field(s),  &
                PYpsi(:,:,s), atoms, elements, grid, parallel) 
             enddo
        endif
        
        if(grid%gamma) then
            call Apply_projectors_RL_gamma(field, atoms, elements, grid, PYpsi, Spsi=Sfield) 
        else
            call Apply_projectors_RL(field, atoms, elements, grid, PYpsi, Spsi=Sfield) 
        endif
            call Apply_projectors_Recip(field, atoms, elements, grid, PYpsi, Spsi=Sfield) 


        call cpu_time(time2)
        !print *, 'Rl time:', time2-time1
        call cpu_time(time1)
        do s=1, size(field)
            grid%work%R =  Sfield(s)%R
            call real_to_recip(grid%work, grids)
            Sfield(s)%G=Sfield(s)%G + grid%work%G
            Sfield(s)%G=Sfield(s)%G*grid%cutwf
            call recip_to_real(Sfield(s), grids)
        enddo
        call cpu_time(time2)
        !print *, 'Loc time:', time2-time1
        if(.not.present(PYpsi_in)) then
            deallocate(PYpsi)
        else
            nullify(PYpsi)
        endif
    end subroutine

    subroutine Apply_H_and_S(field, Hfield, Sfield, grids, potentials, atoms, elements, parallel, calc_R, PYpsi_in)
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct

        use odp_type, only : field_struct
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: real_to_recip, recip_to_real
        use Non_Local_ion, only : Apply_projectors_RL, Apply_projectors_RL_gamma, calculate_projector_overlaps, &
                                  Apply_projectors_Recip
        type(parallel_struct), intent(in) :: parallel
        type(field_struct), intent(inout) :: field(:)
        type(field_struct), intent(inout) :: Hfield(:), Sfield(:)
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        logical, intent(in) :: calc_R
        complex(dp), intent(in), optional, target :: PYPsi_in(:,:,:)
        complex(dp), pointer :: PYPsi(:,:,:)

        integer :: s
        real(dp) :: time1, time2

        if(calc_R) then
            do s=1,size(field)
                call recip_to_real(field(s), grids)
            enddo
        endif

        grid=>grids(field(1)%grid)
        do s=1, size(field)
            Hfield(s)%G=0.0_dp
            Hfield(s)%R=0.0_dp
            Sfield(s)%R=0.0_dp
            Sfield(s)%G=0.0_dp
        enddo

        call cpu_time(time1)

        if(present(PYpsi_in)) then
            PYpsi=>PYpsi_in
        else
            allocate(PYpsi(maxval(elements(:)%n_proj), size(atoms), size(field)))
            PYpsi=0.0_dp
            do s=1, size(field)
                call Calculate_Projector_overlaps(field(s),  &
                PYpsi(:,:,s), atoms, elements, grid, parallel) 
             enddo
        endif
       
        if(grid%gamma) then
            call Apply_projectors_RL_gamma(field, atoms, elements, grid, PYpsi, Vpsi=Hfield, Spsi=Sfield) 
        else
            call Apply_projectors_RL(field, atoms, elements, grid, PYpsi, Vpsi=Hfield, Spsi=Sfield) 
        endif
            call Apply_projectors_Recip(field, atoms, elements, grid, PYpsi, Vpsi=Hfield, Spsi=Sfield) 


        call cpu_time(time2)
        !print *, 'Rl time:', time2-time1
        call cpu_time(time1)
        do s=1, size(field)
            grid%work%R =  Sfield(s)%R
            call real_to_recip(grid%work, grids)
            Sfield(s)%G=(Sfield(s)%G+grid%work%G)*grid%cutwf

            Hfield(s)%R=  Hfield(s)%R + field(s)%R*potentials%total_local%of(s)%R
            grid%work%R =  Hfield(s)%R
            call real_to_recip(grid%work, grids)
            Hfield(s)%G = (Hfield(s)%G + grid%work%G)*grid%cutwf
            Hfield(s)%G = Hfield(s)%G + grid%G2*0.5_dp*field(s)%G
        enddo
        call cpu_time(time2)
        !print *, 'Loc time:', time2-time1
        if(.not.present(PYpsi_in)) then
            deallocate(PYpsi)
        else
            nullify(PYpsi)
        endif
    end subroutine

end module