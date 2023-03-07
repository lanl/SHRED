module Chebychev
    use types, only : dp
    use constants, only: i_, pi

    implicit none

    public :: cheby_exp_and_app, calculate_Chebychev_coefficients_fd, &
              calculate_Chebychev_coefficients_sqrt_fd, calculate_Chebychev_coefficients_propagator, &
              calculate_Chebychev_coefficients_entropy, calculate_Chebychev_coefficients_Gaussian, &
              calculate_Chebychev_coefficients_omfd, calculate_Chebychev_coefficients_FD_Gaussian, &
              calculate_Chebychev_coefficients_root_omfd

    interface cheby_exp_and_app
        module procedure cheby_exp_and_app_real
        module procedure cheby_exp_and_app_complex
    end interface

    contains

    subroutine cheby_exp_and_app_complex(orbitals,  &
        Ebar, deltaE, Mn, Cn, grids, potentials, atoms, elements, parallel, all_PAW)
        use odp_type, only : field_struct, orbital_struct
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use operations_3D, only:  integrate_3D_G
        use Apply_Hamiltonian, only: Apply_H
        use odp_type, only : field_struct
        use odp, only : allocate_field, deallocate_field
        use simulation_type, only : potentials_struct, all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: recip_to_real
        use parallel_mod, only: parallel_task
        complex(dp), intent(inout) :: Cn(:)

        type(parallel_struct), intent(in) :: parallel
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(grid_struct), intent(inout) :: grids(:)
        complex(dp), intent(inout) :: Mn(:)
        real(dp), intent(in) :: Ebar, deltaE
        type(all_PAW_struct), intent(inout) :: all_PAW
        
        if(all_PAW%N_PAW_atoms.lt.1) then 
            call cheby_exp_and_app_NC(orbitals,  &
            Ebar, deltaE, Mn, Cn, grids, potentials, atoms, elements, parallel)
        else
            call cheby_exp_and_app_Gen(orbitals,  &
            Ebar, deltaE, Mn, Cn, grids, potentials, atoms, elements, parallel, all_PAW)
        endif

    end subroutine

    subroutine cheby_exp_and_app_real(orbitals,  &
        Ebar, deltaE, Mn, Cn, grids, potentials, atoms, elements, parallel, all_PAW)
        use odp_type, only : field_struct, orbital_struct
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use operations_3D, only:  integrate_3D_G
        use Apply_Hamiltonian, only: Apply_H
        use odp_type, only : field_struct
        use odp, only : allocate_field, deallocate_field
        use simulation_type, only : potentials_struct, all_PAW_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: recip_to_real
        use parallel_mod, only: parallel_task
        real(dp), intent(inout) :: Cn(:)

        type(parallel_struct), intent(in) :: parallel
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(grid_struct), intent(inout) :: grids(:)
        complex(dp), intent(inout) :: Mn(:)
        real(dp), intent(in) :: Ebar, deltaE
        type(all_PAW_struct), intent(inout) :: all_PAW
        complex(dp) :: Cn_(size(Cn))

        Cn_=Cn
        if(all_PAW%N_PAW_atoms.lt.1) then 
            call cheby_exp_and_app_NC(orbitals,  &
            Ebar, deltaE, Mn, Cn_, grids, potentials, atoms, elements, parallel)
        else
            call cheby_exp_and_app_Gen(orbitals,  &
            Ebar, deltaE, Mn, Cn_, grids, potentials, atoms, elements, parallel, all_PAW)
        endif

    end subroutine


    subroutine cheby_exp_and_app_NC(orbitals,  &
            Ebar, deltaE, Mn, Cn, grids, potentials, atoms, elements, parallel)
        use odp_type, only : field_struct, orbital_struct
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use operations_3D, only:  integrate_3D_G
        use Apply_Hamiltonian, only: Apply_H
        use odp_type, only : field_struct
        use odp, only : allocate_field, deallocate_field
        use simulation_type, only : potentials_struct
        use element_type, only : element_struct
        use atom_type, only : atom_struct
        use fft, only: recip_to_real
        use parallel_mod, only: parallel_task
        use grids_mod, only : allocate_local_fields_G


        complex(dp), intent(inout) :: Cn(:)

        type(parallel_struct), intent(in) :: parallel
        type(potentials_struct), intent(in) :: potentials
        type(atom_struct), intent(inout) :: atoms(:)
        type(element_struct), intent(inout) :: elements(:)
        type(orbital_struct), intent(inout) :: orbitals(:)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer:: grid
        complex(dp), intent(inout) :: Mn(:)
        real(dp), intent(in) :: Ebar, deltaE
        type(field_struct), allocatable :: Qn(:), Qnm1(:)
        complex(dp), allocatable :: Q0(:,:,:,:),Qnm2(:,:,:,:)
        integer :: i,s, n

                !A Chebychev Filter approach will require all orbitals to have same number of components
        allocate(Qn(orbitals(1)%n_spinor))
        allocate(Qnm1(orbitals(1)%n_spinor))
        grid=>grids(orbitals(1)%of(1)%grid)

        call allocate_local_fields_G(Q0, grid,orbitals(1)%n_spinor)
        call allocate_local_fields_G(Qnm2, grid,orbitals(1)%n_spinor)

        do s=1, orbitals(1)%n_spinor
            Qn(s)%grid=orbitals(1)%of(s)%grid
            call allocate_field(Qn(s),grid,parallel)
            Qnm1(s)%grid=orbitals(1)%of(s)%grid
            call allocate_field(Qnm1(s),grid,parallel)
        enddo

        Mn(:)=0.0_dp
        if(size(Mn).ne.size(Cn)) then
            print *, 'Inconsistent size of Cn and Mn Chebychev arrays, stopping'
            stop
        endif
        do i=1, size(orbitals)
            grid=>grids(orbitals(1)%of(1)%grid)
            do n=1, size(Cn)
                if(n.eq.1) then
                    do s=1, orbitals(i)%n_spinor
                        Q0(:,:,:,s)=orbitals(i)%of(s)%G(:,:,:)
                        Qn(s)%G=Q0(:,:,:,s)
                        orbitals(i)%of(s)%G=0.0_dp
                    enddo
                else
                    call Apply_H(Qnm1, Qn, grids, potentials, atoms, elements, parallel, calc_R=.true.)
                    do s=1, orbitals(i)%n_spinor
                        Qn(s)%G=(Qn(s)%G-Ebar*Qnm1(s)%G)/deltaE
                        if(n.gt.2) Qn(s)%G=2.0*Qn(s)%G - Qnm2(:,:,:,s)
                    enddo
                endif
                do s=1, orbitals(i)%n_spinor
                    if(n.gt.1) Qnm2(:,:,:,s)=Qnm1(s)%G
                    Qnm1(s)%G=Qn(s)%G
                    orbitals(i)%of(s)%G=orbitals(i)%of(s)%G+Cn(n)*Qn(s)%G
                    Mn(n)=Mn(n) + orbitals(i)%occ(s)*orbitals(i)%weight(s)* &
                                        integrate_3D_G(Qn(s)%G*conjg(Q0(:,:,:,s)), grid, parallel)
                enddo
            enddo
            do s=1, orbitals(i)%n_spinor
                call recip_to_real(orbitals(i)%of(s), grids)
            enddo
        enddo
        call parallel_task('sum', Mn(:), parallel, 'diff_b')
        
        deallocate(Q0)
        deallocate(Qnm2)
        do s=1, orbitals(1)%n_spinor
            call deallocate_field(Qn(s),grids(Qn(s)%grid))
            call deallocate_field(Qnm1(s),grids(Qnm1(s)%grid))
        enddo
        deallocate(Qn)
        deallocate(Qnm1)

    end subroutine

    subroutine cheby_exp_and_app_Gen(orbitals,  &
        Ebar, deltaE, Mn, Cn, grids, potentials, atoms, elements, parallel, all_PAW)
    use odp_type, only : field_struct, orbital_struct
    use grids_type, only : grid_struct
    use parallel_type, only : parallel_struct
    use operations_3D, only:  integrate_3D_R
    use Apply_Hamiltonian, only: Apply_SinvH
    use Non_Local_ion, only : Apply_S_RR
    use odp_type, only : field_struct
    use odp, only : allocate_field, deallocate_field
    use simulation_type, only : potentials_struct
    use element_type, only : element_struct
    use atom_type, only : atom_struct
    use fft, only: recip_to_real, real_to_recip
    use parallel_mod, only: parallel_task
    use simulation_type, only : all_PAW_struct
    use grids_mod, only : allocate_local_fields_R


    complex(dp), intent(inout) :: Cn(:)

    type(parallel_struct), intent(in) :: parallel
    type(potentials_struct), intent(in) :: potentials
    type(atom_struct), intent(inout) :: atoms(:)
    type(element_struct), intent(inout) :: elements(:)
    type(orbital_struct), intent(inout) :: orbitals(:)
    type(grid_struct), intent(inout), target :: grids(:)
    type(grid_struct), pointer:: grid
    complex(dp), intent(inout) :: Mn(:)
    real(dp), intent(in) :: Ebar, deltaE
    type(field_struct), allocatable :: Qn(:), Qnm1(:)
    complex(dp), allocatable :: Q0(:,:,:,:),Qnm2(:,:,:,:)
    type(all_PAW_struct), intent(inout) :: all_PAW

    integer :: i,s, n

            !A Chebychev Filter approach will require all orbitals to have same number of components
    allocate(Qn(orbitals(1)%n_spinor))
    allocate(Qnm1(orbitals(1)%n_spinor))
    grid=>grids(orbitals(1)%of(1)%grid)

    call allocate_local_fields_R(Q0,grid,orbitals(1)%n_spinor)
    call allocate_local_fields_R(Qnm2,grid,orbitals(1)%n_spinor)

    do s=1, orbitals(1)%n_spinor
        Qn(s)%grid=orbitals(1)%of(s)%grid
        call allocate_field(Qn(s),grid,parallel)
        Qnm1(s)%grid=orbitals(1)%of(s)%grid
        call allocate_field(Qnm1(s),grid,parallel)
    enddo

    Mn(:)=0.0_dp
    if(size(Mn).ne.size(Cn)) then
        print *, 'Inconsistent size of Cn and Mn Chebychev arrays, stopping'
        stop
    endif
    do i=1, size(orbitals)
        grid=>grids(orbitals(1)%of(1)%grid)
        do n=1, size(Cn)
            if(n.eq.1) then
                call Apply_S_RR(orbitals(i)%of, Qnm1, grids, atoms, elements, parallel)
                do s=1, orbitals(i)%n_spinor
                    Qn(s)%R=orbitals(i)%of(s)%R(:,:,:)
                    Q0(:,:,:,s)=Qnm1(s)%R
                    orbitals(i)%of(s)%R=0.0_dp
                    Qnm1(s)%R=0.0_dp
                enddo
            else
                call Apply_SinvH(Qnm1, Qn, grids, potentials, atoms, elements, parallel, all_PAW, calc_G=.true.)!, with_CG=.true.)
                do s=1, orbitals(i)%n_spinor
                    Qn(s)%R=(Qn(s)%R-Ebar*Qnm1(s)%R)/deltaE
                    if(n.gt.2) Qn(s)%R=2.0_dp*Qn(s)%R - Qnm2(:,:,:,s)
                enddo
            endif
            do s=1, orbitals(i)%n_spinor
                if(n.gt.1) Qnm2(:,:,:,s)=Qnm1(s)%R
                Qnm1(s)%R=Qn(s)%R
                orbitals(i)%of(s)%R=orbitals(i)%of(s)%R+Cn(n)*Qn(s)%R
                Mn(n)=Mn(n) + orbitals(i)%occ(s)*orbitals(i)%weight(s)* &
                                    integrate_3D_R(Qn(s)%R*conjg(Q0(:,:,:,s)), grid, parallel)
            enddo
        enddo
        do s=1, orbitals(i)%n_spinor
            call real_to_recip(orbitals(i)%of(s), grids)
            orbitals(i)%of(s)%G=orbitals(i)%of(s)%G*grid%cutwf
            call recip_to_real(orbitals(i)%of(s), grids)
        enddo
    enddo
    call parallel_task('sum', Mn(:), parallel, 'diff_b')

    deallocate(Q0)
    deallocate(Qnm2)
    do s=1, orbitals(1)%n_spinor
        call deallocate_field(Qn(s),grids(Qn(s)%grid))
        call deallocate_field(Qnm1(s),grids(Qnm1(s)%grid))
    enddo
    deallocate(Qn)
    deallocate(Qnm1)

end subroutine

    subroutine calculate_Chebychev_coefficients_fd(Ebar, Delta, Cn, Nc, T_au, mu)
        use fourier, only : fct                                                                       
        real(dp), intent(in) :: Ebar, Delta
        real(dp), intent(inout) :: Cn(:)
        integer, intent(in) :: Nc 
        real(dp), intent(in) ::T_au, mu
        real(dp) ::  expon
 
        real(dp) :: x_k, F_k(2*Nc)                                                       
        integer :: i ,j                                                                                 
        Cn=0.0_dp                                                                                     
        do i=1,Nc                                                                              
            x_k=cos(pi*(i*1.0_dp-0.5_dp)/Nc)*Delta+Ebar   
            expon=(x_k-mu)/T_au
             if(expon.gt.50_dp) then
                F_k(i)=0.0 
             else
                F_k(i)=1.0_dp/(exp(expon)+1.0_dp)
             endif
        enddo                                                                                    
        do i=1,Nc
            Cn(i)=0.0_dp
        do j=1,Nc
            Cn(i)= Cn(i) + F_k(j)*cos(pi*(i-1)*(j*1.0_dp-0.5_dp)/Nc)
        enddo
            Cn(i)=Cn(i)*2.0_dp/Nc
        enddo
            Cn(1)=Cn(1)*0.5_dp   
    end subroutine

    subroutine calculate_Chebychev_coefficients_omfd(Ebar, Delta, Cn, Nc, T_au, mu)
        use fourier, only : fct                                                                       
        real(dp), intent(in) :: Ebar, Delta
        real(dp), intent(inout) :: Cn(:)
        integer, intent(in) :: Nc 
        real(dp), intent(in) ::T_au, mu
        real(dp) ::  expon
 
        real(dp) :: x_k, F_k(2*Nc)                                                       
        integer :: i ,j                                                                                 
        Cn=0.0_dp                                                                                     
        do i=1,Nc                                                                              
            x_k=cos(pi*(i*1.0_dp-0.5_dp)/Nc)*Delta+Ebar   
            expon=(x_k-mu)/T_au
             if(expon.gt.50_dp) then
                F_k(i)=1.0 
             else
                F_k(i)=1.0_dp-1.0_dp/(exp(expon)+1.0_dp)
             endif
        enddo                                                                                    
        do i=1,Nc
            Cn(i)=0.0_dp
        do j=1,Nc
            Cn(i)= Cn(i) + F_k(j)*cos(pi*(i-1)*(j*1.0_dp-0.5_dp)/Nc)
        enddo
            Cn(i)=Cn(i)*2.0_dp/Nc
        enddo
            Cn(1)=Cn(1)*0.5_dp   
    end subroutine

    subroutine calculate_Chebychev_coefficients_root_omfd(Ebar, Delta, Cn, Nc, T_au, mu)
        use fourier, only : fct                                                                       
        real(dp), intent(in) :: Ebar, Delta
        real(dp), intent(inout) :: Cn(:)
        integer, intent(in) :: Nc 
        real(dp), intent(in) ::T_au, mu
        real(dp) ::  expon
 
        real(dp) :: x_k, F_k(2*Nc)                                                       
        integer :: i ,j                                                                                 
        Cn=0.0_dp                                                                                     
        do i=1,Nc                                                                              
            x_k=cos(pi*(i*1.0_dp-0.5_dp)/Nc)*Delta+Ebar   
            expon=(x_k-mu)/T_au
             if(expon.gt.50_dp) then
                F_k(i)=1.0 
             else
                F_k(i)=sqrt(1.0_dp-1.0_dp/(exp(expon)+1.0_dp))
             endif
        enddo                                                                                    
        do i=1,Nc
            Cn(i)=0.0_dp
        do j=1,Nc
            Cn(i)= Cn(i) + F_k(j)*cos(pi*(i-1)*(j*1.0_dp-0.5_dp)/Nc)
        enddo
            Cn(i)=Cn(i)*2.0_dp/Nc
        enddo
            Cn(1)=Cn(1)*0.5_dp   
    end subroutine
    subroutine calculate_Chebychev_coefficients_sqrt_fd(Ebar, Delta, Cn, Nc, T_au, mu)
        use fourier, only : fct                                                                       
        real(dp), intent(in) :: Ebar, Delta
        real(dp), intent(inout) :: Cn(:)
        integer, intent(in) :: Nc 
        real(dp), intent(in) ::T_au, mu
        real(dp) ::  expon
 
        real(dp) :: x_k, F_k(2*Nc)                                                       
        integer :: i,j 

        Cn=0.0_dp                                                                                     
        do i=1,Nc                                                                              
            x_k=cos(pi*(i*1.0_dp-0.5_dp)/Nc)*Delta+Ebar   
            expon=(x_k-mu)/T_au
             if(expon.gt.50_dp) then
                F_k(i)=0.0 
             else
                F_k(i)=sqrt(1.0_dp/(exp(expon)+1.0_dp))
             endif
        enddo                                                                                    
        do i=1,Nc
            Cn(i)=0.0_dp
        do j=1,Nc
            Cn(i)= Cn(i) + F_k(j)*cos(pi*(i-1)*(j*1.0_dp-0.5_dp)/Nc)
        enddo
            Cn(i)=Cn(i)*2.0_dp/Nc
       enddo
            Cn(1)=Cn(1)*0.5_dp   
            
    end subroutine

    subroutine calculate_Chebychev_coefficients_Gaussian(Ebar, Delta, Cn, Nc, w, a)
        use fourier, only : fct 
        use constants, only : pi                                                                      
        real(dp), intent(in) :: Ebar, Delta
        real(dp), intent(inout) :: Cn(:)
        integer, intent(in) :: Nc 
        real(dp), intent(in) ::w, a
        real(dp) ::  expon
 
        real(dp) :: x_k, F_k(2*Nc)                                                       
        integer :: i,j 

        Cn=0.0_dp                                                                                     
        do i=1,Nc                                                                              
            x_k=cos(pi*(i*1.0_dp-0.5_dp)/Nc)*Delta+Ebar   
            expon=((x_k-w)**2/a)
             if(expon.gt.50_dp) then
                F_k(i)=0.0 
             else
                F_k(i)=exp(-expon)/(sqrt(pi*a))
             endif
        enddo                                                                                    
        do i=1,Nc
            Cn(i)=0.0_dp
        do j=1,Nc
            Cn(i)= Cn(i) + F_k(j)*cos(pi*(i-1)*(j*1.0_dp-0.5_dp)/Nc)
        enddo
            Cn(i)=Cn(i)*2.0_dp/Nc
       enddo
            Cn(1)=Cn(1)*0.5_dp   
            
    end subroutine

    subroutine calculate_Chebychev_coefficients_FD_Gaussian(Ebar, Delta, Cn, Nc, w, a, mu, T_au)
        use fourier, only : fct 
        use constants, only : pi                                                                      
        real(dp), intent(in) :: Ebar, Delta
        real(dp), intent(inout) :: Cn(:)
        integer, intent(in) :: Nc 
        real(dp), intent(in) ::w, a, mu, T_au
        real(dp) ::  expon
 
        real(dp) :: x_k, F_k(2*Nc)                                                       
        integer :: i,j 

        Cn=0.0_dp                                                                                     
        do i=1,Nc                                                                              
            x_k=cos(pi*(i*1.0_dp-0.5_dp)/Nc)*Delta+Ebar   
            expon=((x_k-w)**2/a)
             if(expon.gt.25_dp) then
                F_k(i)=0.0 
             else
                F_k(i)=exp(-expon)/(sqrt(pi*a))
             endif
             expon=(x_k-mu)/T_au
             if(expon.gt.25_dp) then
                F_k(i)=0.0 
             else
                F_k(i)=1.0_dp/(exp(expon)+1.0_dp) * F_k(i)
             endif
        enddo                                                                                    
        do i=1,Nc
            Cn(i)=0.0_dp
        do j=1,Nc
            Cn(i)= Cn(i) + F_k(j)*cos(pi*(i-1)*(j*1.0_dp-0.5_dp)/Nc)
        enddo
            Cn(i)=Cn(i)*2.0_dp/Nc
       enddo
            Cn(1)=Cn(1)*0.5_dp   
            
    end subroutine

    subroutine calculate_Chebychev_coefficients_propagator(Ebar, Delta, Cn, Nc, dt)
        real(dp), intent(in) :: Ebar, Delta
        complex(dp), intent(inout) :: Cn(:)
        integer, intent(in) :: Nc 
        real(dp), intent(in) :: dt
 
        real(dp) :: x_k
        complex(dp) ::  F_k(2*Nc)                                                       
        integer :: i ,j                                                                                 
        Cn=0.0_dp                                                                                     
        do i=1,Nc                                                                              
            x_k=cos(pi*(i*1.0_dp-0.5_dp)/Nc)*Delta+Ebar   
            F_k(i)=exp(-i_*x_k*dt)
        enddo                                                                                    
        do i=1,Nc
            Cn(i)=0.0_dp
        do j=1,Nc
            Cn(i)= Cn(i) + F_k(j)*cos(pi*(i-1)*(j*1.0_dp-0.5_dp)/Nc)
        enddo
            Cn(i)=Cn(i)*2.0_dp/Nc
        enddo
            Cn(1)=Cn(1)*0.5_dp   
    end subroutine

    subroutine calculate_Chebychev_coefficients_entropy(Ebar, Delta, Cn, Nc, T_au, mu)
        use fourier, only : fct                                                                       
        real(dp), intent(in) :: Ebar, Delta
        real(dp), intent(inout) :: Cn(:)
        integer, intent(in) :: Nc 
        real(dp), intent(in) ::T_au, mu
        real(dp) ::  expon, FD
 
        real(dp) :: x_k, F_k(2*Nc)                                                       
        integer :: i,j 
        Cn=0.0_dp                                                                                     
        do i=1,Nc                                                                              
            x_k=cos(pi*(i*1.0_dp-0.5_dp)/Nc)*Delta+Ebar   
            expon=(x_k-mu)/T_au
             if(expon.gt.50_dp) then
                F_k(i)=0.0 
             else
                FD=1.0_dp/(exp(expon)+1.0_dp)
                if(abs(FD-1.0_dp).lt.1E-8) then
                    F_k(i)=0.0
                    else
                    F_k(i)=FD*log(FD)+(1.0_dp-FD)*log(1.0_dp-FD)
                    endif
             endif
        enddo                                                                                    
        do i=1,Nc
            Cn(i)=0.0_dp
        do j=1,Nc
            Cn(i)= Cn(i) + F_k(j)*cos(pi*(i-1)*(j*1.0_dp-0.5_dp)/Nc)
        enddo
            Cn(i)=Cn(i)*2.0_dp/Nc
       enddo
            Cn(1)=Cn(1)*0.5_dp   
            
    end subroutine

    

end module