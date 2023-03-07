module Thomas_Fermi
    use types, only : op, dp
    use constants, only : pi
    implicit none
    public :: Thomas_Fermi_Energy, calc_Thomas_Fermi_Potential, Thomas_Fermi_Entropy, calc_Dynamic_KEDP_and_energy

    contains

    real(dp) elemental function f(y, deriv,sec_deriv)
            ! Function f(y) from Appendix A in [1].
            !
            ! [1] Perrot, F. (1979). Gradient correction to the statistical electronic free
            ! energy at nonzero temperatures: Application to equation-of-state
            ! calculations. Physical Review A, 20(2), 586–594.
            real(dp), intent(in) :: y ! must be positive
            ! if deriv == .true. compute df/dy instead. Default .false.
            logical, intent(in), optional :: deriv, sec_deriv
            real(dp), parameter :: y0 = 3*pi/(4*sqrt(2._dp))
            real(dp), parameter :: c(8) = [-0.8791880215_dp, 0.1989718742_dp, &
                0.1068697043e-2_dp, -0.8812685726e-2_dp, 0.1272183027e-1_dp, &
                -0.9772758583e-2_dp, 0.3820630477e-2_dp, -0.5971217041e-3_dp]
            real(dp), parameter :: d(9) = [0.7862224183_dp, -0.1882979454e1_dp, &
                0.5321952681_dp, 0.2304457955e1_dp, -0.1614280772e2_dp, &
                0.5228431386e2_dp, -0.9592645619e2_dp, 0.9462230172e2_dp, &
                -0.3893753937e2_dp]
            real(dp) :: u, dudy, d2udy2
            integer :: i
            logical :: deriv_, sec_deriv_
            deriv_ = .false.
            sec_deriv_ = .false.
            if (present(deriv)) deriv_ = deriv
            if (present(sec_deriv)) sec_deriv_ = sec_deriv
            if (.not. sec_deriv_) then
            if (.not. deriv_) then
                    ! Desire f(y)
                    if (y <= y0) then
                        f = log(y)
                        do i = 0, 7
                            f = f + c(i+1) * y**i
                        end do
                    else
                        u = y**(2._dp / 3)
                        f = 0
                        do i = 0, 8
                            f = f + d(i+1) / u**(2*i-1)
                        end do
                    ! Note: Few terms in [1] have "y" instead of "u" in them for y > y0, but
                    ! that is obviously a typo.
                    end if
                else 
                    ! Desire df(y)/dy
                    if (y <= y0) then
                        f = 1 / y
                        do i = 1, 7
                            f = f + (i) * c(i+1) * y**(i-1)
                        end do
                    else
                        u = y**(2._dp / 3)
                        dudy = 2._dp/3 / y**(1._dp/3)
                        f = 0
                        do i = 0, 8
                            f = f + (1-2*i) * d(i+1) / u**(2*i)
                        end do
                        f = f * dudy
                    end if
                end if
            else
                ! Desire d^2f(y)/dy^2
                if (y <= y0) then
                    f = -1 / y**2
                    do i = 2, 7
                        f = f + i * (i-1) * c(i+1) * y**(i-2)
                    end do
                else
                    u = y**(2._dp / 3)
                    dudy = 2._dp/3 / y**(1._dp/3)
                    dudy = dudy**2
                    d2udy2 = -2._dp/9 / y**(4._dp/3)
                    f = 0
                    do i = 1, 8
                        f = f - dudy*(1-2*i) * (2*i) *  d(i+1) / u**(2*i+1)
                    end do
                    do i = 0, 8
                        f = f + d2udy2*(1-2*i) * d(i+1) / u**(2*i)
                    end do 
                end if       
            end if
        end function

        function Thomas_Fermi_Energy(density, grids, parallel, system, tf) result (Energy)
            use system_type, only : system_struct
            use odp_type, only : spin_DenPot_struct 
            use parallel_type, only : parallel_struct
            use grids_type, only : grid_struct
            use operations_3D, only: integrate_3D_R
            use tf_type, only: tf_struct

            type(spin_DenPot_struct), intent(inout) :: density
            type(grid_struct), intent(inout), target :: grids(:)
            type(grid_struct), pointer :: grid
            type(system_struct), intent(in) :: system
            type(parallel_struct), intent(inout) :: parallel
            type(tf_struct), intent(in) :: tf
            real(dp) :: dy0dn
            real(dp) :: Energy
            integer :: s 
            grid=>grids(density%of(1)%grid)
    
            dy0dn = pi**2 / sqrt(2._dp) * system%temperature**(-3._dp/2)
            Energy=0.0_dp
            do s=1,size(density%of)
                Energy=Energy+real(integrate_3D_R( &
                    density%n_s*density%of(s)%R*f(dy0dn*real(density%of(s)%R)*density%n_s)*system%temperature, grid, parallel))
            enddo
            Energy=tf%gamma*Energy/density%n_s
        end function

        subroutine calc_Thomas_Fermi_Potential(density, potentials, system, tf)
            use system_type, only : system_struct
            use odp_type, only : spin_DenPot_struct 
            use simulation_type, only: potentials_struct
            use operations_3D, only: integrate_3D_R
            use tf_type, only: tf_struct

            type(potentials_struct), intent(inout) :: potentials
            type(spin_DenPot_struct), intent(inout) :: density
            type(tf_struct), intent(in) :: tf

            type(system_struct), intent(in) :: system
            real(dp) :: dy0dn
            integer :: s 
    
            dy0dn = pi**2 / sqrt(2._dp) * system%temperature**(-3._dp/2)
            do s=1,size(density%of)
                potentials%kinetic_local%of(s)%R=system%temperature*(f(dy0dn*real(density%of(s)%R)*density%n_s) &
                    + real(density%of(s)%R)*density%n_s*dy0dn*f(dy0dn*real(density%of(s)%R)*density%n_s, deriv=.true.))
                    potentials%kinetic_local%of(s)%R=potentials%kinetic_local%of(s)%R*tf%gamma
            enddo
        end subroutine

        function Thomas_Fermi_Entropy(density, grids, parallel, system) result (Entrop)
            use system_type, only : system_struct
            use odp_type, only : spin_DenPot_struct 
            use parallel_type, only : parallel_struct
            use grids_type, only : grid_struct
            use simulation_type, only: energies_struct, potentials_struct
            use operations_3D, only: integrate_3D_R
    
            type(spin_DenPot_struct), intent(inout) :: density
            type(grid_struct), intent(inout), target :: grids(:)
            type(grid_struct), pointer :: grid
            type(system_struct), intent(in) :: system
            type(parallel_struct), intent(inout) :: parallel
            real(dp) :: dy0dn
            real(dp) :: Entrop
            integer :: s 
            grid=>grids(density%of(1)%grid)
    
            dy0dn = pi**2 / sqrt(2._dp) * system%temperature**(-3._dp/2)
            Entrop=0.0_dp
            do s=1,size(density%of)
                Entrop=Entrop-real(integrate_3D_R( &
                    density%n_s*density%of(s)%R*( &
                        f(dy0dn*real(density%n_s*density%of(s)%R)) + &
                        system%temperature*f(dy0dn*real(density%n_s*density%of(s)%R),deriv=.true.)* &
                        dy0dn*density%n_s*density%of(s)%R*(-1.5_dp*system%temperature**(-2.5_dp))) &
                    , grid, parallel))
            enddo
            Entrop=-Entrop*system%temperature/density%n_s
        end function

    subroutine calc_Dynamic_KEDP_and_energy(KEDP, KE, current_density, density, system, grids, parallel)
        use system_type, only : system_struct
        use odp_type, only : spin_DenPot_struct 
        use simulation_type, only: potentials_struct
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use simulation_type, only: energies_struct, potentials_struct
        use fft, only: real_to_recip, recip_to_real
        use operations_3D, only: integrate_3D_R
        use constants, only : i_, pi

        type(spin_DenPot_struct), intent(inout) :: density, current_density(3)
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(system_struct), intent(in) :: system
        type(parallel_struct), intent(in) :: parallel
        type(spin_DenPot_struct), intent(inout) :: KEDP
        real(dp),intent(out) :: KE
        integer :: dir, s
        real(dp) :: min_den

        do s=1, density%n_s
         KEDP%of(s)%G=0.0_dp
        enddo
        do dir=1,3
            do s=1, density%n_s
                grid=>grids(current_density(dir)%of(s)%grid)
                call real_to_recip(current_density(dir)%of(s), grids)
                where(grid%G2>tiny(1.0_dp))
                    KEDP%of(s)%G =  KEDP%of(s)%G + i_*grid%G(:,:,:,dir)*current_density(dir)%of(s)%G &
                                                        / sqrt(grid%G2(:,:,:))*grid%cutden
                elsewhere
                    KEDP%of(s)%G = 0.0_dp
                end where 
            enddo
        enddo
        KE=0.0_dp
        do s=1,density%n_s
            grid=>grids(density%of(s)%grid)
            min_den=0.001_dp*system%nelec(s)/product(grid%box_length)
            KEDP%of(s)%G =KEDP%of(s)%G * density%n_s
            call recip_to_real(KEDP%of(s), grids)
            KEDP%of(s)%R=KEDP%of(s)%R*(pi**3)/2.0_dp &
                /(3.0_dp*(pi**2.0_dp)*density%n_s*(density%of(s)%R+min_den))**(2.0/3.0)
            KE=KE+real(integrate_3D_R( &
                    KEDP%of(s)%R*((0.75**3.6_dp + (1.69271_dp*sqrt(2.0_dp*system%temperature) &
                    /(3.0_dp*(pi**2.0_dp)*density%n_s*(density%of(s)%R(:,:,:)+min_den))**(1.0/3.0))**3.6_dp)**(1._dp/3.6_dp)) &
                    *density%of(s)%R, &
                grid, parallel))
            KEDP%of(s)%R=KEDP%of(s)%R*((1.0_dp + (1.69271_dp*sqrt(2.0_dp*system%temperature) &
            /(3.0_dp*(pi**2.0_dp)*density%n_s*(density%of(s)%R+min_den))**(1.0/3.0))**3.6_dp)**(1._dp/3.6_dp))
            call real_to_recip(KEDP%of(s), grids)
            KEDP%of(s)%G=KEDP%of(s)%G*grid%cutden
            call recip_to_real(KEDP%of(s), grids)
        enddo
        KE=KE/density%n_s
    end subroutine
end module