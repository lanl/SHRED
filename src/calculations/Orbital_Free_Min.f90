module Orbital_Free_Min
    use types, only : dp
    use constants, only : pi
    implicit none
    public :: orbital_free_energy, orbital_free_energy_min

    interface
        real(dp) function func(x)
        import :: dp
        implicit none
        real(dp), intent(in) :: x
        end function
    end interface

    contains
   
    subroutine orbital_free_energy(system, all_PAW, density, grids,  parallel, potentials, xc, energies, tf, &
                                   calc_value, calc_derivative, vW)
        use system_type, only : system_struct
        use odp_type, only : orbital_struct, spin_DenPot_struct 
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use xc_type, only: xc_struct
        use hartree, only : calc_Hartree_potential, Hartree_Energy
        use Local_ion, only : Calc_Local_PP_Energy
        use xc_mod, only : calc_XC
        use simulation_type, only: energies_struct, potentials_struct
        use operations_3D, only: integrate_3D_R, integrate_3D_G
        use fft, only : real_to_recip, recip_to_real
        use odp, only: allocate_orbital, deallocate_orbital
        use tf_type, only: tf_struct
        use simulation_type, only : all_PAW_struct
        use Thomas_Fermi, only : Thomas_Fermi_Energy, calc_Thomas_Fermi_Potential
        implicit none

        type(tf_struct), intent(inout) :: tf
        type(spin_DenPot_struct), intent(inout) :: density
        type(all_PAW_struct), intent(inout) :: all_PAw
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(system_struct), intent(in) :: system
        type(parallel_struct), intent(inout) :: parallel
        type(xc_struct), intent(inout) :: xc(:)
        type(orbital_struct) :: psi, d2psi
        type(potentials_struct), intent(inout) :: potentials
        type(energies_struct), intent(inout) :: energies
        logical, intent(in) :: calc_value, calc_derivative
        ! include von Weizsäcker term, default .false.
        logical, intent(in), optional :: vW
        integer :: i
        logical :: vW_
        

        grid=>grids(density%of(1)%grid)

        if(.not.calc_value .and. .not.calc_derivative) then
            print *, 'Need to ask orbital_free_energy for energy or derivative or both'
            stop
        endif
        call calc_Hartree_potential(potentials%Hartree, density, grids)
        call calc_XC(density, xc, all_PAW, parallel, grids, energies%xc, potentials%xc)
        vW_ = .false.
        if (present(vW)) vW_ = vW
        if (vW_) then
            psi%n_spinor=density%n_s
            call allocate_orbital(psi, grids, density%of(1)%grid, parallel)
            do i= 1, psi%n_spinor
                psi%of(i)%R = sqrt(density%of(i)%R)
                call real_to_recip(psi%of(i), grids)
            enddo
        end if
   
        if (calc_value) then
            energies%hartree=Hartree_Energy(potentials%hartree, density, parallel, grids)
            energies%ion_local=Calc_Local_PP_Energy(potentials%ion_local, density, grids, parallel)
            !xc energies calculated in Calc_XC
            energies%kinetic_local=Thomas_Fermi_Energy(density, grids, parallel, system, tf)
            energies%kinetic=0.0_dp
            do i=1, density%n_s
                if (vW_) energies%kinetic=energies%kinetic+tf%lambda* &
                    integrate_3D_G(abs(psi%of(i)%G)**2*grid%G2, grid, parallel)*0.5_dp
            enddo
            energies%total = energies%hartree + energies%ion_local + energies%xc + energies%kinetic_local + energies%kinetic
        end if
        
        if (calc_derivative) then
            if (vW_) then
                d2psi%n_spinor=density%n_s
                call allocate_orbital(d2psi, grids, density%of(1)%grid, parallel)
            endif
            call calc_Thomas_Fermi_Potential(density, potentials, system, tf)
            do i=1, d2psi%n_spinor
                if (vW_) then
                    d2psi%of(i)%G=-grid%G2*psi%of(i)%G
                    call recip_to_real(d2psi%of(i), grids)
                    potentials%kinetic_local%of(i)%R= potentials%kinetic_local%of(i)%R  &
                    - density%n_s*tf%lambda*d2psi%of(i)%R/(2.0_dp*psi%of(i)%R)
                endif
                
                potentials%total_local_fine%of(i)%R= potentials%kinetic_local%of(i)%R + &
                                                potentials%xc%of(i)%R + &
                                                potentials%ion_local%of%R + &
                                                potentials%hartree%of%R + &
                                                potentials%external_field%of%R
                if(all_PAW%usepotzero.gt.0)  potentials%total_local_fine%of(i)%R = &
                                                potentials%total_local_fine%of(i)%R + sum(all_PAW%vpotzero) 
            enddo
            if (vW_) call deallocate_orbital(d2psi,grids)
        end if           
        if (vW_) call deallocate_orbital(psi,grids)
end subroutine

    subroutine orbital_free_energy_min(system, density, all_PAW, grids,  parallel, potentials, xc, energies, &
            tf)
        use system_type, only : system_struct
        use odp_type, only : orbital_struct, spin_DenPot_struct 
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use xc_type, only: xc_struct
        use tf_type, only: tf_struct
        use hartree, only : calc_Hartree_potential, Hartree_Energy
        use Local_ion, only : Calc_Local_PP_Energy
        use xc_mod, only : calc_XC
        use simulation_type, only: energies_struct, potentials_struct
        use operations_3D, only: integrate_3D_R
        use fft, only: real_to_recip, recip_to_real
        use odp, only: allocate_DenPot, deallocate_DenPot
        use simulation_type, only : all_PAW_struct
        use grids_mod, only : allocate_local_fields_R


        type(spin_DenPot_struct), intent(inout) :: density
        type(grid_struct), intent(inout), target :: grids(:)
        type(grid_struct), pointer :: grid
        type(system_struct), intent(in) :: system
        type(parallel_struct), intent(inout) :: parallel
        type(xc_struct), intent(inout) :: xc(:)
        type(tf_struct), intent(inout) :: tf
        type(spin_DenPot_struct) :: psi, psi_
        real(dp), dimension(:,:,:,:), allocatable :: psi_prev, phi, phi_prime, ksi, ksi_prev, eta, Hpsi
        type(potentials_struct), intent(inout) :: potentials
        type(energies_struct), intent(inout) :: energies
        type(all_PAW_struct), intent(inout) :: all_PAw

        
        integer, parameter :: update_fletcher_reeves = 1
        integer, parameter :: update_polak_ribiere   = 2
        integer :: iter
        real(dp) :: last2, last3, brent_eps, free_energy_, &
                    gamma_d, gamma_n, theta, theta_a, theta_b, theta_c, fa, fb, fc
        real(dp), allocatable :: mu(:), free_energies(:)
        real(dp) :: f2
        real(dp) :: psi_norm
        integer ::  i
        logical :: vW_

        grid=>grids(density%of(1)%grid)

        brent_eps=tf%brent_eps
        vW_=.false.
        if(abs(tf%lambda).gt.tiny(1.0_dp)) vW_=.true.

        last2 = 0
        last3 = 0
        psi%n_s=density%n_s
        psi_%n_s=density%n_s
 
        call allocate_DenPot(psi, grids, density%of(1)%grid, parallel)
        call allocate_DenPot(psi_, grids, density%of(1)%grid, parallel)
        psi%of(:)%grid=density%of(1)%grid
        psi_%of(:)%grid=density%of(1)%grid

        allocate(mu(density%n_s))
        call allocate_local_fields_R(psi_prev,grid,density%n_s)
        call allocate_local_fields_R(phi,grid,density%n_s)
        call allocate_local_fields_R(phi_prime,grid,density%n_s)
        call allocate_local_fields_R(ksi,grid,density%n_s)
        call allocate_local_fields_R(ksi_prev,grid,density%n_s)
        call allocate_local_fields_R(eta,grid,density%n_s)
        call allocate_local_fields_R(Hpsi,grid,density%n_s)


        do i= 1, density%n_s
            psi%of(i)%R = sqrt(density%of(i)%R)
            call real_to_recip(psi%of(i), grids)
            psi%of(i)%G=psi%of(i)%G*grid%cutwf
            call recip_to_real(psi%of(i), grids)
            psi_norm = integrate_3D_R( abs(psi%of(i)%R)**2, grid, parallel)
            if (parallel%myid == 0) print *, "Initial norm of psi:", psi_norm
            psi%of(i)%R = psi%of(i)%R * sqrt( system%nelec(i) / psi_norm)
            psi%of(i)%G = psi%of(i)%G * sqrt( system%nelec(i) / psi_norm)
            psi_norm = integrate_3D_R( abs(psi%of(i)%R)**2, grid, parallel)
            if (parallel%myid == 0) print *, "norm of psi:", psi_norm
            density%of(i)%R = abs(psi%of(i)%R)**2
            call real_to_recip(density%of(i), grids)
        enddo


        ! This returns H[n] = delta F / delta n, we save it to the Hpsi variable to
        ! save space:
        call orbital_free_energy(system,all_PAW, density, grids,  parallel, potentials, xc, energies, tf, &
            calc_value=.true., calc_derivative=.true., vW=vW_)
        ! Hpsi = H[psi] = delta F / delta psi = 2*H[n]*psi, due to d/dpsi = 2 psi d/dn
        do i= 1, density%n_s
            Hpsi(:,:,:,i)= 2.0_dp * real(potentials%total_local_fine%of(i)%R * psi%of(i)%R)
            mu(i) = 0.5_dp / system%nelec(i) * integrate_3D_R(real(psi%of(i)%R) * Hpsi(:,:,:,i), grid,  parallel)
            ksi(:,:,:,i)=2.0_dp*mu(i)*real(psi%of(i)%R) -  Hpsi(:,:,:,i)
            phi(:,:,:,i)=ksi(:,:,:,i)
            phi_prime(:,:,:,i) = phi(:,:,:,i) - 1.0_dp / system%nelec(i) * &
                                integrate_3D_R(real(phi(:,:,:,i) * psi%of(i)%R), grid,  parallel) *real(psi%of(i)%R)
            eta(:,:,:,i) = sqrt( system%nelec(i) / integrate_3D_R( phi_prime(:,:,:,i)**2, grid,  parallel)) &
                             * phi_prime(:,:,:,i)
        enddo

        !print *, "Summary of energies [a.u.]:"
        !print "('    Ts   = ', f14.8)", Ts
        !print "('    Een  = ', f14.8)", Een
        !print "('    Eee  = ', f14.8)", Eee
        !print "('    Exc  = ', f14.8)", Exc
        !print *, "   ---------------------"
        !print "('    Etot = ', f14.8, ' a.u.')", free_energy_
        allocate(free_energies(tf%max_iter))
        gamma_n = 0
        theta = pi/2
        do iter = 1, tf%max_iter
            !Hpsi = Hpsi/(2*psi)
            ! Formula (36) in Jiang & Yang
            !A = integral(L, psi*Hpsi*psi) - integral(L, eta*Hpsi*eta)
            !B = 2*integral(L, eta*Hpsi*psi)
            !print *, "theta? =", 0.5_dp * atan(B/A)
            theta_a = 0
            theta_b = mod(theta, 2*pi)
            if (iter == 1) then
                call bracket(func_fft, theta_a, theta_b, theta_c, fa, fb, fc, 100._dp, 20, verbose=.false.)
                call brent(func_fft, theta_a, theta_b, theta_c, brent_eps, 50, theta, &
                    free_energy_, verbose=.false.)
            else
                call bracket(func_fft, theta_a, theta_b, theta_c, fa, fb, fc, 100._dp, 20, verbose=.false.)
                call parabola_vertex(theta_a, fa, theta_b, fb, theta_c, fc, theta, f2)
            end if
            do i= 1, density%n_s
                psi%of(i)%R = cos(theta)* psi%of(i)%R + sin(theta) * eta(:,:,:,i)
                call real_to_recip(psi%of(i), grids)
                density%of(i)%R = abs(psi%of(i)%R)**2
                call real_to_recip(density%of(i), grids)
            enddo
            ! TODO: We probably don't need to recalculate free_energy_ here:

            call orbital_free_energy(system,all_PAW, density, grids,  parallel, potentials, xc, energies, tf, &
                 calc_value=.true., calc_derivative=.true., vW=vW_)
           if(.false.) then
            if(parallel%myid.eq.0) print *, "Iteration:", iter
            ! print *, "Norm of psi:", integrate_3D_R(real(psi(1)%of%R* psi(1)%of%R), grid,  parallel)
            if(parallel%myid.eq.0) print *, "mu =", mu
            if(parallel%myid.eq.0) print *, "|ksi| =", sqrt(gamma_n)
            if(parallel%myid.eq.0) print *, "theta =", theta
            if(parallel%myid.eq.0) print *, "Summary of energies [a.u.]:"
            if(parallel%myid.eq.0) print "('    Ts   = ', f14.8)", energies%kinetic_local 
            if(parallel%myid.eq.0) print "('    vW   = ', f14.8)", energies%kinetic 
            if(parallel%myid.eq.0) print "('    Een  = ', f14.8)", energies%ion_local
            if(parallel%myid.eq.0) print "('    Eee  = ', f14.8)", energies%hartree
            if(parallel%myid.eq.0) print "('    Exc  = ', f14.8)", energies%xc
            if(parallel%myid.eq.0) print *, "   ---------------------"
           endif
            free_energies(iter) = energies%total / system%n_atoms
            free_energy_=energies%total / system%n_atoms
            if (iter > 1) then
                last2 = maxval(free_energies(iter-1:iter)) - &
                    minval(free_energies(iter-1:iter))
            end if
            if (iter > 2) then
                last3 = maxval(free_energies(iter-2:iter)) - &
                    minval(free_energies(iter-2:iter))
            end if
        !    if (myid == 0) then
        !        print "('# ', i3, ' Etot/atom = ', f18.8, ' eV; last2 = ', es10.2, ' last3 = ',es10.2)", &
        !        iter, free_energy_ * Ha2eV / Natom, last2 * Ha2eV, last3 * Ha2eV
        !    end if
            if (iter > 3) then
                if (last3 < tf%energy_eps) then
                    do i= 1, density%n_s
                        density%of(i)%R = abs(psi%of(i)%R)**2
                        call real_to_recip(density%of(i), grids)
                    enddo
                    call deallocate_DenPot(psi,grids)
                    call deallocate_DenPot(psi_,grids)
                    if(parallel%myid.eq.0) print "('    Ts   = ', f14.8)", energies%kinetic_local 
                    if(parallel%myid.eq.0) print "('    vW   = ', f14.8)", energies%kinetic 
                    if(parallel%myid.eq.0) print "('    Een  = ', f14.8)", energies%ion_local
                    if(parallel%myid.eq.0) print "('    Eee  = ', f14.8)", energies%hartree
                    if(parallel%myid.eq.0) print "('    Exc  = ', f14.8)", energies%xc
                    if(parallel%myid.eq.0) print *, "   ---------------------"
                    return
                end if
            end if
            do i= 1, density%n_s
                Hpsi(:,:,:,i)= 2.0_dp * real(potentials%total_local_fine%of(i)%R) * real(psi%of(i)%R)
                ksi_prev(:,:,:,i)=ksi(:,:,:,i)
                mu(i) = 0.5_dp / system%nelec(i) * integrate_3D_R(real(psi%of(i)%R) * Hpsi(:,:,:,i), grid,  parallel)
                ksi(:,:,:,i)=2.0_dp*mu(i)*real(psi%of(i)%R) - Hpsi(:,:,:,i)
                select case(tf%update_type)
                    case(update_fletcher_reeves) ! Fletcher-Reeves
                        gamma_n = integrate_3D_R(ksi(:,:,:,i)**2, grid,  parallel)
                    case(update_polak_ribiere)   ! Polak-Ribiere
                        gamma_n = max(integrate_3D_R(ksi(:,:,:,i)*(ksi(:,:,:,i)-ksi_prev(:,:,:,i)), grid,  parallel), 0.0_dp)
                    case default
                        print *, "Unknown update type."
                        stop
                end select
                gamma_d = integrate_3D_R(ksi_prev(:,:,:,i)**2, grid,  parallel)
                phi(:,:,:,i)=ksi(:,:,:,i) + gamma_n / gamma_d * phi(:,:,:,i)
                phi_prime(:,:,:,i) = phi(:,:,:,i) - 1.0_dp / system%nelec(i) * &
                                    integrate_3D_R(real(phi(:,:,:,i) * psi%of(i)%R), grid,  parallel)*real(psi%of(i)%R)
                eta(:,:,:,i)= sqrt(system%nelec(i) / integrate_3D_R(phi_prime(:,:,:,i)**2, grid,  parallel)) * phi_prime(:,:,:,i)
            enddo
        end do
        print *, "free_energy_minimization: The maximum number of iterations exceeded."
        call deallocate_DenPot(psi,grids)
        call deallocate_DenPot(psi_,grids)
        stop

        contains
            !Note that there is only one theta, there should be a theta for each spin
            !this should be fixed 
            real(dp) function func_fft(theta) result(energy)
                use fft, only: real_to_recip
                real(dp), intent(in) :: theta

                do i= 1, density%n_s
                    psi_%of(i)%R = cos(theta)* psi%of(i)%R + sin(theta) * eta(:,:,:,i)

                    call real_to_recip(psi_%of(i), grids)
                    density%of(i)%R = abs(psi_%of(i)%R)**2
                    call real_to_recip(density%of(i), grids)
                enddo
                
                call orbital_free_energy(system, all_PAW, density, grids,  parallel, potentials, xc, energies, tf, &
                calc_value=.true., calc_derivative=.false., vW=vW_)
                energy=energies%total
            end function

            subroutine brent(f, xa, xb, xc, tol, maxiter, xmin, fxmin, verbose)
                procedure(func) :: f
                real(dp), intent(in) :: xa, xb, xc, tol
                integer, intent(in) :: maxiter
                real(dp), intent(out) :: xmin, fxmin
                logical, intent(in), optional :: verbose
                real(dp), parameter :: mintol = 1e-3_dp*epsilon(xa), &
                    golden_ratio = (1+sqrt(5._dp))/2, cg = 2._dp - golden_ratio
                real(dp) :: tol1, tol2, tmp1, tmp2, deltax, dx_temp, rat, xmid, &
                        a, b, p, u, v, w, x, fu, fv, fw, fx
                integer :: iter
                logical :: verbose_
                verbose_ = .false.
                if (present(verbose)) verbose_ = verbose
                
                x = xb; w = xb; v = xb
                fx = f(x)
                fw = fx; fv = fx
                if (xa < xc) then
                    a = xa
                    b = xc
                else
                    a = xc
                    b = xa
                end if
                deltax = 0
                rat = 0
                do iter = 1, maxiter
                    if (verbose_) then
                        print "(i2, ':  x = ', f23.12, '     tol = ', es10.2)", iter, x, &
                            (b - a) / 2
                    end if
                    tol1 = tol*abs(x) + mintol
                    tol2 = 2*tol1
                    xmid = 0.5_dp*(a + b)
                    if (abs(x - xmid) < (tol2 - 0.5_dp*(b - a))) then ! check for convergence
                        xmin = x
                        fxmin = fx
                        return
                    end if
                    if (abs(deltax) <= tol1) then
                        if (x >= xmid) then           ! do a golden section step
                            deltax = a - x
                        else
                            deltax = b - x
                        end if
                        rat = cg*deltax
                    else
                        tmp1 = (x - w)*(fx - fv)      ! do a parabolic step
                        tmp2 = (x - v)*(fx - fw)
                        p = (x - v)*tmp2 - (x - w)*tmp1;
                        tmp2 = 2*(tmp2 - tmp1)
                        if (tmp2 > 0) then
                            p = -p
                        end if
                        tmp2 = abs(tmp2)
                        dx_temp = deltax
                        deltax = rat
                        if ((p > tmp2*(a - x)) .and. (p < tmp2*(b - x)) .and. &
                                (abs(p) < abs(0.5_dp*tmp2*dx_temp))) then ! check parabolic fit
                            rat = p / tmp2            ! if parabolic step is useful.
                            u = x + rat
                            if ((u - a) < tol2 .or. (b - u) < tol2) then
                                if (xmid - x >= 0) then
                                    rat = tol1
                                else
                                    rat = -tol1
                                end if
                            end if
                        else
                            if (x >= xmid) then       ! if it's not do a golden section step
                                deltax = a - x
                            else
                                deltax = b - x
                            end if
                            rat = cg*deltax
                        end if
                    end if
                
                    if (abs(rat) < tol1) then         ! update by at least tol1
                        if (rat >= 0) then
                            u = x + tol1
                        else
                            u = x - tol1
                        end if
                    else
                        u = x + rat
                    end if
                    fu = f(u) ! calculate new output value
                
                    if (fu > fx) then                 ! if it's bigger than current
                        if (u < x) then
                            a = u
                        else
                            b = u
                        end if
                        if ((fu <= fw) .or. abs(w - x) < tiny(1._dp)) then
                            v = w; w = u; fv = fw; fw = fu
                        else if ((fu <= fv) .or. abs(v - x) < tiny(1._dp) .or. abs(v - w) < tiny(1._dp)) then
                            v = u; fv = fu
                        end if
                    else
                        if (u >= x) then
                            a = x
                        else
                            b = x
                        end if
                        v = w; w = x; x = u
                        fv = fw; fw = fx; fx = fu
                    end if
                end do
                print *, "brent: The maximum number of iterations exceeded."
                stop
            end subroutine
                
            subroutine bracket(f, xa, xb, xc, fa, fb, fc, grow_limit, maxiter, verbose)
                procedure(func) :: f
                real(dp), intent(inout) :: xa, xb
                real(dp), intent(out) :: xc, fa, fb, fc
                real(dp), intent(in) :: grow_limit
                integer, intent(in) :: maxiter
                logical, intent(in), optional :: verbose
                real(dp), parameter :: golden_ratio = (1+sqrt(5._dp))/2
                real(dp) :: denom, dum, fw, tmp1, tmp2, val, w, wlim
                integer :: iter
                logical :: verbose_
                verbose_ = .false.
                if (present(verbose)) verbose_ = verbose
                
                fa = f(xa)
                fb = f(xb)
                if (fa < fb) then                      ! Switch so fa > fb
                    dum = xa; xa = xb; xb = dum
                    dum = fa; fa = fb; fb = dum
                end if
                xc = xb + golden_ratio*(xb - xa)
                fc = f(xc)
                iter = 0
                do while (fc < fb)
                    tmp1 = (xb - xa)*(fb - fc)
                    tmp2 = (xb - xc)*(fb - fa)
                    val = tmp2 - tmp1
                    if (abs(val) < tiny(1.0_dp)) then
                        denom = 2*tiny(1.0_dp)
                    else
                        denom = 2*val
                    end if
                    w = xb - ((xb - xc)*tmp2 - (xb - xa)*tmp1) / denom
                    wlim = xb + grow_limit*(xc - xb)
                    if (iter > maxiter) then
                        print *, "Too many iterations."
                        stop
                    endif
                    iter = iter + 1
                    if ((w - xc)*(xb - w) > 0) then
                        fw = f(w)
                        if (fw < fc) then
                            xa = xb; xb = w; fa = fb; fb = fw
                            return
                        else if (fw > fb) then
                            xc = w; fc = fw
                            return
                        end if
                        w = xc + golden_ratio*(xc - xb)
                        fw = f(w)
                    else if ((w - wlim)*(wlim - xc) >= 0) then
                        w = wlim
                        fw = f(w)
                    else if ((w - wlim)*(xc - w) > 0) then
                        fw = f(w)
                        if (fw < fc) then
                            xb = xc; xc = w; w = xc + golden_ratio*(xc - xb)
                            fb = fc; fc = fw; fw = f(w)
                        end if
                    else
                        w = xc + golden_ratio*(xc - xb)
                        fw = f(w)
                    end if
                    xa = xb; xb = xc; xc = w
                    fa = fb; fb = fc; fc = fw
                    if (verbose_) then
                        print "(i2, ':  xa = ', f23.12, ' xb = ', f23.12, ' xc = ', f23.12)", &
                            iter, xa, xb, xc
                    end if
                end do
            end subroutine
        
            subroutine parabola_vertex(x1, y1, x2, y2, x3, y3, xv, yv)
                ! Calculates the position (xv, yv) of the parabola vertex.
                ! The parabola is defined by 3 points (x1, y1), (x2, y2), (x3, y3).
                real(dp), intent(in) :: x1, y1, x2, y2, x3, y3
                real(dp), intent(out) :: xv, yv
                real(dp) :: denom, A, B, C
                denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
                A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
                B     = (x3**2 * (y1 - y2) + x2**2 * (y3 - y1) + x1**2 * (y2 - y3)) / denom
                C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + &
                            x1 * x2 * (x1 - x2) * y3) / denom
                xv = -B / (2*A)
                yv = C - B**2 / (4*A)
            end subroutine
    end subroutine

    
        
end module
