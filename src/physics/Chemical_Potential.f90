module Chemical_Potential
    use types, only : dp

    implicit none

    public :: Find_Chemical_Potential


    contains

    subroutine Find_Chemical_Potential(Minimum_eigenvalue, all_eigs, all_occ, all_filter, all_degeneracy, all_weight, &
        stochastic, temperature, chem_potential, system, grids, nelec, parallel, success)
        use constants, only : Ha2eV
        use odp_type, only : orbital_struct
        use grids_type, only : grid_struct
        use parallel_type, only : parallel_struct
        use parallel_mod, only : parallel_task
        use system_type, only : system_struct
        use stochastic_type, only : stochastic_struct
        use Chebychev, only : calculate_Chebychev_coefficients_fd, &
        calculate_Chebychev_coefficients_sqrt_fd


        real(dp) , intent(in) :: all_eigs(:,:), nelec, all_filter(:,:), all_weight(:), temperature, Minimum_eigenvalue(:)
        integer, intent(in) :: all_degeneracy(:,:)
        real(dp) , intent(inout) :: all_occ(:,:), chem_potential

        type(grid_struct), intent(inout) :: grids(:)
        type(parallel_struct), intent(in) :: parallel
        type(system_struct), intent(in) :: system
        type(stochastic_struct), intent(inout) :: stochastic(:)
        logical, intent(inout) :: success(:)

        real(dp) :: ne_det, ne_stoc, mu, mu_high, mu_low, expon
        real(dp) :: mu_nm1, ne_n, ne_nm1, tmp
        integer :: n, iter, Nc_max(system%n_kpoints), power, k, k_loc, success_int(system%n_kpoints)
        real(dp), allocatable :: stoc_nelec(:)

        allocate(stoc_nelec(size(stochastic)))
        do k_loc=1,size(stochastic)
            if(stochastic(k_loc)%do_stoc) then
                Nc_max(k_loc)=5*(floor(stochastic(k_loc)%deltaE/system%temperature)+1)+10
                Nc_max(k_loc)=max(Nc_max(k_loc),100) !lower bound on Nc 
                power=1
                do while(power<Nc_max(k_loc))
                    power=power*2
                enddo
                Nc_max(k_loc)=power
                if(allocated(stochastic(k_loc)%Density_Cn)) deallocate(stochastic(k_loc)%Density_Cn)
                allocate(stochastic(k_loc)%Density_Cn(Nc_max(k_loc)))
            endif
        enddo
        mu=chem_potential
        stochastic(:)%nelec_calc=0.0_dp
        stoc_nelec(:)=0.0_dp

        ne_det=0.0_dp
        ne_stoc=0.0_dp
        mu_high=grids(2)%Ecut
        mu_low=minval(Minimum_eigenvalue(:))-1.0_dp
        iter=0
        ne_n=ne_det+ne_stoc
        do while(abs(ne_n-nelec).gt.1E-6_dp)
            iter=iter+1
            ne_det=0.0_dp
            if(system%n_deterministic.gt.0) then
                all_occ(:,:)=0.0_dp
                if(parallel%myid_spin.eq.0) then
                    do k=1,system%n_kpoints
                    do n=1, system%n_deterministic
                            expon=(all_eigs(n,k)-mu)/temperature
                            if(expon.lt.50_dp.and.expon.gt.-50.0_dp) then
                                all_occ(n,k)=all_filter(n,k)*all_degeneracy(n,k)/(exp(expon)+1.0_dp)
                                ne_det=ne_det+all_occ(n,k)*all_weight(k)
                            else if(expon.lt.-50.0_dp) then
                                all_occ(n,k)=all_filter(n,k)*all_degeneracy(n,k)
                                ne_det=ne_det+all_occ(n,k)*all_weight(k)
                            endif
                    enddo
                    enddo
                endif
                call parallel_task('bcast', ne_det, parallel, 'spin', root=0)
                call parallel_task('bcast', all_occ, parallel, 'spin', root=0)
            endif

            if(stochastic(1)%do_stoc) then
                ne_stoc=0.0_dp
                success_int(:)=1
                if(parallel%myid_k.eq.0) then
                    do k_loc=1,size(stochastic)
                            k=stochastic(k_loc)%k_point
                            call calculate_Chebychev_coefficients_fd(stochastic(k_loc)%Ebar, stochastic(k_loc)%deltaE, &
                                                stochastic(k_loc)%Density_Cn(:), Nc_max(k_loc), system%temperature, mu)
                            stochastic(k_loc)%nelec_calc=0.0_dp
                            do n=1, min(size(stochastic(k_loc)%Mn(:)),Nc_max(k_loc))
                                    stochastic(k_loc)%nelec_calc=stochastic(k_loc)%nelec_calc&
                                        +real(stochastic(k_loc)%Mn(n))*stochastic(k_loc)%Density_Cn(n)
                            enddo
                            if(stochastic(k_loc)%nelec_calc.lt.0.0_dp) then
                                success_int(k)=0
                            endif
                    enddo
                    ne_stoc=sum(stochastic(:)%nelec_calc)
                endif
                call parallel_task('prd', success_int, parallel, 'spin')
                call parallel_task('sum', ne_stoc, parallel, 'spin')

                success(:)=.true.
                do k=1,size(success_int)
                    if(success_int(k).eq.0) success(k)=.false.
                enddo

                if(any(.not.success(:))) exit !Need to exit and go re-compute the stochastic with larger Energy range
            endif

            ne_n=ne_stoc+ne_det
            if(parallel%myid_spin.eq.0) then
                if(iter.gt.100) &
                print *,  iter, 'ne_stoc: ', ne_stoc, 'ne_det: ',  ne_det, &
                      'ne_system: ',  nelec, 'difference: ', ne_n - nelec, 'mu:', mu
                if(iter.gt.1000) then
                    print *, 'Problem with Chemical Potential search'
                    stop
                endif
                if(abs(ne_n-nelec).le.1E-6_dp) then
                    success(:)=.true.
                else
                    success(:)=.false.                
                    !Adjust Chemical Potential
                    if(system%temperature.gt.0.5_dp/Ha2eV) then
                        if(iter.eq.1) then
                            ne_nm1=ne_n
                            mu_nm1=mu
                            mu=mu-(ne_n-nelec)/abs((ne_n-nelec))*min(abs((ne_n-nelec))/nelec,0.1_dp) &
                                *abs((mu-minval(Minimum_eigenvalue(:))))
                        else
                            if(iter.eq.2) then
                                !if 10% T shift moved too far away from root, swap so original mu is mu_n
                                if(abs(ne_n-nelec).gt.abs(ne_nm1-nelec)) then
                                    tmp=mu_nm1
                                    mu_nm1=mu
                                    mu=tmp
                                    tmp=ne_nm1
                                    ne_nm1=ne_n
                                    ne_n=tmp
                                endif
                            endif
                            !Secant Step
                            tmp=mu-(ne_n-nelec)*(mu-mu_nm1)/(ne_n-ne_nm1)
                        !  print *, tmp, mu, (ne_n-nelec), (mu-mu_nm1), (ne_n-ne_nm1)
                            mu_nm1=mu
                            mu=tmp
                            ne_nm1=ne_n
                        endif
                    else
                        if((ne_stoc+ne_det).gt.nelec) then
                                mu_high=mu
                                mu=(mu+mu_low)/2.0_dp
                        else
                                mu_low=mu
                                mu=(mu+mu_high)/2.0_dp
                        endif

                        if(abs((ne_stoc+ne_det)-nelec).gt.1E-3_dp) then
                            if(abs(mu-mu_low).lt.1E-6) then
                                mu_low=mu_low-abs(mu_low) -0.1
                            else if (abs(mu-mu_high).lt.1E-6) then
                                mu_high=mu_high+abs(mu_high) +0.1
                                exit
                            endif
                        endif
                    endif
                endif

                
            endif
            
            call parallel_task('bcast', mu, parallel, 'spin', root=0)
            call parallel_task('bcast', success, parallel, 'spin', root=0)
            stoc_nelec(:)=stochastic(:)%nelec_calc
            call parallel_task('bcast', stoc_nelec, parallel, 'k', root=0)
            stochastic(:)%nelec_calc=stoc_nelec(:)
        enddo
        chem_potential=mu

    end subroutine

end module