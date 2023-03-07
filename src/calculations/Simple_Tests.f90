module Simple_Tests  

use types, only : dp
    use constants, only : pi, i_
    implicit none
    
    
    public :: test_apply_Hamiltonian, test_current

    contains


    subroutine test_apply_Hamiltonian(sim)
        use simulation_type, only : simulation_struct
        use fft, only: real_to_recip, recip_to_real
        use Apply_Hamiltonian, only : Apply_H
        use parallel_mod,only : parallel_wait
        use Non_Local_ion, only : Calculate_deriv_Projector_R_overlaps
        use constants, only : i_
        use Update_Potentials_and_Energies, only : fine_to_coarse_all_potentials
        use operations_3D, only : integrate_3D_G

        type(simulation_struct), intent(inout) :: sim
        integer :: i, j,l
        complex(dp) :: norm

        if(sim%grids(3)%FFT_type.eq.1) then
            sim%orbitals(1,1,1)%of(1)%R=real(sin(sim%grids(3)%R(:,:,:,1))*cos(sim%grids(3)%R(:,:,:,2))&
            *tan(sim%grids(3)%R(:,:,:,3)))
        else
            sim%orbitals(1,1,1)%of(1)%R=real(tan(sim%grids(3)%R(:,:,:,1))*cos(sim%grids(3)%R(:,:,:,2))&
            *sin(sim%grids(3)%R(:,:,:,3)))
        endif


        sim%potentials%total_local_fine%of(1)%R=real(sim%potentials%ion_local%of%R) + &
        sim%potentials%hartree%of%R  + &
        sim%potentials%xc%of(1)%R

        call fine_to_coarse_all_potentials(sim%potentials, sim%grids, sim%parallel, sim%fine_to_coarse)
        !sim%potentials%total_local%of(1)%R=0.0_dp!real(sim%potentials%ion_local%of%R) !+ &
        !sim%potentials%hartree%of%R + &
        !sim%potentials%xc%of(1)%R

        call real_to_recip(sim%orbitals(1,1,1)%of(1), sim%grids)
        sim%orbitals(1,1,1)%of(1)%G=sim%orbitals(1,1,1)%of(1)%G*sim%grids(3)%cutwf

        call recip_to_real(sim%orbitals(1,1,1)%of(1), sim%grids)
        do l=1,1
        call Apply_H(sim%orbitals(1,1,1)%of(:), sim%orbitals(2,1,1)%of(:), &
                    sim%grids,sim%potentials,sim%atoms, sim%elements,sim%parallel, calc_R=.true.)
           ! sim%orbitals(1,1,1)%of(1)%G=sim%orbitals(2,1,1)%of(1)%G   
           ! norm=integrate_3D_G(abs(sim%orbitals(1,1,1)%of(1)%G)**2,sim%grids(3), sim%parallel)
           ! sim%orbitals(1,1,1)%of(1)%G=sim%orbitals(1,1,1)%of(1)%G/sqrt(norm)
        enddo
        call recip_to_real(sim%orbitals(2,1,1)%of(1), sim%grids)
        call parallel_wait(sim%parallel)

                do j = 0 ,sim%parallel%nproc-1
                    do i=1, sim%grids(3)%Nr_local(3)
                        if(sim%parallel%myid.eq.j) print *, i, sim%grids(3)%R(1,1,i,3), &
                                                real(sim%orbitals(2,1,1)%of(1)%R(1,1,i)),&
                                                 sim%orbitals(2,1,1)%of(1)%G(1,1,i)
                                                flush(6)
                    enddo
                    call parallel_wait(sim%parallel)
                enddo
                flush(6)
                call parallel_wait(sim%parallel)
                if(sim%parallel%myid.eq.0) print *, ''; flush(6)
                call parallel_wait(sim%parallel)

                do j = 0 ,sim%parallel%nproc-1
                    if(sim%grids(3)%FFT_type.gt.1) then
                        do i=1, sim%grids(3)%Ng_local(1)
                            if(sim%parallel%myid.eq.j) print *, 'x', i, sim%grids(3)%G(2,2,i,1), &
                            sim%orbitals(2,1,1)%of(1)%G(2,2,i), abs(sim%orbitals(2,1,1)%of(1)%G(2,2,i)),&
                             abs(sim%orbitals(2,1,1)%of(1)%G(2,2,i))/sqrt(2.0_dp)
                        enddo
                        do i=1, sim%grids(3)%Ng_local(3)
                            if(sim%parallel%myid.eq.j) print *, 'z', i, sim%grids(3)%G(i,2,2,3), &
                            sim%orbitals(2,1,1)%of(1)%G(i,2,2), abs(sim%orbitals(2,1,1)%of(1)%G(i,2,2)),&
                             abs(sim%orbitals(2,1,1)%of(1)%G(i,2,2))/sqrt(2.0_dp)
                        enddo
                    else
                        do i=1, sim%grids(3)%Ng_local(3)
                            if(sim%parallel%myid.eq.j) print *, 'z', i, sim%grids(3)%G(2,2,i,3), &
                            sim%orbitals(2,1,1)%of(1)%G(2,2,i), abs(sim%orbitals(2,1,1)%of(1)%G(2,2,i)), &
                            abs(sim%orbitals(2,1,1)%of(1)%G(2,2,i))/sqrt(2.0_dp)
                        enddo
                    endif
                    flush(6)
                    call parallel_wait(sim%parallel)
                enddo
                call parallel_wait(sim%parallel)
        stop
    end subroutine

    subroutine test_grid_transfer(sim)
        use simulation_type, only : simulation_struct
        use grids_mod, only:  Gspace_grid_to_grid_transfer
        use parallel_mod, only : parallel_wait

        type(simulation_struct), target :: sim

        integer ::  n3, n1, k

        sim%potentials%total_local_fine%of(1)%G=sim%grids(1)%G2

       
        call Gspace_grid_to_grid_transfer(sim%grids(1), sim%grids(2), &
         sim%parallel, sim%potentials%total_local_fine%of(1)%G, &
         sim%potentials%total_local%of(1)%G, sim%fine_to_coarse)

        if(sim%grids(2)%FFT_type.gt.1) then
            n1=3;n3=1
        else
            n1=1;n3=3
        endif
        if(sim%parallel%myid_band.eq.0) then

            do k=1, sim%grids(2)%Ng_local(n3)
                print *, '1:',  real(sim%potentials%total_local%of(1)%G(1,1,k)), sim%grids(2)%G2(1,1,k)
                flush(6)
            enddo
            print *; flush(6)

            do k=1, sim%grids(2)%Ng_local(2)
                print *, '2:', real(sim%potentials%total_local%of(1)%G(1,k,1)), sim%grids(2)%G2(1,k,1)
                flush(6)
            enddo
            print *; flush(6)

            do k=1, sim%grids(2)%Ng_local(n1)
                print *, '3:', real(sim%potentials%total_local%of(1)%G(k,1,1)), sim%grids(2)%G2(k,1,1)
                flush(6)
            enddo
            print *; flush(6)
            print *; flush(6)
            print *; flush(6)


         endif
         call parallel_wait(sim%parallel)
         stop
    end subroutine

    subroutine test_grid_to_matrix(sim)
        use simulation_type, only : simulation_struct
        use grids_mod, only:  Gspace_grid_to_grid_transfer
        use parallel_mod, only : parallel_wait

        type(simulation_struct), target :: sim

        integer ::  n3, n1, k

        sim%potentials%total_local_fine%of(1)%G=sim%grids(1)%G2

       
        call Gspace_grid_to_grid_transfer(sim%grids(1), sim%grids(2), &
         sim%parallel, sim%potentials%total_local_fine%of(1)%G, &
         sim%potentials%total_local%of(1)%G, sim%fine_to_coarse)

        if(sim%grids(2)%FFT_type.gt.1) then
            n1=3;n3=1
        else
            n1=1;n3=3
        endif
        if(sim%parallel%myid_band.eq.0) then

            do k=1, sim%grids(2)%Ng_local(n3)
                print *, '1:',  real(sim%potentials%total_local%of(1)%G(1,1,k)), sim%grids(2)%G2(1,1,k)
                flush(6)
            enddo
            print *; flush(6)

            do k=1, sim%grids(2)%Ng_local(2)
                print *, '2:', real(sim%potentials%total_local%of(1)%G(1,k,1)), sim%grids(2)%G2(1,k,1)
                flush(6)
            enddo
            print *; flush(6)

            do k=1, sim%grids(2)%Ng_local(n1)
                print *, '3:', real(sim%potentials%total_local%of(1)%G(k,1,1)), sim%grids(2)%G2(k,1,1)
                flush(6)
            enddo
            print *; flush(6)
            print *; flush(6)
            print *; flush(6)


         endif
         call parallel_wait(sim%parallel)
         stop
    end subroutine

    subroutine test_new_FFT(sim)
        use, intrinsic :: iso_c_binding
        use simulation_type, only : simulation_struct
        use grids_type, only : grid_struct
        use fft, only : real_to_recip, recip_to_real
        use parallel_mod, only : parallel_wait
        type(simulation_struct), target :: sim
        type(grid_struct), pointer :: grid

        real(dp) :: time_1, time_2

        integer :: ix,iy,iz, l, i,j,k
        
        !transposed FFT of slow variavle (y of F(x,y))
       ! plan_fft_xy_y_lot = fftw_plan_magrid%Nr(2)_dft(rank, [grid%Nr(2)], int(loty), &
       ! in_xy, [grid%Nr(2)], nx, 1,                          &
       ! zw, [ncache], int(loty), 1, FFTW_FORWARD, FFTW_MEASURE)

        !if(mod(nx,int(loty)).gt.0) &
        !plan_fft_xy_y_left = fftw_plan_magrid%Nr(2)_dft(rank, [grid%Nr(2)], mod(nx,int(loty)), &
        !in_xy, [grid%Nr(2)], nx, 1,                          &
        !zw, [ncache], mod(nx,int(loty)), 1, FFTW_FORWARD, FFTW_MEASURE)

        !Normal FFT of fast variable (x of F(x,y))
        !plan_fft_xy_x_1_lot = fftw_plan_magrid%Nr(2)_dft(rank, [nx], int(lotx), &
        !in_xy, [nx], 1, nx,                          &
        !zw, [ncache], 1, nx, FFTW_FORWARD, FFTW_MEASURE)

        !if(mod(grid%Nr(2),int(lotx)).gt.0) &
        !plan_fft_xy_x_1_left = fftw_plan_magrid%Nr(2)_dft(rank, [nx], mod(grid%Nr(2),int(lotx)), &
        !in_xy, [nx], 1, nx,                          &
        !zw, [ncache], 1, nx, FFTW_FORWARD, FFTW_MEASURE)
    

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        grid=>sim%grids(3)

        !XY plane waves
        sim%orbitals(1,1,1)%of(1)%R=0.0_dp
        do iz=1, sim%grids(3)%Nr_local(3)
        do iy=1, sim%grids(3)%Nr_local(2)
        do ix=1, sim%grids(3)%Nr_local(1)

                sim%orbitals(1,1,1)%of(1)%R(ix,iy,iz)=&
                sim%orbitals(1,1,1)%of(1)%R(ix,iy,iz)+ &
        (2-1+100*(5-1)+10000*(7-1))*exp(i_*(ix-1)*sim%grids(3)%dR(1)*(2-1)*sim%grids(3)%dG(1))&
                                   *exp(i_*(iy-1)*sim%grids(3)%dR(2)*(5-1)*sim%grids(3)%dG(2))&
                                   *exp(i_*(iz+sim%grids(3)%myxyz_r(3)*sim%grids(3)%Nr_local(3)-1)&
                                        *sim%grids(3)%dR(3)*(7-1)*sim%grids(3)%dG(3))
        enddo;enddo; enddo;


        if(sim%grids(3)%gamma) sim%orbitals(1,1,1)%of(1)%R=(sim%orbitals(1,1,1)%of(1)%R+conjg(sim%orbitals(1,1,1)%of(1)%R))

        if(sim%parallel%myid_band.eq.0) then

            do k=1, sim%grids(3)%Nr_local(3)
                print *, 'z:',  sim%orbitals(1,1,1)%of(1)%R(1,1,k)/60401.0_dp
                flush(6)
            enddo
            print *; flush(6)

            do k=1, sim%grids(3)%Nr_local(2)
                print *, 'y:', sim%orbitals(1,1,1)%of(1)%R(1,k,1)/60401.0_dp
                flush(6)
            enddo
            print *; flush(6)

            do k=1, sim%grids(3)%Nr_local(1)
                print *, 'x:',sim%orbitals(1,1,1)%of(1)%R(k,1,1)/60401.0_dp
                flush(6)
            enddo
            print *; flush(6)
            print *; flush(6)
            print *; flush(6)


        endif
        
        call parallel_wait(sim%parallel)
     
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call cpu_time(time_1)
        !FFT  then transpose
        do l=1,10000
            call real_to_recip(sim%orbitals(1,1,1)%of(1), sim%grids)
            call recip_to_real(sim%orbitals(1,1,1)%of(1), sim%grids)
        enddo

        if(sim%parallel%myid_band.eq.0) then

            do k=1, sim%grids(3)%Nr_local(3)
                print *, 'z:',  sim%orbitals(1,1,1)%of(1)%R(1,1,k)/60401.0_dp
                flush(6)
            enddo
            print *; flush(6)

            do k=1, sim%grids(3)%Nr_local(2)
                print *, 'y:', sim%orbitals(1,1,1)%of(1)%R(1,k,1)/60401.0_dp
                flush(6)
            enddo
            print *; flush(6)

            do k=1, sim%grids(3)%Nr_local(1)
                print *, 'x:',sim%orbitals(1,1,1)%of(1)%R(k,1,1)/60401.0_dp
                flush(6)
            enddo
            print *; flush(6)
            print *; flush(6)
            print *; flush(6)


        endif
        call parallel_wait(sim%parallel)

        call cpu_time(time_2)
        print *, 'FFT Forward new x10000:', time_2-time_1
        call parallel_wait(sim%parallel)

        do k=1, grid%Ng_local(3); do j=1, grid%Ng_local(2); do i=1, grid%Ng_local(1)

            if(grid%FFT_type.gt.1) then
                if(abs(sim%orbitals(1,1,1)%of(1)%G(k,j,i)).gt.1.0_dp) then
                print *,&
                        sim%parallel%myid , k,j,i , sim%orbitals(1,1,1)%of(1)%G(k,j,i)/sqrt(2.0_dp)
                        flush(6)
                endif
                if(k.eq.7.and.j.eq.5.and.i.eq.2) print *, '2:', sim%parallel%myid ,  k,j,i , &
                     sim%orbitals(1,1,1)%of(1)%G(k,j,i)/sqrt(2.0_dp)
            else
                if(abs(sim%orbitals(1,1,1)%of(1)%G(i,j,k)).gt.1.0_dp) then
                    print *,&
                            sim%parallel%myid , i,j,k , sim%orbitals(1,1,1)%of(1)%G(i,j,k)/sqrt(2.0_dp)
                            flush(6)
                    endif
                if(k.eq.7.and.j.eq.5.and.i.eq.2) print *, '2:', sim%parallel%myid , i,j,k , &
                     sim%orbitals(1,1,1)%of(1)%G(i,j,k)/sqrt(2.0_dp)
            endif
        enddo;enddo;enddo
        call parallel_wait(sim%parallel)

        call cpu_time(time_1)
        do l=1,10000
            call recip_to_real(sim%orbitals(1,1,1)%of(1), sim%grids)
            !stop
        enddo
        call cpu_time(time_2)

        call parallel_wait(sim%parallel)

        call cpu_time(time_2)
        print *, 'FFT New Back time x10000:', time_2-time_1
        print *; flush(6)
        call parallel_wait(sim%parallel)



        if(sim%parallel%myid_band.eq.0) then

            do k=1, sim%grids(3)%Nr_local(3)
                print *, 'z:',  sim%orbitals(1,1,1)%of(1)%R(1,1,k)/60401.0_dp
                flush(6)
            enddo
            print *; flush(6)

            do k=1, sim%grids(3)%Nr_local(2)
                print *, 'y:', sim%orbitals(1,1,1)%of(1)%R(1,k,1)/60401.0_dp
                flush(6)
            enddo
            print *; flush(6)

            do k=1, sim%grids(3)%Nr_local(1)
                print *, 'x:',sim%orbitals(1,1,1)%of(1)%R(k,1,1)/60401.0_dp
                flush(6)
            enddo
            print *; flush(6)
            print *; flush(6)
            print *; flush(6)

        endif
        call parallel_wait(sim%parallel)
            
    end subroutine


    subroutine test_distribution(sim)
        use simulation_type, only : simulation_struct
        use fft, only: real_to_recip, recip_to_real
        use Apply_Hamiltonian, only : Apply_H
        use grids_type, only : grid_struct
        use parallel_mod,only : parallel_wait
        use constants, only : i_
        use odp_type, only : field_struct
        use odp, only : allocate_field, deallocate_field
        use FFT2Matrix, only : field_FFT2Matrix, field_Matrix2FFT
        use distribute_field

        type(simulation_struct), intent(inout), target :: sim
        integer :: i, j
        complex(dp), allocatable :: orbital_full(:,:,:)
        complex(dp), allocatable :: orbital_full_G(:,:,:)
        real(dp), allocatable :: big_Z(:,:,:), big_G2(:,:,:)


        sim%orbitals(1,1,1)%of(1)%R=sin(sim%grids(3)%R(:,:,:,1))*cos(sim%grids(3)%R(:,:,:,2))&
        *tan(sim%grids(3)%R(:,:,:,3))
        call real_to_recip(sim%orbitals(1,1,1)%of(1), sim%grids)
        do i=0,sim%parallel%nproc_band-1
            if(i.eq.sim%parallel%myid_band) then
                do j=1, sim%grids(3)%Ng_local(3)
                    print *, sim%parallel%myid_band, real(sim%orbitals(1,1,1)%of(1)%R(2,2,j)), sim%grids(3)%R(2,2,j,3), &
                     sim%orbitals(1,1,1)%of(1)%G(2,2,j), sim%grids(3)%G2(2,2,j)
                enddo
            endif
            flush(6)
            call parallel_wait(sim%parallel)
        enddo
    
        allocate(orbital_full(sim%grids(3)%Ng(1),sim%grids(3)%Ng(2),sim%grids(3)%Ng(3)))
        allocate(orbital_full_G(sim%grids(3)%Ng(1),sim%grids(3)%Ng(2),sim%grids(3)%Ng(3)))
        allocate(big_Z(sim%grids(3)%Ng(1),sim%grids(3)%Ng(2),sim%grids(3)%Ng(3)))
        allocate(big_G2(sim%grids(3)%Ng(1),sim%grids(3)%Ng(2),sim%grids(3)%Ng(3)))

        call collate_R(root=0, dsource=sim%orbitals(1,1,1)%of(1)%R, dtarget=orbital_full, &
                       comm='band', grid=sim%grids(3), parallel=sim%parallel)

        call collate_R(root=0, dsource=sim%orbitals(1,1,1)%of(1)%G, dtarget=orbital_full_G, &
                       comm='band', grid=sim%grids(3), parallel=sim%parallel)

        call collate_R(root=0, dsource=sim%grids(3)%R(:,:,:,3), dtarget=big_Z, &
                       comm='band', grid=sim%grids(3), parallel=sim%parallel)

        call collate_R(root=0, dsource=sim%grids(3)%G2(:,:,:), dtarget=big_G2, &
                       comm='band', grid=sim%grids(3), parallel=sim%parallel)

        if(sim%parallel%myid_band.eq.0) then
            print *, 'Collated'
            do j=1, sim%grids(3)%Ng(3)
                print *, sim%parallel%myid_band, real(orbital_full(2,2,j)), big_Z(2,2,j), &
                         orbital_full_G(2,2,j), big_G2(2,2,j)
            enddo
            orbital_full_G=orbital_full_G*2.0_dp
            orbital_full=orbital_full*2.0_dp
        endif
        flush(6)
        call parallel_wait(sim%parallel)       

        call distribute_R(root=0, dsource=orbital_full, dtarget=sim%orbitals(1,1,1)%of(1)%R, &
        comm='band', grid=sim%grids(3), parallel=sim%parallel)

        call distribute_R(root=0, dsource=orbital_full_G, dtarget=sim%orbitals(1,1,1)%of(1)%G, &
        comm='band', grid=sim%grids(3), parallel=sim%parallel)
        
        if(0.eq.sim%parallel%myid_band) print *, 'Distributed'
        call parallel_wait(sim%parallel)

        do i=0,sim%parallel%nproc_band-1
            if(i.eq.sim%parallel%myid_band) then
                do j=1, sim%grids(3)%Ng_local(3)
                    print *, sim%parallel%myid_band, real(sim%orbitals(1,1,1)%of(1)%R(2,2,j))/2.0_dp, sim%grids(3)%R(2,2,j,3), &
                    sim%orbitals(1,1,1)%of(1)%G(2,2,j)/2.0_dp, sim%grids(3)%G2(2,2,j)               
                enddo
            endif
            flush(6)
            call parallel_wait(sim%parallel)
        enddo
    

        stop
    end subroutine

    subroutine test_current(sim)
        use constants, only : i_
        use simulation_type, only : simulation_struct
        use fft, only: real_to_recip, recip_to_real
        use Apply_Hamiltonian, only : Apply_H
        use parallel_mod,only : parallel_wait
        use operations_3D, only : integrate_3D_R
        use Current, only : Current_Density_Calculation
        use grids_type,only : grid_struct

        type(simulation_struct), intent(inout), target :: sim
        type(grid_struct), pointer :: grid
        complex(dp), allocatable :: total_current(:,:)
        integer :: i, j, g

        sim%orbitals(1,1,1)%of(1)%R=sin(sim%grids(3)%R(:,:,:,1))*cos(sim%grids(3)%R(:,:,:,2))&
        *tan(sim%grids(3)%R(:,:,:,3)) + i_*cos(sim%grids(3)%R(:,:,:,1))*tan(sim%grids(3)%R(:,:,:,2))&
        *sin(sim%grids(3)%R(:,:,:,3))

        sim%potentials%total_local%of(1)%R=0.0_dp

        do g=3,size(sim%grids) !grids 1 and 2 belong to density which doesnt know about Vector Potential
            grid=>sim%grids(g)
            grid%G2(:,:,:) = 0.0_dp
            do i=1,1
                grid%G(:,:,:,i) = grid%G(:,:,:,i) + 0.01_dp !(G->G+A(t))
                grid%G2(:,:,:) = grid%G2(:,:,:) + abs(grid%G(:,:,:,i))**2
            enddo
            grid%A(1)=0.01_dp
        enddo

        call real_to_recip(sim%orbitals(1,1,1)%of(1), sim%grids)
        sim%orbitals(1,1,1)%of(1)%G=sim%orbitals(1,1,1)%of(1)%G*sim%grids(3)%cutwf
        call recip_to_real(sim%orbitals(1,1,1)%of(1), sim%grids)
        sim%orbitals(1,1,1)%occ(1)=1.0_dp
        sim%orbitals(1,1,1)%weight(1)=1.0_dp
        call Current_Density_Calculation(sim%orbitals, sim%parallel, sim%potentials, &
                 sim%atoms, sim%elements, sim%all_PAW, sim%grids, &
                sim%coarse_current_density, sim%current_density, sim%coarse_to_fine)
        allocate(total_current(3, sim%current_density(1)%n_s))
        do i=1,3;do j=1,sim%current_density(i)%n_s
            total_current(i,j)=real(integrate_3D_R(sim%current_density(i)%of(j)%R, &
                sim%grids(1), sim%parallel))/product(sim%grids(1)%box_length)
        enddo;enddo
        if(sim%parallel%myid.eq.0) print *, total_current(:,1)
        call parallel_wait(sim%parallel)
        do j = 0 ,sim%parallel%nproc-1
            do i=1, sim%grids(3)%Ng_local(3)
                if(sim%parallel%myid.eq.j) print *, sim%grids(3)%R(1,1,i,3), &
                sim%current_density(1)%of(1)%R(1,1,i), sim%current_density(2)%of(1)%R(1,1,i), &
                sim%current_density(3)%of(1)%R(1,1,i)
            enddo
            call parallel_wait(sim%parallel)
        enddo
        stop
    end subroutine

end module