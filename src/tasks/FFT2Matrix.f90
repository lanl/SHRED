module FFT2Matrix
    use types, only : dp

    implicit none

    public :: orbitals_FFT2Matrix, orbitals_Matrix2FFT, &
    grids_G2_FFT2Matrix, orbital_FFT2Matrix, orbital_Matrix2FFT, &
    field_FFT2Matrix, field_Matrix2FFT, orbitals_G_mix_FFT2Matrix, orbitals_G_mix_Matrix2FFT

    interface full_to_small_1_2D
        module procedure full_to_small_1_2D_real
        module procedure full_to_small_1_2D_complex
    end interface

    interface full_to_small_11_2D
        module procedure full_to_small_11_2D_real
        module procedure full_to_small_11_2D_complex
    end interface

    interface full_to_small_12_1D
        module procedure full_to_small_12_1D_real
        module procedure full_to_small_12_1D_complex
    end interface

    contains

    subroutine orbitals_FFT2Matrix(orbitals, f_out ,grids, parallel, n_band_group, reduction_option)
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use parallel_mod,only : parallel_task
        use odp_type, only : orbital_struct

        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer :: grid

        type(orbital_struct), intent(inout) :: orbitals(:)
        integer, intent(in) :: n_band_group, reduction_option
        complex(dp), allocatable :: orbital_n(:,:,:)

        integer :: n, s,band_start, band_end, nG3, Ng_1(3)
        complex(dp), intent(inout), contiguous, target :: f_out(:,:)

        grid=>grids(orbitals(1)%of(1)%grid)
        Ng_1=grid%Ng_small
        !if(reduction_option.eq.2) then
        !    Ng_1(:)=grid%Ng_small
        !    if(grid%Ng_local(3).ne.grid%Ng_small(3)) then
        !         print *, 'Wrong reduction option on FFT2Matrix'; stop
        !    endif
        !else
        !    Ng_1(:)=grid%Ng_local
        !endif

        if(grid%FFT_type.gt.1) then
            Ng_1(1)=Ng_1(1)*maxval(orbitals(:)%n_spinor)
            Ng3=Ng_1(1)
            allocate(orbital_n(Ng_1(3),Ng_1(2),Ng_1(1)))
        else
            Ng_1(3)=Ng_1(3)*maxval(orbitals(:)%n_spinor)
            Ng3=Ng_1(3)
            allocate(orbital_n(Ng_1(1),Ng_1(2),Ng_1(3)))
        endif

        band_start=1
        band_end=0
        do n=1,size(orbitals)
            !It is required that all bands within a group have same n_spinor, this is enforced on setup
            band_end=band_end + n_band_group
            grid=>grids(orbitals(n)%of(1)%grid)
            do s=1,orbitals(n)%n_spinor
                if(reduction_option.eq.2) then
                    if(grid%FFT_type.eq.1) then 
                        call full_to_small_1_2D(orbitals(n)%of(s)%G(:,:,:),  &
                                          orbital_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, 1)
                    else if(grid%FFT_type.eq.11) then 
                        call full_to_small_11_2D(orbitals(n)%of(s)%G(:,:,:),  &
                                            orbital_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, 1)
                    else
                        print *, 'Error: reduction option not matching FFT_type', grid%FFT_type
                        stop
                    endif
                else
                    orbital_n(:,:,(1+nG3*(s-1)):(nG3*s))=orbitals(n)%of(s)%G(:,:,:)
                endif
            enddo
            call FFT2Matrix_complex(orbital_n(:,:,:), &
                    f_out(:,band_start:band_end), Ng_1, parallel, n_band_group)
            band_start=band_start+n_band_group
        enddo
        f_out=f_out*sqrt(product(grid%Box_length))
        
    end subroutine

    subroutine orbitals_G_mix_FFT2Matrix(orbitals, f_out ,grids, parallel, n_band_group, reduction_option, k)
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use parallel_mod,only : parallel_task
        use odp_type, only : orbital_struct

        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer :: grid

        type(orbital_struct), intent(inout) :: orbitals(:)
        integer, intent(in) :: n_band_group, reduction_option, k
        complex(dp), allocatable :: orbital_n(:,:,:)

        integer :: n, s,band_start, band_end, nG3, Ng_1(3)
        complex(dp), intent(inout), contiguous, target :: f_out(:,:)

        grid=>grids(orbitals(1)%of(1)%grid)

        Ng_1=grid%Ng_small
        if(grid%FFT_type.gt.1) then
            Ng_1(1)=Ng_1(1)*maxval(orbitals(:)%n_spinor)
            Ng3=Ng_1(1)
            allocate(orbital_n(Ng_1(3),Ng_1(2),Ng_1(1)))
        else
            Ng_1(3)=Ng_1(3)*maxval(orbitals(:)%n_spinor)
            Ng3=Ng_1(3)
            allocate(orbital_n(Ng_1(1),Ng_1(2),Ng_1(3)))
        endif

        band_start=1
        band_end=0
        do n=1,size(orbitals)
            !It is required that all bands within a group have same n_spinor, this is enforced on setup
            band_end=band_end + n_band_group
            grid=>grids(orbitals(n)%of(1)%grid)
            do s=1,orbitals(n)%n_spinor
                if(reduction_option.eq.2) then
                    call full_to_small_1_2D(orbitals(n)%of_G_mix(:,:,:,s,k),  &
                                          orbital_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, 1)
                    if(grid%FFT_type.eq.1) then 
                        call full_to_small_1_2D(orbitals(n)%of_G_mix(:,:,:,s,k),  &
                                          orbital_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, 1)
                    else if(grid%FFT_type.eq.11) then 
                        call full_to_small_11_2D(orbitals(n)%of_G_mix(:,:,:,s,k),  &
                                            orbital_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, 1)
                    else
                        print *, 'Error: reduction option not matching FFT_type', grid%FFT_type 
                        stop
                    endif
                else
                    orbital_n(:,:,(1+nG3*(s-1)):(nG3*s))=orbitals(n)%of_G_mix(:,:,:,s,k)
                endif
            enddo
            call FFT2Matrix_complex(orbital_n(:,:,:), &
                    f_out(:,band_start:band_end), Ng_1, parallel, n_band_group)
            band_start=band_start+n_band_group
        enddo
        f_out=f_out*sqrt(product(grid%Box_length))
        
    end subroutine

    subroutine orbital_FFT2Matrix(orbital, f_out ,grids, parallel,n_band_group, reduction_option)
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use parallel_mod,only : parallel_task
        use odp_type, only : orbital_struct

        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer :: grid

        type(orbital_struct), intent(inout) :: orbital
        integer, intent(in) :: n_band_group, reduction_option

        complex(dp), intent(inout), contiguous, target :: f_out(:,:)
        complex(dp), allocatable :: orbital_n(:,:,:)
        integer :: s,band_start, band_end, nG3, Ng_1(3)

        grid=>grids(orbital%of(1)%grid)

        Ng_1=grid%Ng_small
        if(grid%FFT_type.gt.1) then
            Ng_1(1)=Ng_1(1)*(orbital%n_spinor)
            Ng3=Ng_1(1)
            allocate(orbital_n(Ng_1(3),Ng_1(2),Ng_1(1)))
        else
            Ng_1(3)=Ng_1(3)*(orbital%n_spinor)
            Ng3=Ng_1(3)
            allocate(orbital_n(Ng_1(1),Ng_1(2),Ng_1(3)))
        endif

        band_start=1
        band_end=n_band_group
        !It is required that all bands within a group have same n_spinor, this is enforced on setup
        grid=>grids(orbital%of(1)%grid)
        do s=1,orbital%n_spinor
            if(reduction_option.eq.2) then
                if(grid%FFT_type.eq.1) then 
                    call full_to_small_1_2D(orbital%of(s)%G(:,:,:),  &
                                        orbital_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, 1)
                else if(grid%FFT_type.eq.11) then 
                    call full_to_small_11_2D(orbital%of(s)%G(:,:,:),  &
                                        orbital_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, 1)
                else
                    print *, 'Error: reduction option not matching FFT_type', grid%FFT_type 
                    stop
                endif
            else
                orbital_n(:,:,(1+nG3*(s-1)):(nG3*s))=orbital%of(s)%G(:,:,:)
            endif
        enddo
        call FFT2Matrix_complex(orbital_n(:,:,:), &
                f_out(:,band_start:band_end), Ng_1, parallel, n_band_group)

        f_out=f_out*sqrt(product(grid%Box_length))
    end subroutine

    subroutine grids_G2_FFT2Matrix(grid, f_out, parallel, n_band_groups, reduction_option)
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use parallel_mod,only : parallel_task

        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(inout)  :: grid
        real(dp), intent(inout), contiguous, target :: f_out(:,:)
        integer, intent(in) :: n_band_groups, reduction_option
        real(dp), allocatable :: G2_small(:,:,:), Gp(:,:,:)

        if(reduction_option.eq.2) then
            
            if(grid%FFT_type.gt.1) then 
                allocate(G2_small(grid%Ng_small(3),grid%Ng_small(2),grid%Ng_local(1)))
            else
                allocate(G2_small(grid%Ng_small(1),grid%Ng_small(2),grid%Ng_local(3)))
            endif

            if(grid%ecutsm.gt.0.0_dp) then
                allocate(Gp(grid%Ng_local(1),grid%Ng_local(2),grid%Ng_local(3)))
                Gp(:,:,:)=grid%G2*grid%p(:,:,:)
                if(grid%FFT_type.eq.1) then 
                    call full_to_small_1_2D(Gp(:,:,:), G2_small(:,:,:), grid, 1)
                else if(grid%FFT_type.eq.11) then 
                    call full_to_small_11_2D(Gp(:,:,:), G2_small(:,:,:), grid, 1)
                else
                    print *, 'Error: reduction option not matching FFT_type', grid%FFT_type 
                    stop
                endif
                deallocate(Gp)
            else
                if(grid%FFT_type.eq.1) then 
                    call full_to_small_1_2D(grid%G2(:,:,:), G2_small(:,:,:), grid, 1)
                else if(grid%FFT_type.eq.11) then 
                    call full_to_small_11_2D(grid%G2(:,:,:), G2_small(:,:,:), grid, 1)
                else
                    print *, 'Error: reduction option not matching FFT_type', grid%FFT_type 
                    stop
                endif
            endif

            call FFT2Matrix_real(G2_small(:,:,:), &
                    f_out(:,1:n_band_groups), grid%Ng_small, parallel, n_band_groups, 1.0_dp)
        else
            call FFT2Matrix_real(grid%G2(:,:,:), &
            f_out(:,1:n_band_groups), grid%Ng_local, parallel, n_band_groups, 1.0_dp)
        endif

    end subroutine

    subroutine field_FFT2Matrix(field, f_out, grid, parallel, n_band_group, reduction_option)
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use parallel_mod,only : parallel_task
        use odp_type, only : field_struct

        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(in)  :: grid
        type(field_struct), intent(inout)  :: field(:)

        complex(dp), intent(inout), contiguous, target :: f_out(:,:)
        integer, intent(in) :: n_band_group, reduction_option
        complex(dp), allocatable :: field_n(:,:,:)
        integer :: s, nG3, Ng_1(3)


        Ng_1=grid%Ng_small
        if(grid%FFT_type.gt.1) then
            Ng_1(1)=Ng_1(1)*size(field)
            Ng3=Ng_1(1)
            allocate(field_n(Ng_1(3),Ng_1(2),Ng_1(1)))
        else
            Ng_1(3)=Ng_1(3)*size(field)
            Ng3=Ng_1(3)
            allocate(field_n(Ng_1(1),Ng_1(2),Ng_1(3)))
        endif

        do s=1,size(field)
            if(reduction_option.eq.2) then
                if(grid%FFT_type.eq.1) then 
                    call full_to_small_1_2D(field(s)%G(:,:,:),  &
                                        field_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, 1)
                else if(grid%FFT_type.eq.11) then 
                    call full_to_small_11_2D(field(s)%G(:,:,:),  &
                                        field_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, 1)
                else
                    print *, 'Error: reduction option not matching FFT_type', grid%FFT_type 
                    stop
                endif
            else
                field_n(:,:,(1+nG3*(s-1)):(nG3*s))=field(s)%G(:,:,:)
            endif
        enddo
        
        call FFT2Matrix_complex(field_n(:,:,:), &
                f_out(:,1:n_band_group), Ng_1, parallel, n_band_group)
        f_out=f_out*sqrt(product(grid%Box_length))

        deallocate(field_n)
    end subroutine

    subroutine field_Matrix2FFT(field, f_in ,grid, parallel, n_band_group, reduction_option)
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use parallel_mod,only : parallel_task
        use odp_type, only : field_struct

        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(in):: grid

        type(field_struct), intent(inout)  :: field(:)
        integer, intent(in) :: n_band_group, reduction_option

        complex(dp), intent(in), contiguous, target :: f_in(:,:)
        complex(dp), allocatable :: field_n(:,:,:)
        integer :: s, nG3, Ng_1(3)


        Ng_1=grid%Ng_small
        if(grid%FFT_type.gt.1) then
            Ng_1(1)=Ng_1(1)*size(field)
            Ng3=Ng_1(1)
            allocate(field_n(Ng_1(3),Ng_1(2),Ng_1(1)))
        else
            Ng_1(3)=Ng_1(3)*size(field)
            Ng3=Ng_1(3)
            allocate(field_n(Ng_1(1),Ng_1(2),Ng_1(3)))
        endif

        call Matrix2FFT_complex( &
                f_in(:,1:n_band_group), &
                field_n(:,:,:), Ng_1, parallel, n_band_group)
                field_n(:,:,:)=field_n(:,:,:)/sqrt(product(grid%Box_length))
        do s=1,size(field)
            if(reduction_option.eq.2) then
                if(grid%FFT_type.eq.1) then 
                    call full_to_small_1_2D(field(s)%G(:,:,:), field_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, -1)
                else if(grid%FFT_type.eq.11) then 
                    call full_to_small_11_2D(field(s)%G(:,:,:), field_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, -1)
                else
                    print *, 'Error: reduction option not matching FFT_type', grid%FFT_type 
                    stop
                endif
            else
                field(s)%G(:,:,:)=field_n(:,:,(1+nG3*(s-1)):(nG3*s))
            endif
        enddo

        deallocate(field_n)

    end subroutine

    subroutine orbitals_Matrix2FFT(orbitals, f_in ,grids, parallel, n_band_group, reduction_option)
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use parallel_mod,only : parallel_task
        use odp_type, only : orbital_struct

        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer :: grid

        type(orbital_struct), intent(inout) :: orbitals(:)
        integer, intent(in) :: n_band_group, reduction_option

        integer :: n, s,band_start, band_end, nG3, Ng_1(3)
        complex(dp), intent(in), contiguous, target :: f_in(:,:)
        complex(dp), allocatable :: orbital_n(:,:,:)

        grid=>grids(orbitals(1)%of(1)%grid)

        Ng_1=grid%Ng_small
        if(grid%FFT_type.gt.1) then
            Ng_1(1)=Ng_1(1)*maxval(orbitals(:)%n_spinor)
            Ng3=Ng_1(1)
            allocate(orbital_n(Ng_1(3),Ng_1(2),Ng_1(1)))
        else
            Ng_1(3)=Ng_1(3)*maxval(orbitals(:)%n_spinor)
            Ng3=Ng_1(3)
            allocate(orbital_n(Ng_1(1),Ng_1(2),Ng_1(3)))
        endif

        band_start=1
        band_end=0
        do n=1,size(orbitals)
            !It is required that all bands within a group have same n_s, this is enforced on setup
            band_end=band_end + n_band_group
            grid=>grids(orbitals(n)%of(1)%grid)
            call Matrix2FFT_complex( &
                    f_in(:,band_start:band_end), &
                    orbital_n(:,:,:), Ng_1, parallel, n_band_group)
                    orbital_n(:,:,:)=orbital_n(:,:,:)/sqrt(product(grid%Box_length))

            do s=1, orbitals(n)%n_spinor
                if(reduction_option.eq.2) then
                    if(grid%FFT_type.eq.1) then 
                        call full_to_small_1_2D(orbitals(n)%of(s)%G(:,:,:), orbital_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, -1)
                    else if(grid%FFT_type.eq.11) then 
                        call full_to_small_11_2D(orbitals(n)%of(s)%G(:,:,:), orbital_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, -1)
                    else
                        print *, 'Error: reduction option not matching FFT_type', grid%FFT_type 
                        stop
                    endif
                else
                    orbitals(n)%of(s)%G(:,:,:)=orbital_n(:,:,(1+nG3*(s-1)):(nG3*s))
                endif
            enddo
            band_start=band_start+n_band_group
        enddo

    end subroutine

    subroutine orbitals_G_mix_Matrix2FFT(orbitals, f_in ,grids, parallel, n_band_group, reduction_option, k)
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use parallel_mod,only : parallel_task
        use odp_type, only : orbital_struct

        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer :: grid

        type(orbital_struct), intent(inout) :: orbitals(:)
        integer, intent(in) :: n_band_group, reduction_option, k

        integer :: n, s,band_start, band_end, nG3, Ng_1(3)
        complex(dp), intent(in), contiguous, target :: f_in(:,:)
        complex(dp), allocatable :: orbital_n(:,:,:)

        grid=>grids(orbitals(1)%of(1)%grid)

        Ng_1=grid%Ng_small
        if(grid%FFT_type.gt.1) then
            Ng_1(1)=Ng_1(1)*maxval(orbitals(:)%n_spinor)
            Ng3=Ng_1(1)
            allocate(orbital_n(Ng_1(3),Ng_1(2),Ng_1(1)))
        else
            Ng_1(3)=Ng_1(3)*maxval(orbitals(:)%n_spinor)
            Ng3=Ng_1(3)
            allocate(orbital_n(Ng_1(1),Ng_1(2),Ng_1(3)))
        endif

        band_start=1
        band_end=0
        do n=1,size(orbitals)
            !It is required that all bands within a group have same n_s, this is enforced on setup
            band_end=band_end + n_band_group
            grid=>grids(orbitals(n)%of(1)%grid)
            call Matrix2FFT_complex( &
                    f_in(:,band_start:band_end), &
                    orbital_n(:,:,:), Ng_1, parallel, n_band_group)
                    orbital_n(:,:,:)=orbital_n(:,:,:)/sqrt(product(grid%Box_length))

            do s=1, orbitals(n)%n_spinor
                if(reduction_option.eq.2) then
                    if(grid%FFT_type.eq.1) then 
                        call full_to_small_1_2D(orbitals(n)%of_G_mix(:,:,:,s,k), orbital_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, -1)
                    else if(grid%FFT_type.eq.11) then 
                        call full_to_small_11_2D(orbitals(n)%of_G_mix(:,:,:,s,k), orbital_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, -1)
                    else
                        print *, 'Error: reduction option not matching FFT_type', grid%FFT_type 
                        stop
                    endif
                else
                    orbitals(n)%of_G_mix(:,:,:,s,k)=orbital_n(:,:,(1+nG3*(s-1)):(nG3*s))
                endif
            enddo
            band_start=band_start+n_band_group
        enddo

    end subroutine

    subroutine orbital_Matrix2FFT(orbital, f_in ,grids, parallel, n_band_group, reduction_option)
        use parallel_type, only : parallel_struct
        use grids_type, only : grid_struct
        use parallel_mod,only : parallel_task
        use odp_type, only : orbital_struct

        type(parallel_struct), intent(in) :: parallel
        type(grid_struct), intent(in), target :: grids(:)
        type(grid_struct), pointer :: grid

        type(orbital_struct), intent(inout) :: orbital
        integer, intent(in) :: n_band_group, reduction_option

        complex(dp), intent(in), contiguous, target :: f_in(:,:)
        complex(dp), allocatable :: orbital_n(:,:,:)
        integer :: s,band_start, band_end, nG3, Ng_1(3)

        grid=>grids(orbital%of(1)%grid)

        Ng_1=grid%Ng_small
        if(grid%FFT_type.gt.1) then
            Ng_1(1)=Ng_1(1)*(orbital%n_spinor)
            Ng3=Ng_1(1)
            allocate(orbital_n(Ng_1(3),Ng_1(2),Ng_1(1)))
        else
            Ng_1(3)=Ng_1(3)*(orbital%n_spinor)
            Ng3=Ng_1(3)
            allocate(orbital_n(Ng_1(1),Ng_1(2),Ng_1(3)))
        endif

        band_start=1
        band_end=0
        !It is required that all bands within a group have same n_s, this is enforced on setup
        band_end=band_end + n_band_group
        grid=>grids(orbital%of(1)%grid)
        call Matrix2FFT_complex( &
                f_in(:,1:n_band_group), &
                orbital_n(:,:,:), Ng_1, parallel, n_band_group)
                orbital_n(:,:,:)=orbital_n(:,:,:)/sqrt(product(grid%Box_length))

        do s=1,orbital%n_spinor
            if(reduction_option.eq.2) then
                if(grid%FFT_type.eq.1) then 
                    call full_to_small_1_2D(orbital%of(s)%G(:,:,:), orbital_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, -1)
                else if(grid%FFT_type.eq.11) then 
                    call full_to_small_11_2D(orbital%of(s)%G(:,:,:), orbital_n(:,:,(1+nG3*(s-1)):(nG3*s)), grid, -1)
                else
                    print *, 'Error: reduction option not matching FFT_type', grid%FFT_type 
                    stop
                endif
            else
                orbital%of(s)%G(:,:,:)=orbital_n(:,:,(1+nG3*(s-1)):(nG3*s))
            endif
        enddo
        band_start=band_start+n_band_group

    end subroutine

    subroutine FFT2Matrix_complex(f_in, f_out, Ng_local, parallel, n_band_group)
        use parallel_type, only : parallel_struct
        use parallel_mod,only : parallel_task

            type(parallel_struct), intent(in) :: parallel

            complex(dp), intent(in), contiguous, target :: f_in(:,:,:)
            complex(dp), target, allocatable :: f_tmp(:)
            integer, intent(in) :: n_band_group, Ng_local(3)

            complex(dp), pointer :: buffer1(:), buffer2(:)
            complex(dp), intent(inout), contiguous, target :: f_out(:,:)
            integer :: buffered_grid_size

            buffered_grid_size=product(Ng_local)
            buffer2(1:product(Ng_local))=>f_in(:,:,:)
            
            if(mod(buffered_grid_size,n_band_group).ne.0) then
                buffered_grid_size=(product(Ng_local)/n_band_group)*n_band_group &
                                    + n_band_group
                allocate(f_tmp(buffered_grid_size))
                f_tmp=0.0_dp
                f_tmp(1:product(Ng_local))=buffer2(1:product(Ng_local))
                buffer2(1:buffered_grid_size)=>f_tmp(:)
            endif
            
            if(size(f_out).ne.buffered_grid_size) then
                print *, 'FFT2Matrix Error: Vectorized Array does not contain the required buffer space,'
                print *, size(f_out), buffered_grid_size
                stop
            endif
            buffer1(1:size(f_out))=>f_out(:,:)

            if(n_band_group.gt.1) then
                call parallel_task('alltoall', buffer2(:), parallel, 'diff_b', buffer1(:),  &
                send_count=buffered_grid_size/n_band_group, recv_count=buffered_grid_size/n_band_group)
            else
                buffer1(:)=buffer2(:)
            endif

            if(allocated(f_tmp)) deallocate(f_tmp)
    end subroutine

    subroutine FFT2Matrix_real(f_in, f_out, Ng_local, parallel, n_band_group, set_val)
        use parallel_type, only : parallel_struct
        use parallel_mod,only : parallel_task

            type(parallel_struct), intent(in) :: parallel

            real(dp), intent(in), contiguous, target :: f_in(:,:,:)
            real(dp), target, allocatable :: f_tmp(:)
            integer, intent(in) :: n_band_group, Ng_local(3)
            real(dp), intent(in), optional :: set_val
            real(dp)  :: val

            real(dp), pointer :: buffer1(:), buffer2(:)
            real(dp), intent(inout), contiguous, target :: f_out(:,:)
            integer :: buffered_grid_size


            buffered_grid_size=product(Ng_local)
            buffer2(1:product(Ng_local))=>f_in(:,:,:)
            val=0.0_dp
            if(present(set_val)) val=set_val
            
            if(mod(buffered_grid_size,n_band_group).ne.0) then
                buffered_grid_size=(product(Ng_local)/n_band_group)*n_band_group &
                                    + n_band_group
                allocate(f_tmp(buffered_grid_size))
                f_tmp=val
                f_tmp(1:product(Ng_local))=buffer2(1:product(Ng_local))
                buffer2(1:buffered_grid_size)=>f_tmp(:)
            endif

            if(size(f_out).ne.buffered_grid_size) then
                print *, 'FFT2Matrix Error: Vectorized Array does not contain the required buffer space,'
                print *, size(f_out), buffered_grid_size
                stop
            endif
            buffer1(1:size(f_out))=>f_out(:,:)
            
            if(n_band_group.gt.1) then
                call parallel_task('alltoall', buffer2(:), parallel, 'diff_b', buffer1(:),  &
                send_count=buffered_grid_size/n_band_group, recv_count=buffered_grid_size/n_band_group)
            else
                buffer1(:)=buffer2(:)
            endif
            
            if(allocated(f_tmp)) deallocate(f_tmp)
    end subroutine

    subroutine Matrix2FFT_complex(f_in, f_out, Ng_local, parallel, n_band_groups)
        use parallel_type, only : parallel_struct
        use parallel_mod,only : parallel_task

            type(parallel_struct), intent(in) :: parallel

            complex(dp), intent(in), contiguous, target :: f_in(:,:)
            complex(dp), intent(inout), contiguous, target :: f_out(:,:,:)
            integer, intent(in) :: n_band_groups, Ng_local(3)

            complex(dp), target, allocatable :: f_tmp(:)

            complex(dp), pointer :: buffer1(:), buffer2(:), buffer3(:)
            integer :: buffered_grid_size

            buffered_grid_size=product(Ng_local)
            buffer3(1:product(Ng_local))=>f_out(:,:,:)

            buffer1(1:size(f_in))=>f_in(:,:)

            if(mod(buffered_grid_size,n_band_groups).ne.0) then
                buffered_grid_size=(product(Ng_local)/n_band_groups)*n_band_groups &
                                    + n_band_groups
                allocate(f_tmp(buffered_grid_size))
                buffer2(1:buffered_grid_size)=>f_tmp(:)
            else
                buffer2(1:product(Ng_local))=>f_out(:,:,:)
            endif

            if(size(f_in).ne.buffered_grid_size) then
                print *, 'Matrix2FFT Error: Vectorized Array does not contain the required buffer space,'
                print *, size(f_in), buffered_grid_size
                stop
            endif

            if(n_band_groups.gt.1) then
                call parallel_task('alltoall', buffer1(:), parallel, 'diff_b', buffer2(:),  &
                send_count=buffered_grid_size/n_band_groups, recv_count=buffered_grid_size/n_band_groups)
            else
                buffer2(:)=buffer1(:)
            endif


            if(buffered_grid_size.ne.product(Ng_local)) then
                buffer3(1:product(Ng_local))=buffer2(1:product(Ng_local))
            endif

            if(allocated(f_tmp)) deallocate(f_tmp)
    end subroutine

    subroutine full_to_small_1_2D_real(f_in, f_in_small, grid, forward_or_back)
        use grids_type, only : grid_struct
        type(grid_struct), intent(in) :: grid
        real(dp), intent(inout)  :: f_in(:,:,:)
        real(dp), intent(inout)  :: f_in_small(:,:,:)
        integer, intent(in):: forward_or_back
        integer :: n_positive(2), n_negative(2)
        integer :: is1(2), ie1(2)
        integer :: is2(2), ie2(2)
        integer :: i, j

        n_positive(:)=grid%Ng_small(:2)/2 + 1
        n_negative(:)=grid%Ng_small(:2)/2 + mod(grid%Ng_small(:2),2)-1
        if(forward_or_back.eq.-1) f_in=0.0
        do i=1,2 !1=> positive G's 2=> negative G's
            if(i.eq.1) then
                is1(1)=1
                ie1(1)=n_positive(1)
                is2(1)=is1(1)
                ie2(1)=ie1(1)
            else
                is1(1)=grid%Ng_small(1)-n_negative(1)+1
                is2(1)=grid%Ng_local(1)-n_negative(1)+1
                ie1(1)=grid%Ng_small(1)
                ie2(1)=grid%Ng_local(1)
            endif
            do j=1,2 !1=> positive G's 2=> negative G's
                if(j.eq.1) then
                    is1(2)=1
                    ie1(2)=n_positive(2)
                    is2(2)=is1(2)
                    ie2(2)=ie1(2)
                else
                    is1(2)=grid%Ng_small(2)-n_negative(2)+1
                    is2(2)=grid%Ng_local(2)-n_negative(2)+1
                    ie1(2)=grid%Ng_small(2)
                    ie2(2)=grid%Ng_local(2)
                endif
                if(forward_or_back.eq.1)  then
                    f_in_small(is1(1):ie1(1),is1(2):ie1(2),:) = &
                        f_in(is2(1):ie2(1),is2(2):ie2(2),:)
                else if( forward_or_back.eq.-1 ) then
                    f_in(is2(1):ie2(1),is2(2):ie2(2),:) = &
                f_in_small(is1(1):ie1(1),is1(2):ie1(2),:) 
                else
                    print *, 'forward or back transform of full grid to small grid not specified'
                    stop
                endif
            enddo
        enddo
    end subroutine

    subroutine full_to_small_11_2D_real(f_in, f_in_small, grid, forward_or_back)
        use grids_type, only : grid_struct
        type(grid_struct), intent(in) :: grid
        real(dp), intent(inout)  :: f_in(:,:,:)
        real(dp), intent(inout)  :: f_in_small(:,:,:)
        integer, intent(in):: forward_or_back
        integer :: n_positive(2), n_negative(2)
        integer :: is1(2), ie1(2)
        integer :: is2(2), ie2(2)
        integer :: i, j

        n_positive(1)=grid%Ng_small(3)/2 + 1
        n_negative(1)=grid%Ng_small(3)/2 + mod(grid%Ng_small(3),2)-1

        n_positive(2)=grid%Ng_small(2)/2 + 1
        n_negative(2)=grid%Ng_small(2)/2 + mod(grid%Ng_small(2),2)-1

        if(forward_or_back.eq.-1) f_in=0.0
        do i=1,2 !1=> positive G's 2=> negative G's
            if(i.eq.1) then
                is1(1)=1
                ie1(1)=n_positive(1)
                is2(1)=is1(1)
                ie2(1)=ie1(1)
            else
                is1(1)=grid%Ng_small(3)-n_negative(1)+1
                is2(1)=grid%Ng_local(3)-n_negative(1)+1
                ie1(1)=grid%Ng_small(3)
                ie2(1)=grid%Ng_local(3)
            endif
            do j=1,2 !1=> positive G's 2=> negative G's
                if(j.eq.1) then
                    is1(2)=1
                    ie1(2)=n_positive(2)
                    is2(2)=is1(2)
                    ie2(2)=ie1(2)
                else
                    is1(2)=grid%Ng_small(2)-n_negative(2)+1
                    is2(2)=grid%Ng_local(2)-n_negative(2)+1
                    ie1(2)=grid%Ng_small(2)
                    ie2(2)=grid%Ng_local(2)
                endif
                if(forward_or_back.eq.1)  then
                    f_in_small(is1(1):ie1(1),is1(2):ie1(2),:) = &
                        f_in(is2(1):ie2(1),is2(2):ie2(2),:)
                else if( forward_or_back.eq.-1 ) then
                    f_in(is2(1):ie2(1),is2(2):ie2(2),:) = &
                f_in_small(is1(1):ie1(1),is1(2):ie1(2),:) 
                else
                    print *, 'forward or back transform of full grid to small grid not specified'
                    stop
                endif
            enddo
        enddo
    end subroutine

    subroutine full_to_small_12_1D_real(f_in, f_in_small, grid, forward_or_back)
        use grids_type, only : grid_struct
        type(grid_struct), intent(in) :: grid
        real(dp), intent(inout)  :: f_in(:,:,:)
        real(dp), intent(inout)  :: f_in_small(:,:,:)
        integer, intent(in):: forward_or_back
        integer :: n_positive, n_negative
        integer :: is1, ie1
        integer :: is2, ie2
        integer :: i

        n_positive=grid%Ng_small(3)/2 + 1
        n_negative=grid%Ng_small(3)/2 + mod(grid%Ng_small(3),2)-1


        if(forward_or_back.eq.-1) f_in=0.0
        do i=1,2 !1=> positive G's 2=> negative G's
            if(i.eq.1) then
                is1=1
                ie1=n_positive
                is2=is1
                ie2=ie1
            else
                is1=grid%Ng_small(3)-n_negative+1
                is2=grid%Ng_local(3)-n_negative+1
                ie1=grid%Ng_small(3)
                ie2=grid%Ng_local(3)
            endif

            if(forward_or_back.eq.1)  then
                f_in_small(is1:ie1,:,:) = &
                    f_in(is2:ie2,:,:)
            else if( forward_or_back.eq.-1 ) then
                f_in(is2:ie2,:,:) = &
                    f_in_small(is1:ie1,:,:) 
            else
                print *, 'forward or back transform of full grid to small grid not specified'
                stop
            endif
        enddo
    end subroutine

    subroutine full_to_small_1_2D_complex(f_in, f_in_small, grid, forward_or_back)
        use grids_type, only : grid_struct
        type(grid_struct), intent(in) :: grid
        complex(dp), intent(inout)  :: f_in(:,:,:)
        complex(dp), intent(inout):: f_in_small(:,:,:)
        integer, intent(in):: forward_or_back
        integer :: n_positive(2), n_negative(2)
        integer :: is1(2), ie1(2)
        integer :: is2(2), ie2(2)
        integer :: i, j

        n_positive(:)=grid%Ng_small(:2)/2 + 1
        n_negative(:)=grid%Ng_small(:2)/2 + mod(grid%Ng_small(:2),2)-1
        if(forward_or_back.eq.-1) f_in=0.0
        do i=1,2 !1=> positive G's 2=> negative G's
            if(i.eq.1) then
                is1(1)=1
                ie1(1)=n_positive(1)
                is2(1)=is1(1)
                ie2(1)=ie1(1)
            else
                is1(1)=grid%Ng_small(1)-n_negative(1)+1
                is2(1)=grid%Ng_local(1)-n_negative(1)+1
                ie1(1)=grid%Ng_small(1)
                ie2(1)=grid%Ng_local(1)
            endif
            do j=1,2 !1=> positive G's 2=> negative G's
                if(j.eq.1) then
                    is1(2)=1
                    ie1(2)=n_positive(2)
                    is2(2)=is1(2)
                    ie2(2)=ie1(2)
                else
                    is1(2)=grid%Ng_small(2)-n_negative(2)+1
                    is2(2)=grid%Ng_local(2)-n_negative(2)+1
                    ie1(2)=grid%Ng_small(2)
                    ie2(2)=grid%Ng_local(2)
                endif

                if(forward_or_back.eq.1)  then
                    f_in_small(is1(1):ie1(1),is1(2):ie1(2),:) = &
                        f_in(is2(1):ie2(1),is2(2):ie2(2),:)
                else if( forward_or_back.eq.-1 ) then
                    f_in(is2(1):ie2(1),is2(2):ie2(2),:) = &
                f_in_small(is1(1):ie1(1),is1(2):ie1(2),:) 
                else
                    print *, 'forward or back transform of full grid to small grid not specified'
                    stop
                endif
            enddo
        enddo
    end subroutine

    subroutine full_to_small_11_2D_complex(f_in, f_in_small, grid, forward_or_back)
        use grids_type, only : grid_struct
        type(grid_struct), intent(in) :: grid
        complex(dp), intent(inout)  :: f_in(:,:,:)
        complex(dp), intent(inout)  :: f_in_small(:,:,:)
        integer, intent(in):: forward_or_back
        integer :: n_positive(2), n_negative(2)
        integer :: is1(2), ie1(2)
        integer :: is2(2), ie2(2)
        integer :: i, j

        n_positive(1)=grid%Ng_small(3)/2 + 1
        n_negative(1)=grid%Ng_small(3)/2 + mod(grid%Ng_small(3),2)-1

        n_positive(2)=grid%Ng_small(2)/2 + 1
        n_negative(2)=grid%Ng_small(2)/2 + mod(grid%Ng_small(2),2)-1

        if(forward_or_back.eq.-1) f_in=0.0
        do i=1,2 !1=> positive G's 2=> negative G's
            if(i.eq.1) then
                is1(1)=1
                ie1(1)=n_positive(1)
                is2(1)=is1(1)
                ie2(1)=ie1(1)
            else
                is1(1)=grid%Ng_small(3)-n_negative(1)+1
                is2(1)=grid%Ng_local(3)-n_negative(1)+1
                ie1(1)=grid%Ng_small(3)
                ie2(1)=grid%Ng_local(3)
            endif
            do j=1,2 !1=> positive G's 2=> negative G's
                if(j.eq.1) then
                    is1(2)=1
                    ie1(2)=n_positive(2)
                    is2(2)=is1(2)
                    ie2(2)=ie1(2)
                else
                    is1(2)=grid%Ng_small(2)-n_negative(2)+1
                    is2(2)=grid%Ng_local(2)-n_negative(2)+1
                    ie1(2)=grid%Ng_small(2)
                    ie2(2)=grid%Ng_local(2)
                endif
                if(forward_or_back.eq.1)  then
                    f_in_small(is1(1):ie1(1),is1(2):ie1(2),:) = &
                        f_in(is2(1):ie2(1),is2(2):ie2(2),:)
                else if( forward_or_back.eq.-1 ) then
                    f_in(is2(1):ie2(1),is2(2):ie2(2),:) = &
                f_in_small(is1(1):ie1(1),is1(2):ie1(2),:) 
                else
                    print *, 'forward or back transform of full grid to small grid not specified'
                    stop
                endif
            enddo
        enddo
    end subroutine

    subroutine full_to_small_12_1D_complex(f_in, f_in_small, grid, forward_or_back)
        use grids_type, only : grid_struct
        type(grid_struct), intent(in) :: grid
        complex(dp), intent(inout)  :: f_in(:,:,:)
        complex(dp), intent(inout)  :: f_in_small(:,:,:)
        integer, intent(in):: forward_or_back
        integer :: n_positive, n_negative
        integer :: is1, ie1
        integer :: is2, ie2
        integer :: i

        n_positive=grid%Ng_small(3)/2 + 1
        n_negative=grid%Ng_small(3)/2 + mod(grid%Ng_small(3),2)-1


        if(forward_or_back.eq.-1) f_in=0.0
        do i=1,2 !1=> positive G's 2=> negative G's
            if(i.eq.1) then
                is1=1
                ie1=n_positive
                is2=is1
                ie2=ie1
            else
                is1=grid%Ng_small(3)-n_negative+1
                is2=grid%Ng_local(3)-n_negative+1
                ie1=grid%Ng_small(3)
                ie2=grid%Ng_local(3)
            endif

            if(forward_or_back.eq.1)  then
                f_in_small(is1:ie1,:,:) = &
                    f_in(is2:ie2,:,:)
            else if( forward_or_back.eq.-1 ) then
                f_in(is2:ie2,:,:) = &
                    f_in_small(is1:ie1,:,:) 
            else
                print *, 'forward or back transform of full grid to small grid not specified'
                stop
            endif
        enddo
    end subroutine

end module