module operations_3D
    use types, only: dp
    implicit none

    public :: integrate_3D_R, integrate_3D_G

    interface integrate_3D_R
        module procedure integrate_3D_R_complex
        module procedure integrate_3D_R_real
    end interface
    
    interface integrate_3D_G
        module procedure integrate_3D_G_complex
        module procedure integrate_3D_G_real
    end interface
contains
    function integrate_3D_R_complex(R, grid, parallel) result(integral)
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use parallel_mod, only: parallel_task

        complex(dp), intent(in) :: R(:,:,:)
        type(grid_struct), intent(in) :: grid
        type(parallel_struct), intent(in) :: parallel


        complex(dp) :: integral

        integral=sum(R(:,:,:))*product(grid%dR(:))
        call parallel_task('sum', integral, parallel, 'band')
    end function

    function integrate_3D_R_real(R, grid, parallel) result(integral)
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use parallel_mod, only: parallel_task

        real(dp), intent(in) :: R(:,:,:)
        type(grid_struct), intent(in) :: grid
        type(parallel_struct), intent(in) :: parallel


        real(dp) :: integral

        integral=sum(R(:,:,:))*product(grid%dR(:))

        call parallel_task('sum', integral, parallel, 'band')

    end function

    function integrate_3D_G_complex(G, grid, parallel) result(integral)
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use parallel_mod, only: parallel_task, parallel_wait
        use constants, only: pi
        
        complex(dp), intent(in) :: G(:,:,:)
        type(grid_struct), intent(in) :: grid
        type(parallel_struct), intent(in) :: parallel


        complex(dp) :: integral

        integral=sum(G(:,:,:))*product(grid%box_length(:))

        call parallel_task('sum', integral, parallel, 'band')

        if(grid%gamma) integral=real(integral)

    end function

    function integrate_3D_G_real(G, grid, parallel) result(integral)
        use grids_type, only: grid_struct
        use parallel_type, only: parallel_struct
        use parallel_mod, only: parallel_task, parallel_wait
        use constants, only: pi
        
        real(dp), intent(in) :: G(:,:,:)
        type(grid_struct), intent(in) :: grid
        type(parallel_struct), intent(in) :: parallel


        real(dp) :: integral

        integral=sum(G(:,:,:))*product(grid%box_length(:))

        call parallel_task('sum', integral, parallel, 'band')

    end function

end module