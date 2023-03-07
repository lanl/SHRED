module grids_type
    use types, only : dp
    use, intrinsic :: iso_c_binding
    use odp_type, only : field_struct

    implicit none
    private

    public :: grid_struct, inter_grid_struct

    type grid_struct
    integer  :: Ng(3), Ng_local(3), Ng_small(3)
    integer :: myxyz(3)
    integer :: FFT_type, reduction_option
    integer :: FFT_comm, FFT_comm_yz, FFT_comm_xy
    integer :: nproc_FFT, nproc_FFT_yz, nproc_FFT_xy

    real(dp), allocatable :: G(:,:,:,:), R(:,:,:,:), G2(:,:,:), cutwf(:,:,:), cutden(:,:,:), &
         p(:,:,:), pprime(:,:,:), pdoubleprime(:,:,:)
    real(dp) :: dG(3), dR(3), Ecut, G2cut_den, G2cut_wf, k(3), A(3), ecutsm
    real(dp) :: Box_Length(3)
    type(field_struct) :: work ! a saved space for temporary work on this grid
    logical :: gamma

    !FFTW3 plans & local memory size 
    type(C_ptr) :: plan_forward, plan_reverse
    integer(C_INTPTR_T) :: alloc_local, alloc_local_R, alloc_local_G, alloc_local_R_gamma

    integer  :: Nr(3), Nr_local(3), myxyz_r(3)
    integer :: ncache, lotx, loty, lotz

    complex(C_DOUBLE_COMPLEX), pointer, contiguous :: zw(:)
    type(C_ptr) :: p_zw 

    complex(C_DOUBLE_COMPLEX), pointer, contiguous :: outXy(:,:), ybuff(:,:), &
                   outyX(:,:),outYzX(:,:,:), zbuff(:,:,:),outzY(:,:)

    type(C_ptr) :: p_Xy_out, p_ybuff, p_yX_out, p_YzX_out, p_zbuff, p_zY_out
    type(C_ptr) :: plan_fft_xy_x_lot, plan_fft_xy_x_left, plan_fft_yX_y_lot, plan_fft_yX_y_left, &
                   plan_fft_zY_z_lot, plan_fft_zY_z_left

    type(C_ptr) :: plan_Ifft_xy_x_lot, plan_Ifft_xy_x_left, plan_Ifft_yX_y_lot, plan_Ifft_yX_y_left, &
                   plan_Ifft_zY_z_lot, plan_Ifft_zY_z_left

    complex(C_DOUBLE_COMPLEX), pointer, contiguous :: buff(:,:), outyXz2(:,:,:), outYXz(:,:), outzYX(:,:)
    type(C_ptr) :: p_yXz2_out, p_zYX_out, p_YXz_out, p_buff
    type(C_ptr) :: plan_fft_xyz_x_lot, plan_fft_xyz_x_left, plan_fft_yXz_y_lot, plan_fft_yXz_y_left, &
                   plan_fft_zYX_z_lot, plan_fft_zYX_z_left

    type(C_ptr) :: plan_Ifft_xyz_x_lot, plan_Ifft_xyz_x_left, plan_Ifft_yXz_y_lot, plan_Ifft_yXz_y_left, &
                   plan_Ifft_zYX_z_lot, plan_Ifft_zYX_z_left
    end type

    type inter_grid_struct
        integer, allocatable :: pack_index(:,:,:)
        integer, allocatable :: count_(:,:)
        integer, allocatable :: sdisp(:)
        integer, allocatable :: rdisp(:)
        integer, allocatable :: unpack_index(:,:)
    end type

end module
