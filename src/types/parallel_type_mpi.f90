module parallel_type
    implicit none
    private
    public :: parallel_struct

type parallel_struct
    integer :: myid, nproc, nproc_k, nproc_band, nproc_spin, nproc_diff_b, nproc_diff_k, nproc_diff_s, nproc_space
    
    !Spin is division of all procs by the spin groups (max 2 at the moment)
    !diff_s will have the same k, bands, and space but for the different spins
    integer :: my_spin_group, n_spin_groups, myid_spin
    integer :: my_diff_s_group, n_diff_s_groups, myid_diff_s

    !k is division of the spin group into groups of k_points
    !diff k will have the same bands, spin, and space but for the different k_points
    integer :: my_k_group, n_k_groups, myid_k
    integer :: my_diff_k_group, n_diff_k_groups, myid_diff_k

    !band is division of the k group into groups of bands
    !diff b will have the same k, spin, and space but for the different bands
    integer :: my_band_group, n_band_groups, myid_band
    integer :: my_diff_b_group, n_diff_b_groups, myid_diff_b

    !The grids are spatially decomposed in the FFT setup further dividing the band_gorup into local grids
    !with the same bands, ks, and spin(s)

    !The space group is special in that it is not a sub-divison like band or k, it is a division of the all communicator by myid_band
    !Thus it includes all processors which have different bands, k, and spin groups, but share the same local real grid.
    !It is good for summing up real-space density contributions or dividing pseudopotential overhead work
    integer :: my_space_group, myid_space

    integer :: n_band_groups_sub, myid_sub, nproc_sub, my_sub
    integer :: comm_diff_b_sub, myid_diff_b_sub, nproc_diff_b_sub

    
    !2D FFT parallelism 
    integer:: nproc_fft_xy, nproc_fft_yz, myid_fft_xy, myid_fft_yz, my_fft_xy, my_fft_yz
    integer :: comm_fft_yz, comm_fft_xy

    !processes on same node for shared memeory
    integer:: myid_node, nproc_node
    integer:: myid_diff_node, nproc_diff_node
    integer:: nnodes, my_node, my_node_ordered_id, my_band_node
    integer, allocatable :: nproc_all_node(:)

    integer :: comm_all, comm_bands, comm_k, comm_spin, comm_space, comm_node
    integer :: comm_diff_b, comm_diff_k, comm_diff_s, comm_diff_node, comm_sub

    integer ::  myxyz(3), ierr

    integer :: my_atom_max, my_atom_min, my_natom

    logical :: by_node
end type

end module