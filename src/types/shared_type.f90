module shared_type
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
    use types, only : dp
    use mpi

    implicit none
    
    private

    public  :: shared_struct, allocate_shared, shared_fence

type shared_struct
    integer, allocatable :: arrayshape(:)
    INTEGER(KIND=MPI_ADDRESS_KIND)  :: window_size
    TYPE(C_PTR) ::  baseptr
    integer :: win
end type

interface  allocate_shared
    module procedure allocate_shared_doubles
    module procedure allocate_shared_complexs
    module procedure allocate_shared_doubles2
    module procedure allocate_shared_complexs2
    module procedure allocate_shared_doubles3
    module procedure allocate_shared_complexs3
end interface

contains


subroutine allocate_shared_doubles(array, arrayshape, shared, parallel)
    use parallel_type, only: parallel_struct
    type(shared_struct), intent(inout) :: shared
    type(parallel_struct), intent(in) :: parallel
    integer,intent(in) :: arrayshape(:)
    real(dp), pointer, intent(inout) :: array(:)
    integer(kind=MPI_COUNT_KIND) :: byte_count
    integer :: ierr
    call MPI_Type_size_x(MPI_DOUBLE_PRECISION, byte_count, ierr)

    if(parallel%myid_node.eq.0) then
        shared%window_size= int(product(arrayshape),MPI_ADDRESS_KIND)*byte_count
    else
        shared%window_size= 0
    endif

    call shared_allocate_common(shared, parallel%myid_node, parallel%comm_node, byte_count)

    CALL C_F_POINTER(shared%baseptr, array , arrayshape)


end subroutine

subroutine allocate_shared_complexs(array, arrayshape, shared, parallel)
    use parallel_type, only: parallel_struct
    type(shared_struct), intent(inout) :: shared
    type(parallel_struct), intent(in) :: parallel
    integer,intent(in) :: arrayshape(:)
    complex(dp), pointer, intent(inout) :: array(:)
    integer(kind=MPI_COUNT_KIND) :: byte_count
    integer :: ierr
    call MPI_Type_size_x(MPI_DOUBLE_COMPLEX, byte_count, ierr)

    if(parallel%myid_node.eq.0) then
        shared%window_size= int(product(arrayshape),MPI_ADDRESS_KIND)*byte_count
    else
        shared%window_size= 0
    endif

    call shared_allocate_common(shared, parallel%myid_node, parallel%comm_node, byte_count)

    CALL C_F_POINTER(shared%baseptr, array , arrayshape)


end subroutine

subroutine allocate_shared_doubles2(array, arrayshape, shared, parallel)
    use parallel_type, only: parallel_struct
    type(shared_struct), intent(inout) :: shared
    type(parallel_struct), intent(in) :: parallel
    integer,intent(in) :: arrayshape(:)
    real(dp), pointer, intent(inout) :: array(:,:)
    integer(kind=MPI_COUNT_KIND) :: byte_count
    integer :: ierr
    call MPI_Type_size_x(MPI_DOUBLE_PRECISION, byte_count, ierr)

    if(parallel%myid_node.eq.0) then
        shared%window_size= int(product(arrayshape),MPI_ADDRESS_KIND)*byte_count
    else
        shared%window_size= 0
    endif

    call shared_allocate_common(shared, parallel%myid_node, parallel%comm_node, byte_count)

    CALL C_F_POINTER(shared%baseptr, array , arrayshape)

end subroutine

subroutine allocate_shared_complexs2(array, arrayshape, shared, parallel)
    use parallel_type, only: parallel_struct
    type(shared_struct), intent(inout) :: shared
    type(parallel_struct), intent(in) :: parallel
    integer,intent(in) :: arrayshape(:)
    complex(dp), pointer, intent(inout) :: array(:,:)
    integer(kind=MPI_COUNT_KIND) :: byte_count
    integer :: ierr

    call MPI_Type_size_x(MPI_DOUBLE_COMPLEX, byte_count, ierr)

    if(parallel%myid_node.eq.0) then
        shared%window_size= int(product(arrayshape),MPI_ADDRESS_KIND)*byte_count
    else
        shared%window_size= 0
    endif

    call shared_allocate_common(shared, parallel%myid_node, parallel%comm_node, byte_count)

    CALL C_F_POINTER(shared%baseptr, array , arrayshape)

end subroutine

subroutine allocate_shared_doubles3(array, arrayshape, shared, parallel)
    use parallel_type, only: parallel_struct
    type(shared_struct), intent(inout) :: shared
    type(parallel_struct), intent(in) :: parallel
    integer,intent(in) :: arrayshape(:)
    real(dp), pointer, intent(inout) :: array(:,:,:)
    integer(kind=MPI_COUNT_KIND) :: byte_count
    integer :: ierr
    call MPI_Type_size_x(MPI_DOUBLE_COMPLEX, byte_count, ierr)

    if(parallel%myid_node.eq.0) then
        shared%window_size= int(product(arrayshape),MPI_ADDRESS_KIND)*byte_count
    else
        shared%window_size= 0
    endif

    call shared_allocate_common(shared, parallel%myid_node, parallel%comm_node, byte_count)

    CALL C_F_POINTER(shared%baseptr, array , arrayshape)


end subroutine

subroutine allocate_shared_complexs3(array, arrayshape, shared, parallel)
    use parallel_type, only: parallel_struct
    type(shared_struct), intent(inout) :: shared
    type(parallel_struct), intent(in) :: parallel
    integer,intent(in) :: arrayshape(:)
    complex(dp), pointer, intent(inout) :: array(:,:,:)
    integer(kind=MPI_COUNT_KIND) :: byte_count
    integer :: ierr
    call MPI_Type_size_x(MPI_DOUBLE_COMPLEX, byte_count, ierr)

    if(parallel%myid_node.eq.0) then
        shared%window_size= int(product(arrayshape),MPI_ADDRESS_KIND)*byte_count
    else
        shared%window_size= 0
    endif

    call shared_allocate_common(shared, parallel%myid_node, parallel%comm_node, byte_count)

    CALL C_F_POINTER(shared%baseptr, array , arrayshape)


end subroutine

subroutine shared_allocate_common(shared, myid, comm, byte_count)
    type(shared_struct), intent(inout) :: shared
    integer,intent(in) :: myid, comm
    integer(kind=MPI_COUNT_KIND), intent(in) :: byte_count

    integer :: ierr, target_rank, disp_unit

    disp_unit=int(byte_count)

    CALL MPI_Win_allocate_shared(shared%window_size, disp_unit, MPI_INFO_NULL, &
    comm, shared%baseptr, shared%win, ierr)

    target_rank=0

    if(myid.ne.0) then
    CALL MPI_Win_shared_query(shared%win, target_rank, shared%window_size, &
    disp_unit, shared%baseptr, ierr)
    endif

end subroutine

subroutine shared_fence(shared,parallel)
    use parallel_type, only: parallel_struct
    type(shared_struct), intent(inout) :: shared
    type(parallel_struct), intent(in) :: parallel
    integer :: ierr
    if(parallel%nproc_node.gt.1) CALL MPI_WIN_FENCE(0, shared%win, ierr)
end subroutine



end module