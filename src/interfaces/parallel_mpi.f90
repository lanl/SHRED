module parallel_mod
    use types, only : op, dp
    use parallel_type, only : parallel_struct
    use mpi


    implicit none
    public :: initialize_parallel, finalize_parallel, parallel_task

    interface parallel_task
        module procedure parallel_int
        module procedure parallel_ints
        module procedure parallel_ints2
        module procedure parallel_ints3
        module procedure parallel_ints4
        module procedure parallel_double
        module procedure parallel_doubles
        module procedure parallel_doubles2
        module procedure parallel_doubles3
        module procedure parallel_doubles4
        module procedure parallel_complex
        module procedure parallel_complexs
        module procedure parallel_complexs2
        module procedure parallel_complexs3
        module procedure parallel_complexs4
        module procedure parallel_logical
        module procedure parallel_logicals
        module procedure parallel_character
    end interface

 

    contains    
    
    subroutine parallel_wait(parallel, comm)
        character(len=*), intent(in), optional :: comm
        type(parallel_struct), intent(in) :: parallel
        integer :: comm_send
        integer :: ierr

        if(present(comm)) then
            call comm_check(comm,comm_send, parallel, ierr)
            call mpi_barrier(comm_send, ierr)
        else
            call mpi_barrier(parallel%comm_all, ierr)
        endif

    end subroutine

    subroutine initialize_parallel(parallel, comm_world_in)
        type(parallel_struct), intent(inout) :: parallel
        integer :: ierr, color, u, comm_world
        integer, optional :: comm_world_in
        integer :: n_spin_groups, n_k_groups, n_band_groups, i, nproc_i, j, n_space, my_node_check,nproc_FFT_y
        logical :: by_node
        namelist/parallel_groups/ n_spin_groups, n_k_groups, n_band_groups, by_node, nproc_FFT_y

        if(present(comm_world_in)) then 
            comm_world = comm_world_in
        else
            !Initialize MPI
            call mpi_init(ierr)

            !Initialize threads
            parallel%comm_all  = MPI_COMM_WORLD
            call mpi_comm_rank(parallel%comm_all, parallel%myid, ierr)
            call mpi_comm_size(parallel%comm_all, parallel%nproc, ierr)    
        endif

        n_band_groups=1
        n_k_groups=1
        n_spin_groups=1
        nproc_FFT_y=1
        by_node=.false.

        if(parallel%myid.eq.0) then
            open(newunit=u, file="input", status="old")
            read(u,nml=parallel_groups)
            close(u)
            parallel%n_spin_groups=n_spin_groups
            parallel%n_k_groups=n_k_groups
            parallel%n_band_groups=n_band_groups
            parallel%by_node=by_node
            parallel%nproc_fft_xy=nproc_FFT_y

            if(parallel%n_spin_groups.gt.2) then
                print *, 'Error: There are only 2 spin components to parallelize over (alpha / beta): n_spin_groups:', &
                     parallel%n_spin_groups
                     print *, 'fixing, will continue if number of processors & groups still works out'
                     parallel%n_spin_groups=2
            endif
            if(parallel%n_k_groups.lt.1) then
                print *, 'Error: You cannot have negative or zero k-points or n_k_groups, need at least 1 (setting to 1)'
                parallel%n_k_groups=1
            endif
        endif

        call parallel_task('bcast',parallel%n_spin_groups, parallel, 'all', root=0)
        call parallel_task('bcast',parallel%n_k_groups, parallel, 'all', root=0)
        call parallel_task('bcast',parallel%n_band_groups, parallel, 'all', root=0)
        call parallel_task('bcast',parallel%by_node, parallel, 'all', root=0)
        call parallel_task('bcast',parallel%nproc_fft_xy, parallel, 'all', root=0)

        !Split  for same nodes (shared memory)
        CALL MPI_Comm_split_type(parallel%comm_all, MPI_COMM_TYPE_SHARED, 0, &
        MPI_INFO_NULL, parallel%comm_node, ierr)
        call mpi_comm_rank(parallel%comm_node, parallel%myid_node, ierr)
        call mpi_comm_size(parallel%comm_node, parallel%nproc_node, ierr)

        parallel%my_node=-1
        j=0
        do i=0, parallel%nproc-1
            if(parallel%myid_node.eq.0.and.parallel%myid.eq.i) then
                parallel%my_node=j
                j=j+1
            endif
            call parallel_task('bcast',parallel%my_node, parallel, 'node', root=0)
            call parallel_task('bcast',j, parallel, 'all', root=i)
        enddo

        !create comm for same node_id (not on same node)
        call mpi_comm_split(parallel%comm_all, parallel%myid_node, parallel%myid, parallel%comm_diff_node, ierr)
        call mpi_comm_rank(parallel%comm_diff_node, parallel%myid_diff_node, ierr)
        call mpi_comm_size(parallel%comm_diff_node, parallel%nproc_diff_node, ierr) 
        !only actually # of nodes if myid_node=0 or evenly distibuted processors across nodes
        parallel%nnodes=parallel%nproc_diff_node
        call parallel_task('bcast', parallel%nnodes, parallel, 'node', root=0)
        allocate(parallel%nproc_all_node(parallel%nnodes))
        do i=1, parallel%nnodes
            nproc_i=parallel%nproc_node
            if(parallel%myid_node.eq.0) call parallel_task('bcast', nproc_i, parallel, 'diff_node', root=(i-1))
            call parallel_task('bcast', nproc_i, parallel, 'node', root=0)
            parallel%nproc_all_node(i)= nproc_i
        enddo

        call parallel_wait(parallel)
        parallel%my_node_ordered_id=sum(parallel%nproc_all_node(:parallel%my_node)) + parallel%myid_node

        do i=0,parallel%nproc-1
            if(parallel%myid.eq.i) print *, 'myid:', i, ', mynode_ordered_id ', parallel%my_node_ordered_id, &
                                             'mynode: ', parallel%my_node+1, 'of', parallel%nnodes, 'has ', &
                                             parallel%nproc_all_node(parallel%my_node+1), 'processors'
        enddo

        call parallel_wait(parallel)

        if(parallel%nproc.lt.(parallel%n_band_groups*parallel%n_k_groups*parallel%n_spin_groups)) then
            print *, 'Error: communicator is over-split, not enough threads available: N_proc: ', parallel%nproc, &
                    ' n_band_groups*n_k_groups*n_spin_groups: ', &
                    (parallel%n_band_groups*parallel%n_k_groups*parallel%n_spin_groups)
            stop
        endif

        if(mod(parallel%nproc,(parallel%n_band_groups*parallel%n_k_groups*parallel%n_spin_groups)).ne.0) then
            print *, 'number threads must be divisible by number of total splits: N_proc: ', parallel%nproc, &
                     ' n_band_groups*n_k_groups*n_spin_groups: ', &
                     (parallel%n_band_groups*parallel%n_k_groups*parallel%n_spin_groups)
            stop
        endif
        n_space=parallel%nproc/(parallel%n_band_groups*parallel%n_k_groups*parallel%n_spin_groups)
        if(parallel%myid.eq.0) print *, 'Nproc for spacial decompositions: ', n_space
        if(parallel%by_node.and.mod(parallel%nproc_node,n_space).ne.0) then
            print *, 'Error in MPI setup in by_node, spatial grid seperated across nodes &
                & you need n_band_groups*n_k_groups*n_spin_groups to be evenly divisible by the number of nodes you use'
        endif    
        !Split for spin
        color=mod(parallel%myid,parallel%n_band_groups*parallel%n_k_groups*parallel%n_spin_groups)
        if(parallel%n_spin_groups.gt.1) then
            if(parallel%by_node) then
                parallel%my_spin_group=parallel%my_node_ordered_id/(parallel%nproc/parallel%n_spin_groups)
                call mpi_comm_split(parallel%comm_all, parallel%my_spin_group, &
                                    parallel%my_node_ordered_id, parallel%comm_spin, ierr)
            else
                parallel%my_spin_group=mod(color,parallel%n_spin_groups)
                call mpi_comm_split(parallel%comm_all, parallel%my_spin_group, parallel%myid, &
                parallel%comm_spin, ierr)
            endif
            
            call mpi_comm_rank(parallel%comm_spin, parallel%myid_spin, ierr)
            call mpi_comm_size(parallel%comm_spin, parallel%nproc_spin, ierr)

        else
            parallel%comm_spin=parallel%comm_all
            parallel%my_spin_group=0
            if(parallel%by_node) then
                parallel%myid_spin=parallel%my_node_ordered_id
            else
                parallel%myid_spin=parallel%myid
            endif
            parallel%nproc_spin=parallel%nproc
        endif

        parallel%my_diff_s_group=parallel%myid_spin
        parallel%n_diff_s_groups=parallel%nproc_spin

        call mpi_comm_split(parallel%comm_all, parallel%my_diff_s_group, parallel%my_spin_group, &
            parallel%comm_diff_s, ierr)
        call mpi_comm_rank(parallel%comm_diff_s, parallel%myid_diff_s, ierr)
        call mpi_comm_size(parallel%comm_diff_s, parallel%nproc_diff_s, ierr)


        !Split for K-points
            color=mod(parallel%myid_spin,parallel%n_band_groups*parallel%n_k_groups)
            if(parallel%n_k_groups.gt.1) then
                if(parallel%by_node) then
                    parallel%my_k_group=parallel%myid_spin/(parallel%nproc_spin/parallel%n_k_groups)
                else
                    parallel%my_k_group=mod(color,parallel%n_k_groups)
                endif
                call mpi_comm_split(parallel%comm_spin, parallel%my_k_group, parallel%myid_spin, &
                    parallel%comm_k, ierr)
                call mpi_comm_rank(parallel%comm_k, parallel%myid_k, ierr)
                call mpi_comm_size(parallel%comm_k, parallel%nproc_k, ierr)

            else
                parallel%comm_k=parallel%comm_spin
                parallel%my_k_group=0
                parallel%myid_k=parallel%myid_spin
                parallel%nproc_k=parallel%nproc_spin
            endif

            parallel%my_diff_k_group=parallel%myid_k
            parallel%n_diff_k_groups=parallel%nproc_k

            call mpi_comm_split(parallel%comm_all, parallel%my_diff_k_group, parallel%my_k_group, &
                parallel%comm_diff_k, ierr)
            call mpi_comm_rank(parallel%comm_diff_k, parallel%myid_diff_k, ierr)
            call mpi_comm_size(parallel%comm_diff_k, parallel%nproc_diff_k, ierr)

        !Split for Bands
            color=mod(parallel%myid_k,parallel%n_band_groups)
            if(parallel%n_band_groups.gt.1) then
                if(parallel%by_node) then
                    parallel%my_band_group=parallel%myid_k/(parallel%nproc_k/parallel%n_band_groups)
                else
                    parallel%my_band_group=mod(color,parallel%n_band_groups)
                endif
                call mpi_comm_split(parallel%comm_k, parallel%my_band_group, parallel%myid_k, &
                    parallel%comm_bands, ierr)
                call mpi_comm_rank(parallel%comm_bands, parallel%myid_band, ierr)
                call mpi_comm_size(parallel%comm_bands, parallel%nproc_band, ierr)
                    
            else
                parallel%comm_bands=parallel%comm_k
                parallel%my_band_group=0
                parallel%myid_band=parallel%myid_k
                parallel%nproc_band=parallel%nproc_k
            endif

            parallel%my_diff_b_group=parallel%myid_band
            parallel%n_diff_b_groups=parallel%nproc_band

            call mpi_comm_split(parallel%comm_k, parallel%my_diff_b_group, parallel%my_band_group, &
                parallel%comm_diff_b, ierr)
            call mpi_comm_rank(parallel%comm_diff_b, parallel%myid_diff_b, ierr)
            call mpi_comm_size(parallel%comm_diff_b, parallel%nproc_diff_b, ierr)

        !Split  world for same space (different k's, bands, and spins)
            parallel%my_space_group=parallel%myid_band
            call mpi_comm_split(parallel%comm_all, parallel%my_space_group, parallel%myid, &
                parallel%comm_space, ierr)
            call mpi_comm_rank(parallel%comm_space, parallel%myid_space, ierr)
            call mpi_comm_size(parallel%comm_space, parallel%nproc_space, ierr)

        !Create subgroup (split of comm_k) for only some of the band groups
            if(parallel%n_band_groups.gt.4) then
                parallel%n_band_groups_sub=4
                if(parallel%my_band_group.lt.4) then
                    parallel%my_sub=0
                else
                    parallel%my_sub=1
                endif
                call mpi_comm_split(parallel%comm_k, parallel%my_sub, parallel%myid_k, &
                    parallel%comm_sub, ierr)
                call mpi_comm_rank(parallel%comm_sub, parallel%myid_sub, ierr)
                call mpi_comm_size(parallel%comm_sub, parallel%nproc_sub, ierr)
                call mpi_comm_split(parallel%comm_sub, parallel%my_diff_b_group, parallel%my_band_group, &
                        parallel%comm_diff_b_sub, ierr)
                call mpi_comm_rank(parallel%comm_diff_b_sub, parallel%myid_diff_b_sub, ierr)
                call mpi_comm_size(parallel%comm_diff_b_sub, parallel%nproc_diff_b_sub, ierr)
            else
                parallel%comm_sub=parallel%comm_k
                parallel%myid_sub=parallel%myid_k
                parallel%nproc_sub=parallel%nproc_k
                parallel%my_sub=0
                parallel%n_band_groups_sub=parallel%n_band_groups
                parallel%comm_diff_b_sub=parallel%comm_diff_b
                parallel%myid_diff_b_sub=parallel%myid_diff_b
                parallel%nproc_diff_b=parallel%nproc_diff_b
            endif

        !Create sub-band_groups for extra FFT parallelism over y (real) / x (Recip) coordinates
            color=mod(parallel%myid_band, parallel%nproc_fft_xy) 
            if(parallel%nproc_fft_xy.gt.1) then
                !if multiple proc used for y/x coordinates, nees to create smaller groups for z/y parallel

                !i.e. these have the same x_local, y_local_xy
                parallel%my_fft_yz=mod(color,parallel%nproc_fft_xy)
                call mpi_comm_split(parallel%comm_bands, parallel%my_fft_yz, parallel%myid_band, &
                    parallel%comm_fft_yz, ierr)
                call mpi_comm_rank(parallel%comm_fft_yz, parallel%myid_fft_yz, ierr)
                call mpi_comm_size(parallel%comm_fft_yz, parallel%nproc_fft_yz, ierr)
            else
                parallel%my_fft_yz=0
                parallel%comm_fft_yz=parallel%comm_bands
                parallel%myid_fft_yz=parallel%myid_band
                parallel%nproc_fft_yz=parallel%nproc_band
            endif

            !i.e. these have the same z_local, y_local_yz
            parallel%my_fft_xy=parallel%myid_fft_yz
            call mpi_comm_split(parallel%comm_bands, parallel%my_fft_xy, parallel%myid_band, &
                parallel%comm_fft_xy, ierr)
            call mpi_comm_rank(parallel%comm_fft_xy, parallel%myid_fft_xy, ierr)
            call mpi_comm_size(parallel%comm_fft_xy, parallel%nproc_fft_xy, ierr)

    if(parallel%by_node) then
        my_node_check=parallel%my_node
        call parallel_task('bcast', my_node_check, parallel, 'band', root=0)
        if(parallel%my_node.ne.my_node_check) then
            print *, 'somehow the by_node parallelization did not force all &
                        & comm_band members to be on the same node:', &
            parallel%myid, my_node_check, parallel%my_node
            stop
        endif
    endif

    j=0
    parallel%my_band_node=-1
    if(parallel%by_node) then
        do i=0, parallel%nproc_node-1
            if(parallel%myid_node.eq.i) then
                if(parallel%myid_band.eq.0) then
                    parallel%my_band_node=j
                    j=j+1
                endif
            endif
            call parallel_task('bcast', parallel%my_band_node, parallel, 'band', root=0)
            call parallel_task('bcast', j, parallel, 'node', root=i)
        enddo
    endif
        
    do i=0, parallel%nproc
            if(parallel%my_node_ordered_id.eq.i.and.parallel%by_node) then
                print *, parallel%my_node_ordered_id, parallel%myid_node, parallel%myid_space, parallel%myid_band, &
                         parallel%myid_k, parallel%myid_spin, parallel%my_band_node
            else if(parallel%myid.eq.i) then
                 print *, parallel%myid, parallel%myid_node, parallel%myid_space, parallel%myid_band, &
                          parallel%myid_k, parallel%myid_spin
            endif
            flush(6)
            call parallel_wait(parallel)
    enddo

    do i=0, parallel%nproc
        if(parallel%my_node_ordered_id.eq.i.and.parallel%by_node) then
            print *, parallel%my_node_ordered_id, parallel%myid_node, parallel%myid_space, parallel%myid_band, &
                     parallel%myid_k, parallel%myid_spin
        else if(parallel%myid.eq.i) then
             print *, parallel%myid, parallel%myid_node, parallel%myid_space, parallel%myid_band, &
                      parallel%myid_k, parallel%myid_spin
        endif
        flush(6)
        call parallel_wait(parallel)
enddo
    end subroutine

    subroutine finalize_parallel(parallel)
        type(parallel_struct), intent(in) :: parallel
        integer :: ierr=0
        flush(6)
        call parallel_wait(parallel)
        call mpi_finalize(ierr)
    end subroutine

    subroutine comm_check(comm,comm_send,parallel,ierr)
        type(parallel_struct), intent(in) :: parallel
        character(len=*), intent(in) :: comm
        integer,intent(out) :: comm_send
        integer, intent(out) :: ierr

        ierr=0
        if(comm.eq.'all') then
            comm_send=parallel%comm_all
            if(parallel%nproc.eq.1) ierr=1
        elseif(comm.eq.'band') then
            comm_send=parallel%comm_bands
            if(parallel%nproc_band.eq.1) ierr=1
        elseif(comm.eq.'k') then
            comm_send=parallel%comm_k
            if(parallel%nproc_k.eq.1) ierr=1
        elseif(comm.eq.'spin') then
            comm_send=parallel%comm_spin
            if(parallel%nproc_spin.eq.1) ierr=1
        elseif(comm.eq.'diff_b') then
            comm_send=parallel%comm_diff_b
            if(parallel%nproc_diff_b.eq.1) ierr=1
            elseif(comm.eq.'diff_k') then
            comm_send=parallel%comm_diff_k
            if(parallel%nproc_diff_k.eq.1) ierr=1
        elseif(comm.eq.'diff_s') then
            comm_send=parallel%comm_diff_s
            if(parallel%nproc_diff_s.eq.1) ierr=1
        elseif(comm.eq.'space') then
                comm_send=parallel%comm_space
                if(parallel%nproc_space.eq.1) ierr=1
        else if(comm.eq.'node') then
                comm_send=parallel%comm_node
                if(parallel%nproc_node.eq.1) ierr=1
            else if(comm.eq.'diff_node') then
                comm_send=parallel%comm_diff_node
                if(parallel%nproc_diff_node.eq.1) ierr=1
        else
            if(parallel%myid.eq.0) print *, 'Unrecognized comm: ', comm
        endif
    end subroutine

    subroutine operation_check(operation,operation_do)

        character(len=*), intent(in) :: operation
        integer,intent(out) :: operation_do
        if(trim(operation).eq.'max') then
            operation_do=MPI_MAX
        elseif(trim(operation).eq.'min') then
            operation_do=MPI_MIN
        elseif(trim(operation).eq.'sum') then
            operation_do=MPI_SUM
        elseif(trim(operation).eq.'prd') then
            operation_do=MPI_PROD
        elseif(trim(operation).eq.'land') then
            operation_do=MPI_LAND
        elseif(trim(operation).eq.'lor') then
            operation_do=MPI_LOR
        elseif(trim(operation).eq.'band') then
            operation_do=MPI_BAND
        elseif(trim(operation).eq.'bor') then
            operation_do=MPI_BOR
        elseif(trim(operation).eq.'maxloc') then
            operation_do=MPI_MAXLOC
        elseif(trim(operation).eq.'minloc') then
            operation_do=MPI_MINLOC
        else
            print *, 'Error: unrecognized parallel operation: ', trim(operation)
            stop
        endif
    end subroutine
        
    subroutine parallel_int(operation, send, parallel, comm, recv, root, to, tag, send_count, recv_count)
        type(parallel_struct), intent(in) :: parallel
        integer, intent(inout) :: send
        integer, optional, intent(inout) :: recv
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, red_operation, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        integer :: send_add(1)
        integer :: recv_add(1)
        logical  :: present_recv_add


        datatype=MPI_INTEGER
        send_add(1)=send
        size_send=1
        present_recv_add=present(recv)
        if(present_recv_add) then
            recv_add(1)=recv
            size_recv=1
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                endif
            else if(trim(operation).eq.'alltoall') then
                    print *, 'Error: alltoall an array to send'
            else !reduce
                call operation_check(operation,red_operation)
                if(present_recv_add) then
                    call mpi_allreduce(send_add, recv_add, size_send, datatype, red_operation, comm_send, ierr)
                else
                    call mpi_allreduce(MPI_IN_PLACE, send_add, size_send, datatype, red_operation, comm_send, ierr)
                endif
            endif
        endif

        if(present_recv_add) recv=recv_add(1)
        send=send_add(1)
    end subroutine

    subroutine parallel_ints(operation, send_add, parallel, comm, recv_add, root, to, tag, send_count, recv_count, &
                             send_countv, recv_countv, sdispls, rdispls)
        type(parallel_struct), intent(in) :: parallel
        integer, intent(inout) :: send_add(:)
        integer, optional, intent(inout) :: recv_add(:)
        integer, optional, intent(in) :: send_countv(:), recv_countv(:), sdispls(:), rdispls(:)
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, red_operation, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        logical  :: present_recv_add


        datatype=MPI_INTEGER
        size_send=size(send_add)
        present_recv_add=present(recv_add)
        if(present_recv_add) then
            size_recv=size(recv_add)
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype, root, comm_send, ierr)
                    return
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                    return
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                    return
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                    return
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count,  &
                    datatype, comm_send, ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_count).and.present(rdispls)) then
                    call mpi_allgatherv(send_add, send_count, datatype, recv_add, recv_countv, &
                    rdispls, datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                    print *, '      allgatherv requires recv, send_count, recv_countv, rdispls'
                endif
            else if(trim(operation).eq.'alltoall') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_alltoall(send_add, send_count, datatype, recv_add, recv_count, datatype,comm_send,ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_countv)) then

                     if(present(sdispls).and.present(rdispls)) then
                        call mpi_alltoallv(send_add, send_countv, sdispls, datatype, recv_add, recv_countv,rdispls, &
                         datatype,comm_send,ierr)
                         return
                     else
                        print *, 'Error: alltoallv requires recv, send_count, recv_count'
                     endif
                else
                    print *, 'Error: alltoallv requires sdispls, rdispls'
                endif
            else !reduce
                call operation_check(operation,red_operation)
                if(present_recv_add) then
                    call mpi_allreduce(send_add, recv_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                else
                    call mpi_allreduce(MPI_IN_PLACE, send_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                endif
            endif



        endif

    end subroutine

    subroutine parallel_ints2(operation, send_add, parallel, comm, recv_add, root, to, tag, send_count, recv_count, &
        send_countv, recv_countv, sdispls, rdispls)
        type(parallel_struct), intent(in) :: parallel
        integer, intent(inout) :: send_add(:,:)
        integer, optional, intent(inout) :: recv_add(:,:)
        integer, optional, intent(in) :: send_countv(:), recv_countv(:), sdispls(:), rdispls(:)
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, red_operation, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        logical  :: present_recv_add


        datatype=MPI_INTEGER
        size_send=size(send_add)
        present_recv_add=present(recv_add)
        if(present_recv_add) then
            size_recv=size(recv_add)
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                    return
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                    return
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                    return
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                    return
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_count).and.present(rdispls)) then
                    call mpi_allgatherv(send_add, send_count, datatype, recv_add, recv_countv, &
                    rdispls, datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                    print *, '      allgatherv requires recv, send_count, recv_countv, rdispls'
                endif
            else if(trim(operation).eq.'alltoall') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_alltoall(send_add, send_count, datatype, recv_add, recv_count, datatype,comm_send,ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_countv)) then
                    if(present(sdispls).and.present(rdispls)) then
                       call mpi_alltoallv(send_add, send_countv, sdispls, datatype, recv_add, recv_countv,rdispls, &
                        datatype,comm_send,ierr)
                        return
                    else
                       print *, 'Error: alltoall requires recv, send_count, recv_count'
                    endif
                else
                    print *, 'Error: alltoall requires recv, send_count, recv_count'
                endif
            else !reduce
                call operation_check(operation,red_operation)
                if(present_recv_add) then
                    call mpi_allreduce(send_add, recv_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                else
                    call mpi_allreduce(MPI_IN_PLACE, send_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                endif
            endif

        endif
    end subroutine

    subroutine parallel_ints3(operation, send_add, parallel, comm, recv_add, root, to, tag, send_count, recv_count, &
        send_countv, recv_countv, sdispls, rdispls)
        

        type(parallel_struct), intent(in) :: parallel
        integer, intent(inout) :: send_add(:,:,:)
        integer, optional, intent(inout) :: recv_add(:,:,:)
        integer, optional, intent(in) :: send_countv(:), recv_countv(:), sdispls(:), rdispls(:)
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, red_operation, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        logical  :: present_recv_add


        datatype=MPI_INTEGER
        size_send=size(send_add)
        present_recv_add=present(recv_add)
        if(present_recv_add) then
            size_recv=size(recv_add)
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                    return
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                    return
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                    return
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                    return
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_count).and.present(rdispls)) then
                    call mpi_allgatherv(send_add, send_count, datatype, recv_add, recv_countv, &
                    rdispls, datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                    print *, '      allgatherv requires recv, send_count, recv_countv, rdispls'

                endif
            else if(trim(operation).eq.'alltoall') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_alltoall(send_add, send_count, datatype, recv_add, recv_count, datatype,comm_send,ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_countv)) then
                    if(present(sdispls).and.present(rdispls)) then
                       call mpi_alltoallv(send_add, send_countv, sdispls, datatype, recv_add, recv_countv,rdispls, &
                        datatype,comm_send,ierr)
                        return
                    else
                       print *, 'Error: alltoallv requires recv, send_count, recv_count'
                    endif
                else
                    print *, 'Error: alltoall requires recv, send_count, recv_count'
                endif
            else !reduce
                call operation_check(operation,red_operation)
                if(present_recv_add) then
                    call mpi_allreduce(send_add, recv_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                else
                    call mpi_allreduce(MPI_IN_PLACE, send_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                endif
            endif

        endif
    end subroutine

    subroutine parallel_ints4(operation, send_add, parallel, comm, recv_add, root, to, tag, send_count, recv_count, &
        send_countv, recv_countv, sdispls, rdispls)
        

        type(parallel_struct), intent(in) :: parallel
        integer, intent(inout) :: send_add(:,:,:,:)
        integer, optional, intent(inout) :: recv_add(:,:,:,:)
        integer, optional, intent(in) :: send_countv(:), recv_countv(:), sdispls(:), rdispls(:)
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, red_operation, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        logical  :: present_recv_add


        datatype=MPI_INTEGER
        size_send=size(send_add)
        present_recv_add=present(recv_add)
        if(present_recv_add) then
            size_recv=size(recv_add)
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                    return
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                    return
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                    return
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                    return
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_count).and.present(rdispls)) then
                    call mpi_allgatherv(send_add, send_count, datatype, recv_add, recv_countv, &
                    rdispls, datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                    print *, '      allgatherv requires recv, send_count, recv_countv, rdispls'
                endif
            else if(trim(operation).eq.'alltoall') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_alltoall(send_add, send_count, datatype, recv_add, recv_count, datatype,comm_send,ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_countv)) then
                    if(present(sdispls).and.present(rdispls)) then
                       call mpi_alltoallv(send_add, send_countv, sdispls, datatype, recv_add, recv_countv,rdispls, &
                        datatype,comm_send,ierr)
                        return
                    else
                       print *, 'Error: alltoallv requires recv, send_count, recv_count'
                    endif
                else
                    print *, 'Error: alltoall requires recv, send_count, recv_count'
                endif
            else !reduce
                call operation_check(operation,red_operation)
                if(present_recv_add) then
                    call mpi_allreduce(send_add, recv_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                else
                    call mpi_allreduce(MPI_IN_PLACE, send_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                endif
            endif

        endif
    end subroutine

    subroutine parallel_double(operation, send, parallel, comm, recv, root, to, tag, send_count, recv_count)
        

        type(parallel_struct), intent(in) :: parallel
        real(dp), intent(inout) :: send
        real(dp), optional, intent(inout) :: recv
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, red_operation, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        real(dp) :: send_add(1)
        real(dp) :: recv_add(1)
        logical  :: present_recv_add


        datatype=MPI_DOUBLE_PRECISION
        send_add(1)=send
        size_send=1
        present_recv_add=present(recv)
        if(present_recv_add) then
            recv_add(1)=recv
            size_recv=1
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                endif
            else if(trim(operation).eq.'alltoall') then
                print *, 'Error: alltoall needs an array to send'
            else !reduce
                call operation_check(operation,red_operation)
                if(present_recv_add) then
                    call mpi_allreduce(send_add, recv_add, size_send, datatype, red_operation, comm_send, ierr)
                else
                    call mpi_allreduce(MPI_IN_PLACE, send_add, size_send, datatype, red_operation, comm_send, ierr)
                endif
            endif
        endif

        if(present_recv_add) recv=recv_add(1)
        send=send_add(1)
    end subroutine

    subroutine parallel_doubles(operation, send_add, parallel, comm, recv_add, root, to, tag, send_count, recv_count, &
        send_countv, recv_countv, sdispls, rdispls)
        

        type(parallel_struct), intent(in) :: parallel
        real(dp), intent(inout) :: send_add(:)
        real(dp), optional, intent(inout) :: recv_add(:)
        integer, optional, intent(in) :: send_countv(:), recv_countv(:), sdispls(:), rdispls(:)
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, red_operation, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        logical  :: present_recv_add


        datatype=MPI_DOUBLE_PRECISION
        size_send=size(send_add)
        present_recv_add=present(recv_add)
        if(present_recv_add) then
            size_recv=size(recv_add)
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                    return
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                    return
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                    return
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                    return
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_count).and.present(rdispls)) then
                    call mpi_allgatherv(send_add, send_count, datatype, recv_add, recv_countv, &
                    rdispls, datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                    print *, '      allgatherv requires recv, send_count, recv_countv, rdispls'
                endif
            else if(trim(operation).eq.'alltoall') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_alltoall(send_add, send_count, datatype, recv_add, recv_count, datatype,comm_send,ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_countv)) then
                    if(present(sdispls).and.present(rdispls)) then
                       call mpi_alltoallv(send_add, send_countv, sdispls, datatype, recv_add, recv_countv,rdispls, &
                        datatype,comm_send,ierr)
                        return
                    else
                       print *, 'Error: alltoallv requires recv, send_count, recv_count'
                    endif
                else
                    print *, 'Error: alltoall requires recv, send_count, recv_count'
                endif
            else !reduce
                call operation_check(operation,red_operation)
                if(present_recv_add) then
                    call mpi_allreduce(send_add, recv_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                else
                    call mpi_allreduce(MPI_IN_PLACE, send_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                endif
            endif


        endif

    end subroutine

    subroutine parallel_doubles2(operation, send_add, parallel, comm, recv_add, root, to, tag, send_count, recv_count, &
        send_countv, recv_countv, sdispls, rdispls)
        

        type(parallel_struct), intent(in) :: parallel
        real(dp), intent(inout) :: send_add(:,:)
        real(dp), optional, intent(inout) :: recv_add(:,:)
        integer, optional, intent(in) :: send_countv(:), recv_countv(:), sdispls(:), rdispls(:)
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, red_operation, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        logical  :: present_recv_add


        datatype=MPI_DOUBLE_PRECISION
        size_send=size(send_add)
        present_recv_add=present(recv_add)
        if(present_recv_add) then
            size_recv=size(recv_add)
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                    return
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                    return
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                    return
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                    return
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_count).and.present(rdispls)) then
                    call mpi_allgatherv(send_add, send_count, datatype, recv_add, recv_countv, &
                    rdispls, datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                    print *, '      allgatherv requires recv, send_count, recv_countv, rdispls'
                endif
            else if(trim(operation).eq.'alltoall') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_alltoall(send_add, send_count, datatype, recv_add, recv_count, datatype,comm_send,ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_countv)) then
                    if(present(sdispls).and.present(rdispls)) then
                       call mpi_alltoallv(send_add, send_countv, sdispls, datatype, recv_add, recv_countv,rdispls, &
                        datatype,comm_send,ierr)
                        return
                    else
                       print *, 'Error: allgather requires recv, send_count, recv_count'
                    endif
                else
                    print *, 'Error: alltoall requires recv, send_count, recv_count'
                endif
            else !reduce
                call operation_check(operation,red_operation)
                if(present_recv_add) then
                    call mpi_allreduce(send_add, recv_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                else
                    call mpi_allreduce(MPI_IN_PLACE, send_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                endif
            endif

        endif
    end subroutine

    subroutine parallel_doubles3(operation, send_add, parallel, comm, recv_add, root, to, tag, send_count, recv_count, &
        send_countv, recv_countv, sdispls, rdispls)
        

        type(parallel_struct), intent(in) :: parallel
        real(dp), intent(inout) :: send_add(:,:,:)
        real(dp), optional, intent(inout) :: recv_add(:,:,:)
        integer, optional, intent(in) :: send_countv(:), recv_countv(:), sdispls(:), rdispls(:)
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, red_operation, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        logical  :: present_recv_add


        datatype=MPI_DOUBLE_PRECISION
        size_send=size(send_add)
        present_recv_add=present(recv_add)
        if(present_recv_add) then
            size_recv=size(recv_add)
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                    return
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                    return
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                    return
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                    return
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_count).and.present(rdispls)) then
                    call mpi_allgatherv(send_add, send_count, datatype, recv_add, recv_countv, &
                    rdispls, datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                    print *, '      allgatherv requires recv, send_count, recv_countv, rdispls'
                endif
            else if(trim(operation).eq.'alltoall') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_alltoall(send_add, send_count, datatype, recv_add, recv_count, datatype,comm_send,ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_countv)) then
                    if(present(sdispls).and.present(rdispls)) then
                       call mpi_alltoallv(send_add, send_countv, sdispls, datatype, recv_add, recv_countv,rdispls, &
                        datatype,comm_send,ierr)
                        return
                    else
                       print *, 'Error: alltoall requires recv, send_count, recv_count'
                    endif
                else
                    print *, 'Error: alltoall requires recv, send_count, recv_count'
                endif
            else !reduce
                call operation_check(operation,red_operation)
                if(present_recv_add) then
                    call mpi_allreduce(send_add, recv_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                else
                    call mpi_allreduce(MPI_IN_PLACE, send_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                endif
            endif

        endif
    end subroutine

    subroutine parallel_doubles4(operation, send_add, parallel, comm, recv_add, root, to, tag, send_count, recv_count, &
        send_countv, recv_countv, sdispls, rdispls)
        

        type(parallel_struct), intent(in) :: parallel
        real(dp), intent(inout) :: send_add(:,:,:,:)
        real(dp), optional, intent(inout) :: recv_add(:,:,:,:)
        integer, optional, intent(in) :: send_countv(:), recv_countv(:), sdispls(:), rdispls(:)
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, red_operation, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        logical  :: present_recv_add


        datatype=MPI_DOUBLE_PRECISION
        size_send=size(send_add)
        present_recv_add=present(recv_add)
        if(present_recv_add) then
            size_recv=size(recv_add)
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                    return
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                    return
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                    return
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_count).and.present(rdispls)) then
                    call mpi_allgatherv(send_add, send_count, datatype, recv_add, recv_countv, &
                    rdispls, datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                    print *, '      allgatherv requires recv, send_count, recv_countv, rdispls'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                endif
            else if(trim(operation).eq.'alltoall') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_alltoall(send_add, send_count, datatype, recv_add, recv_count, datatype,comm_send,ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_countv)) then
                    if(present(sdispls).and.present(rdispls)) then
                       call mpi_alltoallv(send_add, send_countv, sdispls, datatype, recv_add, recv_countv,rdispls, &
                        datatype,comm_send,ierr)
                        return
                    else
                       print *, 'Error: allgather requires recv, send_count, recv_count'
                    endif
                else
                    print *, 'Error: alltoall requires recv, send_count, recv_count'
                endif
            else !reduce
                call operation_check(operation,red_operation)
                if(present_recv_add) then
                    call mpi_allreduce(send_add, recv_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                else
                    call mpi_allreduce(MPI_IN_PLACE, send_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                endif
            endif

        endif
    end subroutine

    subroutine parallel_complex(operation, send, parallel, comm, recv, root, to, tag, send_count, recv_count)
        

        type(parallel_struct), intent(in) :: parallel
        complex(dp), intent(inout) :: send
        complex(dp), optional, intent(inout) :: recv
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, red_operation, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        complex(dp) :: send_add(1)
        complex(dp) :: recv_add(1)
        logical  :: present_recv_add


        datatype=MPI_DOUBLE_COMPLEX
        send_add(1)=send
        size_send=1
        present_recv_add=present(recv)
        if(present_recv_add) then
            recv_add(1)=recv
            size_recv=1
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                endif
            else if(trim(operation).eq.'alltoall') then
                print *, 'Error: alltoall an array to send'
            else !reduce
                call operation_check(operation,red_operation)
                if(present_recv_add) then
                    call mpi_allreduce(send_add, recv_add, size_send, datatype, red_operation, comm_send, ierr)
                else
                    call mpi_allreduce(MPI_IN_PLACE, send_add, size_send, datatype, red_operation, comm_send, ierr)
                endif
            endif
        endif

        if(present_recv_add) recv=recv_add(1)
        send=send_add(1)
    end subroutine

    subroutine parallel_complexs(operation, send_add, parallel, comm, recv_add, root, to, tag, send_count, recv_count, &
        send_countv, recv_countv, sdispls, rdispls)
        

        type(parallel_struct), intent(in) :: parallel
        complex(dp), intent(inout) :: send_add(:)
        complex(dp), optional, intent(inout) :: recv_add(:)
        integer, optional, intent(in) :: send_countv(:), recv_countv(:), sdispls(:), rdispls(:)
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, red_operation, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        logical  :: present_recv_add

        datatype=MPI_DOUBLE_COMPLEX
        size_send=size(send_add)
        present_recv_add=present(recv_add)
        if(present(recv_add)) then
            size_recv=size(recv_add)
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                    return
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                    return
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                    return
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                    return
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_count).and.present(rdispls)) then
                    call mpi_allgatherv(send_add, send_count, datatype, recv_add, recv_countv, &
                    rdispls, datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                    print *, '      allgatherv requires recv, send_count, recv_countv, rdispls'
                endif
            else if(trim(operation).eq.'alltoall') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_alltoall(send_add, send_count, datatype, recv_add, recv_count, datatype,comm_send,ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_countv)) then
                    if(present(sdispls).and.present(rdispls)) then
                       call mpi_alltoallv(send_add, send_countv, sdispls, datatype, recv_add, recv_countv,rdispls, &
                        datatype,comm_send,ierr)
                        return
                    else
                       print *, 'Error: alltoall requires recv, send_count, recv_count'
                    endif
                else
                    print *, 'Error: alltoall requires recv, send_count, recv_count'
                endif
            else !reduce
                call operation_check(operation,red_operation)
                if(present_recv_add) then
                    call mpi_allreduce(send_add, recv_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                else
                    call mpi_allreduce(MPI_IN_PLACE, send_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                endif
            endif

        endif
    end subroutine

    subroutine parallel_complexs2(operation, send_add, parallel, comm, recv_add, root, to, tag, send_count, recv_count, &
        send_countv, recv_countv, sdispls, rdispls)
        

        type(parallel_struct), intent(in) :: parallel
        complex(dp), intent(inout) :: send_add(:,:)
        complex(dp), optional, intent(inout) :: recv_add(:,:)
        integer, optional, intent(in) :: send_countv(:), recv_countv(:), sdispls(:), rdispls(:)
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, red_operation, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        logical  :: present_recv_add

        datatype=MPI_DOUBLE_COMPLEX
        size_send=size(send_add)
        present_recv_add=present(recv_add)
        if(present(recv_add)) then
            size_recv=size(recv_add)
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                    return
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                    return
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                    return
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                    return
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_count).and.present(rdispls)) then
                    call mpi_allgatherv(send_add, send_count, datatype, recv_add, recv_countv, &
                    rdispls, datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                    print *, '      allgatherv requires recv, send_count, recv_countv, rdispls'
                endif
            else if(trim(operation).eq.'alltoall') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_alltoall(send_add, send_count, datatype, recv_add, recv_count, datatype,comm_send,ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_countv)) then
                    if(present(sdispls).and.present(rdispls)) then
                       call mpi_alltoallv(send_add, send_countv, sdispls, datatype, recv_add, recv_countv,rdispls, &
                        datatype,comm_send,ierr)
                        return
                    else
                       print *, 'Error: alltoall requires recv, send_count, recv_count'
                    endif
                else
                    print *, 'Error: alltoall requires recv, send_count, recv_count'
                endif
            else !reduce
                call operation_check(operation,red_operation)
                if(present_recv_add) then
                    call mpi_allreduce(send_add, recv_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                else
                    call mpi_allreduce(MPI_IN_PLACE, send_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                endif
            endif

        endif
    end subroutine

    subroutine parallel_complexs3(operation, send_add, parallel, comm, recv_add, root, to, tag, send_count, recv_count, &
        send_countv, recv_countv, sdispls, rdispls)
        

        type(parallel_struct), intent(in) :: parallel
        complex(dp), intent(inout) :: send_add(:,:,:)
        complex(dp), optional, intent(inout) :: recv_add(:,:,:)
        integer, optional, intent(in) :: send_countv(:), recv_countv(:), sdispls(:), rdispls(:)
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, red_operation, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        logical  :: present_recv_add

        datatype=MPI_DOUBLE_COMPLEX
        size_send=size(send_add)
        present_recv_add=present(recv_add)
        if(present(recv_add)) then
            size_recv=size(recv_add)
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                    return
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                    return
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                    return
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                    return
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_count).and.present(rdispls)) then
                    call mpi_allgatherv(send_add, send_count, datatype, recv_add, recv_countv, &
                    rdispls, datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                    print *, '      allgatherv requires recv, send_count, recv_countv, rdispls'
                endif
            else if(trim(operation).eq.'alltoall') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_alltoall(send_add, send_count, datatype, recv_add, recv_count, datatype,comm_send,ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_countv)) then
                    if(present(sdispls).and.present(rdispls)) then
                       call mpi_alltoallv(send_add, send_countv, sdispls, datatype, recv_add, recv_countv,rdispls, &
                        datatype,comm_send,ierr)
                        return
                    else
                       print *, 'Error: alltoallv requires recv, send_count, recv_count'
                    endif
                else
                    print *, 'Error: alltoall requires recv, send_count, recv_count'
                endif
            else !reduce
                call operation_check(operation,red_operation)
                if(present_recv_add) then
                    call mpi_allreduce(send_add, recv_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                else
                    call mpi_allreduce(MPI_IN_PLACE, send_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                endif
            endif

        endif
    end subroutine

    subroutine parallel_complexs4(operation, send_add, parallel, comm, recv_add, root, to, tag, send_count, recv_count, &
        send_countv, recv_countv, sdispls, rdispls)
        

        type(parallel_struct), intent(in) :: parallel
        complex(dp), intent(inout) :: send_add(:,:,:,:)
        complex(dp), optional, intent(inout) :: recv_add(:,:,:,:)
        integer, optional, intent(in) :: send_countv(:), recv_countv(:), sdispls(:), rdispls(:)
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, red_operation, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        logical  :: present_recv_add

        datatype=MPI_DOUBLE_COMPLEX
        size_send=size(send_add)
        present_recv_add=present(recv_add)
        if(present(recv_add)) then
            size_recv=size(recv_add)
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                    return
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                    return
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                    return
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                    return
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_count).and.present(rdispls)) then
                    call mpi_allgatherv(send_add, send_count, datatype, recv_add, recv_countv, &
                    rdispls, datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                    print *, '      allgatherv requires recv, send_count, recv_countv, rdispls'
                endif
            else if(trim(operation).eq.'alltoall') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_alltoall(send_add, send_count, datatype, recv_add, recv_count, datatype,comm_send,ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_countv)) then
                    if(present(sdispls).and.present(rdispls)) then
                       call mpi_alltoallv(send_add, send_countv, sdispls, datatype, recv_add, recv_countv,rdispls, &
                        datatype,comm_send,ierr)
                        return
                    else
                       print *, 'Error: alltoallv requires recv, send_count, recv_count'
                    endif
                else
                    print *, 'Error: alltoall requires recv, send_count, recv_count'
                endif
            else !reduce
                call operation_check(operation,red_operation)
                if(present_recv_add) then
                    call mpi_allreduce(send_add, recv_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                else
                    call mpi_allreduce(MPI_IN_PLACE, send_add, size_send, datatype, red_operation, comm_send, ierr)
                    return
                endif
            endif

        endif
    end subroutine

    subroutine parallel_character(operation, send_add, parallel, comm, recv_add, root, to, tag, send_count, recv_count)
        

        type(parallel_struct), intent(in) :: parallel
        character(len=*), intent(inout) :: send_add
        character(len=*), optional, intent(inout) :: recv_add
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        logical  :: present_recv_add

        datatype=MPI_CHARACTER
        size_send=len(send_add)
        present_recv_add=present(recv_add)
        if(present(recv_add)) then
            size_recv=len(recv_add)
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                    return
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                    return
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                    return
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                    return
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                endif
            else if(trim(operation).eq.'alltoall') then
                print *, 'Error: alltoall requires an array to send'
            else ! no reduce for this type
            endif

        endif
    end subroutine

    subroutine parallel_logical(operation, send_add, parallel, comm, recv_add, root, to, tag, send_count, recv_count)
        
        type(parallel_struct), intent(in) :: parallel
        logical, intent(inout) :: send_add
        logical, optional, intent(inout) :: recv_add
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        logical  :: present_recv_add

        datatype=MPI_LOGICAL
        size_send=1
        present_recv_add=present(recv_add)
        if(present(recv_add)) then
            size_recv=1
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                    return
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                    return
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                    return
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                    return
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                endif
            else if(trim(operation).eq.'alltoall') then
                       print *, 'Error: alltoall requires an array'
            else ! no reduce for this type
            endif

        endif
    end subroutine

    subroutine parallel_logicals(operation, send_add, parallel, comm, recv_add, root, to, tag, send_count, recv_count, &
        send_countv, recv_countv, sdispls, rdispls)
        

        type(parallel_struct), intent(in) :: parallel
        logical, intent(inout) :: send_add(:)
        logical, optional, intent(inout) :: recv_add(:)
        integer, optional, intent(in) :: send_countv(:), recv_countv(:), sdispls(:), rdispls(:)
        character(len=*), intent(in) :: comm
        character(len=*), intent(in) :: operation
        integer, optional, intent(in) :: root, to, tag, send_count, recv_count
        integer :: comm_send, ierr, stat(MPI_STATUS_SIZE)
        integer :: datatype, size_send,size_recv
        logical  :: present_recv_add

        datatype=MPI_LOGICAL
        size_send=size(send_add)
        present_recv_add=present(recv_add)
        if(present(recv_add)) then
            size_recv=size(recv_add)
        endif

        if(.true.) then !just setting a wrap for visual studio code
            call comm_check(comm,comm_send,parallel, ierr)
            if((ierr.eq.1).and.(.not.present_recv_add)) return
            if(trim(operation).eq.'bcast') then
                if(present(root)) then
                    call mpi_bcast(send_add, size_send, datatype,  root, comm_send, ierr)
                    return
                else 
                    print *, 'Error: bcast requires root'
                    stop
                endif
            else if(trim(operation).eq.'send') then
                if(present(to).and.present(tag)) then
                    call mpi_send(send_add, size_send, datatype, to, tag, comm_send, ierr)
                    return
                else 
                    print *, 'Error: send requires to, tag'
                    stop
                endif
            else if(trim(operation).eq.'recv') then
                if(present(root).and.present(tag)) then
                    call mpi_recv(send_add, size_send, datatype, root, tag, comm_send, stat, ierr)
                    return
                else 
                    print *, 'Error: recv requires root, tag'
                    stop
                endif
            else if(trim(operation).eq.'gather') then
                if(present(root).and.present_recv_add) then
                    call mpi_gather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, root, comm_send, ierr)
                    return
                else
                    print *, 'Error: gather requires recv and root'
                endif
            else if(trim(operation).eq.'allgather') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_allgather(send_add, send_count, datatype, recv_add, recv_count, &
                    datatype, comm_send, ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_count).and.present(rdispls)) then
                    call mpi_allgatherv(send_add, send_count, datatype, recv_add, recv_countv, &
                    rdispls, datatype, comm_send, ierr)
                    return
                else
                    print *, 'Error: allgather requires recv, send_count, recv_count'
                    print *, '      allgatherv requires recv, send_count, recv_countv, rdispls'
                endif
            else if(trim(operation).eq.'alltoall') then
                if(present_recv_add.and.present(recv_count).and.present(send_count)) then
                    call mpi_alltoall(send_add, send_count, datatype, recv_add, recv_count, datatype,comm_send,ierr)
                    return
                else if(present_recv_add.and.present(recv_countv).and.present(send_countv)) then
                    if(present(sdispls).and.present(rdispls)) then
                       call mpi_alltoallv(send_add, send_countv, sdispls, datatype, recv_add, recv_countv,rdispls, &
                        datatype,comm_send,ierr)
                        return
                    else
                       print *, 'Error: alltoallv requires recv, send_count, recv_count'
                    endif
                else
                    print *, 'Error: alltoall requires recv, send_count, recv_count'
                endif
            else ! no reduce for this type
            endif

        endif
    end subroutine



end module
