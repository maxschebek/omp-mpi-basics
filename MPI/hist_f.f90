! This program calculates a histogram with 16 bins (the least
! significant four bits) from the results of the standard
! rand_r() random number generator.
! MPI is used for parallelization.
program histo
    use mpi
    use, intrinsic :: iso_c_binding
    implicit none

    interface
        integer(c_int) function rand_r(seed) bind(c)
            import :: c_int
            integer(c_int), intent(in) :: seed
        end function
    end interface

    integer :: i, n
    integer(c_int) :: value, seed
    ! Total histogram
    integer, dimension(0:15) :: hist
    ! Local histogram  sent from the individual processes
    integer, dimension(0:15) :: hist_send
    ! MPI variables
    integer :: size, rank, ierror, status(MPI_STATUS_SIZE)
    ! Total number of random examples
    integer, parameter :: n_total = 2E9
    ! Number of random examples per rank
    integer :: n_local
    ! Timing variables
    real(8) :: t_start, t_end

    call MPI_INIT(ierror)
    ! Get total number of ranks and number of local rank
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

    seed = rank*6

    n_local = n_total/size

    if (rank < mod(n_total, size)) then
        n_local = n_local + 1
    endif

    t_start = MPI_Wtime()

    ! Compute individual histogram on each rank
    do i = 1, n_local
        value = rand_r(seed)
        hist(IAND(value, 15)) = hist(IAND(value, 15)) + 1
    enddo

    ! Collect and add all individual histograms on rank 0
    if (rank == 0) then
        do i = 1, size - 1
            call MPI_Recv(hist_send, 16, MPI_INT, i, 0, &
                          MPI_COMM_WORLD, status, ierror)
            hist = hist + hist_send
        end do
        do i = 0, 15
            print *, "hist (", i, ")", hist(i)
        enddo
        t_end = MPI_Wtime()

        write (*, '(a, f7.4)') 'Time: ', t_end - t_start

    else
        call MPI_Send(hist, 16, MPI_INT, 0, 0, MPI_COMM_WORLD, ierror)
    end if

    call MPI_FINALIZE(ierror)

end program histo
