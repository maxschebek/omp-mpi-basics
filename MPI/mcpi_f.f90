! This program computes pi by the Monte Carlo method.
! MPI is used for parallelization.
program mcpi
    use mpi
    use, intrinsic :: iso_c_binding
    implicit none

    interface
        integer(c_int) function rand_r(seed) bind(c)
            import :: c_int
            integer(c_int), intent(in) :: seed
        end function rand_r
    end interface

    integer(c_int) :: seed
    ! Maximum value produced by rand_r
    integer(c_int), parameter :: RAND_MAX = 2147483647
    real(8) :: x, y, pi
    ! Exact value for pi
    real(8), parameter :: M_PI = 3.14159265358979d0
    ! Run_totaling variable for
    integer :: i
    ! MPI variables
    integer :: ierror, size, rank, status(MPI_STATUS_SIZE)
    ! Local counter if point is in or outside the circle
    integer :: count
    ! Received value of counter from other ranks
    integer :: val
    ! Number of random (x,y) combinations
    integer, parameter :: n_total = 342000000
    ! Number of random (x,y) combinations for one processor
    integer :: n_local
    ! Timing variables
    real(8) :: t_start, t_end

    call MPI_INIT(ierror)
    ! Get total number of ranks and rank of calling processor
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

    count = 0
    seed = rank*6
    n_local = n_total/size

    if (rank < mod(n_total, size)) then
        n_local = n_local + 1
    endif
    t_start = MPI_Wtime()

    do i = 1, n_local
        x = rand_r(seed)/DBLE(RAND_MAX)
        y = rand_r(seed)/DBLE(RAND_MAX)
        if (sqrt(x**2 + y**2) < 1.d0) count = count + 1
    enddo

    ! Collect results from all participating processors on rank 0
    if (rank == 0) then
        do i = 1, size - 1
            call MPI_Recv(val, 1, MPI_INT, i, 0, &
                          MPI_COMM_WORLD, status, ierror)
            count = count + val
        end do
        pi = 4.d0*DBLE(count)/n_total
        t_end = MPI_Wtime()

        write (*, '(a20, f7.4, a20, E15.5)') 'Time: ', t_end - t_start, &
            'Accuracy: ', abs(M_PI - pi)/M_PI
    else
        call MPI_Send(count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, ierror)
    end if

    call MPI_FINALIZE(ierror)
end program mcpi
