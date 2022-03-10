! This program computes pi by the Monte Carlo method.
! OMP is used for parallelization.
program mcpi
!$  use omp_lib
    use, intrinsic :: iso_c_binding
    implicit none

    interface
        integer(c_int) function rand_r(seed) bind(c)
            import :: c_int
            integer(c_int), intent(in) :: seed
        end function rand_r
    end interface

    ! Random seed
    integer(c_int) :: seed
    ! Maximum random value
    integer(c_int), parameter :: RAND_MAX = 2147483647
    ! Cartesian coordinates
    real(8) :: x, y
    ! Pi as otbained by Monte Carlo
    real(8) :: pi
    ! Reference value for pi
    real(8), parameter :: M_PI = 3.14159265358979d0
    ! Counter for values in or outside the uni circle
    integer ::  count
    ! Running index and OMP thread information
    integer :: i, my_thread
    ! MPI error variable
    integer :: ierror
    ! Total number of random examples
    integer, parameter :: n_total = 342000000
    ! Timing variables
    real(8) :: t_start, t_end

    count = 0
    seed = 2
    call cpu_time(t_start)
!$  t_start = OMP_get_wtime()
    ! Each thread computes an individual number of random examples
    ! Addition is done automatically by using reduction
!$omp parallel private(seed, x, y, my_thread)
!$  my_thread = OMP_GET_THREAD_NUM()
!$  seed = my_thread + 1
!$omp do reduction(+:count)
    do i = 1, n_total
        x = rand_r(seed)/DBLE(RAND_MAX)
        y = rand_r(seed)/DBLE(RAND_MAX)
        if (sqrt(x**2 + y**2) < 1.d0) count = count + 1
    enddo
!$omp end do
!$omp end parallel
    pi = 4.d0*DBLE(count)/n_total
    call cpu_time(t_end)
!$  t_end = OMP_get_wtime()

    write (*, '(a20, f7.4, a20, E15.5)') 'Time: ', t_end - t_start, &
        'Accuracy: ', abs(M_PI - pi)/M_PI

end program mcpi
