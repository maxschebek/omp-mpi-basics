! This program calculates a histogram with 16 bins (the least
! significant four bits) from the results of the standard
! rand_r() random number generator.
! OMP is used for parallelization. 
program histo
!$  use omp_lib
    use, intrinsic :: iso_c_binding
    implicit none

    interface
        integer(c_int) function rand_r(seed) bind(c)
            import :: c_int
            integer(c_int), intent(in) :: seed
        end function
    end interface

    integer :: i, n
    ! Random value
    integer(c_int) :: value
    ! Random seed
    integer(c_int) ::  seed
    ! Histogram obtained by using OMP reduction
    integer(c_int), dimension(0:15) :: hist_red
    ! Histogram obtained by using OMP critical
    integer(c_int), dimension(0:15) ::  hist_loc, hist_crit
    ! Timing variables
    real(8) :: t_start, t_end

    seed = 123
    hist_red = 0
    hist_crit = 0
    call cpu_time(t_start)
!$  t_start = OMP_get_wtime()

! ---- PARALLEL LOOP WITH REDUCTION
!$omp parallel private(value, seed)
!$  seed = OMP_GET_THREAD_NUM()*6
!$OMP do reduction(+:hist_red)
    do i = 1, 2000000000
        value = rand_r(seed)
        hist_red(IAND(value, 15)) = hist_red(IAND(value, 15)) + 1
    enddo
!$OMP end do
!$OMP end parallel

! ---- PARALLEL LOOP WITH CRITICAL
!$omp parallel private(value, seed, hist_loc)
!$  seed = OMP_GET_THREAD_NUM()*6
    hist_loc = 0
!$OMP do
    do i = 1, 2000000000
        value = rand_r(seed)
        hist_loc(IAND(value, 15)) = hist_loc(IAND(value, 15)) + 1
    enddo
!$OMP end do
!$OMP critical
    hist_crit = hist_crit + hist_loc
!$OMP end critical
!$OMP end parallel

    call cpu_time(t_end)
!$  t_end = OMP_get_wtime()
    do i = 0, 15
        write (*, *) "hist_red (", i, ")", hist_red(i), "hist_crit (", i, ")", &
            hist_crit(i)
    enddo
    write (*, '(a20, f7.4, a20, E15.5)') 'Time: ', t_end - t_start

end program histo
