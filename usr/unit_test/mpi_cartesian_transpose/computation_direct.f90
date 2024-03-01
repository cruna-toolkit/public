module computation_direct

  use helper
  use io
  USE mpi
  use parallelism
  use parameter
  use transpose_mpi_cartesian

contains

!!!=================================================================================================
  subroutine calc_direct(Q,q0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:), intent(inout) :: Q
    real(kind=rk), dimension(:,:,:,:)  , intent(in)    :: q0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:)  , allocatable   :: q0b
    real(kind=rk), dimension(:,:), allocatable         :: fries_2d
    real(kind=rk)                                      :: transpose_max_err,transpose_max_err_all,transpose_max_err_sum
    integer                                            :: i,l
    integer                                            :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real   , dimension(:)              , allocatable   :: tocs
    real                                               :: tocs_min,tocs_max,tocs_avg,tocs_min_all,tocs_max_all,tocs_avg_all
    integer, parameter                                 :: t_loops = 10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%world_image.eq.1) then
       write(*,*) "============================"
       write(*,*) "TRANSPOSE_MPI_CARTESIAN_TEST"
       write(*,*) "============================"
    end if

    !! allocates
    allocate(q0b(size(q0,1),size(q0,2),size(q0,3),size(q0,4)))

    !! consistency check
    if (params%parallelism%world_image.eq.1) then
       write(*,*) "consistency: (requires suitable inital condition)"
       write(*,*) ""
    end if

    ! i1-pens
    do l = 1,params%equation%nbr_vars
       call ftranspose_i1(fries_2d, q0(:,:,:,l))
       call btranspose_i1(q0b(:,:,:,l),fries_2d)
    end do

    transpose_max_err = maxval(abs(q0b - q0))
    call mpi_reduce(transpose_max_err,transpose_max_err_all,1,PARAMS%PARALLELISM%MPI_TYPE_RK,MPI_MAX,0,params%parallelism%block_comm,ierr)
    call mpi_reduce(transpose_max_err,transpose_max_err_sum,1,PARAMS%PARALLELISM%MPI_TYPE_RK,MPI_SUM,0,params%parallelism%block_comm,ierr)
    if(params%parallelism%block_image.eq.1) then
       write(*,*) "i1-transpose errors (max,sum,block): ",transpose_max_err,transpose_max_err_sum,params%parallelism%world_block
    end if

    ! i2-pens
    do l = 1,params%equation%nbr_vars
       call ftranspose_i2(fries_2d, q0(:,:,:,l))
       call btranspose_i2(q0b(:,:,:,l),fries_2d)
    end do

    transpose_max_err = maxval(abs(q0b - q0))
    call mpi_reduce(transpose_max_err,transpose_max_err_all,1,PARAMS%PARALLELISM%MPI_TYPE_RK,MPI_MAX,0,params%parallelism%block_comm,ierr)
    call mpi_reduce(transpose_max_err,transpose_max_err_sum,1,PARAMS%PARALLELISM%MPI_TYPE_RK,MPI_SUM,0,params%parallelism%block_comm,ierr)
    if(params%parallelism%block_image.eq.1) then
       write(*,*) "i2-transpose errors (max,sum,block): ",transpose_max_err,transpose_max_err_sum,params%parallelism%world_block
    end if

    ! i3-pens
    do l = 1,params%equation%nbr_vars
       call ftranspose_i3(fries_2d, q0(:,:,:,l))
       call btranspose_i3(q0b(:,:,:,l),fries_2d)
    end do

    transpose_max_err = maxval(abs(q0b - q0))
    call mpi_reduce(transpose_max_err,transpose_max_err_all,1,PARAMS%PARALLELISM%MPI_TYPE_RK,MPI_MAX,0,params%parallelism%block_comm,ierr)
    call mpi_reduce(transpose_max_err,transpose_max_err_sum,1,PARAMS%PARALLELISM%MPI_TYPE_RK,MPI_SUM,0,params%parallelism%block_comm,ierr)
    if(params%parallelism%block_image.eq.1) then
       write(*,*) "i3-transpose errors (max,sum,block): ",transpose_max_err,transpose_max_err_sum,params%parallelism%world_block
    end if

    !! scaling/performance
    if (params%parallelism%world_image.eq.1) then
       write(*,*) ""
       write(*,*) "performance: (parameters are hard coded in computation_direct.f90 in unit test case)"
       write(*,*) "  --> loops: ",t_loops
       write(*,*) ""
    end if

    allocate(tocs(t_loops))
    tocs = 0.0

    ! i1-pens
    do i=1,t_loops
       call tic()
       do l = 1,params%equation%nbr_vars
          call ftranspose_i1(fries_2d, q0(:,:,:,l))
          call btranspose_i1(q0b(:,:,:,l),fries_2d)
       end do
       tocs(i) = toc()
    end do

    tocs_min = minval(tocs)
    tocs_max = maxval(tocs)
    tocs_avg = real(sum(tocs)/t_loops)

    call mpi_reduce(tocs_min,tocs_min_all,1,MPI_REAL,MPI_MIN,0,params%parallelism%block_comm,ierr)
    call mpi_reduce(tocs_max,tocs_max_all,1,MPI_REAL,MPI_MAX,0,params%parallelism%block_comm,ierr)
    call mpi_reduce(tocs_avg,tocs_avg_all,1,MPI_REAL,MPI_SUM,0,params%parallelism%block_comm,ierr)
    tocs_avg_all = real(tocs_avg_all/params%parallelism%block_size)

    if(params%parallelism%block_image.eq.1) then
       write(*,*) "i1-transpose performance (min,max,avg,block): ",tocs_min_all,tocs_max_all,tocs_avg_all,params%parallelism%world_block
    end if

    ! i2-pens
    do i=1,t_loops
       call tic()
       do l = 1,params%equation%nbr_vars
          call ftranspose_i2(fries_2d, q0(:,:,:,l))
          call btranspose_i2(q0b(:,:,:,l),fries_2d)
       end do
       tocs(i) = toc()
    end do

    tocs_min = minval(tocs)
    tocs_max = maxval(tocs)
    tocs_avg = real(sum(tocs)/t_loops)

    call mpi_reduce(tocs_min,tocs_min_all,1,MPI_REAL,MPI_MIN,0,params%parallelism%block_comm,ierr)
    call mpi_reduce(tocs_max,tocs_max_all,1,MPI_REAL,MPI_MAX,0,params%parallelism%block_comm,ierr)
    call mpi_reduce(tocs_avg,tocs_avg_all,1,MPI_REAL,MPI_SUM,0,params%parallelism%block_comm,ierr)
    tocs_avg_all = real(tocs_avg_all/params%parallelism%block_size)

    if(params%parallelism%block_image.eq.1) then
       write(*,*) "i2-transpose performance (min,max,avg,block): ",tocs_min_all,tocs_max_all,tocs_avg_all,params%parallelism%world_block
    end if

    ! i3-pens
    do i=1,t_loops
       call tic()
       do l = 1,params%equation%nbr_vars
          call ftranspose_i3(fries_2d, q0(:,:,:,l))
          call btranspose_i3(q0b(:,:,:,l),fries_2d)
       end do
       tocs(i) = toc()
    end do

    tocs_min = minval(tocs)
    tocs_max = maxval(tocs)
    tocs_avg = real(sum(tocs)/t_loops)

    call mpi_reduce(tocs_min,tocs_min_all,1,MPI_REAL,MPI_MIN,0,params%parallelism%block_comm,ierr)
    call mpi_reduce(tocs_max,tocs_max_all,1,MPI_REAL,MPI_MAX,0,params%parallelism%block_comm,ierr)
    call mpi_reduce(tocs_avg,tocs_avg_all,1,MPI_REAL,MPI_SUM,0,params%parallelism%block_comm,ierr)
    tocs_avg_all = real(tocs_avg_all/params%parallelism%block_size)

    if(params%parallelism%block_image.eq.1) then
       write(*,*) "i3-transpose performance (min,max,avg,block): ",tocs_min_all,tocs_max_all,tocs_avg_all,params%parallelism%world_block
    end if

    call mpi_barrier(params%parallelism%world_comm,ierr)
    call stop_cruna

  end subroutine calc_direct
!!!=================================================================================================

end module computation_direct
