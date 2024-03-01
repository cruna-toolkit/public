module computation_direct

  use filter_x
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
    real(kind=rk), dimension(:,:,:,:)  , allocatable   :: q1,q1f
    real(kind=rk), dimension(:,:), allocatable         :: fries_2d
    real(kind=rk)                                      :: filter_max_err,filter_max_err_all,filter_max_err_sum
    integer                                            :: i,l
    integer                                            :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%world_image.eq.1) then
       write(*,*) "============================"
       write(*,*) "FILTER_TEST"
       write(*,*) "============================"
    end if

    !! allocates
    allocate( q1(size(q0,1),size(q0,2),size(q0,3),size(q0,4)))
    allocate(q1f(size(q0,1),size(q0,2),size(q0,3),size(q0,4)))

    !! consistency check
    if (params%parallelism%world_image.eq.1) then
       write(*,*) "consistency: (q(:,:,:,:) = 1.0_rk)"
       write(*,*) ""
    end if

    q1 = 1.0_rk

    ! i1-filter
    q1f = q1    

    call filter_i1(q1f)

    filter_max_err = maxval(abs(q1f - q1))
    call mpi_reduce(filter_max_err,filter_max_err_all,1,PARAMS%PARALLELISM%MPI_TYPE_RK,MPI_MAX,0,params%parallelism%block_comm,ierr)
    call mpi_reduce(filter_max_err,filter_max_err_sum,1,PARAMS%PARALLELISM%MPI_TYPE_RK,MPI_SUM,0,params%parallelism%block_comm,ierr)
    if(params%parallelism%block_image.eq.1) then
       write(*,*) "i1-filter errors (max,sum,block): ",filter_max_err,filter_max_err_sum,params%parallelism%world_block
    end if

    ! i2-filter
    q1f = q1    
    call filter_i2(q1f)

    filter_max_err = maxval(abs(q1f - q1))
    call mpi_reduce(filter_max_err,filter_max_err_all,1,PARAMS%PARALLELISM%MPI_TYPE_RK,MPI_MAX,0,params%parallelism%block_comm,ierr)
    call mpi_reduce(filter_max_err,filter_max_err_sum,1,PARAMS%PARALLELISM%MPI_TYPE_RK,MPI_SUM,0,params%parallelism%block_comm,ierr)
    if(params%parallelism%block_image.eq.1) then
       write(*,*) "i2-filter errors (max,sum,block): ",filter_max_err,filter_max_err_sum,params%parallelism%world_block
    end if

    ! i3-filter
    q1f = q1
    call filter_i3(q1f)

    filter_max_err = maxval(abs(q1f - q1))
    call mpi_reduce(filter_max_err,filter_max_err_all,1,PARAMS%PARALLELISM%MPI_TYPE_RK,MPI_MAX,0,params%parallelism%block_comm,ierr)
    call mpi_reduce(filter_max_err,filter_max_err_sum,1,PARAMS%PARALLELISM%MPI_TYPE_RK,MPI_SUM,0,params%parallelism%block_comm,ierr)
    if(params%parallelism%block_image.eq.1) then
       write(*,*) "i3-filter errors (max,sum,block): ",filter_max_err,filter_max_err_sum,params%parallelism%world_block
    end if

    call mpi_barrier(params%parallelism%world_comm,ierr)
    call stop_cruna

  end subroutine calc_direct
!!!=================================================================================================

end module computation_direct
