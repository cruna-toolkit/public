module computation_direct

  use helper
  use io
  use io_wrapper
  USE mpi
  use parallelism
  use parameter

contains

!!!=================================================================================================
  subroutine calc_direct(Q,q0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:), intent(inout) :: Q
    real(kind=rk), dimension(:,:,:,:)  , intent(in)    :: q0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                            , parameter     :: loops = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:), allocatable   :: field5d,field5d_in
    real(kind=rk)                                      :: max_err,max_err_all,max_err_sum
    real(kind=rk)                                      :: tocs_w(loops),tocs_l(loops), tocs(2)
    integer                                            :: d,i,ierr
    integer                                            :: n(5),nc(5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%world_image.eq.1) then
       write(*,*) "==============="
       write(*,*) "STORE/LOAD TEST"
       write(*,*) "==============="
    end if

    !! allocates
    n(1) = size(Q,1)
    n(2) = size(Q,2)
    n(3) = size(Q,3)
    n(4) = size(Q,4)
    n(5) = size(Q,5)
    nc   = n

    allocate(   field5d(n(1),n(2),n(3),n(4),n(5)))
    allocate(field5d_in(n(1),n(2),n(3),n(4),n(5)))

    ! fill field
    call random_number(field5d)

    ! cycles
    do d = 1,5

       nc        = n
       nc(d+1:5) = 1
       
       if(params%parallelism%world_image.eq.1) then
          write(*,*) ""
          write(*,*) "test store/load of ",trim(num2str(d,'(I1.1)')),"d fields:"
          write(*,*) nc
          write(*,*) ""
       end if

       do i = 1,loops

          ! write field
          call tic()
          call store(field5d(1:nc(1),1:nc(2),1:nc(3),1:nc(4),1:nc(5)),'debug_field')
          tocs_w(i) = toc()

          ! load field
          call tic()
          call load(field5d_in(1:nc(1),1:nc(2),1:nc(3),1:nc(4),1:nc(5)),'debug_field')
          tocs_l(i) = toc()

          ! mean tocs
          tocs(1) = real(sum(tocs_w),rk)/real(loops,rk)
          tocs(2) = real(sum(tocs_l),rk)/real(loops,rk)

          ! compare field
          max_err = maxval(abs(field5d_in(1:nc(1),1:nc(2),1:nc(3),1:nc(4),1:nc(5)) - field5d(1:nc(1),1:nc(2),1:nc(3),1:nc(4),1:nc(5))))
          call mpi_reduce(max_err,max_err_all,1,params%parallelism%mpi_type_rk,mpi_max,0,params%parallelism%block_comm,ierr)
          call mpi_reduce(max_err,max_err_sum,1,params%parallelism%mpi_type_rk,mpi_sum,0,params%parallelism%block_comm,ierr)

          ! output
          if(params%parallelism%world_image.eq.1) then
             write(*,*) "   times in loop ",trim(num2str(i,'(I3.3)'))," (write): ",tocs(1)
             write(*,*) "   times in loop ",trim(num2str(i,'(I3.3)')),"  (load): ",tocs(2)
             write(*,*) "   error in loop ",trim(num2str(i,'(I3.3)')),"   (max): ",max_err_all
             write(*,*) "   error in loop ",trim(num2str(i,'(I3.3)')),"   (sum): ",max_err_sum
          end if

       end do
    end do

    call stop_cruna

  end subroutine calc_direct
!!!=================================================================================================

end module computation_direct
