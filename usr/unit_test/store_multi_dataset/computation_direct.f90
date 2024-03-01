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
    real(kind=rk), dimension(:,:,:,:,:), allocatable   :: field5d,field5d_in
    real(kind=rk)                                      :: max_err,max_err_all,max_err_sum
    integer                                            :: n1,n2,n3,n4,n5,i,ierr
    integer                            , parameter     :: loops = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%world_image.eq.1) then
       write(*,*) "========================"
       write(*,*) "STORE MULTI DATASET TEST"
       write(*,*) "========================"
    end if

    !! allocates
    n1 = size(Q,1)
    n2 = size(Q,2)
    n3 = size(Q,3)
    n4 = size(Q,4)
    n5 = size(Q,5)

    allocate(   field5d(n1,n2,n3,n4,n5))

    ! fill field
    call random_number(field5d)

!!! 5d field
    ! default dset=/data
    !call store(field5d,'debug_field')
    call store(field5d,'debug_field',dset_name_optin='test')
    !call store(field5d,'debug_field',dset_name_optin='test2')

    call stop_cruna

  end subroutine calc_direct
!!!=================================================================================================

end module computation_direct
