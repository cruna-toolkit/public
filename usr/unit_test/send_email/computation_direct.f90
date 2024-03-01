module computation_direct

  use helper
  use parallelism
  use parameter

contains

!!!=================================================================================================
  subroutine calc_direct(Q,q0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:), intent(inout) :: Q
    real(kind=rk), dimension(:,:,:,:)  , intent(in)    :: q0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:)      , allocatable   :: rk_2d_array                     
    integer                                            :: i,ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%world_image.eq.1) then
       write(*,*) "==============="
       write(*,*) "SEND_EMAIL_TEST"
       write(*,*) "==============="
    end if

    ! text mail
    call send_mail("Test 1 - string","Test Text 1")

    ! rk 2d array
    allocate(rk_2d_array(5,2))
    do i = 1,5
       rk_2d_array(i,1) = real(i,rk)
       rk_2d_array(i,2) = real(i,rk)*10
    end do
    call send_mail("Test 2 - rk 2d array",transpose(rk_2d_array))

    call mpi_barrier(params%parallelism%world_comm,ierr)
    call stop_cruna

  end subroutine calc_direct
!!!=================================================================================================

end module computation_direct
