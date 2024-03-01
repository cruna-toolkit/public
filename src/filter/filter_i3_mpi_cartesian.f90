module filter_i3_mpi_cartesian

  use parameter
  use filter_x_core           , only: filter_core  => filter_i3 , filter_i3_name => filter_i3_name , filter_i3_type => filter_i3_type
  use transpose_mpi_cartesian , only: ftranspose => ftranspose_i3, btranspose => btranspose_i3

  interface filter
     module procedure filter_1d
     module procedure filter_2d
     module procedure filter_3d
     module procedure filter_4d
  end interface filter

contains

!!!=================================================================================================
  subroutine filter_1d(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:), intent(inout)        :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), allocatable      :: q_3d
    real(kind=rk), dimension(:,:), allocatable        :: fries_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(q_3d(size(q,1),1,1))

    q_3d(:,1,1) = q

    call ftranspose(fries_2d,q_3d)
    call filter_core(fries_2d)
    call btranspose(q_3d,fries_2d)

    q = q_3d(:,1,1)

  end subroutine filter_1d
!!!=================================================================================================



!!!=================================================================================================
  subroutine filter_2d(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)        :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), allocatable        :: q_3d
    real(kind=rk), dimension(:,:), allocatable          :: fries_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(q_3d(size(q,1),size(q,2),1))

    q_3d(:,:,1) = q

    call ftranspose(fries_2d,q_3d)
    call filter_core(fries_2d)
    call btranspose(q_3d,fries_2d)

    q = q_3d(:,:,1)

  end subroutine filter_2d
!!!=================================================================================================



!!!=================================================================================================
  subroutine filter_3d(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(inout)        :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), allocatable            :: fries_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call ftranspose(fries_2d,q)
    call filter_core(fries_2d)
    call btranspose(q,fries_2d)

  end subroutine filter_3d
!!!=================================================================================================



!!!=================================================================================================
  subroutine filter_4d(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)        :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), allocatable              :: fries_2d
    integer                                                 :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1,size(q,4)
       call ftranspose(fries_2d,q(:,:,:,i))
       call filter_core(fries_2d)
       call btranspose(q(:,:,:,i),fries_2d)
    end do

  end subroutine filter_4d
!!!=================================================================================================

end module filter_i3_mpi_cartesian
