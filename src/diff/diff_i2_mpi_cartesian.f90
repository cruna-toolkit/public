module diff_i2_mpi_cartesian

  use parameter
  use discretisation_x_core   , only: diff_core  => diff_i2 , diff_i2_name => diff_i2_name , diff_i2_type => diff_i2_type
  use transpose_mpi_cartesian , only: ftranspose => ftranspose_i2, btranspose => btranspose_i2

  interface diff
     module procedure diff_1d
     module procedure diff_2d
     module procedure diff_3d
     module procedure diff_4d
  end interface diff

contains

!!!=================================================================================================
  subroutine diff_1d(q_x,q,h)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:), intent(out)          :: q_x
    real(kind=rk), dimension(:), intent(in)           :: q
    real(kind=rk)              , intent(in),optional  :: h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), allocatable      :: q_3d
    real(kind=rk), dimension(:,:), allocatable        :: fries_2d
    real(kind=rk)                                     :: dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(present(h)) then
       dx = h
    else
       dx = 1.0_rk
    end if

    allocate(q_3d(size(q,1),1,1))

    q_3d(:,1,1) = q

    call ftranspose(fries_2d,q_3d)
    call diff_core(fries_2d,dx)
    call btranspose(q_3d,fries_2d)

    q_x = q_3d(:,1,1)

  end subroutine diff_1d
!!!=================================================================================================



!!!=================================================================================================
  subroutine diff_2d(q_x,q,h)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(out)          :: q_x
    real(kind=rk), dimension(:,:), intent(in)           :: q
    real(kind=rk)                , intent(in),optional  :: h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), allocatable        :: q_3d
    real(kind=rk), dimension(:,:), allocatable          :: fries_2d
    real(kind=rk)                                       :: dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(present(h)) then
       dx = h
    else
       dx = 1.0_rk
    end if

    allocate(q_3d(size(q,1),size(q,2),1))

    q_3d(:,:,1) = q

    call ftranspose(fries_2d,q_3d)
    call diff_core(fries_2d,dx)
    call btranspose(q_3d,fries_2d)

    q_x = q_3d(:,:,1)

  end subroutine diff_2d
!!!=================================================================================================



!!!=================================================================================================
  subroutine diff_3d(q_x,q,h)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(out)          :: q_x
    real(kind=rk), dimension(:,:,:), intent(in)           :: q
    real(kind=rk)                  , intent(in),optional  :: h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), allocatable            :: fries_2d
    real(kind=rk)                                         :: dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(present(h)) then
       dx = h
    else
       dx = 1.0_rk
    end if

    call ftranspose(fries_2d,q)
    call diff_core(fries_2d,dx)
    call btranspose(q_x,fries_2d)

  end subroutine diff_3d
!!!=================================================================================================



!!!=================================================================================================
  subroutine diff_4d(q_x,q,h)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(out)          :: q_x
    real(kind=rk), dimension(:,:,:,:), intent(in)           :: q
    real(kind=rk)                    , intent(in),optional  :: h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), allocatable              :: fries_2d
    real(kind=rk)                                           :: dx
    integer                                                 :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(present(h)) then
       dx = h
    else
       dx = 1.0_rk
    end if

    do i = 1,size(q,4)
       call ftranspose(fries_2d,q(:,:,:,i))
       call diff_core(fries_2d,dx)
       call btranspose(q_x(:,:,:,i),fries_2d)
    end do

  end subroutine diff_4d
!!!=================================================================================================


end module diff_i2_mpi_cartesian
