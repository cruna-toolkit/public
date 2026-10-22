module diff2_i1_mpi_cartesian

  use parameter
  use discretisation_x_core   , only: diff2_core => diff2_i1 , diff2_i1_name => diff2_i1_name , diff2_i1_type => diff2_i1_type
  use transpose_mpi_cartesian , only: ftranspose => ftranspose_i1, btranspose => btranspose_i1

  interface diff2
     module procedure diff2_1d
     module procedure diff2_2d
     module procedure diff2_3d
     module procedure diff2_4d
  end interface diff2

contains

!!!=================================================================================================
  subroutine diff2_1d(q_x,q,h)
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
    call diff2_core(fries_2d,dx)
    call btranspose(q_3d,fries_2d)

    q_x = q_3d(:,1,1)

  end subroutine diff2_1d
!!!=================================================================================================



!!!=================================================================================================
  subroutine diff2_2d(q_x,q,h)
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
    call diff2_core(fries_2d,dx)
    call btranspose(q_3d,fries_2d)

    q_x = q_3d(:,:,1)

  end subroutine diff2_2d
!!!=================================================================================================



!!!=================================================================================================
  subroutine diff2_3d(q_x,q,h)
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
    call diff2_core(fries_2d,dx)
    call btranspose(q_x,fries_2d)

  end subroutine diff2_3d
!!!=================================================================================================



!!!=================================================================================================
  subroutine diff2_4d(q_x,q,h)
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
       call diff2_core(fries_2d,dx)
       call btranspose(q_x(:,:,:,i),fries_2d)
    end do

  end subroutine diff2_4d
!!!=================================================================================================


end module diff2_i1_mpi_cartesian
