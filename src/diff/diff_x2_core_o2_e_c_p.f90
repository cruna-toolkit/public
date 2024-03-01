module diff_x2_core_o2_e_c_p

  use parameter

  private

  public :: diff

  character(len=max_length_parameter), parameter, public :: diff_type = 'periodic'
  character(len=max_length_parameter), parameter, public :: diff_name = 'diff_x2_core_o2_e_c_p'

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
    real(kind=rk)                                     :: dx
    integer                                           :: i,n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(present(h)) then
       dx = h
    else
       dx = 1.0_rk
    end if

    n = size(q)

    do i = 2,n-1
       q_x(i) = q(i+1) - q(i-1)
    end do

    q_x(1)= q(2) - q(n  )
    q_x(n)= q(1) - q(n-1)

    q_x = q_x/(2.0_rk*dx)

  end subroutine diff_1d
!!!=================================================================================================



!!!=================================================================================================
  subroutine diff_2d(q_x,q,h)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(out)          :: q_x
    real(kind=rk), dimension(:,:), intent(in)           :: q
    real(kind=rk)                , intent(in),optional  :: h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                       :: dx
    integer                                             :: i,n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(present(h)) then
       dx = h
    else
       dx = 1.0_rk
    end if

    n = size(q,2)

    do i = 2,n-1
       q_x(:,i) = q(:,i+1) - q(:,i-1)
    end do

    q_x(:,1)= q(:,2) - q(:,n  )
    q_x(:,n)= q(:,1) - q(:,n-1)

    q_x = q_x/(2.0_rk*dx)

  end subroutine diff_2d
!!!=================================================================================================



!!!=================================================================================================
  subroutine diff_3d(q_x,q,h)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(out)          :: q_x
    real(kind=rk), dimension(:,:,:), intent(in)           :: q
    real(kind=rk)                  , intent(in),optional  :: h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                         :: dx
    integer                                               :: i,n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(present(h)) then
       dx = h
    else
       dx = 1.0_rk
    end if

    n = size(q,2)

    do i = 2,n-1
       q_x(:,i,:) = q(:,i+1,:) - q(:,i-1,:)
    end do

    q_x(:,1,:)= q(:,2,:) - q(:,n  ,:)
    q_x(:,n,:)= q(:,1,:) - q(:,n-1,:)

    q_x = q_x/(2.0_rk*dx)

  end subroutine diff_3d
!!!=================================================================================================



!!!=================================================================================================
  subroutine diff_4d(q_x,q,h)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(out)          :: q_x
    real(kind=rk), dimension(:,:,:,:), intent(in)           :: q
    real(kind=rk)                    , intent(in),optional  :: h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                           :: dx
    integer                                                 :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(present(h)) then
       dx = h
    else
       dx = 1.0_rk
    end if

    do i = 1,size(q,4)
       call diff_3d(q_x(:,:,:,i),q(:,:,:,i),dx)
    end do

  end subroutine diff_4d
!!!=================================================================================================

end module diff_x2_core_o2_e_c_p
