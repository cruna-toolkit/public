module diff_pen_core_o6_e_c_p

  use parameter

  private

  public :: diff_i1, diff_i2, diff_i3

  character(len=max_length_parameter), parameter, public :: diff_type = 'periodic'
  character(len=max_length_parameter), parameter, public :: diff_name = 'diff_pen_core_o6_e_c_p'

contains

!!!=================================================================================================
  subroutine diff_i1(q,dx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)      :: q
    real(kind=rk)                , intent(in)         :: dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call diff_main(q,dx)

  end subroutine diff_i1
!!!==================================================================================================

!!!=================================================================================================
  subroutine diff_i2(q,dx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
    real(kind=rk)                , intent(in)          :: dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call diff_main(q,dx)

  end subroutine diff_i2
!!!==================================================================================================

!!!=================================================================================================
  subroutine diff_i3(q,dx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
    real(kind=rk)                , intent(in)          :: dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call diff_main(q,dx)

  end subroutine diff_i3
!!!==================================================================================================


!!!=================================================================================================
  subroutine diff_main(q,dx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! 6th order explicit, central stencil (periodic)
!!! Interior: 7-point central  (-1, 9, -45, 0, 45, -9, 1) / 60
!!! Periodic wrap handled explicitly for i=1,2,3 and i=n1-2,n1-1,n1
!!! Source: Fornberg, B. 1988,Generation of finite difference formulas on arbitrarily spaced grids
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(:,:), intent(inout)       :: q
    real(kind=rk)                 , intent(in)          :: dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(:,:), allocatable         :: q_tmp
    integer                                             :: n1,n2
    integer                                             :: i,j

    n1 = size(q,1)
    n2 = size(q,2)

    allocate(q_tmp(n1,n2))

    ! compute A*U (periodic, 6th order)
    ! Left wrap block
    q_tmp(1   ,:) = -1.0_rk*q(n1-2,:) +  9.0_rk*q(n1-1,:) - 45.0_rk*q(n1  ,:) &
                    +45.0_rk*q(2     ,:) -  9.0_rk*q(3   ,:) +  1.0_rk*q(4   ,:)
    q_tmp(2   ,:) = -1.0_rk*q(n1-1,:) +  9.0_rk*q(n1  ,:) - 45.0_rk*q(1   ,:) &
                    +45.0_rk*q(3     ,:) -  9.0_rk*q(4   ,:) +  1.0_rk*q(5   ,:)
    q_tmp(3   ,:) = -1.0_rk*q(n1  ,:) +  9.0_rk*q(1   ,:) - 45.0_rk*q(2   ,:) &
                    +45.0_rk*q(4     ,:) -  9.0_rk*q(5   ,:) +  1.0_rk*q(6   ,:)

    ! Interior
    do j = 1,n2
       do i = 4,(n1-3)
          q_tmp(i,j) = -1.0_rk*q(i-3,j) +  9.0_rk*q(i-2,j) - 45.0_rk*q(i-1,j) &
                        +45.0_rk*q(i+1,j) -  9.0_rk*q(i+2,j) +  1.0_rk*q(i+3,j)
       end do
    end do

    ! Right wrap block (mirror)
    q_tmp(n1-2,:) = -1.0_rk*q(n1-5,:) +  9.0_rk*q(n1-4,:) - 45.0_rk*q(n1-3,:) &
                     +45.0_rk*q(n1-1,:) -  9.0_rk*q(n1  ,:) +  1.0_rk*q(1     ,:)
    q_tmp(n1-1,:) = -1.0_rk*q(n1-4,:) +  9.0_rk*q(n1-3,:) - 45.0_rk*q(n1-2,:) &
                     +45.0_rk*q(n1  ,:) -  9.0_rk*q(1     ,:) +  1.0_rk*q(2     ,:)
    q_tmp(n1  ,:) = -1.0_rk*q(n1-3,:) +  9.0_rk*q(n1-2,:) - 45.0_rk*q(n1-1,:) &
                     +45.0_rk*q(1     ,:) -  9.0_rk*q(2     ,:) +  1.0_rk*q(3     ,:)

    ! dx factor
    q = q_tmp/(dx*60.0_rk)

  end subroutine diff_main
!!!==================================================================================================

end module diff_pen_core_o6_e_c_p
