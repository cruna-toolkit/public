module filter_weighted_pen_core_o2_e_c_p

  use parameter

  private

  public :: init_filter
  public :: filter_i1, filter_i2, filter_i3

  character(len=max_length_parameter), parameter, public :: filter_type = 'periodic'
  character(len=max_length_parameter), parameter, public :: filter_name = 'filter_pen_core_o2_e_c_p'

  real(kind=rk), dimension(2)        , save              :: stencil

contains

!!!=================================================================================================
  subroutine init_filter(strength,weight,q0,penalization)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:)  , intent(out)                     :: strength
    real(kind=rk), dimension(:,:,:)  , intent(out)                     :: weight
    real(kind=rk), dimension(:,:,:,:), intent(in) , optional           :: q0, penalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! strength
    strength = 1.0_rk

    ! weight
    if(present(penalization)) then
       weight = penalization(:,:,:,1) ! phi
    else
       weight = 1.0_rk
    end if

    call init_stencil(stencil)

  end subroutine init_filter
!!!==================================================================================================

!!!=================================================================================================
  subroutine filter_i1(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call filter_main(q)

  end subroutine filter_i1
!!!==================================================================================================

!!!=================================================================================================
  subroutine filter_i2(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call filter_main(q)

  end subroutine filter_i2
!!!==================================================================================================

!!!=================================================================================================
  subroutine filter_i3(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call filter_main(q)

  end subroutine filter_i3
!!!==================================================================================================


!!!=================================================================================================
  subroutine filter_main(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(:,:), intent(inout)       :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(:,:), allocatable         :: q_tmp
    real(kind=rk)                                       :: a0,a1
    integer                                             :: n1,n2
    integer                                             :: i,j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n1 = size(q,1)
    n2 = size(q,2)

    allocate(q_tmp(n1,n2))

    ! unpack stencil for readability
    a0  = stencil(1)
    a1  = stencil(2)

    ! compute A*U
    do j = 1,n2
       do i = 2,(n1-1)
          q_tmp(i,j) = a0*q(i,j) + a1*(q(i+1,j) + q(i-1,j))
       end do
    end do

    ! compute A*U (boundary)
    q_tmp( 1,:) = a0*q(1 ,:) + a1*(q(2,:) + q(n1  ,:)) 
    q_tmp(n1,:) = a0*q(n1,:) + a1*(q(1,:) + q(n1-1,:))

    ! write back
    q = q_tmp

  end subroutine filter_main
!!!==================================================================================================

!!!==================================================================================================
  subroutine init_stencil(stencil) ! output args
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Init. Values from
    ! Reiss, Julius. (2021). Pressure-tight and non-stiff volume penalization for compressible flows.
    ! https://arxiv.org/abs/2103.08144
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(2)    , intent(out)       :: stencil
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                       :: a0,a1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    a0      = -2.0_rk*0.25_rk
    a1      =  1.0_rk*0.25_rk

    stencil = (/a0, a1/)

  end subroutine init_stencil
!!!==================================================================================================

end module filter_weighted_pen_core_o2_e_c_p
