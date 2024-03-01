module filter_weighted_x

  use filter_weighted_i1_mpi_cartesian , only: filter_weighted_i1_local => filter , filter_weighted_i1_name => filter_i1_name , filter_weighted_i1_type => filter_i1_type
  use filter_weighted_i2_mpi_cartesian , only: filter_weighted_i2_local => filter , filter_weighted_i2_name => filter_i2_name , filter_weighted_i2_type => filter_i2_type
  use filter_weighted_i3_mpi_cartesian , only: filter_weighted_i3_local => filter , filter_weighted_i3_name => filter_i3_name , filter_weighted_i3_type => filter_i3_type
  use filter_weighted_x_core           , only: init_filter
  use parameter
  use volume_penalization              , only: penalization

  public                                             :: init_filter_weighted
  public                                             :: filter_weighted_i1
  public                                             :: filter_weighted_i2
  public                                             :: filter_weighted_i3

  real(kind=rk), dimension(:,:,:), allocatable, save :: strength_from_init
  real(kind=rk), dimension(:,:,:), allocatable, save :: weight_from_init

contains

!!!=================================================================================================
  subroutine init_filter_weighted(q0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(in)                       :: q0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! nothing to do

  end subroutine init_filter_weighted
!!!==================================================================================================

!!!=================================================================================================
  subroutine filter_weighted_i1(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)                       :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! nothing to do
    q = q

  end subroutine filter_weighted_i1
!!!==================================================================================================

!!!=================================================================================================
  subroutine filter_weighted_i2(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)                       :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! nothing to do
    q = q

  end subroutine filter_weighted_i2
!!!==================================================================================================

!!!=================================================================================================
  subroutine filter_weighted_i3(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)                       :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! nothing to do
    q = q

  end subroutine filter_weighted_i3
!!!==================================================================================================  

end module filter_weighted_x
