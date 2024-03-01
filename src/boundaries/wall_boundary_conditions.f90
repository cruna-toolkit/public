module wall_boundary_conditions

  use data_geom
  use parameter

  private

  public :: wall_slip_no_thermodynamics ! (100)

contains

!!!=================================================================================================
  subroutine wall_slip_no_thermodynamics(rhs,component)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:),intent(inout)                       :: rhs
    integer,intent(in)                                                   :: component
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    rhs(:,:,:,component) = 0

    ! print*,"wall_slip_no_thermodynamics u is component 2, v is component 3 ... now !!!!"

  end subroutine wall_slip_no_thermodynamics
!!!=================================================================================================

end module wall_boundary_conditions
