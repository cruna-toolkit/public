module data_reference

  use parameter

  real(kind=rk), dimension(:,:,:,:) , allocatable, save, public :: qref

  private

  public init_reference

contains

!!!=================================================================================================
  subroutine init_reference(q0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind = rk), dimension(:,:,:,:), intent(in)                     :: q0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(.not.allocated(qref)) then
       allocate(qref(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
       qref = 0.0_rk
    end if

    qref = q0

  end subroutine init_reference
!!!=================================================================================================

end module data_reference
