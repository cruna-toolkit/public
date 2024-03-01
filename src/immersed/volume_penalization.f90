module volume_penalization

  use parameter

  real(kind=rk), dimension(:,:,:,:), allocatable, save, public :: penalization

  private

  public allocate_volume_penalization
  public apply_volume_penalization
  public backup_volume_penalization
  public calc_grad_volume_penalization
  public init_volume_penalization
  public store_volume_penalization
  public restore_volume_penalization
  public update_volume_penalization

contains

!!!=================================================================================================
  subroutine allocate_volume_penalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a forcing please use over-ride functionality of the makefile
  end subroutine allocate_volume_penalization
!!!=================================================================================================

!!!=================================================================================================
  subroutine apply_volume_penalization(rhs,q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)                    :: rhs
    real(kind=rk), dimension(:,:,:,:), intent(in)                       :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a forcing please use over-ride functionality of the makefile
    ! rhs = rhs
  end subroutine apply_volume_penalization
!!!=================================================================================================

!!!=================================================================================================
  subroutine backup_volume_penalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  end subroutine backup_volume_penalization
!!!=================================================================================================

!!!=================================================================================================
  subroutine calc_grad_volume_penalization(Qs,Q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:),intent(in),optional :: Qs
    real(kind=rk), dimension(:,:,:,:,:),intent(in),optional :: Q      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  end subroutine calc_grad_volume_penalization
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_volume_penalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  end subroutine init_volume_penalization
!!!=================================================================================================

!!!=================================================================================================
  subroutine store_volume_penalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  end subroutine store_volume_penalization
!!!=================================================================================================

!!!=================================================================================================
  subroutine restore_volume_penalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  end subroutine restore_volume_penalization
!!!=================================================================================================

!!!=================================================================================================
  subroutine update_volume_penalization(d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),intent(in),optional  :: d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  end subroutine update_volume_penalization
!!!=================================================================================================

end module volume_penalization
