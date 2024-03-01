module force

  use parameter

  real(kind=rk), dimension(:,:,:,:), allocatable, save, public :: f

  private

  public allocate_force
  public apply_force
  public backup_force
  public calc_grad
  public init_force
  public store_force
  public restore_force
  public update_force

contains

!!!=================================================================================================
  subroutine allocate_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a forcing please use over-ride functionality of the makefile
  end subroutine allocate_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine apply_force(rhs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)                    :: rhs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a forcing please use over-ride functionality of the makefile
    ! rhs = rhs
  end subroutine apply_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine backup_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a forcing please use over-ride functionality of the makefile
  end subroutine backup_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine calc_grad(Qs,Q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:),intent(in),optional :: Qs
    real(kind=rk), dimension(:,:,:,:,:),intent(in),optional :: Q      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a forcing please use over-ride functionality of the makefile
  end subroutine calc_grad
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a forcing please use over-ride functionality of the makefile
  end subroutine init_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine store_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a forcing please use over-ride functionality of the makefile
  end subroutine store_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine restore_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a forcing please use over-ride functionality of the makefile
  end subroutine restore_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine update_force(d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),intent(in),optional  :: d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a forcing please use over-ride functionality of the makefile
  end subroutine update_force
!!!=================================================================================================

end module force
