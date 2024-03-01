module objective

  use helper
  use parallelism
  use parameter

  private

  real(kind=rk),dimension(:,:),allocatable,save, public :: objective_function

  public :: calc_objective
  public :: calc_objective_g


contains


!!!=================================================================================================
  subroutine calc_objective(J,Q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                     ,intent(out)                :: J
    real(kind=rk),dimension(:,:,:,:,:),intent(in)                 :: Q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                       :: loop,maxloop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! allocate objective function
    if (.not.allocated(objective_function)) then
       call get_parameter(maxloop,'opt.maxloop',default = 1)
       allocate(objective_function(maxloop,1))
    end if

    ! do nothing, if you need more please use over-ride functionality of the makefile
    call get_parameter(loop,'opt.loop',default = 1)
    objective_function(loop,:) = -1.0_rk
    J = 0.0_rk

  end subroutine calc_objective
!!!=================================================================================================

!!!=================================================================================================
  subroutine calc_objective_g(g,q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:),intent(out)                :: g
    real(kind=rk),dimension(:,:,:,:),intent(in) , optional      :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    g = 0.0_rk

  end subroutine calc_objective_g
!!!=================================================================================================

end module objective
