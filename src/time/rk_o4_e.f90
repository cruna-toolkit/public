module rk_o4_e

  use governing_equations
  use parameter

  private

  public time_stepper_direct
  public time_stepper_adjoint

  character(len=max_length_parameter), parameter, public :: time_stepper_direct_name  = 'rk_o4_e'
  character(len=max_length_parameter), parameter, public :: time_stepper_adjoint_name = 'rk_o4_e'

contains

!!!=================================================================================================
  subroutine time_stepper_direct(q1,q0,t)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(out)          :: q1
    real(kind=rk), dimension(:,:,:,:), intent(in)           :: q0
    real(kind=rk), intent(in)        , optional             :: t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable          :: k1,k2,k3,k4
    real(kind=rk)                                           :: time
    real(kind=rk)                                           :: dt, dth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(t)) then
       time = t
    else
       time = params%time%t
    end if

    dt  = params%time%dt
    dth = params%time%dt*(0.5_rk)

    ! allocate subsets
    allocate(k1(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    allocate(k2(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    allocate(k3(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    allocate(k4(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))

    ! stepping
    call rhs_direct(k1, q0          , time       )
    call rhs_direct(k2, q0 + dth*k1 , time + dth )
    call rhs_direct(k3, q0 + dth*k2 , time + dth )
    call rhs_direct(k4, q0 + dt *k3 , time + dt  )

    q1 = q0 + dt/6.0_rk * (k1 + 2*k2 + 2*k3 + k4)

    ! deallocate
    deallocate(k1)
    deallocate(k2)
    deallocate(k3)
    deallocate(k4)

  end subroutine time_stepper_direct
!!!=================================================================================================

!!!=================================================================================================
  subroutine time_stepper_adjoint(qs1,qs0,q0,g,t)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(out)          :: qs1
    real(kind=rk), dimension(:,:,:,:), intent(in)           :: qs0
    real(kind=rk), dimension(:,:,:,:), intent(in)           :: q0
    real(kind=rk), dimension(:,:,:,:), intent(in)           :: g
    real(kind=rk), intent(in)        , optional             :: t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable          :: k1,k2,k3,k4
    real(kind=rk)                                           :: time
    real(kind=rk)                                           :: dt, dth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(t)) then
       time = t
    else
       time = params%time%t
    end if

    dt  = -params%time%dt
    dth = -params%time%dt*(0.5_rk)

    ! allocate subsets
    allocate(k1(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    allocate(k2(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    allocate(k3(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    allocate(k4(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))

    ! stepping
    call rhs_adjoint(k1, qs0          , q0, g, time       )
    call rhs_adjoint(k2, qs0 + dth*k1 , q0, g, time + dth )
    call rhs_adjoint(k3, qs0 + dth*k2 , q0, g, time + dth )
    call rhs_adjoint(k4, qs0 + dt *k3 , q0, g, time + dt  )

    qs1 = qs0 + dt/6.0_rk * (k1 + 2*k2 + 2*k3 + k4)

    ! deallocate
    deallocate(k1)
    deallocate(k2)
    deallocate(k3)
    deallocate(k4)

  end subroutine time_stepper_adjoint
!!!=================================================================================================

end module rk_o4_e
