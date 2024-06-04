module linear_acoustic_p_u

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! information/remarks:                                                                             !
! - q_in, qs_in in subroutine calls, this is needed due to intent(in) of q/qs in time stepper,     !
!   thus no modification inside of the rhs is permitted, however no intent works fine, even        !
!   nobody understands why ... (2019-01-02, ML)                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use boundary_conditions
  use data_rhs , only: q0
  use discretisation_x
  use linear_acoustic_p_u_trafos
  use force
  use parameter
  use sponge_layer
  use io

  private

  public right_hand_side_direct_3d
  public right_hand_side_adjoint_3d

  public give_T, give_c, give_r, give_u, give_p

  character(len=max_length_parameter), parameter, public :: viscosity_type = "none"
  character(len=max_length_parameter), parameter, public :: equations_name = "linear_acoustic_p_u"

contains

!!!=================================================================================================
  subroutine right_hand_side_direct_3d(rhs,q_in,t)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(out)   :: rhs
    real(kind=rk), dimension(:,:,:,:), intent(in)    :: q_in
    real(kind=rk)                    , intent(in)    :: t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable   :: q
    real(kind=rk), dimension(:,:,:)  , allocatable   :: a_x, a_y, a_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: p_x, p_y, p_z
    real(kind=rk)                                    :: rho0, p0, c0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(q(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))

    allocate(a_x(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(a_y(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(a_z(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

    allocate(p_x(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(p_y(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(p_z(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))


!!! set q
    q = q_in

!!! boundary conditions 1/2: q
    call set_boundary_condition_q(q)
     
!!! set density of air at T= 293,15K and speed of sound
    rho0 = params%init%rho
    p0   = params%init%p
    c0   = sqrt(params%material%gamma * p0 / rho0)

!!! mass equation -------------------------------------------------------------------------------
    call Dx1(a_x,q(:,:,:,2),params%geom%dx1)
    call Dx2(a_y,q(:,:,:,3),params%geom%dx2)
    call Dx3(a_z,q(:,:,:,4),params%geom%dx3)
    rhs(:,:,:,1) = -rho0 * c0**2 * (a_x + a_y + a_z)


!!! momentum equation (1) ------------------------------------------------------------------------
    call Dx1(p_x,q(:,:,:,1),params%geom%dx1)

    rhs(:,:,:,2) = -p_x/rho0

!!! momentum equation (2) ------------------------------------------------------------------------
    call Dx2(p_y,q(:,:,:,1),params%geom%dx2)

    rhs(:,:,:,3) = -p_y/rho0

!!! momentum equation (3) ------------------------------------------------------------------------
    call Dx3(p_z,q(:,:,:,1),params%geom%dx3)

    rhs(:,:,:,4) = -p_z/rho0

!!! forcing (field)
    !call update_force
    !call apply_force(rhs)

!!! sponge layer
    call sponge(rhs,q,q0)

!!! boundary conditions 2/2: rhs
    call set_boundary_condition_rhs(rhs,q)

    deallocate(q)

    deallocate(a_x)
    deallocate(a_y)
    deallocate(a_z)

    deallocate(p_x)
    deallocate(p_y)
    deallocate(p_z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine right_hand_side_direct_3d
!!!=================================================================================================

!!!=================================================================================================
  subroutine right_hand_side_adjoint_3d(rhs,qs_in,q,g,t)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(out)   :: rhs
    real(kind=rk), dimension(:,:,:,:), intent(in)    :: qs_in
    real(kind=rk), dimension(:,:,:,:), intent(in)    :: q
    real(kind=rk), dimension(:,:,:,:), intent(in)    :: g
    real(kind=rk)                    , intent(in)    :: t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable   :: qs
    real(kind=rk)                                    :: rho0, p0, c0
    real(kind=rk), dimension(:,:,:)  , allocatable   :: a_x,a_y,a_z,p_x,p_y,p_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(qs(params%geom%n1b,params%geom%n2b,params%geom%n3b,5))

    allocate(a_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(a_y(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(a_z(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(p_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(p_y(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(p_z(params%geom%n1b,params%geom%n2b,params%geom%n3b))


!!! set qs
    qs = qs_in

!!! boundary conditions 1/2: qs
    call set_boundary_condition_qs(qs,q)

    rhs = 0.0_rk
    
!!! set density of air at T= 293,15K and speed of sound
    rho0 = params%init%rho
    p0   = params%init%p
    c0   = sqrt(params%material%gamma * p0 / rho0)

!!! mass equation -------------------------------------------------------------------------------
    call Dx1(a_x,qs(:,:,:,2),params%geom%dx1)
    call Dx2(a_y,qs(:,:,:,3),params%geom%dx2)
    call Dx3(a_z,qs(:,:,:,4),params%geom%dx3)
    rhs(:,:,:,1) = -(a_x + a_y + a_z)/rho0


!!! momentum equation (1) ------------------------------------------------------------------------
    call Dx1(p_x,qs(:,:,:,1),params%geom%dx1)

    rhs(:,:,:,2) = -rho0 * c0**2 * p_x

!!! momentum equation (2) ------------------------------------------------------------------------
    call Dx2(p_y,qs(:,:,:,1),params%geom%dx2)

    rhs(:,:,:,3) = -rho0 * c0**2 * p_y

!!! momentum equation (3) ------------------------------------------------------------------------
    call Dx3(p_z,qs(:,:,:,1),params%geom%dx3)

    rhs(:,:,:,4) = -rho0 * c0**2 * p_z

!!! adjoint force
    rhs = rhs + g

!!! adjoint sponge
    call sponge_adjoint(rhs,qs,q)

!!! boundary conditions 2/2: rhss
    call set_boundary_condition_rhss(rhs,qs)

    !call store(rhs,'debug_field_A') 

    deallocate(qs)

    deallocate(a_x)
    deallocate(a_y)
    deallocate(a_z)

    deallocate(p_x)
    deallocate(p_y)
    deallocate(p_z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine right_hand_side_adjoint_3d
!!!=================================================================================================

end module linear_acoustic_p_u
