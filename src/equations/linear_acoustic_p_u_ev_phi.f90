module linear_acoustic_p_u_ev_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! information/remarks:                                                                             !
! - q_in, qs_in in subroutine calls, this is needed due to intent(in) of q/qs in time stepper,     !
!   thus no modification inside of the rhs is permitted, however no intent works fine, even        !
!   nobody understands why ... (2019-01-02, ML)                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use boundary_conditions
  use data_rhs , only: q0
  use discretisation_x
  use linear_acoustic_p_u_ev_phi_trafos
  use force
  use parameter
  use sponge_layer
  use volume_penalization

  private

  public right_hand_side_direct_3d
  public right_hand_side_adjoint_3d

  public give_T, give_c, give_r, give_u, give_p

  character(len=max_length_parameter), parameter, public :: viscosity_type = "none"
  character(len=max_length_parameter), parameter, public :: equations_name = "linear_acoustic_p_u_ev_phi"

contains

!!!=================================================================================================
  subroutine right_hand_side_direct_3d(rhs,q_in,t)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(out)   :: rhs
    real(kind=rk), dimension(:,:,:,:), intent(in)    :: q_in
    real(kind=rk)                    , intent(in)    :: t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable   :: q
    real(kind=rk), dimension(:,:,:,:), allocatable   :: phi_u
    real(kind=rk), dimension(:,:,:)  , allocatable   :: a_x, a_y, a_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: p_x, p_y, p_z
    real(kind=rk)                                    :: rho0, p0, c0
    real(kind=rk), dimension(:,:,:)  , allocatable   :: phi, div_1_phi
    integer                                          :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(q(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))

    allocate(a_x(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(a_y(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(a_z(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

    allocate(p_x(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(p_y(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(p_z(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

    allocate(      phi(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(div_1_phi(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(    phi_u(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))


!!! set q
    q = q_in

!!! boundary conditions 1/2: q
    call set_boundary_condition_q(q)
     
!!! set density of air at T= 293,15K and speed of sound
    rho0 = params%init%rho
    p0   = params%init%p
    c0  = sqrt(params%material%gamma * p0 / rho0)

!!! unpack phi variables
    phi        = penalization(:,:,:,1)
    div_1_phi  = penalization(:,:,:,2)
    do i = 1, params%geom%n3b
       phi_u(:,:,i,1) = phi(:,:,i) * q(:,:,i,2)
       phi_u(:,:,i,2) = phi(:,:,i) * q(:,:,i,3)
       phi_u(:,:,i,3) = phi(:,:,i) * q(:,:,i,4)
    end do


!!! mass equation -------------------------------------------------------------------------------
    call Dx1(a_x,phi_u(:,:,:,1),params%geom%dx1)
    call Dx2(a_y,phi_u(:,:,:,2),params%geom%dx2)
    call Dx3(a_z,phi_u(:,:,:,3),params%geom%dx3)
    rhs(:,:,:,1) = -rho0 * c0**2 * div_1_phi * (a_x + a_y + a_z)


!!! momentum equation (1) ------------------------------------------------------------------------
    call Dx1(p_x,q(:,:,:,1),params%geom%dx1)

    rhs(:,:,:,2) = -p_x/rho0

!!! momentum equation (2) ------------------------------------------------------------------------
    call Dx2(p_y,q(:,:,:,1),params%geom%dx2)

    rhs(:,:,:,3) = -p_y/rho0

!!! momentum equation (3) ------------------------------------------------------------------------
    call Dx3(p_z,q(:,:,:,1),params%geom%dx3)

    rhs(:,:,:,4) = -p_z/rho0


!!! darcy
    call apply_volume_penalization(rhs,q)

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

    deallocate(phi)
    deallocate(div_1_phi)
    deallocate(phi_u)

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
    real(kind=rk), dimension(:,:,:)  , allocatable   :: u_1,u_2,u_3
    real(kind=rk), dimension(:,:,:)  , allocatable   :: a_x,b_x,c_x,d_x,e_x
    real(kind=rk), dimension(:,:,:)  , allocatable   :: gamma_DOT_p_INV_rho_gm1, INV_rho
    real(kind=rk)                                    :: gm1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    rhs = 0.0_rk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine right_hand_side_adjoint_3d
!!!=================================================================================================

end module linear_acoustic_p_u_ev_phi
