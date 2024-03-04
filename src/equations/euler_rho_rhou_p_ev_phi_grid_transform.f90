module euler_rho_rhou_p_ev_phi_grid_transform

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! information/remarks:                                                                           !
  ! - q_in, qs_in in subroutine calls, this is needed due to intent(in) of q/qs in time stepper,   !
  !   thus no modification inside of the rhs is permitted, however no intent works fine, even      !
  !   nobody understands why ... (2019-01-02, ML)                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use boundary_conditions
  use data_geom
  use data_reference , only: qref
  use discretisation_x
  use euler_rho_rhou_p_ev_phi_grid_transform_trafos
  use force
  use parallelism
  use parameter
  use sponge_layer
  use volume_penalization

  private

  public right_hand_side_direct_3d
  public right_hand_side_adjoint_3d

  public give_T, give_c, give_r, give_u, give_p

  character(len=max_length_parameter), parameter, public   :: equations_name = "euler_rho_rhou_p_ev_phi_grid_transform"
  character(len=max_length_parameter), parameter, public   :: viscosity_type = "none"

contains

!!!=================================================================================================
  subroutine right_hand_side_direct_3d(rhs,q_in,t)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(out)   :: rhs
    real(kind=rk), dimension(:,:,:,:), intent(in)    :: q_in
    real(kind=rk)                    , intent(in)    :: t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable   :: u,q
    real(kind=rk), dimension(:,:,:,:), allocatable   :: phi_rho_u
    real(kind=rk), dimension(:,:,:)  , allocatable   :: a_x, a_y, a_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: p_x, p_y, p_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: phi,div_1_phi,phi_p, phi_rho
    integer                                          :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable   :: u_proj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(     q(params%geom%n1b,params%geom%n2b,params%geom%n3b,5   ))
    allocate(     u(params%geom%n1b,params%geom%n2b,params%geom%n3b,3   ))
    allocate(u_proj(params%geom%n1b,params%geom%n2b,params%geom%n3b,3   ))

    allocate(a_x(params%geom%n1b,params%geom%n2b,params%geom%n3b        ))
    allocate(a_y(params%geom%n1b,params%geom%n2b,params%geom%n3b        ))
    allocate(a_z(params%geom%n1b,params%geom%n2b,params%geom%n3b        ))

    allocate(p_x(params%geom%n1b,params%geom%n2b,params%geom%n3b        ))
    allocate(p_y(params%geom%n1b,params%geom%n2b,params%geom%n3b        ))
    allocate(p_z(params%geom%n1b,params%geom%n2b,params%geom%n3b        ))

    allocate(      phi(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(div_1_phi(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(    phi_p(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(  phi_rho(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(phi_rho_u(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))

!!! set q
    q = q_in

!!! boundary conditions 1/2: q
    call set_boundary_condition_q(q)

!!! unpack variables
    phi       = penalization(:,:,:,1)
    div_1_phi = penalization(:,:,:,2)

    do i = 1,params%geom%n3b ! loop improves performance :)
       u(:,:,i,1) = q(:,:,i,2)/q(:,:,i,1)
       u(:,:,i,2) = q(:,:,i,3)/q(:,:,i,1)
       u(:,:,i,3) = q(:,:,i,4)/q(:,:,i,1)

       phi_rho_u(:,:,i,1) = phi(:,:,i)*q(:,:,i,2)
       phi_rho_u(:,:,i,2) = phi(:,:,i)*q(:,:,i,3)
       phi_rho_u(:,:,i,3) = phi(:,:,i)*q(:,:,i,4)

       phi_rho(:,:,i) = phi(:,:,i)*q(:,:,i,1)
       phi_p(:,:,i)   = phi(:,:,i)*q(:,:,i,5)
    end do

!!! computing projected velocity components
    u_proj(:,:,:,1) = X_metric(:,:,:,1,1) * u(:,:,:,1) + X_metric(:,:,:,1,2) * u(:,:,:,2) + X_metric(:,:,:,1,3) * u(:,:,:,3)
    u_proj(:,:,:,2) = X_metric(:,:,:,2,1) * u(:,:,:,1) + X_metric(:,:,:,2,2) * u(:,:,:,2) + X_metric(:,:,:,2,3) * u(:,:,:,3)
    u_proj(:,:,:,3) = X_metric(:,:,:,3,1) * u(:,:,:,1) + X_metric(:,:,:,3,2) * u(:,:,:,2) + X_metric(:,:,:,3,3) * u(:,:,:,3)

!!! computing derivatives of p with respect to physical coordinates
    call Dx1(a_x,q(:,:,:,5),params%geom%dx1)
    call Dx2(a_y,q(:,:,:,5),params%geom%dx2)
    call Dx3(a_z,q(:,:,:,5),params%geom%dx3)
    p_x = X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * a_y + X_metric(:,:,:,3,1) * a_z
    p_y = X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * a_y + X_metric(:,:,:,3,2) * a_z
    p_z = X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * a_y + X_metric(:,:,:,3,3) * a_z

!!! mass equation projected velcocity -----------------------------------------------------------------
    call Dx1(a_x,phi_rho * u_proj(:,:,:,1),params%geom%dx1)
    call Dx2(a_y,phi_rho * u_proj(:,:,:,2),params%geom%dx2)
    call Dx3(a_z,phi_rho * u_proj(:,:,:,3),params%geom%dx3)

    rhs(:,:,:,1) = -a_x - a_y - a_z
    rhs(:,:,:,1) = X_jacobian * div_1_phi * rhs(:,:,:,1)

!!! momentum equation (1) ------------------------------------------------------------------------
    call Dx1(a_x,phi_rho_u(:,:,:,1)*u_proj(:,:,:,1),params%geom%dx1)
    call Dx2(a_y,phi_rho_u(:,:,:,1)*u_proj(:,:,:,2),params%geom%dx2)
    call Dx3(a_z,phi_rho_u(:,:,:,1)*u_proj(:,:,:,3),params%geom%dx3)

    rhs(:,:,:,2) = div_1_phi * (-a_x - a_y - a_z) - p_x
    rhs(:,:,:,2) = X_jacobian * rhs(:,:,:,2)

!!! momentum equation (2) ------------------------------------------------------------------------
    call Dx1(a_x,phi_rho_u(:,:,:,2)*u_proj(:,:,:,1),params%geom%dx1)
    call Dx2(a_y,phi_rho_u(:,:,:,2)*u_proj(:,:,:,2),params%geom%dx2)
    call Dx3(a_z,phi_rho_u(:,:,:,2)*u_proj(:,:,:,3),params%geom%dx3)

    rhs(:,:,:,3) = div_1_phi * (-a_x - a_y - a_z) - p_y
    rhs(:,:,:,3) = X_jacobian * rhs(:,:,:,3)

!!! momentum equation (3) ------------------------------------------------------------------------
    call Dx1(a_x,phi_rho_u(:,:,:,3)*u_proj(:,:,:,1),params%geom%dx1)
    call Dx2(a_y,phi_rho_u(:,:,:,3)*u_proj(:,:,:,2),params%geom%dx2)
    call Dx3(a_z,phi_rho_u(:,:,:,3)*u_proj(:,:,:,3),params%geom%dx3)

    rhs(:,:,:,4) = div_1_phi * (-a_x - a_y - a_z) - p_z
    rhs(:,:,:,4) = X_jacobian * rhs(:,:,:,4)

 !!! energy (pressure) equation -------------------------------------------------------------------
    call Dx1(a_x,phi_p * u_proj(:,:,:,1),params%geom%dx1)
    call Dx2(a_y,phi_p * u_proj(:,:,:,2),params%geom%dx2)
    call Dx3(a_z,phi_p * u_proj(:,:,:,3),params%geom%dx3)

    rhs(:,:,:,5) = -params%material%gamma * div_1_phi * (a_x + a_y + a_z)
    rhs(:,:,:,5) = rhs(:,:,:,5) + (params%material%gamma - 1.0_rk)*(u(:,:,:,1)*p_x(:,:,:) + u(:,:,:,2)*p_y(:,:,:) + u(:,:,:,3)*p_z(:,:,:))
    rhs(:,:,:,5) = rhs(:,:,:,5) * X_jacobian

!!! darcy
    call apply_volume_penalization(rhs,q)

!!! forcing (field)
    call apply_force(rhs)

!!! sponge layer
    call sponge(rhs,q,qref)

!!! boundary conditions 2/2: rhs
    call set_boundary_condition_rhs(rhs,q)

    deallocate(q)
    deallocate(u)
    deallocate(u_proj)

    deallocate(a_x)
    deallocate(a_y)
    deallocate(a_z)

    deallocate(p_x)
    deallocate(p_y)
    deallocate(p_z)

    deallocate(      phi)
    deallocate(div_1_phi)
    deallocate(    phi_p)
    deallocate(  phi_rho)
    deallocate(phi_rho_u)

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

    rhs = 0.0_rk
    write(*,*) "error in euler_rho_rhou_p_ev_phi.f90:right_hand_side_adjoint_3d:not yet implemented --> STOP"
    call stop_cruna

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine right_hand_side_adjoint_3d
!!!=================================================================================================

end module euler_rho_rhou_p_ev_phi_grid_transform
