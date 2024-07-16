module ns_rho_rhou_p_ev_phi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! information/remarks:                                                                           !
  ! - q_in, qs_in in subroutine calls, this is needed due to intent(in) of q/qs in time stepper,   !
  !   thus no modification inside of the rhs is permitted, however no intent works fine, even      !
  !   nobody understands why ... (2019-01-02, ML)                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use boundary_conditions
  use data_reference , only: qref
  use discretisation_x
  use ns_rho_rhou_p_ev_phi_trafos
  use force
  use parallelism
  use parameter
  use sponge_layer
  use volume_penalization

  private

  public right_hand_side_direct_3d
  public right_hand_side_adjoint_3d

  public give_T, give_c, give_r, give_u, give_p

  character(len=max_length_parameter), parameter, public   :: equations_name = "ns_rho_rhou_p_ev_phi"
  character(len=max_length_parameter), parameter, public   :: viscosity_type = "newtonian"

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
    real(kind=rk), dimension(:,:,:)  , allocatable   :: phi,div_1_phi,phi_p
    integer                                          :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:)  , allocatable   :: u_x, u_y, u_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: v_x, v_y, v_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: w_x, w_y, w_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: tau_11, tau_12, tau_13
    real(kind=rk), dimension(:,:,:)  , allocatable   :: tau_22, tau_23
    real(kind=rk), dimension(:,:,:)  , allocatable   :: tau_33
    real(kind=rk), dimension(:,:,:)  , allocatable   :: temp,mu,mu_div,phi_lambda
    real(kind=rk)                                    :: gamma_term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(  q(params%geom%n1b,params%geom%n2b,params%geom%n3b,5))
    allocate(  u(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))

    allocate(a_x(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(a_y(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(a_z(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

    allocate(p_x(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(p_y(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(p_z(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

    allocate(      phi(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(div_1_phi(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(    phi_p(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(phi_rho_u(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))

    ! friction part
    allocate(u_x(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(u_y(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(u_z(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

    allocate(v_x(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(v_y(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(v_z(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

    allocate(w_x(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(w_y(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(w_z(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

    allocate(tau_11(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(tau_12(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(tau_13(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

    allocate(tau_22(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(tau_23(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

    allocate(tau_33(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

    allocate(  temp(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(    mu(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(mu_div(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(phi_lambda(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

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

       phi_p(:,:,i) = phi(:,:,i)*q(:,:,i,5)
    end do

!!! compute derivatives used multiple times
    call Dx1(u_x,u(:,:,:,1),params%geom%dx1)
    call Dx2(u_y,u(:,:,:,1),params%geom%dx2)
    call Dx3(u_z,u(:,:,:,1),params%geom%dx3)

    call Dx1(v_x,u(:,:,:,2),params%geom%dx1)
    call Dx2(v_y,u(:,:,:,2),params%geom%dx2)
    call Dx3(v_z,u(:,:,:,2),params%geom%dx3)

    call Dx1(w_x,u(:,:,:,3),params%geom%dx1)
    call Dx2(w_y,u(:,:,:,3),params%geom%dx2)
    call Dx3(w_z,u(:,:,:,3),params%geom%dx3)

!!! compute friction
    temp   = give_T(q)
    mu     = give_mu(temp)

    mu_div = -2.0_rk/3.0_rk*mu*(u_x + v_y + w_z)

    ! compute only 6 entries due to symmetry
    tau_11 = mu_div + 2*mu*u_x
    tau_22 = mu_div + 2*mu*v_y
    tau_33 = mu_div + 2*mu*w_z

    tau_12 = mu*(u_y + v_x)
    tau_13 = mu*(u_z + w_x)

    tau_23 = mu*(v_z + w_y)

    call Dx1(p_x,q(:,:,:,5),params%geom%dx1) ! dx p
    call Dx2(p_y,q(:,:,:,5),params%geom%dx2) ! dy p
    call Dx3(p_z,q(:,:,:,5),params%geom%dx3) ! dz p

!!! mass equation --------------------------------------------------------------------------------
    call Dx1(a_x,phi_rho_u(:,:,:,1),params%geom%dx1) ! d_x (phi rho u)
    call Dx2(a_y,phi_rho_u(:,:,:,2),params%geom%dx2) ! d_y (phi rho v)
    call Dx3(a_z,phi_rho_u(:,:,:,3),params%geom%dx3) ! d_z (phi rho w)

    rhs(:,:,:,1) = -a_x - a_y - a_z
    rhs(:,:,:,1) = div_1_phi*rhs(:,:,:,1)

!!! momentum equation (1) ------------------------------------------------------------------------
    call Dx1(a_x,phi_rho_u(:,:,:,1)*u(:,:,:,1),params%geom%dx1) ! d_x (phi rho u u)
    call Dx2(a_y,phi_rho_u(:,:,:,2)*u(:,:,:,1),params%geom%dx2) ! d_y (phi rho v u)
    call Dx3(a_z,phi_rho_u(:,:,:,3)*u(:,:,:,1),params%geom%dx3) ! d_z (phi rho w u)
    rhs(:,:,:,2) = -a_x - a_y - a_z

    ! friction
    call Dx1(a_x,phi*tau_11,params%geom%dx1) ! d_x (phi tau_11)
    call Dx2(a_y,phi*tau_12,params%geom%dx2) ! d_y (phi tau_12)
    call Dx3(a_z,phi*tau_13,params%geom%dx3) ! d_z (phi tau_13)
    rhs(:,:,:,2) = rhs(:,:,:,2) + a_x + a_y + a_z

    rhs(:,:,:,2) = div_1_phi*rhs(:,:,:,2) - p_x

!!! momentum equation (2) ------------------------------------------------------------------------
    call Dx1(a_x,phi_rho_u(:,:,:,1)*u(:,:,:,2),params%geom%dx1) ! d_x (phi rho u v)
    call Dx2(a_y,phi_rho_u(:,:,:,2)*u(:,:,:,2),params%geom%dx2) ! d_y (phi rho v v)
    call Dx3(a_z,phi_rho_u(:,:,:,3)*u(:,:,:,2),params%geom%dx3) ! d_z (phi rho w v)
    rhs(:,:,:,3) = -a_x - a_y - a_z

    ! friction
    call Dx1(a_x,phi*tau_12,params%geom%dx1) ! d_x (phi tau_21)
    call Dx2(a_y,phi*tau_22,params%geom%dx2) ! d_y (phi tau_22)
    call Dx3(a_z,phi*tau_23,params%geom%dx3) ! d_z (phi tau_23)
    rhs(:,:,:,3) = rhs(:,:,:,3) + a_x + a_y + a_z

    rhs(:,:,:,3) = div_1_phi*rhs(:,:,:,3) - p_y

!!! momentum equation (3) ------------------------------------------------------------------------
    call Dx1(a_x,phi_rho_u(:,:,:,1)*u(:,:,:,3),params%geom%dx1) ! d_x (phi rho u w)
    call Dx2(a_y,phi_rho_u(:,:,:,2)*u(:,:,:,3),params%geom%dx2) ! d_y (phi rho v w)
    call Dx3(a_z,phi_rho_u(:,:,:,3)*u(:,:,:,3),params%geom%dx3) ! d_z (phi rho w w)
    rhs(:,:,:,4) = -a_x - a_y - a_z

    ! friction
    call Dx1(a_x,phi*tau_13,params%geom%dx1) ! d_x (phi tau_31)
    call Dx2(a_y,phi*tau_23,params%geom%dx2) ! d_y (phi tau_32)
    call Dx3(a_z,phi*tau_33,params%geom%dx3) ! d_z (phi tau_33)
    rhs(:,:,:,4) = rhs(:,:,:,4) + a_x + a_y + a_z

    rhs(:,:,:,4) = div_1_phi*rhs(:,:,:,4) - p_z

!!! energy (pressure) equation -------------------------------------------------------------------
    call Dx1(a_x,phi_p*u(:,:,:,1),params%geom%dx1)
    call Dx2(a_y,phi_p*u(:,:,:,2),params%geom%dx2)
    call Dx3(a_z,phi_p*u(:,:,:,3),params%geom%dx3)

    gamma_term = params%material%gamma / (params%material%gamma - 1.0_rk)  ! gamma /(gamma -1)
    rhs(:,:,:,5) = gamma_term  * (-a_x - a_y - a_z)

    ! heat transfer
    call Dx1(a_x,temp) ! d_x T
    call Dx2(a_y,temp) ! d_y T
    call Dx3(a_z,temp) ! d_z T

    phi_lambda = give_lambda(mu)
    phi_lambda = phi*phi_lambda
    call Dx1(a_x,phi_lambda*a_x) ! d_x (phi lambda d_x T)
    call Dx2(a_y,phi_lambda*a_y) ! d_y (phi lambda d_y T)
    call Dx3(a_z,phi_lambda*a_z) ! d_z (phi lambda d_z T)

    ! heat term
    rhs(:,:,:,5) = rhs(:,:,:,5) + a_x + a_y + a_z ! rhs + d_x (phi lambda d_x T) + d_y (phi lambda d_y T) + d_z (phi lambda d_z T)
    rhs(:,:,:,5) = div_1_phi * rhs(:,:,:,5)

    rhs(:,:,:,5) = rhs(:,:,:,5) + u(:,:,:,1)*p_x(:,:,:) + u(:,:,:,2)*p_y(:,:,:) + u(:,:,:,3)*p_z(:,:,:)

    ! friction term
    rhs(:,:,:,5) = rhs(:,:,:,5) + tau_11*u_x + tau_12*u_y + tau_13*u_z & ! 11 12 13
                                + tau_12*v_x + tau_22*v_y + tau_23*v_z & ! 21 22 23
                                + tau_13*w_x + tau_23*w_y + tau_33*w_z   ! 31 32 33

    gamma_term   = (params%material%gamma - 1.0_rk)          ! (gamma -1)
    rhs(:,:,:,5) = gamma_term * rhs(:,:,:,5)

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

    deallocate(a_x)
    deallocate(a_y)
    deallocate(a_z)

    deallocate(p_x)
    deallocate(p_y)
    deallocate(p_z)

    deallocate(u_x)
    deallocate(u_y)
    deallocate(u_z)

    deallocate(v_x)
    deallocate(v_y)
    deallocate(v_z)

    deallocate(w_x)
    deallocate(w_y)
    deallocate(w_z)

    deallocate(tau_11)
    deallocate(tau_12)
    deallocate(tau_13)

    deallocate(tau_22)
    deallocate(tau_23)

    deallocate(tau_33)

    deallocate(temp)
    deallocate(mu)
    deallocate(mu_div)
    deallocate(phi_lambda)

    deallocate(      phi)
    deallocate(div_1_phi)
    deallocate(    phi_p)
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

end module ns_rho_rhou_p_ev_phi
