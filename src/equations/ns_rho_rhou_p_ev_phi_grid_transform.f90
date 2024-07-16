module ns_rho_rhou_p_ev_phi_grid_transform

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
  use ns_rho_rhou_p_ev_phi_grid_transform_trafos
  use force
  use parallelism
  use parameter
  use sponge_layer
  use volume_penalization

  private

  public right_hand_side_direct_3d
  public right_hand_side_adjoint_3d

  public give_T, give_c, give_r, give_u, give_p

  character(len=max_length_parameter), parameter, public   :: equations_name = "ns_rho_rhou_p_ev_phi_grid_transform"
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
    real(kind=rk), dimension(:,:,:)  , allocatable   :: phi,div_1_phi,phi_p, phi_rho
    integer                                          :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:)  , allocatable   :: u_x, u_y, u_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: v_x, v_y, v_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: w_x, w_y, w_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: tau_11, tau_12, tau_13
    real(kind=rk), dimension(:,:,:)  , allocatable   :: tau_22, tau_23
    real(kind=rk), dimension(:,:,:)  , allocatable   :: tau_33
    real(kind=rk), dimension(:,:,:)  , allocatable   :: temp,mu,mu_div,phi_lambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable   :: u_proj
    real(kind=rk), dimension(:,:,:)  , allocatable   :: temp_x, temp_y, temp_z
    real(kind=rk)                                    :: gamma_term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(     q(params%geom%n1b,params%geom%n2b,params%geom%n3b,5   ))
    allocate(     u(params%geom%n1b,params%geom%n2b,params%geom%n3b,3   ))

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

    allocate(  temp(params%geom%n1b,params%geom%n2b,params%geom%n3b     ))
    allocate(    mu(params%geom%n1b,params%geom%n2b,params%geom%n3b     ))
    allocate(mu_div(params%geom%n1b,params%geom%n2b,params%geom%n3b     ))
    allocate(phi_lambda(params%geom%n1b,params%geom%n2b,params%geom%n3b ))

    ! grid tranformation
    allocate(u_proj(params%geom%n1b,params%geom%n2b,params%geom%n3b,3   ))

    allocate(temp_x(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(temp_y(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))
    allocate(temp_z(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

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

!!! compute derivatives used multiple times (in physical coordinates)
    ! without Jacobian (u_tilde)
    call Dx1(a_x,u(:,:,:,1),params%geom%dx1)
    call Dx2(a_y,u(:,:,:,1),params%geom%dx2)
    call Dx3(a_z,u(:,:,:,1),params%geom%dx3)
    u_x = X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * a_y + X_metric(:,:,:,3,1) * a_z
    u_y = X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * a_y + X_metric(:,:,:,3,2) * a_z
    u_z = X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * a_y + X_metric(:,:,:,3,3) * a_z

    call Dx1(a_x,u(:,:,:,2),params%geom%dx1)
    call Dx2(a_y,u(:,:,:,2),params%geom%dx2)
    call Dx3(a_z,u(:,:,:,2),params%geom%dx3)
    v_x = X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * a_y + X_metric(:,:,:,3,1) * a_z
    v_y = X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * a_y + X_metric(:,:,:,3,2) * a_z
    v_z = X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * a_y + X_metric(:,:,:,3,3) * a_z

    call Dx1(a_x,u(:,:,:,3),params%geom%dx1)
    call Dx2(a_y,u(:,:,:,3),params%geom%dx2)
    call Dx3(a_z,u(:,:,:,3),params%geom%dx3)
    w_x = X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * a_y + X_metric(:,:,:,3,1) * a_z
    w_y = X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * a_y + X_metric(:,:,:,3,2) * a_z
    w_z = X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * a_y + X_metric(:,:,:,3,3) * a_z

!!! computing projected quantities
    ! projected velocity components
    u_proj(:,:,:,1) = X_metric(:,:,:,1,1) * u(:,:,:,1) + X_metric(:,:,:,1,2) * u(:,:,:,2) + X_metric(:,:,:,1,3) * u(:,:,:,3)
    u_proj(:,:,:,2) = X_metric(:,:,:,2,1) * u(:,:,:,1) + X_metric(:,:,:,2,2) * u(:,:,:,2) + X_metric(:,:,:,2,3) * u(:,:,:,3)
    u_proj(:,:,:,3) = X_metric(:,:,:,3,1) * u(:,:,:,1) + X_metric(:,:,:,3,2) * u(:,:,:,2) + X_metric(:,:,:,3,3) * u(:,:,:,3)

    ! computing derivatives of p with respect to physical coordinates (without Jacobian)
    call Dx1(a_x,q(:,:,:,5),params%geom%dx1)
    call Dx2(a_y,q(:,:,:,5),params%geom%dx2)
    call Dx3(a_z,q(:,:,:,5),params%geom%dx3)
    p_x = X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * a_y + X_metric(:,:,:,3,1) * a_z
    p_y = X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * a_y + X_metric(:,:,:,3,2) * a_z
    p_z = X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * a_y + X_metric(:,:,:,3,3) * a_z

!!! compute friction
    temp   = give_T(q)
    mu     = give_mu(temp)

    mu_div = -2.0_rk/3.0_rk*mu*(u_x + v_y + w_z)

    ! compute only 6 entries due to symmetry (first without Jacobian)
    tau_11 = mu_div + 2*mu*u_x
    tau_22 = mu_div + 2*mu*v_y
    tau_33 = mu_div + 2*mu*w_z

    tau_12 = mu*(u_y + v_x)
    tau_13 = mu*(u_z + w_x)

    tau_23 = mu*(v_z + w_y)

    ! compute heat
    phi_lambda = give_lambda(mu) ! get lambda (only)

    ! computing derivatives of T with respect to physical coordinates (with Jacobian and lambda)
    call Dx1(a_x,temp,params%geom%dx1) ! d_xi   (T)
    call Dx2(a_y,temp,params%geom%dx2) ! d_eta  (T)
    call Dx3(a_z,temp,params%geom%dx3) ! d_zeta (T)
    temp_x = X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * a_y + X_metric(:,:,:,3,1) * a_z ! d_x (T) (without Jacobian)
    temp_y = X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * a_y + X_metric(:,:,:,3,2) * a_z ! d_y (T) (without Jacobian)
    temp_z = X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * a_y + X_metric(:,:,:,3,3) * a_z ! d_z (T) (without Jacobian)

    ! phi*lambda*grad(T) in physical coordinates (with Jacobian)
    phi_lambda = phi * phi_lambda * X_jacobian ! phi * lambda * J
    temp_x = phi_lambda * temp_x               ! phi * lambda * d_x (T)
    temp_y = phi_lambda * temp_y               ! phi * lambda * d_y (T)
    temp_z = phi_lambda * temp_z               ! phi * lambda * d_z (T)

!!! mass equation projected velcocity -----------------------------------------------------------------
    call Dx1(a_x,phi_rho * u_proj(:,:,:,1),params%geom%dx1) ! d_xi   (phi rho u_hat)
    call Dx2(a_y,phi_rho * u_proj(:,:,:,2),params%geom%dx2) ! d_eta  (phi rho v_hat)
    call Dx3(a_z,phi_rho * u_proj(:,:,:,3),params%geom%dx3) ! d_zeta (phi rho w_hat)

    rhs(:,:,:,1) = -a_x - a_y - a_z
    rhs(:,:,:,1) = X_jacobian * div_1_phi * rhs(:,:,:,1)

!!! momentum equation (1) ------------------------------------------------------------------------
    call Dx1(a_x,phi_rho_u(:,:,:,1)*u_proj(:,:,:,1),params%geom%dx1) ! d_xi   (phi rho u u_hat)
    call Dx2(a_y,phi_rho_u(:,:,:,1)*u_proj(:,:,:,2),params%geom%dx2) ! d_eta  (phi rho u v_hat)
    call Dx3(a_z,phi_rho_u(:,:,:,1)*u_proj(:,:,:,3),params%geom%dx3) ! d_zeta (phi rho u w_hat)
    rhs(:,:,:,2) = -a_x - a_y - a_z                                  !  -d_xi (phi rho u u_hat) -d_eta (phi rho u v_hat) -d_zeta (phi rho u w_hat)

    ! friction
    a_x = X_metric(:,:,:,1,1) * tau_11 + X_metric(:,:,:,1,2) * tau_12 + X_metric(:,:,:,1,3) * tau_13 ! tau_proj_11
    a_y = X_metric(:,:,:,2,1) * tau_11 + X_metric(:,:,:,2,2) * tau_12 + X_metric(:,:,:,2,3) * tau_13 ! tau_proj_12
    a_z = X_metric(:,:,:,3,1) * tau_11 + X_metric(:,:,:,3,2) * tau_12 + X_metric(:,:,:,3,3) * tau_13 ! tau_proj_13
    call Dx1(a_x,phi*a_x,params%geom%dx1) ! d_xi   (phi tau_proj_11)
    call Dx2(a_y,phi*a_y,params%geom%dx2) ! d_eta  (phi tau_proj_12)
    call Dx3(a_z,phi*a_z,params%geom%dx3) ! d_zeta (phi tau_proj_13)
    rhs(:,:,:,2) = rhs(:,:,:,2) + a_x + a_y + a_z   ! +d_xi(phi tau_proj_11) +d_eta  (phi tau_proj_12) +d_zeta (phi tau_proj_13)

    rhs(:,:,:,2) = div_1_phi * rhs(:,:,:,2) - p_x
    rhs(:,:,:,2) = X_jacobian * rhs(:,:,:,2)

!!! momentum equation (2) ------------------------------------------------------------------------
    call Dx1(a_x,phi_rho_u(:,:,:,2)*u_proj(:,:,:,1),params%geom%dx1) ! d_xi   (phi rho v u_hat)
    call Dx2(a_y,phi_rho_u(:,:,:,2)*u_proj(:,:,:,2),params%geom%dx2) ! d_eta  (phi rho v v_hat)
    call Dx3(a_z,phi_rho_u(:,:,:,2)*u_proj(:,:,:,3),params%geom%dx3) ! d_zeta (phi rho v w_hat)
    rhs(:,:,:,3) = -a_x - a_y - a_z                                  ! -d_xi (phi rho v u_hat) -d_eta (phi rho v v_hat) -d_zeta (phi rho v w_hat)

    ! friction
    a_x = X_metric(:,:,:,1,1) * tau_12 + X_metric(:,:,:,1,2) * tau_22 + X_metric(:,:,:,1,3) * tau_23 ! tau_proj_21
    a_y = X_metric(:,:,:,2,1) * tau_12 + X_metric(:,:,:,2,2) * tau_22 + X_metric(:,:,:,2,3) * tau_23 ! tau_proj_22
    a_z = X_metric(:,:,:,3,1) * tau_12 + X_metric(:,:,:,3,2) * tau_22 + X_metric(:,:,:,3,3) * tau_23 ! tau_proj_23
    call Dx1(a_x,phi*a_x,params%geom%dx1) ! d_xi   (phi tau_proj_21)
    call Dx2(a_y,phi*a_y,params%geom%dx2) ! d_eta  (phi tau_proj_22)
    call Dx3(a_z,phi*a_z,params%geom%dx3) ! d_zeta (phi tau_proj_23)

    rhs(:,:,:,3) = rhs(:,:,:,3) + a_x + a_y + a_z   ! +d_xi(phi tau_proj_21) +d_eta(phi tau_proj_22) +d_zeta(phi tau_proj_23)
    rhs(:,:,:,3) = div_1_phi * rhs(:,:,:,3) - p_y
    rhs(:,:,:,3) = X_jacobian * rhs(:,:,:,3)

!!! momentum equation (3) ------------------------------------------------------------------------
    call Dx1(a_x,phi_rho_u(:,:,:,3)*u_proj(:,:,:,1),params%geom%dx1) ! d_xi   (phi rho w u_hat)
    call Dx2(a_y,phi_rho_u(:,:,:,3)*u_proj(:,:,:,2),params%geom%dx2) ! d_eta  (phi rho w v_hat)
    call Dx3(a_z,phi_rho_u(:,:,:,3)*u_proj(:,:,:,3),params%geom%dx3) ! d_zeta (phi rho w w_hat)
    rhs(:,:,:,4) = -a_x - a_y - a_z                                  ! -d_xi (phi rho w u_hat) -d_eta (phi rho w v_hat) -d_zeta (phi rho w w_hat)

    ! friction
    a_x = X_metric(:,:,:,1,1) * tau_13 + X_metric(:,:,:,1,2) * tau_23 + X_metric(:,:,:,1,3) * tau_33 ! tau_proj_31
    a_y = X_metric(:,:,:,2,1) * tau_13 + X_metric(:,:,:,2,2) * tau_23 + X_metric(:,:,:,2,3) * tau_33 ! tau_proj_32
    a_z = X_metric(:,:,:,3,1) * tau_13 + X_metric(:,:,:,3,2) * tau_23 + X_metric(:,:,:,3,3) * tau_33 ! tau_proj_33
    call Dx1(a_x,phi*a_x,params%geom%dx1) ! d_xi   (phi tau_proj_31)
    call Dx2(a_y,phi*a_y,params%geom%dx2) ! d_eta  (phi tau_proj_32)
    call Dx3(a_z,phi*a_z,params%geom%dx3) ! d_zeta (phi tau_proj_33)

    rhs(:,:,:,4) = rhs(:,:,:,4) + a_x + a_y + a_z ! +d_xi(phi tau_proj_31) + d_eta(phi tau_proj_32) + d_zeta(phi tau_proj_33)
    rhs(:,:,:,4) = div_1_phi * rhs(:,:,:,4) - p_z
    rhs(:,:,:,4) = X_jacobian * rhs(:,:,:,4)

!!! energy (pressure) equation -------------------------------------------------------------------
    call Dx1(a_x,phi_p * u_proj(:,:,:,1),params%geom%dx1) ! d_xi   (phi p u_hat)
    call Dx2(a_y,phi_p * u_proj(:,:,:,2),params%geom%dx2) ! d_eta  (phi p v_hat)
    call Dx3(a_z,phi_p * u_proj(:,:,:,3),params%geom%dx3) ! d_zeta (phi p w_hat)

    gamma_term = params%material%gamma / (params%material%gamma - 1.0_rk)  ! gamma /(gamma -1)
    rhs(:,:,:,5) = gamma_term * (-a_x - a_y - a_z)

    ! heat flux in x
    call Dx1(a_x,temp_x) ! d_xi   phi * lambda * T_x (phi,lambda in temp_x)
    call Dx2(a_y,temp_x) ! d_eta  phi * lambda * T_x (phi,lambda in temp_x)
    call Dx3(a_z,temp_x) ! d_zeta phi * lambda * T_x (phi,lambda in temp_x)
    rhs(:,:,:,5) = rhs(:,:,:,5) + X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * a_y + X_metric(:,:,:,3,1) * a_z

    ! heat flux in y
    call Dx1(a_x,temp_y) ! d_xi   phi * lambda * T_y (phi,lambda in temp_y)
    call Dx2(a_y,temp_y) ! d_eta  phi * lambda * T_y (phi,lambda in temp_y)
    call Dx3(a_z,temp_y) ! d_zeta phi * lambda * T_y (phi,lambda in temp_y)
    rhs(:,:,:,5) = rhs(:,:,:,5) + X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * a_y + X_metric(:,:,:,3,2) * a_z

    ! heat flux in z
    call Dx1(a_x,temp_z) ! d_xi   phi * lambda * T_z (phi,lambda in temp_z)
    call Dx2(a_y,temp_z) ! d_eta  phi * lambda * T_z (phi,lambda in temp_z)
    call Dx3(a_z,temp_z) ! d_zeta phi * lambda * T_z (phi,lambda in temp_z)
    rhs(:,:,:,5) = rhs(:,:,:,5) + X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * a_y + X_metric(:,:,:,3,3) * a_z

    rhs(:,:,:,5) = div_1_phi * rhs(:,:,:,5)

    rhs(:,:,:,5) = rhs(:,:,:,5) + u(:,:,:,1)*p_x(:,:,:) + u(:,:,:,2)*p_y(:,:,:) + u(:,:,:,3)*p_z(:,:,:)

    ! friction term: u_tilde (accounting for Jacobian of stresses later of velocities)
    rhs(:,:,:,5) = rhs(:,:,:,5) + X_jacobian * (tau_11*u_x + tau_12*u_y + tau_13*u_z) ! 11 12 13
    rhs(:,:,:,5) = rhs(:,:,:,5) + X_jacobian * (tau_12*v_x + tau_22*v_y + tau_23*v_z) ! 21 22 23
    rhs(:,:,:,5) = rhs(:,:,:,5) + X_jacobian * (tau_13*w_x + tau_23*w_y + tau_33*w_z) ! 31 32 33

    gamma_term   =  (params%material%gamma - 1.0_rk)          ! (gamma -1)
    rhs(:,:,:,5) = X_jacobian * gamma_term * rhs(:,:,:,5)     ! accounting for Jacobian of velocities

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

    deallocate(u_proj)

    deallocate(temp_x)
    deallocate(temp_y)
    deallocate(temp_z)

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

end module ns_rho_rhou_p_ev_phi_grid_transform
