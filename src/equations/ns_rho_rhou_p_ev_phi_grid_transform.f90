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
!   use sponge_layer
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
    
!!! move volume penalization in time
    call move_volume_penalization(t) 

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
    tau_11 = X_jacobian * (mu_div + 2*mu*u_x)
    tau_22 = X_jacobian * (mu_div + 2*mu*v_y)
    tau_33 = X_jacobian * (mu_div + 2*mu*w_z)

    tau_12 = X_jacobian * (mu*(u_y + v_x))
    tau_13 = X_jacobian * (mu*(u_z + w_x))

    tau_23 = X_jacobian * (mu*(v_z + w_y))


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
    rhs(:,:,:,5) = rhs(:,:,:,5) +  tau_11*u_x + tau_12*u_y + tau_13*u_z ! 11 12 13
    rhs(:,:,:,5) = rhs(:,:,:,5) +  tau_12*v_x + tau_22*v_y + tau_23*v_z ! 21 22 23
    rhs(:,:,:,5) = rhs(:,:,:,5) +  tau_13*w_x + tau_23*w_y + tau_33*w_z ! 31 32 33

    gamma_term   =  (params%material%gamma - 1.0_rk)          ! (gamma -1)
    rhs(:,:,:,5) = X_jacobian * gamma_term * rhs(:,:,:,5)     ! accounting for Jacobian of velocities

!!! darcy
    call apply_volume_penalization(rhs,q)

!!! forcing (field)
    call apply_force(rhs)

!!! sponge layer
!     call sponge(rhs,q,qref)

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
    real(kind=rk), dimension(:,:,:,:), allocatable   :: qs
    real(kind=rk), dimension(:,:,:,:), allocatable   :: qs_x1,qs_x2,qs_x3
    real(kind=rk), dimension(:,:,:)  , allocatable   :: u_1,u_2,u_3
    real(kind=rk), dimension(:,:,:)  , allocatable   :: u_x, u_y, u_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: v_x, v_y, v_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: w_x, w_y, w_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: a_x,b_x,c_x,d_x,e_x
    real(kind=rk), dimension(:,:,:)  , allocatable   :: qs1_x,qs2_x,qs3_x,qs4_x,qs5_x
    real(kind=rk), dimension(:,:,:)  , allocatable   :: qs1_y,qs2_y,qs3_y,qs4_y,qs5_y
    real(kind=rk), dimension(:,:,:)  , allocatable   :: qs1_z,qs2_z,qs3_z,qs4_z,qs5_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: mu,dmu_dT,lambda, temp
    real(kind=rk), dimension(:,:,:)  , allocatable   :: gamma_DOT_p_INV_rho_gm1,INV_rho,p_INV_rho,p_INV_rho_squared
    real(kind=rk), dimension(:,:,:)  , allocatable   :: phi
    real(kind=rk)                                    :: gm1,gamma_INV_gm1
    real(kind=rk)                                    :: W,R,W_INV_R,Cp,Pr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(qs(params%geom%n1b,params%geom%n2b,params%geom%n3b,5))
    allocate(qs_x1(params%geom%n1b,params%geom%n2b,params%geom%n3b,5))
    allocate(qs_x2(params%geom%n1b,params%geom%n2b,params%geom%n3b,5))
    allocate(qs_x3(params%geom%n1b,params%geom%n2b,params%geom%n3b,5))

    allocate(u_1(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(u_2(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(u_3(params%geom%n1b,params%geom%n2b,params%geom%n3b))

    allocate(gamma_DOT_p_INV_rho_gm1(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(      p_INV_rho_squared(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(              p_INV_rho(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(                INV_rho(params%geom%n1b,params%geom%n2b,params%geom%n3b))

    allocate(a_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(b_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(c_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(d_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(e_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))

    allocate(qs1_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(qs2_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(qs3_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(qs4_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(qs5_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))

    allocate(qs1_y(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(qs2_y(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(qs3_y(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(qs4_y(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(qs5_y(params%geom%n1b,params%geom%n2b,params%geom%n3b))

    allocate(qs1_z(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(qs2_z(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(qs3_z(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(qs4_z(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(qs5_z(params%geom%n1b,params%geom%n2b,params%geom%n3b))


    allocate(u_x(params%geom%n1b,params%geom%n2b,params%geom%n3b   ))
    allocate(u_y(params%geom%n1b,params%geom%n2b,params%geom%n3b   ))
    allocate(u_z(params%geom%n1b,params%geom%n2b,params%geom%n3b   ))

    allocate(v_x(params%geom%n1b,params%geom%n2b,params%geom%n3b   ))
    allocate(v_y(params%geom%n1b,params%geom%n2b,params%geom%n3b   ))
    allocate(v_z(params%geom%n1b,params%geom%n2b,params%geom%n3b   ))

    allocate(w_x(params%geom%n1b,params%geom%n2b,params%geom%n3b   ))
    allocate(w_y(params%geom%n1b,params%geom%n2b,params%geom%n3b   ))
    allocate(w_z(params%geom%n1b,params%geom%n2b,params%geom%n3b   ))

    allocate(  temp(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(    mu(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(dmu_dT(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(lambda(params%geom%n1b,params%geom%n2b,params%geom%n3b))

    allocate(        phi(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

!!! set qs
    qs = qs_in

!!! boundary conditions 1/2: qs
    call set_boundary_condition_qs(qs,q)

    gm1                     = params%material%gamma - 1.0_rk
    gamma_DOT_p_INV_rho_gm1 = params%material%gamma*q(:,:,:,5)/(q(:,:,:,1) * gm1)
    gamma_INV_gm1           = params%material%gamma/gm1
    INV_rho                 = 1.0_rk/q(:,:,:,1)
    p_INV_rho               = q(:,:,:,5)/q(:,:,:,1)
    p_INV_rho_squared       = q(:,:,:,5)/q(:,:,:,1)**2

    !!! compute friction
    temp    = give_T(q)
    mu      = give_mu(temp) ! give_mu(give_T(temp)) ?
    lambda  = give_lambda(mu)
    W       = give_W()
    R       = give_R_const()
    Cp      = give_cp()
    Pr      = give_Pr()
    W_INV_R = W/R
    dmu_dT  = give_dmu_dT(temp,mu)

    rhs = 0.0_rk

    ! unpack direct state
    u_1 = q(:,:,:,2)/q(:,:,:,1)
    u_2 = q(:,:,:,3)/q(:,:,:,1)
    u_3 = q(:,:,:,4)/q(:,:,:,1)

    ! unpack variables
    phi       = penalization(:,:,:,1)

!!! D_xi qs
    call Dx1(qs1_x,qs(:,:,:,1),params%geom%dx1)
    call Dx1(qs2_x,qs(:,:,:,2),params%geom%dx1)
    call Dx1(qs3_x,qs(:,:,:,3),params%geom%dx1)
    call Dx1(qs4_x,qs(:,:,:,4),params%geom%dx1)
    call Dx1(qs5_x,qs(:,:,:,5),params%geom%dx1)

!!! D_eta qs
    call Dx2(qs1_y,qs(:,:,:,1),params%geom%dx2)
    call Dx2(qs2_y,qs(:,:,:,2),params%geom%dx2)
    call Dx2(qs3_y,qs(:,:,:,3),params%geom%dx2)
    call Dx2(qs4_y,qs(:,:,:,4),params%geom%dx2)
    call Dx2(qs5_y,qs(:,:,:,5),params%geom%dx2)

!!! D_zeta qs
    call Dx3(qs1_z,qs(:,:,:,1),params%geom%dx3)
    call Dx3(qs2_z,qs(:,:,:,2),params%geom%dx3)
    call Dx3(qs3_z,qs(:,:,:,3),params%geom%dx3)
    call Dx3(qs4_z,qs(:,:,:,4),params%geom%dx3)
    call Dx3(qs5_z,qs(:,:,:,5),params%geom%dx3)



!!! compute derivatives used multiple times
    ! d_x1(u1), d_x2(u1), d_x3(u1)
    call Dx1(a_x,u_1,params%geom%dx1)
    call Dx2(b_x,u_1,params%geom%dx2)
    call Dx3(c_x,u_1,params%geom%dx3)

    u_x =  X_jacobian * (X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x)
    u_y =  X_jacobian * (X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x)
    u_z =  X_jacobian * (X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x)

    ! d_x1(u2), d_x2(u2), d_x3(u2)
    call Dx1(a_x,u_2,params%geom%dx1)
    call Dx2(b_x,u_2,params%geom%dx2)
    call Dx3(c_x,u_2,params%geom%dx3)

    v_x =  X_jacobian * (X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x)
    v_y =  X_jacobian * (X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x)
    v_z =  X_jacobian * (X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x)

    ! d_x1(u3), d_x2(u3), d_x3(u3)
    call Dx1(a_x,u_3,params%geom%dx1)
    call Dx2(b_x,u_3,params%geom%dx2)
    call Dx3(c_x,u_3,params%geom%dx3)

    w_x =  X_jacobian * (X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x)
    W_y =  X_jacobian * (X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x)
    w_z =  X_jacobian * (X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x)

    !  qs_x1 = dx_1(q*) (without Jacobian)
    qs_x1(:,:,:,1) = (X_metric(:,:,:,1,1) * qs1_x + X_metric(:,:,:,2,1) * qs1_y + X_metric(:,:,:,3,1) * qs1_z)
    qs_x1(:,:,:,2) = (X_metric(:,:,:,1,1) * qs2_x + X_metric(:,:,:,2,1) * qs2_y + X_metric(:,:,:,3,1) * qs2_z)
    qs_x1(:,:,:,3) = (X_metric(:,:,:,1,1) * qs3_x + X_metric(:,:,:,2,1) * qs3_y + X_metric(:,:,:,3,1) * qs3_z)
    qs_x1(:,:,:,4) = (X_metric(:,:,:,1,1) * qs4_x + X_metric(:,:,:,2,1) * qs4_y + X_metric(:,:,:,3,1) * qs4_z)
    qs_x1(:,:,:,5) = (X_metric(:,:,:,1,1) * qs5_x + X_metric(:,:,:,2,1) * qs5_y + X_metric(:,:,:,3,1) * qs5_z)

    !  qs_x2 = dx_2(q*) (without Jacobian)
    qs_x2(:,:,:,1) = (X_metric(:,:,:,1,2) * qs1_x + X_metric(:,:,:,2,2) * qs1_y + X_metric(:,:,:,3,2) * qs1_z)
    qs_x2(:,:,:,2) = (X_metric(:,:,:,1,2) * qs2_x + X_metric(:,:,:,2,2) * qs2_y + X_metric(:,:,:,3,2) * qs2_z)
    qs_x2(:,:,:,3) = (X_metric(:,:,:,1,2) * qs3_x + X_metric(:,:,:,2,2) * qs3_y + X_metric(:,:,:,3,2) * qs3_z)
    qs_x2(:,:,:,4) = (X_metric(:,:,:,1,2) * qs4_x + X_metric(:,:,:,2,2) * qs4_y + X_metric(:,:,:,3,2) * qs4_z)
    qs_x2(:,:,:,5) = (X_metric(:,:,:,1,2) * qs5_x + X_metric(:,:,:,2,2) * qs5_y + X_metric(:,:,:,3,2) * qs5_z)

    !  qs_x3 = dx_3(q*) (without Jacobian)
    qs_x3(:,:,:,1) = (X_metric(:,:,:,1,3) * qs1_x + X_metric(:,:,:,2,3) * qs1_y + X_metric(:,:,:,3,3) * qs1_z)
    qs_x3(:,:,:,2) = (X_metric(:,:,:,1,3) * qs2_x + X_metric(:,:,:,2,3) * qs2_y + X_metric(:,:,:,3,3) * qs2_z)
    qs_x3(:,:,:,3) = (X_metric(:,:,:,1,3) * qs3_x + X_metric(:,:,:,2,3) * qs3_y + X_metric(:,:,:,3,3) * qs3_z)
    qs_x3(:,:,:,4) = (X_metric(:,:,:,1,3) * qs4_x + X_metric(:,:,:,2,3) * qs4_y + X_metric(:,:,:,3,3) * qs4_z)
    qs_x3(:,:,:,5) = (X_metric(:,:,:,1,3) * qs5_x + X_metric(:,:,:,2,3) * qs5_y + X_metric(:,:,:,3,3) * qs5_z)

! It is possible to factor out the Jacobian from alomst every term (except H)

    ! Part  1: ===============================
    ! B1)             -B1^T * dx_1(q*)
    ! B2)             -B2^T * dx_2(q*)
    ! B3)             -B3^T * dx_3(q*)
    !
    ! Part  2 ================================
    ! C1)         -D_x1(C1^T * q*)
    ! C1)         -D_x2(C2^T * q*)
    ! C1)         -D_x3(C3^T * q*)
    !
    ! Part  3 ================================
    ! C1*)              +C1* * dx_1(q)
    ! C2*)              +C2* * dx_2(q)
    ! C3*)              +C3* * dx_3(q)
    !
    ! Part  4 ================================
    ! D11)       -D_x1(D11^T * dx_1(q*))
    ! D12)       -D_x2(D12^T * dx_1(q*))
    ! D13)       -D_x3(D13^T * dx_1(q*))
    ! D21)       -D_x1(D21^T * dx_2(q*))
    ! D22)       -D_x2(D22^T * dx_2(q*))
    ! D23)       -D_x3(D23^T * dx_2(q*))
    ! D31)       -D_x1(D31^T * dx_3(q*))
    ! D32)       -D_x2(D32^T * dx_3(q*))
    ! D33)       -D_x3(D33^T * dx_3(q*))
    !
    ! Part  5 ================================
    ! E1)             +E1^T * dx_1(q*)
    ! E2)             +E2^T * dx_2(q*)
    ! E3)             +E3^T * dx_3(q*)
    !
    ! Part  6 ================================
    ! F1) -FII_T*D_x1(F_I^T * dx_1(q*))
    ! F2) -FII_T*D_x2(F_I^T * dx_2(q*))
    ! F3) -FII_T*D_x3(F_I^T * dx_3(q*))
    !
    ! Part  7 ================================
    ! G1)        +D_x1(G1^T * q*)
    ! G1)        +D_x2(G2^T * q*)
    ! G1)        +D_x3(G3^T * q*)
    !
    ! Part  8 ================================
    !                     J * rhs
    !
    ! Part  9 ================================
    ! H)               -H^T * qs
    !
    ! Part 10 ================================
    !               darcy
    !
    ! Part 11 ================================
    !               (A⁻1)^T * rhs








!!! ======================     Part 1    =============================================!
    ! B.1)  -B1^T * dx_1(q*)
    ! B.2)  -B2^T * dx_2(q*)
    ! B.3)  -B3^T * dx_3(q*)

  ! B1) ===============================================================================
    !   -B1^T * dx_1(q*)
    !   -B1^T = phi * [ -u_1,   -u_1*u_1,  -u_1*u_2,  -u_1*u_3,                      0]
    !                 [ -rho, -2*rho*u_1,  -rho*u_2,  -rho*u_3,   -gamma*p/(gamma - 1)]
    !                 [    0,          0,  -rho*u_1,         0,                      0]
    !                 [    0,          0,         0,  -rho*u_1,                      0]
    !                 [    0,          0,         0,         0, -gamma*u_1/(gamma - 1)]

    rhs(:,:,:,1) = rhs(:,:,:,1) + phi * (                 &
                                   - u_1 * qs_x1(:,:,:,1) &
                               - u_1*u_1 * qs_x1(:,:,:,2) &
                               - u_1*u_2 * qs_x1(:,:,:,3) &
                               - u_1*u_3 * qs_x1(:,:,:,4))
    !                                  0 * qs_x1(:,:,:,5)

    rhs(:,:,:,2) = rhs(:,:,:,2) + phi * (                 &
                            - q(:,:,:,1) * qs_x1(:,:,:,1) &
                   - 2.0_rk * q(:,:,:,2) * qs_x1(:,:,:,2) &
                            - q(:,:,:,3) * qs_x1(:,:,:,3) &
                            - q(:,:,:,4) * qs_x1(:,:,:,4) &
             - gamma_INV_gm1* q(:,:,:,5) * qs_x1(:,:,:,5))

    rhs(:,:,:,3) = rhs(:,:,:,3) + phi * (                 &
                            - q(:,:,:,2) * qs_x1(:,:,:,3))

    rhs(:,:,:,4) = rhs(:,:,:,4) + phi * (                 &
                            - q(:,:,:,2) * qs_x1(:,:,:,4))

    rhs(:,:,:,5) = rhs(:,:,:,5) + phi * (                 &
                     - gamma_INV_gm1*u_1 * qs_x1(:,:,:,5))

  ! B2) ============================================================================
    !   -B2^T * dx_2(q*)
    !   -B2^T = phi * [-u_2, -u_1*u_2,   -u_2*u_2, -u_2*u_3,                      0]
    !                 [   0, -rho*u_2,          0,        0,                      0]
    !                 [-rho, -rho*u_1, -2*rho*u_2, -rho*u_3,   -gamma*p/(gamma - 1)]
    !                 [   0,        0,          0, -rho*u_2,                      0]
    !                 [   0,        0,          0,        0, -gamma*u_2/(gamma - 1)]

    rhs(:,:,:,1) = rhs(:,:,:,1) + phi * (                 &
                                   - u_2 * qs_x2(:,:,:,1) &
                               - u_2*u_1 * qs_x2(:,:,:,2) &
                               - u_2*u_2 * qs_x2(:,:,:,3) &
                               - u_2*u_3 * qs_x2(:,:,:,4))
    !                                  0 * qs_x2(:,:,:,5)

    rhs(:,:,:,2) = rhs(:,:,:,2) + phi * (                 &
                            - q(:,:,:,3) * qs_x2(:,:,:,2))

    rhs(:,:,:,3) = rhs(:,:,:,3) + phi * (                 &
                            - q(:,:,:,1) * qs_x2(:,:,:,1) &
                            - q(:,:,:,2) * qs_x2(:,:,:,2) &
                   - 2.0_rk * q(:,:,:,3) * qs_x2(:,:,:,3) &
                            - q(:,:,:,4) * qs_x2(:,:,:,4) &
             - gamma_INV_gm1* q(:,:,:,5) * qs_x2(:,:,:,5))

    rhs(:,:,:,4) = rhs(:,:,:,4) + phi * (                 &
                            - q(:,:,:,3) * qs_x2(:,:,:,4))

    rhs(:,:,:,5) = rhs(:,:,:,5) + phi * (                 &
                     - gamma_INV_gm1*u_2 * qs_x2(:,:,:,5))

  ! B3) ============================================================================
    !   -B3^T * dx_3(q*)
    !   -B3^T = phi * [-u_3, -u_1*u_3, -u_2*u_3,   -u_3*u_3,                      0]
    !                 [   0, -rho*u_3,        0,          0,                      0]
    !                 [   0,        0, -rho*u_3,          0,                      0]
    !                 [-rho, -rho*u_1, -rho*u_2, -2*rho*u_3,   -gamma*p/(gamma - 1)]
    !                 [   0,        0,        0,          0, -gamma*u_3/(gamma - 1)]

    rhs(:,:,:,1) = rhs(:,:,:,1) + phi * (                 &
                                   - u_3 * qs_x3(:,:,:,1) &
                               - u_3*u_1 * qs_x3(:,:,:,2) &
                               - u_3*u_2 * qs_x3(:,:,:,3) &
                               - u_3*u_3 * qs_x3(:,:,:,4))
    !                                  0 * qs_x3(:,:,:,5)

    rhs(:,:,:,2) = rhs(:,:,:,2) + phi * (                 &
                            - q(:,:,:,4) * qs_x3(:,:,:,2))

    rhs(:,:,:,3) = rhs(:,:,:,3) + phi * (                 &
                            - q(:,:,:,4) * qs_x3(:,:,:,3))

    rhs(:,:,:,4) = rhs(:,:,:,4) + phi * (                 &
                            - q(:,:,:,1) * qs_x3(:,:,:,1) &
                            - q(:,:,:,2) * qs_x3(:,:,:,2) &
                            - q(:,:,:,3) * qs_x3(:,:,:,3) &
                   - 2.0_rk * q(:,:,:,4) * qs_x3(:,:,:,4) &
             - gamma_INV_gm1* q(:,:,:,5) * qs_x3(:,:,:,5))

    rhs(:,:,:,5) = rhs(:,:,:,5) + phi * (                 &
                     - gamma_INV_gm1*u_3 * qs_x3(:,:,:,5))

!!! ======================     Part 2    ============================================!
   ! C1)        -D_x1(C1^T * q*)
   ! C2)        -D_x2(C2^T * q*)
   ! C3)        -D_x3(C3^T * q*)

  ! C1) ================================
    ! -D_x1(C1^T * q*)
    !   -C1^T = phi * [0,  0, 0, 0,   0]
    !                 [0,  0, 0, 0,   0]
    !                 [0,  0, 0, 0,   0]
    !                 [0,  0, 0, 0,   0]
    !                 [0, -1, 0, 0, u_1]

  ! call Dx1(a_x,   0*qs(:,:,:,5) ,params%geom%dx1) = 0
  ! call Dx1(b_x,   0*qs(:,:,:,5) ,params%geom%dx1) = 0
  ! call Dx1(c_x,   0*qs(:,:,:,5) ,params%geom%dx1) = 0
  ! call Dx1(d_x,   0*qs(:,:,:,5) ,params%geom%dx1) = 0
  ! call Dx1(e_x, phi * (- qs(:,:,:,2) + u_1*qs(:,:,:,5)),params%geom%dx1) --> grid transformed

  ! call Dx1(a_x, phi * (- qs(:,:,:,2) + u_1 * qs(:,:,:,5)),params%geom%dx1)
  ! call Dx2(b_x, phi * (- qs(:,:,:,2) + u_1 * qs(:,:,:,5)),params%geom%dx2)
  ! call Dx3(c_x, phi * (- qs(:,:,:,2) + u_1 * qs(:,:,:,5)),params%geom%dx3)

            d_x = phi * (- qs(:,:,:,2) + u_1 * qs(:,:,:,5))
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,5) = rhs(:,:,:,5) + X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x

  ! C2) =========================
    ! -D_x2(C2^T * q*)
    !  -C2^T = phi * [0, 0,  0, 0,   0]
    !                [0, 0,  0, 0,   0]
    !                [0, 0,  0, 0,   0]
    !                [0, 0,  0, 0,   0]
    !                [0, 0, -1, 0, u_2]

 ! call Dx2(a_x,   0*qs(:,:,:,5) ,params%geom%dx2) = 0
 ! call Dx2(b_x,   0*qs(:,:,:,5) ,params%geom%dx2) = 0
 ! call Dx2(c_x,   0*qs(:,:,:,5) ,params%geom%dx2) = 0
 ! call Dx2(d_x,   0*qs(:,:,:,5) ,params%geom%dx2) = 0
 ! call Dx2(e_x, phi * (- qs(:,:,:,3) + u_2 * qs(:,:,:,5)),params%geom%dx2)  --> grid transformed

 !  call Dx1(a_x, phi * (- qs(:,:,:,3) + u_2 * qs(:,:,:,5)),params%geom%dx1)
 !  call Dx2(b_x, phi * (- qs(:,:,:,3) + u_2 * qs(:,:,:,5)),params%geom%dx2)
 !  call Dx3(c_x, phi * (- qs(:,:,:,3) + u_2 * qs(:,:,:,5)),params%geom%dx3)

           d_x =  phi * (- qs(:,:,:,3) + u_2 * qs(:,:,:,5))
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)

   rhs(:,:,:,5) = rhs(:,:,:,5) + X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x

  ! C3) ========================
    !   -D_x3(C3^T * q*)
    !   -C3^T = phi * [0, 0, 0,  0,   0]
    !                 [0, 0, 0,  0,   0]
    !                 [0, 0, 0,  0,   0]
    !                 [0, 0, 0,  0,   0]
    !                 [0, 0, 0, -1, u_3]

 ! call Dx3(a_x,   0*qs(:,:,:,5) ,params%geom%dx3) = 0
 ! call Dx3(b_x,   0*qs(:,:,:,5) ,params%geom%dx3) = 0
 ! call Dx3(c_x,   0*qs(:,:,:,5) ,params%geom%dx3) = 0
 ! call Dx3(d_x,   0*qs(:,:,:,5) ,params%geom%dx3) = 0
 ! call Dx3(e_x, phi * (- qs(:,:,:,4) + u_3 * qs(:,:,:,5)),params%geom%dx3)

 ! call Dx1(a_x, phi * (- qs(:,:,:,4) + u_3 * qs(:,:,:,5)),params%geom%dx1)
 ! call Dx2(b_x, phi * (- qs(:,:,:,4) + u_3 * qs(:,:,:,5)),params%geom%dx2)
 ! call Dx3(c_x, phi * (- qs(:,:,:,4) + u_3 * qs(:,:,:,5)),params%geom%dx3)

          d_x =  phi * (- qs(:,:,:,4) + u_3 * qs(:,:,:,5))
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,5) = rhs(:,:,:,5) + X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x

   !!! ======================     Part 8    =========================================!
    ! C1* part (no derivatives of q*)
    ! C1*)              +C1* * Dx_1(q)
    ! C2*)              +C2* * Dx_2(q)
    ! C3*)              +C3* * Dx_3(q)

    ! C1*)
    ! +C1* * Dx_1(q)
    ! +C1* = phi * [0, 0, 0, 0,   0]
    !              [0, 0, 0, 0, -ps]
    !              [0, 0, 0, 0,   0]
    !              [0, 0, 0, 0,   0]
    !              [0, 0, 0, 0,   0]
   call Dx1(a_x, q(:,:,:,5), params%geom%dx1)
   call Dx2(b_x, q(:,:,:,5), params%geom%dx2)
   call Dx3(c_x, q(:,:,:,5), params%geom%dx3)
   rhs(:,:,:,2) = rhs(:,:,:,2) - phi * qs(:,:,:,5) * (X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x)

    ! C2*)
    ! +C2* * Dx_2(q)
    !  C2* = phi * [0, 0, 0, 0,   0]
    !              [0, 0, 0, 0,   0]
    !              [0, 0, 0, 0, -ps]
    !              [0, 0, 0, 0,   0]
    !              [0, 0, 0, 0,   0]
   rhs(:,:,:,3) = rhs(:,:,:,3) - phi * qs(:,:,:,5) * (X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x)

    ! C3*)
    ! +C3* * Dx_3(q)
    !  C3*= phi * [0, 0, 0, 0,   0]
    !             [0, 0, 0, 0,   0]
    !             [0, 0, 0, 0,   0]
    !             [0, 0, 0, 0, -ps]
    !             [0, 0, 0, 0,   0]
   rhs(:,:,:,4) = rhs(:,:,:,4) - phi * qs(:,:,:,5) * (X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x)


!!! ======================     Part 4    =========================================!
    ! D11)       -D_x1(D11^T * dx_1(q*))
    ! D12)       -D_x2(D12^T * dx_1(q*))
    ! D13)       -D_x3(D13^T * dx_1(q*))
    ! D21)       -D_x1(D21^T * dx_2(q*))
    ! D22)       -D_x2(D22^T * dx_2(q*))
    ! D23)       -D_x3(D23^T * dx_2(q*))
    ! D31)       -D_x1(D31^T * dx_3(q*))
    ! D32)       -D_x2(D32^T * dx_3(q*))
    ! D33)       -D_x3(D33^T * dx_3(q*))

  ! ========   dx1(q*) =========================
  ! D11) =======================================
    !  -D_x1(D11^T * dx_1(q*))
    !   -D11^T = phi * [0,       0,   0,   0, 0]
    !                  [0, -4*mu/3,   0,   0, 0]
    !                  [0,       0, -mu,   0, 0]
    !                  [0,       0,   0, -mu, 0]
    !                  [0,       0,   0,   0, 0]

  ! call Dx1(a_x,                        0 * qs_x1(:,:,:,1),params%geom%dx1) = 0
  ! call Dx1(b_x,-4.0_rk/3.0_rk * phi * mu * qs_x1(:,:,:,2),params%geom%dx1) --> grid transformed
  ! call Dx1(c_x,               - phi * mu * qs_x1(:,:,:,3),params%geom%dx1) --> grid transformed
  ! call Dx1(d_x,               - phi * mu * qs_x1(:,:,:,4),params%geom%dx1) --> grid transformed
  ! call Dx1(e_x,                        0 * qs_x1(:,:,:,5),params%geom%dx1) = 0

 ! rhs(:,:,:,1) = rhs(:,:,:,1) + 0
 ! call Dx1(a_x, -4.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x1(:,:,:,2),params%geom%dx1)
 ! call Dx2(b_x, -4.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x1(:,:,:,2),params%geom%dx2)
 ! call Dx3(c_x, -4.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x1(:,:,:,2),params%geom%dx3)
          d_x =  -4.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x1(:,:,:,2)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,2) = rhs(:,:,:,2) +  X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x
 ! call Dx1(a_x,                - phi * mu * qs_x1(:,:,:,3),params%geom%dx1)
 ! call Dx2(b_x,                - phi * mu * qs_x1(:,:,:,3),params%geom%dx2)
 ! call Dx3(c_x,                - phi * mu * qs_x1(:,:,:,3),params%geom%dx3)
          d_x =                 - phi * mu * qs_x1(:,:,:,3)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,3) = rhs(:,:,:,3) + X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x
 ! call Dx1(a_x,                - phi * mu * qs_x1(:,:,:,4),params%geom%dx1)
 ! call Dx2(b_x,                - phi * mu * qs_x1(:,:,:,4),params%geom%dx2)
 ! call Dx3(c_x,                - phi * mu * qs_x1(:,:,:,4),params%geom%dx3)
          d_x =                 - phi * mu * qs_x1(:,:,:,4)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,4) = rhs(:,:,:,4) + X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x
 ! rhs(:,:,:,5) = rhs(:,:,:,5) + 0


  ! D12) ===================================
    !  -D_x2(D12^T * dx_1(q*))
    !  -D12^T = phi * [0,      0,   0, 0, 0]
    !                 [0,      0, -mu, 0, 0]
    !                 [0, 2*mu/3,   0, 0, 0]
    !                 [0,      0,   0, 0, 0]
    !                 [0,      0,   0, 0, 0]

  ! call Dx2(a_x,                        0 * qs_x1(:,:,:,1),params%geom%dx2) = 0
  ! call Dx2(b_x,               - phi * mu * qs_x1(:,:,:,3),params%geom%dx2) --> grid transformed
  ! call Dx2(c_x, 2.0_rk/3.0_rk * phi * mu * qs_x1(:,:,:,2),params%geom%dx2) --> grid transformed
  ! call Dx2(d_x,                        0 * qs_x1(:,:,:,4),params%geom%dx2) = 0
  ! call Dx2(e_x,                        0 * qs_x1(:,:,:,5),params%geom%dx2) = 0

 ! rhs(:,:,:,1) = rhs(:,:,:,1) + 0
 ! call Dx1(a_x,               - phi * mu * X_jacobian * qs_x1(:,:,:,3),params%geom%dx1) ! qs index is 3 not 2 as implied by matrix-vector multiplication
 ! call Dx2(b_x,               - phi * mu * X_jacobian * qs_x1(:,:,:,3),params%geom%dx2) ! qs index is 3 not 2 as implied by matrix-vector multiplication
 ! call Dx3(c_x,               - phi * mu * X_jacobian * qs_x1(:,:,:,3),params%geom%dx3) ! qs index is 3 not 2 as implied by matrix-vector multiplication
          d_x =                - phi * mu * X_jacobian * qs_x1(:,:,:,3)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,2) = rhs(:,:,:,2) + X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x
 ! call Dx1(a_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x1(:,:,:,2),params%geom%dx1) ! qs index is 2 not 3 as implied by matrix-vector multiplication
 ! call Dx2(b_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x1(:,:,:,2),params%geom%dx2) ! qs index is 2 not 3 as implied by matrix-vector multiplication
 ! call Dx3(c_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x1(:,:,:,2),params%geom%dx3) ! qs index is 2 not 3 as implied by matrix-vector multiplication
          d_x =  2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x1(:,:,:,2)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,3) = rhs(:,:,:,3) + X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x
 ! rhs(:,:,:,4) = rhs(:,:,:,4) + 0
 ! rhs(:,:,:,5) = rhs(:,:,:,5) + 0


  ! D13) ====================================
    ! -D_x3(D13^T * dx_1(q*))
    !   -D13^T  = phi * [0,      0, 0,   0, 0]
    !                   [0,      0, 0, -mu, 0]
    !                   [0,      0, 0,   0, 0]
    !                   [0, 2*mu/3, 0,   0, 0]
    !                   [0,      0, 0,   0, 0]

  ! call Dx3(a_x,                        0 * qs_x1(:,:,:,1),params%geom%dx3) = 0
  ! call Dx3(b_x,               - phi * mu * qs_x1(:,:,:,4),params%geom%dx3) --> grid transformed
  ! call Dx3(c_x                         0 * qs_x1(:,:,:,3),params%geom%dx3) = 0
  ! call Dx3(d_x, 2.0_rk/3.0_rk * phi * mu * qs_x1(:,:,:,2),params%geom%dx3) --> grid transformed
  ! call Dx3(e_x,                        0 * qs_x1(:,:,:,5),params%geom%dx3) = 0

 ! rhs(:,:,:,1) = rhs(:,:,:,1) + 0
 ! call Dx1(a_x,               - phi * mu * X_jacobian * qs_x1(:,:,:,4),params%geom%dx1) ! qs index is 4 not 2 as implied by matrix-vector multiplication
 ! call Dx2(b_x,               - phi * mu * X_jacobian * qs_x1(:,:,:,4),params%geom%dx2) ! qs index is 4 not 2 as implied by matrix-vector multiplication
 ! call Dx3(c_x,               - phi * mu * X_jacobian * qs_x1(:,:,:,4),params%geom%dx3) ! qs index is 4 not 2 as implied by matrix-vector multiplication
          d_x =                - phi * mu * X_jacobian * qs_x1(:,:,:,4)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,2) = rhs(:,:,:,2) + X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x
 ! rhs(:,:,:,3) = rhs(:,:,:,3) + 0
 ! call Dx1(a_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x1(:,:,:,2),params%geom%dx1) ! qs index is 2 not 4 as implied by matrix-vector multiplication
 ! call Dx2(b_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x1(:,:,:,2),params%geom%dx2) ! qs index is 2 not 4 as implied by matrix-vector multiplication
 ! call Dx3(c_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x1(:,:,:,2),params%geom%dx3) ! qs index is 2 not 4 as implied by matrix-vector multiplication
          d_x =  2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x1(:,:,:,2)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,4) = rhs(:,:,:,4) + X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x
 ! rhs(:,:,:,5) = rhs(:,:,:,5) + 0

  ! ========   dx2(q*) =========================
  ! D21) =======================================
    !  -D_x1(D21^T * dx_2(q*))
    !     -D21^T = phi * [ 0,   0,      0, 0, 0]
    !                    [ 0,   0, 2/3*mu, 0, 0]
    !                    [ 0, -mu,      0, 0, 0]
    !                    [ 0,   0,      0, 0, 0]
    !                    [ 0,   0,      0, 0, 0]

  ! call Dx1(a_x,                        0 * qs_x2(:,:,:,1),params%geom%dx1) = 0
  ! call Dx1(b_x, 2.0_rk/3.0_rk * phi * mu * qs_x2(:,:,:,3),params%geom%dx1) --> grid transformed
  ! call Dx1(c_x,               - phi * mu * qs_x2(:,:,:,2),params%geom%dx1) --> grid transformed
  ! call Dx1(d_x,                        0 * qs_x2(:,:,:,4),params%geom%dx1) = 0
  ! call Dx1(e_x,                        0 * qs_x2(:,:,:,5),params%geom%dx1) = 0

 ! rhs(:,:,:,1) = rhs(:,:,:,1) + 0
 ! call Dx1(a_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x2(:,:,:,3),params%geom%dx1) ! qs index is 3 not 2 as implied by matrix-vector multiplication
 ! call Dx2(b_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x2(:,:,:,3),params%geom%dx2) ! qs index is 3 not 2 as implied by matrix-vector multiplication
 ! call Dx3(c_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x2(:,:,:,3),params%geom%dx3) ! qs index is 3 not 2 as implied by matrix-vector multiplication
          d_x =  2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x2(:,:,:,3)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,2) = rhs(:,:,:,2) + X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x
 ! call Dx1(a_x,               - phi * mu * X_jacobian * qs_x2(:,:,:,2),params%geom%dx1) ! qs index is 2 not 3 as implied by matrix-vector multiplication
 ! call Dx2(b_x,               - phi * mu * X_jacobian * qs_x2(:,:,:,2),params%geom%dx2) ! qs index is 2 not 3 as implied by matrix-vector multiplication
 ! call Dx3(c_x,               - phi * mu * X_jacobian * qs_x2(:,:,:,2),params%geom%dx3) ! qs index is 2 not 3 as implied by matrix-vector multiplication
          d_x =                - phi * mu * X_jacobian * qs_x2(:,:,:,2)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,3) = rhs(:,:,:,3) + X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x
 ! rhs(:,:,:,4) = rhs(:,:,:,4) + 0
 ! rhs(:,:,:,5) = rhs(:,:,:,5) + 0


 ! D22) ========================================
    !  -D_x2(D22^T * dx_2(q*))
    !   -D22^T = phi * [0,   0,       0,   0, 0]
    !                  [0, -mu,       0,   0, 0]
    !                  [0,   0, -4/3*mu,   0, 0]
    !                  [0,   0,       0, -mu, 0]
    !                  [0,   0,       0,   0, 0]

  ! call Dx2(a_x,                         0 * qs_x2(:,:,:,1),params%geom%dx2) = 0
  ! call Dx2(b_x,                - phi * mu * qs_x2(:,:,:,2),params%geom%dx2) --> grid transformed
  ! call Dx2(c_x, -4.0_rk/3.0_rk * phi * mu * qs_x2(:,:,:,3),params%geom%dx2) --> grid transformed
  ! call Dx2(d_x,                - phi * mu * qs_x2(:,:,:,4),params%geom%dx2) --> grid transformed
  ! call Dx2(e_x,                         0 * qs_x2(:,:,:,5),params%geom%dx2) = 0

 ! rhs(:,:,:,1) = rhs(:,:,:,1) + 0
 ! call Dx1(a_x,                - phi * mu * X_jacobian * qs_x2(:,:,:,2),params%geom%dx1)
 ! call Dx2(b_x,                - phi * mu * X_jacobian * qs_x2(:,:,:,2),params%geom%dx2)
 ! call Dx3(c_x,                - phi * mu * X_jacobian * qs_x2(:,:,:,2),params%geom%dx3)
          d_x =                 - phi * mu * X_jacobian * qs_x2(:,:,:,2)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,2) = rhs(:,:,:,2) + X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x
 ! call Dx1(a_x, -4.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x2(:,:,:,3),params%geom%dx1)
 ! call Dx2(b_x, -4.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x2(:,:,:,3),params%geom%dx2)
 ! call Dx3(c_x, -4.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x2(:,:,:,3),params%geom%dx3)
          d_x =  -4.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x2(:,:,:,3)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,3) = rhs(:,:,:,3) + X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x
 ! call Dx1(a_x,                - phi * mu * X_jacobian * qs_x2(:,:,:,4),params%geom%dx1)
 ! call Dx2(b_x,                - phi * mu * X_jacobian * qs_x2(:,:,:,4),params%geom%dx2)
 ! call Dx3(c_x,                - phi * mu * X_jacobian * qs_x2(:,:,:,4),params%geom%dx3)
          d_x =                 - phi * mu * X_jacobian * qs_x2(:,:,:,4)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,4) = rhs(:,:,:,4) + X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x
 ! rhs(:,:,:,5) = rhs(:,:,:,5) + 0


  ! D23) ===================================
    ! -D_x3(D23^T * dx_2(q*))
    !  -D23^T = phi * [0, 0,      0,   0, 0]
    !                 [0, 0,      0,   0, 0]
    !                 [0, 0,      0, -mu, 0]
    !                 [0, 0, 2/3*mu,   0, 0]
    !                 [0, 0,      0,   0, 0]

  ! call Dx1(a_x,                        0 * qs_x2(:,:,:,1),params%geom%dx1) = 0
  ! call Dx1(b_x,                        0 * qs_x2(:,:,:,3),params%geom%dx1) = 0
  ! call Dx1(c_x,               - phi * mu * qs_x2(:,:,:,4),params%geom%dx1) --> grid transformed
  ! call Dx1(d_x, 2.0_rk/3.0_rk * phi * mu * qs_x2(:,:,:,3),params%geom%dx1) --> grid transformed
  ! call Dx1(e_x,                        0 * qs_x2(:,:,:,5),params%geom%dx1) = 0

 ! rhs(:,:,:,1) = rhs(:,:,:,1) + 0
 ! rhs(:,:,:,2) = rhs(:,:,:,2) + 0
 ! call Dx1(a_x,               - phi * mu * X_jacobian * qs_x2(:,:,:,4),params%geom%dx1) ! qs index is 4 not 3 as implied by matrix-vector multiplication
 ! call Dx2(b_x,               - phi * mu * X_jacobian * qs_x2(:,:,:,4),params%geom%dx2) ! qs index is 4 not 3 as implied by matrix-vector multiplication
 ! call Dx3(c_x,               - phi * mu * X_jacobian * qs_x2(:,:,:,4),params%geom%dx3) ! qs index is 4 not 3 as implied by matrix-vector multiplication
          d_x =                - phi * mu * X_jacobian * qs_x2(:,:,:,4)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,3) = rhs(:,:,:,3) + X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x
 ! call Dx1(a_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x2(:,:,:,3),params%geom%dx1) ! qs index is 3 not 4 as implied by matrix-vector multiplication
 ! call Dx2(b_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x2(:,:,:,3),params%geom%dx2) ! qs index is 3 not 4 as implied by matrix-vector multiplication
 ! call Dx3(c_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x2(:,:,:,3),params%geom%dx3) ! qs index is 3 not 4 as implied by matrix-vector multiplication
          d_x =  2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x2(:,:,:,3)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,4) = rhs(:,:,:,4) + X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x
 ! rhs(:,:,:,5) = rhs(:,:,:,5) + 0


  ! ========   dx3(q*) ========================
  ! D31) ======================================
    !  -D_x1(D31^T * dx_3(q*))
    !     -D31^T = phi * [0,   0, 0,      0, 0]
    !                    [0,   0, 0, 2/3*mu, 0]
    !                    [0,   0, 0,      0, 0]
    !                    [0, -mu, 0,      0, 0]
    !                    [0,   0, 0,      0, 0]

  ! call Dx1(a_x,                        0 * qs_x3(:,:,:,1),params%geom%dx1) = 0
  ! call Dx1(b_x, 2.0_rk/3.0_rk * phi * mu * qs_x3(:,:,:,4),params%geom%dx1) --> grid transformed
  ! call Dx1(c_x,                        0 * qs_x3(:,:,:,4),params%geom%dx1) = 0
  ! call Dx1(d_x,               - phi * mu * qs_x3(:,:,:,2),params%geom%dx1) --> grid transformed
  ! call Dx1(e_x,                        0 * qs_x3(:,:,:,5),params%geom%dx1) = 0

 ! rhs(:,:,:,1) = rhs(:,:,:,1) + 0
 ! call Dx1(a_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x3(:,:,:,4),params%geom%dx1) ! qs index is 4 not 2 as implied by matrix-vector multiplication
 ! call Dx2(b_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x3(:,:,:,4),params%geom%dx2) ! qs index is 4 not 2 as implied by matrix-vector multiplication
 ! call Dx3(c_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x3(:,:,:,4),params%geom%dx3) ! qs index is 4 not 2 as implied by matrix-vector multiplication
          d_x =  2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x3(:,:,:,4)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,2) = rhs(:,:,:,2) + X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x
 ! rhs(:,:,:,3) = rhs(:,:,:,3) + 0
 ! call Dx1(a_x,               - phi * mu * X_jacobian * qs_x3(:,:,:,2),params%geom%dx1) ! qs index is 2 not 4 as implied by matrix-vector multiplication
 ! call Dx2(b_x,               - phi * mu * X_jacobian * qs_x3(:,:,:,2),params%geom%dx2) ! qs index is 2 not 4 as implied by matrix-vector multiplication
 ! call Dx3(c_x,               - phi * mu * X_jacobian * qs_x3(:,:,:,2),params%geom%dx3) ! qs index is 2 not 4 as implied by matrix-vector multiplication
          d_x =                - phi * mu * X_jacobian * qs_x3(:,:,:,2)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,4) = rhs(:,:,:,4) + X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x
 ! rhs(:,:,:,5) = rhs(:,:,:,5) + 0

 ! D32) ====================================
    ! -D_x2(D32^T * dx_3(q*))
    !  -D32^T = phi * [0, 0,   0,      0, 0]
    !                 [0, 0,   0,      0, 0]
    !                 [0, 0,   0, 2/3*mu, 0]
    !                 [0, 0, -mu,      0, 0]
    !                 [0, 0,   0,      0, 0]

  ! call Dx2(a_x,                        0 * qs_x3(:,:,:,1),params%geom%dx2) = 0
  ! call Dx2(b_x,                        0 * qs_x3(:,:,:,3),params%geom%dx2) = 0
  ! call Dx2(c_x, 2.0_rk/3.0_rk * phi * mu * qs_x3(:,:,:,4),params%geom%dx2) --> grid transformed
  ! call Dx2(d_x,               - phi * mu * qs_x3(:,:,:,3),params%geom%dx2) --> grid transformed
  ! call Dx2(e_x,                        0 * qs_x3(:,:,:,5),params%geom%dx2) = 0

 ! rhs(:,:,:,1) = rhs(:,:,:,1) + 0
 ! rhs(:,:,:,2) = rhs(:,:,:,2) + 0
 ! call Dx1(a_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x3(:,:,:,4),params%geom%dx1) ! qs index is 4 not 3 as implied by matrix-vector multiplication
 ! call Dx2(b_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x3(:,:,:,4),params%geom%dx2) ! qs index is 4 not 3 as implied by matrix-vector multiplication
 ! call Dx3(c_x, 2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x3(:,:,:,4),params%geom%dx3) ! qs index is 4 not 3 as implied by matrix-vector multiplication
          d_x =  2.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x3(:,:,:,4)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,3) = rhs(:,:,:,3) + X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x
 ! call Dx1(a_x,               - phi * mu * X_jacobian * qs_x3(:,:,:,3),params%geom%dx1) ! qs index is 3 not 4 as implied by matrix-vector multiplication
 ! call Dx2(b_x,               - phi * mu * X_jacobian * qs_x3(:,:,:,3),params%geom%dx2) ! qs index is 3 not 4 as implied by matrix-vector multiplication
 ! call Dx3(c_x,               - phi * mu * X_jacobian * qs_x3(:,:,:,3),params%geom%dx3) ! qs index is 3 not 4 as implied by matrix-vector multiplication
          d_x =                - phi * mu * X_jacobian * qs_x3(:,:,:,3)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,4) = rhs(:,:,:,4) + X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x
 ! rhs(:,:,:,5) = rhs(:,:,:,5) + 0

 ! D33) ======================================
   !  -D_x3(D33^T * dx_3(q*))
   !   -D33^T = phi * [0,   0,   0,       0, 0]
   !                  [0, -mu,   0,       0, 0]
   !                  [0,   0, -mu,       0, 0]
   !                  [0,   0,   0, -4/3*mu, 0]
   !                  [0,   0,   0,       0, 0]

  ! call Dx3(a_x,                          0 * qs_x3(:,:,:,1),params%geom%dx3) = 0
  ! call Dx3(b_x,                  -phi * mu * qs_x3(:,:,:,2),params%geom%dx3) --> grid transformed
  ! call Dx3(c_x,                  -phi * mu * qs_x3(:,:,:,3),params%geom%dx3) --> grid transformed
  ! call Dx3(d_x, - 4.0_rk/3.0_rk * phi * mu * qs_x3(:,:,:,4),params%geom%dx3) --> grid transformed
  ! call Dx3(e_x,                          0 * qs_x3(:,:,:,5),params%geom%dx3) = 0

 ! rhs(:,:,:,1) = rhs(:,:,:,1) + 0
 ! call Dx1(a_x,                  - phi * mu * X_jacobian * qs_x3(:,:,:,2),params%geom%dx1)
 ! call Dx2(b_x,                  - phi * mu * X_jacobian * qs_x3(:,:,:,2),params%geom%dx2)
 ! call Dx3(c_x,                  - phi * mu * X_jacobian * qs_x3(:,:,:,2),params%geom%dx3)
          d_x =                   - phi * mu * X_jacobian * qs_x3(:,:,:,2)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,2) = rhs(:,:,:,2) + X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x
 ! call Dx1(a_x,                 - phi * mu * X_jacobian * qs_x3(:,:,:,3),params%geom%dx1)
 ! call Dx2(b_x,                 - phi * mu * X_jacobian * qs_x3(:,:,:,3),params%geom%dx2)
 ! call Dx3(c_x,                 - phi * mu * X_jacobian * qs_x3(:,:,:,3),params%geom%dx3)
          d_x =                  - phi * mu * X_jacobian * qs_x3(:,:,:,3)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,3) = rhs(:,:,:,3) + X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x
 ! call Dx1(a_x, - 4.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x3(:,:,:,4),params%geom%dx1)
 ! call Dx2(b_x, - 4.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x3(:,:,:,4),params%geom%dx2)
 ! call Dx3(c_x, - 4.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x3(:,:,:,4),params%geom%dx3)
          d_x =  - 4.0_rk/3.0_rk * phi * mu * X_jacobian * qs_x3(:,:,:,4)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,4) = rhs(:,:,:,4) + X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x
 ! rhs(:,:,:,5) = rhs(:,:,:,5) + 0

!!! ======================     Part 5    =========================================!
    ! E1)             +E1^T * dx_1(q*)
    ! E2)             +E2^T * dx_2(q*)
    ! E3)             +E3^T * dx_3(q*)
 ! E1)
    ! +E1^T * dx_1(q*)
    !  E1^T = W/R*d_T(mu)* phi *
    !   [0, -p/rho**2*(4/3*dx_1(u_1) - 2/3*dx_2(u_2) - 2/3*dx_3(u_3)), -p/rho**2*(dx_2(u_1) + dx_1(u_2)), -p/rho**2*(dx_3(u_1) + dx_1(u_3)), -C_p/Pr*p/(rho**2)*dx_1(W/R*p/rho)]
    !   [0, 0, 0, 0, 0]
    !   [0, 0, 0, 0, 0]
    !   [0, 0, 0, 0, 0]
    !   [0, (4*dx_1(u_1)/3 - 2/3*dx_2(u_2) - 2/3)*dx_3(u_3)/rho, (dx_2(u_1) + dx_1(u_2))/rho, (dx_3(u_1) + dx_1(u_3))/rho, C_p/Pr*1/rho*dx_1(W/R*p/rho))]

 ! call Dx1(a_x, W_INV_R * p_INV_rho, params%geom%dx1)  --> grid transformed
 ! call Dx1(a_x, W_INV_R * p_INV_rho,params%geom%dx1)
 ! call Dx2(b_x, W_INV_R * p_INV_rho,params%geom%dx2)
 ! call Dx3(c_x, W_INV_R * p_INV_rho,params%geom%dx3)
          d_x =  W_INV_R * p_INV_rho
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   e_x = X_jacobian * (X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x)
   rhs(:,:,:,1) = rhs(:,:,:,1)  + W_INV_R * phi * dmu_dT * (                                    & ! 0 * qs_x(:,:,:,1)
         - p_INV_rho_squared * 2.0_rk/3.0_rk * (2.0_rk * u_x - v_y - w_z) * qs_x1(:,:,:,2) &
                                        - p_INV_rho_squared * (u_y + v_x) * qs_x1(:,:,:,3) &
                                        - p_INV_rho_squared * (u_z + w_x) * qs_x1(:,:,:,4) &
                                        - Cp/Pr * p_INV_rho_squared * e_x * qs_x1(:,:,:,5))

   rhs(:,:,:,5) = rhs(:,:,:,5) + W_INV_R * phi * dmu_dT * (                                     & ! 0 * qs_x(:,:,:,1)
                   + INV_rho * 2.0_rk/3.0_rk * (2.0_rk * u_x - v_y - w_z) * qs_x1(:,:,:,2) &
                                                  + INV_rho * (u_y + v_x) * qs_x1(:,:,:,3) &
                                                  + INV_rho * (u_z + w_x) * qs_x1(:,:,:,4) &
                                                  + Cp/Pr * INV_rho * e_x * qs_x1(:,:,:,5))
  ! E2)
    ! +E2^T * dx_2(q*)
    !  E2^T = W/R * d_T(mu) * phi *
    !   [0, -p/rho**2*(dx_2(u_1) + dx_1(u_2)), -p/rho**2 * 2/3*(-dx_1(u_1) + 2*dx_2(u_2) - dx_3(u_3)), -p/rho**2*(dx_3(u_2) + dx_2(u_3)), -C_p/Pr*p/rho**2*dx_2(W*p/(R*rho))]
    !   [0, 0, 0, 0, 0]
    !   [0, 0, 0, 0, 0]
    !   [0, 0, 0, 0, 0]
    !   [0, (dx_2(u_1) + dx_1(u_2))/rho, (-2*dx_1(u_1)/3 + 4*dx_2(u_2)/3 - 2*dx_3(u_3)/3)/rho, (dx_3(u_2) + dx_2(u_3))/rho, C_p/Pr*1//rho*dx_2(W*p/(R*rho))]

  ! call Dx2(a_x, W_INV_R * p_INV_rho, params%geom%dx2)  --> grid transformed
    e_x = X_jacobian * (X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x)
    rhs(:,:,:,1) = rhs(:,:,:,1)  + W_INV_R * phi * dmu_dT * (                                                    & ! 0 * qs_x(:,:,:,1)
                                                          - p_INV_rho_squared * (u_y + v_x) * qs_x2(:,:,:,2) &
                         - p_INV_rho_squared * 2.0_rk/3.0_rk * (- u_x + 2.0_rk * v_y - w_z) * qs_x2(:,:,:,3) &
                                                          - p_INV_rho_squared * (v_z + w_y) * qs_x2(:,:,:,4) &
                                                          - Cp/Pr * p_INV_rho_squared * e_x * qs_x2(:,:,:,5))

    rhs(:,:,:,5) = rhs(:,:,:,5) + W_INV_R * phi * dmu_dT * (                                                     & ! 0 * qs_x(:,:,:,1)
                                                                    + INV_rho * (u_y + v_x) * qs_x2(:,:,:,2) &
                                   + INV_rho * 2.0_rk/3.0_rk * (- u_x + 2.0_rk * v_y - w_z) * qs_x2(:,:,:,3) &
                                                                    + INV_rho * (v_z + w_y) * qs_x2(:,:,:,4) &
                                                                    + Cp/Pr * INV_rho * e_x * qs_x2(:,:,:,5))


 ! E3)
    ! +E3^T * dx_3(q*)
    !  E3^T = W/R * d_T(mu) * phi *
    !   [0, -p/rho**2*(dx_3(u_1) + dx_1(u_3)), -p/rho**2*(dx_3(u_2) + dx_2(u_3)), -p/rho**2*2/3*(-dx_1(u_1) - dx_2(u_2) + 2*dx_3(u_3)), -C_p*p*dx_3(W*p/(R*rho))/(Pr*rho**2)],
    !   [0, 0, 0, 0, 0]
    !   [0, 0, 0, 0, 0]
    !   [0, 0, 0, 0, 0]
    !   [0, (dx_3(u_1) + dx_1(u_3))/rho, (dx_3(u_2) + dx_2(u_3))/rho, (-2*dx_1(u_1)/3 - 2*dx_2(u_2)/3 + 4*dx_3(u_3)/3)/rho, C_p*dx_3(W*p/(R*rho))/(Pr*rho)]

    ! call Dx3(a_x, W_INV_R * p_INV_rho, params%geom%dx3)  --> grid transformed
    e_x = X_jacobian * (X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x)
    rhs(:,:,:,1) = rhs(:,:,:,1)  + W_INV_R * phi * dmu_dT * (                                                    & ! 0 * qs_x(:,:,:,1)
                                                          - p_INV_rho_squared * (u_z + w_x) * qs_x3(:,:,:,2) &
                                                          - p_INV_rho_squared * (v_z + w_y) * qs_x3(:,:,:,3) &
                         - p_INV_rho_squared * 2.0_rk/3.0_rk * (- u_x - v_y + 2.0_rk * w_z) * qs_x3(:,:,:,4) &
                                                          - Cp/Pr * p_INV_rho_squared * e_x * qs_x3(:,:,:,5))

    rhs(:,:,:,5) = rhs(:,:,:,5) + W_INV_R * phi * dmu_dT * (                                                     & ! 0 * qs_x(:,:,:,1)
                                                                    + INV_rho * (u_z + w_x) * qs_x3(:,:,:,2) &
                                                                    + INV_rho * (v_z + w_y) * qs_x3(:,:,:,3) &
                                   + INV_rho * 2.0_rk/3.0_rk * (- u_x - v_y + 2.0_rk * w_z) * qs_x3(:,:,:,4) &
                                                                    + Cp/Pr * INV_rho * e_x * qs_x3(:,:,:,5))

!!! ======================     Part 6    =========================================!
    ! F1) -FII_T*D_x1(F_I^T * dx_1(q*))
    ! F2) -FII_T*D_x2(F_I^T * dx_2(q*))
    ! F3) -FII_T*D_x3(F_I^T * dx_3(q*))
  ! F1)
    ! -FII_T*D_x1(F_I^T * dx_1(q*))
    ! -F_II^T = [0, 0, 0, 0, W*p/(R*rho**2)]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,      -W/(R*rho)]

    !  F_I^T = phi *  [0, 0, 0, 0,      0]
    !                 [0, 0, 0, 0,      0]
    !                 [0, 0, 0, 0,      0]
    !                 [0, 0, 0, 0,      0]
    !                 [0, 0, 0, 0, lambda]

 ! call Dx1(e_x,              phi * lambda * qs_x1(:,:,:,5),params%geom%dx1)  --> grid transformed
 ! call Dx1(a_x, phi * lambda * X_jacobian * qs_x1(:,:,:,5),params%geom%dx1)
 ! call Dx2(b_x, phi * lambda * X_jacobian * qs_x1(:,:,:,5),params%geom%dx2)
 ! call Dx3(c_x, phi * lambda * X_jacobian * qs_x1(:,:,:,5),params%geom%dx3)
          d_x = phi * lambda * X_jacobian * qs_x1(:,:,:,5)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   e_x = X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x
    rhs(:,:,:,1) = rhs(:,:,:,1) + W_INV_R * p_INV_rho_squared * e_x
    rhs(:,:,:,5) = rhs(:,:,:,5)           - W_INV_R * INV_rho * e_x


 ! F2)
    ! -FII_T*D_x2(F_I^T * dx_2(q*))
    ! -F_II^T = [0, 0, 0, 0, W*p/(R*rho**2)]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,      -W/(R*rho)]

    !  F_I^T = phi * [0, 0, 0, 0,      0]
    !                [0, 0, 0, 0,      0]
    !                [0, 0, 0, 0,      0]
    !                [0, 0, 0, 0,      0]
    !                [0, 0, 0, 0, lambda]

 ! call Dx2(e_x,              phi * lambda * qs_x2(:,:,:,5),params%geom%dx2)  --> grid transformed
 ! call Dx1(a_x, phi * lambda * X_jacobian * qs_x2(:,:,:,5),params%geom%dx1)
 ! call Dx2(b_x, phi * lambda * X_jacobian * qs_x2(:,:,:,5),params%geom%dx2)
 ! call Dx3(c_x, phi * lambda * X_jacobian * qs_x2(:,:,:,5),params%geom%dx3)
          d_x =  phi * lambda * X_jacobian * qs_x2(:,:,:,5)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
    e_x = X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x
    rhs(:,:,:,1) = rhs(:,:,:,1) + W_INV_R * p_INV_rho_squared * e_x
    rhs(:,:,:,5) = rhs(:,:,:,5)           - W_INV_R * INV_rho * e_x


 ! F3)
    ! -FII_T*D_x3(F_I^T * dx_3(q*))
    ! -F_II^T = [0, 0, 0, 0, W*p/(R*rho**2)]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,      -W/(R*rho)]

    !  F_I^T =  phi * [0, 0, 0, 0,      0]
    !                 [0, 0, 0, 0,      0]
    !                 [0, 0, 0, 0,      0]
    !                 [0, 0, 0, 0,      0]
    !                 [0, 0, 0, 0, lambda]

 ! call Dx3(e_x,              lambda * qs_x3(:,:,:,5),params%geom%dx3)  --> grid transformed
 ! call Dx1(a_x, lambda * X_jacobian * qs_x3(:,:,:,5),params%geom%dx1)
 ! call Dx2(b_x, lambda * X_jacobian * qs_x3(:,:,:,5),params%geom%dx2)
 ! call Dx3(c_x, lambda * X_jacobian * qs_x3(:,:,:,5),params%geom%dx3)
          d_x =  lambda * X_jacobian * qs_x3(:,:,:,5)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
    e_x = X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x
    rhs(:,:,:,1) = rhs(:,:,:,1) + W_INV_R * p_INV_rho_squared * e_x
    rhs(:,:,:,5) = rhs(:,:,:,5)           - W_INV_R * INV_rho * e_x


!!! ======================     Part 7    =========================================!
    ! G1)        +D_x1(G1^T * q*)
    ! G1)        +D_x2(G2^T * q*)
    ! G1)        +D_x3(G3^T * q*)
  ! G1)
    ! +D_x1(G1^T * q*)
    ! G1^T = phi * mu * [0, 0, 0, 0,                                                      0]
    !                   [0, 0, 0, 0, -4/3*(dx_1(u_1) + dx_2(u_2) + dx_3(u_3)) + 2*dx_1(u_1)]
    !                   [0, 0, 0, 0,                           2 * ( dx_2(u_1) + dx_1(u_2))]
    !                   [0, 0, 0, 0,                           2 * ( dx_3(u_1) + dx_1(u_3))]
    !                   [0, 0, 0, 0,                                                      0]

  ! call Dx1(a_x,                                                                     0 * qs(:,:,:,5),params%geom%dx1) = 0
  ! call Dx1(b_x, 2.0_rk * phi * mu * (-2.0_rk/3.0_rk*(u_x + v_y + w_z) + 2.0_rk * u_x) * qs(:,:,:,5),params%geom%dx1)  --> grid transformed
  ! call Dx1(c_x,                                       2.0_rk * phi * mu * (u_y + v_x) * qs(:,:,:,5),params%geom%dx1)  --> grid transformed
  ! call Dx1(d_x,                                       2.0_rk * phi * mu * (u_z + w_x) * qs(:,:,:,5),params%geom%dx1)  --> grid transformed
  ! call Dx1(e_x,                                                                     0 * qs(:,:,:,5),params%geom%dx1) = 0

 ! rhs(:,:,:,1) = rhs(:,:,:,1) + 0
 ! call Dx1(a_x, 2.0_rk * phi * mu * (-2.0_rk/3.0_rk*(u_x + v_y + w_z) + 2.0_rk * u_x) * qs(:,:,:,5),params%geom%dx1)
 ! call Dx2(b_x, 2.0_rk * phi * mu * (-2.0_rk/3.0_rk*(u_x + v_y + w_z) + 2.0_rk * u_x) * qs(:,:,:,5),params%geom%dx2)
 ! call Dx3(c_x, 2.0_rk * phi * mu * (-2.0_rk/3.0_rk*(u_x + v_y + w_z) + 2.0_rk * u_x) * qs(:,:,:,5),params%geom%dx3)
          d_x =  2.0_rk * phi * mu * (-2.0_rk/3.0_rk*(u_x + v_y + w_z) + 2.0_rk * u_x) * qs(:,:,:,5)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,2) = rhs(:,:,:,2) + X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x
 ! call Dx1(a_x,                                       2.0_rk * phi * mu * (u_y + v_x) * qs(:,:,:,5),params%geom%dx1)
 ! call Dx2(b_x,                                       2.0_rk * phi * mu * (u_y + v_x) * qs(:,:,:,5),params%geom%dx2)
 ! call Dx3(c_x,                                       2.0_rk * phi * mu * (u_y + v_x) * qs(:,:,:,5),params%geom%dx3)
          d_x =                                        2.0_rk * phi * mu * (u_y + v_x) * qs(:,:,:,5)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,3) = rhs(:,:,:,3) + X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x
 ! call Dx1(a_x,                                       2.0_rk * phi * mu * (u_z + w_x) * qs(:,:,:,5),params%geom%dx1)
 ! call Dx2(b_x,                                       2.0_rk * phi * mu * (u_z + w_x) * qs(:,:,:,5),params%geom%dx2)
 ! call Dx3(c_x,                                       2.0_rk * phi * mu * (u_z + w_x) * qs(:,:,:,5),params%geom%dx3)
          d_x =                                        2.0_rk * phi * mu * (u_z + w_x) * qs(:,:,:,5)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,4) = rhs(:,:,:,4) + X_metric(:,:,:,1,1) * a_x + X_metric(:,:,:,2,1) * b_x + X_metric(:,:,:,3,1) * c_x
 ! rhs(:,:,:,5) = rhs(:,:,:,5) + 0


 ! G2)
    !+D_x2(G2^T * q*)
    ! G2^T = phi * mu * [0, 0, 0, 0,                                                      0]
    !                   [0, 0, 0, 0,                            2 *( dx_2(u_1) + dx_1(u_2))]
    !                   [0, 0, 0, 0, -4/3*(dx_1(u_1) + dx_2(u_2) + dx_3(u_3)) + 4*dx_2(u_2)]
    !                   [0, 0, 0, 0,                            2 *( dx_3(u_2) + dx_2(u_3))]
    !                   [0, 0, 0, 0,                                                      0]
  ! call Dx2(a_x,                                                            0 * qs(:,:,:,5),params%geom%dx2) = 0
  ! call Dx2(b_x,                              2.0_rk * phi * mu * (u_y + v_x) * qs(:,:,:,5),params%geom%dx2)  --> grid transformed
  ! call Dx2(c_x, phi * mu * (-4.0_rk/3.0_rk*(u_x + v_y + w_z) + 4.0_rk * v_y) * qs(:,:,:,5),params%geom%dx2)  --> grid transformed
  ! call Dx2(d_x,                              2.0_rk * phi * mu * (v_z + w_y) * qs(:,:,:,5),params%geom%dx2)  --> grid transformed
  ! call Dx2(e_x,                                                            0 * qs(:,:,:,5),params%geom%dx2) = 0

 ! rhs(:,:,:,1) = rhs(:,:,:,1) + 0
 ! call Dx1(a_x,                              2.0_rk * phi * mu * (u_y + v_x) * qs(:,:,:,5),params%geom%dx1)
 ! call Dx2(b_x,                              2.0_rk * phi * mu * (u_y + v_x) * qs(:,:,:,5),params%geom%dx2)
 ! call Dx3(c_x,                              2.0_rk * phi * mu * (u_y + v_x) * qs(:,:,:,5),params%geom%dx3)
          d_x =                               2.0_rk * phi * mu * (u_y + v_x) * qs(:,:,:,5)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,2) = rhs(:,:,:,2) + X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x
 ! call Dx1(a_x, phi * mu * (-4.0_rk/3.0_rk*(u_x + v_y + w_z) + 4.0_rk * v_y) * qs(:,:,:,5),params%geom%dx1)
 ! call Dx2(b_x, phi * mu * (-4.0_rk/3.0_rk*(u_x + v_y + w_z) + 4.0_rk * v_y) * qs(:,:,:,5),params%geom%dx2)
 ! call Dx3(c_x, phi * mu * (-4.0_rk/3.0_rk*(u_x + v_y + w_z) + 4.0_rk * v_y) * qs(:,:,:,5),params%geom%dx3)
          d_x =  phi * mu * (-4.0_rk/3.0_rk*(u_x + v_y + w_z) + 4.0_rk * v_y) * qs(:,:,:,5)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,3) = rhs(:,:,:,3) + X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x
 ! call Dx1(a_x,                              2.0_rk * phi * mu * (v_z + w_y) * qs(:,:,:,5),params%geom%dx1)
 ! call Dx2(b_x,                              2.0_rk * phi * mu * (v_z + w_y) * qs(:,:,:,5),params%geom%dx2)
 ! call Dx3(c_x,                              2.0_rk * phi * mu * (v_z + w_y) * qs(:,:,:,5),params%geom%dx3)
          d_x =                               2.0_rk * phi * mu * (v_z + w_y) * qs(:,:,:,5)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,4) = rhs(:,:,:,4) + X_metric(:,:,:,1,2) * a_x + X_metric(:,:,:,2,2) * b_x + X_metric(:,:,:,3,2) * c_x
 ! rhs(:,:,:,5) = rhs(:,:,:,5) + 0

 ! G3)
    !+D_x3(G3^T * q*)
    ! G3^T = phi * mu * [0, 0, 0, 0,                                                      0]
    !                   [0, 0, 0, 0,                              2*(dx_3(u_1) + dx_1(u_3))]
    !                   [0, 0, 0, 0,                              2*(dx_3(u_2) + dx_2(u_3))]
    !                   [0, 0, 0, 0, -4/3*(dx_1(u_1) + dx_2(u_2) + dx_3(u_3)) + 4*dx_3(u_3)]
    !                   [0, 0, 0, 0,                                                      0]
  ! call Dx3(a_x,                                                            0 * qs(:,:,:,5),params%geom%dx2)
  ! call Dx3(b_x,                              2.0_rk * phi * mu * (u_z + w_x) * qs(:,:,:,5),params%geom%dx3)
  ! call Dx3(c_x,                              2.0_rk * phi * mu * (v_z + w_y) * qs(:,:,:,5),params%geom%dx3)
  ! call Dx3(d_x, phi * mu * (-4.0_rk/3.0_rk*(u_x + v_y + w_z) + 4.0_rk * w_z) * qs(:,:,:,5),params%geom%dx3)
  ! call Dx3(e_x,                                                            0 * qs(:,:,:,5),params%geom%dx3)

 ! rhs(:,:,:,1) = rhs(:,:,:,1) + 0
 ! call Dx1(a_x,                              2.0_rk * phi * mu * (u_z + w_x) * qs(:,:,:,5),params%geom%dx1)
 ! call Dx2(b_x,                              2.0_rk * phi * mu * (u_z + w_x) * qs(:,:,:,5),params%geom%dx2)
 ! call Dx3(c_x,                              2.0_rk * phi * mu * (u_z + w_x) * qs(:,:,:,5),params%geom%dx3)
          d_x =                               2.0_rk * phi * mu * (u_z + w_x) * qs(:,:,:,5)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,2) = rhs(:,:,:,2) + X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x
 ! call Dx1(a_x,                              2.0_rk * phi * mu * (v_z + w_y) * qs(:,:,:,5),params%geom%dx1)
 ! call Dx2(b_x,                              2.0_rk * phi * mu * (v_z + w_y) * qs(:,:,:,5),params%geom%dx2)
 ! call Dx3(c_x,                              2.0_rk * phi * mu * (v_z + w_y) * qs(:,:,:,5),params%geom%dx3)
          d_x =                               2.0_rk * phi * mu * (v_z + w_y) * qs(:,:,:,5)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,3) = rhs(:,:,:,3) + X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x
 ! call Dx1(a_x, phi * mu * (-4.0_rk/3.0_rk*(u_x + v_y + w_z) + 4.0_rk * w_z) * qs(:,:,:,5),params%geom%dx1)
 ! call Dx2(b_x, phi * mu * (-4.0_rk/3.0_rk*(u_x + v_y + w_z) + 4.0_rk * w_z) * qs(:,:,:,5),params%geom%dx2)
 ! call Dx3(c_x, phi * mu * (-4.0_rk/3.0_rk*(u_x + v_y + w_z) + 4.0_rk * w_z) * qs(:,:,:,5),params%geom%dx3)
          d_x =  phi * mu * (-4.0_rk/3.0_rk*(u_x + v_y + w_z) + 4.0_rk * w_z) * qs(:,:,:,5)
   call Dx1(a_x, d_x, params%geom%dx1)
   call Dx2(b_x, d_x, params%geom%dx2)
   call Dx3(c_x, d_x, params%geom%dx3)
   rhs(:,:,:,4) = rhs(:,:,:,4) + X_metric(:,:,:,1,3) * a_x + X_metric(:,:,:,2,3) * b_x + X_metric(:,:,:,3,3) * c_x
 ! rhs(:,:,:,5) = rhs(:,:,:,5) + 0

!!! ======================     Part 8    =========================================!
 ! Acounting for Jacobian
     rhs(:,:,:,1) =  X_jacobian * rhs(:,:,:,1)
     rhs(:,:,:,2) =  X_jacobian * rhs(:,:,:,2)
     rhs(:,:,:,3) =  X_jacobian * rhs(:,:,:,3)
     rhs(:,:,:,4) =  X_jacobian * rhs(:,:,:,4)
     rhs(:,:,:,5) =  X_jacobian * rhs(:,:,:,5)


!!! ======================     Part 9    =========================================!

   !  -H^T * qs
    !
    ! e_x =  W/R*d_T(mu)*phi * ( 4/3*(dx_1(u_1)*dx_1(u_1) + dx_2(u_2)*dx_2(u_2) + dx_3(u_3)*dx_3(u_3)
    !                               - dx_1(u_1)*dx_2(u_2) - dx_1(u_1)*dx_3(u_3) - dx_2(u_2)*dx_3(u_3))
    !                           + 4* (dx_2(u_1)*dx_1(u_2) + dx_3(u_1)*dx_1(u_3) + dx_3(u_2)*dx_2(u_3)))
    !  h15 =  +p/rho**2 * e_x
    !
    !  h55 =   -1/rho *  e_x
    !  -H^T = [0, 0, 0, 0, h15]
    !         [0, 0, 0, 0,   0]
    !         [0, 0, 0, 0,   0]
    !         [0, 0, 0, 0,   0]
    !         [0, 0, 0, 0, h55]
    e_x =  W_INV_R * dmu_dT * phi * (                                      &
              4.0_rk/3.0_rk * (  u_x * u_x + v_y * v_y + w_z * w_z   &
                               - u_x * v_y - u_x * w_z - v_y * w_z)  &
                   + 4.0_rk * (  u_y * v_x + u_z * w_x + v_z * w_x))

    rhs(:,:,:,1) = rhs(:,:,:,1) + p_INV_rho_squared * e_x * qs(:,:,:,5)
    rhs(:,:,:,5) = rhs(:,:,:,5)           - INV_rho * e_x * qs(:,:,:,5)
!!! ======================     Part 10   =========================================!
    ! Darcy term (accounting for (A⁻1)^T follows at the end)
    rhs(:,:,:,1) = rhs(:,:,:,1) + penalization(:,:,:,3)*qs(:,:,:,1) ! addtional darcy type for density
    rhs(:,:,:,2) = rhs(:,:,:,2) + penalization(:,:,:,4)*qs(:,:,:,2) ! classical darcy in u_1
    rhs(:,:,:,3) = rhs(:,:,:,3) + penalization(:,:,:,5)*qs(:,:,:,3) ! classical darcy in u_2
    rhs(:,:,:,4) = rhs(:,:,:,4) + penalization(:,:,:,6)*qs(:,:,:,4) ! classical darcy in u_3
    rhs(:,:,:,5) = rhs(:,:,:,5) + penalization(:,:,:,7)*qs(:,:,:,5) ! addtional darcy type for pressure

!!! ======================     Part 11   =========================================!

    !    (A⁻1)^T * rhs
    !    (A⁻1)^T =  [1, -u_1/rho, -u_2/rho, -u_3/rho,         0]
    !               [0,    1/rho,        0,        0,         0]
    !               [0,        0,    1/rho,        0,         0]
    !               [0,        0,        0,    1/rho,         0]
    !               [0,        0,        0,        0, gamma - 1]

    a_x          =                       rhs(:,:,:,1) &
                        - u_1 *INV_rho * rhs(:,:,:,2) &
                        - u_2 *INV_rho * rhs(:,:,:,3) &
                        - u_3 *INV_rho * rhs(:,:,:,4)
    !                                0 * rhs(:,:,:,5)

    b_x          =             INV_rho * rhs(:,:,:,2)
    !                                0 * rhs(:,:,:,3)
    !                                0 * rhs(:,:,:,4)
    !                                0 * rhs(:,:,:,5)

    c_x          =             INV_rho * rhs(:,:,:,3)
    !                                0 * rhs(:,:,:,4)
    !                                0 * rhs(:,:,:,5)

    d_x          =             INV_rho * rhs(:,:,:,4)
    !                                0 * rhs(:,:,:,5)

    e_x          =                 gm1 * rhs(:,:,:,5)


    rhs(:,:,:,1) = a_x
    rhs(:,:,:,2) = b_x
    rhs(:,:,:,3) = c_x
    rhs(:,:,:,4) = d_x
    rhs(:,:,:,5) = e_x

!!! adjoint force
! Note: Although this technically represents -g, the sign convention is handled
!       within calc_objective_g to maintain consistency with existing user cases.
!       The multiplication by (A^T)^(-1) is an additional step already applied within g,
!       so this line effectively corresponds to rhs - (A^T)^(-1) * g.
    rhs = rhs + g

!!! adjoint sponge
!     call sponge_adjoint(rhs,qs,q)

!!! boundary conditions 2/2: rhss
    call set_boundary_condition_rhss(rhs,qs,q)

    deallocate(qs)
    deallocate(qs_x1)
    deallocate(qs_x2)
    deallocate(qs_x3)

    deallocate(u_1)
    deallocate(u_2)
    deallocate(u_3)

    deallocate(gamma_DOT_p_INV_rho_gm1)
    deallocate(      p_INV_rho_squared)
    deallocate(              p_INV_rho)
    deallocate(                INV_rho)

    deallocate(a_x)
    deallocate(b_x)
    deallocate(c_x)
    deallocate(d_x)
    deallocate(e_x)

    deallocate(qs1_x)
    deallocate(qs2_x)
    deallocate(qs3_x)
    deallocate(qs4_x)
    deallocate(qs5_x)

    deallocate(qs1_y)
    deallocate(qs2_y)
    deallocate(qs3_y)
    deallocate(qs4_y)
    deallocate(qs5_y)

    deallocate(qs1_z)
    deallocate(qs2_z)
    deallocate(qs3_z)
    deallocate(qs4_z)
    deallocate(qs5_z)

    deallocate(u_x)
    deallocate(u_y)
    deallocate(u_z)

    deallocate(v_x)
    deallocate(v_y)
    deallocate(v_z)

    deallocate(w_x)
    deallocate(w_y)
    deallocate(w_z)

    deallocate(  temp)
    deallocate(    mu)
    deallocate(dmu_dT)
    deallocate(lambda)

    deallocate(   phi)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine right_hand_side_adjoint_3d
!!!=================================================================================================

end module ns_rho_rhou_p_ev_phi_grid_transform
