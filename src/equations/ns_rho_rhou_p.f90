module ns_rho_rhou_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! information/remarks:                                                                             !
  ! - q_in, qs_in in subroutine calls, this is needed due to intent(in) of q/qs in time stepper,     !
  !   thus no modification inside of the rhs is permitted, however no intent works fine, even        !
  !   nobody understands why ... (2019-01-02, ML)                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use boundary_conditions
  use data_reference , only: qref
  use discretisation_x
  use ns_rho_rhou_p_trafos
  use force
  use parameter
!   use sponge_layer

  private

  public right_hand_side_direct_3d
  public right_hand_side_adjoint_3d

  public give_T, give_c, give_r, give_u, give_p

  character(len=max_length_parameter), parameter, public   :: equations_name = "ns_rho_rhou_p"
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
    real(kind=rk), dimension(:,:,:)  , allocatable   :: a_x, a_y, a_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: p_x, p_y, p_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: u_x, u_y, u_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: v_x, v_y, v_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: w_x, w_y, w_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: tau_11, tau_12, tau_13
    real(kind=rk), dimension(:,:,:)  , allocatable   :: tau_22, tau_23
    real(kind=rk), dimension(:,:,:)  , allocatable   :: tau_33
    real(kind=rk), dimension(:,:,:)  , allocatable   :: temp,mu,mu_div,lambda
    integer                                          :: i
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
    allocate(lambda(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

!!! set q
    q = q_in

!!! boundary conditions 1/2: q
    call set_boundary_condition_q(q)

!!! unpack variables
    do i = 1,params%geom%n3b
       u(:,:,i,1) = q(:,:,i,2)/q(:,:,i,1)
       u(:,:,i,2) = q(:,:,i,3)/q(:,:,i,1)
       u(:,:,i,3) = q(:,:,i,4)/q(:,:,i,1)
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

!!! mass equation --------------------------------------------------------------------------------
    call Dx1(a_x,q(:,:,:,2),params%geom%dx1)
    call Dx2(a_y,q(:,:,:,3),params%geom%dx2)
    call Dx3(a_z,q(:,:,:,4),params%geom%dx3)

    rhs(:,:,:,1) = -a_x - a_y - a_z

!!! momentum equation (1) ------------------------------------------------------------------------
    call Dx1(a_x,q(:,:,:,2)*u(:,:,:,1),params%geom%dx1)
    call Dx2(a_y,q(:,:,:,3)*u(:,:,:,1),params%geom%dx2)
    call Dx3(a_z,q(:,:,:,4)*u(:,:,:,1),params%geom%dx3)
    call Dx1(p_x,q(:,:,:,5)           ,params%geom%dx1)

    rhs(:,:,:,2) = -a_x - a_y - a_z - p_x

    ! friction
    call Dx1(a_x,tau_11,params%geom%dx1) ! d_x tau_11
    call Dx2(a_y,tau_12,params%geom%dx2) ! d_y tau_12
    call Dx3(a_z,tau_13,params%geom%dx3) ! d_z tau_13
    rhs(:,:,:,2) = rhs(:,:,:,2) + a_x + a_y + a_z

!!! momentum equation (2) ------------------------------------------------------------------------
    call Dx1(a_x,q(:,:,:,2)*u(:,:,:,2),params%geom%dx1)
    call Dx2(a_y,q(:,:,:,3)*u(:,:,:,2),params%geom%dx2)
    call Dx3(a_z,q(:,:,:,4)*u(:,:,:,2),params%geom%dx3)
    call Dx2(p_y,q(:,:,:,5)           ,params%geom%dx2)

    rhs(:,:,:,3) = -a_x - a_y - a_z - p_y

    ! friction
    call Dx1(a_x,tau_12,params%geom%dx1) ! d_x tau_21
    call Dx2(a_y,tau_22,params%geom%dx2) ! d_y tau_22
    call Dx3(a_z,tau_23,params%geom%dx3) ! d_z tau_23
    rhs(:,:,:,3) = rhs(:,:,:,3) + a_x + a_y + a_z

!!! momentum equation (3) ------------------------------------------------------------------------
    call Dx1(a_x,q(:,:,:,2)*u(:,:,:,3),params%geom%dx1)
    call Dx2(a_y,q(:,:,:,3)*u(:,:,:,3),params%geom%dx2)
    call Dx3(a_z,q(:,:,:,4)*u(:,:,:,3),params%geom%dx3)
    call Dx3(p_z,q(:,:,:,5)           ,params%geom%dx3)

    rhs(:,:,:,4) = -a_x - a_y - a_z - p_z

    ! friction
    call Dx1(a_x,tau_13,params%geom%dx1) ! d_x tau_31
    call Dx2(a_y,tau_23,params%geom%dx2) ! d_y tau_32
    call Dx3(a_z,tau_33,params%geom%dx3) ! d_z tau_33
    rhs(:,:,:,4) = rhs(:,:,:,4) + a_x + a_y + a_z

!!! energy (pressure) equation -------------------------------------------------------------------
    call Dx1(a_x,q(:,:,:,5)*u(:,:,:,1),params%geom%dx1)
    call Dx2(a_y,q(:,:,:,5)*u(:,:,:,2),params%geom%dx2)
    call Dx3(a_z,q(:,:,:,5)*u(:,:,:,3),params%geom%dx3)

    gamma_term = params%material%gamma/(params%material%gamma - 1.0_rk) ! gamma /(gamma -1)
    rhs(:,:,:,5) = gamma_term*(-a_x - a_y - a_z)
    rhs(:,:,:,5) = rhs(:,:,:,5) + u(:,:,:,1)*p_x(:,:,:) + u(:,:,:,2)*p_y(:,:,:) + u(:,:,:,3)*p_z(:,:,:)

    ! heat transfer
    call Dx1(a_x,temp) ! d_x T
    call Dx2(a_y,temp) ! d_y T
    call Dx3(a_z,temp) ! d_z T

    lambda = give_lambda(mu)

    call Dx1(a_x,lambda*a_x) ! d_x (lambda d_x T)
    call Dx2(a_y,lambda*a_y) ! d_y (lambda d_y T)
    call Dx3(a_z,lambda*a_z) ! d_z (lambda d_z T)

    ! heat term
    rhs(:,:,:,5) = rhs(:,:,:,5) + a_x + a_y + a_z


    ! friction term
    rhs(:,:,:,5) = rhs(:,:,:,5) + tau_11*u_x + tau_12*u_y + tau_13*u_z & ! 11 12 13
                                + tau_12*v_x + tau_22*v_y + tau_23*v_z & ! 21 22 23
                                + tau_13*w_x + tau_23*w_y + tau_33*w_z   ! 31 32 33

    gamma_term = (params%material%gamma - 1.0_rk)   ! (gamma -1)
    rhs(:,:,:,5) = gamma_term*rhs(:,:,:,5)

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
    real(kind=rk), dimension(:,:,:,:), allocatable   :: qs, qs_x
    real(kind=rk), dimension(:,:,:)  , allocatable   :: u_1,u_2,u_3
    real(kind=rk), dimension(:,:,:)  , allocatable   :: u_x, u_y, u_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: v_x, v_y, v_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: w_x, w_y, w_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: a_x,b_x,c_x,d_x,e_x
    real(kind=rk), dimension(:,:,:)  , allocatable   :: mu,dmu_dT,lambda, temp
    real(kind=rk), dimension(:,:,:)  , allocatable   :: gamma_DOT_p_INV_rho_gm1,INV_rho,p_INV_rho,p_INV_rho_squared
    real(kind=rk)                                    :: gm1,gamma_INV_gm1
    real(kind=rk)                                    :: W,R,W_INV_R,Cp,Pr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(qs(params%geom%n1b,params%geom%n2b,params%geom%n3b,5))
    allocate(qs_x(params%geom%n1b,params%geom%n2b,params%geom%n3b,5))

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

    !!! compute derivatives used multiple times
    call Dx1(u_x,u_1,params%geom%dx1)
    call Dx2(u_y,u_1,params%geom%dx2)
    call Dx3(u_z,u_1,params%geom%dx3)

    call Dx1(v_x,u_2,params%geom%dx1)
    call Dx2(v_y,u_2,params%geom%dx2)
    call Dx3(v_z,u_2,params%geom%dx3)

    call Dx1(w_x,u_3,params%geom%dx1)
    call Dx2(w_y,u_3,params%geom%dx2)
    call Dx3(w_z,u_3,params%geom%dx3)

!!! Sequence: (corresponding to ???)
    !      A). D_x1 part
    !      B). D_x2 part
    !      C). D_x3 part
    !      D). ...
    !      E). (A^-1)^T*rhs

!!! A) D_x1 part
    ! A.1)             -B1^T * dx_1(q*)
    ! A.2)        -D_x1(C1^T * q*)
    ! A.3)       -D_x1(D11^T * dx_1(q*))
    ! A.4)       -D_x2(D12^T * dx_1(q*))
    ! A.5)       -D_x3(D13^T * dx_1(q*))
    ! A.6)             +E1^T * dx_1(q*)
    ! A.7) -FII_T*D_x1(F_I^T * dx_1(q*))
    ! A.8)        +D_x1(G1^T * q*)


    !  qs_x = dx_1(q*)
    call Dx1(qs_x(:,:,:,1),qs(:,:,:,1),params%geom%dx1)
    call Dx1(qs_x(:,:,:,2),qs(:,:,:,2),params%geom%dx1)
    call Dx1(qs_x(:,:,:,3),qs(:,:,:,3),params%geom%dx1)
    call Dx1(qs_x(:,:,:,4),qs(:,:,:,4),params%geom%dx1)
    call Dx1(qs_x(:,:,:,5),qs(:,:,:,5),params%geom%dx1)


  ! A.1)
    !   -B1^T * dx_1(q*)
    !   -B1^T = [ -u_1,   -u_1*u_1,  -u_1*u_2,  -u_1*u_3,                      0]
    !           [ -rho, -2*rho*u_1,  -rho*u_2,  -rho*u_3,   -gamma*p/(gamma - 1)]
    !           [    0,          0,  -rho*u_1,         0,                      0]
    !           [    0,          0,         0,  -rho*u_1,                      0]
    !           [    0,          0,         0,         0, -gamma*u_1/(gamma - 1)]

    rhs(:,:,:,1) = rhs(:,:,:,1)                          &
                                   - u_1 * qs_x(:,:,:,1) &
                               - u_1*u_1 * qs_x(:,:,:,2) &
                               - u_1*u_2 * qs_x(:,:,:,3) &
                               - u_1*u_3 * qs_x(:,:,:,4)
    !                                  0 * qs_x(:,:,:,5)

    rhs(:,:,:,2) = rhs(:,:,:,2)                          &
                            - q(:,:,:,1) * qs_x(:,:,:,1) &
                   - 2.0_rk * q(:,:,:,2) * qs_x(:,:,:,2) &
                            - q(:,:,:,3) * qs_x(:,:,:,3) &
                            - q(:,:,:,4) * qs_x(:,:,:,4) &
             - gamma_INV_gm1* q(:,:,:,5) * qs_x(:,:,:,5)

    rhs(:,:,:,3) = rhs(:,:,:,3)                          &
                            - q(:,:,:,2) * qs_x(:,:,:,3)

    rhs(:,:,:,4) = rhs(:,:,:,4)                          &
                            - q(:,:,:,2) * qs_x(:,:,:,4)

    rhs(:,:,:,5) = rhs(:,:,:,5)                          &
                     - gamma_INV_gm1*u_1 * qs_x(:,:,:,5)

  ! A.2)
    ! -D_x1(C1^T * q*)
    !   -C1^T = [0,  0, 0, 0,   0]
    !           [0,  0, 0, 0,   0]
    !           [0,  0, 0, 0,   0]
    !           [0,  0, 0, 0,   0]
    !           [0, -1, 0, 0, u_1]

  ! call Dx1(a_x,   0*qs(:,:,:,5) ,params%geom%dx1)
  ! call Dx1(b_x,   0*qs(:,:,:,5) ,params%geom%dx1)
  ! call Dx1(c_x,   0*qs(:,:,:,5) ,params%geom%dx1)
  ! call Dx1(d_x,   0*qs(:,:,:,5) ,params%geom%dx1)
   call Dx1(e_x,  - qs(:,:,:,2) + u_1*qs(:,:,:,5) ,params%geom%dx1)
   rhs(:,:,:,5) = rhs(:,:,:,5) + e_x

  ! A.3)
    !  -D_x1(D11^T * dx_1(q*))
    !   -D11^T = [0,       0,   0,   0, 0]
    !            [0, -4*mu/3,   0,   0, 0]
    !            [0,       0, -mu,   0, 0]
    !            [0,       0,   0, -mu, 0]
    !            [0,       0,   0,   0, 0]

  ! call Dx1(a_x,                  0 * qs_x(:,:,:,1),params%geom%dx1)
    call Dx1(b_x,-4.0_rk/3.0_rk * mu * qs_x(:,:,:,2),params%geom%dx1)
    call Dx1(c_x,               - mu * qs_x(:,:,:,3),params%geom%dx1)
    call Dx1(d_x,               - mu * qs_x(:,:,:,4),params%geom%dx1)
  ! call Dx1(e_x,                  0 * qs_x(:,:,:,5),params%geom%dx1)

  ! rhs(:,:,:,1) = rhs(:,:,:,1) + a_x
    rhs(:,:,:,2) = rhs(:,:,:,2) + b_x
    rhs(:,:,:,3) = rhs(:,:,:,3) + c_x
    rhs(:,:,:,4) = rhs(:,:,:,4) + d_x
  ! rhs(:,:,:,5) = rhs(:,:,:,5) + e_x


  ! A.4)
    !  -D_x2(D12^T * dx_1(q*))
    !  -D12^T = [0,      0,   0, 0, 0]
    !           [0,      0, -mu, 0, 0]
    !           [0, 2*mu/3,   0, 0, 0]
    !           [0,      0,   0, 0, 0]
    !           [0,      0,   0, 0, 0]

  ! call Dx2(a_x,                  0 * qs_x(:,:,:,1),params%geom%dx2)
    call Dx2(b_x,               - mu * qs_x(:,:,:,3),params%geom%dx2) ! qs index is 3 not 2 as implied by matrix-vector multiplication
    call Dx2(c_x, 2.0_rk/3.0_rk * mu * qs_x(:,:,:,2),params%geom%dx2) ! qs index is 2 not 3 as implied by matrix-vector multiplication
  ! call Dx2(d_x,                  0 * qs_x(:,:,:,4),params%geom%dx2)
  ! call Dx2(e_x,                  0 * qs_x(:,:,:,5),params%geom%dx2)

    rhs(:,:,:,2) = rhs(:,:,:,2) + b_x
    rhs(:,:,:,3) = rhs(:,:,:,3) + c_x


  ! A.5)
    ! -D_x3(D13^T * dx_1(q*))
    !   -D13^T  = [0,      0, 0,   0, 0]
    !             [0,      0, 0, -mu, 0]
    !             [0,      0, 0,   0, 0]
    !             [0, 2*mu/3, 0,   0, 0]
    !             [0,      0, 0,   0, 0]

  ! call Dx3(a_x,                  0 * qs_x(:,:,:,1),params%geom%dx3)
    call Dx3(b_x,               - mu * qs_x(:,:,:,4),params%geom%dx3) ! qs index is 4 not 2 as implied by matrix-vector multiplication
  ! call Dx3(c_x                   0 * qs_x(:,:,:,3),params%geom%dx3)
    call Dx3(d_x, 2.0_rk/3.0_rk * mu * qs_x(:,:,:,2),params%geom%dx3) ! qs index is 2 not 4 as implied by matrix-vector multiplication
  ! call Dx3(e_x,                  0 * qs_x(:,:,:,5),params%geom%dx3)

    rhs(:,:,:,2) = rhs(:,:,:,2) + b_x
    rhs(:,:,:,4) = rhs(:,:,:,4) + d_x

  ! A.6)
    ! +E1^T * dx_1(q*)
    !  E1^T = W/R*d_T(mu)*
    !   [0, -p/rho**2*(4/3*dx_1(u_1) - 2/3*dx_2(u_2) - 2/3*dx_3(u_3)), -p/rho**2*(dx_2(u_1) + dx_1(u_2)), -p/rho**2*(dx_3(u_1) + dx_1(u_3)), -C_p/Pr*p/(rho**2)*dx_1(W/R*p/rho)]
    !   [0, 0, 0, 0, 0]
    !   [0, 0, 0, 0, 0]
    !   [0, 0, 0, 0, 0]
    !   [0, (4*dx_1(u_1)/3 - 2/3*dx_2(u_2) - 2/3)*dx_3(u_3)/rho, (dx_2(u_1) + dx_1(u_2))/rho, (dx_3(u_1) + dx_1(u_3))/rho, C_p/Pr*1/rho*dx_1(W/R*p/rho))]

    call Dx1(a_x, W_INV_R * p_INV_rho, params%geom%dx1)
    rhs(:,:,:,1) = rhs(:,:,:,1)  + W_INV_R *  dmu_dT * (                                  & ! 0 * qs_x(:,:,:,1)
         - p_INV_rho_squared * 2.0_rk/3.0_rk * (2.0_rk * u_x - v_y - w_z) * qs_x(:,:,:,2) &
                                        - p_INV_rho_squared * (u_y + v_x) * qs_x(:,:,:,3) &
                                        - p_INV_rho_squared * (u_z + w_x) * qs_x(:,:,:,4) &
                                        - Cp/Pr * p_INV_rho_squared * a_x * qs_x(:,:,:,5))

    rhs(:,:,:,5) = rhs(:,:,:,5) + W_INV_R *  dmu_dT * (                                   & ! 0 * qs_x(:,:,:,1)
                   + INV_rho * 2.0_rk/3.0_rk * (2.0_rk * u_x - v_y - w_z) * qs_x(:,:,:,2) &
                                                  + INV_rho * (u_y + v_x) * qs_x(:,:,:,3) &
                                                  + INV_rho * (u_z + w_x) * qs_x(:,:,:,4) &
                                                  + Cp/Pr * INV_rho * a_x * qs_x(:,:,:,5))

  ! A.7)
    ! -FII_T*D_x1(F_I^T * dx_1(q*))
    ! -F_II^T = [0, 0, 0, 0, W*p/(R*rho**2)]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,      -W/(R*rho)]

    !  F_I^T =  [0, 0, 0, 0,      0]
    !           [0, 0, 0, 0,      0]
    !           [0, 0, 0, 0,      0]
    !           [0, 0, 0, 0,      0]
    !           [0, 0, 0, 0, lambda]

    call Dx1(e_x, lambda * qs_x(:,:,:,5), params%geom%dx1)
    rhs(:,:,:,1) = rhs(:,:,:,1) + W_INV_R * p_INV_rho_squared * e_x
    rhs(:,:,:,5) = rhs(:,:,:,5)           - W_INV_R * INV_rho * e_x

  ! A.8)
    ! +D_x1(G1^T * q*)
    ! G1^T = mu * [0, 0, 0, 0,                                                      0]
    !             [0, 0, 0, 0, -4/3*(dx_1(u_1) + dx_2(u_2) + dx_3(u_3)) + 2*dx_1(u_1)]
    !             [0, 0, 0, 0,                           2 * ( dx_2(u_1) + dx_1(u_2))]
    !             [0, 0, 0, 0,                           2 * ( dx_3(u_1) + dx_1(u_3))]
    !             [0, 0, 0, 0,                                                      0]

  ! call Dx1(a_x,                                                               0 * qs(:,:,:,5),params%geom%dx1)
    call Dx1(b_x, 2.0_rk * mu * (-2.0_rk/3.0_rk*(u_x + v_y + w_z) + 2.0_rk * u_x) * qs(:,:,:,5),params%geom%dx1)
    call Dx1(c_x,                                       2.0_rk * mu * (u_y + v_x) * qs(:,:,:,5),params%geom%dx1)
    call Dx1(d_x,                                       2.0_rk * mu * (u_z + w_x) * qs(:,:,:,5),params%geom%dx1)
  ! call Dx1(e_x,                                                               0 * qs(:,:,:,5),params%geom%dx1)

  ! rhs(:,:,:,1) = rhs(:,:,:,1) + 0
    rhs(:,:,:,2) = rhs(:,:,:,2) + b_x
    rhs(:,:,:,3) = rhs(:,:,:,3) + c_x
    rhs(:,:,:,4) = rhs(:,:,:,4) + d_x
  ! rhs(:,:,:,5) = rhs(:,:,:,5) + 0

!!! B) D_x2 part
    ! B.1)             -B2^T * dx_2(q*)
    ! B.2)        -D_x2(C2^T * q*)
    ! B.3)       -D_x1(D21^T * dx_2(q*))
    ! B.4)       -D_x2(D22^T * dx_2(q*))
    ! B.5)       -D_x3(D23^T * dx_2(q*))
    ! B.6)             +E2^T * dx_2(q*)
    ! B.7) -FII_T*D_x2(F_I^T * dx_2(q*))´
    ! B.8)        +D_x2(G2^T * q*)


    !  qs_x = dx_2(q*)
    call Dx2(qs_x(:,:,:,1),qs(:,:,:,1),params%geom%dx2)
    call Dx2(qs_x(:,:,:,2),qs(:,:,:,2),params%geom%dx2)
    call Dx2(qs_x(:,:,:,3),qs(:,:,:,3),params%geom%dx2)
    call Dx2(qs_x(:,:,:,4),qs(:,:,:,4),params%geom%dx2)
    call Dx2(qs_x(:,:,:,5),qs(:,:,:,5),params%geom%dx2)

  ! B.1)
    !   -B2^T * dx_2(q*)
    !   -B2^T = [-u_2, -u_1*u_2,   -u_2*u_2, -u_2*u_3,                      0]
    !           [   0, -rho*u_2,          0,        0,                      0]
    !           [-rho, -rho*u_1, -2*rho*u_2, -rho*u_3,   -gamma*p/(gamma - 1)]
    !           [   0,        0,          0, -rho*u_2,                      0]
    !           [   0,        0,          0,        0, -gamma*u_2/(gamma - 1)]

    rhs(:,:,:,1) = rhs(:,:,:,1)                          &
                                   - u_2 * qs_x(:,:,:,1) &
                               - u_2*u_1 * qs_x(:,:,:,2) &
                               - u_2*u_2 * qs_x(:,:,:,3) &
                               - u_2*u_3 * qs_x(:,:,:,4)
    !                                  0 * qs_x(:,:,:,5)

    rhs(:,:,:,2) = rhs(:,:,:,2)                          &
                            - q(:,:,:,3) * qs_x(:,:,:,2)

    rhs(:,:,:,3) = rhs(:,:,:,3)                          &
                            - q(:,:,:,1) * qs_x(:,:,:,1) &
                            - q(:,:,:,2) * qs_x(:,:,:,2) &
                   - 2.0_rk * q(:,:,:,3) * qs_x(:,:,:,3) &
                            - q(:,:,:,4) * qs_x(:,:,:,4) &
             - gamma_INV_gm1* q(:,:,:,5) * qs_x(:,:,:,5)

    rhs(:,:,:,4) = rhs(:,:,:,4)                          &
                            - q(:,:,:,3) * qs_x(:,:,:,4)

    rhs(:,:,:,5) = rhs(:,:,:,5)                          &
                     - gamma_INV_gm1*u_2 * qs_x(:,:,:,5)


  ! B.2)
    ! -D_x2(C2^T * q*)
    !  -C2^T = [0, 0,  0, 0,   0]
    !          [0, 0,  0, 0,   0]
    !          [0, 0,  0, 0,   0]
    !          [0, 0,  0, 0,   0]
    !          [0, 0, -1, 0, u_2]

 ! call Dx2(a_x,   0*qs(:,:,:,5) ,params%geom%dx2)
 ! call Dx2(b_x,   0*qs(:,:,:,5) ,params%geom%dx2)
 ! call Dx2(c_x,   0*qs(:,:,:,5) ,params%geom%dx2)
 ! call Dx2(d_x,   0*qs(:,:,:,5) ,params%geom%dx2)
   call Dx2(e_x, - qs(:,:,:,3) + u_2*qs(:,:,:,5) ,params%geom%dx2)
   rhs(:,:,:,5) = rhs(:,:,:,5) + e_x


  ! B.3)
    !  -D_x1(D21^T * dx_2(q*))
    !     -D21^T = [ 0,   0,      0, 0, 0]
    !              [ 0,   0, 2/3*mu, 0, 0]
    !              [ 0, -mu,      0, 0, 0]
    !              [ 0,   0,      0, 0, 0]
    !              [ 0,   0,      0, 0, 0]

  ! call Dx1(a_x,                  0 * qs_x(:,:,:,1),params%geom%dx1)
    call Dx1(b_x, 2.0_rk/3.0_rk * mu * qs_x(:,:,:,3),params%geom%dx1) ! qs index is 3 not 2 as implied by matrix-vector multiplication
    call Dx1(c_x,               - mu * qs_x(:,:,:,2),params%geom%dx1) ! qs index is 2 not 3 as implied by matrix-vector multiplication
  ! call Dx1(d_x,                  0 * qs_x(:,:,:,4),params%geom%dx1)
  ! call Dx1(e_x,                  0 * qs_x(:,:,:,5),params%geom%dx1)

  ! rhs(:,:,:,1) = rhs(:,:,:,1) + a_x
    rhs(:,:,:,2) = rhs(:,:,:,2) + b_x
    rhs(:,:,:,3) = rhs(:,:,:,3) + c_x
  ! rhs(:,:,:,4) = rhs(:,:,:,4) + d_x
  ! rhs(:,:,:,5) = rhs(:,:,:,5) + e_x

 ! B.4)
    !  -D_x2(D22^T * dx_2(q*))
    !   -D22^T = [0,   0,       0,   0, 0]
    !            [0, -mu,       0,   0, 0]
    !            [0,   0, -4/3*mu,   0, 0]
    !            [0,   0,       0, -mu, 0]
    !            [0,   0,       0,   0, 0]

  ! call Dx2(a_x,                   0 * qs_x(:,:,:,1),params%geom%dx2)
    call Dx2(b_x,                - mu * qs_x(:,:,:,2),params%geom%dx2)
    call Dx2(c_x, -4.0_rk/3.0_rk * mu * qs_x(:,:,:,3),params%geom%dx2)
    call Dx2(d_x,                - mu * qs_x(:,:,:,4),params%geom%dx2)
  ! call Dx2(e_x,                   0 * qs_x(:,:,:,5),params%geom%dx2)

  ! rhs(:,:,:,1) = rhs(:,:,:,1) + a_x
    rhs(:,:,:,2) = rhs(:,:,:,2) + b_x
    rhs(:,:,:,3) = rhs(:,:,:,3) + c_x
    rhs(:,:,:,4) = rhs(:,:,:,4) + d_x
  ! rhs(:,:,:,5) = rhs(:,:,:,5) + e_x


  ! B.5)
    ! -D_x3(D23^T * dx_2(q*))
    !  -D23^T = [0, 0,      0,   0, 0]
    !           [0, 0,      0,   0, 0]
    !           [0, 0,      0, -mu, 0]
    !           [0, 0, 2/3*mu,   0, 0]
    !           [0, 0,      0,   0, 0]

  ! call Dx1(a_x,                  0 * qs_x(:,:,:,1),params%geom%dx1)
  ! call Dx1(b_x,                  0 * qs_x(:,:,:,3),params%geom%dx1)
    call Dx1(c_x,               - mu * qs_x(:,:,:,4),params%geom%dx1) ! qs index is 4 not 3 as implied by matrix-vector multiplication
    call Dx1(d_x, 2.0_rk/3.0_rk * mu * qs_x(:,:,:,3),params%geom%dx1) ! qs index is 3 not 4 as implied by matrix-vector multiplication
  ! call Dx1(e_x,                  0 * qs_x(:,:,:,5),params%geom%dx1)

  ! rhs(:,:,:,1) = rhs(:,:,:,1) + a_x
  ! rhs(:,:,:,2) = rhs(:,:,:,2) + b_x
    rhs(:,:,:,3) = rhs(:,:,:,3) + c_x
    rhs(:,:,:,4) = rhs(:,:,:,4) + d_x
  ! rhs(:,:,:,5) = rhs(:,:,:,5) + e_x


 ! B.6)
    ! +E2^T * dx_2(q*)
    !  E2^T = W/R*d_T(mu)*
    !   [0, -p/rho**2*(dx_2(u_1) + dx_1(u_2)), -p/rho**2 * 2/3*(-dx_1(u_1) + 2*dx_2(u_2) - dx_3(u_3)), -p/rho**2*(dx_3(u_2) + dx_2(u_3)), -C_p/Pr*p/rho**2*dx_2(W*p/(R*rho))]
    !   [0, 0, 0, 0, 0]
    !   [0, 0, 0, 0, 0]
    !   [0, 0, 0, 0, 0]
    !   [0, (dx_2(u_1) + dx_1(u_2))/rho, (-2*dx_1(u_1)/3 + 4*dx_2(u_2)/3 - 2*dx_3(u_3)/3)/rho, (dx_3(u_2) + dx_2(u_3))/rho, C_p/Pr*1//rho*dx_2(W*p/(R*rho))]

    call Dx2(a_x, W_INV_R * p_INV_rho, params%geom%dx2)
    rhs(:,:,:,1) = rhs(:,:,:,1)  + W_INV_R *  dmu_dT * (                                                    & ! 0 * qs_x(:,:,:,1)
                                                          - p_INV_rho_squared * (u_y + v_x) * qs_x(:,:,:,2) &
                         - p_INV_rho_squared * 2.0_rk/3.0_rk * (- u_x + 2.0_rk * v_y - w_z) * qs_x(:,:,:,3) &
                                                          - p_INV_rho_squared * (v_z + w_y) * qs_x(:,:,:,4) &
                                                          - Cp/Pr * p_INV_rho_squared * a_x * qs_x(:,:,:,5))

    rhs(:,:,:,5) = rhs(:,:,:,5) + W_INV_R *  dmu_dT * (                                                     & ! 0 * qs_x(:,:,:,1)
                                                                    + INV_rho * (u_y + v_x) * qs_x(:,:,:,2) &
                                   + INV_rho * 2.0_rk/3.0_rk * (- u_x + 2.0_rk * v_y - w_z) * qs_x(:,:,:,3) &
                                                                    + INV_rho * (v_z + w_y) * qs_x(:,:,:,4) &
                                                                    + Cp/Pr * INV_rho * a_x * qs_x(:,:,:,5))

 ! B.7)
    ! -FII_T*D_x2(F_I^T * dx_2(q*))
    ! -F_II^T = [0, 0, 0, 0, W*p/(R*rho**2)]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,      -W/(R*rho)]

    !  F_I^T =  [0, 0, 0, 0,      0]
    !           [0, 0, 0, 0,      0]
    !           [0, 0, 0, 0,      0]
    !           [0, 0, 0, 0,      0]
    !           [0, 0, 0, 0, lambda]

    call Dx2(e_x, lambda * qs_x(:,:,:,5), params%geom%dx2)
    rhs(:,:,:,1) = rhs(:,:,:,1) + W_INV_R * p_INV_rho_squared * e_x
    rhs(:,:,:,5) = rhs(:,:,:,5)           - W_INV_R * INV_rho * e_x

 ! B.8)
    !+D_x2(G2^T * q*)
    ! G2^T = mu* [0, 0, 0, 0,                                                      0]
    !            [0, 0, 0, 0,                            2 *( dx_2(u_1) + dx_1(u_2))]
    !            [0, 0, 0, 0, -4/3*(dx_1(u_1) + dx_2(u_2) + dx_3(u_3)) + 4*dx_2(u_2)]
    !            [0, 0, 0, 0,                            2 *( dx_3(u_2) + dx_2(u_3))]
    !            [0, 0, 0, 0,                                                      0]
  ! call Dx2(a_x,                                                      0 * qs(:,:,:,5),params%geom%dx2)
    call Dx2(b_x,                              2.0_rk * mu * (u_y + v_x) * qs(:,:,:,5),params%geom%dx2)
    call Dx2(c_x, mu * (-4.0_rk/3.0_rk*(u_x + v_y + w_z) + 4.0_rk * v_y) * qs(:,:,:,5),params%geom%dx2)
    call Dx2(d_x,                              2.0_rk * mu * (v_z + w_y) * qs(:,:,:,5),params%geom%dx2)
  ! call Dx2(e_x,                                                      0 * qs(:,:,:,5),params%geom%dx2)

  ! rhs(:,:,:,1) = rhs(:,:,:,1) + 0
    rhs(:,:,:,2) = rhs(:,:,:,2) + b_x
    rhs(:,:,:,3) = rhs(:,:,:,3) + c_x
    rhs(:,:,:,4) = rhs(:,:,:,4) + d_x
  ! rhs(:,:,:,5) = rhs(:,:,:,5) + 0

!!! C) D_x3 part

    ! C.1)             -B3^T * dx_3(q*)
    ! C.2)        -D_x2(C2^T * q*)
    ! C.3)       -D_x1(D21^T * dx_2(q*))
    ! C.4)       -D_x2(D22^T * dx_2(q*))
    ! C.5)       -D_x3(D23^T * dx_2(q*))
    ! C.6)             +E2^T * dx_2(q*)
    ! C.7) -FII_T*D_x2(F_I^T * dx_2(q*))´
    ! C.8)        +D_x2(G2^T * dx_2(q*))


    !  qs_x = dx_3(q*)
    call Dx3(qs_x(:,:,:,1),qs(:,:,:,1),params%geom%dx3)
    call Dx3(qs_x(:,:,:,2),qs(:,:,:,2),params%geom%dx3)
    call Dx3(qs_x(:,:,:,3),qs(:,:,:,3),params%geom%dx3)
    call Dx3(qs_x(:,:,:,4),qs(:,:,:,4),params%geom%dx3)
    call Dx3(qs_x(:,:,:,5),qs(:,:,:,5),params%geom%dx3)

  ! C.1)
    !   -B3^T * dx_3(q*)
    !   -B3^T = [-u_3, -u_1*u_3, -u_2*u_3,   -u_3*u_3,                      0]
    !           [   0, -rho*u_3,        0,          0,                      0]
    !           [   0,        0, -rho*u_3,          0,                      0]
    !           [-rho, -rho*u_1, -rho*u_2, -2*rho*u_3,   -gamma*p/(gamma - 1)]
    !           [   0,        0,        0,          0, -gamma*u_3/(gamma - 1)]

    rhs(:,:,:,1) = rhs(:,:,:,1)                          &
                                   - u_3 * qs_x(:,:,:,1) &
                               - u_3*u_1 * qs_x(:,:,:,2) &
                               - u_3*u_2 * qs_x(:,:,:,3) &
                               - u_3*u_3 * qs_x(:,:,:,4)
    !                                  0 * qs_x(:,:,:,5)

    rhs(:,:,:,2) = rhs(:,:,:,2)                          &
                            - q(:,:,:,4) * qs_x(:,:,:,2)

    rhs(:,:,:,3) = rhs(:,:,:,3)                          &
                            - q(:,:,:,4) * qs_x(:,:,:,3)

    rhs(:,:,:,4) = rhs(:,:,:,4)                           &
                            - q(:,:,:,1) * qs_x(:,:,:,1) &
                            - q(:,:,:,2) * qs_x(:,:,:,2) &
                            - q(:,:,:,3) * qs_x(:,:,:,3) &
                   - 2.0_rk * q(:,:,:,4) * qs_x(:,:,:,4) &
             - gamma_INV_gm1* q(:,:,:,5) * qs_x(:,:,:,5)

    rhs(:,:,:,5) = rhs(:,:,:,5)                          &
                     - gamma_INV_gm1*u_3 * qs_x(:,:,:,5)



  ! C.2)
    !   -D_x3(C3^T * q*)
    !   -C3^T =[0, 0, 0,  0,   0]
    !          [0, 0, 0,  0,   0]
    !          [0, 0, 0,  0,   0]
    !          [0, 0, 0,  0,   0]
    !          [0, 0, 0, -1, u_3]

 ! call Dx3(a_x,   0*qs(:,:,:,5) ,params%geom%dx3)
 ! call Dx3(b_x,   0*qs(:,:,:,5) ,params%geom%dx3)
 ! call Dx3(c_x,   0*qs(:,:,:,5) ,params%geom%dx3)
 ! call Dx3(d_x,   0*qs(:,:,:,5) ,params%geom%dx3)
   call Dx3(e_x,  - qs(:,:,:,4) + u_3*qs(:,:,:,5) ,params%geom%dx3)
   rhs(:,:,:,5) = rhs(:,:,:,5) + e_x

  ! C.3)
    !  -D_x1(D31^T * dx_3(q*))
    !     -D31^T = [0,   0, 0,      0, 0]
    !              [0,   0, 0, 2/3*mu, 0]
    !              [0,   0, 0,      0, 0]
    !              [0, -mu, 0,      0, 0]
    !              [0,   0, 0,      0, 0]

  ! call Dx1(a_x,                  0 * qs_x(:,:,:,1),params%geom%dx1)
    call Dx1(b_x, 2.0_rk/3.0_rk * mu * qs_x(:,:,:,4),params%geom%dx1) ! qs index is 4 not 2 as implied by matrix-vector multiplication
  ! call Dx1(c_x,                  0 * qs_x(:,:,:,4),params%geom%dx1)
    call Dx1(d_x,               - mu * qs_x(:,:,:,2),params%geom%dx1) ! qs index is 2 not 4 as implied by matrix-vector multiplication
  ! call Dx1(e_x,                  0 * qs_x(:,:,:,5),params%geom%dx1)

  ! rhs(:,:,:,1) = rhs(:,:,:,1) + a_x
    rhs(:,:,:,2) = rhs(:,:,:,2) + b_x
  ! rhs(:,:,:,3) = rhs(:,:,:,3) + c_x
    rhs(:,:,:,4) = rhs(:,:,:,4) + d_x
  ! rhs(:,:,:,5) = rhs(:,:,:,5) + e_x

 ! C.4)
    ! -D_x2(D32^T * dx_3(q*))
    !  -D32^T = [0, 0,   0,      0, 0]
    !           [0, 0,   0,      0, 0]
    !           [0, 0,   0, 2/3*mu, 0]
    !           [0, 0, -mu,      0, 0]
    !           [0, 0,   0,      0, 0]

  ! call Dx2(a_x,                  0 * qs_x(:,:,:,1),params%geom%dx2)
  ! call Dx2(b_x,                  0 * qs_x(:,:,:,3),params%geom%dx2)
    call Dx2(c_x, 2.0_rk/3.0_rk * mu * qs_x(:,:,:,4),params%geom%dx2) ! qs index is 4 not 3 as implied by matrix-vector multiplication
    call Dx2(d_x,               - mu * qs_x(:,:,:,3),params%geom%dx2) ! qs index is 3 not 4 as implied by matrix-vector multiplication
  ! call Dx2(e_x,                  0 * qs_x(:,:,:,5),params%geom%dx2)

  ! rhs(:,:,:,1) = rhs(:,:,:,1) + a_x
  ! rhs(:,:,:,2) = rhs(:,:,:,2) + b_x
    rhs(:,:,:,3) = rhs(:,:,:,3) + c_x
    rhs(:,:,:,4) = rhs(:,:,:,4) + d_x
  ! rhs(:,:,:,5) = rhs(:,:,:,5) + e_x

 ! C.5)
   !  -D_x3(D33^T * dx_3(q*))
   !   -D33^T = [0,   0,   0,       0, 0]
   !            [0, -mu,   0,       0, 0]
   !            [0,   0, -mu,       0, 0]
   !            [0,   0,   0, -4/3*mu, 0]
   !            [0,   0,   0,       0, 0]

  ! call Dx3(a_x,                     0 * qs_x(:,:,:,1),params%geom%dx3)
    call Dx3(b_x,                  - mu * qs_x(:,:,:,2),params%geom%dx3)
    call Dx3(c_x,                  - mu * qs_x(:,:,:,3),params%geom%dx3)
    call Dx3(d_x, - 4.0_rk/3.0_rk *  mu * qs_x(:,:,:,4),params%geom%dx3)
  ! call Dx3(e_x,                     0 * qs_x(:,:,:,5),params%geom%dx3)

  ! rhs(:,:,:,1) = rhs(:,:,:,1) + a_x
    rhs(:,:,:,2) = rhs(:,:,:,2) + b_x
    rhs(:,:,:,3) = rhs(:,:,:,3) + c_x
    rhs(:,:,:,4) = rhs(:,:,:,4) + d_x
  ! rhs(:,:,:,5) = rhs(:,:,:,5) + e_x

 ! C.6)
    ! +E3^T * dx_3(q*)
    !  E3^T = W/R*d_T(mu)*
    !   [0, -p/rho**2*(dx_3(u_1) + dx_1(u_3)), -p/rho**2*(dx_3(u_2) + dx_2(u_3)), -p/rho**2*2/3*(-dx_1(u_1) - dx_2(u_2) + 2*dx_3(u_3)), -C_p*p*dx_3(W*p/(R*rho))/(Pr*rho**2)],
    !   [0, 0, 0, 0, 0]
    !   [0, 0, 0, 0, 0]
    !   [0, 0, 0, 0, 0]
    !   [0, (dx_3(u_1) + dx_1(u_3))/rho, (dx_3(u_2) + dx_2(u_3))/rho, (-2*dx_1(u_1)/3 - 2*dx_2(u_2)/3 + 4*dx_3(u_3)/3)/rho, C_p*dx_3(W*p/(R*rho))/(Pr*rho)]

    call Dx3(a_x, W_INV_R * p_INV_rho, params%geom%dx3)
    rhs(:,:,:,1) = rhs(:,:,:,1)  + W_INV_R *  dmu_dT * (                                                    & ! 0 * qs_x(:,:,:,1)
                                                          - p_INV_rho_squared * (u_z + w_x) * qs_x(:,:,:,2) &
                                                          - p_INV_rho_squared * (v_z + w_y) * qs_x(:,:,:,3) &
                         - p_INV_rho_squared * 2.0_rk/3.0_rk * (- u_x - v_y + 2.0_rk * w_z) * qs_x(:,:,:,4) &
                                                          - Cp/Pr * p_INV_rho_squared * a_x * qs_x(:,:,:,5))

    rhs(:,:,:,5) = rhs(:,:,:,5) + W_INV_R *  dmu_dT * (                                                     & ! 0 * qs_x(:,:,:,1)
                                                                    + INV_rho * (u_z + w_x) * qs_x(:,:,:,2) &
                                                                    + INV_rho * (v_z + w_y) * qs_x(:,:,:,3) &
                                   + INV_rho * 2.0_rk/3.0_rk * (- u_x - v_y + 2.0_rk * w_z) * qs_x(:,:,:,4) &
                                                                    + Cp/Pr * INV_rho * a_x * qs_x(:,:,:,5))

 ! C.7)
    ! -FII_T*D_x3(F_I^T * dx_3(q*))
    ! -F_II^T = [0, 0, 0, 0, W*p/(R*rho**2)]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,               0]
    !           [0, 0, 0, 0,      -W/(R*rho)]

    !  F_I^T =  [0, 0, 0, 0,      0]
    !           [0, 0, 0, 0,      0]
    !           [0, 0, 0, 0,      0]
    !           [0, 0, 0, 0,      0]
    !           [0, 0, 0, 0, lambda]

    call Dx3(e_x, lambda * qs_x(:,:,:,5), params%geom%dx3)
    rhs(:,:,:,1) = rhs(:,:,:,1) + W_INV_R * p_INV_rho_squared * e_x
    rhs(:,:,:,5) = rhs(:,:,:,5)           - W_INV_R * INV_rho * e_x

 ! C.8)
    !+D_x3(G3^T * q*)
    ! G3^T = mu*[0, 0, 0, 0,                                                      0]
    !           [0, 0, 0, 0,                              2*(dx_3(u_1) + dx_1(u_3))]
    !           [0, 0, 0, 0,                              2*(dx_3(u_2) + dx_2(u_3))]
    !           [0, 0, 0, 0, -4/3*(dx_1(u_1) + dx_2(u_2) + dx_3(u_3)) + 4*dx_3(u_3)]
    !           [0, 0, 0, 0,                                                      0]
  ! call Dx3(a_x,                                                      0 * qs(:,:,:,5),params%geom%dx2)
    call Dx3(b_x,                              2.0_rk * mu * (u_z + w_x) * qs(:,:,:,5),params%geom%dx3)
    call Dx3(c_x,                              2.0_rk * mu * (v_z + w_y) * qs(:,:,:,5),params%geom%dx3)
    call Dx3(d_x, mu * (-4.0_rk/3.0_rk*(u_x + v_y + w_z) + 4.0_rk * w_z) * qs(:,:,:,5),params%geom%dx3)
  ! call Dx3(e_x,                                                               0 * qs(:,:,:,5),params%geom%dx3)

  ! rhs(:,:,:,1) = rhs(:,:,:,1) + 0
    rhs(:,:,:,2) = rhs(:,:,:,2) + b_x
    rhs(:,:,:,3) = rhs(:,:,:,3) + c_x
    rhs(:,:,:,4) = rhs(:,:,:,4) + d_x
  ! rhs(:,:,:,5) = rhs(:,:,:,5) + 0

   ! D part (no derivatives of q*)
    ! D.1)              +C1* * Dx_1(q)
    ! D.2)              +C2* * Dx_2(q)
    ! D.3)              +C3* * Dx_3(q)


    ! D.1)
    ! +C1* * Dx_1(q)
    ! +C1* =[0, 0, 0, 0,   0]
    !       [0, 0, 0, 0, -ps]
    !       [0, 0, 0, 0,   0]
    !       [0, 0, 0, 0,   0]
    !       [0, 0, 0, 0,   0]

    ! D.2)
    ! +C2* * Dx_2(q)
    !  C2* = [0, 0, 0, 0,   0]
    !        [0, 0, 0, 0,   0]
    !        [0, 0, 0, 0, -ps]
    !        [0, 0, 0, 0,   0]
    !        [0, 0, 0, 0,   0]

    ! D.3)
    ! +C3* * Dx_3(q)
    !  C3*= [0, 0, 0, 0,   0]
    !       [0, 0, 0, 0,   0]
    !       [0, 0, 0, 0,   0]
    !       [0, 0, 0, 0, -ps]
    !       [0, 0, 0, 0,   0]
   call Dx1(b_x, q(:,:,:,5), params%geom%dx1)
   call Dx2(c_x, q(:,:,:,5), params%geom%dx2)
   call Dx3(d_x, q(:,:,:,5), params%geom%dx3)
    rhs(:,:,:,2) = rhs(:,:,:,2) - qs(:,:,:,5) * b_x
    rhs(:,:,:,3) = rhs(:,:,:,3) - qs(:,:,:,5) * c_x
    rhs(:,:,:,4) = rhs(:,:,:,4) - qs(:,:,:,5) * d_x


   ! D.4) -H^T * qs
    ! -H^T =  W/R*d_T(mu)*
    !  h15 =  +p/rho**2*
    !                  4/3*(dx_1(u_1)*dx_1(u_1) + dx_2(u_2)*dx_2(u_2) + dx_3(u_3)*dx_3(u_3)
    !                     - dx_1(u_1)*dx_2(u_2) - dx_1(u_1)*dx_3(u_3) - dx_2(u_2)*dx_3(u_3))
    !                 + 4* (dx_2(u_1)*dx_1(u_2) + dx_3(u_1)*dx_1(u_3) + dx_3(u_2)*dx_2(u_3))

    !  h55 =   -1/rho   *(4/3*dx_1(u_1) - 2/3*dx_2(u_2) - 2/3*dx_3(u_3))
    !  [0, 0, 0, 0, h15]
    !  [0, 0, 0, 0,   0]
    !  [0, 0, 0, 0,   0]
    !  [0, 0, 0, 0,   0]
    !  [0, 0, 0, 0, h55]
    e_x =  W_INV_R * dmu_dT * (                                      &
              4.0_rk/3.0_rk * (  u_x * u_x + v_y * v_y + w_z * w_z   &
                               - u_x * v_y - u_x * w_z - v_y * w_z)  &
                   + 4.0_rk * (  u_y * v_x + u_z * w_x + v_z * w_x))

    rhs(:,:,:,1) = rhs(:,:,:,1) + p_INV_rho_squared * e_x * qs(:,:,:,5)
    rhs(:,:,:,5) = rhs(:,:,:,5)           - INV_rho * e_x * qs(:,:,:,5)



    ! E) (A⁻1)^T * rhs
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
    deallocate(qs_x)

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine right_hand_side_adjoint_3d
!!!=================================================================================================

end module ns_rho_rhou_p
