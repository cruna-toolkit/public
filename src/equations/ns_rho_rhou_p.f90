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
  use sponge_layer

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
    allocate(qs(params%geom%n1b,params%geom%n2b,params%geom%n3b,5))

    allocate(u_1(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(u_2(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(u_3(params%geom%n1b,params%geom%n2b,params%geom%n3b))

    allocate(gamma_DOT_p_INV_rho_gm1(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(                INV_rho(params%geom%n1b,params%geom%n2b,params%geom%n3b))

    allocate(a_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(b_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(c_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(d_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(e_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))

!!! set qs
    qs = qs_in

!!! boundary conditions 1/2: qs
    call set_boundary_condition_qs(qs,q)

    gm1                     = params%material%gamma - 1.0_rk
    gamma_DOT_p_INV_rho_gm1 = params%material%gamma*q(:,:,:,5)/(q(:,:,:,1) * gm1)
    INV_rho                 = 1.0_rk/q(:,:,:,1)

    rhs = 0.0_rk

    ! unpack direct state
    u_1 = q(:,:,:,2)/q(:,:,:,1)
    u_2 = q(:,:,:,3)/q(:,:,:,1)
    u_3 = q(:,:,:,4)/q(:,:,:,1)

!!! D_x1 part
    call Dx1(a_x,qs(:,:,:,1),params%geom%dx1)
    call Dx1(b_x,qs(:,:,:,2),params%geom%dx1)
    call Dx1(c_x,qs(:,:,:,3),params%geom%dx1)
    call Dx1(d_x,qs(:,:,:,4),params%geom%dx1)
    call Dx1(e_x,qs(:,:,:,5),params%geom%dx1)

    !-(A^T)^(-1)*B1^T * qs_x1
    !   [0,   -u_1**2, -u_1*u_2, -u_1*u_3, -gamma*p*u_1/(rho*(gamma - 1))]
    !   [1,     2*u_1,      u_2,      u_3,      gamma*p/(rho*(gamma - 1))]
    !   [0,         0,      u_1,        0,                              0]
    !   [0,         0,        0,      u_1,                              0]
    !   [0, gamma - 1,        0,        0,                      gamma*u_1]

    rhs(:,:,:,1) = rhs(:,:,:,1)                             &
         +(u_1**2 )                                   * b_x &
         +(u_1*u_2)                                   * c_x &
         +(u_1*u_3)                                   * d_x &
         +gamma_DOT_p_INV_rho_gm1 * u_1               * e_x

    rhs(:,:,:,2) = rhs(:,:,:,2)                             &
         -                                              a_x &
         - 2.0_rk*u_1                                 * b_x &
         -        u_2                                 * c_x &
         -        u_3                                 * d_x &
         - gamma_DOT_p_INV_rho_gm1                    * e_x

    rhs(:,:,:,3) = rhs(:,:,:,3)                             &
         - u_1                                        * c_x

    rhs(:,:,:,4) = rhs(:,:,:,4)                             &
         - u_1                                        * d_x

    rhs(:,:,:,5) = rhs(:,:,:,5)                             &
         - gm1                                        * b_x &
         - params%material%gamma*u_1                  * e_x

    !-(A^T)^(-1) * (C1^T * qs)_x1
    !
    !   [1, -u_1/rho, -u_2/rho, -u_3/rho,         0]
    !   [0,    1/rho,        0,        0,         0]
    !   [0,        0,    1/rho,        0,         0]
    !   [0,        0,        0,    1/rho,         0]
    !   [0,        0,        0,        0, gamma - 1]
    !
    !   [0, 0, 0, 0,    0]
    !   [0, 0, 0, 0,    0]
    !   [0, 0, 0, 0,    0]
    !   [0, 0, 0, 0,    0]
    !   [0, 0, 0, 0, -u_1]

    call Dx1(a_x,-u_1*qs(:,:,:,5),params%geom%dx1)
    rhs(:,:,:,5) = rhs(:,:,:,5)                             &
         -gm1*a_x

    !+(A^T)^(-1) * C1~ * c_x1
    !   [0, 0, 0, 0, ps*u_1/rho]
    !   [0, 0, 0, 0,    -ps/rho]
    !   [0, 0, 0, 0,          0]
    !   [0, 0, 0, 0,          0]
    !   [0, 0, 0, 0,          0]

    call Dx1(a_x,q(:,:,:,5),params%geom%dx1)
    rhs(:,:,:,1) = rhs(:,:,:,1)                             &
         + qs(:,:,:,5)*u_1*INV_rho                     * a_x

    rhs(:,:,:,2) = rhs(:,:,:,2)                             &
         - qs(:,:,:,5)*INV_rho                         * a_x

!!! D_x2 part
    !-(A^T)^(-1)*B2^T = X2
    !   [0, -u_1*u_2,   -u_2**2, -u_2*u_3, -gamma*p*u_2/(rho*(gamma - 1))]
    !   [0,      u_2,         0,        0,                              0]
    !   [1,      u_1,     2*u_2,      u_3,      gamma*p/(rho*(gamma - 1))]
    !   [0,        0,         0,      u_2,                              0]
    !   [0,        0, gamma - 1,        0,                      gamma*u_2]

    call Dx2(a_x,qs(:,:,:,1),params%geom%dx2)
    call Dx2(b_x,qs(:,:,:,2),params%geom%dx2)
    call Dx2(c_x,qs(:,:,:,3),params%geom%dx2)
    call Dx2(d_x,qs(:,:,:,4),params%geom%dx2)
    call Dx2(e_x,qs(:,:,:,5),params%geom%dx2)

    rhs(:,:,:,1) = rhs(:,:,:,1)                             &
         +(u_2*u_1)                                   * b_x &
         +(u_2**2)                                    * c_x &
         +(u_2*u_3)                                   * d_x &
         +gamma_DOT_p_INV_rho_gm1 * u_2               * e_x

    rhs(:,:,:,2) = rhs(:,:,:,2)                             &
         - u_2                                        * b_x

    rhs(:,:,:,3) = rhs(:,:,:,3)                             &
         -                                              a_x &
         -        u_1                                 * b_x &
         -      2*u_2                                 * c_x &
         -        u_3                                 * d_x &
         - gamma_DOT_p_INV_rho_gm1                    * e_x


    rhs(:,:,:,4) = rhs(:,:,:,4)                             &
         - u_2                                        * d_x

    rhs(:,:,:,5) = rhs(:,:,:,5)                             &
         - gm1                                        * c_x &
         - params%material%gamma*u_2                  * e_x

    !-(A^T)^(-1) * (C1^T * qs)_x2
    !
    !   [1, -u_1/rho, -u_2/rho, -u_3/rho,         0]
    !   [0,    1/rho,        0,        0,         0]
    !   [0,        0,    1/rho,        0,         0]
    !   [0,        0,        0,    1/rho,         0]
    !   [0,        0,        0,        0, gamma - 1]
    !
    !   [0, 0, 0, 0,    0]
    !   [0, 0, 0, 0,    0]
    !   [0, 0, 0, 0,    0]
    !   [0, 0, 0, 0,    0]
    !   [0, 0, 0, 0, -u_2]

    call Dx2(a_x,-u_2*qs(:,:,:,5),params%geom%dx2)
    rhs(:,:,:,5) = rhs(:,:,:,5)                             &
         -gm1*a_x

    !+(A^T)^(-1) * C2~ * c_x2
    !   [0, 0, 0, 0, ps*u_2/rho]
    !   [0, 0, 0, 0,          0]
    !   [0, 0, 0, 0,    -ps/rho]
    !   [0, 0, 0, 0,          0]
    !   [0, 0, 0, 0,          0]

    call Dx2(a_x,q(:,:,:,5),params%geom%dx2)
    rhs(:,:,:,1) = rhs(:,:,:,1)                             &
         + qs(:,:,:,5)*u_2*INV_rho                     * a_x

    rhs(:,:,:,3) = rhs(:,:,:,3)                             &
         - qs(:,:,:,5)*INV_rho                         * a_x

!!! D_x3 part
    !-(A^T)^(-1)*B3^T = X3
    !   [0, -u_1*u_3, -u_2*u_3,   -u_3**2, -gamma*p*u_3/(rho*(gamma - 1))]
    !   [0,      u_3,        0,         0,                              0]
    !   [0,        0,      u_3,         0,                              0]
    !   [1,      u_1,      u_2,     2*u_3,      gamma*p/(rho*(gamma - 1))]
    !   [0,        0,        0, gamma - 1,                      gamma*u_3]

    call Dx3(a_x,qs(:,:,:,1),params%geom%dx3)
    call Dx3(b_x,qs(:,:,:,2),params%geom%dx3)
    call Dx3(c_x,qs(:,:,:,3),params%geom%dx3)
    call Dx3(d_x,qs(:,:,:,4),params%geom%dx3)
    call Dx3(e_x,qs(:,:,:,5),params%geom%dx3)

    rhs(:,:,:,1) = rhs(:,:,:,1)                             &
         +(u_3*u_1)                                   * b_x &
         +(u_3*u_2)                                   * c_x &
         +(u_3**2)                                    * d_x &
         +gamma_DOT_p_INV_rho_gm1 * u_3               * e_x

    rhs(:,:,:,2) = rhs(:,:,:,2)                             &
         - u_3                                        * b_x

    rhs(:,:,:,3) = rhs(:,:,:,3)                             &
         - u_3                                        * c_x

    rhs(:,:,:,4) = rhs(:,:,:,4)                             &
         -                                              a_x &
         -        u_1                                 * b_x &
         -        u_2                                 * c_x &
         -      2*u_3                                 * d_x &
         - gamma_DOT_p_INV_rho_gm1                    * e_x

    rhs(:,:,:,5) = rhs(:,:,:,5)                             &
         - gm1                                        * d_x &
         - params%material%gamma*u_3                  * e_x

    !-(A^T)^(-1) * (C1^T * qs)_x3
    !
    !   [1, -u_1/rho, -u_2/rho, -u_3/rho,         0]
    !   [0,    1/rho,        0,        0,         0]
    !   [0,        0,    1/rho,        0,         0]
    !   [0,        0,        0,    1/rho,         0]
    !   [0,        0,        0,        0, gamma - 1]
    !
    !   [0, 0, 0, 0,    0]
    !   [0, 0, 0, 0,    0]
    !   [0, 0, 0, 0,    0]
    !   [0, 0, 0, 0,    0]
    !   [0, 0, 0, 0, -u_3]

    call Dx3(a_x,-u_3*qs(:,:,:,5),params%geom%dx3)
    rhs(:,:,:,5) = rhs(:,:,:,5)                             &
         -gm1*a_x

    !+(A^T)^(-1) * C3~ * c_x3
    !   [0, 0, 0, 0, ps*u_3/rho]
    !   [0, 0, 0, 0,          0]
    !   [0, 0, 0, 0,          0]
    !   [0, 0, 0, 0,    -ps/rho]
    !   [0, 0, 0, 0,          0]

    call Dx3(a_x,q(:,:,:,5),params%geom%dx3)
    rhs(:,:,:,1) = rhs(:,:,:,1)                             &
         + qs(:,:,:,5)*u_3*INV_rho                     * a_x

    rhs(:,:,:,4) = rhs(:,:,:,4)                             &
         - qs(:,:,:,5)*INV_rho                         * a_x

!!! adjoint force
    rhs = rhs + g

!!! adjoint sponge
    call sponge_adjoint(rhs,qs,q)

!!! boundary conditions 2/2: rhss
    call set_boundary_condition_rhss(rhs,qs)

    deallocate(qs)

    deallocate(u_1)
    deallocate(u_2)
    deallocate(u_3)

    deallocate(gamma_DOT_p_INV_rho_gm1)
    deallocate(                INV_rho)

    deallocate(a_x)
    deallocate(b_x)
    deallocate(c_x)
    deallocate(d_x)
    deallocate(e_x)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine right_hand_side_adjoint_3d
!!!=================================================================================================

end module ns_rho_rhou_p
