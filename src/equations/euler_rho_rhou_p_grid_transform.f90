module euler_rho_rhou_p_grid_transform

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! information/remarks:                                                                             !
! - q_in, qs_in in subroutine calls, this is needed due to intent(in) of q/qs in time stepper,     !
!   thus no modification inside of the rhs is permitted, however no intent works fine, even        !
!   nobody understands why ... (2019-01-02, ML)                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use boundary_conditions
  use data_geom
  use data_reference , only: qref
  use discretisation_x
  use euler_rho_rhou_p_grid_transform_trafos
  use force
  use io
  use parameter
  use sponge_layer

  private

  public right_hand_side_direct_3d
  public right_hand_side_adjoint_3d

  public give_T, give_c, give_r, give_u, give_p

  character(len=max_length_parameter), parameter, public   :: equations_name = "euler_rho_rhou_p_grid_transform"
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
    real(kind=rk), dimension(:,:,:)  , allocatable   :: a_x, a_y, a_z
    real(kind=rk), dimension(:,:,:)  , allocatable   :: p_x, p_y, p_z
    integer                                          :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable   :: u_proj
    real(kind=rk)                                    :: gamma_term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    allocate(  q(params%geom%n1b,params%geom%n2b,params%geom%n3b,5))
    allocate(  u(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))
    allocate(  u_proj(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))

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

!!! unpack variables
    do i = 1,params%geom%n3b
       u(:,:,i,1) = q(:,:,i,2)/q(:,:,i,1)
       u(:,:,i,2) = q(:,:,i,3)/q(:,:,i,1)
       u(:,:,i,3) = q(:,:,i,4)/q(:,:,i,1)
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
    call Dx1(a_x,q(:,:,:,1)*u_proj(:,:,:,1),params%geom%dx1)
    call Dx2(a_y,q(:,:,:,1)*u_proj(:,:,:,2),params%geom%dx2)
    call Dx3(a_z,q(:,:,:,1)*u_proj(:,:,:,3),params%geom%dx3)

    rhs(:,:,:,1) = -a_x - a_y - a_z
    rhs(:,:,:,1) = X_jacobian * rhs(:,:,:,1)

!!! momentum equation (1) ------------------------------------------------------------------------
    call Dx1(a_x,q(:,:,:,1)*u(:,:,:,1)*u_proj(:,:,:,1),params%geom%dx1)
    call Dx2(a_y,q(:,:,:,1)*u(:,:,:,1)*u_proj(:,:,:,2),params%geom%dx2)
    call Dx3(a_z,q(:,:,:,1)*u(:,:,:,1)*u_proj(:,:,:,3),params%geom%dx3)

    rhs(:,:,:,2) = -a_x - a_y - a_z - p_x
    rhs(:,:,:,2) = X_jacobian * rhs(:,:,:,2)

!!! momentum equation (2) ------------------------------------------------------------------------
    call Dx1(a_x,q(:,:,:,1)*u(:,:,:,2)*u_proj(:,:,:,1),params%geom%dx1)
    call Dx2(a_y,q(:,:,:,1)*u(:,:,:,2)*u_proj(:,:,:,2),params%geom%dx2)
    call Dx3(a_z,q(:,:,:,1)*u(:,:,:,2)*u_proj(:,:,:,3),params%geom%dx3)

    rhs(:,:,:,3) = -a_x - a_y - a_z - p_y
    rhs(:,:,:,3) = X_jacobian * rhs(:,:,:,3)

!!! momentum equation (3) ------------------------------------------------------------------------
    call Dx1(a_x,q(:,:,:,1)*u(:,:,:,3)*u_proj(:,:,:,1),params%geom%dx1)
    call Dx2(a_y,q(:,:,:,1)*u(:,:,:,3)*u_proj(:,:,:,2),params%geom%dx2)
    call Dx3(a_z,q(:,:,:,1)*u(:,:,:,3)*u_proj(:,:,:,3),params%geom%dx3)

    rhs(:,:,:,4) = -a_x - a_y - a_z - p_z
    rhs(:,:,:,4) = X_jacobian * rhs(:,:,:,4)

!!! energy (pressure) equation -------------------------------------------------------------------
    call Dx1(a_x,q(:,:,:,5)*u_proj(:,:,:,1),params%geom%dx1)
    call Dx2(a_y,q(:,:,:,5)*u_proj(:,:,:,2),params%geom%dx2)
    call Dx3(a_z,q(:,:,:,5)*u_proj(:,:,:,3),params%geom%dx3)

    gamma_term = params%material%gamma / (params%material%gamma - 1.0_rk)  ! gamma /(gamma -1)
    rhs(:,:,:,5) = gamma_term*(-a_x - a_y - a_z)

    rhs(:,:,:,5) = rhs(:,:,:,5) + u(:,:,:,1)*p_x(:,:,:) + u(:,:,:,2)*p_y(:,:,:) + u(:,:,:,3)*p_z(:,:,:)

    gamma_term   = (params%material%gamma - 1.0_rk)          ! (gamma -1)
    rhs(:,:,:,5) = X_jacobian * gamma_term * rhs(:,:,:,5)

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

end module euler_rho_rhou_p_grid_transform