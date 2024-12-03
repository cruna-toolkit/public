module grid_metric

  use data_geom
  use discretisation_x
  use helper
  use parameter
  use io

  private

  public :: calc_metric

contains

!!!=================================================================================================
  subroutine calc_metric(X_calc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , dimension(:,:,:,:), intent(in)  :: X_calc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , dimension(:,:,:,:), allocatable :: X_temp
    real(kind=rk)   , dimension(:,:,:)  , allocatable :: dx1dxi1,dx1dxi2,dx1dxi3
    real(kind=rk)   , dimension(:,:,:)  , allocatable :: dx2dxi1,dx2dxi2,dx2dxi3
    real(kind=rk)   , dimension(:,:,:)  , allocatable :: dx3dxi1,dx3dxi2,dx3dxi3
    real(kind=rk)   , dimension(:,:,:)  , allocatable :: temp1, temp2
    real(kind=rk)   , dimension(3)                    :: force_periodic = 0.0_rk
    real(kind=rk)                                     :: c1, c2, c3                   ! correction used to make pseudo periodic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! grid factors, using  rewriting of
    ! Miguel R Visbal and Datta V Gaitonde. On the use of higher-order
    ! finite-difference schemes on curvilinear and deforming meshes. Journal
    ! of Computational Physics, 181(1):155 – 185, 2002.

    ! remove periodicity if artifical periodic to avoid jump at end
    ! ATTENTION: does not work for real periodic grids (maybe --> ask JR)

    ! general remark:
    ! two nasty  problems arise
    ! for grids which are asumed to be continued in a periodic manner the grid derivatives jump at the end.
    ! The force periodic switch compensates it, see second homework in CFD2
    ! for full three D transformed grids a rewriting of the metric factors is needed, to hava a zero gradient of a constant function
    ! see CFD2 script or Visbal Gaitonde

    if (trim(diff_i1_type).eq.'periodic') then
       force_periodic(1) = 1.0_rk
    end if

    if (trim(diff_i2_type).eq.'periodic') then
       force_periodic(2) = 1.0_rk
    end if

    if (trim(diff_i3_type).eq.'periodic') then
       force_periodic(3) = 1.0_rk
    end if

    ! for readability define correction factors
    c1 = force_periodic(1)
    c2 = force_periodic(2)
    c3 = force_periodic(3)

    ! save c1,c2,c3 in params.geom
    call set_parameter(c1,'geom.c1')
    call set_parameter(c2,'geom.c2')
    call set_parameter(c3,'geom.c3')

    call store(X_calc,'geometry_calc_space')
    ! define corrected grid --> X_temp
    allocate(X_temp(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))
    X_temp(:,:,:,1)  = X(:,:,:,1) - c1 * X_calc(:,:,:,1)                         ! this is needed as otherwise a jump happens at the grid boundary if assuming periodicity
    X_temp(:,:,:,2)  = X(:,:,:,2) - c2 * X_calc(:,:,:,2)
    X_temp(:,:,:,3)  = X(:,:,:,3) - c3 * X_calc(:,:,:,3)

    allocate(dx1dxi1(params%geom%n1b,params%geom%n2b,params%geom%n3b))       ! all  derivatves of the physical grid
    allocate(dx1dxi2(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(dx1dxi3(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(dx2dxi1(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(dx2dxi2(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(dx2dxi3(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(dx3dxi1(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(dx3dxi2(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(dx3dxi3(params%geom%n1b,params%geom%n2b,params%geom%n3b))

    allocate(temp1(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(temp2(params%geom%n1b,params%geom%n2b,params%geom%n3b))

    ! dxi = dx, since the grid created by cartesian_2_cartesian_strech_shift is used as the calculation space
    call Dx1(dx1dxi1,X_temp(:,:,:,1), params%geom%dx1)     ! first local base vector, only periodic part if force_periodic = .true.
    call Dx1(dx2dxi1,X_temp(:,:,:,2), params%geom%dx1)
    call Dx1(dx3dxi1,X_temp(:,:,:,3), params%geom%dx1)

    call Dx2(dx1dxi2,X_temp(:,:,:,1), params%geom%dx2)     ! second local base vector, only periodic part if force_periodic = .true.
    call Dx2(dx2dxi2,X_temp(:,:,:,2), params%geom%dx2)
    call Dx2(dx3dxi2,X_temp(:,:,:,3), params%geom%dx2)

    call Dx3(dx1dxi3,X_temp(:,:,:,1), params%geom%dx3)     ! third  local base vector, only periodic part if force_periodic = .true.
    call Dx3(dx2dxi3,X_temp(:,:,:,2), params%geom%dx3)
    call Dx3(dx3dxi3,X_temp(:,:,:,3), params%geom%dx3)

    X_metric(:,:,:,1,1) = dx1dxi1               ! dx1dxi1 = dxdxi
    X_metric(:,:,:,1,2) = dx1dxi2               ! dx1dxi2 = dxdeta
    X_metric(:,:,:,1,3) = dx1dxi3               ! dx1dxi3 = dxdzeta
    X_metric(:,:,:,2,1) = dx2dxi1               ! dx2dxi1 = dydxi
    X_metric(:,:,:,2,2) = dx2dxi2               ! dx2dxi2 = dydeta
    X_metric(:,:,:,2,3) = dx2dxi3               ! dx2dxi3 = dydzeta
    X_metric(:,:,:,3,1) = dx3dxi1               ! dx3dxi1 = dzdxi
    X_metric(:,:,:,3,2) = dx3dxi2               ! dx3dxi2 = dzdeta
    X_metric(:,:,:,3,3) = dx3dxi3               ! dx3dxi3 = dzdzeta

    ! for a central difference scheme  the error of the metric factors can be reduce by using the following computation
    ! e.g. y_eta*z_zeta - y_zeta*z_eta = (y_eta*z)_zeta - (y_zeta*z)_eta = y_eta*z_zeta - y_zeta*z_eta + (y_eta_zeta*z) - (y_zeta_eta*z)
    ! finally in case of periodicity these terms are adjusted with the correction factors c1, c2, c3
    ! (y_eta + c2)(z_zeta + c3) - y_zeta*z_eta = y_eta*z_zeta - y_zeta*z_eta + c2*z_zeta + c3*y_eta + c2*c3

    !params.geom.m11      = dx2dxi2.*dx3dxi3  -  dx3dxi2.*dx2dxi3
    !params.geom.m11      = dx2dxi2.*dx3dxi3  -  dx3dxi2.*dx2dxi3 +c2*dx3dxi3+c3*dx2dxi2 +c2*c3
    !(y_eta*z)_zeta - (y_zeta*z)_eta + c2*dx3dxi3 + c3*dx2dxi2 + c2*c3
    call Dx3(temp1, dx2dxi2 * X_temp(:,:,:,3), params%geom%dx3)    ! D_zeta(X_temp(:,:,:,3).*dx2dxi2)
    call Dx2(temp2, dx2dxi3 * X_temp(:,:,:,3), params%geom%dx2)    ! D_eta(X_temp(:,:,:,3).*dx2dxi3)
    X_metric(:,:,:,1,1)   =  temp1  -  temp2  + c2*dx3dxi3 + c3*dx2dxi2 + c2*c3


    !params.geom.m12      = dx3dxi2.*dx1dxi3  -  dx1dxi2.*dx3dxi3
    !params.geom.m12      = dx3dxi2.*dx1dxi3  -  dx1dxi2.*dx3dxi3  - c3*dx1dxi2
    !(x_zeta*z)_eta-(x_eta*z)_zeta - c3*dx1dxi2
    call Dx2(temp1, dx1dxi3 * X_temp(:,:,:,3), params%geom%dx2)    ! D_eta(X_temp(:,:,:,3).*dx1dxi3)
    call Dx3(temp2, dx1dxi2 * X_temp(:,:,:,3), params%geom%dx3)    ! D_zeta(X_temp(:,:,:,3).*dx1dxi2)
    X_metric(:,:,:,1,2)   =  temp1  -  temp2  - c3*dx1dxi2


    !params.geom.m13      = dx1dxi2.*dx2dxi3  -  dx2dxi2.*dx1dxi3
    !params.geom.m13      = dx1dxi2.*dx2dxi3  -  dx2dxi2.*dx1dxi3  - c2*dx1dxi3
    !(x_eta*y)_zeta - (x_zeta*y)_eta - c2*dx1dxi3
    call Dx3(temp1, dx1dxi2 * X_temp(:,:,:,2), params%geom%dx3)    ! D_zeta(X_temp(:,:,:,2).*dx1dxi2)
    call Dx2(temp2, dx1dxi3 * X_temp(:,:,:,2), params%geom%dx2)    ! D_eta(X_temp(:,:,:,2)*dx1dxi3)
    X_metric(:,:,:,1,3)   =  temp1 -  temp2   - c2*dx1dxi3

    !params.geom.m21      = dx2dxi3.*dx3dxi1  -  dx3dxi3.*dx2dxi1
    !params.geom.m21      = dx2dxi3.*dx3dxi1  -  dx3dxi3.*dx2dxi1 -c3*dx2dxi1
    !(y_zeta*z)_xi - (y_xi*z)_zeta - c3*dx2dxi1
    call Dx1(temp1, dx2dxi3 * X_temp(:,:,:,3), params%geom%dx1)    ! D_xi(X_temp(:,:,:,2).*dx3dxi1)
    call Dx3(temp2, dx2dxi1 * X_temp(:,:,:,3), params%geom%dx3)    ! D_zeta(dx3dxi3.*X_temp(:,:,:,2))
    X_metric(:,:,:,2,1)   = temp1  -  temp2   - c3*dx2dxi1


    !params.geom.m22      = dx3dxi3.*dx1dxi1  -  dx1dxi3.*dx3dxi1
    !params.geom.m22      = dx3dxi3.*dx1dxi1  + dx3dxi3 *c1  + c3 dx1dxi1 + c1 c3 -  dx1dxi3.*dx3dxi1
    !(x_xi*z)_zeta - (x_zeta*z)_xi + c1*dx3dxi3   + c3*dx1dxi1 + c1*c3
    call Dx3(temp1, dx1dxi1 * X_temp(:,:,:,3), params%geom%dx3)    ! D_zeta(X_temp(:,:,:,3).*dx1dxi1)
    call Dx1(temp2, dx1dxi3 * X_temp(:,:,:,3), params%geom%dx1)    ! D_xi(dx1dxi3.*X_temp(:,:,:,3))
    X_metric(:,:,:,2,2)   = temp1  -  temp2   + c1*dx3dxi3  + c3*dx1dxi1 + c1*c3


    !params.geom.m23      = dx1dxi3.*dx2dxi1  -  dx2dxi3.*dx1dxi1
    !params.geom.m23      = dx1dxi3.*dx2dxi1  -  dx2dxi3.*dx1dxi1 - c1*dx2dxi3
    !(x_zeta*y)_xi - (x_xi*y)_zeta - c1*dx2dxi3
    call Dx1(temp1, dx1dxi3 * X_temp(:,:,:,2), params%geom%dx1)   ! D_xi(dx1dxi3.*X_temp(:,:,:,2))
    call Dx3(temp2, dx1dxi1 * X_temp(:,:,:,2), params%geom%dx3)   ! D_zeta(X_temp(:,:,:,2).*dx1dxi1)
    X_metric(:,:,:,2,3)   = temp1  -  temp2  - c1*dx2dxi3

    !params.geom.m31      = dx2dxi1.*dx3dxi2  -  dx2dxi2.*dx3dxi1
    !params.geom.m31      = dx2dxi1.*dx3dxi2  -  dx2dxi2.*dx3dxi1   -c2*dx3dxi1
    !(y_xi*z)_eta - (y_eta*z)_xi
    call Dx2(temp1, X_temp(:,:,:,3) * dx2dxi1, params%geom%dx2)   ! D_eta(X_temp(:,:,:,3).*dx2dxi1)
    call Dx1(temp2, X_temp(:,:,:,3) * dx2dxi2, params%geom%dx1)   ! D_xi(X_temp(:,:,:,3).*dx2dxi2)
    X_metric(:,:,:,3,1)   = temp1  -  temp2  - c2*dx3dxi1

    !params.geom.m32      = dx3dxi1.*dx1dxi2  -  dx3dxi2.*dx1dxi1
    !params.geom.m32      = dx3dxi1.*dx1dxi2  -  dx3dxi2.*dx1dxi1 -c1*dx3dxi2
    !(x_eta*z)_xi - (x_xi*z)_eta
    call Dx1(temp1, dx1dxi2 * X_temp(:,:,:,3), params%geom%dx1)   ! D_xi(X_temp(:,:,:,3).*dx1dxi2)
    call Dx2(temp2, dx1dxi1 * X_temp(:,:,:,3), params%geom%dx2)   ! D_eta(X_temp(:,:,:,3).*dx1dxi1)
    X_metric(:,:,:,3,2)   = temp1  -  temp2   - c1*dx3dxi2

    !params.geom.m33      = dx1dxi1.*dx2dxi2  -  dx1dxi2.*dx2dxi1
    !params.geom.m33      = dx1dxi1.*dx2dxi2  -  dx1dxi2.*dx2dxi1 + c1*dx2dxi2 + dx1dxi1*c2 + c1*c2
    !(x_xi*y)_eta - (x_eta*y)_xi
    call Dx2(temp1, dx1dxi1 * X_temp(:,:,:,2), params%geom%dx2)   ! D_eta(dx1dxi1.*X_temp(:,:,:,2))
    call Dx1(temp2, dx1dxi2 * X_temp(:,:,:,2), params%geom%dx1)   ! D_xi(dx1dxi2.*X_temp(:,:,:,2))
    X_metric(:,:,:,3,3)   = temp1   -  temp2  + c1*dx2dxi2 + dx1dxi1*c2 + c1*c2

    !X_metric is the metric without sclaing with the Jacobian
    ! Definition of Jacobian follows Visbal et al.
    X_jacobian = X_metric(:,:,:,3,1)*dx1dxi3 + &
         + X_metric(:,:,:,3,2)*dx2dxi3 + &
         + X_metric(:,:,:,3,3)*(dx3dxi3+c3)
    X_jacobian = 1.0_rk/X_jacobian


  end subroutine calc_metric
!!!=================================================================================================

end module grid_metric
