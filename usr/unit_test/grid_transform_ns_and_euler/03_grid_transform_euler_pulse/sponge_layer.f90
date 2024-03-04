module sponge_layer

  use data_geom
  use io
  use parameter

  character(len=max_length_parameter), parameter, public :: sponge_direct_name  = 'quadratic sponge layer (case)'
  character(len=max_length_parameter), parameter, public :: sponge_adjoint_name = 'quadratic sponge layer (case)'

  private

  public sponge
  public sponge_adjoint

contains

!!!=================================================================================================
  subroutine sponge(rhs,q,q0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)             :: rhs
    real(kind=rk), dimension(:,:,:,:), intent(in), optional      :: q,q0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:), allocatable                :: X_calc                                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), allocatable, save           :: sponge_function
    real(kind=rk), dimension(:,:,:), allocatable                 :: sigma
    real(kind=rk)                                                :: len,amp
    real(kind=rk)                                                :: x_min,x_max
    real(kind=rk)                                                :: x10,x11,x20,x21,x30,x31
    real(kind=rk), save                                          :: rho0,u10,u20,u30,p0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(.not.allocated(sponge_function)) then
       allocate(sponge_function(params%geom%n1b,params%geom%n2b,params%geom%n3b))
       allocate(       X_calc(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))
       allocate(          sigma(params%geom%n1b,params%geom%n2b,params%geom%n3b))

       sponge_function = 0.0_rk

       call get_parameter(len,'sponge.length')
       call get_parameter(amp,'sponge.amplitude')
       call get_parameter(x10,'geom.x10')
       call get_parameter(x11,'geom.x11')
       call get_parameter(x20,'geom.x20')
       call get_parameter(x21,'geom.x21')
       call get_parameter(x30,'geom.x30')
       call get_parameter(x31,'geom.x31')
       
       call load(X_calc,'geometry_calc_space')
       ! xm
!       do k = 1,params%geom%n3b
!            do j = 1,params%geom%n2b
!                do i = 1,params%geom%n1b
!                x_min = X(1,j,k,1) + len
!                x_max = X(1,j,k,1)
!                sigma(i,j,k) = (((X(i,j,k,1) - x_min)/(x_max - x_min))**2) - 1.0_rk
!                where(X(i,j,k,1) > x_min ) sigma = -1.0_rk
!                sponge_function(i,j,k) = abs(sigma(i,j,k))
!                end do
!            end do
!       end do
       
       ! xm
       x_min = x10 + len
       x_max = x10
       sigma = (((X_calc(:,:,:,1) - x_min)/(x_max - x_min))**2) - 1.0_rk
       where(X_calc(:,:,:,1) > x_min ) sigma = -1.0_rk
       sponge_function = abs(sigma)
       

       ! xp
       
       x_min  = x11 - len
       x_max  = x11   
       sigma = (((X_calc(:,:,:,1) - x_min)/(x_max - x_min))**2) - 1.0_rk
       where(X_calc(:,:,:,1) < x_min ) sigma = -1.0_rk
       sponge_function = sponge_function * abs(sigma)

       ! ym
       x_min = x20 + len
       x_max = x20
       sigma = (((X_calc(:,:,:,2) - x_min)/(x_max - x_min))**2) - 1.0_rk
       where(X_calc(:,:,:,2) > x_min ) sigma = -1.0_rk
       sponge_function = sponge_function * abs(sigma)

       ! yp
       x_min  = x21 - len
       x_max  = x21   
       sigma = (((X_calc(:,:,:,2) - x_min)/(x_max - x_min))**2) - 1.0_rk
       where(X_calc(:,:,:,2) < x_min ) sigma = -1.0_rk
       sponge_function = sponge_function * abs(sigma)

       ! zm
       x_min = x30 + len
       x_max = x30
       sigma = (((X_calc(:,:,:,3) - x_min)/(x_max - x_min))**2) - 1.0_rk
       where(X_calc(:,:,:,3) > x_min ) sigma = -1.0_rk
       sponge_function = sponge_function * abs(sigma)

       ! zp
       x_min  = x31 - len
       x_max  = x31   
       sigma = (((X_calc(:,:,:,3) - x_min)/(x_max - x_min))**2) - 1.0_rk
       where(X_calc(:,:,:,3) < x_min ) sigma = -1.0_rk
       sponge_function = sponge_function * abs(sigma)

       ! modify sigma
       sponge_function = -1.0_rk * (sponge_function - 1.0_rk)

       ! set amplitude (equal for all quantities)
       sponge_function = amp*sponge_function

       ! get reference values
       call get_parameter(rho0,'init.rho')
       call get_parameter(u10 ,'init.u1' )
       call get_parameter(u20 ,'init.u2' )
       call get_parameter(u30 ,'init.u3' )
       call get_parameter(p0  ,'init.p'  )

       ! store sponge
       call store(sponge_function,'sponge_function')

    end if

    rhs(:,:,:,1) = rhs(:,:,:,1) - sponge_function*(q(:,:,:,1) - rho0    )
    rhs(:,:,:,2) = rhs(:,:,:,2) - sponge_function*(q(:,:,:,2) - rho0*u10)
    rhs(:,:,:,3) = rhs(:,:,:,3) - sponge_function*(q(:,:,:,3) - rho0*u20)
    rhs(:,:,:,4) = rhs(:,:,:,4) - sponge_function*(q(:,:,:,4) - rho0*u30)
    rhs(:,:,:,5) = rhs(:,:,:,5) - sponge_function*(q(:,:,:,5) - p0      )

  end subroutine sponge
!!!=================================================================================================

!!!=================================================================================================
  subroutine sponge_adjoint(rhs,qs,qs0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)             :: rhs
    real(kind=rk), dimension(:,:,:,:), intent(in), optional      :: qs,qs0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rhs = rhs

  end subroutine sponge_adjoint
!!!=================================================================================================

end module sponge_layer
