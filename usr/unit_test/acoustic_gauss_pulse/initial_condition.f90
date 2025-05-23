module initial_condition

  use data_geom
  use parameter

  private

  public :: init_direct
  public :: init_adjoint

contains

!!!=================================================================================================
  subroutine init_direct(q0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)                :: q0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(3)                                     :: center
    real(kind=rk)                                                   :: rho0, u0, v0, w0, p0 
    real(kind=rk)                                                   :: radius, amplitude
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_parameter(rho0       ,'init.rho')
    call get_parameter(u0         ,'init.u1' )
    call get_parameter(v0         ,'init.u2' )
    call get_parameter(w0         ,'init.u3' )
    call get_parameter(p0         ,'init.p'  )
    call get_parameter(center(1)  ,'init.center', 1 )
    call get_parameter(center(2)  ,'init.center', 2 )
    call get_parameter(center(3)  ,'init.center', 3 )
    call get_parameter(radius     ,'init.radius')
    call get_parameter(amplitude  ,'init.amplitude')

    !print*,"xmin xmax",minval(X(:,:,:,1)),maxval(X(:,:,:,1)),params%parallelism%world_image
    !print*,"dx dy dz",params%geom%dx1,params%geom%dx2,params%geom%dx3,params%parallelism%world_image

    q0(:,:,:,2) = u0
    q0(:,:,:,3) = v0
    q0(:,:,:,4) = w0
    q0(:,:,:,5) = p0 + amplitude                 &
         *exp(-(center(1)-X(:,:,:,1))**2/radius) &
         *exp(-(center(2)-X(:,:,:,2))**2/radius) &
         *exp(-(center(3)-X(:,:,:,3))**2/radius)
    q0(:,:,:,1) = rho0 * (q0(:,:,:,5)/p0)**(1/params%material%gamma)

  end subroutine init_direct
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_adjoint(qs0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)                :: qs0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(3)                                     :: center
    real(kind=rk)                                                   :: radius
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_parameter(center(1)  ,'init.center', 1 )
    call get_parameter(center(2)  ,'init.center', 2 )
    call get_parameter(center(3)  ,'init.center', 3 )
    call get_parameter(radius     ,'init.radius')

    qs0(:,:,:,1) = 0.0_rk
    qs0(:,:,:,2) = 0.0_rk
    qs0(:,:,:,3) = 0.0_rk
    qs0(:,:,:,4) = 0.0_rk
    qs0(:,:,:,5) = 0.1_rk                        &
         *exp(-(center(1)-X(:,:,:,1))**2/radius) &
         *exp(-(center(2)-X(:,:,:,2))**2/radius) &
         *exp(-(center(3)-X(:,:,:,3))**2/radius)
    qs0(:,:,:,1) = (qs0(:,:,:,5)/1._rk)**(1/params%material%gamma)  
    qs0 = 0

  end subroutine init_adjoint
!!!=================================================================================================

end module initial_condition
