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
    real(kind=rk)                                                   :: rho0, u0, v0, w0, p0, dp
    real(kind=rk)                                                   :: gamma
    real(kind=rk)                                                   :: x11,fwidth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_parameter(rho0,'init.rho')
    call get_parameter(u0  ,'init.u1' )
    call get_parameter(v0  ,'init.u2' )
    call get_parameter(w0  ,'init.u3' )

    call get_parameter(p0  ,'init.p'  )
    call get_parameter(dp  ,'init.dp' )
    call get_parameter(x11 ,'geom.x11')
    call get_parameter(fwidth ,'init.fwidth_dp')

    call get_parameter(rho0,'init.rho')
    call get_parameter(gamma,'material.gamma', default = 1.4_rk)

    q0(:,:,:,2) =  u0 		! defintion u
    q0(:,:,:,3) =  0.0_rk	! defintion v	
    q0(:,:,:,4) =  0.0_rk	! defintion w
    q0(:,:,:,5) =  p0 + (dp*(0.5_rk*(tanh(1.0_rk/fwidth*(X(:,:,:,1) - (x11/2)))   - tanh(1.0_rk/fwidth*(X(:,:,:,1) - 1e4)))) -dp/2)	! defintion p

    q0(:,:,:,1) =  rho0*(q0(:,:,:,5)/p0)**(1.0/gamma)
    q0(:,:,:,2) =  q0(:,:,:,2) * q0(:,:,:,1) 			! rho * u	
    q0(:,:,:,3) =  q0(:,:,:,3) * q0(:,:,:,1) 			! rho * v
    q0(:,:,:,4) =  q0(:,:,:,4) * q0(:,:,:,1) 			! rho * w
    

  end subroutine init_direct
!!!=================================================================================================


!!!=================================================================================================
  subroutine init_adjoint(qs0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)                :: qs0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(3)                                     :: center
    real(kind=rk)                                                   :: rho0, u0, v0, w0, p0 
    real(kind=rk)                                                   :: sigma, fwhm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_parameter(center(1)  ,'init.center', 1 )
    call get_parameter(center(2)  ,'init.center', 2 )
    call get_parameter(center(3)  ,'init.center', 3 )
    call get_parameter(fwhm       ,'init.fwhm'      )

    sigma = fwhm/2.3548_rk*params%geom%dx1

    qs0 = 0.0_rk

    qs0(:,:,:,5) = 0.0_rk*exp(-((X(:,:,:,1) - center(1))**2 + (X(:,:,:,2) - center(2))**2 + (X(:,:,:,3) - center(3))**2)/(2.0_rk*sigma**2))

  end subroutine init_adjoint
!!!=================================================================================================


end module initial_condition
