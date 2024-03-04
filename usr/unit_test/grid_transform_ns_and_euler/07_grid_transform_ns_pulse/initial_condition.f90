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
    real(kind=rk)                                                   :: x10, x20, x30                     
    real(kind=rk)                                                   :: x11, x21, x31
    real(kind=rk)                                                   :: L1 , L2 , L3
    real(kind=rk)                                                   :: sigma, fwhm, gamma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_parameter(rho0,'init.rho')
    call get_parameter(u0  ,'init.u1' )
    call get_parameter(v0  ,'init.u2' )
    call get_parameter(w0  ,'init.u3' )
    call get_parameter(p0  ,'init.p'  )

    call get_parameter(x10,'geom.x10',default = 0.0_rk)    
    call get_parameter(x20,'geom.x20',default = 0.0_rk)
    call get_parameter(x30,'geom.x30',default = 0.0_rk)

    call get_parameter(x11,'geom.x11',default = 1.0_rk)
    call get_parameter(x21,'geom.x21',default = 1.0_rk)
    call get_parameter(x31,'geom.x31',default = 1.0_rk)

    call get_parameter(gamma,'material.gamma', default = 1.4_rk)

    L1 = (x11 - x10)/2
    L2 = (x21 - x20)/2
    L3 = (x31 - x30)/2

    call get_parameter(center(1)  ,'init.center', 1 , default = L1    )
    call get_parameter(center(2)  ,'init.center', 2 , default = L2    )
    call get_parameter(center(3)  ,'init.center', 3 , default = L3    )
    call get_parameter(fwhm       ,'init.fwhm'      , default = 2.0_rk)

    sigma = fwhm/2.3548_rk*params%geom%dx1

    q0(:,:,:,2) = q0(:,:,:,1)*u0
    q0(:,:,:,3) = q0(:,:,:,1)*v0
    q0(:,:,:,4) = q0(:,:,:,1)*w0
    q0(:,:,:,5) = p0 + 0.1*exp(-((X(:,:,:,1) - center(1))**2 + (X(:,:,:,2) - center(2))**2 + (X(:,:,:,3) - center(3))**2)/(2.0_rk*sigma**2))
    q0(:,:,:,1) = rho0*(q0(:,:,:,5)/p0)**(1.0/gamma)    
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
