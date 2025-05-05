module shape_functions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Collection of various phi shapes                                                                  !
! The length of the region where phi>0 is defined by sponge.length of parameter.dat.                !
! Max(phi) is 1 unless defined by sponge.amplitude in parameter.dat.                                !
! sponge.amplitude = 1.0 (recommended for nonreflecting zonal characteristic BC acting on q or qs)  !
! sponge.amplitude = 1.0e3 (recommended for classic sponge acting on rhs or rhss)                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use data_geom , only: X
  use io , only: store
  use parameter , only: rk, params, get_parameter

  private

  public :: init_quadratic_shape_fct, init_cos_shape_fct

  private :: quad_shape_fct, cos_shape_fct
  private :: init_general_shape_fct

contains

!!!=================================================================================================
  subroutine init_quadratic_shape_fct(phi_shape)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), allocatable, intent(out)    :: phi_shape
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call init_general_shape_fct(phi_shape, quad_shape_fct)

  end subroutine init_quadratic_shape_fct
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_cos_shape_fct(phi_shape)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), allocatable, intent(out)    :: phi_shape
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call init_general_shape_fct(phi_shape, cos_shape_fct)

  end subroutine init_cos_shape_fct
!!!=================================================================================================

!!!=================================================================================================
  function quad_shape_fct(x,x0,len0) result(y)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(in)                 :: x
    real(kind=rk), intent(in)                                   :: x0, len0
    real(kind=rk), dimension(size(x,1),size(x,2),size(x,3))     :: y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    y = 1.0_rk - ((x-x0)/len0)**2
  end function quad_shape_fct
!!!=================================================================================================

!!!=================================================================================================
  function cos_shape_fct(x,x0,len0) result(y)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(in)                 :: x
    real(kind=rk), intent(in)                                   :: x0, len0
    real(kind=rk), dimension(size(x,1),size(x,2),size(x,3))     :: y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), parameter                                    :: pi=4._rk*atan(1.0_rk)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    y = 0.5_rk + 0.5_rk*cos((x-x0)*(pi/len0))
  end function cos_shape_fct
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_general_shape_fct(phi_shape, fct_pointer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), allocatable, intent(out)    :: phi_shape
    interface
        function fct_pointer(x,x0,len0) result(y)
            import                                                  :: rk
            real(kind=rk), dimension(:,:,:), intent(in)             :: x
            real(kind=rk), intent(in)                               :: x0, len0
            real(kind=rk), dimension(size(x,1),size(x,2),size(x,3)) :: y
        end function fct_pointer
    end interface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), allocatable                 :: sigma
    real(kind=rk)                                                :: x10,x11,x20,x21,x30,x31
    real(kind=rk)                                                :: x0,len1,len2,len3,amplitude
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(phi_shape(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    allocate(sigma(params%geom%n1b,params%geom%n2b,params%geom%n3b))

    ! get block dimensions from parameter.dat
    call get_parameter(x10,'geom.x10',default = 0.0_rk)
    call get_parameter(x20,'geom.x20',default = 0.0_rk)
    call get_parameter(x30,'geom.x30',default = 0.0_rk)
    call get_parameter(x11,'geom.x11',default = 1.0_rk)
    call get_parameter(x21,'geom.x21',default = 1.0_rk)
    call get_parameter(x31,'geom.x31',default = 1.0_rk)

    ! get sponge dimensions from parameter.dat
    call get_parameter(len1,'sponge.length',default=0.3_rk,subindex=1) ! for 1-4kHz acoustic waves
    call get_parameter(len2,'sponge.length',default=0.3_rk,subindex=2) ! for 1-4kHz acoustic waves
    call get_parameter(len3,'sponge.length',default=0.3_rk,subindex=3) ! for 1-4kHz acoustic waves

    call get_parameter(amplitude,'sponge.amplitude',default=1.0_rk) ! default

    ! xm
    x0 = x10 + len1
    sigma = fct_pointer(X(:,:,:,1), x0, len1)
    where(X(:,:,:,1) > x0 ) sigma = 1.0_rk
    phi_shape = sigma

    ! xp
    x0  = x11 - len1
    sigma = fct_pointer(X(:,:,:,1), x0, len1)
    where(X(:,:,:,1) < x0 ) sigma = 1.0_rk
    phi_shape = phi_shape * sigma

    ! ym
    x0 = x20 + len2
    sigma = fct_pointer(X(:,:,:,2), x0, len2)
    where(X(:,:,:,2) > x0 ) sigma = 1.0_rk
    phi_shape = phi_shape * sigma

    ! yp
    x0  = x21 - len2
    sigma = fct_pointer(X(:,:,:,2), x0, len2)
    where(X(:,:,:,2) < x0 ) sigma = 1.0_rk
    phi_shape = phi_shape * sigma

    ! zm
    x0 = x30 + len3
    sigma = fct_pointer(X(:,:,:,3), x0, len3)
    where(X(:,:,:,3) > x0 ) sigma = 1.0_rk
    phi_shape = phi_shape * sigma

    ! zp
    x0  = x31 - len3
    sigma = fct_pointer(X(:,:,:,3), x0, len3)
    where(X(:,:,:,3) < x0 ) sigma = 1.0_rk
    phi_shape = phi_shape * sigma
    deallocate(sigma)

    ! invert phi_shape to the range [0,1]. 1 is where the boundary is located. 0 is inner point.
    phi_shape = 1.0_rk - phi_shape

    ! pre-multiply amplitude (equal for all quantities)
    phi_shape = amplitude*phi_shape

    ! store sponge
    if (params%io%verbosity.ge.2) call store(phi_shape,'shape_fct')

  end subroutine init_general_shape_fct
!!!=================================================================================================

end module shape_functions