module sponge_layer_core
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Collection of various sponge routines                                                             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use euler_rho_rhou_p_trafos , only: give_r,give_u
  use parameter

  private

  ! legacy sponge routines - quadratic shape fct
  public sponge_rhs ! recommended classic direct sponge
  public sponge_rhss ! recommended classic adjoint sponge
  public sponge_q
  public sponge_qs

  ! empty/dummy routine
  public empty_sponge
  public empty_sponge_inrhs_warning ! for compatibility of old direct sponge calls by rhs routines

  interface sponge_rhss
     module procedure sponge_rhss_correct
     module procedure empty_sponge_rhss_warning
     ! Empty version needed for compatibility with old rhs routines
     ! because their set_boundary_condition_rhss calls (boundary_conditions.f90) lack q argument.
  end interface sponge_rhss

  private sponge_rhss_correct
  private empty_sponge_rhss_warning

contains

!!!=================================================================================================
  subroutine sponge_rhs(rhs,sponge_shape,q,qref)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)             :: rhs
    real(kind=rk), dimension(:,:,:),   intent(in)                :: sponge_shape
    real(kind=rk), dimension(:,:,:,:), intent(in)                :: q
    real(kind=rk), dimension(:,:,:,:), intent(in), optional      :: qref
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                      :: i
    real(kind=rk), save                                          :: rho0=0.0_rk, u10, u20, u30, p0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (present(qref)) then
      do i = 1, size(rhs,4)
        rhs(:,:,:,i) = rhs(:,:,:,i) - sponge_shape*(q(:,:,:,i) - qref(:,:,:,i))
      end do

    else
      if (rho0 == 0.0_rk) then ! get reference values only once
        call get_parameter(rho0,'init.rho')
        call get_parameter(u10 ,'init.u1' )
        call get_parameter(u20 ,'init.u2' )
        call get_parameter(u30 ,'init.u3' )
        call get_parameter(p0  ,'init.p'  )
      end if
      rhs(:,:,:,1) = rhs(:,:,:,1) - sponge_shape*(q(:,:,:,1) - rho0    )
      rhs(:,:,:,2) = rhs(:,:,:,2) - sponge_shape*(q(:,:,:,2) - rho0*u10)
      rhs(:,:,:,3) = rhs(:,:,:,3) - sponge_shape*(q(:,:,:,3) - rho0*u20)
      rhs(:,:,:,4) = rhs(:,:,:,4) - sponge_shape*(q(:,:,:,4) - rho0*u30)
      rhs(:,:,:,5) = rhs(:,:,:,5) - sponge_shape*(q(:,:,:,5) - p0      )

    end if

  end subroutine sponge_rhs
!!!=================================================================================================

!!!=================================================================================================
  subroutine sponge_q(q,sponge_shape,qref)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)             :: q
    real(kind=rk), dimension(:,:,:),   intent(in)                :: sponge_shape
    real(kind=rk), dimension(:,:,:,:), intent(in), optional      :: qref
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                      :: i
    real(kind=rk), save                                          :: rho0=0.0_rk, u10, u20, u30, p0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (present(qref)) then

      do i = 1, params%equation%nbr_vars
        q(:,:,:,i) = q(:,:,:,i) - sponge_shape*(q(:,:,:,i) - qref(:,:,:,i))
      end do

    else

      if (rho0 == 0.0_rk) then ! get reference values only once
        call get_parameter(rho0,'init.rho')
        call get_parameter(u10 ,'init.u1' )
        call get_parameter(u20 ,'init.u2' )
        call get_parameter(u30 ,'init.u3' )
        call get_parameter(p0  ,'init.p'  )
      end if
      q(:,:,:,1) = q(:,:,:,1) - sponge_shape*(q(:,:,:,1) - rho0    )
      q(:,:,:,2) = q(:,:,:,2) - sponge_shape*(q(:,:,:,2) - rho0*u10)
      q(:,:,:,3) = q(:,:,:,3) - sponge_shape*(q(:,:,:,3) - rho0*u20)
      q(:,:,:,4) = q(:,:,:,4) - sponge_shape*(q(:,:,:,4) - rho0*u30)
      q(:,:,:,5) = q(:,:,:,5) - sponge_shape*(q(:,:,:,5) - p0      )

    end if

  end subroutine sponge_q
!!!=================================================================================================

!!!=================================================================================================
  subroutine sponge_rhss_correct(rhss,sponge_shape,qs,q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)             :: rhss
    real(kind=rk), dimension(:,:,:),   intent(in)                :: sponge_shape
    real(kind=rk), dimension(:,:,:,:), intent(in)                :: qs,q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable, save         :: rhss_sponge, ui
    real(kind=rk), dimension(:,:,:)  , allocatable, save         :: rho
    integer                                                      :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(.not.allocated(rhss_sponge)) then
      allocate(rhss_sponge(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
      allocate(        rho(params%geom%n1b,params%geom%n2b,params%geom%n3b                         ))
      allocate(         ui(params%geom%n1b,params%geom%n2b,params%geom%n3b,3                       ))
    end if

    rho = give_r(q)
    ui  = give_u(q)

    ! This trafo. is based on the adjoint of the Euler eqs incl. a sponge.
    ! This ensures that the direct and adjoint sponge are consistent.

    !-(A^T)^(-1) * C1^T
    !(A.transpose()).inv()
    ![1, -u_1/rho, -u_2/rho, -u_3/rho,         0]
    ![0,    1/rho,        0,        0,         0]
    ![0,        0,    1/rho,        0,         0]
    ![0,        0,        0,    1/rho,         0]
    ![0,        0,        0,        0, gamma - 1]

    rhss_sponge(:,:,:,1) = 1.0_rk*qs(:,:,:,1) &
        -ui(:,:,:,1)/rho*qs(:,:,:,2) &
        -ui(:,:,:,2)/rho*qs(:,:,:,3) &
        -ui(:,:,:,3)/rho*qs(:,:,:,4)
    rhss_sponge(:,:,:,2) = 1.0_rk/rho*qs(:,:,:,2)
    rhss_sponge(:,:,:,3) = 1.0_rk/rho*qs(:,:,:,3)
    rhss_sponge(:,:,:,4) = 1.0_rk/rho*qs(:,:,:,4)
    rhss_sponge(:,:,:,5) = (params%material%gamma - 1.0_rk)*qs(:,:,:,5)

    do i = 1, params%equation%nbr_vars
      rhss(:,:,:,i) = rhss(:,:,:,i) + sponge_shape*rhss_sponge(:,:,:,i)
    end do

  end subroutine sponge_rhss_correct
!!!=================================================================================================

!!!=================================================================================================
  subroutine sponge_qs(qs,sponge_shape)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)             :: qs
    real(kind=rk), dimension(:,:,:),   intent(in)                :: sponge_shape
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                      :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1, params%equation%nbr_vars
      qs(:,:,:,i) = qs(:,:,:,i) - sponge_shape*qs(:,:,:,i)
    end do

  end subroutine sponge_qs
!!!=================================================================================================

!!!=================================================================================================
  subroutine empty_sponge(y,sponge_shape,a_opt,b_opt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout), optional   :: y
    real(kind=rk), dimension(:,:,:),   intent(in), optional      :: sponge_shape
    real(kind=rk), dimension(:,:,:,:), intent(in), optional      :: a_opt,b_opt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine empty_sponge
!!!=================================================================================================

!!!=================================================================================================
  subroutine empty_sponge_inrhs_warning(y,a_opt,b_opt) ! Only needed for compatibility with old rhs routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)             :: y
    real(kind=rk), dimension(:,:,:,:), intent(in), optional      :: a_opt,b_opt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    y = y
    print*,"WARNING in sponge_layer_core.f90:empty_sponge_inrhs_warning: deprecated sponge/_adjoint is called in rhs."
    print*,"FIX: remove all direct sponge calls in RHS and use sponge BC types from boundary_conditions.f90/.dat instead."

  end subroutine empty_sponge_inrhs_warning
!!!=================================================================================================

!!!=================================================================================================
  subroutine empty_sponge_rhss_warning ! Only needed for compatibility with old rhs routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    print*,"WARNING in sponge_layer_core.f90:empty_sponge_rhss_warning: BC type 400 sponge_rhss is called without q."
    print*,"FIX: call set_boundary_condition_rhss in RHS with q as third argument."

  end subroutine empty_sponge_rhss_warning
!!!=================================================================================================

end module sponge_layer_core