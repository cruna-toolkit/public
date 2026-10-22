module sponge_layer_core
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Collection of various sponge routines                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use parameter
  use sponge_layer_trafos

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
     module procedure sponge_rhss_consistent
     module procedure sponge_rhss_inconsistent
     module procedure empty_sponge_rhss_warning
     ! Empty version needed for compatibility with old rhs routines
     ! because their set_boundary_condition_rhss calls (boundary_conditions.f90) lack q argument.
  end interface sponge_rhss

  private sponge_rhss_consistent
  private sponge_rhss_inconsistent
  private empty_sponge_rhss_warning ! to be removed by LS :)

contains

!!!=================================================================================================
  subroutine sponge_rhs(rhs,sponge_shape,q,qref)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)             :: rhs
    real(kind=rk), dimension(:,:,:),   intent(in)                :: sponge_shape
    real(kind=rk), dimension(:,:,:,:), intent(in)                :: q
    real(kind=rk), dimension(:,:,:,:), intent(in)                :: qref
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                      :: i
    real(kind=rk), save                                          :: rho0=0.0_rk, u10, u20, u30, p0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1, size(rhs,4)
       rhs(:,:,:,i) = rhs(:,:,:,i) - sponge_shape*(q(:,:,:,i) - qref(:,:,:,i))
    end do

  end subroutine sponge_rhs
!!!=================================================================================================

!!!=================================================================================================
  subroutine sponge_q(q,sponge_shape,qref)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)             :: q
    real(kind=rk), dimension(:,:,:),   intent(in)                :: sponge_shape
    real(kind=rk), dimension(:,:,:,:), intent(in)                :: qref
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                      :: i
    real(kind=rk), save                                          :: rho0=0.0_rk, u10, u20, u30, p0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1, params%equation%nbr_vars
       q(:,:,:,i) = q(:,:,:,i) - sponge_shape*(q(:,:,:,i) - qref(:,:,:,i))
    end do

  end subroutine sponge_q
!!!=================================================================================================

!!!=================================================================================================
  subroutine sponge_rhss_consistent(rhss,sponge_shape,qs,q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)             :: rhss
    real(kind=rk), dimension(:,:,:),   intent(in)                :: sponge_shape
    real(kind=rk), dimension(:,:,:,:), intent(in)                :: qs,q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable, save         :: rhss_sponge ! simply ATm1*qs
    integer                                                      :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(.not.allocated(rhss_sponge)) then
       allocate(rhss_sponge(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    end if

    rhss_sponge = 0.0_rk

    do i = 1, params%equation%nbr_vars
       rhss_sponge(:,:,:,i) = rhss_sponge(:,:,:,i) + sponge_shape*qs(:,:,:,i) ! sponge_shape*(qs-0.0_rk) 
    end do

    ! This ensures that the direct and adjoint sponge are consistent.
    !-(A^T)^(-1)*rhss_sponge
    !(A.transpose()).inv()
    ![1, -u_1/rho, -u_2/rho, -u_3/rho,         0]
    ![0,    1/rho,        0,        0,         0]
    ![0,        0,    1/rho,        0,         0]
    ![0,        0,        0,    1/rho,         0]
    ![0,        0,        0,        0, gamma - 1]
    call apply_matrix_ATm1(rhss_sponge,qs,q)

    do i = 1, params%equation%nbr_vars
       rhss(:,:,:,i) = rhss(:,:,:,i) + rhss_sponge(:,:,:,i)
    end do

  end subroutine sponge_rhss_consistent
!!!=================================================================================================


!!!=================================================================================================
  subroutine sponge_rhss_inconsistent(rhss,sponge_shape,qs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)             :: rhss
    real(kind=rk), dimension(:,:,:),   intent(in)                :: sponge_shape
    real(kind=rk), dimension(:,:,:,:), intent(in)                :: qs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                      :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1, params%equation%nbr_vars
       rhss(:,:,:,i) = rhss(:,:,:,i) + sponge_shape*qs(:,:,:,i)  ! sponge_shape*(qs-0.0_rk)
    end do

  end subroutine sponge_rhss_inconsistent
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
       qs(:,:,:,i) = qs(:,:,:,i) - sponge_shape*qs(:,:,:,i)      ! sponge_shape*(qs-0.0_rk)
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
