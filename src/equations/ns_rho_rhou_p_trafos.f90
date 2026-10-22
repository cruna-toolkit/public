module ns_rho_rhou_p_trafos

  use parameter

  private

  public give_T, give_c, give_r, give_u, give_p
  public give_W, give_mu, give_cp, give_lambda
  public give_dmu_dT, give_R_const, give_Pr
  public apply_matrix_ATm1

contains

!!!=================================================================================================
  function give_W() result(W)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                           :: W
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), save                                     :: W_const
    logical      , save                                     :: calc_const = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(calc_const.eqv..true.) then
       ! get mol weight
       ! default: air (source: Baehr, H. D., "Thermodynamik, 7.edition" (1989) (page 207)
       call get_parameter(W_const,'material.W',default=0.0289647_rk)

       ! save state
       calc_const = .false.
    end if

    W = W_const

  end function give_W
!!!=================================================================================================

!!!=================================================================================================
  function give_R_const() result(R)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                           :: R
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), save                                     :: R_const
    logical      , save                                     :: calc_const = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(calc_const.eqv..true.) then
       call get_parameter(R_const, "material.R", default = 8.3144598_rk)

       ! save state
       calc_const = .false.
    end if

    R = R_const

  end function give_R_const
!!!=================================================================================================

!!!=================================================================================================
  function give_Pr() result(Pr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                           :: Pr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), save                                     :: Pr_const
    logical      , save                                     :: calc_const = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(calc_const.eqv..true.) then
       ! get Prandtl-number
       call get_parameter(Pr_const,'material.Pr',default = 0.71_rk)

       ! save state
       calc_const = .false.
    end if

    Pr = Pr_const

  end function give_Pr
!!!=================================================================================================


!!!=================================================================================================
  function give_T(q) result(T)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(in)           :: q
    real(kind=rk), dimension(size(q,1),size(q,2),size(q,3)) :: T
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                           :: W
    real(kind=rk), save                                     :: W_INV_R_const
    logical      , save                                     :: calc_const = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(calc_const.eqv..true.) then
       ! get mol weight
       W = give_W()

       ! compute W_INV_R
       W_INV_R_const = W/params%material%R

       ! save state
       calc_const = .false.
    end if

    T = W_INV_R_const*q(:,:,:,5)/q(:,:,:,1)

  end function give_T
!!!=================================================================================================

!!!=================================================================================================
  function give_c(q) result(c)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(in)           :: q
    real(kind=rk), dimension(size(q,1),size(q,2),size(q,3)) :: c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                           :: gamma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    gamma = params%material%gamma

    c = sqrt(gamma*q(:,:,:,5)/q(:,:,:,1))

  end function give_c
!!!=================================================================================================

!!!=================================================================================================
  function give_r(q) result(rho)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(in)           :: q
    real(kind=rk), dimension(size(q,1),size(q,2),size(q,3)) :: rho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    rho = q(:,:,:,1)

  end function give_r
!!!=================================================================================================

!!!=================================================================================================
  function give_p(q) result(p)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(in)           :: q
    real(kind=rk), dimension(size(q,1),size(q,2),size(q,3)) :: p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    p = q(:,:,:,5)

  end function give_p
!!!=================================================================================================

!!!=================================================================================================
  function give_u(q) result(u)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(in)               :: q
    real(kind=rk), dimension(size(q,1),size(q,2),size(q,3),3)   :: u
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                     :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1,3
       u(:,:,:,i) = q(:,:,:,i+1)/q(:,:,:,1)
    end do

  end function give_u
!!!=================================================================================================

!!!=================================================================================================
  function give_mu(T) result(mu)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:)   ,intent(in)           :: T
    real(kind=rk), dimension(size(T,1),size(T,2),size(T,3)) :: mu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), parameter                                :: S1 = 110.4_rk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), save                                     :: mu0 = 0.0_rk, T0 = 0.0_rk
    logical      , save                                     :: calc_const = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(calc_const.eqv..true.) then
       ! get reference temperature
       ! default: air (source: F. M. White, Viscous Fluid Flow, 2nd ed., McGraw-Hill, (1991), page 29)
       call get_parameter( T0,'material.T0' ,default = 273.0_rk)
       call get_parameter(mu0,'material.mu0')

       ! save state
       calc_const = .false.
    end if

    mu = mu0*(T/T0)**(1.5_rk)*((T0 + S1)/(T + S1))

  end function give_mu

!!!=================================================================================================
  function give_dmu_dT(T,mu) result(dmu_dT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:)   ,intent(in)           :: T
    real(kind=rk), dimension(:,:,:)   ,intent(in)           :: mu
    real(kind=rk), dimension(size(T,1),size(T,2),size(T,3)) :: dmu_dT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), parameter                                :: S1 = 110.4_rk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), save                                     :: mu0 = 0.0_rk, T0 = 0.0_rk
    logical      , save                                     :: calc_const = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(calc_const.eqv..true.) then
       ! get reference temperature
       ! default: air (source: F. M. White, Viscous Fluid Flow, 2nd ed., McGraw-Hill, (1991), page 29)
       call get_parameter( T0,'material.T0' ,default = 273.0_rk)
       call get_parameter(mu0,'material.mu0')

       ! save state
       calc_const = .false.
    end if

    dmu_dT = mu0 * (T/T0)**0.5_rk * ((3.0_rk*S1 + T)*(S1 + T0)/(2.0_rk*(T + S1)**2.0_rk))

  end function give_dmu_dT



!!!=================================================================================================
  function give_lambda(mu) result(lambda)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:)    , intent(in)            :: mu
    real(kind=rk), dimension(size(mu,1),size(mu,2),size(mu,3)) :: lambda
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                              :: W
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), save                                        :: cp = 0.0_rk, Pr = 0.0_rk
    logical      , save                                        :: calc_const = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(calc_const.eqv..true.) then
       ! get Prandtl-number
       call get_parameter(Pr,'material.Pr',default = 0.71_rk)

       ! get cp
       W  = give_W()
       cp = params%material%gamma*(params%material%R/W)/(params%material%gamma - 1.0_rk)

       ! save state
       calc_const = .false.
    end if

    lambda = mu*cp/Pr

  end function give_lambda

!!!=================================================================================================
  function give_cp() result(cp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                           :: cp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), save                                     :: W, cp_const
    logical      , save                                     :: calc_const = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(calc_const.eqv..true.) then
       ! calc constant cp
       W        = give_W()
       cp_const = params%material%gamma*(params%material%R/W)/(params%material%gamma - 1.0_rk)

       ! save state
       calc_const = .false.
    end if

    cp = cp_const

  end function give_cp

!!!=================================================================================================
  subroutine apply_matrix_ATm1(vec_out,vec_in,q) ! vec: typically qs (adjoint state) or g (drive of the adjoint)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(out)              :: vec_out    
    real(kind=rk), dimension(:,:,:,:), intent(in)               :: vec_in
    real(kind=rk), dimension(:,:,:,:), intent(in)               :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable              :: ui
    real(kind=rk), dimension(:,:,:)  , allocatable              :: rho
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(ui(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))
    ui  = give_u(q)
    
    allocate(rho(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    rho = give_r(q)

    !   [1, -u_1/rho, -u_2/rho, -u_3/rho,         0]
    !   [0,    1/rho,        0,        0,         0]
    !   [0,        0,    1/rho,        0,         0]
    !   [0,        0,        0,    1/rho,         0]
    !   [0,        0,        0,        0, gamma - 1]

    vec_out(:,:,:,1) = 1.0_rk*vec_in(:,:,:,1) &
         -ui(:,:,:,1)/rho*vec_in(:,:,:,2) &
         -ui(:,:,:,2)/rho*vec_in(:,:,:,3) &
         -ui(:,:,:,3)/rho*vec_in(:,:,:,4)
    vec_out(:,:,:,2) = 1.0_rk/rho*vec_in(:,:,:,2)
    vec_out(:,:,:,3) = 1.0_rk/rho*vec_in(:,:,:,3)
    vec_out(:,:,:,4) = 1.0_rk/rho*vec_in(:,:,:,4)
    vec_out(:,:,:,5) = (params%material%gamma - 1.0_rk)*vec_in(:,:,:,5)

  end subroutine apply_matrix_ATm1
!!!=================================================================================================

end module ns_rho_rhou_p_trafos
