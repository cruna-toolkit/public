module nonreflecting

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                                                   !
  ! Please note                                                                                       !
  ! Because general trafo. fcts. between primitive and rhs variables are not implemented              !
  ! this module is only valid for a euler_rho_rhou_p rhs:                                             !
  ! - euler_rho_rhou_p_trafos are called here                                                         !
  ! - back-trafos. explictly for a euler_rho_rhou_p rhs are hardcoded below                           !
  !                                                                                                   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use data_geom
  use data_rhs
  use parameter
  use euler_rho_rhou_p_trafos , only: give_r, give_u, give_p, give_c

  private

  ! interface wrapper
  public :: characteristic_q, characteristic_qs, characteristic_rhs

  public :: characteristic_q_k1d, characteristic_q_k4d
  public :: characteristic_rhs_k1d, characteristic_rhs_k4d
  public :: characteristic_qs_k1d, characteristic_qs_k4d
  public :: characteristic_zero_grad_init

  interface characteristic_q
     module procedure characteristic_q_k1d
     module procedure characteristic_q_k4d
  end interface characteristic_q

  interface characteristic_rhs
     module procedure characteristic_rhs_k1d
     module procedure characteristic_rhs_k4d
  end interface characteristic_rhs

  interface characteristic_qs
     module procedure characteristic_qs_k1d
     module procedure characteristic_qs_k4d
  end interface characteristic_qs

contains

!!!=================================================================================================
  subroutine characteristic_q_k1d(q,aim,kvector,shape)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),intent(inout)                       :: q              ! [rho,rho u,rho v,rho w,p]
    real(kind=rk), dimension(:,:,:,:),intent(in)                          :: aim            ! [rho,rho u,rho v,rho w,p]
    real(kind=rk), dimension(3)      ,intent(in)                          :: kvector
    real(kind=rk), dimension(:,:,:)  ,intent(in), optional                :: shape
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),allocatable                         :: dq,u
    real(kind=rk)                                                         :: kx,ky,kz       ! -normalvector
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: c, r, p
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: R1,R2,R3,R4,R5 ! riemann invariants / waves
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: lamS,lamP,lamM
    integer                                                               :: n1,n2,n3,n4,i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n1 = size(q,1)
    n2 = size(q,2)
    n3 = size(q,3)
    n4 = size(q,4)

    allocate(dq(n1,n2,n3,n4))
    allocate(u(n1,n2,n3,3),c(n1,n2,n3),r(n1,n2,n3),p(n1,n2,n3))
    allocate(R1(n1,n2,n3),R2(n1,n2,n3),R3(n1,n2,n3),R4(n1,n2,n3),R5(n1,n2,n3))
    allocate(lamS(n1,n2,n3),lamP(n1,n2,n3),lamM(n1,n2,n3))

    ! trafo to primitive variables rho, u, v, w, p
    r = give_r(q)
    u = give_u(q)
    p = give_p(q)
    c = give_c(q)

    ! q  : last value of internal rhs calculation
    ! aim: can be any reference
    ! dq is total change required by aim
    dq(:,:,:,1)   = give_r(aim) - r
    dq(:,:,:,2:4) = give_u(aim) - u
    dq(:,:,:,5)   = give_p(aim) - p

    ! Apply shape if provided
    if (present(shape)) then
      do i = 1, n4
        dq(:,:,:,i) = dq(:,:,:,i) * shape
      end do
    end if

    ! negative boundary normal vector (points inside the domain)
    kx  = -kvector(1)
    ky  = -kvector(2)
    kz  = -kvector(3)

    ! transfrom from primitive variables rho, u, v, w, p to characteristic waves
    R1 = kx*dq(:,:,:,1)                  + kz*dq(:,:,:,3) - ky*dq(:,:,:,4) - kx*dq(:,:,:,5)/(c**2)
    R2 = ky*dq(:,:,:,1) - kz*dq(:,:,:,2)                  + kx*dq(:,:,:,4) - ky*dq(:,:,:,5)/(c**2)
    R3 = kz*dq(:,:,:,1) + ky*dq(:,:,:,2) - kx*dq(:,:,:,3)                  - kz*dq(:,:,:,5)/(c**2)
    R4 =                  kx*dq(:,:,:,2) + ky*dq(:,:,:,3) + kz*dq(:,:,:,4) +    dq(:,:,:,5)/(r*c)
    R5 =                - kx*dq(:,:,:,2) - ky*dq(:,:,:,3) - kz*dq(:,:,:,4) +    dq(:,:,:,5)/(r*c)

    ! characteristic velocities (eigenvalues) multiplied with a negative boundary normal vector kx, ky, kz
    lamS = kx*u(:,:,:,1) + ky*u(:,:,:,2) + kz*u(:,:,:,3)
    lamP = lamS + c
    lamM = lamS - c

    ! remove outgoing characteristic waves with negative: eigenvalue*negative normal vector
    ! the incoming wave part i.e. the allowed change is left
    where(lamS < 0._rk) R1 = 0._rk
    where(lamS < 0._rk) R2 = 0._rk
    where(lamS < 0._rk) R3 = 0._rk
    where(lamP < 0._rk) R4 = 0._rk
    where(lamM < 0._rk) R5 = 0._rk

    ! backtrafo from characteristic waves to primitive variables rho, u, v, w, p
    ! only the incoming wave part of dq is left
    dq(:,:,:,1) =  kx*R1 + ky*R2 + kz*R3 + r/(2.0_rk*c)*R4 + r/(2.0_rk*c)*R5
    dq(:,:,:,2) =        - kz*R2 + ky*R3 +    0.5_rk*kx*R4 -    0.5_rk*kx*R5
    dq(:,:,:,3) =  kz*R1         - kx*R3 +    0.5_rk*ky*R4 -    0.5_rk*ky*R5
    dq(:,:,:,4) = -ky*R1 + kx*R2         +    0.5_rk*kz*R4 -    0.5_rk*kz*R5
    dq(:,:,:,5) =                        +   0.5_rk*r*c*R4 +   0.5_rk*r*c*R5

    ! add allowed change to last value of internal rhs calculation q
    !q = q + dq
    r = r + dq(:,:,:,1)
    u = u + dq(:,:,:,2:4)
    p = p + dq(:,:,:,5)

    ! backtrafo
    q(:,:,:,1) = r
    q(:,:,:,2) = r*u(:,:,:,1)
    q(:,:,:,3) = r*u(:,:,:,2)
    q(:,:,:,4) = r*u(:,:,:,3)
    q(:,:,:,5) = p

  end subroutine characteristic_q_k1d
!!!=================================================================================================

!!!=================================================================================================
  subroutine characteristic_q_k4d(q,aim,kvector,shape)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),intent(inout)                       :: q              ! [rho,rho u,rho v,rho w,p]
    real(kind=rk), dimension(:,:,:,:),intent(in)                          :: aim            ! [rho,rho u,rho v,rho w,p]
    real(kind=rk), dimension(:,:,:,:),intent(in)                          :: kvector
    real(kind=rk), dimension(:,:,:)  ,intent(in), optional                :: shape
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),allocatable                         :: dq,u
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: kx,ky,kz      ! -normalvector
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: c, r, p
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: R1,R2,R3,R4,R5 ! riemann invariants / waves
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: lamS,lamP,lamM
    integer                                                               :: n1,n2,n3,n4,i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n1 = size(q,1)
    n2 = size(q,2)
    n3 = size(q,3)
    n4 = size(q,4)

    allocate(dq(n1,n2,n3,n4))
    allocate(u(n1,n2,n3,3),c(n1,n2,n3),r(n1,n2,n3),p(n1,n2,n3))
    allocate(kx(n1,n2,n3),ky(n1,n2,n3),kz(n1,n2,n3))
    allocate(R1(n1,n2,n3),R2(n1,n2,n3),R3(n1,n2,n3),R4(n1,n2,n3),R5(n1,n2,n3))
    allocate(lamS(n1,n2,n3),lamP(n1,n2,n3),lamM(n1,n2,n3))

    ! trafo to primitive variables rho, u, v, w, p
    r = give_r(q)
    u = give_u(q)
    p = give_p(q)
    c = give_c(q)

    ! q  : last value of internal rhs calculation
    ! aim: can be any reference
    ! dq is total change required by aim
    dq(:,:,:,1)   = give_r(aim) - r
    dq(:,:,:,2:4) = give_u(aim) - u
    dq(:,:,:,5)   = give_p(aim) - p

    ! Apply shape if provided
    if (present(shape)) then
      do i = 1, n4
        dq(:,:,:,i) = dq(:,:,:,i) * shape
      end do
    end if

    ! negative boundary normal vector (points inside the domain)
    kx  = -kvector(:,:,:,1)
    ky  = -kvector(:,:,:,2)
    kz  = -kvector(:,:,:,3)

    ! transfrom from primitive variables rho, u, v, w, p to characteristic waves
    R1 = kx*dq(:,:,:,1)                  + kz*dq(:,:,:,3) - ky*dq(:,:,:,4) - kx*dq(:,:,:,5)/(c**2)
    R2 = ky*dq(:,:,:,1) - kz*dq(:,:,:,2)                  + kx*dq(:,:,:,4) - ky*dq(:,:,:,5)/(c**2)
    R3 = kz*dq(:,:,:,1) + ky*dq(:,:,:,2) - kx*dq(:,:,:,3)                  - kz*dq(:,:,:,5)/(c**2)
    R4 =                  kx*dq(:,:,:,2) + ky*dq(:,:,:,3) + kz*dq(:,:,:,4) +    dq(:,:,:,5)/(r*c)
    R5 =                - kx*dq(:,:,:,2) - ky*dq(:,:,:,3) - kz*dq(:,:,:,4) +    dq(:,:,:,5)/(r*c)

    ! characteristic velocities (eigenvalues) multiplied with a negative boundary normal vector kx, ky, kz
    lamS = kx*u(:,:,:,1) + ky*u(:,:,:,2) + kz*u(:,:,:,3)
    lamP = lamS + c
    lamM = lamS - c

    ! remove outgoing characteristic waves with negative: eigenvalue*negative normal vector
    ! the incoming wave part i.e. the allowed change is left
    where(lamS < 0._rk) R1 = 0._rk
    where(lamS < 0._rk) R2 = 0._rk
    where(lamS < 0._rk) R3 = 0._rk
    where(lamP < 0._rk) R4 = 0._rk
    where(lamM < 0._rk) R5 = 0._rk

    ! backtrafo from characteristic waves to primitive variables rho, u, v, w, p
    ! only the incoming wave part of dq is left
    dq(:,:,:,1) =  kx*R1 + ky*R2 + kz*R3 + r/(2.0_rk*c)*R4 + r/(2.0_rk*c)*R5
    dq(:,:,:,2) =        - kz*R2 + ky*R3 +    0.5_rk*kx*R4 -    0.5_rk*kx*R5
    dq(:,:,:,3) =  kz*R1         - kx*R3 +    0.5_rk*ky*R4 -    0.5_rk*ky*R5
    dq(:,:,:,4) = -ky*R1 + kx*R2         +    0.5_rk*kz*R4 -    0.5_rk*kz*R5
    dq(:,:,:,5) =                        +   0.5_rk*r*c*R4 +   0.5_rk*r*c*R5

    ! add allowed change to last value of internal rhs calculation q
    !q = q + dq
    r = r + dq(:,:,:,1)
    u = u + dq(:,:,:,2:4)
    p = p + dq(:,:,:,5)

    ! backtrafo
    q(:,:,:,1) = r
    q(:,:,:,2) = r*u(:,:,:,1)
    q(:,:,:,3) = r*u(:,:,:,2)
    q(:,:,:,4) = r*u(:,:,:,3)
    q(:,:,:,5) = p

  end subroutine characteristic_q_k4d
!!!=================================================================================================

!!!=================================================================================================
  subroutine characteristic_rhs_k1d(rhs,q_in,aim,kvector,shape)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),intent(inout)                       :: rhs            ! [rho,rho u,rho v,rho w,p]
    real(kind=rk), dimension(:,:,:,:),intent(in)                          :: q_in, aim      ! [rho,rho u,rho v,rho w,p]
    real(kind=rk), dimension(3)      ,intent(in)                          :: kvector
    real(kind=rk), dimension(:,:,:)  ,intent(in), optional                :: shape
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),allocatable                         :: dq,u
    real(kind=rk)                                                         :: kx,ky,kz       ! -normalvector
    real(kind=rk)                                                         :: dt, inv_dt
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: c, r, p, tmp
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: R1,R2,R3,R4,R5 ! riemann invariants / waves
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: lamS,lamP,lamM
    integer                                                               :: n1,n2,n3,n4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    dt = params%time%dt
    inv_dt = 1.0_rk/dt

    n1 = size(rhs,1)
    n2 = size(rhs,2)
    n3 = size(rhs,3)
    n4 = size(rhs,4)

    allocate(dq(n1,n2,n3,n4))
    allocate(u(n1,n2,n3,3),c(n1,n2,n3),r(n1,n2,n3),p(n1,n2,n3),tmp(n1,n2,n3))
    allocate(R1(n1,n2,n3),R2(n1,n2,n3),R3(n1,n2,n3),R4(n1,n2,n3),R5(n1,n2,n3))
    allocate(lamS(n1,n2,n3),lamP(n1,n2,n3),lamM(n1,n2,n3))

    ! trafo. to primitive variables rho, u, v, w, p
    r = give_r(q_in)
    u = give_u(q_in)
    c = give_c(q_in)

    ! q  : last value of internal rhs calculation
    ! aim: can be any reference
    ! dq is total change required by aim
    dq(:,:,:,1)   = give_r(aim) - r
    dq(:,:,:,2:4) = give_u(aim) - u
    dq(:,:,:,5)   = give_p(aim) - give_p(q_in)

    ! trafo. rhs to primitive variables and combine inner rhs with needed change (due to outer aim)
    ! 
    tmp = dt/r
    dq(:,:,:,1) = dq(:,:,:,1) + rhs(:,:,:,1)*dt
    dq(:,:,:,2) = dq(:,:,:,2) + rhs(:,:,:,2)*tmp
    dq(:,:,:,3) = dq(:,:,:,3) + rhs(:,:,:,3)*tmp
    dq(:,:,:,4) = dq(:,:,:,4) + rhs(:,:,:,4)*tmp
    dq(:,:,:,5) = dq(:,:,:,5) + rhs(:,:,:,5)*dt

    ! negative boundary normal vector (points inside the domain)
    kx  = -kvector(1)
    ky  = -kvector(2)
    kz  = -kvector(3)

    ! transfrom from primitive variables rho, u, v, w, p to characteristic waves
    R1 = kx*dq(:,:,:,1)                  + kz*dq(:,:,:,3) - ky*dq(:,:,:,4) - kx*dq(:,:,:,5)/(c**2)
    R2 = ky*dq(:,:,:,1) - kz*dq(:,:,:,2)                  + kx*dq(:,:,:,4) - ky*dq(:,:,:,5)/(c**2)
    R3 = kz*dq(:,:,:,1) + ky*dq(:,:,:,2) - kx*dq(:,:,:,3)                  - kz*dq(:,:,:,5)/(c**2)
    R4 =                  kx*dq(:,:,:,2) + ky*dq(:,:,:,3) + kz*dq(:,:,:,4) +    dq(:,:,:,5)/(r*c)
    R5 =                - kx*dq(:,:,:,2) - ky*dq(:,:,:,3) - kz*dq(:,:,:,4) +    dq(:,:,:,5)/(r*c)

    ! characteristic velocities (eigenvalues) multiplied with a negative boundary normal vector kx, ky, kz
    lamS = kx*u(:,:,:,1) + ky*u(:,:,:,2) + kz*u(:,:,:,3)
    lamP = lamS + c
    lamM = lamS - c

    ! remove outgoing characteristic waves with negative: eigenvalue*negative normal vector
    ! the incoming wave part i.e. the allowed change is left
    where(lamS < 0._rk) R1 = 0._rk
    where(lamS < 0._rk) R2 = 0._rk
    where(lamS < 0._rk) R3 = 0._rk
    where(lamP < 0._rk) R4 = 0._rk
    where(lamM < 0._rk) R5 = 0._rk

    ! backtrafo from characteristic waves to primitive variables rho, u, v, w, p
    ! only the incoming wave part of dq is left
    dq(:,:,:,1) =  kx*R1 + ky*R2 + kz*R3 + r/(2.0_rk*c)*R4 + r/(2.0_rk*c)*R5
    dq(:,:,:,2) =        - kz*R2 + ky*R3 +    0.5_rk*kx*R4 -    0.5_rk*kx*R5
    dq(:,:,:,3) =  kz*R1         - kx*R3 +    0.5_rk*ky*R4 -    0.5_rk*ky*R5
    dq(:,:,:,4) = -ky*R1 + kx*R2         +    0.5_rk*kz*R4 -    0.5_rk*kz*R5
    dq(:,:,:,5) =                        +   0.5_rk*r*c*R4 +   0.5_rk*r*c*R5

    ! add allowed change to last value of internal rhs calculation q
    r = r + dq(:,:,:,1)

    ! backtrafo and apply shape if provided
    if (present(shape)) then
      tmp = shape*inv_dt
      rhs(:,:,:,1) = (1._rk-shape)*rhs(:,:,:,1) + dq(:,:,:,1)*tmp
      rhs(:,:,:,5) = (1._rk-shape)*rhs(:,:,:,5) + dq(:,:,:,5)*tmp
      tmp = tmp*r
      rhs(:,:,:,2) = (1._rk-shape)*rhs(:,:,:,2) + dq(:,:,:,2)*tmp
      rhs(:,:,:,3) = (1._rk-shape)*rhs(:,:,:,3) + dq(:,:,:,3)*tmp
      rhs(:,:,:,4) = (1._rk-shape)*rhs(:,:,:,4) + dq(:,:,:,4)*tmp
    else
      rhs(:,:,:,1) = dq(:,:,:,1)*inv_dt
      tmp = inv_dt*r
      rhs(:,:,:,2) = dq(:,:,:,2)*tmp
      rhs(:,:,:,3) = dq(:,:,:,3)*tmp
      rhs(:,:,:,4) = dq(:,:,:,4)*tmp
      rhs(:,:,:,5) = dq(:,:,:,5)*inv_dt
    end if

  end subroutine characteristic_rhs_k1d
!!!=================================================================================================

!!!=================================================================================================
  subroutine characteristic_rhs_k4d(rhs,q_in,aim,kvector,shape)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),intent(inout)                       :: rhs      ! [rho,rho u,rho v,rho w,p]
    real(kind=rk), dimension(:,:,:,:),intent(in)                          :: q_in,aim ! [rho,rho u,rho v,rho w,p]
    real(kind=rk), dimension(:,:,:,:),intent(in)                          :: kvector
    real(kind=rk), dimension(:,:,:)  ,intent(in), optional                :: shape
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),allocatable                         :: dq,u
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: kx,ky,kz ! -normalvector
    real(kind=rk)                                                         :: dt, inv_dt
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: c, r, p, tmp
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: R1,R2,R3,R4,R5 ! riemann invariants / waves
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: lamS,lamP,lamM
    integer                                                               :: n1,n2,n3,n4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    dt = params%time%dt
    inv_dt = 1.0_rk/dt

    n1 = size(rhs,1)
    n2 = size(rhs,2)
    n3 = size(rhs,3)
    n4 = size(rhs,4)

    allocate(dq(n1,n2,n3,n4))
    allocate(u(n1,n2,n3,3),c(n1,n2,n3),r(n1,n2,n3),p(n1,n2,n3),tmp(n1,n2,n3))
    allocate(kx(n1,n2,n3),ky(n1,n2,n3),kz(n1,n2,n3))
    allocate(R1(n1,n2,n3),R2(n1,n2,n3),R3(n1,n2,n3),R4(n1,n2,n3),R5(n1,n2,n3))
    allocate(lamS(n1,n2,n3),lamP(n1,n2,n3),lamM(n1,n2,n3))

    ! trafo. to primitive variables rho, u, v, w, p
    r = give_r(q_in)
    u = give_u(q_in)
    c = give_c(q_in)

    ! q  : last value of internal rhs calculation
    ! aim: can be any reference
    ! dq is total change required by aim
    dq(:,:,:,1)   = give_r(aim) - r
    dq(:,:,:,2:4) = give_u(aim) - u
    dq(:,:,:,5)   = give_p(aim) - give_p(q_in)

    ! trafo. rhs to primitive variables and combine inner rhs with needed change (due to outer aim)
    tmp = dt/r
    dq(:,:,:,1) = dq(:,:,:,1) + rhs(:,:,:,1)*dt
    dq(:,:,:,2) = dq(:,:,:,2) + rhs(:,:,:,2)*tmp
    dq(:,:,:,3) = dq(:,:,:,3) + rhs(:,:,:,3)*tmp
    dq(:,:,:,4) = dq(:,:,:,4) + rhs(:,:,:,4)*tmp
    dq(:,:,:,5) = dq(:,:,:,5) + rhs(:,:,:,5)*dt

    ! negative boundary normal vector (points inside the domain)
    kx  = -kvector(:,:,:,1)
    ky  = -kvector(:,:,:,2)
    kz  = -kvector(:,:,:,3)

    ! transfrom from primitive variables rho, u, v, w, p to characteristic waves
    R1 = kx*dq(:,:,:,1)                  + kz*dq(:,:,:,3) - ky*dq(:,:,:,4) - kx*dq(:,:,:,5)/(c**2)
    R2 = ky*dq(:,:,:,1) - kz*dq(:,:,:,2)                  + kx*dq(:,:,:,4) - ky*dq(:,:,:,5)/(c**2)
    R3 = kz*dq(:,:,:,1) + ky*dq(:,:,:,2) - kx*dq(:,:,:,3)                  - kz*dq(:,:,:,5)/(c**2)
    R4 =                  kx*dq(:,:,:,2) + ky*dq(:,:,:,3) + kz*dq(:,:,:,4) +    dq(:,:,:,5)/(r*c)
    R5 =                - kx*dq(:,:,:,2) - ky*dq(:,:,:,3) - kz*dq(:,:,:,4) +    dq(:,:,:,5)/(r*c)

    ! characteristic velocities (eigenvalues) multiplied with a negative boundary normal vector kx, ky, kz
    lamS = kx*u(:,:,:,1) + ky*u(:,:,:,2) + kz*u(:,:,:,3)
    lamP = lamS + c
    lamM = lamS - c

    ! remove outgoing characteristic waves with negative: eigenvalue*negative normal vector
    ! the incoming wave part i.e. the allowed change is left
    where(lamS < 0._rk) R1 = 0._rk
    where(lamS < 0._rk) R2 = 0._rk
    where(lamS < 0._rk) R3 = 0._rk
    where(lamP < 0._rk) R4 = 0._rk
    where(lamM < 0._rk) R5 = 0._rk

    ! backtrafo from characteristic waves to primitive variables rho, u, v, w, p
    ! only the incoming wave part of dq is left
    dq(:,:,:,1) =  kx*R1 + ky*R2 + kz*R3 + r/(2.0_rk*c)*R4 + r/(2.0_rk*c)*R5
    dq(:,:,:,2) =        - kz*R2 + ky*R3 +    0.5_rk*kx*R4 -    0.5_rk*kx*R5
    dq(:,:,:,3) =  kz*R1         - kx*R3 +    0.5_rk*ky*R4 -    0.5_rk*ky*R5
    dq(:,:,:,4) = -ky*R1 + kx*R2         +    0.5_rk*kz*R4 -    0.5_rk*kz*R5
    dq(:,:,:,5) =                        +   0.5_rk*r*c*R4 +   0.5_rk*r*c*R5

    ! add allowed change to last value of internal rhs calculation q
    r = r + dq(:,:,:,1)

    ! backtrafo and apply shape if provided
    if (present(shape)) then
      tmp = shape*inv_dt
      rhs(:,:,:,1) = (1._rk-shape)*rhs(:,:,:,1) + dq(:,:,:,1)*tmp
      rhs(:,:,:,5) = (1._rk-shape)*rhs(:,:,:,5) + dq(:,:,:,5)*tmp
      tmp = tmp*r
      rhs(:,:,:,2) = (1._rk-shape)*rhs(:,:,:,2) + dq(:,:,:,2)*tmp
      rhs(:,:,:,3) = (1._rk-shape)*rhs(:,:,:,3) + dq(:,:,:,3)*tmp
      rhs(:,:,:,4) = (1._rk-shape)*rhs(:,:,:,4) + dq(:,:,:,4)*tmp
    else
      rhs(:,:,:,1) = dq(:,:,:,1)*inv_dt
      tmp = inv_dt*r
      rhs(:,:,:,2) = dq(:,:,:,2)*tmp
      rhs(:,:,:,3) = dq(:,:,:,3)*tmp
      rhs(:,:,:,4) = dq(:,:,:,4)*tmp
      rhs(:,:,:,5) = dq(:,:,:,5)*inv_dt
    end if

  end subroutine characteristic_rhs_k4d
!!!=================================================================================================

!!!=================================================================================================
  subroutine characteristic_qs_k1d(qs,q,kvector,shape)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),intent(inout)                       :: qs ! adjoint [rho,u,v,w,p]
    real(kind=rk), dimension(:,:,:,:),intent(in)                          :: q  ! [rho,rho u,rho v,rho w,p]
    real(kind=rk), dimension(3)      ,intent(in)                          :: kvector
    real(kind=rk), dimension(:,:,:)  ,intent(in), optional                :: shape
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),allocatable                         :: ds
    real(kind=rk)                                                         :: kx,ky,kz
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: c, g_c, u, v, w, ku, kus
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: ruus, rups, tmp
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: drs, dus, dvs, dws, dps
    logical,       dimension(:,:,:)  ,allocatable                         :: Sincoming, Pincoming, Mincoming
    real(kind=rk)                                                         :: gm1
    integer                                                               :: n1,n2,n3,n4,i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n1 = size(qs,1)
    n2 = size(qs,2)
    n3 = size(qs,3)
    n4 = size(qs,4)

    allocate(ds(n1,n2,n3,n4))
    allocate(c(n1,n2,n3),g_c(n1,n2,n3))
    allocate(u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3))
    allocate(ku(n1,n2,n3),kus(n1,n2,n3))
    allocate(ruus(n1,n2,n3),rups(n1,n2,n3),tmp(n1,n2,n3))
    allocate(drs(n1,n2,n3),dus(n1,n2,n3),dvs(n1,n2,n3),dws(n1,n2,n3),dps(n1,n2,n3))

    allocate(Sincoming(n1,n2,n3),Pincoming(n1,n2,n3),Mincoming(n1,n2,n3))

    ! qs : last value of internal rhs calculation
    ! aim: zero if adjoint
    ! dq is total change required by aim = 0
    ds = - qs

    ! apply shape if provided
    if (present(shape)) then
      do i = 1, n4
        ds(:,:,:,i) = ds(:,:,:,i) * shape
      end do
    end if

    ! boundary normal vector (points outside the domain, sign switched by ML because adjoint, LS needs to check)
    kx  = -kvector(1)
    ky  = -kvector(2)
    kz  = -kvector(3)

    gm1 = params%material%gamma - 1_rk

    c    = sqrt(params%material%gamma*(q(:,:,:,5)/q(:,:,:,1)))
    g_c  = gm1/c
    u    = q(:,:,:,2)/q(:,:,:,1)
    v    = q(:,:,:,3)/q(:,:,:,1)
    w    = q(:,:,:,4)/q(:,:,:,1)
    ku   =               kx*u           + ky*v           + kz*w
    kus  =               kx*ds(:,:,:,2) + ky*ds(:,:,:,3) + kz*ds(:,:,:,4)
    ruus = (ds(:,:,:,1) + u*ds(:,:,:,2) +  v*ds(:,:,:,3) +  w*ds(:,:,:,4))/c
    rups = ruus + ds(:,:,:,5)/g_c

    ! check if incoming (external reference has impact, internal contribution is zero)
    Sincoming = ku     < 0._rk
    Pincoming = ku + c < 0._rk
    Mincoming = ku - c < 0._rk

    drs  = 0._rk
    dus  = 0._rk
    dvs  = 0._rk
    dws  = 0._rk
    dps  = 0._rk

    ! u wave (entropie)
    where(Sincoming) drs = ds(:,:,:,1) + ku*kus
    where(Sincoming) dus = ds(:,:,:,2) - kx*kus
    where(Sincoming) dvs = ds(:,:,:,3) - ky*kus
    where(Sincoming) dws = ds(:,:,:,4) - kz*kus
    where(Sincoming) dps =           - g_c*ruus

    ! u + c wave (sound)
    where(Pincoming) tmp = 0.5*(kus+rups)
    where(Pincoming) drs  = drs -  ku*tmp
    where(Pincoming) dus  = dus +  kx*tmp
    where(Pincoming) dvs  = dvs +  ky*tmp
    where(Pincoming) dws  = dws +  kz*tmp
    where(Pincoming) dps  = dps + g_c*tmp

    ! u - c wave (sound)
    where(Mincoming) tmp = 0.5*(kus-rups)
    where(Mincoming) drs  = drs -  ku*tmp
    where(Mincoming) dus  = dus +  kx*tmp
    where(Mincoming) dvs  = dvs +  ky*tmp
    where(Mincoming) dws  = dws +  kz*tmp
    where(Mincoming) dps  = dps - g_c*tmp

    ! add allowed change to last value of internal rhs calculation q
    qs(:,:,:,1) = qs(:,:,:,1) + drs
    qs(:,:,:,2) = qs(:,:,:,2) + dus
    qs(:,:,:,3) = qs(:,:,:,3) + dvs
    qs(:,:,:,4) = qs(:,:,:,4) + dws
    qs(:,:,:,5) = qs(:,:,:,5) + dps

  end subroutine characteristic_qs_k1d
!!!=================================================================================================

!!!=================================================================================================
  subroutine characteristic_qs_k4d(qs,q,kvector,shape)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),intent(inout)                       :: qs ! adjoint [rho,u,v,w,p]
    real(kind=rk), dimension(:,:,:,:),intent(in)                          :: q  ! [rho,rho u,rho v,rho w,p]
    real(kind=rk), dimension(:,:,:,:),intent(in)                          :: kvector
    real(kind=rk), dimension(:,:,:)  ,intent(in), optional                :: shape
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),allocatable                         :: ds
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: kx,ky,kz      ! normalvector
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: c, g_c, u, v, w, ku, kus
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: ruus, rups, tmp
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: drs, dus, dvs, dws, dps
    logical,       dimension(:,:,:)  ,allocatable                         :: Sincoming, Pincoming, Mincoming
    real(kind=rk)                                                         :: gm1
    integer                                                               :: n1,n2,n3,n4,i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n1 = size(qs,1)
    n2 = size(qs,2)
    n3 = size(qs,3)
    n4 = size(qs,4)

    allocate(ds(n1,n2,n3,n4))
    allocate(kx(n1,n2,n3),ky(n1,n2,n3),kz(n1,n2,n3))

    allocate(c(n1,n2,n3),g_c(n1,n2,n3))
    allocate(u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3))
    allocate(ku(n1,n2,n3),kus(n1,n2,n3))
    allocate(ruus(n1,n2,n3),rups(n1,n2,n3),tmp(n1,n2,n3))
    allocate(drs(n1,n2,n3),dus(n1,n2,n3),dvs(n1,n2,n3),dws(n1,n2,n3),dps(n1,n2,n3))

    allocate(Sincoming(n1,n2,n3),Pincoming(n1,n2,n3),Mincoming(n1,n2,n3))

    ! qs : last value of internal rhs calculation
    ! aim: zero if adjoint
    ! dq is total change required by aim = 0
    ds = - qs

    ! apply shape if provided
    if (present(shape)) then
      do i = 1, n4
        ds(:,:,:,i) = ds(:,:,:,i) * shape
      end do
    end if

    ! boundary normal vector (points outside the domain, sign switched by ML because adjoint, LS needs to check)
    kx  = -kvector(:,:,:,1)
    ky  = -kvector(:,:,:,2)
    kz  = -kvector(:,:,:,3)

    gm1 = params%material%gamma - 1_rk

    c    = sqrt(params%material%gamma*(q(:,:,:,5)/q(:,:,:,1)))
    g_c  = gm1/c
    u    = q(:,:,:,2)/q(:,:,:,1)
    v    = q(:,:,:,3)/q(:,:,:,1)
    w    = q(:,:,:,4)/q(:,:,:,1)
    ku   =               kx*u           + ky*v           + kz*w
    kus  =               kx*ds(:,:,:,2) + ky*ds(:,:,:,3) + kz*ds(:,:,:,4)
    ruus = (ds(:,:,:,1) + u*ds(:,:,:,2) +  v*ds(:,:,:,3) +  w*ds(:,:,:,4))/c
    rups = ruus + ds(:,:,:,5)/g_c

    ! check if incoming (external reference has impact, internal contribution is zero)
    Sincoming = ku     < 0._rk
    Pincoming = ku + c < 0._rk
    Mincoming = ku - c < 0._rk

    drs  = 0._rk
    dus  = 0._rk
    dvs  = 0._rk
    dws  = 0._rk
    dps  = 0._rk

    ! u wave (entropie)
    where(Sincoming) drs = ds(:,:,:,1) + ku*kus
    where(Sincoming) dus = ds(:,:,:,2) - kx*kus
    where(Sincoming) dvs = ds(:,:,:,3) - ky*kus
    where(Sincoming) dws = ds(:,:,:,4) - kz*kus
    where(Sincoming) dps =           - g_c*ruus

    ! u + c wave (sound)
    where(Pincoming) tmp = 0.5*(kus+rups)
    where(Pincoming) drs  = drs -  ku*tmp
    where(Pincoming) dus  = dus +  kx*tmp
    where(Pincoming) dvs  = dvs +  ky*tmp
    where(Pincoming) dws  = dws +  kz*tmp
    where(Pincoming) dps  = dps + g_c*tmp

    ! u - c wave (sound)
    where(Mincoming) tmp = 0.5*(kus-rups)
    where(Mincoming) drs  = drs -  ku*tmp
    where(Mincoming) dus  = dus +  kx*tmp
    where(Mincoming) dvs  = dvs +  ky*tmp
    where(Mincoming) dws  = dws +  kz*tmp
    where(Mincoming) dps  = dps - g_c*tmp

    ! add allowed change to last value of internal rhs calculation q
    qs(:,:,:,1) = qs(:,:,:,1) + drs
    qs(:,:,:,2) = qs(:,:,:,2) + dus
    qs(:,:,:,3) = qs(:,:,:,3) + dvs
    qs(:,:,:,4) = qs(:,:,:,4) + dws
    qs(:,:,:,5) = qs(:,:,:,5) + dps

  end subroutine characteristic_qs_k4d
!!!=================================================================================================

!!!=================================================================================================
  subroutine characteristic_zero_grad_init(is2,ie2,js2,je2,ks2,ke2,boundary_conditions_entry)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer,intent(out)                                                  :: is2,ie2,js2,je2,ks2,ke2
    real(kind=rk),dimension(:)      ,intent(in)                          :: boundary_conditions_entry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(3)                                           :: kvector
    integer,dimension(3)                                                 :: shift
    integer                                                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! boundary normal vector pointing outwards of the domain
    kvector = boundary_conditions_entry(8:10)

    ! normalize k-vector (ml)
    kvector(1) = kvector(1)/sqrt(sum(kvector**2))
    kvector(2) = kvector(2)/sqrt(sum(kvector**2))
    kvector(3) = kvector(3)/sqrt(sum(kvector**2))

    ! relative index of inner neighbour
    do i = 1,3
       if (kvector(i)>0) then
          shift(i) = -1
       elseif (kvector(i)<0) then
          shift(i) = 1
       else
          shift(i) = 0
       end if
    end do

    ! absolute index of inner neighbour i.e. aim value index of zero gradient boundary condition
    is2 = nint(boundary_conditions_entry(2))+shift(1)
    ie2 = nint(boundary_conditions_entry(3))+shift(1)
    js2 = nint(boundary_conditions_entry(4))+shift(2)
    je2 = nint(boundary_conditions_entry(5))+shift(2)
    ks2 = nint(boundary_conditions_entry(6))+shift(3)
    ke2 = nint(boundary_conditions_entry(7))+shift(3)

    ! security checks that image contains needed neighbour
    if (is2 <1) then
       write(*,*) "error in nonreflecting.f90:characteristic_zero_grad_init: inner neighbour is2 not part of proc:",params%parallelism%world_image
       stop
    endif
    if (js2 <1) then
       write(*,*) "error in nonreflecting.f90:characteristic_zero_grad_init: inner neighbour js2 not part of proc:",params%parallelism%world_image
       stop
    endif
    if (ks2 <1) then
       write(*,*) "error in nonreflecting.f90:characteristic_zero_grad_init: inner neighbour ks2 not part of proc:",params%parallelism%world_image
       stop
    endif
    if (ie2 >params%geom%n1b) then
       write(*,*) "error in nonreflecting.f90:characteristic_zero_grad_init: inner neighbour ie2 not part of proc:",params%parallelism%world_image
       stop
    endif
    if (je2 >params%geom%n2b) then
       write(*,*) "error in nonreflecting.f90:characteristic_zero_grad_init: inner neighbour je2 not part of proc:",params%parallelism%world_image
       stop
    endif
    if (ke2 >params%geom%n3b) then
       write(*,*) "error in nonreflecting.f90:characteristic_zero_grad_init: inner neighbour ke2 not part of proc:",params%parallelism%world_image
       stop
    endif

  end subroutine characteristic_zero_grad_init
!!!=================================================================================================

! !!!=================================================================================================
!   subroutine characteristic_q0_init(q0_loc_bnd,boundary_conditions_entry)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     real(kind=rk),dimension(:,:,:,:),intent(out)                       :: q0_loc_bnd
!     real(kind=rk),dimension(:)      ,intent(in)                        :: boundary_conditions_entry
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     integer                                                            :: is,ie,js,je,ks,ke
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     is = nint(boundary_conditions_entry(2))
!     ie = nint(boundary_conditions_entry(3))
!     js = nint(boundary_conditions_entry(4))
!     je = nint(boundary_conditions_entry(5))
!     ks = nint(boundary_conditions_entry(6))
!     ke = nint(boundary_conditions_entry(7))

!     q0_loc_bnd = q0(is:ie,js:je,ks:ke,:)

!   end subroutine characteristic_q0_init
! !!!=================================================================================================

end module nonreflecting
