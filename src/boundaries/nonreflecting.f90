module nonreflecting

  use data_geom
  use data_rhs
  use parameter
  use euler_rho_rhou_p_trafos

  private

  public :: characteristic_q, characteristic_qs
  public :: characteristic_rhs ! (200)
  public :: characteristic_zero_grad_init

contains

!!!=================================================================================================
  subroutine characteristic_rhs(rhs,q,aim,kvector)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),intent(inout)                       :: rhs
    real(kind=rk), dimension(:,:,:,:),intent(in)                          :: q,aim          ! [rho,rho u,rho v,rho w,p]
    real(kind=rk), dimension(3)      ,intent(in)                          :: kvector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),allocatable                         :: du
    !real(kind=rk), dimension(:,:,:)  ,allocatable                         :: kx,ky,kz       ! -normalvector 
    real(kind=rk)                                                         :: kx,ky,kz       ! -normalvector
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: R1,R2,R3,R4,R5 ! riemann invariants / waves
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: lamS,lamP,lamM,c
    integer                                                               :: n1,n2,n3,n4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n1 = size(rhs,1)
    n2 = size(rhs,2)
    n3 = size(rhs,3)
    n4 = size(rhs,4)

    allocate(du(n1,n2,n3,n4))
    !allocate(kx(n1,n2,n3),ky(n1,n2,n3),kz(n1,n2,n3))
    allocate(c(n1,n2,n3),R1(n1,n2,n3),R2(n1,n2,n3),R3(n1,n2,n3),R4(n1,n2,n3),R5(n1,n2,n3))
    allocate(lamS(n1,n2,n3),lamP(n1,n2,n3),lamM(n1,n2,n3))

    ! q  : last value of internal rhs calculation
    ! aim: can be any reference
    ! du is total change required by aim
    du = aim - q

    ! negative boundary normal vector (points inside the domain)
    kx  = -kvector(1)
    ky  = -kvector(2)
    kz  = -kvector(3)

    ! speed of sound
    c  = sqrt(params%material%gamma*(q(:,:,:,5)/q(:,:,:,1)))

    ! transfrom from primitive variables rho, u, v, w, p to characteristic waves
    R1 = kx*du(:,:,:,1)                  + kz*du(:,:,:,3) - ky*du(:,:,:,4) -      kx*du(:,:,:,5)/(c**2)
    R2 = ky*du(:,:,:,1) - kz*du(:,:,:,2)                  + kx*du(:,:,:,4) -      ky*du(:,:,:,5)/(c**2)
    R3 = kz*du(:,:,:,1) + ky*du(:,:,:,2) - kx*du(:,:,:,3)                  -      kz*du(:,:,:,5)/(c**2)
    R4 =                  kx*du(:,:,:,2) + ky*du(:,:,:,3) + kz*du(:,:,:,4) + du(:,:,:,5)/(q(:,:,:,1)*c)
    R5 =                - kx*du(:,:,:,2) - ky*du(:,:,:,3) - kz*du(:,:,:,4) + du(:,:,:,5)/(q(:,:,:,1)*c)

    ! characteristic velocities (eigenvalues) multiplied with a negative boundary normal vector kx, ky, kz
    lamS = kx*q(:,:,:,2) + ky*q(:,:,:,3) + kz*q(:,:,:,4) 
    lamP = kx*q(:,:,:,2) + ky*q(:,:,:,3) + kz*q(:,:,:,4)  + c
    lamM = kx*q(:,:,:,2) + ky*q(:,:,:,3) + kz*q(:,:,:,4)  - c

    ! remove outgoing characteristic waves with negative: eigenvalue*negative normal vector
    ! the incoming wave part i.e. the allowed change is left
    where(lamS < 0._rk) R1 = 0._rk
    where(lamS < 0._rk) R2 = 0._rk
    where(lamS < 0._rk) R3 = 0._rk
    where(lamP < 0._rk) R4 = 0._rk
    where(lamM < 0._rk) R5 = 0._rk

    ! backtrafo from characteristic waves to primitive variables rho, u, v, w, p
    ! only the incoming wave part of du is left
    du(:,:,:,1) =  kx*R1 + ky*R2 + kz*R3 + q(:,:,:,1)/(2.0_rk*c)*R4 + q(:,:,:,1)/(2.0_rk*c)*R5
    du(:,:,:,2) =        - kz*R2 + ky*R3 +             0.5_rk*kx*R4 -             0.5_rk*kx*R5
    du(:,:,:,3) =  kz*R1         - kx*R3 +             0.5_rk*ky*R4 -             0.5_rk*ky*R5
    du(:,:,:,4) = -ky*R1 + kx*R2         +             0.5_rk*kz*R4 -             0.5_rk*kz*R5
    du(:,:,:,5) =                        +   0.5_rk*q(:,:,:,1)*c*R4 +   0.5_rk*q(:,:,:,1)*c*R5

    ! add allowed change to last value of internal rhs calculation q
    rhs = rhs + du/params%time%dt ! full dt (not substep) is more stable according to julius

  end subroutine characteristic_rhs
!!!=================================================================================================

!!!=================================================================================================
  subroutine characteristic_q(q,aim,kvector)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),intent(inout)                       :: q              ! [rho,rho u,rho v,rho w,p]
    real(kind=rk), dimension(:,:,:,:),intent(in)                          :: aim            ! [rho,rho u,rho v,rho w,p]
    real(kind=rk), dimension(3)      ,intent(in)                          :: kvector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),allocatable                         :: du,u
    !real(kind=rk), dimension(:,:,:)  ,allocatable                         :: kx,ky,kz      ! -normalvector 
    real(kind=rk)                                                         :: kx,ky,kz       ! -normalvector
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: c, r, p
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: R1,R2,R3,R4,R5 ! riemann invariants / waves
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: lamS,lamP,lamM
    integer                                                               :: n1,n2,n3,n4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n1 = size(q,1)
    n2 = size(q,2)
    n3 = size(q,3)
    n4 = size(q,4)

    allocate(du(n1,n2,n3,n4))
    allocate(u(n1,n2,n3,3),c(n1,n2,n3),r(n1,n2,n3),p(n1,n2,n3))
    !allocate(kx(n1,n2,n3),ky(n1,n2,n3),kz(n1,n2,n3))
    allocate(R1(n1,n2,n3),R2(n1,n2,n3),R3(n1,n2,n3),R4(n1,n2,n3),R5(n1,n2,n3))
    allocate(lamS(n1,n2,n3),lamP(n1,n2,n3),lamM(n1,n2,n3))

    r = give_r(q)
    u = give_u(q)
    p = give_p(q)
    c = give_c(q)

    ! q  : last value of internal rhs calculation
    ! aim: can be any reference
    ! du is total change required by aim
    du(:,:,:,1)   = give_r(aim) - r
    du(:,:,:,2:4) = give_u(aim) - u
    du(:,:,:,5)   = give_p(aim) - p

    ! negative boundary normal vector (points inside the domain)
    kx  = -kvector(1)
    ky  = -kvector(2)
    kz  = -kvector(3)

    ! transfrom from primitive variables rho, u, v, w, p to characteristic waves
    R1 = kx*du(:,:,:,1)                  + kz*du(:,:,:,3) - ky*du(:,:,:,4) - kx*du(:,:,:,5)/(c**2)
    R2 = ky*du(:,:,:,1) - kz*du(:,:,:,2)                  + kx*du(:,:,:,4) - ky*du(:,:,:,5)/(c**2)
    R3 = kz*du(:,:,:,1) + ky*du(:,:,:,2) - kx*du(:,:,:,3)                  - kz*du(:,:,:,5)/(c**2)
    R4 =                  kx*du(:,:,:,2) + ky*du(:,:,:,3) + kz*du(:,:,:,4) +    du(:,:,:,5)/(r*c)
    R5 =                - kx*du(:,:,:,2) - ky*du(:,:,:,3) - kz*du(:,:,:,4) +    du(:,:,:,5)/(r*c)

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
    ! only the incoming wave part of du is left
    du(:,:,:,1) =  kx*R1 + ky*R2 + kz*R3 + r/(2.0_rk*c)*R4 + r/(2.0_rk*c)*R5
    du(:,:,:,2) =        - kz*R2 + ky*R3 +    0.5_rk*kx*R4 -    0.5_rk*kx*R5
    du(:,:,:,3) =  kz*R1         - kx*R3 +    0.5_rk*ky*R4 -    0.5_rk*ky*R5
    du(:,:,:,4) = -ky*R1 + kx*R2         +    0.5_rk*kz*R4 -    0.5_rk*kz*R5
    du(:,:,:,5) =                        +   0.5_rk*r*c*R4 +   0.5_rk*r*c*R5

    ! add allowed change to last value of internal rhs calculation q
    !q = q + du
    r = r + du(:,:,:,1)
    u = u + du(:,:,:,2:4)
    p = p + du(:,:,:,5)

    ! backtrafo
    q(:,:,:,1) = r
    q(:,:,:,2) = r*u(:,:,:,1)
    q(:,:,:,3) = r*u(:,:,:,2)
    q(:,:,:,4) = r*u(:,:,:,3)
    q(:,:,:,5) = p

  end subroutine characteristic_q
!!!=================================================================================================

!!!=================================================================================================
  subroutine characteristic_qs(qs,q,kvector)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),intent(inout)                       :: qs
    real(kind=rk), dimension(:,:,:,:),intent(in)                          :: q ! [rho,rho u,rho v,rho w,p]
    real(kind=rk), dimension(3)      ,intent(in)                          :: kvector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),allocatable                         :: ds
    !real(kind=rk), dimension(:,:,:)  ,allocatable                         :: kx,ky,kz      ! normalvector 
    real(kind=rk)                                                         :: kx,ky,kz
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: c, g_c, u, v, w, ku, kus
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: ruus, rups, tmp
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: drs, dus, dvs, dws, dps
    logical,       dimension(:,:,:)  ,allocatable                         :: Sincoming, Pincoming, Mincoming
    real(kind=rk)                                                         :: gm1
    integer                                                               :: n1,n2,n3,n4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n1 = size(qs,1)
    n2 = size(qs,2)
    n3 = size(qs,3)
    n4 = size(qs,4)

    allocate(ds(n1,n2,n3,n4))
    !allocate(kx(n1,n2,n3),ky(n1,n2,n3),kz(n1,n2,n3))

    allocate(c(n1,n2,n3),g_c(n1,n2,n3))
    allocate(u(n1,n2,n3),v(n1,n2,n3),w(n1,n2,n3))
    allocate(ku(n1,n2,n3),kus(n1,n2,n3))
    allocate(ruus(n1,n2,n3),rups(n1,n2,n3),tmp(n1,n2,n3))
    allocate(drs(n1,n2,n3),dus(n1,n2,n3),dvs(n1,n2,n3),dws(n1,n2,n3),dps(n1,n2,n3))

    allocate(Sincoming(n1,n2,n3),Pincoming(n1,n2,n3),Mincoming(n1,n2,n3))

    ! qs : last value of internal rhs calculation
    ! aim: zero if adjoint
    ! du is total change required by aim = 0
    ds = - qs

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

  end subroutine characteristic_qs
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
