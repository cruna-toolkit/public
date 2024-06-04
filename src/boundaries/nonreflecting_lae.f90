module nonreflecting_lae

  use data_geom
  use data_rhs
  use parameter
  use linear_acoustic_p_u_trafos

  private

  public :: characteristic_q_lae, characteristic_qs_lae

contains

!!!=================================================================================================
  subroutine characteristic_q_lae(q,aim,kvector)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),intent(inout)                       :: q              ! [p,u,v,w]
    real(kind=rk), dimension(:,:,:,:),intent(in)                          :: aim            ! [p,u,v,w]
    real(kind=rk), dimension(3)      ,intent(in)                          :: kvector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),allocatable                         :: du 
    real(kind=rk)                                                         :: kx,ky,kz       ! -normalvector
    real(kind=rk)                                                         :: rho0, p0, c
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: R1,R2,R3,R4 ! riemann invariants / waves
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: lamS,lamP,lamM
    integer                                                               :: n1,n2,n3,n4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    n1 = size(q,1)
    n2 = size(q,2)
    n3 = size(q,3)
    n4 = size(q,4)

    allocate(du(n1,n2,n3,n4))
    allocate(R1(n1,n2,n3),R2(n1,n2,n3),R3(n1,n2,n3),R4(n1,n2,n3))
    allocate(lamS(n1,n2,n3),lamP(n1,n2,n3),lamM(n1,n2,n3))

    rho0 = params%init%rho
    p0   = params%init%p

    ! negative boundary normal vector (points inside the domain)
    kx  = -kvector(1)
    ky  = -kvector(2)
    kz  = -kvector(3)
   
    ! q  : last value of internal rhs calculation
    ! aim: can be any reference
    ! du is total change required by aim
    du = aim - q

    ! speed of sound
    c  = sqrt(params%material%gamma*(p0/rho0))

    ! characteristic velocities (eigenvalues) multiplied with a negative boundary normal vector kx, ky, kz
    lamS =  0
    lamP =  c
    lamM = -c

    if (abs(kx) > 0.01_rk) then
       ! transfrom from primitive variables rho, u, v, w, p to characteristic waves
       R1 =                                - (kx*ky)    *du(:,:,:,2) + (kx**2+kz**2)*du(:,:,:,3) - (kz*ky)      *du(:,:,:,4)
       R2 =                                - (kx*kz)    *du(:,:,:,2) - (ky*kz)      *du(:,:,:,3) + (kx**2+ky**2)*du(:,:,:,4)
       R3 = -(0.5_rk/(rho0*c))*du(:,:,:,1) + (0.5_rk*kx)*du(:,:,:,2) + (0.5_rk*ky)  *du(:,:,:,3) + (0.5_rk*kz)  *du(:,:,:,4)
       R4 =  (0.5_rk/(rho0*c))*du(:,:,:,1) + (0.5_rk*kx)*du(:,:,:,2) + (0.5_rk*ky)  *du(:,:,:,3) + (0.5_rk*kz)  *du(:,:,:,4)

       ! [0               , -kx*ky/K^2    , +(kx^2+ky^2/K^2, -kz*ky/K^2      ]
       ! [0               , -kx*kz/K^2    , -ky*kz/K^2     , +(kx^2+ky^2)/K^2]
       ! [-kz/(2*rho0*K*c), +ky*kz/(2*K^2), +ky*kz/(2*K^2) , +kz^2/(2*K^2)   ]
       ! [+kz/(2*rho0*K*c), +kx*kz/(2*K^2), +ky*kz/(2*K^2) , +kz^2/(2*K^2)   ]

       ! remove outgoing characteristic waves with negative: eigenvalue*negative normal vector
       ! the incoming wave part i.e. the allowed change is left
       where(lamS < 0._rk) R1 = 0._rk
       where(lamS < 0._rk) R2 = 0._rk
       where(lamM < 0._rk) R3 = 0._rk
       where(lamP < 0._rk) R4 = 0._rk

       ! backtrafo from characteristic waves to primitive variables rho, u, v, w, p
       ! only the incoming wave part of du is left
       du(:,:,:,1) =                          - (c*rho0)*R3 + (c*rho0)*R4
       du(:,:,:,2) = -(ky/kx)*R1 - (kz/kx)*R2 + kx      *R3 + kx      *R4
       du(:,:,:,3) =          R1              + ky      *R3 + ky      *R4
       du(:,:,:,4) =                       R2 + kz      *R3 + kz      *R4
       
       ! [0     , 0     , -K*c*rho0/kz, K*c*rho0/kz]
       ! [-ky/kx, -kz/kx, kx/kz       , kx/kz      ]
       ! [1     , 0     , ky/kz       , ky/kz      ]
       ! [0     , 1     , 1           , 1          ]

    elseif (abs(ky) > 0.01_rk) then

       ! transfrom from primitive variables rho, u, v, w, p to characteristic waves
       R1 =                                  (ky**2+kz**2)*du(:,:,:,2) - (ky*kx)    *du(:,:,:,3) - (kz*kx)      *du(:,:,:,4)
       R2 =                                - (kx*kz)      *du(:,:,:,2) - (ky*kz)    *du(:,:,:,3) + (ky**2+kx**2)*du(:,:,:,4)
       R3 = -(0.5_rk/(rho0*c))*du(:,:,:,1) + (0.5_rk*kx)  *du(:,:,:,2) + (0.5_rk*ky)*du(:,:,:,3) + (0.5_rk*kz)  *du(:,:,:,4)
       R4 =  (0.5_rk/(rho0*c))*du(:,:,:,1) + (0.5_rk*kx)  *du(:,:,:,2) + (0.5_rk*ky)*du(:,:,:,3) + (0.5_rk*kz)  *du(:,:,:,4)

       ! remove outgoing characteristic waves with negative: eigenvalue*negative normal vector
       ! the incoming wave part i.e. the allowed change is left
       where(lamS < 0._rk) R1 = 0._rk
       where(lamS < 0._rk) R2 = 0._rk
       where(lamM < 0._rk) R3 = 0._rk
       where(lamP < 0._rk) R4 = 0._rk

       ! backtrafo from characteristic waves to primitive variables rho, u, v, w, p
       ! only the incoming wave part of du is left
       du(:,:,:,1) =                          - (c*rho0)*R3 + (c*rho0)*R4
       du(:,:,:,2) =          R1              + kx      *R3 + kx      *R4
       du(:,:,:,3) = -(kx/ky)*R1 - (kz/ky)*R2 + ky      *R3 + ky      *R4
       du(:,:,:,4) =                       R2 + kz      *R3 + kz      *R4
       
    elseif (abs(kz) > 0.01_rk) then
       ! transfrom from primitive variables rho, u, v, w, p to characteristic waves
       R1 =                                  (ky**2+kz**2)*du(:,:,:,2) - (kx*ky)      *du(:,:,:,3) - (kx*kz)    *du(:,:,:,4)
       R2 =                                - (kx*ky)      *du(:,:,:,2) + (kx**2+kz**2)*du(:,:,:,3) - (ky*kz)    *du(:,:,:,4)
       R3 = -(0.5_rk/(rho0*c))*du(:,:,:,1) + (0.5_rk*kx)  *du(:,:,:,2) + (0.5_rk*ky)  *du(:,:,:,3) + (0.5_rk*kz)*du(:,:,:,4)
       R4 =  (0.5_rk/(rho0*c))*du(:,:,:,1) + (0.5_rk*kx)  *du(:,:,:,2) + (0.5_rk*ky)  *du(:,:,:,3) + (0.5_rk*kz)*du(:,:,:,4)

       ! remove outgoing characteristic waves with negative: eigenvalue*negative normal vector
       ! the incoming wave part i.e. the allowed change is left
       where(lamS < 0._rk) R1 = 0._rk
       where(lamS < 0._rk) R2 = 0._rk
       where(lamM < 0._rk) R3 = 0._rk
       where(lamP < 0._rk) R4 = 0._rk

       ! backtrafo from characteristic waves to primitive variables p, u, v, w
       ! only the incoming wave part of du is left
       du(:,:,:,1) =                          - (c*rho0)*R3 + (c*rho0)*R4
       du(:,:,:,2) =          R1              + kx      *R3 + kx      *R4
       du(:,:,:,3) =                       R2 + ky      *R3 + ky      *R4
       du(:,:,:,4) = -(kx/kz)*R1 - (ky/kz)*R2 + kz      *R3 + kz      *R4
    
    else
       write(*,*) "kx, ky, kz were zero, please code"
    end if

    ! add allowed change to last value of internal rhs calculation q
    q  = q + du

  end subroutine characteristic_q_lae
!!!=================================================================================================

!!!=================================================================================================
  subroutine characteristic_qs_lae(qs,q,kvector)

    real(kind=rk), dimension(:,:,:,:),intent(inout)                       :: qs              ! [p,u,v,w]
    real(kind=rk), dimension(:,:,:,:),intent(in)                          :: q               ! [p,u,v,w]
    real(kind=rk), dimension(3)      ,intent(in)                          :: kvector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),allocatable                         :: du,ds 
    real(kind=rk)                                                         :: kx,ky,kz       ! -normalvector
    real(kind=rk)                                                         :: rho0, p0, c
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: R1,R2,R3,R4 ! riemann invariants / waves
    real(kind=rk), dimension(:,:,:)  ,allocatable                         :: lamS,lamP,lamM
    integer                                                               :: n1,n2,n3,n4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n1 = size(q,1)
    n2 = size(q,2)
    n3 = size(q,3)
    n4 = size(q,4)

    allocate(du(n1,n2,n3,n4),ds(n1,n2,n3,n4))
    allocate(R1(n1,n2,n3),R2(n1,n2,n3),R3(n1,n2,n3),R4(n1,n2,n3))
    allocate(lamS(n1,n2,n3),lamP(n1,n2,n3),lamM(n1,n2,n3))

    rho0 = params%init%rho
    p0   = params%init%p

    ! negative boundary normal vector (points inside the domain)
    kx  = -kvector(1)
    ky  = -kvector(2)
    kz  = -kvector(3)

    ! qs  : last value of internal rhs adjoint calculation
    ! aim: reference, here zero
    ! ds is total change required by aim
    ! ds = aim - q
    ds = - qs

    ! speed of sound
    c  = sqrt(params%material%gamma*(p0/rho0))

    ! characteristic velocities (eigenvalues) multiplied with a negative boundary normal vector kx, ky, kz
    lamS =  0
    lamP =  c
    lamM = -c

    if (abs(kx) > 0.01_rk) then
       ! transfrom from primitive variables rho, u, v, w, p to characteristic waves
       R1 =                              - (kx*ky)    *ds(:,:,:,2) + (kx**2+kz**2)*ds(:,:,:,3) - (kz*ky)      *ds(:,:,:,4)
       R2 =                              - (kx*kz)    *ds(:,:,:,2) - (ky*kz)      *ds(:,:,:,3) + (kx**2+ky**2)*ds(:,:,:,4)
       R3 = -(0.5_rk*rho0*c)*ds(:,:,:,1) + (0.5_rk*kx)*ds(:,:,:,2) + (0.5_rk*ky)  *ds(:,:,:,3) + (0.5_rk*kz)  *ds(:,:,:,4)
       R4 =  (0.5_rk*rho0*c)*ds(:,:,:,1) + (0.5_rk*kx)*ds(:,:,:,2) + (0.5_rk*ky)  *ds(:,:,:,3) + (0.5_rk*kz)  *ds(:,:,:,4)

       ! [0               , -kx*ky/K^2    , +(kx^2+ky^2/K^2, -kz*ky/K^2      ]
       ! [0               , -kx*kz/K^2    , -ky*kz/K^2     , +(kx^2+ky^2)/K^2]
       ! [-kz/(2*rho0*K*c), +ky*kz/(2*K^2), +ky*kz/(2*K^2) , +kz^2/(2*K^2)   ]
       ! [+kz/(2*rho0*K*c), +kx*kz/(2*K^2), +ky*kz/(2*K^2) , +kz^2/(2*K^2)   ]

       ! remove outgoing characteristic waves with positive: positive*negative normal vector
       ! the incoming wave part i.e. the allowed change is left
       where(lamS > 0._rk) R1 = 0._rk
       where(lamS > 0._rk) R2 = 0._rk
       where(lamM > 0._rk) R3 = 0._rk
       where(lamP > 0._rk) R4 = 0._rk
 ! backtrafo from characteristic waves to primitive variables rho, u, v, w, p
       ! only the incoming wave part of du is left
       ds(:,:,:,1) =                          - (1.0_rk/(c*rho0))*R3 + (1.0_rk/(c*rho0))*R4
       ds(:,:,:,2) = -(ky/kx)*R1 - (kz/kx)*R2 + kx               *R3 + kx               *R4
       ds(:,:,:,3) =          R1              + ky               *R3 + ky               *R4
       ds(:,:,:,4) =                       R2 + kz               *R3 + kz               *R4

       ! [0     , 0     , -K*c*rho0/kz, K*c*rho0/kz]
       ! [-ky/kx, -kz/kx, kx/kz       , kx/kz      ]
       ! [1     , 0     , ky/kz       , ky/kz      ]
       ! [0     , 1     , 1           , 1          ]

    elseif (abs(ky) > 0.01_rk) then

       ! transfrom from primitive variables rho, u, v, w, p to characteristic waves
       R1 =                                (ky**2+kz**2)*ds(:,:,:,2) - (ky*kx)    *ds(:,:,:,3) - (kz*kx)      *ds(:,:,:,4)
       R2 =                              - (kx*kz)      *ds(:,:,:,2) - (ky*kz)    *ds(:,:,:,3) + (ky**2+kx**2)*ds(:,:,:,4)
       R3 = -(0.5_rk*rho0*c)*ds(:,:,:,1) + (0.5_rk*kx)  *ds(:,:,:,2) + (0.5_rk*ky)*ds(:,:,:,3) + (0.5_rk*kz)  *ds(:,:,:,4)
       R4 =  (0.5_rk*rho0*c)*ds(:,:,:,1) + (0.5_rk*kx)  *ds(:,:,:,2) + (0.5_rk*ky)*ds(:,:,:,3) + (0.5_rk*kz)  *ds(:,:,:,4)

       ! remove outgoing characteristic waves with positive: eigenvalue*positive normal vector
       ! the incoming wave part i.e. the allowed change is left
       where(lamS > 0._rk) R1 = 0._rk
       where(lamS > 0._rk) R2 = 0._rk
       where(lamM > 0._rk) R3 = 0._rk
       where(lamP > 0._rk) R4 = 0._rk

       ! backtrafo from characteristic waves to primitive variables rho, u, v, w, p
       ! only the incoming wave part of du is left
       ds(:,:,:,1) =                          - (1.0_rk/(c*rho0))*R3 + (1.0_rk/(c*rho0))*R4
       ds(:,:,:,2) =          R1              + kx               *R3 + kx               *R4
       ds(:,:,:,3) = -(kx/ky)*R1 - (kz/ky)*R2 + ky               *R3 + ky               *R4
       ds(:,:,:,4) =                       R2 + kz               *R3 + kz               *R4

    elseif (abs(kz) > 0.01_rk) then
       ! transfrom from primitive variables rho, u, v, w, p to characteristic waves
       R1 =                                (ky**2+kz**2)*ds(:,:,:,2) - (kx*ky)      *ds(:,:,:,3) - (kx*kz)    *ds(:,:,:,4)
       R2 =                              - (kx*ky)      *ds(:,:,:,2) + (kx**2+kz**2)*ds(:,:,:,3) - (ky*kz)    *ds(:,:,:,4)
       R3 = -(0.5_rk*rho0*c)*ds(:,:,:,1) + (0.5_rk*kx)  *ds(:,:,:,2) + (0.5_rk*ky)  *ds(:,:,:,3) + (0.5_rk*kz)*ds(:,:,:,4)
       R4 =  (0.5_rk*rho0*c)*ds(:,:,:,1) + (0.5_rk*kx)  *ds(:,:,:,2) + (0.5_rk*ky)  *ds(:,:,:,3) + (0.5_rk*kz)*ds(:,:,:,4)

       ! remove outgoing characteristic waves with positive: eigenvalue*positive normal vector
       ! the incoming wave part i.e. the allowed change is left
       where(lamS > 0._rk) R1 = 0._rk
       where(lamS > 0._rk) R2 = 0._rk
       where(lamM > 0._rk) R3 = 0._rk
       where(lamP > 0._rk) R4 = 0._rk

       ! backtrafo from characteristic waves to primitive variables p, u, v, w
       ! only the incoming wave part of ds is left
       ds(:,:,:,1) =                          - (1.0_rk/(c*rho0))*R3 + (1.0_rk/(c*rho0))*R4
       ds(:,:,:,2) =          R1              + kx               *R3 + kx               *R4
       ds(:,:,:,3) =                       R2 + ky               *R3 + ky               *R4
       ds(:,:,:,4) = -(kx/kz)*R1 - (ky/kz)*R2 + kz               *R3 + kz               *R4

    else
       write(*,*) "kx, ky, kz were zero, please code"
    end if

    ! add allowed change to last value of internal rhs calculation qs
    qs  = qs + ds

  end subroutine characteristic_qs_lae

!!!=================================================================================================

end module nonreflecting_lae
