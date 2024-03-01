module diff_pen_core_o6_i_c_p

  use parameter

  private

  public :: diff_i1, diff_i2, diff_i3

  character(len=max_length_parameter), parameter, public :: diff_type = 'periodic'
  character(len=max_length_parameter), parameter, public :: diff_name = 'diff_pen_core_o6_i_c_p'

  real(kind=rk), dimension(6)        , save              :: stencil      ! right explicit stencil A

contains

!!!=================================================================================================
  subroutine diff_i1(q,dx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
    real(kind=rk)                , intent(in)          :: dx
    real(kind=rk), dimension(:)  , allocatable, save   :: dl,dc,du,du2   ! left implicit stencil Btri
    real(kind=rk),                              save   :: fbrat          ! left implicit stencil B non-tridiagonal corrections
    real(kind=rk), dimension(:)  , allocatable, save   :: fbnorm         ! "
    integer      , dimension(:)  , allocatable, save   :: ipiv           ! left implicit stencil Btri
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(dl)) then
       ! Init. right explicit stencil A
       call init_right_explicit_diff_part(stencil)
       ! Init. left implicit stencil B
       call init_left_implicit_diff_part(dl,dc,du,du2,fbrat,fbnorm,ipiv,size(q,1))
    end if

    call diff_main(q,dx,stencil,dl,dc,du,du2,fbrat,fbnorm,ipiv)

  end subroutine diff_i1
!!!==================================================================================================

!!!=================================================================================================
  subroutine diff_i2(q,dx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
    real(kind=rk)                , intent(in)          :: dx
    real(kind=rk), dimension(:)  , allocatable, save   :: dl,dc,du,du2   ! left implicit stencil Btri
    real(kind=rk),                              save   :: fbrat          ! left implicit stencil B non-tridiagonal corrections
    real(kind=rk), dimension(:)  , allocatable, save   :: fbnorm         ! "
    integer      , dimension(:)  , allocatable, save   :: ipiv           ! left implicit stencil Btri
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(dl)) then
       ! Init. right explicit stencil A
       call init_right_explicit_diff_part(stencil)
       ! Init. left implicit stencil B
       call init_left_implicit_diff_part(dl,dc,du,du2,fbrat,fbnorm,ipiv,size(q,1))
    end if

    call diff_main(q,dx,stencil,dl,dc,du,du2,fbrat,fbnorm,ipiv)

  end subroutine diff_i2
!!!==================================================================================================

!!!=================================================================================================
  subroutine diff_i3(q,dx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
    real(kind=rk)                , intent(in)          :: dx
    real(kind=rk), dimension(:)  , allocatable, save   :: dl,dc,du,du2   ! left implicit stencil Btri
    real(kind=rk),                              save   :: fbrat          ! left implicit stencil B non-tridiagonal corrections
    real(kind=rk), dimension(:)  , allocatable, save   :: fbnorm         ! "
    integer      , dimension(:)  , allocatable, save   :: ipiv           ! left implicit stencil Btri
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(dl)) then
       ! Init. right explicit stencil A
       call init_right_explicit_diff_part(stencil)
       ! Init. left implicit stencil B
       call init_left_implicit_diff_part(dl,dc,du,du2,fbrat,fbnorm,ipiv,size(q,1))
    end if

    call diff_main(q,dx,stencil,dl,dc,du,du2,fbrat,fbnorm,ipiv)

  end subroutine diff_i3
!!!==================================================================================================


!!!=================================================================================================
  subroutine diff_main(q,dx,stencil,dl,dc,du,du2,fbrat,fbnorm,ipiv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! solve BU' = AU = F by a BLAS implementation
!!! i.e. implicit derivation along first dimension
!!! B is tridiagonal
!!! U(') are matrices
!!! stencil values are provided by Lele1992
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(:,:), intent(inout)       :: q
    real(kind=rk)                 , intent(in)          :: dx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(:)  , intent(in)          :: stencil        ! right explicit stencil A
    real(kind=rk) , dimension(:)  , intent(in)          :: dl,dc,du,du2   ! left implicit stencil  Btri
    real(kind=rk) ,                 intent(in)          :: fbrat          ! left implicit stencil  B non-tridiagonal corrections
    real(kind=rk) , dimension(:)  , intent(in)          :: fbnorm         ! "
    integer       , dimension(:)  , intent(in)          :: ipiv           ! left implicit stencil  Btri
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(:,:), allocatable         :: q_tmp
    real(kind=rk)                                       :: a,b
    integer                                             :: n1,n2
    integer                                             :: i,j,info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n1 = size(q,1)
    n2 = size(q,2)

    allocate(q_tmp(n1,n2))

    a  = stencil(1)
    b  = stencil(2)

    ! compute A*U (central)
    ! inner:    sixth order (incl. implicit B matrix)
    ! boundary: forth order (incl. implicit B matrix)
    q_tmp(1   ,:) = a*(q(1+1   ,:) - q(n1    ,:)) + b*(q(1+2   ,:) - q(n1-1  ,:))
    q_tmp(2   ,:) = a*(q(2+1   ,:) - q(2-1   ,:)) + b*(q(2+2   ,:) - q(n1    ,:))
    do j = 1,n2
       do i = 3,(n1-2)
          q_tmp(i,j) = a*(q(i+1,j) - q(i-1,j)) + b*(q(i+2,j) - q(i-2,j))
       end do
    end do
    q_tmp(n1-1,:) = a*(q(n1-1+1,:) - q(n1-1-1,:)) + b*(q(1     ,:) - q(n1-1-2,:))
    q_tmp(n1  ,:) = a*(q(1     ,:) - q(n1  -1,:)) + b*(q(2     ,:) - q(n1  -2,:))

    ! dx factor
    q = q_tmp/dx

    ! solve Btri U' = F
    call dgttrs('no transpose', n1, n2, dl, dc, du, du2, ipiv, q, n1, info)
    if(info.ne.0) then
       write(*,*) "error in diff_pen_core_o6_i_c_o.f90:diff_main:error with id: ",info
    end if

    ! correct non-tridiagonal part of B
    do i = 1,n2
       q(:,i) = q(:,i) - fbnorm*(q(1,i) + fbrat*q(n1,i))
    end do

  end subroutine diff_main
!!!==================================================================================================



!!!==================================================================================================
  subroutine init_right_explicit_diff_part(stencil) ! output args
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Init. right explicit stencil A for o6_i_c_o
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(2)  , intent(out)         :: stencil
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                       :: a,b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    a       = 14.00_rk / 9.0_rk / 2.0_rk
    b       =  1.00_rk / 9.0_rk / 4.0_rk
    stencil = (/ a, b /)

  end subroutine init_right_explicit_diff_part
!!!==================================================================================================



!!!==================================================================================================
  subroutine init_left_implicit_diff_part(dl,dc,du,du2,fbrat,fbnorm,ipiv,n1) ! output args
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! LU decomposition for dgttrf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:), allocatable, intent(out) :: dl,dc,du,du2
    real(kind=rk),                            intent(out) :: fbrat
    real(kind=rk), dimension(:), allocatable, intent(out) :: fbnorm
    integer      , dimension(:), allocatable, intent(out) :: ipiv
    integer                                 , intent(in)  :: n1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:), allocatable              :: fa,fb
    real(kind=rk)                                         :: alpha
    integer                                               :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! allocates
    allocate(   dl(n1-1 ))
    allocate(   dc(n1   ))
    allocate(   du(n1-1 ))
    allocate(  du2(n1-2 ))
    allocate( ipiv(n1   ))
    allocate(   fa(n1   ))
    allocate(   fb(n1   ))
    allocate(fbnorm(n1  ))

    ! setup full B (pure central)
    alpha    = 1.0_rk / 3.0_rk
    dl       = alpha
    dc       = 1.0_rk
    du       = alpha
    ! setup non-tridiagonal B part as vectors (central)
    fa       = 0._rk
    fb       = 0._rk
    ! setup non-tridiagonal B part as vectors (boundaries)
    fa(1 )   = -dc(1) ! free choice, this one is recommended for unknown reasons
    fa(n1)   = +alpha
    fb(1 )   = +1._rk
    fb(n1)   = -alpha/dc(1)
    ! setup corrected tridiagonal Btri
    dc(n1)   = dc(n1) + (alpha*alpha/dc(1))
    dc(1)    = dc(1)  + dc(1)

    ! init LU decomposition by dgttrf
    call dgttrf(n1, dl, dc, du, du2, ipiv, info)
    if(info.ne.0) then
       write(*,*) "error in diff_pen_core_o6_i_c_o.f90:init_left_implicit_diff_part:error with id: ",info
    end if

    ! solve Btri y = fa
    call dgttrs('no transpose', n1,  1, dl, dc, du, du2, ipiv, fa, n1, info)
    if(info.ne.0) then
       write(*,*) "error in diff_pen_core_o6_i_c_o.f90:init_left_implicit_diff_part:error with id: ",info
    end if

    ! init prefactors for normalization
    fbnorm = fa*fb(1)/(1._rk+sum(fb*fa))
    fbrat  = fb(n1)/fb(1)

  end subroutine init_left_implicit_diff_part
!!!==================================================================================================


end module diff_pen_core_o6_i_c_p
