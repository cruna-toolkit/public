module filter_pen_core_o2_i_c_p

  use parameter

  private

  public :: filter_i1, filter_i2, filter_i3

  character(len=max_length_parameter), parameter, public :: filter_type = 'periodic'
  character(len=max_length_parameter), parameter, public :: filter_name = 'filter_pen_core_o2_i_c_p'

  real(kind=rk), dimension(2)        , save              :: stencil      ! right explicit stencil A

contains

!!!=================================================================================================
  subroutine filter_i1(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
    real(kind=rk), dimension(:)  , allocatable, save   :: dl,dc,du,du2   ! left implicit stencil B
    real(kind=rk),                              save   :: fbrat          ! left implicit stencil B non-tridiagonal corrections
    real(kind=rk), dimension(:)  , allocatable, save   :: fbnorm         ! "
    integer      , dimension(:)  , allocatable, save   :: ipiv           ! left implicit stencil B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(dl)) then
       ! Init. right explicit stencil A
       call init_right_explicit_filter_part(stencil)
       ! Init. left implicit stencil B
       call init_left_implicit_filter_part(dl,dc,du,du2,fbrat,fbnorm,ipiv,size(q,1))
    end if

    call filter_main(q,stencil,dl,dc,du,du2,fbrat,fbnorm,ipiv)

  end subroutine filter_i1
!!!==================================================================================================

!!!=================================================================================================
  subroutine filter_i2(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
    real(kind=rk), dimension(:)  , allocatable, save   :: dl,dc,du,du2   ! left implicit stencil B
    real(kind=rk),                              save   :: fbrat          ! left implicit stencil B non-tridiagonal corrections
    real(kind=rk), dimension(:)  , allocatable, save   :: fbnorm         ! "
    integer      , dimension(:)  , allocatable, save   :: ipiv           ! left implicit stencil B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(dl)) then
       ! Init. right explicit stencil A
       call init_right_explicit_filter_part(stencil)
       ! Init. left implicit stencil B
       call init_left_implicit_filter_part(dl,dc,du,du2,fbrat,fbnorm,ipiv,size(q,1))
    end if

    call filter_main(q,stencil,dl,dc,du,du2,fbrat,fbnorm,ipiv)

  end subroutine filter_i2
!!!==================================================================================================

!!!=================================================================================================
  subroutine filter_i3(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
    real(kind=rk), dimension(:)  , allocatable, save   :: dl,dc,du,du2   ! left implicit stencil  B
    real(kind=rk),                              save   :: fbrat          ! left implicit stencil B non-tridiagonal corrections
    real(kind=rk), dimension(:)  , allocatable, save   :: fbnorm         ! "
    integer      , dimension(:)  , allocatable, save   :: ipiv           ! left implicit stencil  B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(dl)) then
       ! Init. right explicit stencil A
       call init_right_explicit_filter_part(stencil)
       ! Init. left implicit stencil B
       call init_left_implicit_filter_part(dl,dc,du,du2,fbrat,fbnorm,ipiv,size(q,1))
    end if

    call filter_main(q,stencil,dl,dc,du,du2,fbrat,fbnorm,ipiv)

  end subroutine filter_i3
!!!==================================================================================================


!!!=================================================================================================
  subroutine filter_main(q,stencil,dl,dc,du,du2,fbrat,fbnorm,ipiv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! solve BU' = AU = F by a BLAS implementation
!!! i.e. implicit filter along first dimension
!!! B is tridiagonal
!!! stencil values are provided by GaitondeVisbal2000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(:,:), intent(inout)       :: q
    real(kind=rk) , dimension(:)  , intent(in)          :: stencil                                 ! right explicit stencil A
    real(kind=rk) , dimension(:)  , intent(in)          :: dl,dc,du,du2                            ! left implicit stencil  B
    real(kind=rk) ,                 intent(in)          :: fbrat                                   ! left implicit stencil  B non-tridiagonal corrections
    real(kind=rk) , dimension(:)  , intent(in)          :: fbnorm                                  ! "
    integer       , dimension(:)  , intent(in)          :: ipiv                                    ! left implicit stencil  B 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(:,:), allocatable         :: q_tmp
    real(kind=rk)                                       :: a0,a1
    integer                                             :: n1,n2
    integer                                             :: i,j,info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n1 = size(q,1)
    n2 = size(q,2)

    allocate(q_tmp(n1,n2))

    ! unpack stencil for readability
    a0  = stencil(1)
    a1  = stencil(2)

    ! compute A*U (interior)
    do j = 1,n2
       do i = 2,(n1-1)
          q_tmp(i,j) = a0*q(i,j) + a1*(q(i+1,j) + q(i-1,j)) ! the factor 1/2 is already considered in stencil
       end do
    end do

    ! compute A*U (boundary)
    q_tmp( 1,:) = a0*q(1,:) + a1*(q(2,:) + q(n1,:))

    q_tmp(n1  ,:) = a0*q(n1  ,:) + a1*(q( 1  ,:) + q(n1-1,:))

    ! write back
    q = q_tmp

    ! solve BU'
    call dgttrs('no transpose', n1, n2, dl, dc, du, du2, ipiv, q, n1, info)
    if(info.ne.0) then
       write(*,*) "error in filter_pen_core_o6_i_c_p.f90:filter_main:error with id: ",info
    end if

    ! correct non-tridiagonal part of B
    do i = 1,n2
       q(:,i) = q(:,i) - fbnorm*(q(1,i) + fbrat*q(n1,i))
    end do

  end subroutine filter_main
!!!==================================================================================================



!!!==================================================================================================
  subroutine init_right_explicit_filter_part(stencil) ! output args
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Init. right explicit stencil A for o2_i_c_p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(2)    , intent(out)       :: stencil
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                       :: a0,a1, alpha
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! get free parameter alpha
    ! choice of alpha   : alpha = 0     --> explicit filter
    !                   : alpha = 0.5   --> no filtering
    !                   : alpha = 0.35  --> default value
    call get_parameter(alpha,'filter.alpha',default = 0.35_rk)

    ! fixed stencil, see table 1
    a0       = 1.0_rk/2.0_rk +  alpha

    a1       = 1.0_rk/2.0_rk + alpha
    a1       = a1/2.0_rk

    stencil  = (/ a0, a1 /)

  end subroutine init_right_explicit_filter_part
!!!==================================================================================================

!!!==================================================================================================
  subroutine init_left_implicit_filter_part(dl,dc,du,du2,fbrat,fbnorm,ipiv,n1) ! output args
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

    ! get free parameter alpha
    ! choice of alpha   : alpha = 0     --> explicit filter
    !                   : alpha = 0.5   --> no filtering
    !                   : alpha = 0.35  --> default value
    call get_parameter(alpha,'filter.alpha',default = 0.35_rk)

    ! allocates
    allocate(   dl(n1-1))
    allocate(    dc(n1 ))
    allocate(   du(n1-1))
    allocate(  du2(n1-2))
    allocate( ipiv(n1  ))
    allocate(   fa(n1  ))
    allocate(   fb(n1  ))
    allocate(fbnorm(n1 ))

    ! pre-stuff for matrix B
    ! fill B (central)
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
    dc(1)    = dc(1)  + dc(1)    ! setup non-tridiagonal B part as vectors (central)

    ! init LU decomposition by dgttrf
    call dgttrf(n1, dl, dc, du, du2, ipiv, info)
    if(info.ne.0) then
       write(*,*) "error in filter_pen_core_o6_i_c_p.f90:init_left_implicit_diff_part:error with id: ",info
    end if

    ! solve Btri y = fa
    call dgttrs('no transpose', n1,  1, dl, dc, du, du2, ipiv, fa, n1, info)
    if(info.ne.0) then
       write(*,*) "error in filter_pen_core_o6_i_c_p.f90:init_left_implicit_diff_part:error with id: ",info
    end if

    ! init prefactors for normalization
    fbnorm = fa*fb(1)/(1._rk+sum(fb*fa))
    fbrat  = fb(n1)/fb(1)

  end subroutine init_left_implicit_filter_part
!!!==================================================================================================

end module filter_pen_core_o2_i_c_p
