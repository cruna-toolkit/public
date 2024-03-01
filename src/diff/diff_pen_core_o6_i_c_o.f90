module diff_pen_core_o6_i_c_o

  use parameter

  private

  public :: diff_i1, diff_i2, diff_i3

  character(len=max_length_parameter), parameter, public :: diff_type = 'open'
  character(len=max_length_parameter), parameter, public :: diff_name = 'diff_pen_core_o6_i_c_o'

  real(kind=rk), dimension(6)        , save              :: stencil      ! right explicit stencil A

contains

!!!=================================================================================================
  subroutine diff_i1(q,dx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)      :: q
    real(kind=rk)                , intent(in)         :: dx
    real(kind=rk), dimension(:) , allocatable, save   :: dl,dc,du,du2   ! left implicit stencil  B
    integer      , dimension(:) , allocatable, save   :: ipiv           ! left implicit stencil  B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(dl)) then
       ! Init. right explicit stencil A
       call init_right_explicit_diff_part(stencil)
       ! Init. left implicit stencil B
       call init_left_implicit_diff_part(dl,dc,du,du2,ipiv,size(q,1))
    end if

    call diff_main(q,dx,stencil,dl,dc,du,du2,ipiv)

  end subroutine diff_i1
!!!==================================================================================================

!!!=================================================================================================
  subroutine diff_i2(q,dx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
    real(kind=rk)                , intent(in)          :: dx
    real(kind=rk) , dimension(:) , allocatable, save   :: dl,dc,du,du2   ! left implicit stencil  B
    integer       , dimension(:) , allocatable, save   :: ipiv           ! left implicit stencil  B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(dl)) then
       ! Init. right explicit stencil A
       call init_right_explicit_diff_part(stencil)
       ! Init. left implicit stencil B
       call init_left_implicit_diff_part(dl,dc,du,du2,ipiv,size(q,1))
    end if

    call diff_main(q,dx,stencil,dl,dc,du,du2,ipiv)

  end subroutine diff_i2
!!!==================================================================================================

!!!=================================================================================================
  subroutine diff_i3(q,dx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
    real(kind=rk)                , intent(in)          :: dx
    real(kind=rk) , dimension(:) , allocatable, save   :: dl,dc,du,du2   ! left implicit stencil  B
    integer       , dimension(:) , allocatable, save   :: ipiv           ! left implicit stencil  B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(dl)) then
       ! Init. right explicit stencil A
       call init_right_explicit_diff_part(stencil)
       ! Init. left implicit stencil B
       call init_left_implicit_diff_part(dl,dc,du,du2,ipiv,size(q,1))
    end if

    call diff_main(q,dx,stencil,dl,dc,du,du2,ipiv)

  end subroutine diff_i3
!!!==================================================================================================


!!!=================================================================================================
  subroutine diff_main(q,dx,stencil,dl,dc,du,du2,ipiv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! solve BU' = AU = F by a BLAS implementation
!!! i.e. implicit derivation along first dimension
!!! B is tridiagonal
!!! U(') are matrices
!!! stencil values are provided by Lele1992
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(:,:), intent(inout)       :: q
    real(kind=rk)                 , intent(in)          :: dx
    real(kind=rk) , dimension(:)  , intent(in)          :: stencil        ! right explicit stencil A
    real(kind=rk) , dimension(:)  , intent(in)          :: dl,dc,du,du2   ! left implicit stencil  B
    integer       , dimension(:)  , intent(in)          :: ipiv           ! left implicit stencil  B 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(:,:), allocatable         :: q_tmp
    real(kind=rk)                                       :: a,b,c,d,e,f
    integer                                             :: n1,n2
    integer                                             :: i,j,info

    n1 = size(q,1)
    n2 = size(q,2)

    allocate(q_tmp(n1,n2))

    a  = stencil(1)
    b  = stencil(2)
    c  = stencil(3)
    d  = stencil(4)
    e  = stencil(5)
    f  = stencil(6)

    ! compute A*U (central)
    ! inner:    sixth order (incl. implicit B matrix)
    ! boundary: forth order (incl. implicit B matrix)
    q_tmp(1   ,:) = -c*q(  1,:) + d*(q(2   ,:)  +    q(   3,:)) - e*q(   4,:)
    q_tmp(2   ,:) = -f*q(  1,:) + f* q(3   ,:)
    do j = 1,n2
       do i = 3,(n1-2)
          q_tmp(i,j) = a*(q(i+1,j) - q(i -1,j)) + b*(q(i +2,j) - q(i -2,j))
       end do
    end do
    q_tmp(n1-1,:) =  f*q(n1 ,:) - f* q(n1-2,:)
    q_tmp(n1  ,:) =  c*q(n1 ,:) - d*(q(n1-1,:)  +    q(n1-2,:)) + e*q(n1-3,:)

    ! dx factor
    q = q_tmp/dx

    ! solve BU'
    call dgttrs('no transpose', n1, n2, dl, dc, du, du2, ipiv, q, n1, info)
    if(info.ne.0) then
       write(*,*) "error in diff_pen_core_o6_i_c_o.f90:diff_main:error with id: ",info
    end if

  end subroutine diff_main
!!!==================================================================================================



!!!==================================================================================================
  subroutine init_right_explicit_diff_part(stencil) ! output args
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Init. right explicit stencil A for o6_i_c_o
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(6)  , intent(out)         :: stencil
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                       :: a,b,c,d,e,f

    a       = 14.00_rk / 9.0_rk / 2.0_rk
    b       =  1.00_rk / 9.0_rk / 4.0_rk
    c       = 17.00_rk / 6.0_rk
    d       =  1.50_rk
    e       =  1.00_rk / 6.0_rk
    f       =  0.75_rk
    stencil = (/ a, b, c, d, e, f /)

  end subroutine init_right_explicit_diff_part
!!!==================================================================================================



!!!==================================================================================================
  subroutine init_left_implicit_diff_part(dl,dc,du,du2,ipiv,n1) ! output args
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! LU decomposition for dgttrf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:), allocatable, intent(out) :: dl,dc,du,du2
    integer      , dimension(:), allocatable, intent(out) :: ipiv
    integer                                 , intent(in)  :: n1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                         :: alpha
    integer                                               :: info

    ! allocates
    allocate(   dl(n1-1 ))
    allocate(    dc(n1   ))
    allocate(   du(n1-1 ))
    allocate(  du2(n1-2 ))
    allocate( ipiv(n1   ))

    ! pre-stuff for matrix B
    ! fill B (central)
    alpha    = 1.0_rk  / 3.0_rk
    dl       = alpha
    dc       = 1.0_rk
    du       = alpha
    ! fill B (boundaries)
    du(1)    = 3.00_rk
    du(2)    = 0.25_rk
    du(n1-1) = 0.25_rk
    dl(1)    = 0.25_rk
    dl(n1-1) = 3.00_rk
    dl(n1-2) = 0.25_rk

    call dgttrf(n1, dl, dc, du, du2, ipiv, info)
    if(info.ne.0) then
       write(*,*) "error in diff_pen_core_o6_i_c_o.f90:init_left_implicit_diff_part:error with id: ",info
    end if

  end subroutine init_left_implicit_diff_part
!!!==================================================================================================


end module diff_pen_core_o6_i_c_o
