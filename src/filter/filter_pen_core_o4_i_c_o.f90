module filter_pen_core_o4_i_c_o

  use parameter

  private

  public :: filter_i1, filter_i2, filter_i3

  character(len=max_length_parameter), parameter, public :: filter_type = 'open'
  character(len=max_length_parameter), parameter, public :: filter_name = 'filter_pen_core_o4_i_c_o'

  real(kind=rk), dimension(3)        , save              :: stencil      ! right explicit stencil A
  real(kind=rk), dimension(1,5)      , save              :: coefficients_for_boundary

contains

!!!=================================================================================================
  subroutine filter_i1(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
    real(kind=rk), dimension(:)  , allocatable, save   :: dl,dc,du,du2   ! left implicit stencil  B
    integer      , dimension(:)  , allocatable, save   :: ipiv           ! left implicit stencil  B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(dl)) then
       ! Init. right explicit stencil A
       call init_right_explicit_filter_part(stencil,coefficients_for_boundary)
       ! Init. left implicit stencil B
       call init_left_implicit_filter_part(dl,dc,du,du2,ipiv,size(q,1))
    end if

    call filter_main(q,stencil,coefficients_for_boundary,dl,dc,du,du2,ipiv)

  end subroutine filter_i1
!!!==================================================================================================

!!!=================================================================================================
  subroutine filter_i2(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
    real(kind=rk), dimension(:)  , allocatable, save   :: dl,dc,du,du2   ! left implicit stencil  B
    integer      , dimension(:)  , allocatable, save   :: ipiv           ! left implicit stencil  B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(dl)) then
       ! Init. right explicit stencil A
       call init_right_explicit_filter_part(stencil,coefficients_for_boundary)
       ! Init. left implicit stencil B
       call init_left_implicit_filter_part(dl,dc,du,du2,ipiv,size(q,1))
    end if

    call filter_main(q,stencil,coefficients_for_boundary,dl,dc,du,du2,ipiv)

  end subroutine filter_i2
!!!==================================================================================================

!!!=================================================================================================
  subroutine filter_i3(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(inout)       :: q
    real(kind=rk), dimension(:)  , allocatable, save   :: dl,dc,du,du2   ! left implicit stencil  B
    integer      , dimension(:)  , allocatable, save   :: ipiv           ! left implicit stencil  B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(dl)) then
       ! Init. right explicit stencil A
       call init_right_explicit_filter_part(stencil,coefficients_for_boundary)
       ! Init. left implicit stencil B
       call init_left_implicit_filter_part(dl,dc,du,du2,ipiv,size(q,1))
    end if

    call filter_main(q,stencil,coefficients_for_boundary,dl,dc,du,du2,ipiv)

  end subroutine filter_i3
!!!==================================================================================================


!!!=================================================================================================
  subroutine filter_main(q,stencil,coefficients_for_boundary,dl,dc,du,du2,ipiv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! solve BU' = AU = F by a BLAS implementation
!!! i.e. implicit filter along first dimension
!!! B is tridiagonal
!!! stencil values are provided by GaitondeVisbal2000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(:,:), intent(inout)       :: q
    real(kind=rk) , dimension(:,:), intent(in)          :: coefficients_for_boundary               ! right explicit stencil A (boundary points)
    real(kind=rk) , dimension(:)  , intent(in)          :: stencil                                 ! right explicit stencil A
    real(kind=rk) , dimension(:)  , intent(in)          :: dl,dc,du,du2                            ! left implicit stencil  B
    integer       , dimension(:)  , intent(in)          :: ipiv                                    ! left implicit stencil  B 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(:,:), allocatable         :: q_tmp
    real(kind=rk)                                       :: a0,a1,a2,a3
    integer                                             :: n1,n2
    integer                                             :: i,j,info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    n1 = size(q,1)
    n2 = size(q,2)

    allocate(q_tmp(n1,n2))

    ! unpack stencil for readability
    a0  = stencil(1)
    a1  = stencil(2)
    a2  = stencil(3)

    ! compute A*U (interior)
    do j = 1,n2
       do i = 3,(n1-2)
          q_tmp(i,j) = a0*q(i,j) + a1*(q(i+1,j) + q(i-1,j)) + a2*(q(i+2,j) + q(i-2,j)) ! the factor 1/2 is already considered in stencil
       end do
    end do

    ! compute A*U (boundary)
    q_tmp( 1,:) = q( 1,:) ! do nothing at the boundary points
    q_tmp(n1,:) = q(n1,:)

    q_tmp(   2,:) = 0.0_rk
    q_tmp(n1-1,:) = 0.0_rk
    do i = 1,5                ! point 2
       q_tmp(   2,:) = q_tmp(   2,:) + coefficients_for_boundary(1,i)*q(i         ,:)
       q_tmp(n1-1,:) = q_tmp(n1-1,:) + coefficients_for_boundary(1,i)*q(n1 - i + 1,:)
    end do


    ! write back
    q = q_tmp

    ! solve BU'
    call dgttrs('no transpose', n1, n2, dl, dc, du, du2, ipiv, q, n1, info)
    if(info.ne.0) then
       write(*,*) "error in filter_pen_core_o6_i_c_o.f90:filter_main:error with id: ",info
    end if

  end subroutine filter_main
!!!==================================================================================================



!!!==================================================================================================
  subroutine init_right_explicit_filter_part(stencil,coefficients_for_boundary) ! output args
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Init. right explicit stencil A for o4_i_c_o
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(1,5)  , intent(out)       :: coefficients_for_boundary
    real(kind=rk) , dimension(3)    , intent(out)       :: stencil
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                       :: a0,a1,a2,a3, alpha
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! get free parameter alpha
    ! choice of alpha   : alpha = 0     --> explicit filter
    !                   : alpha = 0.5   --> no filtering
    !                   : alpha = 0.35  --> default value
    call get_parameter(alpha,'filter.alpha',default = 0.35_rk)

    ! fixed stencil, see table 1
    a0       = 5.0_rk/8.0_rk +  3.0_rk*alpha/ 4.0_rk

    a1       = 1.0_rk/2.0_rk + alpha
    a1       = a1/2.0_rk

    a2       = -1.0_rk/8.0_rk +  alpha/ 4.0_rk
    a2       = a2/2.0_rk


    stencil  = (/ a0, a1, a2 /)

    ! point 2 (non-per case)
    coefficients_for_boundary(1,1) =  1.0_rk/16.0_rk + 7.0_rk*alpha/8.0_rk
    coefficients_for_boundary(1,2) =  3.0_rk/4.0_rk  + 1.0_rk*alpha/2.0_rk
    coefficients_for_boundary(1,3) =  3.0_rk/8.0_rk  + 1.0_rk*alpha/4.0_rk
    coefficients_for_boundary(1,4) = -1.0_rk/4.0_rk  + 1.0_rk*alpha/2.0_rk
    coefficients_for_boundary(1,5) =  1.0_rk/16.0_rk - 1.0_rk*alpha/8.0_rk

  end subroutine init_right_explicit_filter_part
!!!==================================================================================================

!!!==================================================================================================
  subroutine init_left_implicit_filter_part(dl,dc,du,du2,ipiv,n1) ! output args
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! LU decomposition for dgttrf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:), allocatable, intent(out) :: dl,dc,du,du2
    integer      , dimension(:), allocatable, intent(out) :: ipiv
    integer                                 , intent(in)  :: n1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

    ! pre-stuff for matrix B
    ! fill B (central)
    dl       = alpha
    dc       = 1.0_rk
    du       = alpha
    ! fill B (boundaries are unchanged, not filtered)
    du(1)    = 0.0_rk  
    dl(n1-1) = 0.0_rk

    call dgttrf(n1, dl, dc, du, du2, ipiv, info)
    if(info.ne.0) then
       write(*,*) "error in filter_pen_core_o6_i_c_o.f90:init_left_implicit_filter_part:error with id: ",info
    end if

  end subroutine init_left_implicit_filter_part
!!!==================================================================================================


end module filter_pen_core_o4_i_c_o
