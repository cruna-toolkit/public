module filter_pen_core_o10_i_c_o

  use parameter

  private

  public :: filter_i1, filter_i2, filter_i3

  character(len=max_length_parameter), parameter, public  :: filter_type = 'open'
  character(len=max_length_parameter), parameter, public  :: filter_name = 'filter_pen_core_o10_i_c_o'

  real(kind=rk), dimension(6)        , save               :: stencil      ! right explicit stencil A
  real(kind=rk), dimension(4,11)      , save              :: coefficients_for_boundary

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
    real(kind=rk)                                       :: a0,a1,a2,a3,a4,a5
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
    a3  = stencil(4)
    a4  = stencil(5)
    a5  = stencil(6)

    ! compute A*U (interior)
    do j = 1,n2
       do i = 6,(n1-5)
          q_tmp(i,j) = a0*q(i,j) + a1*(q(i+1,j) + q(i-1,j)) + a2*(q(i+2,j) + q(i-2,j)) + a3*(q(i+3,j) + q(i-3,j)) + a4*(q(i+4,j) + q(i-4,j)) + a5*(q(i+5,j) + q(i-5,j)) ! the factor 1/2 is already considered in stencil
       end do
    end do

    ! compute A*U (boundary)
    q_tmp( 1,:) = q( 1,:) ! do nothing at the boundary points
    q_tmp(n1,:) = q(n1,:)

    q_tmp(   2,:) = 0.0_rk
    q_tmp(n1-1,:) = 0.0_rk
    do i = 1,11                ! point 2
       q_tmp(   2,:) = q_tmp(   2,:) + coefficients_for_boundary(1,i)*q(i         ,:)
       q_tmp(n1-1,:) = q_tmp(n1-1,:) + coefficients_for_boundary(1,i)*q(n1 - i + 1,:)
    end do

    q_tmp(   3,:) = 0.0_rk
    q_tmp(n1-2,:) = 0.0_rk
    do i = 1,11                ! point 3
       q_tmp(   3,:) = q_tmp(   3,:) + coefficients_for_boundary(2,i)*q(i         ,:)
       q_tmp(n1-2,:) = q_tmp(n1-2,:) + coefficients_for_boundary(2,i)*q(n1 - i + 1,:)
    end do
    
    q_tmp(   4,:) = 0.0_rk
    q_tmp(n1-3,:) = 0.0_rk
    do i = 1,11                ! point 4
       q_tmp(   4,:) = q_tmp(   4,:) + coefficients_for_boundary(3,i)*q(i         ,:)
       q_tmp(n1-3,:) = q_tmp(n1-3,:) + coefficients_for_boundary(3,i)*q(n1 - i + 1,:)
    end do
    
    q_tmp(   5,:) = 0.0_rk
    q_tmp(n1-4,:) = 0.0_rk
    do i = 1,11                ! point 5
       q_tmp(   5,:) = q_tmp(   5,:) + coefficients_for_boundary(4,i)*q(i         ,:)
       q_tmp(n1-4,:) = q_tmp(n1-4,:) + coefficients_for_boundary(4,i)*q(n1 - i + 1,:)
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
    ! Init. right explicit stencil A for o6_i_c_o
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk) , dimension(4,11)  , intent(out)       :: coefficients_for_boundary
    real(kind=rk) , dimension(6)    , intent(out)        :: stencil
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                        :: a0,a1,a2,a3,a4,a5, alpha
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! get free parameter alpha
    ! choice of alpha   : alpha = 0     --> explicit filter
    !                   : alpha = 0.5   --> no filtering
    !                   : alpha = 0.35  --> default value
    call get_parameter(alpha,'filter.alpha',default = 0.35_rk)

    ! fixed stencil, see table 1
    a0       = 193.0_rk/256.0_rk +  126.0_rk*alpha/ 256.0_rk

    a1       = 105.0_rk/256.0_rk + 302.0_rk*alpha/256.0_rk
    a1       = a1/2.0_rk

    a2       = -15.0_rk/64.0_rk +  30.0_rk*alpha/ 64.0_rk
    a2       = a2/2.0_rk

    a3       = 45.0_rk/512.0_rk -  90.0_rk*alpha/512.0_rk
    a3       = a3/2.0_rk
    
    a4       = -5.0_rk/256.0_rk + 10.0_rk*alpha/256.0_rk
    a4       = a4/2.0_rk
    
    a5       = 1.0_rk/512.0_rk  - 2.0_rk*alpha/512.0_rk
    a5       = a5/2.0_rk

    stencil  = (/ a0, a1, a2, a3, a4, a5 /)

    ! point 2 (non-per case)
    coefficients_for_boundary(1,1)  =   1.0_rk/1024.0_rk + 1022.0_rk*alpha/1024.0_rk
    coefficients_for_boundary(1,2)  = 507.0_rk/512.0_rk  +  10.0_rk*alpha/512.0_rk
    coefficients_for_boundary(1,3)  =  45.0_rk/1024.0_rk + 934.0_rk*alpha/1024.0_rk
    coefficients_for_boundary(1,4)  = -15.0_rk/128.0_rk  +  30.0_rk*alpha/128.0_rk
    coefficients_for_boundary(1,5)  = 105.0_rk/512.0_rk  - 210.0_rk*alpha/512.0_rk
    coefficients_for_boundary(1,6)  = -63.0_rk/256.0_rk  + 126.0_rk*alpha/256.0_rk
    coefficients_for_boundary(1,7)  =  105.0_rk/512.0_rk - 210.0_rk*alpha/512.0_rk
    coefficients_for_boundary(1,8)  = -15.0_rk/128.0_rk  + 30.0_rk*alpha/128.0_rk
    coefficients_for_boundary(1,9)  =  45.0_rk/1024.0_rk - 90.0_rk*alpha/1024.0_rk
    coefficients_for_boundary(1,10) = -5.0_rk/512.0_rk   + 10.0_rk*alpha/512.0_rk
    coefficients_for_boundary(1,11) =  1.0_rk/1024.0_rk  -  2.0_rk*alpha/1024.0_rk

    ! point 3 (non-per case)
    coefficients_for_boundary(2,1)  = -1.0_rk/1024.0_rk   +  2.0_rk*alpha/1024.0_rk
    coefficients_for_boundary(2,2)  =  5.0_rk/512.0_rk    + 502.0_rk*alpha/512.0_rk
    coefficients_for_boundary(2,3)  = 979.0_rk/1024.0_rk  +  90.0_rk*alpha/1024.0_rk
    coefficients_for_boundary(2,4)  =  15.0_rk/128.0_rk   +  98.0_rk*alpha/128.0_rk
    coefficients_for_boundary(2,5)  = -105.0_rk/512.0_rk  + 210.0_rk*alpha/512.0_rk
    coefficients_for_boundary(2,6)  =  63.0_rk/256.0_rk   - 126.0_rk*alpha/256.0_rk
    coefficients_for_boundary(2,7)  = -105.0_rk/512.0_rk  + 210.0_rk*alpha/512.0_rk
    coefficients_for_boundary(2,8)  =  15.0_rk/128.0_rk   -  30.0_rk*alpha/128.0_rk
    coefficients_for_boundary(2,9)  = -45.0_rk/1024.0_rk  +  90.0_rk*alpha/1024.0_rk
    coefficients_for_boundary(2,10) =  5.0_rk/512.0_rk    -  10.0_rk*alpha/512.0_rk
    coefficients_for_boundary(2,11) = -1.0_rk/1024.0_rk   +  2.0_rk*alpha/1024.0_rk
    
    ! point 4 (non-per case)
    coefficients_for_boundary(3,1)  =  1.0_rk/1024.0_rk  -  2.0_rk*alpha/1024.0_rk
    coefficients_for_boundary(3,2)  = -5.0_rk/512.0_rk   +  10.0_rk*alpha/512.0_rk
    coefficients_for_boundary(3,3)  = 45.0_rk/1024.0_rk  + 934.0_rk*alpha/1024.0_rk
    coefficients_for_boundary(3,4)  = 113.0_rk/128.0_rk  +  30.0_rk*alpha/128.0_rk
    coefficients_for_boundary(3,5)  = 105.0_rk/512.0_rk  + 302.0_rk*alpha/512.0_rk
    coefficients_for_boundary(3,6)  = -63.0_rk/256.0_rk  +  126.0_rk*alpha/256.0_rk
    coefficients_for_boundary(3,7)  =  105.0_rk/512.0_rk - 210.0_rk*alpha/512.0_rk
    coefficients_for_boundary(3,8)  = -15.0_rk/128.0_rk  +  30.0_rk*alpha/128.0_rk
    coefficients_for_boundary(3,9)  =  45.0_rk/1024.0_rk -  90.0_rk*alpha/1024.0_rk
    coefficients_for_boundary(3,10) = -5.0_rk/512.0_rk   +  10.0_rk*alpha/512.0_rk
    coefficients_for_boundary(3,11) =  1.0_rk/1024.0_rk  -  2.0_rk*alpha/1024.0_rk
    
    ! point 5 (non-per case)
    coefficients_for_boundary(4,1)  = -1.0_rk/1024.0_rk   +  2.0_rk*alpha/1024.0_rk
    coefficients_for_boundary(4,2)  =  5.0_rk/512.0_rk    -  10.0_rk*alpha/512.0_rk
    coefficients_for_boundary(4,3)  = -45.0_rk/1024.0_rk  + 90.0_rk*alpha/1024.0_rk
    coefficients_for_boundary(4,4)  =  15.0_rk/128.0_rk   + 98.0_rk*alpha/128.0_rk
    coefficients_for_boundary(4,5)  = 407.0_rk/512.0_rk   + 210.0_rk*alpha/512.0_rk
    coefficients_for_boundary(4,6)  =  63.0_rk/256.0_rk   +  130.0_rk*alpha/256.0_rk
    coefficients_for_boundary(4,7)  = -105.0_rk/512.0_rk  + 210.0_rk*alpha/512.0_rk
    coefficients_for_boundary(4,8)  =  15.0_rk/128.0_rk   - 30.0_rk*alpha/128.0_rk
    coefficients_for_boundary(4,9)  = -45.0_rk/1024.0_rk  + 90.0_rk*alpha/1024.0_rk
    coefficients_for_boundary(4,10) =  5.0_rk/512.0_rk    - 10.0_rk*alpha/512.0_rk
    coefficients_for_boundary(4,11) = -1.0_rk/1024.0_rk   +  2.0_rk*alpha/1024.0_rk

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


end module filter_pen_core_o10_i_c_o
