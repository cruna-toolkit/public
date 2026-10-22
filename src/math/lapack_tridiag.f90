module lapack_tridiag

  USE iso_fortran_env , only: real32, real64

  implicit none

  private

  public :: gttrf, gttrs


  interface gttrf
     module procedure gttrf_32, gttrf_64
  end interface gttrf

  interface gttrs
     module procedure gttrs_32_1d, gttrs_32_2d
     module procedure gttrs_64_1d, gttrs_64_2d
  end interface gttrs


  interface

     subroutine sgttrf(n,dl,d,du,du2,ipiv,info)
       import :: real32
       integer                                , intent(in)    :: n
       real(kind=real32), dimension(*)        , intent(inout) :: dl,d,du
       real(kind=real32), dimension(*)        , intent(out)   :: du2
       integer          , dimension(*)        , intent(out)   :: ipiv
       integer                                , intent(out)   :: info
     end subroutine sgttrf

     subroutine dgttrf(n,dl,d,du,du2,ipiv,info)
       import :: real64
       integer                                , intent(in)    :: n
       real(kind=real64), dimension(*)        , intent(inout) :: dl,d,du
       real(kind=real64), dimension(*)        , intent(out)   :: du2
       integer          , dimension(*)        , intent(out)   :: ipiv
       integer                                , intent(out)   :: info
     end subroutine dgttrf


     subroutine sgttrs(trans,n,nrhs,dl,d,du,du2,ipiv,b,ldb,info)
       import :: real32
       character(len=1)                       , intent(in)    :: trans
       integer                                , intent(in)    :: n,nrhs,ldb
       real(kind=real32), dimension(*)        , intent(in)    :: dl,d,du,du2
       integer          , dimension(*)        , intent(in)    :: ipiv
       real(kind=real32), dimension(ldb,*)    , intent(inout) :: b
       integer                                , intent(out)   :: info
     end subroutine sgttrs

     subroutine dgttrs(trans,n,nrhs,dl,d,du,du2,ipiv,b,ldb,info)
       import :: real64
       character(len=1)                       , intent(in)    :: trans
       integer                                , intent(in)    :: n,nrhs,ldb
       real(kind=real64), dimension(*)        , intent(in)    :: dl,d,du,du2
       integer          , dimension(*)        , intent(in)    :: ipiv
       real(kind=real64), dimension(ldb,*)    , intent(inout) :: b
       integer                                , intent(out)   :: info
     end subroutine dgttrs

  end interface


contains


!!!=================================================================================================
  subroutine gttrf_32(n,dl,d,du,du2,ipiv,info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                , intent(in)      :: n
    real(kind=real32), dimension(:)        , intent(inout)   :: dl,d,du
    real(kind=real32), dimension(:)        , intent(out)     :: du2
    integer          , dimension(:)        , intent(out)     :: ipiv
    integer                                , intent(out)     :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call sgttrf(n,dl,d,du,du2,ipiv,info)

  end subroutine gttrf_32
!!!==================================================================================================


!!!=================================================================================================
  subroutine gttrf_64(n,dl,d,du,du2,ipiv,info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                , intent(in)      :: n
    real(kind=real64), dimension(:)        , intent(inout)   :: dl,d,du
    real(kind=real64), dimension(:)        , intent(out)     :: du2
    integer          , dimension(:)        , intent(out)     :: ipiv
    integer                                , intent(out)     :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call dgttrf(n,dl,d,du,du2,ipiv,info)

  end subroutine gttrf_64
!!!==================================================================================================


!!!=================================================================================================
  subroutine gttrs_32_1d(trans,n,nrhs,dl,d,du,du2,ipiv,b,ldb,info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=1)                       , intent(in)      :: trans
    integer                                , intent(in)      :: n,nrhs,ldb
    real(kind=real32), dimension(:)        , intent(in)      :: dl,d,du,du2
    integer          , dimension(:)        , intent(in)      :: ipiv
    real(kind=real32), dimension(:), target, intent(inout)   :: b
    integer                                , intent(out)     :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=real32), dimension(:,:), pointer               :: b_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (nrhs.ne.1) then
       info = -3
       return
    end if

    b_2d(1:ldb,1:1) => b(1:ldb)

    call sgttrs(trans,n,nrhs,dl,d,du,du2,ipiv,b_2d,ldb,info)

  end subroutine gttrs_32_1d
!!!==================================================================================================


!!!=================================================================================================
  subroutine gttrs_32_2d(trans,n,nrhs,dl,d,du,du2,ipiv,b,ldb,info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=1)                       , intent(in)      :: trans
    integer                                , intent(in)      :: n,nrhs,ldb
    real(kind=real32), dimension(:)        , intent(in)      :: dl,d,du,du2
    integer          , dimension(:)        , intent(in)      :: ipiv
    real(kind=real32), dimension(:,:)      , intent(inout)   :: b
    integer                                , intent(out)     :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call sgttrs(trans,n,nrhs,dl,d,du,du2,ipiv,b,ldb,info)

  end subroutine gttrs_32_2d
!!!==================================================================================================


!!!=================================================================================================
  subroutine gttrs_64_1d(trans,n,nrhs,dl,d,du,du2,ipiv,b,ldb,info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=1)                       , intent(in)      :: trans
    integer                                , intent(in)      :: n,nrhs,ldb
    real(kind=real64), dimension(:)        , intent(in)      :: dl,d,du,du2
    integer          , dimension(:)        , intent(in)      :: ipiv
    real(kind=real64), dimension(:), target, intent(inout)   :: b
    integer                                , intent(out)     :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=real64), dimension(:,:), pointer               :: b_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (nrhs.ne.1) then
       info = -3
       return
    end if

    b_2d(1:ldb,1:1) => b(1:ldb)

    call dgttrs(trans,n,nrhs,dl,d,du,du2,ipiv,b_2d,ldb,info)

  end subroutine gttrs_64_1d
!!!==================================================================================================


!!!=================================================================================================
  subroutine gttrs_64_2d(trans,n,nrhs,dl,d,du,du2,ipiv,b,ldb,info)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=1)                       , intent(in)      :: trans
    integer                                , intent(in)      :: n,nrhs,ldb
    real(kind=real64), dimension(:)        , intent(in)      :: dl,d,du,du2
    integer          , dimension(:)        , intent(in)      :: ipiv
    real(kind=real64), dimension(:,:)      , intent(inout)   :: b
    integer                                , intent(out)     :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call dgttrs(trans,n,nrhs,dl,d,du,du2,ipiv,b,ldb,info)

  end subroutine gttrs_64_2d
!!!==================================================================================================


end module lapack_tridiag
