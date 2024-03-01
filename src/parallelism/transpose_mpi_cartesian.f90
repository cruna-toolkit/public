module transpose_mpi_cartesian

  USE mpi
  use parameter

  private

  public :: ftranspose_i1 , btranspose_i1 , ftranspose_i2 , btranspose_i2 , ftranspose_i3 , btranspose_i3

  interface ftranspose_i1
     !module procedure ftranspose_i1_3d_rk
     module procedure ftranspose_i1_wrapper
  end interface ftranspose_i1

  interface btranspose_i1
     !module procedure btranspose_i1_3d_rk
     module procedure btranspose_i1_wrapper
  end interface btranspose_i1

  interface ftranspose_i2
     !module procedure ftranspose_i2_3d_rk
     module procedure ftranspose_i2_wrapper
  end interface ftranspose_i2

  interface btranspose_i2
     !module procedure btranspose_i2_3d_rk
     module procedure btranspose_i2_wrapper
  end interface btranspose_i2

  interface ftranspose_i3
     module procedure ftranspose_i3_wrapper
     !module procedure ftranspose_i3_3d_rk
  end interface ftranspose_i3

  interface btranspose_i3
     module procedure btranspose_i3_wrapper
     !module procedure btranspose_i3_3d_rk
  end interface btranspose_i3

contains

!!!==================================================================================================
  subroutine ftranspose_i1_3d_rk(fries_2d,field_3d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:)  , intent(out), allocatable :: fries_2d
    real(kind=rk), dimension(:,:,:), intent(in)               :: field_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, save, dimension(:,:)  , allocatable   :: pen_fries_mapping
    integer, save, dimension(:)    , allocatable   :: send_counts, recv_counts
    integer, save, dimension(:)    , allocatable   :: send_displs, recv_displs
    integer, save                                  :: pen_fries_nbr, pen_len, pen_comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:),save,allocatable   :: field_vector, field_vector_t
    integer                                        :: i,pen_fries_nbr_times_pen_len
    integer                                        :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(send_counts)) then
       ! here the resort tupples and some dimension information for mpi_alltoallv are constructed ONLY once
       ! this is different for different pencils
       ! save     variables are reused in each transpose call
       ! not save variables are just needed once for creation of the resort tupples and mpi_alltoallv args
       call init_i1(send_counts,recv_counts,send_displs,recv_displs,pen_fries_mapping,pen_fries_nbr,pen_len,pen_comm)

       ! allocate temporary vectors
       allocate(field_vector(  sum(send_counts)))
       allocate(field_vector_t(sum(recv_counts)))

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! alltoallv super magic (forth) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! pack field vector for i1 pencils and fries, different for different pens
    !write(*,*) size(field_3d,1),size(field_3d,2),size(field_3d,3)
    !write(*,*) size(field_vector),size(field_3d)
    field_vector = reshape(field_3d,shape(field_vector))

    ! l = 1
    ! ! loop different fries i1-i2-i3-order
    ! do k = 1,size(field_3d,3)
    !    do j = 1,size(field_3d,2)
    !       ! fries lengthwise
    !       do i = 1,size(field_3d,1)
    !          field_vector(l) = field_3d(i,j,k)
    !          l = l + 1
    !       end do
    !    end do
    ! end do

    ! transpose (the magic)
    call mpi_alltoallv(field_vector, send_counts, send_displs, params%parallelism%mpi_type_rk, &   ! send args
         field_vector_t, recv_counts, recv_displs, params%parallelism%mpi_type_rk, &               ! recv args
         pen_comm, ierr)                                                                           ! mpi stuff

    ! allocate temporary array
    allocate(fries_2d(pen_len,pen_fries_nbr))

    ! unpack transposed vector via pen fries mapping tupples
    pen_fries_nbr_times_pen_len = pen_fries_nbr*pen_len
    do i = 1,pen_fries_nbr_times_pen_len
       fries_2d(pen_fries_mapping(1,i),pen_fries_mapping(2,i)) = field_vector_t(i)
    end do

  end subroutine ftranspose_i1_3d_rk
!!!==================================================================================================

!!!==================================================================================================
  subroutine f0transpose_i1_3d_rk(fries_2d,field_3d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:)  , intent(out), allocatable :: fries_2d
    real(kind=rk), dimension(:,:,:), intent(in)               :: field_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !integer                                                   :: j,k,l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(fries_2d(size(field_3d,1),size(field_3d,2)*size(field_3d,3)))
    ! l = 0
    ! do k = 1,size(field_3d,3)
    !    do j = 1,size(field_3d,2)
    !       l = l + 1
    !       fries_2d(:,l) = field_3d(:,j,k)
    !    end do
    ! end do
    fries_2d = reshape(field_3d,shape(fries_2d))

  end subroutine f0transpose_i1_3d_rk
!!!==================================================================================================

!!!==================================================================================================
  subroutine ftranspose_i1_wrapper(fries_2d,field_3d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:)  , intent(out), allocatable :: fries_2d
    real(kind=rk), dimension(:,:,:), intent(in)               :: field_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%i1.eq.1) then
       call f0transpose_i1_3d_rk(fries_2d,field_3d)
    else
       call ftranspose_i1_3d_rk(fries_2d,field_3d)
    end if

  end subroutine ftranspose_i1_wrapper
!!!==================================================================================================

!!!==================================================================================================
  subroutine btranspose_i1_3d_rk(field_3d,fries_2d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(out)          :: field_3d
    real(kind=rk), dimension(:,:)  , intent(in)           :: fries_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, save, dimension(:,:)  , allocatable          :: pen_fries_mapping
    integer, save, dimension(:)    , allocatable          :: send_counts, recv_counts
    integer, save, dimension(:)    , allocatable          :: send_displs, recv_displs
    integer, save                                         :: pen_fries_nbr, pen_len, pen_comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:),save,allocatable          :: field_vector, field_vector_t
    integer                                               :: i,pen_fries_nbr_times_pen_len
    integer                                               :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(send_counts)) then
       ! here the resort tupples and some dimension information for mpi_alltoallv are constructed ONLY once
       ! this is different for different pencils
       !     save variables are reused in each transpose call
       ! not save variables are just needed once for creation of the resort tupples and mpi_alltoallv args
       call init_i1(send_counts,recv_counts,send_displs,recv_displs,pen_fries_mapping,pen_fries_nbr,pen_len,pen_comm)

       ! allocate temporary vectors
       allocate(field_vector(  sum(send_counts)))
       allocate(field_vector_t(sum(recv_counts)))
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! alltoallv super magic (back) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! pack transposed vector via pen fries mapping tupples
    pen_fries_nbr_times_pen_len = pen_fries_nbr*pen_len
    do i = 1,pen_fries_nbr_times_pen_len
       field_vector_t(i) = fries_2d(pen_fries_mapping(1,i),pen_fries_mapping(2,i))
    end do

    ! transpose (the magic)
    call mpi_alltoallv(field_vector_t, recv_counts, recv_displs, params%parallelism%mpi_type_rk, &  ! send args
         field_vector, send_counts, send_displs, params%parallelism%mpi_type_rk, &                  ! recv args
         pen_comm, ierr)

    ! unpack vector for i1 pencils and fries, different for different pens
    field_3d = reshape(field_vector, shape(field_3d))

    ! l = 1
    ! ! loop different fries i1-i2-i3-order
    ! do k = 1,size(field_3d,3)
    !    do j = 1,size(field_3d,2)
    !       ! frie lengthwise
    !       do i = 1,size(field_3d,1)
    !          field_3d(i,j,k) = field_vector(l)
    !          l = l + 1
    !       end do
    !    end do
    ! end do

  end subroutine btranspose_i1_3d_rk
!!!==================================================================================================

!!!==================================================================================================
  subroutine b0transpose_i1_3d_rk(field_3d,fries_2d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(out)          :: field_3d
    real(kind=rk), dimension(:,:)  , intent(in)           :: fries_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    integer                                               :: j,k,l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! l = 0
    ! do k = 1,size(field_3d,3)
    !    do j = 1,size(field_3d,2)
    !       l = l + 1
    !       field_3d(:,j,k) = fries_2d(:,l)
    !    end do
    ! end do
    field_3d = reshape(fries_2d, shape(field_3d))

  end subroutine b0transpose_i1_3d_rk
!!!==================================================================================================

!!!==================================================================================================
  subroutine btranspose_i1_wrapper(field_3d,fries_2d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(out)          :: field_3d
    real(kind=rk), dimension(:,:)  , intent(in)           :: fries_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%i1.eq.1) then
       call b0transpose_i1_3d_rk(field_3d,fries_2d)
    else
       call btranspose_i1_3d_rk(field_3d,fries_2d)
    end if

  end subroutine btranspose_i1_wrapper
!!!==================================================================================================

!!!==================================================================================================
  subroutine ftranspose_i2_3d_rk(fries_2d,field_3d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:)  , intent(out), allocatable :: fries_2d
    real(kind=rk), dimension(:,:,:), intent(in)               :: field_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, save, dimension(:,:)  , allocatable              :: pen_fries_mapping
    integer, save, dimension(:)    , allocatable              :: send_counts, recv_counts
    integer, save, dimension(:)    , allocatable              :: send_displs, recv_displs
    integer, save                                             :: pen_fries_nbr, pen_len, pen_comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !real(kind=rk), dimension(:,:)   ,allocatable              :: field_2d
    real(kind=rk), dimension(:),save,allocatable              :: field_vector, field_vector_t
    integer                                                   :: s1,s2,s3
    integer                                                   :: i,pen_fries_nbr_times_pen_len
    integer                                                   :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(send_counts)) then
       ! here the resort tupples and some dimension information for mpi_alltoallv are constructed ONLY once
       ! this is different for different pencils
       ! save     variables are reused in each transpose call
       ! not save variables are just needed once for creation of the resort tupples and mpi_alltoallv args
       call init_i2(send_counts,recv_counts,send_displs,recv_displs,pen_fries_mapping,pen_fries_nbr,pen_len,pen_comm)

       ! allocate temporary vectors
       allocate(field_vector(  sum(send_counts)))
       allocate(field_vector_t(sum(recv_counts)))
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! alltoallv super magic (forth) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    s1 = size(field_3d,1)
    s2 = size(field_3d,2)
    s3 = size(field_3d,3)

    ! allocate(field_2d(s1,s2))
    ! l = 1
    ! ! loop different fries i1-i2-i3-order
    ! do k = 1,s3
    !    field_2d = field_3d(:,:,k)
    !    do i = 1,s1
    !       ! fries lengthwise
    !       field_vector(l:l+s2-1) = field_2d(i,:)
    !       l = l + s2
    !    end do
    ! end do
    
    ! pack field vector for i1 pencils and fries, different for different pens
    field_vector = reshape( reshape(field_3d,(/s2,s1,s3/),order =(/2,1,3/)), shape(field_vector))

    ! l = 1
    ! ! loop different fries i1-i2-i3-order
    ! do k = 1,size(field_3d,3)
    !    do i = 1,size(field_3d,1)
    !       ! fries lengthwise
    !       do j = 1,size(field_3d,2)
    !          field_vector(l) = field_3d(i,j,k)
    !          l = l + 1
    !       end do
    !    end do
    ! end do

    ! transpose (the magic)
    call mpi_alltoallv(field_vector, send_counts, send_displs, params%parallelism%mpi_type_rk, &   ! send args
         field_vector_t, recv_counts, recv_displs, params%parallelism%mpi_type_rk, &               ! recv args
         pen_comm, ierr)  

    ! allocate temporary array
    allocate(fries_2d(pen_len,pen_fries_nbr))                                              ! mpi stuff

    ! unpack transposed vector via pen fries mapping tupples
    pen_fries_nbr_times_pen_len = pen_fries_nbr*pen_len
    do i = 1,pen_fries_nbr_times_pen_len
       fries_2d(pen_fries_mapping(1,i),pen_fries_mapping(2,i)) = field_vector_t(i)
    end do

  end subroutine ftranspose_i2_3d_rk
!!!==================================================================================================

!!!==================================================================================================
  subroutine f0transpose_i2_3d_rk(fries_2d,field_3d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:)  , intent(out), allocatable :: fries_2d
    real(kind=rk), dimension(:,:,:), intent(in)               :: field_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                   :: s1,s2,s3 !i,k,l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    s1 = size(field_3d,1)
    s2 = size(field_3d,2)
    s3 = size(field_3d,3)
    
    allocate(fries_2d(s2,s1*s3))
    
    fries_2d = reshape( reshape(field_3d,(/s2,s1,s3/),order =(/2,1,3/)), shape(fries_2d))

    ! l = 0
    ! do k = 1,size(field_3d,3)
    !    do i = 1,size(field_3d,1)
    !       l = l + 1
    !       fries_2d(:,l) = field_3d(i,:,k)
    !    end do
    ! end do

  end subroutine f0transpose_i2_3d_rk
!!!==================================================================================================

!!!==================================================================================================
  subroutine ftranspose_i2_wrapper(fries_2d,field_3d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:)  , intent(out), allocatable :: fries_2d
    real(kind=rk), dimension(:,:,:), intent(in)               :: field_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%i2.eq.1) then
       call f0transpose_i2_3d_rk(fries_2d,field_3d)
    else
       call ftranspose_i2_3d_rk(fries_2d,field_3d)
    end if

  end subroutine ftranspose_i2_wrapper
!!!==================================================================================================

!!!==================================================================================================
  subroutine btranspose_i2_3d_rk(field_3d,fries_2d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(out)          :: field_3d
    real(kind=rk), dimension(:,:)  , intent(in)           :: fries_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, save, dimension(:,:)  , allocatable          :: pen_fries_mapping
    integer, save, dimension(:)    , allocatable          :: send_counts, recv_counts
    integer, save, dimension(:)    , allocatable          :: send_displs, recv_displs
    integer, save                                         :: pen_fries_nbr, pen_len, pen_comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !real(kind=rk), dimension(:,:)   , allocatable         :: field_2d
    real(kind=rk), dimension(:),save, allocatable         :: field_vector, field_vector_t
    !integer                                               :: s1,s2,s3
    integer                                               :: i,pen_fries_nbr_times_pen_len
    integer                                               :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(send_counts)) then
       ! here the resort tupples and some dimension information for mpi_alltoallv are constructed ONLY once
       ! this is different for different pencils
       ! save     variables are reused in each transpose call
       ! not save variables are just needed once for creation of the resort tupples and mpi_alltoallv args
       call init_i2(send_counts,recv_counts,send_displs,recv_displs,pen_fries_mapping,pen_fries_nbr,pen_len,pen_comm)

       ! allocate temporary vectors
       allocate(field_vector(  sum(send_counts)))
       allocate(field_vector_t(sum(recv_counts)))
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! alltoallv super magic (back) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! pack transposed vector via pen fries mapping tupples
    pen_fries_nbr_times_pen_len = pen_fries_nbr*pen_len
    do i = 1,pen_fries_nbr_times_pen_len
       field_vector_t(i) = fries_2d(pen_fries_mapping(1,i),pen_fries_mapping(2,i))
    end do

    ! transpose (the magic)
    call mpi_alltoallv(field_vector_t, recv_counts, recv_displs, params%parallelism%mpi_type_rk, &  ! send args
         field_vector, send_counts, send_displs, params%parallelism%mpi_type_rk, &                  ! recv args
         pen_comm, ierr)

    ! unpack vector for i1 pencils and fries, different for different pens
    field_3d = reshape(field_vector, shape(field_3d), order = (/2,1,3/))
    
    ! s1 = size(field_3d,1)
    ! s2 = size(field_3d,2)
    ! s3 = size(field_3d,3)    
    ! allocate(field_2d(s1,s2))
    ! l = 1
    ! ! loop different fries i1-i2-i3-order
    ! do k = 1,s3
    !    do i = 1,s1
    !       ! frie lengthwise
    !       field_2d(i,:) = field_vector(l:l+s2-1)
    !       l = l + s2
    !    end do
    !    field_3d(:,:,k) = field_2d
    ! end do

    ! l = 1
    ! ! loop different fries i1-i2-i3-order
    ! do k = 1,size(field_3d,3)
    !    do i = 1,size(field_3d,1)
    !       ! frie lengthwise
    !       do j = 1,size(field_3d,2)
    !          field_3d(i,j,k) = field_vector(l)
    !          l = l + 1
    !       end do
    !    end do
    ! end do

  end subroutine btranspose_i2_3d_rk
!!!==================================================================================================

!!!==================================================================================================
  subroutine b0transpose_i2_3d_rk(field_3d,fries_2d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(out)          :: field_3d
    real(kind=rk), dimension(:,:)  , intent(in)           :: fries_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !integer                                               :: i,k,l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! l = 0
    ! do k = 1,size(field_3d,3)
    !    do i = 1,size(field_3d,1)
    !       l = l + 1
    !       field_3d(i,:,k) = fries_2d(:,l)
    !    end do
    ! end do

    field_3d = reshape(fries_2d, shape(field_3d), order = (/2,1,3/))

  end subroutine b0transpose_i2_3d_rk
!!!==================================================================================================

!!!==================================================================================================
  subroutine btranspose_i2_wrapper(field_3d,fries_2d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(out)          :: field_3d
    real(kind=rk), dimension(:,:)  , intent(in)           :: fries_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%i2.eq.1) then
       call b0transpose_i2_3d_rk(field_3d,fries_2d)
    else
       call btranspose_i2_3d_rk(field_3d,fries_2d)
    end if

  end subroutine btranspose_i2_wrapper
!!!==================================================================================================

!!!==================================================================================================
  subroutine ftranspose_i3_3d_rk(fries_2d,field_3d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:)  , intent(out), allocatable :: fries_2d
    real(kind=rk), dimension(:,:,:), intent(in)               :: field_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, save, dimension(:,:)  , allocatable              :: pen_fries_mapping
    integer, save, dimension(:)    , allocatable              :: send_counts, recv_counts
    integer, save, dimension(:)    , allocatable              :: send_displs, recv_displs
    integer, save                                             :: pen_fries_nbr, pen_len, pen_comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:),save, allocatable             :: field_vector, field_vector_t
    integer                                                   :: s1,s2,s3
    integer                                                   :: i,pen_fries_nbr_times_pen_len
    integer                                                   :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(send_counts)) then
       ! here the resort tupples and some dimension information for mpi_alltoallv are constructed ONLY once
       ! this is different for different pencils
       ! save     variables are reused in each transpose call
       ! not save variables are just needed once for creation of the resort tupples and mpi_alltoallv args
       call init_i3(send_counts,recv_counts,send_displs,recv_displs,pen_fries_mapping,pen_fries_nbr,pen_len,pen_comm)

       ! allocate temporary vectors
       allocate(field_vector(  sum(send_counts)))
       allocate(field_vector_t(sum(recv_counts)))
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! alltoallv super magic (forth) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    s1 = size(field_3d,1)
    s2 = size(field_3d,2)
    s3 = size(field_3d,3)
    
    ! pack field vector for i1 pencils and fries, different for different pens
    field_vector = reshape( reshape(field_3d,(/s3,s1,s2/),order =(/2,3,1/)), shape(field_vector))

    ! allocate(field_2d(s1,s3))

    ! l = 1
    ! ! loop different fries i1-i2-i3-order
    ! do j = 1,s2
    !    field_2d = field_3d(:,j,:)
    !    do i = 1,s1
    !       ! fries lengthwise
    !       field_vector(l:l+s3-1) = field_2d(i,:)
    !       l = l + s3
    !    end do
    ! end do

    ! ! loop different fries i1-i2-i3-order
    ! do j = 1,size(field_3d,2)
    !    do i = 1,size(field_3d,1)
    !       ! fries lengthwise
    !       do k = 1,size(field_3d,3)
    !          field_vector(l) = field_3d(i,j,k)
    !          l = l + 1
    !       end do
    !    end do
    ! end do

    ! transpose (the magic)
    call mpi_alltoallv(field_vector, send_counts, send_displs, params%parallelism%mpi_type_rk, &   ! send args
         field_vector_t, recv_counts, recv_displs, params%parallelism%mpi_type_rk, &               ! recv args
         pen_comm, ierr)                                                       ! mpi stuff

    ! allocate temporary array
    allocate(fries_2d(pen_len,pen_fries_nbr))

    ! unpack transposed vector via pen fries mapping tupples
    pen_fries_nbr_times_pen_len = pen_fries_nbr*pen_len
    do i = 1,pen_fries_nbr_times_pen_len
       fries_2d(pen_fries_mapping(1,i),pen_fries_mapping(2,i)) = field_vector_t(i)
    end do

  end subroutine ftranspose_i3_3d_rk
!!!==================================================================================================

!!!==================================================================================================
  subroutine f0transpose_i3_3d_rk(fries_2d,field_3d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:)  , intent(out), allocatable :: fries_2d
    real(kind=rk), dimension(:,:,:), intent(in)               :: field_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                   :: s1,s2,s3 ! i,j,l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    s1 = size(field_3d,1)
    s2 = size(field_3d,2)
    s3 = size(field_3d,3)
    
    allocate(fries_2d(s3,s1*s2))

    !allocate(tmp(s3,s1,s2))
    !tmp = reshape(field_3d,(/s3,s1,s2/),order =(/2,3,1/))
    !fries_2d = reshape(tmp, shape(fries_2d))

    ! l = 0
    ! do j = 1,size(field_3d,2)
    !    do i = 1,size(field_3d,1)
    !       l = l + 1
    !       fries_2d(:,l) = field_3d(i,j,:)
    !    end do
    ! end do
        
    fries_2d = reshape( reshape(field_3d,(/s3,s1,s2/),order=(/2,3,1/)), shape(fries_2d))

  end subroutine f0transpose_i3_3d_rk
!!!==================================================================================================

!!!==================================================================================================
  subroutine ftranspose_i3_wrapper(fries_2d,field_3d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:)  , intent(out), allocatable :: fries_2d
    real(kind=rk), dimension(:,:,:), intent(in)               :: field_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%i3.eq.1) then
       call f0transpose_i3_3d_rk(fries_2d,field_3d)
    else
       call ftranspose_i3_3d_rk(fries_2d,field_3d)
    end if

  end subroutine ftranspose_i3_wrapper
!!!==================================================================================================

!!!==================================================================================================
  subroutine btranspose_i3_3d_rk(field_3d,fries_2d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(out)          :: field_3d
    real(kind=rk), dimension(:,:)  , intent(in)           :: fries_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, save, dimension(:,:)  , allocatable          :: pen_fries_mapping
    integer, save, dimension(:)    , allocatable          :: send_counts, recv_counts
    integer, save, dimension(:)    , allocatable          :: send_displs, recv_displs
    integer, save                                         :: pen_fries_nbr, pen_len, pen_comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !real(kind=rk), dimension(:,:)   , allocatable         :: field_2d
    real(kind=rk), dimension(:),save,allocatable          :: field_vector, field_vector_t
    integer                                               :: s1,s2,s3
    integer                                               :: i,pen_fries_nbr_times_pen_len
    integer                                               :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(send_counts)) then
       ! here the resort tupples and some dimension information for mpi_alltoallv are constructed ONLY once
       ! this is different for different pencils
       ! save     variables are reused in each transpose call
       ! not save variables are just needed once for creation of the resort tupples and mpi_alltoallv args
       call init_i3(send_counts,recv_counts,send_displs,recv_displs,pen_fries_mapping,pen_fries_nbr,pen_len,pen_comm)

       ! allocate temporary vectors
       allocate(field_vector(  sum(send_counts)))
       allocate(field_vector_t(sum(recv_counts)))
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! alltoallv super magic (back) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    s1 = size(field_3d,1)
    s2 = size(field_3d,2)
    s3 = size(field_3d,3)
    
    ! pack transposed vector via pen fries mapping tupples
    pen_fries_nbr_times_pen_len = pen_fries_nbr*pen_len
    do i = 1,pen_fries_nbr_times_pen_len
       field_vector_t(i) = fries_2d(pen_fries_mapping(1,i),pen_fries_mapping(2,i))
    end do

    ! transpose (the magic)
    call mpi_alltoallv(field_vector_t, recv_counts, recv_displs, params%parallelism%mpi_type_rk, &  ! send args
         field_vector, send_counts, send_displs, params%parallelism%mpi_type_rk, &                  ! recv args
         pen_comm, ierr)

    ! unpack vector for i1 pencils and fries, different for different pens
    ! allocate(field_2d(s1,s3))

    ! l = 1
    ! ! loop different fries i1-i2-i3-order
    ! do j = 1,s2
    !    do i = 1,s1
    !       ! frie lengthwise
    !       field_2d(i,:) = field_vector(l:l+s3-1)
    !       l = l + s3
    !    end do
    !    field_3d(:,j,:) = field_2d
    ! end do

    field_3d = reshape( reshape(field_vector,(/s3,s1,s2/)), shape(field_3d), order = (/3,1,2/) )

    ! ! loop different fries i1-i2-i3-order
    ! do j = 1,size(field_3d,2)
    !   do i = 1,size(field_3d,1)
    !      ! frie lengthwise
    !      do k = 1,size(field_3d,3)
    !         field_3d(i,j,k) = field_vector(l)
    !         l = l + 1
    !      end do
    !   end do
    ! end do

  end subroutine btranspose_i3_3d_rk
!!!==================================================================================================

!!!==================================================================================================
  subroutine b0transpose_i3_3d_rk(field_3d,fries_2d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(out)          :: field_3d
    real(kind=rk), dimension(:,:)  , intent(in)           :: fries_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                               :: s1,s2,s3 ! i,j,l
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    s1 = size(field_3d,1)
    s2 = size(field_3d,2)
    s3 = size(field_3d,3)

    ! l = 0
    ! do j = 1,size(field_3d,2)
    !    do i = 1,size(field_3d,1)
    !       l = l + 1
    !       field_3d(i,j,:) = fries_2d(:,l)
    !    end do
    ! end do

    !field_3d = reshape( transpose(fries_2d), shape(field_3d) )
    
    !field_3d = reshape( reshape(fries_2d,(/s1*s2,s3/),order=(/2,1/)), shape(field_3d) )
    ! inner reshape is transpose
    
    field_3d = reshape( reshape(fries_2d,(/s3,s1,s2/)), shape(field_3d), order = (/3,1,2/) )
    ! this form is direct transferable to btranspose_i3_3d_rk
    
  end subroutine b0transpose_i3_3d_rk
!!!==================================================================================================

!!!==================================================================================================
  subroutine btranspose_i3_wrapper(field_3d,fries_2d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(out)          :: field_3d
    real(kind=rk), dimension(:,:)  , intent(in)           :: fries_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%i3.eq.1) then
       call b0transpose_i3_3d_rk(field_3d,fries_2d)
    else
       call btranspose_i3_3d_rk(field_3d,fries_2d)
    end if

  end subroutine btranspose_i3_wrapper
!!!==================================================================================================

!!!==================================================================================================
  subroutine init_transpose(send_counts,recv_counts,send_displs,recv_displs,pen_fries_mapping, &    ! output args
       total_fries_nr,pen_comm,pen_rank,pen_image_nr,pen_len,image_frie_len)                        ! input  args
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:)  , intent(out), allocatable   :: send_counts, recv_counts
    integer, dimension(:)  , intent(out), allocatable   :: send_displs, recv_displs
    integer, dimension(:,:), intent(out), allocatable   :: pen_fries_mapping
    integer                , intent(in)                 :: total_fries_nr
    integer                , intent(in)                 :: pen_comm,pen_rank
    integer                , intent(in)                 :: pen_image_nr,pen_len
    integer                , intent(in)                 :: image_frie_len
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:)  , allocatable                :: image_frie_lens
    integer, dimension(:)  , allocatable                :: fries_per_pen
    integer                                             :: pen_fries_nbr
    integer                                             :: ierr
    integer                                             :: i,l,n,offset,frie
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! here the resort tupples and some dimension information for mpi_alltoallv are constructed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! nXb of all images
    allocate(image_frie_lens(pen_image_nr))
    call mpi_allgather(image_frie_len,1,mpi_integer,image_frie_lens,1,mpi_integer,pen_comm,ierr)

    allocate(fries_per_pen(pen_image_nr))                                  ! array for all procs in pen
    fries_per_pen = floor(real(total_fries_nr) / real(pen_image_nr))
    do i = pen_image_nr,pen_image_nr-mod(total_fries_nr,pen_image_nr)+1,-1 ! distribute the rest
       fries_per_pen(i) = fries_per_pen(i) + 1
    end do
    pen_fries_nbr = fries_per_pen(pen_rank+1);

    !! build 2d pen fries mapping tupples to unpack recv vector from mpi_alltoallv ------------------
    allocate(pen_fries_mapping(2,pen_fries_nbr*pen_len))
    l        = 0
    offset   = 0
    do i = 1,pen_image_nr             ! loop images
       do frie = 1,pen_fries_nbr      ! loop fries
          do n = 1,image_frie_lens(i) ! loop pts
             l = l + 1
             pen_fries_mapping(1,l) = n + offset
             pen_fries_mapping(2,l) = frie
             if (n + offset.gt.pen_len)  then
                write(*,*) "error in transpose_mpi_cartesian.f90:init_transpose:elements in pen .ne. elements in tupple mapping"
             end if
          end do
       end do
       offset = offset + image_frie_lens(i);
    end do

    !! share alltoallv args
    allocate(send_counts(pen_image_nr))
    allocate(recv_counts(pen_image_nr))
    allocate(send_displs(pen_image_nr))
    allocate(recv_displs(pen_image_nr))

    ! build send_count
    send_counts = fries_per_pen * image_frie_len

    ! build recv_counts
    call mpi_alltoall(send_counts,1,mpi_integer,recv_counts,1,mpi_integer,pen_comm,ierr)

    ! build send_displs
    send_displs(1) = 0
    do i = 2,pen_image_nr
       send_displs(i) = sum(send_counts(1:i-1))
    end do

    ! build recv_displs
    recv_displs(1) = 0
    do i = 2,pen_image_nr
       recv_displs(i) = sum(recv_counts(1:i-1))
    end do

  end subroutine init_transpose
!!!==================================================================================================

!!!==================================================================================================
  subroutine init_i1(send_counts,recv_counts,send_displs,recv_displs,pen_fries_mapping,pen_fries_nbr,pen_len,pen_comm) ! output args
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:)  , intent(out), allocatable   :: send_counts, recv_counts
    integer, dimension(:)  , intent(out), allocatable   :: send_displs, recv_displs
    integer, dimension(:,:), intent(out), allocatable   :: pen_fries_mapping
    integer                , intent(out)                :: pen_fries_nbr, pen_len, pen_comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                             :: pen_rank, pen_image_nr
    integer                                             :: total_fries_nr, image_frie_len
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! here the resort tupples and some dimension information for mpi_alltoallv are constructed for the i1 pencil

    !! pen mpi cartesian stuff
    pen_comm        = params%parallelism%i1pen_comm
    pen_rank        = params%parallelism%i1pen_rank
    pen_image_nr    = params%parallelism%i1pen_size

    !! pen dimension stuff (support variables)
    pen_len         = params%geom%n1
    total_fries_nr  = params%geom%n2b*params%geom%n3b ! pts per i1 frie in i2/i3 plane

    !! image dimension stuff
    image_frie_len  = params%geom%n1b

    ! get mpi_alltoallv args and resort tupples
    call init_transpose(send_counts,recv_counts,send_displs,recv_displs,pen_fries_mapping, &     ! output args
         total_fries_nr,pen_comm,pen_rank,pen_image_nr,pen_len,image_frie_len)                   ! input  args  

    pen_fries_nbr = size(pen_fries_mapping,2)/pen_len

    call set_parameter(pen_fries_nbr,'parallelism.i1_pen_fries_nbr')
    call set_parameter(pen_len,'parallelism.i1_pen_len')

  end subroutine init_i1
!!!==================================================================================================

!!!==================================================================================================
  subroutine init_i2(send_counts,recv_counts,send_displs,recv_displs,pen_fries_mapping,pen_fries_nbr,pen_len,pen_comm) ! output args
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:)  , intent(out), allocatable   :: send_counts, recv_counts
    integer, dimension(:)  , intent(out), allocatable   :: send_displs, recv_displs
    integer, dimension(:,:), intent(out), allocatable   :: pen_fries_mapping
    integer                , intent(out)                :: pen_fries_nbr, pen_len, pen_comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                             :: pen_rank, pen_image_nr
    integer                                             :: total_fries_nr, image_frie_len
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! here the resort tupples and some dimension information for mpi_alltoallv are constructed for the i2 pencil

    !! pen mpi cartesian stuff
    pen_comm        = params%parallelism%i2pen_comm
    pen_rank        = params%parallelism%i2pen_rank
    pen_image_nr    = params%parallelism%i2pen_size

    !! pen dimension stuff (support variables)
    pen_len         = params%geom%n2
    total_fries_nr  = params%geom%n1b*params%geom%n3b ! pts per i2 frie in i1/i3 plane

    !! image dimension stuff
    image_frie_len  = params%geom%n2b

    ! get mpi_alltoallv args and resort tupples
    call init_transpose(send_counts,recv_counts,send_displs,recv_displs,pen_fries_mapping, &     ! output args
         total_fries_nr,pen_comm,pen_rank,pen_image_nr,pen_len,image_frie_len)                   ! input  args  

    pen_fries_nbr = size(pen_fries_mapping,2)/pen_len

    call set_parameter(pen_fries_nbr,'parallelism.i2_pen_fries_nbr')
    call set_parameter(pen_len,'parallelism.i2_pen_len')

  end subroutine init_i2
!!!==================================================================================================

!!!==================================================================================================
  subroutine init_i3(send_counts,recv_counts,send_displs,recv_displs,pen_fries_mapping,pen_fries_nbr,pen_len,pen_comm) ! output args
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:)  , intent(out), allocatable   :: send_counts, recv_counts
    integer, dimension(:)  , intent(out), allocatable   :: send_displs, recv_displs
    integer, dimension(:,:), intent(out), allocatable   :: pen_fries_mapping
    integer                , intent(out)                :: pen_fries_nbr, pen_len, pen_comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                             :: pen_rank, pen_image_nr
    integer                                             :: total_fries_nr, image_frie_len
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! here the resort tupples and some dimension information for mpi_alltoallv are constructed for the i3 pencil

    !! pen mpi cartesian stuff
    pen_comm        = params%parallelism%i3pen_comm
    pen_rank        = params%parallelism%i3pen_rank
    pen_image_nr    = params%parallelism%i3pen_size

    !! pen dimension stuff (support variables)
    pen_len         = params%geom%n3
    total_fries_nr  = params%geom%n1b*params%geom%n2b ! pts per i1 frie in i1/i2 plane

    !! image dimension stuff
    image_frie_len  = params%geom%n3b

    ! get mpi_alltoallv args and resort tupples
    call init_transpose(send_counts,recv_counts,send_displs,recv_displs,pen_fries_mapping, &     ! output args
         total_fries_nr,pen_comm,pen_rank,pen_image_nr,pen_len,image_frie_len)                   ! input  args  

    pen_fries_nbr = size(pen_fries_mapping,2)/pen_len

    call set_parameter(pen_fries_nbr,'parallelism.i3_pen_fries_nbr')
    call set_parameter(pen_len,'parallelism.i3_pen_len')

  end subroutine init_i3
!!!==================================================================================================

end module transpose_mpi_cartesian
