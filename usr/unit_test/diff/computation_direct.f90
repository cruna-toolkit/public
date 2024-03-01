module computation_direct

  use helper
  USE mpi
  use parallelism
  use parameter
  use data_geom
  use discretisation_x

contains

!!!=================================================================================================
  subroutine calc_direct(Q,q0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:), intent(inout) :: Q
    real(kind=rk), dimension(:,:,:,:)  , intent(in)    :: q0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:)  , allocatable   :: y
    real(kind=rk), dimension(:,:,:)    , allocatable   :: y_x
    real(kind=rk), dimension(3)                        :: tot_errs
    real(kind=rk)                                      :: loc_err,tot_err
    real(kind=rk)                                      :: x10, x20, x30
    real(kind=rk)                                      :: x11, x21, x31
    real(kind=rk)                                      :: l1, l2, l3
    integer                                            :: i,j,k,l
    integer                                            :: n1,n2,n3
    integer                                            :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                      :: pi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%world_image.eq.1) then
       write(*,*) "============================"
       write(*,*) "DIFF_TEST"
       write(*,*) "============================"
       write(*,*) " .. works only on uniform grids"
    end if

    ! allocates and constants
    allocate(  y(params%geom%n1b,params%geom%n2b,params%geom%n3b,5))
    allocate(y_x(params%geom%n1b,params%geom%n2b,params%geom%n3b  ))

    pi = 4.0*atan(1.0)

    ! fill y(x)
    if(params%parallelism%world_image.eq.1) then
       write(*,*) "y(x) (components)"
       write(*,*) "1: constant zero"
       write(*,*) "2: constant one"
       write(*,*) "3: sin(x1)"
       write(*,*) "4: sin(x2)"
       write(*,*) "5: sin(x3)"
       write(*,*) ""
    end if

    ! get physical block dimensions
    call get_parameter(x10,'geom.x10',default = 0.0_rk)
    call get_parameter(x20,'geom.x20',default = 0.0_rk)
    call get_parameter(x30,'geom.x30',default = 0.0_rk)

    call get_parameter(x11,'geom.x11',default = 1.0_rk)
    call get_parameter(x21,'geom.x21',default = 1.0_rk)
    call get_parameter(x31,'geom.x31',default = 1.0_rk)

    l1 = x11 - x10
    l2 = x21 - x20
    l3 = x31 - x30

    ! set components
    y(:,:,:,1) = 0.0_rk
    y(:,:,:,2) = 1.0_rk
    y(:,:,:,3) = sin((X(:,:,:,1)))
    y(:,:,:,4) = sin((X(:,:,:,2)))
    y(:,:,:,5) = sin((X(:,:,:,3)))

!!! component 1
    call Dx1(y_x,y(:,:,:,1),params%geom%dx1)
    loc_err = sum(abs(y_x - 0.0_rk))/(params%geom%n1b * params%geom%n2b * params%geom%n3b) 
    call mpi_reduce(loc_err,tot_err,1,mpi_double,mpi_sum,0,params%parallelism%block_comm,ierr)
    tot_errs(1) = tot_err

    call Dx2(y_x,y(:,:,:,1),params%geom%dx1)
    loc_err = sum(abs(y_x - 0.0_rk))/(params%geom%n1b * params%geom%n2b * params%geom%n3b)
    call mpi_reduce(loc_err,tot_err,1,mpi_double,mpi_sum,0,params%parallelism%block_comm,ierr)
    tot_errs(2) = tot_err

    call Dx3(y_x,y(:,:,:,1),params%geom%dx1)
    loc_err = sum(abs(y_x - 0.0_rk))/(params%geom%n1b * params%geom%n2b * params%geom%n3b)
    call mpi_reduce(loc_err,tot_err,1,mpi_double,mpi_sum,0,params%parallelism%block_comm,ierr)
    tot_errs(3) = tot_err

    if (params%parallelism%world_image.eq.1) then
       write(*,*) "c1, tot_err:",tot_errs
    end if

!!! component 2
    call Dx1(y_x,y(:,:,:,2),params%geom%dx1)
    loc_err = sum(abs(y_x - 0.0_rk))/(params%geom%n1b * params%geom%n2b * params%geom%n3b)
    call mpi_reduce(loc_err,tot_err,1,mpi_double,mpi_sum,0,params%parallelism%block_comm,ierr)
    tot_errs(1) = tot_err

    call Dx2(y_x,y(:,:,:,2),params%geom%dx1)
    loc_err = sum(abs(y_x - 0.0_rk))/(params%geom%n1b * params%geom%n2b * params%geom%n3b)
    call mpi_reduce(loc_err,tot_err,1,mpi_double,mpi_sum,0,params%parallelism%block_comm,ierr)
    tot_errs(2) = tot_err

    call Dx3(y_x,y(:,:,:,2),params%geom%dx1)
    loc_err = sum(abs(y_x - 0.0_rk))/(params%geom%n1b * params%geom%n2b * params%geom%n3b)
    call mpi_reduce(loc_err,tot_err,1,mpi_double,mpi_sum,0,params%parallelism%block_comm,ierr)
    tot_errs(3) = tot_err

    if (params%parallelism%world_image.eq.1) then
       write(*,*) "c2, tot_err:",tot_errs
    end if

!!! component 3
    call Dx1(y_x,y(:,:,:,3),params%geom%dx1)
    loc_err = sum(abs(y_x - cos(X(:,:,:,1))))/(params%geom%n1b * params%geom%n2b * params%geom%n3b)
    call mpi_reduce(loc_err,tot_err,1,mpi_double,mpi_sum,0,params%parallelism%block_comm,ierr)

    if (params%parallelism%world_image.eq.1) then
       write(*,*) "c3, Dx1, tot_err:",tot_err
    end if

!!! component 4
    call Dx2(y_x,y(:,:,:,4),params%geom%dx2)
    loc_err = sum(abs(y_x - cos(X(:,:,:,2))))/(params%geom%n1b * params%geom%n2b * params%geom%n3b)
    call mpi_reduce(loc_err,tot_err,1,mpi_double,mpi_sum,0,params%parallelism%block_comm,ierr)

    if (params%parallelism%world_image.eq.1) then
       write(*,*) "c4, Dx2, tot_err:",tot_err
    end if

!!! component 5
    call Dx3(y_x,y(:,:,:,5),params%geom%dx3)
    loc_err = sum(abs(y_x - cos(X(:,:,:,3))))/(params%geom%n1b * params%geom%n2b * params%geom%n3b)
    call mpi_reduce(loc_err,tot_err,1,mpi_double,mpi_sum,0,params%parallelism%block_comm,ierr)

    if (params%parallelism%world_image.eq.1) then
       write(*,*) "c5, Dx3, tot_err:",tot_err
    end if

    call stop_cruna

  end subroutine calc_direct
!!!=================================================================================================

end module computation_direct
