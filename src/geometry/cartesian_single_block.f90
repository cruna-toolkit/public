module cartesian_single_block

  USE mpi
  use data_geom
  use grid_metric
  use grid_transform  
  use helper
  use io
  use parallelism
  use parameter

  private

  public :: init_geometry
  public :: calc_spatial_distance

  interface calc_spatial_distance
     module procedure calc_spatial_distance_rk_3d
  end interface calc_spatial_distance

contains

!!!=================================================================================================
  subroutine calc_spatial_distance_rk_3d(distance,X,point)
    real(kind=rk)   , dimension(:,:,:)  , intent(out) :: distance
    real(kind=rk)   , dimension(:,:,:,:), intent(in)  :: X
    real(kind=rk)   , dimension(3)      , intent(in)  :: point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    distance = sqrt( (X(:,:,:,1)-point(1))**2 + (X(:,:,:,2)-point(2))**2 + (X(:,:,:,3)-point(3))**2 )

  end subroutine calc_spatial_distance_rk_3d
!!!=================================================================================================


!!!=================================================================================================
  subroutine init_geometry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , dimension(:,:,:,:), allocatable :: X_calc_space 
    integer         , dimension(:)      , allocatable :: x1i, x2i, x3i
    integer         , dimension(3)                    :: mpi_cart_coord
    integer                                           :: x10i,x20i,x30i
    integer                                           :: x11i,x21i,x31i
    integer                                           :: i,s
    integer                                           :: ierr
    character(len=max_length_parameter)               :: grid_transform_type
    logical                                           :: grid_deformed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
!!! create local global block indices
    if (params%parallelism%world_size.eq.1) then
       ! no parallelism (or just one process)
       mpi_cart_coord = (/0,0,0/)

       x10i = 1
       x11i = params%geom%n1

       x20i = 1
       x21i = params%geom%n2

       x30i = 1
       x31i = params%geom%n3

    else
       call mpi_cart_coords(params%parallelism%block_comm,params%parallelism%block_image-1,size(mpi_cart_coord),mpi_cart_coord,ierr)

       if (ierr.ne.0) then
          write(*,*) "error in cartesian_single_block.f90:init_geometry:mpi_cart_coords:",ierr
          call stop_cruna
       end if

       !! create local cartesian grid, use for calculation of different grids
       ! x1-direction
       s = 0
       do i = 1,mpi_cart_coord(1) + 1
          s = s + floor(real(params%geom%n1,rk) / real(params%parallelism%i1,rk))                         ! seperate into equal parts
          if (mod(params%geom%n1,params%parallelism%i1).ne.0) then                                        ! distribute extra points, analogous to mpi_cartesian
             if ((i).gt.(params%parallelism%i1 - mod(params%geom%n1,params%parallelism%i1)) ) then        ! the idea here is to compute the number of grid points
                s = s + 1                                                                                 ! of all images in front of the actual and and to sum them up
             end if
          end if
       end do
       x10i = s - params%geom%n1b + 1
       x11i = s

       ! x2-direction
       s = 0
       do i = 1,mpi_cart_coord(2) + 1
          s = s + floor(real(params%geom%n2,rk) / real(params%parallelism%i2,rk))       
          if (mod(params%geom%n2,params%parallelism%i2).ne.0) then                      
             if ((i).gt.(params%parallelism%i2 - mod(params%geom%n2,params%parallelism%i2)) ) then
                s = s + 1
             end if
          end if
       end do
       x20i = s - params%geom%n2b + 1
       x21i = s

       ! x3-direction
       s = 0
       do i = 1,mpi_cart_coord(3) + 1
          s = s + floor(real(params%geom%n3,rk) / real(params%parallelism%i3,rk))       
          if (mod(params%geom%n3,params%parallelism%i3).ne.0) then                      
             if ((i).gt.(params%parallelism%i3 - mod(params%geom%n3,params%parallelism%i3)) ) then
                s = s + 1
             end if
          end if
       end do
       x30i = s - params%geom%n3b + 1
       x31i = s
    end if

!     write(*,*) x10i,x11i,x20i,x21i,x30i,x31i,params%parallelism%block_image

!!! create and store index grid (needed for e.g. block reconstruction)
    allocate(x1i(params%geom%n1b))
    allocate(x2i(params%geom%n2b))
    allocate(x3i(params%geom%n3b))

    call linspace(x1i,x10i,x11i)
    call linspace(x2i,x20i,x21i)
    call linspace(x3i,x30i,x31i)

    call ndgrid(Xi,x1i,x2i,x3i)

    call store(Xi,'geometry_indices')

    call set_parameter(x10i,'geom.xi10i')
    call set_parameter(x11i,'geom.xi11i')
    call set_parameter(x20i,'geom.xi20i')
    call set_parameter(x21i,'geom.xi21i')
    call set_parameter(x30i,'geom.xi30i')
    call set_parameter(x31i,'geom.xi31i')

!     write(*,*) x1i,x2i,x3i,params%parallelism%block_image

!!! compute the physical grid based on the obtained indices
    call get_parameter(grid_transform_type,'geom.grid_transform_type',default = 'cartesian_2_cartesian_strech_shift')
    call get_parameter(grid_deformed      ,'geom.grid_deformed'      ,default = .false.)

   select case (trim(grid_transform_type))    


    case ('cartesian_2_cartesian_strech_shift')
       call cartesian_2_cartesian_strech_shift(X,x10i,x20i,x30i)

    case ('cartesian_2_distorted_wave')
       if (grid_deformed.eqv..false.) then
          write(*,*) "error in cartesian_single_block.f90: using grid transform (cartesian_2_distorted_wave) without setting geom.grid_deformed = .true."
          stop
       end if

       call cartesian_2_cartesian_strech_shift(X,x10i,x20i,x30i)

       allocate( X_calc_space(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))
       X_calc_space = X     ! X_temp are is the calculation space. It is needed to remove the non periodic part in calc_metric  
       call distort_wave(X)  

    case ('cartesian_2_distorted_local_refinement')
       if(params%parallelism%block_image==1) write(*,*) 'trying new distortion (refinement)'
       if (grid_deformed.eqv..false.) then
          write(*,*) "error in cartesian_single_block.f90: using grid transform (cartesian_2_distorted_local_refinement) without setting geom.grid_deformed = .true."
          stop
       end if

       call cartesian_2_cartesian_strech_shift(X,x10i,x20i,x30i)

       allocate( X_calc_space(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))
       X_calc_space = X     ! X_temp are is the calculation space. It is needed to remove the non periodic part in calc_metric  
       call distort_local_refinement(X) 
       
    case ('cartesian_2_distorted_user_defined')
       if(params%parallelism%block_image==1) write(*,*) 'using user defined grid'
       if (grid_deformed.eqv..false.) then
          write(*,*) "error in cartesian_single_block.f90: using grid transform (cartesian_2_distorted_local_refinement_ndgrid) without setting geom.grid_deformed = .true."
          stop
       end if

       call cartesian_2_cartesian_strech_shift(X,x10i,x20i,x30i)

       allocate( X_calc_space(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))
       X_calc_space = X     ! X_temp are is the calculation space. It is needed to remove the non periodic part in calc_metric  
       call distort_user_defined(X,x10i,x20i,x30i)       
    case default
       write(*,*) "error in cartesian_single_block.f90:init_geometry:unknown grid_transform_type:",trim(grid_transform_type)
       call stop_cruna

    end select

    ! store physical geometry
    call store(X,'geometry')

    ! call metric if needed    
    if (grid_deformed.eqv..true.) then
       if ((params%io%verbosity.ge.1).and.(params%parallelism%world_image.eq.1)) then
          write(*,*) " calc metric coefficients" 
       end if
       call calc_metric(X_calc_space)                                  ! actually X is the important quantity, but it is global
       call store(X_metric  ,'geometry_metrics')
       call store(X_jacobian,'geometry_jacobian')
    end if

  end subroutine init_geometry
!!!=================================================================================================

end module cartesian_single_block
