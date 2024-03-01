module grid_transform

  use discretisation_x
  use helper
  use parameter
  use io
  use mpi_cartesian

  private

  public :: cartesian_2_cartesian_strech_shift
  public :: distort_wave
  public :: distort_local_refinement
  public :: distort_user_defined

contains


!!!=================================================================================================
  subroutine cartesian_2_cartesian_strech_shift(X,x10i,x20i,x30i)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:), intent(inout) :: X                                 ! local block grid   
    integer                         , intent(in)    :: x10i,x20i,x30i                    ! start indices of imaginary block grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:),allocatable          :: x1 , x2 , x3                      ! grid vectors for nd grid
    real(kind=rk)                                   :: x10, x20, x30                     ! global min
    real(kind=rk)                                   :: x11, x21, x31                     ! global max 
    real(kind=rk)                                   :: dx1, dx2, dx3                     ! dx 
    integer                                         :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! get physical block dimensions
    call get_parameter(x10,'geom.x10',default = 0.0_rk)
    call get_parameter(x20,'geom.x20',default = 0.0_rk)
    call get_parameter(x30,'geom.x30',default = 0.0_rk)

    call get_parameter(x11,'geom.x11',default = 1.0_rk)
    call get_parameter(x21,'geom.x21',default = 1.0_rk)
    call get_parameter(x31,'geom.x31',default = 1.0_rk)

    ! compute block dx
    if (trim(diff_i1_type).eq.'periodic') then
       dx1 = (x11-x10)/real(params%geom%n1 - 0,rk)
    else
       dx1 = (x11-x10)/real(params%geom%n1 - 1,rk)
    end if

    if (trim(diff_i2_type).eq.'periodic') then
       dx2 = (x21-x20)/real(params%geom%n2 - 0,rk)
    else
       dx2 = (x21-x20)/real(params%geom%n2 - 1,rk)
    end if

    if (trim(diff_i3_type).eq.'periodic') then
       dx3 = (x31-x30)/real(params%geom%n3 - 0,rk)
    else
       dx3 = (x31-x30)/real(params%geom%n3 - 1,rk)
    end if

    params%geom%dx1 = dx1
    params%geom%dx2 = dx2
    params%geom%dx3 = dx3

    call set_parameter(dx1,'geom.dx1')
    call set_parameter(dx2,'geom.dx2')
    call set_parameter(dx3,'geom.dx3')

    ! allocate grid vectors
    allocate(x1(params%geom%n1b))
    allocate(x2(params%geom%n2b))
    allocate(x3(params%geom%n3b))

    x1 = 0._rk
    x2 = 0._rk
    x3 = 0._rk

    ! write(*,*) params%geom%n1b,params%geom%n3b,params%geom%n3b,params%parallelism%block_image

    if (params%geom%n1.gt.1) then                               ! check for collapsed (global) grid
       do i = 1,params%geom%n1b                                 ! loop over all grid points
          x1(i) = real(i - 1 + x10i - 1,rk) * dx1
          !            ^----------^ global grid index        
       end do
       ! else
       !    x1 = x10                                            ! use minimum value for collapsed grid
    end if

    ! origin
    x1 = x1 + x10

    if (params%geom%n2.gt.1) then
       do i = 1,params%geom%n2b
          x2(i) = real(i - 1 + x20i - 1,rk) * dx2
          !            ^----------^ global grid index   
       end do
       ! else
       !    x2 = x20
    end if

    ! origin
    x2 = x2 + x20

    if (params%geom%n3.gt.1) then
       do i = 1,params%geom%n3b
          x3(i) = real(i - 1 + x30i - 1,rk) * dx3
          !            ^----------^ global grid index   
       end do
       ! else
       !    x3 = x30
    end if

    ! origin
    x3 = x3 + x30

    ! debug
    ! write(*,*) 'x1 start dx end',x1(1),dx1,x1(params%geom%n1b), params%parallelism%block_image
    ! write(*,*) 'x2 start dy end',x2(1),dx2,x2(params%geom%n2b), params%parallelism%block_image
    ! write(*,*) 'x3 start dz end',x3(1),dx3,x3(params%geom%n3b), params%parallelism%block_image

    ! create the full grid
    call ndgrid(X,x1,x2,x3)
  end subroutine cartesian_2_cartesian_strech_shift
!!!==========================================================================================================


!!!=================================================================================================
  subroutine distort_wave(X)   
    ! a test grid which looks like a flag in wind, no deeper technical meaning 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:), intent(inout) :: X                                 ! local block grid      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                      
    real(kind=rk), dimension(:,:,:,:), allocatable  :: X_temp                            ! 
    real(kind=rk)                                   :: alpha   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_parameter(alpha,'geom.distort_alpha')

    allocate(X_temp(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))
    X_temp     = X
    if (alpha.ne.0) then
       X(:,:,:,1) = X_temp(:,:,:,1) + alpha*sin(X_temp(:,:,:,1) + X_temp(:,:,:,2) + X_temp(:,:,:,3))              ! here X_1 is overwritten 
       X(:,:,:,2) = X_temp(:,:,:,2) + alpha*sin(X_temp(:,:,:,1) + X_temp(:,:,:,2) + X_temp(:,:,:,3))    
       X(:,:,:,3) = X_temp(:,:,:,3) + alpha*sin(X_temp(:,:,:,1) + X_temp(:,:,:,2) + X_temp(:,:,:,3)) 
    end if

  end subroutine distort_wave
!!!=================================================================================================

!!!=================================================================================================
  subroutine distort_local_refinement(X)   
    ! 2D refinement for bluff body simulations (can be adjusted for 3D)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:), intent(inout) :: X                                 ! local block grid                            ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable  :: X_temp   
    real(kind=rk)                                   :: alpha                     
    real(kind=rk)                                   :: delta 
    real(kind=rk)                                   :: sigma 
    real(kind=rk)                                   :: x0,y0
    real(kind=rk)                                   :: x10, x20, x30                     ! global min
    real(kind=rk)                                   :: x11, x21, x31                     ! global max 
    real(kind=rk)                                   :: L1, L2, L3                        ! global length
    real(kind=rk)                                   :: x1_min, x2_min, x3_min            ! global new min
    real(kind=rk)                                   :: x1_max, x2_max, x3_max            ! global new max
    integer                                         :: n1, n2, n3                        ! global number of grid points
    integer, dimension(3), save                     :: point_idx = (/0,0,0/)
    integer                                         :: image     = 0                     ! core
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_parameter(alpha,'geom.distort_alpha')
    call get_parameter(delta,'geom.distort_delta')
    call get_parameter(sigma,'geom.distort_sigma')
    call get_parameter(x0   ,'geom.distort_x0')
    call get_parameter(y0   ,'geom.distort_y0')

    allocate(X_temp(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))
    X_temp = X
    ! compute x_tilde (see Reiss2022 "Pressure-Tight and Non-stiff Volume Penalization for Compressible Flows")
    ! length correction (considering L  +2*alpha)

    ! get physical block dimensions
    call get_parameter(x10,'geom.x10',default = 0.0_rk)
    call get_parameter(x20,'geom.x20',default = 0.0_rk)
    call get_parameter(x30,'geom.x30',default = 0.0_rk)
    call get_parameter(x11,'geom.x11',default = 1.0_rk)
    call get_parameter(x21,'geom.x21',default = 1.0_rk)
    call get_parameter(x31,'geom.x31',default = 1.0_rk)

    call get_parameter(n1,'geom.n1',default = 1)
    call get_parameter(n2,'geom.n2',default = 1)
    call get_parameter(n3,'geom.n3',default = 1)

    ! compute physical (global) domain length
    L1 = (x11-x10) 
    L2 = (x21-x20) 
    L3 = (x31-x30) 

    ! compute length correction
    X_temp(:,:,:,1) = X_temp(:,:,:,1) /  L1 * ( L1+ 2 * alpha) ! /Lxi   * (Lxi + 2*alpha)
    X_temp(:,:,:,2) = X_temp(:,:,:,2) /  L2 * ( L2+ 2 * alpha) ! /Leta  * (Leta + 2*alpha)
    !     X_temp(:,:,:,3) = X_temp(:,:,:,3) /  L3 * ( L3 + 2 * alpha) ! currently no z deformation

    ! compute refinement
    X(:,:,:,1) = X_temp(:,:,:,1) - alpha * (tanh( (X_temp(:,:,:,1) - (x0 + alpha + delta))/ sigma) + 1)
    X(:,:,:,2) = X_temp(:,:,:,2) - alpha * (tanh( (X_temp(:,:,:,2) - (y0 + alpha)        )/ sigma)    ) -alpha
    X(:,:,:,3) = X_temp(:,:,:,3)

    ! get global minimum and maximum values
    call get_index(point_idx,image,x1_min,X(:,:,:,1),'min')
    call get_index(point_idx,image,x2_min,X(:,:,:,2),'min')

    call get_index(point_idx,image,x1_max,X(:,:,:,1),'max')
    call get_index(point_idx,image,x2_max,X(:,:,:,2),'max')

    ! correct for original block length (account for periodicity)
    if (trim(diff_i1_type).eq.'periodic') then
       X(:,:,:,1) =  (X(:,:,:,1) - x1_min)/(x1_max - x1_min) * L1 * ((n1-1)/n1)   
    else
       X(:,:,:,1) =  (X(:,:,:,1) - x1_min)/(x1_max - x1_min) * L1  
    end if

    if (trim(diff_i2_type).eq.'periodic') then
       X(:,:,:,2) =  (X(:,:,:,2) - x2_min)/(x2_max - x2_min) * L2 * ((n2-1)/n2)       
    else
       X(:,:,:,2) =  (X(:,:,:,2) - x2_min)/(x2_max - x2_min) * L2 
    end if

    !     if (trim(diff_i3_type).eq.'periodic') then
    !        X(:,:,:,3) =  (X(:,:,:,3) - X(1,1,1,3))/(X(1,1,size(X,3),3) - X(1,1,1,3)) * (x31-x30) * ((size(X,3)-1)/size(X,3))        
    !     else
    !        X(:,:,:,3) =  (X(:,:,:,3) - X(1,1,1,3))/(X(1,1,size(X,3),3) - X(1,1,1,3)) * (x31-x30) 
    !     end if

  end subroutine distort_local_refinement
!!!=================================================================================================


!!!=================================================================================================
  subroutine distort_user_defined(X,x10i,x20i,x30i)   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout) :: X                                 ! local block grid       
    integer                         , intent(in)     :: x10i,x20i,x30i                    ! start indices of imaginary block grid! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:),allocatable          :: x1 , x2 , x3                      ! grid vectors for nd grid (local)
    real(kind=rk), dimension(:), allocatable         :: x_vec, y_vec, z_vec               ! grid vectors for nd grid (global)
    integer                                          :: dim_x_vec(1), dim_y_vec(1), dim_z_vec(1)
    character(len=max_length_fname)                  :: filename
    character(len=max_length_fname)                  :: data_directory
    integer                                          :: image     = 1                     ! core
    integer                                          :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! add working path if present
    call get_parameter(data_directory,'geom.data_directory',default = ".")
    data_directory = trim(data_directory) 
    
    
    call get_parameter(filename,'geom.filename',default = ".")
    filename = trim(filename) 
    
    ! first core read grid vector file and gets dimensions
    filename = trim(data_directory)//trim(filename) 
!     if (params%parallelism%world_image.eq.image) write(*,*) 'filename new', filename
    if (params%parallelism%world_image.eq.image) then
       call get_dset_dim(dim_x_vec, 'dummy', fname_optin=  filename, dset_name_optin='x')
       call get_dset_dim(dim_y_vec, 'dummy', fname_optin=  filename, dset_name_optin='y')
       call get_dset_dim(dim_z_vec, 'dummy', fname_optin=  filename, dset_name_optin='z')
    end if  
    
    ! broad cast dimensions and allocate memory
    call bcast(dim_x_vec,image-1 , params%parallelism%block_comm)
    call bcast(dim_y_vec,image-1 , params%parallelism%block_comm)
    call bcast(dim_z_vec,image-1 , params%parallelism%block_comm)
    
    allocate(x_vec(dim_x_vec(1)))
    allocate(y_vec(dim_y_vec(1)))
    allocate(z_vec(dim_z_vec(1)))
    
    
    if (params%parallelism%world_image.eq.image) then
       call load(x_vec, 'dummy', fname_optin= filename,dset_name_optin='x')
       call load(y_vec, 'dummy', fname_optin= filename,dset_name_optin='y')
       call load(z_vec, 'dummy', fname_optin= filename,dset_name_optin='z')
    end if
        
    
    ! broad cast dimensions and allocate memory
    call bcast(x_vec,image-1 , params%parallelism%block_comm)
    call bcast(y_vec,image-1 , params%parallelism%block_comm)
    call bcast(z_vec,image-1 , params%parallelism%block_comm)
    
       
    ! allocate grid vectors
    allocate(x1(params%geom%n1b))
    allocate(x2(params%geom%n2b))
    allocate(x3(params%geom%n3b))
   
    x1 = 0._rk
    x2 = 0._rk
    x3 = 0._rk


    if (params%geom%n1.gt.1) then                               ! check for collapsed (global) grid
       do i = 1,params%geom%n1b                                 ! loop over all grid points
          x1(i) = x_vec(i - 1 + x10i) 
          !             ^----------^ global grid index      
       end do
    end if

    if (params%geom%n2.gt.1) then
       do i = 1,params%geom%n2b
          x2(i) = y_vec(i - 1 + x20i) 
          !             ^----------^ global grid index   

       end do
    end if

    if (params%geom%n3.gt.1) then
       do i = 1,params%geom%n3b
          x3(i) = z_vec(i - 1 + x30i) 
          !             ^----------^ global grid index  
       end do
    end if    
    
   
    ! create the full grid
    call ndgrid(X,x1,x2,x3)


  end subroutine distort_user_defined
!!!=================================================================================================




end module grid_transform
