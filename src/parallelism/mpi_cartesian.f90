module mpi_cartesian

  use boundary_conditions
  USE mpi
  use parameter

  private

  public  :: allreduce
  public  :: bcast
  public  :: end_parallelism
  public  :: init_parallelism
  public  :: init_topology
  public  :: spread_parameter
  public  :: spread_boundary_conditions_array
  public  :: stop_cruna
  private :: init_parameter_array_block
  private :: spread_parameter_array_world
  private :: spread_parameter_array_block
  private :: spread_boundary_conditions_array_world  

  interface allreduce
     module procedure allreduce_ik
     module procedure allreduce_rk
     module procedure allreduce_rk_loc
  end interface allreduce

  interface bcast
     module procedure bcast_ik
     module procedure bcast_ik_1d
     module procedure bcast_rk
     module procedure bcast_rk_1d
  end interface bcast

contains 

!!!=================================================================================================
  subroutine init_parallelism
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                               :: ierr, comm_rank, comm_size
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! type
    params%parallelism%type = "mpi_cartesian"

    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world,comm_rank,ierr)
    call mpi_comm_size(mpi_comm_world,comm_size,ierr)

    ! set spmd block and image
    params%parallelism%world_comm  = mpi_comm_world
    params%parallelism%world_size  = comm_size
    params%parallelism%world_image = comm_rank + 1
    params%parallelism%world_block = 1

    call set_parameter(params%parallelism%world_comm  ,'parallelism.world_comm' )
    call set_parameter(params%parallelism%world_size  ,'parallelism.world_size' )
    call set_parameter(params%parallelism%world_image ,'parallelism.world_image')
    call set_parameter(params%parallelism%world_block ,'parallelism.world_block')

    params%parallelism%block_comm  = mpi_comm_world
    params%parallelism%block_size  = comm_size
    params%parallelism%block_image = comm_rank + 1

    call set_parameter(params%parallelism%block_comm  ,'parallelism.block_comm' )
    call set_parameter(params%parallelism%block_size  ,'parallelism.block_size' )
    call set_parameter(params%parallelism%block_image ,'parallelism.block_image')

    ! types and stuff
    params%parallelism%mpi_type_ik    = MPI_INTEGER
    select case (rk)
    case(4)
       params%parallelism%mpi_type_rk  = MPI_REAL
       params%parallelism%mpi_type_2rk = MPI_2REAL               ! needed for minloc,maxloc
       params%parallelism%mpi_type_ck  = MPI_COMPLEX
    case(8)
       params%parallelism%mpi_type_rk  = MPI_DOUBLE_PRECISION
       params%parallelism%mpi_type_2rk = MPI_2DOUBLE_PRECISION   ! needed for minloc,maxloc
       params%parallelism%mpi_type_ck  = MPI_DOUBLE_COMPLEX
    case default
       if (params%parallelism%world_image.eq.1) then
          write(*,*) "error in mpi_cartesian.f90:init_parallelism:real rk not implemented:",rk
       end if
    end select
    call set_parameter(params%parallelism%mpi_type_rk,'parallelim.mpi_type_rk')

  end subroutine init_parallelism
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_topology
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(3)                                 :: dims
    integer, dimension(3)                                 :: mpi_cart_coord
    logical, dimension(3)                                 :: periodic
    integer                                               :: ierr, new_comm, comm_rank, comm_size
    integer                                               :: i1 , i2 , i3 
    integer                                               :: n1b, n2b, n3b 
    logical, parameter                                    :: reorder = .true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! create block_comm
    ! unpack i1-i3 from struct
    i1 = params%parallelism%i1
    i2 = params%parallelism%i2
    i3 = params%parallelism%i3

    ! check for correct number of processes
    if (i1*i2*i3 .ne. params%parallelism%world_size) then
       if (params%parallelism%world_image.eq.1) then
          write(*,*) "error in mpi_cartesian.f90:init_topology:n1*n2*n3 do not match total number of processes:",i1,i2,i3,params%parallelism%world_size
       end if

       call stop_cruna
    end if

    ! set parameters for mpi cart create
    dims     = (/i1,i2,i3/)
    periodic = .false.

    ! create cartesian comm
    call mpi_cart_create(params%parallelism%block_comm,size(dims),dims,periodic,reorder,new_comm,ierr)

    if (ierr.ne.0) then
       write(*,*) "error in mpi_cartesian.f90:init_topology:mpi_cart_create:",ierr
       call stop_cruna
    end if

    ! update parallelism block parameter (needed as reorder is true)
    call mpi_comm_size(new_comm,comm_size,ierr)
    call mpi_comm_rank(new_comm,comm_rank,ierr)

    params%parallelism%block_comm  = new_comm
    params%parallelism%block_size  = comm_size
    params%parallelism%block_image = comm_rank + 1

    call set_parameter(params%parallelism%block_comm  ,'parallelism.block_comm' )
    call set_parameter(params%parallelism%block_size  ,'parallelism.block_size' )
    call set_parameter(params%parallelism%block_image ,'parallelism.block_image')

!!! split directions, calculate number of grid points in each direction
    call mpi_cart_coords(params%parallelism%block_comm,params%parallelism%block_image-1,size(mpi_cart_coord),mpi_cart_coord,ierr)

    if (ierr.ne.0) then
       write(*,*) "error in mpi_cartesian.f90:init_topology:mpi_cart_coord:",ierr
       call stop_cruna
    end if

    ! write cart coords to struct and list (needed by e.g. boundary conditions)
    params%parallelism%i1b = mpi_cart_coord(1) + 1
    params%parallelism%i2b = mpi_cart_coord(2) + 1
    params%parallelism%i3b = mpi_cart_coord(3) + 1

    call set_parameter(params%parallelism%i1b,'parallelism.i1b')
    call set_parameter(params%parallelism%i2b,'parallelism.i2b')
    call set_parameter(params%parallelism%i3b,'parallelism.i3b')

    ! split x1-direction ---------------------------------------------------------------------------
    n1b = floor(real(params%geom%n1,rk) / real(params%parallelism%i1,rk))                           ! seperate into equal parts

    if (mod(params%geom%n1,params%parallelism%i1).ne.0) then                                        ! distribute extra points
       if (params%parallelism%block_image.eq.1) then
          write(*,*) "warning in mpi_cartesian.f90:init_topology:uneven grid distribution (i1), performance losses are possible" 
       end if

       if ((mpi_cart_coord(1)+1).gt.(params%parallelism%i1 - mod(params%geom%n1,params%parallelism%i1)) ) then
          n1b = n1b + 1
       end if
    end if

    if (n1b.lt.1) then                                                                              ! check for over parallelisation
       write(*,*) "error in mpi_cartesian.f90:init_topology:geom.n1 < parallelism.i1" 
       call stop_cruna
    end if

    ! split x2-direction ---------------------------------------------------------------------------
    n2b = floor(real(params%geom%n2,rk) / real(params%parallelism%i2,rk))

    if (mod(params%geom%n2,params%parallelism%i2).ne.0) then                                        
       if (params%parallelism%block_image.eq.1) then
          write(*,*) "warning in mpi_cartesian.f90:init_topology:uneven grid distribution (i2), performance losses are possible" 
       end if

       if ((mpi_cart_coord(2)+1).gt.(params%parallelism%i2 - mod(params%geom%n2,params%parallelism%i2)) ) then
          n2b = n2b + 1
       end if
    end if

    if (n2b.lt.1) then 
       write(*,*) "error in mpi_cartesian.f90:init_topology:geom.n2 < parallelism.i2" 
       call stop_cruna
    end if

    ! split x3-direction ---------------------------------------------------------------------------
    n3b = floor(real(params%geom%n3,rk) / real(params%parallelism%i3,rk))

    if (mod(params%geom%n3,params%parallelism%i3).ne.0) then                                        
       if (params%parallelism%block_image.eq.1) then
          write(*,*) "warning in mpi_cartesian.f90:init_topology:uneven grid distribution (i3), performance losses are possible" 
       end if

       if ((mpi_cart_coord(3)+1).gt.(params%parallelism%i3 - mod(params%geom%n3,params%parallelism%i3)) ) then
          n3b = n3b + 1
       end if
    end if

    if (n3b.lt.1) then 
       write(*,*) "error in mpi_cartesian.f90:init_topology:geom.n3 < parallelism.i3" 
       call stop_cruna
    end if

    ! write(*,*) "nXl:",n1b,n2b,n3b,params%parallelism%block_image,mpi_cart_coord

    ! write local block resolutions to struct and parameter list
    params%geom%n1b = n1b
    params%geom%n2b = n2b
    params%geom%n3b = n3b

    call set_parameter(n1b,'geom.n1b')
    call set_parameter(n2b,'geom.n2b')
    call set_parameter(n3b,'geom.n3b')

!!! create pencil communicators (for tranpose operations)
    ! i1 pencil comm
    call mpi_cart_sub(params%parallelism%block_comm,(/.true.,.false.,.false./),new_comm,ierr)
    call mpi_comm_size(new_comm,comm_size,ierr)
    call mpi_comm_rank(new_comm,comm_rank,ierr)

    params%parallelism%i1pen_comm = new_comm
    params%parallelism%i1pen_size = comm_size
    params%parallelism%i1pen_rank = comm_rank

    call set_parameter(params%parallelism%i1pen_comm,'parallelism.i1pen_comm')
    call set_parameter(params%parallelism%i1pen_size,'parallelism.i1pen_size')
    call set_parameter(params%parallelism%i1pen_rank,'parallelism.i1pen_rank')

    ! i2 pencil comm
    call mpi_cart_sub(params%parallelism%block_comm,(/.false.,.true.,.false./),new_comm,ierr)
    call mpi_comm_size(new_comm,comm_size,ierr)
    call mpi_comm_rank(new_comm,comm_rank,ierr)

    params%parallelism%i2pen_comm = new_comm
    params%parallelism%i2pen_size = comm_size
    params%parallelism%i2pen_rank = comm_rank

    call set_parameter(params%parallelism%i2pen_comm,'parallelism.i2pen_comm')
    call set_parameter(params%parallelism%i2pen_size,'parallelism.i2pen_size')
    call set_parameter(params%parallelism%i2pen_rank,'parallelism.i2pen_rank')

    ! i3 pencil comm
    call mpi_cart_sub(params%parallelism%block_comm,(/.false.,.false.,.true./),new_comm,ierr)
    call mpi_comm_size(new_comm,comm_size,ierr)
    call mpi_comm_rank(new_comm,comm_rank,ierr)

    params%parallelism%i3pen_comm = new_comm
    params%parallelism%i3pen_size = comm_size
    params%parallelism%i3pen_rank = comm_rank

    call set_parameter(params%parallelism%i3pen_comm,'parallelism.i3pen_comm')
    call set_parameter(params%parallelism%i3pen_size,'parallelism.i3pen_size')    
    call set_parameter(params%parallelism%i3pen_rank,'parallelism.i3pen_rank')    

  end subroutine init_topology
!!!=================================================================================================


!!!=================================================================================================
  subroutine spread_parameter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! spread the parameter array to all images (world)
    call spread_parameter_array_world

    ! read parameter files for each block
    call init_parameter_array_block

    ! spread the updated parameter array to all images (block)
    call spread_parameter_array_block

  end subroutine spread_parameter
!!!=================================================================================================

!!!=================================================================================================
  subroutine spread_boundary_conditions_array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! spread the parameter array to all images (world)
    call spread_boundary_conditions_array_world

  end subroutine spread_boundary_conditions_array
!!!=================================================================================================

!!!=================================================================================================
  subroutine spread_parameter_array_world
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                               :: ierr, size_parameter_array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! calculate parameter array size
    size_parameter_array = size(parameter_array) * max_length_parameter

    ! broadcast parameter_array
    call mpi_bcast(parameter_array, size_parameter_array, mpi_character, 0, params%parallelism%world_comm, ierr)

    if (ierr.ne.0) then
       write(*,*) "error in mpi_cartesian.f90:spread_parameter_array_world:mpi_bcast" 
       call stop_cruna
    end if

    ! add/create parameterlist entries
    !    just to have a parallelism documentation in all output files
    call set_parameter(params%parallelism%type        ,'parallelism.type'       )

    call set_parameter(params%parallelism%world_comm  ,'parallelism.world_comm' )
    call set_parameter(params%parallelism%world_size  ,'parallelism.world_size' )
    call set_parameter(params%parallelism%world_image ,'parallelism.world_image')
    call set_parameter(params%parallelism%world_block ,'parallelism.world_block')

    call set_parameter(params%parallelism%block_comm  ,'parallelism.block_comm' )
    call set_parameter(params%parallelism%block_size  ,'parallelism.block_size' )
    call set_parameter(params%parallelism%block_image ,'parallelism.block_image')

  end subroutine spread_parameter_array_world
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_parameter_array_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! nothing to do here
  end subroutine init_parameter_array_block
!!!=================================================================================================

!!!=================================================================================================
  subroutine spread_parameter_array_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! nothing to do here
  end subroutine spread_parameter_array_block
!!!=================================================================================================

!!!=================================================================================================
  subroutine spread_boundary_conditions_array_world
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                               :: s1,s2
    integer                                               :: ierr, size_boundary_conditions_array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! broadcast size of boundary conditions array
    s1 = size(boundary_conditions_array,1)
    s2 = size(boundary_conditions_array,2)
    call mpi_bcast(s1, 1, mpi_integer, 0, params%parallelism%world_comm, ierr)
    call mpi_bcast(s2, 1, mpi_integer, 0, params%parallelism%world_comm, ierr)

    ! allocate boundary_conditions_array (all but rank 0)
    if (params%parallelism%world_image.ne.1) then
       allocate(boundary_conditions_array(s1,s2))
       boundary_conditions_array = 0.0_rk
    end if

    ! calculate parameter array size
    size_boundary_conditions_array = size(boundary_conditions_array)

    ! broadcast parameter_array
    call mpi_bcast(boundary_conditions_array, size_boundary_conditions_array, params%parallelism%mpi_type_rk, 0, params%parallelism%world_comm, ierr)

    if (ierr.ne.0) then
       write(*,*) "error in mpi_cartesian.f90:spread_boundary_conditions_array_world:mpi_bcast" 
       call stop_cruna
    end if

  end subroutine spread_boundary_conditions_array_world
!!!=================================================================================================

!!!=================================================================================================
  subroutine stop_cruna
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                               :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call mpi_barrier(params%parallelism%world_comm,ierr)
    call mpi_finalize(ierr)
    stop

  end subroutine stop_cruna
!!!=================================================================================================

!!!=================================================================================================
  subroutine end_parallelism
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                               :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call mpi_barrier(params%parallelism%world_comm,ierr)
    call mpi_finalize(ierr)

  end subroutine end_parallelism
!!!=================================================================================================

!!!=================================================================================================
  subroutine allreduce_rk(val_out, val_in, operation, comm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , intent(inout) :: val_out
    real(kind=rk)   , intent(inout) :: val_in
    character(len=*), intent(in)    :: operation
    integer         , intent(in)    :: comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                   :: val_in_buffer
    integer                         :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    val_in_buffer = val_in          ! to avoid aliasing of mpi_buffers

    select case (trim(adjustl(operation)))

    case ('max')
       call mpi_allreduce(val_in_buffer,val_out,1,params%parallelism%mpi_type_rk,mpi_max,comm,ierr)

    case ('min')
       call mpi_allreduce(val_in_buffer,val_out,1,params%parallelism%mpi_type_rk,mpi_min,comm,ierr)

    case ('sum')
       call mpi_allreduce(val_in_buffer,val_out,1,params%parallelism%mpi_type_rk,mpi_sum,comm,ierr)

    case default
       write(*,*) "error in mpi_cartesian.f90:allreduce_rk:unknown operation, stop" 
       call stop_cruna

    end select

  end subroutine allreduce_rk
!!!=================================================================================================

!!!=================================================================================================
  subroutine allreduce_rk_loc(val_out, val_in, operation, comm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , intent(out) :: val_out(2)
    real(kind=rk)   , intent(in)  :: val_in(2)
    character(len=*), intent(in)  :: operation
    integer         , intent(in)  :: comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                       :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    select case (trim(adjustl(operation)))

    case ('maxloc')
!       call mpi_allreduce(val_in,val_out,2,params%parallelism%mpi_type_rk,mpi_maxloc,comm,ierr)
       call mpi_allreduce(val_in,val_out,2,params%parallelism%mpi_type_2rk,mpi_maxloc,comm,ierr)

    case ('minloc')
!       call mpi_allreduce(val_in,val_out,2,params%parallelism%mpi_type_rk,mpi_minloc,comm,ierr)
       call mpi_allreduce(val_in,val_out,2,params%parallelism%mpi_type_2rk,mpi_minloc,comm,ierr)

    case default
       write(*,*) "error in mpi_cartesian.f90:allreduce_rk_loc:unknown operation, stop" 
       call stop_cruna

    end select

  end subroutine allreduce_rk_loc
!!!=================================================================================================


!!!=================================================================================================
  subroutine allreduce_ik(val_out, val_in, operation, comm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer         , intent(inout) :: val_out
    integer         , intent(inout) :: val_in
    character(len=*), intent(in)    :: operation
    integer         , intent(in)    :: comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                         :: val_in_buffer
    integer                         :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    val_in_buffer = val_in          ! to avoid aliasing of mpi_buffers

    select case (trim(adjustl(operation)))

    case ('max')
       call mpi_allreduce(val_in_buffer,val_out,1,params%parallelism%mpi_type_ik,mpi_max,comm,ierr)

    case ('min')
       call mpi_allreduce(val_in_buffer,val_out,1,params%parallelism%mpi_type_ik,mpi_min,comm,ierr)

    case ('sum')
       call mpi_allreduce(val_in_buffer,val_out,1,params%parallelism%mpi_type_ik,mpi_sum,comm,ierr)

    case default
       write(*,*) "error in mpi_cartesian.f90:allreduce_ik:unknown operation, stop" 
       call stop_cruna

    end select

  end subroutine allreduce_ik
!!!=================================================================================================


!!!=================================================================================================
  subroutine bcast_ik(val, root, comm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer         , intent(inout)  :: val
    integer         , intent(in)     :: root
    integer         , intent(in)     :: comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                          :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call mpi_bcast(val,1,params%parallelism%mpi_type_ik,root,comm,ierr)

  end subroutine bcast_ik
!!!=================================================================================================


!!!=================================================================================================
  subroutine bcast_ik_1d(val, root, comm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:)         , intent(inout)  :: val
    integer         , intent(in)                   :: root
    integer         , intent(in)                   :: comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                        :: ierr
    integer                                        :: n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    n = size(val)
    call mpi_bcast(val,n,params%parallelism%mpi_type_ik,root,comm,ierr)

  end subroutine bcast_ik_1d
!!!=================================================================================================

!!!=================================================================================================
  subroutine bcast_rk(val, root, comm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , intent(inout)  :: val
    integer         , intent(in)     :: root
    integer         , intent(in)     :: comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                          :: ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call mpi_bcast(val,1,params%parallelism%mpi_type_rk,root,comm,ierr)

  end subroutine bcast_rk
!!!=================================================================================================


!!!=================================================================================================
  subroutine bcast_rk_1d(val, root, comm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:)  , intent(inout)  :: val
    integer                      , intent(in)     :: root
    integer                      , intent(in)     :: comm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                       :: ierr
    integer                                       :: n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    n = size(val)
    call mpi_bcast(val,n,params%parallelism%mpi_type_rk,root,comm,ierr)

  end subroutine bcast_rk_1d
!!!=================================================================================================


end module mpi_cartesian
