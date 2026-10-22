!> Central parameter storage, conversion, and runtime access module.
!!
!! This module defines the project-wide parameter precision, structured runtime
!! parameter types, the backing character-based parameter table, and generic
!! interfaces for reading and updating parameter values.
!!
!! Parameters may originate from the parameter file, command-line selection of that
!! file, defaults supplied by callers, or direct updates through `set_parameter`.
!! The structured `params` object provides typed access to frequently used settings,
!! while `parameter_array` retains the complete textual representation.
module parameter

  USE iso_fortran_env , only: real32, real64
  use parameter_file_io
  use precision , only: rk

  private

!!! general parameter
  ! parameter definitions
  integer, parameter, public       :: max_length_parameter = 384
  integer, parameter, public       :: max_nbr_parameters   = 512

  integer, parameter, public       :: max_length_fname     = 384

  ! variables types
  ! Project-wide floating-point precision used throughout the parameter system.
  ! real32 (single) or real64 (double, default)
  public                           :: rk ! to be set in precision.f90 

!!! parameter array & struct
  ! array
  ! table: first column stores names; remaining columns store values/subindices.
  character(len=max_length_parameter), dimension(max_nbr_parameters,5), save, public :: parameter_array

  ! struct
  !> Equation-related runtime parameters.
  !!
  !! Stores the number of dependent variables used by the active equation system.
  type equation
     integer :: nbr_vars
  end type equation

  !> Geometry and local-grid parameters.
  !!
  !! Contains grid spacings, global grid sizes, and block-local grid sizes for the
  !! three spatial coordinate directions.
  type geom
     real(kind=rk) :: dx1,dx2,dx3
     integer       :: n1 ,n2 ,n3
     integer       :: n1b,n2b,n3b
  end type geom

  !> Initial-condition parameters.
  !!
  !! Stores primitive initial-state quantities.
  type init
     real(kind=rk) :: rho,u1,u2,u3,p
  end type init

  !> Input/output scheduling and backend parameters.
  !!
  !! Contains output frequencies, optional start/end limits, verbosity, and MPI-I/O
  !! configuration.
  type io
     integer :: sfreq, afreq, dfreq, asfreq, dsfreq, aslfreq, dslfreq, gfreq, ffreq
     integer :: afreq_start, dfreq_start, asfreq_start, dsfreq_start, aslfreq_start, dslfreq_start, gfreq_start, ffreq_start
     integer :: afreq_end  , dfreq_end  , asfreq_end  , dsfreq_end  , aslfreq_end  , dslfreq_end  , gfreq_end  , ffreq_end
     integer :: verbosity
     logical :: mpiio, mpiio_single_precision
  end type io

  !> Material-property parameters.
  !!
  !! Stores thermodynamic material constants.
  type material
     real(kind=rk) :: gamma
     real(kind=rk) :: R
  end type material

  !> Parallel decomposition and communicator parameters.
  !!
  !! Stores communicator metadata, MPI datatype handles, and decomposition indices.
  type parallelism
     character(len=max_length_parameter) :: type
     integer                             :: world_comm, world_size, world_image, world_block, block_comm, block_size, block_image
     integer                             :: i1pen_comm, i1pen_size, i1pen_rank, i2pen_comm, i2pen_size, i2pen_rank, i3pen_comm, i3pen_size, i3pen_rank
     integer                             :: mpi_type_rk, mpi_type_2rk, mpi_type_ik, mpi_type_ck
     integer                             :: i1 ,i2 ,i3
     integer                             :: i1b,i2b,i3b
  end type parallelism

  !> Time-integration parameters.
  !!
  !! Stores timestep size, physical time, subset state, and timestep counters.
  type time
     real(kind=rk)     :: dt,t
     integer           :: ns, subsets
     integer           :: nt, steps
  end type time

  !> Aggregated typed runtime parameter structure.
  !!
  !! Groups the individual parameter categories into the global `params` instance.
  type params_struct
     type(equation)    :: equation
     type(geom)        :: geom
     type(init)        :: init
     type(io)          :: io
     type(material)    :: material
     type(parallelism) :: parallelism
     type(time)        :: time
  end type params_struct

  ! Global typed parameter structure.
  type(params_struct), save, public :: params

!!! methods of the module
  public  :: init_parameter_array
  public  :: init_parameter_struct
  public  :: create_parameter_list
  public  :: get_parameter
  public  :: set_parameter

  private :: get_parameter_index

!!! interfaces of the methods
  !> Generic typed parameter retrieval interface.
  !!
  !! Resolves to character, integer, logical, or `REAL(kind=rk)` implementations.
  interface get_parameter
     module procedure get_parameter_c
     module procedure get_parameter_i
     module procedure get_parameter_l
     module procedure get_parameter_rk
  end interface get_parameter

  !> Generic typed parameter update interface.
  !!
  !! Resolves to character, integer, logical, or `REAL(kind=rk)` implementations.
  interface set_parameter
     module procedure set_parameter_c
     module procedure set_parameter_i
     module procedure set_parameter_l
     module procedure set_parameter_rk
  end interface set_parameter

  ! All routines below are module procedures.
contains

!!!=================================================================================================
  !> Populate the typed runtime parameter structure from the parameter table.
  !!
  !! Reads parallelism, equation, geometry, initial-condition, I/O, material, and
  !! time-integration settings into the global `params` structure.
  subroutine init_parameter_struct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! parallelism
    !    type, world_comm, world_image, world_block, block_comm, block_image, mpi_type_rk are already set in init_parallelism
    ! Populate decomposition and communicator-related runtime parameters.
    call get_parameter(params%parallelism%i1 ,"parallelism.i1" , default = 1)
    call get_parameter(params%parallelism%i2 ,"parallelism.i2" , default = 1)
    call get_parameter(params%parallelism%i3 ,"parallelism.i3" , default = 1)

    call get_parameter(params%parallelism%i1b,"parallelism.i1b", default = 1)
    call get_parameter(params%parallelism%i2b,"parallelism.i2b", default = 1)
    call get_parameter(params%parallelism%i3b,"parallelism.i3b", default = 1)

    call get_parameter(params%parallelism%i1pen_comm,"parallelism.i1pen_comm",default = 0)
    call get_parameter(params%parallelism%i2pen_comm,"parallelism.i2pen_comm",default = 0)
    call get_parameter(params%parallelism%i3pen_comm,"parallelism.i3pen_comm",default = 0)

    call get_parameter(params%parallelism%i1pen_size,"parallelism.i1pen_size",default = 1)
    call get_parameter(params%parallelism%i2pen_size,"parallelism.i2pen_size",default = 1)
    call get_parameter(params%parallelism%i3pen_size,"parallelism.i3pen_size",default = 1)

    call get_parameter(params%parallelism%i1pen_rank,"parallelism.i1pen_rank",default = 1)
    call get_parameter(params%parallelism%i2pen_rank,"parallelism.i2pen_rank",default = 1)
    call get_parameter(params%parallelism%i3pen_rank,"parallelism.i3pen_rank",default = 1)

    ! equations
    ! Populate equation configuration.
    call get_parameter(params%equation%nbr_vars, "equation.nbr_vars")

    ! geometry
    ! Populate global and block-local geometry parameters.
    call get_parameter(params%geom%n1 , "geom.n1")
    call get_parameter(params%geom%n2 , "geom.n2")
    call get_parameter(params%geom%n3 , "geom.n3")

    call get_parameter(params%geom%n1b, "geom.n1b", default = 0)
    call get_parameter(params%geom%n2b, "geom.n2b", default = 0)
    call get_parameter(params%geom%n3b, "geom.n3b", default = 0)

    call get_parameter(params%geom%dx1, "geom.dx1", default = 0.0_rk)
    call get_parameter(params%geom%dx2, "geom.dx2", default = 0.0_rk)
    call get_parameter(params%geom%dx3, "geom.dx3", default = 0.0_rk)

    ! init
    ! Populate initial-condition parameters.
    call get_parameter(params%init%rho             , "init.rho"                                           )
    call get_parameter(params%init%p               , "init.p"                                             )
    call get_parameter(params%init%u1              , "init.u1"               , default = 0.0_rk           )
    call get_parameter(params%init%u2              , "init.u2"               , default = 0.0_rk           )
    call get_parameter(params%init%u3              , "init.u3"               , default = 0.0_rk           )

    ! io
    ! Populate output frequencies and I/O options.
    call get_parameter(params%io%sfreq  , "io.sfreq"  , default = 1      )
    call get_parameter(params%io%afreq  , "io.afreq"  , default = huge(1))
    call get_parameter(params%io%asfreq , "io.asfreq" , default = huge(1))
    call get_parameter(params%io%aslfreq, "io.aslfreq", default = huge(1))
    call get_parameter(params%io%dfreq  , "io.dfreq"  , default = huge(1))
    call get_parameter(params%io%dsfreq , "io.dsfreq" , default = huge(1))
    call get_parameter(params%io%dslfreq, "io.dslfreq", default = huge(1))
    call get_parameter(params%io%ffreq  , "io.ffreq"  , default = huge(1))
    call get_parameter(params%io%gfreq  , "io.gfreq"  , default = huge(1))

    call get_parameter(params%io%mpiio                   , "io.mpiio"                   , default = .false.)
    call get_parameter(params%io%mpiio_single_precision  , "io.mpiio_single_precision"  , default = .false.)

    ! io skip freq, when to start
    ! Populate optional start indices for periodic output actions.
    call get_parameter(params%io%afreq_start  , "io.afreq"  , 2, default = 1)
    call get_parameter(params%io%asfreq_start , "io.asfreq" , 2, default = 1)
    call get_parameter(params%io%aslfreq_start, "io.aslfreq", 2, default = 1)
    call get_parameter(params%io%dfreq_start  , "io.dfreq"  , 2, default = 1)
    call get_parameter(params%io%dsfreq_start , "io.dsfreq" , 2, default = 1)
    call get_parameter(params%io%dslfreq_start, "io.dslfreq", 2, default = 1)
    call get_parameter(params%io%ffreq_start  , "io.ffreq"  , 2, default = 1)
    call get_parameter(params%io%gfreq_start  , "io.gfreq"  , 2, default = 1)

    ! Populate optional end indices for periodic output actions.
    call get_parameter(params%io%afreq_end  , "io.afreq"  , 3, huge(1))
    call get_parameter(params%io%asfreq_end , "io.asfreq" , 3, huge(1))
    call get_parameter(params%io%aslfreq_end, "io.aslfreq", 3, huge(1))
    call get_parameter(params%io%dfreq_end  , "io.dfreq"  , 3, huge(1))
    call get_parameter(params%io%dsfreq_end , "io.dsfreq" , 3, huge(1))
    call get_parameter(params%io%dslfreq_end, "io.dslfreq", 3, huge(1))
    call get_parameter(params%io%ffreq_end  , "io.ffreq"  , 3, huge(1))
    call get_parameter(params%io%gfreq_end  , "io.gfreq"  , 3, huge(1))

    call get_parameter(params%io%verbosity, "io.verbosity", default = 0)

    ! material
    ! Populate material properties.
    call get_parameter(params%material%gamma, "material.gamma",default = 1.4_rk)
    call get_parameter(params%material%R    , "material.R", default = 8.3144598_rk)

    ! time
    ! Populate time-integration settings.
    call get_parameter(params%time%steps  , "time.steps"                     )
    call get_parameter(params%time%subsets, "time.subsets", default = 1      )
    call get_parameter(params%time%t      , "time.t"      , default = 0.0_rk )
    call get_parameter(params%time%dt     , "time.dt"                        )

  end subroutine init_parameter_struct
!!!=================================================================================================




!!!=================================================================================================
  !> Initialize the textual parameter table from the selected parameter file.
  !!
  !! Parses the optional `-p <file>` command-line argument, falls back to
  !! `parameter.dat`, reads the parameter file on the root image, transfers entries
  !! into `parameter_array`, and records precision metadata.
  subroutine init_parameter_array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_parameter), dimension(:,:), allocatable     :: parameter_list
    real(kind=rk)                                                        :: prec_rk = 0
    integer                                                              :: prec_i  = 0
    integer                                                              :: i,j
    integer                                                              :: a,arg_count
    character(len=max_length_parameter)                                  :: arg_val
    character(len=max_length_fname)                                      :: parameter_fname
    logical                                                              :: p_arg_found = .false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! eval command line arguments for user-defined parameter-file
    ! just "-p" is evaluated and implemented as all other parameters should be set within the parameter file
    ! so do not implement other command line argument stuff
    ! Inspect command-line arguments for a parameter-file override.
    arg_count = command_argument_count()
    if (arg_count.gt.0) then

       if (mod(arg_count,2).ne.0) then ! wrong format or syntax, must be "-a <arg>"
          write(*,*) "error in parameter.f90:init_parameter_array:command line arguments syntax error"
          stop
       end if

       ! look for "-p <arg>"
       ! Scan option/value pairs for the supported `-p` argument.
       do a = 1,arg_count,2
          call get_command_argument(a,arg_val)
          if (trim(arg_val).eq."-p") then
             p_arg_found = .true.

             call get_command_argument(a+1,arg_val)
             parameter_fname = trim(arg_val)
          end if
       end do
    end if

    ! set default parameter fname
    ! Fall back to the conventional parameter filename.
    if (p_arg_found.eqv..false.) then
       parameter_fname = "parameter.dat"
    end if

    ! init parameter array
    ! Reset the persistent parameter table before importing file contents.
    parameter_array = ""

    ! just done for image 1
    if (params%parallelism%world_image.eq.1) then
       ! allocate
       allocate(parameter_list(max_nbr_parameters,5))

       ! read parameter.dat
       ! Read the selected parameter file into a temporary table.
       call read_parameter_file(parameter_list,parameter_fname)

       ! transfer to parameter_array
       ! Transfer populated parameter entries and subindices into the persistent table.
       do i = 1,max_nbr_parameters
          if (len_trim(parameter_list(i,1)).gt.0) then
             ! debug:
             ! write (*,*) "write ",trim(parameter_list(i,1))," to parameter array ",trim(parameter_list(i,2))

             do j = 1,5 ! parameter name subindex=0 AND arguments are set subindex>0
                call set_parameter(trim(parameter_list(i,j)),trim(parameter_list(i,1)),subindex = j-1)
             end do
          end if
       end do

    end if

    ! debug:
    ! write(*,*) "array: ",parameter_array(1,2)

    ! add precision stuff to parameter list
    ! Record precision metadata in the parameter table.
    call set_parameter(rk           ,'precision.rk')
    call set_parameter(tiny(prec_rk),'precision.tiny_rk')
    call set_parameter(huge(prec_rk),'precision.huge_rk')
    call set_parameter(huge(prec_i) ,'precision.huge_int')

  end subroutine init_parameter_array
!!!=================================================================================================





!!!=================================================================================================
  !> Store a character parameter value.
  !!
  !! **Arguments**
  !!
  !! * `val` - Character value to store.
  !! * `name` - Parameter name used as lookup key.
  !! * `subindex` - Optional logical subindex.
  !! * `struct_val` - Optional structured field updated with the same value.
  subroutine set_parameter_c(val, name, subindex, struct_val)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*), intent(in)              :: val
    character(len=*), intent(in)              :: name
    integer         , intent(in)   , optional :: subindex
    character(len=*), intent(inout), optional :: struct_val
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                   :: i1,i2,i3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional subindex
    if (present(subindex)) then
       i3 = subindex + 1
    else ! first/default is 2
       i3 = 2
    end if

    ! if parameter part of params structure struct_val is given
    if (present(struct_val)) then
       struct_val = val
    end if

    ! debug:
    ! write(*,*) "add: ", trim(name),trim(val),i3 

    ! check for already present parameter
    call get_parameter_index(i1,trim(name))

    if (i1.eq.0) then
       !! new parameter

       do i2 = 1,max_nbr_parameters
          ! check for free index
          if (trim(parameter_array(i2,1)).eq."") then
             parameter_array(i2,1 ) = trim(name)
             parameter_array(i2,i3) = trim(val )
             exit
          end if
       end do

    else
       !! already present parameter
       parameter_array(i1,1 ) = trim(name)
       parameter_array(i1,i3) = trim(val )
    end if

  end subroutine set_parameter_c
!!!=================================================================================================

!!!=================================================================================================
  !> Store an integer parameter value.
  !!
  !! **Arguments**
  !!
  !! * `val` - Integer value to store.
  !! * `name` - Parameter name used as lookup key.
  !! * `subindex` - Optional logical subindex.
  !! * `struct_val` - Optional structured field updated with the same value.
  subroutine set_parameter_i(val, name, subindex, struct_val)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer         , intent(in)              :: val
    character(len=*), intent(in)              :: name
    integer         , intent(in)   , optional :: subindex
    integer         , intent(inout), optional :: struct_val
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_parameter)       :: val_c
    integer                                   :: i1,i2,i3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional subindex
    if (present(subindex)) then
       i3 = subindex + 1
    else ! first/default is 2
       i3 = 2
    end if

    ! if parameter part of params structure struct_val is given
    if (present(struct_val)) then
       struct_val = val
    end if

    write(unit=val_c,fmt=*) val

    ! check for already present parameter
    call get_parameter_index(i1,trim(name))

    if (i1.eq.0) then
       !! new parameter

       do i2 = 1,max_nbr_parameters
          ! check for free index
          if (trim(parameter_array(i2,1)).eq."") then
             parameter_array(i2,1 ) = trim(name)
             parameter_array(i2,i3) = trim(val_c)
             exit
          end if
       end do

    else
       !! already present parameter
       parameter_array(i1,1 ) = trim(name)
       parameter_array(i1,i3) = trim(val_c)
    end if

  end subroutine set_parameter_i
!!!=================================================================================================

!!!=================================================================================================
  !> Store a logical parameter value.
  !!
  !! **Arguments**
  !!
  !! * `val` - Logical value to store.
  !! * `name` - Parameter name used as lookup key.
  !! * `subindex` - Optional logical subindex.
  !! * `struct_val` - Optional structured field updated with the same value.
  subroutine set_parameter_l(val, name, subindex, struct_val)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical         , intent(in)              :: val
    character(len=*), intent(in)              :: name
    integer         , intent(in)   , optional :: subindex
    logical         , intent(inout), optional :: struct_val
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_parameter)       :: val_c
    integer                                   :: i1,i2,i3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional subindex
    if (present(subindex)) then
       i3 = subindex + 1
    else ! first/default is 2
       i3 = 2
    end if

    ! if parameter part of params structure struct_val is given
    if (present(struct_val)) then
       struct_val = val
    end if

    write(unit=val_c,fmt=*) val

    ! check for already present parameter
    call get_parameter_index(i1,trim(name))

    if (i1.eq.0) then
       !! new parameter

       do i2 = 1,max_nbr_parameters
          ! check for free index
          if (trim(parameter_array(i2,1)).eq."") then
             parameter_array(i2,1 ) = trim(name)
             parameter_array(i2,i3) = trim(val_c)
             exit
          end if
       end do

    else
       !! already present parameter
       parameter_array(i1,1 ) = trim(name)
       parameter_array(i1,i3) = trim(val_c)
    end if

  end subroutine set_parameter_l
!!!=================================================================================================  



!!!=================================================================================================
  !> Store a `REAL(kind=rk)` parameter value.
  !!
  !! **Arguments**
  !!
  !! * `val` - Real value to store.
  !! * `name` - Parameter name used as lookup key.
  !! * `subindex` - Optional logical subindex.
  !! * `struct_val` - Optional structured field updated with the same value.
  subroutine set_parameter_rk(val, name, subindex, struct_val)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , intent(in)              :: val
    character(len=*), intent(in)              :: name
    integer         , intent(in)   , optional :: subindex
    real(kind=rk)   , intent(inout), optional :: struct_val
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_parameter)       :: val_c
    integer                                   :: i1,i2,i3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional subindex
    if (present(subindex)) then
       i3 = subindex + 1
    else ! first/default is 2
       i3 = 2
    end if

    ! if parameter part of params structure struct_val is given
    if (present(struct_val)) then
       struct_val = val
    end if

    write(unit=val_c,fmt=*) val

    ! check for already present parameter
    call get_parameter_index(i1,trim(name))

    if (i1.eq.0) then
       !! new parameter

       do i2 = 1,max_nbr_parameters
          ! check for free index
          if (trim(parameter_array(i2,1)).eq."") then
             parameter_array(i2,1 ) = trim(name)
             parameter_array(i2,i3) = trim(val_c)
             exit
          end if
       end do

    else
       !! already present parameter
       parameter_array(i1,1 ) = trim(name)
       parameter_array(i1,i3) = trim(val_c)
    end if

  end subroutine set_parameter_rk
!!!=================================================================================================





!!!=================================================================================================
  !> Retrieve a character parameter value.
  !!
  !! Uses the selected subindex when present and applies the optional default policy
  !! when the requested entry is missing or empty.
  !!
  !! **Arguments**
  !!
  !! * `val` - Retrieved character value.
  !! * `name` - Parameter name used as lookup key.
  !! * `subindex` - Optional logical subindex.
  !! * `default` - Optional fallback value.
  !! * `set_default` - Optional flag controlling persistence of the fallback.
  subroutine get_parameter_c(val, name, subindex, default, set_default)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_parameter), intent(out)          :: val
    character(len=*)                   , intent(in)           :: name
    integer                            , intent(in), optional :: subindex
    character(len=*)                   , intent(in), optional :: default
    logical                            , intent(in), optional :: set_default
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_parameter)                       :: transfer_dummy
    integer                                                   :: i1,i3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (present(subindex)) then
       i3 = subindex + 1
    else ! first/default is 2
       i3 = 2
    end if

    ! find index
    call get_parameter_index(i1,trim(name))

    if (i1.ne.0) then
       ! transfer to dummy (as read requires an integer or a character variable as unit)
       transfer_dummy = "'" // trim(parameter_array(i1,i3)) // "'"

       ! check for empty transfer_dummy, which could happen if subindex's are used, e.g.:
       ! a parameter is created with set_default for a certain index, so all other subindex's are ""
       if (len_trim(adjustl(transfer_dummy)).eq.0) then
          ! transfer dummy is empty as parameter entry(@ index) is empty
          ! check for default value
          if (present(default)) then
             val = default

             if (present(set_default)) then
                if (set_default.eqv..true.) then
                   call set_parameter(val,trim(name),i3-1)
                end if
             else
                call set_parameter(val,trim(name),i3-1)
             end if
          else
             ! no default value given, stop
             write(*,*) "params: parameter '",trim(name),"' found, but empty at index:",i3-1,", no default value given"
             write(*,*) ">> STOP <<"
             stop
          end if
       else
          ! read val to output --> this should happen most of the time (!)
          read(unit=transfer_dummy,fmt=*) val
       end if
    else
       if (present(default)) then
          val = default

          if (present(set_default)) then
             if (set_default.eqv..true.) then
                call set_parameter(val,trim(name),i3-1)
             end if
          else
             call set_parameter(val,trim(name),i3-1)
          end if
       else
          write(*,*) "params: parameter '",trim(name),"' not found, no default value given"
          write(*,*) ">> STOP <<"
          stop
       end if
    end if

  end subroutine get_parameter_c
!!!=================================================================================================

!!!=================================================================================================
  !> Retrieve an integer parameter value.
  !!
  !! Uses the selected subindex when present and applies the optional default policy
  !! when the requested entry is missing or empty.
  !!
  !! **Arguments**
  !!
  !! * `val` - Retrieved integer value.
  !! * `name` - Parameter name used as lookup key.
  !! * `subindex` - Optional logical subindex.
  !! * `default` - Optional fallback value.
  !! * `set_default` - Optional flag controlling persistence of the fallback.
  subroutine get_parameter_i(val, name, subindex, default, set_default)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer         , intent(out)             :: val
    character(len=*), intent(in)              :: name
    integer         , intent(in), optional    :: subindex
    integer         , intent(in), optional    :: default
    logical         , intent(in), optional    :: set_default
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_parameter)       :: transfer_dummy
    integer                                   :: i1,i3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (present(subindex)) then
       i3 = subindex + 1
    else ! first/default is 2
       i3 = 2
    end if

    ! find index
    call get_parameter_index(i1,trim(name))

    if (i1.ne.0) then
       ! transfer to dummy (as read requires an integer or a character variable as unit)
       transfer_dummy = trim(parameter_array(i1,i3))

       ! check for empty transfer_dummy, which could happen if subindex's are used, e.g.:
       ! a parameter is created with set_default for a certain index, so all other subindex's are ""
       if (len_trim(adjustl(transfer_dummy)).eq.0) then
          ! transfer dummy is empty as parameter entry(@ index) is empty
          ! check for default value
          if (present(default)) then
             val = default

             if (present(set_default)) then
                if (set_default.eqv..true.) then
                   call set_parameter(val,trim(name),i3-1)
                end if
             else
                call set_parameter(val,trim(name),i3-1)
             end if
          else
             ! no default value given, stop
             write(*,*) "params: parameter '",trim(name),"' found, but empty at index:",i3-1,", no default value given"
             write(*,*) ">> STOP <<"
             stop
          end if

       else
          ! read val to output --> this should happen most of the time (!)
          read(unit=transfer_dummy,fmt=*) val
       end if
    else
       if (present(default)) then
          val = default

          if (present(set_default)) then
             if (set_default.eqv..true.) then
                call set_parameter(val,trim(name),i3-1)
             end if
          else
             call set_parameter(val,trim(name),i3-1)
          end if
       else
          write(*,*) "params: parameter '",trim(name),"' not found, no default value given"
          write(*,*) ">> STOP <<"
          stop
       end if
    end if

  end subroutine get_parameter_i
!!!=================================================================================================




!!!=================================================================================================
  !> Retrieve a logical parameter value.
  !!
  !! Uses the selected subindex when present and applies the optional default policy
  !! when the requested entry is missing or empty.
  !!
  !! **Arguments**
  !!
  !! * `val` - Retrieved logical value.
  !! * `name` - Parameter name used as lookup key.
  !! * `subindex` - Optional logical subindex.
  !! * `default` - Optional fallback value.
  !! * `set_default` - Optional flag controlling persistence of the fallback.
  subroutine get_parameter_l(val, name, subindex, default, set_default)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical                            , intent(out)          :: val
    character(len=*)                   , intent(in)           :: name
    integer                            , intent(in), optional :: subindex
    logical                            , intent(in), optional :: default
    logical                            , intent(in), optional :: set_default
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_parameter)                       :: transfer_dummy
    integer                                                   :: i1,i3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (present(subindex)) then
       i3 = subindex + 1
    else ! first/default is 2
       i3 = 2
    end if

    ! find index
    call get_parameter_index(i1,trim(name))

    if (i1.ne.0) then
       ! transfer to dummy (as read requires an integer or a character variable as unit)
       transfer_dummy = "" // trim(parameter_array(i1,i3)) // ""

       ! check for empty transfer_dummy, which could happen if subindex's are used, e.g.:
       ! a parameter is created with set_default for a certain index, so all other subindex's are ""
       if (len_trim(adjustl(transfer_dummy)).eq.0) then
          ! transfer dummy is empty as parameter entry(@ index) is empty
          ! check for default value
          if (present(default)) then
             val = default

             if (present(set_default)) then
                if (set_default.eqv..true.) then
                   call set_parameter(val,trim(name),i3-1)
                end if
             else
                call set_parameter(val,trim(name),i3-1)
             end if
          else
             ! no default value given, stop
             write(*,*) "params: parameter '",trim(name),"' found, but empty at index:",i3-1,", no default value given"
             write(*,*) ">> STOP <<"
             stop
          end if
       else
          ! read val to output --> this should happen most of the time (!)
          read(unit=transfer_dummy,fmt=*) val
       end if
    else
       if (present(default)) then
          val = default

          if (present(set_default)) then
             if (set_default.eqv..true.) then
                call set_parameter(val,trim(name),i3-1)
             end if
          else
             call set_parameter(val,trim(name),i3-1)
          end if
       else
          write(*,*) "params: parameter '",trim(name),"' not found, no default value given"
          write(*,*) ">> STOP <<"
          stop
       end if
    end if

  end subroutine get_parameter_l
!!!=================================================================================================



!!!=================================================================================================
  !> Retrieve a `REAL(kind=rk)` parameter value.
  !!
  !! Uses the selected subindex when present and applies the optional default policy
  !! when the requested entry is missing or empty.
  !!
  !! **Arguments**
  !!
  !! * `val` - Retrieved real value.
  !! * `name` - Parameter name used as lookup key.
  !! * `subindex` - Optional logical subindex.
  !! * `default` - Optional fallback value.
  !! * `set_default` - Optional flag controlling persistence of the fallback.
  subroutine get_parameter_rk(val, name, subindex, default, set_default)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , intent(out)             :: val
    character(len=*), intent(in)              :: name
    integer         , intent(in), optional    :: subindex
    real(kind=rk)   , intent(in), optional    :: default
    logical         , intent(in), optional    :: set_default
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_parameter)       :: transfer_dummy
    integer                                   :: i1,i3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (present(subindex)) then
       i3 = subindex + 1
    else ! first/default is 2
       i3 = 2
    end if

    ! find index
    call get_parameter_index(i1,trim(name))

    if (i1.ne.0) then
       ! transfer to dummy (as read requires an integer or a character variable as unit)
       transfer_dummy = trim(parameter_array(i1,i3))
       ! check for empty transfer_dummy, which could happen if subindex's are used, e.g.:
       ! a parameter is created with set_default for a certain index, so all other subindex's are ""
       if (len_trim(adjustl(transfer_dummy)).eq.0) then
          ! transfer dummy is empty as parameter entry(@ index) is empty
          ! check for default value

          if (present(default)) then
             val = default

             if (present(set_default)) then
                if (set_default.eqv..true.) then
                   call set_parameter(val,trim(name),i3-1)
                end if
             else
                call set_parameter(val,trim(name),i3-1)
             end if
          else
             ! no default value given, stop
             write(*,*) "params: parameter '",trim(name),"' found, but empty at index:",i3-1,", no default value given"
             write(*,*) ">> STOP <<"
             stop
          end if

       else
          ! read val to output --> this should happen most of the time (!)
          read(unit=transfer_dummy,fmt=*) val
       end if
    else
       if (present(default)) then
          val = default

          if (present(set_default)) then
             if (set_default.eqv..true.) then
                call set_parameter(val,trim(name),i3-1)
             end if
          else
             call set_parameter(val,trim(name),i3-1)
          end if
       else
          write(*,*) "params: parameter '",trim(name),"' not found, no default value given"
          write(*,*) ">> STOP <<"
          stop
       end if
    end if

  end subroutine get_parameter_rk
!!!=================================================================================================





!!!=================================================================================================
  !> Locate a parameter in the backing parameter table.
  !!
  !! **Arguments**
  !!
  !! * `index` - Matching table row, or zero when the parameter is not found.
  !! * `name` - Parameter name to search for.
  subroutine get_parameter_index(index,name)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer         , intent(out) :: index
    character(len=*), intent(in)  :: name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                       :: i
    logical                       :: found
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! returns zero if parameter is not found

    ! Zero acts as the sentinel for a parameter that was not found.
    index = 0
    found = .false.

    ! Search linearly through all possible parameter rows.
    do i=1,max_nbr_parameters
       ! debug
       ! write(*,*) i,trim(name),trim(parameter_array(i,1))
       if (trim(parameter_array(i,1)).eq.trim(name)) then
          found = .true.
          index = i
          exit
       end if

    end do

  end subroutine get_parameter_index
!!!=================================================================================================


!!!=================================================================================================
  !> Serialize the parameter table into restart-compatible strings.
  !!
  !! **Arguments**
  !!
  !! * `parameter_list` - Output array containing serialized parameter entries.
  subroutine create_parameter_list(parameter_list)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! general constants for the parameter list (matches to file setup, used for restart)
    character(len=1), parameter                                   :: parameter_equal           = "="
    character(len=1), parameter                                   :: parameter_seperator       = ";"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_parameter*5), dimension(max_nbr_parameters), intent(out) :: parameter_list
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                                           :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Initialize serialized output entries as empty.
    parameter_list = ""

    ! Serialize each populated parameter row.
    do i = 1,max_nbr_parameters
       if(trim(parameter_array(i,1)).ne."") then
          parameter_list(i) = trim(parameter_array(i,1)) // parameter_equal &
               // trim(parameter_array(i,2)) // parameter_seperator &
               // trim(parameter_array(i,3)) // parameter_seperator &
               // trim(parameter_array(i,4)) // parameter_seperator &
               // trim(parameter_array(i,5))
       end if
    end do

  end subroutine create_parameter_list
!!!=================================================================================================

end module parameter
