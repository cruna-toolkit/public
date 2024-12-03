module parameter

  use parameter_file_io

  private

!!! general parameter
  ! parameter definitions
  integer, parameter, public       :: max_length_parameter = 384
  integer, parameter, public       :: max_nbr_parameters   = 512

  integer, parameter, public       :: max_length_fname     = 384

  ! variables types
  integer, parameter, public       :: rk = selected_real_kind(15,307) ! double:(15,307); quad:(33, 4931) 

!!! parameter array & struct
  ! array
  character(len=max_length_parameter), dimension(max_nbr_parameters,5), save, public :: parameter_array

  ! struct
  type equation
     integer :: nbr_vars
  end type equation

  type geom
     real(kind=rk) :: dx1,dx2,dx3
     integer       :: n1 ,n2 ,n3
     integer       :: n1b,n2b,n3b
  end type geom

  type init
     real(kind=rk) :: rho,u1,u2,u3,p
  end type init

  type io
     integer :: sfreq, afreq, dfreq, asfreq, dsfreq, gfreq, ffreq
     integer :: verbosity
  end type io

  type material
     real(kind=rk) :: gamma
     real(kind=rk) :: R
  end type material

  type parallelism
     character(len=max_length_parameter) :: type
     integer                             :: world_comm, world_size, world_image, world_block, block_comm, block_size, block_image
     integer                             :: i1pen_comm, i1pen_size, i1pen_rank, i2pen_comm, i2pen_size, i2pen_rank, i3pen_comm, i3pen_size, i3pen_rank
     integer                             :: mpi_type_rk, mpi_type_2rk, mpi_type_ik, mpi_type_ck
     integer                             :: i1 ,i2 ,i3
     integer                             :: i1b,i2b,i3b
  end type parallelism

  type time
     real(kind=rk)     :: dt,t
     integer           :: ns, subsets
     integer           :: nt, steps
  end type time

  type params_struct
     type(equation)    :: equation
     type(geom)        :: geom
     type(init)        :: init
     type(io)          :: io
     type(material)    :: material
     type(parallelism) :: parallelism
     type(time)        :: time
  end type params_struct

  type(params_struct), save, public :: params

!!! methods of the module
  public  :: init_parameter_array
  public  :: init_parameter_struct
  public  :: create_parameter_list
  public  :: get_parameter
  public  :: set_parameter

  private :: get_parameter_index

!!! interfaces of the methods
  interface get_parameter
     module procedure get_parameter_c
     module procedure get_parameter_i
     module procedure get_parameter_l
     module procedure get_parameter_rk
  end interface get_parameter

  interface set_parameter
     module procedure set_parameter_c
     module procedure set_parameter_i
     module procedure set_parameter_l
     module procedure set_parameter_rk
  end interface set_parameter

contains 

!!!=================================================================================================
  subroutine init_parameter_struct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! parallelism
    !    type, world_comm, world_image, world_block, block_comm, block_image, mpi_type_rk are already set in init_parallelism
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
    call get_parameter(params%equation%nbr_vars, "equation.nbr_vars")

    ! geometry
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
    call get_parameter(params%init%rho             , "init.rho"                                           )
    call get_parameter(params%init%p               , "init.p"                                             )
    call get_parameter(params%init%u1              , "init.u1"               , default = 0.0_rk           )
    call get_parameter(params%init%u2              , "init.u2"               , default = 0.0_rk           )
    call get_parameter(params%init%u3              , "init.u3"               , default = 0.0_rk           )

    ! io
    call get_parameter(params%io%sfreq , "io.sfreq" , default = 1      )
    call get_parameter(params%io%afreq , "io.afreq" , default = huge(1))
    call get_parameter(params%io%asfreq, "io.asfreq", default = huge(1))
    call get_parameter(params%io%dfreq , "io.dfreq" , default = huge(1))
    call get_parameter(params%io%dsfreq, "io.dsfreq", default = huge(1))
    call get_parameter(params%io%ffreq , "io.ffreq" , default = huge(1))
    call get_parameter(params%io%gfreq , "io.gfreq" , default = huge(1))

    call get_parameter(params%io%verbosity, "io.verbosity", default = 0)

    ! material
    call get_parameter(params%material%gamma, "material.gamma")
    call get_parameter(params%material%R    , "material.R", default = 8.3144598_rk)

    ! time
    call get_parameter(params%time%steps  , "time.steps"                     )
    call get_parameter(params%time%subsets, "time.subsets", default = 1      )
    call get_parameter(params%time%t      , "time.t"      , default = 0.0_rk )
    call get_parameter(params%time%dt     , "time.dt"                        )

  end subroutine init_parameter_struct
!!!=================================================================================================




!!!=================================================================================================
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
    arg_count = command_argument_count()
    if (arg_count.gt.0) then

       if (mod(arg_count,2).ne.0) then ! wrong format or syntax, must be "-a <arg>"
          write(*,*) "error in parameter.f90:init_parameter_array:command line arguments syntax error"
          stop
       end if

       ! look for "-p <arg>"
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
    if (p_arg_found.eqv..false.) then
       parameter_fname = "parameter.dat"
    end if

    ! init parameter array
    parameter_array = ""

    ! just done for image 1
    if (params%parallelism%world_image.eq.1) then
       ! allocate
       allocate(parameter_list(max_nbr_parameters,5))

       ! read parameter.dat
       call read_parameter_file(parameter_list,parameter_fname)

       ! transfer to parameter_array
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
    call set_parameter(rk           ,'precision.rk')
    call set_parameter(tiny(prec_rk),'precision.tiny_rk')
    call set_parameter(huge(prec_rk),'precision.huge_rk')
    call set_parameter(huge(prec_i) ,'precision.huge_int')

  end subroutine init_parameter_array
!!!=================================================================================================





!!!=================================================================================================
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
  subroutine get_parameter_index(index,name)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer         , intent(out) :: index
    character(len=*), intent(in)  :: name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                       :: i
    logical                       :: found
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! returns zero if parameter is not found

    index = 0
    found = .false.

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

    parameter_list = ""

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
