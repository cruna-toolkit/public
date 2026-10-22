!> General-purpose helper routines used throughout CRUNA.
!!
!! This module provides coordinate conversion, global extrema lookup, evenly spaced
!! vector generation, structured grid construction, numeric-to-string conversion,
!! timing utilities, mail helpers, shell-command execution, and a text progress bar.
!!
!! Generic interfaces expose type-specific implementations where appropriate.
module helper

  USE iso_fortran_env , only: real32, real64
  use parallelism , only: allreduce, bcast, stop_cruna
  use parameter

  private

  public cart2pol
  public command
  public get_index
  public linspace
  public meshgrid
  public ndgrid
  public num2str
  public progressbar
  public tic
  public toc
  public send_mail

  !! Persistent clock value used by `tic` and `toc` when no explicit timestamp is supplied.
  integer, save, private :: tic_persistent
  !! Fixed output length used by all `num2str` implementations.
  integer, parameter     :: num2str_len = 32

  !> Generic global-extremum index interface.
  !!
  !! Resolves to the `REAL(kind=rk)` implementation.
  interface get_index
     module procedure get_index_rk
  end interface get_index

  !> Generic numeric-to-string conversion interface.
  !!
  !! Resolves the concrete implementation from the numeric argument type.
  interface num2str
     module procedure num2str_i
     module procedure num2str_r32
     module procedure num2str_r64
  end interface num2str

  !> Generic evenly spaced vector generation interface.
  !!
  !! Supports integer and `REAL(kind=rk)` output arrays.
  interface linspace
     module procedure linspace_i
     module procedure linspace_rk
  end interface linspace

  !> Generic three-dimensional ND-grid construction interface.
  !!
  !! Supports integer and `REAL(kind=rk)` coordinate vectors.
  interface ndgrid
     module procedure ndgrid_i
     module procedure ndgrid_rk
  end interface ndgrid

  !> Generic mail helper interface.
  !!
  !! Supports text messages and two-dimensional `REAL(kind=rk)` arrays.
  interface send_mail
     module procedure send_mail_c
     module procedure send_mail_rk_2d
  end interface send_mail

  ! All routines below are module procedures.
contains

!!!=================================================================================================
  !> Convert Cartesian coordinates to polar coordinates.
  !!
  !! Computes the radial distance and the quadrant-aware polar angle of `(x1,x2)`.
  !!
  !! **Arguments**
  !!
  !! * `rho` - Radial distance from the origin.
  !! * `theta` - Polar angle returned by `atan2(x2,x1)`.
  !! * `x1` - First Cartesian coordinate.
  !! * `x2` - Second Cartesian coordinate.
  subroutine cart2pol(rho,theta,x1,x2)
    real(kind=rk),intent(out)              :: rho,theta
    real(kind=rk),intent(in)               :: x1,x2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Euclidean distance from the origin.
    rho   = sqrt(x1**2+x2**2)
    ! `atan2` preserves the correct quadrant of the polar angle.
    theta = atan2(x2,x1)

  end subroutine cart2pol
!!!=================================================================================================


!!!=================================================================================================
  !> Determine the global minimum or maximum of a distributed three-dimensional field.
  !!
  !! The local extremum and owning block image are combined across `block_comm`.
  !! The index of the selected extremum is then broadcast from the owning image so
  !! every participant receives the same value, image, and index tuple.
  !!
  !! **Arguments**
  !!
  !! * `idx` - Indices of the selected extremum within the owning local field.
  !! * `image` - Block image containing the selected global extremum.
  !! * `val` - Selected global minimum or maximum value.
  !! * `field` - Local three-dimensional field searched on each image.
  !! * `type` - Selection mode; supported values are `max` and `min`.
  subroutine get_index_rk(idx,image,val,field,type)
    integer         , dimension(:)    , intent(inout)            :: idx
    integer                           , intent(inout)            :: image
    real(kind=rk)                     , intent(inout)            :: val
    real(kind=rk)   , dimension(:,:,:), intent(in)               :: field
    character(len=*)                  , intent(in)               :: type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                                :: val_in(2),val_out(2)
    integer                                                      :: i,idx_bc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Normalize the selector string before choosing minimum or maximum handling.
    select case (trim(adjustl(type)))

    case ('max')
       ! value and image (rank)
       val_in(1) = maxval(field)
       val_in(2) = real(params%parallelism%block_image,kind=rk)

       ! Combine local maxima and owning images across the block communicator.
       call allreduce(val_out,val_in,'maxloc',params%parallelism%block_comm)

       val   = val_out(1)
       image = int(val_out(2))

       ! location (indices)
       idx = maxloc(field)
       do i = 1,size(idx)
          idx_bc = idx(i)
          ! Broadcast each index component from the image owning the global extremum.
          call bcast(idx_bc,image-1,params%parallelism%block_comm)
          idx(i) = idx_bc
       end do

    case ('min')
       ! value and image (rank)
       val_in(1) = minval(field)
       val_in(2) = real(params%parallelism%block_image,kind=rk)

       ! Combine local minima and owning images across the block communicator.
       call allreduce(val_out,val_in,'minloc',params%parallelism%block_comm)

       val   = val_out(1)
       image = int(val_out(2))

       ! location (indices)
       idx = minloc(field)
       do i = 1,size(idx)
          idx_bc = idx(i)
          ! Broadcast each index component from the image owning the global extremum.
          call bcast(idx_bc,image-1,params%parallelism%block_comm)
          idx(i) = idx_bc
       end do

    case default
       write(*,*) "error in helper.f90:get_index_rk:unknown type, stop"
       call stop_cruna

    end select

  end subroutine get_index_rk
!!!=================================================================================================


!!!=================================================================================================
  !> Generate an evenly spaced integer sequence between two endpoints.
  !!
  !! The number of values is determined by `size(x)`. For more than one output
  !! element, the integer step is computed as `(x1-x0)/(n-1)`.
  !!
  !! **Arguments**
  !!
  !! * `x` - Output vector containing the generated integer sequence.
  !! * `x0` - First value of the sequence.
  !! * `x1` - Last nominal endpoint used to determine the integer step.
  subroutine linspace_i(x,x0,x1)
    integer,dimension(:),intent(out) :: x
    integer,intent(in)               :: x0,x1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                          :: step
    integer                          :: i,n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Derive the requested number of output values from the caller-provided array.
    n = size(x)

    if (n.eq.1) then
       x    = x0
    else
       step = (x1-x0)/(n-1)
       x    = (/(((i-1)*step + x0),i=1,n)/)
    end if

  end subroutine linspace_i
!!!=================================================================================================


!!!=================================================================================================
  !> Generate an evenly spaced `REAL(kind=rk)` sequence between two endpoints.
  !!
  !! The number of values is determined by `size(x)`. For more than one output
  !! element, the spacing is `(x1-x0)/(n-1)`.
  !!
  !! **Arguments**
  !!
  !! * `x` - Output vector containing the generated real sequence.
  !! * `x0` - First value of the sequence.
  !! * `x1` - Last value of the sequence.
  subroutine linspace_rk(x,x0,x1)
    real(kind=rk),dimension(:),intent(out) :: x
    real(kind=rk),intent(in)               :: x0,x1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                          :: step
    integer                                :: i,n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Derive the requested number of output values from the caller-provided array.
    n = size(x)

    if (n.eq.1) then
       x    = x0
    else
       step = (x1-x0)/(n-1)
       x    = (/((real(i-1,rk)*step + x0),i=1,n)/)
    end if

  end subroutine linspace_rk
!!!=================================================================================================


!!!=================================================================================================
  !> Construct a three-dimensional meshgrid from coordinate vectors.
  !!
  !! The last dimension of `X` stores the coordinate component. The first two
  !! components are replicated across the third spatial direction; the third
  !! component is replicated across the first two spatial directions.
  !!
  !! **Arguments**
  !!
  !! * `X` - Output coordinate grid with component index in the fourth dimension.
  !! * `x1v` - Coordinate vector for the first direction.
  !! * `x2v` - Coordinate vector for the second direction.
  !! * `x3v` - Coordinate vector for the third direction.
  subroutine meshgrid(X,x1v,x2v,x3v)
    real(kind=rk), dimension(:,:,:,:), intent(out)  :: X
    real(kind=rk), dimension(:)      , intent(in)   :: x1v, x2v, x3v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer           :: sX1, sX2, sX3, i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Determine the coordinate-vector lengths used to replicate each component.
    sX1 = size(x1v)
    sX2 = size(x2v)
    sX3 = size(x3v)

    do i = 1,sX3
       X(:,:,i,1) = spread( x1v, 1, sX2 )
       X(:,:,i,2) = spread( x2v, 2, sX1 )
    enddo

    do i = 1,sX1
       X(i,:,:,3) = spread( x3v, 1, sX2)
    enddo

  end subroutine meshgrid
!!!=================================================================================================


!!!=================================================================================================
  !> Construct an integer three-dimensional ND grid.
  !!
  !! Each coordinate vector is expanded along the other two directions and stored
  !! in the corresponding component of the fourth dimension of `X`.
  !!
  !! **Arguments**
  !!
  !! * `X` - Output integer coordinate grid.
  !! * `x1v` - Coordinate vector for the first direction.
  !! * `x2v` - Coordinate vector for the second direction.
  !! * `x3v` - Coordinate vector for the third direction.
  subroutine ndgrid_i(X,x1v,x2v,x3v)
    integer, dimension(:,:,:,:), intent(out)  :: X
    integer, dimension(:)      , intent(in)   :: x1v, x2v, x3v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                   :: sX1, sX2, sX3, i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Determine the coordinate-vector lengths used to replicate each component.
    sX1 = size(x1v)
    sX2 = size(x2v)
    sX3 = size(x3v)

    ! write(*,*) sX1,sX2,sX3,params%parallelism%block_image

    do i = 1,sX1
       X(i,:,:,1) = x1v(i)
    enddo

    do i = 1,sX2
       X(:,i,:,2) = x2v(i)
    enddo

    do i = 1,sX3
       X(:,:,i,3) = x3v(i)
    enddo

  end subroutine ndgrid_i
!!!=================================================================================================


!!!=================================================================================================
  !> Construct a `REAL(kind=rk)` three-dimensional ND grid.
  !!
  !! Each coordinate vector is expanded along the other two directions and stored
  !! in the corresponding component of the fourth dimension of `X`.
  !!
  !! **Arguments**
  !!
  !! * `X` - Output real coordinate grid.
  !! * `x1v` - Coordinate vector for the first direction.
  !! * `x2v` - Coordinate vector for the second direction.
  !! * `x3v` - Coordinate vector for the third direction.
  subroutine ndgrid_rk(X,x1v,x2v,x3v)
    real(kind=rk), dimension(:,:,:,:), intent(out)  :: X
    real(kind=rk), dimension(:)      , intent(in)   :: x1v, x2v, x3v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                         :: sX1, sX2, sX3, i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Determine the coordinate-vector lengths used to replicate each component.
    sX1 = size(x1v)
    sX2 = size(x2v)
    sX3 = size(x3v)

    ! write(*,*) sX1,sX2,sX3,params%parallelism%block_image

    do i = 1,sX1
       X(i,:,:,1) = x1v(i)
    enddo

    do i = 1,sX2
       X(:,i,:,2) = x2v(i)
    enddo

    do i = 1,sX3
       X(:,:,i,3) = x3v(i)
    enddo

  end subroutine ndgrid_rk
!!!=================================================================================================


!!!=================================================================================================
  !> Convert an integer value to a left-adjusted character string.
  !!
  !! When `form` is present it is used as the internal `WRITE` format; otherwise
  !! list-directed formatting is used. Conversion failure returns `***`.
  !!
  !! **Arguments**
  !!
  !! * `num` - Integer value to convert.
  !! * `form` - Optional Fortran format specification.
  !!
  !! **Returns**
  !!
  !! Left-adjusted character representation of `num`.
  character(len=num2str_len) function num2str_i(num,form)
    integer                   , intent(in)     :: num
    character(len=*), optional, intent(in)     :: form
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                    :: io_stat
    character(len=num2str_len)                 :: str
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(present(form))then
       write(str,form,iostat=io_stat) num
    else
       write(str,*   ,iostat=io_stat) num
    end if

    if (io_stat.eq.0) then
       num2str_i = adjustl(str)
    else
       num2str_i = '***'
    end if


  end function num2str_i
!!!=================================================================================================


!!!=================================================================================================
  !> Convert a default `REAL(kind=real32)` value to a left-adjusted character string.
  !!
  !! When `form` is present it is used as the internal `WRITE` format; otherwise
  !! list-directed formatting is used. Conversion failure returns `***`.
  !!
  !! **Arguments**
  !!
  !! * `num` - Real(kind=real32) value to convert.
  !! * `form` - Optional Fortran format specification.
  !!
  !! **Returns**
  !!
  !! Left-adjusted character representation of `num`.
  character(len=num2str_len) function num2str_r32(num,form)
    real(kind=real32)         , intent(in)     :: num
    character(len=*), optional, intent(in)     :: form
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                    :: io_stat
    character(len=num2str_len)                 :: str
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(present(form))then
       write(str,form,iostat=io_stat) num
    else
       write(str,*   ,iostat=io_stat) num
    end if

    if (io_stat.eq.0) then
       num2str_r32 = adjustl(str)
    else
       num2str_r32 = '***'
    end if

  end function num2str_r32
!!!=================================================================================================


!!!=================================================================================================
  !> Convert a `REAL(kind=real64)` value to a left-adjusted character string.
  !!
  !! When `form` is present it is used as the internal `WRITE` format; otherwise
  !! list-directed formatting is used. Conversion failure returns `***`.
  !!
  !! **Arguments**
  !!
  !! * `num` - `REAL(kind=real64)` value to convert.
  !! * `form` - Optional Fortran format specification.
  !!
  !! **Returns**
  !!
  !! Left-adjusted character representation of `num`.
  character(len=num2str_len) function num2str_r64(num,form)
    real(kind=real64)          , intent(in)     :: num
    character(len=*) , optional, intent(in)     :: form
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                     :: io_stat
    character(len=num2str_len)                  :: str
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(present(form))then
       write(str,form,iostat=io_stat) num
    else
       write(str,*   ,iostat=io_stat) num
    end if

    if (io_stat.eq.0) then
       num2str_r64 = adjustl(str)
    else
       num2str_r64 = '***'
    end if

  end function num2str_r64
!!!=================================================================================================


!!!=================================================================================================
  !> Start a wall-clock timing interval using the Fortran system clock.
  !!
  !! The current counter value is stored internally and may optionally be returned
  !! so multiple independent timing intervals can be maintained by the caller.
  !!
  !! **Arguments**
  !!
  !! * `tic_out` - Optional returned system-clock counter value.
  subroutine tic(tic_out)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, optional, intent(out) :: tic_out
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Capture the current system-clock counter.
    call system_clock(tic_persistent)

    if(present(tic_out)) then
       tic_out = tic_persistent
    end if

  end subroutine tic
!!!=================================================================================================


!!!=================================================================================================
  !> Return elapsed wall-clock time in seconds.
  !!
  !! If `tic_in` is present it is used as the starting counter; otherwise the most
  !! recent counter stored by `tic` is used.
  !!
  !! **Arguments**
  !!
  !! * `tic_in` - Optional system-clock counter marking the start of the interval.
  !!
  !! **Returns**
  !!
  !! Elapsed wall-clock time in seconds.
  real function toc(tic_in)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, optional, intent(in) :: tic_in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                       :: tic_val, toc_val, tic_rate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(present(tic_in)) then
       tic_val = tic_in
    else
       tic_val = tic_persistent
    end if

    ! Read the current counter and the number of clock ticks per second.
    call system_clock(toc_val, tic_rate)

    ! Convert elapsed clock ticks to seconds.
    toc = real(toc_val - tic_val)/real(tic_rate)

  end function toc
!!!=================================================================================================


!!!=================================================================================================
  !> Prepare a text message for delivery through the system mail command.
  !!
  !! Only `world_image == 1` performs the mail handling. The recipient is read from
  !! `io.mail_recv`; if no recipient is configured, an error message is printed.
  !! The actual `execute_command_line` call is currently commented out in the code.
  !!
  !! **Arguments**
  !!
  !! * `mail_subj` - Subject text appended to the CRUNA host identifier.
  !! * `mail_text` - Message body.
  subroutine send_mail_c(mail_subj,mail_text)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*), intent(in)           :: mail_subj,mail_text
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=  4096                )  :: mail_text_out
    character(len=  max_length_parameter)  :: mail_recv
    character(len=3*max_length_parameter)  :: system_cmd
    ! integer                                :: estat, cstat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=128)                     :: host_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%world_image.eq.1) then                                                  ! just image 1 send mails

       ! Read the configured recipient; the sentinel `no_recv` disables mail preparation.
       call get_parameter(mail_recv,'io.mail_recv',default='no_recv')

       if (mail_recv.ne.'no_recv') then                                                            ! check for recv

          ! str2str
          mail_text_out = mail_text

          ! host name
          call hostnm( host_name )

          ! system command
          system_cmd = "echo '" // trim(mail_text) //  "' | mail -s '" // "cruna@" // trim(host_name) // ":" // trim(mail_subj) // "' " // trim(mail_recv) // " >> /dev/null"
          !call execute_command_line(system_cmd,.false.,estat,cstat)

       else
          write(*,*) "error in helper.f90:send_mail:no 'io.mail_recv' given"
       end if
    end if

  end subroutine send_mail_c
!!!=================================================================================================


!!!=================================================================================================
  !> Prepare a two-dimensional `REAL(kind=rk)` array for delivery by system mail.
  !!
  !! The array is first converted to text using an internal `WRITE`. Only
  !! `world_image == 1` performs the mail handling. The actual command execution
  !! remains disabled by the commented `execute_command_line` call.
  !!
  !! **Arguments**
  !!
  !! * `mail_subj` - Subject text appended to the CRUNA host identifier.
  !! * `mail_text` - Two-dimensional real array converted to the message body.
  subroutine send_mail_rk_2d(mail_subj,mail_text)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*), intent(in)                     :: mail_subj
    real(kind=rk)   , intent(in), dimension(:,:)     :: mail_text
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=  4096                )            :: mail_text_out
    character(len=  max_length_parameter)            :: mail_recv
    character(len=3*max_length_parameter)            :: system_cmd
    ! integer                                          :: estat, cstat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=128)                               :: host_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%world_image.eq.1) then                                                  ! just image 1 send mails

       ! Read the configured recipient; the sentinel `no_recv` disables mail preparation.
       call get_parameter(mail_recv,'io.mail_recv',default='no_recv')

       if (mail_recv.ne.'no_recv') then                                                            ! check for recv

          ! real2str
          write(mail_text_out,*) mail_text

          ! host name
          call hostnm( host_name )

          ! system command
          system_cmd = "echo '" // trim(mail_text_out) //  "' | mail -s '" // "cruna@" // trim(host_name) // ":" // trim(mail_subj) // "' " // trim(mail_recv) !// " >> /dev/null"
          !call execute_command_line(system_cmd,.false.,estat,cstat)

       else
          write(*,*) "error in helper.f90:send_mail:no 'io.mail_recv' given"
       end if
    end if

  end subroutine send_mail_rk_2d
!!!=================================================================================================

!!!=================================================================================================
  !> Execute a shell command asynchronously.
  !!
  !! An optional integer selects a predefined command prefix. Invalid prefix indices
  !! return without executing anything. Without `pre_cmd`, `cmd` is executed directly.
  !!
  !! **Arguments**
  !!
  !! * `cmd` - Command text or argument appended to the selected prefix.
  !! * `pre_cmd` - Optional index into the internal list of allowed command prefixes.
  subroutine command(cmd,pre_cmd)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*), intent(in)           :: cmd
    integer,          intent(in), optional :: pre_cmd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=:), allocatable :: full_cmd
    character(len=32), parameter  :: command_list(*) = [ character(len=32) ::  &
         "echo",               & ! 1 - echo
         "vmtouch -t",         & ! 2 - pre-read file using vmtouch
         "cat" ]                 ! 3 - pre-read file using cat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (present(pre_cmd)) then

       ! Reject prefix indices outside the predefined command list.
       if (pre_cmd < 1 .or. pre_cmd > size(command_list)) return

       full_cmd = trim(command_list(pre_cmd)) // " " // trim(cmd)

    else

       full_cmd = trim(cmd)

    end if

    ! Launch the command without waiting for completion.
    call execute_command_line(full_cmd,wait=.false.)

  end subroutine command
!!!=================================================================================================

!!!=================================================================================================
  !> Write a simple single-line text progress bar.
  !!
  !! The supplied fraction controls both the number of `=` characters and the
  !! displayed percentage. Output uses a carriage return and no line advance so
  !! repeated calls update the same terminal line.
  !!
  !! **Arguments**
  !!
  !! * `p` - Progress fraction used to construct the bar and percentage.
  subroutine progressbar(p)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real, intent(in)  :: p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=40) :: progbar
    character(len=25) :: prog
    integer           :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Start with an empty visible progress segment.
    prog = ""

    ! Fill the visible segment proportionally to the supplied progress fraction.
    do i = 1,ceiling(p*real(len(prog)))
       prog = trim(prog) // "="
    end do

    progbar = adjustl("[" // prog // "] " // trim(num2str(ceiling(p*100),'(I3.3)')) // "%" )
    ! Rewrite the current terminal line without advancing to a new line.
    write(*,'(a1,a1,a40)',advance='no') char(13)," ", progbar

  end subroutine progressbar
!!!=================================================================================================

end module helper
