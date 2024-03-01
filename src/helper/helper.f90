module helper

  use parallelism
  use parameter

  private

  public cart2pol
  public get_index
  public linspace
  public meshgrid
  public ndgrid
  public num2str
  public progressbar
  public tic
  public toc
  public send_mail

  integer, save, private :: tic_persistent
  integer, parameter     :: num2str_len = 32

  interface get_index
     module procedure get_index_rk
  end interface get_index

  interface num2str
     module procedure num2str_i
     module procedure num2str_r
     module procedure num2str_rk
  end interface num2str

  interface linspace 
     module procedure linspace_i
     module procedure linspace_rk
  end interface linspace

  interface ndgrid
     module procedure ndgrid_i
     module procedure ndgrid_rk
  end interface ndgrid

  interface send_mail
     module procedure send_mail_c
     module procedure send_mail_rk_2d
  end interface send_mail

contains

!!!=================================================================================================
  subroutine cart2pol(rho,theta,x1,x2)
    real(kind=rk),intent(out)              :: rho,theta
    real(kind=rk),intent(in)               :: x1,x2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   rho   = sqrt(x1**2+x2**2)
   theta = atan2(x2,x1)

  end subroutine cart2pol
!!!=================================================================================================


!!!=================================================================================================
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

    select case (trim(adjustl(type)))

    case ('max')
       ! value and image (rank)
       val_in(1) = maxval(field)
       val_in(2) = real(params%parallelism%block_image,kind=rk)

       call allreduce(val_out,val_in,'maxloc',params%parallelism%block_comm)

       val   = val_out(1)
       image = int(val_out(2))

       ! location (indices)
       idx = maxloc(field)
       do i = 1,size(idx)
          idx_bc = idx(i)
          call bcast(idx_bc,image-1,params%parallelism%block_comm)
          idx(i) = idx_bc          
       end do

    case ('min')
       ! value and image (rank)
       val_in(1) = minval(field)
       val_in(2) = real(params%parallelism%block_image,kind=rk)

       call allreduce(val_out,val_in,'minloc',params%parallelism%block_comm)

       val   = val_out(1)
       image = int(val_out(2))

       ! location (indices)
       idx = minloc(field)
       do i = 1,size(idx)
          idx_bc = idx(i)
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
  subroutine linspace_i(x,x0,x1)
    integer,dimension(:),intent(out) :: x
    integer,intent(in)               :: x0,x1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                          :: step
    integer                          :: i,n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  subroutine linspace_rk(x,x0,x1)
    real(kind=rk),dimension(:),intent(out) :: x
    real(kind=rk),intent(in)               :: x0,x1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                          :: step
    integer                                :: i,n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  subroutine meshgrid(X,x1v,x2v,x3v)
    real(kind=rk), dimension(:,:,:,:), intent(out)  :: X
    real(kind=rk), dimension(:)      , intent(in)   :: x1v, x2v, x3v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer           :: sX1, sX2, sX3, i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  subroutine ndgrid_i(X,x1v,x2v,x3v)
    integer, dimension(:,:,:,:), intent(out)  :: X
    integer, dimension(:)      , intent(in)   :: x1v, x2v, x3v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                   :: sX1, sX2, sX3, i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  subroutine ndgrid_rk(X,x1v,x2v,x3v)
    real(kind=rk), dimension(:,:,:,:), intent(out)  :: X
    real(kind=rk), dimension(:)      , intent(in)   :: x1v, x2v, x3v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                         :: sX1, sX2, sX3, i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
  character(len=num2str_len) function num2str_r(num,form)
    real            , intent(in)               :: num
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
       num2str_r = adjustl(str)
    else
       num2str_r = '***'
    end if

  end function num2str_r
!!!=================================================================================================


!!!=================================================================================================
  character(len=num2str_len) function num2str_rk(num,form)
    real(kind=rk)    , intent(in)               :: num
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
       num2str_rk = adjustl(str)
    else
       num2str_rk = '***'
    end if

  end function num2str_rk
!!!=================================================================================================


!!!=================================================================================================
  subroutine tic(tic_out)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    integer, optional, intent(out) :: tic_out
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    call system_clock(tic_persistent)

    if(present(tic_out)) then
       tic_out = tic_persistent
    end if

  end subroutine tic
!!!=================================================================================================


!!!=================================================================================================
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

    call system_clock(toc_val, tic_rate)

    toc = real(toc_val - tic_val)/real(tic_rate)

  end function toc
!!!=================================================================================================


!!!=================================================================================================
  subroutine send_mail_c(mail_subj,mail_text)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*), intent(in)           :: mail_subj,mail_text
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=  4096                )  :: mail_text_out
    character(len=  max_length_parameter)  :: mail_recv
    character(len=3*max_length_parameter)  :: system_cmd
    integer                                :: estat, cstat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=128)                     :: host_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%world_image.eq.1) then                                                  ! just image 1 send mails

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
          write(*,*) "error in helper.f90:send_mail:no 'params.io.mail_recv' given"
       end if
    end if

  end subroutine send_mail_c
!!!=================================================================================================


!!!=================================================================================================
  subroutine send_mail_rk_2d(mail_subj,mail_text)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*), intent(in)                     :: mail_subj
    real(kind=rk)   , intent(in), dimension(:,:)     :: mail_text
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=  4096                )            :: mail_text_out
    character(len=  max_length_parameter)            :: mail_recv
    character(len=3*max_length_parameter)            :: system_cmd
    integer                                          :: estat, cstat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=128)                               :: host_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%parallelism%world_image.eq.1) then                                                  ! just image 1 send mails

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
          write(*,*) "error in helper.f90:send_mail:no 'params.io.mail_recv' given"
       end if
    end if

  end subroutine send_mail_rk_2d
!!!=================================================================================================

!!!=================================================================================================
  subroutine progressbar(p)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real, intent(in)  :: p
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=40) :: progbar
    character(len=25) :: prog
    integer           :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    prog = ""

    do i = 1,ceiling(p*real(len(prog)))
       prog = trim(prog) // "="
    end do

    progbar = adjustl("[" // prog // "] " // trim(num2str(ceiling(p*100),'(I3.3)')) // "%" )
    write(*,'(a1,a1,a40)',advance='no') char(13)," ", progbar

  end subroutine progressbar
!!!=================================================================================================

end module helper
