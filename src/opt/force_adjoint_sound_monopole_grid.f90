module force_adjoint_sound_monopole_grid

  use data_geom
  use io
  USE mpi
  use parameter
  use tub_ak_data

  real(kind=rk), dimension(:,:,:,:), allocatable, save, public :: f,f_backup,f_grad
  real(kind=rk), dimension(:,:,:)  , allocatable, save, public :: theta_x
  real(kind=rk), dimension(:,:)    , allocatable, save, public :: theta_t

  private

  public allocate_force
  public apply_force
  public backup_force
  public calc_grad
  public init_force
  public store_force
  public restore_force
  public update_force
  public init_theta_t
  public init_theta_x

contains

!!!=================================================================================================
  subroutine allocate_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call tub_ak_read_params

    allocate(     f(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%time%steps))
    allocate(f_grad(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%time%steps))

  end subroutine allocate_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine apply_force(rhs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)                    :: rhs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    rhs(:,:,:,5) = rhs(:,:,:,5) + f(:,:,:,params%time%nt)

  end subroutine apply_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine backup_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                               :: s,s_backup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (.not.allocated(f_backup)) then
       allocate(f_backup(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%time%steps))
    end if

    if (params%time%subsets.eq.1) then
       f_backup = f
    else
       s_backup = params%time%ns

       do s = 1,params%time%subsets
          if(params%parallelism%block_image.eq.1) then
             if (params%io%verbosity.ge.1) then
                write(*,*) " backup force of subset", s
             end if
          end if
          
          call set_parameter(s,'time.ns',struct_val=params%time%ns)
          
          call load(f     ,'force_subsets_cache')
          call store(f    ,'force_backup_subsets_cache')
       end do

       call set_parameter(s_backup,'time.ns',struct_val=params%time%ns)
    end if

  end subroutine backup_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine calc_grad(Qs,Q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:),intent(inout)       :: Qs     ! (inout) is required for multiple subsets
    real(kind=rk), dimension(:,:,:,:,:),intent(in),optional :: Q      ! if regularization is applied
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                 :: s,s_backup,loop,gfreq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(params%parallelism%block_image.eq.1) then
       ! all sub block own information of all speaker, see above
       write(*,*) "write gradient"
    end if

    if (params%time%subsets.eq.1) then
       f_grad = -Qs(:,:,:,5,:)
       call get_parameter(loop,'opt.loop')
       call get_parameter(gfreq,'io.gfreq',default=huge(1))
       if (mod(loop,gfreq).eq.0) then
          call store(f_grad,'gradient_loop')
       end if
    else

       s_backup = params%time%ns

       do s = 1,params%time%subsets
          if(params%parallelism%block_image.eq.1) then
             if (params%io%verbosity.ge.1) then
                write(*,*) " store force gradient of subset", s
             end if
          end if

          call set_parameter(s,'time.ns',struct_val=params%time%ns)

          call load(Qs,'data_adjoint_subsets_cache')
          f_grad = -Qs(:,:,:,5,:)

          call store(f_grad,'force_grad_subsets_cache')
          call store(f_grad,'gradient_loop') 
       end do

       call set_parameter(s_backup,'time.ns',struct_val=params%time%ns)
    end if

  end subroutine calc_grad
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                               :: s,s_backup
    character(len=max_length_parameter)                   :: fname    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! init f with all zero
    f              = 0.0_rk
    f_grad         = 0.0_rk

    ! check for present pre-defined force
    call get_parameter(fname,'opt.init_force_file',default = '-')

    ! handle multiple subsets
    if (params%time%subsets.eq.1) then

       if (trim(fname).ne.'-') then
          if (params%parallelism%block_image.eq.1) then
             write(*,*) "warning in force_adjoint_sound_monopole_array.f90:init_force:load present force from file:", fname
          end if

          ! load single pre-defined force from file, e.g. from optimization run
          call load(f, 'dummy', fname_optin = trim(fname))
       end if

    else

       if (trim(fname).ne.'-') then
          if (params%parallelism%block_image.eq.1) then
             write(*,*) "warning in force_adjoint_sound_monopole_array.f90:init_force:load present force from file not implemented for multiple subsets:stop"
             stop
          end if
       end if

       s_backup = params%time%ns

       if(params%parallelism%block_image.eq.1) then
          write(*,*) "write initial force and force_grad for all subsets"
       end if

       do s = 1,params%time%subsets
          if(params%parallelism%block_image.eq.1) then
             if (params%io%verbosity.ge.1) then
                write(*,*) " store force and force_grad of subset", s
             end if
          end if
          
          call set_parameter(s,'time.ns',struct_val=params%time%ns)

          call store(f     ,'force_subsets_cache')
          call store(f_grad,'force_grad_subsets_cache')
       end do
       call set_parameter(s_backup,'time.ns',struct_val=params%time%ns)

    end if

  end subroutine init_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine store_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                               :: s,s_backup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(params%parallelism%block_image.eq.1) then
       ! all sub block own information of all speaker, see above
       write(*,*) "store force"
    end if

    if (params%time%subsets.eq.1) then
       call store(f,'force_loop')
    else
       s_backup = params%time%ns

       do s = 1,params%time%subsets
          call set_parameter(s,'time.ns',struct_val=params%time%ns)

          if(params%parallelism%block_image.eq.1) then
             if (params%io%verbosity.ge.1) then
                write(*,*) " store force of subset", s
             end if
          end if

          call load(f     ,'force_subsets_cache')
          call store(f    ,'force_loop')
       end do

       call set_parameter(s_backup,'time.ns',struct_val=params%time%ns)
    end if

  end subroutine store_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine restore_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                               :: s,s_backup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%time%subsets.eq.1) then
       f = f_backup
    else
       s_backup = params%time%ns

       do s = 1,params%time%subsets
          
          if(params%parallelism%block_image.eq.1) then
             if (params%io%verbosity.ge.1) then
                write(*,*) " restore force of subset", s
             end if
          end if

          call set_parameter(s,'time.ns',struct_val=params%time%ns)

          call load(f     ,'force_backup_subsets_cache')
          call store(f    ,'force_subsets_cache')
       end do

       call set_parameter(s_backup,'time.ns',struct_val=params%time%ns)
    end if

  end subroutine restore_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine update_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                         :: alpha,alpha0
    integer                                               :: nt,s,s_backup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! allocate/init theta_t
    if (.not.allocated(theta_t)) then
       call init_theta_t
    end if

    ! allocate/init theta_x
    if (.not.allocated(theta_x)) then
       call init_theta_x
    end if

    ! get opt value
    call get_parameter(alpha0,'opt.alpha0'                 )
    call get_parameter(alpha ,'opt.alpha' ,default = alpha0)

    if (params%time%subsets.eq.1) then
       do nt = 1,params%time%steps
          f(:,:,:,nt) = f(:,:,:,nt) + alpha*f_grad(:,:,:,nt)*theta_x*theta_t(nt,params%time%ns)
       end do
    else
       s_backup = params%time%ns
       do s=1,params%time%subsets
          if(params%parallelism%block_image.eq.1) then
             if (params%io%verbosity.ge.1) then
                write(*,*) " update force of subset", s
             end if
          end if

          call set_parameter(s,'time.ns',struct_val=params%time%ns)

          call load(f     ,'force_subsets_cache')
          call load(f_grad,'force_grad_subsets_cache')

          do nt = 1,params%time%steps
             f(:,:,:,nt) = f(:,:,:,nt) + alpha*f_grad(:,:,:,nt)*theta_x*theta_t(nt,params%time%ns)
          end do

          call store(f,'force_subsets_cache')
       end do
       call set_parameter(s_backup,'time.ns',struct_val=params%time%ns)
    end if

  end subroutine update_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_theta_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:), allocatable                           :: theta_t_vector
    integer                                                      :: discarded_steps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! allocates
    allocate(theta_t(params%time%steps, params%time%subsets))
    allocate(theta_t_vector(params%time%steps * params%time%subsets))

    theta_t        = 1.0_rk
    theta_t_vector = 1.0_rk

    ! discard time steps
    call get_parameter(discarded_steps,'opt.sigma_t.discarded_steps',default = 0)

    if (discarded_steps.gt.0) then
       !theta_t_vector(1:discarded_steps) = 0.0_rk ! see sketch below
       theta_t_vector((params%time%steps*params%time%subsets) - discarded_steps + 1 : (params%time%steps*params%time%subsets)) = 0.0_rk
    end if

    ! reshape to array
    theta_t = reshape(theta_t_vector,(/params%time%steps,params%time%subsets/))

    call store(theta_t,'theta_t_function')

    !    characteristic acoustic path
    !    signal ////////////////////|
    !    is    //////////////////// |
    !    zero ////////////////////  |
    !    here//// slope is c ////   |
    !   ^   ////////////////////    |
    ! x |  ////////////////////     |
    !   | ////////////////////      |
    !   |111111111111111111110000000|  <- forcing mask (avoids supersonic effects)
    !   ------>              |------|     implies that signal has trailing zeros
    !       t                extend/c


  end subroutine init_theta_t
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_theta_x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), allocatable                 :: src_pos_distance
    integer                                                      :: src_radius_dx1fac, src_flank_width_dx1fac
    real(kind=rk)                                                :: src_radius, src_flank_width
    real(kind=rk)                                                :: pi
    integer                                                      :: s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! constants
    pi = 4.0_rk*atan(1.0_rk)

    ! parameter
    call get_parameter(src_radius_dx1fac     ,'opt.src_radius_dx1fac'     ,default = 10)
    call get_parameter(src_flank_width_dx1fac,'opt.src_flank_width_dx1fac',default =  5)

    src_radius      = src_radius_dx1fac     *params%geom%dx1
    src_flank_width = src_flank_width_dx1fac*params%geom%dx1
    ! flank width is defined as x-length between 0 and 1 of theta_x in linear approximation:
    ! df/dx = slope = Dx[0.5*erf(sqrt(pi)*x/len)]= 1/len -> since df=1 it follows dx = len 

    ! allocates
    allocate(theta_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    theta_x = 0.0_rk

    allocate(src_pos_distance(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    src_pos_distance = 0.0_rk

    do s = 1,speaker_nbr
       ! compute distance
       src_pos_distance = sqrt((X(:,:,:,1) - speaker_pos(s,1))**2 + (X(:,:,:,2) - speaker_pos(s,2))**2 + (X(:,:,:,3) - speaker_pos(s,3))**2)

       ! compute theta_x
       theta_x = theta_x + (0.5_rk - 0.5_rk*erf(sqrt(pi)/src_flank_width*(src_pos_distance - src_radius)))
    end do

    call store(theta_x,'theta_x_function')

  end subroutine init_theta_x
!!!=================================================================================================

end module force_adjoint_sound_monopole_grid





