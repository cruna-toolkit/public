module force_adjoint_monopole_synthesis

  use data_geom
  use io
  USE mpi
  use parameter
  use tub_ak_data

  real(kind=rk), dimension(:,:,:,:), allocatable, save, public :: f
  real(kind=rk), dimension(:,:)    , allocatable, save, public :: f_speaker, f_speaker_backup
  real(kind=rk), dimension(:,:)    , allocatable, save, public :: f_speaker_grad
  real(kind=rk), dimension(:)      , allocatable, save, public :: theta_t

  private

  public allocate_force
  public apply_force
  public backup_force
  public calc_grad
  public init_force
  public store_force
  public restore_force
  public update_force
  public update_speaker
  public init_theta_t

contains

!!!=================================================================================================
  subroutine allocate_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call tub_ak_read_params

    allocate(f(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%time%steps))

    allocate(     f_speaker(params%time%steps,speaker_nbr))
    allocate(f_speaker_grad(params%time%steps,speaker_nbr))

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
    if (.not.allocated(f_speaker_backup)) then
       allocate(f_speaker_backup(params%time%steps,speaker_nbr))
    end if

    f_speaker_backup = f_speaker

  end subroutine backup_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine calc_grad(Qs,Q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:),intent(in)          :: Qs
    real(kind=rk), dimension(:,:,:,:,:),intent(in),optional :: Q      ! if regularization is applied
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:), allocatable                :: f_speaker_grad_transfer
    integer                                                 :: s, ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(f_speaker_grad_transfer(size(f_speaker_grad,1)))

    do s = 1,speaker_nbr

       ! only the block in which the speaker is located evaluates the adjoint/gradient
       if (speaker_idx(s,7).ne.0) then
          f_speaker_grad(:,s) = -Qs(speaker_idx(s,4),speaker_idx(s,5),speaker_idx(s,6),5,:)
       else
          f_speaker_grad(:,s) = 0.0_rk
       end if

       ! mpi_reduce to update f_speaker_grad on all blocks
       call mpi_allreduce(f_speaker_grad(:,s),f_speaker_grad_transfer,size(f_speaker_grad(:,s)), mpi_double, mpi_sum, params%parallelism%block_comm, ierr)
       f_speaker_grad(:,s) = f_speaker_grad_transfer

       ! now each block should have the temporal gradient of all speakers
       ! even if they are located in other blocks

    end do

    if(params%parallelism%block_image.eq.1) then
       ! all sub block own information of all speaker, see above
       write(*,*) "write gradient"
       call store(f_speaker_grad,'gradient_loop')
    end if

  end subroutine calc_grad
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                               :: i,use_only_speaker_nbr
    character(len=max_length_parameter)                   :: fname
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! init f with all zero
    f              = 0.0_rk
    f_speaker      = 0.0_rk
    f_speaker_grad = 0.0_rk

    ! load present force from file
    call get_parameter(fname,'opt.init_force_file',default = '-')

    if (trim(fname).ne.'-') then                                        ! parameter is not the default value, there is something to do

       if (params%parallelism%block_image.eq.1) then
          write(*,*) "warning in force_adjoint_sound_monopole_array.f90:init_force:load present force from file:", fname
       end if

       call load(f_speaker, 'dummy', fname_optin = trim(fname))         ! load present force, e.g. from optimization run

       call get_parameter(use_only_speaker_nbr,'opt.use_only_speaker_nbr',default = 0)
       if(use_only_speaker_nbr.ne.0) then                               ! use only one out of x speaker signals
          do i = 1,speaker_nbr                                          ! in order to test phase delay between speakers
             if (i.ne.use_only_speaker_nbr) then                        ! not beautiful but it works
                f_speaker(:,i) = 0.0_rk
             end if
          end do
       end if

       call update_force                                                ! use update_force routine to compute and set f
       !                                                                ! ATTENTION: f_speaker_grad must be zero, see lines above
    end if

  end subroutine init_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine store_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(params%parallelism%block_image.eq.1) then
       ! all sub block own information of all speaker, see above
       write(*,*) "write force"
       call store(f_speaker,'force_loop')
    end if
  end subroutine store_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine restore_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    f_speaker = f_speaker_backup

  end subroutine restore_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine update_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), allocatable          :: fx
    real(kind=rk)                                         :: alpha,alpha0
    real(kind=rk)                                         :: pi,sigma_fac,sigma
    real(kind=rk)                                         :: fx_alpha,fx_norm
    integer                                               :: s,nt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! allocate/init theta_t
    if (.not.allocated(theta_t)) then
       call init_theta_t
    end if

    ! allocate fx/def of sigma/pi
    allocate(fx(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    pi       = 4.0_rk*atan(1.0_rk)

    call get_parameter(sigma_fac,'opt.force_gauss_var_fac',default = 2.0_rk)
    sigma    = sigma_fac*params%geom%dx1
    fx_alpha = 0.5_rk/(sigma**2)
    fx_norm  = (2.0_rk*pi)**(-1.5_rk) * sigma**(-3_rk)

    ! reset current f
    f = 0.0_rk

    ! get opt value
    call get_parameter(alpha0,'opt.alpha0'                 )
    call get_parameter(alpha ,'opt.alpha' ,default = alpha0)

    do s = 1,speaker_nbr

       ! update speaker signal
       do nt = 1,params%time%steps
          f_speaker(nt,s) = f_speaker(nt,s) + alpha*f_speaker_grad(nt,s)*theta_t(nt)
       end do

       fx = fx_norm*exp(-fx_alpha*((X(:,:,:,1) - speaker_pos(s,1))**2 + (X(:,:,:,2) - speaker_pos(s,2))**2 + (X(:,:,:,3) - speaker_pos(s,3))**2))

       do nt = 1,params%time%steps
          f(:,:,:,nt) = f(:,:,:,nt) + fx*f_speaker(nt,s)
       end do
    end do

    !call store(fx,'debug_field')

  end subroutine update_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_theta_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                      :: discarded_steps
    integer                                                      :: nt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(theta_t(params%time%steps))
    theta_t = 1.0_rk

    call get_parameter(discarded_steps,'opt.sigma_t.discarded_steps',default = 0)

    do nt = 1,discarded_steps
       theta_t(nt) = 0.0_rk
    end do

    do  nt = params%time%steps-discarded_steps+1,params%time%steps
       theta_t(nt) = 0.0_rk
    end do

    call store(theta_t,'theta_t_function')

  end subroutine init_theta_t
!!!=================================================================================================

!!!=================================================================================================
  subroutine calc_time_av_sensitivity(time_av_sensitivity,Qs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:)    , intent(out)            :: time_av_sensitivity
    real(kind=rk), dimension(:,:,:,:,:), intent(in)             :: Qs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    time_av_sensitivity = sum(abs(Qs(:,:,:,5,:)),4)

  end subroutine calc_time_av_sensitivity
!!!=================================================================================================

!!!=================================================================================================
  subroutine update_speaker(Qs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:), intent(in)             :: Qs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:)    , allocatable            :: time_av_sensitivity
    real(kind=rk)                                               :: time_av_sensitivity_max
    real(kind=8) , dimension(2)                                 :: max_on_rank, max_all ! arrays 4 allreduce
    integer      , dimension(3)                                 :: time_av_sensitivity_max_idxs
    integer                                                     :: ierr,s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_parameter(s,'opt.active_monopoles',default=0) ! should be zero in the first run
    
    if (s.lt.speaker_nbr) then
       
       ! compute time_av_sensitivities
       allocate(time_av_sensitivity(params%geom%n1b,params%geom%n2b,params%geom%n3b))
       time_av_sensitivity = 0.0_rk
       call calc_time_av_sensitivity(time_av_sensitivity,Qs)

       ! get max sensitivity (on each rank)
       time_av_sensitivity_max = maxval(time_av_sensitivity)

       ! get max sensitivity and rank holding it
       max_on_rank(1) = time_av_sensitivity_max
       max_on_rank(2) = params%parallelism%block_image

       call mpi_allreduce(max_on_rank,max_all,1,mpi_2double_precision,mpi_maxloc,params%parallelism%block_comm,ierr)
       ! now max_all(1) should contain the overall max sensitivity
       !     max_all(2) should contain the rank holding the overall max sensitivity

       ! increase monopole counter
       s = s + 1
       call set_parameter(s,'opt.active_monopoles')

       ! update speaker positions
       time_av_sensitivity_max_idxs = -42
       if (int(max_all(2)).eq.params%parallelism%block_image) then
          ! this happens only on the rank which holds the current maxima
          time_av_sensitivity_max_idxs = maxloc(time_av_sensitivity)

          ! get current speaker (count)
          call set_parameter(s,'opt.active_monopoles')

          ! enable the speaker
          speaker_idx(s,7) = params%parallelism%block_image

          ! update speaker locate indices
          speaker_idx(s,4) = time_av_sensitivity_max_idxs(1)
          speaker_idx(s,5) = time_av_sensitivity_max_idxs(2)
          speaker_idx(s,6) = time_av_sensitivity_max_idxs(3)

          ! update speaker position
          speaker_pos(s,1) = X(time_av_sensitivity_max_idxs(1),time_av_sensitivity_max_idxs(2),time_av_sensitivity_max_idxs(3),1)
          speaker_pos(s,2) = X(time_av_sensitivity_max_idxs(1),time_av_sensitivity_max_idxs(2),time_av_sensitivity_max_idxs(3),2)
          speaker_pos(s,3) = X(time_av_sensitivity_max_idxs(1),time_av_sensitivity_max_idxs(2),time_av_sensitivity_max_idxs(3),3)

          if (params%io%verbosity.ge.1) then
             write(*,*) ""
             write(*,*) " maximum time avg. sensitivity:", max_all(1)
             write(*,*) " on image                     :", int(max_all(2))
             write(*,*) ""
             write(*,*) " -number  : ",s
             write(*,*) " -position: ",speaker_pos(s,1:3)
             write(*,*) " -indices : ",speaker_idx(s,4:6)
             write(*,*) ""
          end if
       end if

       call store(speaker_pos,'speaker_positions_loop')

    else
       
       if (params%parallelism%world_image.eq.1) then
          write(*,*) ""
          write(*,*) "warning in force_adjoint_monopole_synthesis.f90:update_speaker:can not add more speakers, maximum defined in p_target:",speaker_nbr
          write(*,*) ""
       end if
       
    end if

  end subroutine update_speaker
!!!=================================================================================================

end module force_adjoint_monopole_synthesis





