module force_adjoint_sound_monopole_array

  use data_geom
  use io
  USE mpi
  use parameter
  use tub_ak_data

  real(kind=rk), dimension(:,:,:,:), allocatable, save, public :: f
  real(kind=rk), dimension(:,:,:)  , allocatable, save, public :: f_speaker, f_speaker_backup
  real(kind=rk), dimension(:,:,:)  , allocatable, save, public :: f_speaker_grad
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

contains

!!!=================================================================================================
  subroutine allocate_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call tub_ak_read_params

    allocate(f(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%time%steps))

    allocate(     f_speaker(params%time%steps,params%time%subsets,speaker_nbr))
    allocate(f_speaker_grad(params%time%steps,params%time%subsets,speaker_nbr))

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
       allocate(f_speaker_backup(params%time%steps,params%time%subsets,speaker_nbr))
    end if

    f_speaker_backup = f_speaker

  end subroutine backup_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine calc_grad(Qs,Q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:),intent(inout)       :: Qs     ! (inout) required for multiple subsets
    real(kind=rk), dimension(:,:,:,:,:),intent(in),optional :: Q      ! if regularization is applied
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:)         , allocatable       :: f_speaker_grad_transfer
    integer                                                 :: ns, ns_backup, s, ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(f_speaker_grad_transfer(size(f_speaker_grad,1)))

    ns_backup = params%time%ns
    do ns = 1,params%time%subsets

       if (params%time%subsets.gt.1) then
          if(params%parallelism%block_image.eq.1) then
             if (params%io%verbosity.ge.1) then
                write(*,*) " load data adjoint of subset", s
             end if
          end if

          call set_parameter(ns,'time.ns',struct_val=params%time%ns)
          call load(Qs,'data_adjoint_subsets_cache')
       end if

       do s = 1,speaker_nbr

          ! only the block in which the speaker is located evaluates the adjoint/gradient
          if (speaker_idx(s,7).ne.0) then
             f_speaker_grad(:,ns,s) = -Qs(speaker_idx(s,4),speaker_idx(s,5),speaker_idx(s,6),5,:)
          else
             f_speaker_grad(:,ns,s) = 0.0_rk
          end if

          ! mpi_reduce to update f_speaker_grad on all blocks
          call mpi_allreduce(f_speaker_grad(:,ns,s),f_speaker_grad_transfer,size(f_speaker_grad(:,ns,s)), mpi_double, mpi_sum, params%parallelism%block_comm, ierr)
          f_speaker_grad(:,ns,s) = f_speaker_grad_transfer
          ! now each block should have the temporal gradient of all speakers
          ! even if they are located in other blocks
       end do
    end do

    call set_parameter(ns_backup,'time.ns',struct_val=params%time%ns)

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
                f_speaker(:,:,i) = 0.0_rk
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
    integer                                               :: s,nt,ns,ns_backup
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

    ns_backup = params%time%ns
    do ns = 1,params%time%subsets
       do s = 1,speaker_nbr

          ! update speaker signal
          do nt = 1,params%time%steps
             f_speaker(nt,ns,s) = f_speaker(nt,ns,s) + alpha*f_speaker_grad(nt,ns,s)*theta_t(nt,ns)
          end do

          fx = fx_norm*exp(-fx_alpha*((X(:,:,:,1) - speaker_pos(s,1))**2 + (X(:,:,:,2) - speaker_pos(s,2))**2 + (X(:,:,:,3) - speaker_pos(s,3))**2))

          do nt = 1,params%time%steps
             f(:,:,:,nt) = f(:,:,:,nt) + fx*f_speaker(nt,ns,s)
          end do
       end do

       if (params%time%subsets.gt.1) then
          if(params%parallelism%block_image.eq.1) then
             if (params%io%verbosity.ge.1) then
                write(*,*) " store force of subset", s
             end if
          end if

          call set_parameter(ns,'time.ns',struct_val=params%time%ns)
          call store(f,'force_subsets_cache')
       end if
    end do
    call set_parameter(ns_backup,'time.ns',struct_val=params%time%ns)

    !call store(fx,'debug_field')

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
       theta_t_vector(1:discarded_steps) = 0.0_rk
       theta_t_vector((params%time%steps*params%time%subsets) - discarded_steps + 1 : (params%time%steps*params%time%subsets)) = 0.0_rk
    end if

    ! reshape to array
    theta_t = reshape(theta_t_vector,(/params%time%steps,params%time%subsets/))

    call store(theta_t,'theta_t_function')

  end subroutine init_theta_t
!!!=================================================================================================

end module force_adjoint_sound_monopole_array





