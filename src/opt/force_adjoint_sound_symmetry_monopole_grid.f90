module force_adjoint_sound_symmetry_monopole_grid

  use data_geom
  use io
  USE mpi
  use parameter
  use tub_ak_data

  real(kind=rk), dimension(:,:,:,:), allocatable, save, public :: f,f_backup,f_grad
  real(kind=rk), dimension(:,:,:)  , allocatable, save, public :: theta_x
  real(kind=rk), dimension(:,:,:)  , allocatable, save, public :: sym_box,sym_box_tmp
  real(kind=rk), dimension(:)      , allocatable, save, public :: theta_t
  integer      , dimension(6)      , save, public              :: idx_SinB, idx_BinS

  private

  public allocate_force
  public apply_force
  public backup_force
  public calc_grad
  public init_force
  public store_force
  public restore_force
  public update_force
  public symmetrice_f_grad
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
    if (.not.allocated(f_backup)) then
       allocate(f_backup(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%time%steps))
    end if

    f_backup = f

  end subroutine backup_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine calc_grad(Qs,Q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:),intent(in)          :: Qs
    real(kind=rk), dimension(:,:,:,:,:),intent(in),optional :: Q      ! if regularization is applied
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                 :: loop,gfreq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    f_grad = -Qs(:,:,:,5,:)

    call symmetrice_f_grad

    if(params%parallelism%block_image.eq.1) then
       ! all sub block own information of all speaker, see above
       write(*,*) "write gradient"
    end if

    call get_parameter(loop,'opt.loop')
    call get_parameter(gfreq,'io.gfreq',default=huge(1))
    if (mod(loop,gfreq).eq.0) then
       call store(f_grad,'gradient_loop')
    end if

  end subroutine calc_grad
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_parameter)                   :: fname
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! init f with all zero
    f              = 0.0_rk
    f_grad         = 0.0_rk

    ! load present force from file
    call get_parameter(fname,'opt.init_force_file',default = '-')

    if (trim(fname).ne.'-') then                                        ! parameter is not the default value, there is something to do

       if (params%parallelism%block_image.eq.1) then
          write(*,*) "warning in force_adjoint_sound_monopole_array.f90:init_force:load present force from file:", fname
       end if

       call load(f, 'dummy', fname_optin = trim(fname))                 ! load present force, e.g. from optimization run

    end if

    ! allocate/init theta_x
    if (.not.allocated(theta_x)) then
       call init_theta_x
    end if

  end subroutine init_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine store_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(params%parallelism%block_image.eq.1) then
       ! all sub block own information of all speaker, see above
       write(*,*) "write force"
    end if

    call store(f,'force_loop')
    
  end subroutine store_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine restore_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    f = f_backup

  end subroutine restore_force
!!!=================================================================================================

!!!=================================================================================================
  subroutine update_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                         :: alpha,alpha0
    integer                                               :: nt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! allocate/init theta_t
    if (.not.allocated(theta_t)) then
       call init_theta_t
    end if

    ! allocate/init theta_x
    !if (.not.allocated(theta_x)) then
    !   call init_theta_x
    !end if

    ! get opt value
    call get_parameter(alpha0,'opt.alpha0'                 )
    call get_parameter(alpha ,'opt.alpha' ,default = alpha0)

    do nt = 1,params%time%steps
       f(:,:,:,nt) = f(:,:,:,nt) + alpha*f_grad(:,:,:,nt)*theta_x*theta_t(nt)
    end do

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

    !do nt = 1,discarded_steps
    !   theta_t(nt) = 0.0_rk
    !end do

    do  nt = params%time%steps-discarded_steps+1,params%time%steps
       theta_t(nt) = 0.0_rk
    end do

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

    call store(theta_t,'theta_t_function')

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
!    integer                                                      :: ierr
    integer                                                      :: sym_box_n, sym_box_nhalf
    integer                                                      :: sym_box_gis,sym_box_gie,sym_box_gjs,sym_box_gje,sym_box_gks,sym_box_gke
    integer                                                      :: common_gis,common_gie,common_gjs,common_gje,common_gks,common_gke
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! constants
    pi = 4.0_rk*atan(1.0_rk)

    ! parameter
    call get_parameter(src_radius_dx1fac     ,'opt.src_radius_dx1fac'     ,default = 10)
    call get_parameter(src_flank_width_dx1fac,'opt.src_flank_width_dx1fac',default =  5)

    src_radius      = src_radius_dx1fac     *params%geom%dx1
    src_flank_width = src_flank_width_dx1fac*params%geom%dx1

    ! allocate
    allocate(theta_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    theta_x = 0.0_rk

    allocate(src_pos_distance(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    src_pos_distance = 0.0_rk

    ! flank width is defined as x-length between 0 and 1 of theta_x in linear approximation:
    ! df/dx = slope = Dx[0.5*erf(sqrt(pi)*x/len)]= 1/len -> since df=1 it follows dx = len 

    do s = 1,speaker_nbr
       ! compute distance
       src_pos_distance = sqrt((X(:,:,:,1) - speaker_pos(s,1))**2 + (X(:,:,:,2) - speaker_pos(s,2))**2 + (X(:,:,:,3) - speaker_pos(s,3))**2)

       ! compute theta_x
       theta_x = theta_x + (0.5_rk - 0.5_rk*erf(sqrt(pi)/src_flank_width*(src_pos_distance - src_radius)))
    end do

    call store(theta_x,'theta_x_function')

    ! create sym_box below:
    ! relevant part of f_grad is copied to "symmetrization box" (and back)
    ! to enforce spatial symmetries of the forcing
    ! all sizes of sym_box are always odd & speaker center is in its middle

    ! sym_box size
    sym_box_n = ceiling(2.0_rk*(src_radius+3.0_rk*src_flank_width)/params%geom%dx1)
    if (mod(sym_box_n,2).eq.0) then ! if even
       sym_box_n = sym_box_n+1
    end if
    sym_box_nhalf = (sym_box_n-1)/2

    ! allocate
    allocate(sym_box(sym_box_n,sym_box_n,sym_box_n))
    allocate(sym_box_tmp(sym_box_n,sym_box_n,sym_box_n))

    ! sym_box indices (global basis)
    sym_box_gis = speaker_idx(1,1)-sym_box_nhalf
    sym_box_gie = speaker_idx(1,1)+sym_box_nhalf
    sym_box_gjs = speaker_idx(1,2)-sym_box_nhalf
    sym_box_gje = speaker_idx(1,2)+sym_box_nhalf
    sym_box_gks = speaker_idx(1,3)-sym_box_nhalf
    sym_box_gke = speaker_idx(1,3)+sym_box_nhalf

    ! common indices of sym_box and local MPI block (global basis)
    common_gis = max(sym_box_gis,Xi(1,1,1,1))
    common_gie = min(sym_box_gie,Xi(params%geom%n1b,1,1,1))
    common_gjs = max(sym_box_gjs,Xi(1,1,1,2))
    common_gje = min(sym_box_gje,Xi(1,params%geom%n2b,1,2))
    common_gks = max(sym_box_gks,Xi(1,1,1,3))
    common_gke = min(sym_box_gke,Xi(1,1,params%geom%n3b,3))

    ! common indices (local/block basis)
    idx_SinB(1) = common_gis + 1 - Xi(1,1,1,1)
    idx_SinB(2) = common_gie + 1 - Xi(1,1,1,1)
    idx_SinB(3) = common_gjs + 1 - Xi(1,1,1,2)
    idx_SinB(4) = common_gje + 1 - Xi(1,1,1,2)
    idx_SinB(5) = common_gks + 1 - Xi(1,1,1,3)
    idx_SinB(6) = common_gke + 1 - Xi(1,1,1,3)

    ! common indices (sym_box basis)
    idx_BinS(1) = common_gis + 1 - sym_box_gis
    idx_BinS(2) = common_gie + 1 - sym_box_gis
    idx_BinS(3) = common_gjs + 1 - sym_box_gjs
    idx_BinS(4) = common_gje + 1 - sym_box_gjs
    idx_BinS(5) = common_gks + 1 - sym_box_gks
    idx_BinS(6) = common_gke + 1 - sym_box_gks

    ! debug / check

    !print*,'spe',params%parallelism%block_image,speaker_idx(1,1:3),speaker_pos(1,1),params%geom%dx1
    !print*,'box',params%parallelism%block_image,sym_box_gis,sym_box_gie,sym_box_gjs,sym_box_gje,sym_box_gks,sym_box_gke
    !print*,'c_g',params%parallelism%block_image,common_gis,common_gie,common_gjs,common_gje,common_gks,common_gke
    !print*,'c_l',params%parallelism%block_image,idx_SinB
    !print*,'c_s',params%parallelism%block_image,idx_BinS

    ! fill sym_box
    !sym_box(idx_BinS(1):idx_BinS(2),idx_BinS(3):idx_BinS(4),idx_BinS(5):idx_BinS(6)) = &
    !     theta_x(idx_SinB(1):idx_SinB(2),idx_SinB(3):idx_SinB(4),idx_SinB(5):idx_SinB(6))

    !print*,'a sym_box x',params%parallelism%block_image,sym_box(24:28,26,26)
    !print*,'a sym_box y',params%parallelism%block_image,sym_box(26,24:28,26)
    !print*,'a sym_box z',params%parallelism%block_image,sym_box(26,26,24:28)

    ! sum sym_box content of all images
    !call mpi_allreduce(sym_box, sym_box_tmp, size(sym_box), mpi_double, mpi_sum, params%parallelism%block_comm, ierr)
    !sym_box = sym_box_tmp

    !print*,'b sym_box x',params%parallelism%block_image,sym_box(24:28,26,26)
    !print*,'b sym_box y',params%parallelism%block_image,sym_box(26,24:28,26)
    !print*,'b sym_box z',params%parallelism%block_image,sym_box(26,26,24:28)

    ! force x symmetry
    !sym_box_tmp = sym_box(sym_box_n:1:-1,:,:)    ! flip x
    !sym_box     = 0.5_rk*(sym_box + sym_box_tmp) ! mean

    ! force y symmetry
    !sym_box_tmp = sym_box(:,sym_box_n:1:-1,:)    ! flip y
    !sym_box     = 0.5_rk*(sym_box + sym_box_tmp) ! mean

    !print*,'c sym_box x',params%parallelism%block_image,sym_box(24:28,26,26)
    !print*,'c sym_box y',params%parallelism%block_image,sym_box(26,24:28,26)
    !print*,'c sym_box z',params%parallelism%block_image,sym_box(26,26,24:28)

  end subroutine init_theta_x
!!!=================================================================================================

!!!=================================================================================================
  subroutine symmetrice_f_grad
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                               :: ierr
    integer                                               :: sym_box_n,sym_box_nhalf,nt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! allocate/init theta_t
    if (.not.allocated(sym_box)) then
       call init_theta_x
    end if

    sym_box_n     = size(sym_box,1)
    sym_box_nhalf = (sym_box_n-1)/2 ! sym_box_n is always odd

    do nt = 1,params%time%steps

       ! clear sym_box
       sym_box     = 0.0_rk
       sym_box_tmp = 0.0_rk

       ! fill sym_box with local f_grad
       sym_box(idx_BinS(1):idx_BinS(2),idx_BinS(3):idx_BinS(4),idx_BinS(5):idx_BinS(6)) = &
            f_grad(idx_SinB(1):idx_SinB(2),idx_SinB(3):idx_SinB(4),idx_SinB(5):idx_SinB(6),nt)

       ! sum sym_box content of all images
       call mpi_allreduce(sym_box, sym_box_tmp, size(sym_box), mpi_double, mpi_sum, params%parallelism%block_comm, ierr)
       sym_box = sym_box_tmp

       ! force x symmetry (incl. positive x & y quarter only)
       ! i.e., assume theta_x "lives" in positive x and positive y quarter
       sym_box(1:sym_box_nhalf,:,:) = sym_box(sym_box_n:sym_box_n-sym_box_nhalf+1:-1,:,:) ! flip x
       sym_box(:,1:sym_box_nhalf,:) = sym_box(:,sym_box_n:sym_box_n-sym_box_nhalf+1:-1,:) ! flip y

       ! force x symmetry (incl. all quadrants)
       !sym_box_tmp = sym_box(sym_box_n:1:-1,:,:)    ! flip x
       !sym_box     = 0.5_rk*(sym_box + sym_box_tmp) ! mean

       ! force y symmetry (incl. all quadrants)
       !sym_box_tmp = sym_box(:,sym_box_n:1:-1,:)    ! flip y
       !sym_box     = 0.5_rk*(sym_box + sym_box_tmp) ! mean

       ! overwrite local f_grad with symmetriced form
       f_grad(idx_SinB(1):idx_SinB(2),idx_SinB(3):idx_SinB(4),idx_SinB(5):idx_SinB(6),nt) = &
            sym_box(idx_BinS(1):idx_BinS(2),idx_BinS(3):idx_BinS(4),idx_BinS(5):idx_BinS(6))

    end do

  end subroutine symmetrice_f_grad
!!!=================================================================================================

end module force_adjoint_sound_symmetry_monopole_grid





