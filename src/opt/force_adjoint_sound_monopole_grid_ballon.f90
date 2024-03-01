module force_adjoint_sound_monopole_grid_ballon

  use data_geom
  use io
  USE mpi
  use parameter
  use tub_ak_data

  real(kind=rk), dimension(:,:,:,:,:), allocatable, save, public :: ballon,b_grad
  real(kind=rk), dimension(:,:,:,:)  , allocatable, save, public :: f,f_backup,f_grad
  real(kind=rk), dimension(:,:,:)    , allocatable, save, public :: theta_x
  real(kind=rk), dimension(:,:)      , allocatable, save, public :: theta_t
  real(kind=rk), dimension(:,:)      , allocatable, save, public :: drive,d_grad
  integer      ,                                    save, public :: freq_n
  logical      ,                                    save, public :: in_training_mode
  
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

    ! check if in ballon training mode
    call get_parameter(in_training_mode,'opt.in_training_mode',default=.true.)
    ! set number of frequency components, i.e., number of ballons
    call get_parameter(freq_n,'opt.freq_n',default=1)
    
    allocate(ballon(2,freq_n,params%geom%n1b,params%geom%n2b,params%geom%n3b))
    if (in_training_mode) then
       allocate(b_grad(2,freq_n,params%geom%n1b,params%geom%n2b,params%geom%n3b))
    end if
 
    allocate( drive(2,freq_n))
    allocate(d_grad(2,freq_n))

    ! 2 stands for a cos/sin coefficient -> in a future version complex ballons & FFT might more efficient 

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
       !f_grad = -Qs(:,:,:,5,:)
       if (in_training_mode) then
          call calc_b_grad(-Qs(:,:,:,5,:))
       else
          call calc_d_grad(-Qs(:,:,:,5,:))
       end if
       ! partial grad is written here instead of (non-existent) f_grad
       call get_parameter(loop,'opt.loop')
       call get_parameter(gfreq,'io.gfreq',default=huge(1))
       if (mod(loop,gfreq).eq.0) then
          ! call store(f_grad,'gradient_loop')
          if (in_training_mode) then
             call store(b_grad,'ballon_grad_loop')
          else
             call store(d_grad,'drive_grad_loop')
          end if
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
  subroutine calc_b_grad(ps)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),intent(in)            :: ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),allocatable,save      :: tmp
    real(kind=rk), dimension(:,:)    ,allocatable           :: omega_t
    real(kind=rk), dimension(:)      ,allocatable           :: oneomega_allt
    integer                                                 :: i1,i2,i3,nf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! b_grad is manuell FT of qs at distinct frequencies

    ! allocate/init theta_t
    ! if (.not.allocated(theta_t)) then
    !    call init_theta_t
    ! end if
    
    if (.not.allocated(tmp)) then
       allocate(tmp(params%time%steps,params%geom%n1b,params%geom%n2b,params%geom%n3b))
    end if

    ! shift time to first dimension
    tmp = reshape(ps, shape(tmp), order=(/2,3,4,1/))

    ! calc. ballon gradient / manuell FT
    call calc_omega_t(omega_t)
    allocate(oneomega_allt(size(omega_t,2)))
    do i3 = 1,params%geom%n3b
       do i2 = 1,params%geom%n2b
          do i1 = 1,params%geom%n1b
             do nf = 1,freq_n
                oneomega_allt = omega_t(nf,:)
                ! sum over all times
                b_grad(1,nf,i1,i2,i3) = sum( tmp(:,i1,i2,i3) * cos(oneomega_allt) )
                b_grad(2,nf,i1,i2,i3) = sum( tmp(:,i1,i2,i3) * sin(oneomega_allt) )
                ! b_grad(2,nf,i1,i2,i3) = sum( tmp(:,i1,i2,i3) * sin(oneomega_allt) * theta_t(:,params%time%ns) )
                ! line above: example incl. theta_t,
                ! however for FT current tophat window is bad -> better implement hann before use
             end do
          end do
       end do
    end do
    b_grad = b_grad/real(params%time%steps,rk) ! FT normalization
    
  end subroutine calc_b_grad
!!!=================================================================================================

!!!=================================================================================================
  subroutine calc_d_grad(ps)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),intent(in)            :: ps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:)  ,allocatable,save      :: tmp
    real(kind=rk), dimension(:,:)    ,allocatable           :: omega_t
    real(kind=rk), dimension(:)      ,allocatable           :: oneomega_allt
    integer                                                 :: i1,i2,i3,nf,nt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! d_grad is ballon-projection and manuell FT of qs (at distinct frequencies)

    ! allocate/init theta_t
    ! if (.not.allocated(theta_t)) then
    !    call init_theta_t
    ! end if
    
    if (.not.allocated(tmp)) then
       allocate(tmp(2,freq_n,params%time%steps))
    end if
    tmp = 0.0_rk

    !print*,'ext b',params%parallelism%block_image,maxval(reshape(ballon,(/size(ballon),1/))),minval(reshape(ballon,(/size(ballon),1/)))

    ! project all qs xyz-space onto ballons  
    do nt = 1,params%time%steps
       do i3 = 1,params%geom%n3b
          do i2 = 1,params%geom%n2b
             do i1 = 1,params%geom%n1b
                do nf = 1,freq_n
                   tmp(1,nf,nt) = tmp(1,nf,nt) + ps(i1,i2,i3,nt)*ballon(1,nf,i1,i2,i3)
                   tmp(2,nf,nt) = tmp(2,nf,nt) + ps(i1,i2,i3,nt)*ballon(2,nf,i1,i2,i3)
                end do
             end do
          end do
       end do
    end do

    !print*,'ext t',params%parallelism%block_image,maxval(reshape(tmp,(/size(tmp),1/))),minval(reshape(tmp,(/size(tmp),1/)))

    ! calc. drive gradient / manuell FT
    call calc_omega_t(omega_t)
    allocate(oneomega_allt(size(omega_t,2)))
    do i3 = 1,params%geom%n3b
       do i2 = 1,params%geom%n2b
          do i1 = 1,params%geom%n1b
             do nf = 1,freq_n
                oneomega_allt = omega_t(nf,:)
                ! sum over all times
                d_grad(1,nf) = sum( tmp(1,nf,:) * cos(oneomega_allt) )
                d_grad(2,nf) = sum( tmp(2,nf,:) * sin(oneomega_allt) )
                !d_grad(2,nf) = sum( tmp(2,nf,:) * sin(oneomega_allt) * theta_t(:,params%time%ns) )
                ! line above: example incl. theta_t,
                ! however for FT current tophat window is bad -> better implement hann before use
             end do
          end do
       end do
    end do
    d_grad = d_grad/real(params%time%steps,rk) ! FT normalization
    
  end subroutine calc_d_grad
!!!=================================================================================================

!!!=================================================================================================
  subroutine calc_omega_t(omega_t)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), allocatable, intent(out) :: omega_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:),   allocatable              :: omega_vec
    real(kind=rk)                                           :: twopi,frequency
    integer                                                 :: nf, nt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    twopi = 8.0_rk*atan(1.0_rk)

    ! get frequency vector
    ! in a future version possible build here dense FFT like freq. vector
    allocate(omega_vec(freq_n))
    do nf = 1,freq_n
       call get_parameter(frequency,'opt.freq_list',subindex=nf)
       omega_vec(nf) = twopi*frequency
    end do
      
    allocate(omega_t(freq_n,params%time%steps))
    do nt = 1,params%time%steps
       do nf = 1,freq_n
          ! warning subset support missing here
          omega_t(nf,nt) = omega_vec(nf) * real(nt,rk)*params%time%dt
       end do
    end do
    
  end subroutine calc_omega_t
!!!=================================================================================================
  
!!!=================================================================================================
  subroutine init_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                               :: s,s_backup
    integer                                               :: ballon_input_loop
    character(len=max_length_parameter)                   :: fname
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! init f with all zero
    f              = 0.0_rk
    f_grad         = 0.0_rk

    ! init. ballon driving function (unused if training mode)
    drive  = 0.0_rk
    d_grad = 0.0_rk
    ! if training mode init (and build) own ballons else import ballons
    if (in_training_mode) then
       ballon = 0.0_rk
       b_grad = 0.0_rk
    else ! production mode -> determining driving function
       if(params%parallelism%block_image.eq.1) then
          write(*,*) "load ballons"
       end if
       !call get_parameter(fname,'opt.init_ballon_file')      ! set file for ballon import
       !call load(ballon, 'dummy', fname_optin = trim(fname)) ! import ballon
       call get_parameter(ballon_input_loop,'opt.ballon_input_loop',default = 1)
       call set_parameter(ballon_input_loop,'opt.loop') ! set loop for ballon import
       call load(ballon,'ballon_loop')                  ! import ballon (has no time dependence)
       call set_parameter(1,'opt.loop')                 ! reset first loop nr
       if(params%parallelism%block_image.eq.1) then
          write(*,*) size(ballon,2), "ballons loaded"
          write(*,*) ""
       end if
    end if

    ! check for present pre-defined force (feature not used by lewin)
    call get_parameter(fname,'opt.init_force_file',default = '-')

    ! handle multiple subsets
    if (params%time%subsets.eq.1) then

       if (trim(fname).ne.'-') then
          if (params%parallelism%block_image.eq.1) then
             write(*,*) "warning in force_adjoint_sound_monopole_array.f90:init_force:load present force from file:", fname
          end if

          ! load pre-defined relevant force fraction from file, e.g. from optimization run
          if (in_training_mode) then
             call load(ballon, 'dummy', fname_optin = trim(fname))
          else
             call load(drive, 'dummy', fname_optin = trim(fname))
          end if
          
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
       write(*,*) "store force, and ballons or drive"
    end if

    if (params%time%subsets.eq.1) then
       call store(f,'force_loop')
       if (in_training_mode) then
          call store(ballon,'ballon_loop')
       else
          call store(drive,'drive_loop')
       end if
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
    ! force holds relevant degrees of freedom only - not anymore entire space-time

    if (params%time%subsets.eq.1) then
       ! do nt = 1,params%time%steps
       !    f(:,:,:,nt) = f(:,:,:,nt) + alpha*f_grad(:,:,:,nt)*theta_x*theta_t(nt,params%time%ns)
       ! end do
       if (in_training_mode) then
          call update_ballon_and_f
       else
          call update_drive_and_f
       end if
    else
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
  subroutine update_ballon_and_f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), allocatable              :: omega_t
    real(kind=rk)                                           :: alpha,alpha0
    integer                                                 :: i1,i2,i3,nf,nt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! update (train) ballon and calc. new f via manuell iFT

    ! allocate/init theta_t
    ! if (.not.allocated(theta_t)) then
    !    call init_theta_t
    ! end if
    
    ! allocate/init theta_x
    if (.not.allocated(theta_x)) then
       call init_theta_x
    end if
    
    ! get opt value
    call get_parameter(alpha0,'opt.alpha0'                 )
    call get_parameter(alpha ,'opt.alpha' ,default = alpha0)

    ! update ballon (using present alpha)
    do nf = 1,freq_n
       ballon(1,nf,:,:,:) = ballon(1,nf,:,:,:) + alpha*b_grad(1,nf,:,:,:)*theta_x
       ballon(2,nf,:,:,:) = ballon(2,nf,:,:,:) + alpha*b_grad(2,nf,:,:,:)*theta_x
    end do

    ! calc. new f / manuell iFT @ distinct frequencies
    call calc_omega_t(omega_t)
    do nt = 1,params%time%steps
       do i3 = 1,params%geom%n3b
          do i2 = 1,params%geom%n2b
             do i1 = 1,params%geom%n1b
                ! sum over all frequencies
                f(i1,i2,i3,nt) = sum( ballon(1,:,i1,i2,i3) * cos(omega_t(:,nt)) ) &
                               + sum( ballon(2,:,i1,i2,i3) * sin(omega_t(:,nt)) )
                !f(i1,i2,i3,nt) = f(i1,i2,i3,nt) * theta_t(nt,params%time%ns)
                ! the line above: example including theta_t,
                ! however for FT current tophat window is bad -> better implement hann before use
             end do
          end do
       end do
    end do
    
  end subroutine update_ballon_and_f
!!!=================================================================================================

!!!=================================================================================================
  subroutine update_drive_and_f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), allocatable              :: omega_t
    real(kind=rk)                                           :: alpha,alpha0
    integer                                                 :: i1,i2,i3,nf,nt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! update driving function and calc. new f via manuell iFT  
    ! theta_x is included in ballon already

    ! allocate/init theta_t
    ! if (.not.allocated(theta_t)) then
    !    call init_theta_t
    ! end if
    
    ! get opt value
    call get_parameter(alpha0,'opt.alpha0'                 )
    call get_parameter(alpha ,'opt.alpha' ,default = alpha0)

    ! update drive (using present alpha)
    do nf = 1,freq_n
       drive(1,nf) = drive(1,nf) + alpha*d_grad(1,nf)
       drive(2,nf) = drive(2,nf) + alpha*d_grad(2,nf)
    end do

    ! calc. new f / manuell iFT @ distinct frequencies
    call calc_omega_t(omega_t)
    do nt = 1,params%time%steps
       do i3 = 1,params%geom%n3b
          do i2 = 1,params%geom%n2b
             do i1 = 1,params%geom%n1b
                ! sum over all frequencies
                f(i1,i2,i3,nt) = sum( drive(1,:) * ballon(1,:,i1,i2,i3) * cos(omega_t(:,nt)) ) &
                               + sum( drive(2,:) * ballon(2,:,i1,i2,i3) * sin(omega_t(:,nt)) )
                !f(i1,i2,i3,nt) = f(i1,i2,i3,nt) * theta_t(nt,params%time%ns)
                ! the line above: example including theta_t,
                ! however for FT current tophat window is bad -> better implement hann before use
             end do
          end do
       end do
    end do
    
  end subroutine update_drive_and_f
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
    real(kind=rk)                                                :: src_radius, src_flank_width
    real(kind=rk)                                                :: pi
    integer                                                      :: s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! constants
    pi = 4.0_rk*atan(1.0_rk)

    ! parameter
    call get_parameter(src_radius     ,'opt.src_radius'     ,default = 0.10_rk)
    call get_parameter(src_flank_width,'opt.src_flank_width',default = 0.04_rk)
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

end module force_adjoint_sound_monopole_grid_ballon


