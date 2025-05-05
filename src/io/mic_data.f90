module mic_data

  use data_geom
  use helper
  use io
  use parameter

  private

  real(kind=rk), public, dimension(:,:,:,:), allocatable, save :: mic_distance
  real(kind=rk), public, dimension(:,:)    , allocatable, save :: mic_pos
  real(kind=rk), public, dimension(:,:)    , allocatable, save :: mic_sig
  integer      , public                                 , save :: mic_nbr

  public  :: give_opt_target
  private :: read_params_mic

contains

!!!=================================================================================================
  subroutine give_opt_target(opt_target)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(out)          :: opt_target
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), save                                   :: mic_offset
    integer                                               :: i,j,k,m
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.allocated(mic_pos)) then
       ! read mic params
       call read_params_mic

       ! get offset stuff
       call get_parameter(mic_offset,'opt.mic_offset',default = 5.0_rk)
       mic_offset = mic_offset * Params%geom%dx1
    end if

    ! create p_target
    opt_target = 0.0_rk

    ! set mic data
    do m = 1,mic_nbr
       do i = 1,params%geom%n1b
          do j = 1,params%geom%n2b
             do k = 1,params%geom%n3b
                if (mic_distance(i,j,k,m).le.mic_offset) then
                   opt_target(i,j,k) = mic_sig(params%time%nt + (params%time%ns-1)*params%time%steps, m)
                end if
             end do
          end do
       end do
    end do

  end subroutine give_opt_target
!!!=================================================================================================

!!!=================================================================================================
  subroutine read_params_mic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                         :: dx1,dx2,dx3
    real(kind=rk), dimension(1)                           :: mic_samplerate, mic_pos_case
    real(kind=rk)                                         :: mic_pos_norm
    real(kind=rk)                                         :: delta = 1e-12
    integer                                               :: i,j,k,m
    integer      , dimension(2)                           :: mic_sig_dims,mic_pos_dims
    character(len=max_length_fname)                       :: fname
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_parameter(fname,'opt.p_target_file')

    ! read parameter
    call get_parameter(fname,'opt.p_target_file')

    call get_dset_dim(mic_sig_dims, 'dummy', fname_optin = trim(adjustl(fname)), dset_name_optin = 'mic_sig')
    call get_dset_dim(mic_pos_dims, 'dummy', fname_optin = trim(adjustl(fname)), dset_name_optin = 'mic_pos')

    call load(mic_samplerate,'dummy', fname_optin = trim(adjustl(fname)), dset_name_optin = 'mic_samplerate')
    call load(mic_pos_case  ,'dummy', fname_optin = trim(adjustl(fname)), dset_name_optin = 'mic_pos_case'  )

    mic_nbr = mic_sig_dims(2)

    ! read data
    allocate(mic_pos(mic_pos_dims(1),mic_pos_dims(2)))
    call load(mic_pos,'dummy', fname_optin = trim(adjustl(fname)), dset_name_optin = 'mic_pos')

    ! norm mic pos data
    call get_parameter(mic_pos_norm,'opt.mic_pos_norm',default = 1.0_rk)
    mic_pos = mic_pos / mic_pos_norm

    ! output
    if (params%parallelism%world_image.eq.1) then

       write(*,*) "target: tub_rar_data read data from file ", trim(fname)
       write(*,*) "      : number of microphones         : " , mic_nbr
       write(*,*) "      : mic setup case  (always check): " , mic_pos_case

       if (params%io%verbosity.ge.1) then
          write(*,*) ""
          write(*,*) " number of mics: ", mic_nbr
          write(*,*) " samples       : ", mic_sig_dims(1)
          write(*,*) " sample rate   : ", mic_samplerate
          write(*,*) " mic pos norm  : ", mic_pos_norm
          write(*,*) " mic positions : "
          do i = 1,mic_pos_dims(1)
             write(*,*) " #",num2str(i,'(I4.4)'),mic_pos(i,:)
          end do
          write(*,*) ""
       end if

       ! check parameter
       if (mic_sig_dims(1).lt.params%time%steps) then
          write(*,*) "warning in tub_rar_data.f90:tub_rar_read_data:signal length shorter than time.steps:fill with zeros",mic_sig_dims(1),params%time%steps
       end if

       if (mic_sig_dims(2).ne.mic_pos_dims(1)) then
          write(*,*) "warning in tub_rar_data.f90:tub_rar_read_data:dimension mismatch between number of mics in signal and position data",mic_sig_dims(2),mic_pos_dims(1)
       end if

       if (.not.(((1.0_rk/mic_samplerate(1) - delta).lt.params%time%dt).and.((1.0_rk/mic_samplerate(1) + delta).gt.params%time%dt))) then
          write(*,*) "warning in tub_rar_data.f90:tub_rar_read_data:sampling frequencies do not match",mic_samplerate,params%time%dt
       end if

    end if

    ! mic distance
    allocate(mic_distance(params%geom%n1b,params%geom%n2b,params%geom%n3b,mic_nbr))
    do m = 1,mic_nbr
       do k = 1,params%geom%n3b
          do j = 1,params%geom%n2b
             do i = 1,params%geom%n1b
                dx1 = X(i,j,k,1) - mic_pos(m,1)
                dx2 = X(i,j,k,2) - mic_pos(m,2)
                dx3 = X(i,j,k,3) - mic_pos(m,3)

                mic_distance(i,j,k,m) = sqrt(dx1**2 + dx2**2 + dx3**2)
             end do
          end do
       end do
    end do

    ! read data (allocate)
    allocate(mic_sig(mic_sig_dims(1),mic_sig_dims(2)))

    ! read mic data
    mic_sig = 0.0_rk
    call load(mic_sig,'dummy', fname_optin = trim(adjustl(fname)), dset_name_optin = 'mic_sig')

  end subroutine read_params_mic
!!!=================================================================================================

end module mic_data
