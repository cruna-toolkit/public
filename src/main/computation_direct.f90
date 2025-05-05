module computation_direct

  use discretisation_t
  use filter_x
  use force
  use helper
  use io
  use io_wrapper
  use parameter
  use time_series
  use boundary_conditions

contains

!!!=================================================================================================
  subroutine calc_direct(Q,q0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:), intent(inout) :: Q
    real(kind=rk), dimension(:,:,:,:)  , intent(in)    :: q0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:)  ,allocatable    :: qn,q1
    real                                               :: toc_rhs
    integer                                            :: nt
    integer                                            :: s
    integer                                            :: fi1freq,fi2freq,fi3freq
    logical                                            :: const_direct = .false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! get filter freqs
    call get_parameter(fi1freq,'filter.di1freq',default = huge(1))
    call get_parameter(fi2freq,'filter.di2freq',default = huge(1))
    call get_parameter(fi3freq,'filter.di3freq',default = huge(1))

    ! store initial condition
    call set_parameter(1,'time.ns',struct_val=params%time%ns)
    call set_parameter(0,'time.nt',struct_val=params%time%nt)
    call store(q0,'data_direct_snapshot')

    ! init
    allocate(q1(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    allocate(qn(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    qn           = q0
    Q(:,:,:,:,1) = q0 ! this sets Q=q0 in case of size(Q,5).eq.1 (linear adjoint)

    call get_parameter(const_direct,'opt.const_direct',default = .false.)
    if(const_direct.eqv..true.) then
       if(params%parallelism%block_image.eq.1) then
          write(*,*) ""
          write(*,*) "skip direct computation due to params.opt.const_direct = .true. "
          write(*,*) ""
       end if
       return
    endif

    ! loop subsets
    do s = 1,params%time%subsets
       !params%time%ns = s
       call set_parameter(s,'time.ns',struct_val=params%time%ns)

       ! load subset data force
       if(params%time%subsets.gt.1) then
          ! load actual force (actual subset is set above)
          call load_force_subsets_cache ! formerly load(f,'force_subsets_cache')
       end if

       ! loop timesteps
       do nt = 1,params%time%steps
          !params%time%nt = nt
          call set_parameter(nt,'time.nt',struct_val=params%time%nt)

          ! increment time
          !! params%time%t = params%time%t + params%time%dt ! this is not robust in case of optimisation
          params%time%t = ((s-1)*params%time%steps + nt) * params%time%dt

          call tic()

          ! call time-steper
          call time_stepper_direct(q1,qn,params%time%t)
          qn = q1

          ! rhs toc
          toc_rhs = toc()

          ! boundary handling
          call set_boundary_condition_q(qn,outside_rhs=.true.)

          ! filter
          if (mod(nt,fi1freq).eq.0) then
             call filter_i1(qn)
          end if

          if (mod(nt,fi2freq).eq.0) then
             call filter_i2(qn)
          end if

          if (mod(nt,fi3freq).eq.0) then
             call filter_i3(qn)
          end if

          ! save to memory
          if (size(Q,5).gt.1) then
             Q(:,:,:,:,nt) = qn;
          end if

          ! snapshot
          if (mod(nt,params%io%dsfreq).eq.0) then
             call store(qn,'data_direct_snapshot')
          end if

          ! time_series
          call sample(qn)
          call probe(qn)

          ! screen
          if(params%parallelism%world_image.eq.1) then
             if (mod(nt,params%io%sfreq).eq.0) then
                write(*,*) "computation (direct) of timestep ", trim(num2str(params%time%nt,'(I7.7)'))," in sub-set ", trim(num2str(params%time%ns,'(I7.7)')), " T/nt: ",trim(num2str(toc(),'(F00.3)'))," (",trim(num2str(toc_rhs,'(F0.3)')),")"
             end if
          end if

       end do

       if (params%time%subsets.gt.1) then
          if(params%parallelism%block_image.eq.1) then
             if (params%io%verbosity.ge.1) then
                write(*,*) " store data direct of subset", s
             end if
          end if

          ! store subset data
          call store(Q,'data_direct_subsets_cache')
       end if

    end do

    ! clean
    deallocate(q1)

  end subroutine calc_direct
!!!=================================================================================================

end module computation_direct
