module computation_adjoint

  use boundary_conditions
  use discretisation_t
  use filter_x
  use helper
  use io
  use io_wrapper
  use objective
  use parameter

contains

!!!=================================================================================================
  subroutine calc_adjoint(Qs,qs0,Q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:), intent(inout) :: Qs
    real(kind=rk), dimension(:,:,:,:,:), intent(inout) :: Q     ! inout is required in case of multiple subsets load
    real(kind=rk), dimension(:,:,:,:)  , intent(in)    :: qs0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:)  , allocatable   :: qsn,qs1,q0,g
    real                                               :: toc_rhs
    integer                                            :: nt
    integer                                            :: s
    integer                                            :: fi1freq,fi2freq,fi3freq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! get filter freqs
    call get_parameter(fi1freq,'filter.ai1freq',default = huge(1))
    call get_parameter(fi2freq,'filter.ai2freq',default = huge(1))
    call get_parameter(fi3freq,'filter.ai3freq',default = huge(1))

    ! store initial condition
    call set_parameter(params%time%steps + 1,'time.nt',struct_val=params%time%nt)
    call store(qs0,'data_adjoint_snapshot')

    ! init
    allocate(qs1(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    allocate(qsn(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    allocate( q0(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    allocate(  g(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    qsn = qs0

    ! loop subsets
    do s = params%time%subsets,1,-1
       call set_parameter(s,'time.ns',struct_val=params%time%ns)

       ! load subset data direct
       if(params%time%subsets.gt.1) then
          if(params%parallelism%block_image.eq.1) then
             if (params%io%verbosity.ge.1) then
                write(*,*) " load data direct of subset", s
             end if
          end if

          ! load actual data direct (actual subset is set above)
          call load(Q,'data_direct_subsets_cache')
       end if

       ! loop timesteps
       do nt = params%time%steps,1,-1
          call set_parameter(nt,'time.nt',struct_val=params%time%nt)

          ! increment time
          ! params%time%t = params%time%t - params%time%dt ! this is not reboust in case of optimisation
          params%time%t = ((s-1)*params%time%steps + nt) * params%time%dt

          ! load direct solution and compute g
          if (size(Q,5).gt.1) then
             q0 = Q(:,:,:,:,nt)
          else
             q0 = Q(:,:,:,:,1)
          end if
          
          call calc_objective_g(g,q0)

          call tic()

          ! call time-steper
          call time_stepper_adjoint(qs1,qsn,q0,g,params%time%t)
          qsn = qs1

          ! rhs toc
          toc_rhs = toc()

          ! boundary handling
          call set_boundary_condition_qs(qsn,q0)

          ! filter
          if (mod(nt,fi1freq).eq.0) then
             call filter_i1(qsn)
          end if

          if (mod(nt,fi2freq).eq.0) then
             call filter_i2(qsn)
          end if

          if (mod(nt,fi3freq).eq.0) then
             call filter_i3(qsn)
          end if

          ! save to memory
          Qs(:,:,:,:,nt) = qsn;

          ! snapshot
          if (mod(nt,params%io%asfreq).eq.0) then
             call store(qsn,'data_adjoint_snapshot')
          end if

          ! screen
          if(params%parallelism%world_image.eq.1) then
             if (mod(nt,params%io%sfreq).eq.0) then
                write(*,*) "computation (adjoint) of timestep ", trim(num2str(params%time%nt,'(I6.6)'))," in sub-set ", trim(num2str(params%time%ns,'(I6.6)')), " T/nt: ", &
                     trim(num2str(toc(),'(F0.3)'))," (",trim(num2str(toc_rhs,'(F0.3)')),")"
             end if
          end if

       end do

       if (params%time%subsets.gt.1) then
          if(params%parallelism%block_image.eq.1) then
             if (params%io%verbosity.ge.1) then
                write(*,*) " store data adjoint of subset", s
             end if
          end if
          
          ! store adjoint subset data
          call store(Qs,'data_adjoint_subsets_cache')
       end if

    end do

    ! clean
    deallocate(qs1)
    deallocate(qsn)
    deallocate( q0)
    deallocate(  g)

  end subroutine calc_adjoint
!!!=================================================================================================

end module computation_adjoint
