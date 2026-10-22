!> Adjoint computation driver.
!!
!! This module advances the adjoint state backward through all configured subsets and
!! timesteps. It coordinates direct-state loading, objective-gradient evaluation,
!! adjoint time integration, boundary-condition handling, filtering, snapshot output,
!! time-series sampling, and progress reporting.
!!
!! The implementation relies on the shared runtime parameter structure `params` for
!! geometry, time integration, filtering, parallelism, and I/O configuration.
module computation_adjoint

  use boundary_conditions
  use discretisation_t
  use filter_x
  use helper
  use io
  use io_helper
  use io_wrapper
  use objective
  use parameter
  use time_series

  ! All routines below are module procedures.
contains

!!!=================================================================================================
  !> Compute the adjoint solution backward in time.
  !!
  !! The routine starts from the supplied terminal adjoint state `qs0` and marches
  !! backward over all configured subsets and timesteps. For each timestep, it obtains
  !! the corresponding direct solution, computes the objective contribution `g`,
  !! advances the adjoint state by one step, applies boundary conditions and optional
  !! filters, stores the result in memory, and triggers configured output operations.
  !!
  !! When multiple subsets are used, the direct solution for the current subset is
  !! loaded from cache and the next subset is pre-read to overlap I/O with computation.
  !! Adjoint subset data are written back after each subset has been processed.
  !!
  !! **Arguments**
  !!
  !! * `Qs` - Five-dimensional storage for the computed adjoint solution over time.
  !!   The current adjoint field is written into the fifth dimension at each timestep.
  !! * `qs0` - Initial state for the backward adjoint integration, i.e. the terminal
  !!   adjoint condition supplied by the caller.
  !! * `Q` - Direct solution data used to evaluate the objective contribution and the
  !!   adjoint time step. The array is `intent(inout)` because multiple-subset runs
  !!   reload its contents from cached direct-solution data.
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
    ! Read the filter application frequencies; a missing setting effectively disables filtering.
    call get_parameter(fi1freq,'filter.ai1freq',default = huge(1))
    call get_parameter(fi2freq,'filter.ai2freq',default = huge(1))
    call get_parameter(fi3freq,'filter.ai3freq',default = huge(1))

    ! store initial condition
    ! Move the global timestep parameter to the terminal state before writing the initial snapshot.
    call set_parameter(params%time%steps + 1,'time.nt',struct_val=params%time%nt)
    ! Persist the terminal adjoint condition as the first adjoint snapshot.
    call store(qs0,'data_adjoint_snapshot')

    ! init
    ! Allocate working fields using the local block geometry and number of equation variables.
    allocate(qs1(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    allocate(qsn(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    allocate( q0(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    allocate(  g(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
    ! Initialize the current adjoint state with the supplied terminal condition.
    qsn = qs0

    ! loop subsets
    ! Traverse subsets in reverse order because the adjoint integration proceeds backward in time.
    do s = params%time%subsets,1,-1
       ! Publish the active subset index through the global parameter structure.
       call set_parameter(s,'time.ns',struct_val=params%time%ns)

       ! load subset data direct
       if(params%time%subsets.gt.1) then
          if(params%parallelism%block_image.eq.1) then
             if (params%io%verbosity.ge.1) then
                write(*,*) " load data direct of subset", s
             end if
          end if

          ! load actual data direct (actual subset is set above)
          ! Load the direct solution corresponding to the current subset into the reusable buffer.
          call load(Q,'data_direct_subsets_cache')

          ! pre-read next subset
          if (s.gt.1) then
             ! Start pre-reading the next subset needed by the backward traversal.
             call pre_read('data_direct_subsets_cache',ns=s-1)
          end if

       end if

       ! loop timesteps
       ! Process timesteps of the active subset in reverse chronological order.
       do nt = params%time%steps,1,-1
          ! Publish the active timestep through the global parameter structure.
          call set_parameter(nt,'time.nt',struct_val=params%time%nt)

          ! increment time
          ! params%time%t = params%time%t - params%time%dt ! this is not reboust in case of optimisation
          ! Reconstruct physical time from subset index and local timestep instead of decrementing cumulatively.
          params%time%t = ((s-1)*params%time%steps + nt) * params%time%dt

          ! load direct solution and compute g
          ! Select the direct state for the current timestep; single-slice buffers reuse index 1.
          if (size(Q,5).gt.1) then
             q0 = Q(:,:,:,:,nt)
          else
             q0 = Q(:,:,:,:,1)
          end if

          ! Evaluate the objective forcing/source term required by the adjoint equations.
          call calc_objective_g(g,q0)

          ! Start timing the adjoint right-hand-side/time-step computation.
          call tic()

          ! call time-steper
          ! Advance the adjoint solution by one backward step using the direct state and objective term.
          call time_stepper_adjoint(qs1,qsn,q0,g,params%time%t)
          ! Promote the newly computed adjoint state to the current state for the next iteration.
          qsn = qs1

          ! rhs toc
          ! Store the elapsed time of the adjoint time-step calculation for diagnostics.
          toc_rhs = toc()

          ! boundary handling
          ! Apply adjoint boundary conditions after the time-step update.
          call set_boundary_condition_qs(qsn,q0,outside_rhs=.true.)

          ! filter
          ! Apply the first directional filter at its configured timestep frequency.
          if (mod(nt,fi1freq).eq.0) then
             call filter_i1(qsn)
          end if

          ! Apply the second directional filter at its configured timestep frequency.
          if (mod(nt,fi2freq).eq.0) then
             call filter_i2(qsn)
          end if

          ! Apply the third directional filter at its configured timestep frequency.
          if (mod(nt,fi3freq).eq.0) then
             call filter_i3(qsn)
          end if

          ! save to memory
          ! Store the current adjoint field in the in-memory time-history array.
          Qs(:,:,:,:,nt) = qsn;

          ! snapshot
          ! Write regular adjoint snapshots according to the configured snapshot frequency.
          if (mod(nt,params%io%asfreq).eq.0) then
             call store(qsn,'data_adjoint_snapshot')
          end if

          ! snapshot with loop number
          ! Write loop-indexed adjoint snapshots according to the configured loop snapshot frequency.
          if (mod(nt,params%io%aslfreq).eq.0) then
             call store(qsn,'data_adjoint_snapshot_loop')
          end if

          ! time_series
          ! Sample configured adjoint time-series quantities from the current state.
          call sample_adj(qsn)
          ! Evaluate configured adjoint probes for the current state.
          call probe_adj(qsn)

          ! screen
          ! Restrict progress output to the global root image to avoid duplicate terminal messages.
          if(params%parallelism%world_image.eq.1) then
             if (mod(nt,params%io%sfreq).eq.0) then
                write(*,*) "computation (adjoint) of timestep ", trim(num2str(params%time%nt,'(I7.7)'))," in sub-set ", trim(num2str(params%time%ns,'(I7.7)')), " T/nt: ", &
                     trim(num2str(toc(),'(F0.3)'))," (",trim(num2str(toc_rhs,'(F0.3)')),")"
             end if
          end if

       end do

       ! Multi-subset runs persist the completed adjoint subset before moving to the preceding subset.
       if (params%time%subsets.gt.1) then
          if(params%parallelism%block_image.eq.1) then
             if (params%io%verbosity.ge.1) then
                write(*,*) " store data adjoint of subset", s
             end if
          end if

          ! store adjoint subset data
          ! Cache the computed adjoint data for the completed subset.
          call store(Qs,'data_adjoint_subsets_cache')
       end if

    end do

    ! clean
    ! Release all temporary working arrays allocated by this routine.
    deallocate(qs1)
    deallocate(qsn)
    deallocate( q0)
    deallocate(  g)

  end subroutine calc_adjoint
!!!=================================================================================================

end module computation_adjoint

! module computation_adjoint

!   use boundary_conditions
!   use discretisation_t
!   use filter_x
!   use helper
!   use io
!   use io_helper
!   use io_wrapper
!   use objective
!   use parameter
!   use time_series

! contains

! !!!=================================================================================================
!   subroutine calc_adjoint(Qs,qs0,Q)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     real(kind=rk), dimension(:,:,:,:,:), intent(inout) :: Qs
!     real(kind=rk), dimension(:,:,:,:,:), intent(inout) :: Q     ! inout is required in case of multiple subsets load
!     real(kind=rk), dimension(:,:,:,:)  , intent(in)    :: qs0
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     real(kind=rk), dimension(:,:,:,:)  , allocatable   :: qsn,qs1,q0,g
!     real                                               :: toc_rhs
!     integer                                            :: nt
!     integer                                            :: s
!     integer                                            :: fi1freq,fi2freq,fi3freq
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ! get filter freqs
!     call get_parameter(fi1freq,'filter.ai1freq',default = huge(1))
!     call get_parameter(fi2freq,'filter.ai2freq',default = huge(1))
!     call get_parameter(fi3freq,'filter.ai3freq',default = huge(1))

!     ! store initial condition
!     call set_parameter(params%time%steps + 1,'time.nt',struct_val=params%time%nt)
!     call store(qs0,'data_adjoint_snapshot')

!     ! init
!     allocate(qs1(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
!     allocate(qsn(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
!     allocate( q0(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
!     allocate(  g(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
!     qsn = qs0

!     ! loop subsets
!     do s = params%time%subsets,1,-1
!        call set_parameter(s,'time.ns',struct_val=params%time%ns)

!        ! load subset data direct
!        if(params%time%subsets.gt.1) then
!           if(params%parallelism%block_image.eq.1) then
!              if (params%io%verbosity.ge.1) then
!                 write(*,*) " load data direct of subset", s
!              end if
!           end if

!           ! load actual data direct (actual subset is set above)
!           call load(Q,'data_direct_subsets_cache')

!           ! pre-read next subset
!           if (s.gt.1) then
!              call pre_read('data_direct_subsets_cache',ns=s-1)
!           end if

!        end if

!        ! loop timesteps
!        do nt = params%time%steps,1,-1
!           call set_parameter(nt,'time.nt',struct_val=params%time%nt)

!           ! increment time
!           ! params%time%t = params%time%t - params%time%dt ! this is not reboust in case of optimisation
!           params%time%t = ((s-1)*params%time%steps + nt) * params%time%dt

!           ! load direct solution and compute g
!           if (size(Q,5).gt.1) then
!              q0 = Q(:,:,:,:,nt)
!           else
!              q0 = Q(:,:,:,:,1)
!           end if

!           call calc_objective_g(g,q0)

!           call tic()

!           ! call time-steper
!           call time_stepper_adjoint(qs1,qsn,q0,g,params%time%t)
!           qsn = qs1

!           ! rhs toc
!           toc_rhs = toc()

!           ! boundary handling
!           call set_boundary_condition_qs(qsn,q0,outside_rhs=.true.)

!           ! filter
!           if (mod(nt,fi1freq).eq.0) then
!              call filter_i1(qsn)
!           end if

!           if (mod(nt,fi2freq).eq.0) then
!              call filter_i2(qsn)
!           end if

!           if (mod(nt,fi3freq).eq.0) then
!              call filter_i3(qsn)
!           end if

!           ! save to memory
!           Qs(:,:,:,:,nt) = qsn;

!           ! snapshot
!           if (mod(nt,params%io%asfreq).eq.0) then
!              call store(qsn,'data_adjoint_snapshot')
!           end if

!           ! snapshot with loop number
!           if (mod(nt,params%io%aslfreq).eq.0) then
!              call store(qsn,'data_adjoint_snapshot_loop')
!           end if

!           ! time_series
!           call sample_adj(qsn)
!           call probe_adj(qsn)

!           ! screen
!           if(params%parallelism%world_image.eq.1) then
!              if (mod(nt,params%io%sfreq).eq.0) then
!                 write(*,*) "computation (adjoint) of timestep ", trim(num2str(params%time%nt,'(I7.7)'))," in sub-set ", trim(num2str(params%time%ns,'(I7.7)')), " T/nt: ", &
!                      trim(num2str(toc(),'(F0.3)'))," (",trim(num2str(toc_rhs,'(F0.3)')),")"
!              end if
!           end if

!        end do

!        if (params%time%subsets.gt.1) then
!           if(params%parallelism%block_image.eq.1) then
!              if (params%io%verbosity.ge.1) then
!                 write(*,*) " store data adjoint of subset", s
!              end if
!           end if

!           ! store adjoint subset data
!           call store(Qs,'data_adjoint_subsets_cache')
!        end if

!     end do

!     ! clean
!     deallocate(qs1)
!     deallocate(qsn)
!     deallocate( q0)
!     deallocate(  g)

!   end subroutine calc_adjoint
! !!!=================================================================================================

! end module computation_adjoint
