!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CRUNA - compressible reactive unsteady navier-stokes adjoint - code                               !
! ML (2017)                                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! please: read README.md                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cruna
  ! general
  use boundary_conditions
  use check_state
  use computation_direct
  use computation_adjoint
  use data_geom
  use data_reference
  use data_rhs
  use display
  use helper
  use io_helper
  use io_wrapper
  use line_search
  use parameter
  use volume_penalization

  ! case dependent
  use force
  use geometry
  use initial_condition
  use objective
  use parallelism

  ! main routine vars
  real(kind=rk)      , dimension(2)                  :: loop_start = (/1.0_rk,1.0_rk/)
  real(kind=rk)                                      :: J
  integer                                            :: loop,starting_loop,maxloop, ls_start,loop_exit_freq
  integer                                            :: afreq, dfreq, ffreq !, pxy_lfreq, pxy_zidx
  integer                                            :: stop_after_direct_loop,stop_after_adjoint_loop
  logical                                            :: stop_after_adjoint, stop_after_direct
  logical                                            :: restart_loop, file_exist
  character(len=max_length_fname)                    :: restartfile_loop
!!! START !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call init_parallelism                   ! set image, block, communicators, spmd_type
  call display_cruna_start                ! display logo and general information

!!! PARAMETER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call init_parameter_array               ! read parameter.dat
  call spread_parameter                   ! spread parameter corresponding to parallelism
  call init_parameter_struct              ! fill params% struct

!!! TOPOLOGY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call init_topology                      ! initializes mpi topologies

!!! ALLOCATES & GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call allocate_geometry                  ! allocates the memory in data_geom
  call init_geometry                      ! init the geometry

  call allocate_data_direct               ! allocates the memory in data_rhs
  call allocate_data_adjoint              ! allocates the memory in data_rhs

  call allocate_force                     ! allocates  force f (optimization stuff)
  call init_force                         ! initialize force f (optimization stuff)

  call allocate_volume_penalization       ! allocates volume penalization
  call init_volume_penalization           ! initialize volume penalization

!!! BOUNDARY CONDITION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call init_boundary_conditions_array     ! read boundary condition file
  call spread_boundary_conditions_array   ! spread boundary conditions corresponding to parallelism
  call modify_boundary_conditions_array   ! modify array to local indices and clean up unn. conditions

!!! INIT CHECK
  call init_direct(q0)
  call init_reference(q0)
  call check_cfl(q0)
  call check_c(q0)
  call check_Re

!!! THE MAIN ACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call get_parameter(maxloop,'opt.maxloop',default = 1)

!!! RESTART PROCEDURE (only for completed loops)

  !===============================================================================
  ! NOTE:
  ! The restart procedure implemented here requires additional steps from the user
  ! to function correctly. Since the volume penalization, objective, and force
  ! terms are user-defined, users must ensure compatibility with their specific
  ! implementations. In particular, modifications may be required in the following
  ! files:
  !   - force.f90
  !   - volume_penalization.f90
  !   - objective.f90
  !
  ! Please review and adapt these components as necessary to ensure proper
  ! operation of the restart functionality for your case.
  !===============================================================================

  call get_parameter(restart_loop,'init.restart_loop',default = .false.)
  if(restart_loop.eqv..true.) then

     if(params%parallelism%block_image.eq.1) then
        call get_fname(restartfile_loop,'restartfile_loop')
        restartfile_loop = trim(restartfile_loop) // ".h5"
        call inquire_file(restartfile_loop,file_exist)

        if(file_exist.eqv..true.) then
           call load(loop_start, 'restartfile_loop')
           call get_parameter(loop_start(2) ,'opt.alpha0') ! safety check, if user as defined opt.alpha0 in parameter.dat
           call set_parameter(loop_start(2) ,'opt.alpha0') ! set alpha0
           write(*,*)
           write(*,*) "WARNING restart loop file <",trim(restartfile_loop),"> start from loop: ",trim(num2str(int(loop_start(1)),'(I4.4)'))
           write(*,*) "WARNING restart loop file <",trim(restartfile_loop),"> set alpha0 to: "  ,loop_start(2)
           write(*,*)
        else
           write(*,*)
           write(*,*) "WARNING restart loop file <",trim(restartfile_loop),"> not found: start from loop: 1"
           write(*,*)

        end if
     end if

     ! distribute loop_start
     call allreduce(loop_start(1),loop_start(1),'max',params%parallelism%world_comm)
     call allreduce(loop_start(2),loop_start(2),'max',params%parallelism%world_comm)
  end if
  starting_loop = int(loop_start(1))

  ! --- Only start simulation simulation is not completed
  if (starting_loop <= maxloop) then
     ! Proceed with main simulation
     ! start main loop
     do loop = starting_loop, maxloop
        call set_parameter(loop,'opt.loop')

        ! DIRECT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(params%parallelism%world_image.eq.1) then
           write(*,*) ''
           write(*,*) 'Starting direct compuation of loop: ', trim(num2str(loop,'(I4.4)'))
           write(*,*) ''
        end if

        call init_direct(q0)
        call calc_direct(Q,q0)

        ! store full Q
        call get_parameter(dfreq,'io.dfreq',default=huge(1))
        if (mod(loop,dfreq).eq.0) then
           call store_state(Q,'data_direct_loop')
        end if

        ! exit condition (direct)
        call get_parameter(stop_after_direct,'opt.stop_after_direct_pre_J',default=.false.)
        if (stop_after_direct.eqv..true.) then
           write(*,*) 'stop_after_direct --> stop cruna (loop ',trim(num2str(loop,'(I4.4)')),')'
           call stop_cruna
        end if

        ! objective
        call calc_objective(J,Q)
        if(params%parallelism%world_image.eq.1) then
           call store_J
           write(*,*) 'current relative objective function value (J): ', trim(num2str(objective_function(loop,1)/objective_function(1,1),'(F0.6)'))
        end if

        ! exit condition (direct)
        call get_parameter(stop_after_direct,'opt.stop_after_direct',default=.false.)
        if (stop_after_direct.eqv..true.) then
           call get_parameter(stop_after_direct_loop,'opt.stop_after_direct',2,default=loop)
           if (loop.eq.stop_after_direct_loop) then
              if(params%parallelism%world_image.eq.1) then
                 write(*,*) 'stop_after_direct --> stop cruna (loop ',trim(num2str(loop,'(I4.4)')),')'
              end if
              call stop_cruna
           end if
        end if

        ! store forcing
        call get_parameter(ffreq,'io.ffreq',default=huge(1))
        if (mod(loop,ffreq).eq.0) then
           call store_force
        end if

        ! ADJOINT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(params%parallelism%world_image.eq.1) then
           write(*,*) ''
           write(*,*) 'Starting adjoint compuation of loop: ', trim(num2str(loop,'(I4.4)'))
           write(*,*) ''
        end if

        call init_adjoint(qs0)
        call calc_adjoint(Qs,qs0,Q)

        ! store full Qs
        call get_parameter(afreq,'io.afreq',default=1)
        if ((mod(loop,afreq).eq.0).or.(loop.eq.1)) then
           call store_state(Qs,'data_adjoint_loop')
        end if

        ! exit condition (adjoint)
        call get_parameter(stop_after_adjoint,'opt.stop_after_adjoint',default = .false.)
        if (stop_after_adjoint.eqv..true.) then
           call get_parameter(stop_after_adjoint_loop,'opt.stop_after_adjoint',2,default = loop)
           if (loop.eq.stop_after_adjoint_loop) then
              if(params%parallelism%world_image.eq.1) then
                 write(*,*) 'stop_after_adjoint --> stop cruna (loop ',trim(num2str(loop,'(I4.4)')),')'
              end if
           end if
           call stop_cruna
        end if

        ! GRADIENT/DIRECTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call calc_grad(Qs,Q) ! includes store grad

        ! LINE-SEARCH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call get_parameter(ls_start,'opt.ls_start',default = huge(1))

        if (loop.ge.ls_start) then
           if(params%parallelism%world_image.eq.1) then
              write(*,*) ''
              write(*,*) 'starting line-search of loop: ', trim(num2str(loop,'(I4.4)'))
              write(*,*) ''
           end if

           call ls_quadratic_fit ! includes two direct computations

        end if

        ! UPDATE FORCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call update_force


        ! UPDATE restart loop (if enabled)
        if(restart_loop.eqv..true.) then
           if(params%parallelism%world_image.eq.1) then
              write(*,*) ''
              write(*,*) 'udpate loop restart counter to: ', trim(num2str(loop+1,'(I4.4)'))
              write(*,*) ''
              call get_parameter(loop_start(2),'opt.alpha0') ! load current alpha0
              loop_start(1) = real(loop+1,rk)                ! overwrite with last loop
              call store(loop_start, 'restartfile_loop', file_overwrite_optin=.true.)
           end if
        end if

        !-----------------------------------------------------------------------
        ! Stop CRUNA every 'loop_exit_freq' iterations (disabled if 0),
        ! but NEVER at the (re)start loop.
        !-----------------------------------------------------------------------
        call get_parameter(loop_exit_freq, 'run.loop_exit_freq', default=0)
        if (loop_exit_freq > 0) then
           if (mod(loop, loop_exit_freq) .eq. 0) then
              if (params%parallelism%world_image .eq. 1) then
                 write(*,*) 'Loop exit frequency reached (loop = ', trim(num2str(loop,'(I4.4)')), ') --> stopping CRUNA.'
              end if
              call stop_cruna
           end if
        end if
     end do
  else
     ! simulation already finished
     if (params%parallelism%block_image.eq.1) write(*,*) 'Maximum number of loops reached, simulation is completed --> stopping CRUNA.'
     call stop_cruna
  end if

!!! END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call end_parallelism

end program cruna
