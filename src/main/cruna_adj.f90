!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CRUNA_ADJ - compressible reactive unsteady navier-stokes adjoint - code                           !
! ML (202+)                                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! please: read README.md                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cruna_adj
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
  use time_series
  use volume_penalization

  ! case dependent
  use force
  use geometry
  use initial_condition
  use objective
  use parallelism

  ! variables
  real(kind=rk), dimension(:,:,:,:)  ,allocatable    :: qsn,qs1,g
  real                                               :: toc_rhs
  integer      , dimension(2)                        :: ns_nt_start
  integer                                            :: nt,nt0,nt_total
  integer                                            :: s,s0
  integer                                            :: fi1freq ,fi2freq ,fi3freq
  character(len=max_length_fname)                    :: restartfile
  logical                                            :: restart, file_exisit

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

  call allocate_volume_penalization       ! allocates volume penalization
  call init_volume_penalization           ! initialize volume penalization

  allocate( q0(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
  allocate(qs0(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
  allocate(qs1(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
  allocate(qsn(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
  allocate(  g(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))

!!! BOUNDARY CONDITION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call init_boundary_conditions_array     ! read boundary condition file
  call spread_boundary_conditions_array   ! spread boundary conditions corresponding to parallelism
  call modify_boundary_conditions_array   ! modify array to local indices and clean up unn. conditions

!!! INIT CHECK
  call init_direct(q0)
  call init_reference(q0)
  call init_adjoint(qs0)

  call check_cfl(q0)
  call check_c(q0)
  call check_Re

  ! get filter freqs
  call get_parameter(fi1freq,'filter.ai1freq',default = huge(1))
  call get_parameter(fi2freq,'filter.ai2freq',default = huge(1))
  call get_parameter(fi3freq,'filter.ai3freq',default = huge(1))

  ! store initial condition
  call set_parameter(params%time%steps + 1,'time.nt',struct_val=params%time%nt)
  call store(qs0,'data_adjoint_snapshot')

  qsn = qs0 ! normal start, maybe overwritten by restart procedure

!!! RESTART PROCEDURE
  call get_parameter(ns_nt_start(1),'time.subsets', default = 1)
  call get_parameter(ns_nt_start(2),'time.steps'               )

  call get_parameter(restart,'init.restart',default = .false.)
  if(restart.eqv..true.) then

     if(params%parallelism%block_image.eq.1) then

        call get_fname(restartfile,'restartfile')
        restartfile = trim(restartfile) // ".h5"
        call inquire_file(restartfile,file_exisit)

        if(file_exisit.eqv..true.) then
           call load(ns_nt_start, 'restartfile')

           write(*,*)
           write(*,*) "WARNING restart file <",trim(restartfile),"> loading ns:",ns_nt_start(1)," nt:",ns_nt_start(2)
           write(*,*)
        else
           write(*,*)
           write(*,*) "WARNING restart file <",trim(restartfile),"> not found: start from 0"
           write(*,*)
        end if
     end if

     ! distribute ns_nt_start
     call allreduce(ns_nt_start(1),ns_nt_start(1),'min',params%parallelism%world_comm)
     call allreduce(ns_nt_start(2),ns_nt_start(2),'min',params%parallelism%world_comm)

     call load(qsn,'data_adjoint_snapshot', ns_optin = ns_nt_start(1), nt_optin = ns_nt_start(2))

  end if

  ! get start subset s0 und start time-step nt0:
  nt_total = (ns_nt_start(1)-1)*params%time%steps + ns_nt_start(2) - 1                                   ! computing overall time-step based on subset und time-step, minus 1 (to avoid recomputing), nt_total: nt with subsets = 1
  s0       = nt_total/params%time%steps + 1                                                              ! gives s0 (using integer arithmetic)
  nt0      = nt_total - (s0-1)*params%time%steps

  call set_parameter(s0 ,'init.s0' )
  call set_parameter(nt0,'init.nt0')

!!! THE MAIN ACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! loop subsets
  do s = s0,1,-1
     call set_parameter(s,'time.ns',struct_val=params%time%ns)

     ! loop timesteps
     do nt = nt0,1,-1
        call set_parameter(nt,'time.nt',struct_val=params%time%nt)

        ! increment time
        ! params%time%t = params%time%t - params%time%dt ! this is not reboust
        params%time%t = ((s-1)*params%time%steps + nt) * params%time%dt

        ! load direct solution and compute g
        call calc_objective_g(g,q0)

        call tic()

        ! call time-steper
        call time_stepper_adjoint(qs1,qsn,q0,g,params%time%t)
        qsn = qs1

        ! rhs toc
        toc_rhs = toc()

        ! boundary handling
        call set_boundary_condition_qs(qsn,q0,outside_rhs=.true.)

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

        ! snapshot
        if (mod(nt,params%io%asfreq).eq.0) then
           call store(qsn,'data_adjoint_snapshot')
        end if

        ! time_series
        call sample_adj(qsn)
        call probe_adj(qsn)

        ! screen
        if(params%parallelism%world_image.eq.1) then
           if (mod(nt,params%io%sfreq).eq.0) then
              write(*,*) "computation (adjoint) of timestep ", trim(num2str(params%time%nt,'(I7.7)'))," in sub-set ", trim(num2str(params%time%ns,'(I7.7)')), " T/nt: ", &
                   trim(num2str(toc(),'(F0.3)'))," (",trim(num2str(toc_rhs,'(F0.3)')),")"
           end if
        end if

     end do

  end do

  if(params%parallelism%world_image.eq.1) then
     write(*,*) '(adjoint) computation finished --> stop cruna_adj'
  end if

!!! END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call end_parallelism

end program cruna_adj
