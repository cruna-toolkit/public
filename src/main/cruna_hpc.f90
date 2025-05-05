!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CRUNA_HPC - compressible reactive unsteady navier-stokes hpc - code                               !
! ML (202+)                                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! please: read README.md                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cruna_hpc
  ! general
  use boundary_conditions
  use check_state
  use data_geom
  use data_reference
  use data_rhs
  use discretisation_t
  use display
  use helper
  use filter_x
  use filter_weighted_x
  use io
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
  real(kind=rk), dimension(:,:,:,:)  ,allocatable    :: qn,q1
  real(kind=rk)                                      :: dq_max_loc,dq_max_glo,dq_max_tol
  real(kind=rk)                                      :: cputime_glo,walltime,cputime_loc
  real                                               :: toc_rhs,cputime_loc_r
  integer      , dimension(2)                        :: ns_nt_start = (/1,0/)
  integer                                            :: nt,nt0,nt_total
  integer                                            :: s,s0
  integer                                            :: cfreq
  integer                                            :: fi1freq ,fi2freq ,fi3freq
  integer                                            :: fwi1freq,fwi2freq,fwi3freq
  character(len=max_length_fname)                    :: restartfile
  logical                                            :: restart, file_exisit

!!! START !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call init_parallelism                   ! set image, block, communicators, spmd_type
  call display_cruna_hpc_start            ! display logo and general information

!!! PARAMETER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call init_parameter_array               ! read parameter.dat
  call spread_parameter                   ! spread parameter corresponding to parallelism
  call init_parameter_struct              ! fill params% struct

!!! TOPOLOGY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call init_topology                      ! initializes mpi topologies

!!! ALLOCATES & GEOMETRY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call allocate_geometry                  ! allocates the memory in data_geom
  call init_geometry                      ! init the geometry

  call allocate_force                     ! allocates  force f (optimization stuff)
  call init_force                         ! initialize force f (optimization stuff)

  call allocate_volume_penalization       ! allocates volume penalization
  call init_volume_penalization           ! initialize volume penalization

  allocate(q0(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
  allocate(q1(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))
  allocate(qn(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))

!!! BOUNDARY CONDITION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call init_boundary_conditions_array     ! read boundary condition file
  call spread_boundary_conditions_array   ! spread boundary conditions corresponding to parallelism
  call modify_boundary_conditions_array   ! modify array to local indices and clean up unn. conditions

!!! GET FREQS
  call get_parameter(cfreq     ,'run.cfreq'   ,default = huge(1))
  call get_parameter(dq_max_tol,'run.ctol'    ,default = tiny(1.0_rk))

  call get_parameter(walltime  ,'run.walltime',default = huge(1.0_rk)) ! walltime in hours
  walltime = walltime*3600_rk

  call get_parameter(fi1freq,'filter.di1freq',default = huge(1))
  call get_parameter(fi2freq,'filter.di2freq',default = huge(1))
  call get_parameter(fi3freq,'filter.di3freq',default = huge(1))

  call get_parameter(fwi1freq,'filter.dwi1freq',default = huge(1))
  call get_parameter(fwi2freq,'filter.dwi2freq',default = huge(1))
  call get_parameter(fwi3freq,'filter.dwi3freq',default = huge(1))

!!! INIT & CHECK
  call init_direct(q0)
  call init_reference(q0)
  call init_filter_weighted(q0)

  ! store initial condition
  call set_parameter(1,'time.ns',struct_val=params%time%ns)
  call set_parameter(0,'time.nt',struct_val=params%time%nt)
  call store(q0,'data_direct_snapshot')

  call check_cfl(q0)
  call check_c(q0)
  call check_Re

  qn = q0 ! normal start, maybe overwritten by restart procedure

!!! RESTART PROCEDURE
  call get_parameter(restart,'init.restart',default = .false.)
  if(restart.eqv..true.) then

     if(params%parallelism%block_image.eq.1) then

        call get_fname(restartfile,'restartfile')
        restartfile = trim(restartfile) // ".h5"
        call inquire_file(restartfile,file_exisit)

        if(file_exisit.eqv..true.) then
           call load(ns_nt_start, 'restartfile')

           write(*,*)
           write(*,*) "WARNING restart file <",trim(restartfile),"> loading ns: ",trim(num2str(ns_nt_start(1),'(I7.7)'))," nt: ",trim(num2str(ns_nt_start(2),'(I7.7)'))
           write(*,*)
        else
           write(*,*)
           write(*,*) "WARNING restart file <",trim(restartfile),"> not found: start from 0"
           write(*,*)
        end if
     end if

     ! distribute ns_nt_start
     call allreduce(ns_nt_start(1),ns_nt_start(1),'max',params%parallelism%world_comm)
     call allreduce(ns_nt_start(2),ns_nt_start(2),'max',params%parallelism%world_comm)

     call load(qn,'data_direct_snapshot', ns_optin = ns_nt_start(1), nt_optin = ns_nt_start(2))

  end if

  ! get start subset s0 und start time-step nt0:
  nt_total = (ns_nt_start(1)-1)*params%time%steps + ns_nt_start(2) + 1                                   ! computing overall time-step based on subset und time-step, add 1 (to avoid recomputing), nt_total: nt with subsets = 1
  s0       = nt_total/params%time%steps + 1                                                              ! gives s0 (using integer arithmetic)
  nt0      = nt_total - (s0-1)*params%time%steps

  call set_parameter(s0 ,'init.s0' )
  call set_parameter(nt0,'init.nt0')

!!! THE MAIN ACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! loop subsets
  do s = s0,params%time%subsets

     !params%time%ns = s
     call set_parameter(s,'time.ns',struct_val=params%time%ns)

     ! load subset data force
     if(params%time%subsets.gt.1) then
        if(params%parallelism%block_image.eq.1) then
           if (params%io%verbosity.ge.1) then
              write(*,*) " load force of subset", s
           end if
        end if

        ! load actual force (actual subset is set above)
        call load(f,'force_subsets_cache')
     end if

     ! loop timesteps
     do nt = nt0,params%time%steps
        !params%time%nt = nt
        call set_parameter(nt,'time.nt',struct_val=params%time%nt)

        ! increment time
        !! params%time%t = params%time%t + params%time%dt ! this is not robust in case of optimisation
        params%time%t = ((s-1)*params%time%steps + nt) * params%time%dt

        call tic()

        ! call time-steper
        call time_stepper_direct(q1,qn,params%time%t)

        ! convergence check
        if (mod(nt,cfreq).eq.0) then
           dq_max_loc = real(maxval(abs((qn-q1)/params%time%dt)),rk)

           call allreduce(dq_max_glo,dq_max_loc,'max',params%parallelism%world_comm)

           if(params%parallelism%world_image.eq.1) then
              write(*,*) 'residual:',trim(num2str(dq_max_glo,'(E13.8)')) ,' (',trim(num2str(dq_max_tol,'(E13.8)')),')'
           end if

           if(dq_max_glo.le.dq_max_tol) then
              call store(q1,'data_direct_snapshot')
              if(params%parallelism%block_image.eq.1) then
                 ns_nt_start(1) = s
                 ns_nt_start(2) = nt
                 call store(ns_nt_start, 'restartfile', file_overwrite_optin=.true.)
              end if

              if(params%parallelism%world_image.eq.1) then
                 write(*,*) 'residual tolerance reached:',trim(num2str(dq_max_glo,'(E13.8)')),' <= ',trim(num2str(dq_max_tol,'(E13.8)'))
                 write(*,*) "last timestep ", trim(num2str(params%time%nt,'(I7.7)'))," in sub-set ", trim(num2str(params%time%ns,'(I7.7)'))
              end if
              call stop_cruna
           end if
        end if

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

        ! filter_weighted
        if (mod(nt,fwi1freq).eq.0) then
           call filter_weighted_i1(qn)
        end if

        if (mod(nt,fwi2freq).eq.0) then
           call filter_weighted_i2(qn)
        end if

        if (mod(nt,fwi3freq).eq.0) then
           call filter_weighted_i3(qn)
        end if

        ! time_series
        call sample(qn)
        call probe(qn)

        ! wall time
        call cpu_time(cputime_loc_r)
        cputime_loc = real(cputime_loc_r,rk)                                            ! transfer from real (cpu_time intrinsic) to real_rk (cruna)
        call allreduce(cputime_glo,cputime_loc,'max',params%parallelism%world_comm)     ! allreduce not defined for real, just real_rk

        if (cputime_glo.gt.walltime) then

           call store(qn,'data_direct_snapshot')
           if(params%parallelism%block_image.eq.1) then
              ns_nt_start(1) = s
              ns_nt_start(2) = nt
              call store(ns_nt_start, 'restartfile', file_overwrite_optin=.true.)
           end if

           if(params%parallelism%world_image.eq.1) then
              write(*,*) 'cruna walltime reached: ',trim(num2str(cputime_glo/3600_rk,'(F10.7)'))," (",trim(num2str(walltime/3600_rk,'(F10.7)')) ,") hrs"
              write(*,*) "last timestep ", trim(num2str(params%time%nt,'(I7.7)'))," in sub-set ", trim(num2str(params%time%ns,'(I7.7)'))
           end if

           call stop_cruna
        end if

        ! snapshot
        if (mod(nt,params%io%dsfreq).eq.0) then
           call store(qn,'data_direct_snapshot')
           if(params%parallelism%block_image.eq.1) then
              ns_nt_start(1) = s
              ns_nt_start(2) = nt
              call store(ns_nt_start, 'restartfile', file_overwrite_optin=.true.)
           end if
        end if

        ! screen
        if(params%parallelism%world_image.eq.1) then
           if (mod(nt,params%io%sfreq).eq.0) then
              write(*,*) "computation (direct) of timestep ", trim(num2str(params%time%nt,'(I7.7)'))," in sub-set ", trim(num2str(params%time%ns,'(I7.7)')), " T/nt: ",trim(num2str(toc(),'(F00.3)'))," (",trim(num2str(toc_rhs,'(F0.3)')),")"
           end if
        end if

     end do

  end do

  ! store last timestep as snapshot
  call store(qn,'data_direct_snapshot')

  ! clean
  deallocate(q1)

  if(params%parallelism%world_image.eq.1) then
     write(*,*) 'computation finished --> stop cruna_hpc'
  end if

!!! END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call end_parallelism

end program cruna_hpc
