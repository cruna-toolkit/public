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
  real(kind=rk) :: J
  integer       :: loop , maxloop, ls_start
  integer       :: afreq, dfreq, ffreq, pxy_lfreq, pxy_zidx
  logical       :: stop_after_adjoint, stop_after_direct

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

  do loop = 1, maxloop
     call set_parameter(loop,'opt.loop')

     ! DIRECT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(params%parallelism%world_image.eq.1) then
        write(*,*) ''
        write(*,*) 'Starting direct compuation of loop: ', trim(num2str(loop,'(I3.3)'))
        write(*,*) ''
     end if

     call init_direct(q0)
     call calc_direct(Q,q0)

     ! store full Q
     call get_parameter(dfreq,'io.dfreq',default=huge(1))
     if (mod(loop,dfreq).eq.0) then
        call store_state(Q,'data_direct_loop')
     end if

     ! store Q parts
     call get_parameter(pxy_lfreq,'io.pxy_lfreq',default=huge(1))
     if (mod(loop,pxy_lfreq).eq.0) then
        call get_parameter(pxy_zidx ,'io.pxy_zidx',default=1)
        call store(Q(:,:,pxy_zidx,5,:),'pxy_direct_loop')
     end if

     ! exit condition (direct)
     call get_parameter(stop_after_direct,'opt.stop_after_direct',default=.false.)
     if (stop_after_direct.eqv..true.) then
        if(params%parallelism%world_image.eq.1) then  
           write(*,*) 'stop_after_direct --> stop cruna'
        end if
        call stop_cruna
     end if

     ! objective
     call calc_objective(J,Q)
     if(params%parallelism%world_image.eq.1) then  
        call store_J        
        write(*,*) 'current relative objective function value (J): ', trim(num2str(objective_function(loop,1)/objective_function(1,1),'(F0.6)'))
     end if

     ! store forcing
     call get_parameter(ffreq,'io.ffreq',default=huge(1))
     if (mod(loop,ffreq).eq.0) then
        call store_force
     end if

     ! return conditions
     if (loop.eq.maxloop) then
        if(params%parallelism%world_image.eq.1) then  
           write(*,*) 'maximum number of loops reached --> stop cruna'
        end if
        call stop_cruna
     end if

     ! ADJOINT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(params%parallelism%world_image.eq.1) then
        write(*,*) ''
        write(*,*) 'Starting adjoint compuation of loop: ', trim(num2str(loop,'(I3.3)'))
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
        if(params%parallelism%world_image.eq.1) then  
           write(*,*) 'stop_after_adjoint --> stop cruna'
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
           write(*,*) 'starting line-search of loop: ', trim(num2str(loop,'(I3.3)'))
           write(*,*) ''
        end if

        call ls_quadratic_fit ! includes two direct computations

     end if

     ! UPDATE FORCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call update_force

  end do

!!! END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call end_parallelism

end program cruna
