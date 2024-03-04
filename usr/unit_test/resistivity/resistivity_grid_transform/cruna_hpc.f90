!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CRUNA - compressible reactive unsteady navier-stokes adjoint - code                               !
! ML (2017)                                                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! please: read README.txt                                                                           !
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
  real                                               :: toc_rhs
  integer                                            :: nt
  integer                                            :: s
  integer                                            :: fi1freq ,fi2freq ,fi3freq
  integer                                            :: fwi1freq,fwi2freq,fwi3freq
  real(kind=rk)                                      :: u, dp, chi, resistivity, length_porous
  real(kind=rk)                                      :: p_inlet, p_outlet ! data position (with offset)
  real(kind=rk)                                      :: p_start, p_end    ! actual position of porous material
  real(kind=rk), dimension(3)                        :: inlet_pos, outlet_pos
  integer, dimension(3), save                        :: inlet_idx  = (/0.0_rk,0.0_rk,0.0_rk/)
  integer, dimension(3), save                        :: outlet_idx = (/0.0_rk,0.0_rk,0.0_rk/)
  real(kind=rk), dimension(:,:,:), allocatable       :: distance
  integer                                            :: image_inlet  = 0.0_rk
  integer                                            :: image_outlet  = 0.0_rk
  real(kind=rk)                                      :: val    = 0.0_rk
  real(kind=rk)                                      :: x21, x31
  real(kind=rk)                                      :: pos_data_inlet, pos_data_outlet
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

!!! INIT & CHECK
  call init_direct(q0)
  call init_reference(q0)
  call init_filter_weighted(q0)

  ! store initial condition
  call set_parameter(1,'time.ns',struct_val=params%time%ns)
  call set_parameter(0,'time.nt',struct_val=params%time%nt)
  call store(q0,'data_direct_snapshot')
  qn = q0

  call check_cfl(q0)
  call check_c(q0)
  call check_Re

!!! GET FREQS
  call get_parameter(fi1freq,'filter.di1freq',default = huge(1))
  call get_parameter(fi2freq,'filter.di2freq',default = huge(1))
  call get_parameter(fi3freq,'filter.di3freq',default = huge(1))

  call get_parameter(fwi1freq,'filter.dwi1freq',default = huge(1))
  call get_parameter(fwi2freq,'filter.dwi2freq',default = huge(1))
  call get_parameter(fwi3freq,'filter.dwi3freq',default = huge(1))

!!! THE MAIN ACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! loop subsets
  do s = 1,params%time%subsets

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
        call set_boundary_condition_q(qn)

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

        ! snapshot
        if (mod(nt,params%io%dsfreq).eq.0) then
           call store(qn,'data_direct_snapshot')
        end if

        ! screen
        if(params%parallelism%world_image.eq.1) then
           if (mod(nt,params%io%sfreq).eq.0) then
              write(*,*) "computation (direct) of timestep ", trim(num2str(params%time%nt,'(I6.6)'))," in sub-set ", trim(num2str(params%time%ns,'(I6.6)')), " T/nt: ",trim(num2str(toc(),'(F00.3)'))," (",trim(num2str(toc_rhs,'(F0.3)')),")"
           end if
        end if

     end do

  end do

  !!! define inlet and outlet position (before and after the porous material)
 call get_parameter(x21,'geom.x21')
 call get_parameter(x31,'geom.x31')
 call get_parameter(pos_data_inlet, 'geom.pos_data_inlet')
 call get_parameter(pos_data_outlet,'geom.pos_data_outet')
 allocate(distance(params%geom%n1b,params%geom%n2b,params%geom%n3b))
 inlet_pos(1) = pos_data_inlet
 inlet_pos(2) = x21/2
 inlet_pos(3) = x31/2
 outlet_pos(1) = pos_data_outlet
 outlet_pos(2) = x21/2
 outlet_pos(3) = x31/2

 ! get length porous material
 call get_parameter(p_start, 'geom.p_start')
 call get_parameter(p_end ,'geom.p_end ')
 length_porous = p_end - p_start

 ! get chi / analytical resistivtiy
 call get_parameter(chi,'geom.darcy_amp')

 !!! get core that contains the chosen grid point
 call calc_spatial_distance(distance,X,inlet_pos)
 call get_index(inlet_idx, image_inlet, val,distance,'min')
 call calc_spatial_distance(distance,X,outlet_pos)
 call get_index(outlet_idx,image_outlet,val,distance,'min')

 ! get u
  u = 0.0_rk
 if (image_inlet.eq.params%parallelism%block_image) then
    ! compute u at chosen point (indices have been determined in advance for grid: 128x32x32)
        u  = qn(inlet_idx(1),inlet_idx(2),inlet_idx(3),2)/qn(inlet_idx(1),inlet_idx(2),inlet_idx(3),1)
 end if
 call bcast(u, image_inlet-1, params%parallelism%block_comm)

 ! get p inlet
 p_inlet= 0.0_rk
 if (image_inlet.eq.params%parallelism%block_image) then
    ! compute u at chosen point (indices have been determined in advance for grid: 128x32x32)
        p_inlet  = qn(inlet_idx(1),inlet_idx(2),inlet_idx(3),5)
 end if
 call bcast(p_inlet, image_inlet-1, params%parallelism%block_comm)


 ! get p outlet
 p_outlet= 0.0_rk
 if (image_outlet.eq.params%parallelism%block_image) then
    ! compute u at chosen point (indices have been determined in advance for grid: 128x32x32)
        p_outlet  = qn(outlet_idx(1),outlet_idx(2),outlet_idx(3),5)
 end if
 call bcast(p_outlet, image_outlet-1, params%parallelism%block_comm)


 ! compute dp
  dp = p_inlet - p_outlet

 !  compute resisitvity
 resistivity = dp/length_porous/u

 if (params%parallelism%world_image.eq.1) then
    write(*,*) 'computed resisitvity : ',resistivity
    write(*,*) 'given resisitvity : ',chi
    write(*,*) 'relative difference (in %): ',(resistivity-chi)/chi *100
    if ((resistivity-chi)/chi  < 0.01) then
        write(*,*) 'Unit Test Resisitvity PASSED'
    else
        write(*,*) 'Unit Test Resisitvity FAILED'
        stop
    end if
 end if


  ! clean
  deallocate(q1)

!!! END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call end_parallelism

end program cruna_hpc
