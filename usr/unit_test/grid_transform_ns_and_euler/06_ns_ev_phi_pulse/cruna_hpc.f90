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
  use discretisation_x
  use display
  use evaluating_unit_test
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
  
  ! unit test variables
  real(kind=rk), dimension(:,:,:), allocatable    :: tmp,I, tmp_I, Ident 
  real(kind=rk)                                   :: I1, I2, I3 

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

        ! time_series
        call sample(qn)
        call probe(qn)

        ! screen
        if(params%parallelism%world_image.eq.1) then
           if (mod(nt,params%io%sfreq).eq.0) then
              write(*,*) "computation (direct) of timestep ", trim(num2str(params%time%nt,'(I6.6)'))," in sub-set ", trim(num2str(params%time%ns,'(I6.6)')), " T/nt: ",trim(num2str(toc(),'(F00.3)'))," (",trim(num2str(toc_rhs,'(F0.3)')),")"
           end if
        end if

     end do     

  end do

  ! clean
  deallocate(q1)
  write(*,*) "run simulation"
  call snapshot_comparison()
!!! END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call end_parallelism
  
    !!
  

end program cruna_hpc
