module objective_adjoint_sound_monopole_array

  use data_geom
  use force
  use helper
  use io
  use parallelism
  use parameter
  use tub_ak_data

  private

  real(kind=rk),dimension(:,:,:,:),allocatable,save, public :: opt_target
  real(kind=rk),dimension(:,:,:)  ,allocatable,save         :: sigma_x
  real(kind=rk),dimension(:,:)    ,allocatable,save, public :: objective_function
  real(kind=rk),dimension(:,:)    ,allocatable,save, public :: sigma_t

  public :: calc_objective
  public :: calc_objective_g

contains

!!!=================================================================================================
  subroutine calc_objective(J,Q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                     ,intent(out)                :: J
    real(kind=rk),dimension(:,:,:,:,:),intent(inout)              :: Q      ! (inout) is required for mutiple subsets
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                                 :: J_loc
    real(kind=rk)                                                 :: p_init
    integer                                                       :: loop,maxloop  
    integer                                                       :: nt,ns,ns_backup
    integer                                                       :: i1,i2,i3, in_line_search
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! get p_init for offset
    call get_parameter(p_init,'init.p')

    ! allocate/init opt_target
    if (.not.allocated(opt_target)) then
       allocate(opt_target(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%time%steps))

       if (params%time%subsets.eq.1) then
          call load(opt_target,'opt_target')
          opt_target = opt_target + p_init
       end if
    end if

    ! allocate/init sigma_x
    if (.not.allocated(sigma_x)) then
       call init_sigma_x
    end if

    ! allocate/init sigma_t
    if (.not.allocated(sigma_t)) then
       call init_sigma_t
    end if

    ! allocate objective function
    if (.not.allocated(objective_function)) then
       call get_parameter(maxloop,'opt.maxloop',default = 1)
       allocate(objective_function(maxloop,1))
       objective_function = -1.0_rk
    end if

    J_loc     = 0.0_rk

    ns_backup = params%time%ns
    do ns = 1,params%time%subsets

       if (params%time%subsets.gt.1) then
          if(params%parallelism%block_image.eq.1) then
             if (params%io%verbosity.ge.1) then
                write(*,*) " read data direct and opt_target of subset", ns
             end if
          end if

          call set_parameter(ns,'time.ns',struct_val=params%time%ns)
          call load(Q          ,'data_direct_subsets_cache')

          call load(opt_target ,'opt_target'               )
          opt_target = opt_target + p_init          
       end if

       do nt = 1,params%time%steps
          do i3 = 1,params%geom%n3b
             do i2 = 1,params%geom%n2b
                do i1 = 1,params%geom%n1b
                   J_loc = J_loc + ((Q(i1,i2,i3,5,nt)-opt_target(i1,i2,i3,nt))**2)*sigma_x(i1,i2,i3)*sigma_t(nt,ns)
                end do
             end do
          end do
       end do

    end do
    call set_parameter(ns_backup,'time.ns',struct_val=params%time%ns)

    J_loc = J_loc*params%time%dt*params%geom%dx1*params%geom%dx2*params%geom%dx3

    if (params%io%verbosity.ge.1) then
       write(*,*) " J (block):",J_loc,params%parallelism%block_image
    end if

    call allreduce(J, J_loc, 'sum', params%parallelism%block_comm)

    call get_parameter(in_line_search,'opt.in_line_search',default = 0)
    if (in_line_search.eq.0) then
       call get_parameter(loop,'opt.loop',default = 1)
       objective_function(loop,1) = J
    end if

  end subroutine calc_objective
!!!=================================================================================================

!!!=================================================================================================
  subroutine calc_objective_g(g,q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:),intent(out)                :: g
    real(kind=rk),dimension(:,:,:,:),intent(in)                 :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                               :: p_init
    integer,save                                                :: ns_persistent = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle subsets by load of current opt_target
    if (params%time%subsets.gt.1) then
       if (params%time%ns.ne.ns_persistent) then
          if(params%parallelism%block_image.eq.1) then
             if (params%io%verbosity.ge.1) then
                write(*,*) " read opt_target of subset", params%time%ns," for us in calc_objective_g"
             end if
          end if

          ns_persistent = params%time%ns
          call load(opt_target ,'opt_target')

          ! get p_init for offset
          call get_parameter(p_init,'init.p')
          opt_target = opt_target + p_init
       end if
    end if

    !(A.transpose()).inv()
    ![1, -u_1/rho, -u_2/rho, -u_3/rho,         0]
    ![0,    1/rho,        0,        0,         0]
    ![0,        0,    1/rho,        0,         0]
    ![0,        0,        0,    1/rho,         0]
    ![0,        0,        0,        0, gamma - 1]

    g          = 0.0_rk
    g(:,:,:,5) = -(params%material%gamma - 1.0_rk)*(q(:,:,:,5) - opt_target(:,:,:,params%time%nt))*sigma_x*sigma_t(params%time%nt,params%time%ns)

  end subroutine calc_objective_g
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_sigma_x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), allocatable                             :: r
    real(kind=rk)                                                            :: x10,x11,x20,x21,x30,x31,len,rad,rad_inner,rad_outer
    real(kind=rk), dimension(3)                                              :: default_origin,origin
    real(kind=rk)                                                            :: pi
    character(len=max_length_parameter)                                      :: sigma_x_type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(sigma_x(params%geom%n1b,params%geom%n2b,params%geom%n3b))
    sigma_x = 1.0_rk

    call get_parameter(sigma_x_type,'opt.sigma_x.type',default='rectangular')

    select case (trim(sigma_x_type))

    case('rectangular')

       call get_parameter(x10,'opt.sigma_x.x10')
       call get_parameter(x11,'opt.sigma_x.x11')
       call get_parameter(x20,'opt.sigma_x.x20')
       call get_parameter(x21,'opt.sigma_x.x21')
       call get_parameter(x30,'opt.sigma_x.x30')
       call get_parameter(x31,'opt.sigma_x.x31')
       call get_parameter(len,'opt.sigma_x.len') ! flank width

       pi = 4.0_rk*atan(1.0_rk)
       ! df/dx = slope = Dx[0.5*erf(sqrt(pi)*x/len)]= 1/len -> since df=1 it follows dx = len

       sigma_x = 0.5_rk*(erf(sqrt(pi)/len*(X(:,:,:,1) - x10)) - erf(sqrt(pi)/len*(X(:,:,:,1) - x11))) &
            * 0.5_rk*(erf(sqrt(pi)/len*(X(:,:,:,2) - x20)) - erf(sqrt(pi)/len*(X(:,:,:,2) - x21))) &
            * 0.5_rk*(erf(sqrt(pi)/len*(X(:,:,:,3) - x30)) - erf(sqrt(pi)/len*(X(:,:,:,3) - x31)))

    case ('circular','sphere')
       ! full sphere

       allocate(r(params%geom%n1b,params%geom%n2b,params%geom%n3b))

       call get_parameter(x10,'opt.sigma_x.x10')
       call get_parameter(x20,'opt.sigma_x.x20')
       call get_parameter(x30,'opt.sigma_x.x30')
       call get_parameter(len,'opt.sigma_x.len') ! flank width
       call get_parameter(rad,'opt.sigma_x.rad')

       pi = 4.0_rk*atan(1.0_rk)
       ! df/dx = slope = Dx[0.5*erf(sqrt(pi)*x/len)]= 1/len -> since df=1 it follows dx = len

       r = sqrt((X(:,:,:,1)-x10)**2 + (X(:,:,:,2)-x20)**2 + (X(:,:,:,3)-x30)**2)

       sigma_x = 0.5_rk - 0.5_rk*erf(sqrt(pi)/len*(r - rad))

    case ('circular_2','shell')
       ! spherical shell

       allocate(r(params%geom%n1b,params%geom%n2b,params%geom%n3b))

       call get_parameter(x10,'opt.sigma_x.x10') ! shell origin
       call get_parameter(x20,'opt.sigma_x.x20')
       call get_parameter(x30,'opt.sigma_x.x30')
       call get_parameter(len,'opt.sigma_x.len') ! flank width
       call get_parameter(rad_inner,'opt.sigma_x.rad_inner')
       call get_parameter(rad_outer,'opt.sigma_x.rad_outer')

       pi = 4.0_rk*atan(1.0_rk)
       ! df/dx = slope = Dx[0.5*erf(sqrt(pi)*x/len)]= 1/len -> since df=1 it follows dx = len 

       r = sqrt((X(:,:,:,1)-x10)**2 + (X(:,:,:,2)-x20)**2 + (X(:,:,:,3)-x30)**2)

       sigma_x = - 0.5_rk*erf(sqrt(pi)/len*(r - rad_outer)) + 0.5_rk*erf(sqrt(pi)/len*(r - rad_inner))

    case('swiss_cheese')
       ! box with spherical hole

       allocate(r(params%geom%n1b,params%geom%n2b,params%geom%n3b))

       call get_parameter(x10,'opt.sigma_x.x10') ! outer box extents
       call get_parameter(x11,'opt.sigma_x.x11')
       call get_parameter(x20,'opt.sigma_x.x20')
       call get_parameter(x21,'opt.sigma_x.x21')
       call get_parameter(x30,'opt.sigma_x.x30')
       call get_parameter(x31,'opt.sigma_x.x31')

       call get_parameter(len,'opt.sigma_x.len') ! error function flank width

       call get_parameter(rad,'opt.sigma_x.rad') ! inner hole radius

       ! default hole origin is center of outer box
       default_origin(1) = 0.5_rk*(x10+x11)
       default_origin(2) = 0.5_rk*(x20+x21)
       default_origin(3) = 0.5_rk*(x30+x31)
       call get_parameter(origin(1),'opt.sigma_x.hole_orig',subindex=1,default=default_origin(1))
       call get_parameter(origin(2),'opt.sigma_x.hole_orig',subindex=2,default=default_origin(2))
       call get_parameter(origin(3),'opt.sigma_x.hole_orig',subindex=3,default=default_origin(3))

       pi = 4.0_rk*atan(1.0_rk)

       ! df/dx = slope = Dx[0.5*erf(sqrt(pi)*x/len)]= 1/len -> since df=1 it follows dx = len

       ! set outer rectangular box
       sigma_x = 0.5_rk*(erf(sqrt(pi)/len*(X(:,:,:,1) - x10)) - erf(sqrt(pi)/len*(X(:,:,:,1) - x11))) &
               * 0.5_rk*(erf(sqrt(pi)/len*(X(:,:,:,2) - x20)) - erf(sqrt(pi)/len*(X(:,:,:,2) - x21))) &
               * 0.5_rk*(erf(sqrt(pi)/len*(X(:,:,:,3) - x30)) - erf(sqrt(pi)/len*(X(:,:,:,3) - x31)))

       ! dig inner spherical hole
       r = sqrt((X(:,:,:,1)-origin(1))**2 + (X(:,:,:,2)-origin(2))**2 + (X(:,:,:,3)-origin(3))**2)
       !sigma_x = sigma_x - ( 0.5_rk - 0.5_rk*erf(sqrt(pi)/len*(r - rad)) )
       sigma_x = sigma_x * ( 0.5_rk + 0.5_rk*erf(sqrt(pi)/len*(r - rad)) )

    case default

       write(*,*) "error in objective_adjoint_sound_monopole_array.f90:init_sigma_x:unknown opt.sigma_x.type (rectangular/circular/...)"
       stop

    end select

    call store(sigma_x,'sigma_x_function')

  end subroutine init_sigma_x
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_sigma_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:), allocatable                           :: sigma_t_vector
    integer                                                      :: discarded_steps
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! allocates
    allocate(sigma_t(params%time%steps, params%time%subsets))
    allocate(sigma_t_vector(params%time%steps * params%time%subsets))

    sigma_t        = 1.0_rk
    sigma_t_vector = 1.0_rk

    ! discard time steps
    !call get_parameter(discarded_steps,'opt.sigma_t.discarded_steps',default = 0)
    !
    !if (discarded_steps.gt.0) then
    ! sigma_t_vector(1:discarded_steps) = 0.0_rk
    ! sigma_t_vector((params%time%steps*params%time%subsets) - discarded_steps + 1 : (params%time%steps*params%time%subsets)) = 0.0_rk
    !end if

    !    characteristic acoustic path    mask of adjoint forcing "g"
    !    signal ////////////////////|    111111111111111111111111111|
    !    is    //////////////////// |    111111111111111111111111111|
    !    zero ////////////////////  |    111111111111111111111111111|
    !    here//// slope is c ////   |    111111111111111111111111111|
    !   ^   ////////////////////    |   ^111111111111111111111111111|
    ! x |  ////////////////////     | x |111111111111111111111111111|
    !   | ////////////////////      |   |111111111111111111111111111|
    !   ------>              |------|   ------>              |------|
    !       t                extend/c       t                extend/c

    ! reshape to array
    sigma_t = reshape(sigma_t_vector,(/params%time%steps,params%time%subsets/))

    call store(sigma_t,'sigma_t_function')

  end subroutine init_sigma_t
!!!=================================================================================================

end module objective_adjoint_sound_monopole_array
