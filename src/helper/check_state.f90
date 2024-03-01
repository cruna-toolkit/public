module check_state

  use cartesian_single_block
  use data_geom
  use discretisation_x
  use governing_equations  
  use helper
  use io
  use parameter
  use parallelism
  use mpi_cartesian
  use transpose_mpi_cartesian

  private

  public :: check_cfl
  public :: check_c
  public :: check_Re

contains 

!!!=================================================================================================
  subroutine check_cfl(q,echo)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(in)  :: q
    logical      , optional                        :: echo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable :: u
    real(kind=rk)                                  :: dx_min, dx_min_loc
    real(kind=rk)                                  :: u_max , u_max_loc
    real(kind=rk)                                  :: cfl
    logical                                        :: echo_cfl
    real(kind=rk)                                  :: x10, x20, x30                     ! global min
    real(kind=rk)                                  :: x11, x21, x31                     ! global max
    real(kind=rk)                                  :: x1_min, x2_min, x3_min            ! global min (distorted)
    real(kind=rk), dimension(:,:,:), allocatable   :: deltaX1, deltaX2, deltaX3
    real(kind=rk), dimension(:,:), allocatable     :: dx_fries_2d
    real(kind=rk), dimension(:), allocatable       :: X_end
    real(kind=rk)                                  :: Lx, Ly, Lz
    integer, dimension(3), save                    :: point_idx = (/0,0,0/)                                                                                   
    integer                                        :: image = 0
    logical                                        :: grid_deformed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.present(echo)) then
       echo_cfl = .true.
    else
       echo_cfl = echo
    end if

    ! calc dx_min
    call get_parameter(grid_deformed,'geom.grid_deformed',default = .false.)
    if (grid_deformed.eqv..true.) then

       if(params%parallelism%block_image==1) then
          write(*,*) 'CFL computing for distorted grid:'
       end if
       ! compute dx, dy, dz for distorted grid
       call get_parameter(x10,'geom.x10',default = 0.0_rk)
       call get_parameter(x20,'geom.x20',default = 0.0_rk)
       call get_parameter(x30,'geom.x30',default = 0.0_rk)

       call get_parameter(x11,'geom.x11',default = 1.0_rk)
       call get_parameter(x21,'geom.x21',default = 1.0_rk)
       call get_parameter(x31,'geom.x31',default = 1.0_rk)

       Lx = x11 - x10
       Ly = x21 - x20
       Lz = x31 - x30

       allocate(deltaX1(size(X,1),size(X,2),size(X,3)))
       allocate(deltaX2(size(X,1),size(X,2),size(X,3)))
       allocate(deltaX3(size(X,1),size(X,2),size(X,3)))

       ! get dx
       call ftranspose_i1(dx_fries_2d,X(:,:,:,1))
       allocate(X_end(size(dx_fries_2d,2)))
       X_end = dx_fries_2d(size(dx_fries_2d,1),:)
       dx_fries_2d(1:(size(dx_fries_2d,1)-1),:) = dx_fries_2d(2:(size(dx_fries_2d,1)),:) - dx_fries_2d(1:(size(dx_fries_2d,1)-1),:) ! x(2:end)-x(1:end-1)
       if (trim(diff_i1_type).eq.'periodic') then
          dx_fries_2d(size(dx_fries_2d,1),:) = Lx - X_end      
       end if
       call btranspose_i1(deltaX1,dx_fries_2d)  
       deallocate(X_end)

       ! get dy
       call ftranspose_i2(dx_fries_2d,X(:,:,:,2))
       allocate(X_end(size(dx_fries_2d,2)))
       X_end = dx_fries_2d(size(dx_fries_2d,1),:)
       dx_fries_2d(1:(size(dx_fries_2d,1)-1),:) = dx_fries_2d(2:(size(dx_fries_2d,1)),:) - dx_fries_2d(1:(size(dx_fries_2d,1)-1),:) ! y(2:end)-y(1:end-1)
       if (trim(diff_i2_type).eq.'periodic') then
          dx_fries_2d(size(dx_fries_2d,1),:) = Ly - X_end      
       end if
       deallocate(X_end)
       call btranspose_i2(deltaX2,dx_fries_2d)  

       ! get dz
       call ftranspose_i3(dx_fries_2d,X(:,:,:,3))
       allocate(X_end(size(dx_fries_2d,2)))
       X_end = dx_fries_2d(size(dx_fries_2d,1),:)
       dx_fries_2d(1:(size(dx_fries_2d,1)-1),:) = dx_fries_2d(2:(size(dx_fries_2d,1)),:) - dx_fries_2d(1:(size(dx_fries_2d,1)-1),:) ! z(2:end)-z(1:end-1)
       if (trim(diff_i2_type).eq.'periodic') then
          dx_fries_2d(size(dx_fries_2d,1),:) = Lz - X_end      
       end if
       deallocate(X_end)
       call btranspose_i3(deltaX3,dx_fries_2d) 

       ! compute minimal distance
       call get_index(point_idx,image,x1_min,deltaX1,'min')
       call get_index(point_idx,image,x2_min,deltaX2,'min')
       call get_index(point_idx,image,x3_min,deltaX3,'min')  

       if(params%parallelism%block_image==1) then
          write(*,*) 'dx1 (min) =', num2str(real(x1_min),'(F0.4)')
          write(*,*) 'dx2 (min) =', num2str(real(x2_min),'(F0.4)')
          write(*,*) 'dx3 (min) =', num2str(real(x3_min),'(F0.4)')        
       end if

       dx_min_loc = minval( (/ x1_min, x2_min, x3_min /) )
       call allreduce(dx_min, dx_min_loc, 'min', params%parallelism%world_comm)
       if (params%parallelism%world_image.eq.1) then
          write(*,*) 'dx (min) = ',num2str(real(dx_min),'(F0.4)')
       end if
    else
       dx_min_loc = minval( (/ params%geom%dx1, params%geom%dx2, params%geom%dx3 /) )
       call allreduce(dx_min, dx_min_loc, 'min', params%parallelism%world_comm)
    end if
    ! calc u_max + c
    allocate(u(size(q,1),size(q,2),size(q,3),3))
    u = give_u(q)

    u_max_loc = maxval(sqrt(u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2) + give_c(q))
    call allreduce(u_max, u_max_loc, 'max', params%parallelism%world_comm)

    ! calc cfl
    cfl = u_max*params%time%dt/dx_min

    if (echo_cfl.eqv..true.) then
       if (params%parallelism%world_image.eq.1) then
          write(*,*) "cfl-number (max)    : ", num2str(real(cfl),'(F0.3)')
       end if
    end if

  end subroutine check_cfl
!!!=================================================================================================

!!!=================================================================================================
  subroutine check_c(q,echo)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(in)  :: q
    logical      , optional                        :: echo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:)  , allocatable :: c
    real(kind=rk)                                  :: c_max,c_max_loc
    logical                                        :: echo_c_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.present(echo)) then
       echo_c_init = .true.
    else
       echo_c_init = echo
    end if

    ! calc c
    allocate(c(size(q,1),size(q,2),size(q,3)))
    c = give_c(q)

    c_max_loc = maxval(c)
    call allreduce(c_max, c_max_loc, 'max', params%parallelism%world_comm)

    if (echo_c_init.eqv..true.) then
       if (params%parallelism%world_image.eq.1) then
          write(*,*) "speed of sound (max): ", num2str(real(c_max),'(F0.3)')
       end if
    end if

  end subroutine check_c
!!!=================================================================================================

!!!=================================================================================================
  subroutine check_Re(echo)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical      , optional                        :: echo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                  :: u0,v0,w0,rho0,mu0,um,uc
    real(kind=rk)                                  :: x10, x20, x30                     
    real(kind=rk)                                  :: x11, x21, x31
    real(kind=rk)                                  :: L1 , L2 , L3, Lm, Lc
    logical                                        :: echo_Re_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (.not.present(echo)) then
       echo_Re_init = .true.
    else
       echo_Re_init = echo
    end if

    select case (trim(viscosity_type)) ! from governing equations

    case("none")

       if (echo_Re_init.eqv..true.) then
          if (params%parallelism%world_image.eq.1) then
             write(*,*) "Reynolds-number: oo (no friction)"
          end if
       end if

    case("newtonian")
       ! velocity
       call get_parameter(u0,'init.u1')
       call get_parameter(v0,'init.u2')
       call get_parameter(w0,'init.u3')

       um = sqrt(u0**2.0_rk + v0**2.0_rk + w0**2.0_rk)

       call get_parameter(uc,'reference.velocity',default = um)

       ! characteristic length
       call get_parameter(x10,'geom.x10',default = 0.0_rk)
       call get_parameter(x20,'geom.x20',default = 0.0_rk)
       call get_parameter(x30,'geom.x30',default = 0.0_rk)

       call get_parameter(x11,'geom.x11',default = 1.0_rk)
       call get_parameter(x21,'geom.x21',default = 1.0_rk)
       call get_parameter(x31,'geom.x31',default = 1.0_rk)

       L1 = (x11 - x10)/2.0_rk
       L2 = (x21 - x20)/2.0_rk
       L3 = (x31 - x30)/2.0_rk
       Lm = maxval((/L1,L2,L3/))

       call get_parameter(Lc,'reference.length',default=Lm)

       ! dynamic viscosity
       call get_parameter( mu0,'material.mu0',default = 0.0_rk)
       call get_parameter(rho0,'init.rho'    )

       if (echo_Re_init.eqv..true.) then
          if (params%parallelism%world_image.eq.1) then
             write(*,*) "Reynolds-number: ", num2str(real(rho0*uc*Lc/mu0),'(F0.8)')
             write(*,*) "          - mu0: ", num2str(real(mu0),'(F0.8)')
             write(*,*) "          -  uc: ", num2str(real(uc) ,'(F0.8)')
             write(*,*) "          -  Lc: ", num2str(real(Lc) ,'(F0.8)')
          end if
       end if


    case default
       write(*,*) "error in check_state.f90:check_Re:unknown viscosity type in governing equations"
       stop

    end select

  end subroutine check_Re
!!!=================================================================================================

end module check_state
