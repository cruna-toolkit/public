module time_series

  use cartesian_single_block
  use data_geom
  use helper
  use io
  use mpi_cartesian
  use parameter

  
  private

  public :: sample
  public :: probe

  real(kind=rk), dimension(:,:)    , allocatable, save :: kinetic_energy

contains

!!!=================================================================================================
  subroutine sample(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(in), optional                      :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                                                :: dx1,dx2,dx3
    real(kind=rk)                                 				 :: x10, x20, x30                     ! global min
    real(kind=rk)                                  				 :: x11, x21, x31                     ! global max
    real(kind=rk)                                                                :: pi
    integer                                                                      :: i,j,k
    integer                                                                      :: n1,n2,n3
    real(kind=rk)                                                                :: rho0
    real(kind=rk), dimension(:,:,:,:), allocatable                               :: u
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! compute and save kinectic energy for every time_step

        
!         write(*,*) 'params%time%nt', params%time%nt
!         write(*,*) 'params%time%steps', params%time%steps    
!         write(*,*) 'shape q', shape(q)    
        
        
        ! allocate u and kinetic_energy
        if(.not.allocated(kinetic_energy)) then 
            allocate(kinetic_energy(1,params%time%steps))
        end if
        allocate(  u(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))
!         write(*,*) 'shape kinetic_energy', shape(kinetic_energy)    
!         write(*,*) 'muh'
        ! compute u
        do i = 1,params%geom%n3b
            u(:,:,i,1) = q(:,:,i,2)/q(:,:,i,1)
            u(:,:,i,2) = q(:,:,i,3)/q(:,:,i,1)
            u(:,:,i,3) = q(:,:,i,4)/q(:,:,i,1)
        end do
                    
        kinetic_energy(:,params%time%nt) =sum(q(:,:,:,1)*(u(:,:,:,1)**2+u(:,:,:,2)**2 +u(:,:,:,3)**2))

!         if(params%parallelism%block_image==1) write(*,*) 'nt = ', params%time%nt   ,'  kinetic_energy(:,params%time%nt)' , kinetic_energy(:,params%time%nt)         
        call allreduce(kinetic_energy(1,params%time%nt), kinetic_energy(1,params%time%nt), 'sum', params%parallelism%block_comm) 		
!         if (params%parallelism%block_image==1) write(*,*) 'nt = ', params%time%nt   ,'  kinetic_energy(:,params%time%nt)' , kinetic_energy(:,params%time%nt)

    if (1.eq.params%parallelism%block_image) then
	call get_parameter(rho0,'init.rho')
	call get_parameter(n1,'geom.n1')
	call get_parameter(n2,'geom.n2')
	call get_parameter(n3,'geom.n3')
	call get_parameter(dx1,'geom.dx1')
	call get_parameter(dx2,'geom.dx2')
        call get_parameter(dx3,'geom.dx3') 
	call get_parameter(x10,'geom.x10',default = 0.0_rk)
	call get_parameter(x20,'geom.x20',default = 0.0_rk)
	call get_parameter(x30,'geom.x30',default = 0.0_rk)

	call get_parameter(x11,'geom.x11',default = 1.0_rk)
	call get_parameter(x21,'geom.x21',default = 1.0_rk)
	call get_parameter(x31,'geom.x31',default = 1.0_rk)
	kinetic_energy(1,params%time%nt) = kinetic_energy(1,params%time%nt)*dx1*dx2*dx3/rho0/2 /((x11-x10)*(x21-x20)*(x31-x30))
! 	write(*,*) 'kin E(t)',kinetic_energy(1,params%time%nt)
        ! store at last time-step or at mod == dsfreq  
        if((params%time%nt.eq.params%time%steps).or.(mod(params%time%nt,params%io%dsfreq).eq.0)) then
            write(*,*) 'start storing'
                   
            call store(kinetic_energy,'sample_time',dset_name_optin="kinetic_energy",file_overwrite_optin=.true.)
        end if

    end if

  end subroutine sample
!!!=================================================================================================

!!!=================================================================================================
  subroutine probe(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(in), optional                      :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a probe func please use over-ride functionality of the makefile
  end subroutine probe
!!!=================================================================================================

end module time_series
