module data_geom

  use parameter

  real(kind=rk), dimension(:,:,:,:,:), allocatable, public :: X_metric
  real(kind=rk), dimension(:,:,:,:)  , allocatable, public :: X
  integer      , dimension(:,:,:,:)  , allocatable, public :: Xi
  real(kind=rk), dimension(:,:,:)    , allocatable, public :: X_jacobian
 
  private

  public allocate_geometry

contains

!!!=================================================================================================
  subroutine allocate_geometry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical                                                  :: grid_deformed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate( X(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))
    allocate(Xi(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))

    call get_parameter(grid_deformed,'geom.grid_deformed',default = .false.)
    if (grid_deformed.eqv..true.) then
        if ((params%io%verbosity.ge.1).and.(params%parallelism%world_image.eq.1)) then
          write(*,*) " allocate metric coefficients" 
        end if
        allocate(  X_metric(params%geom%n1b,params%geom%n2b,params%geom%n3b,3,3))
        allocate(X_jacobian(params%geom%n1b,params%geom%n2b,params%geom%n3b    ))        
    end if 
  end subroutine allocate_geometry
!!!=================================================================================================

end module data_geom
