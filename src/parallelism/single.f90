module single

  use parameter

  private

  public  :: end_parallelism
  public  :: init_parallelism
  public  :: init_topology
  public  :: spread_parameter
  public  :: spread_boundary_conditions_array
  public  :: stop_cruna
  private :: init_parameter_array_block
  private :: spread_parameter_array_world
  private :: spread_parameter_array_block
  private :: spread_boundary_conditions_array_world

contains 

!!!=================================================================================================
  subroutine init_parallelism
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! type
    params%parallelism%type = "single"

    ! set spmd block and image
    params%parallelism%world_comm  = 0
    params%parallelism%world_size  = 1
    params%parallelism%world_image = 0
    params%parallelism%world_block = 0

    params%parallelism%block_comm  = 0
    params%parallelism%block_size  = 1
    params%parallelism%block_image = 0

    ! create parameterlist entries
    call set_parameter(params%parallelism%world_comm  ,'parallelism.world_comm' )
    call set_parameter(params%parallelism%world_size  ,'parallelism.world_size' )
    call set_parameter(params%parallelism%world_image ,'parallelism.world_image')
    call set_parameter(params%parallelism%world_block ,'parallelism.world_block')

    call set_parameter(params%parallelism%block_comm  ,'parallelism.block_comm' )
    call set_parameter(params%parallelism%block_size  ,'parallelism.block_size' )
    call set_parameter(params%parallelism%block_image ,'parallelism.block_image')

  end subroutine init_parallelism
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_topology
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    params%geom%n1b = params%geom%n1
    params%geom%n2b = params%geom%n2
    params%geom%n3b = params%geom%n3

    call set_parameter(params%geom%n1,'params.geom.n1b')
    call set_parameter(params%geom%n2,'params.geom.n2b')
    call set_parameter(params%geom%n3,'params.geom.n3b')

  end subroutine init_topology
!!!=================================================================================================

!!!=================================================================================================
  subroutine spread_parameter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! spread the parameter array to all images (world)
    call spread_parameter_array_world

    ! read parameter files for each block
    call init_parameter_array_block

    ! spread the updated parameter array to all images (block)
    call spread_parameter_array_block

  end subroutine spread_parameter
!!!=================================================================================================

!!!=================================================================================================
  subroutine spread_boundary_conditions_array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! spread the parameter array to all images (world)
    call spread_boundary_conditions_array_world

  end subroutine spread_boundary_conditions_array
!!!=================================================================================================

!!!=================================================================================================
  subroutine spread_parameter_array_world
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! nothing to do here
  end subroutine spread_parameter_array_world
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_parameter_array_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! nothing to do here
  end subroutine init_parameter_array_block
!!!=================================================================================================

!!!=================================================================================================
  subroutine spread_parameter_array_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! nothing to do here
  end subroutine spread_parameter_array_block
!!!=================================================================================================

!!!=================================================================================================
  subroutine spread_boundary_conditions_array_world
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! nothing to do here
  end subroutine spread_boundary_conditions_array_world
!!!=================================================================================================

!!!=================================================================================================
  subroutine stop_cruna
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    stop
  end subroutine stop_cruna
!!!=================================================================================================

!!!=================================================================================================
  subroutine end_parallelism
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! nothing to do here
  end subroutine end_parallelism
!!!=================================================================================================

end module single
