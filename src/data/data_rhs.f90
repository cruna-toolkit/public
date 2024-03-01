module data_rhs

  use parameter

  real(kind=rk), dimension(:,:,:,:,:), allocatable, public :: Q , Qs
  real(kind=rk), dimension(:,:,:,:)  , allocatable, public :: q0, qs0

  private

  public allocate_data_direct
  public allocate_data_adjoint

contains

!!!=================================================================================================
  subroutine allocate_data_direct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical :: const_direct = .false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_parameter(const_direct,'params.opt.const_direct',default = .false.)
    if(const_direct.eqv..true.) then
       allocate( Q(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars,                1))
    else
       allocate( Q(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars,params%time%steps))
    endif

    allocate(q0(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))

  end subroutine allocate_data_direct
!!!=================================================================================================

!!!=================================================================================================
  subroutine allocate_data_adjoint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate( Qs(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars,params%time%steps))
    
    allocate(qs0(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars))

  end subroutine allocate_data_adjoint
!!!=================================================================================================

end module data_rhs
