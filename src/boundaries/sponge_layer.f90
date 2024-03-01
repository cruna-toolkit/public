module sponge_layer

  use parameter

  character(len=max_length_parameter), parameter, public :: sponge_direct_name  = 'empty sponge (default)'
  character(len=max_length_parameter), parameter, public :: sponge_adjoint_name = 'empty sponge (default)'

  private

  public sponge
  public sponge_adjoint

contains

!!!=================================================================================================
  subroutine sponge(rhs,q,q0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)             :: rhs
    real(kind=rk), dimension(:,:,:,:), intent(in), optional      :: q,q0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rhs = rhs

  end subroutine sponge
!!!=================================================================================================

!!!=================================================================================================
  subroutine sponge_adjoint(rhs,qs,q0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)             :: rhs
    real(kind=rk), dimension(:,:,:,:), intent(in), optional      :: qs,q0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rhs = rhs

  end subroutine sponge_adjoint
!!!=================================================================================================

end module sponge_layer
