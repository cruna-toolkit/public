module sponge_layer
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Select active sponge routines both direct and ajoint                                              !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use parameter

  use sponge_layer_core , only: sponge_rhs, sponge_rhss, sponge_q, sponge_qs
  ! Attention: sponge_rhss requires q as argument which is only provided by set_boundary_condition_rhss calls of new rhs routines.

  !!! For compatibility with old approach where sponge is called directly in rhs:
  ! Should always link to empty_sponge here, because new approach applies sponge only via set_boundary_condition_*
  ! TODO in rhs: remove all sponge calls and add q argument to all set_boundary_condition_rhss calls,
  ! then this section can be removed..
  use sponge_layer_core , only: sponge => empty_sponge_inrhs_warning ! called inside rhs (not any more supported)
  use sponge_layer_core , only: sponge_adjoint => empty_sponge_inrhs_warning ! called inside rhs (not any more supported)

  character(len=max_length_parameter), parameter, public :: sponge_direct_name  = 'sponge_rhs,sponge_q'
  character(len=max_length_parameter), parameter, public :: sponge_adjoint_name = 'sponge_rhss,sponge_qs'

end module sponge_layer