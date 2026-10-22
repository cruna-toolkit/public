module precision

  USE iso_fortran_env , only: real32, real64

  implicit none

  private

  integer         , parameter, public :: rk        = real64 ! << change here only
  character(len=6), parameter, public :: prec_type = merge('real32','real64',rk.eq.real32)

end module precision
