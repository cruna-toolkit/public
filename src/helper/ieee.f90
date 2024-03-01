module ieee

  use parameter

  implicit none

  public  :: is_nan
  public  :: is_inf

  private :: is_nan_rk, is_nan_ck
  private :: is_inf_rk, is_inf_ck

  interface is_nan
     module procedure is_nan_rk, is_nan_ck
  end interface

  interface is_inf
     module procedure is_inf_rk, is_inf_ck
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical function is_nan_rk(x)
    real(kind=rk), intent(in)    :: x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    is_nan_rk = isnan(x)
  end function is_nan_rk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical function is_nan_ck(x)
    complex(kind=rk), intent(in)    :: x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                   :: xabs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    xabs = abs(x) ! real or imag part is nan
    is_nan_ck = isnan(xabs)
  end function is_nan_ck

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical function is_inf_rk(x)
    real(kind=rk), intent(in)    :: x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(x.gt.huge(x)) then
       is_inf_rk = .true.
    else
       is_inf_rk = .false.
    end if
  end function is_inf_rk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  logical function is_inf_ck(x)
    complex(kind=rk), intent(in)    :: x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                   :: xabs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    xabs = abs(x)
    if(xabs.gt.huge(xabs)) then
       is_inf_ck = .true.
    else
       is_inf_ck = .false.
    end if
  end function is_inf_ck

end module ieee
