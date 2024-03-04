module volume_penalization_stl

  use io
  use parameter
  use stl

  real(kind=rk), dimension(:,:,:,:), allocatable, save, public :: penalization

  private

  public allocate_volume_penalization
  public apply_volume_penalization
  public backup_volume_penalization
  public calc_grad_volume_penalization
  public init_volume_penalization
  public store_volume_penalization
  public restore_volume_penalization
  public update_volume_penalization

contains

!!!=================================================================================================
  subroutine allocate_volume_penalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(penalization(params%geom%n1b,params%geom%n2b,params%geom%n3b,7))

  end subroutine allocate_volume_penalization
!!!=================================================================================================

!!!=================================================================================================
  subroutine apply_volume_penalization(rhs,q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)                    :: rhs
    real(kind=rk), dimension(:,:,:,:), intent(in)                       :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a forcing please use over-ride functionality of the makefile
    ! rhs = rhs
  end subroutine apply_volume_penalization
!!!=================================================================================================

!!!=================================================================================================
  subroutine backup_volume_penalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  end subroutine backup_volume_penalization
!!!=================================================================================================

!!!=================================================================================================
  subroutine calc_grad_volume_penalization(Qs,Q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:),intent(in),optional :: Qs
    real(kind=rk), dimension(:,:,:,:,:),intent(in),optional :: Q      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  end subroutine calc_grad_volume_penalization
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_volume_penalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)                                           :: phi, fwidth
    integer                                                 :: i,j,k
    character(len=max_length_fname)                         :: fname
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    penalization(:,:,:,1) = 1.0_rk ! effective volume phi
    penalization(:,:,:,2) = 1.0_rk ! 1/phi
    penalization(:,:,:,3) = 0.0_rk ! darcy rho
    penalization(:,:,:,4) = 0.0_rk ! darcy rho u
    penalization(:,:,:,5) = 0.0_rk ! darcy rho w
    penalization(:,:,:,6) = 0.0_rk ! darcy rho v
    penalization(:,:,:,7) = 0.0_rk ! signed dist. function

    call get_parameter(fname,'stl.file')
    call read_stl_file(fname)
    call calc_signed_distance(penalization(:,:,:,7))

    call get_parameter(phi      ,'penalization.phi'       ,default = 1.0_rk)
    call get_parameter(fwidth   ,'penalization.fwidth_dx'                  )


    fwidth = fwidth * params%geom%dx1


    do k=1,params%geom%n3b
       do j=1,params%geom%n2b
          do i=1,params%geom%n1b
             penalization(i,j,k,1) = 1.0_rk - (1.0_rk-phi)*(0.5_rk*(tanh(penalization(i,j,k,7)/fwidth) + 1.0_rk))
             penalization(i,j,k,2) = 1.0_rk/penalization(i,j,k,1)

          end do
       end do
    end do

    call store(penalization,'penalization')

  end subroutine init_volume_penalization
!!!=================================================================================================

!!!=================================================================================================
  subroutine store_volume_penalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  end subroutine store_volume_penalization
!!!=================================================================================================

!!!=================================================================================================
  subroutine restore_volume_penalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  end subroutine restore_volume_penalization
!!!=================================================================================================

!!!=================================================================================================
  subroutine update_volume_penalization(d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:),intent(in),optional  :: d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  end subroutine update_volume_penalization
!!!=================================================================================================

end module volume_penalization_stl
