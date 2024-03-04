module volume_penalization_classic

  use data_geom
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

   ! darcy
    rhs(:,:,:,2) = rhs(:,:,:,2) - penalization(:,:,:,4)*(q(:,:,:,2)/q(:,:,:,1))
    rhs(:,:,:,3) = rhs(:,:,:,3) - penalization(:,:,:,5)*(q(:,:,:,3)/q(:,:,:,1))
    rhs(:,:,:,4) = rhs(:,:,:,4) - penalization(:,:,:,6)*(q(:,:,:,4)/q(:,:,:,1))

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
    real(kind=rk)                                           :: phi
    real(kind=rk)                                           :: p_start,p_end
    real(kind=rk)                                           :: fwidth, darcy_amp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    call get_parameter(phi,'penalization.phi',default = 1.0_rk)

    penalization(:,:,:,1) = 1.0_rk ! effective volume phi
    penalization(:,:,:,2) = 1.0_rk ! 1/phi
    penalization(:,:,:,3) = 0.0_rk ! darcy rho
    penalization(:,:,:,4) = 0.0_rk ! darcy rho u
    penalization(:,:,:,5) = 0.0_rk ! darcy rho w
    penalization(:,:,:,6) = 0.0_rk ! darcy rho v
    penalization(:,:,:,7) = 0.0_rk ! signed dist. function
    
    
    penalization(:,:,:,1) = phi                       ! effective volume phi
    penalization(:,:,:,2) = 1.0_rk/phi ! 1/phi
    
    call get_parameter(p_start,'geom.p_start')
    call get_parameter(p_end  ,'geom.p_end')

    call get_parameter(fwidth   ,'geom.fwidth_dx')
    call get_parameter(darcy_amp,'geom.darcy_amp')

    ! darcy (material)
     penalization(:,:,:,4) = &
         0.5_rk*(tanh(1.0_rk/fwidth*(X(:,:,:,3) - (-1e4)))   - tanh(1.0_rk/fwidth*(X(:,:,:,3) - 1e4))) &
       * 0.5_rk*(tanh(1.0_rk/fwidth*(X(:,:,:,2) - (-1e4)))   - tanh(1.0_rk/fwidth*(X(:,:,:,2) - 1e4))) &
       * 0.5_rk*(tanh(1.0_rk/fwidth*(X(:,:,:,1) - p_start )) - tanh(1.0_rk/fwidth*(X(:,:,:,1) - p_end   )))

!          0.5_rk*(tanh(1.0_rk/fwidth*(X(:,:,:,3) - (-1e4)))   - tanh(1.0_rk/fwidth*(X(:,:,:,3) - 1e4))) &
!        * 0.5_rk*(tanh(1.0_rk/fwidth*(X(:,:,:,2) - (-1e4)))   - tanh(1.0_rk/fwidth*(X(:,:,:,2) - 1e4))) &
!        * 0.5_rk*(tanh(1.0_rk/fwidth*(X(:,:,:,1) - p_start )) - tanh(1.0_rk/fwidth*(X(:,:,:,1) - p_end   )))

    penalization(:,:,:,4) = darcy_amp*penalization(:,:,:,4)
    penalization(:,:,:,5) =           penalization(:,:,:,4)
    penalization(:,:,:,6) =           penalization(:,:,:,4)
    

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

end module volume_penalization_classic
