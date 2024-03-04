module volume_penalization

  use volume_penalization_classic
!  use volume_penalization_stl

  !   use data_geom
  !   use io
  !   use parameter
  !   use stl

  !   real(kind=rk), dimension(:,:,:,:), allocatable, save, public :: penalization

  !   private

  !   public allocate_volume_penalization
  !   public apply_volume_penalization
  !   public backup_volume_penalization
  !   public calc_grad_volume_penalization
  !   public init_volume_penalization
  !   public store_volume_penalization
  !   public restore_volume_penalization
  !   public update_volume_penalization

  ! contains

  ! !!!=================================================================================================
  !   subroutine allocate_volume_penalization
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !     allocate(penalization(params%geom%n1b,params%geom%n2b,params%geom%n3b,7))

  !   end subroutine allocate_volume_penalization
  ! !!!=================================================================================================

  ! !!!=================================================================================================
  !   subroutine apply_volume_penalization(rhs,q)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     real(kind=rk), dimension(:,:,:,:), intent(inout)                    :: rhs
  !     real(kind=rk), dimension(:,:,:,:), intent(in)                       :: q
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !     ! do nothing, if you need a forcing please use over-ride functionality of the makefile
  !     ! rhs = rhs

  !   end subroutine apply_volume_penalization
  ! !!!=================================================================================================

  ! !!!=================================================================================================
  !   subroutine backup_volume_penalization
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  !   end subroutine backup_volume_penalization
  ! !!!=================================================================================================

  ! !!!=================================================================================================
  !   subroutine calc_grad_volume_penalization(Qs,Q)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     real(kind=rk), dimension(:,:,:,:,:),intent(in),optional :: Qs
  !     real(kind=rk), dimension(:,:,:,:,:),intent(in),optional :: Q      
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  !   end subroutine calc_grad_volume_penalization
  ! !!!=================================================================================================

  ! !!!=================================================================================================
  !   subroutine init_volume_penalization
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     real(kind=rk), dimension(:,:,:,:), allocatable          :: Xrot
  !     real(kind=rk)                                           :: phi,x10,x11,x20,x21,x30,x31,len
  !     real(kind=rk)                                           :: xt,yt,zt
  !     real(kind=rk)                                           :: ax,ay,az,cos_a,sin_a                ! angles
  !     real(kind=rk)                                           :: pi
  !     integer                                                 :: i,j,k
  !     character(len=max_length_fname)                         :: fname
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !     allocate(Xrot(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))
  !     Xrot = X

  !     ! parameter
  !     call get_parameter(ax,'room.rot_angles',1,default = 0.0_rk)
  !     call get_parameter(ay,'room.rot_angles',2,default = 0.0_rk)
  !     call get_parameter(az,'room.rot_angles',3,default = 0.0_rk)

  !     call get_parameter(x10,'room.x10')
  !     call get_parameter(x11,'room.x11')
  !     call get_parameter(x20,'room.x20')
  !     call get_parameter(x21,'room.x21')
  !     call get_parameter(x30,'room.x30')
  !     call get_parameter(x31,'room.x31')

  !     call get_parameter(phi,'penalization.phi',default = 1.0_rk)
  !     call get_parameter(len,'room.len') ! flank width

  !     ! deg2rad
  !     pi = 4.0_rk*atan(1.0_rk)
  !     ax = -ax * pi / 180.0_rk
  !     ay = -ay * pi / 180.0_rk
  !     az = -az * pi / 180.0_rk

  !     ! shift grid to rotation point (room center)
  !     Xrot(:,:,:,1) = Xrot(:,:,:,1) - (x11 + x10)/2.0_rk
  !     Xrot(:,:,:,2) = Xrot(:,:,:,2) - (x21 + x20)/2.0_rk
  !     Xrot(:,:,:,3) = Xrot(:,:,:,3) - (x31 + x30)/2.0_rk

  !     ! rotate (x,x1)
  !     cos_a = cos(ax)
  !     sin_a = sin(ax)
  !     do k = 1,params%geom%n3b
  !        do j = 1,params%geom%n2b
  !           do i = 1,params%geom%n1b
  !              xt = 1.0_rk*Xrot(i,j,k,1) + 0.0_rk*Xrot(i,j,k,2) + 0.0_rk*Xrot(i,j,k,3)
  !              yt = 0.0_rk*Xrot(i,j,k,1) + cos_a *Xrot(i,j,k,2) - sin_a *Xrot(i,j,k,3)
  !              zt = 0.0_rk*Xrot(i,j,k,1) + sin_a *Xrot(i,j,k,2) + cos_a *Xrot(i,j,k,3)
  !              Xrot(i,j,k,1) = xt
  !              Xrot(i,j,k,2) = yt
  !              Xrot(i,j,k,3) = zt
  !           end do
  !        end do
  !     end do

  !     ! rotate (y,x2)
  !     cos_a = cos(ay)
  !     sin_a = sin(ay)    
  !     do k = 1,params%geom%n3b
  !        do j = 1,params%geom%n2b
  !           do i = 1,params%geom%n1b
  !              xt =  cos_a *Xrot(i,j,k,1) + 0.0_rk*Xrot(i,j,k,2) + sin_a *Xrot(i,j,k,3)
  !              yt =  0.0_rk*Xrot(i,j,k,1) + 1.0_rk*Xrot(i,j,k,2) + 0.0_rk*Xrot(i,j,k,3)
  !              zt = -sin_a *Xrot(i,j,k,1) + 0.0_rk*Xrot(i,j,k,2) + cos_a *Xrot(i,j,k,3)
  !              Xrot(i,j,k,1) = xt
  !              Xrot(i,j,k,2) = yt
  !              Xrot(i,j,k,3) = zt
  !           end do
  !        end do
  !     end do

  !     ! rotate (z,x3)
  !     cos_a = cos(az)
  !     sin_a = sin(az)
  !     do k = 1,params%geom%n3b
  !        do j = 1,params%geom%n2b
  !           do i = 1,params%geom%n1b
  !              xt = cos_a *Xrot(i,j,k,1) - sin_a *Xrot(i,j,k,2) + 0.0_rk*Xrot(i,j,k,3)
  !              yt = sin_a *Xrot(i,j,k,1) + cos_a *Xrot(i,j,k,2) + 0.0_rk*Xrot(i,j,k,3)
  !              zt = 0.0_rk*Xrot(i,j,k,1) + 0.0_rk*Xrot(i,j,k,2) + 1.0_rk*Xrot(i,j,k,3)             
  !              Xrot(i,j,k,1) = xt
  !              Xrot(i,j,k,2) = yt
  !              Xrot(i,j,k,3) = zt
  !           end do
  !        end do
  !     end do

  !     ! UNshift grid to rotation point (room center)
  !     Xrot(:,:,:,1) = Xrot(:,:,:,1) + (x11 + x10)/2.0_rk
  !     Xrot(:,:,:,2) = Xrot(:,:,:,2) + (x21 + x20)/2.0_rk
  !     Xrot(:,:,:,3) = Xrot(:,:,:,3) + (x31 + x30)/2.0_rk

  !     penalization(:,:,:,1) = 1.0_rk ! effective volume phi
  !     penalization(:,:,:,2) = 1.0_rk ! 1/phi
  !     penalization(:,:,:,3) = 0.0_rk ! darcy rho
  !     penalization(:,:,:,4) = 0.0_rk ! darcy rho u
  !     penalization(:,:,:,5) = 0.0_rk ! darcy rho w
  !     penalization(:,:,:,6) = 0.0_rk ! darcy rho v
  !     penalization(:,:,:,7) = 0.0_rk ! signed dist. function

  !     len = len * params%geom%dx1

  !     !    penalization(:,:,:,1) = 0.5_rk*(erf(sqrt(pi)/len*(X(:,:,:,1) - x10)) - erf(sqrt(pi)/len*(X(:,:,:,1) - x11))) &
  !     !                          * 0.5_rk*(erf(sqrt(pi)/len*(X(:,:,:,2) - x20)) - erf(sqrt(pi)/len*(X(:,:,:,2) - x21))) &
  !     !                          * 0.5_rk*(erf(sqrt(pi)/len*(X(:,:,:,3) - x30)) - erf(sqrt(pi)/len*(X(:,:,:,3) - x31)))

  !     penalization(:,:,:,1) = 0.5_rk*(tanh(1.0_rk/len*(Xrot(:,:,:,1) - x10)) - tanh(1.0_rk/len*(Xrot(:,:,:,1) - x11))) &
  !                           * 0.5_rk*(tanh(1.0_rk/len*(Xrot(:,:,:,2) - x20)) - tanh(1.0_rk/len*(Xrot(:,:,:,2) - x21))) &
  !                           * 0.5_rk*(tanh(1.0_rk/len*(Xrot(:,:,:,3) - x30)) - tanh(1.0_rk/len*(Xrot(:,:,:,3) - x31)))    

  !     do k=1,params%geom%n3b
  !        do j=1,params%geom%n2b
  !           do i=1,params%geom%n1b
  !              if (penalization(i,j,k,1).lt.(phi)) then
  !                 penalization(i,j,k,1) = phi
  !              end if
  !              penalization(i,j,k,2) = 1.0_rk / penalization(i,j,k,1)
  !           end do
  !        end do
  !     end do

  !     call store(penalization,'penalization')

  !   end subroutine init_volume_penalization
  ! !!!=================================================================================================

  ! !!!=================================================================================================
  !   subroutine store_volume_penalization
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  !   end subroutine store_volume_penalization
  ! !!!=================================================================================================

  ! !!!=================================================================================================
  !   subroutine restore_volume_penalization
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  !   end subroutine restore_volume_penalization
  ! !!!=================================================================================================

  ! !!!=================================================================================================
  !   subroutine update_volume_penalization(d)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     real(kind=rk), dimension(:,:,:,:),intent(in),optional  :: d
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     ! do nothing, if you need a volume_penalization please use over-ride functionality of the makefile
  !   end subroutine update_volume_penalization
  ! !!!=================================================================================================

end module volume_penalization
