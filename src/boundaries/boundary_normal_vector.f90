module boundary_normal_vector
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Provides 4D boundary normal vector fields - holding k-vectors values everywhere.                  !
  ! This subroutine is usefull if a single k-vector of the boundary_conditions_array(i,8:10)          !
  ! (see boundary_conditions.f90) is not sufficient.                                                  !
  !                                                                                                   !
  ! LS (2025)                                                                                         !
  !                                                                                                   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use data_geom , only: X
  use parameter , only: rk, params, get_parameter
  use io        , only: store

  private

  public init_orthogonal_and_diagonal_kvector
  ! only works if rectilinear/cuboid grid - no skew distortions allowed; includes:
  ! a) near BC surfaces: 0°,90°,180°,270° k-vectors (orthogonal/rectangular to BC surface)
  ! b) near BC edge (in "2D" diagonal extension where different k-vectors of type a) meet): diagonal k-vectors
  ! c) near BC corner (in "3D" diagonal extension where different k-vectors of type b) meet): diagonal k-vectors

  public init_orthogonal_kvector
  ! only works if rectilinear/cuboid grid - no skew distortions allowed; oncly includes:
  ! a) near BC surfaces: 0°,90°,180°,270° k-vectors (orthogonal/rectangular to BC surface)
  ! b) near BC edge (treaded like one of its neighouring surfaces)
  ! c) near BC corner (treaded like one of its neighouring surfaces)

  public init_hedgehog_kvector
  ! hedgehog (Igel in German), pointing circular/outwards from geom. center

  private normalize_kvector

contains

!!!=================================================================================================
  subroutine init_orthogonal_and_diagonal_kvector(kvector)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable, intent(out)  :: kvector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(3)                                  :: length, center, tmp_abs_kvec
    real(kind=rk)                                                :: x10,x11,x20,x21,x30,x31
    real(kind=rk)                                                :: half_norm_dx, max_abs_value
    integer                                                      :: i,j,k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(kvector(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))

    ! get global block dimensions from parameter.dat
    call get_parameter(x10,'geom.x10',default = 0.0_rk)
    call get_parameter(x20,'geom.x20',default = 0.0_rk)
    call get_parameter(x30,'geom.x30',default = 0.0_rk)
    call get_parameter(x11,'geom.x11',default = 1.0_rk)
    call get_parameter(x21,'geom.x21',default = 1.0_rk)
    call get_parameter(x31,'geom.x31',default = 1.0_rk)

    length(1) = x11 - x10
    length(2) = x21 - x20
    length(3) = x31 - x30

    ! calculate geomtric center
    center(1) = x10 + 0.5_rk*length(1)
    center(2) = x20 + 0.5_rk*length(2)
    center(3) = x30 + 0.5_rk*length(3)

    ! build hedgehog
    kvector(:,:,:,1) = X(:,:,:,1) - center(1)
    kvector(:,:,:,2) = X(:,:,:,2) - center(2)
    kvector(:,:,:,3) = X(:,:,:,3) - center(3)

    ! convert hedgehog to orthogonal (90°) and diagonal (45°) k-vectors near edge or corner
    kvector(:,:,:,1) = kvector(:,:,:,1)/length(1)
    kvector(:,:,:,2) = kvector(:,:,:,2)/length(2)
    kvector(:,:,:,3) = kvector(:,:,:,3)/length(3)
    half_norm_dx = 0.5_rk*min(params%geom%dx1/length(1), params%geom%dx2/length(2), params%geom%dx3/length(3))
    do k = 1, size(kvector, 3)
      do j = 1, size(kvector, 2)
        do i = 1, size(kvector, 1)
            tmp_abs_kvec = abs(kvector(i,j,k,:))
            max_abs_value = maxval(tmp_abs_kvec)
            ! set largest component to 1.0
            ! set other components - if equal length to largest - to 1.0 and otherwise to 0.0
            if (max_abs_value == tmp_abs_kvec(1)) then

              kvector(i,j,k,1) = sign(1.0_rk, kvector(i,j,k,1))
              if (max_abs_value-tmp_abs_kvec(2)<half_norm_dx) then ! diagonal element
                kvector(i,j,k,2) = sign(1.0_rk, kvector(i,j,k,2))
              else
                kvector(i,j,k,2) = 0.0_rk
              end if
              if (max_abs_value-tmp_abs_kvec(3)<half_norm_dx) then ! diagonal element
                kvector(i,j,k,3) = sign(1.0_rk, kvector(i,j,k,3))
              else
                kvector(i,j,k,3) = 0.0_rk
              end if

            elseif (max_abs_value == tmp_abs_kvec(2)) then

              if (max_abs_value-tmp_abs_kvec(1)<half_norm_dx) then ! diagonal element
                kvector(i,j,k,1) = sign(1.0_rk, kvector(i,j,k,1))
              else
                kvector(i,j,k,1) = 0.0_rk
              end if
              kvector(i,j,k,2) = sign(1.0_rk, kvector(i,j,k,2))
              if (max_abs_value-tmp_abs_kvec(3)<half_norm_dx) then ! diagonal element
                kvector(i,j,k,3) = sign(1.0_rk, kvector(i,j,k,3))
              else
                kvector(i,j,k,3) = 0.0_rk
              end if

            else ! max_abs_value == tmp_abs_kvec(3)

              if (max_abs_value-tmp_abs_kvec(1)<half_norm_dx) then ! diagonal element
                kvector(i,j,k,1) = sign(1.0_rk, kvector(i,j,k,1))
              else
                kvector(i,j,k,1) = 0.0_rk
              end if
              if (max_abs_value-tmp_abs_kvec(2)<half_norm_dx) then ! diagonal element
                kvector(i,j,k,2) = sign(1.0_rk, kvector(i,j,k,2))
              else
                kvector(i,j,k,2) = 0.0_rk
              end if
              kvector(i,j,k,3) = sign(1.0_rk, kvector(i,j,k,3))

            end if
        end do
      end do
    end do

    call normalize_kvector(kvector)

    if (params%io%verbosity.ge.2) call store(kvector,'kvector')

  end subroutine init_orthogonal_and_diagonal_kvector
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_orthogonal_kvector(kvector)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable, intent(out)  :: kvector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(3)                                  :: length, center
    real(kind=rk)                                                :: x10,x11,x20,x21,x30,x31
    real(kind=rk)                                                :: max_abs_value
    integer                                                      :: i,j,k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(kvector(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))

    call get_parameter(x10,'geom.x10',default = 0.0_rk)
    call get_parameter(x20,'geom.x20',default = 0.0_rk)
    call get_parameter(x30,'geom.x30',default = 0.0_rk)
    call get_parameter(x11,'geom.x11',default = 1.0_rk)
    call get_parameter(x21,'geom.x21',default = 1.0_rk)
    call get_parameter(x31,'geom.x31',default = 1.0_rk)

    length(1) = x11 - x10
    length(2) = x21 - x20
    length(3) = x31 - x30

    ! calculate geomtric center
    center(1) = x10 + 0.5_rk*length(1)
    center(2) = x20 + 0.5_rk*length(2)
    center(3) = x30 + 0.5_rk*length(3)

    ! build hedgehog
    kvector(:,:,:,1) = X(:,:,:,1) - center(1)
    kvector(:,:,:,2) = X(:,:,:,2) - center(2)
    kvector(:,:,:,3) = X(:,:,:,3) - center(3)

    ! convert hedgehog to strictly orthogonal (90°) k-vector
    ! (only works if rectilinear/cuboid grid)
    kvector(:,:,:,1) = kvector(:,:,:,1)/length(1)
    kvector(:,:,:,2) = kvector(:,:,:,2)/length(2)
    kvector(:,:,:,3) = kvector(:,:,:,3)/length(3)
    do k = 1, size(kvector, 3)
       do j = 1, size(kvector, 2)
          do i = 1, size(kvector, 1)
             max_abs_value = maxval(abs(kvector(i,j,k,:)))
             if (max_abs_value == abs(kvector(i,j,k,1))) then
                kvector(i,j,k,1) = sign(1.0_rk, kvector(i,j,k,1))
                kvector(i,j,k,2) = 0.0_rk
                kvector(i,j,k,3) = 0.0_rk
             elseif (max_abs_value == abs(kvector(i,j,k,2))) then
                kvector(i,j,k,1) = 0.0_rk
                kvector(i,j,k,2) = sign(1.0_rk, kvector(i,j,k,2))
                kvector(i,j,k,3) = 0.0_rk
             else
                kvector(i,j,k,1) = 0.0_rk
                kvector(i,j,k,2) = 0.0_rk
                kvector(i,j,k,3) = sign(1.0_rk, kvector(i,j,k,3))
             end if
          end do
       end do
    end do

    if (params%io%verbosity.ge.2) call store(kvector,'kvector')

  end subroutine init_orthogonal_kvector
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_hedgehog_kvector(kvector)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), allocatable, intent(out)  :: kvector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(3)                                  :: length, center
    real(kind=rk)                                                :: x10,x11,x20,x21,x30,x31
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(kvector(params%geom%n1b,params%geom%n2b,params%geom%n3b,3))

    ! get global block dimensions from parameter.dat
    call get_parameter(x10,'geom.x10',default = 0.0_rk)
    call get_parameter(x20,'geom.x20',default = 0.0_rk)
    call get_parameter(x30,'geom.x30',default = 0.0_rk)
    call get_parameter(x11,'geom.x11',default = 1.0_rk)
    call get_parameter(x21,'geom.x21',default = 1.0_rk)
    call get_parameter(x31,'geom.x31',default = 1.0_rk)

    length(1) = x11 - x10
    length(2) = x21 - x20
    length(3) = x31 - x30

    ! calculate geomtric center
    center(1) = x10 + 0.5_rk*length(1)
    center(2) = x20 + 0.5_rk*length(2)
    center(3) = x30 + 0.5_rk*length(3)

    ! build hedgehog a.k.a. circular kvector field (pointing outwards from geom. center)
    kvector(:,:,:,1) = X(:,:,:,1) - center(1)
    kvector(:,:,:,2) = X(:,:,:,2) - center(2)
    kvector(:,:,:,3) = X(:,:,:,3) - center(3)

    call normalize_kvector(kvector)

    if (params%io%verbosity.ge.2) call store(kvector,'kvector')

  end subroutine init_hedgehog_kvector
!!!=================================================================================================

!!!=================================================================================================
  subroutine normalize_kvector(kvector)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)            :: kvector
    real(kind=rk), dimension(:,:,:), allocatable                :: r_distance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(r_distance(size(kvector, 1), size(kvector, 2), size(kvector, 3)))

    ! Compute magnitude
    r_distance = sqrt(kvector(:,:,:,1)**2 + kvector(:,:,:,2)**2 + kvector(:,:,:,3)**2)

    ! Normalize
    kvector(:,:,:,1) = kvector(:,:,:,1) / r_distance
    kvector(:,:,:,2) = kvector(:,:,:,2) / r_distance
    kvector(:,:,:,3) = kvector(:,:,:,3) / r_distance

    deallocate(r_distance)

  end subroutine normalize_kvector
!!!=================================================================================================

end module boundary_normal_vector