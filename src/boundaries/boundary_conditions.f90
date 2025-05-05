module boundary_conditions

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                                                   !
  ! Please note:                                                                                      !
  ! The subroutines *_rhs or *_rhss act inside rhs after the (interior) rhs is calculated.            !
  ! => next timestep fully includes their effect.                                                     !
  !                                                                                                   !
  ! *_q or *_qs routines act on q inside rhs before the rhs calc. and outside on the next timestep.   !
  ! Optionally, some *_q / *_qs types can be called selectively only inside and/or only outside:      !
  ! Per default, outside_rhs_only and inside_rhs_only are false => next timestep is fully BC conform. !
  ! But, e.g., if inside_rhs_only is set true => next timestep only partially includes the BC effect. !
  !                                                                                                   !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use boundary_conditions_file_io
  use data_geom
  use data_reference , only: qref
  use parallelism , only: bcast, stop_cruna
  use parameter
  use wall_boundary_conditions
  use nonreflecting
  use nonreflecting_lae
  use sponge_layer , only: sponge_q, sponge_qs, sponge_rhs, sponge_rhss
  use shape_functions , only: init_quadratic_shape_fct, init_cos_shape_fct
  use boundary_normal_vector , only: init_orthogonal_and_diagonal_kvector, init_orthogonal_kvector, init_hedgehog_kvector

  private

  ! all fields which are initialized once and kept via save attribute
  real(kind=rk), dimension(:,:), allocatable, save, public :: boundary_conditions_array
  real(kind=rk), dimension(:,:,:), allocatable, save       :: shape
  real(kind=rk), dimension(:,:,:,:), allocatable, save     :: kvector
  logical, save                                            :: outside_rhs_only=.false., inside_rhs_only=.false.

  ! applied inside the rhs to q input and outside the rhs to the new timestep
  public set_boundary_condition_q , set_boundary_condition_qs

  ! applied inside the rhs to the newly calc. rhs itself => no need to apply to new timestep
  ! (it is already inherited by rhs)
  public set_boundary_condition_rhs, set_boundary_condition_rhss

  ! this block is always called together once as init. in computation_direct, computation_adjoint, cruna_hpc, and cruna_adj
  public init_boundary_conditions_array ! if rank=0 alloc. and init.
  public spread_boundary_conditions_array ! if rank\=0 alloc. and bcast from rank=0 to others: todo join with above
  public modify_boundary_conditions_array ! index fixes (global>local,... )
  private init_boundary_conditions ! called at end of modify_boundary_conditions_array

contains

!!!=================================================================================================
  subroutine set_boundary_condition_q(q,outside_rhs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:),intent(inout)                       :: q
    logical,optional,intent(in)                                          :: outside_rhs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical                                                              :: inside_rhs
    integer                                                              :: i,is,ie,js,je,ks,ke
    integer                                                              :: is2,ie2,js2,je2,ks2,ke2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1,size(boundary_conditions_array,1)
      !if(params%parallelism%block_image.eq.1) print*,"BCnr",i,":",boundary_conditions_array(i,1)

      is = nint(boundary_conditions_array(i,2))
      ie = nint(boundary_conditions_array(i,3))
      js = nint(boundary_conditions_array(i,4))
      je = nint(boundary_conditions_array(i,5))
      ks = nint(boundary_conditions_array(i,6))
      ke = nint(boundary_conditions_array(i,7))

      select case(nint(boundary_conditions_array(i,1)))
      case (0)
        ! empty, do nothing

      case (100)
        ! wall, slip, no care about thermodynamics
        ! This is a rhs boundary condition only. See content there.

      case (200)
        ! nonreflecting, characteristic decomposition, target values are inner neighbours i.e. zero gradient
        call characteristic_zero_grad_init(is2,ie2,js2,je2,ks2,ke2,boundary_conditions_array(i,:))
        call characteristic_q(q(is:ie,js:je,ks:ke,:),q(is2:ie2,js2:je2,ks2:ke2,:),boundary_conditions_array(i,8:10))

      case (201)
        ! nonreflecting, characteristic decomposition, target values are initial conditions
        call characteristic_q_k1d(q(is:ie,js:je,ks:ke,:),qref(is:ie,js:je,ks:ke,:),boundary_conditions_array(i,8:10))

      case (202)
        ! same as 201 exept kvector - to be consistent with 203
        call characteristic_q_k4d(q(is:ie,js:je,ks:ke,:),qref(is:ie,js:je,ks:ke,:),kvector(is:ie,js:je,ks:ke,:))

      case (203)
        ! nonreflecting zonal characteristic acting on q, target values are initial conditions
        ! attention: this type is applied to the whole domain and computationally demanding.

        if (present(outside_rhs)) then
          inside_rhs=.not.outside_rhs
        else
          inside_rhs=.true. ! default
        end if

        if (inside_rhs) then
          if (.not.outside_rhs_only) call characteristic_q_k4d(q,qref,kvector,shape)
        else ! outside rhs
          if (.not.inside_rhs_only) call characteristic_q_k4d(q,qref,kvector,shape)
        end if

      case (204)
        ! same as 201 exept rhs as inout argument (instead of q). See content there.

      case (205)
        ! same as 202 exept rhs as inout argument (instead of q). See content there.

      case (206)
        ! same as 203 exept rhs as inout argument (instead of q). See content there.

      case (301)
        ! nonreflecting for LAE, characteristic decomposition, target values are zero
        call characteristic_q_lae(q(is:ie,js:je,ks:ke,:),qref(is:ie,js:je,ks:ke,:),boundary_conditions_array(i,8:10))

      case (400)
        ! classic sponge acting on rhs (see content there).

      case (401)
        ! sponge acting on q
        call sponge_q(q,shape,qref)

      case default
        write(*,*) "error in boundary_conditions.f90:set_boundary_condition_q:boundary type:",boundary_conditions_array(i,1)," not implemented @proc:",params%parallelism%world_image
        call stop_cruna

      end select
    end do

  end subroutine set_boundary_condition_q
!!!=================================================================================================

!!!=================================================================================================
  subroutine set_boundary_condition_qs(qs,q,outside_rhs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:),intent(inout)                       :: qs
    real(kind=rk),dimension(:,:,:,:),intent(in)                          :: q ! from direct calc.
    logical,optional,intent(in)                                          :: outside_rhs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical                                                              :: inside_rhs=.true. ! default
    integer                                                              :: i,is,ie,js,je,ks,ke
    integer                                                              :: is2,ie2,js2,je2,ks2,ke2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1,size(boundary_conditions_array,1)

      is = nint(boundary_conditions_array(i,2))
      ie = nint(boundary_conditions_array(i,3))
      js = nint(boundary_conditions_array(i,4))
      je = nint(boundary_conditions_array(i,5))
      ks = nint(boundary_conditions_array(i,6))
      ke = nint(boundary_conditions_array(i,7))

      select case(nint(boundary_conditions_array(i,1)))
      case (0)
        ! empty

      case (100)
        ! This is a rhss boundary condition only.

      case (200)
        ! nonreflecting, characteristic decomposition, target values are inner neighbours i.e. zero gradient
        call characteristic_zero_grad_init(is2,ie2,js2,je2,ks2,ke2,boundary_conditions_array(i,:))
        call characteristic_qs(qs(is:ie,js:je,ks:ke,:),q(is2:ie2,js2:je2,ks2:ke2,:),boundary_conditions_array(i,8:10))

      case (201)
        ! nonreflecting, characteristic decomposition, target values are zero
        call characteristic_qs(qs(is:ie,js:je,ks:ke,:),q(is:ie,js:je,ks:ke,:),boundary_conditions_array(i,8:10))

      case (202)
        ! nonreflecting, characteristic decomposition, target values are zero
        !call characteristic_qs(qs(is:ie,js:je,ks:ke,:),q(is:ie,js:je,ks:ke,:),boundary_conditions_array(i,8:10))
        call characteristic_qs(qs(is:ie,js:je,ks:ke,:),q(is:ie,js:je,ks:ke,:),kvector(is:ie,js:je,ks:ke,:))

      case (203)
        ! characteristic nonrefecting sponge acting on qs, target values are zero (no argument needed here)
        ! attention: this type is applied to the whole domain and computationally demanding.

        if (present(outside_rhs)) then
          inside_rhs=.not.outside_rhs
        else
          inside_rhs=.true. ! default
        end if

        if (inside_rhs) then
          if (.not.outside_rhs_only) call characteristic_qs_k4d(qs,q,kvector,shape)
        else ! outside rhs
          if (.not.inside_rhs_only) call characteristic_qs_k4d(qs,q,kvector,shape)
        end if

      case (204)
        ! same as 201 exept rhss as inout argument (instead of qs). See content there.

      case (205)
        ! same as 202 exept rhss as inout argument (instead of qs). See content there.

      case (206)
        ! same as 203 exept rhss as inout argument (instead of qs). See content there.

      case (301)
        ! nonreflecting, characteristic decomposition, target values are initial conditions
        call characteristic_qs_lae(qs(is:ie,js:je,ks:ke,:),q(is:ie,js:je,ks:ke,:),boundary_conditions_array(i,8:10))

      case (400)
        ! This is a rhss boundary condition only.

      case (401)
        ! sponge acting on qs
        call sponge_qs(qs,shape)

      case default
        write(*,*) "error in boundary_conditions.f90:set_boundary_condition_qs:boundary type:",boundary_conditions_array(i,1)," not implemented @proc:",params%parallelism%world_image
        call stop_cruna

      end select
    end do

  end subroutine set_boundary_condition_qs
!!!=================================================================================================

!!!=================================================================================================
  subroutine set_boundary_condition_rhs(rhs,q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:),intent(inout)                       :: rhs
    real(kind=rk),dimension(:,:,:,:),intent(in)                          :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                              :: i,is,ie,js,je,ks,ke
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1,size(boundary_conditions_array,1)

      is = nint(boundary_conditions_array(i,2))
      ie = nint(boundary_conditions_array(i,3))
      js = nint(boundary_conditions_array(i,4))
      je = nint(boundary_conditions_array(i,5))
      ks = nint(boundary_conditions_array(i,6))
      ke = nint(boundary_conditions_array(i,7))

      select case(nint(boundary_conditions_array(i,1)))
      case (0)
        ! empty

      case (100)
        ! wall, slip, no care about thermodynamics (100)
        call wall_slip_no_thermodynamics(rhs(is:ie,js:je,ks:ke,:),nint(boundary_conditions_array(i,8)))

      case (200)
        ! "q" version of boundary condition active, rhs version below inactive:
        ! nonreflecting, characteristic decomposition, target values are inner neighbors i.e. zero gradient (200)
        ! call characteristic_zero_grad_init(is2,ie2,js2,je2,ks2,ke2,boundary_conditions_array(i,:))
        ! call characteristic_rhs(rhs(is:ie,js:je,ks:ke,:),q(is:ie,js:je,ks:ke,:),q(is2:ie2,js2:je2,ks2:ke2,:),boundary_conditions_array(i,8:10))

      case (201)
        ! This is a "q" boundary condition only.

      case (202)
        ! This is a "q" boundary condition only.

      case (203)
        ! This is a "q" boundary condition only.

      case (204)
        ! nonreflecting, characteristic decomposition, target values are initial conditions
        call characteristic_rhs_k1d(rhs(is:ie,js:je,ks:ke,:),q(is:ie,js:je,ks:ke,:),qref(is:ie,js:je,ks:ke,:),boundary_conditions_array(i,8:10))

      case (205)
        ! same as 204 exept kvector - to be consistent with 206
        call characteristic_rhs_k4d(rhs(is:ie,js:je,ks:ke,:),q(is:ie,js:je,ks:ke,:),qref(is:ie,js:je,ks:ke,:),kvector(is:ie,js:je,ks:ke,:))

      case (206)
        ! nonreflecting zonal characteristic acting on q, target values are initial conditions
        ! attention: this type is applied to the whole domain and computationally demanding.
        call characteristic_rhs_k4d(rhs,q,qref,kvector,shape)

      case (301)
        ! This is a "q" boundary condition only.

      case (400)
        ! classic sponge acting on rhs
        call sponge_rhs(rhs,shape,q,qref)

      case (401)
        ! This is a "q" boundary condition only.

      case default
        write(*,*) "error in boundary_conditions.f90:set_boundary_condition_rhs:boundary type:",boundary_conditions_array(i,1)," not implemented @proc:",params%parallelism%world_image
        call stop_cruna

      end select
    end do

  end subroutine set_boundary_condition_rhs
!!!=================================================================================================

!!!=================================================================================================
  subroutine set_boundary_condition_rhss(rhss,qs,q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:),intent(inout)                       :: rhss
    real(kind=rk),dimension(:,:,:,:),intent(in)                          :: qs
    real(kind=rk),dimension(:,:,:,:),intent(in),optional                 :: q ! from direct calc.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                              :: i,is,ie,js,je,ks,ke
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1,size(boundary_conditions_array,1)

      is = nint(boundary_conditions_array(i,2))
      ie = nint(boundary_conditions_array(i,3))
      js = nint(boundary_conditions_array(i,4))
      je = nint(boundary_conditions_array(i,5))
      ks = nint(boundary_conditions_array(i,6))
      ke = nint(boundary_conditions_array(i,7))

      select case(nint(boundary_conditions_array(i,1)))
      case (0)
        ! empty

      case (100)
        ! wall, slip, no care about thermodynamics (100)
        call wall_slip_no_thermodynamics(rhss(is:ie,js:je,ks:ke,:),nint(boundary_conditions_array(i,8)))

      case (200)
        ! This is a "qs" boundary condition only.

      case (201)
        ! This is a "qs" boundary condition only.

      case (202)
        ! This is a "qs" boundary condition only.

      case (203)
        ! This is a "qs" boundary condition only.

      !case (204)
        ! nonreflecting, characteristic decomposition, target values are initial conditions
        !call characteristic_rhss_k1d(rhss(is:ie,js:je,ks:ke,:),qs(is:ie,js:je,ks:ke,:),q(is:ie,js:je,ks:ke,:),boundary_conditions_array(i,8:10))
        ! subroutine called above not yet implemented in nonreflecting.f90

      !case (205)
        ! same as 204 except kvector - to be consistent with 206
        !call characteristic_rhss_k4d(rhss(is:ie,js:je,ks:ke,:),qs(is:ie,js:je,ks:ke,:),q(is:ie,js:je,ks:ke,:),kvector(is:ie,js:je,ks:ke,:))
        ! subroutine called above not yet implemented in nonreflecting.f90

      !case (206)
        ! nonreflecting zonal characteristic acting on q, target values are initial conditions
        ! attention: this type is applied to the whole domain and computationally demanding.
        !call characteristic_rhss_k4d(rhss,qs,q,kvector,shape)
        ! subroutine called above not yet implemented in nonreflecting.f90

      case (301)
        ! This is a "qs" boundary condition only.

      case (400)
        ! sponge acting on rhss, reference values are zero (no argument needed here)
        if(present(q)) then
          call sponge_rhss(rhss,shape,qs,q)
        else ! for compatibility with old rhs routines
          call sponge_rhss
        end if

      case (401)
        ! This is a "qs" boundary condition only.

      case default
        write(*,*) "error in boundary_conditions.f90:set_boundary_condition_rhss:boundary type:",boundary_conditions_array(i,1)," not implemented @proc:",params%parallelism%world_image
        call stop_cruna

      end select
    end do

  end subroutine set_boundary_condition_rhss
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_boundary_conditions
    ! Only active BC are allocated here.
    ! Empty case below indicate that implementor decided consciously that no allocation is needed.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_parameter)                                  :: type
    integer                                                              :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    interface
      subroutine init_shape_interface(shape)
        import                                                           :: rk
        real(kind=rk), dimension(:,:,:), allocatable, intent(out)        :: shape
      end subroutine init_shape_interface
    end interface
    interface
      subroutine init_kvector_interface(kvector)
        import                                                           :: rk
        real(kind=rk), dimension(:,:,:,:), allocatable, intent(out)      :: kvector
      end subroutine init_kvector_interface
    end interface
    procedure(init_shape_interface), pointer                             :: init_shape => null()
    procedure(init_kvector_interface), pointer                           :: init_kvector => null()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_parameter(type,'sponge.shape',default="quadratic")
    if (trim(type).eq."quadratic") then
      init_shape => init_quadratic_shape_fct
    else if (trim(type).eq."cos") then
      init_shape => init_cos_shape_fct
    else
      print*, "error in boundary_conditions.f90:init_boundary_conditions:shape_fct:",type," not implemented."
      call stop_cruna
    end if

    call get_parameter(type,'sponge.kvector',default="diagonal")
    if (trim(type).eq."diagonal") then
      init_kvector => init_orthogonal_and_diagonal_kvector
    else if (trim(type).eq."orthogonal") then
      init_kvector => init_orthogonal_kvector
    else if (trim(type).eq."hedgehog") then
      init_kvector => init_hedgehog_kvector
    else
      print*, "error in boundary_conditions.f90:init_boundary_conditions:kvector_fct:",type," not implemented."
      call stop_cruna
    end if

    do i = 1,size(boundary_conditions_array,1)
      select case(nint(boundary_conditions_array(i,1)))
      case (0)
      case (100)
      case (200)
      case (201)
      case (202)
        if(.not.allocated(kvector)) call init_kvector(kvector)
      case (203)
        if(.not.allocated(shape)) call init_shape(shape)
        if(.not.allocated(kvector)) call init_kvector(kvector)
        outside_rhs_only = (nint(boundary_conditions_array(i,8))==-1)
        inside_rhs_only = (nint(boundary_conditions_array(i,8))==1)
      case (204)
      case (205)
        if(.not.allocated(kvector)) call init_kvector(kvector)
      case (206)
        if(.not.allocated(shape)) call init_shape(shape)
        if(.not.allocated(kvector)) call init_kvector(kvector)
      case (301)
      case (400)
        if(.not.allocated(shape)) call init_shape(shape)
      case (401)
        if(.not.allocated(shape)) call init_shape(shape)
      case default
        write(*,*) "error in boundary_conditions.f90:init_boundary_conditions:boundary type:",boundary_conditions_array(i,1)," not implemented @proc:",params%parallelism%world_image
        call stop_cruna
      end select

    end do

  end subroutine init_boundary_conditions
!!!=================================================================================================

!!!=================================================================================================
  subroutine init_boundary_conditions_array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_fname)                       :: boundary_conditions_file_name
    integer                                               :: i,j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call get_parameter(boundary_conditions_file_name,'bounds.boundary_conditions_file',default='boundary_conditions.dat')

    ! just done for image 1
    if (params%parallelism%world_image.eq.1) then
       ! read boundary conditions file
       call read_boundary_conditions_file(boundary_conditions_array,boundary_conditions_file_name)

!!! replace inf entries for flexibility
       do j = 2,3
          do i = 1,size(boundary_conditions_array,1)
             if ((boundary_conditions_array(i,j) - 1).eq.boundary_conditions_array(i,j)) then
                boundary_conditions_array(i,j) = params%geom%n1
             end if
          end do
       end do

       do j = 4,5
          do i = 1,size(boundary_conditions_array,1)
             if ((boundary_conditions_array(i,j) - 1).eq.boundary_conditions_array(i,j)) then
                boundary_conditions_array(i,j) = params%geom%n2
             end if
          end do
       end do

       do j = 6,7
          do i = 1,size(boundary_conditions_array,1)
             if ((boundary_conditions_array(i,j) - 1).eq.boundary_conditions_array(i,j)) then
                boundary_conditions_array(i,j) = params%geom%n3
             end if
          end do
       end do

!!! echo boundary conditions
       if (params%io%verbosity.ge.1 .and. params%parallelism%world_image.eq.1) then
          write(*,*) ""
          write(*,*) ' boundary conditions: global list'
          do i = 1,size(boundary_conditions_array,1)
             write(*,"(A,I3.3,A,99F8.3)") "  #",i,":",boundary_conditions_array(i,:)
          end do
          write(*,*) ""
       end if

    end if

  end subroutine init_boundary_conditions_array
!!!=================================================================================================

!!!=================================================================================================
  subroutine spread_boundary_conditions_array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                               :: s1,s2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! broadcast size of boundary conditions array
    s1 = size(boundary_conditions_array,1)
    s2 = size(boundary_conditions_array,2)
    call bcast(s1,0,params%parallelism%world_comm)
    call bcast(s2,0,params%parallelism%world_comm)

    ! allocate boundary_conditions_array (all but rank 0)
    if (params%parallelism%world_image.ne.1) then
       allocate(boundary_conditions_array(s1,s2))
       boundary_conditions_array = 0.0_rk
    end if

    ! broadcast parameter_array
    call bcast(boundary_conditions_array,0,params%parallelism%world_comm)

  end subroutine spread_boundary_conditions_array
!!!=================================================================================================

!!!=================================================================================================
  subroutine modify_boundary_conditions_array ! and init_boundary_conditions
! convert global indices to rank-local indices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer,dimension(3)                                             :: shift
    integer                                                          :: i,j,ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    shift = 0

!!! find the shift (Xi are global indices)
    shift(1) = Xi(1,1,1,1) - 1
    shift(2) = Xi(1,1,1,2) - 1
    shift(3) = Xi(1,1,1,3) - 1

    ! write(*,*) "shift: ",shift,params%parallelism%world_image

!!! do the shift
    do i = 1,size(boundary_conditions_array,1)
       ! i1
       boundary_conditions_array(i,2:3) = boundary_conditions_array(i,2:3) - shift(1)
       ! i2
       boundary_conditions_array(i,4:5) = boundary_conditions_array(i,4:5) - shift(2)
       ! i3
       boundary_conditions_array(i,6:7) = boundary_conditions_array(i,6:7) - shift(3)
    end do

!!! set not members to zero (might be improved be recreating/reallocation of boundary_conditions_array)
    do j = 1,size(boundary_conditions_array,1)
       ! i1
       if (boundary_conditions_array(j,2)>params%geom%n1b.or.boundary_conditions_array(j,3)<1) then
          boundary_conditions_array(j,1)=0
       end if
       ! i2
       if (boundary_conditions_array(j,4)>params%geom%n2b.or.boundary_conditions_array(j,5)<1) then
          boundary_conditions_array(j,1)=0
       end if
       ! i3
       if (boundary_conditions_array(j,6)>params%geom%n3b.or.boundary_conditions_array(j,7)<1) then
          boundary_conditions_array(j,1)=0
       end if
    end do

!!! adapt member bounds (remove outer indices: negative or larger than image points )
    do j = 1,size(boundary_conditions_array,1)
       if (boundary_conditions_array(j,1).ne.0) then
          boundary_conditions_array(j,2) = maxval((/1.0_rk                  , boundary_conditions_array(j,2)/))
          boundary_conditions_array(j,4) = maxval((/1.0_rk                  , boundary_conditions_array(j,4)/))
          boundary_conditions_array(j,6) = maxval((/1.0_rk                  , boundary_conditions_array(j,6)/))
          boundary_conditions_array(j,3) = minval((/real(params%geom%n1b,rk), boundary_conditions_array(j,3)/))
          boundary_conditions_array(j,5) = minval((/real(params%geom%n2b,rk), boundary_conditions_array(j,5)/))
          boundary_conditions_array(j,7) = minval((/real(params%geom%n3b,rk), boundary_conditions_array(j,7)/))
       end if
    end do

!!! echo all boundary conditions
    if (params%io%verbosity.ge.2) then
       call mpi_barrier(params%parallelism%world_comm,ierr)
       do i = 1,params%parallelism%world_size
          call mpi_barrier(params%parallelism%world_comm,ierr)
          if (params%parallelism%world_image.eq.i) then
             write(*,*) 'boundary conditions: list of image ',params%parallelism%world_image
             do j = 1,size(boundary_conditions_array,1)
                write(*,"(A,I3.3,A,99F8.3)") " #",j,":",boundary_conditions_array(j,:)
             end do
          end if
          call mpi_barrier(params%parallelism%world_comm,ierr)
       enddo
    end if

    call init_boundary_conditions

  end subroutine modify_boundary_conditions_array
!!!=================================================================================================

end module boundary_conditions
