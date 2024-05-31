module boundary_conditions

  use boundary_conditions_file_io
  use data_geom
  use data_rhs , only: q0
  use parameter
  use wall_boundary_conditions
  use nonreflecting
  use nonreflecting_lae

  private

  real(kind=rk), dimension(:,:), allocatable, save, public :: boundary_conditions_array

  public init_boundary_conditions_array
  public set_boundary_condition_q , set_boundary_condition_rhs
  public set_boundary_condition_qs , set_boundary_condition_rhss
  public modify_boundary_conditions_array

contains

!!!=================================================================================================
  subroutine set_boundary_condition_q(q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:),intent(inout)                       :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical, save                                                        :: do_all_initialization=.true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
          ! This is a rhs boundary condition only.

       case (200) 
          ! nonreflecting, characteristic decomposition, target values are inner neighbours i.e. zero gradient
          call characteristic_zero_grad_init(is2,ie2,js2,je2,ks2,ke2,boundary_conditions_array(i,:))
          call characteristic_q(q(is:ie,js:je,ks:ke,:),q(is2:ie2,js2:je2,ks2:ke2,:),boundary_conditions_array(i,8:10))

       case (201) 
          ! nonreflecting, characteristic decomposition, target values are initial conditions
          call characteristic_q(q(is:ie,js:je,ks:ke,:),q0(is:ie,js:je,ks:ke,:),boundary_conditions_array(i,8:10))

       case (301)
          ! nonreflecting for LAE, characteristic decomposition, target values are zero
          call characteristic_q_lae(q(is:ie,js:je,ks:ke,:),q0(is:ie,js:je,ks:ke,:),boundary_conditions_array(i,8:10))

       case default
          write(*,*) "error in boundary_conditions.f90:set_boundary_condition_rhs:boundary type:",boundary_conditions_array(i,1)," not implemented @proc:",params%parallelism%world_image
          stop

       end select

       if (do_all_initialization) then
          do_all_initialization = .false.
       end if

    end do

  end subroutine set_boundary_condition_q
!!!=================================================================================================

!!!=================================================================================================
  subroutine set_boundary_condition_qs(qs,q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:),intent(inout)                       :: qs
    real(kind=rk),dimension(:,:,:,:),intent(in)                          :: q ! from direct calc.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical, save                                                        :: do_all_initialization=.true.
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
          ! This is a rhss boundary condition only.

       case (200) 
          ! nonreflecting, characteristic decomposition, target values are initial conditions
          call characteristic_qs(qs(is:ie,js:je,ks:ke,:),q(is:ie,js:je,ks:ke,:),boundary_conditions_array(i,8:10))
          
       case (201) 
          ! nonreflecting, characteristic decomposition, target values are initial conditions
          call characteristic_qs(qs(is:ie,js:je,ks:ke,:),q(is:ie,js:je,ks:ke,:),boundary_conditions_array(i,8:10))

       case (301)
          ! nonreflecting, characteristic decomposition, target values are initial conditions
          call characteristic_qs_lae(qs(is:ie,js:je,ks:ke,:),q(is:ie,js:je,ks:ke,:),boundary_conditions_array(i,8:10))

       case default
          write(*,*) "error in boundary_conditions.f90:set_boundary_condition_rhs:boundary type:",boundary_conditions_array(i,1)," not implemented @proc:",params%parallelism%world_image
          stop

       end select

       if (do_all_initialization) then
          do_all_initialization = .false.
       end if

    end do

  end subroutine set_boundary_condition_qs
!!!=================================================================================================

!!!=================================================================================================
  subroutine set_boundary_condition_rhs(rhs,q)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:),intent(inout)                       :: rhs
    real(kind=rk),dimension(:,:,:,:),intent(in)                          :: q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical, save                                                        :: do_all_initialization=.true.
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
          ! nonreflecting, characteristic decomposition, target values are inner neighbours i.e. zero gradient (200)
          ! call characteristic_zero_grad_init(is2,ie2,js2,je2,ks2,ke2,boundary_conditions_array(i,:))
          ! call characteristic_rhs(rhs(is:ie,js:je,ks:ke,:),q(is:ie,js:je,ks:ke,:),q(is2:ie2,js2:je2,ks2:ke2,:),boundary_conditions_array(i,8:10))

       case (201)
          ! This is a "q" boundary condition only.
       
       case (301)
          ! This is a "q" boundary condition only.

       case default
          write(*,*) "error in boundary_conditions.f90:set_boundary_condition_rhs:boundary type:",boundary_conditions_array(i,1)," not implemented @proc:",params%parallelism%world_image
          stop

       end select

       if (do_all_initialization) then
          do_all_initialization = .false.
       end if

    end do

  end subroutine set_boundary_condition_rhs
!!!=================================================================================================

!!!=================================================================================================
  subroutine set_boundary_condition_rhss(rhs,qs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:),intent(inout)                       :: rhs
    real(kind=rk),dimension(:,:,:,:),intent(in)                          :: qs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical, save                                                        :: do_all_initialization=.true.
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
          ! This is a "qs" boundary condition only.

       case (201)
          ! This is a "qs" boundary condition only.
      
       case (301)
          ! This is a "q" boundary condition only.

       case default
          write(*,*) "error in boundary_conditions.f90:set_boundary_condition_rhs:boundary type:",boundary_conditions_array(i,1)," not implemented @proc:",params%parallelism%world_image
          stop

       end select

       if (do_all_initialization) then
          do_all_initialization = .false.
       end if

    end do

  end subroutine set_boundary_condition_rhss
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
  subroutine modify_boundary_conditions_array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer,dimension(3)                                             :: shift
    integer                                                          :: i,j,ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    shift = 0

!!! find the shift
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

  end subroutine modify_boundary_conditions_array
!!!=================================================================================================


end module boundary_conditions
