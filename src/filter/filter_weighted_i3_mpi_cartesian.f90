module filter_weighted_i3_mpi_cartesian

  use parameter
  use filter_weighted_x_core  , only: filter_core  => filter_i3 , filter_i3_name => filter_i3_name , filter_i3_type => filter_i3_type , init_filter => init_filter
  use transpose_mpi_cartesian , only: ftranspose => ftranspose_i3, btranspose => btranspose_i3
  USE mpi

  interface filter
     module procedure filter_4d
  end interface filter

contains

!!!=================================================================================================
  subroutine filter_4d(q,weight,strength)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(inout)           :: q
    real(kind=rk), dimension(:,:,:)  , intent(in)              :: weight
    real(kind=rk), dimension(:,:,:)  , intent(in)              :: strength
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), allocatable                 :: fries_2d, q_tilde
    real(kind=rk), dimension(:,:), allocatable                 :: strength_2d
    real(kind=rk), dimension(:,:), allocatable                 :: weight_2d
    integer                                                    :: i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call ftranspose(strength_2d,strength)
    call ftranspose(weight_2d  ,weight  )

    do i = 1,size(q,4)
       call ftranspose(fries_2d,q(:,:,:,i))

       ! backup unfiltered vector
       allocate(q_tilde(size(fries_2d,1),size(fries_2d,2)))
       q_tilde = fries_2d

       ! 1 st call
       call filter_core(fries_2d)
       ! strengt
       fries_2d = strength_2d*fries_2d
       ! 2nd call
       call filter_core(fries_2d)

       ! diff
       fries_2d = q_tilde - fries_2d

       ! back transpose
       call btranspose(q(:,:,:,i),fries_2d)

       deallocate(q_tilde)
    end do

  end subroutine filter_4d
!!!=================================================================================================

end module filter_weighted_i3_mpi_cartesian
