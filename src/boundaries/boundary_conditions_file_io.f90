module boundary_conditions_file_io

  use parameter

  private

  public :: read_boundary_conditions_file

contains 

!!!=================================================================================================
  subroutine read_boundary_conditions_file(boundary_conditions,boundary_conditions_file_name)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! general constants for the boundary_conditions file
    integer         , parameter                                             :: boundary_conditions_max_file_length = 2048
    integer         , parameter                                             :: boundary_conditions_max_number_sub  = 11
    character(len=1), parameter                                             :: boundary_conditions_comment         = "%"
    character(len=1), parameter                                             :: boundary_conditions_seperator       = ";"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , dimension(:,:), allocatable, intent(inout)            :: boundary_conditions
    character(len=*)                             , intent(in)               :: boundary_conditions_file_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=512), dimension(boundary_conditions_max_file_length)      :: boundary_conditions_strings                      
    character(len=512)                                                      :: boundary_conditions_string
    integer                                                                 :: io_unit, io_stat
    integer                                                                 :: d,i,l,p,s
    integer                                                                 :: no_lines, no_boundary_conditions
    logical                                                                 :: boundary_conditions_file_exist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! CHECK EXISTENCE OF BOUNDARY_CONDITIONS_FILE
    inquire(file=trim(boundary_conditions_file_name), exist=boundary_conditions_file_exist)
    if (boundary_conditions_file_exist.neqv..true.) then
       write(*,*) "bounds: WARNING boundary conditions file <",trim(boundary_conditions_file_name),"> not found"
       allocate(boundary_conditions(0,boundary_conditions_max_number_sub))
       boundary_conditions = 0.0_rk
       return
    end if

!!! READ BOUNDARY_CONDITIONS FILE
    ! open boundary_conditions file
    open(newunit=io_unit, file=trim(boundary_conditions_file_name), iostat=io_stat, status='old')

    ! read into boundary_conditions_strings (loop)
    l = 0
    do 
       l = l + 1

       read(unit=io_unit, fmt="(A)", iostat=io_stat) boundary_conditions_string

       if (io_stat.ne.0) then
          exit ! end of file
       end if

       if (l.gt.boundary_conditions_max_file_length) then 
          write(*,*)
          write(*,"(A,A,A,I4.4,A)") "error in boundary_conditions_file_io.f90:read_boundary_conditions_file: <",trim(boundary_conditions_file_name),"> too long (max. ",boundary_conditions_max_file_length," ll)"
          stop
       end if

       boundary_conditions_strings(l) = boundary_conditions_string
    enddo
    no_lines = l - 1

    !debug:
    !write(*,"(A,A,A,I4.4,A,I4.4,A)") " bounds: read ", trim(boundary_conditions_file_name), " with ", no_lines," (max. ", size(boundary_conditions_strings), ") lines"

    ! close boundary_conditions file
    close(io_unit, iostat=io_stat)

!!! COUNT BOUNDARY CONDITIONS / ALLOCATE FIELD
    p = 0
    do l=1,no_lines
       ! check for comment lines and empty lines
       boundary_conditions_string = trim(adjustl(boundary_conditions_strings(l)))

       if ((boundary_conditions_string(1:1).ne.boundary_conditions_comment).and.(len_trim(boundary_conditions_string).gt.0)) then
          ! boundary condition present
          p = p + 1 
       endif
    end do
    no_boundary_conditions = p

    allocate(boundary_conditions(no_boundary_conditions,boundary_conditions_max_number_sub))
    boundary_conditions = 0.0_rk

    ! debug:
    ! write(*,*) "number of bound conds",no_boundary_conditions

!!! UNPACK BOUNDARY_CONDITIONS
    p = 0
    do l=1,no_lines
       ! check for comment lines and empty lines
       boundary_conditions_string = trim(adjustl(boundary_conditions_strings(l)))

       if ((boundary_conditions_string(1:1).ne.boundary_conditions_comment).and.(len_trim(boundary_conditions_string).gt.0)) then
          ! continue for non-comment entries
          p = p + 1 

          ! remove trailing comment
          i = scan(boundary_conditions_string,boundary_conditions_comment,.false.) 
          if (i.gt.0) then ! comment present
             boundary_conditions_string = trim(adjustl(boundary_conditions_string(1:i-1)))
          end if

          ! add closing seperator
          boundary_conditions_string = trim(boundary_conditions_string)//boundary_conditions_seperator

          ! extract boundary_conditions value
          d = scan(boundary_conditions_string,boundary_conditions_seperator,.false.)

          ! extract (sub-)values
          s = 0
          do while (s.le.size(boundary_conditions,2))
             s = s + 1

             read(boundary_conditions_string(1:d-1),*) boundary_conditions(p,s)
             boundary_conditions_string = trim(adjustl(boundary_conditions_string(d+1:len(boundary_conditions_string))))

             if (len_trim(boundary_conditions_string).eq.0) then ! no more sub boundary_conditionss
                exit
             end if

             d = scan(boundary_conditions_string,boundary_conditions_seperator,.false.)
          end do
       endif
    end do

    write(*,"(A,I3.3,A,A,A,I4.4,A)") " bounds: read ", no_boundary_conditions, " condition(s) in ", trim(boundary_conditions_file_name), " (", no_lines, " ll)"

  end subroutine read_boundary_conditions_file
!!!=================================================================================================

end module boundary_conditions_file_io
