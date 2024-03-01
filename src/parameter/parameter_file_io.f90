module parameter_file_io

  private

  public :: read_parameter_file

contains 

!!!=================================================================================================
  subroutine read_parameter_file(parameter_list,parameter_file_name)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! general constants for the parameter file
    integer         , parameter                                   :: parameter_max_file_length = 2048
    character(len=1), parameter                                   :: parameter_comment         = "%"
    character(len=1), parameter                                   :: parameter_equal           = "="
    character(len=1), parameter                                   :: parameter_seperator       = ";"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*), dimension(:,:), allocatable, intent(inout)  :: parameter_list
    character(len=*)                             , intent(in)     :: parameter_file_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=512), dimension(parameter_max_file_length)      :: parameter_strings                      
    character(len=512)                                            :: parameter_string
    integer                                                       :: io_unit, io_stat
    integer                                                       :: d,i,l,p,s,q
    integer                                                       :: no_lines, no_parameters
    logical                                                       :: parameter_file_exist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    parameter_list = ""

!!! CHECK EXISTENCE OF PARAMETER_FILE
    inquire(file=trim(parameter_file_name), exist=parameter_file_exist)
    if (parameter_file_exist.neqv..true.) then
       write(*,*) "error in read_parameter_file_io.f90:read_parameter_file:parameter file <",trim(parameter_file_name),"> not found"
       stop
    end if

!!! READ PARAMETER FILE
    ! open parameter file
    open(newunit=io_unit, file=trim(parameter_file_name), iostat=io_stat, status='old')

    ! read into parameter_strings (loop)
    l = 0
    do 
       l = l + 1

       read(unit=io_unit, fmt="(A)", iostat=io_stat) parameter_string

       if (io_stat.ne.0) then
          exit ! end of file
       end if

       if (l.gt.parameter_max_file_length) then 
          write(*,*)
          write(*,"(A,A,A,I4.4,A)") "error in read_parameter_file_io.f90:read_parameter_file: <",trim(parameter_file_name),"> too long (max. ",parameter_max_file_length," ll)"
          stop
       end if

       parameter_strings(l) = parameter_string
    enddo
    no_lines = l - 1

    !debug:
    !write(*,"(A,A,A,I4.4,A,I4.4,A)") " params: read ", trim(parameter_file_name), " with ", no_lines," (max. ", size(parameter_strings), ") lines"

    ! close parameter file
    close(io_unit, iostat=io_stat)

!!! UNPACK PARAMETER
    p = 0
    do l=1,no_lines
       ! check for comment lines and empty lines
       parameter_string = trim(adjustl(parameter_strings(l)))

       if ((parameter_string(1:1).ne.parameter_comment).and.(len_trim(parameter_string).gt.0)) then
          ! continue for non-comment entries
          p = p + 1 

          ! remove trailing comment
          i = scan(parameter_string,parameter_comment,.false.) 
          if (i.gt.0) then ! comment present
             parameter_string = trim(adjustl(parameter_string(1:i-1)))
          end if

          ! extract parameter name
          i = scan(parameter_string,parameter_equal,.false.) 
          if (i.eq.0) then ! no equal present
             write(*,*)
             write(*,*) "error in read_parameter_file.f90:read_parameter_file_io:syntax error in <", parameter_file_name, ">: ", trim(adjustl(parameter_strings(l)))
             stop
          end if
          parameter_list(p,1) = trim(parameter_string(1:i-1))
          parameter_string    = trim(adjustl(parameter_string(i+1:len_trim(parameter_string))))

          ! check for additional equals
          i = scan(parameter_string,parameter_equal,.false.) 
          if (i.ne.0) then ! no equal present
             write(*,*)
             write(*,*) "error in read_parameter_file.f90:read_parameter_file_io:syntax error in <", parameter_file_name, ">: ", trim(adjustl(parameter_strings(l)))
             stop
          end if

          ! check for double entries
          do q = 1,p-1
             if (trim(adjustl(parameter_list(p,1))).eq.trim(adjustl(parameter_list(q,1)))) then
                write(*,*)
                write(*,*) "error in read_parameter_file.f90:read_parameter_file_io:parameter has been defined more than once:",trim(adjustl(parameter_list(p,1)))
                stop
             end if
          end do

          ! add closing seperator
          parameter_string = trim(parameter_string)//parameter_seperator

          ! extract parameter value
          d = scan(parameter_string,parameter_seperator,.false.)

          ! extract (sub-)values
          s = 1
          do while (s.le.size(parameter_list,2))
             s = s + 1
             parameter_list(p,s) = trim(parameter_string(1:d-1))
             parameter_string    = trim(adjustl(parameter_string(d+1:len(parameter_string))))

             if (len_trim(parameter_string).eq.0) then ! no more sub parameters
                exit
             end if

             d = scan(parameter_string,parameter_seperator,.false.)
          end do

          !debug:
          ! write(*,*) "  > ",trim(parameter_list(p,1)),"=",trim(parameter_list(p,2)),";",trim(parameter_list(p,3)),";",trim(parameter_list(p,4)),";",trim(parameter_list(p,5))

       endif
    end do
    no_parameters = p

    write(*,"(A,I3.3,A,A,A,I4.4,A)") " params: read ", no_parameters, " parameter(s) in ", trim(parameter_file_name), " (", no_lines, " ll)"

  end subroutine read_parameter_file
!!!=================================================================================================

end module parameter_file_io
