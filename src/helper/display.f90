module display

  use discretisation_t    , only: time_stepper_direct_name, time_stepper_adjoint_name
  use discretisation_x    , only: diff_i1_name, diff_i2_name, diff_i3_name
  use filter_x            , only: filter_i1_name, filter_i2_name, filter_i3_name
  use filter_weighted_x   , only: filter_weighted_i1_name, filter_weighted_i2_name, filter_weighted_i3_name
  use governing_equations , only: equations_name
  use helper              , only: num2str
  use parameter
  use reporting
  use sponge_layer        , only: sponge_direct_name, sponge_adjoint_name

  private

  public :: display_cruna_start
  public :: display_cruna_hpc_start
  public :: display_cruna_cdps_start

contains

!!!=================================================================================================
  subroutine display_cruna_start

    if (params%parallelism%world_image.eq.1) then
       write(*,*) ""
       call display_logo()
       write(*,*) ""
       call display_compile_info
       write(*,*) ""
       call display_date_time()
       write(*,*) ""
       call display_host_name()
       call display_login_name
       call display_case_info
       write(*,*) ""
       call display_parallelism
       write(*,*) ""
       call display_discretisation
       write(*,*) ""
       call display_equations
       write(*,*) ""
       call display_filter
       write(*,*) ""
       call display_sponge
       write(*,*) ""
    end if

  end subroutine display_cruna_start
!!!=================================================================================================

!!!=================================================================================================
  subroutine display_cruna_hpc_start

    if (params%parallelism%world_image.eq.1) then
       write(*,*) ""
       call display_logo_hpc()
       write(*,*) ""
       call display_compile_info
       write(*,*) ""
       call display_date_time()
       write(*,*) ""
       call display_host_name()
       call display_login_name
       call display_case_info
       write(*,*) ""
       call display_parallelism
       write(*,*) ""
       call display_discretisation
       write(*,*) ""
       call display_equations
       write(*,*) ""
       call display_filter
       write(*,*) ""
       call display_sponge
       write(*,*) ""
    end if

  end subroutine display_cruna_hpc_start
!!!=================================================================================================

!!!=================================================================================================
  subroutine display_cruna_cdps_start

    if (params%parallelism%world_image.eq.1) then
       write(*,*) ""
       call display_logo_cdps()
       write(*,*) ""
       call display_compile_info
       write(*,*) ""
       call display_date_time()
       write(*,*) ""
       call display_host_name()
       call display_login_name
       call display_case_info
       write(*,*) ""
       call display_parallelism
       write(*,*) ""
       call display_discretisation
       write(*,*) ""
    end if

  end subroutine display_cruna_cdps_start
!!!=================================================================================================

!!!=================================================================================================
  subroutine display_logo

    write(*,*) "                                                      "
    write(*,*) "              ________  __  ___  _____                "
    write(*,*) "             / ___/ _ \/ / / / |/ / _ |               "
    write(*,*) "            / /__/ , _/ /_/ /    / __ |               "
    write(*,*) "            \___/_/|_|\____/_/|_/_/ |_|               "
    write(*,*) "                                                      "
    write(*,*) "(Compressible Reactive Unsteady Navier-Stokes Adjoint)"
    write(*,*) "                                                      "

  end subroutine display_logo
!!!=================================================================================================

!!!=================================================================================================
  subroutine display_logo_hpc

    write(*,*) "                                                      "
    write(*,*) "     ________  __  ___  _____      __ _____  _____    "
    write(*,*) "    / ___/ _ \/ / / / |/ / _ |    / // / _ \/ ___/    "
    write(*,*) "   / /__/ , _/ /_/ /    / __ |   / _  / ___/ /__      "
    write(*,*) "   \___/_/|_|\____/_/|_/_/ |_|__/_//_/_/   \___/      "
    write(*,*) "                                                      "
    write(*,*) "(Compressible Reactive Unsteady Navier-Stokes Adjoint)"
    write(*,*) " --> forward high performance computing               "
    write(*,*) "                                                      "

  end subroutine display_logo_hpc
!!!=================================================================================================

!!!=================================================================================================
  subroutine display_logo_cdps

    write(*,*) "                                                      "
    write(*,*) "   ________  __  ___  _____     ________  ___  ____   "
    write(*,*) "  / ___/ _ \/ / / / |/ / _ |   / ___/ _ \/ _ \/ __/   "
    write(*,*) " / /__/ , _/ /_/ /    / __ |  / /__/ // / ___/\ \     "
    write(*,*) " \___/_/|_|\____/_/|_/_/ |_|  \___/____/_/  /___/     "
    write(*,*) "                                                      "
    write(*,*) "(Compressible Reactive Unsteady Navier-Stokes Adjoint)"
    write(*,*) " --> complex directivity source point model           "
    write(*,*) "                                                      "

  end subroutine display_logo_cdps
!!!=================================================================================================

!!!=================================================================================================
  subroutine display_date_time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer,dimension(8) :: date_time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call date_and_time(VALUES=date_time)

    write(*,"(A,I4.4,A,I2.2,A,I2.2)") " date: ", date_time(1), "-", date_time(2), "-", date_time(3)
    write(*,"(A,I2.2,A,I2.2,A,I2.2)") " time: ", date_time(5), ":", date_time(6), ":", date_time(7)

  end subroutine display_date_time
!!!=================================================================================================

!!!=================================================================================================
  subroutine display_host_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character*8 host_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call hostnm( host_name )

    write(*,*) "host: ",host_name
  end subroutine display_host_name
!!!=================================================================================================

!!!=================================================================================================
  subroutine display_case_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "case: ",trim(cruna_report_case)
  end subroutine display_case_info
!!!=================================================================================================

!!!=================================================================================================
  subroutine display_login_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character*128 login_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call getlog( login_name )

    write(*,*) "user: ",trim(login_name)
  end subroutine display_login_name
!!!=================================================================================================

!!!=================================================================================================
  subroutine display_equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(*,*) "governing equations: ",trim(equations_name)

  end subroutine display_equations
!!!=================================================================================================

!!!=================================================================================================
  subroutine display_discretisation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(*,*) "discretisation time: ",trim(time_stepper_direct_name)
    write(*,*) "discretisation i1  : ",trim(diff_i1_name)
    write(*,*) "discretisation i2  : ",trim(diff_i2_name)
    write(*,*) "discretisation i3  : ",trim(diff_i3_name)

  end subroutine display_discretisation
!!!=================================================================================================

!!!=================================================================================================
  subroutine display_sponge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    character(len=max_length_parameter) :: sponge_shape_type

    call get_parameter(sponge_shape_type,'sponge.shape',default="quadratic")
    ! originates from subroutine init_boundary_conditions of boundary_conditions.f90
    sponge_shape_type = " ("//trim(sponge_shape_type)//" shape)"

    write(*,*) "sponge (direct) : ",trim(sponge_direct_name)//trim(sponge_shape_type)
    write(*,*) "sponge (adjoint): ",trim(sponge_adjoint_name)//trim(sponge_shape_type)

  end subroutine display_sponge
!!!=================================================================================================

!!!=================================================================================================
  subroutine display_filter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(*,*) "filter i1: "            ,trim(filter_i1_name)
    write(*,*) "filter i2: "            ,trim(filter_i2_name)
    write(*,*) "filter i3: "            ,trim(filter_i3_name)
    write(*,*) "filter i1 (weighted): " ,trim(filter_weighted_i1_name)
    write(*,*) "filter i2 (weighted): " ,trim(filter_weighted_i2_name)
    write(*,*) "filter i3 (weighted): " ,trim(filter_weighted_i3_name)


  end subroutine display_filter
!!!=================================================================================================

!!!=================================================================================================
  subroutine display_compile_info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(*,*) "compile time case: ",trim(cruna_report_case)
    write(*,*) "compile time comp: ",trim(cruna_report_comp)
    write(*,*) "compile time copt: ",trim(cruna_report_copt)
    write(*,*) "compile time ctar: ",trim(cruna_report_ctar)
    write(*,*) "compile time date: ",trim(cruna_report_time)
    write(*,*) "compile time host: ",trim(cruna_report_host)
    write(*,*) "compile time lopt: ",trim(cruna_report_lopt)
    write(*,*) "compile time path: ",trim(cruna_report_path)
    write(*,*) "compile time root: ",trim(cruna_report_root)
    write(*,*) "compile time user: ",trim(cruna_report_user)

  end subroutine display_compile_info
!!!=================================================================================================

!!!=================================================================================================
  subroutine display_parallelism
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(*,*) "parallelism type : ",trim(params%parallelism%type)
    write(*,*) "parallelism cores: ",trim(num2str(params%parallelism%world_size,'(I5.5)'))

  end subroutine display_parallelism
!!!=================================================================================================

end module display
