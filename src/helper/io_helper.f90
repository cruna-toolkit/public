!> Utility module for project-specific file-name generation and file inquiries.
!!
!! The module centralizes construction of file names used by the application.
!! File names are assembled from a type-dependent base name and optional
!! identifiers such as optimization loop, MPI world block, block image, subset,
!! and timestep. The module also provides a convenience routine for checking
!! whether a file exists.
!!
!! Naming details depend on the global parameter state imported from `parameter`
!! and on helper functionality imported from `helper`.
!!
!! Only `get_fname` and `inquire_file` are part of the public module interface.

module io_helper

  use helper
  use parallelism
  use parameter

  private

  public get_fname
  public inquire_file
  public pre_read

  !> Integer format used for optimization-loop identifiers (`I4.4`).
  character(len=6), save :: fmt_loop        = '(I4.4)'
  !> Integer format used for MPI world-block identifiers (`I2.2`).
  character(len=6), save :: fmt_world_block = '(I2.2)'
  !> Integer format used for block-image identifiers (`I5.5`).
  character(len=6), save :: fmt_block_image = '(I5.5)'
  !> Integer format used for subset identifiers (`I2.2`).
  character(len=6), save :: fmt_subset      = '(I2.2)'
  !> Integer format used for timestep identifiers (`I7.7`).
  character(len=6), save :: fmt_tstep       = '(I7.7)'

contains

!!!=================================================================================================
  !> Construct the file name for a specified application data type.
  !!
  !! The routine generates a project-specific file name according to `type`.
  !! Depending on the selected type, the resulting name may contain identifiers
  !! for the optimization loop, MPI world block, block image, subset, and/or
  !! timestep. Components that do not apply to a particular file type remain
  !! empty and are therefore omitted from the final name.
  !!
  !! The general non-MPI-I/O naming scheme is:
  !!
  !!     <workdirectory>/<general><block><block_image><subset><timestep><loop>
  !!
  !! For effective MPI-I/O, the block-image component is omitted.
  !!
  !! If `nt` is not present, `params%time%nt` supplies the timestep. If `ns` is
  !! not present, `params%time%ns` supplies the subset index. File types that
  !! depend on the current optimization loop obtain it from parameter `opt.loop`.
  !!
  !! MPI-I/O is enabled only for file types that request it and only when the
  !! global switch `params%io%mpiio` is also enabled. The optional `mpiio`
  !! argument returns this effective state.
  !!
  !! The working directory is read from parameter `io.workdirectory` and defaults
  !! to `.`. A slash is appended before the file-name components are assembled.
  !!
  !! @param[out] fname Generated file name including the working directory.
  !! @param[in] type File-type identifier selecting the naming convention.
  !! @param[in] nt Optional timestep. Defaults to `params%time%nt`.
  !! @param[in] ns Optional subset index. Defaults to `params%time%ns`.
  !! @param[out] mpiio Optional flag indicating whether MPI-I/O is effectively enabled.
  !!
  !! @warning An unsupported `type` prints an error message and terminates execution
  !! with `stop`.
  subroutine get_fname(fname,type,nt,ns,mpiio)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_fname),          intent(out)           :: fname
    character(len=*)               ,          intent(in)            :: type
    integer                        ,          intent(in) , optional :: nt
    integer                        ,          intent(in) , optional :: ns
    logical                        ,          intent(out), optional :: mpiio
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                         :: tstep
    integer                                                         :: subset, loop
    character(len=1)               , parameter                      :: separator = "_"
    character(len=128)                                              :: fn_general, fn_loop, fn_block, fn_block_image, fn_subset, fn_timestep
    character(len=max_length_parameter)                             :: fn_wdir
    logical                                                         :: mpi_io
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! handle optional arguments
    if(present(nt))then
       tstep = nt
    else
       tstep = params%time%nt
    endif

    if(present(ns))then
       subset = ns
    else
       subset = params%time%ns
    endif

    ! no mpi write by defaul, only for selected file types, see below
    mpi_io = .false.

    ! set all string empty
    fn_general     = ""
    fn_loop        = ""
    fn_block       = ""
    fn_block_image = ""
    fn_subset      = ""
    fn_timestep    = ""
    fn_wdir        = ""

    ! select type
    !    defines the general file name and more ...
    select case (trim(type))

    case ("alpha_J_list")
       fn_general     = "alpha_J_list_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty
       ! timestep is empty

    case ("ballon_grad_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "ballon_grad_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset is empty since ballons live in frequency-domain (all freq.-components are written)
       ! timestep is empty since ballons live in frequency-domain (all freq.-components are written)

    case ("ballon_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "ballon_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset is empty since ballons live in frequency-domain (all freq.-components are written)
       ! timestep is empty since ballons live in frequency-domain (all freq.-components are written)

    case ("data_adjoint")
       fn_general     = "data_adjoint_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("data_adjoint_subsets_cache")
       fn_general     = "data_adjoint_subsets_cache"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("data_adjoint_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "data_adjoint_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("data_adjoint_snapshot")
       fn_general     = "data_adjoint_snapshot_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       fn_timestep    = separator // trim(num2str(tstep                         ,fmt_tstep      ))

       mpi_io = .true.

    case ("data_adjoint_snapshot_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "data_adjoint_snapshot_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       fn_timestep    = separator // trim(num2str(tstep                         ,fmt_tstep      ))

       mpi_io = .true.

    case ("data_cdps")
       fn_general     = "data_cdps_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty
       ! timestep is empty as full subsets will be written


    case ("data_direct")
       fn_general     = "data_direct_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("data_direct_subsets_cache")
       fn_general     = "data_direct_subsets_cache_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("data_direct_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "data_direct_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("pxy_direct_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "pxy_direct_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

    case ("data_direct_snapshot")
       fn_general     = "data_direct_snapshot_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       fn_timestep    = separator // trim(num2str(tstep                         ,fmt_tstep      ))

       mpi_io = .true.

    case ("data_direct_snapshot_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "data_direct_snapshot_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       fn_timestep    = separator // trim(num2str(tstep                         ,fmt_tstep      ))

       mpi_io = .true.

    case ("debug_field")
       fn_general     = "debug_field_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty
       ! timestep is empty

    case ("debug_field_A")
       fn_general     = "debug_field_A_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       fn_timestep    = separator // trim(num2str(tstep                         ,fmt_tstep      ))

    case ("debug_field_B")
       fn_general     = "debug_field_B_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       fn_timestep    = separator // trim(num2str(tstep                         ,fmt_tstep      ))

    case ("debug_field_C")
       fn_general     = "debug_field_C_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       fn_timestep    = separator // trim(num2str(tstep                         ,fmt_tstep      ))

    case ("debug_field_D")
       fn_general     = "debug_field_D_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

    case ("drive_grad_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "drive_grad_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset is empty since driving fct. is frequency dependent (all freq.-components are written)
       ! timestep is empty since driving fct. is frequency dependent

    case ("drive_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "drive_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset is empty since driving fct. is frequency dependent (all freq.-components are written)
       ! timestep is empty since driving fct. is frequency dependent

    case ("force_backup_subsets_cache")
       fn_general     = "force_backup_subsets_cache_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("force_grad_subsets_cache")
       fn_general     = "force_grad_subsets_cache_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("force_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "force_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("force_subsets_cache")
       fn_general     = "force_subsets_cache_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("geometry")
       fn_general     = "geometry_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant

       mpi_io = .true.

    case ("geometry_calc_space")
       fn_general     = "geometry_calc_space_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant

       mpi_io = .true.

    case ("geometry_indices")
       fn_general     = "geometry_indices_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant

       mpi_io = .true.

    case ("geometry_jacobian")
       fn_general     = "geometry_jacobian_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant

       mpi_io = .true.

    case ("geometry_metrics")
       fn_general     = "geometry_metric_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant

       mpi_io = .true.

    case ("gradient_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "gradient_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("gradient_force_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "gradient_force_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("gradient_penalization_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "gradient_penalization_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is contant
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("init_data_direct_snapshot")
       fn_general     = "init_data_direct_snapshot_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       fn_timestep    = separator // trim(num2str(tstep                         ,fmt_tstep      ))

       mpi_io = .true.

    case ("J_shift_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "J_shift_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       ! world block is empty as J is global
       ! block image is empty as J is global
       ! subset   is empty as ts_shift_function_loop is function of time
       ! timestep is empty as ts_shift_function_loop is function of time

    case ("J_shift_sample_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "J_shift_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       ! world block is empty as J is global
       ! block image is empty as J is global
       ! subset   is empty as ts_shift_function_loop is function of time
       ! timestep is empty as ts_shift_function_loop is function of time

    case ("kvector")
       fn_general     = "kvector_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is constant
       ! timestep is empty as geometry is constant       

    case ("noise")
       fn_general     = "noise_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

    case ("objective")
       fn_general     = "objective"
       ! world block is empty as J is global
       ! block image is empty as J is global
       ! subset   is empty as objective is contant
       ! timestep is empty as objective is contant

    case ("objective_ascii")
       fn_general     = "objective.dat"
       ! world block is empty as J is global
       ! block image is empty as J is global
       ! subset   is empty as objective_ascii is contant
       ! timestep is empty as objective_ascii is contant

    case ("opt_target")
       fn_general     = "opt_target_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("opt_target_all")
       fn_general     = "opt_target_all"
       ! world block is empty
       ! block image is empty
       ! subset   is empty
       ! timestep is empty

       mpi_io = .true.

    case ("penalization")
       fn_general     = "penalization_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("penalization_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "penalization_loop_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! timestep is empty as full subsets will be written

       mpi_io = .true.

    case ("penalization_time")
       fn_general     = "penalization_time_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       fn_timestep    = separator // trim(num2str(tstep                         ,fmt_tstep      ))

       mpi_io = .true.

    case ("ref_data_direct_snapshot")
       fn_general     = "ref_data_direct_snapshot_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       fn_timestep    = separator // trim(num2str(tstep                         ,fmt_tstep      ))

       mpi_io = .true.

    case ("restartfile")
       fn_general     = "restartfile"
       ! fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       ! fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty
       ! timestep is empty as full subsets will be written

    case ("restartfile_loop")
       fn_general     = "restartfile_loop"
       ! fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       ! fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty
       ! timestep is empty as full subsets will be written

    case ("sample")
       fn_general     = "sample_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant
       
   case ("sample_loop")
       fn_general     = "sample_loop_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant

    case ("sample_time")
       fn_general     = "sample_time_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       fn_timestep    = separator // trim(num2str(tstep                         ,fmt_tstep      ))
       
    case ("speaker_positions_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "speaker_positions_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty
       ! timestep is empty

    case ("sponge_shape_function")
       fn_general     = "sponge_shape_function_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is constant
       ! timestep is empty as geometry is constant

       mpi_io = .true.

    case ("sigma_x_function")
       fn_general     = "sigma_x_function_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is constant
       ! timestep is empty as geometry is constant

       mpi_io = .true.

    case ("sigma_t_function")
       fn_general     = "sigma_t_function_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as sigma_t is function of time

    case ("theta_x_function")
       fn_general     = "theta_x_function_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as theta_x is constant
       ! timestep is empty as theta_x is constant

       mpi_io = .true.

    case ("theta_x_function_penalization")
       fn_general     = "theta_x_function_penalization_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as theta_x is constant
       ! timestep is empty as theta_x is constant

       mpi_io = .true.

    case ("theta_t_function")
       fn_general     = "theta_t_function_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as theta_t is function of time

    case default
       write(*,*) "error in io_helper.f90:get_fname:type '",trim(type),"' not found, can not create file name"
       write(*,*) ">> STOP <<"
       stop

    end select

    ! globally enable/disable parallel HDF5 I/O, default = .false.
    mpi_io = (mpi_io).and.(params%io%mpiio)

    ! add working path if present
    call get_parameter(fn_wdir,'io.workdirectory',default = ".")
    fn_wdir = trim(fn_wdir) // "/"

    if (mpi_io.eqv..false.) then
       ! assemble the file name
       fname = trim(fn_wdir)     // &
            trim(fn_general)     // &
            trim(fn_block)       // &
            trim(fn_block_image) // &
            trim(fn_subset)      // &
            trim(fn_timestep)    // &
            trim(fn_loop)
    else
       ! assemble the file name
       fname = trim(fn_wdir)     // &
            trim(fn_general)     // &
            trim(fn_block)       // &
                                ! trim(fn_block_image) // &
            trim(fn_subset)      // &
            trim(fn_timestep)    // &
            trim(fn_loop)
    end if

    ! handle optional arguments
    if(present(mpiio))then
       mpiio = mpi_io
    end if

  end subroutine get_fname
!!!=================================================================================================

!!!=================================================================================================
  !> Test whether a file exists.
  !!
  !! This routine wraps the Fortran `inquire` statement and checks the file named
  !! by `fname`. Leading and trailing blanks are removed from the file name before
  !! the inquiry is performed.
  !!
  !! @param[in] fname Name or path of the file to test.
  !! @param[out] file_exists `.true.` if the file exists, otherwise `.false.`.
  subroutine inquire_file(fname,file_exists)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*)               ,          intent(in)          :: fname
    logical                        ,          intent(out)         :: file_exists
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    inquire(file=trim(fname), exist=file_exists )

  end subroutine inquire_file
!!!=================================================================================================

!!!=================================================================================================
  !> Pre-read a file into the Linux page cache asynchronously.
  !!
  !! Pre-reading can be enabled globally using parameter `io.prefetch`. The file name is
  !! generated using `get_fname`. Files using MPI-I/O are not pre-read. For non-MPI-I/O
  !! files, each MPI process pre-reads its rank-local file independently.
  !!
  !! On the first invocation, MPI world rank zero checks whether `vmtouch` is available.
  !! If unavailable, `cat` is used instead. Errors during pre-reading are silently ignored.
  !!
  !! **Arguments**
  !!
  !! * `type` - File-type identifier passed to `get_fname`.
  !! * `nt` - Optional timestep passed to `get_fname`.
  !! * `ns` - Optional subset index passed to `get_fname`.
  subroutine pre_read(type,nt,ns)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*), intent(in)           :: type
    integer,          intent(in), optional :: nt
    integer,          intent(in), optional :: ns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_fname) :: fname
    integer                         :: exitstat,cmdstat
    integer                         :: use_vmtouch_int
    logical                         :: mpiio
    logical, save                   :: initialized = .false.
    logical, save                   :: prefetch    = .false.
    logical, save                   :: use_vmtouch = .false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Initialize pre-fetching and test once whether vmtouch is available.
    if (.not.initialized) then

       call get_parameter(prefetch,'io.prefetch',default = .false.)

       if (prefetch) then

          use_vmtouch_int = 0

          if (params%parallelism%world_image == 1) then

             exitstat = 1

             call execute_command_line("command -v vmtouch >/dev/null 2>&1", &
                  wait=.true.,exitstat=exitstat)

             if (exitstat == 0) then
                use_vmtouch_int = 1
                write(*,'(A)') "pre_read: vmtouch found"
             else
                write(*,'(A)') "pre_read: vmtouch not found, using cat"
             end if

          end if

          call bcast(use_vmtouch_int,0,params%parallelism%world_comm)

          use_vmtouch = (use_vmtouch_int == 1)

       end if

       initialized = .true.

    end if

    ! Return if pre-fetching is disabled globally.
    if (.not.prefetch) return

    ! Generate rank-local file name and determine effective MPI-I/O state.
    call get_fname(fname,type,nt,ns,mpiio)

    ! Do not pre-read shared MPI-I/O files.
    if (mpiio) return

    ! Launch asynchronous rank-local pre-read and silently ignore all errors.
    if (use_vmtouch) then
       call command('"' // trim(fname) // '" >/dev/null 2>&1',2)
    else
       call command('"' // trim(fname) // '" >/dev/null 2>&1',3)
    end if

  end subroutine pre_read
!!!=================================================================================================


end module io_helper
