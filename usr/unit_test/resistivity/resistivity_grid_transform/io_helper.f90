module io_helper

  use helper
  use parallelism
  use parameter

  private

  public get_fname

  character(len=6), save :: fmt_loop        = '(I3.3)'
  character(len=6), save :: fmt_world_block = '(I2.2)'
  character(len=6), save :: fmt_block_image = '(I5.5)'
  character(len=6), save :: fmt_subset      = '(I2.2)'
  character(len=6), save :: fmt_tstep       = '(I6.6)'

contains

!!!=================================================================================================
  subroutine get_fname(fname,type,nt,ns)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_fname),          intent(out)          :: fname
    character(len=*)               ,          intent(in)           :: type
    integer                        ,          intent(in), optional :: nt
    integer                        ,          intent(in), optional :: ns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                        :: tstep                                                   
    integer                                                        :: subset, loop
    character(len=1)               , parameter                     :: separator = "_"
    character(len=128)                                             :: fn_general, fn_loop, fn_block, fn_block_image, fn_subset, fn_timestep
    character(len=max_length_parameter)                            :: fn_wdir
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

    case ("data_adjoint")
       fn_general     = "data_adjoint_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

    case ("data_adjoint_subsets_cache")
       fn_general     = "data_adjoint_subsets_cache"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written       

    case ("data_adjoint_loop")
       call get_parameter(loop,'opt.loop')

       fn_general     = "data_adjoint_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

    case ("data_adjoint_snapshot")
       fn_general     = "data_adjoint_snapshot_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       fn_timestep    = separator // trim(num2str(tstep                         ,fmt_tstep      ))

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

    case ("data_direct_subsets_cache")
       fn_general     = "data_direct_subsets_cache"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written       

    case ("data_direct_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "data_direct_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

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

    case ("alpha_J_list")
       fn_general     = "alpha_J_list_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty 
       ! timestep is empty 

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

    case ("force_backup_subsets_cache")
       fn_general     = "force_backup_subsets_cache"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

    case ("force_grad_subsets_cache")
       fn_general     = "force_grad_subsets_cache"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

    case ("force_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "force_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

    case ("force_subsets_cache")
       fn_general     = "force_subsets_cache"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

    case ("ballon_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "ballon_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset is empty since ballons live in frequency-domain (all freq.-components are written)
       ! timestep is empty since ballons live in frequency-domain (all freq.-components are written)

    case ("drive_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "drive_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset is empty since driving fct. is frequency dependent (all freq.-components are written)
       ! timestep is empty since driving fct. is frequency dependent

    case ("geometry")
       fn_general     = "geometry_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant

     case ("geometry_calc_space")
       fn_general     = "geometry_calc_space_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant   

    case ("geometry_indices")
       fn_general     = "geometry_indices_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant
       
    case ("geometry_metrics")
       fn_general     = "geometry_metric_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant
       
    case ("geometry_jacobian")
       fn_general     = "geometry_jacobian_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant


    case ("gradient_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "gradient_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

    case ("ballon_grad_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "ballon_grad_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset is empty since ballons live in frequency-domain (all freq.-components are written)
       ! timestep is empty since ballons live in frequency-domain (all freq.-components are written)

    case ("drive_grad_loop")
       call get_parameter(loop,'opt.loop')
       fn_general     = "drive_grad_"
       fn_loop        = separator // "loop_" // trim(num2str(loop,fmt_loop))
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset is empty since driving fct. is frequency dependent (all freq.-components are written)
       ! timestep is empty since driving fct. is frequency dependent

    case ("penalization")
       fn_general     = "penalization_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty
       ! timestep is empty as full subsets will be written       

    case ("objective")
       fn_general     = "objective"
       ! world block is empty as J is global
       ! block image is empty as J is global
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant

    case ("objective_ascii")
       fn_general     = "objective.dat"
       ! world block is empty as J is global
       ! block image is empty as J is global
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant

    case ("opt_target")
       fn_general     = "opt_target_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as full subsets will be written

    case ("opt_target_all")
       fn_general     = "opt_target_all"
       ! world block is empty
       ! block image is empty
       ! subset   is empty
       ! timestep is empty

    case ("sample")
       fn_general     = "sample_"
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

    case ("sponge_function")
       fn_general     = "sponge_function_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant

    case ("sigma_x_function")
       fn_general     = "sigma_x_function_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       ! subset   is empty as geometry is contant
       ! timestep is empty as geometry is contant

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
       ! subset   is empty as theta_x is contant
       ! timestep is empty as theta_x is contant

    case ("theta_t_function")
       fn_general     = "theta_t_function_"
       fn_block       = separator // trim(num2str(params%parallelism%world_block,fmt_world_block))
       fn_block_image = separator // trim(num2str(params%parallelism%block_image,fmt_block_image))
       fn_subset      = separator // trim(num2str(subset                        ,fmt_subset     ))
       ! timestep is empty as theta_t is function of time

    case default
       write(*,*) "error in helper.f90:get_fname:type '",trim(type),"' not found, can not create file name"
       write(*,*) ">> STOP <<"
       stop

    end select

    ! add working path if present
    call get_parameter(fn_wdir,'io.workdirectory',default = ".")
    fn_wdir = trim(fn_wdir) // "/"

    ! assemble the file name
    fname = trim(fn_wdir)     // &
         trim(fn_general)     // &
         trim(fn_block)       // &
         trim(fn_block_image) // &
         trim(fn_subset)      // &
         trim(fn_timestep)    // &
         trim(fn_loop)       

  end subroutine get_fname
!!!=================================================================================================

end module io_helper
