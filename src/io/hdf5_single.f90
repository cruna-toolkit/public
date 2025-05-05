module hdf5_single

  USE hdf5
  use helper
  use io_helper
  use parameter

  private

  public get_dset_dim
  public load
  public load_sub
  public store

  interface load
     module procedure load_5d_rk
     module procedure load_4d_rk
     module procedure load_4d_i
     module procedure load_3d_rk
     module procedure load_2d_rk
     module procedure load_1d_rk
     module procedure load_1d_i
  end interface load

  interface load_sub
     module procedure load_sub_4d_rk
  end interface load_sub

  interface store
     module procedure store_1d_rk
     module procedure store_1d_i
     module procedure store_2d_rk
     module procedure store_3d_rk
     module procedure store_4d_i
     module procedure store_4d_rk
     module procedure store_5d_rk
  end interface store

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! s t o r e !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!=================================================================================================
  subroutine store_1d_rk(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin, file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , dimension(:)      , intent(in)                     :: data
    character(len=*),                     intent(in)                     :: type
    character(len=*),                     intent(in), optional           :: fname_optin
    integer         ,                     intent(in), optional           :: nt_optin
    integer         ,                     intent(in), optional           :: ns_optin
    character(len=*),                     intent(in), optional           :: dset_name_optin
    logical         ,                     intent(in), optional           :: file_overwrite_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:), allocatable                          :: dims
    integer(hid_t)                                                       :: file_id, c_type, dset_id, dspace_id 
    integer(size_t)                                                      :: s_dim 
    integer                                                              :: iostat
    logical                                                              :: file_exist, write_params = .true., file_overwrite = .false.

    character(len=max_length_parameter*5), dimension(max_nbr_parameters) :: parameter_list
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)

       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5") 
    endif

    if(present(file_overwrite_optin))then
       file_overwrite = file_overwrite_optin
    end if

    ! open hdf5
    call h5open_f(iostat)

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)

       ! inquire file
       inquire(file = trim(fname),exist = file_exist, iostat = iostat)

       ! open exisiting file or create new one
       if((file_exist.eqv..true.).and.(file_overwrite.eqv..false.)) then
          call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)
          ! do not write params (again)
          write_params = .false.
       else
          call h5fcreate_f_error_control(trim(fname), h5f_acc_trunc_f, file_id, iostat)
       end if
    else ! the defaul way
       ! default dset name = 'data'
       dset_name = 'data'

       ! create new h5 file, deletes ALL old stuff, fastest way of writing, the way to go!
       ! params are written as write_params = .true. by default
       call h5fcreate_f_error_control(trim(fname), h5f_acc_trunc_f, file_id, iostat)
    end if

    ! - data - 

    ! determine size
    allocate(dims(1))
    dims(1) = size(data,1)
    ! create dataspace
    call h5screate_simple_f(1, dims, dspace_id, iostat)
    ! create dataset
    call h5dcreate_f_error_control(file_id, trim(dset_name),  h5t_native_double, dspace_id, dset_id, iostat)
    ! write dataset
    call h5dwrite_f(dset_id, h5t_native_double, data, dims, iostat)
    ! close dataset and dataspace
    call h5dclose_f(dset_id,iostat)
    call h5sclose_f(dspace_id,iostat)
    deallocate(dims)

    ! - parameter - 
    if(write_params.eqv..true.) then
       ! convert parameter to char
       call create_parameter_list(parameter_list)
       ! determine sizes
       s_dim = max_length_parameter*5 ! chars per row
       allocate(dims(1))
       dims  = max_nbr_parameters     ! row number
       ! create datatype (for entire row)
       call h5tcopy_f(h5t_fortran_s1,c_type,iostat)
       call h5tset_size_f(c_type,s_dim,iostat)
       ! create dataspace
       call h5screate_simple_f(1,dims,dspace_id,iostat)
       ! create dataset
       call h5dcreate_f_error_control(file_id,"params",c_type,dspace_id,dset_id,iostat)
       ! write dataset
       call h5dwrite_f(dset_id,c_type,parameter_list,dims,iostat)
       ! close datatype, dataset and dataspace
       call h5tclose_f(c_type,iostat)
       call h5dclose_f(dset_id,iostat)
       call h5sclose_f(dspace_id,iostat)

       deallocate(dims)
    end if

    ! close hdf5 file
    call h5fclose_f(file_id, iostat)
    ! close hdf5
    call h5close_f(iostat)

  end subroutine store_1d_rk
!!!=================================================================================================

!!!=================================================================================================
  subroutine store_1d_i(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin, file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer         , dimension(:)      , intent(in)                     :: data
    character(len=*),                     intent(in)                     :: type
    character(len=*),                     intent(in), optional           :: fname_optin
    integer         ,                     intent(in), optional           :: nt_optin
    integer         ,                     intent(in), optional           :: ns_optin
    character(len=*),                     intent(in), optional           :: dset_name_optin
    logical         ,                     intent(in), optional           :: file_overwrite_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:), allocatable                          :: dims
    integer(hid_t)                                                       :: file_id, c_type, dset_id, dspace_id 
    integer(size_t)                                                      :: s_dim 
    integer                                                              :: iostat
    logical                                                              :: file_exist, write_params = .true., file_overwrite = .false.

    character(len=max_length_parameter*5), dimension(max_nbr_parameters) :: parameter_list
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)

       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5") 
    endif

    if(present(file_overwrite_optin))then
       file_overwrite = file_overwrite_optin
    end if

    ! open hdf5
    call h5open_f(iostat)

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)

       ! inquire file
       inquire(file = trim(fname),exist = file_exist, iostat = iostat)

       ! open exisiting file or create new one
       if((file_exist.eqv..true.).and.(file_overwrite.eqv..false.)) then
          call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)
          ! do not write params (again)
          write_params = .false.
       else
          call h5fcreate_f_error_control(trim(fname), h5f_acc_trunc_f, file_id, iostat)
       end if
    else ! the defaul way
       ! default dset name = 'data'
       dset_name = 'data'

       ! create new h5 file, deletes ALL old stuff, fastest way of writing, the way to go!
       ! params are written as write_params = .true. by default
       call h5fcreate_f_error_control(trim(fname), h5f_acc_trunc_f, file_id, iostat)
    end if

    ! - data - 

    ! determine size
    allocate(dims(1))
    dims(1) = size(data,1)
    ! create dataspace
    call h5screate_simple_f(1, dims, dspace_id, iostat)
    ! create dataset
    ! call h5dcreate_f_error_control(file_id, trim(dset_name),  h5t_native_double, dspace_id, dset_id, iostat)
    call h5dcreate_f_error_control(file_id, trim(dset_name),  h5t_native_integer, dspace_id, dset_id, iostat)
    ! write dataset
    call h5dwrite_f(dset_id, h5t_native_integer, data, dims, iostat)
    ! close dataset and dataspace
    call h5dclose_f(dset_id,iostat)
    call h5sclose_f(dspace_id,iostat)
    deallocate(dims)

    ! - parameter - 
    if(write_params.eqv..true.) then
       ! convert parameter to char
       call create_parameter_list(parameter_list)
       ! determine sizes
       s_dim = max_length_parameter*5 ! chars per row
       allocate(dims(1))
       dims  = max_nbr_parameters     ! row number
       ! create datatype (for entire row)
       call h5tcopy_f(h5t_fortran_s1,c_type,iostat)
       call h5tset_size_f(c_type,s_dim,iostat)
       ! create dataspace
       call h5screate_simple_f(1,dims,dspace_id,iostat)
       ! create dataset
       call h5dcreate_f_error_control(file_id,"params",c_type,dspace_id,dset_id,iostat)
       ! write dataset
       call h5dwrite_f(dset_id,c_type,parameter_list,dims,iostat)
       ! close datatype, dataset and dataspace
       call h5tclose_f(c_type,iostat)
       call h5dclose_f(dset_id,iostat)
       call h5sclose_f(dspace_id,iostat)

       deallocate(dims)
    end if

    ! close hdf5 file
    call h5fclose_f(file_id, iostat)
    ! close hdf5
    call h5close_f(iostat)

  end subroutine store_1d_i
!!!=================================================================================================

!!!=================================================================================================
  subroutine store_2d_rk(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin, file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , dimension(:,:)    , intent(in)                     :: data
    character(len=*),                     intent(in)                     :: type
    character(len=*),                     intent(in), optional           :: fname_optin
    integer         ,                     intent(in), optional           :: nt_optin
    integer         ,                     intent(in), optional           :: ns_optin
    character(len=*),                     intent(in), optional           :: dset_name_optin
    logical         ,                     intent(in), optional           :: file_overwrite_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:), allocatable                          :: dims
    integer(hid_t)                                                       :: file_id, c_type, dset_id, dspace_id 
    integer(size_t)                                                      :: s_dim 
    integer                                                              :: iostat
    logical                                                              :: file_exist, write_params = .true., file_overwrite = .false.

    character(len=max_length_parameter*5), dimension(max_nbr_parameters) :: parameter_list
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)

       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5")
    endif

    if(present(file_overwrite_optin))then
       file_overwrite = file_overwrite_optin
    end if

    ! open hdf5
    call h5open_f(iostat)

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)

       ! inquire file
       inquire(file = trim(fname),exist = file_exist, iostat = iostat)

       ! open exisiting file or create new one
       if((file_exist.eqv..true.).and.(file_overwrite.eqv..false.)) then
          call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)
          ! do not write params (again)
          write_params = .false.
       else
          call h5fcreate_f_error_control(trim(fname), h5f_acc_trunc_f, file_id, iostat)
       end if
    else ! the defaul way
       ! default dset name = 'data'
       dset_name = 'data'

       ! create new h5 file, deletes ALL old stuff, fastest way of writing, the way to go!
       ! params are written as write_params = .true. by default
       call h5fcreate_f_error_control(trim(fname), h5f_acc_trunc_f, file_id, iostat)
    end if

    ! - data - 
    allocate(dims(1:2))
    dims(1) = size(data,1)
    dims(2) = size(data,2)

    call h5screate_simple_f(2, dims, dspace_id, iostat)
    call h5dcreate_f_error_control(file_id, trim(dset_name),  h5t_native_double, dspace_id, dset_id, iostat)
    call h5dwrite_f(dset_id, h5t_native_double, data, dims, iostat)

    ! close dataset and dataspace
    call h5dclose_f(dset_id,iostat)
    call h5sclose_f(dspace_id,iostat)

    deallocate(dims)

    ! - parameter -
    if(write_params.eqv..true.) then
       call create_parameter_list(parameter_list)
       s_dim = max_length_parameter*5

       allocate(dims(1))
       dims  = max_nbr_parameters

       call h5tcopy_f(h5t_fortran_s1,c_type,iostat)
       call h5tset_size_f(c_type,s_dim,iostat)
       call h5screate_simple_f(1,dims,dspace_id,iostat)
       call h5dcreate_f_error_control(file_id,"params",c_type,dspace_id,dset_id,iostat)
       call h5dwrite_f(dset_id,c_type,parameter_list,dims,iostat)
       call h5tclose_f(c_type,iostat)
       call h5dclose_f(dset_id,iostat)
       call h5sclose_f(dspace_id,iostat)

       deallocate(dims)
    end if

    call h5fclose_f(file_id, iostat)
    call h5close_f(iostat)

  end subroutine store_2d_rk
!!!=================================================================================================

!!!=================================================================================================
  subroutine store_3d_rk(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin, file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , dimension(:,:,:)  , intent(in)                     :: data
    character(len=*),                     intent(in)                     :: type
    character(len=*),                     intent(in), optional           :: fname_optin
    integer         ,                     intent(in), optional           :: nt_optin
    integer         ,                     intent(in), optional           :: ns_optin
    character(len=*),                     intent(in), optional           :: dset_name_optin
    logical         ,                     intent(in), optional           :: file_overwrite_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:), allocatable                          :: dims
    integer(hid_t)                                                       :: file_id, c_type, dset_id, dspace_id 
    integer(size_t)                                                      :: s_dim 
    integer                                                              :: iostat
    logical                                                              :: file_exist, write_params = .true., file_overwrite = .false.

    character(len=max_length_parameter*5), dimension(max_nbr_parameters) :: parameter_list
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)

       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5")
    endif

    if(present(file_overwrite_optin))then
       file_overwrite = file_overwrite_optin
    end if

    ! open hdf5
    call h5open_f(iostat)

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)

       ! inquire file
       inquire(file = trim(fname),exist = file_exist, iostat = iostat)

       ! open exisiting file or create new one
       if((file_exist.eqv..true.).and.(file_overwrite.eqv..false.)) then
          call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)
          ! do not write params (again)
          write_params = .false.
       else
          call h5fcreate_f_error_control(trim(fname), h5f_acc_trunc_f, file_id, iostat)
       end if
    else ! the defaul way
       ! default dset name = 'data'
       dset_name = 'data'

       ! create new h5 file, deletes ALL old stuff, fastest way of writing, the way to go!
       ! params are written as write_params = .true. by default
       call h5fcreate_f_error_control(trim(fname), h5f_acc_trunc_f, file_id, iostat)
    end if

    ! - data - 
    allocate(dims(1:3))
    dims(1) = size(data,1)
    dims(2) = size(data,2)
    dims(3) = size(data,3)

    call h5screate_simple_f(3, dims, dspace_id, iostat)
    call h5dcreate_f_error_control(file_id, trim(dset_name),  h5t_native_double, dspace_id, dset_id, iostat)
    call h5dwrite_f(dset_id, h5t_native_double, data, dims, iostat)
    call h5dclose_f(dset_id,iostat)
    call h5sclose_f(dspace_id,iostat)

    deallocate(dims)

    ! - parameter -
    if(write_params.eqv..true.) then
       call create_parameter_list(parameter_list)
       s_dim = max_length_parameter*5

       allocate(dims(1))
       dims  = max_nbr_parameters

       call h5tcopy_f(h5t_fortran_s1,c_type,iostat)
       call h5tset_size_f(c_type,s_dim,iostat)
       call h5screate_simple_f(1,dims,dspace_id,iostat)
       call h5dcreate_f_error_control(file_id,"params",c_type,dspace_id,dset_id,iostat)
       call h5dwrite_f(dset_id,c_type,parameter_list,dims,iostat)
       call h5tclose_f(c_type,iostat)
       call h5dclose_f(dset_id,iostat)
       call h5sclose_f(dspace_id,iostat)

       deallocate(dims)
    end if

    call h5fclose_f(file_id, iostat)
    call h5close_f(iostat)

  end subroutine store_3d_rk
!!!=================================================================================================


!!!=================================================================================================
  subroutine store_4d_i(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin, file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer         , dimension(:,:,:,:), intent(in)                     :: data
    character(len=*),                     intent(in)                     :: type
    character(len=*),                     intent(in), optional           :: fname_optin
    integer         ,                     intent(in), optional           :: nt_optin
    integer         ,                     intent(in), optional           :: ns_optin
    character(len=*),                     intent(in), optional           :: dset_name_optin
    logical         ,                     intent(in), optional           :: file_overwrite_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:), allocatable                          :: dims
    integer(hid_t)                                                       :: file_id, c_type, dset_id, dspace_id 
    integer(size_t)                                                      :: s_dim 
    integer                                                              :: iostat
    logical                                                              :: file_exist, write_params = .true., file_overwrite = .false.

    character(len=max_length_parameter*5), dimension(max_nbr_parameters) :: parameter_list
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)
       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5")
    endif

    if(present(file_overwrite_optin))then
       file_overwrite = file_overwrite_optin
    end if

    ! open hdf5
    call h5open_f(iostat)

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)

       ! inquire file
       inquire(file = trim(fname),exist = file_exist, iostat = iostat)

       ! open exisiting file or create new one
       if((file_exist.eqv..true.).and.(file_overwrite.eqv..false.)) then
          call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)
          ! do not write params (again)
          write_params = .false.
       else
          call h5fcreate_f_error_control(trim(fname), h5f_acc_trunc_f, file_id, iostat)
       end if
    else ! the defaul way
       ! default dset name = 'data'
       dset_name = 'data'

       ! create new h5 file, deletes ALL old stuff, fastest way of writing, the way to go!
       ! params are written as write_params = .true. by default
       call h5fcreate_f_error_control(trim(fname), h5f_acc_trunc_f, file_id, iostat)
    end if

    ! - data - 
    allocate(dims(1:4))
    dims(1) = size(data,1)
    dims(2) = size(data,2)
    dims(3) = size(data,3)
    dims(4) = size(data,4)

    call h5screate_simple_f(4, dims, dspace_id, iostat)
    call h5dcreate_f_error_control(file_id, trim(dset_name),  h5t_native_integer, dspace_id, dset_id, iostat)
    call h5dwrite_f(dset_id, h5t_native_integer, data, dims, iostat)
    call h5dclose_f(dset_id,iostat)
    call h5sclose_f(dspace_id,iostat)

    deallocate(dims)

    ! - parameter -
    if(write_params.eqv..true.) then
       call create_parameter_list(parameter_list)
       s_dim = max_length_parameter*5

       allocate(dims(1))
       dims  = max_nbr_parameters

       call h5tcopy_f(h5t_fortran_s1,c_type,iostat)
       call h5tset_size_f(c_type,s_dim,iostat)
       call h5screate_simple_f(1,dims,dspace_id,iostat)
       call h5dcreate_f_error_control(file_id,"params",c_type,dspace_id,dset_id,iostat)
       call h5dwrite_f(dset_id,c_type,parameter_list,dims,iostat)
       call h5tclose_f(c_type,iostat)
       call h5dclose_f(dset_id,iostat)
       call h5sclose_f(dspace_id,iostat)

       deallocate(dims)
    end if

    call h5fclose_f(file_id, iostat)
    call h5close_f(iostat)

  end subroutine store_4d_i
!!!=================================================================================================

!!!=================================================================================================
  subroutine store_4d_rk(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin, file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , dimension(:,:,:,:), intent(in)                     :: data
    character(len=*),                     intent(in)                     :: type
    character(len=*),                     intent(in), optional           :: fname_optin
    integer         ,                     intent(in), optional           :: nt_optin
    integer         ,                     intent(in), optional           :: ns_optin
    character(len=*),                     intent(in), optional           :: dset_name_optin
    logical         ,                     intent(in), optional           :: file_overwrite_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:), allocatable                          :: dims
    integer(hid_t)                                                       :: file_id, c_type, dset_id, dspace_id 
    integer(size_t)                                                      :: s_dim 
    integer                                                              :: iostat
    logical                                                              :: file_exist, write_params = .true., file_overwrite = .false.

    character(len=max_length_parameter*5), dimension(max_nbr_parameters) :: parameter_list
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)

       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5")
    endif

    if(present(file_overwrite_optin))then
       file_overwrite = file_overwrite_optin
    end if

    ! open hdf5
    call h5open_f(iostat)

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)

       ! inquire file
       inquire(file = trim(fname),exist = file_exist, iostat = iostat)

       ! open exisiting file or create new one
       if((file_exist.eqv..true.).and.(file_overwrite.eqv..false.)) then
          call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)
          ! do not write params (again)
          write_params = .false.
       else
          call h5fcreate_f_error_control(trim(fname), h5f_acc_trunc_f, file_id, iostat)
       end if
    else ! the defaul way
       ! default dset name = 'data'
       dset_name = 'data'

       ! create new h5 file, deletes ALL old stuff, fastest way of writing, the way to go!
       ! params are written as write_params = .true. by default
       call h5fcreate_f_error_control(trim(fname), h5f_acc_trunc_f, file_id, iostat)
    end if

    ! - data - 
    allocate(dims(1:4))
    dims(1) = size(data,1)
    dims(2) = size(data,2)
    dims(3) = size(data,3)
    dims(4) = size(data,4)

    call h5screate_simple_f(4, dims, dspace_id, iostat)
    call h5dcreate_f_error_control(file_id, trim(dset_name),  h5t_native_double, dspace_id, dset_id, iostat)
    call h5dwrite_f(dset_id, h5t_native_double, data, dims, iostat)
    call h5dclose_f(dset_id,iostat)
    call h5sclose_f(dspace_id,iostat)

    deallocate(dims)

    ! - parameter -
    if(write_params.eqv..true.) then
       call create_parameter_list(parameter_list)
       s_dim = max_length_parameter*5

       allocate(dims(1))
       dims  = max_nbr_parameters

       call h5tcopy_f(h5t_fortran_s1,c_type,iostat)
       call h5tset_size_f(c_type,s_dim,iostat)
       call h5screate_simple_f(1,dims,dspace_id,iostat)
       call h5dcreate_f_error_control(file_id,"params",c_type,dspace_id,dset_id,iostat)
       call h5dwrite_f(dset_id,c_type,parameter_list,dims,iostat)
       call h5tclose_f(c_type,iostat)
       call h5dclose_f(dset_id,iostat)
       call h5sclose_f(dspace_id,iostat)

       deallocate(dims)
    end if

    call h5fclose_f(file_id, iostat)
    call h5close_f(iostat)

  end subroutine store_4d_rk
!!!=================================================================================================


!!!=================================================================================================
  subroutine store_5d_rk(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin, file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , dimension(:,:,:,:,:), intent(in)                   :: data
    character(len=*),                       intent(in)                   :: type
    character(len=*),                       intent(in), optional         :: fname_optin
    integer         ,                       intent(in), optional         :: nt_optin
    integer         ,                       intent(in), optional         :: ns_optin
    character(len=*),                       intent(in), optional         :: dset_name_optin
    logical         ,                     intent(in), optional           :: file_overwrite_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:), allocatable                          :: dims
    integer(hid_t)                                                       :: file_id, c_type, dset_id, dspace_id 
    integer(size_t)                                                      :: s_dim 
    integer                                                              :: iostat
    logical                                                              :: file_exist, write_params = .true., file_overwrite = .false.

    character(len=max_length_parameter*5), dimension(max_nbr_parameters) :: parameter_list
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)

       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5")
    endif

    if(present(file_overwrite_optin))then
       file_overwrite = file_overwrite_optin
    end if

    ! open hdf5
    call h5open_f(iostat)

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)

       ! inquire file
       inquire(file = trim(fname),exist = file_exist, iostat = iostat)

       ! open exisiting file or create new one
       if((file_exist.eqv..true.).and.(file_overwrite.eqv..false.)) then
          call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)
          ! do not write params (again)
          write_params = .false.
       else
          call h5fcreate_f_error_control(trim(fname), h5f_acc_trunc_f, file_id, iostat)
       end if
    else ! the defaul way
       ! default dset name = 'data'
       dset_name = 'data'

       ! create new h5 file, deletes ALL old stuff, fastest way of writing, the way to go!
       ! params are written as write_params = .true. by default
       call h5fcreate_f_error_control(trim(fname), h5f_acc_trunc_f, file_id, iostat)
    end if

    ! - data - 
    allocate(dims(1:5))
    dims(1) = size(data,1)
    dims(2) = size(data,2)
    dims(3) = size(data,3)
    dims(4) = size(data,4)
    dims(5) = size(data,5)

    call h5screate_simple_f(5, dims, dspace_id, iostat)
    call h5dcreate_f_error_control(file_id, trim(dset_name),  h5t_native_double, dspace_id, dset_id, iostat)
    call h5dwrite_f(dset_id, h5t_native_double, data, dims, iostat)
    call h5dclose_f(dset_id,iostat)
    call h5sclose_f(dspace_id,iostat)

    deallocate(dims)

    ! - parameter -
    if(write_params.eqv..true.) then
       call create_parameter_list(parameter_list)
       s_dim = max_length_parameter*5

       allocate(dims(1))
       dims  = max_nbr_parameters

       call h5tcopy_f(h5t_fortran_s1,c_type,iostat)
       call h5tset_size_f(c_type,s_dim,iostat)
       call h5screate_simple_f(1,dims,dspace_id,iostat)
       call h5dcreate_f_error_control(file_id,"params",c_type,dspace_id,dset_id,iostat)
       call h5dwrite_f(dset_id,c_type,parameter_list,dims,iostat)
       call h5tclose_f(c_type,iostat)
       call h5dclose_f(dset_id,iostat)
       call h5sclose_f(dspace_id,iostat)

       deallocate(dims)
    end if

    call h5fclose_f(file_id, iostat)
    call h5close_f(iostat)

  end subroutine store_5d_rk
!!!=================================================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! l o a d !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!=================================================================================================
  subroutine load_5d_rk(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , dimension(:,:,:,:,:), intent(out)                  :: data
    character(len=*),                       intent(in)                   :: type
    character(len=*),                       intent(in), optional         :: fname_optin
    integer         ,                       intent(in), optional         :: nt_optin
    integer         ,                       intent(in), optional         :: ns_optin
    character(len=*),                       intent(in), optional         :: dset_name_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:), allocatable                          :: dims
    integer(hid_t)                                                       :: file_id, dset_id
    integer                                                              :: iostat
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)

       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5")
    endif

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! open hdf5
    call h5open_f(iostat)

    ! open hdf5 file
    call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)

    ! - data - 
    allocate(dims(1:5))
    dims(1) = size(data,1)
    dims(2) = size(data,2)
    dims(3) = size(data,3)
    dims(4) = size(data,4)
    dims(5) = size(data,5)

    call h5dopen_f_error_control(file_id, trim(dset_name), dset_id, iostat)
    call h5dread_f(dset_id, h5t_native_double, data, dims, iostat)
    call h5dclose_f(dset_id,iostat)

    deallocate(dims)

    ! close hdf5 (file)
    call h5fclose_f(file_id, iostat)
    call h5close_f(iostat)

  end subroutine load_5d_rk
!!!=================================================================================================


!!!=================================================================================================
  subroutine load_4d_rk(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , dimension(:,:,:,:)  , intent(out)                  :: data
    character(len=*),                       intent(in)                   :: type
    character(len=*),                       intent(in), optional         :: fname_optin
    integer         ,                       intent(in), optional         :: nt_optin
    integer         ,                       intent(in), optional         :: ns_optin
    character(len=*),                       intent(in), optional         :: dset_name_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:), allocatable                          :: dims
    integer(hid_t)                                                       :: file_id, dset_id
    integer                                                              :: iostat
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)

       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5")
    endif

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! open hdf5
    call h5open_f(iostat)

    ! open hdf5 file
    call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)

    ! - data - 
    allocate(dims(1:4))
    dims(1) = size(data,1)
    dims(2) = size(data,2)
    dims(3) = size(data,3)
    dims(4) = size(data,4)

    call h5dopen_f_error_control(file_id, trim(dset_name), dset_id, iostat)
    call h5dread_f(dset_id, h5t_native_double, data, dims, iostat)
    call h5dclose_f(dset_id,iostat)

    deallocate(dims)

    ! close hdf5 (file)
    call h5fclose_f(file_id, iostat)
    call h5close_f(iostat)

  end subroutine load_4d_rk
!!!=================================================================================================

!!!=================================================================================================
  subroutine load_sub_4d_rk(data, type, lbounds, ubounds, fname_optin, nt_optin, ns_optin, dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , dimension(:,:,:,:)  , intent(out)                  :: data
    character(len=*),                       intent(in)                   :: type
    integer         , dimension(4)        , intent(in)                   :: lbounds
    integer         , dimension(4)        , intent(in)                   :: ubounds
    character(len=*),                       intent(in), optional         :: fname_optin
    integer         ,                       intent(in), optional         :: nt_optin
    integer         ,                       intent(in), optional         :: ns_optin
    character(len=*),                       intent(in), optional         :: dset_name_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t)                                                       :: file_id, dset_id
    integer(hsize_t), dimension(4)                                       :: dims_file, dims_max_file
    integer(hsize_t), dimension(4)                                       :: dims_local
    integer(hid_t)                                                       :: filespace     
    integer(hid_t)                                                       :: memspace      
    integer(hsize_t), dimension(4)                                       :: count = 1
    integer(hssize_t), dimension(4)                                      :: offset
    integer(hsize_t), dimension(4)                                       :: stride = 1
    integer                                                              :: iostat
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)

       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5")
    endif

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! open hdf5
    call h5open_f(iostat)

    ! open hdf5 file
    call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)

    ! get dims of data (file)
    call h5dopen_f_error_control(file_id, dset_name, dset_id, iostat)
    call h5dget_space_f(dset_id, filespace, iostat)
    call h5sget_simple_extent_dims_f( filespace, dims_file, dims_max_file, iostat)

    ! get dims of data (memory)
    dims_local(1) = size(data,1)
    dims_local(2) = size(data,2)
    dims_local(3) = size(data,3)
    dims_local(4) = size(data,4)

    if (sum(dims_local - (ubounds - lbounds + 1)).ne.0) then
       write(*,*) "warning in hdf5_single.f90:load_sub_4d_rk:dimension mismatch",dims_local,(ubounds - lbounds + 1)
    end if

    ! create hyperslap, memspace
    offset = lbounds - 1

    call h5sselect_hyperslab_f(filespace, h5s_select_set_f, offset, count, iostat, stride, dims_local)
    call h5screate_simple_f(4, dims_local, memspace, iostat)

    ! read data
    call h5dread_f( dset_id, h5t_native_double, data, dims_local, iostat, mem_space_id = memspace, file_space_id = filespace)

    ! close hdf5 (all)
    call h5sclose_f(memspace , iostat)
    call h5sclose_f(filespace, iostat)
    call h5dclose_f(dset_id  , iostat)
    call h5fclose_f(file_id  , iostat)
    call h5close_f(iostat)

  end subroutine load_sub_4d_rk
!!!=================================================================================================

!!!=================================================================================================
  subroutine load_4d_i(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer         , dimension(:,:,:,:)  , intent(out)                  :: data
    character(len=*),                       intent(in)                   :: type
    character(len=*),                       intent(in), optional         :: fname_optin
    integer         ,                       intent(in), optional         :: nt_optin
    integer         ,                       intent(in), optional         :: ns_optin
    character(len=*),                       intent(in), optional         :: dset_name_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:), allocatable                          :: dims
    integer(hid_t)                                                       :: file_id, dset_id
    integer                                                              :: iostat
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)

       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5")
    endif

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! open hdf5
    call h5open_f(iostat)

    ! open hdf5 file
    call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)

    ! - data - 
    allocate(dims(1:4))
    dims(1) = size(data,1)
    dims(2) = size(data,2)
    dims(3) = size(data,3)
    dims(4) = size(data,4)

    call h5dopen_f_error_control(file_id, trim(dset_name), dset_id, iostat)
    call h5dread_f(dset_id, h5t_native_integer, data, dims, iostat)
    call h5dclose_f(dset_id,iostat)

    deallocate(dims)

    ! close hdf5 (file)
    call h5fclose_f(file_id, iostat)
    call h5close_f(iostat)

  end subroutine load_4d_i
!!!=================================================================================================

!!!=================================================================================================
  subroutine load_3d_rk(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , dimension(:,:,:)    , intent(out)                  :: data
    character(len=*),                       intent(in)                   :: type
    character(len=*),                       intent(in), optional         :: fname_optin
    integer         ,                       intent(in), optional         :: nt_optin
    integer         ,                       intent(in), optional         :: ns_optin
    character(len=*),                       intent(in), optional         :: dset_name_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:), allocatable                          :: dims
    integer(hid_t)                                                       :: file_id, dset_id
    integer                                                              :: iostat
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)

       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5")
    endif

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! open hdf5
    call h5open_f(iostat)

    ! open hdf5 file
    call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)

    ! - data - 
    allocate(dims(1:3))
    dims(1) = size(data,1)
    dims(2) = size(data,2)
    dims(3) = size(data,3)

    call h5dopen_f_error_control(file_id, trim(dset_name), dset_id, iostat)
    call h5dread_f(dset_id, h5t_native_double, data, dims, iostat)
    call h5dclose_f(dset_id,iostat)

    deallocate(dims)

    ! close hdf5 (file)
    call h5fclose_f(file_id, iostat)
    call h5close_f(iostat)

  end subroutine load_3d_rk
!!!=================================================================================================

!!!=================================================================================================
  subroutine load_2d_rk(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , dimension(:,:)      , intent(out)                  :: data
    character(len=*),                       intent(in)                   :: type
    character(len=*),                       intent(in), optional         :: fname_optin
    integer         ,                       intent(in), optional         :: nt_optin
    integer         ,                       intent(in), optional         :: ns_optin
    character(len=*),                       intent(in), optional         :: dset_name_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:), allocatable                          :: dims
    integer(hid_t)                                                       :: file_id, dset_id
    integer                                                              :: iostat
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)

       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5")
    endif

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! open hdf5
    call h5open_f(iostat)

    ! open hdf5 file
    call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)

    ! - data - 
    allocate(dims(1:2))
    dims(1) = size(data,1)
    dims(2) = size(data,2)

    call h5dopen_f_error_control(file_id, trim(dset_name), dset_id, iostat)
    call h5dread_f(dset_id, h5t_native_double, data, dims, iostat)
    call h5dclose_f(dset_id,iostat)

    deallocate(dims)

    ! close hdf5 (file)
    call h5fclose_f(file_id, iostat)
    call h5close_f(iostat)

  end subroutine load_2d_rk
!!!=================================================================================================


!!!=================================================================================================
  subroutine load_1d_rk(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk)   , dimension(:)        , intent(out)                  :: data
    character(len=*),                       intent(in)                   :: type
    character(len=*),                       intent(in), optional         :: fname_optin
    integer         ,                       intent(in), optional         :: nt_optin
    integer         ,                       intent(in), optional         :: ns_optin
    character(len=*),                       intent(in), optional         :: dset_name_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:), allocatable                          :: dims
    integer(hid_t)                                                       :: file_id, dset_id
    integer                                                              :: iostat
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)

       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5")
    endif

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! open hdf5
    call h5open_f(iostat)

    ! open hdf5 file
    call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)

    ! - data - 
    allocate(dims(1))
    dims(1) = size(data,1)

    call h5dopen_f_error_control(file_id, trim(dset_name), dset_id, iostat)
    call h5dread_f(dset_id, h5t_native_double, data, dims, iostat)
    call h5dclose_f(dset_id,iostat)

    deallocate(dims)

    ! close hdf5 (file)
    call h5fclose_f(file_id, iostat)
    call h5close_f(iostat)

  end subroutine load_1d_rk
!!!=================================================================================================

!!!=================================================================================================
  subroutine load_1d_i(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer         , dimension(:)        , intent(out)                  :: data
    character(len=*),                       intent(in)                   :: type
    character(len=*),                       intent(in), optional         :: fname_optin
    integer         ,                       intent(in), optional         :: nt_optin
    integer         ,                       intent(in), optional         :: ns_optin
    character(len=*),                       intent(in), optional         :: dset_name_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:), allocatable                          :: dims
    integer(hid_t)                                                       :: file_id, dset_id
    integer                                                              :: iostat
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)

       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5")
    endif

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! open hdf5
    call h5open_f(iostat)

    ! open hdf5 file
    call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)

    ! - data - 
    allocate(dims(1))
    dims(1) = size(data,1)

    call h5dopen_f_error_control(file_id, trim(dset_name), dset_id, iostat)
    ! call h5dread_f(dset_id, h5t_native_double, data, dims, iostat)
    call h5dread_f(dset_id, h5t_native_integer, data, dims, iostat)
    call h5dclose_f(dset_id,iostat)

    deallocate(dims)

    ! close hdf5 (file)
    call h5fclose_f(file_id, iostat)
    call h5close_f(iostat)

  end subroutine load_1d_i
!!!=================================================================================================



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! m i s c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!=================================================================================================
  subroutine get_dset_dim(dims, type, fname_optin, nt_optin, ns_optin, dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer         , dimension(:)        , intent(out)                  :: dims
    character(len=*),                       intent(in)                   :: type
    character(len=*),                       intent(in), optional         :: fname_optin
    integer         ,                       intent(in), optional         :: nt_optin
    integer         ,                       intent(in), optional         :: ns_optin
    character(len=*),                       intent(in), optional         :: dset_name_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:),allocatable                           :: dims_file, dims_max_file
    integer(hid_t)                                                       :: file_id, file_space, dset_id
    integer                                                              :: iostat,n1
    character(len=max_length_fname)                                      :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    dims=0

    ! handle optional arguments
    if(present(fname_optin))then
       ! use optional fname only
       fname = trim(fname_optin)
    else
       if ((present(nt_optin)).and.(present(ns_optin))) then
          ! construct fname using optional timestep and subset
          call get_fname(fname,type,nt_optin,ns_optin)

       else if (present(nt_optin)) then
          ! construct fname using optional timestep
          call get_fname(fname,type,nt = nt_optin)

       else if (present(ns_optin)) then
          ! construct fname using optional subset
          call get_fname(fname,type,ns = ns_optin)

       else           
          ! construct fname (the default way for all io operations)
          call get_fname(fname,type)

       end if

       ! add filename extension
       fname = trim(trim(fname)//".h5")
    endif

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! allocate dims
    n1 = size(dims)
    allocate(    dims_file(n1))
    allocate(dims_max_file(n1))

    ! open hdf5
    call h5open_f(iostat)

    ! open hdf5 file
    call h5fopen_f_error_control(trim(fname), h5f_acc_rdwr_f, file_id, iostat)

    ! get dims of data
    call h5dopen_f_error_control(file_id, dset_name, dset_id, iostat)
    call h5dget_space_f(dset_id, file_space, iostat)
    call h5sget_simple_extent_dims_f( file_space, dims_file, dims_max_file, iostat) ! who knows the difference between size and max size ???

    ! close hdf5 (file)
    call h5dclose_f(dset_id, iostat)
    call h5fclose_f(file_id, iostat)
    call h5close_f(iostat)

    ! set output
    dims = int(dims_file)

  end subroutine get_dset_dim
!!!=================================================================================================

!!!=================================================================================================  
  subroutine h5fopen_f_error_control(fname, file_access_flag, file_id, iostat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*)                                  , intent(in)       :: fname
    integer                                           , intent(in)       :: file_access_flag
    integer(hid_t)                                    , intent(inout)    :: file_id
    integer                                           , intent(inout)    :: iostat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call h5fopen_f(trim(fname), file_access_flag, file_id, iostat)

    if (iostat.ne.0) then
       write(*,*) "error in hdf5_single.f90:h5fopen_f_error_control:iostat:",iostat,":",trim(fname)
       stop
    end if

  end subroutine h5fopen_f_error_control
!!!=================================================================================================

!!!=================================================================================================  
  subroutine h5fcreate_f_error_control(fname, file_access_flag, file_id, iostat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*)                                  , intent(in)       :: fname
    integer                                           , intent(in)       :: file_access_flag
    integer(hid_t)                                    , intent(inout)    :: file_id
    integer                                           , intent(inout)    :: iostat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call h5fcreate_f(trim(fname), file_access_flag, file_id, iostat)

    if (iostat.ne.0) then
       write(*,*) "error in hdf5_single.f90:h5fcreate_f_error_control:iostat:",iostat,":",trim(fname)
       stop
    end if

  end subroutine h5fcreate_f_error_control
!!!=================================================================================================

!!!=================================================================================================  
  subroutine h5dopen_f_error_control(file_id, dset_name, dset_id, iostat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t)                                    , intent(in)       :: file_id
    character(len=*)                                  , intent(in)       :: dset_name
    integer(hid_t)                                    , intent(inout)    :: dset_id
    integer                                           , intent(inout)    :: iostat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call h5dopen_f(file_id, dset_name, dset_id, iostat)

    if (iostat.ne.0) then
       write(*,*) "error in hdf5_single.f90:h5dopen_f_error_control:iostat:",iostat,":",trim(dset_name)
       stop
    end if

  end subroutine h5dopen_f_error_control
!!!=================================================================================================

!!!=================================================================================================  
  subroutine h5dcreate_f_error_control(file_id, dset_name, data_type, dspace_id, dset_id, iostat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t)                                    , intent(in)       :: file_id
    character(len=*)                                  , intent(in)       :: dset_name
    integer(hid_t)                                    , intent(in)       :: data_type
    integer(hid_t)                                    , intent(inout)    :: dspace_id, dset_id
    integer                                           , intent(inout)    :: iostat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call h5dcreate_f(file_id, trim(dset_name), data_type, dspace_id, dset_id, iostat)

    if (iostat.ne.0) then
       write(*,*) "error in hdf5_single.f90:h5dcreate_f_error_control:iostat:",iostat,":",trim(dset_name)
       stop
    end if

  end subroutine h5dcreate_f_error_control
!!!=================================================================================================    

end module hdf5_single
