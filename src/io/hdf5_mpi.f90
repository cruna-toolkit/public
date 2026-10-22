!> Parallel HDF5 I/O backend with transparent serial fallback.
!!
!! This module exposes the same public `store`, `load`, `load_sub`, and `get_dset_dim`
!! interfaces as `hdf5_single`. File types selected for MPI-I/O are handled through a
!! shared-file Parallel-HDF5 backend; all other operations are delegated to the legacy
!! serial/per-rank implementation.
!!
!! Explicit filenames preserve the legacy behavior and bypass automatic MPI-I/O policy
!! selection. For automatically generated filenames, `io_helper:get_fname` determines
!! whether a file type uses the shared Parallel-HDF5 path.
!!
!! REAL(kind=rk) memory datatypes always follow `rk`. For newly created parallel datasets,
!! `params%io%mpiio_single_precision` may request IEEE 32-bit file storage while HDF5 performs
!! conversion between the in-memory and on-disk representations. Reads always derive the memory
!! datatype from `rk`; the stored datatype is converted by HDF5 when required. Integer datasets
!! and the character dataset `params` are not affected by this precision switch.
!!
!! Parallel HDF5 operations use `params%parallelism%block_comm`; participating ranks must enter
!! collective metadata operations consistently. In particular, `params%io%mpiio_single_precision`
!! must have the same value on all ranks in the communicator. Requested collective dataset I/O may
!! still be executed internally by HDF5 through a non-collective low-level path when datatype
!! conversion or other transformations require it.
module hdf5_mpi

  ! HDF5 provides the Fortran bindings and HDF5 datatype/property-list constants.
  ! MPI is required for communicator collectives, error propagation, and the MPIO file driver.
  USE hdf5
  USE mpi
  use helper
  use io_helper
  use parameter

  ! Reuse the complete serial/per-rank implementation.
  ! Import the legacy implementation under explicit local aliases.  The public generics below
  ! can therefore expose the same API while choosing at runtime between per-rank and shared-file I/O.
  use hdf5_single , only:                  &
       store_single        => store,       &
       load_single         => load,        &
       load_sub_single     => load_sub,    &
       get_dset_dim_single => get_dset_dim

  !> Name identifying this I/O backend.
  character(len=max_length_parameter), parameter, public   :: io_type_name = "hdf5_mpi"

  ! Keep implementation helpers private.  Callers see only the same four generic entry points
  ! that are provided by hdf5_single.
  private

  public :: get_dset_dim
  public :: load
  public :: load_sub
  public :: store

  !> Generic load interface.
  !!
  !! Resolves the concrete load wrapper from array rank and element type. The wrapper
  !! then selects either the shared Parallel-HDF5 backend or `hdf5_single`.
  interface load
     module procedure load_5d_rk
     module procedure load_4d_rk
     module procedure load_4d_i
     module procedure load_3d_rk
     module procedure load_2d_rk
     module procedure load_1d_rk
     module procedure load_1d_i
  end interface load

  !> Generic explicit-bounds subregion load interface.
  !!
  !! Loads a caller-selected hyperslab while preserving the serial/parallel dispatch policy.
  interface load_sub
     module procedure load_sub_4d_rk
  end interface load_sub

  !> Generic store interface.
  !!
  !! Resolves the concrete store wrapper from array rank and element type and dispatches
  !! to Parallel HDF5 only for file types configured for MPI-I/O.
  interface store
     module procedure store_1d_rk
     module procedure store_1d_i
     module procedure store_2d_rk
     module procedure store_3d_rk
     module procedure store_4d_i
     module procedure store_4d_rk
     module procedure store_5d_rk
  end interface store


  ! All routines below are module procedures.  Helpers are private; only the generic interfaces
  ! declared above are visible outside this module.
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! dispatch / filename policy helpers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!=================================================================================================
  !> Determine whether the requested operation shall use the shared-file Parallel-HDF5 backend.
  !!
  !! Backward compatibility with hdf5_single is important here: historically an explicitly supplied
  !! fname_optin bypassed get_fname completely. Consequently type could be empty, a dummy string, or
  !! otherwise not present in io_helper:get_fname without affecting the I/O operation. To preserve
  !! this behavior, every operation with fname_optin is classified as serial/per-rank I/O and
  !! forwarded unchanged to hdf5_single. In that case type is deliberately not evaluated at all.
  !! Only when fname_optin is absent is get_fname used as the authoritative source for the MPI-I/O
  !! policy. The optional timestep/subset values are forwarded so backend selection follows exactly
  !! the same logical file-type policy as ordinary filename generation.
  !!
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param mpi_io Returns whether the shared-file Parallel-HDF5 backend is selected.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  subroutine get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*), intent(in)           :: type
    logical,          intent(out)          :: mpi_io
    integer,          intent(in), optional :: nt_optin
    integer,          intent(in), optional :: ns_optin
    character(len=*), intent(in), optional :: fname_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_fname)        :: fname_dummy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Preserve the exact legacy semantics of hdf5_single for explicitly named files.
    ! No call to get_fname is made here, so type may safely be '', 'dummy', or any other placeholder.
    if(present(fname_optin))then
       mpi_io = .false.
       return
    end if

    ! Without an explicit filename, get_fname is the authoritative source for the MPI-I/O policy.
    ! Forward exactly the optional timestep/subset arguments supplied by the caller.
    if(present(nt_optin).and.present(ns_optin))then
       call get_fname(fname_dummy,type,nt=nt_optin,ns=ns_optin,mpiio=mpi_io)
    else if(present(nt_optin))then
       call get_fname(fname_dummy,type,nt=nt_optin,mpiio=mpi_io)
    else if(present(ns_optin))then
       call get_fname(fname_dummy,type,ns=ns_optin,mpiio=mpi_io)
    else
       call get_fname(fname_dummy,type,mpiio=mpi_io)
    end if

  end subroutine get_mpiio_policy
!!!=================================================================================================
  !> Construct the filename used by the shared-file Parallel-HDF5 backend.
  !!
  !! get_fname also returns the MPI-I/O policy. Reaching this routine for a non-MPI-I/O type
  !! indicates an internal dispatch error and aborts block_comm. The helper still accepts
  !! fname_optin for internal completeness. Public dispatchers preserve legacy behavior by routing
  !! explicitly named files to hdf5_single before this routine is reached.
  !!
  !! @param fname Filename selected for the shared-file Parallel-HDF5 backend.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  subroutine get_parallel_fname(fname,type,fname_optin,nt_optin,ns_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_fname), intent(out)          :: fname
    character(len=*),                intent(in)           :: type
    character(len=*),                intent(in), optional :: fname_optin
    integer,                         intent(in), optional :: nt_optin
    integer,                         intent(in), optional :: ns_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=max_length_fname)                       :: fname_generated
    logical                                               :: mpi_io
    integer                                               :: mpierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Build the canonical filename exactly as io_helper would do for the current timestep/subset
    ! and retrieve the MPI-I/O policy in the same call.
    if(present(nt_optin).and.present(ns_optin))then
       call get_fname(fname_generated,type,nt=nt_optin,ns=ns_optin,mpiio=mpi_io)
    else if(present(nt_optin))then
       call get_fname(fname_generated,type,nt=nt_optin,mpiio=mpi_io)
    else if(present(ns_optin))then
       call get_fname(fname_generated,type,ns=ns_optin,mpiio=mpi_io)
    else
       call get_fname(fname_generated,type,mpiio=mpi_io)
    end if

    ! This helper must only be reached after a dispatcher selected the parallel backend.
    ! A false policy here therefore indicates a programming/dispatch consistency error.
    if(.not.mpi_io)then
       if(params%parallelism%block_image.eq.1)then
          write(*,*) "internal error in hdf5_mpi:get_parallel_fname:"
          write(*,*) "parallel backend called for non-MPI-I/O type: ",trim(type)
       end if
       call mpi_abort(params%parallelism%block_comm,1,mpierr)
    end if

    ! If this internal helper is called with an explicit filename, use it verbatim. Under normal public
    ! dispatch this branch is not reached because fname_optin is handled by the legacy serial path first.
    if(present(fname_optin))then
       fname = trim(fname_optin)
    else
       fname = trim(fname_generated)//".h5"
    end if

  end subroutine get_parallel_fname
!!!=================================================================================================
  !> Protect 1D/2D wrappers from accidentally using an MPI-I/O file type.
  !!
  !! These ranks are intentionally handled by hdf5_single, so a type that requests a shared MPI file
  !! would otherwise make several ranks access the same serial filename.
  !!
  !! @param where Name of the calling context used in diagnostics.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param mpi_io Returns whether the shared-file Parallel-HDF5 backend is selected.
  subroutine reject_mpiio_for_serial_rank(where,type,mpi_io)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*), intent(in) :: where
    character(len=*), intent(in) :: type
    logical,          intent(in) :: mpi_io
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                      :: mpierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! A shared MPI filename must never be used by hdf5_single from multiple ranks.  Abort early
    ! rather than allowing concurrent serial access to the same file.
    if(mpi_io)then

       if(params%parallelism%block_image.eq.1)then
          write(*,*) "error in hdf5_mpi.f90:",trim(where), &
               ": MPI-I/O type used with 1D/2D serial backend: ",trim(type)
       end if

       call mpi_abort(params%parallelism%block_comm,1,mpierr)

    end if

  end subroutine reject_mpiio_for_serial_rank
!!!=================================================================================================
  !> Return the global 1-based spatial bounds owned by the current MPI rank.
  !!
  !! The geometry indices are constant after initialization, therefore they are read from the
  !! parameter system once and cached locally in each process.
  !!
  !! @param lbounds Global 1-based lower spatial bounds owned by the current MPI rank.
  !! @param ubounds Global 1-based upper spatial bounds owned by the current MPI rank.
  subroutine get_parallel_spatial_bounds(lbounds,ubounds)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(3), intent(out) :: lbounds
    integer, dimension(3), intent(out) :: ubounds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(3), save        :: lbounds_cache
    integer, dimension(3), save        :: ubounds_cache
    logical, save                      :: initialized = .false.

    integer                            :: x10i,x11i
    integer                            :: x20i,x21i
    integer                            :: x30i,x31i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! The decomposition does not change during a run.  Avoid repeated parameter lookups by
    ! initializing a process-local cache on the first call.
    if(.not.initialized)then

       call get_parameter(x10i,'geom.xi10i')
       call get_parameter(x11i,'geom.xi11i')
       call get_parameter(x20i,'geom.xi20i')
       call get_parameter(x21i,'geom.xi21i')
       call get_parameter(x30i,'geom.xi30i')
       call get_parameter(x31i,'geom.xi31i')

       lbounds_cache = (/x10i,x20i,x30i/)
       ubounds_cache = (/x11i,x21i,x31i/)

       initialized = .true.

    end if

    ! Return copies of the cached global 1-based indices to the caller.
    lbounds = lbounds_cache
    ubounds = ubounds_cache

  end subroutine get_parallel_spatial_bounds
!!!=================================================================================================
  !> Map the project-wide REAL kind rk to the matching native HDF5 memory datatype.
  !!
  !! This keeps the HDF5 interface independent of whether rk is single or double precision.
  !!
  !! @param mem_type HDF5 datatype describing the REAL(kind=rk) representation in memory.
  subroutine get_hdf5_real_memory_type(mem_type)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t), intent(out) :: mem_type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! The memory datatype must describe the actual Fortran array passed to HDF5.
    ! Do not hard-code H5T_NATIVE_DOUBLE here: rk is the single source of truth for
    ! floating-point precision inside CRUNA. h5kind_to_type maps this Fortran KIND
    ! to the matching native HDF5 REAL datatype.
    !
    ! Examples:
    !   rk = single precision -> mem_type represents a native 32-bit REAL
    !   rk = double precision -> mem_type represents a native 64-bit REAL
    !
    ! This datatype is used as the memory datatype in h5dwrite_f/h5dread_f. It does
    ! not by itself determine the datatype stored in the file.
    mem_type = h5kind_to_type(rk,h5_real_kind)

  end subroutine get_hdf5_real_memory_type
!!!=================================================================================================
  !> Select the HDF5 memory datatype and the datatype physically stored in a new REAL dataset.
  !!
  !! Memory always follows rk. File precision follows rk unless mpiio_single_precision requests IEEE
  !! float32 storage, in which case HDF5 performs the conversion.
  !!
  !! @param mem_type HDF5 datatype describing the REAL(kind=rk) representation in memory.
  !! @param file_type HDF5 datatype used for the dataset representation in the file.
  subroutine get_hdf5_real_store_types(mem_type,file_type)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t), intent(out) :: mem_type
    integer(hid_t), intent(out) :: file_type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! First determine the datatype of the actual Fortran data in memory from rk.
    call get_hdf5_real_memory_type(mem_type)

    ! By default the HDF5 file precision follows rk exactly. Using the same HDF5
    ! datatype for memory and file also means that no precision conversion is needed.
    file_type = mem_type

    ! Optional storage reduction for MPI-I/O datasets:
    !
    ! If requested, only the FILE datatype is changed to IEEE float32. mem_type still
    ! describes REAL(kind=rk). HDF5 therefore performs the precision conversion during
    ! the transfer. No temporary single-precision Fortran array is required here.
    !
    ! Example for rk = double precision:
    !   mem_type  -> native float64
    !   file_type -> IEEE float32
    !   HDF5      -> converts float64 -> float32 during h5dwrite_f
    !
    ! This parameter has to be identical on all ranks of block_comm because dataset
    ! creation and its datatype are collective metadata operations.
    if(params%io%mpiio_single_precision)then
       file_type = h5t_ieee_f32le
    end if

  end subroutine get_hdf5_real_store_types
!!!=================================================================================================
  !> Central fatal-error path for MPI-HDF5 operations.
  !!
  !! Only block_image=1 prints the generic error message; all ranks in block_comm are then
  !! terminated together with MPI_ABORT.
  !!
  !! @param where Name of the calling context used in diagnostics.
  !! @param iostat Error/status code associated with the failed MPI-HDF5 operation.
  subroutine mpi_hdf5_abort(where,iostat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*), intent(in) :: where
    integer,          intent(in) :: iostat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                      :: mpierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Avoid duplicate diagnostics from every rank.  MPI_ABORT below still terminates all ranks
    ! belonging to the block communicator.
    if(params%parallelism%block_image.eq.1)then
       write(*,*) "error in hdf5_mpi.f90:",trim(where),":iostat:",iostat
    end if

    call mpi_abort(params%parallelism%block_comm,1,mpierr)

  end subroutine mpi_hdf5_abort
!!!=================================================================================================
  !> Verify collectively that each local array shape matches the bounds assigned to that rank.
  !!
  !! The lowest block_image reporting an error prints the diagnostics. Fatal checks abort the
  !! communicator; non-fatal checks are used where only a warning is desired.
  !!
  !! @param where Name of the calling context used in diagnostics.
  !! @param dims_local Local array extents for the current MPI rank.
  !! @param lbounds Lower bounds of the selected region, using global 1-based indices.
  !! @param ubounds Upper bounds of the selected region, using global 1-based indices.
  !! @param fatal Optional flag selecting fatal abort behavior; defaults to true.
  subroutine check_parallel_dimensions(where,dims_local,lbounds,ubounds,fatal)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*),      intent(in)           :: where
    integer, dimension(:), intent(in)           :: dims_local
    integer, dimension(:), intent(in)           :: lbounds
    integer, dimension(:), intent(in)           :: ubounds
    logical,               intent(in), optional :: fatal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                     :: error_image,error_image_local
    integer                                     :: mpierr
    logical                                     :: fatal_error
    logical                                     :: local_error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Dimension mismatches are fatal by default.  load_sub uses fatal=.false. where the historical
    ! behavior is to warn before proceeding with the caller-provided output buffer.
    fatal_error = .true.
    if(present(fatal)) fatal_error = fatal

    ! For an inclusive 1-based interval [lbounds,ubounds], the expected local extent is
    ! ubounds-lbounds+1 in every dimension.
    local_error = any(dims_local.ne.(ubounds-lbounds+1))

    ! Encode 'no local error' as a very large rank index.  MPI_MIN then identifies one deterministic
    ! failing rank without printing the same diagnostic on every process.
    error_image_local = huge(0)
    if(local_error) error_image_local = params%parallelism%block_image

    call mpi_allreduce(error_image_local,error_image,1,mpi_integer,mpi_min, &
         params%parallelism%block_comm,mpierr)
    if(mpierr.ne.mpi_success) call mpi_hdf5_abort("check_parallel_dimensions:mpi_allreduce",mpierr)

    if(error_image.ne.huge(0))then

       if(params%parallelism%block_image.eq.error_image)then
          if(fatal_error)then
             write(*,*) "error in hdf5_mpi.f90:",trim(where),":dimension mismatch"
          else
             write(*,*) "warning in hdf5_mpi.f90:",trim(where),":dimension mismatch"
          end if
          write(*,*) "block_image:",params%parallelism%block_image
          write(*,*) "shape(data):",dims_local
          write(*,*) "lbounds:",lbounds
          write(*,*) "ubounds:",ubounds
       end if

       if(fatal_error) call mpi_abort(params%parallelism%block_comm,1,mpierr)

    end if

  end subroutine check_parallel_dimensions
!!!=================================================================================================
  !> Verify collectively that every selected hyperslab lies inside the dataset extent.
  !!
  !! Negative offsets and selections extending beyond dims_file are treated as fatal consistency
  !! errors.
  !!
  !! @param where Name of the calling context used in diagnostics.
  !! @param offset Zero-based HDF5 hyperslab offset for the current MPI rank.
  !! @param dims_local Local array extents for the current MPI rank.
  !! @param dims_file Global HDF5 dataset extents.
  subroutine check_parallel_hyperslab(where,offset,dims_local,dims_file)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*),                intent(in) :: where
    integer(hssize_t), dimension(:), intent(in) :: offset
    integer(hsize_t),  dimension(:), intent(in) :: dims_local
    integer(hsize_t),  dimension(:), intent(in) :: dims_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                       :: error_image,error_image_local
    integer                                       :: mpierr
    logical                                       :: local_error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! HDF5 offsets are zero based.  Reject negative offsets first, then check that
    ! offset+dims_local does not extend beyond the dataset extent.
    local_error = any(offset.lt.0_hssize_t)
    if(.not.local_error)then
       local_error = any(int(offset,kind=hsize_t)+dims_local.gt.dims_file)
    end if

    ! Select one failing rank collectively so exactly one process emits the detailed diagnostics.
    error_image_local = huge(0)
    if(local_error) error_image_local = params%parallelism%block_image

    call mpi_allreduce(error_image_local,error_image,1,mpi_integer,mpi_min, &
         params%parallelism%block_comm,mpierr)
    if(mpierr.ne.mpi_success) call mpi_hdf5_abort("check_parallel_hyperslab:mpi_allreduce",mpierr)

    if(error_image.ne.huge(0))then

       if(params%parallelism%block_image.eq.error_image)then
          write(*,*) "error in hdf5_mpi.f90:",trim(where),":hyperslab outside dataset"
          write(*,*) "block_image:",params%parallelism%block_image
          write(*,*) "dims_file:",dims_file
          write(*,*) "offset:",offset
          write(*,*) "dims_local:",dims_local
       end if

       call mpi_abort(params%parallelism%block_comm,1,mpierr)

    end if

  end subroutine check_parallel_hyperslab
!!!=================================================================================================
  !> Create the HDF5 file-access property list used for Parallel HDF5.
  !!
  !! The MPIO VFD is attached to params%parallelism%block_comm with MPI_INFO_NULL, so all file
  !! open/create calls using this FAPL operate collectively on that communicator.
  !!
  !! @param fapl_id Created HDF5 file-access property-list identifier configured for MPIO.
  subroutine create_parallel_fapl(fapl_id)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t), intent(out) :: fapl_id
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                     :: iostat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Start from a file-access property list and replace the default serial file driver with MPIO.
    call h5pcreate_f(h5p_file_access_f,fapl_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("create_parallel_fapl:h5pcreate_f",iostat)

    ! block_comm defines the set of processes that collectively participate in this shared file.
    call h5pset_fapl_mpio_f(fapl_id,params%parallelism%block_comm,mpi_info_null,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("create_parallel_fapl:h5pset_fapl_mpio_f",iostat)

  end subroutine create_parallel_fapl
!!!=================================================================================================
  !> Open an existing HDF5 file collectively for read-only access using the MPIO file driver.
  !!
  !! The temporary file-access property list is closed immediately after H5Fopen because the open
  !! file keeps the required access configuration.
  !!
  !! @param fname Resolved filename used for the HDF5 operation.
  !! @param file_id HDF5 file identifier.
  subroutine open_parallel_file_read(fname,file_id)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*), intent(in)  :: fname
    integer(hid_t),   intent(out) :: file_id
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t)                :: fapl_id
    integer                       :: iostat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! The FAPL carries the communicator and selects the MPI-IO virtual file driver.
    call create_parallel_fapl(fapl_id)

    call h5fopen_f(trim(fname),h5f_acc_rdonly_f,file_id,iostat,access_prp=fapl_id)
    if(iostat.ne.0) call mpi_hdf5_abort("open_parallel_file_read:h5fopen_f",iostat)

    call h5pclose_f(fapl_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("open_parallel_file_read:h5pclose_f",iostat)

  end subroutine open_parallel_file_read
!!!=================================================================================================
  !> Collectively open an existing shared file or create/truncate a new one for output.
  !!
  !! When a named dataset is appended, block rank 0 first checks whether the file exists and
  !! broadcasts the result so every rank performs the same HDF5 metadata operation. write_params
  !! tells the caller whether the file is new and therefore still needs the canonical params
  !! dataset.
  !!
  !! @param fname Resolved filename used for the HDF5 operation.
  !! @param append_dataset Whether the operation may append a named dataset to an existing file.
  !! @param file_overwrite Whether an existing target file may be truncated/overwritten.
  !! @param file_id HDF5 file identifier.
  !! @param write_params Returns true when the caller must create the canonical params dataset.
  subroutine open_parallel_file_store(fname,append_dataset,file_overwrite,file_id,write_params)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    character(len=*), intent(in)  :: fname
    logical,          intent(in)  :: append_dataset
    logical,          intent(in)  :: file_overwrite
    integer(hid_t),   intent(out) :: file_id
    logical,          intent(out) :: write_params
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t)                :: fapl_id
    integer                       :: iostat
    integer                       :: mpierr
    integer                       :: file_exist_i
    logical                       :: file_exist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Every rank constructs an equivalent MPIO file-access property list before opening/creating.
    call create_parallel_fapl(fapl_id)

    ! Represent the existence result as INTEGER because it is broadcast with MPI_BCAST.
    file_exist_i = 0
    file_exist   = .false.

    ! A supplied dataset name means the caller may be adding a second dataset to an existing file.
    ! In that case preserve the file unless file_overwrite explicitly requests truncation.
    if(append_dataset)then

       ! Only one process queries the filesystem.  The result is broadcast so all ranks take
       ! exactly the same collective HDF5 open/create branch.
       if(params%parallelism%block_image.eq.1)then
          inquire(file=trim(fname),exist=file_exist,iostat=iostat)
          if(iostat.ne.0) call mpi_hdf5_abort("open_parallel_file_store:inquire",iostat)

          if(file_exist)then
             file_exist_i = 1
          else
             ! Represent the existence result as INTEGER because it is broadcast with MPI_BCAST.
             file_exist_i = 0
          end if
       end if

       ! Root in block_comm is MPI rank 0; block_image is rank+1, hence block_image==1 above.
       call mpi_bcast(file_exist_i,1,mpi_integer,0,params%parallelism%block_comm,mpierr)
       if(mpierr.ne.mpi_success) call mpi_hdf5_abort("open_parallel_file_store:mpi_bcast",mpierr)

       file_exist = (file_exist_i.eq.1)

       ! Existing file + no overwrite: open read/write and keep the existing params dataset.
       ! Otherwise create a fresh file and request params to be written by the caller.
       if(file_exist.and.(.not.file_overwrite))then
          call h5fopen_f(trim(fname),h5f_acc_rdwr_f,file_id,iostat,access_prp=fapl_id)
          write_params = .false.
       else
          call h5fcreate_f(trim(fname),h5f_acc_trunc_f,file_id,iostat,access_prp=fapl_id)
          write_params = .true.
       end if

    else

       call h5fcreate_f(trim(fname),h5f_acc_trunc_f,file_id,iostat,access_prp=fapl_id)
       write_params = .true.

    end if

    if(iostat.ne.0) call mpi_hdf5_abort("open_parallel_file_store:h5fopen/create_f",iostat)

    call h5pclose_f(fapl_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("open_parallel_file_store:h5pclose_f",iostat)

  end subroutine open_parallel_file_store
!!!=================================================================================================
  !> Create the dataset-transfer property list used by parallel reads and writes.
  !!
  !! H5FD_MPIO_COLLECTIVE requests collective dataset I/O for transfers that HDF5 can execute
  !! collectively.
  !!
  !! @param dxpl_id Created HDF5 dataset-transfer property-list identifier configured for collective
  !! MPIO.
  subroutine create_collective_dxpl(dxpl_id)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t), intent(out) :: dxpl_id
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                     :: iostat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Dataset-transfer properties control the actual H5Dread/H5Dwrite transfer mode.
    call h5pcreate_f(h5p_dataset_xfer_f,dxpl_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("create_collective_dxpl:h5pcreate_f",iostat)

    ! Request collective MPI-I/O.  HDF5 may internally fall back for operations requiring
    ! transformations such as datatype conversion, but this remains the desired transfer mode.
    call h5pset_dxpl_mpio_f(dxpl_id,h5fd_mpio_collective_f,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("create_collective_dxpl:h5pset_dxpl_mpio_f",iostat)

  end subroutine create_collective_dxpl
!!!=================================================================================================
  !> Create the character dataset params in a newly created shared file.
  !!
  !! Dataset creation is performed collectively. Only block_image=1 writes the canonical parameter
  !! list, followed by a barrier before collective metadata cleanup.
  !!
  !! @param file_id HDF5 file identifier.
  subroutine write_parallel_params(file_id)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t), intent(in)                                           :: file_id
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t)                                                       :: c_type,dset_id,dspace_id
    integer(size_t)                                                      :: s_dim
    integer(hsize_t), dimension(1)                                       :: dims
    integer                                                              :: iostat,mpierr
    character(len=max_length_parameter*5), dimension(max_nbr_parameters) :: parameter_list
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! The serial backend stores the parameter list as fixed-length Fortran strings.  Reuse the
    ! same dimensions and string width so files from both backends keep the same metadata layout.
    s_dim   = max_length_parameter*5
    dims(1) = max_nbr_parameters

    ! Create a fixed-length character datatype matching one element of parameter_list.
    call h5tcopy_f(h5t_fortran_s1,c_type,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("write_parallel_params:h5tcopy_f",iostat)

    call h5tset_size_f(c_type,s_dim,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("write_parallel_params:h5tset_size_f",iostat)

    ! One-dimensional dataspace: one string entry per possible parameter.
    call h5screate_simple_f(1,dims,dspace_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("write_parallel_params:h5screate_simple_f",iostat)

    ! Dataset creation is a collective metadata operation.
    ! Dataset creation modifies shared metadata and is therefore executed by every process.
    call h5dcreate_f(file_id,"params",c_type,dspace_id,dset_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("write_parallel_params:h5dcreate_f",iostat)

    ! Store one canonical parameter list. The parameter array contains rank-local
    ! parallelism entries, therefore only block rank 0 writes the contents.
    ! Only one canonical list is needed.  Other ranks participate in metadata operations but do
    ! not issue a redundant independent data write.
    if(params%parallelism%block_image.eq.1)then
       call create_parameter_list(parameter_list)
       call h5dwrite_f(dset_id,c_type,parameter_list,dims,iostat)
       if(iostat.ne.0) call mpi_hdf5_abort("write_parallel_params:h5dwrite_f",iostat)
    end if

    ! Ensure the independent root write is complete before all ranks close shared metadata objects.
    call mpi_barrier(params%parallelism%block_comm,mpierr)
    if(mpierr.ne.mpi_success) call mpi_hdf5_abort("write_parallel_params:mpi_barrier",mpierr)

    call h5dclose_f(dset_id,iostat)
    call h5sclose_f(dspace_id,iostat)
    call h5tclose_f(c_type,iostat)

  end subroutine write_parallel_params
!!!=================================================================================================
  !> Public 1D REAL store wrapper.
  !!
  !! 1D data are not spatially decomposed by this backend, so the operation is forwarded to
  !! hdf5_single after rejecting MPI-I/O-only file types.
  !!
  !! @param data Array to store.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  !! @param file_overwrite_optin Optional override controlling whether an existing file is
  !! overwritten.
  subroutine store_1d_rk(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:), intent(in)           :: data
    character(len=*)           , intent(in)           :: type
    character(len=*)           , intent(in), optional :: fname_optin
    integer                    , intent(in), optional :: nt_optin,ns_optin
    character(len=*)           , intent(in), optional :: dset_name_optin
    logical                    , intent(in), optional :: file_overwrite_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical                                           :: mpi_io
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Determine the backend policy; fname_optin preserves the legacy serial path without evaluating type.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)
    ! Rank-1/2 arrays have no spatial parallel implementation in this module.
    call reject_mpiio_for_serial_rank("store_1d_rk",type,mpi_io)
    ! Preserve the original serial/per-rank semantics and filename handling.
    call store_single(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)

  end subroutine store_1d_rk
!!!=================================================================================================
  !> Public 1D INTEGER store wrapper.
  !!
  !! The data are written by hdf5_single; the policy guard prevents shared-file MPI-I/O types from
  !! reaching a serial/per-rank implementation.
  !!
  !! @param data Array to store.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  !! @param file_overwrite_optin Optional override controlling whether an existing file is
  !! overwritten.
  subroutine store_1d_i(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:), intent(in)           :: data
    character(len=*)     , intent(in)           :: type
    character(len=*)     , intent(in), optional :: fname_optin
    integer              , intent(in), optional :: nt_optin,ns_optin
    character(len=*)     , intent(in), optional :: dset_name_optin
    logical              , intent(in), optional :: file_overwrite_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    logical                                     :: mpi_io
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Determine the backend policy; fname_optin preserves the legacy serial path without evaluating type.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)
    ! Rank-1/2 arrays have no spatial parallel implementation in this module.
    call reject_mpiio_for_serial_rank("store_1d_i",type,mpi_io)
    ! Preserve the original serial/per-rank semantics and filename handling.
    call store_single(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)

  end subroutine store_1d_i
!!!=================================================================================================
  !> Public 2D REAL store wrapper.
  !!
  !! 2D data remain on the serial/per-rank path and are forwarded unchanged to hdf5_single after the
  !! MPI-I/O policy check.
  !!
  !! @param data Array to store.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  !! @param file_overwrite_optin Optional override controlling whether an existing file is
  !! overwritten.
  subroutine store_2d_rk(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(in) :: data
    character(len=*), intent(in) :: type
    character(len=*), intent(in), optional :: fname_optin
    integer, intent(in), optional :: nt_optin,ns_optin
    character(len=*), intent(in), optional :: dset_name_optin
    logical, intent(in), optional :: file_overwrite_optin
    logical :: mpi_io

    ! Determine the backend policy; fname_optin preserves the legacy serial path without evaluating type.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)
    ! Rank-1/2 arrays have no spatial parallel implementation in this module.
    call reject_mpiio_for_serial_rank("store_2d_rk",type,mpi_io)
    ! Preserve the original serial/per-rank semantics and filename handling.
    call store_single(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
  end subroutine store_2d_rk
!!!=================================================================================================
  !> Public 3D REAL store dispatcher.
  !!
  !! Spatial MPI-I/O types use store_3d_rk_parallel; all other types preserve the legacy hdf5_single
  !! behavior.
  !!
  !! @param data Array to store.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  !! @param file_overwrite_optin Optional override controlling whether an existing file is
  !! overwritten.
  subroutine store_3d_rk(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(in) :: data
    character(len=*), intent(in) :: type
    character(len=*), intent(in), optional :: fname_optin
    integer, intent(in), optional :: nt_optin,ns_optin
    character(len=*), intent(in), optional :: dset_name_optin
    logical, intent(in), optional :: file_overwrite_optin
    logical :: mpi_io

    ! Without fname_optin, backend selection is based on the logical file type, not merely on array rank.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)

    ! Shared-file types use Parallel HDF5; all other types retain the legacy hdf5_single path.
    if(mpi_io)then
       call store_3d_rk_parallel(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
    else
       call store_single(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
    end if
  end subroutine store_3d_rk
!!!=================================================================================================
  !> Public 4D INTEGER store dispatcher.
  !!
  !! The first three dimensions are considered spatial on the parallel path; non-MPI-I/O types are
  !! delegated to hdf5_single.
  !!
  !! @param data Array to store.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  !! @param file_overwrite_optin Optional override controlling whether an existing file is
  !! overwritten.
  subroutine store_4d_i(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:,:,:,:), intent(in) :: data
    character(len=*), intent(in) :: type
    character(len=*), intent(in), optional :: fname_optin
    integer, intent(in), optional :: nt_optin,ns_optin
    character(len=*), intent(in), optional :: dset_name_optin
    logical, intent(in), optional :: file_overwrite_optin
    logical :: mpi_io

    ! Without fname_optin, backend selection is based on the logical file type, not merely on array rank.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)

    ! Shared-file types use Parallel HDF5; all other types retain the legacy hdf5_single path.
    if(mpi_io)then
       call store_4d_i_parallel(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
    else
       call store_single(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
    end if
  end subroutine store_4d_i
!!!=================================================================================================
  !> Public 4D REAL store dispatcher.
  !!
  !! The first three dimensions are spatial on the parallel path. The fourth dimension is stored in
  !! full by every rank and is not decomposed.
  !!
  !! @param data Array to store.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  !! @param file_overwrite_optin Optional override controlling whether an existing file is
  !! overwritten.
  subroutine store_4d_rk(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(in) :: data
    character(len=*), intent(in) :: type
    character(len=*), intent(in), optional :: fname_optin
    integer, intent(in), optional :: nt_optin,ns_optin
    character(len=*), intent(in), optional :: dset_name_optin
    logical, intent(in), optional :: file_overwrite_optin
    logical :: mpi_io

    ! Without fname_optin, backend selection is based on the logical file type, not merely on array rank.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)

    ! Shared-file types use Parallel HDF5; all other types retain the legacy hdf5_single path.
    if(mpi_io)then
       call store_4d_rk_parallel(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
    else
       call store_single(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
    end if
  end subroutine store_4d_rk
!!!=================================================================================================
  !> Public 5D REAL store dispatcher.
  !!
  !! The first three dimensions are spatially decomposed; dimensions four and five are local full
  !! dimensions shared by the rank-local buffer.
  !!
  !! @param data Array to store.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  !! @param file_overwrite_optin Optional override controlling whether an existing file is
  !! overwritten.
  subroutine store_5d_rk(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:), intent(in) :: data
    character(len=*), intent(in) :: type
    character(len=*), intent(in), optional :: fname_optin
    integer, intent(in), optional :: nt_optin,ns_optin
    character(len=*), intent(in), optional :: dset_name_optin
    logical, intent(in), optional :: file_overwrite_optin
    logical :: mpi_io

    ! Without fname_optin, backend selection is based on the logical file type, not merely on array rank.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)

    ! Shared-file types use Parallel HDF5; all other types retain the legacy hdf5_single path.
    if(mpi_io)then
       call store_5d_rk_parallel(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
    else
       call store_single(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin,file_overwrite_optin)
    end if
  end subroutine store_5d_rk
!!!=================================================================================================
  !> Write a spatially decomposed 3D REAL(kind=rk) field into one shared HDF5 dataset.
  !!
  !! Each rank selects its global 3D hyperslab, while dataset dimensions are obtained from the
  !! maximum global upper bounds across block_comm. Dataset storage precision follows rk unless
  !! mpiio_single_precision requests float32.
  !!
  !! @param data Rank-local array written to the shared HDF5 dataset.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  !! @param file_overwrite_optin Optional override controlling whether an existing file is
  !! overwritten.
  subroutine store_3d_rk_parallel(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin, file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(in) :: data
    character(len=*), intent(in) :: type
    character(len=*), intent(in), optional :: fname_optin
    integer, intent(in), optional :: nt_optin
    integer, intent(in), optional :: ns_optin
    character(len=*), intent(in), optional :: dset_name_optin
    logical, intent(in), optional :: file_overwrite_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t) :: file_id,dset_id,filespace,memspace,dxpl_id
    integer(hid_t) :: mem_type,file_type
    integer(hsize_t), dimension(3) :: dims_local,dims_file,count,stride
    integer(hssize_t), dimension(3) :: offset
    integer, dimension(3) :: lbounds,ubounds,dims_local_i,dims_file_i
    integer, dimension(3) :: lbounds_space,ubounds_space
    integer :: iostat,mpierr
    logical :: write_params,file_overwrite,append_dataset
    character(len=max_length_fname) :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Match hdf5_single: unless requested explicitly, a named dataset is appended to an existing
    ! file while the default unnamed dataset creates/truncates its target file.
    file_overwrite = .false.

    ! Resolve the shared filename.  This also verifies that type is actually configured for MPI-I/O.
    call get_parallel_fname(fname,type,fname_optin,nt_optin,ns_optin)

    ! Optional caller override for the append/truncate behavior.
    if(present(file_overwrite_optin)) file_overwrite = file_overwrite_optin

    ! A custom dataset name denotes the multi-dataset/append use case used by hdf5_single as well.
    append_dataset = present(dset_name_optin)
    if(append_dataset)then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! Obtain the cached global 1-based indices of the spatial cells owned by this rank.
    call get_parallel_spatial_bounds(lbounds_space,ubounds_space)

    ! A 3D field is purely spatial, so its complete local bounds are the decomposition bounds.
    lbounds = lbounds_space
    ubounds = ubounds_space

    ! Record the actual local buffer shape and verify it against the selected global interval.
    dims_local_i = shape(data)

    ! This collective validation catches inconsistent geometry bounds before entering HDF5.
    call check_parallel_dimensions("store_3d_rk_parallel",dims_local_i,lbounds,ubounds)

    ! The decomposition bounds are global 1-based indices.  Their component-wise maximum therefore
    ! gives the required global dataset extent for the spatial dimensions (and identical full
    ! extents for non-spatial dimensions).
    call mpi_allreduce(ubounds,dims_file_i,3,mpi_integer,mpi_max, &
         params%parallelism%block_comm,mpierr)
    if(mpierr.ne.mpi_success) call mpi_hdf5_abort("store_3d_rk_parallel:mpi_allreduce",mpierr)

    ! Convert Fortran default INTEGER extents to HDF5 size kinds and translate 1-based lower bounds
    ! to the zero-based offsets required by HDF5 hyperslabs.
    dims_local = int(dims_local_i,kind=hsize_t)
    dims_file  = int(dims_file_i,kind=hsize_t)
    offset     = int(lbounds-1,kind=hssize_t)
    count      = 1_hsize_t
    stride     = 1_hsize_t

    ! Initialize the HDF5 Fortran interface before creating/opening HDF5 objects.
    call h5open_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_3d_rk_parallel:h5open_f",iostat)

    ! Determine the datatype of the REAL(kind=rk) array in memory and the
    ! datatype to be stored in the HDF5 file. Normally both follow rk. If
    ! mpiio_single_precision is enabled, only file_type becomes IEEE float32
    ! and HDF5 performs the conversion from mem_type during the write.
    call get_hdf5_real_store_types(mem_type,file_type)

    ! Open or create the shared file collectively with the MPIO file driver.
    call open_parallel_file_store(fname,append_dataset,file_overwrite,file_id,write_params)

    ! Create the global dataset dataspace.  Every rank uses the same global dims_file values.
    call h5screate_simple_f(3,dims_file,filespace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_3d_rk_parallel:h5screate_simple_f(filespace)",iostat)

    ! Dataset creation is collective for a file opened with the MPIO driver.
    ! file_type determines the precision physically stored in the HDF5 file.
    ! Create the dataset collectively.  The datatype argument defines the representation stored
    ! in the file; each rank subsequently selects only its own spatial hyperslab.
    call h5dcreate_f(file_id,trim(dset_name),file_type,filespace,dset_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_3d_rk_parallel:h5dcreate_f",iostat)

    ! Select the contiguous rank-owned region in the global file dataspace.
    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offset,count,iostat,stride,dims_local)
    if(iostat.ne.0) call mpi_hdf5_abort("store_3d_rk_parallel:h5sselect_hyperslab_f",iostat)

    ! Describe the contiguous local Fortran buffer independently of the global file dataspace.
    call h5screate_simple_f(3,dims_local,memspace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_3d_rk_parallel:h5screate_simple_f(memspace)",iostat)

    ! Attach the requested collective MPI-I/O transfer mode to the following dataset operation.
    call create_collective_dxpl(dxpl_id)

    ! mem_type describes the REAL(kind=rk) buffer. If file_type differs, HDF5
    ! performs the floating-point conversion as part of the dataset transfer.
    ! Write only this rank's selected hyperslab from its local memory buffer.
    call h5dwrite_f(dset_id,mem_type,data,dims_local,iostat, &
         mem_space_id=memspace,file_space_id=filespace,xfer_prp=dxpl_id)
    if(iostat.ne.0) call mpi_hdf5_abort("store_3d_rk_parallel:h5dwrite_f",iostat)

    ! Release transfer and dataspace objects after the dataset transfer is complete.
    call h5pclose_f(dxpl_id,iostat)
    call h5sclose_f(memspace,iostat)
    call h5sclose_f(filespace,iostat)
    call h5dclose_f(dset_id,iostat)

    ! New files receive the canonical params dataset once; existing append targets keep theirs.
    if(write_params) call write_parallel_params(file_id)

    ! Close the collectively opened file before shutting down the HDF5 Fortran interface.
    call h5fclose_f(file_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_3d_rk_parallel:h5fclose_f",iostat)

    ! Release HDF5 library-level resources associated with this operation.
    call h5close_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_3d_rk_parallel:h5close_f",iostat)

  end subroutine store_3d_rk_parallel
!!!=================================================================================================
  !> Write a spatially decomposed 4D INTEGER field into one shared HDF5 dataset.
  !!
  !! The first three dimensions use the geometry decomposition; the fourth dimension is stored
  !! completely and is not decomposed.
  !!
  !! @param data Rank-local array written to the shared HDF5 dataset.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  !! @param file_overwrite_optin Optional override controlling whether an existing file is
  !! overwritten.
  subroutine store_4d_i_parallel(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin, file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:,:,:,:), intent(in) :: data
    character(len=*), intent(in) :: type
    character(len=*), intent(in), optional :: fname_optin
    integer, intent(in), optional :: nt_optin
    integer, intent(in), optional :: ns_optin
    character(len=*), intent(in), optional :: dset_name_optin
    logical, intent(in), optional :: file_overwrite_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t) :: file_id,dset_id,filespace,memspace,dxpl_id
    integer(hsize_t), dimension(4) :: dims_local,dims_file,count,stride
    integer(hssize_t), dimension(4) :: offset
    integer, dimension(4) :: lbounds,ubounds,dims_local_i,dims_file_i
    integer, dimension(3) :: lbounds_space,ubounds_space
    integer :: iostat,mpierr
    logical :: write_params,file_overwrite,append_dataset
    character(len=max_length_fname) :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Match hdf5_single: unless requested explicitly, a named dataset is appended to an existing
    ! file while the default unnamed dataset creates/truncates its target file.
    file_overwrite = .false.

    ! Resolve the shared filename.  This also verifies that type is actually configured for MPI-I/O.
    call get_parallel_fname(fname,type,fname_optin,nt_optin,ns_optin)

    ! Optional caller override for the append/truncate behavior.
    if(present(file_overwrite_optin)) file_overwrite = file_overwrite_optin

    ! A custom dataset name denotes the multi-dataset/append use case used by hdf5_single as well.
    append_dataset = present(dset_name_optin)
    if(append_dataset)then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! Obtain the cached global 1-based indices of the spatial cells owned by this rank.
    call get_parallel_spatial_bounds(lbounds_space,ubounds_space)

    ! Dimensions 1:3 are spatially decomposed.  Dimension 4 starts at one and is present in full
    ! in every rank-local buffer.
    lbounds = (/lbounds_space,1/)
    ubounds = (/ubounds_space,size(data,4)/)

    ! Record the actual local buffer shape and verify it against the selected global interval.
    dims_local_i = shape(data)

    ! This collective validation catches inconsistent geometry bounds before entering HDF5.
    call check_parallel_dimensions("store_4d_i_parallel",dims_local_i,lbounds,ubounds)

    ! The decomposition bounds are global 1-based indices.  Their component-wise maximum therefore
    ! gives the required global dataset extent for the spatial dimensions (and identical full
    ! extents for non-spatial dimensions).
    call mpi_allreduce(ubounds,dims_file_i,4,mpi_integer,mpi_max, &
         params%parallelism%block_comm,mpierr)
    if(mpierr.ne.mpi_success) call mpi_hdf5_abort("store_4d_i_parallel:mpi_allreduce",mpierr)

    ! Convert Fortran default INTEGER extents to HDF5 size kinds and translate 1-based lower bounds
    ! to the zero-based offsets required by HDF5 hyperslabs.
    dims_local = int(dims_local_i,kind=hsize_t)
    dims_file  = int(dims_file_i,kind=hsize_t)
    offset     = int(lbounds-1,kind=hssize_t)
    count      = 1_hsize_t
    stride     = 1_hsize_t

    ! Initialize the HDF5 Fortran interface before creating/opening HDF5 objects.
    call h5open_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_i_parallel:h5open_f",iostat)

    ! Open or create the shared file collectively with the MPIO file driver.
    call open_parallel_file_store(fname,append_dataset,file_overwrite,file_id,write_params)

    ! Create the global dataset dataspace.  Every rank uses the same global dims_file values.
    call h5screate_simple_f(4,dims_file,filespace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_i_parallel:h5screate_simple_f(filespace)",iostat)

    ! Dataset creation is collective for a file opened with the MPIO driver.
    ! Create the dataset collectively.  The datatype argument defines the representation stored
    ! in the file; each rank subsequently selects only its own spatial hyperslab.
    call h5dcreate_f(file_id,trim(dset_name),h5t_native_integer,filespace,dset_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_i_parallel:h5dcreate_f",iostat)

    ! Select the contiguous rank-owned region in the global file dataspace.
    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offset,count,iostat,stride,dims_local)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_i_parallel:h5sselect_hyperslab_f",iostat)

    ! Describe the contiguous local Fortran buffer independently of the global file dataspace.
    call h5screate_simple_f(4,dims_local,memspace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_i_parallel:h5screate_simple_f(memspace)",iostat)

    ! Attach the requested collective MPI-I/O transfer mode to the following dataset operation.
    call create_collective_dxpl(dxpl_id)

    ! Write only this rank's selected hyperslab from its local memory buffer.
    call h5dwrite_f(dset_id,h5t_native_integer,data,dims_local,iostat, &
         mem_space_id=memspace,file_space_id=filespace,xfer_prp=dxpl_id)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_i_parallel:h5dwrite_f",iostat)

    ! Release transfer and dataspace objects after the dataset transfer is complete.
    call h5pclose_f(dxpl_id,iostat)
    call h5sclose_f(memspace,iostat)
    call h5sclose_f(filespace,iostat)
    call h5dclose_f(dset_id,iostat)

    ! New files receive the canonical params dataset once; existing append targets keep theirs.
    if(write_params) call write_parallel_params(file_id)

    ! Close the collectively opened file before shutting down the HDF5 Fortran interface.
    call h5fclose_f(file_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_i_parallel:h5fclose_f",iostat)

    ! Release HDF5 library-level resources associated with this operation.
    call h5close_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_i_parallel:h5close_f",iostat)

  end subroutine store_4d_i_parallel
!!!=================================================================================================
  !> Write a spatially decomposed 4D REAL(kind=rk) field into one shared HDF5 dataset.
  !!
  !! The first three dimensions use the geometry decomposition and dimension four is taken directly
  !! from size(data,4). REAL file precision follows rk unless the optional MPI-I/O single-precision
  !! storage policy is enabled.
  !!
  !! @param data Rank-local array written to the shared HDF5 dataset.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  !! @param file_overwrite_optin Optional override controlling whether an existing file is
  !! overwritten.
  subroutine store_4d_rk_parallel(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin, file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(in) :: data
    character(len=*), intent(in) :: type
    character(len=*), intent(in), optional :: fname_optin
    integer, intent(in), optional :: nt_optin
    integer, intent(in), optional :: ns_optin
    character(len=*), intent(in), optional :: dset_name_optin
    logical, intent(in), optional :: file_overwrite_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t) :: file_id,dset_id,filespace,memspace,dxpl_id
    integer(hid_t) :: mem_type,file_type
    integer(hsize_t), dimension(4) :: dims_local,dims_file,count,stride
    integer(hssize_t), dimension(4) :: offset
    integer, dimension(4) :: lbounds,ubounds,dims_local_i,dims_file_i
    integer, dimension(3) :: lbounds_space,ubounds_space
    integer :: iostat,mpierr
    logical :: write_params,file_overwrite,append_dataset
    character(len=max_length_fname) :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Match hdf5_single: unless requested explicitly, a named dataset is appended to an existing
    ! file while the default unnamed dataset creates/truncates its target file.
    file_overwrite = .false.

    ! Resolve the shared filename.  This also verifies that type is actually configured for MPI-I/O.
    call get_parallel_fname(fname,type,fname_optin,nt_optin,ns_optin)

    ! Optional caller override for the append/truncate behavior.
    if(present(file_overwrite_optin)) file_overwrite = file_overwrite_optin

    ! A custom dataset name denotes the multi-dataset/append use case used by hdf5_single as well.
    append_dataset = present(dset_name_optin)
    if(append_dataset)then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! Obtain the cached global 1-based indices of the spatial cells owned by this rank.
    call get_parallel_spatial_bounds(lbounds_space,ubounds_space)

    ! Dimensions 1:3 are spatially decomposed.  Dimension 4 starts at one and is present in full
    ! in every rank-local buffer.
    lbounds = (/lbounds_space,1/)
    ubounds = (/ubounds_space,size(data,4)/)

    ! Record the actual local buffer shape and verify it against the selected global interval.
    dims_local_i = shape(data)

    ! This collective validation catches inconsistent geometry bounds before entering HDF5.
    call check_parallel_dimensions("store_4d_rk_parallel",dims_local_i,lbounds,ubounds)

    ! The decomposition bounds are global 1-based indices.  Their component-wise maximum therefore
    ! gives the required global dataset extent for the spatial dimensions (and identical full
    ! extents for non-spatial dimensions).
    call mpi_allreduce(ubounds,dims_file_i,4,mpi_integer,mpi_max, &
         params%parallelism%block_comm,mpierr)
    if(mpierr.ne.mpi_success) call mpi_hdf5_abort("store_4d_rk_parallel:mpi_allreduce",mpierr)

    ! Convert Fortran default INTEGER extents to HDF5 size kinds and translate 1-based lower bounds
    ! to the zero-based offsets required by HDF5 hyperslabs.
    dims_local = int(dims_local_i,kind=hsize_t)
    dims_file  = int(dims_file_i,kind=hsize_t)
    offset     = int(lbounds-1,kind=hssize_t)
    count      = 1_hsize_t
    stride     = 1_hsize_t

    ! Initialize the HDF5 Fortran interface before creating/opening HDF5 objects.
    call h5open_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_rk_parallel:h5open_f",iostat)

    ! Determine the datatype of the REAL(kind=rk) array in memory and the
    ! datatype to be stored in the HDF5 file. Normally both follow rk. If
    ! mpiio_single_precision is enabled, only file_type becomes IEEE float32
    ! and HDF5 performs the conversion from mem_type during the write.
    call get_hdf5_real_store_types(mem_type,file_type)

    ! Open or create the shared file collectively with the MPIO file driver.
    call open_parallel_file_store(fname,append_dataset,file_overwrite,file_id,write_params)

    ! Create the global dataset dataspace.  Every rank uses the same global dims_file values.
    call h5screate_simple_f(4,dims_file,filespace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_rk_parallel:h5screate_simple_f(filespace)",iostat)

    ! Dataset creation is collective for a file opened with the MPIO driver.
    ! file_type determines the precision physically stored in the HDF5 file.
    ! Create the dataset collectively.  The datatype argument defines the representation stored
    ! in the file; each rank subsequently selects only its own spatial hyperslab.
    call h5dcreate_f(file_id,trim(dset_name),file_type,filespace,dset_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_rk_parallel:h5dcreate_f",iostat)

    ! Select the contiguous rank-owned region in the global file dataspace.
    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offset,count,iostat,stride,dims_local)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_rk_parallel:h5sselect_hyperslab_f",iostat)

    ! Describe the contiguous local Fortran buffer independently of the global file dataspace.
    call h5screate_simple_f(4,dims_local,memspace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_rk_parallel:h5screate_simple_f(memspace)",iostat)

    ! Attach the requested collective MPI-I/O transfer mode to the following dataset operation.
    call create_collective_dxpl(dxpl_id)

    ! mem_type describes the REAL(kind=rk) buffer. If file_type differs, HDF5
    ! performs the floating-point conversion as part of the dataset transfer.
    ! Write only this rank's selected hyperslab from its local memory buffer.
    call h5dwrite_f(dset_id,mem_type,data,dims_local,iostat, &
         mem_space_id=memspace,file_space_id=filespace,xfer_prp=dxpl_id)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_rk_parallel:h5dwrite_f",iostat)

    ! Release transfer and dataspace objects after the dataset transfer is complete.
    call h5pclose_f(dxpl_id,iostat)
    call h5sclose_f(memspace,iostat)
    call h5sclose_f(filespace,iostat)
    call h5dclose_f(dset_id,iostat)

    ! New files receive the canonical params dataset once; existing append targets keep theirs.
    if(write_params) call write_parallel_params(file_id)

    ! Close the collectively opened file before shutting down the HDF5 Fortran interface.
    call h5fclose_f(file_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_rk_parallel:h5fclose_f",iostat)

    ! Release HDF5 library-level resources associated with this operation.
    call h5close_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_4d_rk_parallel:h5close_f",iostat)

  end subroutine store_4d_rk_parallel
!!!=================================================================================================
  !> Write a spatially decomposed 5D REAL(kind=rk) field into one shared HDF5 dataset.
  !!
  !! The first three dimensions are spatially decomposed; dimensions four and five are stored in
  !! full for each spatial hyperslab.
  !!
  !! @param data Rank-local array written to the shared HDF5 dataset.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  !! @param file_overwrite_optin Optional override controlling whether an existing file is
  !! overwritten.
  subroutine store_5d_rk_parallel(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin, file_overwrite_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:), intent(in) :: data
    character(len=*), intent(in) :: type
    character(len=*), intent(in), optional :: fname_optin
    integer, intent(in), optional :: nt_optin
    integer, intent(in), optional :: ns_optin
    character(len=*), intent(in), optional :: dset_name_optin
    logical, intent(in), optional :: file_overwrite_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t) :: file_id,dset_id,filespace,memspace,dxpl_id
    integer(hid_t) :: mem_type,file_type
    integer(hsize_t), dimension(5) :: dims_local,dims_file,count,stride
    integer(hssize_t), dimension(5) :: offset
    integer, dimension(5) :: lbounds,ubounds,dims_local_i,dims_file_i
    integer, dimension(3) :: lbounds_space,ubounds_space
    integer :: iostat,mpierr
    logical :: write_params,file_overwrite,append_dataset
    character(len=max_length_fname) :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Match hdf5_single: unless requested explicitly, a named dataset is appended to an existing
    ! file while the default unnamed dataset creates/truncates its target file.
    file_overwrite = .false.

    ! Resolve the shared filename.  This also verifies that type is actually configured for MPI-I/O.
    call get_parallel_fname(fname,type,fname_optin,nt_optin,ns_optin)

    ! Optional caller override for the append/truncate behavior.
    if(present(file_overwrite_optin)) file_overwrite = file_overwrite_optin

    ! A custom dataset name denotes the multi-dataset/append use case used by hdf5_single as well.
    append_dataset = present(dset_name_optin)
    if(append_dataset)then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! Obtain the cached global 1-based indices of the spatial cells owned by this rank.
    call get_parallel_spatial_bounds(lbounds_space,ubounds_space)

    ! Dimensions 1:3 are spatially decomposed; dimensions 4:5 are non-spatial and start at one.
    lbounds = (/lbounds_space,1,1/)
    ubounds = (/ubounds_space,size(data,4),size(data,5)/)

    ! Record the actual local buffer shape and verify it against the selected global interval.
    dims_local_i = shape(data)

    ! This collective validation catches inconsistent geometry bounds before entering HDF5.
    call check_parallel_dimensions("store_5d_rk_parallel",dims_local_i,lbounds,ubounds)

    ! The decomposition bounds are global 1-based indices.  Their component-wise maximum therefore
    ! gives the required global dataset extent for the spatial dimensions (and identical full
    ! extents for non-spatial dimensions).
    call mpi_allreduce(ubounds,dims_file_i,5,mpi_integer,mpi_max, &
         params%parallelism%block_comm,mpierr)
    if(mpierr.ne.mpi_success) call mpi_hdf5_abort("store_5d_rk_parallel:mpi_allreduce",mpierr)

    ! Convert Fortran default INTEGER extents to HDF5 size kinds and translate 1-based lower bounds
    ! to the zero-based offsets required by HDF5 hyperslabs.
    dims_local = int(dims_local_i,kind=hsize_t)
    dims_file  = int(dims_file_i,kind=hsize_t)
    offset     = int(lbounds-1,kind=hssize_t)
    count      = 1_hsize_t
    stride     = 1_hsize_t

    ! Initialize the HDF5 Fortran interface before creating/opening HDF5 objects.
    call h5open_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_5d_rk_parallel:h5open_f",iostat)

    ! Determine the datatype of the REAL(kind=rk) array in memory and the
    ! datatype to be stored in the HDF5 file. Normally both follow rk. If
    ! mpiio_single_precision is enabled, only file_type becomes IEEE float32
    ! and HDF5 performs the conversion from mem_type during the write.
    call get_hdf5_real_store_types(mem_type,file_type)

    ! Open or create the shared file collectively with the MPIO file driver.
    call open_parallel_file_store(fname,append_dataset,file_overwrite,file_id,write_params)

    ! Create the global dataset dataspace.  Every rank uses the same global dims_file values.
    call h5screate_simple_f(5,dims_file,filespace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_5d_rk_parallel:h5screate_simple_f(filespace)",iostat)

    ! Dataset creation is collective for a file opened with the MPIO driver.
    ! file_type determines the precision physically stored in the HDF5 file.
    ! Create the dataset collectively.  The datatype argument defines the representation stored
    ! in the file; each rank subsequently selects only its own spatial hyperslab.
    call h5dcreate_f(file_id,trim(dset_name),file_type,filespace,dset_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_5d_rk_parallel:h5dcreate_f",iostat)

    ! Select the contiguous rank-owned region in the global file dataspace.
    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offset,count,iostat,stride,dims_local)
    if(iostat.ne.0) call mpi_hdf5_abort("store_5d_rk_parallel:h5sselect_hyperslab_f",iostat)

    ! Describe the contiguous local Fortran buffer independently of the global file dataspace.
    call h5screate_simple_f(5,dims_local,memspace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_5d_rk_parallel:h5screate_simple_f(memspace)",iostat)

    ! Attach the requested collective MPI-I/O transfer mode to the following dataset operation.
    call create_collective_dxpl(dxpl_id)

    ! mem_type describes the REAL(kind=rk) buffer. If file_type differs, HDF5
    ! performs the floating-point conversion as part of the dataset transfer.
    ! Write only this rank's selected hyperslab from its local memory buffer.
    call h5dwrite_f(dset_id,mem_type,data,dims_local,iostat, &
         mem_space_id=memspace,file_space_id=filespace,xfer_prp=dxpl_id)
    if(iostat.ne.0) call mpi_hdf5_abort("store_5d_rk_parallel:h5dwrite_f",iostat)

    ! Release transfer and dataspace objects after the dataset transfer is complete.
    call h5pclose_f(dxpl_id,iostat)
    call h5sclose_f(memspace,iostat)
    call h5sclose_f(filespace,iostat)
    call h5dclose_f(dset_id,iostat)

    ! New files receive the canonical params dataset once; existing append targets keep theirs.
    if(write_params) call write_parallel_params(file_id)

    ! Close the collectively opened file before shutting down the HDF5 Fortran interface.
    call h5fclose_f(file_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_5d_rk_parallel:h5fclose_f",iostat)

    ! Release HDF5 library-level resources associated with this operation.
    call h5close_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("store_5d_rk_parallel:h5close_f",iostat)

  end subroutine store_5d_rk_parallel
!!!=================================================================================================
  !> Public 1D REAL load wrapper.
  !!
  !! 1D data stay on the hdf5_single path; MPI-I/O file types are rejected because they require a
  !! shared spatial backend.
  !!
  !! @param data Array receiving the loaded data.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  subroutine load_1d_rk(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:), intent(out) :: data
    character(len=*), intent(in)             :: type
    character(len=*), intent(in), optional   :: fname_optin
    integer, intent(in), optional            :: nt_optin,ns_optin
    character(len=*), intent(in), optional   :: dset_name_optin
    logical                                  :: mpi_io

    ! Determine the backend policy; fname_optin preserves the legacy serial path without evaluating type.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)
    ! Rank-1/2 arrays have no spatial parallel implementation in this module.
    call reject_mpiio_for_serial_rank("load_1d_rk",type,mpi_io)
    ! Preserve the original serial/per-rank load implementation.
    call load_single(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
  end subroutine load_1d_rk
!!!=================================================================================================
  !> Public 1D INTEGER load wrapper.
  !!
  !! The operation is delegated unchanged to hdf5_single after validating that the selected file
  !! type is not configured for MPI-I/O.
  !!
  !! @param data Array receiving the loaded data.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  subroutine load_1d_i(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:), intent(out)     :: data
    character(len=*), intent(in)           :: type
    character(len=*), intent(in), optional :: fname_optin
    integer, intent(in), optional          :: nt_optin,ns_optin
    character(len=*), intent(in), optional :: dset_name_optin
    logical                                :: mpi_io

    ! Determine the backend policy; fname_optin preserves the legacy serial path without evaluating type.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)
    ! Rank-1/2 arrays have no spatial parallel implementation in this module.
    call reject_mpiio_for_serial_rank("load_1d_i",type,mpi_io)
    ! Preserve the original serial/per-rank load implementation.
    call load_single(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
  end subroutine load_1d_i
!!!=================================================================================================
  !> Public 2D REAL load wrapper.
  !!
  !! 2D data are loaded by hdf5_single; this module only performs the policy guard.
  !!
  !! @param data Array receiving the loaded data.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  subroutine load_2d_rk(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:), intent(out) :: data
    character(len=*), intent(in)               :: type
    character(len=*), intent(in), optional     :: fname_optin
    integer, intent(in), optional              :: nt_optin,ns_optin
    character(len=*), intent(in), optional     :: dset_name_optin
    logical                                    :: mpi_io

    ! Determine the backend policy; fname_optin preserves the legacy serial path without evaluating type.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)
    ! Rank-1/2 arrays have no spatial parallel implementation in this module.
    call reject_mpiio_for_serial_rank("load_2d_rk",type,mpi_io)
    ! Preserve the original serial/per-rank load implementation.
    call load_single(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
  end subroutine load_2d_rk
!!!=================================================================================================
  !> Public 3D REAL load dispatcher.
  !!
  !! Spatial MPI-I/O types read their local hyperslab from the shared file; other types use the
  !! original hdf5_single implementation.
  !!
  !! @param data Array receiving the loaded data.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  subroutine load_3d_rk(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(out) :: data
    character(len=*), intent(in)                 :: type
    character(len=*), intent(in), optional       :: fname_optin
    integer, intent(in), optional                :: nt_optin,ns_optin
    character(len=*), intent(in), optional       :: dset_name_optin
    logical                                      :: mpi_io

    ! Without fname_optin, backend selection is based on the logical file type, not merely on array rank.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)

    ! Shared-file types use Parallel HDF5; all other types retain the legacy hdf5_single path.
    if(mpi_io)then
       call load_3d_rk_parallel(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
    else
       call load_single(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
    end if
  end subroutine load_3d_rk
!!!=================================================================================================
  !> Public 4D INTEGER load dispatcher.
  !!
  !! MPI-I/O types use the shared-file backend with spatial decomposition in the first three
  !! dimensions; other types use hdf5_single.
  !!
  !! @param data Array receiving the loaded data.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  subroutine load_4d_i(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:,:,:,:), intent(out) :: data
    character(len=*), intent(in)             :: type
    character(len=*), intent(in), optional   :: fname_optin
    integer, intent(in), optional            :: nt_optin,ns_optin
    character(len=*), intent(in), optional   :: dset_name_optin
    logical                                  :: mpi_io

    ! Without fname_optin, backend selection is based on the logical file type, not merely on array rank.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)

    ! Shared-file types use Parallel HDF5; all other types retain the legacy hdf5_single path.
    if(mpi_io)then
       call load_4d_i_parallel(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
    else
       call load_single(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
    end if
  end subroutine load_4d_i
!!!=================================================================================================
  !> Public 4D REAL load dispatcher.
  !!
  !! Parallel loads use the first three dimensions as spatial dimensions and retain the complete
  !! fourth dimension of the target array.
  !!
  !! @param data Array receiving the loaded data.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  subroutine load_4d_rk(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(out) :: data
    character(len=*), intent(in)                   :: type
    character(len=*), intent(in), optional         :: fname_optin
    integer, intent(in), optional                  :: nt_optin,ns_optin
    character(len=*), intent(in), optional         :: dset_name_optin
    logical                                        :: mpi_io

    ! Without fname_optin, backend selection is based on the logical file type, not merely on array rank.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)

    ! Shared-file types use Parallel HDF5; all other types retain the legacy hdf5_single path.
    if(mpi_io)then
       call load_4d_rk_parallel(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
    else
       call load_single(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
    end if
  end subroutine load_4d_rk
!!!=================================================================================================
  !> Public 5D REAL load dispatcher.
  !!
  !! Parallel loads read the rank-owned spatial hyperslab and all local entries in dimensions four
  !! and five.
  !!
  !! @param data Array receiving the loaded data.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  subroutine load_5d_rk(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:), intent(out) :: data
    character(len=*), intent(in)                     :: type
    character(len=*), intent(in), optional           :: fname_optin
    integer, intent(in), optional                    :: nt_optin,ns_optin
    character(len=*), intent(in), optional           :: dset_name_optin
    logical                                          :: mpi_io

    ! Without fname_optin, backend selection is based on the logical file type, not merely on array rank.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)

    ! Shared-file types use Parallel HDF5; all other types retain the legacy hdf5_single path.
    if(mpi_io)then
       call load_5d_rk_parallel(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
    else
       call load_single(data,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
    end if
  end subroutine load_5d_rk
!!!=================================================================================================
  !> Read a rank-local 5D REAL(kind=rk) hyperslab from a shared Parallel-HDF5 dataset.
  !!
  !! The stored dataset datatype may differ from rk; HDF5 converts it to the rk-derived memory
  !! datatype during H5Dread.
  !!
  !! @param data Array receiving the rank-local data read from the HDF5 dataset.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  subroutine load_5d_rk_parallel(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:,:), intent(out) :: data
    character(len=*), intent(in)                     :: type
    character(len=*), intent(in), optional           :: fname_optin
    integer, intent(in), optional                    :: nt_optin
    integer, intent(in), optional                    :: ns_optin
    character(len=*), intent(in), optional           :: dset_name_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t)                  :: file_id,dset_id,filespace,memspace,dxpl_id
    integer(hid_t)                  :: mem_type
    integer(hsize_t), dimension(5)  :: dims_local,dims_file,dims_max_file,count,stride
    integer(hssize_t), dimension(5) :: offset
    integer, dimension(5)           :: lbounds,ubounds,dims_local_i
    integer, dimension(3)           :: lbounds_space,ubounds_space
    integer                         :: iostat
    character(len=max_length_fname) :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Resolve the shared filename.  This also verifies that type is actually configured for MPI-I/O.
    call get_parallel_fname(fname,type,fname_optin,nt_optin,ns_optin)

    ! Use the caller-supplied dataset name when present; otherwise the conventional name is data.
    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! Obtain the cached global 1-based indices of the spatial cells owned by this rank.
    call get_parallel_spatial_bounds(lbounds_space,ubounds_space)

    ! Dimensions 1:3 are spatially decomposed; dimensions 4:5 are non-spatial and start at one.
    lbounds = (/lbounds_space,1,1/)
    ubounds = (/ubounds_space,size(data,4),size(data,5)/)

    ! Record the actual local buffer shape and verify it against the selected global interval.
    dims_local_i = shape(data)

    ! This collective validation catches inconsistent geometry bounds before entering HDF5.
    call check_parallel_dimensions("load_5d_rk_parallel",dims_local_i,lbounds,ubounds)

    ! Convert Fortran default INTEGER extents to HDF5 size kinds and translate 1-based lower bounds
    ! to the zero-based offsets required by HDF5 hyperslabs.
    dims_local = int(dims_local_i,kind=hsize_t)
    offset     = int(lbounds-1,kind=hssize_t)
    count      = 1_hsize_t
    stride     = 1_hsize_t

    ! Initialize the HDF5 Fortran interface before creating/opening HDF5 objects.
    call h5open_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_5d_rk_parallel:h5open_f",iostat)

    ! The destination array is REAL(kind=rk), therefore the HDF5 memory datatype
    ! always follows rk. The datatype stored in the file is obtained internally
    ! from the dataset. If the file contains float32 and rk is float64, HDF5
    ! automatically converts float32 -> float64 during h5dread_f.
    call get_hdf5_real_memory_type(mem_type)

    ! Open the shared file collectively in read-only mode with the MPIO file driver.
    call open_parallel_file_read(fname,file_id)

    ! Dataset open and dataspace query are performed by all ranks.
    ! Open the existing dataset on every rank; shared metadata operations remain synchronized.
    call h5dopen_f(file_id,trim(dset_name),dset_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_5d_rk_parallel:h5dopen_f",iostat)

    ! Obtain the dataset's file dataspace so the rank can select its hyperslab.
    call h5dget_space_f(dset_id,filespace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_5d_rk_parallel:h5dget_space_f",iostat)

    ! Read the actual on-disk extent instead of assuming it from the current parameter set.
    call h5sget_simple_extent_dims_f(filespace,dims_file,dims_max_file,iostat)
    if(iostat.lt.0) call mpi_hdf5_abort("load_5d_rk_parallel:h5sget_simple_extent_dims_f",iostat)

    ! Refuse to issue an HDF5 read if any rank's requested region lies outside the stored dataset.
    call check_parallel_hyperslab("load_5d_rk_parallel",offset,dims_local,dims_file)

    ! Select the contiguous rank-owned region in the global file dataspace.
    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offset,count,iostat,stride,dims_local)
    if(iostat.ne.0) call mpi_hdf5_abort("load_5d_rk_parallel:h5sselect_hyperslab_f",iostat)

    ! Describe the contiguous local Fortran buffer independently of the global file dataspace.
    call h5screate_simple_f(5,dims_local,memspace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_5d_rk_parallel:h5screate_simple_f(memspace)",iostat)

    ! Attach the requested collective MPI-I/O transfer mode to the following dataset operation.
    call create_collective_dxpl(dxpl_id)

    ! mem_type describes the REAL(kind=rk) destination buffer. The precision
    ! used when the dataset was written does not need to be known by this code.
    ! Read only this rank's selected hyperslab into the local memory buffer.
    call h5dread_f(dset_id,mem_type,data,dims_local,iostat, &
         mem_space_id=memspace,file_space_id=filespace,xfer_prp=dxpl_id)
    if(iostat.ne.0) call mpi_hdf5_abort("load_5d_rk_parallel:h5dread_f",iostat)

    ! Release transfer and dataspace objects after the dataset transfer is complete.
    call h5pclose_f(dxpl_id,iostat)
    call h5sclose_f(memspace,iostat)
    call h5sclose_f(filespace,iostat)
    call h5dclose_f(dset_id,iostat)

    ! Close the collectively opened file before shutting down the HDF5 Fortran interface.
    call h5fclose_f(file_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_5d_rk_parallel:h5fclose_f",iostat)

    ! Release HDF5 library-level resources associated with this operation.
    call h5close_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_5d_rk_parallel:h5close_f",iostat)

  end subroutine load_5d_rk_parallel
!!!=================================================================================================
  !> Read a rank-local 4D REAL(kind=rk) hyperslab from a shared Parallel-HDF5 dataset.
  !!
  !! The first three dimensions follow the spatial decomposition and the complete size(data,4) range
  !! is read for the rank.
  !!
  !! @param data Array receiving the rank-local data read from the HDF5 dataset.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  subroutine load_4d_rk_parallel(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(out) :: data
    character(len=*), intent(in)                   :: type
    character(len=*), intent(in), optional         :: fname_optin
    integer, intent(in), optional                  :: nt_optin
    integer, intent(in), optional                  :: ns_optin
    character(len=*), intent(in), optional         :: dset_name_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t)                  :: file_id,dset_id,filespace,memspace,dxpl_id
    integer(hid_t)                  :: mem_type
    integer(hsize_t), dimension(4)  :: dims_local,dims_file,dims_max_file,count,stride
    integer(hssize_t), dimension(4) :: offset
    integer, dimension(4)           :: lbounds,ubounds,dims_local_i
    integer, dimension(3)           :: lbounds_space,ubounds_space
    integer                         :: iostat
    character(len=max_length_fname) :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Resolve the shared filename.  This also verifies that type is actually configured for MPI-I/O.
    call get_parallel_fname(fname,type,fname_optin,nt_optin,ns_optin)

    ! Use the caller-supplied dataset name when present; otherwise the conventional name is data.
    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! Obtain the cached global 1-based indices of the spatial cells owned by this rank.
    call get_parallel_spatial_bounds(lbounds_space,ubounds_space)

    ! Dimensions 1:3 are spatially decomposed.  Dimension 4 starts at one and is present in full
    ! in every rank-local buffer.
    lbounds = (/lbounds_space,1/)
    ubounds = (/ubounds_space,size(data,4)/)

    ! Record the actual local buffer shape and verify it against the selected global interval.
    dims_local_i = shape(data)

    ! This collective validation catches inconsistent geometry bounds before entering HDF5.
    call check_parallel_dimensions("load_4d_rk_parallel",dims_local_i,lbounds,ubounds)

    ! Convert Fortran default INTEGER extents to HDF5 size kinds and translate 1-based lower bounds
    ! to the zero-based offsets required by HDF5 hyperslabs.
    dims_local = int(dims_local_i,kind=hsize_t)
    offset     = int(lbounds-1,kind=hssize_t)
    count      = 1_hsize_t
    stride     = 1_hsize_t

    ! Initialize the HDF5 Fortran interface before creating/opening HDF5 objects.
    call h5open_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_rk_parallel:h5open_f",iostat)

    ! The destination array is REAL(kind=rk), therefore the HDF5 memory datatype
    ! always follows rk. The datatype stored in the file is obtained internally
    ! from the dataset. If the file contains float32 and rk is float64, HDF5
    ! automatically converts float32 -> float64 during h5dread_f.
    call get_hdf5_real_memory_type(mem_type)

    ! Open the shared file collectively in read-only mode with the MPIO file driver.
    call open_parallel_file_read(fname,file_id)

    ! Dataset open and dataspace query are performed by all ranks.
    ! Open the existing dataset on every rank; shared metadata operations remain synchronized.
    call h5dopen_f(file_id,trim(dset_name),dset_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_rk_parallel:h5dopen_f",iostat)

    ! Obtain the dataset's file dataspace so the rank can select its hyperslab.
    call h5dget_space_f(dset_id,filespace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_rk_parallel:h5dget_space_f",iostat)

    ! Read the actual on-disk extent instead of assuming it from the current parameter set.
    call h5sget_simple_extent_dims_f(filespace,dims_file,dims_max_file,iostat)
    if(iostat.lt.0) call mpi_hdf5_abort("load_4d_rk_parallel:h5sget_simple_extent_dims_f",iostat)

    ! Refuse to issue an HDF5 read if any rank's requested region lies outside the stored dataset.
    call check_parallel_hyperslab("load_4d_rk_parallel",offset,dims_local,dims_file)

    ! Select the contiguous rank-owned region in the global file dataspace.
    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offset,count,iostat,stride,dims_local)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_rk_parallel:h5sselect_hyperslab_f",iostat)

    ! Describe the contiguous local Fortran buffer independently of the global file dataspace.
    call h5screate_simple_f(4,dims_local,memspace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_rk_parallel:h5screate_simple_f(memspace)",iostat)

    ! Attach the requested collective MPI-I/O transfer mode to the following dataset operation.
    call create_collective_dxpl(dxpl_id)

    ! mem_type describes the REAL(kind=rk) destination buffer. The precision
    ! used when the dataset was written does not need to be known by this code.
    ! Read only this rank's selected hyperslab into the local memory buffer.
    call h5dread_f(dset_id,mem_type,data,dims_local,iostat, &
         mem_space_id=memspace,file_space_id=filespace,xfer_prp=dxpl_id)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_rk_parallel:h5dread_f",iostat)

    ! Release transfer and dataspace objects after the dataset transfer is complete.
    call h5pclose_f(dxpl_id,iostat)
    call h5sclose_f(memspace,iostat)
    call h5sclose_f(filespace,iostat)
    call h5dclose_f(dset_id,iostat)

    ! Close the collectively opened file before shutting down the HDF5 Fortran interface.
    call h5fclose_f(file_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_rk_parallel:h5fclose_f",iostat)

    ! Release HDF5 library-level resources associated with this operation.
    call h5close_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_rk_parallel:h5close_f",iostat)

  end subroutine load_4d_rk_parallel
!!!=================================================================================================
  !> Read a rank-local 4D INTEGER hyperslab from a shared Parallel-HDF5 dataset.
  !!
  !! The dataset uses native INTEGER transfers; only the first three dimensions are spatially
  !! decomposed.
  !!
  !! @param data Array receiving the rank-local data read from the HDF5 dataset.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  subroutine load_4d_i_parallel(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:,:,:,:), intent(out) :: data
    character(len=*), intent(in)             :: type
    character(len=*), intent(in), optional   :: fname_optin
    integer, intent(in), optional            :: nt_optin
    integer, intent(in), optional            :: ns_optin
    character(len=*), intent(in), optional   :: dset_name_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t)                  :: file_id,dset_id,filespace,memspace,dxpl_id
    integer(hsize_t), dimension(4)  :: dims_local,dims_file,dims_max_file,count,stride
    integer(hssize_t), dimension(4) :: offset
    integer, dimension(4)           :: lbounds,ubounds,dims_local_i
    integer, dimension(3)           :: lbounds_space,ubounds_space
    integer                         :: iostat
    character(len=max_length_fname) :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Resolve the shared filename.  This also verifies that type is actually configured for MPI-I/O.
    call get_parallel_fname(fname,type,fname_optin,nt_optin,ns_optin)

    ! Use the caller-supplied dataset name when present; otherwise the conventional name is data.
    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! Obtain the cached global 1-based indices of the spatial cells owned by this rank.
    call get_parallel_spatial_bounds(lbounds_space,ubounds_space)

    ! Dimensions 1:3 are spatially decomposed.  Dimension 4 starts at one and is present in full
    ! in every rank-local buffer.
    lbounds = (/lbounds_space,1/)
    ubounds = (/ubounds_space,size(data,4)/)

    ! Record the actual local buffer shape and verify it against the selected global interval.
    dims_local_i = shape(data)

    ! This collective validation catches inconsistent geometry bounds before entering HDF5.
    call check_parallel_dimensions("load_4d_i_parallel",dims_local_i,lbounds,ubounds)

    ! Convert Fortran default INTEGER extents to HDF5 size kinds and translate 1-based lower bounds
    ! to the zero-based offsets required by HDF5 hyperslabs.
    dims_local = int(dims_local_i,kind=hsize_t)
    offset     = int(lbounds-1,kind=hssize_t)
    count      = 1_hsize_t
    stride     = 1_hsize_t

    ! Initialize the HDF5 Fortran interface before creating/opening HDF5 objects.
    call h5open_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_i_parallel:h5open_f",iostat)

    ! Open the shared file collectively in read-only mode with the MPIO file driver.
    call open_parallel_file_read(fname,file_id)

    ! Dataset open and dataspace query are performed by all ranks.
    ! Open the existing dataset on every rank; shared metadata operations remain synchronized.
    call h5dopen_f(file_id,trim(dset_name),dset_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_i_parallel:h5dopen_f",iostat)

    ! Obtain the dataset's file dataspace so the rank can select its hyperslab.
    call h5dget_space_f(dset_id,filespace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_i_parallel:h5dget_space_f",iostat)

    ! Read the actual on-disk extent instead of assuming it from the current parameter set.
    call h5sget_simple_extent_dims_f(filespace,dims_file,dims_max_file,iostat)
    if(iostat.lt.0) call mpi_hdf5_abort("load_4d_i_parallel:h5sget_simple_extent_dims_f",iostat)

    ! Refuse to issue an HDF5 read if any rank's requested region lies outside the stored dataset.
    call check_parallel_hyperslab("load_4d_i_parallel",offset,dims_local,dims_file)

    ! Select the contiguous rank-owned region in the global file dataspace.
    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offset,count,iostat,stride,dims_local)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_i_parallel:h5sselect_hyperslab_f",iostat)

    ! Describe the contiguous local Fortran buffer independently of the global file dataspace.
    call h5screate_simple_f(4,dims_local,memspace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_i_parallel:h5screate_simple_f(memspace)",iostat)

    ! Attach the requested collective MPI-I/O transfer mode to the following dataset operation.
    call create_collective_dxpl(dxpl_id)

    ! Read only this rank's selected hyperslab into the local memory buffer.
    call h5dread_f(dset_id,h5t_native_integer,data,dims_local,iostat, &
         mem_space_id=memspace,file_space_id=filespace,xfer_prp=dxpl_id)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_i_parallel:h5dread_f",iostat)

    ! Release transfer and dataspace objects after the dataset transfer is complete.
    call h5pclose_f(dxpl_id,iostat)
    call h5sclose_f(memspace,iostat)
    call h5sclose_f(filespace,iostat)
    call h5dclose_f(dset_id,iostat)

    ! Close the collectively opened file before shutting down the HDF5 Fortran interface.
    call h5fclose_f(file_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_i_parallel:h5fclose_f",iostat)

    ! Release HDF5 library-level resources associated with this operation.
    call h5close_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_4d_i_parallel:h5close_f",iostat)

  end subroutine load_4d_i_parallel
!!!=================================================================================================
  !> Read a rank-local 3D REAL(kind=rk) hyperslab from a shared Parallel-HDF5 dataset.
  !!
  !! HDF5 converts the on-disk floating-point datatype to rk if the file was stored at a different
  !! precision.
  !!
  !! @param data Array receiving the rank-local data read from the HDF5 dataset.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  subroutine load_3d_rk_parallel(data, type, fname_optin, nt_optin, ns_optin, dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:), intent(out) :: data
    character(len=*), intent(in)                 :: type
    character(len=*), intent(in), optional       :: fname_optin
    integer, intent(in), optional                :: nt_optin
    integer, intent(in), optional                :: ns_optin
    character(len=*), intent(in), optional       :: dset_name_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t)                  :: file_id,dset_id,filespace,memspace,dxpl_id
    integer(hid_t)                  :: mem_type
    integer(hsize_t), dimension(3)  :: dims_local,dims_file,dims_max_file,count,stride
    integer(hssize_t), dimension(3) :: offset
    integer, dimension(3)           :: lbounds,ubounds,dims_local_i
    integer, dimension(3)           :: lbounds_space,ubounds_space
    integer                         :: iostat
    character(len=max_length_fname) :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Resolve the shared filename.  This also verifies that type is actually configured for MPI-I/O.
    call get_parallel_fname(fname,type,fname_optin,nt_optin,ns_optin)

    ! Use the caller-supplied dataset name when present; otherwise the conventional name is data.
    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! Obtain the cached global 1-based indices of the spatial cells owned by this rank.
    call get_parallel_spatial_bounds(lbounds_space,ubounds_space)

    lbounds = lbounds_space
    ubounds = ubounds_space

    ! Record the actual local buffer shape and verify it against the selected global interval.
    dims_local_i = shape(data)

    ! This collective validation catches inconsistent geometry bounds before entering HDF5.
    call check_parallel_dimensions("load_3d_rk_parallel",dims_local_i,lbounds,ubounds)

    ! Convert Fortran default INTEGER extents to HDF5 size kinds and translate 1-based lower bounds
    ! to the zero-based offsets required by HDF5 hyperslabs.
    dims_local = int(dims_local_i,kind=hsize_t)
    offset     = int(lbounds-1,kind=hssize_t)
    count      = 1_hsize_t
    stride     = 1_hsize_t

    ! Initialize the HDF5 Fortran interface before creating/opening HDF5 objects.
    call h5open_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_3d_rk_parallel:h5open_f",iostat)

    ! The destination array is REAL(kind=rk), therefore the HDF5 memory datatype
    ! always follows rk. The datatype stored in the file is obtained internally
    ! from the dataset. If the file contains float32 and rk is float64, HDF5
    ! automatically converts float32 -> float64 during h5dread_f.
    call get_hdf5_real_memory_type(mem_type)

    ! Open the shared file collectively in read-only mode with the MPIO file driver.
    call open_parallel_file_read(fname,file_id)

    ! Dataset open and dataspace query are performed by all ranks.
    ! Open the existing dataset on every rank; shared metadata operations remain synchronized.
    call h5dopen_f(file_id,trim(dset_name),dset_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_3d_rk_parallel:h5dopen_f",iostat)

    ! Obtain the dataset's file dataspace so the rank can select its hyperslab.
    call h5dget_space_f(dset_id,filespace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_3d_rk_parallel:h5dget_space_f",iostat)

    ! Read the actual on-disk extent instead of assuming it from the current parameter set.
    call h5sget_simple_extent_dims_f(filespace,dims_file,dims_max_file,iostat)
    if(iostat.lt.0) call mpi_hdf5_abort("load_3d_rk_parallel:h5sget_simple_extent_dims_f",iostat)

    ! Refuse to issue an HDF5 read if any rank's requested region lies outside the stored dataset.
    call check_parallel_hyperslab("load_3d_rk_parallel",offset,dims_local,dims_file)

    ! Select the contiguous rank-owned region in the global file dataspace.
    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offset,count,iostat,stride,dims_local)
    if(iostat.ne.0) call mpi_hdf5_abort("load_3d_rk_parallel:h5sselect_hyperslab_f",iostat)

    ! Describe the contiguous local Fortran buffer independently of the global file dataspace.
    call h5screate_simple_f(3,dims_local,memspace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_3d_rk_parallel:h5screate_simple_f(memspace)",iostat)

    ! Attach the requested collective MPI-I/O transfer mode to the following dataset operation.
    call create_collective_dxpl(dxpl_id)

    ! mem_type describes the REAL(kind=rk) destination buffer. The precision
    ! used when the dataset was written does not need to be known by this code.
    ! Read only this rank's selected hyperslab into the local memory buffer.
    call h5dread_f(dset_id,mem_type,data,dims_local,iostat, &
         mem_space_id=memspace,file_space_id=filespace,xfer_prp=dxpl_id)
    if(iostat.ne.0) call mpi_hdf5_abort("load_3d_rk_parallel:h5dread_f",iostat)

    ! Release transfer and dataspace objects after the dataset transfer is complete.
    call h5pclose_f(dxpl_id,iostat)
    call h5sclose_f(memspace,iostat)
    call h5sclose_f(filespace,iostat)
    call h5dclose_f(dset_id,iostat)

    ! Close the collectively opened file before shutting down the HDF5 Fortran interface.
    call h5fclose_f(file_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_3d_rk_parallel:h5fclose_f",iostat)

    ! Release HDF5 library-level resources associated with this operation.
    call h5close_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_3d_rk_parallel:h5close_f",iostat)

  end subroutine load_3d_rk_parallel
!!!=================================================================================================
  !> Load an explicitly specified 4D subregion.
  !!
  !! The caller-provided lbounds/ubounds are preserved exactly. Non-MPI-I/O types are forwarded to
  !! hdf5_single; MPI-I/O types read the requested hyperslab collectively from the shared file.
  !!
  !! @param data Array receiving the rank-local data read from the HDF5 dataset.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param lbounds Caller-supplied lower bounds of the requested 4D subregion.
  !! @param ubounds Caller-supplied upper bounds of the requested 4D subregion.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  subroutine load_sub_4d_rk(data, type, lbounds, ubounds, fname_optin, nt_optin, ns_optin, dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(:,:,:,:), intent(out) :: data
    character(len=*), intent(in)                   :: type
    integer, dimension(4), intent(in)              :: lbounds
    integer, dimension(4), intent(in)              :: ubounds
    character(len=*), intent(in), optional         :: fname_optin
    integer, intent(in), optional                  :: nt_optin
    integer, intent(in), optional                  :: ns_optin
    character(len=*), intent(in), optional         :: dset_name_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hid_t)                  :: file_id,dset_id,filespace,memspace,dxpl_id
    integer(hid_t)                  :: mem_type
    integer(hsize_t), dimension(4)  :: dims_file,dims_max_file,dims_local,count,stride
    integer(hssize_t), dimension(4) :: offset
    integer, dimension(4)           :: dims_local_i
    integer                         :: iostat
    logical                         :: mpi_io
    character(len=max_length_fname) :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Keep the explicit lbounds/ubounds interface of load_sub.
    ! Apply the same policy as ordinary loads; fname_optin preserves the legacy serial path.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)

    ! Serial/per-rank files retain the original hdf5_single explicit-bounds implementation.
    if(.not.mpi_io)then
       call load_sub_single(data,type,lbounds,ubounds,fname_optin,nt_optin,ns_optin,dset_name_optin)
       return
    end if

    ! Only the filename/file access is adapted to the shared MPI-HDF5 file.
    ! MPI-I/O case: resolve the shared filename while preserving the caller's explicit bounds.
    call get_parallel_fname(fname,type,fname_optin,nt_optin,ns_optin)

    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! The output array must normally match the inclusive lbounds/ubounds interval exactly.
    dims_local_i = shape(data)

    ! Preserve the historical warning behavior for a shape mismatch, but perform the check
    ! collectively so only one rank prints the warning.
    call check_parallel_dimensions("load_sub_4d_rk",dims_local_i,lbounds,ubounds,fatal=.false.)

    ! Explicit bounds are 1 based; HDF5 requires the corresponding zero-based file offset.
    dims_local = int(dims_local_i,kind=hsize_t)
    offset     = int(lbounds-1,kind=hssize_t)
    count      = 1_hsize_t
    stride     = 1_hsize_t

    ! Initialize HDF5 and obtain the rk-derived memory datatype used for conversion on read.
    call h5open_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_sub_4d_rk:h5open_f",iostat)

    ! The destination array is REAL(kind=rk), therefore the HDF5 memory datatype
    ! always follows rk. The datatype stored in the file is obtained internally
    ! from the dataset. If the file contains float32 and rk is float64, HDF5
    ! automatically converts float32 -> float64 during h5dread_f.
    call get_hdf5_real_memory_type(mem_type)

    ! Open the shared file collectively through the MPIO driver.
    call open_parallel_file_read(fname,file_id)

    ! Open the requested dataset and query its actual dataspace.
    call h5dopen_f(file_id,trim(dset_name),dset_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_sub_4d_rk:h5dopen_f",iostat)

    call h5dget_space_f(dset_id,filespace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_sub_4d_rk:h5dget_space_f",iostat)

    call h5sget_simple_extent_dims_f(filespace,dims_file,dims_max_file,iostat)
    if(iostat.lt.0) call mpi_hdf5_abort("load_sub_4d_rk:h5sget_simple_extent_dims_f",iostat)

    ! Validate the caller-provided subregion against the stored dataset before reading.
    call check_parallel_hyperslab("load_sub_4d_rk",offset,dims_local,dims_file)

    ! Select exactly the explicit [lbounds,ubounds] region in file coordinates.
    call h5sselect_hyperslab_f(filespace,h5s_select_set_f,offset,count,iostat,stride,dims_local)
    if(iostat.ne.0) call mpi_hdf5_abort("load_sub_4d_rk:h5sselect_hyperslab_f",iostat)

    call h5screate_simple_f(4,dims_local,memspace,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_sub_4d_rk:h5screate_simple_f(memspace)",iostat)

    ! Use the same collective dataset-transfer policy as the ordinary parallel load routines.
    call create_collective_dxpl(dxpl_id)

    ! mem_type describes the REAL(kind=rk) destination buffer. The precision
    ! used when the dataset was written does not need to be known by this code.
    ! HDF5 converts from the stored file datatype to REAL(kind=rk) when necessary.
    call h5dread_f(dset_id,mem_type,data,dims_local,iostat, &
         mem_space_id=memspace,file_space_id=filespace,xfer_prp=dxpl_id)
    if(iostat.ne.0) call mpi_hdf5_abort("load_sub_4d_rk:h5dread_f",iostat)

    call h5pclose_f(dxpl_id,iostat)
    call h5sclose_f(memspace,iostat)
    call h5sclose_f(filespace,iostat)
    call h5dclose_f(dset_id,iostat)

    call h5fclose_f(file_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_sub_4d_rk:h5fclose_f",iostat)

    call h5close_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("load_sub_4d_rk:h5close_f",iostat)

  end subroutine load_sub_4d_rk
!!!=================================================================================================
  !> Public dataset-dimension query dispatcher.
  !!
  !! The query uses the same type-based MPI-I/O policy as load/store so it inspects either the
  !! shared file or the legacy per-rank file consistently.
  !!
  !! @param dims Dataset dimensions returned as default Fortran INTEGER values.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  subroutine get_dset_dim(dims,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:), intent(out)     :: dims
    character(len=*), intent(in)           :: type
    character(len=*), intent(in), optional :: fname_optin
    integer, intent(in), optional          :: nt_optin,ns_optin
    character(len=*), intent(in), optional :: dset_name_optin
    logical                                :: mpi_io

    ! Without fname_optin, backend selection is based on the logical file type, not merely on array rank.
    call get_mpiio_policy(type,mpi_io,nt_optin,ns_optin,fname_optin)

    ! Shared-file types use Parallel HDF5; all other types retain the legacy hdf5_single path.
    if(mpi_io)then
       call get_dset_dim_parallel(dims,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
    else
       call get_dset_dim_single(dims,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
    end if

  end subroutine get_dset_dim
!!!=================================================================================================
  !> Read the dimensions of a dataset from a shared Parallel-HDF5 file.
  !!
  !! All ranks open the same file and dataset through the MPIO driver and return the dataset extent
  !! converted to default INTEGER.
  !!
  !! @param dims Dataset dimensions returned as default Fortran INTEGER values.
  !! @param type Logical file-type identifier used for filename generation and MPI-I/O policy
  !! selection.
  !! @param fname_optin Optional explicit filename. Public dispatchers preserve the legacy
  !! serial/per-rank path when present.
  !! @param nt_optin Optional timestep forwarded to the filename and I/O policy logic.
  !! @param ns_optin Optional subset index forwarded to the filename and I/O policy logic.
  !! @param dset_name_optin Optional dataset name; the backend uses its default dataset name when
  !! absent.
  subroutine get_dset_dim_parallel(dims,type,fname_optin,nt_optin,ns_optin,dset_name_optin)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer, dimension(:), intent(out)     :: dims
    character(len=*), intent(in)           :: type
    character(len=*), intent(in), optional :: fname_optin
    integer, intent(in), optional          :: nt_optin,ns_optin
    character(len=*), intent(in), optional :: dset_name_optin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(hsize_t), dimension(:), allocatable :: dims_file,dims_max_file
    integer(hid_t)                              :: file_id,file_space,dset_id
    integer                                     :: iostat,n1
    character(len=max_length_fname)             :: fname,dset_name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Initialize the output defensively before opening the file.
    dims = 0

    ! Resolve the same shared filename that the parallel load/store routines would use.
    call get_parallel_fname(fname,type,fname_optin,nt_optin,ns_optin)

    ! Query the requested dataset, defaulting to the conventional data dataset.
    if(present(dset_name_optin))then
       dset_name = trim(dset_name_optin)
    else
       dset_name = 'data'
    end if

    ! Allocate HDF5-size arrays with the rank expected by the caller-provided dims vector.
    n1 = size(dims)
    allocate(dims_file(n1),dims_max_file(n1))

    ! Open the shared file and dataset collectively through the MPIO driver.
    call h5open_f(iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("get_dset_dim_parallel:h5open_f",iostat)

    call open_parallel_file_read(fname,file_id)

    call h5dopen_f(file_id,trim(dset_name),dset_id,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("get_dset_dim_parallel:h5dopen_f",iostat)

    call h5dget_space_f(dset_id,file_space,iostat)
    if(iostat.ne.0) call mpi_hdf5_abort("get_dset_dim_parallel:h5dget_space_f",iostat)

    ! Obtain current and maximum extents; only the current extent is returned to the caller.
    call h5sget_simple_extent_dims_f(file_space,dims_file,dims_max_file,iostat)
    if(iostat.lt.0) call mpi_hdf5_abort("get_dset_dim_parallel:h5sget_simple_extent_dims_f",iostat)

    ! Convert HDF5 hsize_t extents back to the default INTEGER interface used by hdf5_single.
    dims = int(dims_file)

    ! Close metadata objects and deallocate the temporary HDF5 extent arrays.
    call h5sclose_f(file_space,iostat)
    call h5dclose_f(dset_id,iostat)
    call h5fclose_f(file_id,iostat)
    call h5close_f(iostat)

    deallocate(dims_file,dims_max_file)

  end subroutine get_dset_dim_parallel
!!!=================================================================================================

end module hdf5_mpi
