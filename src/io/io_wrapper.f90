module io_wrapper

  use helper
  use io
  use io_helper
  use objective
  use parallelism
  use parameter
  use reporting

  private

!  public load_state
  public store_state
  public store_J

contains

!!!=================================================================================================
  subroutine store_state(Q,type)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:,:), intent(in)                :: Q
    character(len=*)                  , intent(in)                :: type
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk),dimension(:,:,:,:,:), allocatable               :: Qt
    real                                                          :: toc_store
    integer                                                       :: s,s_backup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (params%time%subsets.gt.1) then
       !write(*,*) "error in io_wrapper.f90:store_state:storing of multiple subsets not yet implemented"
       !call stop_cruna

       ! store actual subset first, as no load is required
       call tic()
       call store(Q,trim(type))

       ! store other subsets by load first then store (a.k.a. rename)
       allocate(Qt(params%geom%n1b,params%geom%n2b,params%geom%n3b,params%equation%nbr_vars,params%time%steps))

       s_backup = params%time%ns

       do s = 1,params%time%subsets
          params%time%ns = s

          if (trim(type).eq.'data_direct_loop') then
             call  load(Qt,'data_direct_subsets_cache')
             call store(Qt,'data_direct_loop')
          end if

          if (trim(type).eq.'data_adjoint_loop') then
             call  load(Qt,'data_adjoint_subsets_cache')
             call store(Qt,'data_adjoint_loop')
          end if
          
       end do

       toc_store = toc()
       params%time%ns = s_backup

    else
       ! store state direct or adjoint
       call tic()
       call store(Q,trim(type))
       toc_store = toc()

    end if

    if (params%parallelism%world_image.eq.1) then
       write(*,*) "store state (",trim(type),") T/io: ",trim(num2str(toc_store,'(F0.2)')), &
            " (",trim(num2str(real(size(Q))*real(params%time%subsets)*real(params%parallelism%world_size)*real(rk)/real(1024)/real(1024)/toc_store,'(F0.2)'))," MB/s)"
    end if

  end subroutine store_state
!!!=================================================================================================

! !!!=================================================================================================
!   subroutine load_state(Q,type)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     real(kind=rk),dimension(:,:,:,:,:), intent(out)               :: Q
!     character(len=*)                  , intent(in)                :: type
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     real                                                          :: toc_store
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     if (params%time%subsets.gt.1) then
!        write(*,*) "error in io_wrapper.f90:load_subsets:loading of mutliple subsets not yet implemented"
!        call stop_cruna
!     else
!        ! load state direct or adjoint
!        call tic()
!        call load(Q,trim(type))
!        toc_store = toc()

!        if (params%parallelism%world_image.eq.1) then
!           write(*,*) "load state (",trim(type),") T/io: ",trim(num2str(toc_store,'(F0.2)')), &
!                " (",trim(num2str(real(size(Q))*real(params%parallelism%world_size)*real(rk)/real(1024)/real(1024)/toc_store,'(F0.2)'))," MB/s)"
!        end if

!     end if

!   end subroutine load_state
! !!!=================================================================================================

!!!=================================================================================================
  subroutine store_J
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer                                                        :: loop
    integer                                                        :: io_unit,io_stat
    character(len=max_length_fname)                                :: fname
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! store progress of J in hdf5
    call store(objective_function,'objective')

!!! store progress of J in ASCII (no one can read h5 file on a console using cat)
    call get_parameter(loop,'opt.loop',default = 1)
    call get_fname(fname,'objective_ascii')

    if (loop.eq.1) then
       ! open output file
       open(newunit=io_unit, file=trim(fname), iostat=io_stat, status='unknown') 

       ! write reporting stuff
       write(io_unit,*) "progress of objective for case: ",cruna_report_case
       close(io_unit)
    end if

    ! open output file
    open(newunit=io_unit, file=trim(fname), iostat=io_stat, status='old',position='append')

    ! write to file
    write(io_unit,*) loop,objective_function(loop,:)

    ! close file
    close(io_unit)

  end subroutine store_J
!!!=================================================================================================

end module io_wrapper
