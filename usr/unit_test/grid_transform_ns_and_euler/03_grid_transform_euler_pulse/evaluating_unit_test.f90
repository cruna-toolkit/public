module evaluating_unit_test
 
  use helper
  use io
  use io_helper
  use parameter

  include 'fftw3.f'

  private

  public :: snapshot_comparison

contains



  subroutine snapshot_comparison()

    integer                                        :: dim_ref_data(4), counter, n, m, o,j
    real(kind=rk), dimension(:,:,:,:), allocatable :: reference_data, simulation_data
    real(kind=rk), dimension(:,:), allocatable     :: difference
    real(kind=rk)                                  :: max_difference, av_difference    
    integer, dimension(5)                          :: snapshots
    character(len=max_length_parameter)            :: filename
    character(len=max_length_parameter)            :: work_directory,reference_directory
    character(len=6)                               :: snapshot
    
        
    ! add reference data path if present
    call get_parameter(reference_directory,'io.referencedirectory',default = ".")
    reference_directory = trim(reference_directory) 
   
    ! add working path if present
    call get_parameter(work_directory,'io.workdirectory',default = ".")
    work_directory = trim(work_directory) 
        
    snapshots = [0, 10, 100, 200, 250]    
   
    filename = 'data_direct_snapshot__01_00001_01_000000.h5'
    call get_dset_dim(dim_ref_data, 'dummy', fname_optin= trim(reference_directory) // trim(filename),dset_name_optin='data')
    allocate(reference_data(dim_ref_data(1), dim_ref_data(2), dim_ref_data(3), dim_ref_data(4)))
    allocate(simulation_data(dim_ref_data(1), dim_ref_data(2), dim_ref_data(3), dim_ref_data(4)))
    allocate(difference(dim_ref_data(1)*dim_ref_data(2)*dim_ref_data(3),dim_ref_data(4)))
    
    do j = 1,size(snapshots)
      
      ! Define snapshot (to load)
      write(snapshot, '(I6.6)')  snapshots(j)    
      filename   = 'data_direct_snapshot__01_00001_01_'// snapshot //'.h5'
      write(*,*) 'iteration: ', j , ' filename: ', filename
      ! Load reference simulation data
      
      call load(reference_data,'dummy',fname_optin=  trim(reference_directory) // trim(filename) ,dset_name_optin='data')
      
      ! Load cruna simulation data
      
      call load(simulation_data,'dummy',fname_optin=  trim(work_directory) // trim(filename),dset_name_optin='data')
      
      counter = 1
      do n = 1, dim_ref_data(1)
         do m = 1, dim_ref_data(2)
           do o = 1, dim_ref_data(3)
                difference(counter,:) = abs(reference_data(n,m,o,:) - simulation_data(n,m,o,:))
                counter = counter + 1
           end do
         end do
      end do
      max_difference = maxval(difference)
      av_difference = sum(difference) / size(difference)

      write(*,*) '############ Unit Test Simulation test sinusiodial grid ############'
      write(*,*) '############  Compare ',j,', Snapshot: ',snapshots(j), '      ############'
      write(*,*) 'Maximum Difference between Reference and CRUNA simulation is: ', max_difference
      write(*,*) 'Average Difference between Reference and CRUNA simulation is: ', av_difference
      write(*,*) '############      Difference: rho, u,v,w, p    ############'
      write(*,*) 'Maximum Difference rho: ', maxval(difference(:,1)) 
      write(*,*) 'Maximum Difference   u: ', maxval(difference(:,2))
      write(*,*) 'Maximum Difference   v: ', maxval(difference(:,3))
      write(*,*) 'Maximum Difference   w: ', maxval(difference(:,4))
      write(*,*) 'Maximum Difference   p: ', maxval(difference(:,5))
      write(*,*) 'Average Difference rho: ', sum(difference(:,1)) / size(difference(:,1))
      write(*,*) 'Average Difference   u: ', sum(difference(:,2)) / size(difference(:,1))
      write(*,*) 'Average Difference   v: ', sum(difference(:,3)) / size(difference(:,1))
      write(*,*) 'Average Difference   w: ', sum(difference(:,4)) / size(difference(:,1))
      write(*,*) 'Average Difference   p: ', sum(difference(:,5)) / size(difference(:,1))
      if (av_difference < 1.0E-009) then
        write(*,*) 'Unit Test Compare snapshot: ',snapshots(j),' PASSED'
      else
        write(*,*) 'Unit Test Compare snapshot ',snapshots(j),' FAILED'
        stop
      end if    

    end do
  end subroutine snapshot_comparison

end module evaluating_unit_test
