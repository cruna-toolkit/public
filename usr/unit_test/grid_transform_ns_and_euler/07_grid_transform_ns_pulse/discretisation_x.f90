module discretisation_x
  ! controls parallelization

  use diff_i1_mpi_cartesian , only: Dx1 => diff , diff_i1_name => diff_i1_name , diff_i1_type => diff_i1_type
  use diff_i2_mpi_cartesian , only: Dx2 => diff , diff_i2_name => diff_i2_name , diff_i2_type => diff_i2_type
  use diff_i3_mpi_cartesian , only: Dx3 => diff , diff_i3_name => diff_i3_name , diff_i3_type => diff_i3_type

end module discretisation_x
