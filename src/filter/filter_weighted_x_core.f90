module filter_weighted_x_core

  use filter_weighted_pen_core_empty    , only: filter_i1 , filter_i1_type => filter_type, filter_i1_name => filter_name
  use filter_weighted_pen_core_empty    , only: filter_i2 , filter_i2_type => filter_type, filter_i2_name => filter_name
  use filter_weighted_pen_core_empty    , only: filter_i3 , filter_i3_type => filter_type, filter_i3_name => filter_name

  ! for init (init must be compatible with ALL pen_cores above)
  use filter_weighted_pen_core_empty    , only: init_filter => init_filter

end module filter_weighted_x_core
