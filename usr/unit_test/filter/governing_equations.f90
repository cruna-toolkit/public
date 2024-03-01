module governing_equations

  !use rhs_zero , rhs_direct  => right_hand_side_3d
  use euler_rho_rhou_p , rhs_direct => right_hand_side_direct_3d , rhs_Adjoint => right_hand_side_adjoint_3d

end module governing_equations
