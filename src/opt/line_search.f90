module line_search

  use computation_direct
  use data_rhs
  use force
  use helper
  use parameter

  ! case dependent
  use force
  use geometry
  use initial_condition
  use objective
  use parallelism

  private
  
  public :: ls_quadratic_fit

contains

!!!=================================================================================================
  subroutine ls_quadratic_fit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(kind=rk), dimension(3,3)                                       :: M
    real(kind=rk), dimension(3)                                         :: Js,alphas,polynom_coefficients
    real(kind=rk), dimension(3)                                         :: ipiv
    real(kind=rk)                                                       :: alpha,alpha0,J
    integer                                                             :: loop,ls
    integer                                                             :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! backup old optima (must be code in force, ...)
    if (params%io%verbosity.ge.1) then
       if (params%parallelism%world_image.eq.1) then
          write(*,*) "line-search: backup old optima"
       end if
    end if

    call backup_force

!!! init line-search
    call set_parameter(1     ,'opt.in_line_search'             )
    call get_parameter(loop  ,'opt.loop'          , default = 1)
    call get_parameter(alpha0,'opt.alpha0'                     )

    alphas = (/ 0.0_rk, 0.5_rk, 1.0_rk /)
    alphas = alphas * alpha0

    do ls = 1,3
       M(ls,:) = (/ alphas(ls)**2, alphas(ls), 1.0_rk /)
    end do

    Js(1) = objective_function(loop,1)

!!! the line-search
    do ls = 1,2

       ! add half alpha0 each search
       call set_parameter(alphas(ls+1) ,'opt.alpha')

       ! update the force
       call update_force

       ! direct computation
       if(params%parallelism%world_image.eq.1) then
          write(*,*) 'line-search: starting direct compuation (', trim(num2str(ls,'(I1.1)')),'/2)'
       end if

       call init_direct(q0)
       call calc_direct(Q,q0)

       ! calc objective
       call calc_objective(J,Q)
       Js(ls+1) = J

    end do

!!! calc opt alpha
    polynom_coefficients = Js
    call dgesv( 3, 1, M, 3, ipiv, polynom_coefficients, 3, info )

    if(info.gt.0) then

       if (params%parallelism%world_image.eq.1) then
          write(*,*) "line-search: *** WARNING *** :coefficient matrix is (close to) singular:use alpha0 as new alpha"
       end if
       call set_parameter(alpha0,'opt.alpha')

    else

       alpha = -polynom_coefficients(2)/(2*polynom_coefficients(1))

       if (params%parallelism%world_image.eq.1) then
          write(*,*) "line-search: alpha0, new alpha: ", alpha0,alpha
          write(*,*) "line-search: objectives 0,1,2"                , Js
          write(*,*) "line-search: norm. objectives 0,1,2"          , Js/objective_function(1,1)
          write(*,*) "line-search: polynom_coefficients ax^2+bx+c: ", polynom_coefficients
       end if

       ! check for large alpha
       if(alpha.gt.(10.0_rk * alpha0)) then

          alpha = 10.0_rk * alpha0
          
          if (params%parallelism%world_image.eq.1) then
             write(*,*) "line-search: *** WARNING *** :alpha too large:alpha is reduced to 10x alpha0", alpha, alpha0
          end if

       end if

    end if

    call set_parameter(alpha,'opt.alpha')

!!! restore force
    if (params%io%verbosity.ge.1) then
       if (params%parallelism%world_image.eq.1) then
          write(*,*) "line-search: restore old optima"
       end if
    end if

    call restore_force

!!! exit line-search
    call set_parameter(0,'opt.in_line_search')

  end subroutine ls_quadratic_fit
!!!=================================================================================================

end module line_search
