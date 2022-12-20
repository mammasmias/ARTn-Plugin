
module debug
  !> @brief
  !!   This module contains all the routine used to debug the code

  !> @Note
  !!   - info_field
  !!   - report_atom_prop()
  !!   - compute_curve        ...Useless
  !!   - fire2_integration    ...Boh

  use units, only : DP
  implicit none

  ! 
  ! Curvature (DEBUG thing)
  REAL(DP), allocatable :: f0(:)
  REAL(DP) :: rcurv



 contains



  !............................................................
  subroutine info_field( u0, n, v, txt )
    !> @brief 
    !!   compute and print some information on the field (arrays rank 2)
    !!   Can be upgrade for efficiency
    !
    !> @param[in] u0    output unit channel
    !> @param[in] n     size of 2nd rank of array
    !> @param[in] v     array of rank 2
    !> @param[in] txt   comment on the array
    ! 
    use units, only : DP
    implicit none
 
    integer, intent(in) :: u0, n
    real(DP), intent(in) :: v(3,n)
    character(*), intent(in) :: txt
 
    integer :: i
    real(DP) :: lmax, vtot
    real(DP), external :: dsum
 
    write(u0,'("****** INFO FIELD:",x,a)') trim(txt)
 
    lmax = 0.0_DP
    do i = 1,n
       lmax = max( lmax, norm2(v(:,i)) )
    enddo
    call sum_force( v, n, vtot )
    vtot = dsum( 3*n, v )
 
    write(u0,'("* total field: ",g12.4," | max norm: ",g12.4/)') vtot, lmax

  end subroutine info_field




  !............................................................
  subroutine report_atom_prop( filename, comment, nat, order, prop, prop1, prop2 )
  !subroutine report_atom_prop( filename, comment, nat, prop, prop1, prop2 )
    !> @brief
    !!   write in file the atomic properties we want in more than the position
    !
    !> @param[in] filename    name of the output file
    !! @param[in] comment     text for the comment line in xyz file
    !! @param[in] nat     
    !
    use units,       only : DP, unconvert_length
    use artn_params, only : tau_step
    implicit none 
 
    character(*), intent(in) :: filename, comment
    integer, intent(in) :: nat, order(nat)
    real(DP), intent(in) :: prop(3,nat)
 
    real(DP), intent(in) :: prop1(3,nat), prop2(3,nat)
    !real(DP), intent(in), optional :: prop1(3,nat), prop2(3,nat)
 
    integer :: ii, i, u0
 
    logical, save :: first_time = .true.
    !logical :: opt1, opt2
  
    !opt1 = present(prop1)
    !opt2 = ( present(prop1).AND.present(prop2) )
  
    if( first_time )then !! Open a new file
      open( newunit=u0, file=trim(filename), status='replace', action='write')
      first_time = .false.
    else
      open( newunit=u0, file=trim(filename), status='old', action='write', access='append')
    endif
  
    write(u0,*) nat
    write(u0,'(5x,"STEP: ", a)') comment
  
    !! We print in Engine ORDER

    !if( opt2 )then
       do i = 1,nat
          write(u0,15) order(i), tau_step(:,order(i)), prop(:,i), prop1(:,i), prop2(:,i)
       enddo
       15 format(i0,x,*(f12.4))
    !else if( opt1 )then
    !  do i = 1,nat
    !     write(u0,*) i, tau_step(:,i), prop(:,i), prop1(:,i)
    !  enddo
    !else
    !   do i=1,nat
    !      write(u0,*) i, tau_step(:,i), prop(:,i)
    !   enddo
    !endif

    close( u0 )
  
  end subroutine report_atom_prop




  !............................................................
  subroutine compute_curve( iter, n, r, f )

    use units, only : DP, convert_length
    use artn_params, only : tau_step
    implicit none
 
    integer, intent(in) :: iter, n
    REAL(DP), intent(in) :: f(*), r(*)
 
    real(DP) :: df(n), dr(n),nf
 
    IF( iter == 0 )THEN
 
      !r0 = [ convert_length(r(1:n)) ]
      f0 = [ f(1:n) ]
      rcurv = 0.0
 
    ELSE
 
      rcurv = ( MAXVAL(ABS(f(1:n))) - MAXVAL(ABS(f0(1:n))) ) / MAXVAL(ABS(f0(1:n)))
      !dr = convert_length(r(1:n)) - r0(1:n)
      !curv = [ df(1:n) / dr(1:n) ]
      !rcurv =  MAXVAL(df/dnrm(f0)) 
 
      f0 = [ f(1:n) ]
      !r0 = [ convert_length(r(1:n)) ]
 
    ENDIF

  end subroutine compute_curve



  subroutine fire2_integration( istep, n, f, v, dt, alpha, delaystep, rmax )
    !> @brief 
    !!    Fire intgration following FIRE in lammps
 
    use units, only : DP
    use artn_params, only : u => iunartout, filout
    implicit none
 
    integer, intent( IN ) :: istep, n, delaystep
    REAL(DP), intent( in ) :: f(3,n), v(3,n), dt, alpha
    real(DP), intent( out ) :: rmax
 
    integer :: last_neg, step_start, nmov, i
 
    real(DP) :: vdotv, vdotf, fdotf, a1, a2, dtv, dtf, dtfm, ftm2v, mass, &
                dmax, vmax, dtgrow, dtshrink, dtmin, dtmax,  &
                alfa, alpha0, alpha_shrink
    real(DP) :: rsum
    real(DP) :: x(3,n), vloc(3,n), dr2(n)
    logical :: flagv0, halfback
 
    real(DP), external :: ddot, dsum
 
    !! -- Parameters
    halfback = .true.
    ftm2v = 9648.53_DP
    mass = 1.0_DP
    dmax = 0.5_DP
    dtgrow = 1.1
    dtshrink = 0.5
    dtmin = 6.0e-6
    dtmax = 0.006
    alpha0 = 0.0
    alpha_shrink = 0.9
    !! --------------------
 
 
    if( istep == 1 )step_start = istep
    vdotv = ddot( 3*n, v, 1, v, 1 )
    vdotf = ddot( 3*n, v, 1, f, 1 )
    fdotf = ddot( 3*n, f, 1, f, 1 )
 
    dtv = dt
 
    flagv0 = .false.
    if( vdotf > 0.0_DP )then
      a1 = 1. - alpha
      if( fdotf < 1.0e-8 )then; a2 = 0.0
      else; a2 = alpha * sqrt( vdotv / fdotf )
      endif
      !! delaystep 
      if( istep - last_neg > delaystep )then
        dtv = MIN(dt*dtgrow, dtmax)
        alfa = alpha * alpha_shrink
      endif
      call dcopy( 3*n, v, 1, vloc, 1 )
    else
      last_neg = istep
      !! delaystep
      if( .not.(istep - step_start < delaystep) )then
        alfa = alpha0
        if( dt*dtshrink >= dtmin )dtv = dt*dtshrink
      endif
      if( halfback ) x = x - 0.5 * dtv * v
      flagv0 = .true.
      vloc = 0.0
    endif
 
    !! now we work with vloc
 
    dtfm = dtv * ftm2v / mass
    if( flagv0 ) vloc = f * dtfm
 
    ! ...Rescale dtv
    !dtv = dt 
    vmax = MAXVAL( ABS(vloc) )
    if( dtv*vmax > dmax )dtv = dmax / vmax
 
    if( flagv0 )vloc = 0.0_DP
 
 
    ! ...Euler integration
 
    dtf = dtv * ftm2v
    dtfm = dtf / mass
    vloc = vloc + dtfm * f
    if( vdotf > 0.0_DP )vloc = a1 * vloc + a2 * f
    x = x + dtv * vloc
 
 
    ! ...Analysis
 
    !write(*,'("******* FIRE INTEGRATION ANALYSIS::",x,i0)') istep
    !write(*,'("* MAX displacment",x,f10.4)') MAXVAL(ABS(x))
    !write(*,'("* sum displcmeemt",x,f10.4)') sqrt( dsum(3*n,x) )
 
    rmax = 0.0
    dr2 = 0.0
    do i = 1,n
       dr2(i) = dot_product( x(:,i), x(:,i))
       rmax = MAX(rmax, dr2(i))
    enddo
    nmov = 0
    rsum = 0.0
    do i = 1,n
       if( dr2(i) > rmax*(1.-0.05)**2 )then
         nmov = nmov + 1
         rsum = rsum + sqrt(dr2(i))
       endif
    enddo
    rmax = sqrt( rmax )
 
    OPEN ( UNIT = u, FILE = filout, FORM = 'formatted', STATUS = 'old', POSITION = 'append' )
 
    write(u,'(5x,"|> ",i3," FIRE integration: rmax",x,g10.4,x,"| sum_d",x,g10.4,x,"| ",i0,x,g10.4)') &
       istep,  rmax, sqrt( dsum(3*n,x) ), nmov, rsum
 
    close( unit=u, STATUS='KEEP' )


  end subroutine fire2_integration
  

end module debug
