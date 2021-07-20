


SUBROUTINE move_mode( nat, force, vel, etot, nsteppos, dt_curr, alpha, alpha_init, dt_init, disp )
  !
  ! translate specified move to appropriate force and set FIRE parameters accordingly  
  !
  !use iso_c_binding, only : c_char
  USE artn_params, ONLY: DP, AMU_RY, iperp, push0 => push, push=>eigenvec, dlanc, MOVE
  USE UNITS
  !
  IMPLICIT NONE

  ! -- Arguments
  INTEGER, INTENT(IN), value                       :: nat

  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: force
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: vel

  REAL(DP), INTENT(IN)                      :: alpha_init, dt_init
  REAL(DP), INTENT(INOUT)                   :: etot, alpha, dt_curr
  INTEGER,  INTENT(INOUT)                   :: nsteppos
  
  INTEGER, INTENT(IN)                       :: disp
  ! 
  ! -- Local Variable
  REAL(DP) :: dt0, dt
  REAL(DP), EXTERNAL               :: ddot,dnrm2

  !
  ! do things depending on mode of the move
  ! NOTE force units of Ry/a.u. are assumed ... 
  !


  ! .. Convert the force & time
  force = convert_force( force )
  dt = convert_time( dt_curr )
  dt0 = convert_time( dt_init )   !%! Finally we don't touch di_init


  SELECT CASE( MOVE(disp) )

  CASE( 'init' )
     !
     etot = 1.D0
     vel(:,:) = 0.D0
     alpha = 0.0_DP
     !dt_curr = dt_init
     dt = dt0
     nsteppos = 0
     !force(:,:) = push0(:,:)*amu_ry/dt_curr**2
     force(:,:) = push0(:,:)*amu_ry/dt**2

  CASE( 'perp' )
     !
     !IF( iperp .eq. 0 ) THEN
     IF( iperp - 1 .eq. 0 ) THEN  !%! Because I increment iperp before to enter in move_mode
        ! for the first step forget previous velocity (prevent P < 0)
        etot = 0.D0
        vel(:,:) = 0.D0
        alpha = alpha_init
        !dt_curr = dt_init
        dt = dt0
        nsteppos = 5
     ELSE
        ! subtract the components that are parallel
        vel(:,:) = vel(:,:) - ddot( 3*nat, vel, 1, push, 1 )*push(:,:) / ddot( 3*nat, push(:,:), 1, push(:,:), 1 )
        !%! vel = 0.D0
        
     ENDIF
        !
  CASE( 'lanc' )
     !
     ! set the velocity and acceleration and alpha of previous step to move correctly
     !
     etot = 0.D0
     vel(:,:) = 0.D0
     !dt_curr = dt_init
     dt = dt0
     alpha = 0.D0
     nsteppos = 0
     ! the step performed should be like this now translate it into the correct force
     !force(:,:) = force(:,:)*dlanc*amu_ry/dt_curr**2
     force(:,:) = force(:,:)*dlanc*amu_ry/dt**2
     !force(:,:) = eigenvec(:,:)*dlanc*amu_ry/dt**2   ! Should be like that
     !
  CASE( 'eign' )
     !
     etot = 0.D0
     vel(:,:) = 0.D0
     alpha = 0.0_DP
     !dt_curr = dt_init
     dt = dt0
     nsteppos = 0
     !force(:,:) = force(:,:)*amu_ry/dt_curr**2
     force(:,:) = force(:,:)*amu_ry/dt**2
     !force(:,:) = eigenvec(:,:)*amu_ry/dt**2   ! Should be like that
     !

  CASE( 'relx' )
     !forc_thr = 10D-8    !! QE dependent
     alpha = alpha_init
     !dt_curr = dt_init
     dt = dt0

  CASE default
     write(*,*) 'Problem with move_mode!'

  END SELECT


  ! ...Unconvert the force & time
  dt_curr = unconvert_time( dt )
  force = unconvert_force( force )


END SUBROUTINE move_mode



