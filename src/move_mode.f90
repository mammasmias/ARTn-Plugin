
!> @author
!!  Matic Poberjnik,
!!  Miha Gunde
!!  Nicolas Salles


SUBROUTINE move_mode( nat, force, vel, etot, nsteppos, dt_curr, alpha, alpha_init, dt_init, disp, displ_vec )
  !
  !> @breif
  !!   translate specified move to appropriate force and set FIRE parameters accordingly  
  !
  !> @param [in]    nat		Size of list: Number of atoms
  !> @param [inout] force	List of force on atoms
  !> @param [inout] vel		List of atomic velicity 
  !> @param [in]    alpha_init	Initial Value of alpha parameter of FIRE algorithm
  !> @param [in]    dt_init     Initial Value of dt parameter of FIRE algorithm
  !> @param [inout] etot	Actual energy total of the system
  !> @param [inout] alpha	Value of alpha paramter of FIRE algorithm
  !> @param [inout] dt_curr	Value of dt paramter of FIRE algorithm
  !> @param [inout] nsteppos	??
  !> @param [in]    disp	Kind of actual displacement 
  !> @param [in]    displ_vec	Displacement field (unit lemgth/force/hessian ) 
  !
  USE artn_params, ONLY:  iperp, push0 => push, push=>eigenvec, dlanc, MOVE
  USE UNITS
  !
  IMPLICIT NONE

  ! -- Arguments
  INTEGER, INTENT(IN), value                       :: nat

  REAL(DP), DIMENSION(3,nat), INTENT(IN)    :: displ_vec
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
  integer :: i
  !
  ! do things depending on mode of the move
  ! NOTE force units of Ry/a.u. are assumed ... 
  !


  ! .. Convert the force & time
  force = convert_force( displ_vec )
  dt = convert_time( dt_curr )
  dt0 = convert_time( dt_init )   !%! Finally we don't touch dt_init


  SELECT CASE( MOVE(disp) )

  CASE( 'init' )
     !
     etot = 1.D0
     vel(:,:) = 0.D0
     alpha = 0.0_DP
     dt = dt0
     nsteppos = 0

     ! ...Displ_vec should be a Length
     !force(:,:) = push0(:,:)*amu_ry/dt_curr**2
     !force(:,:) = force(:,:)*amu_ry/dt**2
     force(:,:) = displ_vec(:,:)*amu_ry/dt**2

     !do i = 1,nat
     !print*, MOVE(disp), force(:,i)
     !enddo

  CASE( 'perp' )
     !

     ! ...Displ_vec is fperp
     force = displ_vec

     IF( iperp - 1 .eq. 0 ) THEN  !%! Because I increment iperp before to enter in move_mode
        ! for the first step forget previous velocity (prevent P < 0)
        etot = 0.D0
        vel(:,:) = 0.D0
        alpha = alpha_init
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
     !force(:,:) = force(:,:)*dlanc*amu_ry/dt**2
     force(:,:) = displ_vec(:,:)*dlanc*amu_ry/dt**2

     !do i = 1,nat
     !print*, MOVE(disp), push(:,i)
     !enddo

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
     !force(:,:) = force(:,:)*amu_ry/dt**2
     force(:,:) = displ_vec(:,:)*amu_ry/dt**2
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



