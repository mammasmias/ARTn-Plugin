
!> @author
!!  Matic Poberjnik,
!!  Miha Gunde
!!  Nicolas Salles

SUBROUTINE move_mode( nat, order, force, vel, etot, nsteppos, dt_curr, alpha, alpha_init, dt_init, disp, displ_vec )
  !
  !> @breif
  !!   translate specified move to appropriate force and set FIRE parameters accordingly  
  !
  !> @param [in]    nat		Size of list: Number of atoms
  !> @param [in]    order	Order of engine atoms list
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
  USE artn_params, ONLY:  lbasin, iperp, irelax, push, eigenvec, lanczos_disp, MOVE !, istep
  USE UNITS
  !
  IMPLICIT NONE

  ! -- Arguments
  INTEGER, INTENT(IN), value                :: nat
  INTEGER, INTENT(IN)                       :: order(nat)

  REAL(DP), DIMENSION(3,nat), INTENT(IN)    :: displ_vec
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: force
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: vel

  REAL(DP), INTENT(IN)                      :: alpha_init, dt_init
  REAL(DP), INTENT(INOUT)                   :: etot, alpha, dt_curr
  INTEGER,  INTENT(INOUT)                   :: nsteppos
  
  INTEGER, INTENT(IN)                       :: disp
  ! 
  ! -- Local Variable
  REAL(DP) :: dt0, dt, tmp0, tmp1, dr(3,nat)
  REAL(DP), EXTERNAL               :: ddot,dnrm2
  integer :: u0
  !
  ! do things depending on mode of the move
  ! NOTE force units of Ry/a.u. are assumed ... 
  !


  ! .. Convert the force & time
  !force = convert_force( displ_vec )
  dt = convert_time( dt_curr )
  dt0 = convert_time( dt_init )   !%! Finally we don't touch dt_init

  u0 = 73
  !if( istep == 1 )write(u0,'("# -- FIRE PARAMETERS"/" istep  |  alpha   |   dt   |  nsteppos")')
  !10 format(x,i0,2x,a,2(x,f15.5),x,i0)
  !open(newunit=u0,file="fire_params.dat", )


  SELECT CASE( MOVE(disp) )

  CASE( 'init' )
     !
     etot = 1.D0
     vel(:,:) = 0.D0
     alpha = 0.0_DP
     dt = dt0
     nsteppos = 0

     ! ...Displ_vec should be a Length
     force(:,:) = displ_vec(:,order(:))*amu_ry/dt**2

     !write(u0,10) istep, MOVE(disp), alpha, dt, nsteppos

  CASE( 'perp' )
     !

     ! ...Displ_vec is fperp
     force(:,:) = displ_vec(:,order(:))

     IF( iperp - 1 .eq. 0 ) THEN  !%! Because I increment iperp before to enter in move_mode
        ! for the first step forget previous velocity (prevent P < 0)
        etot = 0.D0
        vel(:,:) = 0.D0
        alpha = alpha_init
        dt = dt0
        nsteppos = 5

     ELSE
        ! subtract the components that are parallel
        if( lbasin )then
          tmp0 = ddot( 3*nat, vel(:,:), 1, push(:,order(:)), 1 )
          tmp1 = ddot( 3*nat, push(:,:), 1, push(:,:), 1 )  !> Don't need to be ordered
          vel(:,:) = vel(:,:) - tmp0 / tmp1 * push(:,order(:)) 
        else
          tmp0 = ddot( 3*nat, vel(:,:), 1, eigenvec(:,order(:)), 1 )
          tmp1 = ddot( 3*nat, eigenvec(:,:), 1, eigenvec(:,:), 1 )  !> Don't need to be ordered
          vel(:,:) = vel(:,:) - tmp0 / tmp1 * eigenvec(:,order(:)) 
        endif
        
     ENDIF
     !write(u0,10) istep, MOVE(disp), alpha, dt, nsteppos


     ! ...FIRE integration anticipation
     !call FIRE2_integration( iperp, nat, unconvert_force(force), vel, unconvert_time(dt), alpha, nsteppos, tmp0 )
     !if( tmp0 < 1.0e-2 )then
     !   force = force * (1. + 0.02)
     !   call FIRE2_integration( iperp, nat, unconvert_force(force), vel, unconvert_time(dt), alpha, nsteppos, tmp0 )
     !endif


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
     force(:,:) = displ_vec(:,order(:))*lanczos_disp*amu_ry/dt**2

     !write(u0,10) istep, MOVE(disp), alpha, dt, nsteppos

     !
  CASE( 'eign' )
     !
     etot = 0.D0
     vel(:,:) = 0.D0
     alpha = 0.0_DP
     !dt_curr = dt_init
     dt = dt0
     nsteppos = 0
     force(:,:) = displ_vec(:,order(:))*amu_ry/dt**2

     !write(u0,10) istep, MOVE(disp), alpha, dt, nsteppos
     !!

  CASE( 'relx' )
     !forc_thr = 10D-8    !! QE dependent
     if( irelax == 1 )then
       alpha = alpha_init
       !dt_curr = dt_init
       dt = dt0
     endif
     force(:,:) = displ_vec(:,order(:))

     !write(u0,10) istep, MOVE(disp), alpha, dt, nsteppos

  CASE default
     write(*,'(5x,"|> No parameter convertion in move_mode:",x,a)') MOVE(disp)

  END SELECT


  ! ...Unconvert the force & time
  dt_curr = unconvert_time( dt )
  force = unconvert_force( force )


END SUBROUTINE move_mode



