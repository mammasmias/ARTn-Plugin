
!> @author
!!  Matic Poberznik
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
  USE artn_params, ONLY:  lbasin, iperp, irelax, push, &
                          eigenvec, MOVE , &
                          prev_disp, iunartout, filout

  USE UNITS, Only: DP, convert_time, unconvert_time, &
                   unconvert_force, MASS

  !use debug, only: report_atom_prop
  !
  IMPLICIT NONE
  !
  ! -- Arguments
  INTEGER, value,             INTENT(IN)    :: nat
  INTEGER,                    INTENT(IN)    :: order(nat)
  REAL(DP), DIMENSION(3,nat), INTENT(IN)    :: displ_vec
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: force
  REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: vel
  REAL(DP),                   INTENT(IN)    :: alpha_init, dt_init
  REAL(DP),                   INTENT(INOUT) :: etot, alpha, dt_curr
  INTEGER,                    INTENT(INOUT) :: nsteppos
  INTEGER,                    INTENT(IN)    :: disp
  ! 
  ! -- Local Variables
  REAL(DP)                                  :: dt0, dt, tmp0, tmp1 !, dr(3,nat)
  REAL(DP), EXTERNAL                        :: ddot,dnrm2, dsum
  INTEGER                                   :: u0
  !character(256) :: ctmp
  !
  ! do things depending on mode of the move
  ! NOTE force units of Ry/a.u. are assumed ... 
  !
  ! .. Convert the force & time
  !force = convert_force( displ_vec )
  dt  = convert_time( dt_curr )
  dt0 = convert_time( dt_init )   !%! Finally we don't touch dt_init
  u0  = 73

  !if( istep == 1 )write(u0,'("# -- FIRE PARAMETERS"/" istep  |  alpha   |   dt   |  nsteppos")')
  !10 format(x,i0,2x,a,2(x,f15.5),x,i0)
  !open(newunit=u0,file="fire_params.dat", )
  !


  ! ...Save actuall displacement
  prev_disp = disp


  SELECT CASE( MOVE(disp) )
  !
  CASE( 'init' )
     !
     etot     = 1.D0
     vel(:,:) = 0.D0
     alpha    = 0.0_DP
     dt       = dt0
     nsteppos = 0
     !
     ! ...Displ_vec should be a Length
     force(:,:) = displ_vec(:,order(:))*Mass/dt**2
     !
  CASE( 'perp' )
     !
     ! ...Displ_vec is fperp
     force(:,:) = displ_vec(:,order(:))
     !
     IF( iperp - 1 .eq. 0 ) THEN  !%! Because I increment iperp before to enter in move_mode
        ! for the first step forget previous velocity (prevent P < 0)
        etot     = 0.D0
        vel(:,:) = 0.D0
        alpha    = alpha_init
        dt       = dt0
        nsteppos = 5

        !
     ELSE
        ! 
        ! subtract the components that are parallel
        IF( lbasin ) THEN
          tmp0     = ddot( 3*nat, vel(:,:), 1, push(:,order(:)), 1 )
          tmp1     = ddot( 3*nat, push(:,:), 1, push(:,:), 1 )          !> Don't need to be ordered
          vel(:,:) = vel(:,:) - tmp0 / tmp1 * push(:,order(:)) 
        ELSE
          tmp0     = ddot( 3*nat, vel(:,:)     , 1, eigenvec(:,order(:)), 1 )
          tmp1     = ddot( 3*nat, eigenvec(:,:), 1, eigenvec(:,:), 1 )  !> Don't need to be ordered
          vel(:,:) = vel(:,:) - tmp0 / tmp1 * eigenvec(:,order(:)) 
        ENDIF
        !  
     ENDIF
     !
     !write(u0,10) istep, MOVE(disp), alpha, dt, nsteppos

     ! ...Print info in file Atom_step_eigen.xyz
     !write(ctmp,'(i0)') istep
     !ctmp = trim(ctmp)//" - "//trim(MOVE(disp))
     !ctmp = trim(ctmp)//" Column: Pos - eigenvec - vel - force"
     !call report_atom_prop( "Atom_step_eigen.xyz", ctmp, nat, order,  eigenvec(:,order(:)), vel/Mass*dt, force/Mass*dt**2 )
     !call report_atom_prop( istep, disp, nat, "Atom_step_eigen.xyz", &
     !   " Column: Pos - eigenvec - vel - force", eigenvec, vel, force )
     !
     !
  ! CASE( 'lanc' )
  !    !
  !    ! ... set the velocity and acceleration and alpha of previous step to move correctly
  !    etot     = 0.D0
  !    vel(:,:) = 0.D0
  !    dt       = dt0
  !    alpha    = 0.D0
  !    nsteppos = 0
  !    !
  !    ! ... the step performed should be like this now translate it into the correct force
  !    ! force(:,:) = displ_vec(:,order(:))*lanczos_disp*Mass/dt**2
  !    force(:,:) = displ_vec(:,order(:))*Mass/dt**2
     !
  CASE( 'eign', 'over', 'smth', 'lanc' )
     !
     etot       = 0.D0
     vel(:,:)   = 0.D0
     alpha      = 0.0_DP
     dt         = dt0
     nsteppos   = 0
     force(:,:) = displ_vec(:,order(:))*Mass/dt**2

     ! ...Print info in file Atom_step_eigen.xyz
        !write(ctmp,'(i0)') istep
        !ctmp = trim(ctmp)//" - "//trim(MOVE(disp))
        !ctmp = trim(ctmp)//" Column: Pos - eigenvec - vel - force"
        !call report_atom_prop( "Atom_step_eigen.xyz", ctmp, nat, order, eigenvec(:,order(:)), vel/Mass*dt, force/Mass*dt**2 )
     !
  CASE( 'relx' )
     !forc_thr = 10D-8    !! QE dependent
     IF( irelax == 1 ) THEN
       alpha    = alpha_init
       dt       = dt0
     ENDIF
     !
     ! ... We reaload because it is unconverted at this place
     force(:,:) = displ_vec(:,order(:))
     !
  CASE default
     ! 
     write(*,'(5x,"|> No parameter conversion in move_mode:",x,a)') MOVE(disp)
     !  
  END SELECT

  !
  ! ...Unconvert the force & time
  dt_curr = unconvert_time( dt )
  force = unconvert_force( force )
  

END SUBROUTINE move_mode





