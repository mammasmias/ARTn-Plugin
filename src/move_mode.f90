


!SUBROUTINE move_mode(nat, force, vel, etot, nsteppos, dt_curr, alpha, alpha_init, dt_init, cmode )
SUBROUTINE move_mode( nat, force, vel, etot, nsteppos, dt_curr, alpha, alpha_init, dt_init, disp )
  !
  ! translate specified move to appropriate force and set FIRE parameters accordingly  
  !
  use iso_c_binding, only : c_char
  USE artn_params, ONLY: DP, AMU_RY, iperp, push0 => push, push=>eigenvec, dlanc, MOVE
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
  !CHARACTER(LEN=1,KIND=c_har), INTENT(IN)  :: cmode(:)
  ! 
  ! -- Local Variable
  !CHARACTER(LEN=:), allocatable             :: mode
  !CHARACTER(LEN=:), allocatable, external    :: ctrim
  REAL(DP), EXTERNAL               :: ddot,dnrm2

  real(DP) :: vdtf
  integer :: i, atmov;
  !
  ! do things depending on mode of the move
  ! NOTE force units of Ry/a.u. are assumed ... 
  !

  ! .. Convert C_Char to Fortran String
  !print*, " * ARTn::MOVE_MODE::", cmode
  !mode = ctrim( cmode )


  print*, " * ARTn::MOVE_MODE::mode->", MOVE(disp)
  print*, " * ARTn::MOVE_MODE::iperp->", iperp
  print*, " * ARTn::MOVE_MODE::arg: ", nat, etot, amu_ry, dt_curr, dt_init

!  do i = 1,10
!     print*, " * force", i,force(:,i)
!  enddo
     atmov = 242
     print*, " * force", atmov,force(:,atmov)


  !SELECT CASE( TRIM(mode) )
  SELECT CASE( MOVE(disp) )

  CASE( 'init' )
     !
     etot = 1.D0
     vel(:,:) = 0.D0
     alpha = 0.0_DP
     dt_curr = dt_init
     nsteppos = 0
     force(:,:) = push0(:,:)*amu_ry/dt_curr**2

  CASE( 'perp' )
     !
     !IF( iperp .eq. 0 ) THEN
     IF( iperp - 1 .eq. 0 ) THEN  !%! Because I increment iperp before to enter in move_mode
        ! for the first step forget previous velocity (prevent P < 0)
        etot = 0.D0
        vel(:,:) = 0.D0
        alpha = alpha_init
        dt_curr = dt_init
     ELSE
        ! subtract the components that are parallel
        vel(:,:) = vel(:,:) - ddot( 3*nat, vel, 1, push, 1 )*push(:,:) / ddot( 3*nat, push(:,:), 1, push(:,:), 1 )
     ENDIF
        !
  CASE( 'lanc' )
     !
     ! set the velocity and acceleration and alpha of previous step to move correctly
     !
     etot = 0.D0
     vel(:,:) = 0.D0
     dt_curr = dt_init
     alpha = 0.D0
     nsteppos = 0
     ! the step performed should be like this now translate it into the correct force
     force(:,:) = force(:,:)*dlanc*amu_ry/dt_curr**2
     !force(:,:) = eigenvec(:,:)*dlanc*amu_ry/dt_curr**2   ! Should be like that
     !
  CASE( 'eign' )
     !
     etot = 0.D0
     vel(:,:) = 0.D0
     alpha = 0.0_DP
     dt_curr = dt_init
     nsteppos = 0
     force(:,:) = force(:,:)*amu_ry/dt_curr**2
     !force(:,:) = eigenvec(:,:)*amu_ry/dt_curr**2   ! Should be like that
     !

  CASE( 'relx' )
     !forc_thr = 10D-8    !! QE dependent
     alpha = alpha_init
     dt_curr = dt_init

  CASE default
     write(*,*) 'Problem with move_mode!'

  END SELECT


  print*, " MOVE_MODE::END ", alpha, dt_curr, nsteppos, dt_curr, amu_ry, amu_ry/dt_curr**2
  print*, " * force", atmov, force(:,atmov)
  print*, " * push", atmov, push0(:,atmov)
  print*, " * dx", atmov, force(:,atmov)*dt_curr**2/amu_ry
  print*, " * dx", atmov, force(:,atmov)*dt_curr**2/amu_ry * 0.529

  ! ... Print some value
  !print*, " * MOVE_MODE::END::Force", MAXVAL(force)
  !print*, " * MOVE_MODE::END::Veloc", MAXVAL(vel)
  !vdtf = 0.0
  !do i = 1,nat
  !   vdtf = vdtf + dot_product( vel(:,i), force(:,i) )
  !enddo
  !print*, " * MOVE_MODE::END::VdotF", vdtf




END SUBROUTINE move_mode



