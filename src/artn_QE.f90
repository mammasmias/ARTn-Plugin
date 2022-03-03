! 
!> @author
!!   Matic Poberznik,
!!   Miha Gunde, 
!
!> @brief 
!!   Main ARTn plugin subroutine:
!> @details
!!   modifies the input force to perform the ARTn algorithm 
!------------------------------------------------------------------------------
SUBROUTINE artn_QE( force, etot, epsf_qe, nat, ityp, atm, tau, at, alat, istep, if_pos,   &
                    vel, dt_init, fire_alpha_init, lconv, prefix_qe, tmp_dir_qe )
  !----------------------------------------------------------------------------
  !
  USE units, ONLY : DP
  USE artn_params, ONLY: forc_thr 
  !
  !  Interface Quantum ESPRESSO/ARTn:
  !  We convert/compute/adapt some variables 
  !
  ! 
  IMPLICIT NONE
  REAL(DP),           INTENT(INOUT) :: force(3,nat)     !> force calculated by the engine
  REAL(DP),           INTENT(INOUT) :: vel(3,nat)       !> velocity of previous FIRE step
  REAL(DP),           INTENT(INOUT) :: tau(3,nat)       !> atomic positions (needed for output only)
  REAL(DP),           INTENT(INOUT) :: epsf_qe          !> force convergence threshold of the engine
  REAL(DP),           INTENT(IN) ::    etot             !> total energy in current step
  REAL(DP),           INTENT(IN) ::    dt_init          !> default time step in FIRE  
  REAL(DP),           INTENT(IN) ::    fire_alpha_init  !> initial value of alpha in FIRE 
  REAL(DP),           INTENT(IN) ::    alat             !> lattice parameter of QE
  REAL(DP),           INTENT(IN) ::    at(3,3)          !> lattice parameters in alat units 
  INTEGER,            INTENT(IN) ::    nat              !> number of atoms
  INTEGER,            INTENT(IN) ::    ityp(nat)        !> atom types
  INTEGER,            INTENT(IN) ::    istep            !> current step
  INTEGER,            INTENT(IN) ::    if_pos(3,nat)    !> coordinates fixed by engine 
  CHARACTER(LEN=3),   INTENT(IN) :: atm(*)              !> name of atom corresponding to ityp
  CHARACTER(LEN=255), INTENT(IN) :: tmp_dir_qe          !> scratch directory of engine 
  CHARACTER(LEN=255), INTENT(IN) :: prefix_qe           !> prefix for scratch files of engine 
  LOGICAL,            INTENT(OUT) :: lconv              !> flag for controlling convergence 
  !  
  REAL(DP)                  :: box(3,3)
  REAL(DP)                  :: pos(3,nat)
  REAL(DP)                  :: etot_fire, dt_curr, alpha
  REAL(DP)                  :: displ_vec(3,nat)
  INTEGER                   :: nsteppos, order(nat)

  LOGICAL                   :: file_exists
  CHARACTER(len=256)        :: filnam
  INTEGER                   :: ios, i, disp


  !------------------------------------------------------------------------------------------------------------
  interface
    SUBROUTINE artn( force, etot, nat, ityp, atm, tau, order, at, if_pos, disp, displ_vec, lconv )
      USE units, ONLY: DP
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: force(3,nat)     ! force calculated by the engine
      REAL(DP), INTENT(INOUT) :: tau(3,nat)       ! atomic positions (needed for output only)
      REAL(DP), INTENT(OUT) :: displ_vec(3,nat)   ! displacement vector communicated to move mode
      REAL(DP), INTENT(IN) ::    etot             ! total energy in current step
      REAL(DP), INTENT(IN) ::    at(3,3)          ! lattice parameters in alat units 
      INTEGER,  INTENT(IN), value ::    nat       ! number of atoms
      INTEGER,  INTENT(IN) ::    order(nat)       ! Engine order of atom
      INTEGER,  INTENT(IN) ::    ityp(nat)        ! atom types
      INTEGER,  INTENT(IN) ::    if_pos(3,nat)    ! coordinates fixed by engine 
      CHARACTER(LEN=3),   INTENT(IN) :: atm(*)    ! name of atom corresponding to ityp
      INTEGER,          INTENT(OUT) :: disp
      LOGICAL,          INTENT(OUT) :: lconv  
    END SUBROUTINE artn
    SUBROUTINE move_mode(nat, order, force, vel, etot, nsteppos, dt_curr, alpha, alpha_init, dt_init, disp, displ_vec )
      use units, only : DP
      USE artn_params, ONLY: iperp, push0 => push, push=>eigenvec, dlanc, move
      IMPLICIT NONE
      INTEGER, INTENT(IN), value                :: nat
      REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: force
      REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: vel
      REAL(DP), DIMENSION(3,nat), INTENT(IN)    :: displ_vec 
      REAL(DP), INTENT(IN)                      :: alpha_init, dt_init
      REAL(DP), INTENT(INOUT)                   :: etot, alpha, dt_curr
      INTEGER,  INTENT(INOUT)                   :: nsteppos
      INTEGER, INTENT(IN)                       :: disp, order(nat)
    END SUBROUTINE move_mode 
  end interface
  !------------------------------------------------------------------------------------------------------------



  print*, " * IN ARTn_QE::", nat
  !print*, " * ARTn_QE::CONV:: ", epsf_qe


  ! ...Convert the length in angstrom
  box = at * alat
  pos = tau * alat  

  !print*, " * BOX(c,l)::"
  !do i = 1, 3
  !   print*, box(:,i)
  !enddo

  do i = 1,nat
     order(i) = i
  enddo


  ! ...Launch ARTn
  call artn( force, etot, nat, ityp, atm, pos, order, box, if_pos, disp, displ_vec, lconv )


  ! ...Compare the Threshold
  if( epsf_qe < forc_thr )then
    write( *,* ) "WARNING:: QE force threshold is lower than ARTn", epsf_qe, forc_thr
    epsf_qe = forc_thr  
  endif
  
  ! ...Change the position to QE
  tau = pos / alat



  ! ...Read the Fire parameters
  filnam = trim(tmp_dir_qe) // '/' // trim(prefix_qe) // '.' //'fire'
  INQUIRE( file = filnam, exist = file_exists )
  OPEN( unit = 4, file = filnam, form = 'formatted', status = 'unknown', iostat = ios)
  !
  IF (file_exists ) THEN
     ! if file exists read the data, otherwise just close it 
     READ( UNIT = 4, FMT = * ) etot_fire, nsteppos, dt_curr, alpha
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
  ELSE
     CLOSE( UNIT = 4, STATUS = 'DELETE')
  ENDIF


  !print*, " * ARTn_QE::PRE_MOVEMODE "

  ! ...Convert the dR given by ARTn to forces
  !call move_mode( nat, force, vel, etot_fire, nsteppos, dt_curr, alpha, fire_alpha_init, dt_init, disp, displ_vec )
  call move_mode( nat, order, force, vel, etot_fire, nsteppos, dt_curr, alpha, fire_alpha_init, dt_init, disp, displ_vec )



  ! ...Clean ARTn 
  IF( lconv )call clean_artn()


  !
  ! write the FIRE parameters to its scratch file
  ! 
  OPEN( unit = 4, file = filnam, form = 'formatted', status = 'unknown', iostat = ios)
  WRITE( UNIT = 4, FMT = * ) etot_fire, nsteppos, dt_curr, alpha
  !
  CLOSE( UNIT = 4, STATUS = 'KEEP' )
  !
  !


END SUBROUTINE artn_QE 



