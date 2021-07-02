! 
! 
! Main ARTn plugin subroutine:
!        modifies the input force to perform the ARTn algorithm 
!------------------------------------------------------------------------------
!SUBROUTINE artn_QE( force, etot, forc_conv_thr_qe, nat, ityp, atm, tau, at, alat, istep, if_pos, vel, dt_init, fire_alpha_init, lconv, prefix_qe, tmp_dir_qe )
SUBROUTINE artn_QE( force, etot, epsf_qe, nat, ityp, atm, tau, at, alat, istep, if_pos, vel, dt_init, fire_alpha_init, lconv, prefix_qe, tmp_dir_qe )
  !----------------------------------------------------------------------------
  !
  use iso_c_binding, only : c_char, c_null_char
  USE artn_params, ONLY: DP, lartn, convcrit_final 
  !
  !  Interface Quantum ESPRESSO/ARTn:
  !  We convert/compute/adapt some variable 
  !
  ! 
  IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: force(3,nat)     ! force calculated by the engine
  REAL(DP), INTENT(INOUT) :: vel(3,nat)       ! velocity of previous FIRE step
  REAL(DP), INTENT(INOUT) :: tau(3,nat)       ! atomic positions (needed for output only)
  REAL(DP), INTENT(INOUT) :: epsf_qe          ! force convergence threshold of the engine
  REAL(DP), INTENT(IN) ::    etot             ! total energy in current step
  REAL(DP), INTENT(IN) ::    dt_init          ! default time step in FIRE  
  REAL(DP), INTENT(IN) ::    fire_alpha_init  ! initial value of alpha in FIRE 
  REAL(DP), INTENT(IN) ::    alat             ! lattice parameter of QE
  REAL(DP), INTENT(IN) ::    at(3,3)          ! lattice parameters in alat units 
  INTEGER,  INTENT(IN) ::    nat              ! number of atoms
  INTEGER,  INTENT(IN) ::    ityp(nat)        ! atom types
  INTEGER,  INTENT(IN) ::    istep            ! current step
  INTEGER,  INTENT(IN) ::    if_pos(3,nat)    ! coordinates fixed by engine 
  CHARACTER(LEN=3),   INTENT(IN) :: atm(*)    ! name of atom corresponding to ityp
  CHARACTER(LEN=255), INTENT(IN) :: tmp_dir_qe   ! scratch directory of engine 
  CHARACTER(LEN=255), INTENT(IN) :: prefix_qe    ! prefix for scratch files of engine 
  LOGICAL, INTENT(OUT) :: lconv               ! flag for controlling convergence 
  !  
  !integer, save             :: step = 0
  !character(len=4)          :: move
  character(len=1,kind=c_char), allocatable :: cmove(:)
  real(dp)                  :: box(3,3)
  real(dp)                  :: pos(3,nat)
  real(dp)                  :: etot_fire, dt_curr, alpha
  INTEGER                   :: nsteppos

  LOGICAL                   :: file_exists
  CHARACTER(len=256)        :: filnam
  INTEGER                   :: ios, i
  !REAL(DP) :: dlanc      ! dR in Lanczos
  !REAL(DP) :: eigvec(3,nat)   !


  !------------------------------------------------------------------------------------------------------------
  interface
    SUBROUTINE artn( force, etot, nat, ityp, atm, tau, at, if_pos, cmove, lconv )
      use iso_c_binding, only : c_char
      USE artn_params, ONLY: DP
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: force(3,nat)     ! force calculated by the engine
      REAL(DP), INTENT(INOUT) :: tau(3,nat)       ! atomic positions (needed for output only)

      REAL(DP), INTENT(IN) ::    etot             ! total energy in current step
      REAL(DP), INTENT(IN) ::    at(3,3)          ! lattice parameters in alat units 
      INTEGER,  INTENT(IN) ::    nat              ! number of atoms
      INTEGER,  INTENT(IN) ::    ityp(nat)        ! atom types
      INTEGER,  INTENT(IN) ::    if_pos(3,nat)    ! coordinates fixed by engine 
      CHARACTER(LEN=3),   INTENT(IN) :: atm(*)    ! name of atom corresponding to ityp

      CHARACTER(LEN=1,kind=c_char), allocatable, INTENT(OUT) :: cmove(:)       ! Stage for move_mode
      LOGICAL,          INTENT(OUT) :: lconv  
    END SUBROUTINE artn
    SUBROUTINE move_mode(nat, force, vel, etot, nsteppos, dt_curr, alpha, alpha_init, dt_init, cmode )
      use iso_c_binding, only : c_char
      USE artn_params, ONLY: DP, AMU_RY, iperp, push0 => push, push=>eigenvec, dlanc !,  &
      IMPLICIT NONE
      INTEGER, INTENT(IN)                       :: nat
      REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: force
      REAL(DP), DIMENSION(3,nat), INTENT(INOUT) :: vel
      REAL(DP), INTENT(IN)                      :: alpha_init, dt_init
      REAL(DP), INTENT(INOUT)                   :: etot, alpha, dt_curr
      INTEGER,  INTENT(INOUT)                   :: nsteppos
      CHARACTER(LEN=1,KIND=c_char), INTENT(IN)  :: cmode(:)
    END SUBROUTINE move_mode 
  end interface
  !------------------------------------------------------------------------------------------------------------


  !%! Return if artn.in does not exist
  if( .not.lartn )return


  print*, " * IN ARTn_QE::"
  print*, " * ARTn_QE::CONV:: ", epsf_qe


  ! ...Convert the length in angstrom
  box = at * alat
  pos = tau * alat  

  print*, " * BOX(c,l)::"
  do i = 1, 3
     print*, box(:,i)
  enddo




  ! ...Launch ARTn
  call artn( force, etot, nat, ityp, atm, pos, box, if_pos, cmove, lconv )


  ! ...Compare the Threshold
  if( epsf_qe < convcrit_final )then
    write( *,* ) "WARNING:: QE force threshold is lower than ARTn", epsf_qe, convcrit_final
    epsf_qe = convcrit_final  
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


  !print*, " * ARTn_QE::cmove ", cmove, size(cmove), (cmove(5)==c_null_char)

  ! ...Convert the dR given by ARTn to forces
  call move_mode( nat, force, vel, etot_fire, nsteppos, dt_curr, alpha, fire_alpha_init, dt_init, cmove )


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



