! 
! 
! Main ARTn plugin subroutine:
!        modifies the input force to perform the ARTn algorithm 
!------------------------------------------------------------------------------
SUBROUTINE artn_QE( force, etot, forc_conv_thr_qe, nat, ityp, atm, tau, at, alat, istep, if_pos, vel, dt_init, fire_alpha_init, lconv, prefix_qe, tmp_dir_qe )
  !----------------------------------------------------------------------------
  !
  USE artn_params, ONLY: DP, lartn 
  !
  !  Interface Quantum ESPRESSO/ARTn:
  !  We convert/compute/adapt some variable 
  !
  ! 
  IMPLICIT NONE
  REAL(DP), INTENT(INOUT) :: force(3,nat)     ! force calculated by the engine
  REAL(DP), INTENT(INOUT) :: vel(3,nat)       ! velocity of previous FIRE step
  REAL(DP), INTENT(INOUT) :: tau(3,nat)       ! atomic positions (needed for output only)
  REAL(DP), INTENT(INOUT) :: forc_conv_thr_qe ! force convergence threshold of the engine
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
  integer, save             :: step = 0
  character(len=4)          :: move
  real(dp)                  :: box(3,3)
  real(dp)                  :: pos(3,nat)
  real(dp)                  :: etot_fire, dt_curr, alpha
  INTEGER                   :: nsteppos

  LOGICAL                   :: file_exists
  CHARACTER(len=256)        :: filnam
  INTEGER                   :: ios
  !REAL(DP) :: dlanc      ! dR in Lanczos
  !REAL(DP) :: eigvec(3,nat)   !



  !%! Return if artn.in does not exist
  if( .not.lartn )return


  print*, " * IN ARTn_QE::"


  ! ...Convert the length in angstrom
  box = at * alat
  !allocate( pos(3,nat) )
  pos = tau * alat  

  print*, " * BOX::", box


  ! ...Launch ARTn
  !call artn( force, etot, forc_conv_thr_qe, nat, ityp, atm, tau, at, alat, istep, if_pos, vel, dt, fire_alpha_init, lconv,prefix, tmp_dir )
  !          in     in    in   in    in  inout in   in      out     out      out   out   out
  !call artn( force, etot, nat, ityp, atm, pos, box, if_pos, dlanc, eigenvec, iperp, move, lconv )
  call artn( force, etot, nat, ityp, atm, pos, box, if_pos, move, lconv )

  
  tau = pos / alat

  ! ...We give output dir once at the beginning 
! if( step == 0 )then

!   prefix = prefix_qe
!   tmp_dir = tmp_dir_qe

!   !alpha_init = fire_alpha_init
!   !dt_init = dt

!   step = step + 1
! endif



  ! ...Read the Fire parameters
  !filnam = trim(tmpdir) // '/' // trim(prfx) // '.' //'fire'
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



  ! ...Convert the dR given by ARTn to forces
  call move_mode( nat, force, vel, etot_fire, nsteppos, dt_curr, alpha, fire_alpha_init, dt_init, move, forc_conv_thr_qe )


  !
  ! write the FIRE parameters to its scratch file
  ! 
  OPEN( unit = 4, file = filnam, form = 'formatted', status = 'unknown', iostat = ios)
  WRITE( UNIT = 4, FMT = * ) etot_fire, nsteppos, dt_curr, alpha
  !
  CLOSE( UNIT = 4, STATUS = 'KEEP' )
  !
  !


  !deallocate( pos )

END SUBROUTINE artn_QE 



