! 
! 
! Main ARTn plugin subroutine:
!        modifies the input force to perform the ARTn algorithm 
!------------------------------------------------------------------------------
SUBROUTINE artn_QE( force, etot, forc_conv_thr_qe, nat, ityp, atm, tau, at, alat, istep, if_pos, vel, dt, fire_alpha_init, lconv, prefix, tmp_dir )
  !----------------------------------------------------------------------------
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
  REAL(DP), INTENT(IN) ::    dt               ! default time step in FIRE  
  REAL(DP), INTENT(IN) ::    fire_alpha_init  ! initial value of alpha in FIRE 
  REAL(DP), INTENT(IN) ::    alat             ! lattice parameter of QE
  REAL(DP), INTENT(IN) ::    at(3,3)          ! lattice parameters in alat units 
  INTEGER,  INTENT(IN) ::    nat              ! number of atoms
  INTEGER,  INTENT(IN) ::    ityp(nat)        ! atom types
  INTEGER,  INTENT(IN) ::    istep            ! current step
  INTEGER,  INTENT(IN) ::    if_pos(3,nat)    ! coordinates fixed by engine 
  CHARACTER(LEN=3),   INTENT(IN) :: atm(*)    ! name of atom corresponding to ityp
  CHARACTER(LEN=255), INTENT(IN) :: tmp_dir   ! scratch directory of engine 
  CHARACTER(LEN=255), INTENT(IN) :: prefix    ! prefix for scratch files of engine 
  LOGICAL, INTENT(OUT) :: lconv               ! flag for controlling convergence 
  !
  real(dp) :: box(3,3)
  real(dp), allocatable :: pos(:,:)

  ! ...Convert the length in angstrom
  box = at * alat
  allocate( pos(3,nat) )
  pos = tau * alat  

  ! ...

  ! ...Launch ARTn
  !call artn( force, etot, forc_conv_thr_qe, nat, ityp, atm, tau, at, alat, istep, if_pos, vel, dt, fire_alpha_init, lconv,prefix, tmp_dir )
  call artn( force, etot, nat, ityp, atm, pos, box, if_pos, lconv )

  ! ...Convert the dR given by ARTn to forces
  call move_atom( vel, dt, fire_alpha_init, forc_conv_thr_qe )



END SUBROUTINE artn_QE 
