!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_ext_forces()
  !----------------------------------------------------------------------------
  !
  !
  USE mp,               ONLY : mp_bcast
  USE mp_images,        ONLY : intra_image_comm
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE kinds,            ONLY : DP
  !
  USE plugin_flags
  ! modifications start here
  ! we take some stuff from QE
  USE ions_base,     ONLY : nat, tau, if_pos,ntyp => nsp , ityp, atm, amass
  USE cell_base,     ONLY : alat, at
  USE force_mod,     ONLY : force
  USE ener,          ONLY : etot 
  USE relax,         ONLY : epsf, epse
  USE control_flags, ONLY : istep, conv_ions
  USE dynamics_module, ONLY : vel, dt, fire_alpha_init
  USE io_files,      ONLY : prefix,tmp_dir
  !
  IMPLICIT NONE
  ! 
  LOGICAL :: lconv
  REAL(DP), SAVE :: eps(2)

  ! ...ARTn Flag
  ! usage: ./pw.x -artn < input_qe
  IF( .not.use_partn )RETURN
  !
  ! ...ARTn convergence flag 
  lconv = .false. 
  !

  if( istep == 0 )then
    eps = [ epse, epsf ]
    !epsf = 1.0D-8   !! Modify it in artn_QE()
  endif

  IF ( ionode .and. use_partn ) THEN
     CALL artn_QE( force, etot, eps(2), nat,ntyp, ityp, atm, tau, at, alat, istep, if_pos, vel, dt, fire_alpha_init, &
          lconv, prefix, tmp_dir ) 
  ENDIF
  IF ( ionode .and. lconv .and. use_partn ) THEN
     WRITE (*,*) "ARTn calculation converged, stopping" 
     conv_ions = .true.
     epsf = eps(2)
     epse = 1. !eps(1)
     print*, " * PLUGIN_EXT_FORCE:: ", epsf, epse
     !CALL laxlib_end()
     !CALL stop_run( 0 )
     !CALL do_stop( 1 )
     !STOP 1
  END IF
  
END SUBROUTINE plugin_ext_forces
