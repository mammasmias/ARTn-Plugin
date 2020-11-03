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
  USE ions_base,     ONLY : nat, tau, if_pos, ityp, amass
  USE cell_base,     ONLY : alat, at
  USE force_mod,     ONLY : force
  USE control_flags, ONLY : istep
  USE dynamics_module, ONLY : vel, acc, dt, fire_alpha_init
  USE io_files,      ONLY : prefix,seqopn,tmp_dir
  !
  IMPLICIT NONE

  IF ( ionode ) THEN
     CALL artn(force,nat,alat,istep,if_pos,vel,acc,dt,fire_alpha_init,prefix,tmp_dir) 
  ENDIF
END SUBROUTINE plugin_ext_forces
