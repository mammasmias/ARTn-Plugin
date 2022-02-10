

SUBROUTINE clean_artn()
  use artn_params, only : lrelax, linit, lbasin, lperp, &
           llanczos, leigen, lsaddle, lbackward, lend, &
           iartn, istep, iinit, iperp, ilanc, ieigen,   &
           irelax, iover, fpush_factor
  implicit none

  write(*,*) " !> CLEANING PROCEDURE"

  lrelax = .false.
  linit = .true.
  lbasin = .true.
  lperp = .false.
  llanczos = .false.
  leigen = .false.
  lsaddle = .false.

  
  !lpush_final = .false.
  !lrestart = .false.
  !lmove_nextmin = .false.

  !> Internal param
  lbackward = .true.
  fpush_factor = 1.0

  lend = .false.
  !verbose = 0
  !
  iartn = 0
  !istep = 0
  iinit = 0
  iperp = 0
  ilanc = 0
  ieigen = 0
  !ismooth = 1
  !if_pos_ct = 0
  irelax = 0
  iover = 0



END SUBROUTINE clean_artn


