

SUBROUTINE clean_artn()
  use units, only : DP
  use artn_params, only : lrelax, linit, lbasin, lperp, &
           llanczos, leigen, lsaddle, lbackward, lend, &
           iartn, istep, iinit, iperp, ilanc, ieigen,   &
           irelax, iover, istep, fpush_factor, lowest_eigval,  &
           artn_resume, old_lanczos_vec
  implicit none

  write(*,'(5x,"!> CLEANING PROCEDURE")')

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
  istep = 0
  iinit = 0
  iperp = 0
  ilanc = 0
  ieigen = 0
  !ismooth = 1
  !if_pos_ct = 0
  irelax = 0
  iover = 0

  
  lowest_eigval = 0.0_DP
  artn_resume = ""

  if( allocated(old_lanczos_vec) )old_lanczos_vec = 0.0_DP


END SUBROUTINE clean_artn


