

SUBROUTINE clean_artn()
  use units, only : DP
  use artn_params, only : lrelax, linit, lbasin, lperp, &
           llanczos, leigen, lsaddle, lbackward, lend,  &
           iartn, istep, iinit, iperp, ilanc, ieigen, nlanc, ifails,  &
           irelax, iover, istep, fpush_factor, lowest_eigval, nperp, nperp_list,  &
           artn_resume, old_lanczos_vec, H, Vmat, lanc_mat_size,  &
           iunartout, filout, nperp_step
  implicit none

  integer :: ios


  ! ...Fails if finished before it converged
  IF( .NOT.lend ) ifails = ifails + 1

  ! ...Write in output log
  WRITE(*,'(5x,"!> CLEANING ARTn | Fail:",x,i0)') ifails
  OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'old', POSITION = 'append', IOSTAT = ios )
    WRITE(iunartout,'(5x,"!> CLEANING ARTn | Fail:",x,i0/5x,*(a)//)') ifails, repeat("-",50)

  lrelax = .false.
  linit = .true.
  lbasin = .true.
  lperp = .false.
  llanczos = .false.
  leigen = .false.
  lsaddle = .false.
  lend = .false.

  !> Internal param
  lbackward = .true.
  fpush_factor = 1.0

  lend = .false.
  !
  iartn = 0
  istep = 0
  iinit = 0
  iperp = 0
  ilanc = 0
  ieigen = 0
  irelax = 0
  iover = 0

  ! ...Return the initial value of nperp
  nperp = nperp_list(1)
  nperp_step = 1
  !write(*,'(5x,"Reinitialize NPERP:",x,i0)') nperp

  
  lowest_eigval = 0.0_DP
  artn_resume = ""

  if( allocated(old_lanczos_vec) )old_lanczos_vec = 0.0_DP

  ! read the ARTn input file
  !OPEN( UNIT = iunartin, FILE = filin, FORM = 'formatted', STATUS = 'unknown', IOSTAT = ios)
  !READ( NML = artn_parameters, UNIT = iunartin)
  !CLOSE ( UNIT = iunartin, STATUS = 'KEEP')

  nlanc = lanc_mat_size

  H = 0.0_DP
  Vmat = 0.0_DP

  CLOSE ( UNIT = iunartout, STATUS = 'KEEP')

END SUBROUTINE clean_artn


