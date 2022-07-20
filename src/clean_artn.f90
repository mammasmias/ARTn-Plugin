

SUBROUTINE clean_artn()
  !> @brief
  !!   Clean and end the ARTn research to be ready for another or to stop
  !
  use units,       only : DP
  use artn_params, only : lrelax, linit, lbasin, lperp,                 &
           llanczos, leigen, lpush_over, lbackward, lend,               &
           iartn, istep, iinit, iperp, ilanc, ieigen, nlanc, ifails,    &
           irelax, iover, istep, fpush_factor, lowest_eigval,           &
           artn_resume, old_lanczos_vec, H, Vmat, lanczos_max_size,     &
           iunartout, filout, nperp_step, old_lowest_eigval, prev_disp, &
           error_message
  implicit none

  integer :: ios


  ! ...Fails if finished before it converged
  IF( .NOT.lend )then
    ifails = ifails + 1
    error_message = 'ARTn RESEARCH STOP BEFORE THE END'
    call write_fail_report( iunartout, prev_disp, lowest_eigval )
  ENDIF

  ! ...Write in output log
  WRITE(*,'(5x,"!> CLEANING ARTn | Fail:",x,i0)') ifails
  OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'old', POSITION = 'append', IOSTAT = ios )
    WRITE(iunartout,'(5x,"!> CLEANING ARTn | Fail:",x,i0/5x,*(a))') ifails, repeat("-",50)

  lrelax = .false.
  linit = .true.
  lbasin = .true.
  lperp = .false.
  llanczos = .false.
  leigen = .false.
  lpush_over = .false.
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
  call nperp_limitation_step( -1 )
  !nperp = nperp_list(1)
  !nperp_step = 1
  !write(*,'(5x,"Reinitialize NPERP:",x,i0)') nperp

  
  old_lowest_eigval = huge( old_lowest_eigval ) 
  lowest_eigval = 0.0_DP
  artn_resume = ""

  if( allocated(old_lanczos_vec) )old_lanczos_vec = 0.0_DP

  ! read the ARTn input file
  !OPEN( UNIT = iunartin, FILE = filin, FORM = 'formatted', STATUS = 'unknown', IOSTAT = ios)
  !READ( NML = artn_parameters, UNIT = iunartin)
  !CLOSE ( UNIT = iunartin, STATUS = 'KEEP')

  nlanc = lanczos_max_size

  H = 0.0_DP
  Vmat = 0.0_DP

  WRITE(iunartout,'(/)')
  CLOSE ( UNIT = iunartout, STATUS = 'KEEP')

END SUBROUTINE clean_artn
