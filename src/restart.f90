
!> @authors
!!   Matic Poberznik
!!   Miha gunde
!!   Nicolas Salles

!
!---------------------------------------------------------------------------
SUBROUTINE write_restart( filnres )
  !
  !> @breif 
  !!   Subroutine that writes the minimum parameters required for restart of a calculation
  !!   to a file
  !
  !> @param[in]    filnres   restart filename
  ! 
  !> @note 
  !!   LOGICAL FLAGS: linit, lperp, leigen, llanczos, lsaddle, lrelax
  !!   COUNTERS : istep, iinit, ilanc, ieigen, ismooth, nlanc
  !!   ARRAYS:
  !!   - LANCZOS: eigenvec, H, Vmat, force_old, lowest_eigval
  !!   - (INIT/SADDLE/STEP): etot, tau, force (, current_step_size, fpush_factor)
  !
  use artn_params, only : linit, lperp, leigen, llanczos, lpush_over, lrelax, &
                          iartn, istep, iinit, ieigen, iperp, ilanc, irelax, ismooth,   &
                          ninit, neigen, nlanc, lanczos_max_size, nperp, nmin, nsaddle, &
                          etot_init, &
                          etot_step, tau_step, force_step, current_step_size, fpush_factor, &    !> Actual step
                          eigenvec, H, Vmat, force_old, lowest_eigval, &
                          etot_saddle, tau_saddle, iunartres

  implicit none

  CHARACTER(LEN=255), INTENT(IN) :: filnres
  INTEGER :: ios

  !print*, "WRITE_RESTART::", istep
  OPEN( UNIT = iunartres, FILE = filnres, ACTION="WRITE", FORM = 'formatted', STATUS = 'REPLACE', IOSTAT = ios)

  WRITE ( iunartres, * ) linit, lperp, leigen, llanczos, lpush_over, lrelax, &
       iartn, istep, iinit, ieigen, iperp, ilanc, irelax, ismooth,   &
       ninit, neigen, nlanc, lanczos_max_size, nperp, nmin, nsaddle, &
       etot_init, &
       etot_step, tau_step, force_step, current_step_size, fpush_factor, &    !> Actual step
       eigenvec, H, Vmat, force_old, lowest_eigval, &
       etot_saddle, tau_saddle

  CLOSE ( UNIT = iunartres, STATUS = 'KEEP')

END SUBROUTINE write_restart



!
!---------------------------------------------------------------------------
SUBROUTINE read_restart( filnres, nat, order, ityp, ierr )
  !
  !> @brief
  !!   Subroutine that reads the restart file, if a restart is requested
  !
  !> @param[in]   filnres     restart filename
  !> @param[in]   nat         number of atoms
  !> @param[in]   order       index order of atoms
  !> @param[in]   itype       type of atoms
  !> @param[out]  ierr        flag for the error
  !
  !> @note
  !!   LOGICAL FLAGS: lpush_init, lperp, leigen, llanczos, lsaddle, lrelax
  !!   COUNTERS : istep, ipush, ilanc, ieigen, ismooth, nlanc
  !!   ARRAYS: eigenvec, H, Vmat
  !
  !> @warning
  !!   - Maybe doing a verification on the nat parameter in case...
  !!   - order argument is not used in the routine
  !!   - ityp is not used in the routine
  !
  use artn_params, only : linit, lperp, leigen, llanczos, lpush_over, lrelax, &
                          iartn, istep, iinit, ieigen, iperp, ilanc, irelax, ismooth,   &
                          ninit, neigen, nlanc, lanczos_max_size, nperp, nmin, nsaddle, &
                          etot_init, &
                          etot_step, tau_step, force_step, current_step_size, fpush_factor, &    !> Actual step
                          eigenvec, H, Vmat, force_old, lowest_eigval, &
                          etot_saddle, tau_saddle, &
                          tau_init, initpfname, struc_format_out, elements, iunartres, iunartout, &
                          lat, push
  implicit none

  ! -- Arguments
  CHARACTER (LEN=255), INTENT(IN) :: filnres
  INTEGER, INTENT( IN ) :: nat
  INTEGER, intent( inout ) :: order(nat), ityp(nat)   !> We change them or use them

  ! -- Local variable
  LOGICAL :: file_exists
  INTEGER :: ios
  LOGICAL, intent( out ) :: ierr
  CHARACTER(LEN=255) :: fname


  ! ...Verify if the file exist

  INQUIRE( file = filnres, exist = file_exists )
  ierr = .NOT.file_exists

  IF ( file_exists ) THEN

     OPEN( UNIT = iunartres, FILE = filnres, ACTION="READ", FORM = 'formatted', STATUS = 'old', IOSTAT = ios)

     READ( iunartres, * ) linit, lperp, leigen, llanczos, lpush_over, lrelax, &
       iartn, istep, iinit, ieigen, iperp, ilanc, irelax, ismooth,   &
       ninit, neigen, nlanc, lanczos_max_size, nperp, nmin, nsaddle, &
       etot_init, &
       etot_step, tau_step, force_step, current_step_size, fpush_factor, &   !> Actual step
       eigenvec, H, Vmat, force_old, lowest_eigval, &
       etot_saddle, tau_saddle

     CLOSE ( UNIT = iunartres, STATUS = 'KEEP')

     !> Maybe initialize de_back if needed


     ! ...Read the initial configuration => push, tau_init

     print*, "* RESTART:: init_structure file exist: ", trim(initpfname)

     if( .not.allocated(tau_init) )allocate( tau_init, source=tau_step)
     SELECT CASE( struc_format_out )
       CASE( 'xsf' )
         fname = TRIM(initpfname)//"."//TRIM(struc_format_out)
         CALL read_xsf( lat, nat, tau_init, order, elements, ityp, push, fname )

       CASE( 'xyz' )
         fname = TRIM(initpfname)//"."//TRIM(struc_format_out)
         CALL read_xyz( lat, nat, tau_init, order, elements, ityp, push, fname )
     END SELECT

  ELSE

     WRITE(iunartout,*) "ARTn: restart file does not exist, exiting ..."

  ENDIF


END SUBROUTINE read_restart

