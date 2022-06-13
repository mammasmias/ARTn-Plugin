SUBROUTINE start_guess( idum, nat, order, force, push, eigenvec )
  !
  !> @brief
  !!    Initialize the push and eigenvec arrays following the mode keyword
  !!
  !! MIHA
  !! use force input as mask for push_ids when calling push_init for eigenvec.
  !! Why? To not generate initial lanczos vec for fixed atoms.
  !
  !> @param[in]  idum       seed for random number
  !> @param[in]  nat        number of point
  !! @param[in]  order      atom order of engine
  !! @param[in]  push
  !! @param[in]  eigenvec
  !
  USE units,       ONLY : DP
  USE artn_params, ONLY : push_mode, push_step_size, add_const, dist_thr,             &
                          lat, tau_step, eigen_step_size, push_guess, eigenvec_guess, &
                          push_ids, iunartout, filout, verbose
  USE tools,       ONLY : read_guess
  !
  IMPLICIT NONE
  ! 
  ! Arguments
  INTEGER,  INTENT(IN)  :: nat, idum
  INTEGER,  INTENT(IN)  :: order(nat)
  REAL(DP), INTENT(IN)  :: force(3,nat)
  REAL(DP), INTENT(OUT) :: push(3,nat)
  REAL(DP), INTENT(OUT) :: eigenvec(3,nat)
  !
  ! Local variables
  INTEGER               :: mask(nat)
  INTEGER               :: i, j
  !
  IF( verbose >2 ) OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', ACCESS = 'append', STATUS = 'unknown' )
  !
  SELECT CASE( TRIM(push_mode) )
    !
    CASE( 'all', 'list', 'rad' )
       ! 
       IF( verbose>2 ) WRITE(iunartout,'(5x,"|> First PUSH vectors almost RANDOM")')
       CALL push_init( nat, tau_step, order, lat, idum, push_ids, dist_thr, add_const, push_step_size, push, push_mode)
       !
    CASE( 'file' )
       ! 
       IF( verbose >2 ) WRITE(iunartout,'(5x,"|> PUSH vectors read in file",x,a)') TRIM(push_guess)
       CALL read_guess( idum, nat, push, push_guess )
       !
  END SELECT
  !
  ! ...Define EIGENVEC:
  IF( LEN_TRIM(eigenvec_guess) /= 0 ) THEN
    !! read the file
    IF( verbose>2 ) WRITE(iunartout,'(5x,"|> First EIGEN vectors read in file",x,a)') TRIM(eigenvec_guess)
    CALL read_guess( idum, nat, eigenvec, eigenvec_guess )
  ELSE
    !! random
    IF( verbose>2 ) WRITE(iunartout,'(5x,"|> First EIGEN vectors RANDOM")')
    add_const = 0
    !! set up the mask according to input forces
    mask(:) = 0
    j = 1
    DO i = 1, nat
       !! Atoms with norm of the force > 1e-16 are most probably not fixed by the engine, use the array 'mask'
       !! as list of push_ids to generate the initial eigenvec components.
       !! This avoids generating components on atoms that are fixed.
       IF( NORM2(force(:,i)) .GT. 1e-16 ) THEN
          mask(j) = i
          j = j + 1
       ENDIF
       !
    ENDDO
    call push_init( nat, tau_step, order, lat, idum, mask, dist_thr, add_const, eigen_step_size, eigenvec, 'list')
    !
  ENDIF
  !
  IF( verbose>2 ) CLOSE(UNIT=iunartout, STATUS='KEEP')
  !
END SUBROUTINE start_guess
