
!> @author
!!  Matic Poberznik
!!  Miha Gunde
!!  Nicolas Salles


!SUBROUTINE start_guess( idum, nat, order, force, push, eigenvec )
SUBROUTINE start_guess( idum, nat, order, push, eigenvec )
  !
  !> @brief
  !!    Initialize the push and eigenvec arrays following the mode keyword
  !!
  !! MIHA <= Move in push_init
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
  !REAL(DP), INTENT(IN)  :: force(3,nat)
  REAL(DP), INTENT(OUT) :: push(3,nat)
  REAL(DP), INTENT(OUT) :: eigenvec(3,nat)
  !
  ! Local variables
  INTEGER               :: mask(nat)
  !
  IF( verbose >1 ) OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', ACCESS = 'append', STATUS = 'OLD' )
  !
  SELECT CASE( TRIM(push_mode) )
    !
    CASE( 'all', 'list', 'rad' )
       ! 
       IF( verbose>1 ) WRITE(iunartout,'(5x,"|> First PUSH vectors almost RANDOM")')
       CALL push_init( nat, tau_step, order, lat, idum, push_ids, dist_thr, add_const, push_step_size, push, push_mode)
       !
    CASE( 'file' )
       ! 
       IF( verbose >1 ) WRITE(iunartout,'(5x,"|> PUSH vectors read in file",x,a)') TRIM(push_guess)
       CALL read_guess( idum, nat, push, push_guess )
       !
  END SELECT
  !
  ! ...Define EIGENVEC:
  IF( LEN_TRIM(eigenvec_guess) /= 0 ) THEN

    !! read the file
    IF( verbose>1 ) WRITE(iunartout,'(5x,"|> First EIGEN vectors read in file",x,a)') TRIM(eigenvec_guess)
    CALL read_guess( idum, nat, eigenvec, eigenvec_guess )

  ELSE

    !! random
    IF( verbose>1 ) WRITE(iunartout,'(5x,"|> First EIGEN vectors RANDOM")')
    add_const = 0
    !! Replace Mask on norm(force) by keyword 'list_force'. 
    !! keyword 'bias_force' = orient the randomness on the actual atomic forces
    call push_init( nat, tau_step, order, lat, idum, mask, dist_thr, add_const, eigen_step_size, eigenvec, 'list_force')
    !call push_init( nat, tau_step, order, lat, idum, mask, dist_thr, add_const, eigen_step_size, eigenvec, 'bias_force' )

  ENDIF
  !
  IF( verbose>1 ) CLOSE(UNIT=iunartout, STATUS='KEEP')
  !
END SUBROUTINE start_guess


