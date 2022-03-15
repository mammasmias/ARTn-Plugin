

!SUBROUTINE start_guess( idum, nat, order, mask, push, eigenvec )
SUBROUTINE start_guess( idum, nat, order, push, eigenvec )
  !
  !> @brief
  !!    Initialize the push and eigenvec arrays following the 
  !!    mode keyword
  !
  !> @param[in]  idum       seed for random number 
  !> @param[in]  nat        number of point 
  !! @param[in]  order      atom order of engine  
  !! @param[in]  push     
  !! @param[in]  eigenvec     
  !
  use units, only : DP
  use artn_params, only : push_mode, push_step_size, add_const, dist_thr, &
                          lat, tau_step, eigen_step_size, push_guess, eigenvec_guess, &
                          push_ids, iunartout, filout
  implicit none

  integer, intent(in) :: nat, idum
  integer, intent(in) :: order(nat)
  !integer, intent(in) :: mask(nat)
  real(DP), intent(out) :: push(3,nat)
  real(DP), intent(out) :: eigenvec(3,nat)

  logical :: verb

  verb = .true.
  !verb = .false.

  if( verb )OPEN ( UNIT = iunartout, FILE = filout, FORM = 'formatted', ACCESS = 'append', STATUS = 'unknown' )

  ! ...Define PUSH:
  !print*, "PUSH_MODE:", trim(push_mode)

  select case( trim(push_mode) )

    case( 'all', 'list', 'rad' )
      if( verb )write(iunartout,'(5x,"|> First PUSH vectors almost RANDOM")') 
      CALL push_init( nat, tau_step, order, lat, idum, push_ids, dist_thr, add_const, push_step_size, push, push_mode)

    case( 'file' )
      if( verb )write(iunartout,'(5x,"|> PUSH vectors read in file",x,a)') trim(push_guess)
      CALL read_guess( idum, nat, push, push_guess )

  end select
  

  ! ...Define EIGENVEC:
  !print*, "EIGENVEC_GUESS:", trim(eigenvec_guess), len_trim(eigenvec_guess)

  if( len_trim(eigenvec_guess) /= 0 )then
    !! read the file
    if( verb )write(iunartout,'(5x,"|> First EIGEN vectors read in file",x,a)') trim(eigenvec_guess)
    call read_guess( idum, nat, eigenvec, eigenvec_guess )
  else
    !! random 
    if( verb )write(iunartout,'(5x,"|> First EIGEN vectors RANDOM")') 
    add_const = 0
    call push_init( nat, tau_step, order, lat, idum, push_ids, dist_thr, add_const, eigen_step_size, eigenvec, 'all')
  endif


  if( verb )close(unit=iunartout, STATUS='KEEP')


END SUBROUTINE start_guess



