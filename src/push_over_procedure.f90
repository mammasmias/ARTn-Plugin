

SUBROUTINE Push_Over_Procedure( iover, nat, pos, v0, push_factor, order, displ_vec, lstop )
  !
  !> @Brief
  !!    Perform the push over the saddle point
  !
  !> @param[inout] iover         iterator of push_over
  !> @param[in]    nat           number of atoms
  !> @param[out]   pos           atomic position
  !> @param[in]    v0            Vector defining the push over
  !> @param[in]    push_factor   +/- 1 depending the sens of the push 
  !> @param[out]   order         atoms engine order 
  !> @param[out]   displ_vec     displacement vector
  !> @param[out]   lstop         flag to stop the computation
  !
  use units, only : DP
  use artn_params, only : eigen_step_size, push_over, &
                   etot_step, etot_saddle, frelax_ene_thr, tau_saddle, tau_init, &
                   iunartout, OVER, warning, flag_false
  implicit none

  integer, intent(inout) :: iover
  integer, intent(in)    :: nat, order(nat)
  real(DP), intent(in)   :: push_factor, v0(3,nat)
  real(DP), intent(out)  :: pos(3,nat), displ_vec(3,nat)
  logical, intent(out) :: lstop

  real(DP) :: coeff

  lstop = .false.

  iover = iover + 1

  IF( iover > 1 )THEN
    pos(:,:) = tau_saddle(:,order(:))  ! no convertion needed
    !push_over = push_over * 0.80      !! Replace by ternary operator to don't change the value
  ENDIF

  ! Decrease the push_over factor from 0.8^(iover-1)
  coeff = push_over * 0.8**(iover-1)  !merge( 1.0, 0.8**real(iover-1), iover == 1)
  displ_vec(:,:) = push_factor * v0(:,:) * eigen_step_size * coeff

  print*, "PUSH_OVER:", iover, coeff, etot_step - etot_saddle


  ! ** WARNING **
  if( iover > 4 ) &
       CALL WARNING( iunartout, "PUSH_OVER_PROCEDURE()",&
       "Too many push over at saddle point: frelax_ene_thr can be too big or push_over", &
        [ etot_step - etot_saddle, frelax_ene_thr, coeff ] )

  ! ** ERROR **
  IF( iover > 10 )THEN
    call write_fail_report( iunartout, OVER, etot_step )
    call flag_false()

    ! ...Return to Starting configuration
    pos(:,:) = tau_init(:,order(:))
    displ_vec = 0.0_DP
    lstop = .true.
    !return
  ENDIF


END SUBROUTINE Push_Over_Procedure

!
!! SAVE from ARTn() Body
!

    !      !>>>>>>>>>>>>>>>>>>>>>> push_over_procedure()
    !      !! Idea: Push over first time and if does not work return to the saddle 
    !      !!  and do a smaller push. Doing that one or two times and stop the research
    !      !
    !      iover = iover + 1

    !      IF( iover > 1 )THEN
    !        tau(:,:) = tau_saddle(:,order(:))  ! no convertion needed
    !        !push_over = push_over * 0.80      !! Replace by ternary operator to don't change the value
    !      ENDIF

    !      !
    !      displ_vec(:,:) = fpush_factor*eigenvec(:,:)*eigen_step_size * push_over * merge( 1.0, 0.8**real(iover-1), iover == 1)


    !      ! ** WARNING **
    !      if( iover > 4 ) &
    !           CALL WARNING( iunartout, "PUSH_OVER_PROCEDURE()",&
    !           "Too many push over at saddle point: frelax_ene_thr can be too big or push_over", &
    !            [ etot_step, etot_saddle, etot_step - etot_saddle, frelax_ene_thr,   &
    !              push_over*merge( 1.0, 0.8**real(iover-1), iover == 1) ] )

    !      ! ** ERROR **
    !      IF( iover > 10 )THEN
    !        call write_fail_report( iunartout, OVER, etot_step )
    !        !call clean_artn()  !! Over Push_Over
    !        lrelax = .false.
    !        linit = .false.
    !        lbasin = .false.
    !        lperp = .false.
    !        llanczos = .false.
    !        leigen = .false.
    !        lsaddle = .false.

    !        ! ...Return to Starting configuration
    !        tau(:,:) = tau_init(:,order(:))
    !        displ_vec = 0.0_DP
    !        lconv = .true.
    !        !return
    !      ENDIF
    !      !<<<<<<<<<<<<<<<<<<<<<<
