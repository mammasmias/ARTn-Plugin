
!> @author
!!   Matic Poberznik,
!!   Miha Gunde

SUBROUTINE perpforce( force, if_pos, push, fperp, fpara, nat )
!SUBROUTINE splitfield( force, if_pos, push, fperp, fpara, nat )
  !
  !> @brief subroutine that subtracts parallel components to push from force
  !
  !> @param[in]     force	Field input
  !! @param[in]	    if_pos	Constrain in field
  !! @param[in]	    push	Parallel field reference
  !! @param[out]    fperp	Perpendicular force field following Push field
  !! @param[out]    fpara	Parallel force field following Push field
  !! @param[in]	    nat		number of point in the field 
  !
  USE units, only : DP
  USE artn_params, ONLY : iunartout, filout
  IMPLICIT NONE

  ! -- ARGUMENTS
  INTEGER,  INTENT(IN)     :: nat
  REAL(DP), INTENT(IN)     :: push(3,nat)
  REAL(DP), INTENT(IN)     :: force(3,nat)
  REAL(DP), INTENT(OUT)    :: fpara(3,nat)
  REAL(DP), INTENT(OUT)    :: fperp(3,nat)
  INTEGER,  INTENT(IN)     :: if_pos(3,nat)

  ! -- LOCAL VARIABLE
  REAL(DP) :: push_norm(3,nat), a, b
  REAL(DP), EXTERNAL :: ddot,dnrm2
  INTEGER  :: na

  ! calculate components parallel to the push
  !fpara(:,:) = ddot(3*nat,force(:,:),1,push(:,:),1) / ddot(3*nat,push(:,:),1,push(:,:),1) * push(:,:)
  a = ddot(3*nat,force(:,:),1,push(:,:),1)
  b = ddot(3*nat,push(:,:),1,push(:,:),1)
  fpara(:,:) = a / b * push(:,:)

  ! subtract them
  fperp(:,:) = force(:,:) - fpara(:,:)

  ! apply constraints
  IF ( ANY(if_pos(:,:) == 0)  ) fperp(:,:) = fperp(:,:)*if_pos(:,:) 

  !open( unit=iunartout, file=filout, FORM = 'formatted', ACCESS = 'append', STATUS = 'unknown' )
  !write(iunartout,*)"DEBUG::PERPF-coeff", a, b, a/b
  !write(iunartout,*)"DEBUG::PERPF-MAX", maxval(abs(force)), maxval(abs(fpara)), maxval(abs(fperp))
  !close( unit=iunartout, status='KEEP')

END SUBROUTINE perpforce
!END SUBROUTINE splitfield

