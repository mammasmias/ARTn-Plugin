
!> @author
!!   Matic Poberznik,
!!   Miha Gunde

SUBROUTINE perpforce( force, if_pos, push, fperp, fpara, nat )
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
  REAL(DP) :: a, b
  REAL(DP), EXTERNAL :: ddot,dnrm2

  ! calculate components parallel to the push
  !fpara(:,:) = ddot(3*nat,force(:,:),1,push(:,:),1) / ddot(3*nat,push(:,:),1,push(:,:),1) * push(:,:)
  a = ddot(3*nat,force(:,:),1,push(:,:),1)
  b = ddot(3*nat,push(:,:),1,push(:,:),1)
  fpara(:,:) = a / b * push(:,:)

  ! subtract them
  fperp(:,:) = force(:,:) - fpara(:,:)

  ! apply constraints
  IF ( ANY(if_pos(:,:) == 0)  ) fperp(:,:) = fperp(:,:)*if_pos(:,:) 


END SUBROUTINE perpforce


!! N. Salles
!subroutine splitfield( n, field, mask, fref, fperp, fpara )
subroutine field_split( n, field, mask, fref, fperp, fpara )
  !
  !> @brief Extract the parallel and perpendicular component of field 
  !!   followig a reference field (fref) according to a mask.
  !!   (Generalization of perpforce)
  !
  !> @param[in]     nat         number of point in the field 
  !! @param[in]     field       Field input
  !! @param[in]     mask        Constrain in field
  !! @param[in]     fref        Parallel direction field reference
  !! @param[out]    fperp       Perpendicular force field following dir field
  !! @param[out]    fpara       Parallel force field following dir field
  !
  USE units, only : DP
  USE artn_params, ONLY : iunartout, filout
  IMPLICIT NONE

  ! -- ARGUMENTS
  INTEGER,  INTENT(IN)     :: n
  REAL(DP), INTENT(IN)     :: field(*)
  REAL(DP), INTENT(IN)     :: fref(*)
  INTEGER,  INTENT(IN)     :: mask(*)
  REAL(DP), INTENT(OUT)    :: fpara(*)
  REAL(DP), INTENT(OUT)    :: fperp(*)

  ! -- LOCAL VARIABLE
  REAL(DP) :: a, b
  REAL(DP), EXTERNAL :: ddot,dnrm2

  ! calculate components parallel to the dir
  a = ddot( n, field, 1, fref, 1 )
  b = ddot( n, fref, 1, fref, 1 )
  fpara(1:n) = a / b * fref(1:n)

  ! subtract them
  fperp(1:n) = field(1:n) - fpara(1:n)

  ! apply constraints
  IF( ANY(mask(1:n) == 0) )then
     fperp(1:n) = fperp(1:n)*mask(1:n)
     fpara(1:n) = fpara(1:n)*mask(1:n)
  endif


end subroutine field_split







