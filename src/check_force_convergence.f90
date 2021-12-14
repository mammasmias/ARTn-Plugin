SUBROUTINE check_force_convergence(nat,force, if_pos, fperp, fpara, lforc_conv, lsaddle_conv) 
  !
  ! A subroutine that checks the force convergence of a particular step in the artn algorithm 
  !
  USE units
  USE artn_params, ONLY : leigen, lperp, lrelax, init_forc_thr, forc_thr, fpara_thr, push
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: force(3,nat)
  REAL(DP), INTENT(IN) :: fperp(3,nat)
  REAL(DP), INTENT(IN) :: fpara(3,nat)
  INTEGER, INTENT(IN) :: if_pos(3,nat)
  INTEGER, INTENT(IN) :: nat
  REAL(DP) :: fperp_thr
  LOGICAL, INTENT(OUT) :: lforc_conv, lsaddle_conv
  !
  lforc_conv = .false.
  lsaddle_conv = .false.
  !
  !
  IF ( lperp ) THEN 
     !
     IF ( leigen ) THEN
        !
        ! We are outside of the basin check both perpendicular and total force convergence
        ! 
        IF (MAXVAL( ABS(force*if_pos)) < forc_thr  ) lsaddle_conv = .true.
        !
        ! check whether the fperp criterion should be tightened
        ! 
        IF ( MAXVAL(ABS(fpara)) <= fpara_thr ) THEN
           fperp_thr = forc_thr
        ELSE
           fperp_thr = init_forc_thr
        ENDIF
        ! 
        ! check perpendicular force convergence 
        ! 
        IF ((MAXVAL( ABS(fperp))) < fperp_thr  ) lforc_conv = .true.
     ELSE
        fperp_thr = init_forc_thr 
        IF ((MAXVAL( ABS(fperp))) < fperp_thr  ) lforc_conv = .true.
     ENDIF
     ! 
  ELSE IF ( lrelax ) THEN  
     !
     write (*,*) "Debug conv check at lrelax", lforc_conv, forc_thr
     IF ((MAXVAL( ABS(force*if_pos))) < forc_thr  ) lforc_conv = .true.
     !
  ENDIF
     
ENDSUBROUTINE check_force_convergence
