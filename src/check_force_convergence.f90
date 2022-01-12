SUBROUTINE check_force_convergence( nat, force, if_pos, fperp, fpara, lforc_conv, lsaddle_conv) 
  !
  !> @breif 
  !!   A subroutine that checks the force convergence of a particular step in the artn algorithm 
  !
  !> @param [in]   nat             Size of list: number of atoms
  !> @param [in]   force           Force field
  !> @param [in]   if_pos          List of atom move or not
  !> @param [in]   fperp           Perpendicular Force Field
  !> @param [in]   fpara           Parallel Force Field
  !> @param [out]  lforc_conv      Force Convergence Flag
  !> @param [out]  lsaddke_conv    Saddle-point Convergence Flag
  !
  USE units
  USE artn_params, ONLY : linit, lbasin, leigen, llanczos, lperp, lrelax, &
                          iperp, nperp, &
                          init_forc_thr, forc_thr, fpara_thr, push
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: force(3,nat)
  REAL(DP), INTENT(IN) :: fperp(3,nat)
  REAL(DP), INTENT(IN) :: fpara(3,nat)
  INTEGER, INTENT(IN) :: if_pos(3,nat)
  INTEGER, INTENT(IN) :: nat
  REAL(DP) :: fperp_thr
  LOGICAL, INTENT(OUT) :: lforc_conv, lsaddle_conv
  !
  LOGICAL :: C1, C2

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
        IF (MAXVAL( ABS(force*if_pos)) < forc_thr  ) THEN
           lsaddle_conv = .true.
        ENDIF
        !
        ! check whether the fperp criterion should be tightened
        ! 
        IF ( MAXVAL(ABS(fpara)) <= fpara_thr ) THEN
           fperp_thr = forc_thr
           nperp = 0    !
        ELSE
           fperp_thr = init_forc_thr
        ENDIF
        ! 
        ! check perpendicular force convergence 
        ! 
        C1 = ( MAXVAL( ABS(fperp)) < fperp_thr ) ! check on the force field
        C2 = ( nperp > 0.AND.iperp >= nperp )    ! check on the perp-relax iteration

        IF( C1 .OR. C2  ) THEN
           lperp = .false.
           llanczos = .true.
           leigen = .false. 
        ENDIF

     ELSE
        !
        ! In the Basin after perp-relax we always return to init mode
        !

        fperp_thr = init_forc_thr 

        ! ...Do INIT until fperp is > fperp_thr
        IF ((MAXVAL( ABS(fperp))) < fperp_thr  ) THEN
           lperp = .false.
           linit = .true.

        ! ...Max Iteration of Perp-Relax
        ! OPTION::nperp = 0 means no nperp constrain
        ELSEIF( nperp > 0.AND.iperp >= nperp )THEN
           lperp = .false.
           linit = .true.

        ENDIF

     ENDIF
     ! 
  ELSE IF ( lrelax ) THEN  
     !
     IF ((MAXVAL( ABS(force*if_pos))) < forc_thr  ) lforc_conv = .true.
     !
  ENDIF
     
ENDSUBROUTINE check_force_convergence
