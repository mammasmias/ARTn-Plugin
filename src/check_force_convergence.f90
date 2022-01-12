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
                          iperp, nperp, istep, INIT, PERP, EIGN, LANC, RELX, &
                          init_forc_thr, forc_thr, fpara_thr, push, &
                          lowest_eigval, iunartout, etot_step
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
  LOGICAL, parameter :: ArtnStep = .true.
  integer :: ios

  lforc_conv = .false.
  lsaddle_conv = .false.
  !

  OPEN ( UNIT = iunartout, FILE = 'artn.out', FORM = 'formatted', ACCESS = 'append', STATUS = 'unknown', IOSTAT = ios )

  !
  IF ( lperp ) THEN 
     !
     IF ( leigen ) THEN
        !
        ! We are outside of the basin check both perpendicular and total force convergence
        ! 
        IF (MAXVAL( ABS(force*if_pos)) < forc_thr  ) THEN
           lsaddle_conv = .true.
           CALL write_report( etot_step, force, fperp, fpara, lowest_eigval, EIGN, if_pos, istep, nat,  iunartout, ArtnStep )
        ENDIF
        !
        ! check whether the fperp criterion should be tightened
        ! 
        IF ( MAXVAL(ABS(fpara)) <= fpara_thr ) THEN
           fperp_thr = forc_thr
           nperp = 0    ! Remove the Perp-Relax Iteration constrain
        ELSE
           fperp_thr = init_forc_thr
        ENDIF
        ! 
        ! check perpendicular force convergence 
        ! 
        C1 = ( MAXVAL( ABS(fperp)) < fperp_thr ) ! check on the fperp field
        C2 = ( nperp > 0.AND.iperp >= nperp )    ! check on the perp-relax iteration

        IF( C1 .OR. C2  ) THEN
           lperp = .false.
           llanczos = .true.
           leigen = .false. 
           !if( C2 )  &    ! Because if C1 it already write_report before => If f < f_thr then fperp < f_thr
            CALL write_report( etot_step, force, fperp, fpara, lowest_eigval, EIGN, if_pos, istep, nat,  iunartout, ArtnStep )
        ENDIF

     ELSE ! ...IN  BASIN
        !
        ! In the Basin after perp-relax we always return to init mode
        !

        fperp_thr = init_forc_thr 

        ! ...Do INIT until fperp is > fperp_thr
        ! && Max Iteration of Perp-Relax
        !    nperp = 0 means no nperp constrain

        C1 = ( MAXVAL( ABS(fperp)) < fperp_thr ) ! check on the fperp field
        C2 = ( nperp > 0.AND.iperp >= nperp )    ! check on the perp-relax iteration

        IF( C1 .OR. C2 ) THEN
           lperp = .false.
           linit = .true.
           CALL write_report( etot_step, force, fperp, fpara, lowest_eigval, INIT, if_pos, istep, nat,  iunartout, ArtnStep )
        ENDIF

     ENDIF
     ! 
  ELSE IF ( lrelax ) THEN  
     !
     IF( MAXVAL( ABS(force*if_pos)) < forc_thr  )then
       lforc_conv = .true.
       CALL write_report( etot_step, force, fperp, fpara, lowest_eigval, RELX, if_pos, istep, nat,  iunartout, ArtnStep )
     ENDIF
     !
  ENDIF


  close( iunartout )
     
ENDSUBROUTINE check_force_convergence
