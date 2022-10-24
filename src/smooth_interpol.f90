
!> @author
!!  Matic Poberjnik,
!!  Miha Gunde
!!  Nicolas Salles
!
SUBROUTINE smooth_interpol( ismooth, nsmooth, nat, v0, v1, v2 )
  !
  !> @brief 
  !!   Return a smooth interpolation v2 betwwen 2 field v1 and v2:
  !!   v2 is a linear combination between v1 and v2
  !!   v2= v1 when ismooth = 0       -> done in init
  !!   v2= V2 when ismooth = nsmooth -> done in eigen 
  !
  !> @param[in]      ismooth  actual smooth step 
  !> @param[in]      nsmooth  max smooth step 
  !> @param[in]      nat      number of atom
  !> @param[in]      v0       the actual orientation
  !> @param[in]      v1       the direction we come
  !> @param[inout]   v2       the direction we go
  !
  USE units,       ONLY : DP
  USE artn_params, ONLY : iunartout,  dot_field, filout, verbose
  IMPLICIT NONE
  !
  INTEGER,  INTENT( INOUT ) :: ismooth   ! degree of interpolation
  INTEGER,  INTENT( IN )    :: nsmooth   ! number of degree of interpolation
  INTEGER,  INTENT( IN )    :: nat       ! number of points in 3D field
  REAL(DP), INTENT( IN )    :: v0(3,nat) ! Actuel field
  REAL(DP), INTENT( INOUT ) :: v1(3,nat) ! Orientation field 1
  REAL(DP), INTENT( IN )    :: v2(3,nat) ! Orientation field 2
  !
  ! Local variables
  REAL(DP)                  :: smoothing_factor, f_orient
  REAL(DP), external        :: ddot
  INTEGER                   :: ios

  ! save variable 
  REAL(DP), allocatable, save :: Vi(:,:), Vf(:,:)

  !
  ! ...At the first step we save the last push and the direction
  !      we want to go smoothly
  if( ismooth == 1 )then
    Vi = v1
    Vf = v2
  endif

  smoothing_factor = 1.0_DP*ismooth / (nsmooth+1)

  !
  ! ...Define the actual oriention from the final direction
  f_orient = ddot( 3*nat, v0, 1, vf, 1 )

  !
  ! ...Interpolation between Vi (initial) and Vf (final)
  v1 = (1.0_DP - smoothing_factor) * Vi    &
        - SIGN(1.0_DP,f_orient) * smoothing_factor * Vf

  !
  ! ...Info Output                
  IF (verbose>1) THEN              
    OPEN  (UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'unknown', POSITION='append', IOSTAT=ios)
    WRITE (iunartout,'(5x,a23,x,i2,a1,i2,x,a7,x,f15.9)')&
           "|> Smooth interpolation", ismooth,"/",nsmooth, "factor=",smoothing_factor
    CLOSE(iunartout)
  ENDIF
  !
END SUBROUTINE smooth_interpol



