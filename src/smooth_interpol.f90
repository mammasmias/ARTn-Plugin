SUBROUTINE smooth_interpol( ismooth, nsmooth, nat, v0, v1, v2 )
  !
  !> @breif 
  !  Return a smooth interpolation v2 betwwen 2 field v1 and v2:
  !  v2 is a linear combination between v0 and v1
  !  v2= v0 when ismooth = 0       -> done in init
  !  v2= V1 when ismooth = nsmooth -> done in eigen 
  !
  USE units,       ONLY : DP
  USE artn_params, ONLY : iunartout,  dot_field,filout, verbose
  !
  INTEGER,  INTENT( INOUT ) :: ismooth   ! degree of interpolation
  INTEGER,  INTENT( IN )    :: nsmooth   ! number of degree of interpolation
  INTEGER,  INTENT( IN )    :: nat       ! number of points in 3D field
  REAL(DP), INTENT( IN )    :: v0(3,nat) ! Actuel field
  REAL(DP), INTENT( IN )    :: v1(3,nat) ! Orientation field 1
  REAL(DP), INTENT( INOUT ) :: v2(3,nat) ! Orientation field 2
  !
  ! Local variables
  REAL(DP)                  :: smoothing_factor, f_orient
  REAL(DP), external        :: ddot
  INTEGER                   :: ios
  !
  smoothing_factor = 1.0_DP*ismooth / (nsmooth+1)
  f_orient         = ddot(3*nat, v0(:,:), 1, v2(:,:), 1)
  v2(:,:)          = (1.0_DP - smoothing_factor) * v1(:,:) &
                     - SIGN(1.0_DP,f_orient) * smoothing_factor * v2(:,:)
  !               
  IF (verbose>1) THEN              
    OPEN  (UNIT = iunartout, FILE = filout, FORM = 'formatted', STATUS = 'unknown', POSITION='append', IOSTAT=ios)
    WRITE (iunartout,'(5x,a23,x,i2,a1,i2,x,a7,x,f15.9)')&
           "|> Smooth interpolation", ismooth,"/",nsmooth, "factor=",smoothing_factor
    CLOSE(iunartout)
  ENDIF
  !
END SUBROUTINE smooth_interpol
