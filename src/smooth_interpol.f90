


SUBROUTINE smooth_interpol( ismooth, nsmooth, nat, v0, v1, v2 )

  !> @breif 
  !!   Return a smooth interpolation v2 betwwen 2 field v1 and v2
  !!   according the orientation of field v0
  !
  !> @detail 
  !!   When the nsmooth are done this step serve to align in opposite send  
  !!   the field v2 with the direction of the field v0
  !
  !> @param[inout]  ismooth    degree of interpolation
  !! @param[in]     nsmooth    number of degree of interpolation
  !! @param[in]     nat        number of point in 3D field
  !! @param[in]     v0         Actual field
  !! @param[inout]  v1         orientation field 1
  !! @param[inout]  v2         orientation field 2

  use units, only : DP
  use artn_params, only : iunartout, force_step, eigenvec, dot_field

  integer, intent( inout )  :: ismooth
  integer, intent( in )     :: nsmooth, nat
  real(DP), intent( in )    :: v0(3,nat)
  real(DP), intent( in )    :: v1(3,nat)
  real(DP), intent( inout ) :: v2(3,nat)

  real(DP) ::  smoothing_factor, f_orient, vtmp(3,nat)
  REAL(DP), external :: ddot
  !
 
  smoothing_factor = 1.0_DP*ismooth/nsmooth
  !
  f_orient = ddot(3*nat, v0(:,:), 1, v2(:,:), 1)

  if( ismooth <= nsmooth )then
    vtmp(:,:) = (1.0_DP - smoothing_factor) * v1(:,:) &
            - SIGN(1.0_DP,f_orient) * smoothing_factor * v2(:,:)
  else
    vtmp = - SIGN(1.0_DP,f_orient) * v2(:,:)
  endif
  v2 = vtmp
  !write(iunartout,*)"DEBUG::SMOOTH_INTERPOL-3", ismooth, nsmooth, smoothing_factor, f_orient, &
  !       (1.0_DP - smoothing_factor), - SIGN(1.0_DP,f_orient)
  !
  !IF( ismooth < nsmooth) ismooth = ismooth + 1
  if( ismooth > nsmooth )return
  ismooth = ismooth + 1


END SUBROUTINE smooth_interpol

!
!! SAVE from ARTN() body
!
     !>>>>>>>>>>>>>>>
  !  IF( nsmooth > 0 )THEN

  !    smoothing_factor = 1.0_DP*ismooth/nsmooth
  !    !!
  !    fpara_tot = ddot(3*nat, force_step(:,:), 1, eigenvec(:,:), 1)

  !    eigenvec(:,:) = (1.0_DP - smoothing_factor)*push(:,:) &
  !         -SIGN(1.0_DP,fpara_tot)*smoothing_factor*eigenvec(:,:)
  !    !
  !    IF( ismooth < nsmooth)then
  !      ismooth = ismooth + 1
  !    ELSE
  !      ! ...Save the eigenvector
  !      push(:,:) = eigenvec(:,:)
  !      write(iunartout,'(x,"DEBUG::EIGEN::Overwrite push = eigenvec")')
  !    ENDIF

  !  ELSE
  !    ! ...Save the eigenvector
  !     push(:,:) = eigenvec(:,:)
  !  ENDIF
     !<<<<<<<<<<<<<<<


