


!Apply_Smooth_interpol( ismooth, push, eigenvec )
SUBROUTINE smooth_interpol( ismooth, nat, v0, v1, v2, f_orient )

!> @breif Return a smooth interpolation v2 betwwen 2 field v1 and v2
!! according the orientation of field v0

!> @param[inout]  ismooth    degree of interpolation
!! @param[in]     nat        number of point in 3D field
!! @param[in]     v0         Actual field
!! @param[in]     v1         orientation field 1
!! @param[inout]  v2         orientation field 2
!! @param[inout]  f_orient   Oritentation of v0 on v2

  use units
  use artn_params, only : nsmooth

  integer, intent( inout ) :: ismooth
  real(DP), intent( in ) :: v0(3,nat), v1(3,nat)
  real(DP), intent( inout ) :: v2(3,nat)

  real(DP) ::  smoothing_factor, f_orient

  smoothing_factor = 1.0_DP*ismooth/nsmooth
  !
  f_orient = ddot(3*nat, v0(:,:), 1, v2(:,:), 1)

  v2(:,:) = (1.0_DP - smoothing_factor) * v1(:,:) &
          - SIGN(1.0_DP,f_orient) * smoothing_factor * v2(:,:)
  !
  IF ( ismooth < nsmooth) ismooth = ismooth + 1


END SUBROUTINE smooth_interpol
