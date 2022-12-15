
!> @authors
!!  Matic Poberznik,
!!  Miha Gunde,
!!  Nicolas Salles

  !> @brief 
  !!   give the parameters IPERP
  !> @return  iperp
integer function get_iperp()
  USE artn_params, only : iperp
  get_iperp = iperp
end function get_iperp
  

  !> @brief 
  !!   give the parameters PERP
  !> @return  PERP
integer function get_perp()
  USE artn_params, only : perp
  get_perp = perp
end function get_perp


  !> @brief 
  !!   give the parameters RELX
  !> @return RELX
integer function get_relx()
  USE artn_params, only : relx
  get_relx = relx
end function get_relx


  !> @brief 
  !!   give the parameters IRELX
  !> @return  iperp
integer function get_irelx()
  USE artn_params, only : irelax
  get_irelx = irelax
end function get_irelx
