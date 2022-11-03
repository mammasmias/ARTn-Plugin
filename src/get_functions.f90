
!> @authors
!!  Matic Poberznik
!!  Miha Gunde
!!  Nicolas Salles

integer function get_iperp()
  !> @brief 
  !!   give the parameters IPERP
  !> @return  iperp
  USE artn_params, only : iperp
  get_iperp = iperp
end function get_iperp
  

integer function get_perp()
  !> @brief 
  !!   give the parameters PERP
  !> @return  PERP
  USE artn_params, only : perp
  get_perp = perp
end function get_perp


integer function get_relx()
  !> @brief 
  !!   give the parameters RELX
  !> @return RELX
  USE artn_params, only : relx
  get_relx = relx
end function get_relx

integer function get_irelx()
  !> @brief 
  !!   give the parameters IRELX
  !> @return  iperp
  USE artn_params, only : irelax
  get_irelx = irelax
end function get_irelx
