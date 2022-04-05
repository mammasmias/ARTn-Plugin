


SUBROUTINE nperp_limitation( increment )
  !> @brief
  !!   increment in the list only if the actual perp-relax finished

  use artn_params, only: nperp_list, nperp_step, nperp
  implicit none

  integer, intent( in ) :: increment

  integer :: n

  if( increment == 1 )then
    !! Save the basin nperp
    if( nperp_step == 1 )nperp_list(1) = nperp
    !! increment the limitation list
    nperp_step = nperp_step + increment
  endif

  n = size(nperp_list)
  !n = 5 !! By default

  !select case( nperp_step )
  !  case(:n-1); nperp = nperp_list( nperp_step )
  !  case(n:); nperp = nperp_list( n )
  !end select

  !nperp = merge( nperp_list( nperp_step ), nperp_list( n ), nperp_step < n )

  if( nperp_step < n )then
    nperp = nperp_list( nperp_step )
  else 
    nperp = nperp_list( n )
  endif

  if( increment == -1 )then
    nperp_step = 1
    nperp = nperp_list( nperp_step )
  endif

end SUBROUTINE nperp_limitation



SUBROUTINE init_nperp_limitation( flag )

  use artn_params, only: nperp_list, nperp_step, nperp, &
                         def_nperp_list
  implicit none

  logical, intent( in ) :: flag

  !! User says use nperp_limitation
  IF( flag )THEN

    !! but no defines the limitation
    IF( .not.allocated(nperp_list) )THEN
      allocate( nperp_list, source=def_nperp_list )

    !! define just one limitation
    ELSE
      if( size(nperp_list) == 1 )then
        !! and also define nperp
        if( nperp /= -1 )then
          nperp_list = [ nperp, nperp_list(1), -1 ]
        else
          nperp_list = [ def_nperp_list(1), nperp_list(1), -1 ]
        endif
      endif
    ENDIF

  !! User does not use the nperp_limitation 
  ELSE
    !! We still define nperp_list but at nothing
    !nperp_list = [ 0, 0 ]
    nperp_list = [ 4, -1 ]

    !! But he still define is own nperp in basin
    if( nperp /= -1 )nperp_list(1) = nperp
  ENDIF

  !! Define nperp
  nperp = nperp_list(1)


end SUBROUTINE init_nperp_limitation


