


SUBROUTINE nperp_limitation_step( increment )
  !> @brief
  !!   increment in the list only if the actual perp-relax finished
  !
  !> @param[in]  increment   command {-1,0,1} allows to show what it does
  !
  use artn_params, only: nperp_limitation, nperp_step, nperp
  implicit none

  integer, intent( in ) :: increment

  integer :: n

  if( increment == 1 )then
    !! Save the basin nperp
    if( nperp_step == 1 )nperp_limitation(1) = nperp
    !! increment the limitation list
    nperp_step = nperp_step + increment
  endif

  n = size(nperp_limitation)
  !n = 5 !! By default

  !select case( nperp_step )
  !  case(:n-1); nperp = nperp_limitation( nperp_step )
  !  case(n:); nperp = nperp_limitation( n )
  !end select

  !nperp = merge( nperp_limitation( nperp_step ), nperp_limitation( n ), nperp_step < n )

  if( nperp_step < n )then
    nperp = nperp_limitation( nperp_step )
  else 
    nperp = nperp_limitation( n )
  endif

  if( increment == -1 )then
    nperp_step = 1
    nperp = nperp_limitation( nperp_step )
  endif

end SUBROUTINE nperp_limitation_step





!...........................................................................................
SUBROUTINE nperp_limitation_init( flag )

  use artn_params, only: nperp_limitation, nperp_step, nperp, &
                         def_nperp_limitation
  implicit none

  logical, intent( in ) :: flag

  integer :: i , n


  !! User says use nperp_limitation
  IF( flag )THEN

    !! but no defines the limitation
    !IF( .not.allocated(nperp_limitation) )THEN
    IF( ALL(nperp_limitation == -2) )THEN
      deallocate( nperp_limitation )
      allocate( nperp_limitation, source=def_nperp_limitation )

    !! define just one limitation
    ELSE
      n = 0
      do i = 1,size(nperp_limitation)
         if(nperp_limitation(i) > -2)n = n + 1 
      enddo
      !print*, "LENGTH:", n
      !! and also define nperp
      if( nperp /= -1 )then
        nperp_limitation = [ nperp, nperp_limitation(1:n), -1 ]
      else
        nperp_limitation = [ def_nperp_limitation(1), nperp_limitation(1:n), -1 ]
      endif
        
    ENDIF

  !! User does not use the nperp_limitation 
  ELSE
    !! We still define nperp_limitation but at nothing
    !nperp_limitation = [ 0, 0 ]
    nperp_limitation = [ -1, -1 ]

    !! But he still define is own nperp in basin
    if( nperp /= -1 )nperp_limitation(1) = nperp
  ENDIF

  !! Define nperp
  nperp = nperp_limitation(1)

  write(*,'(5x,"|> NPERP_LIMITATION:: Actual nperp",x,i0,/5x,"|> NPERP_LIMITATION::List:",*(x,i0))') &
      nperp, nperp_limitation(:)


end SUBROUTINE nperp_limitation_init

