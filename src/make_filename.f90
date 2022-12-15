
!> @authors
!!  Matic Poberznik
!!  Miha Gune
!!  Nicolas Salles

!> @brief
!!    build a filename from the prefix and the number n
!
!> @param[out]    f        filename
!> @param[in]     prefix   prefix for filename
!> @param[inout]  n        integer for the file name
!
SUBROUTINE make_filename( f, prefix, n )
  !
  implicit none
  character(*), intent(out) :: f
  character(*), intent(in) :: prefix
  integer,      intent(inout) :: n

  character(len=4) :: ctmp
  character(len=256) :: fcounter
  integer :: o0

  n = n + 1
  write( ctmp, '(I0.4)') n
  f = trim(prefix)//trim(ctmp)

  fcounter = trim(prefix)//"counter"
  open( newunit=o0, file=fcounter, action="write" )
    write( o0, '(x,"counter:",2(x,a))' ) trim(prefix), trim(ctmp)
  close( o0 )

END SUBROUTINE make_filename
