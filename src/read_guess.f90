

MODULE TOOLS

  use units, only : DP
  implicit none

 CONTAINS

  ! .............................................................................
  subroutine read_line(fd, line, end_of_file)
    !--------------------
    ! read a line, makes possible to use # for comment lines, skips empty lines, 
    !  is pretty much a copy from QE.
    !---------
    ! fd ==> file descriptor
    ! line ==> what it reads
    implicit none
    integer, intent(in) :: fd
    integer             :: ios
    character(len=256), intent(out) :: line
    logical, optional, intent(out) :: end_of_file
    logical :: tend

    !print*, "in read_line", fd
   
    tend = .false.
    101 read(unit=fd,fmt='(A256)',END=111, iostat=ios) line
    if (ios /= 0) then
       print*, " Reading Problem..."; stop; endif
    if(line == ' ' .or. line(1:1) == '#') go to 101
    go to 105
    111     tend = .true.
    !print*,"read_line", line
    go to 105
    105   continue

    if( present(end_of_file)) then
      end_of_file = tend
    endif
  end subroutine read_line




  !............................................................
  integer function fparse(instrg, FS, args )result( nargs )
    !
    implicit none
 
    ! -- ARGUMENT
    CHARACTER(len=*),              intent( in ) :: instrg
    character(len=1),              intent( in ) :: FS
    CHARACTER(len=:), allocatable, intent( inout ) :: args(:)
 
    ! -- LOCAL VAR
    character(len=:), allocatable :: str
    character(len=25) :: mot
    integer :: idx,leng
 
    ! +++ Copy in local variable the input_string
    str = adjustl(instrg)
    nargs = 0
    !
    ! +++ Repeat for each field
    do
    !   +++ Verification the length of sentence
        leng = len_TRIM( str )
        if( leng == 0 )exit
 
    !   +++ Find the Field Separator
        idx = SCAN( str, FS )
 
    !   +++ extract the word
        if( idx == 0 )then
          mot = trim(str)
        else
          mot = str( :idx-1 )
        endif
 
    !   +++ Add the word in args
        nargs = nargs + 1
        if( nargs == 1 )then
          args = [ mot ]
        else
          args = [ args(:), mot ]
        endif
 
    !   +++ cut the word
        if( idx == 0 )exit
        !str = trim(str(idx+1:))
        str = adjustl(str(idx+1:))
    !
    enddo

  end function fparse


  !......................................................
  elemental FUNCTION is_numeric(string)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: string
    LOGICAL :: is_numeric
    REAL :: x
    INTEGER :: e,n
    CHARACTER(len=12) :: fmt

    n = LEN_TRIM(string)
    WRITE(fmt,'("(F",I0,".0)")') n
    READ(string,fmt,IOSTAT=e) x
    is_numeric = (e == 0)
  END FUNCTION is_numeric


  subroutine random_displacement( idum, id, vec )

    use units, only : DP
    use artn_params, only : ran3
    implicit none

    integer, intent(in) :: id, idum
    real(DP), intent(inout ) :: vec(3)

    real(DP) :: dr
    real(DP), external :: dnrm2

    RDM:DO
       vec(:) = (/ 0.5_DP - ran3(idum), 0.5_DP - ran3(idum), 0.5_DP - ran3(idum) /)
       dr = dnrm2( 3, vec, 1 )
       IF ( dr < 0.25_DP ) RETURN
    ENDDO RDM

  end subroutine random_displacement


  subroutine neigh_random_displacement( idum, nat, id, rcut, vec )

    use units, only : DP, unconvert_length
    use artn_params, only : lat, tau_step, push_ids
    implicit none

    integer, intent( in ) :: idum, id, nat
    real(DP), intent( in ) :: rcut
    real(DP), intent( out ) :: vec(3,nat)

    integer :: na
    real(DP) :: x0(3), dr(3), d, rc
    real(DP), external :: dnrm2

    !
    ! -- WARNING : The position and lattice are stil in engine units
    !       we unconvert the rcut which should be in 
    !
    rc = unconvert_length( rcut )

    x0 = tau_step(:,id)
    DO na = 1,nat
       IF( id == na)cycle
       !IF( ANY(push_ids == na) )cycle
       dr(:) = tau_step(:,na) - x0(:)

       CALL pbc( dr, lat)
       d = dnrm2(3,dr,1) 
       IF( d <= rc )THEN
         ! found an atom within dist_thr 
         call random_displacement( idum, na, vec(:,na)) 
         !print*, id, na, d, "neigh random disp:", vec(:,na)
       ENDIF
    ENDDO


  end subroutine neigh_random_displacement


END MODULE TOOLS










!.....................................................................................................
SUBROUTINE READ_GUESS( idum, nat, vec, filename )

  use units,       only : DP, unconvert_length
  use artn_params, only : warning, iunartout, dist_thr, push_ids
  use tools
  implicit none

  integer,      intent( in ) :: nat, idum
  REAL(DP),     intent( out ) :: vec(3,nat)
  character(*), intent( in ) :: filename

  character(len=256) :: line
  character(:), allocatable :: words(:)
  integer :: i, n, u0, nwords, idx, j
  logical :: ok, neiglist

  !PRINT*, "   ** ENTER IN READ_GUESS()"

  ! ...Look at the file

  inquire( file=filename, exist=ok )
  if( .not.ok )CALL warning( iunartout, 'READ_GUESS','File does not find')



  !  ...Initialization

  if( allocated(push_ids) )deallocate(push_ids)
  neiglist = .false.
  if( dist_thr > 0.0e-8 ) neiglist = .true.
  !print*, "DIST_THR", dist_thr, unconvert_length( dist_thr )



  ! ...Read file

  OPEN( newunit=u0, file=filename, ACTION="READ" )

  READ(u0,*) n
  allocate(push_ids(n))
  READ(u0,*)

  do i = 1, n

     idx = 0
     !read(u0,*) line
     call read_line( u0, line )

     !print*, "line", trim(line)
     nwords = fparse( trim(line), " ", words )
     !print*, nwords, "fparse", ("|",j,words(j)," ",j=1,nwords)

     select case( nwords )

       case( 1 )
         IF( is_numeric(words(1)) )read(words(1),*) idx
         push_ids(i) = idx
         call random_displacement( idum, idx, vec(:,idx) )
         !print*, idx, "random disp:", vec(:,idx)

       case( 2: )
         IF( is_numeric(words(1)) )then
           read(words(1),*) idx
         else
           call warning( iunartout, 'READ_GUESS', 'index  proposed are not valid', words )
         endif
         push_ids(i) = idx

         !print*, "   ** push_ids", idx
         do j = 2,4
            !print*, j, "is num", is_numeric(words(j))
            IF( is_numeric(trim(words(j))) )then
              read(words(j),*) vec(j-1,idx)
              !print*, "read", j, vec(j-1,idx)
            else
              call warning( iunartout, 'READ_GUESS', 'Displacement propose are not valid', words )
            endif
         enddo   
         !print*, idx, "constrain disp:", vec(:,idx)

       case default
         call warning( iunartout, 'READ_GUESS', 'Empty line' )
         exit

     end select

     if( neiglist )call neigh_random_displacement( idum, nat, idx, dist_thr, vec )

  enddo

  CLOSE( u0 )

END SUBROUTINE READ_GUESS






