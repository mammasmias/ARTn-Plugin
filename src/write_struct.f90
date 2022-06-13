
!> @author
!!   Matic Poberznik,
!!   Miha Gunde,
!!   Nicolas Salles


SUBROUTINE write_struct( lat, nat, tau, order, atm, ityp, force, ener, fscale, ounit, form, fname )
  !
  !> @brief
  !!   A subroutine that writes the structure to a file (based on xsf_struct of QE)
  !
  !> @param [in]  nat       number of atoms
  !> @param [in]  ityp      atom type
  !> @param [in]  order     atom type
  !> @param [in]  atm       contains information on atomic types
  !> @param [in]  ounit     output fortran unit
  !> @param [in]  tau       atomic positions
  !> @param [in]  at        lattice parameters in alat units
  !> @param [in]  force     list of atomic forces
  !> @param [in]  ene       energy of the current structure, in engine units
  !> @param [in]  fscale    factor for scaling the force
  !> @param [in]  form      format of the structure file (default xsf)
  !> @param [in]  fname     file name
  !
  !USE artn_params, ONLY: DP
  USE UNITS
  IMPLICIT NONE
  ! -- Arguments
  INTEGER,          INTENT(IN) :: nat            !> number of atoms
  INTEGER,          INTENT(IN) :: ityp(nat)      !> atom type
  INTEGER,          INTENT(IN) :: order(nat)     !> atom type
  CHARACTER(LEN=3), INTENT(IN) :: atm(*)         !> contains information on atomic types
  INTEGER,          INTENT(IN) :: ounit          !> output fortran unit
  REAL(DP),         INTENT(IN) :: tau(3,nat)     !> atomic positions
  REAL(DP),         INTENT(IN) :: lat(3,3)        !> lattice parameters in alat units
  REAL(DP),         INTENT(IN) :: force(3,nat)   !> list of atomic forces
  REAL(DP),         INTENT(IN) :: ener           !> energy of the structure, in engine units
  REAL(DP),         INTENT(IN) :: fscale         !> factor for scaling the force
  CHARACTER(LEN=3), INTENT(IN) :: form           !> format of the structure file (default xsf)
  !CHARACTER(LEN=255), INTENT(IN) :: fname        !> file name
  CHARACTER(*), INTENT(IN) :: fname        !> file name
  !
  ! -- Local Variables
  INTEGER ::  ios
  CHARACTER(:), ALLOCATABLE :: output

  ! ... Open the file with the good extention
  output = TRIM(fname)//"."//TRIM(form)
  OPEN ( UNIT = ounit, FILE = output, FORM = 'formatted',  STATUS = 'unknown', IOSTAT = ios )


  ! ... Select the format of the file
  SELECT CASE( form )

    CASE( 'xsf' )
      CALL write_xsf( lat, nat, tau, order, atm, ityp, force*fscale, ounit )

    CASE( 'xyz')
      CALL write_xyz( lat, nat, tau, order, atm, ityp, force*fscale, ounit, ener )

    CASE DEFAULT
      WRITE (ounit,*) " ** LIB::ARTn::WRITE_STRUC::Specified structure format not supported"

  END SELECT


  ! ... Close the file
  CLOSE (UNIT = ounit , STATUS = 'KEEP')
END SUBROUTINE write_struct


! .......................................................................................... XSF
SUBROUTINE write_xsf( lat, nat, tau, order, atm, ityp, force, ounit )

  USE UNITS
  !USE artn_params, ONLY: DP
  IMPLICIT NONE
  ! -- ARGUMENTS
  INTEGER,            INTENT(IN) :: nat            !> number of atoms
  INTEGER,            INTENT(IN) :: ityp(nat)      !> atom type
  INTEGER,            INTENT(IN) :: order(nat)     !> atom type
  CHARACTER(LEN=3),   INTENT(IN) :: atm(*)         !> contains information on atomic types
  INTEGER,            INTENT(IN) :: ounit          !> output fortran unit
  REAL(DP),           INTENT(IN) :: tau(3,nat)     !> atomic positions
  REAL(DP),           INTENT(IN) :: lat(3,3)        !> lattice parameters in alat units
  REAL(DP),           INTENT(IN) :: force(3,nat)   !> forces
  ! -- LOCAL VARIABLES
  INTEGER :: na, iloc
  REAL(DP) :: at_angs(3,3)

  !at_angs = unconvert_length( lat )  !!LAT is not converted

  WRITE(ounit,*) 'CRYSTAL'
  WRITE(ounit,*) 'PRIMVEC'
  !WRITE(ounit,'(2(3F15.9/),3f15.9)') at_angs
  WRITE(ounit,'(2(3F15.9/),3f15.9)') lat*B2A
  WRITE(ounit,*) 'PRIMCOORD'
  WRITE(ounit,*) nat, 1

  DO na=1,nat
     ! convert positions are in Engine units length
     ! -> And the force???
     iloc = order(na)
     WRITE(ounit,'(a3,3x,6f15.9)') atm(ityp(iloc)), tau(:,iloc)*B2A , unconvert_force( force(:,iloc) )
     !WRITE(ounit,'(a3,3x,6f15.9)') atm(ityp(iloc)), unconvert_length( tau(:,iloc) ), unconvert_force( force(:,iloc) )
     !WRITE(ounit,'(a3,3x,6f15.9)') ityp(iloc), unconvert_length( tau(:,iloc) ), unconvert_force( force(:,iloc) )
  ENDDO

END SUBROUTINE write_xsf


SUBROUTINE read_xsf( lat, nat, tau, order, atm, ityp, force, fname )

  USE UNITS
  implicit none

  ! -- ARGUMENTS
  INTEGER,            INTENT(IN) :: nat            !> number of atoms
  INTEGER,            INTENT(IN) :: ityp(nat)      !> atom type
  INTEGER,            INTENT(IN) :: order(nat)     !> atom type
  CHARACTER(LEN=3),   INTENT(OUT) :: atm(*)         !> contains information on atomic types
  REAL(DP),           INTENT(OUT) :: tau(3,nat)     !> atomic positions
  REAL(DP),           INTENT(OUT) :: lat(3,3)        !> lattice parameters in alat units
  REAL(DP),           INTENT(OUT) :: force(3,nat)   !> forces
  CHARACTER(*),       INTENT(IN) :: fname           !> file name
  ! -- LOCAL VARIABLES
  INTEGER :: na, u0, ios, iloc
  REAL(DP) :: at_angs(3,3)

  OPEN( newunit=u0, file=fname)

    READ( u0,* )
    READ( u0,* )
    READ( u0,* ) lat(:,1)
    READ( u0,* ) lat(:,2)
    READ( u0,* ) lat(:,3)
    READ( u0,* )
    READ( u0,* ) na, ios

    IF( na /= nat )print*, "* PROBLEM IN READ_XSF:: Different number of atoms", nat, na

    DO na=1,nat
       iloc = order(na)
       READ( u0,* ) atm(ityp(iloc)), tau(:,iloc), force(:,iloc)
    ENDDO
    lat   = lat*1/B2A
    tau   = tau*1/B2A
    force = convert_force( force )

  CLOSE( u0 )

END SUBROUTINE read_xsf




! .......................................................................................... XYZ
SUBROUTINE write_xyz( at, nat, tau, order, atm, ityp, f, ounit, ener )

  USE UNITS
  !USE artn_params, ONLY: DP
  IMPLICIT NONE
  ! -- ARGUMENTS
  INTEGER,            INTENT(IN) :: nat            !> number of atoms
  INTEGER,            INTENT(IN) :: ityp(nat)      !> atom type
  INTEGER,            INTENT(IN) :: order(nat)     !> atom type
  CHARACTER(LEN=3),   INTENT(IN) :: atm(*)         !> contains information on atomic types
  INTEGER,            INTENT(IN) :: ounit          !> output fortran unit
  REAL(DP),           INTENT(IN) :: tau(3,nat)     !> atomic positions
  REAL(DP),           INTENT(IN) :: at(3,3)        !> lattice parameters in alat units
  REAL(DP),           INTENT(IN) :: f(3,nat)       !> forces
  REAL(DP),           INTENT(IN) :: ener
  ! -- LOCAL VARIABLES
  INTEGER :: na, iloc, ios

  WRITE(ounit,*) nat
  WRITE(ounit,fmt=11) 'Lattice="',at(:,:),'"', &
       ' properties=species:S:1:pos:R:3:force:R:3:id:I:1',' energy:',ener
  11 format(a,x,9(f10.4,x),a,a,a,f15.9)

  DO na=1,nat
     ! convert positions are in Engine units length
     ! -> And the force???
     iloc = order(na)
     !WRITE( ounit, fmt='(a3,3x,6f15.9)', IOSTAT=ios ) atm(ityp(iloc)), unconvert_length( tau(:,iloc) ), unconvert_force( f(:,iloc) )
     !WRITE( ounit, fmt=10, IOSTAT=ios ) iloc, ityp(na), unconvert_length( tau(:,iloc) ), unconvert_force( f(:,iloc) )
     WRITE( ounit, fmt=10, IOSTAT=ios ) ityp(na), tau(:,iloc) , unconvert_force( f(:,iloc) ), iloc
     10 format(i2,3x,6f15.9,x,i0)
  ENDDO

END SUBROUTINE write_xyz


SUBROUTINE read_xyz( lat, nat, tau, order, atm, ityp, force, fname )

  USE UNITS
  implicit none

  ! -- ARGUMENTS
  INTEGER,            INTENT(IN) :: nat            !> number of atoms
  INTEGER,            INTENT(IN) :: ityp(nat)      !> atom type
  INTEGER,            INTENT(IN) :: order(nat)     !> atom type
  CHARACTER(LEN=3),   INTENT(OUT) :: atm(*)         !> contains information on atomic types
  REAL(DP),           INTENT(OUT) :: tau(3,nat)     !> atomic positions
  REAL(DP),           INTENT(OUT) :: lat(3,3)        !> lattice parameters in alat units
  REAL(DP),           INTENT(OUT) :: force(3,nat)   !> forces
  CHARACTER(*),       INTENT(IN) :: fname           !> file name
  ! -- LOCAL VARIABLES
  INTEGER :: na, u0, ios, iloc, i
  !REAL(DP) :: at_angs(3,3)

  OPEN( newunit=u0, file=fname)

    READ( u0,* ) na
    READ( u0,* )

    IF( na /= nat )print*, "* PROBLEM IN READ_XYZ:: Different number of atoms", nat, na

    DO na=1,nat
       iloc = order(na)
       !READ( u0,* ) atm(ityp(iloc)), tau(:,iloc), force(:,iloc)
       READ( u0,* ) i, ios, tau(:,iloc), force(:,iloc)
    ENDDO
    force = convert_force( force )

  CLOSE( u0 )

END SUBROUTINE read_xyz









