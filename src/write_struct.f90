
!> @author
!!   Matic Poberznik,
!!   Miha Gunde,
!!   Nicolas Salles


SUBROUTINE write_struct( at, nat, tau, order, atm, ityp, force, fscale, ounit, form, fname)
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
  REAL(DP),         INTENT(IN) :: at(3,3)        !> lattice parameters in alat units 
  REAL(DP),         INTENT(IN) :: force(3,nat)   !> list of atomic forces
  REAL(DP),         INTENT(IN) :: fscale         !> factor for scaling the force 
  CHARACTER(LEN=3), INTENT(IN) :: form           !> format of the structure file (default xsf)
  CHARACTER(LEN=255), INTENT(IN) :: fname        !> file name                                        
  ! 
  ! -- Local Variables
  INTEGER :: i, j, na, ios, iloc 
  CHARACTER(:), ALLOCATABLE :: output

  ! ... Open the file with the good extention
  output = TRIM(fname)//"."//TRIM(form)
  OPEN ( UNIT = ounit, FILE = output, FORM = 'formatted',  STATUS = 'unknown', IOSTAT = ios )


  ! ... Select the format of the file
  SELECT CASE( form )

    CASE( 'xsf' )
      CALL write_xsf( at, nat, tau, order, atm, ityp, force*fscale, ounit )

    CASE( 'xyz')
      CALL write_xyz( at, nat, tau, order, atm, ityp, force*fscale, ounit )

    CASE DEFAULT
      WRITE (ounit,*) " ** LIB::ARTn::WRITE_STRUC::Specified structure format not supported"

  END SELECT


  ! ... Close the file
  CLOSE (UNIT = ounit , STATUS = 'KEEP')
END SUBROUTINE write_struct


! .......................................................................................... XSF
SUBROUTINE write_xsf( at, nat, tau, order, atm, ityp, force, ounit )

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
  REAL(DP),           INTENT(IN) :: force(3,nat)   !> forces
  ! -- LOCAL VARIABLES
  INTEGER :: na, u0, ios, iloc
  REAL(DP) :: at_angs (3,3)

  at_angs = unconvert_length( at )

  WRITE(ounit,*) 'CRYSTAL'
  WRITE(ounit,*) 'PRIMVEC'
  WRITE(ounit,'(2(3F15.9/),3f15.9)') at_angs
  WRITE(ounit,*) 'PRIMCOORD'
  WRITE(ounit,*) nat, 1

  DO na=1,nat
     ! convert positions are in Engine units length
     ! -> And the force???
     iloc = order(na)
     WRITE(ounit,'(a3,3x,6f15.9)') atm(ityp(iloc)), unconvert_length( tau(:,iloc) ), force(:,iloc)
     !Print*, na, atm(ityp(iloc)), unconvert_length( tau(:,iloc) ), unconvert_force( force(:,iloc) )
     !Print*, na, atm(ityp(iloc)), tau(:,iloc), force(:,iloc)
  ENDDO

END SUBROUTINE write_xsf




! .......................................................................................... XYZ
SUBROUTINE write_xyz( at, nat, tau, order, atm, ityp, f, ounit )

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
  ! -- LOCAL VARIABLES
  INTEGER :: na, u0, ios, iloc

  WRITE(ounit,*) nat
  WRITE(ounit,*) ""

  DO na=1,nat
     ! convert positions are in Engine units length
     ! -> And the force???
     iloc = order(na)
     WRITE( ounit, fmt='(a3,3x,6f15.9)', IOSTAT=ios ) atm(ityp(iloc)), unconvert_length( tau(:,iloc) ), unconvert_force( f(:,iloc) )
  ENDDO

END SUBROUTINE write_xyz











