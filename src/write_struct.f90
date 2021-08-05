!SUBROUTINE write_struct(alat, at, nat, tau, atm, ityp, force, fscale, ounit, form, fname)
SUBROUTINE write_struct( at, nat, tau, order, atm, ityp, force, fscale, ounit, form, fname)
  !
  ! A subroutine that writes the structure to a file (based on xsf_struct of QE)  
  !
  USE artn_params, ONLY: DP,B2A 
  USE UNITS
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: nat             ! number of atoms 
  INTEGER, INTENT(IN) :: ityp(nat)       ! atom type
  INTEGER, INTENT(IN) :: order(nat)       ! atom type
  CHARACTER(LEN=3), INTENT(IN) :: atm(*) ! contains information on atomic types 
  INTEGER, INTENT(IN) :: ounit           ! output fortran unit 
  !REAL(DP), INTENT(IN) :: alat           ! alat of QE   
  REAL(DP), INTENT(IN) :: tau(3,nat)     ! atomic positions 
  REAL(DP), INTENT(IN) :: at(3,3)        ! lattice parameters in alat units 
  REAL(DP), INTENT(IN) :: force(3,nat)   ! forces
  REAL(DP), INTENT(IN) :: fscale         ! factor for scaling the force 
  CHARACTER(LEN=3), INTENT(IN) :: form   ! format of the structure file (default xsf)
  CHARACTER(LEN=255), INTENT(IN) :: fname  ! file name 
  ! 
  !
  INTEGER :: i, j, na, ios, iloc 
  REAL(DP) :: at_angs (3,3)
  OPEN ( UNIT = ounit, FILE = fname, FORM = 'formatted',  STATUS = 'unknown', IOSTAT = ios )
  IF ( form == 'xsf' ) THEN 
     
     !DO i=1,3
     !   DO j=1,3
     !      !at_angs(j,i) = at(j,i)*B2A
     !      at_angs(j,i) = unconvert_length( at(j,i) )
     !   ENDDO
     !ENDDO
     at_angs = unconvert_length( at )

     WRITE(ounit,*) 'CRYSTAL'
     WRITE(ounit,*) 'PRIMVEC'
     WRITE(ounit,'(2(3F15.9/),3f15.9)') at_angs
     WRITE(ounit,*) 'PRIMCOORD'
     WRITE(ounit,*) nat, 1
     
     DO na=1,nat
        ! positions are in Angstroms
        ! -> And the force???
        iloc = order(na)
        !WRITE(ounit,'(a3,3x,6f15.9)') atm(ityp(iloc)), tau(:,iloc)*B2A, force(:,iloc)*fscale
        WRITE(ounit,'(a3,3x,6f15.9)') atm(ityp(iloc)), unconvert_length( tau(:,iloc) ), force(:,iloc)*fscale

     ENDDO
  ELSE
     WRITE (ounit,*) "Specified structure format not supported"
  ENDIF
  CLOSE (UNIT = ounit , STATUS = 'KEEP')
END SUBROUTINE write_struct
