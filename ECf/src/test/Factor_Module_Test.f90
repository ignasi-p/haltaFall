MODULE FACTOR_Module
USE CHEM, ONLY : z
Integer, PARAMETER :: dp = kind(0.d0) ! double precision
! The model to calculate activity coefficients (ionic strength effects):
!   < 0 for ideal solutions (all activity coefficients = 1)
!   = 0 Davies eqn.
!   = 1 SIT (Specific Ion interaction "Theory")
!   = 2 Simplified HKF (Helson, Kirkham and Flowers)
Integer :: activityCoeffsModel=-1
! The ionic strength, either provided by the user
! or calculated (if the supplied ionic strength is negative).
Real(dp) :: ionicStr = -1.D0
! Sum of molalities of all species in the aqueous solution. This is
! only calculated if the ionic strength (ionicStrength) is negative
! (the ionic strength has to be calculated) and if the activity
! coefficients are also calculated.
Real(dp) :: sumM = 0.D0
! the calculated osmotic coefficient of water
Real(dp) :: phi = 1.D0
Real(dp) :: electricBalance = 0.d0

Real(dp) :: Temperature=-1000.D0, Pressure=-1.D0

! Used internally: true if ionicStr is negative at the start
! of the calculation; false if ionicStr is zero or positive at
! the start of the calculation
Logical, PRIVATE :: calculateIonicStr = .true.
!Real(dp), PRIVATE :: NaN = -1._dp

!PRIVATE :: ErrStop

SAVE

CONTAINS

! !-----------------------------------------------------------------------------
  ! SUBROUTINE FACTOR_MEM_ALLOC (nIon) ! allocates memory for the arrays
    ! USE CHEM
    ! Implicit NONE
    ! Integer, INTENT(IN) :: nIon
    ! Integer :: nIon2, status = 0
    ! NaN = sqrt(m1)
! !    if(isNaN(NaN)) print *, "NaN is a NaN"
! !  nIon = number of soluble species
    ! if (nIon <= 0) then
        ! PRINT *, "Error in 'FACTOR_Module': nIon =",nIon," (must be >0)"
        ! Call ErrStop
    ! endif
    ! nIon2 = nIon+2
    ! if (.not. ALLOCATED(Z)) then
        ! ALLOCATE(Z(nIon), STAT=status)
    ! end if
    ! if (status /= 0) then
        ! PRINT *, "Something went wrong when trying to allocate arrays in 'FACTOR_Module'"
        ! Call ErrStop
    ! end if
    ! Z(1) = huge(0)
    ! activityCoeffsModel=-1
    ! phi= 1.D0
    ! sumM = 0.D0
  ! END SUBROUTINE FACTOR_MEM_ALLOC

! !-----------------------------------------------------------------------------
  ! SUBROUTINE FACTOR_MEM_FREE  ! deallocates the array memory
    ! Implicit NONE
    ! Integer :: status = 0
    ! if (ALLOCATED(Z)) then
        ! DEALLOCATE(Z, STAT=status)
    ! end if
    ! if (status /= 0) then
        ! PRINT *, "Something went wrong when trying to deallocate arrays in 'FACTOR_Module'"
        ! Call ErrStop
    ! end if
  ! END SUBROUTINE FACTOR_MEM_FREE

! !----------------------------------------------------------------------------
! SUBROUTINE ErrStop ! (PRIVATE)
    ! Implicit NONE
    ! Integer :: ierr, i
    ! WRITE (*, '("Press Enter to continue ...")', Advance='No')
    ! READ (*,'(I5)',iostat=ierr) i
    ! STOP 1
! END SUBROUTINE ErrStop

! ----------------------------------------------------------------------------
SUBROUTINE FactorStart (NIons,Conc,lnF)
USE CHEM, ONLY : DBG, IOUT
Implicit NONE
Integer, INTENT(IN) :: NIons
Real(Kind(1.D0)), INTENT (IN) :: Conc(NIons)
Real(Kind(1.D0)), INTENT (OUT) :: lnF(NIons)
Integer :: I, IUT
IUT = IOUT;  if(IUT <= 0) IUT = 6
If(DBG >= 3) Write(IUT,'("factorStart(",I3,"), activityCoeffsModel =",I3)') nIons,activityCoeffsModel
SumM = 0.D0; electricBalance = 0.D0
If(activityCoeffsModel >= 0) Then
    if(ionicStr > 0.D0) then
        ionicStr = min(ionicStr,200.D0)
    else
        ionicStr = -1.D0
    endif
EndIf
calculateIonicStr = ionicStr < 0.D0
If(abs(ionicStr) < 1.D-35 .or. activityCoeffsModel < 0) Then
    DO I=1,NIons
        lnF(I)=0.D0
    EndDo
EndIf
If(DBG >= 3) Write(IUT,'("factorStart() ut, calculateIonicStr = ",L1)') calculateIonicStr
END SUBROUTINE FactorStart

! ----------------------------------------------------------------------------
SUBROUTINE FACTOR (NIons,Conc,lnF)
! Must be provided to calculate activity factors for aqueous species.
USE CHEM
Implicit NONE
Integer, INTENT(IN) :: NIons
Real(Kind(1.D0)), INTENT (IN) :: Conc(NIons)
Real(Kind(1.D0)), INTENT (OUT) :: lnF(NIons)
Real(Kind(1.D0)), PARAMETER :: ln10=2.3025850929940D0
Real(Kind(1.D0)) :: temp, down, phiDH, rootI, SumPrd
!  lnF()= natural log of activity coefficient
Integer :: I
IF (abs(ionicStr) < 1.D-35 .or. activityCoeffsModel /= 0) RETURN
electricBalance=0.D0
IF(calculateIonicStr) THEN
    !Calculate Electrical Balance
    !    if electricBalance > 0,  then  Cl-  are added
    !    if electricBalance < 0,  then  Na+  are added
    !  to maintain an electrically neutral solution
    DO I=1,NIons
        IF (.not.NOLL(I)) electricBalance = electricBalance + Z(I)*Conc(I)
    END DO
    !Calculate Ionic Strength.
    ionicStr=0.D0
    IF(ABS(electricBalance) > 1.D-15) ionicStr = ABS(electricBalance)
    DO I=1,NIons
        IF (.not.NOLL(I)) ionicStr = ionicStr + Z(I)*Z(I)*Conc(I)
    END DO
    ionicStr=0.5D0 * ionicStr
ELSE
    ionicStr = max(ionicStr,0.D0)
ENDIF
rootI=SQRT(ionicStr)
down = 1.0D0 + rootI
temp= - 0.51002D0 * ((rootI/down) - 0.3D0*ionicStr)
DO I=1,NIons
   if(NOLL(I)) then
        lnF(I) = 0.D0
        cycle
   endif
   if(Z(I) /= 0) then
        lnF(I)= ln10 * Z(I)*Z(I)*temp
   else
        lnF(I)= ln10 * 0.1D0*ionicStr
   end if
   IF(lnF(I) > 20.D0) lnF(I)= 20.D0
   IF(lnF(I) < -20.D0) lnF(I)=-20.D0
END DO
!Calculate osmotic coeff. and activity of water
IF(JWATER <= 0) RETURN
phiDH= -0.69685D0 * (down-2.D0*LOG(down)-(1.D0/down))
SumM=0.D0
IF(abs(electricBalance) > 1.D-15) SumM=abs(electricBalance)
SumPrd= 0.D0
DO I=1,NIons
   IF (NOLL(I)) CYCLE
   SumM = SumM + Conc(I)
   IF(Z(I) == 0) CYCLE
   SumPrd=SumPrd + Z(I)*Z(I)*0.15D0 * Conc(I) * ionicStr
END DO
IF(abs(electricBalance) > 1.D-15) SumPrd = SumPrd + 0.15D0 * abs(electricBalance) * ionicStr
! Calculate Osmotic Coeff.
phi= 1.D0
IF (SumM > 1.D-15) phi= 1.D0 + ( phiDH + ln10*SumPrd )/SumM
!Calculate Water Activity.
TEMP = -phi * SumM / 55.5102D0
lnF(JWATER) = max(-3.D0,min(TEMP,0.D0))
RETURN
END SUBROUTINE FACTOR

END MODULE FACTOR_Module
