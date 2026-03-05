!---- SIT ----------------------------------------------------------------
MODULE SIT
! This module contains the parameters for the specific ion interaction
! equations, also known as the specific ion-interaction theory (SIT).
! Default values are assigned, and data are read from a text file.
!
! The parameters (epsilon, or "eps") for two ions "i" and "j" are symmetrical,
! that is, eps(i,j) = eps(j,i).  Therefore, a lower triangular matrix is
! enough to store the eps-values.  In this module the eps-values are stored
! in a one-dimensional array (vector) equivalent to the lower triangular
! matrix, and two procedures (getEps and setEps) are used to manage the
! data vector as if it was a lower triangular matrix.
!
! This version also includes eps-values involving neutral species.
!
! Note that the calculation of activity coefficients, using the
! epsilon values in this module, is performed in subroutine 'factor'
! in 'FACTOR_Module'.
!
USE IO, ONLY : ErrStop
USE CHEM
Implicit NONE
!Character (len = 1), PARAMETER :: nl = new_line('a')
!Integer, PARAMETER :: dp = kind(0.d0) ! double precision
Integer :: epsSize = -1, nIon2 = -1
Real(dp), ALLOCATABLE, DIMENSION(:) :: eps

PRIVATE
PUBLIC :: SIT_MEM_ALLOC, SIT_MEM_FREE, getEps, setEps, isEpsAssigned, printEps, readSITdata, setDefaultSITvalues

SAVE

CONTAINS

!-----------------------------------------------------------------------------
SUBROUTINE SIT_MEM_ALLOC (nIon) ! allocates memory for the arrays
! nIon = number of aqueous species in the chemical system
    Implicit NONE
    Integer, INTENT(IN) :: nIon
    Integer :: status = 0, i

    if(nIon <= 0) then
        PRINT *, "Error in module 'SIT': nIon =",nIon," (must be >0)"
        Call ErrStop
    endif

    nIon2 = nIon+2

    ! the two ions added for electroneutrality: Na+ and Cl-
    z(nIon+1) = +1; z(nIon+2) = -1;

    if (ALLOCATED(eps)) Return

    epsSize = nIon2*(nIon2+1)/2   ! Number of non-zero elements in a lower triangular matrix of size nIon2 by nIon2
    ALLOCATE(eps(epsSize), STAT=status)
    if (status /= 0) then
        PRINT *, "Something went wrong when trying to allocate arrays 'eps' in module 'SIT'"
        Call ErrStop
    end if
    Do i = 1, epsSize
        eps(i) = (1234567.8d0) ! flag to indicate "no value assigned"
    EndDo

END SUBROUTINE SIT_MEM_ALLOC

!-----------------------------------------------------------------------------
SUBROUTINE SIT_MEM_FREE  ! deallocates the array memory
    Implicit NONE
    Integer :: status = 0
    if (ALLOCATED(eps)) then
        DEALLOCATE(eps, STAT=status)
    end if
    if (status /= 0) then
        PRINT *, "Something went wrong when trying to deallocate array 'eps' in module 'SIT'"
        Call ErrStop
    end if
END SUBROUTINE SIT_MEM_FREE

!----------------------------------------------------------------------------
SUBROUTINE printEps (UNIT)
! Prints information (to file UNIT) on how activity coefficients are calculated
! if UNIT > 0 printout to file only
! if UNIT = 0 printout to terminal only
! if UNIT < 0, then printout is to both the terminal and to abs(UNIT)
Implicit NONE
Integer, INTENT(IN)  :: UNIT
Integer :: OUT
Integer :: i,j, zi, zj
Real(dp) :: w
Character (len=9) :: str
Logical notPossible, assigned

OUT = abs(UNIT)
if(OUT == 0) OUT = 6

if(.not. ALLOCATED(eps) .or. nIon2 <=0) then
    PRINT *, "Error in 'SIT.printEps': array ""eps"" is not allocated."
    Call ErrStop
end if

write(OUT,'("List of ""epsilon"" values:")')
! --- print a title line with names
write(OUT,'(10x,100i10)') (j,j=1,nIon2)
write(OUT,'(10x)',advance='no')
Do i = 1,nIon
    str = ident(i); str = adjustR(str); write(OUT,'(1x,a9)',advance='no') str
EndDo
str = "Na+"; str = adjustR(str); write(OUT,'(1x,a9)',advance='no') str
str = "Cl-"; str = adjustR(str); write(OUT,'(1x,a9)') str

! --- table body:
DO i = 1,nIon2
    !write(*,'(i10)',advance='no') i
    if(i<=nIon) then
        str = ident(i); str = adjustR(str); write(OUT,'(i2,1x,a9)',advance='no') i,str
    else if(i==nIon+1) then
        str = "Na+"; str = adjustR(str); write(OUT,'(i2,1x,a9)',advance='no') i,str
    else if(i==nIon2) then
        str = "Cl-"; str = adjustR(str); write(OUT,'(i2,1x,a9)',advance='no') i,str
    endif
    if(i <= nIon)   zi = z(i)
    if(i == nIon+1) zi = +1
    if(i == nIon+2) zi = -1
    DO j = 1,nIon2
        IF(j <= i) THEN       ! Non-zero indices of the matrix
            if(j <= nIon)   zj = z(j)
            if(j == nIon+1) zj = +1
            if(j == nIon+2) zj = -1
            notPossible = .false.
            if(i <= nIon) then
                if(i == jWater .or. NOLL(i)) notPossible = .true.
            endif
            if(j <= nIon) then
                if(j == jWater .or. NOLL(j)) notPossible = .true.
            endif
            if((zi*zj) > 0) notPossible = .true.
            if(notPossible) then
                write(OUT,'(a10)',advance="no") "--  "
            else
                assigned = isEpsAssigned(i,j)
                if(assigned) then
                    w = getEps(i,j)
                    write(OUT,'(f10.4)',advance="no") w
                else
                    write(OUT,'(a10)',advance="no") "*   "
                endif
            endif
        ELSE ! if(j > i)
            write(OUT,'(a10)',advance="no") ':   '
        END IF
    END DO
    write(OUT,*)
END DO
Write(OUT,'("Note that eps[i,j]=eps[j,i], and esp[i,j]=0 if z[i]*z[j]>0",/, &
            "(equal charge sign), and  eps[i,i] is used only if z[i]=0.",/, &
            "Unassigned eps[i,j] values are denoted with a ""*"".",/, &
            "The last two rows/columns, Na+ and Cl-, are used only",/, &
            "if the aqueous solution is not electrically neutral.")')

RETURN
END SUBROUTINE printEps

!----------------------------------------------------------------------------
REAL(dp) FUNCTION getEps(i,j) ! i = row; j = column
USE IO
Implicit NONE
Integer, Intent(IN) :: i,j
Integer :: m,jc, row,col

!N = nint((-1.d0+sqrt(1.d0+8.d0*epsSize))/2.d0) ! size of the equivalent square array (in this case N = nIon2)
!epsSize = N*(N+1)/2   = Number of non-zero elements in a lower triangular matrix of size N by N
!epsSize = real(size(eps), dp)
!write(*,'("epsSize=",i3," N=",i3," nIon2=",i3)') epsSize,N,nIon2

if (.not. ALLOCATED(eps) .or. nIon2 <=0) then
    PRINT *, "Error in 'SIT.getEps': array eps is not allocated."
    Call ErrStop
end if
if(i<0 .or. i>nIon2 .or. j<0 .or. j>nIon2) then
    PRINT *, "Error in 'SIT.getEps': i=",i,"  j=,",j,nl,"both must be >0 and <=",nIon2
    Call ErrStop
endif

if(j>i) then
    row = j
    col = i
else
    row = i
    col = j
endif

m = 1
DO jc = 1, col-1
    m = m + nIon2 - jc
EndDo

getEps = eps(m+row-1)
if(getEps >= 1234567.d0) getEps = 0.d0 ! if not assigned, return zero

RETURN
END FUNCTION getEps

!----------------------------------------------------------------------------
Logical FUNCTION isEpsAssigned(i,j) ! i = row; j = column
! Find out if a value for eps(i,j) has been previously assigned
USE IO
Implicit NONE
Integer, Intent(IN) :: i,j
Integer :: m,jc, row,col
Real(dp):: w

!N = nint((-1.d0+sqrt(1.d0+8.d0*epsSize))/2.d0) ! size of the equivalent square array (in this case N = nIon2)
!epsSize = N*(N+1)/2   = Number of non-zero elements in a lower triangular matrix of size N by N
!epsSize = real(size(eps), dp)
!write(*,'("epsSize=",i3," N=",i3," nIon2=",i3)') epsSize,N,nIon2

if (.not. ALLOCATED(eps) .or. nIon2 <=0) then
    PRINT *, "Error in 'SIT.isEpsAssigned': array eps is not allocated."
    Call ErrStop
end if
if(i<0 .or. i>nIon2 .or. j<0 .or. j>nIon2) then
    PRINT *, "Error in 'SIT.isEpsAssigned': i=",i,"  j=,",j,nl,"both must be >0 and <=",nIon2
    Call ErrStop
endif

if(j>i) then
    row = j
    col = i
else
    row = i
    col = j
endif

m = 1
DO jc = 1, col-1
    m = m + nIon2 - jc
EndDo

isEpsAssigned = .true.
w = eps(m+row-1)
if(w >= 1234567.d0) isEpsAssigned = .false.

RETURN
END FUNCTION isEpsAssigned

!----------------------------------------------------------------------------
SUBROUTINE setEps(i,j,w) ! i = row; j = column
USE IO
Implicit NONE
Integer, Intent(IN) :: i,j
Real(dp), Intent(IN) :: w
Integer :: m,jc, row,col

!N = nint((-1.d0+sqrt(1.d0+8.d0*epsSize))/2.d0) ! size of the equivalent square array (in this case N = nIon2)
!epsSize = N*(N+1)/2   = Number of non-zero elements in a lower triangular matrix of size N by N
!epsSize = real(size(eps), dp)
!write(*,'("epsSize=",i3," N=",i3," nIon2=",i3)') epsSize,N,nIon2

if (.not. ALLOCATED(eps) .or. nIon2 <=0) then
    PRINT *, "Error in 'SIT.setEps': array eps is not allocated."
    Call ErrStop
end if
if(i<0 .or. i>nIon2 .or. j<0 .or. j>nIon2) then
    PRINT *, "Error in 'SIT.setEps': i=",i,"  j=,",j,nl,"both must be >0 and <=",nIon2
    Call ErrStop
endif

if(j>i) then
    row = j
    col = i
else
    row = i
    col = j
endif

m = 1
DO jc = 1, col-1
    m = m + nIon2 - jc
EndDo
eps(m+row-1) = w
RETURN
END SUBROUTINE setEps

!---- SIT ----------------------------------------------------------------
SUBROUTINE readSITdata (SITpath, IUT)
! Reads parameters for the SIT activity coefficients model
! from file SITfile='SIT-coefficients.dta' in folder SITpath.
! Errors are printed to the terminal.
! Messages are printed to IUT if IUT > 0.
USE READIR
USE IO, ONLY : UpCase, Max_Path
USE CHEM, ONLY : bareNameOf, chargeOf
Implicit NONE
!-------------------
Character (LEN=*), INTENT(IN)  :: SITpath
Integer, INTENT(IN)  :: IUT
!------------------- ! MXID=30 ! max length for names
Character (LEN=MXID) :: cation, anion, identNa="Na", identCl="Cl", &
                        bareCation,bareAnion,identI,identI2, nameU, &
                        first, second
Character (LEN=Max_Path) :: SITfile="SIT-coefficients.dta"
Logical :: ok, firstLine, setDefaultSITepsValues, foundAll
!Logical :: dbg
Integer :: i,i2,ierr, inFile, ZCat, Zan, Z1, Z2
Real(dp) :: w0

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! get the memory needed for the SIT coefficients
call SIT_MEM_ALLOC(nIon)

! --- get the name of the file containing the SIT coefficients
i=Len_Trim(SITpath)
If(i > 0) Then
    If(SITpath(i:i) /= "\") Then   !add "\" at the end of SITpath
        if(i >= len(SITpath)) then ! no space left to add a "\"
            write(*,'("Error in ''readSITdata'': SITpath = """,A,"""",/,"   (must end with ''\'')")') trim(SITpath)
            Call ErrStop
        endif
        if((Len_Trim(SITfile)+i+1) > len(SITfile)) then
            write(*,'("Error in ''readSITdata'': SITpath is too long.")')
            Call ErrStop
        endif
        SITfile = SITpath(1:i) //"\"// Trim(SITfile)
    Else ! SITpath ends with "\"
        SITfile = SITpath(1:i) // Trim(SITfile)
    EndIf
EndIf

! --- check that everything seems ok
if (.not. ALLOCATED(Z)) then
    write(*,'("Error  in ''readSITdata'': array ""z"" has not been allocated yet")')
    Call ErrStop
endif
if(Z(1) > 100) then
    write(*,'("Error  in ''readSITdata'': array ""z"" does not contain data")')
    Call ErrStop
endif
if (.not. ALLOCATED(IDENT)) then
    write(*,'("Error  in ''readSITdata'': array IDENT has not been allocated yet")')
    Call ErrStop
endif
ok = .false.
Do i = 1,nIon
    if(len_trim(IDENT(i))<=0) Cycle
    ok = .true.
    Exit
EndDo
if(.not.ok) then
    write(*,'("Error  in ''readSITdata'': array IDENT is empty",/,"it should contain the names of chemical species.")')
    Call ErrStop
endif

! --- Open the SIT file
WRITE(*,7002) Trim(SITfile)
WRITE(IUT,7002) Trim(SITfile)
7002 FORMAT("Reading SIT data from file """,A,"""")

inFile = INFL
if(.not.firstInput) inFile = INFL2
OPEN (UNIT=inFile, FILE=SITfile, STATUS='OLD',IOSTAT=IERR)
if(IERR /= 0) then
    Write(*,'("Error  in ''readSITdata'':",/,"can not open file """,A,"""")') Trim(SITfile)
    Call ErrStop
endif

!7017 FORMAT(5x,'eps(',A,',',A,')=',F8.4,:,' +',F8.4,' log(I)')
7017 FORMAT(3x,'eps(',A,',',A,')=',F8.4,:,' ',a)

setDefaultSITepsValues = .true.; firstLine = .true.

! --- Search the SIT-file

! --- Read the specific interaction parameters (epsilon)
! --- Search the SIT-file (compare names such as "Fe+2" and "Fe 2+")
DO WHILE(.true.)

    ! -- get the first species
    nowReading = "First name"
    call READA(first)
    nameU = first; call UPCASE (nameU)
    if(firstLine .and. trim(nameU) == "NODEFAULTS") then
        Write(*,7018) Trim(SITfile)
        7018 FORMAT('"NoDefaults" in file: ',a)
        setDefaultSITepsValues = .false.
        CYCLE ! DO While(.true.)
    endif
    firstLine = .false.
    if(nameU == 'END') EXIT ! DO While(.true.)

    ! -- get the second species
    nowReading = "Second name"
    CALL READA(second)

    ! Here we assume that the file contains: First_species, Second_species, epsilon_value
    ! The first and second species can be either a cation, an anion, or a neutral species.
    ! For purposes of searching the chemical system we assign the cation/neutral as the
    ! "cation" species, and the anion/neutral as the "anion" species.
    ! Note that both the "cation" and the "anion" can be neutral species
    Z1 = chargeOf(first)
    Z2 = chargeOf(second)
    ! Possible cases:
    !  z1  z2  cation
    !  +   -    <- (first)
    !  +   0    <-
    !  0   -    <-
    !  0   +    ->
    !  -   +    ->
    !  -   0    ->
    cation = first;  anion = second  ! default
    ZCat = Z1;  ZAn = Z2
    if(Z1 >= 0) then ! first is cation or neutral
        if(Z2 <= 0) then
            cation = first; ZCat = Z1;   anion = second; ZAn = Z2
        else ! Z2 > 0 (Z1 should be =0)
            cation = second; ZCat = Z2;   anion = first; ZAn = Z1
        endif
    else ! Z1 < 0
        anion = first; ZAn = Z1;   cation = second; ZCat = Z2
    endif

    write(nowReading,'("SIT eps(",a,",",a,")")') trim(first),trim(second)
    CALL READR(w0)
    if(abs(w0) > 1000.d0) then
        Write(*,'("Error  in ''readSITdata'': eps(",a,",",a,")=",f8.4,"  (must be >-1000 and < +1000)")') &
                    trim(first),trim(second),w0
        call ErrStop
    endif

    if((Z1*Z2) > 0) then
        Write(*,'("Error  in ''readSITdata'': eps(",a,",",a,")=",f8.4,"  both species have the same charge sign.")') &
                    trim(first),trim(second),w0
        call ErrStop
    endif
    if(ZCat < 0 .OR. ZAn > 0) then
        Write(*,'("Error  in ''readSITdata'': eps(",a,",",a,")=",f8.4,/,' // &
                 '"  can not figure out which species is the cation and which is the anion.")') &
                    trim(first),trim(second),w0
        call ErrStop
    endif

    ! Compare names such as "Fe+2" and "Fe 2+" when searching
    bareCation = bareNameOf(cation)
    bareAnion = bareNameOf(anion);

    !dbg = .true.
    !if(dbg) write(IUT,7017)  Trim(cation),Trim(anion),w0,"<--"

    foundAll = .false.
    IF(Zcat > 0) THEN
      !-- is this pair of cation-anion present in the chemical system?
      DO i=1,nIon
        if(NOLL(i) .or. Z(i) == 0) Cycle
        identI = bareNameOf(IDENT(i))
        !if(dbg) write(IUT,'("i=",i3," ",a," z=",i3," bareCation=",a)')  i,trim(IDENT(i)),Z(i),bareCation
        !-- is this species "i" in the chemical system = the "anion" in the SIT-file?
        If(Z(i) < 0) Then
          if(identI /= bareAnion .or. Z(i) /= Zan) Cycle ! this "i" is not the "anion", get next "i"
          ! "i" is the "anion", check if the "cation" is Na+
          if(bareCation /= identNa .or. Zcat /= +1) Cycle
          ! "i" is the "anion" and the "cation" is Na+
          !if(dbg) write(IUT,7017)  Trim(cation),Trim(anion),w0,", z(i) <0 and cation = Na+"
          call setEps((nIon+1),i, w0)
          Cycle ! done with this "i"
        EndIf ! (Z(i) < 0)

        !-- This "i" has Z(i) >= 0: is this "i" the "cation"?
        if(bareCation /= identI .or. Zcat /= Z(i)) Cycle ! this "i" is not the "cation"
        ! yes, it is the "cation", is Cl- the "anion"?
        ! if the "anion" is not Cl-: find out if the "anion" is also in the chemical system.

        foundAll = .false.
        ! "i" is the cation, look in the chemical system for all anions (iterate "i2")
        do i2=1,nIon
          if(NOLL(i2).or.Z(i2) > 0) cycle ! look only for "i2"=anion ("i" is cation)
          identI2 = bareNameOf(IDENT(i2));
          if(identI2 /= bareAnion .or. Z(i2) /= Zan) cycle ! "i2" is not the anion
          ! "i" is the cation and "i2" is the anion, get value and quit search
          call setEps(i,i2,w0)
          write(IUT,7017)  Trim(cation),Trim(anion),w0
          ! if the cation is Na+: set eps for the electroneutrality cation
          if(bareCation == identNa .and. Zcat == 1) then
                !if(dbg) write(IUT,7017)  Trim(cation),Trim(anion),w0,", z(i2) <0 and cation = Na+"
                call setEps((nIon+1),i2,w0)
          endif
          foundAll = .true.
          exit ! i2
        enddo !i2
        If(bareAnion == identCl .and. Zan == -1) then
            !if(dbg) write(IUT,7017)  Trim(cation),Trim(anion),w0,", z(i) >0 and anion = Cl-"
            call setEps(i,nIon2,w0) ! "i" is the cation and the anion is Cl-
        endif
        if(foundAll) Exit ! "i"
      END DO ! "i"
      ! look for epsilon(Na+,Cl-)
      if(bareCation == identNa .and. zCat == 1 .and. &
         bareAnion == identCl .and.  zAn == -1) then
            !if(dbg) write(IUT,7017)  Trim(cation),Trim(anion),w0,", cation = Na+, anion = Cl-"
            call setEps((nIon+1), nIon2, w0);
      endif
    ENDIF ! (Zcat > 0)

    ! Is this pair of cation-anion
    ! in fact a cation-neutral or neutral-anion or neutral-neutral couple
    ! and present in the chemical system?
    IF(.not. foundAll .and. (Zcat == 0 .or. Zan == 0)) THEN
      !if(dbg) write(IUT,'("bareCation=",a," bareAnion=",a," check neutral")')  trim(bareCation),trim(bareAnion)
      DO i=1,nIon
        if(NOLL(i) .or. Z(i) /= 0) Cycle
        identI = bareNameOf(IDENT(i))
        !if(dbg) write(IUT,'("i=",i3," ",a," z=",i3," bareCation=",a," check neutral")')  i,trim(IDENT(i)),Z(i),trim(bareCation)
        if(identI /= bareCation .or. Z(i) /= Zcat) cycle ! this "i" is not the cation, get next "i"
        ! yes, it is the "cation", is Cl- the "anion"?
        ! if the "anion" is not Cl-: find out if the "anion" is also in the chemical system.

        foundAll = .false.
        ! "i" is the cation, look in the chemical system for all anions (iterate "i2")
        do i2=1,nIon
          if(NOLL(i2).or.Z(i2) > 0) cycle ! look only for "i2"=anion ("i" is a cation)
          identI2 = bareNameOf(IDENT(i2));
          if(identI2 /= bareAnion .or. Z(i2) /= Zan) cycle ! "i2" is not the anion
          ! "i" is the cation and "i2" is the anion, get value and quit search
          call setEps(i,i2,w0)
          write(IUT,7017)  Trim(cation),Trim(anion),w0
          ! if the cation is Na+: set eps for the electroneutrality cation
          if(bareCation == identNa .and. Zcat == 1) then
                !if(dbg) write(IUT,7017)  Trim(cation),Trim(anion),w0,", z(i2) <0 and cation = Na+"
                call setEps((nIon+1),i2,w0)
          endif
          foundAll = .true.
          exit ! i2
        enddo !i2
        If(bareAnion == identCl .and. Zan == -1) then
          !if(dbg) write(IUT,7017)  Trim(cation),Trim(anion),w0,", z(i) >0 and anion = Cl-"
          call setEps(i,nIon2,w0) ! "i" is the cation and the anion is Cl-
        endif
        if(foundAll) Exit ! "i"
      END DO ! "i"
    ENDIF ! (.not. foundAll .and. (Zcat == 0 .or. Zan == 0))

END DO ! WHILE(.true.)

CLOSE (UNIT=inFile)
Write(IUT,'("Finished reading SIT data.")')

! --- Set Default values for the specific ion interaction coefficients
If(setDefaultSITepsValues) call setDefaultSITvalues (IUT)

RETURN

END SUBROUTINE readSITdata

!-----------------------------------------------------------------------------
SUBROUTINE setDefaultSITvalues (IUT)
!Sets default values for unassigned specific ion interaction coefficients
! if IUT > 0 printout to to both the terminal and to IUT
! if IUT <= 0 printout to terminal only
USE CHEM, ONLY : bareNameOf, MXID, nIon, IDENT, Z
Implicit NONE
Integer, INTENT(IN) :: IUT
Character (len = 4) :: identClO4 = "ClO4"
Character (len = 2) :: identCl = "Cl",identNa = "Na"
Character (LEN=MXID) :: species1, species2
Real(dp) :: w, zj, zi
Integer :: i,j,zc,za
Logical :: ClO41, Cl1, Na1, ClO42, Cl2, Na2, fnd, assigned


if(IUT>0 .and. IUT /= 6) Write(IUT,100)
Write(*,100)
100 Format("Setting ""empty"" SIT-coefficients to default values (Hummel 2009):",/5x, &
                     "eps(M^z+,ClO4-) = 0.2*z",/5X, &
                     "eps(M^z+,Cl-) = 0.1*z - 0.05",/5X, &
                     "eps(Na+,X^y-) = 0.05*y",/2X, &
                     "and the remaining unassigned values are set to:",/5X, &
                     "eps(M^z+, X^y-) = 0.05*y + 0.15*(z-1)",/5X, &
                     "eps(i,j)= 0  (if: Zi*Zj >= zero)",/5x, &
                     "eps(N^0,M/X) = eps(N^0,N^0) = 0.05  (N = neutral molecule, M/X= any ion)")

!Set the default values according to Hummel (2009) for interactions with ClO4, Cl- and Na+
!For eps(Na+,Cl-) and eps(Na+,ClO4-) the default values for Cl- and ClO4-
! are used (instead of the default for Na+)
Do i = 2, nIon+2
    ClO41 = .false.; Cl1 = .false.; Na1 = .false.;
    if(i<=nIon) then
        species1 = bareNameOf(IDENT(i));
        if(trim(species1) == identClO4 .and. z(i) == -1) then
            ClO41 = .true.;
        else if(trim(species1) == identCl .and. z(i) == -1) then
            Cl1 = .true.;
        else if(trim(species1) == identNa .and. z(i) == +1) then
            Na1 = .true.;
        endif
    else if(i==nIon+1) then    ! species1 = "Na";
        Na1 = .true.;
    else if(i==nIon+2) then  ! species1 = "Cl";
        Cl1 = .true.;
    endif
    Do j = 1, (i-1)
        if(z(i)*z(j) >=0) Cycle
        ClO42 = .false.; Cl2 = .false.; Na2 = .false.;
        if(j <= nIon) then
            species2 = bareNameOf(IDENT(j));
            if(trim(species2) == identClO4 .and. z(j) == -1) then
                ClO42 = .true.;
            else if(trim(species2) == identCl .and. z(j) == -1) then
                Cl2 = .true.;
            else if(trim(species2) == identNa .and. z(j) == +1) then
                Na2 = .true.;
            endif
        else if(j==nIon+1) then  ! species2 = "Na";
              Na2 = .true.;
        endif
        ! note that if ClO41 is true, then ClO42 can not be true, etc
        fnd = .false.
        zj = real(z(j),dp); zi = real(z(i),dp)
        if(ClO41) then
            w = (0.2d0 * zj); fnd = .true.
        else if(ClO42) then
            w = (0.2d0 * zi); fnd = .true.
        else if(Cl1) then
            w = (0.1d0 * zj)-0.05d0; fnd = .true.
        else if(Cl2) then
            w = (0.1d0 * zi)-0.05d0; fnd = .true.
        else if(Na1) then
            w = (0.05d0 * zj); fnd = .true.
        else if(Na2) then
            w = (0.05d0 * zi)
        endif
        assigned = isEpsAssigned(i,j)
        if(fnd .and. .not.assigned) call setEps(i,j,w)
    EndDo ! j
EndDo ! i

! other epsilon defaults
Do i = 2, nIon+2
    if(NOLL(i) .or. z(i)==0) Cycle 
    Do j=1, (i-1)
        if(z(i)*z(j) >=0) Cycle
        ! old default values: eps(M^+z, X^-y) = 0.15 + 0.15 (z + y)
        ! new default values: eps(M^z+, X^y-) = 0.05 y + 0.15 (z-1)
        if(z(i)>0) then
            zc = z(i); za = z(j);
        else
            zc = z(j); za = z(i);
        endif
        w = 0.05d0*real(za,dp) + 0.15d0*real((zc-1),dp)
        assigned = isEpsAssigned(i,j)
        if(.not.assigned) call setEps(i,j,w)
    EndDo ! j
EndDo ! i

! neutral species
w = 0.05d0
Do i=1,nIon
    if(NOLL(i).or.Z(i) /= 0) Cycle !i
    assigned = isEpsAssigned(i,i)
    if(.not.assigned) call setEps(i,i,w)
    do j=1,nIon
        if(.not.NOLL(j) .and. Z(j) /= 0) then
            assigned = isEpsAssigned(i,j)
            if(.not.assigned) call setEps(i,j,w)
        endif
    enddo
    call setEps(i,nIon+1,w)
    call setEps(i,nIon2,w)
EndDo
Write(IUT,'("SIT-coefficients assigned.")')

RETURN
END SUBROUTINE setDefaultSITvalues

!-----------------------------------------------------------------------------
!SUBROUTINE NoAq(FNAME,NAME)
!! Return in NAME the text given in FNAME
!! without any ending "(aq)"
!USE IO, ONLY : UpCase
!Implicit NONE
!Character (LEN=*), INTENT(IN) :: FNAME
!Character (LEN=*), INTENT(OUT) :: NAME
!Character (LEN=4) :: CDum
!Integer :: NChr
!Nchr=Len_Trim(FNAME)
!NAME=FNAME
!if(NChr<=4) Return
!   CDum=FNAME(Nchr-3:Nchr)
!   Call UpCase(CDum)
!   if(Cdum=='(AQ)') NAME=FNAME(1:Nchr-4)
!Return
!End SUBROUTINE NoAq

END MODULE SIT