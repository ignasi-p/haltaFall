!    ================================
!      EC: EQUILIBRIUM CALCULATIONS
!    ================================
! ================================================================
!  Program to perform Chemical Equilibrium Calculations
!  It uses the HALTAFALL algorithm
!  Developed by: Ignasi Puigdomenech
!                The Royal Institute of Technology
!                Dep. Inorganic Chemistry
!                S-100 44 STOCKHOLM
!                SWEDEN
! ================================================================
! Variable-Vector indexes run as follows:
!     KH(i), TOT(i), and
!     IDENTC(i)   for i=1...NA
!     LBETA(i)    for i=1...NX+MSOL
!     Z(i)        for i=1...NA+NX+2
!     A(i,j)      for i=1...NA,  and j=1...NX+MSOL
!     logA(i), C(i), IDENT(i), and
!     NOLL(i)     for i=1...NA+NX+MSOL
!     logF(i)     for i=1...NA+NX
!     BT(i,j)     for i=1...NA,  and j=1...NP (=number of points)
! - - - - - - - - - - -

PROGRAM EC
USE IO, ONLY : Max_Path, Quit, ErrStop, GetExt, TimDat, UpCase, getParent
USE READIR, ONLY : INFL
USE CHEM
USE HALTAF_Module
USE FACTOR_Module, ONLY : Temperature, Pressure, factorPrint, p1Sat
USE SIT
Implicit NONE

!                    **************************
!                              YYYY-MM-DD
Character (LEN=12) :: Version='2026-03-04'
!                    **************************

Character (LEN=Max_Path) :: FILIN, FILUT, FileRES, FILPATH, ACPATH
Integer, PARAMETER :: IUT = 21

Integer ::  J,NP, ierr
Character (LEN=400) :: TEXT
Character (LEN=80) :: dashLINE
Character (LEN=11) :: DateStr, TimeStr
! keep terminal open (pause) after program execution
Logical :: keep = .false.
! default relative tolerance for HaltaFall
Real(dp), PARAMETER :: TOL_HALTA_DEF = 1.D-5
! the default tolerance for activity coefficients
Real(dp), PARAMETER :: TOL_HALTA_DEF_LOGF = 1.D-3
Real(dp), PARAMETER :: ln10 = log(10.d0)
! the tolerances may be superseeded by the command-line arguments
! that are read in routine "FileOPEN"
Real(dp) :: tolHalta0 = TOL_HALTA_DEF
Real(dp) :: tolLogF0 = TOL_HALTA_DEF_LOGF
Logical :: calcFailed, tableFirstTime = .true.
Integer :: NPKT, dbg0
Logical :: debug
Integer, Allocatable :: HOW(:,:)
Real(dp), Allocatable ::  BT(:,:)
Real(dp) :: start, finish
Logical :: POINTS
Integer :: activityCoeffsModel0
Real(dp) :: maxTotConc, tC, pBar
! the separator between data columns in the output table
Character(len=1) :: sep = ","

SAVE

WRITE(*,1) Trim(Version)
1 FORMAT('.',41('-'),'.', &
    /,'| EQUILIBRIUM CALCULATIONS',T43,'|', &
    /,'| Modified HaltaFall, vers.',A,T43,'|', &
    /,'| (by I.Puigdomenech)',T43,'|', &
    /,'|',41('_'),'|')

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- Open input and output files, read command-line arguments
CALL FileOPEN
! at this point the input file (with full path name FILIN)
! is opened on INFL and the output file to IUT
! FileRES is later opened at subroutine "print"

CALL TimDat (DateStr,TimeStr)

WRITE(IUT,2) Trim(DateStr),Trim(Version),Trim(TimeStr)
2 FORMAT('EQUILIBRIUM CALCULATIONS',T62,'Date: ',A,/,  &
'Modified HaltaFall vers.',A,T62,'Time: ',A,/, &
'(by I.Puigdomenech)',/)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- Read Input file - (allocate memory for the arrays)
Call ReadChemSystem (tC, pBar)

IF(tC > -273.) THEN !###
    IF(Temperature <= -273.) Temperature = tC;
    IF(ABS(tC - Temperature) > 0.01) THEN
        WRITE(*,3002) Temperature,tC, tC
        3002 FORMAT('? Warning: temperature mismatch.',/, &
            '  entered in command line as ',F6.2,', but in the input file it is given as ',F6.2,/, &
            '  The value ',F6.2,' will be used.')
        Temperature = tC;
    END IF
END IF
IF(Temperature > 0.) THEN ! ####
    Pressure = max(pBar,p1Sat(Temperature)); ! set the pressure to max(1,pSat,pBar)
ENDIF

! --- print some info ---
WRITE(IUT,1051) Trim(FILPATH),Trim(FILIN), Trim(FILUT), Trim(FileRES)
1051 Format (/,"Input, output and table files:",/4x,"path: """,A,"""",/3x,"input: """,A,"""", &
             /2x,"output: """,A,"""",/3x,"table: """,A,"""")
If(dbg /= 0) Write(IUT,'("Debug printout level: ",i2)') dbg
if((abs(tolHalta0-TOL_HALTA_DEF)/TOL_HALTA_DEF) > 1.d-3) Write(IUT,'("Relative tolerance in HaltaFall = ",1PE9.2)') tolHalta0
If(activityCoeffsModel >= 0 .and. activityCoeffsModel <= 2) then
    if(activityCoeffsModel == 0) then
        TEXT = "Davies eqn."
    else if(activityCoeffsModel == 1) then
        TEXT = "SIT"
    else
        TEXT = "simpl. HKF"
    endif
    Write(IUT,'("Activity Coeffs.:  method =",i3," (",A,")")') activityCoeffsModel, trim(TEXT)
    if(Temperature > 0.) Write(IUT,'(24x,"t = ",F6.2,"""C, p =",F8.2," bar")') Temperature, Pressure
    if(abs(ionicStr) > 1.e-12) Then
        if(ionicStr<0) Then
           Write(IUT,'(11x,"Ionic strength = ",F6.2," (calculated at each point)")') ionicStr
        Else
           Write(IUT,'(11x,"Ionic strength = ",F6.2," mol/kg")') ionicStr
        EndIf
    EndIf
    if((abs(tolLogF0-TOL_HALTA_DEF_LOGF)/TOL_HALTA_DEF_LOGF) > 1.d-3) &
        Write(IUT,'(18x,"tolLogF = ",1PE9.2)') tolLogF0
    if(activityCoeffsModel == 1 .and. len_trim(ACPath) > 0) &
        Write(IUT,'("   path to file ""SIT-coefficients.dta""",/8x,"""",A,"""")') trim(ACPath)
    TEXT = ""
EndIf
WRITE(IUT,*)

Call ReadConcentrations

! --- Read and print out heading (title)
Write(dashLINE,'(a,39(" -"),a)') nl,nl
Write(IUT,'(a)') dashLINE
j = 0
DO WHILE (.true.)
    Read(INFL,'(a)',iostat=ierr) TEXT
    if(ierr /= 0) EXIT
    if(len_trim(TEXT) <= 0 .or. trim(TEXT) <= " ") cycle
    Write(IUT,'(a)') trim(TEXT)
    j = j+1
END DO
CLOSE (UNIT=INFL)

if(j > 0) Write(IUT,'(a)') dashLINE

! --- print the chemical system
!if(dbg > 0) then
    Call printChemSystem (IUT)
    WRITE(IUT,'(a)') dashLINE
!endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! --- Prepare the HaltaFall calculations
!     Activity coefficient calculations:
If(activityCoeffsModel == 1) then !SIT
    ! Allocate memory for epsilon-values
    ! and read them from a file.
    ! Note that this uses the READIR module,
    ! so unit INFL must be free for reading
    call readSITdata (ACPath, IUT) !###
EndIf

! --- (End of input) ---

! --- print act.coeffs. model description,
Call factorPrint (IUT)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!           Begin the calculations
WRITE(*,'("The calculations begin now!")')
call CPU_time(start)

! allocate memory for working arrays

call HALTAF_MEM_ALLOC

CONT=.false.
IOUT = IUT ! debug printout in HaltaFall
! IOUT = 6 ! un-comment if you wish debug print-out to the console
dbg0 = dbg

activityCoeffsModel0 = activityCoeffsModel;

NP=0
DO WHILE (.true.)

    NP = NP + 1
    if(NPKT > 1) write(*,'(A,"starting calculations for point ",I0)', Advance='No') CHAR(13),NP
    WRITE(IUT,'("*** POINT ",i0,"  ***")') NP
    if(POINTS) CONT=.false. ! independent list of points
    DO J=1,NA
        IF(POINTS) KH(J)=HOW(J,NP) ! independent list of points 
        if(KH(J) == 1) then
            TOT(J)=BT(J,NP)
            if(.not.CONT) then
                LOGA(J)=-10.D0
                IF(BT(J,NP) > 0.D0) LOGA(J)= log10(BT(J,NP)) - 3.D0
            endif
        endif
        if(KH(J) == 2) LOGA(J)=BT(J,NP)
    END DO

    maxTotConc = 0.d0;
    do j =1, Na
        if(kh(j) == 1) maxTotConc = max(abs(tot(j)), maxTotConc);
    enddo

    activityCoeffsModel = activityCoeffsModel0

    ! - - - - - - - - - - - - - - - - - - - - - - - -
    !     Do the Chemical Equilibrium Calculation

    TOL = tolHalta0; tolLogF = tolLogF0
    ! ---- debugging ---- ####
    if(NP == 32) then
        !dbg=7
    else
        dbg = dbg0
    endif
    ! ---- ---- ---- ---- ####

    CALL HaltaCalc

    ! check for errors from HaltaFall
    if(errFlags > 0) then
        TEXT = errFlagsGetMessages()
        write(IUT,'(a,/,a)') errFlagsToString(),trim(TEXT)
        !   1 = The numerical solution is uncertain (round-off errors)
        !   2 = Too many iterations when solving the mass balance equations
        !   3 = Failed to find a satisfactory combination of solids
        !   4 = Too many iterations trying to find the solids at equilibrium
        !   5 = Some aqueous concentration(s) too large (>UNREASONABLE_CONC): uncertain activity coefficients
        !   6 = Activity factors did not converge
        !   7 = Calculation interrupted by the user
        calcFailed = (isErrFlagSet(2) .or. isErrFlagSet(3) .or. isErrFlagSet(4) .or. isErrFlagSet(6))
        if(calcFailed) then
            Cont = .false.
            if(isErrFlagSet(6) .and. maxTotConc > 0.3d0 .and. activityCoeffsModel == 1) then
                !if(dbg > 2) &
                write(IUT,'("(performing a calculation with the ""simpl.HKF"" model for act. coeffs.)")')
                activityCoeffsModel = 2;
                CALL HaltaCalc
                write(IUT,'("(performing following calculations with the ""SIT"" model for act. coeffs.)")')
                activityCoeffsModel = activityCoeffsModel0
                CALL HaltaCalc
                calcFailed = (isErrFlagSet(2) .or. isErrFlagSet(3) .or. isErrFlagSet(4) .or. isErrFlagSet(6))
                Cont = .false.
            endif
            do while (.true.)
                if(.not.calcFailed .or. TOL <= 1d-12) exit ! do while (true)
                tol = tol * 0.1D0 ! decrease tolerance and try again
                TEXT = errFlagsGetMessages()
                Write(*,'(/,a,/,"  decreasing tolerance to: ",1PE9.2," and trying again.")') trim(TEXT),TOL
                Write(IUT,'(/,a,/,"  decreasing tolerance to: ",1PE9.2," and trying again.")') trim(TEXT),TOL
                CALL HaltaCalc
                calcFailed = (isErrFlagSet(2) .or. isErrFlagSet(3) .or. isErrFlagSet(4) .or. isErrFlagSet(6))
            enddo
            Cont = .false.
            if(errFlags > 0) then ! success?
                TEXT = errFlagsGetMessages()
                write(*,'(/,a)') trim(TEXT)
                write(IUT,'(/,a,/,a)') errFlagsToString(),trim(TEXT)
                ! ---- debugging ---- ####
                if(NP == 99999) then
                    call factorPrint (IUT)
                    !tol = tolHalta0; tolLogF = tolLogF0
                    !dbg=5
                    !CALL HaltaCalc
                    !dbg = dbg0
                    !call ErrStop
                endif
                ! ---- ---- ---- ---- ####
            endif ! if(errFlags > 0)
            if(NP < NPKT .and. abs((tol-tolHalta0)/tolHalta0) > 0.001) then
                Write(*,'(/,"Restoring tolerance to: ",1PE9.2," for next calculation.")') tolHalta0
                Write(IUT,'(/,"Restoring tolerance to: ",1PE9.2," for next calculation.")') tolHalta0
            endif
            tol = tolHalta0; tolLogF = tolLogF0
        endif ! if(calcFailed)
    endif ! if(errFlags > 0)

    ! - - - - - - - - - - - - - - - - - - - - - - - -
    !            Print out the results

    ! Print results.
    !if(NPKT > 1) write(*,'(A,"end of calculations for point ",I0,"  ")', Advance='No') CHAR(13),NP

    call printConcs (IUT)

    call tablePRINT

    ! WRITE(IUT,'(a)') dashLINE

    ! - - - - - - - - - - - - - - - - - - - - - - - -
    !            Next point
    If(NP >= NPKT) EXIT ! Do While (true)

END DO ! WHILE (true)

if(NPKT > 1) write(*,'(A,"                                     ",A)', Advance='No') CHAR(13),CHAR(13)

call CPU_time(finish)
write(*,'("Calculated ",i0," points, wall-clock time = ",F8.3)') NP,(finish-start)

WRITE(IUT,'("All Done")')
CLOSE (UNIT=IUT)
WRITE(*,'("All Done")')

if(keep) call Quit ! pause

CONTAINS

!-----------------------------------------------------------------------------
SUBROUTINE ReadConcentrations
! Reads information on the concentrations for each chemical component.
!
! The format is the same as for Medusa/Hydra and Spana/Database, except that
! there is no diagram type information to read. Instead, the number of
! calculations is given by the user.
! This routine reads first the number of points (chemical compositions)
! to calculate (NPKT).
!  - If NPKT is negative (POINTS=true), then information for each of
!    the calculation points is specified in succession.
!  - If NPKT is positive, then calculation points between two extreme
!    compositions are evaluated.
USE READIR
Implicit NONE
Real(dp) :: BC, STEP, W
Character (LEN=3) :: CCONC
Integer :: j, n, NP, status
Logical :: pointsOK
! - - - - - - - - - - - - - - - - - - - - - - - -
!   READ THE CONCENTRATIONS OF EACH COMPONENT
!
WRITE (*,'("Reading Input Disk File:   concentration data...  ",A)',Advance='No') CHAR(13)
CALL READI(NPKT)

if(NPKT == 0) then
    WRITE(*,4000)
    WRITE(IUT,4000)
    4000 FORMAT('? Error:  number of points to calculate is zero.')
    Call ErrStop
endif

POINTS=.false.
If(NPKT < 0) Then
    POINTS=.true.
    NPKT=-NPKT
EndIf
pointsOK = (NPKT <= 1)

! allocate memory for the arrays
    if (NA <= 0) then
        Write(*,'("Error in ''ReadConcentrations'': NA =",I0)') NA
        STOP 1
    endif
    status = 0
    if (.not. ALLOCATED(HOW)) then
        ALLOCATE(HOW(NA,NPKT), STAT=status)
        ALLOCATE(BT(NA,NPKT), STAT=status)
    end if
    if (status /= 0) then
        Write(*,'("Something went wrong when trying to allocate arrays in ''ReadConcentrations''")')
        STOP 1
    end if

NP=0
! Input Concentrations
DO WHILE (.true.)
    IF(POINTS) NP = NP+1
    DO J = 1, NA
        CALL READA (CCONC)
        CALL UPCASE (CCONC)
        KH(J)=-1
        IF(CCONC == 'T')    KH(J)=1
        IF(CCONC == 'TV')   KH(J)=2
        IF(CCONC == 'LTV')  KH(J)=3
        IF(CCONC == 'LA')   KH(J)=4
        IF(CCONC == 'LAV')  KH(J)=5
        IF(KH(J) == -1) Then
            WRITE(*,4001) CCONC, IDENTC(J)
            WRITE(IUT,4001)  CCONC, IDENTC(J)
            4001 FORMAT('? Error:  "',A,'"  is NOT any of: T,TV,LTV,TLA,LA',  &
            ' or LAV',/2X,'Input ERROR reading concentration for ',A)
            Call ErrStop
        EndIf
        If(J == jWater .and. KH(J) <= 3) Then
            WRITE(*,4002) J
            WRITE(IUT,4002) J
            4002 FORMAT('? You cannot give a total concentration for water (component',i2,') !',/2X,  &
                          'The calculations are made for 1000 g of water')
            Call ErrStop
        EndIf
        If(KH(J) == 1 .or. KH(J) ==4) Then  ! 'T' or 'LA'
            ! constant value (Total concentration or log(Activity))
            CALL READR (BC)
            If(POINTS) Then
                BT(J,NP)=BC
            Else
                DO N = 1, NPKT
                    BT(J,N) = BC
                END DO
            EndIf
            CYCLE !J
        Else If(KH(J) == 2 .or. KH(J) == 3 .or. KH(J) == 5) Then ! 'TV', 'LTV' or 'LAV'
            ! varied Total concentration, log (Total concentration) or log(Activity)
            pointsOK = .true.
            If(NPKT == 1) Then
                WRITE(*,4003) IDENT(J)
                WRITE(IUT,4003)  IDENT(J)
                4003 FORMAT('? You cannot vary the concentration of ',A,/2x,  &
                'with only one point to calculate ')
                Call ErrStop
            EndIf
            IF(POINTS) Then
                WRITE(*,4004)IDENT(J)
                WRITE(IUT,4004) IDENT(J)
                4004 FORMAT('? You cannot vary the concentration of ',A,/2x,  &
                'if you give negative NPKT')
                Call ErrStop
            EndIf
            CALL READR(BT(J,1))
            CALL READR(W)
            STEP=(W-BT(J,1))/(NPKT-1)
            DO N = 2, NPKT
                BT(J,N) = BT(J,N-1) + STEP
            END DO
            IF(KH(J) /= 3) CYCLE !J
            ! for 'LTV':
            DO N=1,NPKT
                BT(J,N)=EXP(ln10*BT(J,N))
            END DO
        EndIf
    END DO !J

    ! If KH =1  Tot.Conc. is given. The Mass Balance Eq. has to be solved.
    !    KH =2  log(activity) is given. The Tot.Conc. will be calculated.
    DO J=1,NA
        If(KH(J) > 3) Then
            KH(J)=2
            IF(POINTS) HOW(J,NP)=KH(J)
        Else
            KH(J)=1
            IF(POINTS) HOW(J,NP)=KH(J)
        EndIf
    END DO

    If(.not.POINTS .or. NP >= NPKT) EXIT

END DO ! WHILE (.true.)

    if(.not.POINTS .and. .not.pointsOK) then
        WRITE(*,4005)
        WRITE(IUT,4005)
        4005 FORMAT(/,'Error: the number of points is > 0',/2x,  &
        'but no concentration is varied?')
        Call ErrStop
    endif

WRITE (*,'("                                                  ",A)',Advance='No') CHAR(13)
RETURN
END SUBROUTINE ReadConcentrations

!------------------------------------------------------------------------------
SUBROUTINE tablePrint
! This routine is called after each point has been calculated.
! The first time it opens the output table file ("name.csv") and writes a
! heading on it.  Every time a point has been calculated, and the chemical
! equilibrium composition is known, this routine is called.
! The columns of data are separated by the separator "sep".
Implicit NONE
! The chemical composition is printed in a single line, each value in a field
! that is 12 characters wide. If 'truncate' is true, long names are cut-off
! and a "plus" is appended.
Logical :: truncate = .false.
Character(len=1) :: plus = "&"
Integer, PARAMETER :: IPRN=25 ! output unit
Logical :: activityCoeffs = .false.
Integer :: i, j,k,n, Hplus=-1, eMinus = -1
Real(dp) :: pH, pe, Eh, w
Character(len=MXID+5) :: conc_i, solub_i, log_i, press, nmn
Integer :: ioErr

SAVE

If(tableFirstTime) Then
    ! ------------------------------------------------------------------
    !  First time: open the table file (possibly already exists) 
    ! ------------------------------------------------------------------
    OPEN (Unit=IPRN,File=FileRES,Status='UnKnown',encoding='utf-8') ! ,Access='append'

    !  do things that need to be done only once

    ! find out which species are H+ and e-
    do i=1,nIon
      if(ident(i) == 'H+' .or. ident(i) == 'H +') Hplus = i
      if(ident(i) == 'e-' .or. ident(i) == 'e -' .or. ident(i) == 'E-' .or. ident(i) == 'E -') eMinus = i
    end do
    activityCoeffs = (activityCoeffsModel >=0 .and. activityCoeffsModel <= 2)
    ! -------------------------------
    !   write a heading on the file
    ! -------------------------------
    ioErr = 0
    write (IPRN,'(" ")', Advance='No',IOSTAT=ioErr)
    if(ioErr /= 0) then
        write(*,'(/,"ERROR: can not write to file: """,A,"""",/,"   is the file locked?")') trim(FileRES)
        call ErrStop
    endif
    if(activityCoeffs .and. temperature > -273.d0) &
        write (IPRN,'(" Act.Coef.",a," Temp.(''C)",a,"Press.(bar)",a)', Advance='No') sep,sep,sep
    if(Hplus > 0)  write (IPRN,'("     pH    ",a)', Advance='No') sep
    if(eMinus > 0) then
        write (IPRN,'("     pe    ",a)', Advance='No') sep
        if(temperature > -273.d0) write (IPRN,'("    Eh_mV  ",a)', Advance='No') sep
    endif
    if(activityCoeffs) write (IPRN,'(" Ionic.Str.",a," El.Balance",a," Sum_m     ",a)', Advance='No') sep,sep,sep
    if(jWater > 0) then
        write (IPRN,'(" logA(H2O) ",a)', Advance='No') sep
        if(activityCoeffs) write (IPRN,'("    phi    ",a)', Advance='No') sep
    endif
    Do i = 1,NA
        if(i == jWater) Cycle
        nmn = identc(i)
        if(truncate .and. len_trim(nmn) > 6) nmn = identc(i)(1:5) // plus
        Write(conc_i, '("Tot(",a,")")') trim(nmn)
        Write(solub_i,'("Sol(",a,")")') trim(nmn)
        Write(log_i,  '("lgA(",a,")")') trim(nmn)
        j = max(11,len_trim(conc_i))
        k = max(11,len_trim(solub_i))
        n = max(11,len_trim(log_i))
        write (IPRN,'(6a)', Advance='No') conc_i(1:j),sep, solub_i(1:k),sep, log_i(1:n),sep
    EndDo
    Do i = 1,nIon
        if(i == jWater .or. i == eMinus) Cycle
        nmn = ident(i)
        if(truncate .and. len_trim(ident(i)) > 8) nmn = ident(i)(1:7) // plus
        Write(conc_i,'("C(",a,")")') trim(nmn)
        j = max(11,len_trim(conc_i))
        write (IPRN,'(2a)', Advance='No') conc_i(1:j),sep
        if(activityCoeffs) then
            if(truncate .and. len_trim(ident(i)) > 6) nmn = ident(i)(1:5) // plus
            Write(log_i,'("lgF(",a,")")') trim(nmn)
            n = max(11,len_trim(log_i))
            write (IPRN,'(2a)', Advance='No') log_i(1:n), sep
        endif
    EndDo
    if(MSOL > 0) Then
        Do i = nIon+1,MS
            nmn = ident(i)
            if(truncate .and. len_trim(ident(i)) > 8) nmn = ident(i)(1:7) // plus
            Write(conc_i,'("C(",a,")")') trim(nmn)
            if(truncate .and. len_trim(ident(i)) > 6) nmn = ident(i)(1:5) // plus
            Write(log_i,'("lgA(",a,")")') trim(nmn)
            j = max(11,len_trim(conc_i))
            n = max(11,len_trim(log_i))
            write (IPRN,'(4a)', Advance='No') conc_i(1:j),sep, log_i(1:n),sep
        EndDo
    EndIf
    write (IPRN,*)
    tableFirstTime = .false.
EndIf ! tableFirstTime?

! ---------------------------------------------------------
!    write one line with data for the calculated point
! ---------------------------------------------------------
write (IPRN,'(" ")', Advance='No')
! print ionic strength, concentrations and activity coefficients
if(activityCoeffs  .and. temperature > -273.d0) then
    nmn = ""
    if(activityCoeffsModel ==0) nmn = "Davies eqn" // sep
    if(activityCoeffsModel ==1) nmn = "    SIT   " // sep
    if(activityCoeffsModel ==2) nmn = " simpl.HKF" // sep
    Write(press,'(2x,F8.4,1x)') pressure
    if(pressure >= 10.d0) Write(press,'(2x,F8.2,1x)') pressure
    if(pressure >= 100.d0) Write(press,'(2x,F8.1,1x,a)') pressure
    write (IPRN,'(A,1x,F8.2,1x,4A)',Advance='No') nmn(1:11), Temperature, sep, press(1:11), sep
endif
if(Hplus >0) then
    pH = -logA(Hplus)
    write (IPRN,'(2x,F8.4,1x,a)',Advance='No') pH, sep
endif
if(eMinus>0) then
    pe = -logA(eMinus)
    write (IPRN,'(2x,F8.4,1x,a)',Advance='No') pe, sep
    if(temperature > -273.d0) then
        Eh = 1000.d0 * pe * (8.31446261815324d0 * (Temperature + 273.15d0) * ln10) /  96485.3321233100184d0
        write (IPRN,'(2x,F8.2,1x,a)',Advance='No') Eh, sep
    endif
endif
if(activityCoeffs) then
    if(sumM > 99.9d0 .or. sumM < 0.01d0) then
        write(conc_i,'(1PE11.5)') sumM
    else
        write(conc_i,'(F11.7)') sumM
    endif
    write (IPRN,'(1PG11.4,a,E11.4,a,2a)',Advance='No') ionicStrCalc,sep, electricBalance,sep, conc_i(1:11),sep
endif
if(jWater > 0) then
    write (IPRN,'(F11.6,a)',Advance='No') logA(jWater), sep
    if(activityCoeffs) write (IPRN,'(F11.6,a)',Advance='No') phi, sep
endif
Do i = 1,Na
    if(i == jWater) Cycle
    w = max(-999.d0,min(999.d0,logA(i)))
    write (IPRN,'(1PE11.4,a,E11.4,a,0P,F11.6,a)',Advance='No') tot(i),sep,solub(i),sep, w,sep
EndDo
Do i = 1,nIon
    if(i == jWater .or. i == eMinus) Cycle
    write (IPRN,'(1PE11.4,a)',Advance='No') C(i),sep
    if(activityCoeffs) then
        w = max(-999.d0,min(999.d0,logF(i)))
        write (IPRN,'(F11.6,a)',Advance='No') w,sep
    endif
EndDo
If(MSOL > 0) Then
    Do i = nIon+1,MS
        w = max(-999.d0,min(999.d0,logA(i)))
        write (IPRN,'(1PE11.4,a,0P,F11.6,a)',Advance='No') C(i),sep, w,sep
    EndDo
EndIf
write (IPRN,*)
RETURN
END SUBROUTINE tablePrint

!-----------------------------------------------------------------------------
SUBROUTINE FileOPEN
!  Open input (FILIN) and output (FILUT) files,
!  reads command-line parameters
!USE IO, ONLY: UpCase, ###
Implicit NONE

LOGICAL :: YesNo
CHARACTER (len = 3) :: EXT
CHARACTER (len = 20) :: TXT
Integer :: IOERR

dbg=0
debug = .false.
FILIN= 'EC.dat'; FILUT= 'EC.out';  FileRES='EC.csv'

! get file names from command-line (if any)
Call GetCmd
if(FILIN <= ' ') then
    Write(*,'("Input file name not given in command line...")')
    Call Quit
endif

! - - - - - - - - - - - - - - - - - - - - -
! Open the disk file with the input data
Inquire (File=FILIN,Exist=YesNo)
IF(.not.YesNo) Then
    FILUT = FILIN
    Call GetExt (FILUT,EXT)
    call upCase (EXT)
    if(EXT /= 'DAT') then
        FILUT(Len_Trim(FILUT)+1:) = '.dat'
        Inquire (File=FILUT,Exist=YesNo)
        if(.not.YesNo) then
            Write(*,1021) Trim(FILIN), Trim(FILUT)
            1021 FORMAT('? can not find file: ',/5x,'"',A,'"',:/,  &
            '? can not find file: "',A,'"  either.')
            Call ErrStop
            endif
        else ! it was ".dat"
            Write(*,1021) Trim(FILIN)
            Call ErrStop
        endif
    FILIN=FILUT
EndIf
OPEN (UNIT=INFL,FILE=FILIN,STATUS='OLD',IOSTAT=IOERR)
IF (IOERR /= 0) Then
    WRITE(*,1031) IOERR,Trim(FILIN)
    1031 FORMAT('? Error Nbr.',I6,/2X,  &
    'Could NOT open disk file: ',/5x,'"',A,'"')
    Call ErrStop
EndIf

! Get FILPATH
! getParent: Extract from "fName" the path for the file, which is returned in "PATH"
! while the file name without the path is returned in "fName"
FILUT=FILIN
call getParent(FILUT,FILPATH)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! open the output file
FILUT=FILIN
Call GetExt (FILUT,EXT)
FileRES=Trim(FILUT) // '.csv'
FILUT= Trim(FILUT) // '.out'
OPEN(UNIT=IUT,STATUS='UnKnown',FILE=FILUT,IOSTAT=IOERR)
IF(IOERR /= 0) Then
    WRITE(*,1031) IOERR, Trim(FILUT)
    Call ErrStop
EndIf

! - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! print some info
WRITE(*,1051) Trim(FILPATH),Trim(FILIN), Trim(FILUT), Trim(FileRES)
1051 Format (/,"Input, output and table files:",/4x,"path: """,A,"""",/3x,"input: """,A,"""", &
             /2x,"output: """,A,"""",/3x,"table: """,A,"""")
If(dbg /= 0) Write(*,'("Debug printout level: ",i2)') dbg
if((abs(tolHalta0-TOL_HALTA_DEF)/TOL_HALTA_DEF) > 1.d-3) Write(*,'("Relative tolerance in HaltaFall = ",1PE9.2)') tolHalta0
If(activityCoeffsModel >= 0 .and. activityCoeffsModel <= 2) then
    if(activityCoeffsModel == 0) then
        TXT = "Davies eqn."
    else if(activityCoeffsModel == 1) then
        TXT = "SIT"
    else
        TXT = "simpl. HKF"
    endif
    Write(*,'("Activity Coeffs.:  method =",i3," (",A,")")') activityCoeffsModel, trim(TXT)
    if(Temperature > 0.) Write(*,'(24x,"t = ",F6.2,"""C, p =",F8.2," bar")') Temperature, Pressure
    if(abs(ionicStr) > 1.e-12) then
        if(ionicStr>0) Write(*,'(11x,"Ionic strength = ",F6.2," mol/kg")') ionicStr
        if(ionicStr<0) Write(*,'(11x,"Ionic strength = ",F6.2," (calc. at each point)")') ionicStr
    endif
    if((abs(tolLogF0-TOL_HALTA_DEF_LOGF)/TOL_HALTA_DEF_LOGF) > 1.d-3) &
        Write(*,'(18x,"tolLogF = ",1PE9.2)') tolLogF0
    if(activityCoeffsModel == 1 .and. len_trim(ACPath) > 0) &
        Write(*,'("   path to file ""SIT-coefficients.dta""",/8x,"""",A,"""")') trim(ACPath)
EndIf
if(sep /= achar(9)) then
    Write(*,'("Column separator in table file: """,a,"""")') sep
else
    Write(*,'("Column separator in table file: ""tab""")')
endif
WRITE(*,*)

RETURN
END SUBROUTINE FileOPEN

!-----------------------------------------------------------------------------
SUBROUTINE GetCmd
! gets data from the command-line
  USE IO, ONLY : UpCase
  Implicit NONE
  Character (LEN=Max_Path) :: txt, txtU
  Character (LEN=20) :: columnSep
  Logical :: OK
  Integer :: count, i, j, ierr
  FILIN = ''
  Temperature = -273.d0 !###
  ionicStr = 0.d0
  activityCoeffsModel = -1
  ACPath = ''
  count = command_argument_count() ! Minimalist Gnu for Windows
  if(count < 1) then
    WRITE (*,91)
    WRITE (*,92)
    call Quit
  endif
  Do i = 1, count
    call get_command_argument(i, txt)
    OK = .false.
    if(txt(1:1).ne.'-') then
        if(FILIN.le.' ') then
            FILIN = txt
            OK = .true.
            cycle
        else
            write(*,'(" ? Unknown command """,a,"""")') trim(txt)
            stop
        endif
    endif
    txt=txt(2:)
    if(txt(1:1).eq.'?') then ! -?
        WRITE (*,91)
        WRITE (*,95)
        WRITE (*,96)
        call Quit
    endif
    txtU = txt
    call UpCase (txtU)
    if(txtU(1:3).eq.'DBG') then
        debug = .true.
        OK = .true.
    endif
    if(txtU(1:2).eq.'D=') then
        txtU = txtU(3:)
        j=max(1,len_trim(txtU)+1)
        txtU(j:) = ','
        read(txtU(1:len_trim(txtU)),*,iostat=ierr) DBG
        if(ierr.eq.0 .and. DBG>=0 .and. DBG<=10) OK = .true.
        if(debug) write(*,'(3x,"-",a,"   debug printout =",i3)') trim(txt),dbg
    endif
    if(txtU(1:4) == 'KEEP') then
        OK = .true.
        keep = .true.
    endif    
    if(txtU(1:2).eq.'T=') then
        txtU = txtU(3:)
        j=Len_Trim(txtU)+1
        if(j<38 .and. index(txtU,'.').le.0 .and. index(txtU,'E').le.0) then
              txtU(j:) = '.,'
            else
              txtU(j:) = ','
            endif
        read(txtU(1:len_trim(txtU)),*,iostat=ierr) Temperature
        if(ierr.eq.0 .and. Temperature >= 0.d0 .and. Temperature <= 1000.d0) OK = .true.
        if(debug) write(*,'(3x,"-",a,"   temperature =",F8.3)') trim(txt),Temperature
    endif
    if(txtU(1:2).eq.'I=') then
        txtU = txtU(3:)
        j=Len_Trim(txtU)+1
        if(j<38 .and. index(txtU,'.').le.0 .and. index(txtU,'E').le.0) then
              txtU(j:) = '.,'
            else
              txtU(j:) = ','
            endif
        read(txtU(1:len_trim(txtU)),*,iostat=ierr) ionicStr
        if(ierr.eq.0 .and. (abs(ionicStr+1.d0) < 0.001d0 .or. ionicStr >= 0.d0)) OK = .true.
        if(debug) then
        if(ionicStr>=0)write(*,'(3x,"-",a,"   ionic strength = ",f9.4)') trim(txt),ionicStr
        if(ionicStr<0) write(*,'(3x,"-",a,"   ionic strength = ",f9.4" (calc. at each point)")') trim(txt),ionicStr
        endif
    endif
    if(txtU(1:2).eq.'M=') then
        txtU = txtU(3:)
        j=max(1,Len_Trim(txtU)+1)
        txtU(j:) = ','
        read(txtU(1:len_trim(txtU)),*,iostat=ierr) activityCoeffsModel
        if(ierr.eq.0 .and. activityCoeffsModel>=-1 .and. activityCoeffsModel<=2) OK = .true.
        if(debug) write(*,'(3x,"-",a,"   activityCoeffsModel = ",i2)') trim(txt),activityCoeffsModel
    endif
    if(txtU(1:5).eq.'PATH=') then
        ACPath = ''
        if(len_trim(txtU)>5) read(txt(6:len_trim(txt)),'(a)') ACPath     ! -path=directory
        if(len_trim(ACPath)>0) OK = .true.
        if(debug) write(*,'(3x,"-",a)') trim(txt)
    endif
    if(txtU(1:4).eq.'TOL=') then
        txtU = txtU(5:)
        j=Len_Trim(txtU)+1
        if(j<38 .and. index(txtU,'.').le.0 .and. index(txtU,'E').le.0) then
              txtU(j:) = '.,'
            else
              txtU(j:) = ','
            endif
        read(txtU(1:len_trim(txtU)),*,iostat=ierr) tolHalta0
        if(ierr.eq.0 .and. tolHalta0 >= 1.d-10 .and. tolHalta0 <= 0.01) OK = .true.
        if(debug) write(*,'(3x,"-",a,"   relative tolerance =",1PE9.2)') trim(txt),tolHalta0
    endif
    if(txtU(1:5).eq.'TOLA=') then
        txtU = txtU(6:)
        j=Len_Trim(txtU)+1
        if(j<38 .and. index(txtU,'.').le.0 .and. index(txtU,'E').le.0) then
              txtU(j:) = '.,'
            else
              txtU(j:) = ','
            endif
        read(txtU(1:len_trim(txtU)),*,iostat=ierr) tolLogF0
        if(ierr.eq.0 .and. tolLogF0 >= 1.d-10 .and. tolLogF0 <= 0.01) OK = .true.
        if(debug) write(*,'(3x,"-",a,"   tolerance in logF =",1PE9.2)') trim(txt),tolLogF0
    endif
    if(txtU(1:5).eq.'TBLS=') then
        OK = .true.
        columnSep = " ." ! make sure that if an empty string is given, then sep = " "
        if(len_trim(txtU)>5) read(txt(6:len_trim(txt)),'(a)') columnSep
        sep = columnSep(1:1) ! the separator between data columns in the output table
        if(columnSep.eq."\t") sep = achar(9)
    endif
    if(.not.OK) then
      write(*,'("? Unknown command ""-",a,"""")') trim(txt)
      write (*,95)
      call ErrStop
    endif
  End Do
  if(FILIN.le.' ') then
    WRITE (*,91)
    WRITE (*,92)
  endif
  Return
   91 FORMAT(/,'Usage:  EC  Input_File_Name  [-command:value]',/, &
         'Enclose Input_File_Name with double quotes ("file name") it it contains blank space.')
   92 FORMAT('For a list of possible commands type:  EC  -?')
   95 FORMAT(/,'Possible commands are:',/3x, &
   '-?    (prints command-line and input-file help)',/3x, &
   '-dbg  (debug printout)',/3x, &
   '-d=value  (debug printout level in HaltaFall:',/14x,  &
        '>=1  print only errors.',/14x, &
        '>=2  print results (conc. and activity) for all species.',/14x, &
        '>=3  as above, and print data base (equilibrium constants) used.',/14x,  &
        '>=4  as above, and print messages from subroutine HaltaFall.)',/3x, &
   '-m=value  (method used to calculate activity coefficients:',/14x, &
        '-1:ideal solution;  0:Davies Eq.;  1:SIT;  2=simpl.HKF;',/14x, &
        'The temperature is needed)',/3x, &
   '-i=value  (ionic strength; if equal "-1" it will be calculated)',/3x, &
   '-path=directory   (directory containing file "SIT-coefficients.dta",',/22x, &
        'with the specific ion-interaction parameters',/22x, &
        'If not given, the input file path is used)',/3x, &
   '-t=value  (temperature in Celsius; it will be ignored if not needed)',/3x, &
   '-tol=value  (relative tolerance in HaltaFall)',/3x, &
   '-tolA=value  (tolerance in log of activity coefficients;',/17x,'it will be ignored if not needed)',/3x, &
   '-tbls=s  ("s" is a single character used as column separator in',/13x, &
              'the output table, e.g. "," or ";" or "\t" for tab, default ",")',/3x, &
   '-keep  (pause after program execution)',/, &
   'Enclose path-names with double quotes if they contain blank space.',/, &
   'For example:  EC "Fe exp-1" -m=2 -t=25 -path="E:\calcs data" -d=1',/)
   96 FORMAT('Input File: the first part, the chemical system, has the same',/, &
    'format as for the Medusa or Spana software.',/,'The last part is either of:',/2x, &
    'a) N-Points (a positive number of points to calculate)',/5x, &
        'then, for each chemical component, the concentrations in the',/5x, &
        'same format as Medusa/Spana:',/7x, &
        '(T or LA), value (the total conc or the log[activity])',/7x, &
        '(TV, LTV or LAV), min-value, max-value',/25x, &
            '(either the total conc, or log[total conc]',/26x, &
            'or log[activity] are varied between the',/26x, &
            '"min" and "max" values)',/5x, &
        'For example, for the system "H+,CO3-2", the input could be',/5x, &
        '"11, LAV,-5,-12, T,0.01," to calculate eleven points between',/5x, &
        'pH 5 and 12 (both included) with total carbonate = 0.01 mol/kg.',/2x, &
    'b) -N (a negative number, "N" points will be calculated)',/5x, &
        'then, for each point to calculate, the concentration for each',/5x, &
        'chemical component:',/7x, &
            '(T or LA), value (the total conc. or the log[activity])',/5x, &
        'For example, for the system "H+,CO3-2", the input could be',/5x, &
        '"-2, LA,-8,T,0.01, LA,-9,T,0.01," to calculate two points at',/5x, &
        'pH 8 and 9 with total carbonate = 0.01 mol/kg.')
END SUBROUTINE GetCmd

END PROGRAM EC
