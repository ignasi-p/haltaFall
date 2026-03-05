MODULE CHEM
!---------------------
! Data defining a chemical equilibrium system.
! To run Haltafall:
!   - define the values of NA, MS, MSOL
!   - call CHEM_MEM_ALLOC to allocate memory for the arrays in module CHEM
!   - call HALTAF_MEM_ALLOC to allocate memory for the arrays in HaltaFall
!   - make sure there is a module called FACTOR (and any dependent modules)
!   - define a method for activity coefficients calculations
!   - define an initial value for the ionic strength (and temperature etc)
!   - set CONT=.false.
!   - call HaltaCalc
!---------------------
USE IO, ONLY : UpCase
Implicit NONE
Integer, PARAMETER, PUBLIC :: dp = kind(0.d0) ! double precision
Character (len = 1), PARAMETER, PUBLIC :: nl = new_line('a')
! Input variables defining the size of the chemical system:
!  NA = number of components
!  MS = total number of species (components (soluble and solid) + all reaction products (aqueous + solids)
!  SOLIDC = the number of solid components (used only for printout)
!  MSOL = number of solids (components + reaction products).  Note: If some component is a solid phase,
!         then you must add a new solid reaction product for each solid component.  That is, increase MS
!         (total nbr of species) and MSOL (nbr solids) with SOLIDC:
!               MSOL = (solid products)+SOLIDC
!               MS = (Tot. nbr. species)+SOLIDC
!         You have to add the new solid reaction products "i" as follows:
!               DO i = 1, SOLIDC
!                   j = (MS-NA-SOLIDC)+i
!                   k = (NA-SOLIDC)+i
!                   LBETA[j] = 0;
!                   NOLL[k] = .true.
!                   do n = 1,NA
!                       a(j,n)=0.D0
!                       if(n==k) a(j,n)=1.D0
!                   end do
!               END DO
Integer, PUBLIC :: NA, MS, MSOL, SOLIDC
! Input variable used to regulate the debug printout in HaltaFall and Factor
! If larger than zero, different levels of printout information is produced
! depending on DBG (1 to 7).  Information is printed in unit = IOUT.
! But some error messages are printed to the terminal only.
! DBG =0 do not output anything
!  =1 output errors, but no debug information
!  =2 errors and results
!  =3 errors, results and input
!  =4 errors, results, input and debug for procedure fasta()
!  =5 errors, results, input and debug for activity coefficients
! >=6 errors, results, input and full debug print-out
!  Default = 1 (report errors only)
Integer, PUBLIC :: DBG = 0
! Debug output is controlled by DBG and written to any file opened on IOUT.
! But some ERROR messages are printed to the terminal only.
! Note: in Fortran, unit=6 is terminal output, and unit=5 is terminal input;
! setting IOUT=6 will write output to the terminal.
Integer, PUBLIC :: IOUT = 0
! Input array:
! KH(i) =1  Then the total concentration is given to HALTA when called in array TOT(i),
!     and the tolerance is given in TOL (e.g., to 1E-5).  The mass-balance equation
!     has to be solved in subroutine HALTA.
! KH(i) =2  Then log(activity) is given to HALTA in array logA(i).  The total concentration will
!     be calculated by subroutine HALTA.  Except that for H2O (i = jWater) the value of
!     logA(i) may be calculated through activity coefficient calculations.
Integer, PUBLIC, Allocatable :: KH(:)
!  CONT= input variable that must be set =.FALSE. (before calling GASOL)
!        every time that the user changes the value of any of the variables
!        (or array or vector element) in COMMON /HLTF0/
Logical, PUBLIC :: CONT = .false.
!  NOLL(i) = input array: True if the concentration of species "i" is zero ("e-", etc)
!            False if this is a normal species       (i=1...MS)
Logical, PUBLIC, Allocatable :: NOLL(:)
!  LBETA(i) = Input array: log10(Beta) for reaction product "i"
!             (Beta=global equilibrium constant of formation, i=1...MS-NA)
!             (solid phases are those with: i = MS-MSOL+1 ... MS)
Real(dp), PUBLIC, Allocatable :: LBETA(:)
!  A(i,j) = Input array with the formula units (stoichiometry) for reaction product "i" and component "j" (j=1...NA,i=1...MS-NA)
Real(dp), PUBLIC, Allocatable :: A(:,:)
! errFlags (output variable)
!   1 = The numerical solution is uncertain (round-off errors)
!   2 = Too many iterations when solving the mass balance equations
!   3 = Failed to find a satisfactory combination of solids
!   4 = Too many iterations trying to find the solids at equilibrium
!   5 = Some aqueous concentration(s) too large (>UNREASONABLE_CONC): uncertain activity coefficients
!   6 = Activity factors did not converge
!   7 = Calculation interrupted by the user
Integer, PUBLIC :: errFlags
! Input/output array:
! TOT(i)  = total concentration for component "i"       (i=1...NA)
!          (needed on input only if KH(i)=1, calculated on output if KH(i)=2)
Real(dp), PUBLIC, Allocatable :: TOT(:)
! Input/output array:
! logA(i) = log10(activity) for species "i"             (i=1...MS)
!          (needed on input for components only if KH(i)=2,
!           calculated on output for all species (components and reaction products))
Real(dp), PUBLIC, Allocatable :: logA(:)
!  SOLUB(i)= Output array with the solubility for component "i"   (i=1...NA)
Real(dp), PUBLIC, Allocatable :: SOLUB(:)
!  C(i)    = Output array with the concentration for each species "i"   (i=1...MS)
Real(dp), PUBLIC, Allocatable :: C(:)
!  LOGF(i) = Output array with the log10(activity coefficient) for species "i"  (i=1...MS)
Real(dp), PUBLIC, Allocatable :: LOGF(:)
! Input variable with the relative maximum tolerance to solve the mass balance equations
! (used only for components with KH[i]=1). For example: 1e-5.
! Must be between 1e-9 and 1e-2.
! If the total concentration of a component is zero, then a small absolute
! tolerance is used, which depends on the total concentrations for the other
! components, for example = 1e-8 x TOL.
Real(dp), PUBLIC ::  TOL
! Input variable with the requested tolerance (log-10 scale) when calculating
! activity coefficients through iterations.
Real(dp), PUBLIC ::  tolLogF !###
! Output value (log-10 scale) of the lowest achieved variation when calculating
! activity coefficients through iterations (the variation before two iterations
! before the activity coefficients start to diverge).
Real(dp), PUBLIC ::  tolLogFAchieved = -1.d0 !###
! Input variable indicating what species is H2O (if any). See also description of "KH".
Integer, PUBLIC :: jWater
! The model to calculate activity coefficients (ionic strength effects):
!   < 0 for ideal solutions (all activity coefficients = 1)
!   = 0 Davies eqn.
!   = 1 SIT (Specific Ion interaction "Theory")
!   = 2 Simplified HKF (Helson, Kirkham and Flowers)
Integer, PUBLIC :: activityCoeffsModel=-1
! The input ionic strength. This is the input value provided by the user before
! the start of the "HaltaFall" calculation. If this value is negative (e.g. =-1)
! then the ionic strength will be calculated in "Factor". When making a diagram,
! this variable is set equal to "ionicStrengthInput" before each calculation step.
Real(dp), PUBLIC :: ionicStr = -1.D0
! The calculated ionic strength. If the input ionic strength value by the user
! (ChemConcs.ionicStr) before the calculation is negative (e.g. =-1) then
! it is calculated in "Factor" as "ionicStreCalc". Otherwise
! ionicStrCalc = ionicStr.
Real(dp), PUBLIC :: ionicStrCalc = -1.D0
! Sum of molalities of all species in the aqueous solution. This is
! only calculated if the ionic strength (ionicStrength) is negative
! (the ionic strength has to be calculated) and if the activity
! coefficients are also calculated.
Real(dp), PUBLIC :: sumM = 0.D0
Real(dp), PUBLIC :: electricBalance = 0.d0
! The osmotic coefficient of water
Real(dp), PUBLIC :: phi = 1.d0

Integer, PUBLIC, PARAMETER :: MXID=30 ! max length for names
! Input arrays with the names of the chemical species
!   IDENTC(i) for component "i"   (i=1...NA)
!   IDENT(i) for all speccies "i"    (i=1...MS)
Character(LEN=MXID), PUBLIC, Allocatable :: IDENTC(:), IDENT(:)

! The electric charges of all aqueous species, plus the electric charge of
! two fictive species (Na+ and Cl-) that are used to ensure electrically neutral
! aqueous solutions (electroneutrality) when calculating the ionic strength and
! activity coefficients
Integer, PUBLIC, Allocatable :: Z(:)

Integer, PUBLIC :: NX, nIon

Integer, PRIVATE :: nIon2 ! = nIon+2
Real(dp), PARAMETER :: UNREASONABLE_CONC = 200.D0

PRIVATE
PUBLIC CHEM_MEM_ALLOC,CHEM_MEM_FREE, ReadChemSystem, printChemSystem, chargeOf, bareNameOf, &
        errFlagSet,errFlagClear,isErrFlagSet,errFlagsGetMessages,errFlagsToString, UNREASONABLE_CONC

SAVE

CONTAINS

!----------------------------------------------------------------------------
SUBROUTINE CHEM_MEM_ALLOC  ! allocates memory for the arrays
    Implicit NONE
    Integer :: MXB, I, status = 0
 ! NA = number of components
    if (NA <= 0) then
        PRINT *, "Error in module 'CHEM': NA =",NA," (must be >0)"
        Call ErrStop
    endif
 ! MS = total number of species
 !      (components (soluble or solid) + soluble complexes + solid complexes)
    if (MS <= 0) then
        PRINT *, "Error in module 'CHEM': MS =",MS," (must be >0)"
        Call ErrStop
    endif
 ! MSOL = number of solids (components + complexes)
    if (MSOL < 0) then
        PRINT *, "Error in module 'CHEM': MSOL =",MSOL," (must be >=0)"
        Call ErrStop
    endif

    NX = MS - NA - MSOL
    nIon = NA + NX
    if (nIon <= 0) then
        PRINT *, "Error in module 'CHEM': nIon =",nIon," (must be >0)"
        Call ErrStop
    endif

    nIon2 = nIon+2
    MXB = MS - NA      ! nbr of non-component species
    if (.not. ALLOCATED(LBETA)) then
      DO I = 1,1
        ALLOCATE(LBETA(MXB), STAT=status)
        if (status /= 0) exit
        ALLOCATE(A(MXB,NA), STAT=status)
        if (status /= 0) exit
        ALLOCATE(KH(NA),  TOT(NA), SOLUB(NA), STAT=status)
        if (status /= 0) exit
        ALLOCATE(logA(MS), C(MS), LOGF(nIon), NOLL(MS), STAT=status)
        if (status /= 0) exit
        ALLOCATE(IDENTC(NA), IDENT(MS), STAT=status)
        if (status /= 0) exit
        ALLOCATE(Z(nIon2), STAT=status)
      EndDo
    end if
    if (status /= 0) then
        PRINT *, "Something went wrong when trying to allocate arrays in module 'CHEM' (LBeta, etc)"
        Call ErrStop
    end if

    Do i = 1,NA
        TOT(i) =0.d0
        KH(i) = 1
        logA(i)=0.d0
        IDENTC(i) = ""
    EndDo
    Do i = 1,MS
        IDENT(i) = ""
    EndDo
    Do i = 1,nIon
        logF(i) = 0.d0
    End Do
    jWater = 0
    lBeta(1) = huge(1.d0)  ! flag to check that the user provides data
    Z(1) = huge(0)
    C(1) = huge(1.d0) ! flag to check that calculations have been performed

END SUBROUTINE CHEM_MEM_ALLOC

!----------------------------------------------------------------------------
SUBROUTINE CHEM_MEM_FREE  ! deallocates the array memory
    Implicit NONE
    Integer :: I, status = 0
    if (ALLOCATED(LBETA)) then
      DO I = 1,1
        DEALLOCATE(LBETA, STAT=status)
        if (status /= 0) exit
        DEALLOCATE(A, STAT=status)
        if (status /= 0) exit
        DEALLOCATE(KH, TOT, SOLUB, STAT=status)
        if (status /= 0) exit
        DEALLOCATE(logA, C, LOGF, NOLL, STAT=status)
        if (status /= 0) exit
        DEALLOCATE(IDENTC, IDENT, STAT=status)
        if (status /= 0) exit
        DEALLOCATE(Z, STAT=status)
      EndDo
    end if
    if (status /= 0) then
        PRINT *, "Something went wrong when trying to deallocate arrays in module 'CHEM' (LBeta, etc)"
        Call ErrStop
    end if

END SUBROUTINE CHEM_MEM_FREE

!-----------------------------------------------------------------------------
SUBROUTINE ReadChemSystem (Temperature, Pressure)
! READ INPUT DATA
! reads a text data file containing a chemical system,
! the data file is usually created with either
! the Medusa-Hydra or Spana-Database software
!
! The input file has to be already opened.
! Error messages are printed to the terminal.
!**************************************************************
USE IO, ONLY : Max_Path
USE READIR
Implicit NONE
Real(dp), INTENT(OUT) :: Temperature
Real(dp), INTENT(OUT) :: Pressure
Real(dp) :: tC, pBar, ZSUM
Integer :: i, k, j, j2, j3, ix, m, n, NX, &
           NrSolids, NrSolubComp, inFile
Logical :: warn, aqueousSystem
!Character (len = 1), PARAMETER :: nl = new_line('a')
Character (len = 10):: nbr, nbrr
Character (len = 7) :: species
Character (len = Max_Path) :: FILIN
Character (len = MXID) :: IDENT2

inFile = INFL
if(.not.firstInput) inFile = INFL2
INQUIRE (UNIT=inFile,NAME=FILIN) ! get the name of the input file, for error reporting

WRITE (*,'("Reading Input Disk File:   the chemical system... ",A)',Advance='No') CHAR(13)
 ! Read temperature and pressure
 CALL ReadTP (tC,pBar)
 IF(tC > -273.) Then
    Temperature = tC
    Pressure = pBar
 Else
    Temperature = -273.d0
    Pressure = -1.d0
 EndIf
 !print *,"temperature =",Temperature," pressure =",Pressure

 ! ---- read the number of species (chemical components, number of reactions, etc)
 
 nowReading = 'Nbr of chemical components (Na)'
 CALL READI(NA)         ! NA = nr of components
 IF (NA < 1) THEN
   Write(*,'("? Error: Number of components =",I0,"  (must be > 0)",/,"  while reading data file: ",A)') NA,FILIN
   Call ErrStop
 ENDIF
 ! NX = nr of soluble complexes
 nowReading = 'Nbr of complexes (Nx)'
 CALL READI(NX)
 IF (NX < 0) THEN
   Write(*,'("? Error: Number of soluble complexes =",I0,"  (must be >= 0)",/,"  while reading data file: ",A)') NX,FILIN
   Call ErrStop
 ENDIF
 nowReading = 'Nbr of solid reaction products (nrSolids)' 
 CALL READI(NrSolids)
 IF (NrSolids < 0) THEN
   Write(*,'("? Error: Number of solids =",I0,"  (must be >= 0)",/,"  while reading data file: ",A)') NrSolids,FILIN
   Call ErrStop
 ENDIF
 IF ((NX+NrSolids) < 0) THEN
   Write(*,'("? Error: Number of solid reaction products = ",I0,/, &
             "         Number of soluble complexes = ",I0,/, &
             "  (the total number of reactions (soluble + solids) must be >= 0)",/, &
             "  while reading data file: ",A)') NrSolids, NX,FILIN
   Call ErrStop
 ENDIF
 nowReading = 'Nbr of solid chemical components (solidC)' 
 CALL READI(SOLIDC)     ! SOLIDC = nr of solid components
 IF (SOLIDC < 0 .or. SOLIDC > NA) THEN
   Write(*,'("? Error: Number of solid chemical components =",I0,"  must be >= 0 and <=",I0,/, &
             "  (the total number of components)",/,"  while reading data file: ",A)') SOLIDC,NA,FILIN
   Call ErrStop
 ENDIF

! MSOL= number of solids (components + reaction products)
!       because in Halta the solid components are assumed to be
!       fictive soluble components with "zero" concentration
!       (setting noll[] = true) and the corresponding solid complex
!       is added to the list of solids.
MSOL = NrSolids + SOLIDC
! MS = total nr of species (components (soluble or solid) + soluble complexes + solid reaction products)
MS = NA + NX + MSOL
! nx = ms - na - mSol;
! nIon = total number of soluble species in the aqueous solution:
!      all components + soluble complexes
nIon = NA + NX  ! = MS - MSOL
NrSolubComp = NA ! NrSolubComp = number of soluble components

call CHEM_MEM_ALLOC  ! allocate memory for the arrays

! ---- Read names for the components:
jWater=0
DO i=1,NA
    write(nowReading,'("Name of component nr",i4)') I
    CALL READA(IDENTC(i))
    IDENT2=IDENTC(i)
    CALL UPCASE (IDENT2)
    if(IDENT2 == "H2O" .or. IDENT2 == "H2O(L)") then
        jWater=i
        NOLL(i) = .true.
    endif
END DO

! ---- Read data for the reactions
aqueousSystem = .false.
I = 0 ! "I" loops through components, soluble complexes and solid reaction products
Do While (.true.)
    I = I + 1
    if(I > MS) Exit
    IF(I <= NA) THEN ! For components
        IDENT(I) = IDENTC(I)
    ELSE IF (I > MS - SOLIDC) THEN
        IDENT(I) = IDENTC(NA - (MS -I))
    ELSE             ! For reactions
        write(nbr,'(i10)') I
        species = 'solid'
        if(I <= nIon) species = 'soluble'
        nowReading = "Name of "//species//" reaction product "//trim(adjustL(nbr))
        CALL READA (IDENT(I))
        nowReading = 'Formation equilibrium constant for '//species//' reaction product "'//IDENT(I)//'"'
        IX = I-NA
        CALL READR (LBETA(IX))
        DO  J=1,NA
            write(nbrr,'(i10)') J
            nowReading = 'Stoichiometric coeff. nbr '//trim(adjustL(nbrr))//' for '//species &
                            //' reaction product "'//IDENT(I)//'"'
            CALL READR (A(IX,J))
        END DO
    ENDIF

    NOLL(I)=.FALSE.
    IF(I > (NA-SOLIDC) .AND. I <= NA) NOLL(I)=.true.
    IF(IDENT(I)(1:1) == "*") NOLL(I)=.TRUE.
    IDENT2=IDENT(I)
    CALL UPCASE (IDENT2)
    if(IDENT2 == "E-" .or. IDENT2 == "E -" .or. IDENT2 == "H2O" .or. IDENT2 == "H2O(L)") then
        NOLL(I)=.true.
        aqueousSystem = .true.
    endif
    J = len_trim(IDENT2) 
    if(J > 4) then
        J3 = J-3
        if(IDENT2(J3:J) == '(AQ)') aqueousSystem = .true.
    endif
    if(I <= (NA-SOLIDC) .or. (I > NA .and. I <= nIon)) then
        J = len_trim(IDENT2) 
        if(J > 3) then
            j2 = J-2
            if(IDENT2(j2:J) /= '(G)') then
                aqueousSystem = .true.
            endif
        endif
    endif

    IF(I <= nIon) THEN ! Soluble Species:  get the charge
        Z(I)=0
        LOGF(I)=0.D0
        IF(I <= (NA-SOLIDC) .or. I > NA) Z(I) = chargeOf(IDENT2)
        if(Z(I) /= 0) aqueousSystem = .true.
    END IF

  End Do ! while (true)

  ! for HALTA - solid components:
  !      the data for the extra fictive solid reaction products
  !      equilibrium contstant of formation = 1.
  !      and set stoichiometric coefficients for components
  IF(SOLIDC > 0) THEN
    DO m = 1,SOLIDC
        j = (MS - NA - SOLIDC) + m ! j= ((Ms-Na)-solidC)+1 ...(Ms-Na)
        k = (NA - SOLIDC) + m      ! k= (Na-solidC)+1...Na
        NOLL(k) = .true. ! solid components are not aqueous species
        LBETA(j) = 0.d0  ! equilibrium contstant of formation = 1.
        Do n = 1, NA
            A(j,n) = 0.d0
            if(n == k) A(j,n) = 1.d0
        End Do
        n = j + NA ! n= (Ms-solidC)+1 ...(Ms)
        IDENT(n) = IDENTC(k)
        if(IDENT(n)(1:1) == "*") then
            IDENT(n) =  IDENT(n)(2:)
            IDENTC(k) = IDENTC(k)(2:)
            NOLL(n) = .true.
        end if
    END DO
  END IF

  IF(aqueousSystem) THEN
    Do I = 1, nIon
        IDENT2=IDENT(I)
        CALL UPCASE (IDENT2)
        J = len_trim(IDENT2) 
        if(J > 3) then
            if(IDENT2((J-2):J) == '(G)') NOLL(I) = .true.
        endif
    End Do
  ENDIF

! -----------------------
!Check Electric Charges
Z(nIon+1)=1
Z(nIon+2)=-1
!IF (activityCoeffsModel >= 0) THEN
    !Check that all reactions are charge balanced
    warn=.false.
    DO i=NA+1,MS-SOLIDC
        IX=i-NA
        ZSUM=0.
        if(i <= nIon) ZSUM=REAL(-Z(i),dp)
        DO j=1,NA-SOLIDC
          ZSUM = ZSUM + A(IX,j)*REAL(Z(j),dp)
        END DO
        IF(abs(ZSUM)<=5.E-4) CYCLE
        IF(.not.warn) WRITE(*,'("*** Warning:")')
        warn=.true.
        WRITE(*,'(5X,"the formation of ",A," is not charge balanced")') IDENT(I)
    END DO
    IF(warn) WRITE(*,'(" -------")')
!END IF

RETURN

END SUBROUTINE ReadChemSystem

!-----------------------------------------------------------------------------
FUNCTION bareNameOf(speciesName) RESULT(species)
  ! Get the name without the electric charge for a soluble species.
  ! For "A+B " it returns "A+B".
  ! For "Al+++" and "CO3--" it returns "Al" and "CO3".
  ! For "B+-", "B+-+", "B++--" and "B-+++" it returns
  !     "B+",  "B+",   "B++" and   "B-"
  ! For "+++" it returns "+"
  ! Remove also "(aq)" at the end of a name. For a solid, liquid or gas,
  ! remove phase specification such as "(s)", "(g)" etc
  USE IO, ONLY : UpCase
  Implicit NONE
  Character(len = *), INTENT(IN) :: speciesName
  Character(len = :), Allocatable :: species
  Character(len = :), Allocatable :: sName, sNameU
  Integer :: L, ik,ij, ij1, ij2,ij3, ij0, Name_Length
  Logical :: hasSign
  sName = trim(speciesName)
  L = len_trim(sName)
  if(L <= 1) then
    species = sName;  RETURN
  endif
  sNameU = sName
  call UpCase(sNameU)
  L = len_trim(sNameU)
  if(L > 3) then
    if(sNameU(L-2:L) == "(S)" .or. sNameU(L-2:L) == "(C)" .or. sNameU(L-2:L) == "(A)" &
          .or. sNameU(L-2:L) == "(L)" .or. sNameU(L-2:L) == "(G)") then
                species = sName(1:L-3);  RETURN
    endif
  endif
  if(L > 4) then
      if(sNameU(L-3:L) == "(CR)" .or. sNameU(L-3:L) == "(AM)" &
          .or. sNameU(L-3:L) == "(AQ)") then
                species = sName(1:L-4); RETURN
      endif
  endif
  if(L > 5) then
      if(sNameU(L-4:L) == "(VIT)" .or. sNameU(L-4:L) == "(PPT)") then
                species = sName(1:L-5); RETURN
      endif
  endif ! if len >5
  !-- if the last symbol is a letter: then there is no charge
  if(sNameU(L:L) >= "A" .and. sNameU(L:L) <= "Z") then
    species = sName;  RETURN
  endif

  hasSign = .false.
  ij = L
  Do ik = L,1,-1 ! get the last position of either "+" or "-"
    if(sName(ik:ik) == '+' .or. sName(ik:ik) == '-') then
        hasSign = .true.; ij = ik; Exit
    endif
  EndDo
  if(.not.hasSign) then  ! no electrical sign: exit
    species = sName;  RETURN
  endif

  !There is a '+' or '-':  is it located at the end of the name?
  if(ij <= 1 .or. ij < (L-2)) then
    species = sName;  RETURN
  endif

  !Get the length of the name without the electrical charge
  Name_Length = ij-1   ! in case of charge as in "OH-" or "Fe(OH)2+"
                       ! or in the case "CO3-2" and "Y +11"
  IF(ij == L) THEN ! the sign is the last character
    ij1=ij-1            ! note that IJ>1 and IJ1>=1
    IF(sName(ij1:ij1) >= '0' .and. sName(ij1:ij1) <= '9') THEN
        ! there is a number before the sign: as in "CO3 2-" or "Fe(OH)2+"
        ! is there a space before the charge?
        If(ij >= 4) Then ! the name must be, at least, as long as in "X 2-"
            ij2=ij-2
            If(sName(ij2:ij2) == ' ') Then
                ! charge as in "CO3 2-" or "FeOH 2+"
                Name_Length = ij-3
            EndIf
            If(ij >= 5) Then ! if the name is, at least, as long as in "X 12-"
                ij3=ij-3
                if(sName(ij3:ij3) == ' ') then
                    ! charge as in "Y 10-" or "X 10+"
                    Name_Length = ij-4
                endif
            EndIf !IJ>5
        EndIf !IJ>4
    ELSE !sName(ij1:ij1) is not a number
        ! there is not a number before the sign: as in "Cl-" or "Fe++"
        ! is the previous character equal to this? as in "Fe++"
        if(sName(ij1:ij1) == sName(ij:ij) .and. L >= 3) then
            ij0 = 1
            Do ik=ij1,1,-1 ! loop backwards
                if(sName(ik:ik) == sName(ij:ij)) Cycle
                ij0 = ik; Exit
            EndDo
            Name_Length = ij0
        endif ! the previous character equal to this, as in "Fe++"
    END IF !sName(ij1:ij1) >= '0' and sName(ij1:ij1) <= '9' ??
  ELSE   ! the sign is not the last character, such as "A+B" or "Fe+2"
    ! are the following characters numbers?
    Do ik = ij+1,L
        if(sName(ik:ik) >= '0' .and. sName(ik:ik) <= '9') Cycle !OK
        ! not OK, the last "+" or "-" character is followed by non-numbers
        Name_Length = L; Exit
    EndDo
  END IF ! ij = L
  species = trim(sName(1:Name_Length))
  RETURN
END FUNCTION bareNameOf

!-----------------------------------------------------------------------------
FUNCTION chargeOf(speciesName) RESULT(z)
! Returns the charge of a species. For "H+" and "Na +" returns +1;
! for "CO3-2" and "SO4 2-" it returns -2; for "X+10" and "X 10+" returns +10.
! For "A+B " it returns zero. For "Al+++" and "CO3--" it +3 and -2.
! For "B+-", "B+-+", "B++--" and "B-+++" it returns -1,+1,-2, and +3, respectively
! For "+++" it returns +3, and for "A-B-2" it returns -2.
Implicit NONE
Character(len = *), INTENT(IN) :: speciesName
Integer :: z
Character(len = :), Allocatable :: sName
Integer :: L, N, IK,IJ, IJ1,IJ2,IJ3
Logical :: hasSign, ok

  Z=0
  sName = trim(speciesName)
  L = len_trim(sName)
  if(L <= 1) RETURN

  !-- if the last symbol is a letter: then there is no charge
  if((sName(L:L) >= "A" .and. sName(L:L) <= "Z") .or. &
     (sName(L:L) >= "a" .and. sName(L:L) <= "z"))  RETURN

  hasSign = .false.
  ij = L
  Do ik = L,1,-1 ! get the last position of either "+" or "-"
    if(sName(ik:ik) == '+' .or. sName(ik:ik) == '-') then
        hasSign = .true.; ij = ik; Exit
    endif
  EndDo
  if(.not.hasSign) RETURN  ! no electrical sign

  !There is a '+' or '-':  is it located at the end of the name?
  if(ij <= 1 .or. ij < (L-2)) RETURN

  Z=1
  IF(sName(ij:ij) == '-') Z=-1

  IF(ij == L) THEN ! the sign is the last character
    ij1=ij-1            ! note that ij>1 and ij1>=1
    IF(sName(ij1:ij1) >= '0' .and. sName(ij1:ij1) <= '9') THEN
        ! there is a number before the sign: as in "CO3 2-" or "Fe(OH)2+"
        ! is there a space before the charge?
        If(ij >= 4) Then ! the name must be, at least, as long as in "X 2-"
            ij2=ij-2
            If(sName(ij2:ij2) == ' ') Then
                ! charge as in "CO3 2-" or "FeOH 2+"
                READ(sName(ij1:ij1),'(i1)') N
                Z = Z*N;  RETURN
            EndIf
            If(ij >= 5) Then ! if the name is, at least, as long as in "X 12-"
                ij3=ij-3
                if(sName(ij3:ij3) == ' ') then
                    ! charge as in "Y 10-" or "X 10+"
                    READ(sName(ij2:ij1),'(i2)') N
                    Z = Z*N;  RETURN
                endif
            EndIf !IJ>5
        EndIf !IJ>4
    ELSE !sName(ij1:ij1) is not a number
        ! there is not a number before the sign: as in "Cl-" or "Fe++"
        ! is the previous character equal to this? as in "Fe++"
        if(sName(ij1:ij1) == sName(ij:ij) .and. L >= 3) then
            ij2 = 1
            Do ik=ij1,1,-1 ! loop backwards
                if(sName(ik:ik) /= sName(ij:ij)) Exit
                ij2 = ij2 + 1
            EndDo
            Z = Z * ij2;  RETURN
        endif ! the previous character equal to this, as in "Fe++"
    END IF !sName(ij1:ij1) >= '0' and sName(ij1:ij1) <= '9' ??
  ELSE   ! the sign is not the last character, such as "A+B" or "Fe+2"
    ! are the following characters numbers?
    ok = .true.
    Do ik = ij+1,L
        if(sName(ik:ik) >= '0' .and. sName(ik:ik) <= '9') Cycle !OK
        ! not OK, the last "+" or "-" character is followed by non-numbers
        ok = .false.; Exit
    EndDo
    if(ok) then
        READ(sName(ij+1:L),'(i2)') N
        Z = Z*N;  RETURN
    endif
  END IF ! ij = L
  RETURN



   ! ChargeNumber = '  '
   ! IF(IJ == L .AND. IJ > 1) THEN ! charge given as "+", " 2-", or " 10+" (the last character)
      ! IJ1=IJ-1
      ! IF(IDENT(IJ1:IJ1) < '0' .OR. IDENT(IJ1:IJ1) > '9') RETURN
      ! IF(IJ > 2) THEN
         ! IJ2=IJ-2
         ! IF(IDENT(IJ2:IJ2) == ' ') THEN
            ! ChargeNumber = IDENT(IJ1:IJ1)
            ! GO TO 1000
            ! END IF
         ! IF(IJ > 3) THEN
            ! IJ3=IJ-3
            ! IF(IDENT(IJ3:IJ3) == ' ') THEN
               ! ChargeNumber = IDENT(IJ2:IJ1)
               ! GO TO 1000
               ! END IF
            ! END IF
         ! END IF
      ! RETURN
   ! ELSE ! charge given as "-2", or "+11" (not the last character)
      ! LMAX=LEN(IDENT)
      ! IF(IJ > (LMAX-2)) RETURN ! the "+/-" sign must be at least the third position from the end
      ! IMX=IJ+1
      ! IF(IJ <= LMAX-2) IMX=IJ+2
      ! ChargeNumber = IDENT(IJ+1:IMX)
   ! END IF

! 1000 CONTINUE
   ! IF(ChargeNumber(1:1) <= '0' .OR. ChargeNumber(1:1) > '9') ChargeNumber(1:1)=' '
   ! IF(ChargeNumber(2:2) < '0' .OR. ChargeNumber(2:2) > '9') ChargeNumber(2:2)=' '
   ! IF(ChargeNumber(2:2) == ' ') THEN
      ! ChargeNumber(2:2)=ChargeNumber(1:1)
      ! ChargeNumber(1:1)=' '
      ! END IF
   ! IF(ChargeNumber == '  ') RETURN
   ! READ(ChargeNumber,'(i2)') N
   ! Z = Z*N
END FUNCTION chargeOf

!-----------------------------------------------------------------------------
SUBROUTINE ReadTP (Temperature, Press)
! Reads temperature and pressure from the first line of the input file.
! The input file must already be opened in unit "inFile" (="INFL" or "INFL2").
! NOTE: a rewind is performed on the file.  Reading after this subroutine
!       will start from the beginning of the file.
! The first line in the input file must end with " /" (to start a comment)
! followed by "t=xxxx" and "p=xxxx".
! On output Temperature and Pressure:
!  - contain positive (or zero) real values of the
!    temperature and/or pressure found, or
!  - If an error occurs, Temperature = -273 and Press = -1,
!    for example, if the input file is not opened,
!    or if there is no temperature and/or pressure given
!    in the first line of the input file.
USE READIR, ONLY : INFL, INFL2, firstInput
Implicit NONE
Real(Kind(1.D0)), INTENT(OUT)       :: Temperature
Real(Kind(1.D0)), INTENT(OUT)       :: Press

Integer, PARAMETER :: length=200
Character :: cdum*(length),cxxx*(length),q*1
Integer :: I, inFile, ierr, it,ip, last
Real(Kind(1.D0)) :: w
!     this function is true if the character "q" is not part of a number (but it can be a space)
Logical :: nonnum
nonnum(q) = ((q > '9' .and. q /= 'E' .and. q /= 'D') .or.  &
              (q < '0' .and.  &
               q /= '-' .and. q /= '+' .and. q /= ' ' .and. q /= '.'))
Temperature = -273.
Press = -1.
inFile = INFL
if(.not.firstInput) inFile = INFL2
rewind (unit=inFile, iostat=ierr)
if(ierr /= 0) go to 9999
!  read the first line
50 read (inFile,'(a)') cdum
   if(cdum <= ' ') go to 50 ! go to first non-empty line
   i = index (cdum,'/')     ! find '/'
   if (i <= 0 .or. i >= (length-2)) go to 9999  ! not good
   !  is there a "t=" or a "p=" ?
   cxxx = cdum(i+1:)
   cdum = cxxx
   Call UpCase (cdum)
   it = max(index(cdum,'T='),index(cdum,'T ='))
   ip = max(index(cdum,'P='),index(cdum,'P ='))
   if(it > 0) it=it+index(cdum(it:),'=')
   if(ip > 0) ip=ip+index(cdum(ip:),'=')
!  read the temperature
   if(it > 0 .and. it <= length) then
      cxxx = cdum(it:)
      last = len_trim(cxxx)
      if (last >= 1) then
         Do i=1,last      ! find the end of the number
            q = cxxx(i:i)
            if(nonnum(q)) exit
         EndDo
         if(i > 1) then
            i = i-1
            w = -1.
            read(cxxx(1:i),*,iostat=ierr) w
            if(ierr.eq.0 .and. w >= 0 .and. w <= 1000) Temperature = w
         endif
      endif !(last >= 1)
   endif !(it > 0 .and. it <= length)
!  read the pressure
   if(ip > 0 .and. ip <= length) then
      cxxx = cdum(ip:)
      last = len_trim(cxxx)
      if (last >= 1) then
         Do i=1,last      ! find the end of the number
            q = cxxx(i:i)
            if(nonnum(q)) exit
         EndDo
         if(i > 1) then
            i = i-1
            read(cxxx(1:i),*,iostat=ierr) w
            if(ierr.eq.0 .and. w > 0 .and. w <= 1000000) Press = w
         endif
      endif !(last >= 1)
   endif !(ip > 0 .and. ip <= length)
9999 Rewind (unit=inFile, iostat=ierr)
!print *,"temperature =",Temperature," pressure =",Press
Return
END SUBROUTINE ReadTP

!----------------------------------------------------------------------------
SUBROUTINE printChemSystem (UNIT) ! PUBLIC
! prints a description of the chemical system
! if UNIT > 0 printout to file only
! if UNIT <= 0 printout to terminal only
  Implicit NONE
  Integer, INTENT(IN)  :: UNIT
  Character (Len=40), PARAMETER :: &
      LINE = '----------------------------------------'
  Character (len = 10):: nbr1, nbr2
  Character (len = MXID) :: IDENT2
  Character (len = 35):: txt
  Integer :: i,j,js, k, UNT

  if (.not. ALLOCATED(LBETA)) then
    PRINT *, "Memory not allocated for arrays in 'printChemSystem'"
    Call ErrStop
  end if

  UNT = UNIT
  if(UNT <=0) UNT = 6

  write(UNT,'(a)') LINE
  write(nbr1,'(i10)') NA
  nbr1 = adjustL(nbr1)
  write(UNT,'(a)') "Chemical system:"
  if(SOLIDC == 1) then
    txt = ", of which the last one is solid"
  else if(SOLIDC >1) then
    write(nbr2,'(i10)') SOLIDC
    nbr2 = adjustL(nbr2)
    txt = ", of which the last "//trim(nbr2)//" are solid"
  else
    txt = ""
  endif
  write(UNT,'(a)')"   nr components = "//trim(nbr1)//trim(txt)
  !NX = MS - NA - MSOL ! nbr of soluble reaction products
  write(nbr1,'(i10)') NX
  nbr1 = adjustL(nbr1)
  write(nbr2,'(i10)') MSOL
  nbr2 = adjustL(nbr2)
  write(UNT,'(a)')"   reaction products: soluble  = "//trim(nbr1)//", solid = "//trim(nbr2)
  write(UNT,'(a)') "Component nr,  name,  and  noll:"
  Do I = 1, NA
    if(len_trim(IDENTC(I)) <= 0) then
      IDENT2 = '""'
    else
      IDENT2 = IDENTC(I)
    endif
    write(UNT,'(i3,1x,a20,1x,L1)') I,IDENT2,NOLL(I)
  End Do
  If(NX > 0 .or. MSOL > 0) Then
    write(UNT,'(a)') 'Reaction products: name, logBeta, noll, a[]='
    Do I = 1, NX
      J = I + NA
      if(len_trim(IDENT(J)) <= 0) then
        IDENT2 = '""'
      else
        IDENT2 = IDENT(J)
      endif
      write(UNT,'(i3,1x,"(",i3,")",1x,a20,1x,F10.5,1x,L1,1x,999(F8.3))') &
                            J,I,IDENT2,LBETA(I),NOLL(J),(A(I,K), K=1,NA)
    End Do
    js = 0
    Do I = NX+1, (NX+MSOL)
      J = I + NA
      if(len_trim(IDENT(J)) <= 0) then
        IDENT2 = '""'
      else
        IDENT2 = IDENT(J)
      endif
      js = js+1
      write(UNT,'(i3,1x,"(",i3,")",1x,a20,1x,F10.5,1x,L1,1x,999(F8.3))') &
                            J,JS,IDENT2,LBETA(I),NOLL(J),(A(I,K), K=1,NA)
    End Do  
  EndIf ! NX>0 or. MSOL>0
  write(UNT,'("Components:")')
  write(UNT,'("    kh()=",15i4,10(:/6x,15i4))') (kh(j),j=1,na)
  Write(UNT,'("   tot()=",1p,10g14.6,40(:/9x,10g14.6))') (tot(j),j=1,na)
  Write(UNT,'("  logA()=",10f14.6,40(:/9x,10f14.6))') (logA(j),j=1,na)
  write(UNT,'(a)') LINE
  RETURN
END SUBROUTINE printChemSystem

!----------------------------------------------------------------------------
SUBROUTINE errFlagClear (I)
  Implicit NONE
  Integer, INTENT(IN)  :: I
  if(I==1) then
    if((IAND(errFlags,1) == 1)) &
    errFlags = IEOR(errFlags ,1) ! bitwise exclusive OR
  else if(I==2) then
    if((IAND(errFlags,2) == 2)) &
    errFlags = IEOR(errFlags ,2)
  else if(I==3) then
    if(IAND(errFlags,4) == 4) &
    errFlags = IEOR(errFlags ,4)
  else if(I==4) then
    if(IAND(errFlags,8) == 8) &
    errFlags = IEOR(errFlags ,8)
  else if(I==5) then
    if(IAND(errFlags,16) == 16) &
    errFlags = IEOR(errFlags ,16)
  else if(I==6) then
    if(IAND(errFlags,32) == 32) &
    errFlags = IEOR(errFlags ,32)
  else if(I==7) then
    if(IAND(errFlags,64) == 64) &
    errFlags = IEOR(errFlags ,64)
  endif
END SUBROUTINE errFlagClear

!-----------------------------------------------------------------------------
FUNCTION isErrFlagSet (I) RESULT(isFlag)
  Implicit NONE
  Integer, INTENT(IN)  :: I
  Logical              :: isFlag
  isFlag = .false.
  if(I==1) then
    if(IAND(errFlags,1) == 1) isFlag = .true.
  else if(I==2) then
    if(IAND(errFlags,2) == 2) isFlag = .true.
  else if(I==3) then
    if(IAND(errFlags,4) == 4) isFlag = .true.
  else if(I==4) then
    if(IAND(errFlags,8) == 8) isFlag = .true.
  else if(I==5) then
    if(IAND(errFlags,16) == 16) isFlag = .true.
  else if(I==6) then
    if(IAND(errFlags,32) == 32) isFlag = .true.
  else if(I==7) then
    if(IAND(errFlags,64) == 64) isFlag = .true.
  endif
END FUNCTION isErrFlagSet

!-----------------------------------------------------------------------------
SUBROUTINE errFlagSet (I)
  Implicit NONE
  Integer, INTENT(IN)  :: I
  if(I==1) then
    errFlags = IOR(errFlags ,1) ! bitwise Boolean inclusive-OR
  else if(I==2) then
    errFlags = IOR(errFlags ,2)
  else if(I==3) then
    errFlags = IOR(errFlags ,4)
  else if(I==4) then
    errFlags = IOR(errFlags ,8)
  else if(I==5) then
    errFlags = IOR(errFlags ,16)
  else if(I==6) then
    errFlags = IOR(errFlags ,32)
  else if(I==7) then
    errFlags = IOR(errFlags ,64)
  endif
END SUBROUTINE errFlagSet

!-----------------------------------------------------------------------------
FUNCTION errFlagsToString() RESULT(txt) ! txt (len=26)
  Implicit NONE
  Character(len=26) :: txt
  txt = 'errFlags (1 to 7): '
  if(IAND(errFlags,1) == 1) then
    txt = trim(txt)//" 1"
  else
    txt = trim(txt)//" 0"
  endif
  if(IAND(errFlags,2) == 2) then
    txt = trim(txt)//"1"
  else
    txt = trim(txt)//"0"
  endif
  if(IAND(errFlags,4) == 4) then
    txt = trim(txt)//"1"
  else
    txt = trim(txt)//"0"
  endif
  if(IAND(errFlags,8) == 8) then
    txt = trim(txt)//"1"
  else
    txt = trim(txt)//"0"
  endif
  if(IAND(errFlags,16) == 16) then
    txt = trim(txt)//"1"
  else
    txt = trim(txt)//"0"
  endif
  if(IAND(errFlags,32) == 32) then
    txt = trim(txt)//"1"
  else
    txt = trim(txt)//"0"
  endif
  if(IAND(errFlags,64) == 64) then
    txt = trim(txt)//"1"
  else
    txt = trim(txt)//"0"
  endif
END FUNCTION errFlagsToString

!-----------------------------------------------------------------------------
FUNCTION errFlagsGetMessages() RESULT(txt) ! txt (len=400)
  Implicit NONE
  Character (len=400) :: txt
  Character (len=4) :: intTxt
  Character (len=20) :: doubleTxt
  txt = ''
  if(IAND(errFlags,1) == 1) then
    txt = 'The numerical solution is uncertain (round-off errors).'
  endif
  if(IAND(errFlags,2) == 2) then
    if(len_trim(txt)>0) txt = trim(txt)//nl
    txt = trim(txt)//'Too many iterations when solving the mass balance equations.'
  endif
  if(IAND(errFlags,4) == 4) then
    if(len_trim(txt)>0) txt = trim(txt)//nl
    txt = trim(txt)//'Failed to find a satisfactory combination of solids.'
  endif
  if(IAND(errFlags,8) == 8) then
    if(len_trim(txt)>0) txt = trim(txt)//nl
    txt = trim(txt)//'Too many iterations trying to find the solids at equilibrium.'
  endif
  if(IAND(errFlags,16) == 16) then
    if(len_trim(txt)>0) txt = trim(txt)//nl
    write(intTxt,'(i0)') nint(UNREASONABLE_CONC)
    txt = trim(txt)//'Some aqueous concentration(s) too large (>'//trim(intTxt)//'): uncertain activity coefficients.'
  endif
  if(IAND(errFlags,32) == 32) then
    if(len_trim(txt)>0) txt = trim(txt)//nl
    if(abs(tolLogFAchieved)> 1.d0) then
        write(doubleTxt,'(1PE20.4)') tolLogFAchieved
    else
        write(doubleTxt,'(F20.8)') tolLogFAchieved
    endif
    txt = trim(txt)//'Activity factors did not converge. Uncertainty (log10): '//trim(adjustL(doubleTxt))
  endif
  if(IAND(errFlags,64) == 64) then
    if(len_trim(txt)>0) txt = trim(txt)//nl
    txt = trim(txt)//'Calculation interrupted by the user.'
  endif
END FUNCTION errFlagsGetMessages

!----------------------------------------------------------------------------
SUBROUTINE ErrStop ! (PRIVATE)
    Implicit NONE
    Integer :: ierr, i
    WRITE (*, '("Press Enter to continue ...")', Advance='No')
    READ (*,'(I5)',iostat=ierr) i
    STOP 1
END SUBROUTINE ErrStop

END MODULE CHEM
