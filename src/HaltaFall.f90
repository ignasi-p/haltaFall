! For IBM-PC          Version: Jan-1991 (With Activity Coefficients Correction)
!                                This version with a new solid phase selection.
!     ***********
!      HALTAFALL
!     ***********
!
!  "HaltaFall" ("Concentrations and precipitates" in Swedish) calculates the
!  equilibrium composition of a chemical system with fluid and solid phases.
!  From the HALTAFALL Program, Published in:
!  -  N.Ingri, W.Kakolowicz, L.G.Sillen and B.Warnqvist "High-Speed Computers
!     as a Supplement to Graphical Methods - V. HALTAFALL, a General Program
!     for Calculating the Composition of Equilibrium Mixtures".
!     Talanta, Vol.14, 1967, pp.1261-1286
!  -  N.Ingri, W.Kakolowicz, L.G.Sillen and B.Warnqvist "Errata"
!     Talanta, Vol.15, 1968, pp.xi-xii
!  -  R.Ekelund, L.G.Sillen and O.Wahlberg "Fortran editions of Haltafall and
!     Letagrop"  Acta Chemica Scandinavica, Vol.24, 1970, p.3073
!  -  B.Warnqvist and N.Ingri "The HALTAFALL program - some corrections, and
!     comments on recent experience"
!     Talanta, Vol.18, 1971, pp.457-458
!-----------------------------------------------------------------------
!  Subroutines derived by
!     Ignasi Puigdomenech
!     Dept. of Inorganic Chemistry
!     The Royal Institute Of Technology
!     100 44 Stockholm, Sweden
!-----------------------------------------------------------------------
MODULE HALTAF_Module
!-----------------------------------------------------------------------
! See input/output data
! description in module CHEM
!-----------------------------------------------------------------------
! However, in this version the supersaturated solid phases are picked-up
! one at a time, instead of all at once.
!
! The input data is supplied through the CHEM module.
! The equilibrium composition is calculated by method "haltaCalc" and stored
! in the arrays of the CHEM module.
! Note that if the contents of the chemical system is changed _after_
! "haltaCalc" unpredictable results and errors may occur.
! "haltaCalc" is also associated with the module "Factor"
! that calculates the activity coefficients.
!-----------------------------------------------------------------------
USE CHEM
! Debug output is controlled by DBG and written to any file opened on IOUT.
! But some ERROR messages are printed to the terminal only.
! Note: in Fortran, unit=6 is terminal output, and unit=5 is terminal input;
! setting IOUT=6 will write output to the terminal.
Implicit NONE

! all variables are PRIVATE

! -- parameters in alphabetical order:

Real(dp), PARAMETER :: ln10 = log(10.d0)
! maximum number of iterations in procedure kille
Integer, PARAMETER :: ITER_MAX =120 ! it is perhaps ok with 40
! maximum number of "InFall" iterations in procedure fasta()
Integer, PARAMETER ::  ITER_FASTA_MAX = 100
! the value of logA when noCalc is true
Real(dp), PARAMETER :: NOCALC_LOGA = -99999.999D0
! If true, only one solid is allowed to precipitate at each call of fasta();
! if false, all supersaturated solids are allowed to precipitate. In many
! cases 10-20 solids are found to be supersaturated in a system with less
! than 5 chemical components.  This leads to long iterations trying to select
! the few solids that precipitate and to reject the 15 or more solids that
! do not precipitate.
Logical, PARAMETER :: ONLY_ONE_SOLID_AT_A_TIME = .true.
! starting value of steg[] in the iterations of procedure kille
! when dealing with a new chemical system, that is, when cont = false
Real(dp), PARAMETER :: STEG0 = 0.1D0   ! starting value of steg[] in the iterations of procedure kille
! when dealing with the same chemical system as last calculation, that is,
! when cont = true. It must not be smaller than STEGBYT.
Real(dp), PARAMETER :: STEG0_CONT = 0.01D0
! Used in procedure kille(): When the stepwise changes of x
! are smaller than STEGBYT the procedure switches to
! the chord method. Should be small enough to avoid rounding errors.
Real(dp), PARAMETER :: STEGBYT = 0.005D0
! default tolerance when solving mass balance equations
Real(dp), PARAMETER :: TOL_HALTA_DEF = 1.d-4

! -- booleans in alphabetical order:

! false after the first call to haltaCalc()
Logical :: firstTime = .true.
! if true the ongoing calculations will be stopped at the earliest opportunity
Logical :: panic
! true  if matrix "ruta" in procedure fasta() has come out singular
Logical :: singFall

! -- integers in alphabetical order:

! component for which tot[ivar] is being tested and lnA[ivar] adjusted
Integer :: IVAR
!flag to indicate degrees of success or failure from some procedures
Integer ::  INDIK
! iterations counter when calling procedure fasta(), that is,
! how many different sets of mass balance equations have been solved
Integer :: iterFasta =0
! The routine Fasta now cycles through chunks of code in the original code,
! using the value of nextFall as a guide in each case. The meanings
! of nextFall on return to routine Fasta are as follows:
!   nextFall   Next Routine       Line No. in original HaltaFall
!     0        Exit to PROV,          9000, 10000, 25000
!              indik set to ready
!     1        fallProv               14000
!     2        beFall                 16000
!     3        anFall                 24000
!     4        utFall                 17000
!     5        sing                   18000
!     6        fUtt                   20000
Integer :: nextFall
! number of solids present at equilibrium
Integer :: NFALL
! nbr of solids indicated at INFALL in procedure fasta()
Integer :: NFSPAR
! number of solids systematically eliminated at some stage while singFall
! is true in procedure fasta()
Integer :: NUT
! nr of mass balance equations for tot[ia] to be tested (lnA[ivar] to be varied)
! in absence of solids at equilibrium
Integer :: NVA
! nbr of mass balance equations for tot[ia] to be tested (lnA[ivar] to be varied)
! in the presence of solids at equilibrium (nvaf = nva - nfall)
Integer :: NVAF
! For some reason the calculations go faster if they are performed in two steps,
! first with a small tolerance (1.E-3) and then with the user-requested tolerance.
! This is accomplished with the loop "loopTol".
! This variable is either 1 or 2 (first and second loops).
Integer :: tolLoop

! -- doubles in alphabetical order:

! The minimum value of the tolY[] array.
! In some calculations where the tolerance is high, a solid
! might be found both supersaturated and with zero or slightly negative
! concentration, and this leads to a loop. To avoid this situation, if
! the negative concentration of the solid is (in absolute value) less
! than this tolerance, the concentration is assumed to be zero.
Real(dp) :: tolFasta
! independent variable in y(x)=y0 (in procedure kille)
Real(dp) :: X
! dependent variable in y(x)=y0 (in procedure kille)
Real(dp) :: Y
! value for y aimed at in y(x)=y0 (in procedure kille)
Real(dp) :: Y0

! -- boolean arrays in alphabetical order:

! true if lnA[] is calculated from a solubility product (eqn.(14))
Logical, Allocatable :: BER(:)
! true if the solid is assumed to be present at equilibrium
Logical, Allocatable :: FALL(:)
! true if the component is assumed to occur in one or more solids at equilibrium
Logical, Allocatable :: FALLA(:)
! used to accelerate iterations in kille()
Logical, Allocatable :: jumpChord(:)
! mono[] = true for components that only take part in mononuclear soluble
! complexes (a[ix][ia]=0 or +1 for all ix)
Logical, Allocatable :: MONO(:)
! true if any of the stoichiometric coefficients involving this component are negative (<0)
Logical, Allocatable :: NEG(:)
! true if it has been requested to solve the mass balance equation for this component
! (that is, kh=1) but the calculation is not possible as seen from the
! stoichiometric coefficients
Logical, Allocatable :: noCalc(:)
! true if any of the stoichiometric coefficients involving this component are positive (>0)
Logical, Allocatable :: POS(:)
! true if the two components are independent as seen from the stoichiometric coefficients
! (no complex or solid contains both components)
Logical, Allocatable :: OBER(:,:)

! -- integer arrays in alphabetical order:

! control number used in procedures kille and totBer to catch numerical
! rounding errors:  If two practivally equal lnA[ivar] are found (x1 and x2),
! that correspond to calculated values (y1 and y2) higher and lower,
! respectively, than the given tot[ivar] (y0), then x (=lnA[iva]) can not be
! further adjusted, because x1 and x2 are equal. This occurs often if STEGBYT
! is too large.
Integer, Allocatable :: catchRoundingErrors(:)
! fut[i] = ifSpar number of the i'th solid to be eliminated at some stage in
! systematic variation during 'singFall' in procedure fasta().
Integer, Allocatable :: Fut(:)
! component numbers for which tot[ia] is to be calculated
! by means of the solutbility product (1 to nfall)
Integer, Allocatable :: IBE(:)
! iber[j] = j'th of iva numbers to be an ibe (ibe[j] = iva[iber[j]], j= 1 to nfall)
Integer, Allocatable :: Iber(:)
! numbers of the solids present at equilibrium (1 to nfall)
Integer, Allocatable :: IFALL(:)
! nbr of the solids indicated as possible after INFALL in procedure fasta(),
! (needed when 'singFall'
Integer, Allocatable :: IFspar(:)
! (used in procedure INVERT)
Integer, Allocatable :: indexInv(:,:)
! (used in procedure INVERT)
Integer, Allocatable :: iPivot(:)
! control number used in procedures kille and totBer to catch too many iterations
Integer, Allocatable :: ITER(:)
! components numbers for which the mass balance equation has to be solved
! (lnA to be varied) in absence of solids (1 to nva)
Integer, Allocatable :: IVA(:)
! ivaBra[ivar]= the component number to be tested after ivar,
! if the mass balance for tot[ivar] is satisfied
Integer, Allocatable :: IVABRA(:)
! components numbers for which the mass balance equation has to be solved
! (lnA to be varied) when some solids are present at equilibrium (1 to nvaf)
Integer, Allocatable :: IVAF(:)
! ivaNov[ivar]= the component number to be tested after ivar,
! if the mass balance for tot[ivar] is not satisfied
Integer, Allocatable :: IVANOV(:)
! control number for solving the equations in procedure kille
Integer, Allocatable :: KARL(:)
! number of other components which are independent of this one
Integer, Allocatable :: NOBER(:)

! -- double arrays in alphabetical order:

! scale factor for solid phases, used when evaluating supersaturation
Real(dp), Allocatable :: FSCAL(:)
! natural logarithm of the activity of the components
Real(dp), Allocatable :: LnA(:)
! part of ln(c[]) independent of lnA[ivar] (eqn. 4, procedure lnaBas)
Real(dp), Allocatable :: LNBA(:)
! natural logarithm of the formation equilibrium constant of a complex
Real(dp), Allocatable :: LNBETA(:)
! ln(activity coeff.): natural logarithm of the single-ion activity coefficients
Real(dp), Allocatable :: LnG(:)
! natural logarithm of the equilibrium constant for the dissolution of a solid
Real(dp), Allocatable :: LnKF(:)
! term in lnKf', reduced lnKf (eqn 14a, procedure lnaBer)
Real(dp), Allocatable :: LnKMI(:)
! (used in procedure INVERT)
Real(dp), Allocatable :: Pivot(:)
! matrix combining stoichiometric coefficients and "rut1" (eqn.16b)
! at ANFALL in procedure fasta()
Real(dp), Allocatable :: PVA(:,:)
! matrix combining stoichiometric coefficients and "ruta" (eqns 16a, 16b)
! at ANFALL in procedure fasta()
Real(dp), Allocatable :: RUT1(:,:)
! matrix with stoichiometric coefficients for solid phases for the components
! that are "be", that is, those components for which the lnA[] are calculated
! from the solubility products (eqns 11m', 14, 15) at UTFALL in procedure fasta().
! Note that ruta must be inverted.
Real(dp), Allocatable :: RUTA(:,:)
! step for adjusting lnA in procedure kille
Real(dp), Allocatable :: STEG(:)
! tolerance when solving the mass balance equations for tot[]
Real(dp), Allocatable :: TOLY(:)
! totBe[i] = term for component ibe[i] (eqn. 15) used for calculating
! "cf" after FallProv in procedure fasta().
Real(dp), Allocatable :: TotBe(:)
! term in Tot[ia] - C[ia] (or Tot[ia] if Noll[ia]) where ia=ibe[i], i=1 to Nfall (procedure LnABer2)
Real(dp), Allocatable :: TotMi(:)
! term of reduced total concentration in the presence of solids at
! equilibrium (eqn. 16a, procedure lnaBer)
Real(dp), Allocatable :: TOTVA(:)
! value for x below the right value (in procedure kille)
Real(dp), Allocatable :: X1(:)
! value for x above the right value (in procedure kille)
Real(dp), Allocatable :: X2(:)
! the value of x obtained during the last iteration (in procedure kille)
Real(dp), Allocatable :: xOld(:)
! value for y corresponding to x1 (in procedure kille)
Real(dp), Allocatable :: Y1(:)
! value for y corresponding to x2 (in procedure kille)
Real(dp), Allocatable :: Y2(:)

! ---- activity coefficients stuff:

! Used internally (value assigned in HaltaFall):
! If it is .true. activity coefficients (ionic strength effects)
! will be calculated by HaltaFall using using the provided instance of Factor.
! If it is .false. then the calculations in HaltaFall are made for ideal
! solutions (all activity coefficients = 1) and Factor is never called by
! HaltaFall during the iterations.
Logical :: actCoefCalc
! iteration counter when performing activity coefficient calculations
Integer :: iterAc =0
! default tolerance for the log10 of the activity coefficients
Real(dp), PARAMETER :: TOL_HALTA_DEF_LOGF = 0.001
! default maximum number of iterations when calculating activity coefficients
Integer, PARAMETER ::  ITERAC_MAX0 = 200 ! in most cases it is ok with 100 !##
! default maximum number of iterations when calculating activity coefficients
Integer, PARAMETER ::  ITERAC_MAX_SIT = 1000 !##
! maximum number of iterations when calculating activity coefficients
Integer ::  iterAC_MAX
!The default value for the tolerance (in natural log scale) when iterating activity coefficients (lnG[]).
!See "tolLogF" in CHEM_Module
Real(dp), PARAMETER :: TOL_LNG = 0.0023026D0 ! =0.001 in log10 scale
! The tolerance (in log-10 scale) set in the HaltaFall
! when iterating activity coefficients. Initially set to
! "ChemConcs.tolLogF" (the user-requested value).
! Default tolLogF0 = TOL_LNG/ln10 = 0.001.
! For systems where the concentration for ionic species (C[i])
! is large, then tolLnG = (tolLogF0*max(C[i]*z[i]))*ln10.
! For example, if the highest concentration is for Ca+2,
! and C[i] = 12 mol/L, then (tolLnG/ln10) = 0.001 * 12 * 2
! =0.024. But (tolLnG/ln10) is always less than 0.1
!See "tolLogF" in CHEM_Module
Real(dp) :: tolLogF0
! The tolerance (in natural log scale) being used when iterating activity coefficients (lnG)
Real(dp) :: tolLnG
! the ln(activity coeff.) from the previous iteration of activity coefficient
! calculations using procedure factor.
Real(dp), Allocatable :: OldLnG(:)
!In the SIT model or if the ionic strength has to be calculated,
!activity coefficients have a tendency to oscillate between iterations.
!OldOldLnG is true if
! - lnG for the previous-to-the-previous iteration
!   was larger than the lnG of the previous iteration
Logical(dp), Allocatable :: largerOldLnG(:,:)

Logical :: firstProv

PRIVATE

PUBLIC :: HaltaCalc, HaltaCancel, HALTAF_MEM_ALLOC, HALTAF_MEM_FREE, printConcs

SAVE

CONTAINS

!----------------------------------------------------------------------------
SUBROUTINE HALTAF_MEM_ALLOC  ! allocates memory for the arrays
  Integer :: I, status = 0
 ! NA = number of components
    if (NA <= 0) then
        Write(*,'("Error in ''HALTAF_Module'': NA =",I0," (must be >0)")') NA
        Call ErrStop
    endif
 ! MS = total number of species
 !      (components (soluble or solid) + soluble complexes + solid complexes)
    if (MS <= 0) then
        Write(*,'("Error in ''HALTAF_Module'': MS =",I0," (must be >0)")') MS
        Call ErrStop
    endif
 ! MSOL = number of solids (components + complexes)
    if (MSOL < 0) then
        Write(*,'("Error in ''HALTAF_Module'': MSOL =",I0," (must be >=0)")') MSOL
        Call ErrStop
    endif

    firstTime = .true.
    CONT = .false.

    if (.not. ALLOCATED(X1)) then
      DO I = 1,1
        ALLOCATE(X1(NA),X2(NA),Y1(NA),Y2(NA),STEG(NA),KARL(NA),ITER(NA),catchRoundingErrors(NA),xOld(NA), STAT=status)
        if(status /= 0) exit
        ALLOCATE(jumpChord(NA),NOBER(NA),OBER(NA,NA), STAT=status)
      EndDo
    end if
    if (status /= 0) then
        Write(*,'("Something went wrong when trying to allocate arrays in ''HALTAF_Module'' (X1, etc)")')
        Call ErrStop
    end if
    Do i = 1,NA
        X1(i)=0.d0; X2(i)=0.d0; Y1(i)=0.d0; Y2(i) = 0.d0
    EndDo

    if (.not. ALLOCATED(LNBETA)) then
      DO I = 1,1
        ALLOCATE(LNBETA(NX),LnA(MS),LNBA(NX),TOLY(NA),oldLnG(nIon),LnG(nIon), STAT=status)
        if(status /= 0) exit
        ALLOCATE(IVA(NA+1),IVABRA(NA+1),IVANOV(NA), STAT=status)
        if(status /= 0) exit
        ALLOCATE(MONO(NA), STAT=status)
        if(status /= 0) exit
        ALLOCATE(noCalc(NA),POS(NA),NEG(NA),largerOldLnG(4,nIon), STAT=status)
      EndDo
    end if
    if (status /= 0) then
        Write(*,'("Something went wrong when trying to allocate arrays in ''HALTAF_Module'' (LnBeta, etc)")')
        Call ErrStop
    end if
    Do i = 1,nIon
        LnA(i) = 0.d0
        LnG(i) = 0.d0
        oldLnG(i) = 0.d0
    End Do

    if(MSOL > 0) then
        if (.not. ALLOCATED(TotMi)) then
        Do I=1,1
            ALLOCATE(TotMi(NA),TotBe(NA),Pivot(MSOL), &
                 Fut(NA+1),Iber(NA+1),IFspar(NA),iPivot(MSOL),indexInv(MSOL,2), STAT=status)
            if(status /= 0) exit
            ALLOCATE(BER(NA),FALLA(NA),IBE(NA),IFALL(NA),IVAF(NA), STAT=status)
            if(status /= 0) exit
            ALLOCATE(FALL(MSOL),LNKF(MSOL),LNKMI(MSOL),TOTVA(NA),FSCAL(MSOL), STAT=status)
            if(status /= 0) exit
            ALLOCATE(PVA(NA,NX), STAT=status)
            if(status /= 0) exit
            ALLOCATE(RUTA(MSOL,MSOL),RUT1(NA,MSOL), STAT=status)
        EndDo
        end if
        if (status /= 0) then
          Write(*,'("Something went wrong when trying to allocate arrays in ''HALTAF_Module'' (TotMi, etc)")')
          Call ErrStop
        end if
    endif ! (MSOL > 0)

END SUBROUTINE HALTAF_MEM_ALLOC

!----------------------------------------------------------------------------
SUBROUTINE HALTAF_MEM_FREE  ! deallocates the array memory
    Integer :: I, status = 0

    firstTime = .true.

    if (ALLOCATED(X1)) then
      DO I = 1,1
        DEALLOCATE(X1,X2,Y1,Y2,STEG,KARL,ITER,catchRoundingErrors,xOld, STAT=status)
        if(status /= 0) exit
        DEALLOCATE(jumpChord,NOBER,OBER, STAT=status)
      EndDo
    end if
    if (status /= 0) then
        Write(*,'("Something went wrong when trying to deallocate arrays in ''HALTAF_Module'' (X1, etc)")')
        Call ErrStop
    end if

    if (ALLOCATED(LNBETA)) then
      DO I = 1,1
        DEALLOCATE(LNBETA,LnA,LNBA,TOLY,oldLnG,LnG, STAT=status)
        if(status /= 0) exit
        DEALLOCATE(IVA,IVABRA,IVANOV, STAT=status)
        if(status /= 0) exit
        DEALLOCATE(MONO, STAT=status)
        if(status /= 0) exit
        DEALLOCATE(noCalc,POS,NEG,largerOldLnG, STAT=status)
      EndDo
    end if
    if (status /= 0) then
        Write(*,'("Something went wrong when trying to deallocate arrays in ''HALTAF_Module'' (LnBeta, etc)")')
        Call ErrStop
    end if

    if(MSOL > 0) then
        if (ALLOCATED(TotMi)) then
        Do I=1,1
            DEALLOCATE(TotMi,TotBe,Pivot,Fut,Iber,IFspar,iPivot,indexInv, STAT=status)
            if(status /= 0) exit
            DEALLOCATE(BER,FALLA,IBE,IFALL,IVAF, STAT=status)
            if(status /= 0) exit
            DEALLOCATE(FALL,LNKF,LNKMI,TOTVA,FSCAL, STAT=status)
            if(status /= 0) exit
            DEALLOCATE(PVA, STAT=status)
            if(status /= 0) exit
            DEALLOCATE(RUTA,RUT1, STAT=status)
        EndDo
        end if
        if (status /= 0) then
          Write(*,'("Something went wrong when trying to deallocate arrays in ''HALTAF_Module'' (TotMi, etc)")')
          Call ErrStop
        end if
    endif ! (MSOL > 0)

END SUBROUTINE HALTAF_MEM_FREE

!----------------------------------------------------------------------------
SUBROUTINE HaltaCalc
! Calculates the equilibrium composition of a Chemical System.
! Error reporting:  the variable "errFlags" is set.
! Debug output is controlled by DBG and written to IOUT.
! if IOUT > 0 printout to file only
! if IOUT <= 0 printout to terminal only
! But some error messages are printed to the terminal only.
! See also: Factor, haltaCancel
!
! VARIABLES ----------------------------------------------------------------
! Ivar= nr of component for which Tot(Ivar) is being tested and LnA(Ivar)
!       adjusted (varied).
!---------------------------------------------------------------------------
!  Factor= a subroutine supplied by the user.
!  Logf(i)= log of activity coeff. of:
!  lnG(i)= natural log of activity coeff. of:
!  LnA(i)=  natural log of activity of:
!  logA(i)= log of activity of:
!  C(i)= concentrations of:
!                                    components            1 < i <= NA
!                                    soluble complexes  NA+1 < i <= NA+NX
!                                    solids          NA+NX+1 < i <= NA+NX+MSOL
!---------------------------------------------------------------------------
Integer :: ia,n,RVA
Integer :: MAXtolSPAN
Real(dp) :: tolD, logTol, tolStep, w, minTot
Integer :: tolSpan, ivaIndex
Character(len=100) :: txt
! For some reason the calculations go faster if they are performed in two steps,
! first with a large tolerance (e.g. 1.E-3) and then with the user-requested
! tolerance. This is accomplished with "loopTol".
! To remove this loop, set TOL0 < 1.E-10 in next line  */
Real(dp), PARAMETER :: TOL0 = 1.D-2 !TOL0 = 1.D-2 !  log10(TOL0) = -2

If(DBG >= 3) &
    Write(IOUT,'("--- haltaCalc(): Starting calculation with new concentrations.",/,"    debug level = ",i0,"  cont = ",L1)')  &
                   dbg,cont

panic = .false.

DO ia=1,NA ! check values of kh[]
    If(KH(ia) == 1 .or. KH(ia) == 2) Cycle
    Write(IOUT,'("Error in HaltaCalc: KH(",i0,")=",i0,", allowed values are 1 or 2")') ia, KH(ia)
    Call ErrStop
    RETURN
END DO

Call HaltaInit

firstTime = .false.
iterFasta=0

if(abs(TOL) < 1.d-20) TOL = TOL_HALTA_DEF
if(abs(tolLogF) < 1.d-20) tolLogF = TOL_LNG / ln10

!-------------------
tolLoop = 0;
DO WHILE (.true.) ! loopTol:
!-------------------

    tolD = min(max(abs(TOL),1.d-14),1.d-2) ! something between 1.E-14 and 1.E-2
    tolLogF0 = min(max(abs(tolLogF),1.d-10),1.d-2) ! something between 1.E-10 and 1.E-2
    iterAC_MAX = ITERAC_MAX0
    if(activityCoeffsModel == 1) iterAC_MAX = ITERAC_MAX_SIT ! for SIT
    if(tolLoop == 0) then
        if(tolD < TOL0)then
            tolD = TOL0
            ! tolLogF0 should not be much lower than "tol" for the mass balance
            w = tolD * 0.1d0
            tolLogF0 = max(tolLogF0, w + w * 1.d-10)
        else
            tolLoop = tolLoop +1
        endif
        iterAC_MAX = ITERAC_MAX0 / 2
    endif
    tolLoop = tolLoop +1
    errFlags = 0 ! clear all error flags

    If(DBG >= 3 .and. TOL0 > 1.D-10) Then
        if(tolLoop ==1) then
                txt = " (preliminary calculation with higher tolerance)"
            else if(tolLoop ==2) then
                txt = " (with user-requested tolerance)"
            else
                txt = ""
            endif
        Write(IOUT,'(/,"- - - - -  Tolerance loop: ",I0," (of 2),  ITERAC_MAX=",I0," activityCoeffsModel=",I0,A,/)') &
                tolLoop, iterAC_MAX, activityCoeffsModel, trim(txt)
    EndIf

! NYA1  Get LnA() and Toly()

    Do ia=1,NA
        lnA(ia) = ln10 * logA(ia)
    EndDo

    !Smaller tolerances are set to components solved in "inner" iteration loops.
    !  The maximum decrease in the tolerance is "MAXtolSPAN" (in log-units).
    !  For example, if tol=1e-5, the tolerance for the component in the inner
    !  loop might be 1e-7, depending on the number of equations to solve (nva).

    MAXtolSPAN = 2 ! log10-units
    if(NVA <= 3) MAXtolSPAN = 1
    logTol = log10(tolD)
    If(DBG >= 3) &
        Write(IOUT,'("Max tolerance to solve mass-balance eqns. (log) = ",f6.3,", max tolerance span (log) = ",i0)') &
                logTol,MAXtolSPAN
    ! Calculate tolStep (tolerance step) from tolSpan, minTot:
    !    minTot will contain the smallest absolute value of non-zero total concentration
    minTot = huge(1.d0)
    tolStep = 0 ! the step size in log-units for the change in tol between components
    If(nva > 1) Then
        ! find a value between 1 and 3
        ! if nva = 2  then tolSpan = 1
        ! if nva = 3  then tolSpan = 2
        ! if nva >= 4  then tolSpan = 3
        tolSpan = min(MAXtolSPAN, nva-1) !a value between 1 and MAXtolSPAN
        ! if nva = 2  then tolSteg = 1
        ! if nva = 3  then tolSteg = 1
        ! if nva = 4  then tolSteg = 1
        ! if nva = 5  then tolSteg = 0.75
        ! if nva = 6  then tolSteg = 0.6  etc
        tolStep = dble(tolSpan)/dble(nva-1)
        ! get "minTot": the smallest and largest absolute value non-zero total concentrations
        do ia = 1,NA
            if(KH(ia)==1 .and. abs(TOT(ia))>1.d-50 .and. .not.noCalc(ia)) &
                minTot = min(abs(TOT(ia)),minTot)
        enddo
    End If ! NVA >1

    ! minTot value is used to calculate tolerance when a tot.conc.=0 (zero);
    ! use 1e-6 as minimum concentration if the smallest non-zero total concentration is too high
    minTot = min(minTot, 1.d-6)
    If(DBG >= 3) Write(IOUT,'("Mass-balance tolerances (relative):")')
    tolFasta = huge(1.d0)
    Do ia=1,NA
      If(KH(ia)==1 .and. .not.noCalc(ia)) Then
          w = logTol
          If(nva > 1) Then
            ivaIndex = 1
            Do ivar=1,nva
                If(iva(ivar) == ia) Then
                  ivaIndex = ivar
                  Exit
                End If
            End Do
            ivaIndex = nva - ivaIndex
            w = logTol - tolStep * dble(ivaIndex)
         End If
         w = 10.D0**w
         If(DBG >= 3) Write(IOUT,'(1P,G13.4)', Advance='No') w
         ! if the total concentration is zero, use a minimum value
         TOLY(ia)= max(2.d-14,(w * max(abs(TOT(ia)),minTot)))
         tolFasta = min(tolFasta,TOLY(ia))
      Else ! kh = 2
         If(DBG >= 3) Write(IOUT,'("  NaN")', Advance='No')
         TOLY(ia) = huge(1.d0)
      End If ! kh?
    End Do ! ia

    If(DBG >= 3) Then !  ---- debug ----
        Write(IOUT,'(/,"Mass-balance tolerances (absolute):")')
        do ia=1,NA
            w = TOLY(ia)
            if(w > 1.d+100) w = 0.d0
            Write(IOUT,'(1P,G13.4)', Advance='No') w
        enddo
        Write(IOUT,*)
        if(activityCoeffsModel >=0 .and. activityCoeffsModel <=2) &
            Write(IOUT,'("Activity coeffs. tolerance (initial): ",1pE11.3," (log-10 scale)")') tolLogF0
        Write(IOUT,'("Continuation run: ",L1,", STEG0=",F7.4,", STEG0_CONT=",F7.4,", STEGBYT=",F7.4)') &
                    CONT,STEG0,STEG0_CONT,STEGBYT
        write(IOUT,'("Components:")')
        write(IOUT,'("    kh()=",15i4,10(:/6x,15i4))') (kh(ia),ia=1,na)
        Write(IOUT,'("   tot()=",1p,10g14.6,40(:/9x,10g14.6))') (tot(ia),ia=1,na)
        Write(IOUT,'("  logA()=",1p,10g14.6,40(:/9x,10g14.6))') (logA(ia),ia=1,na)
    EndIf

    ! NVA=0   no mass balance equation needs to be solved
    If(NVA <= 0) Then
        if(DBG >= 4) Write(IOUT,'("NVA = 0")')
        iterFasta=0
        IVAR=1
        Call lnaBas
        iterAc = 0
        call errFlagClear(6)
        Do While (.true.)
            Call Cber
            if(.not.actCoefCalc) exit
            if(actCoeffs()) exit !ok!
            if(isErrFlagSet(6)) exit !activity factors did not converge
        End Do
        Call Fasta
        Call Nog
        if(DBG >= 4) Write(IOUT,'("HaltaCalc returns",/)')
        Return
    End If ! if NVA=0

    iterAc = 0 ! iteration counter when performing activity coefficient calculations
    firstProv = .true.
    singFall=.false.

! SLINGOR:  ("loops" in Swedish)

10000 Continue
    ! Either:
    ! 1- calculation starts (a first guess is made for the unknown lnA[])
    ! 2- after fasta()-ANFALL when a different set of solids is found,
    !    this is the starting point to find the new lnA[] with the new set of solids
    ! 3- after activity coefficient corrections, find the new lnA[]
    if(panic) Go To 90000 !Nog
    If(DBG >= 4) Then
        Write(IOUT,'("--- HaltaCalc at Slingor, nva=",i0,", nfall=",i0,", nvaf=",i0)') nva,nfall,nvaf
        If(DBG >= 4) Then
            Call printArrays (.true.,.false.) ! print iva(), ivaBra(), ivaNov()
            ! print ifall() and ivaf()
            !                                ifall_,  iber_,  fall_, fallA_,  ber_,   ivaf_,  ifSpar_,fut_
            if(nvaf>0) Call printArraysFasta(.true., .true.,.false.,.false.,.true., .true., .false.,.false.)
        End If
    End If

    DO rva=1,NVA
        ia=IVA(rva)
        KARL(ia)=1
        if(.not.CONT) then
            STEG(ia)=STEG0
        else
            STEG(ia)=STEG0_CONT
        endif
        iter(ia) = 0
        catchRoundingErrors(ia) = 0
        jumpChord(ia) = .true.
    End Do

    IVAR=IVA(1)

    If(DBG >= 4) Write(IOUT,'("SLINGOR, new ivar=",i0)') ivar

    If(NFALL > 0)  Call lnaBer
    Call lnaBas
    !tjat = true;
    If(MSOL >0) Then
      If(BER(IVAR))  Then ! calculate lnA from a solubility product
        If(DBG >= 4) Then
          Write(IOUT,'("   ber(",i0,")=true (lnA calc from a solubility product),  x=",F15.10,"  tjat=false")') ivar,x
          Write(IOUT,'("   lnA(",i0,")=",F22.16)') ivar,lnA(ivar)
        EndIf
        Call Cber
        If(NVAF > 0) Call lnaBer2
        ! tjat = .false.
        ! if(actCoefCalc) dummy = actCoeffs()
        GO TO 12000 !Prov (tjat = false)
      Else ! not BER(IVAR)
        ! tjat = .true.
        if(NX == 0 .and. NOLL(ivar)) then
            lnA(ivar) = 0.d0 !1.D0
            ! tjat = .false.
            ! if(actCoefCalc) dummy = actCoeffs()
            GO TO 12000 !Prov (tjat = false)
        endif
        If(DBG >= 4) Write(IOUT,'("   ber(",i0,")=false (lnA NOT from a solubility product),  x=",F15.10,"  tjat=true")') ivar,x
      End If
    End If ! MSOL >0

! TJAT  ( = "repeat" or "nagg" in Swedish)

11000 Continue ! Solve the mass balance equation for component: IVAR
        if(panic) Go To 90000 !Nog
        If(DBG >= 4) write(IOUT,'("TJAT; ivar=",i0,"  tjat=true")') ivar

    ! IF(tjat)
        LnA(IVAR)=X
        If(MSOL > 0) Then
            If(FALLA(IVAR))  Then
                Call lnaBer
                Call lnaBas
            End If
        End If
        If(IVAR /= IVANOV(IVAR))  Then
            ! Clear "too Many Iterations" for "inner" loops if they are "dependent" with this ivar
            Do rva = 1, NVA
                if(iva(rva) /= ivar) Cycle
                if(rva == 1) Exit ! do not clear anything for the 1st iva
                do n = 1,rva
                    ia = iva(n)
                    if(ia /= ivar .and. .not.ober(ivar,ia)) then
                        if(DBG >= 6) &
                           Write(IOUT,'("   setting  iter[",i0,"]=0  and  catchRoundingErrors[",i0,"]=0")') ia,ia
                        iter(ia) = 0
                        catchRoundingErrors(ia) = 0
                        jumpChord(ia) = .true.
                    endif
                enddo
                Exit !do
            End Do

            IVAR = IVANOV(IVAR)

            If(DBG >= 6) Write(IOUT,'("   Tjat, new ivar=",i0)') ivar
            Call lnaBas
            If(MSOL > 0) Then
                If(BER(IVAR)) Then
                    Call Cber
                    If(NVAF > 0) Call lnaBer2
                    GO TO 12000 !Prov
                End If
            End If
        End If ! IVAR /= IVANOV(IVAR)

        Call Cber
        If(NVAF > 0) Call lnaBer2

        Call totBer ! Returns:
            !  indik=1 not ok but it is a component only involved in
            !     mononuclear reactions and there are no solids present
            !  indik=2 ok (the Y is equal to Y0 within the tolerance)
            !  indik=3 not ok and it is either not a mono component or
            !     there are solids present
            !  indik=4 if too many iterations (iter(ivar) > ITER_MAX+1)
        If(INDIK==1) Go To 11000 !Tjat
        If(INDIK==2 .or. INDIK==4) Go To 12000 !Prov
        If(INDIK==3) Then
            Call Kille ! ADA
            GO TO 11000 !Tjat
        End If
    !End If tjat

! PROV  ( = "test" in Swedish)

12000 Continue ! The mass balance for component: "IVAR" was satisfied
    if(panic) Go To 90000 !Nog
    If(DBG >= 4) Then
      Write(IOUT,'("PROV; old ivar=",i0,", new ivar=",i0,1x,A,"  firstProv=",L1)') ivar,ivaBra(ivar),errFlagsToString(),firstProv
      Call printArrays (.false.,.true.) ! print lnA
    End If

    IVAR=IVABRA(IVAR)

    If(IVAR <= 0) Then
        ! This was the last ivar: check activity coefficients and the solid phases
        ! Calculate activity coefficients
        If(firstProv .and. actCoefCalc) Then
            If(.not. actCoeffs()) Then !not ok (act.coeffs. changed)
                singFall = .false. ! solid phases might change as well...
                Go To 10000 ! Slingor
            EndIf
        EndIf
        ! ok, act.coeffs. not changed
        firstProv = .false.
        !if(isErrFlagSet(6)) iterAc = 0 !###

        ! check for "too many iterations" of the outer loop
        if(iter(IVA(NVA)) > ITER_MAX+1) Go To 90000 !Nog

        ! Solid phases
        Call Fasta ! Returns:
            !  indik=2 if a new set of solids is chosen
            !  indik=3 if ok: either no solids precipitate in the system,
            !      or no solubility products are exceeded and no solid
            !      assumed present has zero or negative concentration
            !  indik=4 if too may solid phase combinations have been tested
            !  indik=1 if "tji" (all combinations of solids give unsoluble equations,
            !      i.e. negative determinants)

        if(panic) Go To 90000 !Nog

        If(DBG >= 4) Then
                Write(IOUT,'("haltaCalc at PROV after Fasta() INDIK=",I0,1x,A)') INDIK,errFlagsToString()
                Call printArrays(.false.,.true.) ! print lnA
        End If
        If(INDIK /= 2) Then ! "ok" (or can not go on) in fasta()
            If(INDIK==4) Go To 90000 ! too may solid phase iterations: NOG
            if(.not.actCoefCalc) Go To 90000 !Nog
            ! Calculate activity coefficients
            if(actCoeffs()) Go To 90000 !ok, act.coeffs. not changed: NOG
            !not ok (act.coeffs. changed)
            singFall = .false. !solid phases might also change
            Go To 10000 ! Slingor
        Else If(INDIK == 2) Then ! not "ok" in fasta()
            ! new set of solids: clear the "too many iterations" flag
            iterAc = 0
            Call errFlagClear(1) ! the numerical solution is uncertain
            Call errFlagClear(2) ! too many iterations when solving the mass balance equations
            Call errFlagClear(5) ! some aqueous concentration(s) are too large (>20): uncertain act. coeffs
            Call errFlagClear(6) ! activity factors did not converge
            Go To 10000 !Slingor
        EndIf

    End If ! IVAR <= 0

    If(MSOL > 0) Then
        ! ber()=true if lnA is calculated from a solid
        If(BER(IVAR)) Go To 12000 !Prov
    End If

    Call lnaBas

    Call totBer ! Returns:
        !  indik=1 not ok but it is a component only involved in
        !     mononuclear reactions and there are no solids present
        !  indik=2 ok (the Y is equal to Y0 within the tolerance)
        !  indik=3 not ok and it is either not a mono component or
        !     there are solids present
        !  indik=4 if too many iterations (iter(ivar) > ITER_MAX+1)
    If(INDIK==1) Go To 11000 !Tjat
    If(INDIK==2) Go To 12000 !Prov
    If(INDIK==3) Then
        Call Kille ! ADA
        GO TO 11000 !Tjat
    End If
    If(INDIK==4) Then
        iter(IVAR)=0
        catchRoundingErrors(IVAR)=0
        jumpChord(IVAR) = .true.
        Go To 12000 !Prov
    End If

! NOG   ( = "enough" in Swedish)

90000 Continue ! This is the end.  Calculate logA(i),Logf(i),Tot(ia), and Solub(ia)
    If(DBG >= 3) Write(IOUT,'("NOG")')

    If(panic) Then
        Call errFlagSet(7)
        If(DBG >= 1) Then
            Write(IOUT,'(" *****  Interrupted by the user!  *****")')
            EXIT ! tolLoop
        End If
    End If

    Call NOG

!-------------------
    IF(tolLoop > 1) EXIT
END DO ! WHILE (.TRUE.) = tolLoop
!-------------------

! 1: the numerical solution is uncertain (round-off errors)
! 2: too many iterations when solving the mass balance equations
! 3: failed to find a satisfactory combination of solids
! 4: too many iterations trying to find the solids at equilibrium
! 5: some aqueous concentration(s) are too large (20) uncertain activity coefficients
! 6: activity factors did not converge.
! 7: calculation interrupted by the user.
do rva = 1,NVA
    if(catchRoundingErrors(iva(rva)) >= 3) call errFlagSet(1)
    if(iter(iva(rva)) > ITER_MAX+1) call errFlagSet(2)
end do

CONT = .true.
if(isErrFlagSet(3) .or. isErrFlagSet(4)) CONT = .false.

! DBG =0 do not output anything
!  =1 output errors, but no debug information
!  =2 errors and results
!  =3 errors, results and input
!  =4 errors, results, input and debug for procedure fasta()
!  =5 errors, results, input and debug for activity coefficients
! >=6 errors, results, input and full debug print-out
!  Default = 1 (report errors only)
If(DBG >= 3) Then
    Write(IOUT,'("---- HaltaCalc returns;  cont = ",L1,", ",A,/,'// &
        '"   solid combinations = ",I0,/,'// &
        '"   activity coefficient iterations: ",I0)') &
            CONT, trim(errFlagsToString()), iterFasta, iterAc
    txt = ""
    do rva = 1,NVA
        if(catchRoundingErrors(iva(rva)) >= 3) then
            if(txt > "") txt = txt//","
            write(txt,'(a,i0)') trim(txt),iva(rva)
        endif
    end do
    if(txt > "") Write(IOUT,'("   round-off errors for component(s) ",A)') trim(txt)
    txt = ""
    do rva = 1,NVA
        if(iter(iva(rva)) > (ITER_MAX+1)) then
            if(txt > "") txt = txt//","
            write(txt,'(a,i0)') trim(txt),iva(rva)
        endif
    end do
    if(txt > "") Write(IOUT,'("   too many iterations for component(s) ",A)') trim(txt)
EndIf
If(DBG >= 1 .and. errFlags > 0) Then
    ! 1: the numerical solution is uncertain (round-off errors)
    ! 2: too many iterations when solving the mass balance equations
    ! 3: failed to find a satisfactory combination of solids
    ! 4: too many iterations trying to find the solids at equilibrium
    ! 5: some aqueous concentration(s) are too large (20) uncertain activity coefficients
    ! 6: activity factors did not converge.
    ! 7: calculation interrupted by the user.
    Write(IOUT,'(a)') trim(errFlagsGetMessages())
EndIf

if(DBG >= 2) Call printConcs(IOUT)

Return

END SUBROUTINE HaltaCalc

!---- HaltaCancel ---------------------------------------------------------
SUBROUTINE HaltaCancel
    panic = .true.
    Return
END SUBROUTINE HaltaCancel

!----------------------------------------------------------------------------
SUBROUTINE printConcs (IUNIT) ! PUBLIC
! Prints calculated concentrations
! if UNIT > 0 printout to file only
! if UNIT <= 0 printout to terminal only
!
  USE FACTOR_Module, ONLY : Temperature
  Implicit NONE
  Integer, INTENT(IN)  :: IUNIT
  Character (Len=40), PARAMETER :: &
      LINE = '----------------------------------------'
  Integer :: i,j,js, OUT
  Logical :: failed
  Character (len = 20) :: IDENT2
  Real(dp) :: W,W2, pe, Eh

  if(.not.ALLOCATED(C) .or. .not.ALLOCATED(LBeta)) then
    Write(*,'("Memory not allocated for arrays in ''printConcs''")')
    Call ErrStop
  end if

  if(LBeta(1) > 1.d+90) then
    Write(*,'("Arrays not initialized in ''printConcs''")')
    Call ErrStop
  end if

  if(C(1) > 1.d+90) then
    Write(*,'("Calculation not performed yet (in ''printConcs'')")')
    Call ErrStop
  end if

  OUT = IUNIT
  if(OUT <=0) OUT = 6

  Eh = -1234567890.D0

  write(OUT,'(a)') LINE
  failed = .false.
  if(errFlags > 0) then
    failed = isErrFlagSet(2)
    if(.not.failed) failed = isErrFlagSet(3)
    if(.not.failed) failed = isErrFlagSet(4)
    if(.not.failed) failed = isErrFlagSet(6)
  endif
  if(failed) then
    write(OUT,'(a)') '--- HaltaFall - Output composition'//nl//'... Failed calculation: '//nl// &
                trim(errFlagsGetMessages())
  else
    write(OUT,'(a)') '--- HaltaFall - Calculated equilibrium composition:'
    !if(errFlags > 0)    write(OUT,'(a)') trim(errFlagsToString())
    if(isErrFlagSet(1)) write(OUT,'(a)') '    (round-off errors - not within given tolerances)'
  endif
                  !123 12345678901234567890 1234567890123 1234567890123 1234567890123 12345678901 12345678901
  write(OUT,'(a)')'Components:        name    tot conc      solubility    conc          log act     log act-cf'
  Do I = 1, NA
    if(len_trim(IDENTC(I)) <= 0) then
      IDENT2 = '""'
    else
      IDENT2 = IDENTC(I)(1:20)
    endif
    IDENT2 = adjustR(IDENT2)
    W = max(-999.99D0,min(999.99D0,(lnA(I)/ln10)))
    W2 = max(-999.99D0,min(999.99D0,(lnG(I)/ln10)))
    write(OUT,'(i3,1x,a20,1x,1PE13.5,2(1x,E13.5),0P,2(1x,F11.5))') I,IDENT2,TOT(I),SOLUB(I),C(I),W,W2
    IDENT2 = IDENTC(I)(1:20)
    if((IDENT2 == 'e-' .or. IDENT2 == 'e -' .or. &
       IDENT2 == 'E-' .or. IDENT2 == 'E -').and. Temperature > -274.D0) then
        pe = -W
        Eh = 1000.d0 * pe * (8.31446261815324D0 * (Temperature + 273.15D0) * ln10) /  96485.3321233100184D0
    endif
  End Do
                  !123 12345 12345678901234567890 1234567890123 12345678901 12345678901
  write(OUT,'(a)')'Reaction products:        name    conc          log act     log act-cf'
  !NX = MS - NA - MSOL ! nbr of soluble reaction products
  Do I = 1, NX
    J = I + NA
    if(len_trim(IDENT(J)) <= 0) then
      IDENT2 = '""'
    else
      IDENT2 = IDENT(J)(1:20)
    endif
    IDENT2 = adjustR(IDENT2)
    W = max(-999.99D0,min(999.99D0,(lnA(J)/ln10)))
    W2 = max(-999.99D0,min(999.99D0,(lnG(J)/ln10)))
    write(OUT,'(i3,1x,"(",i3,")",1x,a20,1x,1PE13.5,0P,2(1x,F11.5))') J,I,IDENT2,C(J),W,W2
    IDENT2 = IDENT(J)(1:20)
    if((IDENT2 == 'e-' .or. IDENT2 == 'e -' .or. &
       IDENT2 == 'E-' .or. IDENT2 == 'E -').and. Temperature > -274.D0) then
        pe = -W
        Eh = 1000.d0 * pe * (8.31446261815324D0 * (Temperature + 273.15D0) * ln10) /  96485.3321233100184D0
    endif
  End Do
  IF (MSOL > 0) THEN
                    !123 12345 12345678901234567890 1234567890123 12345678901 12345678901
    write(OUT,'(a)')'Solid phases:             name    conc          log act'
    js = 0
    Do I = NX+1, (NX+MSOL)
        J = I + NA
        if(len_trim(IDENT(J)) <= 0) then
            IDENT2 = adjustR('""')
        else
            IDENT2 = adjustR(IDENT(J)(1:20))
        endif
        js = js+1
        W = max(-999.99D0,min(999.99D0,logA(J)))
        write(OUT,'(i3,1x,"(",i3,")",1x,a20,1x,1PE13.5,0P,1x,F11.5)') J,JS,IDENT2,C(J),W
    End Do
  END IF
  Write(OUT,'("Component nbrs. in the order they are iterated:",20(/5X,15I4))') (IVA(j),j=1,NVA+1)
  write(OUT,'("Max. relative tolerance:",1pE10.2)') TOL
  write(OUT,'("Tolerances used (absolute)():")')
  do i = 1,NA
    w = tolY(I)
    if(w > 1.d+100) w = 0.d0
    write(OUT,'(1PE10.3)', Advance='No') w
  enddo
  write(OUT,*)
  if(abs(Eh +1234567890.D0) > 0.1d0) then
    write(OUT,'("Redox potential: Eh =",F9.2," mV,  pe =",F11.5)') Eh, pe
  endif

  IF(activityCoeffsModel >= 0) THEN
    write(OUT,'("Tolerance used for act. coeffs.: ",1pG10.3," (in log-10 scale)")') (tolLnG/ln10)
    if(sumM > 0.01 .and. sumM < 99.) then
        write(IDENT2,'(f12.7)') sumM
    else
        write(IDENT2,'(1PE13.5)') sumM
    endif
    if(jWater > 0) write(OUT,'(A,f12.7,A,f12.7,2A)') "logA[H2O]= ",(lnA(jWater)/ln10),", phi=",phi,", sum(m)=",trim(IDENT2)
    Write(OUT,'(a,1PG11.4,a,G10.3)') "I=",ionicStrCalc,"  electric balance = ",electricBalance
  ENDIF

  write(OUT,'("solid combinations tested: ",i0,"  activity coefficient iterations: ",i0)') iterFasta, iterAc
  if(failed) then
    write(OUT,'(a)') '--- HaltaFall - End of output composition (failed calculation).'
  else
    write(OUT,'(a)') '--- HaltaFall - End of equilibrium composition.'
  endIF
  write(OUT,'(a)') LINE
  RETURN
END SUBROUTINE printConcs

!---- Private Procedures --------------------------------------------------

!---- NOG -----------------------------------------------------------------
SUBROUTINE NOG ! Enough: Calculates logA[], logf[], tot[] and solub[]
Implicit NONE
Integer :: lia,lix,liax,liaf
If(DBG >= 6) Write(IOUT,'("NOG in")')
DO lia=1,NA
    logF(lia)=LnG(lia)/ln10 ! activity coeff.
    SOLUB(lia)=0.D0
    If(lia /= jWater) Then
        if(.not.NOLL(lia)) SOLUB(lia)=C(lia)
        DO lix=1,NX
            liax=NA+lix
            if(.not.NOLL(liax)) SOLUB(lia) = SOLUB(lia) + A(lix,lia)*C(liax)
        End Do
    EndIf ! jWater?
    If(KH(lia) == 2) Then !logA given as input
         ! for water (H2O) the activity might have been changed
         ! when calculating activity coefficients
         If(lia == jWater) Then
            logA(lia) = LnA(LIA)/ln10
            TOT(lia) = 0.D0
         Else
            TOT(lia) = SOLUB(lia)
            If(MSOL > 0) Then
              DO lix=1,MSOL
                liax=NION+lix
                liaf=NX+lix
                if(.not.NOLL(liax)) TOT(lia) = TOT(lia) + A(liaf,lia)*C(liax)
              End Do
            EndIf
         EndIf ! jWater?
    Else ! KH = 1 (Total Concentration given as input)
       logA(LIA)=LnA(LIA)/ln10
    EndIf
End Do

DO lix=1,NX
    liax=NA+lix
    logF(liax)=LnG(LIAX)/ln10
    logA(LIAX)=LnA(LIAX)/ln10
End Do

If(MSOL /= 0) Then
   DO lix=1,MSOL
      liax=NION+lix
      logA(liax) = LnA(liax)/ln10
   End Do
End If

If(DBG >= 6) Write(IOUT,'("NOG returns")')
RETURN
END SUBROUTINE NOG

!0--- INIT ----------------------------------------------------------------
SUBROUTINE HaltaInit
! Checks if for some component the total concentration has changed in a way
! that makes the mass-balance equations to have no solution.
USE FACTOR_Module, ONLY : FactorStart
Implicit NONE
Logical :: OK
Real(dp) :: XM
Integer :: ia,LI,LJ,LIX,LIA,LIAX,LIF,LIAF
SAVE
! VARIABLES -----------------------------------------------------------------
!  Nva= nr of equations for Tot(ia) to be tested (LnA(Ivar) to be varied)
!       in absence of solids
!  Mono(ia)=.True. if component 'ia' only takes part in mononuclear
!           complexes (A(ix,ia)=0 or 1 for all ix)
!  Noll(ia)=.True. if the concentration for species 'ia' is not included
!           in Tot(ia)
!  Iva(m)= ia number for m'th Tot(ia) to be tested in absence of solids
!          (m=1 to Nva)
!  Ivanov(Ivar)= ia to be tested after Ivar, if the mass balance for Tot(Ivar)
!                is not satisfied
!  Ivabra(Ivar)= Ivar to be tested after Ivar if the mass balance for
!                Tot(Ivar) is satisfied.
!----------------------------------------------------------------------------
If(DBG >= 3) Write(IOUT,'("haltaInit() in,  CONT=",L1)') CONT

IF(firstTime) THEN
    If(DBG >= 3) Write(IOUT,'("haltaInit: firstTime, setting CONT=false")')
    if (.not.ALLOCATED(C) .or. .not.ALLOCATED(LBeta)) then
      Write(*,'("HaltaCalc: Memory not allocated for arrays")')
      Call ErrStop
    end if
    if (LBeta(1) > 1.d+90) then
      Write(*,'("HaltaCalc: Arrays not initialized")')
      Call ErrStop
    end if
    NFSPAR = 0
    NVAF = 0
    CONT = .false.
    firstTime = .false.

    ! NYKO

    ! The first time HaltaFall is called, a plan is made to solve the mass
    ! balance equations, some variables are initialized, etc.

    If(MSOL > 0)  Then
        DO lif=1,MSOL
            liaf=NX+lif
            FScal(lif)=1.D0
            liax=NA+NX+lif
            If(NOLL(liax)) Cycle
            ! The solids get a scaling factor FSCAL to determine which
            ! of the oversaturated solids will be allowed to precipitate
            Do  lia=1,NA
                FScal(lif) = FScal(lif) + abs(A(liaf,lia))
            End Do
        END DO
    End If

    DO lia=1,NA
      ITER(LIA)=0
      catchRoundingErrors(LIA)=0
      jumpChord(LIA) = .true.
      MONO(lia)=.true.
      POS(lia)=.false.
      if(.not.NOLL(lia)) POS(lia)=.true.
      NEG(lia)=.false.
      do lix=1,NX
         liax=lix+NA
         if(NOLL(liax)) cycle
         if(abs(A(lix,lia)) > 1.D-5 &
            .and. abs(A(lix,lia)-1.D0) > 1.D-5) MONO(lia)=.false.
         if(A(lix,lia) > 0.D0) POS(lia)=.true.
         if(A(lix,lia) < 0.D0) NEG(lia)=.true.
      end do
    END DO

    DO li=1,NA
        DO lj=1,NA
            OBER(li,lj) = (li /= lj)
            IF(OBER(li,lj)) THEN
                do lix=1,NX
                    if(NOLL(NA+lix)) cycle
                    XM = A(lix,li)*A(lix,lj)
                    if(abs(XM) > 1.D-8) OBER(li,lj)=.false.
                end do
                if(MSOL == 0)  cycle
                do lif=1,MSOL
                    if(NOLL(NION+lif)) cycle
                    liaf=NX+lif
                    XM = A(liaf,li)*A(liaf,lj)
                    if(abs(XM) > 1.D-8) OBER(li,lj)=.false.
                end do
            ENDIF ! OBER(li,lj)
        END DO
    END DO

    DO li=1,NA
      NOBER(li)=0
      Do lj=1,NA
         If(OBER(li,lj))  NOBER(li)=NOBER(li)+1
      End Do
    END DO

END IF ! --- firstTime ---

If(DBG >= 3) call printInput()

Do lix=1,NX
    lnBeta(lix) = ln10*lBeta(lix)
End Do

IF(MSOL > 0)  THEN
   DO lif=1,MSOL
      liaf=NX+lif
      lnKF(lif)=-lBeta(liaf)*ln10
   END DO
END IF

actCoefCalc = ((abs(ionicStr)>1.D-15) .and. (activityCoeffsModel >= 0 .and. activityCoeffsModel <= 2))

call FactorStart(NION,C,LnG)

Do lix =1, NION
    if(.not.CONT) then
        lnG(lix) = 0.d0;
    else
        lnG(lix) = logf(lix) * ln10
    endIf
    oldLnG(lix) = lnG(lix);
    do li = 1,4
        largerOldLnG(li,lix) = .false.
    enddo
EndDo

IF(.not. CONT) THEN ! not a continuation run
  ok = .false.
  DO ia=1,NA
      noCalc(ia)=.false.  ! that is, calc = true, calculation is possible
      if(KH(ia) == 2) then  ! logA given as input, calculation not requested
        if(ia == jWater) logA(ia) = 0.d0;
        cycle
      endif
      logA(ia) = -10.d0
      if(TOT(ia) > 0.D0) logA(ia) = LOG10(TOT(ia))-2.D0
      if(POS(ia) .and. NEG(ia)) cycle
      if(POS(ia) .and. TOT(ia) > 0.D0) cycle
      if(NEG(ia) .and. TOT(ia) < 0.D0) cycle
      noCalc(ia)=.true. ! calculation is not possible
      logA(ia)= NOCALC_LOGA
  End Do
ELSE ! a continuation run
   ok = .true.
   DO ia=1,NA
    iter(ia) = 0
    catchRoundingErrors(ia) = 0
    jumpChord(ia) = .true.
    xOld(ia) = huge(1.d0)
    if(KH(ia) == 2) cycle
    if(POS(ia).AND.NEG(ia)) cycle
    if((POS(ia) .and. TOT(ia) > 0.D0) .or. &
       (NEG(ia) .and. TOT(ia) < 0.D0)) then
            if(noCalc(ia)) then !calculation was not possible
              ok = .false.
              noCalc(ia)=.false. ! calculation is possible
              logA(ia)=-10.d0
              if(TOT(ia) > 0.D0) logA(ia) = LOG10(TOT(ia))-2.D0
              !LnA(ia) = logA(ia) * ln10
            endif
    else
        if(.not.noCalc(ia)) ok = .false. !calculation was possible
        noCalc(ia) = .true. ! calculation is not possible
        logA(ia) = NOCALC_LOGA
        !LnA(ia) = NOCALC_LOGA * ln10
    endif
   END DO
END IF

IF(.not. ok) THEN
    If(CONT) Then
        CONT = .false.
        if(DBG >= 3) Write(IOUT,'("setting cont = false")')
    EndIf
    Call HaltaGetIva
ENDIF

! debug printing
If(DBG >= 3) Then
    Write(IOUT,'(" mono()=",20L4,20(:/,7x,20L4))') (mono(ia), ia=1,NA)
    Write(IOUT,'("nober()=",20I4,20(:/,8x,20I4))') (nober(ia), ia=1,NA)
    Do lia = 1, NA
        Write(IOUT,'("ober(",I2,")()=",20L4,20(:/,8x,20L4))') lia,(ober(lia,ia), ia=1,NA)
    End Do
    Call printArrays (.true., .false.) ! print iva, ivaBra, ivaNov
    Write(IOUT,'("noCalc()=",20L4,20(:/,9x,20L4))') (noCalc(ia), ia=1,NA)
    Write(IOUT,'("    kh()=",20I4,20(:/,9x,20I4))') (kh(ia), ia=1,Na)
    Write(IOUT,'("   Tot()=",1p,10(1x,g13.5),20(:/,9x,10(1x,g13.5)))') (tot(ia), ia=1,Na)
    Write(IOUT,'("  logA()=")', Advance="no")
    Do ia=1,Na
        if(logA(ia) > -99999.) then
            Write(IOUT,'(1p,1x,g22.15)',Advance="no") logA(ia)
        else
            Write(IOUT,'(" -99999")',Advance="no")
        endif
    endDo
    Write(IOUT,*)
    Write(IOUT,'("haltaInit() returns")')
EndIf

RETURN

END SUBROUTINE HaltaInit

!---- HaltaGetIva --------------------------------------------------------
SUBROUTINE HaltaGetIva
! Makes a plan made to solve the mass balance equations,
! initializes some variables, etc.
Implicit NONE
Integer :: i, ia, IVA1, IVA2, jf, li, M
If(DBG >= 3) Write(IOUT,'("haltaGetIva() in at NYA")')

! NYA

NFALL=0
IF(MSOL > 0) THEN
    Do ia=1,NA
      BER(ia)=.false.
      FALLA(ia)=.false.
    End Do
    Do jf=1,MSOL
      FALL(jf)=.false.
    End Do
END IF

NVA=0
Do ia = 1,NA
    if(KH(ia) /= 2 .and. .not.noCalc(ia)) then  ! calculation requested and possible
        NVA = NVA +1
        IVA(NVA) = ia
    else
        ivaBra(ia) = -1
        ivaNov(ia) = -1
    end if
End Do
Do li = (NVA+1), (NA+1)
    IVA(li) = -1
EndDo

! PLAN

If(DBG >= 3) Write(IOUT,'("haltaGetIva() at PLAN")')
IF(NVA > 1) THEN
  DO WHILE (.true.)
    M=-1
    DO li=1,(NVA-1)
      IVA1=IVA(li)
      IVA2=IVA(li+1)
      If(NOBER(IVA2) > NOBER(IVA1)) Then
          M=IVA(li)
          IVA(li)=IVA(li+1)
          IVA(li+1)=M
      Else
          if(MONO(IVA2) .and. .not.MONO(IVA1) .and. &
                 NOBER(IVA2)==NOBER(IVA1)) then
              M=IVA(li)
              IVA(li)=IVA(li+1)
              IVA(li+1)=M
          endif
      End If
    END DO
    IF (M == -1) EXIT
  END DO ! while .true.
END IF ! NVA > 1

IF (NVA > 0) THEN
    DO li=1,NVA
        IVAR=IVA(li)
        IVABRA(IVAR)=IVA(li+1)
        i=0
! HOPP ("jump" in Swedish)
        Do While (.true.)
          i=i+1
          If (.not. OBER(IVAR,IVA(i))) Then 
            IVANOV(IVAR)=IVA(i)
            Exit
          EndIf
        End Do !While (true)
    END DO ! li
END IF
IVABRA(NA+1) = -1

If(DBG >= 3) Write(IOUT,'("haltaGetIva() returns")')
RETURN
END SUBROUTINE HaltaGetIva

!1--- KILLE --------------------------------------------------------------
SUBROUTINE Kille
! Used to solve the mass balance equation of component: Ivar, which forms
! polynuclear complexes.  Locates the solution by setting
! steg[ivar]=0.5,1,2,4,8,16,32,64,...  until a pair of x1 and x2 is found
! (when karl[ivar]=2 or 3).  Then decreases steg[ivar] (=steg[ivar]/2,
! when karl[ivar]=4) until  x1 and x2 are separated less than STEGBYT; at
! that point a chord method is used.
! VARIABLES----------------------------------------------------------------
!  X= independent variable in the equation Y0=Y(X)
!  Y= dependent variable
!  Y0= value aimed at in Y0=Y(X)
!  X1(ia) and X2(ia)= values for X below and above the right value
!  Y1(ia) and Y2(ia)= values for Y corresponding to X1 and X2
!--------------------------------------------------------------------------
Implicit NONE
Logical :: YL
Real(dp) :: W,W1
Integer :: LQ
!  Locate the solution by giving Steg(Ivar)=0.5,1,2,4,8,16,32,64,...
!    until a pair of X1 and X2 is found (in Karl(Ivar)=2 or 3). Then
!    decrease Steg(Ivar) (=Steg(Ivar)*0.5, in Karl(Ivar)=4) until
!    X1 and X2 are separated less than STEGBYT
   YL = (Y > Y0)
   If(YL) Then
        X2(IVAR)=X
        Y2(IVAR)=Y
   Else
      X1(IVAR)=X
      Y1(IVAR)=Y
   EndIf
   If(DBG >= 5) Then
        Write(IOUT,'(3x,"kille() in; ivar=",i0,", karl=",i0,", x=",1P,G20.12,", steg=",G11.3,", STEGBYT=",0P,F6.4)') &
                        ivar,karl(ivar),X,steg(ivar),STEGBYT
        Write(IOUT,'(7x,"x1(",i0,")=",1PG20.12,", x2(",i0,")=",G20.12,",    y=",G14.6,", y0=",G14.6)') &
                        ivar,x1(ivar),ivar,x2(ivar),y,y0
        Write(IOUT,'(7x,"y1(",i0,")=",1P,G14.6,", y2(",i0,")=",G14.6)') &
                        ivar,y1(ivar),ivar,y2(ivar)
   EndIf
   LQ=KARL(IVAR)
   SELECT CASE ( LQ )
      CASE (    1)
         GO TO 1301
      CASE (    2)
         GO TO 1302
      CASE (    3)
         GO TO 1303
      CASE (    4)
         GO TO 1304
      END SELECT
1301 If(YL) Then ! Karl=1  the beginning
        KARL(IVAR)=3
        X=X-STEG(IVAR)
     Else
        KARL(IVAR)=2
        X=X+STEG(IVAR)
     EndIf
     GO TO 9000
1302 If(.not.YL) Then ! Karl=2   it was Y<Y0
        STEG(IVAR)=STEG(IVAR)+STEG(IVAR)
        X=X+STEG(IVAR)
     Else
        KARL(IVAR)=4
     EndIf
     GO TO 9000
1303 If(YL) Then ! Karl=3   it was Y>Y0
        STEG(IVAR)=STEG(IVAR)+STEG(IVAR)
        X=X-STEG(IVAR)
     Else
        KARL(IVAR)=4
     EndIf
     GO TO 9000
1304 Continue ! Karl=4   There is both X1 and X2 corresponding to Y<Y0 and Y>Y0
     If(STEG(IVAR) >= STEGBYT) Then
        STEG(IVAR) = 0.5D0 * STEG(IVAR)
        if(YL) then
            X=X-STEG(IVAR)
        else
            X=X+STEG(IVAR)
        endif
        GO TO 9000
     EndIf
!Korda: Aproximating the solution by chord shooting ('secant method')
     if(iter(IVAR) < ITER_MAX/6 .or. .not.jumpChord(IVAR) .or. &
        abs(1.d0-(x2(ivar)/x1(ivar))) < 1.d-9 ) then
        W=Y0-Y1(IVAR)
        W1=X2(IVAR)-X1(IVAR)
        X=X1(IVAR)+ W*W1/(Y2(IVAR)-Y1(IVAR))
        If(DBG >= 5) Write(IOUT,'(7x,"Korda: new x=",1PG20.12,"  delta-x=",G22.14)') x, (x-xOld(ivar))
        jumpChord(IVAR) = .true.
     else ! iter[ivar] >= ITER_MAX/6 and jumpChord[ivar]
        ! --- Sometimes the chord shooting progresses very slowly
        !     (depending on the inner loops)
        X = (X2(IVAR)+X1(IVAR))/2.d0
        If(DBG >= 5) Write(IOUT,'(6x,"---- iter[",i0,"]=",i0,", new x=",1PG20.12,"  delta-x=",G22.14)') &
                                        ivar,iter(ivar),x, (x-xOld(ivar))
        jumpChord(IVAR) = .false.
     endif
     !--- Avoid rounding errors
     If(abs(1.d0-(xOld(ivar)/x)) < 1.d-10) Then ! could perhaps be 1.d-12 or 1.d-10
        If(catchRoundingErrors(IVAR) == 0) Then
            X = X + 1.d-10*abs(X)  ! could perhaps be 1.d-9
        Else If(catchRoundingErrors(IVAR) == 1) Then
            X = X - 2.d-10*abs(X)  ! could perhaps be 2.d-9
        End If
        If(DBG >= 5) Write(IOUT,'(6x"---- catchRoundingErrors[",i0,"]=",i0,", x=",1PG21.13,"  old x=",G21.13)') &
                                                                    ivar,catchRoundingErrors(ivar),x,xOld(ivar)
        !catchRoundingErrors[] may be 0,1,2 or 3
        If(catchRoundingErrors(IVAR) < 3) Then
            catchRoundingErrors(IVAR) = catchRoundingErrors(IVAR) +1
        End If
     End If

9000 Continue
xOld(ivar) = x
If(DBG >= 5) Write(IOUT,'(3x"kille() returns, karl(",i0,")=",i0,"; indik=",i0,"; steg(",i0,")=",0PF7.4," x=",1PG21.13)') &
                ivar,karl(ivar),indik,ivar,STEG(IVAR),X
RETURN
END SUBROUTINE Kille

!2--- FASTA --------------------------------------------------------------
SUBROUTINE Fasta
! Finds out what solids are present in the system at equilibrium.
! Returns INDIK:
! = 2 if a new set of solids is chosen
! = 3 if ok (either no solids precipitate in the system, or no solubility
!   products are exceeded and no solid assumed present has zero or negative
!   concentration)
! = 4 if too may solid phase combinations tested
! = 1 if "tji" (all combinations of solids give unsolvable equations, i.e.
!   negative determinants).
Implicit NONE
! VARIABLES -----------------------------------------------
!  Ber(ia)=.True. if LnA(ia) is calculated by means of the solubility product
!          in subroutine LnaBer
!  Fall(if)=.True. if the solid 'if' is assumed present at equilibrium
!  Falla(ia)=.True. if the component 'ia' is assumed to occur in one or more
!            solids at equilibrium
!  Fut(i)= Ifspar number of i'th solid to be eliminated at some stage in
!          systematic variation during singFall=.True.
!  Ibe(m)= 'ia' number for m'th LnA(ia) to be calculated with the solubility
!         product in LnaBer (m=1 to Nfall)
!  Iber(j)= j'th of 'Iva' numbers to be an 'Ibe' (Ibe(j)=Iva(Iber(j)), (j=1 to
!           Nfall)
!  Ifall(m)= 'if' number of m'th solid present at equilibrium (m=1 to Nfall)
!  Ifspar(m)= if number of m'th solid indicated as possible a FALLPROV
!  Iva(m)= 'ia' number for m'th component to be tested in absence of solids
!         (m=1 to Nva)
!  Ivaf(m)='ia' number for m'th component to be tested in presence of solids
!          (m=1 to Nvaf)
!  Nfall= nr of solids present at equilibrium
!  Nfspar= nr of solids indicated at FALLPROV
!  Nvaf= nr of mass balance equations to be solved in presence of solids
!  Nut= nr of solids systematically eliminated at some stage while singFalll is
!       .True.
!  singFall= .True. if the first matrix Ruta tried has come out singular
! ---------------------------------------------------------
! Changes by I.Puigdomenech in Jan.2011:
! Curiously, the FUTT-HOPPFUT-INFUT-UPPNUT system published in 1967
! does not work when nfSpar=1, because setting nUt=1 means that all
! solids (i.e. the only one) are "striked out". The new changes may be
! found at:
! UPPNUT (if(nfall == 0 .and. nfSpar == 1 .and. nUt == 1)) and at
! INFUT (if(nfall == 1 .and. nfSpar == 1 .and. nUt == 1)).
! ---------------------------------------------------------
! ---------------------------------------------------------
! Master routine for solid phase calculations. The routine cycles through
! chunks of code in the original FASTA routine, using the value of nextFall
! as switch. The meanings of nextFall on return to this routine
! are as follows:
!   nextFall     next routine       line No. in original programme
!       0        Exit to PROV,          9000, 10000, 25000
!                "indik" set ready
!       1        fallProv_InFall        14000
!       2        beFall                 16000
!       3        anFall                 24000
!       4        utFall                 17000
!       5        sing                   18000
!       6        fUtt                   20000
!       7        hoppFut                21000
!       8        inFut                  22000
!       9        uppNut                 23000

nextFall = 1

DO WHILE (nextFall > 0)
    if(panic) Exit

    SELECT CASE (nextFall)
        CASE(1)
            ! Tests for changes in the number of solids calculated with the given lnA[] values
            call fallProv_Infall
        CASE(2)
            ! Initializes iber[]
            call beFall
        CASE(3)
            ! Calculates ibe[], ivaf[], rut1[][] and pva[][] arrays
            call anFall
        CASE(4)
            ! Calculates ruta[][] and inverts the array
            call utFall
        CASE(5)
            call sing_HoppSi
        CASE(6)
            call fUtt
        CASE(7)
            call hoppFut
        CASE(8)
            call inFut
        CASE(9)
            call uppNut
    END SELECT

END DO

RETURN
END SUBROUTINE Fasta

SUBROUTINE fallProv_Infall
! Tests for changes in the number of solids
! calculated with the given lnA[] values.
Implicit NONE
Logical :: BRA, foundOne
Integer :: ia, li, lia, liaf, liax, lif, lix, lj, lq, lqa
Integer :: NyFall, kMax
Real(dp) :: W, zMax
! FALLPROV  Check that the solubility product is not exeeded for any solid
!            phase assumed to be absent. Find Ifall() and Fall() for the new
!            solids appearing
If(DBG >= 4) Then
    Write(IOUT,'("--- Fasta() in at FallProv, nfall =",i0)') NFALL
    ! print ifall
    call printArraysFasta(.true.,.false.,.true.,.false.,.false.,.false.,.false.,.false.)
    Write(IOUT,'("Testing that no additional solid is supersaturated and that the formal",/, &
                "concentration is not negative for any solid assumed to be present.")')
EndIf
BRA=.TRUE.
NyFall=0
! zMax and kMax are used with ONLY_ONE_SOLID_AT_A_TIME
! to find which of the supersaturated solids is to be included
! in the equilibrium system
zMax=-1.D50; kMax=-1
DO lif=1,MSOL
    liax=NION+lif
    If(.not.FALL(lif)) Then
      liaf=NX+lif
      C(LIAX)=0.D0
      W=0.D0
      DO lia=1,NA
        W = W + A(LIAF,LIA)*LnA(LIA)
      End Do
      LnA(LIAX) = W - LNKF(LIF) ! lnA is now the (over)saturation index
      If(W <= LNKF(LIF) .or. NOLL(LIAX) .or. NVA == 0) Cycle
      ! ---- Block added 2013-Jan.
      !      Solids may be oversaturated that can not be allowed to
      !      precipitate. For example:
      !      let us say C(s) is oversaturated at low redox potential,
      !      if HCO3- is the component and its total conc. is given,
      !      then C(s) may precipitate;  but if CO2(g) is the component
      !      and the partial pressure of CO2(g) is given, then C(s)
      !      may NOT be present at equilibrium.
      foundOne = .false.
      ! loop through components for which the total conc. is given
      Do lia=1,NVA
        ! is the stoichiometric coefficient non-zero?
        if(abs(A(liaf,IVA(lia))) >0.00001D0) then
            foundOne = .true.; Exit
        endif
      End Do
      if(.not. foundOne) Cycle !all coefficients zero?
      ! calculate the scaled oversaturation
      W=LnA(liax)/FSCAL(lif)
      If(w < tolFasta) Then
        if(DBG >= 4) &
            Write(IOUT,'("FallProv solid nbr: ",i0,", lnA[",i0,"]=",1PG14.6,'// &
              '", scaling factor = ",G14.6,", lnA/fscal = ",G14.6,/,'// &
              '"       tolerance = ",1PE13.6,"; Accepting the oversaturation of solid ",i0," as zero within the tolerance...")') &
              lif, liax,lnA(liax), fscal(lif), w, tolFasta, lif
        Cycle
      EndIf
      ! ---- Block end
      If(ONLY_ONE_SOLID_AT_A_TIME) Then
        If(W > zMax) Then
            zMax=W
            kMax=lif
        EndIf
      Else
        BRA=.false.
        Fall(lif)=.true.
        NyFall=NyFall+1
        ifall(Nfall+NyFall)=lif
      End If !ONLY_ONE_SOLID_AT_A_TIME
    Else ! if FALL(lif)
      LnA(LIAX)=0.D0
    EndIf ! FALL(lif)?
End Do

If(ONLY_ONE_SOLID_AT_A_TIME .and. kMax > 0) Then
    BRA = .false.
    FALL(kMax)=.TRUE.
    NyFall=1
    IFALL(NFALL+NyFall)=kMax
End If !ONLY_ONE_SOLID_AT_A_TIME

if(.not.BRA .and. DBG >= 4) then
    Write(IOUT,'("FallProv: NyFall=",i2)', Advance='No') NyFall
    if(ONLY_ONE_SOLID_AT_A_TIME) then
        Write(IOUT,'(", new solid nbr. is: ",i0,", lnA[",i0,"]=",1PG14.6)') &
              ifall(nfall+nyfall), (nIon+ifall(nfall+nyfall)), lnA((nIon+ifall(nfall+nyfall)))
    else
        Write(IOUT,*)
        call printArraysFasta(.true.,.false.,.false.,.false.,.false.,.false.,.false.,.false.) ! print IFALL()
    endif
endif
! ---- Block added 2018-March.
! When singFall is true, solids are removed one at a time until
! a set is found that gives a non-singular matrix ruta. When
! testing a reduced set of solids, if a new solid becomes
! supersaturated which was not in the list "ifspar", then the
! new solid is tested by setting singFall = false
If(ONLY_ONE_SOLID_AT_A_TIME .and. kMax > 0) Then
    If(singFall) Then
        foundOne = .false.
        Do lif = 1, NFALL
            if(ifSpar(lif) == kMax) then
                foundOne = .true.
                Exit
            endif
        EndDo
        if(.not.foundOne) then
            if(DBG >= 4) Write(IOUT,'("The new solid is not among the ""ifSpar"" list.  Setting singFall=false.")')
            singFall = .false.
        endif
    EndIf
End If !ONLY_ONE_SOLID_AT_A_TIME and kMax >0
! ---- Block end

If(NFALL+NyFall >0 .and. NFALL+NyFall < NA) Then
    Do lif = NFALL+NyFall+1, NA
        IFALL(lif) = 0
    EndDo
EndIf

If(NVA == 0) THEN
    If(DBG >= 4) Write(IOUT,'("FASTA returns (nva=0)")')
    nextFall = 0
    Return
END IF
! Check that the quantity of solid is not negative for any solid
! assumed to be present.
IF(NFALL > 0) THEN
    DO li=1,NFALL
      ia=IBE(li)
      W=TOT(ia)-C(ia)
      DO lix=1,NX
         liax=NA+lix
         W = W - A(lix,ia)*C(liax)
      End Do
      TOTBE(li)=W
    End Do
    DO li=1,NFALL
      W=0.D0
      DO lj=1,NFALL
         W=W + RUTA(li,lj)*TOTBE(lj)
      End Do
      lq=IFALL(li)
      lqa=NION+lq
      C(lqa)=W
      If(C(lqa) >= 0.D0)  Cycle
      If(DBG >= 4) Write(IOUT,'("Fasta(): for solid ",i0," the concentration is <0, C[",i0,"]=",1PE13.5,", tolerance=",E13.5)') &
                        lq, lqa, w, tolFasta
      If(-C(lqa) > tolFasta) Then
        FALL(lq)=.false.
        C(lqa)=0.D0
        BRA=.false.
      Else
        If(DBG >= 4) Write(IOUT,'("         accepting the concentration of solid ",i0," as zero within the tolerance...")') lq
      EndIf
    End Do
END IF ! NFALL >0

IF(BRA) THEN ! it is OK
    INDIK = 3
    If(DBG >= 4) Then
        Write(IOUT,'("HaltaFall.fasta() returns OK (to nog()); indik=3; nfall=",i0,"; iterFasta=",i0)') nfall, iterFasta
        ! print most arrays
        Call printArraysFasta(.true.,.true.,.true.,.true.,.true.,.false.,.false.,.false.)
    EndIf
    nextFall = 0
    RETURN
ENDIF

! not OK
If(iterFasta > ITER_FASTA_MAX) Then
    INDIK = 4
    call errFlagSet(4) !too many iterations trying to find the solids at equilibrium
    If(DBG >= 4) Write(IOUT,561) ITER_FASTA_MAX
    561 Format("Error in HaltaFall.fasta(): ",i0," different solid phase combinations",/, &
               "were tested and found NOT satisfactory (too many iterations).")
    nextFall = 0
    RETURN
End If

iterFasta = iterFasta +1

! ---------------------------------------------------------
! INFALL:  find new: nfall, ifall[nfall]
! ---------------------------------------------------------
! The LnA(ia) were not consistent with the solid phases. Either NyFall new
! solid phases appeared, or some solids assumed present had negative C(if).
! At first it is assumed at BEFALL that the first Nfall of the Iva(ia) are
! the Ibe() to be calculated. If the determinant of Ruta is found to be zero
! at UTFALL, the Ibe are changed systematically at SING-HOPPSI by means of
! array Iber, until a non-zero determinant is found.
! Then at ANFALL, the arrays Rut1 and Pva are calculated, and the new mass
! balance equations are solved in SLINGOR (routine HALTA).

NFALL=NFALL+NyFall
li = 0
DO While (li < NFALL)
    li = li+1
    if(FALL(IFALL(li))) Cycle
    NFALL = NFALL-1
    If(NFALL /= 0) Then
        DO lj=li,NFALL
            IFALL(lj)=IFALL(lj+1)
        End Do
    EndIf
    li=li-1
EndDo
If(DBG >= 4) Then
    Write(IOUT,'("Fasta at InFall; nfall=",i0,", NyFall=",i0,", nva=",i0,", singFall=",l1)') nfall,NyFall,nva,singFall
    call printArraysFasta(.true.,.false.,.false.,.false.,.false.,.false.,.false.,.false.) !print ifall
    if(singFall) call printArraysFasta(.false.,.false.,.false.,.false.,.false.,.false.,.true.,.false.) !print ifSpar
EndIf

IF(NFALL <= 0) THEN
    nextFall = 3 ! no solid phases, goto ANFALL;
    Return
ENDIF

! NFALL >0
!Check if singFall=.False., and if Nfall <= Nva
If(.not.singFall) Then
    DO li=1,NFALL
         IFSPAR(li)=IFALL(li)
    End Do
    NFSPAR=NFALL
    NUT=0
    If(NFALL > NVA) Then
        NUT = NFALL - NVA -1
        If(DBG >= 4) Then
            Write(IOUT,'("    nfall (=",i0,") is > nva (=",i0,") & not singFall",/,"    nfSpar=",i0)') nfall,nva,nfSpar
            call printArraysFasta(.false.,.false.,.false.,.false.,.false.,.false.,.true.,.false.)
            Write(IOUT,'("   setting singFall=true")')
        EndIf
        singFall = .true.
        nextFall = 9 ! goto UppNut
        Return
    End If ! nfall > nva
End If !not singFalll

If(NFALL >= NVA .and. singFall) Then
    If(DBG >= 4) Write(IOUT,'("    nfall (=",i0,") is >= nva (=",i0,") & singFall = "l1,/,"    nfSpar=",i0," nUt=",i0)') &
                        nfall,nva,singFall,nfSpar,nut
    NFALL = NFSPAR - NUT
    Do lj=1,MSOL  ! added April-2013
        FALL(lj) = .false.
        Do li=1,NFSPAR
            If(IFSPAR(li) == lj) Then
               FALL(lj) = .true.
               Exit
            End If
        End Do
    End Do
    nextFall = 7 ! goto HoppFut
    RETURN
End If

nextFall = 2 ! solids present: goto BEFALL then UTFALL
RETURN
END SUBROUTINE fallProv_Infall

SUBROUTINE beFall
! Initializes iber()
Implicit NONE
Integer :: li
If(DBG >= 4) Write(IOUT,'("Fasta() at BeFall, nfall=",i0)') nfall

Do li = 1, NFALL
    IBER(li) = li
EndDo
IBER(NFALL+1)=0

nextFall = 4 ! goto UTFALL

RETURN
END SUBROUTINE beFall

SUBROUTINE anFall
! Calculates ibe[], ivaf[], rut1[][] and pva[][] arrays.
Implicit NONE
Integer :: lia, li, liaf, lix, lj, lm, lq, lz
Real(dp) :: W
! ---------------------------------------------------------
! ANFALL:
! A set of nfall solids has been found with a non-zero determinant
! for the array "ruta". Some variables are changed, and the arrays
! rut1 and pva are calculated.
If(DBG >= 4) Write(IOUT,'("Fasta() at AnFall, nfall=",i0)') nfall
DO lia=1,NA
    FALLA(lia)=.false.
    BER(lia)=.false.
End Do

IF(NFALL > 0) THEN
    ! Get Falla(ia), Ber(ia) and Ibe(ia) (from Ifall and Fall)
    Do li=1,NFALL
        liaf=NX+IFALL(li)
        If(FALL(IFALL(li))) Then
            Do lia=1,NA
                If(abs(A(liaf,lia)) > 1.D-10)  FALLA(lia)=.true.
            EndDo
        End If
        IBE(li)=IVA(IBER(li))
        BER(IBE(li))=.true.
    EndDo
    ! Find Nvaf and Ivaf()
    NVAF=0
    If(NVA > 0) Then
      DO li=1,NVA
        lia = IVA(li)
        If(BER(lia))  Cycle
        NVAF=NVAF+1
        IVAF(NVAF) = lia
      End Do
      If(DBG >= 4) Then
        Write(IOUT,'("nvaf=",i0)') nvaf
        ! print arrays "ber" and "ivaf"
        Call printArraysFasta(.false.,.false.,.false.,.false.,.true.,.true.,.false.,.false.)
      EndIf

      ! Calculate RUT1(I,J)
      If(NVAF > 0) Then
        DO li=1,NVAF
          DO lj=1,NFALL
            W=0.D0
            DO lm=1,NFALL
               lq=IVAF(li)
               lz=NX+IFALL(lm)
               W = W + A(lz,lq)*RUTA(lm,lj)
            End Do
            RUT1(li,lj)=W
          End Do
        End Do
        ! Calculate PVA(ia,IX)
        DO li=1,NVAF
          DO lix=1,NX
            lia=IVAF(li)
            W=A(lix,lia)
            If(NFALL > 0)  Then
              DO lm=1,NFALL
                lq=IBE(lm)
                W = W - RUT1(li,lm)*A(lix,lq)
              EndDo
            EndIf !NFALL >0
            PVA(lia,lix)=W
          EndDo !lix
        EndDo !li
      EndIf ! NVAF >0
    EndIf ! NVA >0
ENDIF ! NFALL >0

! SLINGOR(IN)

indik = 2

If(DBG >= 4) Then
    Write(IOUT,'("HaltaFall.fasta() returns to Slingor(In); indik =2, nfall=",i0,"; iterFasta=",i0)') nfall,iterFasta
    ! print most arrays
    Call printArraysFasta(.true.,.true.,.true.,.true.,.true., .false.,.false.,.false.)
EndIf

nextFall = 0 ! exit to Slingor
RETURN
END SUBROUTINE anFall

SUBROUTINE utFall
! Calculates ruta[][] and inverts the array.
Implicit NONE
Integer :: ia, li, lj, lf
If(DBG >= 4) Then
    Write(IOUT,'("Fasta() at UtFall, nfall=",i0)') nfall
    ! print array "ifall" and "iber"
    Call printArraysFasta(.true., .true., .false.,.false.,.false.,.false.,.false.,.false.)
EndIf
DO li=1,NFALL
    ia=IVA(IBER(li))
    DO lj=1,NFALL
        lf=NX+IFALL(lj)
        RUTA(li,lj)=A(lf,ia)
    End Do
End Do
If(DBG >= 4 .and. NFALL>=2) Then
    Do Li=1,NFALL
        Write(IOUT,'(3x"RUTA(",I0,",J)= ",1P,100(E12.3))') LI,(ruta(LI,LJ),LJ=1,NFALL)
    EndDo
EndIf
INDIK=0
Call Invert ! sets indik =0 (ok) or =1 (matrix singular)
If(INDIK /= 0) Then
    nextFall = 5 ! matrix singular, goto SING:
Else
    nextFall = 3 ! matrix inverted, goto Anfall
EndIf
RETURN
END SUBROUTINE utFall

SUBROUTINE sing_HoppSi
! the matrix ruta(,) was singular
Implicit NONE
! --------
! SING:
! --------
If(DBG >= 4) Write(IOUT,'("Fasta() at Sing-HoppSi")')
! --------
! HOPPSI:
! --------
! Bump up iber[] indices  one at a time to give all combinations of nfall
! "be" components from the set of nva
If(hopp(NFALL+1,IBER,NVA)) Then
    !all iber[] sets used up: move to FUTT to remove a solid phase
    If(DBG >= 4) Write(IOUT,'("     ruta[][] was singular for all possible iber[] '//'combinations; setting singFall = true")')
    singFall = .true.
    ! note that nfspar and ifspar[] are already set at INFALL
    nextFall = 6 ! goto FUTT
Else
    ! move to UTFALL to try the new iber() set
    nextFall = 4 ! goto UTFALL
EndIf
If(DBG >= 4) Write(IOUT,'("     Sing-HoppSi returns")')
RETURN
END SUBROUTINE sing_HoppSi

SUBROUTINE fUtt
! Routine to handle systematic removal of solid phases in cases where no
! consistent set can be found.
! ---------------------------------------------------------
! FUTT:
! ---------------------------------------------------------
! iber[] sets are exausted. Either:
!   - It has been proved inpossible to pick out a group of nfall
!     components iber[] such that the array "ruta" has a non zero
!     determinant,
!   - Or more solids are indicated than there are logA[] values
!     to vary (nfall > nva).
! One must try systematically combinations of a smaller number
! (nfall-nUt) of solid phases, until a non-zero determinant for
! "ruta" is found. This is done at labels FUTT, INFUT and UPPNUT
! using the array fut[]. To avoid going into a loop, the program
! sets singFall=true and remembers the first set of solids
! (ifspar[nfspar]) which gave only zero determinants.
! ---------------------------------------------------------
! Changes by I.Puigdomenech in 2011-Jan.:
! Curiously, the FUTT-HOPPFUT-INFUT-UPPNUT system published in 1967
! does not work when nfSpar=1, because setting nUt=1 means that all
! solids (i.e. the only one) are "striked out".
! The changes are at UPPNUT and at INFUT.
! ---------------------------------------------------------
Implicit NONE
If(DBG >= 4) Then
    Write(IOUT,'("Fasta() at Futt;  nUt=",i0,", nfall=",i0,", nfspar=",i0)') nUt,nfall,nfSpar
    ! print arrays "ifspar" and "fut"
    Call printArraysFasta(.false.,.false.,.false.,.false.,.false.,.false., .true., .true.)
EndIf
If(nUt == 0) Then
    nextFall = 9 ! goto UPPNUT
Else
    nextFall = 7 ! goto HOPPFUT
EndIf
RETURN
END SUBROUTINE fUtt

SUBROUTINE hoppFut
Implicit NONE
Logical :: hoppKlart
! ---------
! HOPPFUT:
! ---------
! if nUt >0 call HOPP to pick a new fut[] array,
!  if all used, then increment nUt at UPPNUT;
!  else goto INFUT, BEFALL, UTFALL, ANFALL, and back to Slingor
If(DBG >= 4) Write(IOUT,'("Fasta() at HoppFut")')

hoppKlart = hopp(NUT+1,Fut,NFSpar)
! print array "fut"
If(DBG >= 4) Call printArraysFasta(.false.,.false.,.false.,.false.,.false.,.false.,.false., .true.)

If(hoppKlart) Then
    nextFall = 9 !goto UPPNUT
Else
    nextFall = 8 !goto INFUT
EndIf

RETURN
END SUBROUTINE hoppFut

SUBROUTINE inFut
! Select new fall[] and ifall[] for the reduced group of solids.
Implicit NONE
Integer :: j, li
! ---------
! INFUT:
! ---------
! find new fall[] and ifall[] for the reduced group of solids
If(DBG >= 4) Then
    Write(IOUT,'("Fasta() at InFut; nfall=",i0,", nfSpar=",i0,", nUt=",i0)') nfall,nfSpar,nUt
    ! print array "fut"
    Call printArraysFasta(.false.,.false.,.false.,.false.,.false.,.false.,.false., .true.)
EndIF

If(NFSPAR > 1) Then !line added jan.2011
    DO li=1,NFSPAR
        FALL(IFSPAR(li))=.true.
    End Do
    If(NUT > 0)  Then
        DO li=1,NUT
            FALL(IFSPAR(FUT(li)))=.false.
        EndDo
    EndIf
    J=0
    DO li=1,NFSPAR
        If(.not.FALL(IFSPAR(li))) Cycle
        J=J+1
        IFALL(J)=IFSPAR(li)
    End Do
    NFALL = J !added jan.2011
Else If(NFALL == 1 .and. NFSPAR == 1 .and. NUT == 1) Then !block added jan.2011
    FALL(IFSPAR(1))=.false.
    IFALL(1) = IFALL(2)
End If !nfSpar >1
!end of addition jan.2011

If(DBG >= 4) Then
    Write(IOUT,'("        after InFut: nfall=",i0)') nfall
    ! print array "ifall"
    Call printArraysFasta(.true.,.false.,.false.,.false.,.false.,.false.,.false.,.false.)
EndIF

If(NFALL == 0) Then
    nextFall = 3 ! goto ANFALL and then Slingor
Else
    nextFall = 2 ! goto BEFALL, UTFALL, ANFALL, and back to Slingor
EndIf

RETURN
END SUBROUTINE inFut

SUBROUTINE uppNut
! Increments nUt and sets initial fut[] values.
Implicit NONE
Integer :: li, lf, ia
! ---------
! UPPNUT:
! ---------
! Increment nUt and set initial fut[] values for this nUt.
! Exit in confusion if nUt > nfSpar
NUT = NUT+1
NFALL = NFSPAR -NUT
If(DBG >= 4) Write(IOUT,'("Fasta() at UppNut;  nUt=",i0,", nfSpar=",i0,",  nfall=",i0,"(=nfSpar-nUt)")') nUt,nfSpar,nfall

If(NFALL > 0)  Then
    DO li=1,NUT
        FUT(li)=li
    End Do
    FUT(NUT+1)=0
Else If(NFALL == 0 .and. NFSPAR == 1 .and. NUT == 1) Then !block added jan.2011
    NFALL = 1
    FUT(1) = 1
    FUT(2) = 0
Else ! end addition jan.2011
    ! --------
    ! "TJI" (Tji is Romani (Gypsy) and means "not, nothing")
    ! --------
    NFALL = 0; NVAF = 0 ! ---- line added 2013-Jan.
    DO lf=1,MSOL
        FALL(lf)=.false.
    End Do
    DO ia=1,NA
        BER(ia)=.false.
        FALLA(ia)=.false.
    End Do
    If(DBG >= 1) Then ! ---- error ----
      Write(IOUT,752) NfSpar
      752 Format('Error in HaltaFall.fasta(): ', &
      'The',I0,'  solids found give only zero determinants ("Tji").',/, &
      'List of solids (array ifSpar):')
      ! print array "ifSpar"
      Call printArraysFasta(.false.,.false.,.false.,.false.,.false.,.false., .true., .false.)
    EndIf
    Call errFlagSet(3) ! failed to find a satisfactory combination of solids.
    INDIK=1
    If(DBG >= 4) Write(IOUT,753) trim(errFlagsToString()),nvaf,iterFasta
    753 Format("HaltaFall.fasta() returns at 'Tji', indik =1, ",A &
                    ", nvaf = ",i0,", iterFasta=",i0,/, &
                    "Failed to find a satisfactory combination of solids.")
    nextFall = 0;
    RETURN ! NYP (IN)

End If

!print array "fut"
If(DBG >= 4) Call printArraysFasta(.false.,.false.,.false.,.false.,.false.,.false.,.false.,.true.)

nextFall = 8 ! goto INFUT
RETURN
END SUBROUTINE uppNut

LOGICAL FUNCTION hopp(N,lista, max)
! Moves to the next permutation of variables in the array list.
! lista = array to move
! max
! Returns false if a new permutation has been found;
!         true if the last permutation has already been used
Implicit NONE
Integer, INTENT(IN)     :: N
Integer, INTENT(IN OUT) :: lista(N)
Integer, INTENT(IN)     :: max
Integer :: lj, i
If(DBG >= 4) Write(IOUT,'("Fasta() at Hopp, max=",i0,", lista()=",10(/3x,10i0))') max,(lista(i),i=1,N)
i = 0
DO While(.true.)
    i = i+1
    if(lista(i) == max) then
        !If(DBG >= 4) Write(IOUT,'(2(a,i0))') "  i=",i,", lista(i)=",lista(i)
        !If(DBG >= 4) Write(IOUT,'(a,10(/3x,10i0))') "   new lista()=",(lista(lj),lj=1,N)
        hopp = .true.; RETURN
    else
        !If(DBG >= 4) Write(IOUT,'(3(a,i0))') "  i=",i,", lista(i)=",lista(i),", lista(i+1)=",lista(i+1)
        if(lista(i+1) == lista(i)+1) Cycle
        lista(i) = lista(i)+1
        if(i>1) then
            Do lj=1, (i-1)
                lista(lj) = lj
            EndDo
        endif
        If(DBG >= 4) Write(IOUT,'("   new lista()=",10(/3x,10i0))') (lista(lj),lj=1,N)
        hopp = .false.; RETURN
    endif
EndDo ! while true
END FUNCTION Hopp

!4--- INVERT -----------------------------------------------------------
SUBROUTINE Invert
! Matrix inversion with accompanying solution of linear equations
! Inverts matrix "RUTA" of size NFALLxNFALL
! Working space: PIVOT(NFALL),IPIVOT(NFALL),INDX(NFALL,2)
Implicit NONE
Real(dp) :: SWAP,zMax,T
Integer :: I,J,K,L,IROW=-1,ICOLUMN=-1,L1,JROW,JCOLUMN,IPIVOT_1
If(DBG >= 4) Then
    Write(IOUT,'(3x,"invert(RUTA(:),",I0,")")') NFALL
    !Do I=1,NFALL
    !    Write(IOUT,'"RUTA(",I0,",J)= ",1P,100(E12.3))') I,(RUTA(I,J),J=1,NFALL)
    !EndDo
EndIf

!    Initialization
DO J=1,NFALL
    IPIVOT(J)=0
End Do
!    Search for pivot element
DO I=1,NFALL
    zMax=0.D0
    DO  J=1,NFALL
        If(IPIVOT(J) == 1) Cycle
        DO  K=1,NFALL
            IPIVOT_1 = IPIVOT(K)-1
            If(IPIVOT_1 == 0) Cycle
            If(IPIVOT_1 > 0) Then
                INDIK=1
                If(DBG >= 4) Write(IOUT,'(3x,"invert() returns indik =1 (matrix singular)")') 
                RETURN
            EndIf
            If((ABS(zMax)-ABS(RUTA(J,K))) <= 0) Then
               IROW=J
               ICOLUMN=K
               zMax=RUTA(J,K)
               End If
        End Do
    End Do
    IPIVOT(ICOLUMN)=IPIVOT(ICOLUMN)+1
    ! Interchange rows to put pivot element on diagonal
    If(IROW /= ICOLUMN) Then
        DO L=1,NFALL
            SWAP=RUTA(IROW,L)
            RUTA(IROW,L)=RUTA(ICOLUMN,L)
            RUTA(ICOLUMN,L)=SWAP
        End Do
    EndIf
    indexInv(I,1)=IROW
    indexInv(I,2)=ICOLUMN
    PIVOT(I)=RUTA(ICOLUMN,ICOLUMN)
    ! Divide pivot row by pivot element
    If(abs(zMax) < 1.D-10) Then
        INDIK=1
        If(DBG >= 4) Write(IOUT,'(3x,"invert() returns indik =1 (matrix singular)")')
        RETURN
    EndIf
    RUTA(ICOLUMN,ICOLUMN)=1.0D0
    DO L=1,NFALL
         RUTA(ICOLUMN,L)=RUTA(ICOLUMN,L)/PIVOT(I)
    End Do
    ! Reduce non-pivot rows
    DO L1=1,NFALL
        If(L1 == ICOLUMN)Cycle
        T=RUTA(L1,ICOLUMN)
        RUTA(L1,ICOLUMN)=0.D0
        DO  L=1,NFALL
            RUTA(L1,L) = RUTA(L1,L) - RUTA(ICOLUMN,L)*T
        End Do
    End Do
End Do

! Interchange columns
DO I=1,NFALL
    L=NFALL+1-I
    If(indexInv(L,1) == indexInv(L,2)) Cycle
    JROW=indexInv(L,1)
    JCOLUMN=indexInv(L,2)
    DO K=1,NFALL
        SWAP=RUTA(K,JROW)
        RUTA(K,JROW)=RUTA(K,JCOLUMN)
        RUTA(K,JCOLUMN)=SWAP
    End Do
End Do
INDIK = 0
If(DBG >= 4) Write(IOUT,'(3x,"invert() returns indik =0 (ok)")')
RETURN
END SUBROUTINE Invert

!3---- CBER --------------------------------------------------------------
SUBROUTINE Cber
!  1.-  Calculation of the activity of soluble complexes.
!  2.-  Calculation of the concentrations of all soluble species
!       with old values for activity coefficients.
Implicit NONE
Real(dp) :: lnC, q
Integer :: LIA,LIX,LIAX
If(DBG >= 6) Write(IOUT,'(3x,"cBer in; ivar=",i0," lnA(",i0,")=",1PG23.15)') ivar,ivar,lnA(ivar)
!Calculate activities of soluble complexes
    DO lix=1,NX
      liax=NA+lix
      q = A(lix,ivar)
      !if(noCalc(ivar)) q = abs(q) ! ### is this needed? ###
      lnA(liax) = lnBA(lix) + q * lnA(ivar)
    End Do
!Calculate Concentrations:
!   components and soluble complexes
    DO lia=1,NION
      C(lia)=0.d0
      If(NOLL(lia)) Cycle
      lnC = lnA(lia) - lnG(lia)
      lnC = min(lnC, 80.5904784d0) ! max value 1.E+35
      C(lia)=exp(lnC) ! this might produce an underflow warning
    End Do
If(DBG >= 6) Write(IOUT,'(3x,"cBer returns")')
RETURN
END SUBROUTINE Cber

!5--- LNABER - Part 1 ------------------------------------------------------
SUBROUTINE LnABer
!Calculates LnKmi and LnA(Ibe()) when some solid phase is assumed
!   to be present.
Implicit NONE
Real(dp) :: W
Integer :: LI,LJ,LIA,IF,LIAF,ia
If(DBG >= 6) Write(IOUT,'(3x,"LnABer in, nfall=",i0)') nfall
!Calculate LNKMI
   DO li=1,NFALL
      if=IFALL(li)
      liaf = NX + if
      W = lnKf(if)
      DO lia=1,NA
         If(.not.BER(lia)) Then
            W = W - A(liaf,lia) * LnA(lia)
         EndIf
      End Do
      lnKmi(li) = W
   End Do
!Calculate LnA(IBE())
   DO li=1,NFALL
      W=0.D0
      DO lj=1,NFALL
         W = W + lnKmi(lj) * RUTA(lj,li)
      End Do
      ia=IBE(li)
      LnA(ia)=W
   End Do
If(DBG >= 6) Write(IOUT,'(3x,"LnABer returns")')
RETURN
END SUBROUTINE lnABer

!6--- LNABAS -----------------------------------------------------------------
SUBROUTINE lnABas
!Calculates those parts of the activity of the soluble complexes that are
! independent of LnA(Ivar), (the LnA varied) and stores them in Lnba(ix).
! This subroutine is called once every time the value of Ivar is changed.
Implicit NONE
Integer :: lix,li
Real(dp) :: q
If(DBG >= 6) Write(IOUT,'(3x,"LnABas(",i0,") in, lnA(",i0,") =",1PG21.12)') ivar,ivar,lnA(ivar)
   X = lnA(ivar)
   DO lix=1,NX
      lnBA(lix)=lnBeta(lix)
      DO li=1,NA
        If(li /= ivar) Then
            q = A(lix,li)
            !if(noCalc(li)) q = abs(q) ! ### is this needed? ###
            lnBA(lix) = lnBA(lix) + q * LnA(li)
        EndIf
      End Do
   End Do
If(DBG >= 6) Write(IOUT,'(3x,"LnABas returns,  x=",1PG20.12)') X
RETURN
END SUBROUTINE lnABas

!7--- TOTBER ----------------------------------------------------------------
SUBROUTINE totBer
! Calculates Y and compares with Y0
! Calculates the total concentration (Y) of a component given the free
! concentrations of all components (logA[]). The procedure is different
! if there are solid phases present or not.
! Returns "indik":
! =1 not ok but it is a component only involved in
!    mononuclear reactions and there are no solids present
! =2 if ok (the Y is equal to Y0 within the tolerance)
! =3 not ok and it is either not a mono component or there are solids present
! =4 if too many iterations (iter[ivar] larger than ITER_MAX+1)
Implicit NONE
Integer :: lix,liax
Real(dp) :: W
If(DBG >= 5) Write(IOUT,'(3x,"totBer() in: ivar=",i0,", indik=",i2,", nFall=",i2,", x=",1PG20.12)') ivar,indik,nfall,x

Y=C(IVAR)
If(NFALL > 0) Then ! Some solid is assumed to be present
   DO lix=1,NX
      liax=NA+LIX
      Y = Y + PVA(ivar,lix)*C(liax)
   End Do
   Y = Y - TOTVA(IVAR)
   Y0 = 0.D0
   W = abs(Y-Y0)
   if(TOLY(IVAR) > 0.D0 .and. abs(TOTVA(IVAR)) > 1.D0 .and. &
       W < abs(TOLY(IVAR)*TOTVA(IVAR))) W=0.D0
Else ! nfall == 0; No solid phase assumed to be present
   DO lix=1,NX
      liax = NA+lix
      Y = Y + A(lix,IVAR)*C(liax)
   End Do
   Y0 = TOT(IVAR)
   W = abs(Y-Y0)
EndIf !NFALL

!Compare Y with Y0
IF(TOLY(IVAR) < W) THEN ! not OK
    If(catchRoundingErrors(IVAR) >= 3) Then !the solution can not be made better: it is uncertain
      INDIK = 2 ! ok
      KARL(IVAR) = 1
      STEG(IVAR) = STEG0_CONT
      ! iter(ivar)=0
      If(DBG >= 5) Then
        Call prnt
        Write(IOUT,'(3x,"totBer() returns; indik = 2 (ok), but ""rounding errors""; iter(",i0,")=",i0)') ivar,iter(ivar)
      EndIf
      Return
    End If

    If(MONO(IVAR) .and. NFALL == 0) Then !Mononuclear component
        If(Y0 <= 0.D0  .OR.  Y <= 0.D0) Then ! can not make logarithms of negative numbers
            INDIK = 3
            ITER(IVAR) = ITER(IVAR) +1
            If(ITER(IVAR) >= ITER_MAX)  GO TO 111
            If(DBG >= 5) Then
                Call prnt
                Write(IOUT,'(3x,"totBer() returns; indik = 3 (not ok & (y or y0) <0); iter(",i0,")=",i0)') &
                    ivar,iter(ivar)
            EndIf
            Return
        EndIf

        LnA(IVAR)=LnA(IVAR)+log(Y0)-log(Y)
        X=LnA(IVAR)
        INDIK=1 ! go to TJAT
        ITER(IVAR)=ITER(IVAR)+1
        If(ITER(IVAR) >= ITER_MAX)  GO TO 111
        If(DBG >= 5) Then
            Call prnt
            Write(IOUT,'(3x,"totBer() returns; indik = 1 (not ok & mono & not solids), mono(",i0,")=true, iter(",i0,")=",i0)') &
                ivar,ivar,iter(ivar)
        EndIf
    Else ! not mono or nfall /=0
        INDIK = 3
        ITER(IVAR)=ITER(IVAR)+1
        If(ITER(IVAR) >= ITER_MAX)  GO TO 111
        If(DBG >= 5) Then
            Call prnt
            Write(IOUT,'(3x,"totBer() returns; indik = 3 (not ok & not (mono or solids)); iter(",i0,")=",i0)') ivar,iter(ivar)
        EndIf
    EndIf ! mono and nfall =0 ?
ELSE ! It is OK
    INDIK=2 ! go to PROV
    KARL(IVAR)=1
    STEG(IVAR)=STEG0_CONT
    ITER(IVAR)=ITER(IVAR)+1
    If(DBG >= 5) Then
        Call prnt
        Write(IOUT,'(3x,"totBer() returns; indik = 2 (ok); iter(",i0,")=",i0)') ivar,iter(ivar)
    EndIf
ENDIF
RETURN
!Too many iterations done with component Ivar
111 Continue
    If(DBG >= 5) Write(IOUT,'(3x,"--- Too many iterations for ivar=",i0)') ivar

    IF(IVABRA(IVAR) /= 0 & ! only print error message for the "outer" loop component
        .and. DBG >= 4) call printTooManyIterations (w)

    If(ITER(IVAR) > (ITER_MAX+1)) Then
        INDIK = 4 ! break the iterations
        KARL(IVAR) = 1
        if(.not. CONT) then
            STEG(IVAR) = STEG0
        else
            STEG(IVAR) = STEG0_CONT
        endif
    EndIf

    If(DBG >= 5) Write(IOUT,'(3x,"totBer() returns; indik = ",i0,"; iter(",i0,")=",i0)') INDIK,ivar,iter(ivar)

RETURN
END SUBROUTINE totBer

!8--- LNABER - Part 2 ----------------------------------------------------
SUBROUTINE lnABer2
!Calculates Totmi and Totva when a solid phase is assumed
Implicit NONE
Real(dp) :: W
Integer :: li,ia,lj
If(DBG >= 6) Write(IOUT,'("   LnABer2 in, nfall=",i0,", nvaf=",i0)') nfall,nvaf
!Calculate Totmi
    DO li=1,NFALL
      ia=IBE(li)
      TOTMI(li) = TOT(ia)-C(ia)
    End Do
!Calculate TOTVA
    DO li=1,NVAF
      ia=IVAF(LI)
      W=TOT(ia)
      DO lj=1,NFALL
         W = W - RUT1(li,lj)*TOTMI(lj)
      End Do
      TOTVA(ia)=W
    End Do
If(DBG >= 6) Write(IOUT,'("   LnABer2 returns")')
RETURN
END SUBROUTINE lnABer2

!9--- ACTICO -------------------------------------------------------
LOGICAL FUNCTION actCoeffs()
! Calculates activity coefficients (array lnG[]).
! If the activity coeficients (which depend on the composition of the
! fluid) have varied since they were last calculated, then 
! the equilibrium composition must be recalculated (by iterations elsewhere)
USE FACTOR_Module, ONLY : FACTOR, log10aH2O !calculateIonicStr
Implicit NONE
Real(dp) :: W, absDiff, maxAbsDiff, maxConc, concLimit, p, q
Integer :: i,j,zz, iMaxDiff
!Real(dp), PARAMETER :: fL=0.5D0, cLim = fL*fL ! parameters to calculate "f" (fL can be = 0.25 or 0.5)
!Real(dp) :: f
Integer :: osc, oscMax ! oscillations
Logical :: OK
! Especially adapted to avoid some oscillatoty behaviour when the
! activity coefficient model is the SIT. These oscillations are
! found in some concentrated aqueous solutions, like AlCl3, etc, where
! the activity coefficients and the aqueous composition are affecting
! each other in a strong way. The procedure does give slightly slower
! convergence for solutions of NaCl for example, but it seems more
! reliable to converge.
! However, ITERC_MAX must be sufficiently large (perhaps 200?).
!
If(DBG >= 5) Write(IOUT,'("--- actCoeffs() in, iter ",i0)') (iterAc+1)

!  --- get lnG[]= natural log of activity coefficient ---
Call FACTOR (NION,C,LnG)

if(jWater > 0) lnA(jWater) = log10aH2O * ln10  ! for H2O

! --- Decide what tolerances in the lnG[] to use.
!     Check what tolerance the user wishes
tolLnG = tolLogF0 * ln10 ! tolLogF0 is between 1.e-10 and 0.01

maxConc = -1.d0;
! --- Increase the tolerance if large concentrations and get maxConc
!f = 1.d0;
W = -1.d0;
Do i=1,NION
    lnG(i) = max(-92.1034d0,lnG(i)) ! min value logf = -40
    if(i == jWater) cycle
    zz = max(abs(z(i)),1)
    If(z(i) /= 0 .and. .not.NOLL(i)) Then
        maxConc = max(maxConc,C(i))
        W = max(W, C(i)*real(zz,dp))
    EndIf
    ! u = C(i) * real((zz*zz),dp)
    ! if(u > cLim) f = min(f, (fL/sqrt(u)))
EndDo

if(W > 1.d0) tolLnG = W * tolLnG
tolLnG = min(tolLnG, 0.01d0 * ln10)
concLimit = maxConc * 1.d-4

iMaxDiff = -1
maxAbsDiff = -1.d0

DO i=1,NION
    if(NOLL(i) .or. C(i) < concLimit) Cycle ! skip species with low concentrations when deciding if it converged or not
    if(i == jWater) Cycle ! if act. coeffs. converge for other species, they must do so for H2O as well
    absDiff = abs(lnG(i) - oldLnG(i))
    if(maxAbsDiff < absDiff) then
        iMaxDiff = i
        maxAbsDiff = absDiff
    endif
END DO

!print C[], lnG and diff =lnG-oldLnG
If(DBG >= 5) Then
    ! if(f < 0.99) Write(IOUT,'("   f=",f10.5," (fraction of how much change in lnG should be applied)")') f
    Call printLnG(0)
EndIf

If(DBG >= 5 .and. maxAbsDiff > 0.) &
        Write(IOUT,'("Max abs(diff) for """,A,""", abs(diff[",i0,"])=",1PG13.4," (in log-10 scale)")') &
                trim(ident(iMaxDiff)), iMaxDiff,(maxAbsDiff/ln10)

OK = (maxAbsDiff <= tolLnG)

If(iterAc == 0) Then ! keep track of oscillations
    Do i=1,NION
        largerOldLnG(1,i) = (oldLnG(i) > lnG(i))
        do j = 2,4
            largerOldLnG(j,i) = largerOldLnG(1,i)
        enddo
    EndDo
EndIf

If(.not. OK) Then
    ! ---- For the SIT model: Instead of going ahead and use the new activity
    !      coefficients in a new iteration step, when species have "large"
    !      concentrations, we do not apply the full change: it is scaled down
    !      to avoid oscillations. This is done with "p" (and "q"=1+p)
    oscMax = 0 ! the max nbr of oscillations for any species
    Do i=1,NION
        if(NOLL(i) .or. C(i) < concLimit) Cycle ! skip species with low concentrations
        if(i == jWater) Cycle
        osc = 0
        if(largerOldLnG(1,i) .neqv. (oldLnG(i) > lnG(i))) osc = 1
        do j = 1,3
            if(largerOldLnG(j,i) .neqv. largerOldLnG((j+1),i)) osc = osc + 1
        enddo
        oscMax = max(osc, oscMax)
    EndDo
    p = min(8.d0, max(0.1d0, maxConc));
    if(iterAc > 1) p = min(10.d0, p + oscMax * 2.d0);
    q = 1.d0 + p;
    w = 1.d0/q
    If(DBG >= 5) write(IOUT,'(" 1/q=",f10.5," oscMax=",i0," maxConc=",f8.3)') w,oscMax,maxConc

    Do i=1,NION
        if(NOLL(i) .or. i == jWater) Cycle
        !if(DBG >= 5 .and. C(i) > concLimit) then
        !    osc = 0
        !    if(largerOldLnG(1,i) .neqv. (oldLnG(i) > lnG(i))) osc = 1
        !    do j = 1,3
        !        if(largerOldLnG(j,i) .neqv. largerOldLnG((j+1),i)) osc = osc + 1
        !    enddo
        !    write(IOUT,'("lnG(",i0,") = oldLnG + ",F7.4,"(lnG - oldLnG);' // &
        !            '  (oldLnG > lnG)=",5L2,2x,i0," lnG=",F10.5," best=",F10.5," old=",F10.5)') &
        !                i,(1./q),(oldLnG(i) > lnG(i)),(largerOldLnG(j,i), j=1,4),osc
        !endif
        ! lnG(i) = (lnG(i) + p * oldLnG(i)) / q;
        ! Note: this is equivalent to
        !   lnG = (LnG + p * oldLnG) / q;
        !   lnG = (LnG + (q-1) * oldLnG) / q;
        !   lnG = (LnG + q * oldLnG - oldLnG) / q;
        !   lnG = (q/q) * oldLnG + (LnG  - oldLnG) / q;
        !   lnG = oldLnG + (1/q)(LnG  - oldLnG);
        !   lnG = oldLnG + f * (lnG - oldLnG);
        ! with f = 1/q = 1/(p+1)
        lnG(i) = oldLnG(i) + w * (lnG(i) - oldLnG(i));
    EndDo

    if(DBG >= 5 .and. iterAC >0) then
        Write(IOUT,'("New values:")')
        Call printLnG(1);
    endif
    If(DBG >= 5 .and. maxAbsDiff > 0.) &
        Write(IOUT,'("Max abs(diff) for """,A,""", abs(diff[",i0,"])=",1PG13.4," (in log-10 scale)")') &
                trim(ident(iMaxDiff)), iMaxDiff,(maxAbsDiff/ln10)

EndIf ! if not OK

Do i=1,NION
    if(NOLL(i) .or. (i == jWater)) Cycle
    largerOldLnG(4,i) = largerOldLnG(3,i)
    largerOldLnG(3,i) = largerOldLnG(2,i)
    largerOldLnG(2,i) = largerOldLnG(1,i)
    largerOldLnG(1,i) = (oldLnG(i) > lnG(i))
    oldLnG(i) = lnG(i)
EndDo

iterAc = iterAc +1

if(iterAc > iterAC_MAX .or. (maxAbsDiff/ln10) > 100.) then ! too many iterations
    OK = .true.
    call errFlagSet(6) ! activity factors did not converge
endif

If(OK) Then
    ! some aqueous concentration(s) too large?
    Do i=1,nIon
        if(C(i) > UNREASONABLE_CONC) then
            if(DBG >= 5) Write(IOUT,'("note: C(",i0,")=",1PE13.4," >",0p,F7.3," (too large)")') i,C(i),UNREASONABLE_CONC
            Call errFlagSet(5) ! some aqueous concentration(s) are too large: uncertain activity coefficients
            !Exit !Do
        endif
    EndDo
EndIf ! OK?

If(DBG >= 4) Then
    If(OK) Then
        if(isErrFlagSet(6)) then !activity factors did not converge
                Write(IOUT,'("--- actCoeffs() returns OK after too many iterations, iterAc=",i0,"  iterAC_Max=",i0)',Advance="no") &
                iterAc,iterAC_MAX
        else
            Write(IOUT,'("--- actCoeffs() returns OK after ",i0," iterations")',Advance="no") iterAc
        endif
    Else ! not OK
            Write(IOUT,'("--- actCoeffs() returns NOT OK. Iteration ",i0)',Advance="no") iterAc
    EndIf ! OK?
    Write(IOUT,'(", tol=",1pE11.3," (in log-10 scale), I=",g11.3,", el.bal.=",g11.3)') (tolLnG/ln10),ionicStrCalc,electricBalance
End If

actCoeffs = OK
RETURN
END FUNCTION actCoeffs

!----------------------------------------------------------------------------
SUBROUTINE prnt
  Implicit NONE
  if (.not.ALLOCATED(TOLY)) then
    Write(*,'("Memory not allocated for arrays in ''prnt''")')
    Call ErrStop
  end if
  call printArrays(.false.,.true.) ! print lnA
  write(IOUT,'(5x,"x=",1p,G21.14,", y=",E22.14,", y0=",E22.14,", tolY[",i0,"]=",E11.4)') x,y,y0, ivar,tolY(ivar)
END SUBROUTINE prnt

!----------------------------------------------------------------------------
SUBROUTINE printTooManyIterations (w)
  Implicit NONE
  Real(dp), INTENT(IN) :: w
  Integer :: li
  If(ITER(IVAR) >= ITER_MAX) Then
    Write(IOUT,20) IVAR,(IVA(LI),LI=1,NVA+1)
    Write(IOUT,21) (TOT(LI),LI=1,NA)
    Write(IOUT,22) (C(LI),LI=1,NION+MSOL)
    20 Format(/3x,"Error: too many iterations with ivar = ",I0,/3X, &
                  "Component nbrs. in the order they are iterated:",20(/7X,12I4))
    21 Format(3X,"TOT()=",20(/5X,1P5G15.7))
    22 Format(3X,"  C()=",100(/5X,1P5G15.7))
  EndIf ! iter(ivar) >= ITER_MAX

  ! if(iter(ivar) = ITER_MAX+1
  Write(IOUT,25) IVAR,ITER(IVAR)
  25 Format(/3x,"Component: ",I0,", iteration: ",I0)
  Call printArrays(.false.,.true.) ! print lnA()
  Write(IOUT,'(6x,"lnF()=",10f14.6,40(:/12x,10f14.6))') (LnG(LI),LI=1,NA)
  If(NFALL == 0) Then
    Write(IOUT,27) Y, Y0, TOLY(IVAR)
    if(.not.MONO(IVAR)) write(IOUT,28) X1(IVAR),X2(IVAR),Y1(IVAR),Y2(IVAR)
    27 Format(6X,"Tot(Calc)=",1PG17.9,", Tot(Input)=",G17.9,", tolerance = ",G12.5)
    28 Format(12X,"low LnA=",1PG23.16,", high LnA=",G23.16,/,12X,"low Tot(Calc)=",G17.9,2X,"high Tot(Calc)=",G17.9)
  Else ! NFALL > 0
    Write(IOUT,41) NFALL,W, TOLY(IVAR)
    Write(IOUT,42) X1(IVAR),X2(IVAR),Y1(IVAR),Y2(IVAR)
    41 Format(5X,"Nr.Solids=",I0,",  error in Tot.Conc.=",1PG17.9,", tolerance = ",G12.5)
    42 Format(12X,"low LnA=",1PG23.16,", high LnA=",G23.16,/,12X, &
                    "low Err.Tot(Calc)=",G17.9,2X,"high Err.Tot(Calc)=",G17.9)
  End If ! NFALL == 0 ?
END SUBROUTINE printTooManyIterations

!-----------------------------------------------------------------------------
SUBROUTINE printLnG (printDiffs)
  ! prints debug information on activity coefficients
  Implicit NONE
  Integer, INTENT(IN) :: printDiffs
  Integer :: i
  if (.not.ALLOCATED(C)) then
    Write(*,'("Memory not allocated for arrays in ''printLnG''")')
    Call ErrStop
  end if

if(jWater > 0) Write(IOUT,'("lnA[H2O]= ",1pg13.6," phi=",0p,f12.7," sumM=",1pg13.6)') lnA(jWater),phi,sumM
If(printDiffs == 0) Then
    Write(IOUT,'("      C()=",1p,7g13.4,100(:/,10x,7g13.4))') (C(i),i=1,nIon)
    Write(IOUT,'("old lnG()=",1p,7g13.6,100(:/,10x,7g13.6))') (oldLnG(i),i=1,nIon)
EndIf
Write(IOUT,'("    lnG()=",1p,7g13.6,100(:/,10x,7g13.6))') (LnG(i),i=1,nIon)
If(printDiffs == 1) Return
Write(IOUT,'("  diffs()=",7f13.6,100(:/,10x,7f13.6))') ((lnG(i)-oldLnG(i)),i=1,nIon)
Write(IOUT,'("tolerance requested (log10)=",1PE9.2,"  used =",E11.4,"  I=",G11.4,"  electric balance = ",E11.4)') &
                tolLogF0, (tolLnG/ln10), ionicStrCalc, electricBalance
RETURN
END SUBROUTINE printLnG

!----------------------------------------------------------------------------
SUBROUTINE printArrays (print_iva, print_lnA)
  Implicit NONE
  Logical, INTENT(IN) :: print_iva, print_lnA
  Integer :: ia
  if (.not.ALLOCATED(IVA)) then
    Write(*,'("Memory not allocated for arrays in ''printArrays''")')
    Call ErrStop
  end if

  If(print_iva) Then
    Write(IOUT,'("      nva =",I4)') nva
    Write(IOUT,'("     iva()=",20I4,20(:/9x,20I4))') (iva(ia), ia=1,Nva+1)
    Write(IOUT,'("  ivabra()=",20I4,20(:/9x,20I4))') (ivaBra(ia), ia=1,NA+1)
    Write(IOUT,'("  ivaNov()=",20I4,20(:/9x,20I4))') (ivaNov(ia), ia=1,NA)
  EndIf
  If(print_lnA) then
    Write(IOUT,'("     lnA()=")', Advance="no")
    do ia=1,NA
        if(lnA(ia) > -230258) then
            Write(IOUT,'(1p,10g25.15)',Advance="no") lnA(ia)
        else
            Write(IOUT,'(" -230258.")',Advance="no")
        endif
    enddo
    Write(IOUT,*)
  endif
  RETURN
END SUBROUTINE printArrays

!----------------------------------------------------------------------------
SUBROUTINE printInput
  USE FACTOR_Module, ONLY : factorPrint
  Implicit NONE
  if (.not.ALLOCATED(TOLY)) then
    Write(*,'("Memory not allocated for arrays in ''prnt''")')
    Call ErrStop
  end if
  Write(IOUT,'("--- HaltaFall - input data:")')
  call printChemSystem(IOUT)
  if(dbg >= 5) call factorPrint(IOUT)
  Write(IOUT,'("--- HaltaFall - end of input data.",/)')
END SUBROUTINE printInput

!----------------------------------------------------------------------------
SUBROUTINE printArraysFasta (ifall_,iber_,fall_,fallA_,ber_,ivaf_,ifSpar_,fut_)
  Implicit NONE
  Logical, INTENT(IN)  :: ifall_,iber_,fall_,fallA_,ber_,ivaf_,ifSpar_,fut_
  Integer :: i

if(MSOL <= 0) RETURN

if (.not.ALLOCATED(IFALL)) then
    Write(*,'("Memory not allocated for arrays in ''printArraysFasta''")')
    Call ErrStop
end if

if(NFALL >0) then
    if(ifall_) write(IOUT,'("  ifall()=",10i4,100(:/,8x,10i4))') (ifall(i),i=1,NFALL)
    if(iber_)  write(IOUT,'("   iber()=",10i4,100(:/,8x,10i4))') (iber(i),i=1,NFALL)
endif
if(fall_ .and. MSOL>0)  write(IOUT,'("    fall()=",10L2,100(:/,8x,10L2))') (fall(i),i=1,MSOL)
if(fallA_) write(IOUT,'("  fallA()=",10L2,100(:/,8x,10L2))') (fallA(i),i=1,NA)
if(ber_)   write(IOUT,'("    ber()=",10L2,100(:/,8x,10L2))') (ber(i),i=1,NA)
if(ivaf_ .and. NVAF>0)  write(IOUT,'("   ivaf()=",10i4,100(:/,8x,10i4))') (ivaf(i),i=1,NVAF)
if(ifSpar_ .and. NFSPAR>0)write(IOUT,'(" ifSpar()=",10i4,100(:/,8x,10i4))') (ifSpar(i),i=1,NFSPAR)
if(fut_ .and. NUT>0)   write(IOUT,'("    fut()=",10i4,100(:/,8x,10i4))') (fut(i),i=1,NUT)
RETURN
END SUBROUTINE printArraysFasta

!----------------------------------------------------------------------------
SUBROUTINE ErrStop ! (PRIVATE)
    Implicit NONE
    Integer :: ierr, i
    WRITE (*, '("Press Enter to continue ...")', Advance='No')
    READ (*,'(I5)',iostat=ierr) i
    STOP 1
END SUBROUTINE ErrStop

END MODULE HALTAF_Module
