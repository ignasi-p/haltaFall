MODULE FACTOR_Module
! A module to calculate single ion activity coefficients in aqueous solutions.
! The procedure "factor(C, lnf)" is to be called repeatedly by "HaltaFall":
! it should be as fast as possible.
!
! If "ionicStr" in CHEM (the input ionic strength, I) is negative (e.g.=-1)
! at the start of the calculations, then the ionic strength will be
! calculated and "ionicStrCalc" will be updated all the time during the
! "HaltaFall" iterations.
!
! The calculation model is defined in "activityCoeffsModel" in module CHEM.
! The user must define Temperature and Pressure before
! calling HaltaFall (which calls FACTOR) if either:
!  - the ionic strenght must be calculated
!  - activityCoeffsModel is not zero
!
! Three (3) models are implemented: Davies eqn., SIT, and simplified HKF.
!
! Note that there is no check that reactions are charge balanced.
!
! This version implements the Davies eqn., the SIT model and the simplified HKF model
!
USE CHEM, ONLY : dp, nl, DBG, activityCoeffsModel, ionicStr, ionicStrCalc, phi, sumM, electricBalance, jWater
USE IO, ONLY : Max_Path
Implicit NONE

! The log10 of the activity of water, calculated in this module.
Real(dp) :: log10aH2O = 0.d0
! Debye-Hückel parameter
Real(dp) :: Agamma
! Debye-Hückel parameter including the size parameter (=Bgamma x å)
Real(dp) :: rB
! the extended Debye-Hückel parameter in the HKF model
Real(dp) :: bgi
! Used internally: true if ionicStr is negative at the start
! of the calculation; false if ionicStr is zero or positive at
! the start of the calculation
Logical :: calculateIonicStr = .true.
! Temperature (Celsius) and pressure (bar), must be given by the user.
! Note that in this module pressure shoult be either 1 bar
! or at the steam saturated line (pSat)
Real(dp) :: Temperature=-1000.D0, Pressure=-1.D0

! the "first time" some checks are made, and it is decided which species
! require activity coefficients, for example, gases are excluded
Logical,  PRIVATE :: firstTime = .true.

Character (LEN=Max_Path), PRIVATE :: ACPATH = ""

! the temperature given by the user the last time this procedure was executed
Real(dp), PRIVATE :: lastTemperature = -274.d0
! the square root of the ionic strength
Real(dp), PRIVATE :: rootI
! The dielectric constant of water
Real(dp), PRIVATE :: eps_H2O
! The density of water in units of g/cm3
Real(dp), PRIVATE :: rhoH2O
! Debye-Hückel parameters
Real(dp), PRIVATE :: Bgamma
! the g-function in the HKF model, in units of Å
Real(dp), PRIVATE :: g_function
! the coefficient in Davies eqn.
Real(dp), PRIVATE, PARAMETER :: DAVIES = 0.3d0
Real(dp), PRIVATE, PARAMETER :: ln10=2.3025850929940d0
! The cut-off abs value for a log10 of the activity coefficients
Real(dp), PRIVATE, PARAMETER :: MAX_LOG_G = 50.d0 ! ### = 115.13 in ln-units
! The cut-off abs value for log10aH2O
Real(dp), PRIVATE, PARAMETER :: MAX_LOGAH2O = 10.d0; ! = 23.03 in ln-units !##
! The cut-off concentration when calculating the ionic strength etc
! to avoid "unreasonable" results
Real(dp), PRIVATE, PARAMETER :: MAX_CONC = 10000.d0 ! ## 200?

! = 55.508 = 1000/(15.9994+2*1.00794)
Real(dp), PRIVATE, PARAMETER :: molH2Oin1kg = 1000.d0/(15.9994d0+2.d0*1.008d0);

!Integer,  PRIVATE :: activityCoeffsModel0 !###

Logical :: fDbg = .false.

PRIVATE :: ErrStop

SAVE

CONTAINS

! ----------------------------------------------------------------------------
SUBROUTINE factorSetPath (path)
Implicit NONE
CHARACTER (LEN=*), INTENT(IN)  :: path
Write(*,'("Setting path for SIT-file to: """,A,"""")') trim(path)
ACPath = path
END SUBROUTINE factorSetPath

! ----------------------------------------------------------------------------
SUBROUTINE FactorStart (nIons,Conc,lnF)
USE CHEM, ONLY : DBG, IOUT !, NA, KH, TOT !###
USE SIT, ONLY : SIT_MEM_ALLOC !, readSITdata
Implicit NONE
Integer, INTENT(IN) :: nIons
Real(dp), INTENT (IN) :: Conc(NIons)
Real(dp), INTENT (OUT) :: lnF(NIons)
!Real(dp) :: maxTotConc !###
Integer :: i

ionicStrCalc = 0.d0;
sumM = 0.D0;
electricBalance = 0.D0;
phi = 1.D0;
firstTime = .false.

If(DBG >= 3) Write(IOUT,'("factorStart(",I0,"), activityCoeffsModel =",I0)') nIons,activityCoeffsModel

If(activityCoeffsModel >= 0 .and. activityCoeffsModel <= 2) Then
    if(ionicStr > 0.D0) then
        ionicStrCalc = min(ionicStr,200.D0)
    else if(ionicStr < 0.D0) then
        ionicStrCalc = -1.D0
    endif
EndIf
calculateIonicStr = (ionicStr < 0.D0)

If(abs(ionicStr) < 1.D-35 .or. activityCoeffsModel < 0 .or. activityCoeffsModel > 2) Then
    DO I=1,NIons
        lnF(I)=0.D0
    EndDo
    ionicStrCalc = 0.d0
    Return
EndIf

Call getH2Oparameters

! If the ionic strength is fixed, the activity coefficients for models that
! depend only on the ionic strength and on temperature-pressure
! (that is: Davies and HKF) are calculated here
If(.not.calculateIonicStr) Then
    rootI = 0.d0
    if(ionicStr > 0.d0) rootI=SQRT(ionicStr)
    if(activityCoeffsModel == 0) call calcDavies (nIons,Conc,lnF)
    if(activityCoeffsModel == 2) call calcHKF (nIons,Conc,lnF)
EndIf

!activityCoeffsModel0 = activityCoeffsModel !###
! if(activityCoeffsModel == 1 .and. calculateIonicStr) then !SIT
    ! !for SIT model, first do a calculation using the simplyfied HKF model
    ! maxTotConc = 0.D0
    ! Do i = 1,NA
        ! if(KH(i) == 1) maxTotConc = max(maxTotConc, abs(TOT(i)))
    ! EndDo
    ! if(maxTotConc > 0.3D0) then
        ! activityCoeffsModel = 2
        ! if(dbg >=3) write(IOUT,'(/,"(Note: performing an initial calculation with the ""simpl.HKF"" model for act. coeffs.)",/)')
    ! endif
! endif

if(activityCoeffsModel == 1) then !SIT
        ! Create SIT-arrays (if needed)
        call SIT_MEM_ALLOC(NIons)
endif

If(DBG >= 3) Write(IOUT,'("factorStart() ut, calculateIonicStr = ",L1,/, &
                          "                  t=",F8.2,"C, rho=",F10.4,", eps=",F9.2)') &
                          calculateIonicStr, temperature,rhoH2O, eps_H2O
RETURN
END SUBROUTINE FactorStart

!----------------------------------------------------------------------------
SUBROUTINE getH2Oparameters
! Calculates temperature- and pressure-dependent parameters for water:
! rho (the density),
! eps_H2O (the static relative permittivity, or dielectric constant)
! and the Debye-Hückel constants and HKF parameters
! In this case they are for 1 bar (at t < 100 C) and at pSat above 100 C
Implicit NONE
! Has the user changed the temperature for the calculation?
! Note: at the start lastTemperature = -274
logical :: tChanged
Real(dp) :: prSat

tChanged = abs(lastTemperature - temperature) > 0.1d0
if(.not.tChanged) then
    lastTemperature = temperature
    Return
endif

! ---- get the density of water (H2O)
rhoH2O = -1.d0; prSat = -1.d0

if(temperature >= 0.d0 .and. temperature < 373.946D0) then
    if(abs(pressure-1.d0) < 1.d-4 .and. temperature <= 100.d0) then ! if p = 1bar and t < 100C
        rhoH2O = rho1bar(temperature);
    else ! either p not = 1 bar or t > 100 or both
        prSat = 1.d0
        if(temperature > 100.d0) prSat = pSat(temperature);
        if(abs(pressure - prSat) < abs(prSat * 0.01d0)) then
            rhoH2O = rhoSat(temperature)
        else
            write(*,*) '"getH2Oparameters": tC =',real(temperature),real(pressure),nl, &
                '  Pressure must be =',prSat
            Call ErrStop
        endif
    endif
else
    write(*,*) '"getH2Oparameters": tC =',real(temperature),nl, &
       ' (must be >=0 and < 373.946 C)'
    Call ErrStop
endif
if(rhoH2O <= 0.d0) then
    write(*,*) '"getH2Oparameters": tC =',real(temperature),real(pressure),nl, &
       '  Can not calculate the density of water (H2O)'
    Call ErrStop
endif

! ---- get the dielectric constant of water (H2O)
eps_H2O = -1.d0;
if(temperature >= 0.d0 .and. temperature < 373.946D0) then
    if(abs(pressure-1.d0) < 1.d-4 .and. temperature <= 100.d0) then ! if p = 1bar and t < 100C
        eps_H2O = eps1bar(temperature);
    else ! either p > 1 bar or t > 100 or both
        prSat = pSat(temperature);
        if(abs(pressure - prSat) < abs(prSat * 0.01d0)) eps_H2O = epsSat(temperature);
    endif
endif
if(eps_H2O <= 0.d0) then
    write(*,*) '"getH2Oparameters": tC =',real(temperature),real(pressure),nl, &
       '  Can not calculate the dielectric constant of water (H2O)'
    Call ErrStop
endif

Agamma = A_gamma(temperature, rhoH2O, eps_H2O);
Bgamma = B_gamma(temperature, rhoH2O, eps_H2O);

if(activityCoeffsModel == 1) then ! SIT
    rB = 4.56620482473d0 * Bgamma;  ! this gives 1.5 at 25°C
else if(activityCoeffsModel == 2) then ! HKF
    g_function = 1.d+10 * gHKF(temperature, pressure) ! convert to Å
    ! use distance of closest approach for NaCl (HKF-4, Table 3)
    ! Shock et al., 1992 J. Chem. Soc., Faraday Trans., 88, 803–826. doi: 10.1039/FT9928800803
    ! Table 9(b) r_e(Cl-)=1.81, Table 9(a) r_e(Na+)=1.91 (=0.97+0.94)
    ! Eqs.(3)-(5) where kz = 0.94 for cations and zero for anions.
    rB = Bgamma * ( (1.81d0+g_function) + (0.97d0 + (0.94d0+g_function)) );
    bgi = b_gamma_NaCl(temperature, pressure);
endif

lastTemperature = temperature;

END SUBROUTINE getH2Oparameters

!----------------------------------------------------------------------------
SUBROUTINE factorPrint (UNIT)
! Prints information (to file UNIT) on how activity coefficients are calculated
  USE SIT, ONLY : printEps !, readSITdata
  Implicit NONE
  Integer, INTENT(IN)  :: UNIT
  Character (len=54), PARAMETER :: sigma ='where sigma(x) = (3/x^3) {(1+x) - 1/(1+x) - 2 ln(1+x)}'
  Character (len=20) :: I, T, P
  Character (len=40) :: S
  Character (len=10)  :: At, Bt, Dt, bgit, rhoT, epsT
  Character (LEN=80) :: dashLINE

  Write(dashLINE,'(39(" -"))')

  If(DBG >= 3) Write(UNIT,'("factorPrint(",I0,")")') UNIT

  write(UNIT,'(a)') dashLINE

  if(activityCoeffsModel < 0 .or. activityCoeffsModel > 2) then
    write(UNIT,'(/,"Activity coefficient calculations will not be performed.")')
    write(UNIT,'(a)') dashLINE
    RETURN
  endif

  ! Calculate rho, eps_H2O, Agamma, etc
  Call getH2Oparameters;

  write(T,'(F10.2)') Temperature; T = adjustL(T)
  write(P,'(F10.2)') Pressure; P = adjustL(P)

  write(At,'(F10.6)') Agamma; At = adjustL(At)
  write(rhoT,'(F10.6)') rhoH2O; rhoT = adjustL(rhoT)
  write(epsT,'(F10.5)') eps_H2O; epsT = adjustL(epsT)
  if(ionicStr < 0) then
    I="varied"
  else
    write(I,'(G20.2)') ionicStr; I = adjustL(I)
  endif

  if(activityCoeffsModel == 0) then !---- Davies

    write(Dt,'(F8.2)') Davies; Dt = adjustL(Dt)

    write(UNIT,'("Calculation of activity coefficients with I = ",a,"  t =",a," C,  p =",a," bar",/, &
        "rho=",a,", eps=",a,",  using Davies eqn.:",/, &
        "  log f(i) = -",a," Zi^2 ( I^0.5 / (1 + I^0.5) - ",a," I)")') &
        trim(I),trim(T),trim(P),trim(rhoT),trim(epsT),trim(At),trim(Dt)
    if(jWater > 0) then
        write(UNIT,'("The activity of water is calculated according to",/, &
        "  log(a(H2O)) = - phi Sum[m] /(ln(10) 55.508)",/, &
        "the osmotic coefficient ""phi"" is calculated from:",/, &
        "  phi = 1 - (2/3) (ln(10)/Sum[m]) ",a," I^(3/2) sigma(I^0.5)" ,/, &
        "             + (ln(10)/Sum[m]) (",a," x ",a,") I^2",/,a,".")') &
            trim(At),trim(At),trim(Dt),sigma
    endif

  else if(activityCoeffsModel == 1) then !---- SIT

    write(Bt,'(F8.4)') rB; Bt = adjustL(Bt)

    write(UNIT,'("Calculation of activity coefficients with I = ",a," and t =",a," C,  p =",a," bar",/, &
        "rho=",a,", eps=",a,",  using the SIT model:",/, &
        "  log f(i) = -(",a," Zi^2 I^0.5 / (1 + ",a," I^0.5)) + Sum[ eps(i,j) m(j) ]",/, &
        "for all ""j"" with Zj*Zi<=0;  where eps(i,j) is a specific ion interaction",/, &
        "parameter (in general independent of I).")') trim(I),trim(T),trim(P),trim(rhoT),trim(epsT),trim(At),trim(Bt)
    if(jWater > 0) then
        write(UNIT,'("The activity of water is calculated with:",/, &
        "  log(a(H2O)) = - phi Sum[m] /(ln(10) 55.508)",/, &
        "the osmotic coefficient ""phi"" is calculated from:",/, &
        "  phi = 1 - (2/3) (ln(10)/Sum[m]) ",a," I^(3/2) sigma(",a," I^0.5)",/, &
        "          + (ln(10)/Sum[m]) Sum_i[ Sum_j[ eps(i,j) m(i) m(j) ]]",/,a,";",/, &
        "and where ""i"" are cations or neutral species and ""j"" are anions.")') &
             trim(At),trim(Bt),sigma
    endif
    call printEps(UNIT)

  else if(activityCoeffsModel == 2) then !---- HKF

    if(ionicStr < 0.d0) then
        S="I)"
    else
        write(S,'(G20.6)') -log10(1.d0+0.0180153d0*ionicStr); S = adjustL(S)
        S = "I) = " // trim(S)
    endif
    write(Bt,'(F8.3)') rB; Bt = adjustL(Bt)
    write(bgit,'(F8.4)') bgi; bgit = adjustL(bgit)

    write(UNIT,'("Calculation of activity coefficients with I = ",a,",  t =",a," C,  p =",a," bar",/, &
        "rho=",a,", eps=",a,",  using a simplified Helgeson, Kirkham & Flowers model:",/, &
        "  log f(i) = -(",a," Zi^2 I^0.5) / (1 + ",a," I^0.5) + Gamma + (",a," I)",/, &
        "where: Gamma = -log(1+0.0180153 ",a)') trim(I),trim(T),trim(P),trim(rhoT),trim(epsT),trim(At),trim(Bt),trim(bgit),trim(S)
    if(jWater > 0) then
        write(UNIT,'("The activity of water is calculated according to",/, &
        "  log(a(H2O)) = - phi Sum[m] /(ln(10) 55.508)",/, &
        "the osmotic coefficient ""phi"" is calculated from:",/, &
        "  phi = (ln(10) Gamma / (0.0180153 I))",/, &
        "              - (2/3) (ln(10)/Sum[m]) ",a," I^(3/2) sigma(",a," I^0.5)",/, &
        "              + (ln(10) (1/2) ",a," I)",/,a,".")') &
            trim(At),trim(Bt),trim(bgit),sigma
    endif
  endif

write(UNIT,'(a)') dashLINE

RETURN
END SUBROUTINE factorPrint

!-----------------------------------------------------------------------------
REAL(dp) FUNCTION eps1bar(tC) RESULT(diel)
! Returns the static permittivity of water (dielectric constant, unitless)
! at 1 bar, calculated using the equations reported in:
! Pátek J, Hrubý J, Klomfar J, Součková M, Harvey A H, 2009.
! Reference correlations for thermophysical properties of liquid
! water at 0.1MPa. J. Phys. Chem. Reference Data 38, 21–29.
! doi:10.1063/1.3043575
!
! Range of conditions: -16°C to 110°C
!
Implicit NONE
Real(dp), INTENT(IN) :: tC   ! the temperature in degrees Celsius
Real(dp) :: TK,tStar
Character (len = 1), PARAMETER :: nl = new_line('a')
Integer :: i
! Table 7 in Pátek et al (2009)
Real(dp), PARAMETER :: e(*) = (/ -43.7527d0, 299.504d0, -399.364d0, 221.327d0/)
Real(dp), PARAMETER :: f(*) = (/ -0.05d0,    -1.47d0,   -2.11d0,   -2.31d0/)

if(tC < -16.15d0 .or. tC > 110.d0) then
    write(*,*) '"eps1bar(t)": tC =',real(tC),nl,' (must be >=-16 and <= 110 C)'
    Call ErrStop
    endif
TK = tC + 273.15D0
!  Eqn. (11) in Pátek et al (2009)
tStar = tK / 300.d0
diel = 0.d0;
Do i=1, 4
    diel = diel + e(i) * tStar**f(i)
End Do
RETURN ! diel
END FUNCTION eps1bar

!-----------------------------------------------------------------------------
REAL(dp) FUNCTION epsSat(tC) RESULT(diel)
! Returns the static permittivity of water (dielectric constant, unitless)
! at the liquid-vapor saturated pressure, calculated using the
! equations reported in:
! Fernández, D.P., Goodwin, A.R.H., Lemmon, E.W., Levelt Sengers, J.M.H.,
! Williams, R.C., 1997. A formulation for the static permittivity of water and
! steam at temperatures from 238 K to 873 K at pressures up to 1200 MPa,
! including derivatives and Debye–Hückel coefficients.
! Journal of Physical and Chemical Reference Data 26, 1125–1166.
! doi: 10.1063/1.555997
!
! If tC = 0 it returns 87.956, the value at the triple point (0.01 C and 0.00611657 bar).
! Range of conditions: 0.01°C (triple point) to 373.9°C (critical point).
Implicit NONE
Real(dp), INTENT(IN) :: tC   ! the temperature in degrees Celsius
Real(dp) :: T, theta, epsLiq
Character (len = 1), PARAMETER :: nl = new_line('a')
Integer :: i
! Table 8 in Fernandez et al (1997)
Real(dp), PARAMETER :: L(*) = (/ 2.725384249466D0,1.090337041668D0,21.45259836736D0,-47.12759581194D0, &
     4.346002813555D0,237.5561886971D0,-417.7353077397D0,249.3834003133D0/)

if(abs(tC) < 1.D-8) then
    diel=87.95631D0 ! the value at the triple point: 0.01 C and pBar = 0.00611657 bar
    RETURN
endif
if(tC < 0.01D0 .or. tC >= 373.946D0) then
    write(*,*) '"epsSat(t)": tC =',real(tC),nl,' (must be either zero or >=0.01 and < 373.946 C)'
    Call ErrStop
    endif
T = tC + 273.15D0
! Eqns. (35) and (36) in Fernandez et al (1997)
theta = (1.d0-(T/647.096))**(1.d0/3.d0)
epsLiq = 1.D0;
 Do i=1, 8
      epsLiq = epsLiq + L(i) * theta**(i)
      end Do
diel = epsLiq * 5.36058D0
RETURN
END FUNCTION epsSat

!-----------------------------------------------------------------------------
REAL(dp) FUNCTION epsJN(tC,pBar) RESULT(epsH2O)
! Returns the static permittivity of water (dielectric constant, unitless)
! calculated using the equations reported in:
! Johnson J W, Norton D (1991) Critical phenomena in hydrothermal systems:
! state, thermodynamic, electrostatic, and transport properties of H2O in the
! critical region. American Journal of Science, 291, 541–648.
! doi: 10.2475/ajs.291.6.541
! 
! see also eqns. (36) and (40)-(44) and Table 1 in:
! Johnson J W, Oelkers E H, Helgeson H C (1992) SUPCRT92: A software package
! for calculating the standard molal thermodynamic properties of minerals,
! gases, aqueous species, and reactions from 1 to 5000 bar and 0 to 1000°C.
! Computers & Geosciences, 18, 899–947. doi:10.1016/0098-3004(92)90029-Q
!
! If tC = 0 and pBar = 1, it returns the value at 0.01 Celsius and 1 bar.
! Range of conditions: 0 to 12,000 bar, -35 to 600°C.
Implicit NONE
Real(dp), INTENT(IN) :: tC   ! the temperature in degrees Celsius
Real(dp), INTENT(IN) :: pBar ! the pressure in bar
Real(dp) :: pSt, t_C, r, r2, T, T_, T_2, k1,k2,k3,k4
Real(dp), PARAMETER :: Tr = 298.15D0
Character (len = 1), PARAMETER :: nl = new_line('a')
Real(dp), PARAMETER :: a(10) = (/ &
      0.1470333593D+2, &
      0.2128462733D+3, &
     -0.1154445173D+3, &
      0.1955210915D+2, &
     -0.8330347980D+2, &
      0.3213240048D+2, &
     -0.6694098645D+1, &
     -0.3786202045D+2, &
      0.6887359646D+2, &
     -0.2729401652D+2 /)
epsH2O=0.d0
if(tC < -35.d0 .or. tC > 600.001d0) then
    write(*,*) '"epsJN": tC =',real(tC),' (must be >-35 and < 600 C)'
    Call ErrStop
    endif
if(abs(pBar) < 1.D-10 .or. pBar < 0.d0 .or. pBar > 10000.001d0) then
    write(*,*) '"epsJN": pBar =',real(pBar),' (must be >0 and < 10 kbar)'
    Call ErrStop
    endif
! if pressure = 1 bar and temperature = 0, set temperature to 0.01 C
! (water freezes at zero degrees and 1 bar)
if(abs(pBar - 1.d0) < 1.D-5 .and. abs(tC) < 0.01D0) then
    t_C = 0.01D0
else
    t_C = tC
endif
r = -1.d0
if(t_C >= 0.01d0 .and. t_C < 373.946d0) then
    pSt= pSat(t_C)
    if(t_C < 99.6059d0) then
        if(abs(pBar-1.d0) > 1.d-5) then
            write(*,*) '"epsJN": temperature =',real(t_C),' pressure =',real(pBar),nl, &
                       '  (must be = 1 bar at temperatures <= 100 C)'
            Call ErrStop
        endif
    else if(pBar < (pSt*0.985d0)) then
        write(*,*) '"epsJN": t_C =',real(t_C),' pBar =',real(pBar),nl, &
                   'pBar must be >=',real(pSt)
        Call ErrStop
    else if(pBar > (pSt*1.015d0)) then
        write(*,*) '"epsJN": t_C =',real(t_C),' pBar =',real(pBar),' Steam-sat. p=',real(pSt),nl, &
                   ' Calculation outside steam saturation is not implemented yet.'
        Call ErrStop
    endif
    r = rhoSat(t_C)
else
    write(*,*) '"epsJN": t_C =',real(t_C),' pBar =',real(pBar),nl, &
               ' Calculation above the critical point is not implemented yet.'
    Call ErrStop
    !r = rho(t_C, pBar)
endif
if(r < 0.05D0) then
    write (*,*) '"epsJN(t,p)": calculated density =',real(r),' at tC=',real(t_C),' pBar =',real(pBar),nl, &
                ' (density must be >0.05 g/cm3)'
    Call ErrStop
endif
T = t_C + 273.15D0
T_= T/Tr
T_2 = T_*T_
r2 = r*r
k1 = a(1)/T_
k2 = a(2)/T_ + a(3) + a(4)*T_
k3 = a(5)/T_ + a(6)*T_ + a(7)*T_2
k4 = a(8)/T_2 + a(9)/T_ + a(10)
epsH2O = 1.D0 + k1*r + k2*r2 + k3*r2*r + k4*r2*r2
END FUNCTION epsJN

!-----------------------------------------------------------------------------
REAL(dp) FUNCTION p1Sat(tC) RESULT(p)
! Returns 1 bar at tC below 99.6059 C, and the saturation pressure (bar)
! of ordinary water substance at temperatures above that.
!
! Range of conditions: 0.01 to 373.946 °C.
Implicit NONE
Real(dp), INTENT(IN) :: tC   ! input temperature in degrees Celsius (>= 0 and < 373.946)
p = max(1.d0,pSat(tC))
Return
END FUNCTION p1Sat

!-----------------------------------------------------------------------------
REAL(dp) FUNCTION pSat(tC) RESULT(p)
! Returns the saturation pressure (bar) of ordinary water substance, that is,
! the pressure at the vapor–liquid phase boundary, for a given temperature (in
! degrees Celsius). Uses eqn. (2.5) in Wagner, W., Pruß, A., 2002. "The IAPWS
! formulation 1995 for the thermodynamic properties of ordinary water substance
! for general and scientific use. Journal of Physical and Chemical Reference
! Data 31, 387–535. doi: 10.1063/1.1461829.
!
! Note: it returns pressures below 1 bar at tC below 99.6059 C.
! If tC = 0, the pressure (0.00611657 bar) at the triple point (0.01 C) is returned.
! Range of conditions: 0.01 to 373.946 °C.
Implicit NONE
Real(dp), INTENT(IN) :: tC   ! input temperature in degrees Celsius (>= 0 and < 373.946)
Real(dp) :: a1,a2,a3,a4,a5,a6, tK, theta, ln_pSat

if(abs(tC) < 0.0001D0) then
    p = 0.00611657D0
    Return
end if
if(abs(tC-25.D0)<0.01) then
    p = 0.031698246D0
    Return
endif
! Critical temperature = 647.096 (+/- 0.01) K  (= 373.946 C)
if(tC < 0.01D0 .or. tC >= 373.946D0) then
    write(*,*) '"pSat(t)": tC =',real(tC),' must be either zero or >=0.01 and < 373.946 C'
    Call ErrStop
endif
a1 = -7.85951783D0; a2 = 1.84408259D0;  a3 = -11.7866497D0
a4 = 22.6807411D0;  a5 = -15.9618719D0; a6 = 1.80122502D0
tK = tC + 273.15D0
theta = 1.D0 - (tK / 647.096D0)
! CRITICAL_p = 22.064 MPa
ln_pSat = log(22.064D0) + (647.096D0 / tK)*(a1*theta + a2*(theta**1.5D0) + a3*(theta**3) &
             + a4*(theta**3.5D0) + a5*(theta**4) + a6*(theta**7.5))
p = exp(ln_pSat)*10.D0; ! convert MPa to bar
RETURN
END FUNCTION pSat

!-----------------------------------------------------------------------------
REAL(dp) FUNCTION rho1bar(tC) RESULT(rho) ! the density in units of g/cm3
! Returns the density (g/cm3) of ordinary water substance
! at 1 bar, for a given temperature in degrees Celsius.
! Uses eqn. (3) in:
! Pátek J, Hrubý J, Klomfar J, Součková M, Harvey A H, 2009.
! Reference correlations for thermophysical properties of liquid water
! at 0.1MPa. J. Phys. Chem. Reference Data, 38, 21–29. doi:10.1063/1.3043575
! Range of conditions: -20 to 110 °C.
Implicit NONE
Real(dp), INTENT(IN) :: tC  ! tC input temperature in degrees Celsius (>=0 and < 373.946)
Real(dp), PARAMETER :: TR = 10.d0, Ta = 593.d0, Tb = 232.d0
! the gas constant (kJ/(kg K)) for a molar mass for water = 18.015268 g/mol
Real(dp), PARAMETER :: R = 0.46151805d0

Real(dp) :: v, tK, alpha, beta
Integer :: i,j
! from Table 1:
Real(dp), PARAMETER :: n(6) = (/ -1.d0, 4.d0, 5.d0, 7.d0, 8.d0, 9.d0 /)
Real(dp), PARAMETER :: m(6) = (/  1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0 /)
Real(dp), PARAMETER :: a(6) = (/ &
        1.93763157d-2, &
        6.74458446d+3, &
       -2.22521604d+5, &
        1.00231247d+8, &
       -1.63552118d+9, &
        8.32299658d+9 /)
Real(dp), PARAMETER :: b(6) = (/ &
        5.78545292d-3, &
       -1.53195665d-2, &
        3.11337859d-2, &
       -4.23546241d-2, &
        3.38713507d-2, &
       -1.19946761d-2 /)
if(tC < -20.d0 .or. tC > 110.) then
    write(*,*) '"rho1bar(t)": tC =',real(tC),' (must be >-20 and < 110 C)'
    Call ErrStop
endif
tK = tC + 273.15d0
alpha = TR / (Ta - tK)
beta = TR /(tK - Tb)
v = a(1)
do i = 6,10
    j = i-4
    v = v + a(j) * alpha**n(j)
end do
do i = 5,10
    j = i-4
    v = v + b(j) * beta**m(j)
end do

rho = 0.1d0 / (v * R * TR)
RETURN
END FUNCTION rho1bar

!-----------------------------------------------------------------------------
REAL(dp) FUNCTION rhoSat(tC) RESULT(rho) ! the density in units of g/cm3
! Returns the density (g/cm3) of ordinary water substance at the
! vapor–liquid phase boundary, for a given temperature in
! degrees Celsius. Uses eqn. (2.6) in Wagner, W., Pruß, A., 2002. "The IAPWS
! formulation 1995 for the thermodynamic properties of ordinary water substance
! for general and scientific use. Journal of Physical and Chemical Reference
! Data 31, 387–535. DOI: 10.1063/1.1461829.
!
! If tC = 0, the density (0.9997891 g/cm3) at the triple point
! (0.01 C and 0.00611657 bar) is returned.
! Range of conditions: 0.01 to 373.946 °C.
Implicit NONE
Real(dp), INTENT(IN) :: tC  ! tC input temperature in degrees Celsius (>=0 and < 373.946)
Real(dp) :: b1,b2,b3,b4,b5,b6, tK, theta

if(abs(tC) < 0.0001D0) then
    rho = 0.9997891D0 ! the calculated density at 0.01 C and 0.00611657 bar
    Return
end if
if(abs(tC-25.D0)<0.01) then
    rho = 0.9969994D0
    Return
endif
! Critical temperature = 647.096 (+/- 0.01) K  (= 373.946 C)
if(tC < 0.01D0 .or. tC >= 373.946D0) then
    write(*,*) '"rhoSat(t)": tC =',real(tC),' (must be either zero or >=0.01 and < 373.946 C)'
    Call ErrStop
endif
b1 = 1.99274064D0;  b2 = 1.09965342D0;  b3 = -0.510839303D0
b4 = -1.75493479D0; b5 = -45.5170352D0; b6 = -6.74694450D+5
tK = tC + 273.15D0
theta = 1.D0 - (tK / 647.096D0)
! CRITICAL_rho = 322; ! kg/m3
rho = 322.D0 * &
            (1.D0 + b1*(theta**(1.D0/3.D0)) + b2*(theta**(2.D0/3.D0)) + b3*(theta**(5.D0/3.D0)) &
             + b4*(theta**(16.D0/3.D0)) + b5*(theta**(43.D0/3.D0)) + b6*(theta**(110.D0/3.D0)))
rho = rho/1000.D0
RETURN
END FUNCTION rhoSat

!-----------------------------------------------------------------------------
REAL(dp) FUNCTION A_gamma(tC, rho, eps) RESULT(A_g)
! Calculates the Debye-Hückel slope (in units of (kg/mol)^0.5) as defined
! in eqn(30) of:
! Staples, B.R., Nuttall, R.L., 1977.The activity and osmotic coefficients of
! aqueous calcium chloride at 298.15 K.J. Phys. Chem. Ref. Data vol.6, p.385–407.
! See also eqn(2-2) in:
! Hamer, W.J. and Wu, Y.-C., 1972. Osmotic coefficients and mean activity
! coefficients of uni-univalent electrolytes in water at 25 °C.
! J. Phys. Chem. Ref. Data, vol.1, p.1047–1099.
Implicit NONE
Real(dp), INTENT(IN) :: tC  ! the temperature in degrees Celsius
Real(dp), INTENT(IN) :: rho ! the density of water in g/cm3
Real(dp), INTENT(IN) :: eps ! the dielectric constant of water at the given temperature and density (unitless)
  ! The Debye-Hückel slope is:
  ! ((1/ln(10))*(2*pi*Na*rho))^0.5
  !    * (e^2 / (4*pi*eps0*eps*k*T))^(3/2)
  ! where:
  ! Na = 6.02214076e23 1/mol (Avogadro's number)
  ! eps0 = 8.8541878128e-12 F/m (permittivity of vacuum)
  ! e = 1.60217662e-19 C (elementary charge)
  ! k = 1.380649e-23 J/K (Boltzmann's constant)
  ! rho = the density of water in g/cm^3
  ! eps = the dielectric constant of water
  ! T = the temperature in Kelvin
A_g = (1.82481168d+6 * sqrt(rho)) / ((eps * (tC + 273.15d0))**(1.5d0))
RETURN
END FUNCTION A_gamma

!-----------------------------------------------------------------------------
REAL(dp) FUNCTION B_gamma(tC, rho, eps) RESULT(B_g)
! Calculates the Debye-Hückel slope (in units of (kg/mol)^0.5 * Å^-1)
! as defined by eqn(35) in:
! Staples, B.R., Nuttall, R.L., 1977. The activity and osmotic coefficients of
! aqueous calcium chloride at 298.15 K. J. Phys. Chem. Ref. Data, vol.6, p.385–407.
! See also eqn(2-5) in:
! Hamer, W.J. and Wu, Y.-C., 1972. Osmotic coefficients and mean activity
! coefficients of uni-univalent electrolytes in water at 25 °C.
! J. Phys. Chem. Ref. Data, vol.1, p.1047–1099.
Implicit NONE
Real(dp), INTENT(IN) :: tC  ! the temperature in degrees Celsius
Real(dp), INTENT(IN) :: rho ! the density of water in g/cm3
Real(dp), INTENT(IN) :: eps ! the dielectric constant of water at the given temperature and density (unitless)
  ! The Debye-Hückel [B*å] constant is:
  ! ((8*pi * Na * 1000)^(0.5) * e
  !    / (4*pi*eps0*eps*k*T))^(0.5) * (rho^0.5) * å
  ! where:
  ! Na = 6.02214076e23 1/mol (Avogadro's number)
  ! e = 1.602176634e-19 C (elementary charge)
  ! eps0 = 8.8541878128e-12 F/m (permittivity of vacuum)
  ! k = 1.380649e-23 J/K (Boltzmann's constant)
  ! rho = the density of water in g/cm^3
  ! eps = the dielectric constant of water
  ! T = the temperature in Kelvin
  ! å = distance parameter in units of metres
  B_g = 5.0290371d+1 * sqrt(rho/(eps * (tC + 273.15d0)));
RETURN
END FUNCTION B_gamma

!-----------------------------------------------------------------------------
REAL(dp) FUNCTION gHKF(tC, pBar) RESULT(g_HKF)
! HKF (Helgeson-Krikham-Flowers) model:
! g (in units of metre (m)) designates a P/T-dependent solvent function that
! provides for dielectric saturation and the compressibility of the solvent
! at high temperatures and pressures.  With this function the temperature
! and pressure dependence of the effective electrostatic radii of aqueous ions
! can be calculated.  See eqns.(49) to (52) and Tables 2 and 3 in:
! 
! Johnson J W, Oelkers E H, Helgeson H C, (1992) SUPCRT92: A software package
! for calculating the standard molal thermodynamic properties of minerals,
! gases, aqueous species, and reactions from 1 to 5000 bar and 0 to 1000°C.
! Computers & Geosciences, 18, 899–947. doi: 10.1016/0098-3004(92)90029-Q
! 
! as well as eqns. (25), (26), (32), (33) and Table 5 in:
! Shock E L, Oelkers E H, Johnson J W, Sverjensky D A, Helgeson H C, (1992)
! Calculation of the thermodynamic properties of aqueous species at high
! pressures and temperatures. Effective electrostatic radii, dissociation
! constants and standard partial molal properties to 1000°C and 5 kbar.
! J. Chem. Soc., Faraday Transactions, 88, 803–826. doi:10.1039/FT9928800803
!
! Range of conditions: 0 to 5000 bar, 0 to 1000°C, density 0.35 to 1 g/cm3
! and pressures above the vaporization boundary (if temperature is below the
! critical point). 
Implicit NONE
Real(dp), INTENT(IN) :: tC    ! the temperature in degrees Celsius
Real(dp), INTENT(IN) :: pBar  ! the pressure in bar
Character (len = 1), PARAMETER :: nl = new_line('a')
Real(dp) :: pSt, t_C, r, a1,a2,a3, b1,b2,b3, c1,c2,c3, t, p, f, w, ag,bg
if(tC > 1000.d0) then
    write(*,*) '"gHKF": tC =',real(tC),' must be <= 1000 C'
    Call ErrStop
    endif
if(abs(pBar) < 1.d-10 .or. pBar < 0.d0 .or. pBar > 5000.d0) then
    write(*,*) '"gHKF": pBar =',real(pBar),' must be >0 and <= 5 kbar'
    Call ErrStop
    endif
!--------
  if(tC < 100.01) then
    g_HKF = 0.d0
    RETURN
  endif
!--------
! if pressure = 1 bar and temperature = 0, set temperature to 0.01 C
! (water freezes at zero degrees and 1 bar)
if(abs(pBar - 1.d0) < 1.D-5 .and. abs(tC) < 0.01D0) then
    t_C = 0.01D0
else
    t_C = tC
endif
r = -1.d0
pSt = 1.d0
if(t_C >= 0.01d0 .and. t_C < 373.946d0) then
    pSt= pSat(t_C)
    if(t_C < 99.6059d0) then
        if(abs(pBar-1.d0) > 1.d-5) then
            write(*,*) '"gHKF": temperature =',real(t_C),' pressure =',real(pBar),nl, &
                       '  (must be = 1 bar at temperatures <= 100 C)'
            Call ErrStop
        endif
    else if(pBar < (pSt*0.985d0)) then
        write(*,*) '"gHKF": t_C =',real(t_C),' pBar =',real(pBar),nl, &
                   'pBar must be >=',real(pSt)
        Call ErrStop
    else if(pBar > (pSt*1.015d0)) then
        write(*,*) '"gHKF": t_C =',real(t_C),' pBar =',real(pBar),' Steam-sat. p=',real(pSt),nl, &
                   ' Calculation outside steam saturation is not implemented yet.'
        Call ErrStop
    endif
    r = rhoSat(t_C)
else
    write(*,*) '"gHKF": t_C =',real(t_C),' pBar =',real(pBar),nl, &
               ' Calculation above the critical point is not implemented yet.'
    Call ErrStop
    !r = rho(t_C, pBar)
endif
if(r < 0.15d0) then
        write(*,*) '"gHKF": calculated density =',r,' at tC =',tC,' pBar =',pBar,nl, &
                   ' (density must be >0.15 g/cm3)'
        Call ErrStop
endif
!--------
if(r > 1.d0) then
    g_HKF = 0.d0
    RETURN
endif
!--------
! Expressions in Johnson et al (1992)
!   g = a (1 - ?*)^b - f(T,P)     (eq.49)
!   a = a1 + a2 T + a3 T^2         (eq.50)
!   b = b1 + b2 T + b3 T^2         (eq.51)
!   ?* = ?(T,P) / (1 g/cm3)
!   f(T,P) = [ ((T-155)/300)^(4.8)
!           + c1 ((T-155)/300)^16 ]
!         ×[c2 (1000-P)^3 + c3 (1000-P)^4]     (eq.52)
!Table 2 in Johnson et al (1992)
  a1 = -2.037662d0;   b1 =  6.107361d0
  a2 =  5.747000D-3;  b2 = -1.074377D-2
  a3 = -6.557892D-6;  b3 =  1.268348D-5
!Table 3 in Johnson et al (1992)
  c1 = 36.6666d0;  c2 = -1.504956D-10;  c3 = 5.01799D-14
T = tC; P = pBar
if(T < 155.d0 .or. T > 355.d0 .or. P < pSt .or. P > 1000.d0) then
    f = 0.d0
else
    w = (T-155.d0)/300.d0
    f = ((w**4.8d0) + c1*(w**16)) &
          * ((c2*((1000.d0-P)**3))+(c3*(((1000.d0-P)**4))))
endif
ag = (a1+a2*T+a3*T*T)
bg = (b1+b2*T+b3*T*T)
write(*,*) "t=",T," p=",P," ag=",ag," bg=",bg," f=",f," r=",r
g_HKF = 1.d-10*((ag*(((1.d0-r)**bg)))-f); ! Convert Angstrom to metres
RETURN
END FUNCTION gHKF

!-----------------------------------------------------------------------------
REAL(dp) FUNCTION b_gamma_NaCl(tC, pBar) RESULT(b_gamma)
! Returns b_gamma_NaCl (units of (kg/mol)) in the HKF model of
! activity coefficients.
! Range of conditions: < 5,000 bar, 0 to 1000°C.
! It stops outside this range.
! If pBar < 1 a value of pBar = 1 is used.
! If tC = 0 it returns the value at 0.01 Celsius.
! NOTE: values returned at tC <0 are extrapolations outside the valid range.
! It requires a values of "eps_H2O" and of "gHKF"
!
! b_gamma_NaCl is “defined” in Eqn. (173) of Helgeson, Kirkham and Flowers (1981)
! values at 25°C and 1 bar are listed in Tables 5 and 6, while Table 26 gives
! values at temps up to 325°C and Psat (pressures corresponding to the
! liquid-vapor equilibrium) and Table 27 lists values at temperatures up to
! 500°C and pressures up to 5000 bar.
!
! See also eqn.(22) in  Oelkers and Helgeson (1990) where "b?,NaCl" is given in
! Table A-2, and the eqn. (B-12) and parameters in Table A-1.  Note that the
! values at 3 kbar and temps. of 800-1000°C in the Table A-2 are wrong:
! they must be multiplied by 10.  Note also that eqn.(B-13) in that paper
! is wrong.  The same expression for b? may be found as eqns.(2) and (31)-(32)
! in Pokrovskii and Helgeson (1997).
!
! To calculate values in Table A-2 of Oelkers and Helgeson (1990) using eqn (B-12)
! or the equivalent eqns. (31)-(32) in Pokrovskii and Helgeson (1997) one needs
! values of the dielectric constant (e) and of ?, which requires values of the
! g-function (see e.g. eqn.(15) in Pokrovskii and Helgeson 1997).
!
! Values of the relative permittivity of water (dielectric constant e) are
! calculated using the equations of Johnson and Norton (1991).  See also
! Johnson et al (1992) and Shock et al (1992) which lists values of e in Table C2.
!
! The g-function is described in Johnson et al (1992), eqns. (49)-(51) and
! parameters in Tables 2 and 3.  See also Shock et al. (1992), where values of
! g are given in Table 5.
!
! References: 
! Helgeson, Kirkham and Flowers, Amer. J. Sci. 281 (1981) 1249-1516, doi: 10.2475/ajs.281.10.1249
! Oelkers and Helgeson, Geochim. Cosmochim. Acta 54 (1990) 727-738 doi: 10.1016/0016-7037(90)90368-U
! Johnson, Norton. Amer. J. Sci., 291 (1991) 541–648. doi: 10.2475/ajs.291.6.541
! Johnson, Oelkers, Helgeson. Computers & Geosciences, 18 (1992) 899–947. doi: 10.1016/0098-3004(92)90029-Q
! Shock, Oelkers, Johnson, Sverjensky, Helgeson. J. Chem. Soc., Faraday Trans., 88 (1992) 803–826. doi: 10.1039/FT9928800803
! Pokrovskii and Helgeson, Geochim. Cosmochim. Acta 61 (1997) 2175-2183 doi: 10.1016/S0016-7037(97)00070-7
Implicit NONE
Real(dp), INTENT(IN) :: tC     ! the temperature in degrees Celsius
Real(dp), INTENT(IN) :: pBar   ! the pressure in bar
Real(dp) :: epsH2O ! the dielectric constant
Real(dp) :: gf     ! the g-function (in units of Å) in the HKF model
Real(dp) :: p_Bar, t_C, tK, tr, eta, r_c, r_a, a1,a2,a3,a4,a5, c1,c2, omg, bg, bs, &
                    r_eff_c, r_eff_a, omgpt, f1T,f2T,f1P,f2P,f1PT,f2PT, nbg

if(tC < 0.D0 .or. tC > 1000.01D0) then
    write(*,*) """b_gamma_NaCl"": tC =",real(tC)," must be >=0 and < 1000 C"
    Call ErrStop
endif
if(pBar <= 0.D0 .or. pBar > 5000.01D0) then
    write(*,*) """b_gamma_NaCl"": pBar =",real(pBar)," must be >0 and < 5 kbar"
    Call ErrStop
endif
p_Bar = max(1.D0, pBar);
! if temperature = 0, set temperature to 0.01 C (tripple point of water)
if(abs(tC) < 0.001D0) then
    t_C = 0.01D0
else
    t_C = tC
endif
epsH2O = epsJN(t_C, p_Bar);
gf = 1.d+10 * gHKF(t_C, p_Bar);
if(abs(gf)<0.0004d0) gf=0.d0;
tK = t_C + 273.15D0
tr = 298.15D0
eta = 1.66027d5
r_c = 0.97D0; r_a = 1.81D0
a1 = 0.030056D0; a2 = -202.55D0; a3 = -2.9092D0; a4 = 20302D0; a5 = -0.206D0
c1 = -1.50D0;   c2 = 53300.D0
omg = 178650.D0
bg = -174.623D0; bs = 2.164D0
r_eff_c = r_c + (0.94D0+gf)
r_eff_a = r_a + gf
omgpt = eta*( (1.D0/r_eff_c) + (1.D0/r_eff_a) )
f1T = tK*log(tK/tr) - tK + tr
f2T = ((1.D0/(tK-228.D0))-(1.D0/(tr-228.D0))) * ((228.D0-tK)/228.D0) -(tK/(228.D0*228.D0)) &
                * log((tr*(tK-228.D0))/(tK*(tr-228.D0)))
f1P = p_Bar-1.D0
f2P = log((2600.D0+p_Bar)/(2600.D0+1.D0))
f1PT = (p_Bar-1.D0)/(tK-228.D0)
f2PT = (1.D0/(tK-228.D0))*log((2600.D0+p_Bar)/(2600.D0+1.D0))
    !
    !write(*,*) "r_eff_c = ",real(r_eff_c),", r_eff_a = ",real(r_eff_a);
    !write(*,*) "eps = ",real(epsH2O),", g = ",real(gf*10000.),", omgpt = ",real(omgpt);
    !write(*,*) "- bg + bs*(tK-298.15) = ",real(- bg + bs*(tK-tr));
    !write(*,*) "c1*( f1T ) = ",real(c1*( f1T ));
    !write(*,*) "c2*( f2T ) = ",real(c2*( f2T ));
    !write(*,*) "a3*( f1PT ) + a4*( f2PT )= ",real(a3*( f1PT )+ a4*( f2PT ));
    !
nbg = - bg + bs*(tK-tr) &
            - c1*( f1T ) - c2*( f2T ) &
            + a1*( f1P ) + a2*( f2P ) &
            + a3*( f1PT ) + a4*( f2PT ) &
            + a5*(omgpt*((1.D0/epsH2O)-1.D0) - omg*((1.D0/78.244D0)-1.D0) + (-5.799d-5)*omg*(tK-tr));
b_gamma = nbg/(2.D0*log(10.D0)*(1.98720426D0)*tK) ! gas constant in cal/(K mol)
RETURN
END FUNCTION b_gamma_NaCl

!----------------------------------------------------------------------------
SUBROUTINE FACTOR (nIons,Conc,lnF)
! Must be provided to calculate activity factors for aqueous species.
!  nIons = number of species (dimension of Conc and lnF) for which
!       activity coefficients are to be calculated
!  Conc = concentration of each species (size: 1 to nIons)
!  lnF = calculated natural logarith of the activity coefficient
!           corresponding to Conc (size: 1 to nIons)

! Define upper and lower bounds on the concentrations if concen-
! trations far from the equilibrium values can cause exponents
! out of the range allowed.

!---------------------
Implicit NONE
Integer, INTENT(IN) :: nIons
Real(dp), INTENT(IN)  :: Conc(nIons)
Real(dp), INTENT(OUT) :: lnF(nIons)
!Logical, SAVE :: msg !###
!-------------------
!This routine uses the the variables: NOLL, ionicStr, ionicStrCalc, Z
!-------------------
if(firstTime) call factorStart (nIons,Conc,lnF)

! --- Ideal solutions (activity coefficients = 1.):
! If the ionic strength is fixed and equal to zero, or if the model is not chosen,
! then the log of activity coefficients are zero, already assigned in "factorStart"
If(abs(ionicStr) < 1.D-35 .or. activityCoeffsModel < 0 .or. activityCoeffsModel > 2) RETURN;


! --- If the ionic strength is fixed, the activity coefficients for models
! that only depend on the ionic strength and on temperature-pressure
! (that is: Davies and HKF) have already been calculated in "factorStart"
If(.not.calculateIonicStr .and. (activityCoeffsModel == 0 .or. activityCoeffsModel == 2)) RETURN;

! --- At this point, either
!  - ionicStr > 0 (calculateIonicStr = false) and activityCoeffsModel = 1 (SIT model)
!  - ionicStr < 0 (calculateIonicStr = true) and activityCoeffsModel = 0,1, or 2

! Calculate ionicStrCalc, electricBalance and sumM
if(calculateIonicStr) call calcIonicStr (nIons, Conc)

! ionicStrCalc should now be >=0
rootI = 0.d0
if(ionicStrCalc > 0.d0) rootI=SQRT(ionicStrCalc)

!------------------------------------------------
! ---- Calculate activity coefficients
!------------------------------------------------
!           Davies Equation
!           ---------------
IF (activityCoeffsModel == 0) THEN
    call calcDavies (nIons,Conc,lnF)
ENDIF

!------------------------------------------------
!      Specific Ion Interaction Model
!      ------------------------------
IF (activityCoeffsModel == 1) THEN
    ! With the SIT the activity coefficients may change even if the
    ! ionic strength does not change: they depend on the composition.
    ! Hence, the values of lnf[] will NOT be calculated correctly
    ! in "calcSIT" unless the values of C[] are correct!
    call calcSIT (nIons,Conc,lnF)
ENDIF

!------------------------------------------------
!     Simplified HKF Model
!     --------------------
! H.C.Helgeson, D.H.Kirkham and G.C.Flowers, Amer.Jour.Sci. 281 (1981) 1249-1516
!  see also: H.C.Helgeson in "Chemistry and Geochemistry of Solutions at High
!  Temperatures and Pressures" D.T.Rickard and F.E.Wickman (eds), Pergamon
!  Press, Oxford, 1981, pp. 133-177, and E.H.Oelkers and H.C.Helgeson,
!  Geochim.Cosmochi.Acta 54 (1990) 727-738
IF (activityCoeffsModel == 2) THEN
    call calcHKF (nIons,Conc,lnF)
ENDIF

RETURN
END SUBROUTINE FACTOR

!----------------------------------------------------------------------------
SUBROUTINE calcIonicStr (nIons, Conc)
! Calculates 'ionicStrCalc', 'electricBalance' and 'sumM'.
! First the electrical balance is calculated.  To maintain an electrically
! neutral solution: If 'electricBalance' is larger than zero, then an inert
! anion (X-) is added. If 'electricBalance' is less than zero,  then  an
! inert cation (M+) is added to maintain an electrically neutral solution.
USE CHEM, ONLY : NOLL, Z
Implicit NONE
Integer, INTENT(IN) :: nIons
Real(dp), INTENT(IN)  :: Conc(nIons)
Real(dp) :: Ci,zi
Integer :: i

    call calcSumM(nIons, Conc) ! get electricBalance
    ionicStrCalc = 0.d0
    DO I=1,NIons
        if(NOLL(I)) Cycle
        Ci = max(0.d0, min(1.d+35, Conc(I)))
        if(Z(I) == 0) Cycle
        zi = real(Z(I),dp)
        ionicStrCalc = ionicStrCalc + zi*zi*Ci
    END DO
    ionicStrCalc=0.5D0 * (abs(electricBalance) + ionicStrCalc)

RETURN
END SUBROUTINE calcIonicStr

!----------------------------------------------------------------------------
SUBROUTINE calcSumM (nIons, Conc)
! Calculates 'electricBalance' and 'sumM'.
! First the electrical balance is calculated.
! To maintain an electrically neutral solution:
! If 'electricBalance' is larger than zero, then an inert
! anion (X-) is added. If 'electricBalance' is less than zero,  then  an
! inert cation (M+) is added to maintain an electrically neutral solution.
USE CHEM, ONLY : NOLL, Z
Implicit NONE
Integer, INTENT(IN) :: nIons
Real(dp), INTENT(IN)  :: Conc(nIons)
Real(dp) :: Ci,zi
Integer :: i

    sumM = 0.d0
    electricBalance = 0.d0
    DO I=1,NIons
        if(NOLL(I)) Cycle
        Ci = max(0.d0, min(1.d+35, Conc(I)))
        sumM = sumM + Ci
        if(Z(I) == 0) Cycle
        zi = real(Z(I),dp)
        electricBalance = electricBalance + zi*Ci
    END DO
    sumM = sumM + abs(electricBalance)

RETURN
END SUBROUTINE calcSumM

!----------------------------------------------------------------------------
SUBROUTINE calcDavies (nIons,Conc,lnF)
USE CHEM, ONLY : DBG, IOUT, NOLL, Z, jWater
Implicit NONE
Integer, INTENT(IN) :: nIons
Real(dp), INTENT(IN)  :: Conc(nIons)
Real(dp), INTENT(OUT) :: lnF(nIons)
Real(dp) :: sig, w, phiDH, down, logf, zz
Integer :: i
If(DBG >= 3) Write(IOUT,'("calcDavies(",I0,"), activityCoeffsModel =",I0,", calculateIonicStr=",L1,", ionicStrCalc=",1PE8.2)') &
                            nIons,activityCoeffsModel, calculateIonicStr, ionicStrCalc
if(ionicStrCalc <= 0.d0) return
if(.not.calculateIonicStr) call calcSumM(nIons, Conc) ! get sumM
! both the number "1.0" and "DAVIES" are assumed to be temperature-independent
down = 1.0d0 + rootI
w = - Agamma * ((rootI/down) - DAVIES*ionicStrCalc)
DO I=1,nIons
    if(NOLL(I)) then
        lnF(I) = 0.d0
        cycle
    endif
    if(Z(I) /= 0) then
        zz = real((Z(I)*Z(I)), dp)
        logf = zz * w
    else
        zz = 1.d0
        logf = 0.1d0 * ionicStrCalc
    end if
    logf = max(-MAX_LOG_G*zz, min(logf, MAX_LOG_G*zz))
    lnF(I)= ln10 * logf
END DO

IF(jWater <= 0) RETURN

!Calculate osmotic coeff. and activity of water
w = 0.d0
if(ionicStrCalc > 1.d-50) w = (ionicStrCalc*rootI) ! = I^3/2
sig = sigma(rootI)
! phiDH = Debye-Huckel term for phi
!  2 ln(10) / 3 = 1.5350567286627
phiDH= 1.5350567286627d0 * Agamma * w * sig
phi= 1.d0
IF(sumM > 1.D-15) phi= 1.d0 - (phiDH/sumM) + ln10*(Agamma*DAVIES)*(ionicStrCalc*ionicStrCalc)/sumM
! Calculate Water Activity
log10aH2O = -phi * sumM / (ln10 * molH2Oin1kg)
log10aH2O = max(-MAX_LOGAH2O,min(MAX_LOGAH2O,log10aH2O))
! lnF(jWater) = ln10 * log10aH2O ! ##
RETURN

END SUBROUTINE calcDavies

!----------------------------------------------------------------------------
SUBROUTINE calcHKF (nIons,Conc,lnF)
! Calculates activity coefficients with:
!   log f(i) = - (Aγ(T) Z(i)² √I)/(1 + r B(T) √I) + Γ + (bγ(T) I)
! where: Γ = -log(1+0.0180153 Σm).</pre>
! These eqns. are an approximation to the HKF model
! (Helgeson, Kirkham and Flowers),
! see Eqs.121,165-167,297,298 in: H.C.Helgeson, D.H.Kirkham and G.C.Flowers,
! Amer.Jour.Sci., 281 (1981) 1249-1516,  Eqs.22 & 23 in: E.H.Oelkers and
! H.C.Helgeson, Geochim.Cosmochim.Acta, 54 (1990) 727-738, etc.
USE CHEM, ONLY : DBG, IOUT, NOLL, Z, jWater
Implicit NONE
Integer, INTENT(IN) :: nIons
Double Precision, INTENT(IN)  :: Conc(nIons)
Double Precision, INTENT(OUT) :: lnF(nIons)
Real(dp) :: gamma, sig, w, phiDH, down, logf, zz
Integer :: i
!--- Note that "sumM" is not used in this verions of HKF
If(DBG >= 3) Write(IOUT,'("calcHKF(",I0,"), activityCoeffsModel =",I0,", calculateIonicStr=",L1,", ionicStrCalc=",1PE8.2)') &
                            nIons,activityCoeffsModel, calculateIonicStr, ionicStrCalc
if(ionicStrCalc <= 0.d0) return
if(.not.calculateIonicStr) call calcSumM(nIons, Conc) ! get sumM and electricBalance
down = 1.d0 + rB*rootI
w = - Agamma * (rootI/down)
gamma = log10(1.d0+(0.0180153d0*ionicStrCalc)) ! sumM));
DO I=1,nIons
    if(NOLL(I)) then
        lnF(I) = 0.D0
        cycle
    endif
    if(Z(I) /= 0) then
        zz = real(Z(I)*Z(I), dp)
        logf = zz * w - gamma + bgi * ionicStrCalc
    else
        zz = 1.d0
        logf = 0.1d0 * ionicStrCalc
    end if
    logf = max(-MAX_LOG_G*zz, min(logf, MAX_LOG_G*zz))
    lnF(I)= ln10 * logf
END DO

IF(jWater <= 0) RETURN

! Calculate osmotic coeff. and activity of water
w = 0.d0
if(ionicStrCalc > 1.d-50) w = (ionicStrCalc*rootI) ! = I^3/2
sig = sigma(rB*rootI)
! phiDH = Debye-Huckel term for phi
!  2 ln(10) / 3 = 1.5350567286627
phiDH= 1.5350567286627d0 * Agamma * w * sig
phi = 1.d0
IF (ionicStrCalc > 1.D-25) & ! sumM
    phi = (ln10*gamma/(0.0180153d0*ionicStrCalc)) & ! sumM)
    - (phiDH/ionicStrCalc) & ! sumM)
    + (ln10 * bgi * 0.5d0 * ionicStrCalc)
! Calculate Water Activity
log10aH2O = -(phi * ionicStrCalc) & ! sumM
        / (ln10 * molH2Oin1kg)
log10aH2O = max(-MAX_LOGAH2O, min(MAX_LOGAH2O, log10aH2O))
! lnF(jWater) = ln10 * log10aH2O ! ##
RETURN

END SUBROUTINE calcHKF

!----------------------------------------------------------------------------
SUBROUTINE calcSIT (nIons,Conc,lnF)
USE CHEM, ONLY : DBG, IOUT, NOLL, Z, jWater !, DBG, IOUT
USE SIT, ONLY : getEps
Implicit NONE
Integer, INTENT(IN) :: nIons
Real(dp), INTENT(IN)  :: Conc(nIons)
Real(dp), INTENT(OUT) :: lnF(nIons)
Real(dp) :: w, DH, down, eps, sumEpsM, elBal, Ci, Cj, zz, logf, sig, phiDH, sumPrd, sumPrd_i
Integer :: i,j
If(DBG >= 3 .or. fDbg) &
    Write(IOUT,'("calcSIT(",I0,"), activityCoeffsModel =",I0,", calculateIonicStr=",L1,", ionicStrCalc=",1PE8.2)') &
                            nIons,activityCoeffsModel, calculateIonicStr, ionicStrCalc
if(ionicStrCalc <= 0.d0) return
if(.not.calculateIonicStr) call calcSumM(nIons, Conc) ! get electricBalance and sumM
!Calculate Debye-Huckel term
down= 1.d0 + (rB * rootI)
DH= -AGAMMA * rootI / down
elBal = abs(electricBalance)
! --- Calculate the individual ionic activity coefficients (lnF)
! For neutral species this program uses ε(i,MX)*[M] + ε(i,MX)*[X]
! As a consequence:
!   for a MX electrolyte: ε(i,M)*[MX] + ε(i,X)*[MX]
!   for a M2X (or MX2) electrolyte: ε(i,M)*2*[M2X] + ε(i,X)*[M2X]
!   (or ε(i,M)*[MX2] + ε(i,X)*2*[MX2])
! In the SIT-file you must enter ε = ε(i,MX)/2 (or ε = ε(i,M2X)/3)
DO i=1,nIons
    lnF(i)=0.d0
    if(NOLL(i)) Cycle
    sumEpsM = 0.d0
    if(Z(i) == 0) then
        eps = getEps(i,i)
        Ci = max(0.d0, min(MAX_CONC, Conc(i)))
        sumEpsM = eps * Ci
    endif
    if(elBal > 1.d-10) then
        eps = 0.d0
        if(electricBalance < 0.d0 .and. Z(i) <= 0) eps = getEps(nIons+1,i) ! add Na+
        if(electricBalance > 0.d0 .and. Z(i) >= 0) eps = getEps(i,nIons+2) ! add Cl-
        sumEpsM =  sumEpsM + eps * max(0.d0, min(MAX_CONC, elBal))
    endif
    do j=1,nIons
        if(NOLL(j) .or. (j == i) .or. ((Z(i)*Z(j)) > 0)) Cycle !j
        eps = getEps(i,j)
        Cj = max(0.d0, min(MAX_CONC, Conc(j)))
        sumEpsM = sumEpsM + eps * Cj
    enddo
    zz = real((Z(i)*Z(i)),dp)
    logf = zz*DH + sumEpsM
    zz = max(1.d0,zz) ! make sure zz is not zero
    logf = max(-MAX_LOG_G*zz, min(logf, MAX_LOG_G*zz))
    lnF(i) = ln10 * logf
END DO !i

! --- Calculate Osmotic Coeff.
if(jWater <= 0) RETURN

! phiDH = Debye-Huckel term
w = 0.d0
if(ionicStrCalc > 1.d-50) w = ionicStrCalc * rootI; ! = I^3/2
sig = sigma(rB*rootI)
!  2 ln(10) / 3 = 1.5350567286627
phiDH= 1.5350567286627d0 * Agamma * w * sig
! Calculate the sum of ions and the sum of products of conc. times epsilon
sumPrd=0.d0
! loop through cations and neutral species
DO i=1,nIons
    if(NOLL(i)) cycle
    Ci = max(0.d0, min(MAX_CONC, Conc(i)))
    sumPrd_i = 0.d0
    if(Z(i) > 0) then       !--- i = cation ---
        do j=1,nIons
            if(NOLL(j)) cycle
            if(j == i)  cycle
            if(Z(j) > 0) cycle ! take cations x (anions and neutral)
            Cj = max(0.d0, min(MAX_CONC, Conc(j)))
            eps = getEps(i,j)
            sumPrd_i = sumPrd_i + (eps * Cj)
        enddo !j
        ! if electricBalance > 0, then an inert anion (X-) is added.
        ! if electricBalance < 0,  then  an inert cation (M+) is added.
        if(electricBalance > 1.d-10) then
            eps = getEps(i,nIons+2)
            sumPrd_i = sumPrd_i + (eps * electricBalance)
        endif
        sumPrd_i = sumPrd_i * Ci
    else  if(Z(i) < 0) then !--- i = anion ---
        do j=1,nIons
            if(NOLL(j)) cycle
            if(j == i)  cycle
            if(Z(j) /= 0) cycle ! take anions x neutral
            Cj = max(0.d0, min(MAX_CONC, Conc(j)))
            eps = getEps(i,j)
            sumPrd_i = sumPrd_i + (eps * Cj)
        enddo !J
        sumPrd_i = sumPrd_i * Ci
    else                    !--- i = neutral species ---
        do j=1,nIons
            if(NOLL(j)) cycle
            if(Z(j) /= 0) cycle ! take neutral x neutral
            Cj = max(0.d0, min(MAX_CONC, Conc(j)))
            eps = getEps(i,j)
            sumPrd_i = sumPrd_i + (eps * Cj)
        enddo !J
        if(elBal > 1.d-10) then
            ! presumably eps(n,Na+) = eps(n,Cl-)
            eps = getEps(i,nIons+1) ! interaction with Na+ (elec.balance)
            sumPrd_i = sumPrd_i + (eps * elBal)
        endif
        sumPrd_i = 0.5d0 * sumPrd_i * Ci
    endif ! Z(i)?
    sumPrd = sumPrd + sumPrd_i
END DO !i

! The remaining cation is Na+ "added" for electic balance
! If electricBalance > 0, then an inert anion (X-) is added.
! If electricBalance < 0,  then  an inert cation (M+) is added.
if(electricBalance < -1.D-10) then
    do i=1,nIons
        if(Z(i) >= 0 .or. NOLL(i)) Cycle
        eps = getEps(nIons+1,i)
        Ci = max(0.d0, min(MAX_CONC, Conc(i)))
        sumPrd = sumPrd - (eps * electricBalance * Ci)
    enddo !i
endif

phi = 1.d0
if(sumM > 1.D-15) phi = 1.d0 - (phiDH/sumM) + (ln10*sumPrd)/sumM
! Calculate Water Activity.
log10aH2O = - phi * sumM / (ln10 * molH2Oin1kg)
log10aH2O = max(-MAX_LOGAH2O, min(MAX_LOGAH2O, log10aH2O))
! lnF(jWater) = ln10 * log10aH2O ! ##
RETURN

END SUBROUTINE calcSIT

!-----------------------------------------------------------------------------
FUNCTION sigma(x) RESULT(s)
Implicit NONE
Real(dp), INTENT(IN) :: x
Real(dp) :: s
if(x < 0.d0) then
    s = 0.d0
    RETURN
else if(x < 1.d-10) then
    s = 1.d0
    RETURN
endif
s = (3.d0/(x*x*x)) * ((1.d0+x) - 2.d0*log(1.d0+x) - (1.d0/(1.d0+x)));
RETURN
END FUNCTION sigma

!----------------------------------------------------------------------------
SUBROUTINE ErrStop ! (PRIVATE)
    Implicit NONE
    Integer :: ierr, i
    WRITE (*, '("Press Enter to continue ...")', Advance='No')
    READ (*,'(I5)',iostat=ierr) i
    STOP 1
END SUBROUTINE ErrStop

END MODULE FACTOR_Module
