PROGRAM TEST_HaltaF
USE CHEM
USE HALTAF_Module
USE FACTOR_Module, ONLY : activityCoeffsModel, ionicStr, Temperature
Implicit NONE
!------------------------------------------------------------------
! Test of HaltaFall routine.
! Four different tests: 1) a chemical system without solids (only aqueous
! solution); 2) a system with one solid reaction product; 3) a system
! with several possible solid products; and 4) a system where one of the
! chemical components is solid, and there are several possible solid reaction
! products as well.
!
! Equilibrium constants correspond to the Wateq4F database distributed
! with PhreeqC at https://www.usgs.gov/software/phreeqc-version-3
!------------------------------------------------------------------
! @author Ignasi Puigdomenech
Integer :: I,J,K,N,ierr
Real(Kind(1.D0)) :: pH, pe
    IOUT=21;
    OPEN(UNIT=IOUT,FILE='Test_Halta.out')
    WRITE(*,20) "output file: 'Test_Halta.out'"
    20 FORMAT( &
    '========================',/, &
    'Test of subroutine HALTA',/, &
    '========================',:,/,A)
    21 FORMAT('Nr. Components =',I4," of which ",I2," are solid.",/, &
    'Nr. Aqueous Species =',I5,/,'Nr. Solid Phases =',I3,/)
    22 FORMAT(/,'Components: Input Tot.Concs. =',1P,6(:/5X,4G14.6))
    23 FORMAT('RESULTS:',:,10X,' errFlags = ',A)

! -----------------------------------------------------------------------
! Test 1: a chemical system without solids (only aqueous solution)
! -------
    WRITE(*,'("Test 1 (no solids)")')
    WRITE(IOUT,20)
    WRITE(IOUT,'("Test 1 (no solids)")')
    ! 4 chemical components are: H+, CO3-2, Na+, H2O
    ! 11 aqueous species are:     H+, CO3-2, Na+, H2O
    !    CO2, H2CO3, HCO3-, NaCO3-, NaHCO3, NaOH, OH-
    ! No solids
    NA = 4 !chemical components
    MSOL = 0 ! solids
    MS = 11  !species = components + reaction products
    SOLIDC = 0  ! how many components are solid phases
    ! allocate memory arrays
    Call CHEM_MEM_ALLOC
    Call HALTAF_MEM_ALLOC
    DO I=1,MS
      NOLL(I)=.FALSE.
    EndDo
    NOLL(4) = .TRUE. ! the mass balance of H2O is not included
! Data for chemical components
    IDENTC(1) = "H+";  IDENTC(2) = "CO3-2";  IDENTC(3) = "Na+";  IDENTC(4) ="H2O"
    IDENT(1:NA) = IDENTC(1:NA)
    Z(1:NA)= (/ +1, -2, +1, 0 /) ! electric charges
! Data for reaction products
    A(1,1:4) = (/ 2.,1.,0.,-1./); LBETA(1) = 16.681; IDENT(NA+1) = "CO2";    Z(NA+1) =  0
    A(2,1:4) = (/ 2.,1.,0., 0./); LBETA(2) = 13.859; IDENT(NA+2) = "H2CO3";  Z(NA+2) =  0
    A(3,1:4) = (/ 1.,1.,0., 0./); LBETA(3) = 10.329; IDENT(NA+3) = "HCO3-";  Z(NA+3) = -1
    A(4,1:4) = (/ 0.,1.,1., 0./); LBETA(4) =  1.27;  IDENT(NA+4) = "NaCO3-"; Z(NA+4) = -1
    A(5,1:4) = (/ 1.,1.,1., 0./); LBETA(5) = 10.079; IDENT(NA+5) = "NaHCO3"; Z(NA+5) =  0
    A(6,1:4) = (/-1.,0.,1., 1./); LBETA(6) =-14.4;   IDENT(NA+6) = "NaOH";   Z(NA+6) =  0
    A(7,1:4) = (/-1.,0.,0., 1./); LBETA(7) =-14.00;  IDENT(NA+7) = "OH-";    Z(NA+7) = -1
! Concentrations
    DO I = 1, NA
      KH(I)=1       ! only total concentrations given (KH=1)
    End Do
    KH(4) = 2 ! log(activity) given for H2O
    ! total concentrations for the components
    TOT(1)= 0.01 ! H+
    TOT(2)= 0.01 ! CO3-2
    TOT(3)= 0.01 ! Na+
    TOT(4)=0.D0; LOGA(4) = 0. ! H2O
    TOL=1.D-5  ! tolerance when solving mass-balance equation
    JWATER = 4 ! H2O
    DBG=1; activityCoeffsModel = 0; ionicStr = -1.
    WRITE(IOUT,21) NA,SOLIDC,NIon,MSOL
    WRITE(IOUT,'("Find pH of an aqueous solution with 0.01 mol NaHCO3.")')
    WRITE(IOUT,22) (TOT(I),I=1,NA)
    Write(IOUT,'(" tol=",1pg10.3)') tol
    Call printChemSystem (IOUT)
    CONT=.FALSE.
    CALL HaltaCalc
    WRITE(IOUT,23) errFlagsToString()
    pH=-LOGA(1)
    WRITE(IOUT,'(5X,"pH=",F6.3," (should be 8.222)")')  pH
    Call printConcs (IOUT)
Call CHEM_MEM_FREE
Call HALTAF_MEM_FREE

! -----------------------------------------------------------------------
! Test 2: a system with one solid reaction product
! -------
    WRITE(*,'("Test 2 (one solid)")')
    WRITE(IOUT,20)
    WRITE(IOUT,'("Test 2 (one solid)")')
    ! 6 chemical components are: H+, CO3-2, Ca+2, H2O, Na+, Cl-,
    ! 14 aqueous species are:     H+, CO3-2, Ca+2, H2O, Na+, Cl-,
    !    CaCO3, CaHCO3+, CaOH+, CO2, HCO3-, OH-, NaCO3-, NaHCO3
    ! 1 solid is: CaCO3(cr)
    NA = 6 !chemical components
    MSOL = 1 ! solids
    MS = 15  !species = components + reaction products
    SOLIDC = 0  ! how many components are solid phases
    ! allocate memory arrays
    Call CHEM_MEM_ALLOC
    Call HALTAF_MEM_ALLOC
    DO I=1,MS
      NOLL(I)=.FALSE.
    EndDo
    NOLL(4) = .TRUE. ! the mass balance of H2O is not included
! Data for components
    IDENTC(1) = "H+";  IDENTC(2) = "CO3-2";
    IDENTC(3) = "Ca+2";IDENTC(4) = "H2O";
    IDENTC(5) = "Na+"; IDENTC(6) = "Cl-"
    IDENT(1:NA) = IDENTC(1:NA)
    Z(1:NA)= (/ +1, -2, +2, 0, +1, -1 /) ! electric charges
! Data for reaction products
    A(1,1:6) = (/ 2.,1.,0.,-1.,0.,0./); LBETA(1) = 16.681; IDENT(NA+1) = "CO2";     Z(NA+1) =  0
    A(2,1:6) = (/ 1.,1.,0., 0.,0.,0./); LBETA(2) = 10.329; IDENT(NA+2) = "HCO3-";   Z(NA+2) = -1
    A(3,1:6) = (/ 0.,1.,1., 0.,0.,0./); LBETA(3) =  3.224; IDENT(NA+3) = "CaCO3";   Z(NA+3) =  0
    A(4,1:6) = (/ 1.,1.,1., 0.,0.,0./); LBETA(4) = 11.435; IDENT(NA+4) = "CaHCO3+"; Z(NA+4) = +1
    A(5,1:6) = (/-1.,0.,1., 1.,0.,0./); LBETA(5) =-12.78;  IDENT(NA+5) = "CaOH+";   Z(NA+5) = +1
    A(6,1:6) = (/-1.,0.,0., 1.,0.,0./); LBETA(6) =-14.00;  IDENT(NA+6) = "OH-";     Z(NA+6) = -1
    A(7,1:6) = (/ 0.,1.,0., 0.,1.,0./); LBETA(7) =  1.27;  IDENT(NA+7) = "NaCO3-";  Z(NA+7) = -1
    A(8,1:6) = (/ 1.,1.,0., 0.,1.,0./); LBETA(8) = 10.079; IDENT(NA+8) = "NaHCO3";  Z(NA+8) =  0
    A(9,1:6) = (/ 0.,1.,1., 0.,0.,0./); LBETA(9) =  8.48;  IDENT(NA+9) = "CaCO3(cr)";
! Concentrations
    DO I = 1, NA
      KH(I)=1       ! only total concentrations given (KH=1)
    End Do
    KH(4) = 2 ! log activity given for H2O
    ! total concentrations for the components
    TOT(1)= 0.   ! H+
    TOT(2)= 0.01 ! CO3-2
    TOT(3)= 0.01 ! Ca+2
    TOT(4)=0.D0; LOGA(4) = 0. ! H2O
    TOT(5)= 0.02 ! Na+
    TOT(6)= 0.02 ! Cl-
    TOL=1.D-6  ! tolerance when solving mass-balance equation
    JWATER = 4 ! H2O
    DBG=1; activityCoeffsModel = 0; ionicStr = -1.
    WRITE(IOUT,21) NA,SOLIDC,NIon,MSOL
    WRITE(IOUT,'("Find pH and calcite precipitated of an aqueous solution",/, &
    "with 0.01 mol CaCl2 and 0.01 mol Na2CO3")')
    WRITE(IOUT,22) (TOT(I),I=1,NA)
    Write(IOUT,'(" tol=",1pg10.3)') tol
    Call printChemSystem (IOUT)
    CONT=.FALSE.
    CALL HaltaCalc
    WRITE(IOUT,23) errFlagsToString()
    pH=-LOGA(1)
    WRITE(IOUT,'(5X,"pH=",F6.3," (should be 9.932)")')  pH
    WRITE(IOUT,'(5X,"calcite precipitated: ",F8.6," mol (should be 0.00982)")')  C(NION+1)
    Call printConcs (IOUT)
Call CHEM_MEM_FREE
Call HALTAF_MEM_FREE

! -----------------------------------------------------------------------
! Test 3: a system with several possible solid products
! -------
    WRITE(*,'("Test 3 (several solids)")')
    WRITE(IOUT,20)
    WRITE(IOUT,'("Test 3 (several solids)")')
    ! 6 chemical components are: H+, CO3-2, Ca+2, H2O, Na+, Cl-,
    ! 14 aqueous species are:     H+, CO3-2, Ca+2, H2O, Na+, Cl-,
    !    CaCO3, CaHCO3+, CaOH+, CO2, HCO3-, OH-, NaCO3-, NaHCO3
    ! 4 solid are: CaCO3(cr), CaCO3(s), Ca(OH)2(s), CaO(s)
    !              (both calcite and aragonite are included)
    NA = 6 !chemical components
    MSOL = 4 ! solids
    MS = 18  !species = components + reaction products
    SOLIDC = 0  ! how many components are solid phases
    ! allocate memory arrays
    Call CHEM_MEM_ALLOC
    Call HALTAF_MEM_ALLOC
    DO I=1,MS
      NOLL(I)=.FALSE.
    EndDo
    NOLL(4) = .TRUE. ! the mass balance of H2O is not included
! Data for components
    IDENTC(1) = "H+";  IDENTC(2) = "CO3-2";
    IDENTC(3) = "Ca+2";IDENTC(4) = "H2O";
    IDENTC(5) = "Na+"; IDENTC(6) = "Cl-"
    IDENT(1:NA) = IDENTC(1:NA)
    Z(1:NA)= (/ +1, -2, +2, 0, +1, -1 /) ! electric charges
! Data for reaction products
    A(1,1:6) = (/ 2.d0,1.d0,0.d0,-1.d0,0.d0,0.d0/); LBETA(1) = 16.681d0; IDENT(NA+1) = "CO2";     Z(NA+1) =  0
    A(2,1:6) = (/ 1.d0,1.d0,0.d0, 0.d0,0.d0,0.d0/); LBETA(2) = 10.329d0; IDENT(NA+2) = "HCO3-";   Z(NA+2) = -1
    A(3,1:6) = (/ 0.d0,1.d0,1.d0, 0.d0,0.d0,0.d0/); LBETA(3) =  3.224d0;  IDENT(NA+3) = "CaCO3";   Z(NA+3) =  0
    A(4,1:6) = (/ 1.d0,1.d0,1.d0, 0.d0,0.d0,0.d0/); LBETA(4) = 11.435d0; IDENT(NA+4) = "CaHCO3+"; Z(NA+4) = +1
    A(5,1:6) = (/-1.d0,0.d0,1.d0, 1.d0,0.d0,0.d0/); LBETA(5) =-12.78d0;   IDENT(NA+5) = "CaOH+";   Z(NA+5) = +1
    A(6,1:6) = (/-1.d0,0.d0,0.d0, 1.d0,0.d0,0.d0/); LBETA(6) =-14.00d0;  IDENT(NA+6) = "OH-";     Z(NA+6) = -1
    A(7,1:6) = (/ 0.d0,1.d0,0.d0, 0.d0,1.d0,0.d0/); LBETA(7) =  1.27d0;  IDENT(NA+7) = "NaCO3-"; Z(NA+7) = -1
    A(8,1:6) = (/ 1.d0,1.d0,0.d0, 0.d0,1.d0,0.d0/); LBETA(8) = 10.079d0;IDENT(NA+8) = "NaHCO3"; Z(NA+8) =  0
    A(9,1:6) = (/ 0.d0,1.d0,1.d0, 0.d0,0.d0,0.d0/); LBETA(9) =  8.48d0;  IDENT(NA+9) = "Calcite";
    A(10,1:6)= (/ 0.d0,1.d0,1.d0, 0.d0,0.d0,0.d0/); LBETA(10)=  8.336d0;IDENT(NA+10) = "Aragonite";
    A(11,1:6)= (/-2.d0,0.d0,1.d0, 2.d0,0.d0,0.d0/); LBETA(11)=-22.8d0;  IDENT(NA+11) = "Ca(OH)2(s)";
    A(12,1:6)= (/-2.d0,0.d0,1.d0, 1.d0,0.d0,0.d0/); LBETA(12)=-32.58d0; IDENT(NA+12) = "CaO(s)";
! Concentrations
    DO I = 1, NA
      KH(I)=1       ! only total concentrations given (KH=1)
    End Do
    KH(4) = 2 ! log activity given for H2O
    ! total concentrations for the components
    TOT(1)= -0.1d0   ! H+
    TOT(2)=  0.001d0 ! CO3-2
    TOT(3)=  0.01d0  ! Ca+2
    TOT(4)=0.d0; LOGA(4)= 0.d0    ! H2O
    TOT(5)=  0.102d0 ! Na+
    TOT(6)=  0.02d0  ! Cl-
    TOL=1.D-6  ! tolerance when solving mass-balance equation
    JWATER = 4 ! H2O
    DBG=1; activityCoeffsModel = 0; ionicStr = -1.
    WRITE(IOUT,21) NA,SOLIDC,NIon,MSOL
    WRITE(IOUT,'("Find pH and calcite precipitated of an aqueous solution",/, &
    "with 0.01 mol CaCl2, 0.1 NaOH and 0.001 Na2CO3")')
    WRITE(IOUT,22) (TOT(I),I=1,NA)
    Write(IOUT,'(" tol=",1pg10.3)') tol
    Call printChemSystem (IOUT)
    CONT=.false.
    CALL HaltaCalc
    WRITE(IOUT,23) errFlagsToString()
    pH=-LOGA(1)
    WRITE(IOUT,'(5X,"pH=",F6.3," (should be 12.849)")')  pH
    WRITE(IOUT,'(5X,"calcite precipitated: ",F8.6," mol (should be 0.000982)")')  C(NION+1)
    WRITE(IOUT,'(5X,"Ca(OH)2 precipitated: ",F8.6," mol (should be 0.00354)")')  C(NION+3)
    Call printConcs (IOUT)
    !FLUSH(IOUT)
Call CHEM_MEM_FREE
Call HALTAF_MEM_FREE

! -----------------------------------------------------------------------
! Test 4: a system where one of the chemical components is solid,
!         and there are several possible solid reaction products as well
! -------
    WRITE(*,'("Test 4 (one solid component and several solid reaction products)")')
    WRITE(IOUT,20)
    WRITE(IOUT,'("Test 4 (one solid component and several solid reaction products)")')
    ! 4 chemical components: H+, Fe+2, Cl-, Fe(s)
    ! 17 aqueous species:     H+, Fe+2, Cl-
    !    e-, FeOH+, Fe(OH)2, Fe(OH)3-, FeCl+,
    !    Fe+3, FeOH+2, Fe(OH)2+, Fe(OH)3, Fe(OH)4-, Fe2(OH)2+4,
    !    OH-, O2, H2
    ! 3 solid reaction products   Fe(OH)2(s), Fe3O4(s), FeOOH(s)
    ! because a component is a solid, an extra solid is added:  Fe(s)
    NA = 4 !chemical components (both soluble and solids)
    MS = 25  !species = components (soluble+solid) + complexes (soluble+solid)
    SOLIDC = 1  ! how many components are solid phases
    MS = MS + SOLIDC ! add extra solids for each solid component
    MSOL = 4 !solids = 3 solid reaction products + 1 solid component
    ! allocate memory arrays
    Call CHEM_MEM_ALLOC
    Call HALTAF_MEM_ALLOC
    DO I=1,MS
      NOLL(I)=.FALSE.
    EndDo
   ! Note: all soluble complexes must be given first, followed by all solids.
   ! Note: Here the solids are arranged in the following order:
   ! all solid reaction products, followed by all solid components.
! Data for components
    IDENTC(1) = "H+";  IDENTC(2) = "Fe+2";
    IDENTC(3) = "Cl-"; IDENTC(4) = "Fe(s)"
    IDENT(1:NA) = IDENTC(1:NA)
    DO I = 1, SOLIDC !for solid components:
        J = MS -NA - SOLIDC + I ! j= (Ms-Na)-solidC + 1 ...(Ms-Na)
        K = NA - SOLIDC + I     ! k = (Na-solidC)+1...Na
        !create an extra solid complex
        NOLL(K)=.true. ! solid components are not aqueous species
        LBETA(J)=0.      ! equilibrium contstant of formation = 1.
        IDENT(NA+j) = IDENTC(K)
        DO N = 1, NA    ! set stoichiometry for solid components
            A(J,N)=0.D0
            if(N == K) A(J,N)=1.D0
        End Do
    End Do
    Z(1:NA)= (/ +1, +2, -1, 0 /) ! electric charges
! Data for reaction products
    A( 1,1:4) = (/ 0.,-0.5, 0., 0.5/);  LBETA( 1)= 7.98;   IDENT(NA+1)="e-";         Z(NA+1) = -1
    A( 2,1:4) = (/-1., 1.,  0., 0. /);  LBETA( 2)=-9.5;    IDENT(NA+2)="FeOH+";      Z(NA+2) = +1
    A( 3,1:4) = (/-2., 1.,  0., 0. /);  LBETA( 3)=-20.57;  IDENT(NA+3)="Fe(OH)2";    Z(NA+3) =  0
    A( 4,1:4) = (/-3., 1.,  0., 0. /);  LBETA( 4)=-31.;    IDENT(NA+4)="Fe(OH)3-";   Z(NA+4) = -1
    A( 5,1:4) = (/ 0., 1.,  1., 0. /);  LBETA( 5)= 0.14;   IDENT(NA+5)="FeCl+";      Z(NA+5) = +1
    A( 6,1:4) = (/ 0., 1.5, 0.,-0.5/);  LBETA( 6)=-21.;    IDENT(NA+6)="Fe+3";       Z(NA+6) = +3
    A( 7,1:4) = (/-1., 1.5, 0.,-0.5/);  LBETA( 7)=-23.19;  IDENT(NA+7)="FeOH+2";     Z(NA+7) = +2
    A( 8,1:4) = (/-2., 1.5, 0.,-0.5/);  LBETA( 8)=-26.67;  IDENT(NA+8)="Fe(OH)2+";   Z(NA+8) = +1
    A( 9,1:4) = (/-3., 1.5, 0.,-0.5/);  LBETA( 9)=-33.56;  IDENT(NA+9)="Fe(OH)3";    Z(NA+9) =  0
    A(10,1:4) = (/-4., 1.5, 0.,-0.5/);  LBETA(10)=-42.6;   IDENT(NA+10)="Fe(OH)4-";  Z(NA+10)= -1
    A(11,1:4) = (/-2., 3.,  0.,-1. /);  LBETA(11)=-44.95;  IDENT(NA+11)="Fe2(OH)2+4";Z(NA+11)= +4
    A(12,1:4) = (/ 0., 1.5, 1.,-0.5/);  LBETA(12)=-19.52;  IDENT(NA+12)="FeCl+2";    Z(NA+12)= +2
    A(13,1:4) = (/ 0., 1.5, 2.,-0.5/);  LBETA(13)=-18.87;  IDENT(NA+13)="FeCl2+";    Z(NA+13)= +1
    A(14,1:4) = (/ 0., 1.5, 3.,-0.5/);  LBETA(14)=-19.87;  IDENT(NA+14)="FeCl3";     Z(NA+14)=  0
    A(15,1:4) = (/ 0., 1.5, 4.,-0.5/);  LBETA(15)=-21.98;  IDENT(NA+15)="FeCl4-";    Z(NA+15)= -1
    A(16,1:4) = (/-1., 0.,  0., 0. /);  LBETA(16)=-14.00;  IDENT(NA+16)="OH-";       Z(NA+16)= -1
    A(17,1:4) = (/ 2.,-1.,  0., 1. /);  LBETA(17)= 12.815; IDENT(NA+17)="H2";        Z(NA+17)=  0
    A(18,1:4) = (/-4., 2.,  0.,-2. /);  LBETA(18)=-117.784;IDENT(NA+18)="O2";        Z(NA+18)=  0
    A(19,1:4) = (/-2., 1.,  0., 0. /);  LBETA(19)=-12.76;  IDENT(NA+19)="Fe(OH)2(s)"
    A(20,1:4) = (/-8., 4.,  0.,-1. /);  LBETA(20)=-45.737; IDENT(NA+20)="Fe3O4(s)"
    A(21,1:4) = (/-3., 1.5, 0.,-0.5/);  LBETA(21)=-20.0;   IDENT(NA+21)="FeOOH(s)"
    NOLL(5)=.true. ! do not consider the formation of "e-"
! Concentrations
    DO I = 1, NA
      KH(I)=1       ! only total concentrations given (KH=1)
    End Do
    ! total concentrations for the components
    TOT(1)= 1.D-5  ! H+
    TOT(2)= 0.     ! Fe+2
    TOT(3)= 1.D-5  ! Cl-
    TOT(4)= 0.01   ! Fe(s)
    TOL=1.D-6  ! tolerance when solving mass-balance equation
    JWATER = 0 ! no H2O
    DBG=1; activityCoeffsModel = 0; ionicStr = -1.
    Temperature = 25. ! needed to calculate the redox potential (Eh)
    WRITE(IOUT,21) NA,SOLIDC,NIon,MSOL
    WRITE(IOUT,'("Find pH and redox potential and solids that precipitate",/1X, &
      "for an aqueous solution with 1E-5 mol HCl and 0.01 mol Fe(metal)")')
    ! 3Fe(s) + 4H2O = Fe3O4 + 4H2
    WRITE(IOUT,22) (TOT(I),I=1,NA)
    Write(IOUT,'(" tol=",1pg10.3)') tol
    Call printChemSystem (IOUT)
    CONT=.false.
    CALL HaltaCalc
    WRITE(IOUT,23) errFlagsToString()
    pH=-LOGA(1); pe=-LOGA(5)
    WRITE(IOUT,'(5X,"pH=",F6.3," (should be 7.82),  pe=",F7.3," (should be -8.45)")')  pH,pe
    WRITE(IOUT,'(5X,"magnetite precipitated: ",F8.6," mol (should be 0.00333)")')  C(NION+2)
    Call printConcs (IOUT)
Call CHEM_MEM_FREE
Call HALTAF_MEM_FREE

! -----------------------------------------------------------------------
! Test 5: a system with several possible solid products
! -------
    WRITE(*,'("Test 5 (''noCalc'')")')
    WRITE(IOUT,20)
    WRITE(IOUT,'("Test 5 (''noCalc'')")')
    ! 6 chemical components are: H+, CO3-2, Ca+2, H2O, Na+, Cl-,
    ! 14 aqueous species are:     H+, CO3-2, Ca+2, H2O, Na+, Cl-,
    !    CaCO3, CaHCO3+, CaOH+, CO2, HCO3-, OH-, NaCO3-, NaHCO3
    ! 4 solid are: CaCO3(cr), CaCO3(s), Ca(OH)2(s), CaO(s)
    !              (both calcite and aragonite are included)
    NA = 6 !chemical components
    MSOL = 4 ! solids
    MS = 18  !species = components + reaction products
    SOLIDC = 0  ! how many components are solid phases
    ! allocate memory arrays
    Call CHEM_MEM_ALLOC
    Call HALTAF_MEM_ALLOC
    DO I=1,MS
      NOLL(I)=.FALSE.
    EndDo
    NOLL(4) = .TRUE. ! the mass balance of H2O is not included
! Data for components
    IDENTC(1) = "H+";  IDENTC(2) = "CO3-2";
    IDENTC(3) = "Ca+2";IDENTC(4) = "H2O";
    IDENTC(5) = "Na+"; IDENTC(6) = "Cl-"
    IDENT(1:NA) = IDENTC(1:NA)
    Z(1:NA)= (/ +1, -2, +2, 0, +1, -1 /) ! electric charges
! Data for reaction products
    A(1,1:6) = (/ 2.d0,1.d0,0.d0,-1.d0,0.d0,0.d0/); LBETA(1) = 16.681d0; IDENT(NA+1) = "CO2";     Z(NA+1) =  0
    A(2,1:6) = (/ 1.d0,1.d0,0.d0, 0.d0,0.d0,0.d0/); LBETA(2) = 10.329d0; IDENT(NA+2) = "HCO3-";   Z(NA+2) = -1
    A(3,1:6) = (/ 0.d0,1.d0,1.d0, 0.d0,0.d0,0.d0/); LBETA(3) =  3.224d0;  IDENT(NA+3) = "CaCO3";   Z(NA+3) =  0
    A(4,1:6) = (/ 1.d0,1.d0,1.d0, 0.d0,0.d0,0.d0/); LBETA(4) = 11.435d0; IDENT(NA+4) = "CaHCO3+"; Z(NA+4) = +1
    A(5,1:6) = (/-1.d0,0.d0,1.d0, 1.d0,0.d0,0.d0/); LBETA(5) =-12.78d0;   IDENT(NA+5) = "CaOH+";   Z(NA+5) = +1
    A(6,1:6) = (/-1.d0,0.d0,0.d0, 1.d0,0.d0,0.d0/); LBETA(6) =-14.00d0;  IDENT(NA+6) = "OH-";     Z(NA+6) = -1
    A(7,1:6) = (/ 0.d0,1.d0,0.d0, 0.d0,1.d0,0.d0/); LBETA(7) =  1.27d0;  IDENT(NA+7) = "NaCO3-"; Z(NA+7) = -1
    A(8,1:6) = (/ 1.d0,1.d0,0.d0, 0.d0,1.d0,0.d0/); LBETA(8) = 10.079d0;IDENT(NA+8) = "NaHCO3"; Z(NA+8) =  0
    A(9,1:6) = (/ 0.d0,1.d0,1.d0, 0.d0,0.d0,0.d0/); LBETA(9) =  8.48d0;  IDENT(NA+9) = "Calcite";
    A(10,1:6)= (/ 0.d0,1.d0,1.d0, 0.d0,0.d0,0.d0/); LBETA(10)=  8.336d0;IDENT(NA+10) = "Aragonite";
    A(11,1:6)= (/-2.d0,0.d0,1.d0, 2.d0,0.d0,0.d0/); LBETA(11)=-22.8d0;  IDENT(NA+11) = "Ca(OH)2(s)";
    A(12,1:6)= (/-2.d0,0.d0,1.d0, 1.d0,0.d0,0.d0/); LBETA(12)=-32.58d0; IDENT(NA+12) = "CaO(s)";
! Concentrations
    DO I = 1, NA
      KH(I)=1       ! only total concentrations given (KH=1)
    End Do
    KH(4) = 2 ! log activity given for H2O
    ! total concentrations for the components
    TOT(1)=  0.d0   ! H+
    TOT(2)=  0.001d0 ! CO3-2
    TOT(3)= -0.01  ! Ca+2 ! test what happens with impossible concentrations...
    TOT(4)=0.d0; LOGA(4)= 0.d0    ! H2O
    TOT(5)=  0.022d0 ! Na+
    TOT(6)=  0.02d0  ! Cl-
    TOL=1.D-6  ! tolerance when solving mass-balance equation
    JWATER = 4 ! H2O
    DBG=1; activityCoeffsModel = 0; ionicStr = -1.
    WRITE(IOUT,21) NA,SOLIDC,NIon,MSOL
    WRITE(IOUT,'("Find pH of an aqueous solution",/, &
    "with 0.02 mol NaCl and 0.001 Na2CO3")')
    WRITE(IOUT,22) (TOT(I),I=1,NA)
    Write(IOUT,'(" tol=",1pg10.3)') tol
    Call printChemSystem (IOUT)
    CONT=.false.
    CALL HaltaCalc
    WRITE(IOUT,23) errFlagsToString()
    pH=-LOGA(1)
    WRITE(IOUT,'(5X,"pH=",F6.3," (should be 10.418)")')  pH
    Call printConcs (IOUT)
    !FLUSH(IOUT)
Call CHEM_MEM_FREE
Call HALTAF_MEM_FREE

! -----------------------------------------------------------------------
   CLOSE(UNIT=IOUT)
   WRITE (*, '("All done.",/,"Press Enter to continue ...")', Advance='No')
   READ (*,'(I5)',iostat=ierr) I
END PROGRAM TEST_HaltaF
