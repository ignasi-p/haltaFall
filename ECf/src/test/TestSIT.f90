PROGRAM TestSIT
USE IO
USE SIT
USE CHEM
USE FACTOR_Module
Implicit NONE
Integer :: I
Character (len = 0) :: path = ""
!Integer, PARAMETER :: dp = kind(0.d0) ! double precision

    IOUT=IUT; ! IOUT is used in "CHEM", IUT is defined in "IO"
    OPEN(UNIT=IUT,FILE='Test_SIT.out')
    WRITE(*,'(a)') "output file: 'Test_SIT.out'"

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
    DO I = 1, NA
      KH(I)=1       ! only total concentrations given (KH=1)
    End Do
    KH(4) = 2 ! log(activity) given for H2O
    ! total concentrations for the components
    TOT(1)= 0.01 ! H+
    TOT(2)= 0.01 ! CO3-2
    TOT(3)= 0.01 ! Na+
    TOT(4)=0.D0; LOGA(4) = 0. ! H2O
    JWATER = 4 ! H2O
    activityCoeffsModel = 1;
    temperature = 25.d0; pressure = 1.d0
    ionicStr = -1.
    Call printChemSystem (IUT)

call SIT_MEM_ALLOC(nIon)

call readSITdata(path)

!call setEps(1,1,-0.1d0)

!call setEps(3,3,3.d0)

!call setEps(2,3,5.d0)

!call setEps(5,2,-5.d0)

call printEps (IUT)

close(UNIT=IUT)

call Quit

END PROGRAM TestSIT