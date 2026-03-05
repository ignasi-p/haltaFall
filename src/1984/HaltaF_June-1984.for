C     SUBROUTINE HALTA
C     ****************
C
C  Version: June-1984        With Activity Coefficients Correction.
C
C  From the HALTAFALL Program, Published in:
C  N.Ingri, W.Kakolowicz, L.G.Sillen, and B. Warnqvist, Talanta,14(1967)1261
C
C  This version with a new solid phase selection
C -----------------------
C  Subroutines derived by
C     Ignasi Puigdomenech
C     Dept. of Inorganic Chemistry
C     The Royal Institute Of Technology
C     100 44 Stockholm, Sweden,   Tel. 08/7878153
C -----------------------
C  The subroutine is dimensioned to  Na =  12 chemical components
C                                    Nx = 150 soluble complexes
C                                    Nf =  40 solid phases
C  Some Vectors are dimensioned to Na+Nx+Nf = 202
C                                     Na+Nx = 162
C                                     Nx+Nf = 190
C    If you change that, do not forget to change also the dimensions in
C    arrays outside common areas, the check in subroutine HALTA0,
C    and the lines you are now reading.
C For IBM-PC
$NOFLOATCALLS
$STORAGE:2
      SUBROUTINE HALTA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C---------------------
C This is the information that the user of HALTA either
C          - must supply before calling the subroutine, or
C          - may retrieve after HALTA has been called
      LOGICAL CONT,NOLL
      INTEGER HUR
      DOUBLE PRECISION LBETA,P
      COMMON /HLTF0/ CONT,NA,NX,NF,LBETA(190),P(12,190),HUR(12),
     + NOLL(202)
      DOUBLE PRECISION LOGA,TOT,TOL,C,LOGF
      COMMON /HLTF1/ IERR,LOGA(202),TOT(12),TOL(12),C(202),LOGF(162)
C---------------------
C These variables are only used by HALTA, and the user does not need to
C   worry about them.
      DOUBLE PRECISION LN10,LNBETA,LNKF,LNKMI,LNA,LNBA,LNG
      LOGICAL SINFAL,MONO,BER,FALLA,FALL
      COMMON /HLTF3/ LN10,INDIK,NVA,NFALL,NVAF,NUT,SINFAL,LNBETA(150),
     + LNKF(40),LNKMI(40),LNA(202),MONO(12),TOLY(12),LNG(162),FALL(40),
     + TOTVA(12),PVA(12,150),RUTA(40,40),RUT1(40,40),LNBA(150),IBE(12),
     + IFALL(12),IVAF(12),BER(12),FALLA(12),IVA(12),IVABRA(12),
     + IVANOV(12),ITERS,ITERC,IUT,NION,FSCAL(40)
      COMMON /HLTF2/ IVAR,X,Y,Y0,X1(12),X2(12),Y1(12),Y2(12),KARL(12)
     + ,STEG(12),STEP0,ITER(12)
C
C VARIABLES ----------------------------------------------------------------
C Ivar= nr of component for which Tot(Ivar) is being tested and Lna(Ivar)
C       adjusted (varied).
C---------------------------------------------------------------------------
C Some variables introduced for this version (related to activity coeff.):
C  Factor= a subroutine supplied by the user.
C  Logf(i)= log of activity coeff. of:
C  lng(i)= natural log of activity coeff. of:
C  Lna(i)=  natural log of activity of:
C  Loga(i)= log of activity of:
C  C(i)= concentrations of:
C                                    components          0 < i <= NA
C                                    soluble complexes  NA < i <= NA+NX
C                                    solids          NA+NX < i <= NA+NX+NF
C---------------------------------------------------------------------------
      IUT=0
      IF(IERR.GT.0) IUT=IERR
      IERR=0
      CALL HALTA0
      IF(IERR.LE.-4) GO TO 99999
      SINFAL=.FALSE.
      ITERS=0
      ITERC=0
C NYA1  Get Lna(), or Toly()
      DO 900 IA=1,NA
      LQ=HUR(IA)
      GO TO (901,902),LQ
901   TOLY(IA)=DABS(TOL(IA)*TOT(IA))
      IF(TOT(IA).EQ.0.D0)  TOLY(IA)=TOL(IA)*1.D-5
902   LNA(IA)=LN10*LOGA(IA)
900   CONTINUE
      IF(NVA.NE.0) GO TO 10000
C NVA=0   no mass balance equation needs to be solved
      IVAR=1
      CALL HALTA6(IVAR,X)
31    CALL HALTA9 (IVAR)
      GO TO (31,32),INDIK
32    CALL HALTA2
      GO TO 90000
C SLINGOR  Beguin with first Ivar
10000 DO 41 LRVA=1,NVA
      IA=IVA(LRVA)
      KARL(IA)=1
41    STEG(IA)=STEP0
      IVAR=IVA(1)
      IF(NFALL.GT.0)  CALL HALTA5
      CALL HALTA6(IVAR,X)
      IF(.NOT.BER(IVAR))  GO TO 11000
      CALL HALTA3(IVAR)
      IF(NVAF.GT.0) CALL HALTA8
      GO TO 12000
C TJAT  Solve the mass balance equation for component: Ivar
11000 LNA(IVAR)=X
      IF(.NOT.FALLA(IVAR))  GO TO 43
      CALL HALTA5
      CALL HALTA6(IVAR,X)
43    IF(IVAR.EQ.IVANOV(IVAR))  GO TO 44
      IVAR=IVANOV(IVAR)
      CALL HALTA6(IVAR,X)
      IF(.NOT.BER(IVAR)) GO TO 44
      CALL HALTA3(IVAR)
      IF(NVAF.GT.0) CALL HALTA8
      GO TO 12000
44    CALL HALTA3(IVAR)
      IF(NVAF.GT.0) CALL HALTA8
      CALL HALTA7
      GO TO (11000,12000,13000),INDIK
C PROV  The mass balance for component: Ivar was satisfied
12000 IVAR=IVABRA(IVAR)
C If it was the last Ivar, check the solid phases
      IF(IVAR.NE.0) GO TO 45
      CALL HALTA2
      IF(INDIK.EQ.1)  GO TO 9000
      IF(INDIK.EQ.2)  GO TO 10000
      IF(INDIK.EQ.3)  GO TO 9000
45    IF(BER(IVAR))  GO TO 12000
      CALL HALTA6(IVAR,X)
      CALL HALTA7
      GO TO (11000,12000,13000),INDIK
C ADA
13000 CALL HALTA1 (INDIK)
      GO TO 11000
C
C AKTIV  Calculate activity coefficients
C
9000  CALL HALTA9
      GO TO (10000,90000),INDIK
C
C NOG   This is the end.  Calculate Loga(i),Logf(i) and Tot(ia)
C
90000 DO 87 LIA=1,NA
      LOGF(LIA)=LNG(LIA)/LN10
      IF(HUR(LIA).EQ.1)  GO TO 89
      TOT(LIA)=C(LIA)
      DO 86 LIX=1,NX
      LIAX=NA+LIX
86    TOT(LIA)=TOT(LIA)+P(LIA,LIX)*C(LIAX)
      IF(NF.EQ.0) GO TO 87
      DO 88 LIX=1,NF
      LIAX=NION+LIX
      LIAF=NX+LIX
88    TOT(LIA)=TOT(LIA)+P(LIA,LIAF)*C(LIAX)
      GO TO 87
89    LOGA(LIA)=LNA(LIA)/LN10
87    CONTINUE
C
      DO 85 LIX=1,NX
      LIAX=NA+LIX
      LOGF(LIAX)=LNG(LIAX)/LN10
85    LOGA(LIAX)=LNA(LIAX)/LN10
C
      IF(NF.EQ.0) GO TO 83
      DO 84 LIX=1,NF
      LIAX=NION+LIX
84    LOGA(LIAX)=LNA(LIAX)/LN10
83    CONTINUE
C
99999 CONT=.TRUE.
      RETURN
      END
C0 -----------------------------------------------------------------------
      SUBROUTINE HALTA0
C The first time HALTA is called, a plan is made to solve the mass
C balance equations, some variables are initialized, etc.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONT,NOLL
      INTEGER HUR
      DOUBLE PRECISION LBETA,P
      COMMON /HLTF0/ CONT,NA,NX,NF,LBETA(190),P(12,190),HUR(12),
     + NOLL(202)
      DOUBLE PRECISION LOGA,TOT,TOL,C,LOGF
      COMMON /HLTF1/ IERR,LOGA(202),TOT(12),TOL(12),C(202),LOGF(162)
      DOUBLE PRECISION LN10,LNBETA,LNKF,LNKMI,LNA,LNBA,LNG
      LOGICAL SINFAL,MONO,BER,FALLA,FALL
      COMMON /HLTF3/ LN10,INDIK,NVA,NFALL,NVAF,NUT,SINFAL,LNBETA(150),
     + LNKF(40),LNKMI(40),LNA(202),MONO(12),TOLY(12),LNG(162),FALL(40),
     + TOTVA(12),PVA(12,150),RUTA(40,40),RUT1(40,40),LNBA(150),IBE(12),
     + IFALL(12),IVAF(12),BER(12),FALLA(12),IVA(12),IVABRA(12),
     + IVANOV(12),ITERS,ITERC,IUT,NION,FSCAL(40)
      COMMON /HLTF2/ IVAR,X,Y,Y0,X1(12),X2(12),Y1(12),Y2(12),KARL(12)
     + ,STEG(12),STEP0,ITER(12)
      INTEGER NOBER(12)
      LOGICAL OBER(12,12),NOCALC(12),POS(12),NEG(12),OK
C VARIABLES -----------------------------------------------------------------
C  Nva= nr of equations for Tot(ia) to be tested (Lna(Ivar) to be varied)
C       in absence of solids
C  Mono(ia)=.True. if component 'ia' only takes part in mononuclear
C           complexes (P(ia,ix)=0 or 1 for all ix)
C  Noll(ia)=.True. if the concentration for species 'ia' is not included
C           in Tot(ia)
C  Iva(m)= ia number for m'th Tot(ia) to be tested in absence of solids
C          (m=1 to Nva)
C  Ivanov(Ivar)= ia to be tested after Ivar, if the mass balance for Tot(Ivar)
C                is not satisfied
C  Ivabra(Ivar)= Ivar to be tested after Ivar if the mass balance for
C                Tot(Ivar) is satisfied.
C----------------------------------------------------------------------------
      IF(CONT) GO TO 2000
      STEP0=0.1D0
      LN10=DLOG(10.D0)
C     CALL ERRSET(0)
C NYKO
      IF(NA.LE.0.OR.NA.GT.12) GO TO 97
      IF(NX.LE.0.OR.NX.GT.150) GO TO 97
      IF(NF.LT.0.OR.NF.GT.40) GO TO 97
      GO TO 99
97    IF(IUT.GT.0)  WRITE(IUT,98)
98    FORMAT(' ?? Maximum values:',/,11X,'Components: 0 < NA <= 12'
     + ,/,11X,'Soluble complexes: 1 <= NX <= 150',/,11X
     + ,'Solids: 0 <= NF <= 40')
      IERR=-4
      RETURN
99    DO 2 IA=1,NA
      IF(HUR(IA).EQ.1.OR.HUR(IA).EQ.2) GO TO 2
      IERR=-5
      IF(IUT.GT.0)  WRITE(IUT,95) IA
95    FORMAT(' ?? HUR outside allowed values for component:',I5)
      RETURN
2     CONTINUE
      NION=NA+NX
      DO 11 LIX=1,NX
      LNBETA(LIX)=LN10*LBETA(LIX)
11    CONTINUE
      DO 12 LIA=1,NA
      ITER(LIA)=0
      MONO(LIA)=.TRUE.
      POS(LIA)=.FALSE.
      IF(.NOT.NOLL(LIA)) POS(LIA)=.TRUE.
      NEG(LIA)=.FALSE.
      DO 121 LIX=1,NX
      LIAX=LIX+NA
      IF(NOLL(LIAX)) GO TO 121
      IF(DABS(P(LIA,LIX)).GT.0.001.AND.DABS(P(LIA,LIX)-1.D0).GT.0.001)
     +     MONO(LIA)=.FALSE.
      IF(P(LIA,LIX).GT.0.D0) POS(LIA)=.TRUE.
      IF(P(LIA,LIX).LT.0.D0) NEG(LIA)=.TRUE.
121   CONTINUE
12    CONTINUE
      IF(NF.EQ.0)  GO TO 141
      DO 14 LIF=1,NF
      LIAX=NA+NX+LIF
      LIAF=NX+LIF
      LNKF(LIF)=-LBETA(LIAF)*LN10
      FSCAL(LIF)=1.D0
      IF(NOLL(LIAX)) GO TO 14
C  The solids get a scaling factor FSCAL
C  to determine which one of the oversaturated solids will be allowed
C  to precipitate
          DO 142 LIA=1,NA
          FSCAL(LIF)=FSCAL(LIF) + DABS(P(LIA,LIAF))
142       CONTINUE
14    CONTINUE
141   DO 18 LI=1,NA
      DO 18 LJ=1,NA
      IF(LI.NE.LJ)  GO TO 15
      OBER(LI,LJ)=.FALSE.
      GO TO 16
15    OBER(LI,LJ)=.TRUE.
16    DO 17 LIX=1,NX
      IF(NOLL(NA+LIX)) GO TO 17
      M=P(LI,LIX)*P(LJ,LIX)
      IF(M.NE.0)  OBER(LI,LJ)=.FALSE.
17    CONTINUE
      IF(NF.EQ.0)  GO TO 18
      DO 181 LIF=1,NF
      IF(NOLL(NION+LIF)) GO TO 181
      LIAF=NX+LIF
      M=P(LI,LIAF)*P(LJ,LIAF)
      IF(M.NE.0)  OBER(LI,LJ)=.FALSE.
181   CONTINUE
18    CONTINUE
20    DO 19 LI=1,NA
      NOBER(LI)=0
      DO 19 LJ=1,NA
      IF(OBER(LI,LJ))  NOBER(LI)=NOBER(LI)+1
19    CONTINUE
      DO 21 IA=1,NA
      IF(HUR(IA).EQ.2) GO TO 21
      NOCALC(IA)=.FALSE.
      IF(POS(IA).AND.NEG(IA)) GO TO 21
      IF(POS(IA).AND.TOT(IA).GT.0.D0) GO TO 21
      IF(NEG(IA).AND.TOT(IA).LT.0.D0) GO TO 21
      LOGA(IA)=-99999.999D0
      NOCALC(IA)=.TRUE.
21    CONTINUE
C NYA
999   NFALL=0
      DO 992 IA=1,NA
      BER(IA)=.FALSE.
992   FALLA(IA)=.FALSE.
      IF(NF.EQ.0) GO TO 23
      DO 22 IF=1,NF
      FALL(IF)=.FALSE.
22    CONTINUE
23    NVA=0
      IA=0
3000  IA=IA+1
      IF(IA.GT.NA)  GO TO 4000
      IF(HUR(IA).EQ.2.OR.NOCALC(IA))  GO TO 3000
      NVA=NVA+1
      IVA(NVA)=IA
      GO TO 3000
C PLAN
4000  M=0
      IF(NVA.LE.1)  GO TO 261
      LQ=NVA-1
      DO 26 LI=1,LQ
      IVA1=IVA(LI)
      IVA2=IVA(LI+1)
      IF(NOBER(IVA2).LE.NOBER(IVA1))  GO TO 25
29    M=IVA(LI)
      IVA(LI)=IVA(LI+1)
      IVA(LI+1)=M
      GO TO 26
25    IF(MONO(IVA2).AND..NOT.MONO(IVA1).AND.NOBER(IVA2)
     + .EQ.NOBER(IVA1))        GO TO 29
26    CONTINUE
      IF(M.NE.0)  GO TO 4000
261   IVA(NVA+1)=0
      IVABRA(NVA+1)=0
      IF(NVA.EQ.0) RETURN
      DO 27 LI=1,NVA
      IVAR=IVA(LI)
      IVABRA(IVAR)=IVA(LI+1)
      I=0
28    I=I+1
      IVA3=IVA(I)
      IF(OBER(IVAR,IVA3))  GO TO 28
27    IVANOV(IVAR)=IVA(I)
      RETURN
Check if for some component the Tot.Conc. has changed in a way that makes
C the mass-balance equations to have no solution.
2000  OK=.TRUE.
      DO 2001 IA=1,NA
      IF(HUR(IA).EQ.2) GO TO 2001
      IF(POS(IA).AND.NEG(IA)) GO TO 2001
      IF(POS(IA).AND.TOT(IA).GT.0.D0) GO TO 2003
      IF(NEG(IA).AND.TOT(IA).LT.0.D0) GO TO 2003
      LOGA(IA)=-99999.999D0
      IF(.NOT.NOCALC(IA)) OK=.FALSE.
      NOCALC(IA)=.TRUE.
      GO TO 2001
2003  IF(.NOT.NOCALC(IA)) GO TO 2001
      OK=.FALSE.
      NOCALC(IA)=.FALSE.
      LOGA(IA)=-10.D0
      IF(TOT(IA).GT.0.D0) LOGA(IA)=DLOG10(TOT(IA))-2.D0
2001  CONTINUE
      IF(.NOT.OK) GO TO 999
      RETURN
      END
C1 -----------------------------------------------------------------------
      SUBROUTINE HALTA1 (INDIK)
C  Solving the mass balance equation of component: Ivar, which forms
C    polynuclear complexes
C VARIABLES----------------------------------------------------------------
C  X= independent variable in the equation Y0=Y(X)
C  Y= dependent variable
C  Y0= value aimed at in Y0=Y(X)
C  X1(ia) and X2(ia)= values for X below and above the right value
C  Y1(ia) and Y2(ia)= values for Y corresponding to X1 and X2
C--------------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER INDIK
      COMMON /HLTF2/ IVAR,X,Y,Y0,X1(12),X2(12),Y1(12),Y2(12),KARL(12)
     + ,STEG(12),STEP0,ITER(12)
      LOGICAL YL
      DOUBLE PRECISION W,W1,STEGBY
      DATA STEGBY/0.01D0/
C  Locating the solution by giving Steg(Ivar)=0.5,1,2,4,8,16,32,64,...
C    until a pair of X1 and X2 is found (in Karl(Ivar)=2 or 3). Then
C    decreasing Steg(Ivar) (=Steg(Ivar)*0.5, in Karl(Ivar)=4) until
C    X1 and X2 are separated less than Stegby
      YL=Y.GT.Y0
      IF(YL)  GO TO 45
      X1(IVAR)=X
      Y1(IVAR)=Y
      GO TO 46
45    X2(IVAR)=X
      Y2(IVAR)=Y
46    LQ=KARL(IVAR)
      GO TO (1301,1302,1303,1304),LQ
C Karl=1  the beguining
1301  IF(YL)  GO TO 47
      KARL(IVAR)=2
      X=X+STEG(IVAR)
      RETURN
47    KARL(IVAR)=3
      X=X-STEG(IVAR)
      RETURN
C Karl=2   it was Y<Y0
1302  IF(YL) GO TO 1350
      STEG(IVAR)=STEG(IVAR)+STEG(IVAR)
      X=X+STEG(IVAR)
      RETURN
C Karl=3   it was Y>Y0
1303  IF(.NOT.YL) GO TO 1350
      STEG(IVAR)=STEG(IVAR)+STEG(IVAR)
      X=X-STEG(IVAR)
      RETURN
1350  KARL(IVAR)=4
C Karl=4   There is both X1 and X2 corresponding to Y<Y0 and Y>Y0
1304  IF(STEG(IVAR).LT.STEGBY) GO TO 1370
      STEG(IVAR)=STEG(IVAR)*0.5D0
      IF(YL)  GO TO 49
      X=X+STEG(IVAR)
      RETURN
49    X=X-STEG(IVAR)
      RETURN
C
C Aproximating the solution by chord shooting ('secant method')
C
1370  W=Y0-Y1(IVAR)
      W1=X2(IVAR)-X1(IVAR)
      IF(W .LT.1.D-10) GO TO 1371
      IF(W1.LT.1.D-14) INDIK=10
      X=X1(IVAR)+ W*W1/(Y2(IVAR)-Y1(IVAR))
      RETURN
C Some rounding errors might be disturbing.
1371  IF(YL) X=X - 5.D-6
      IF(.NOT.YL) X=X + 5.D-6
      RETURN
      END
C2 -----------------------------------------------------------------------
C               FASTA
      SUBROUTINE HALTA2
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONT,NOLL
      INTEGER HUR
      DOUBLE PRECISION LBETA,P
      COMMON /HLTF0/ CONT,NA,NX,NF,LBETA(190),P(12,190),HUR(12),
     + NOLL(202)
      DOUBLE PRECISION LOGA,TOT,TOL,C,LOGF
      COMMON /HLTF1/ IERR,LOGA(202),TOT(12),TOL(12),C(202),LOGF(162)
      DOUBLE PRECISION LN10,LNBETA,LNKF,LNKMI,LNA,LNBA,LNG
      LOGICAL SINFAL,MONO,BER,FALLA,FALL
      COMMON /HLTF3/ LN10,INDIK,NVA,NFALL,NVAF,NUT,SINFAL,LNBETA(150),
     + LNKF(40),LNKMI(40),LNA(202),MONO(12),TOLY(12),LNG(162),FALL(40),
     + TOTVA(12),PVA(12,150),RUTA(40,40),RUT1(40,40),LNBA(150),IBE(12),
     + IFALL(12),IVAF(12),BER(12),FALLA(12),IVA(12),IVABRA(12),
     + IVANOV(12),ITERS,ITERC,IUT,NION,FSCAL(40)
      DOUBLE PRECISION TOTBE(12),PIVOT(40)
      INTEGER NFSPAR,NYFALL,FUT(12),IBER(12),IFSPAR(12),IPIVOT(40),
     + INDEX(40,2)
      LOGICAL BRA
      DATA ITERSM/20/
C VARIABLES -------------------------------------------------------------------
C  Ber(ia)=.True. if LNA(ia) is calculated by means of the solubility product
C          in subroutine HALTA5 (Lnaber)
C  Fall(if)=.True. if the solid 'if' is assumed present at equilibrium
C  Falla(ia)=.True. if the component 'ia' is assumed to occur in one or more
C            solids at equilibrium
C  Fut(i)= Ifspar number of i'th solid to be eliminated at some stage in
C          systematic variation during Sinfal=.True.
C  Ibe(m)= 'ia' number for m'th LNA(ia) to be calculated with the solubility
C         product in HALTA5 (Lnaber) (m=1 to Nfall)
C  Iber(j)= j'th of 'Iva' numbers to be an 'Ibe' (Ibe(j)=Iva(Iber(j)), (j=1 to
C           Nfall)
C  Ifall(m)= 'if' number of m'th solid present at equilibrium (m=1 to Nfall)
C  Ifspar(m)= if number of m'th solid indicated as possible a FALLPROV
C  Iva(m)= 'ia' number for m'th component to be tested in absence of solids
C         (m=1 to Nva)
C  Ivaf(m)='ia' number for m'th component to be tested in presence of solids
C          (m=1 to Nvaf)
C  Nfall= nr of solids present at equilibrium
C  Nfspar= nr of solids indicated at FALLPROV
C  Nvaf= nr of mass balance equations to be solved in presence of solids
C  Nut= nr of solids systematically eliminated at some stage whyle Sinfall is
C       .True.
C  Sinfal= .True. if the first matrix Ruta tried has come out singular
C -----------------------------------------------------------------------------
C After a set of Lna(ia) has been found (in routine HALTA) that satisfy all
C mass balance equations, two tests are made at FALLPROV
C
C FALLPROV  Check that the solubility product is not exeeded for any solid
C            phase assumed to be absent. Find Ifall() and Fall() for the new
C            solids appearing
      IF(NF.EQ.0)  GO TO 25000
      BRA=.TRUE.
      NYFALL=0
      IMAX=0
      AMAX=-1.D30
      DO 51 LIF=1,NF
      LIAX=NION+LIF
      IF(FALL(LIF)) GO TO 63
      LIAF=NX+LIF
      C(LIAX)=0.D0
      W=0.D0
      DO 50 LIA=1,NA
50    W=W+LNA(LIA)*P(LIA,LIAF)
      LNA(LIAX)=W-LNKF(LIF)
      IF(W.LE.LNKF(LIF).OR.NOLL(LIAX).OR.NVA.EQ.0) GO TO 51
      BRA=.FALSE.
      W=LNA(LIAX)/FSCAL(LIF)
      IF(W.LE.AMAX) GO TO 51
      AMAX=W
      IMAX=LIF
      GO TO 51
63    LNA(LIAX)=0.D0
51    CONTINUE
      IF(IMAX.EQ.0) GO TO 49
      FALL(IMAX)=.TRUE.
      NYFALL=1
      IFALL(NFALL+NYFALL)=IMAX
49    IF(NVA.EQ.0) RETURN
C  Check that the quantity of solid is not negative for any solid
C  assumed to be present.
      IF(NFALL.EQ.0)  GO TO 56
      DO 53 LI=1,NFALL
      IA=IBE(LI)
      W=TOT(IA)
      W=W-C(IA)
      DO 52 LIX=1,NX
      LIAX=NA+LIX
52    W=W-P(IA,LIX)*C(LIAX)
53    TOTBE(LI)=W
      DO 55 LI=1,NFALL
      W=0.D0
      DO 54 LJ=1,NFALL
54    W=W+RUTA(LI,LJ)*TOTBE(LJ)
      LQ=IFALL(LI)
      LQA=NION+LQ
      C(LQA)=W
      IF(C(LQA).GE.0.D0)  GO TO 55
      FALL(LQ)=.FALSE.
      C(LQA)=0.D0
      BRA=.FALSE.
55    CONTINUE
56    IF(BRA)  GO TO 25000
      IF(ITERS.LE.ITERSM) GO TO 15000
      IERR = -3
      IF(IUT.GT.0)  WRITE(IUT,561)
561   FORMAT(' ?? Error:  20 different solid phase ',
     + 'combinations',/,'        were proved and found not ',
     + 'satisfactory')
      GO TO 25000
C-------------------------------------------------------------------------
C The Lna(ia) were not consistent with the solid phases. Either Nyfall new
C solid phases appeared, or some solids assumed present had negative C(if).
C At first it is assumed at BEFALL that the first Nfall of the Iva(ia) are
C the Ibe() to be calculated. If the determinant of Ruta is found to be zero
C at UTFALL, the Ibe are changed systematically at SING-HOPPSI by means of
C array Iber, until a non-zero determinant is found.
C Then at ANFALL, the arrays Rut1 and Pva are calculated, and the new mass
C balance equations are solved in SLINGOR (routine HALTA).
C
C INFALL  find new Nfall, Ifall(Nfall)
15000 NFALL=NFALL+NYFALL
      LI=0
57    LI=LI+1
      IF=IFALL(LI)
      IF(FALL(IF))  GO TO 62
      NFALL=NFALL-1
      IF(NFALL.EQ.0)  GO TO 59
      DO 58 LJ=LI,NFALL
58    IFALL(LJ)=IFALL(LJ+1)
59    LI=LI-1
62    IF(LI.LT.NFALL)  GO TO 57
C
      IF(NFALL.EQ.0)  GO TO 24000
Check if Sinfal=.False., and if Nfall <= Nva
      IF(SINFAL) GO TO 65
        DO 64 LI=1,NFALL
64      IFSPAR(LI)=IFALL(LI)
        NFSPAR=NFALL
        NUT=0
      IF(NFALL.LE.NVA) GO TO 16000
      SINFAL=.TRUE.
      NUT=NFALL-NVA-1
      GO TO 23000
65    IF(NFALL.LT.NVA) GO TO 16000
      I=0
      NFALL=NFSPAR-NUT
      GO TO 21000
C BEFALL  first guess if Iber(j) set
16000 DO 66 LI=1,NFALL
66    IBER(LI)=LI
      IBER(NFALL+1)=0
C UTFALL  calculate RUTA(I,J) and check that it is not singular
17000 DO 67 LI=1,NFALL
      IA=IVA(IBER(LI))
      DO 67 LJ=1,NFALL
      IF=NX+IFALL(LJ)
67    RUTA(LI,LJ)=P(IA,IF)
      INDIK=0
      CALL HALTA4 (NFALL,RUTA,INDIK,PIVOT,IPIVOT,INDEX)
      IF(INDIK.EQ.0)  GO TO 24000
C SING   the matrix RUTA was singular
18000 I=0
      SINFAL=.TRUE.
C HOPPSI  get a new Iber(j) set
19000 I=I+1
      IF(IBER(I).EQ.NVA)  GO TO 20000
      IF(IBER(I+1).EQ.IBER(I)+1)  GO TO 19000
      IBER(I)=IBER(I)+1
      IF(I.EQ.1)  GO TO 17000
      DO 68 LJ=1,I-1
68    IBER(LJ)=LJ
      GO TO 17000
C--------------------------------------------------------------------------
C Either
C  - It has been proved inpossible to pick out a group of Nfall components
C    Ibe() such that the array Ruta has a non zero determinant,
C  - Or more solids are indicated than there are LOGA() values to vary
C    (Nfall > Nva).
C One must try systematically combinations of a smaller number (Nfall-Nut)
C of solid phases, until a non-zero determinant for Ruta is found. This is done
C at labels FUT, INFUT and UPPNUT using the array Fut(). To avoid going into a
C loop, the program sets Sinfal=.True. and remembers the first set of solids
C (Ifspar(Nfspar)) which gave only zero determinants.
C
C FUTT    Iber sets exausted
20000 IF(NUT.EQ.0)  GO TO 23000
      I=0
C HOPP FUT    Get a new Fut() set
21000 I=I+1
      IF(FUT(I).EQ.NFSPAR)  GO TO 23000
      IF(FUT(I+1).EQ.FUT(I)+1)  GO TO 21000
      FUT(I)=FUT(I)+1
      IF(I.EQ.1)  GO TO 22000
      DO 69 LJ=1,I-1
69    FUT(LJ)=LJ
C INFUT  find new Fall(), Ifall() for the reduced group of solids
22000 DO 70 LI=1,NFSPAR
70    FALL(IFSPAR(LI))=.TRUE.
      IF(NUT.EQ.0)  GO TO 72
      DO 71 LI=1,NUT
71    FALL(IFSPAR(FUT(LI)))=.FALSE.
72    J=0
      DO 73 LI=1,NFSPAR
      IF=IFSPAR(LI)
      IF(.NOT.FALL(IF))  GO TO 73
      J=J+1
      IFALL(J)=IF
73    CONTINUE
      GO TO 16000
C UPPNUT  another solid must be taken away. Get a first guess of Fut()
23000 NUT=NUT+1
      NFALL=NFSPAR-NUT
      IF(NFALL.EQ.0)  GO TO 75
      DO 74 LI=1,NUT
74    FUT(LI)=LI
      FUT(NUT+1)=0
      GO TO 22000
75    CONTINUE
C NYP(IN)
      IF(IUT.GT.0) WRITE(IUT,752) NFSPAR,(IFSPAR(LI),LI=1,NFSPAR)
752   FORMAT(' ?? Error:',/,'  The'I3' solid phases that appeared',
     +' give only zero determinants.',/,'  Solids:',15I4,/,9X,15I4)
      DO 78 IF=1,NF
78    FALL(IF)=.FALSE.
      DO 85 IA=1,NA
      BER(IA)=.FALSE.
85    FALLA(IA)=.FALSE.
      IERR=-2
      INDIK=1
      RETURN
C------------------------------------------------------------------------
C A set of Nfall solids has been found with a non-zero determinant for the
C array Ruta. Some variables are changed, and the arrays Rut1 and Pva are
C calculated.
C
C ANFALL
24000 DO 76 LIA=1,NA
      FALLA(LIA)=.FALSE.
76    BER(LIA)=.FALSE.
C  Get Falla(ia), Ber(ia) and Ibe(ia) (from Ifall and Fall)
      IF(NFALL.EQ.0)  GO TO 10000
      DO 77 LI=1,NFALL
      IF=IFALL(LI)
      LIAF=NX+IF
      IF(.NOT.FALL(IF)) GO TO 60
      DO 61 LIA=1,NA
61    IF(P(LIA,LIAF).NE.0.D0)  FALLA(LIA)=.TRUE.
60    CONTINUE
      IBE(LI)=IVA(IBER(LI))
      IA=IBE(LI)
77    BER(IA)=.TRUE.
C  Find Nvaf and Ivaf()
      NVAF=0
      IF(NVA.EQ.0)  GO TO 10000
      DO 79 LI=1,NVA
      IA=IVA(LI)
      IF(BER(IA))  GO TO 79
      NVAF=NVAF+1
      IVAF(NVAF)=IA
79    CONTINUE
      IF(NVAF.EQ.0)  GO TO 10000
Calculate RUT1(I,J)
      DO 81 LI=1,NVAF
      DO 81 LJ=1,NFALL
      W=0.D0
      DO 80 LM=1,NFALL
      LQ=IVAF(LI)
      LZ=NX+IFALL(LM)
80    W=W+P(LQ,LZ)*RUTA(LM,LJ)
81    RUT1(LI,LJ)=W
Calculate PVA(IA,IX)
82    DO 84 LI=1,NVAF
      DO 84 LIX=1,NX
      IA=IVAF(LI)
      W=P(IA,LIX)
      IF(NFALL.EQ.0)  GO TO 84
      DO 83 LM=1,NFALL
      LQ=IBE(LM)
83    W=W-RUT1(LI,LM)*P(LQ,LIX)
84    PVA(IA,LIX)=W
C SLINGOR(IN)
10000 INDIK=2
      RETURN
C NOG(IN)
25000 INDIK=3
      RETURN
      END
C3 -----------------------------------------------------------------------
C            CBER
      SUBROUTINE HALTA3(IVAR)
C  1.-  Calculation of the activity of soluble complexes.
C  2.-  Calculation of the concentrations of all soluble species
C       with old values for activity coefficients.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONT,NOLL
      INTEGER HUR
      DOUBLE PRECISION LBETA,P
      COMMON /HLTF0/ CONT,NA,NX,NF,LBETA(190),P(12,190),HUR(12),
     + NOLL(202)
      DOUBLE PRECISION LOGA,TOT,TOL,C,LOGF
      COMMON /HLTF1/ IERR,LOGA(202),TOT(12),TOL(12),C(202),LOGF(162)
      DOUBLE PRECISION LN10,LNBETA,LNKF,LNKMI,LNA,LNBA,LNG
      LOGICAL SINFAL,MONO,BER,FALLA,FALL
      COMMON /HLTF3/ LN10,INDIK,NVA,NFALL,NVAF,NUT,SINFAL,LNBETA(150),
     + LNKF(40),LNKMI(40),LNA(202),MONO(12),TOLY(12),LNG(162),FALL(40),
     + TOTVA(12),PVA(12,150),RUTA(40,40),RUT1(40,40),LNBA(150),IBE(12),
     + IFALL(12),IVAF(12),BER(12),FALLA(12),IVA(12),IVABRA(12),
     + IVANOV(12),ITERS,ITERC,IUT,NION,FSCAL(40)
      DOUBLE PRECISION LNC
      INTEGER IVAR
Calculate activities of soluble complexes
      DO 1 LIX=1,NX
      LIAX=NA+LIX
1     LNA(LIAX)=LNBA(LIX)+P(IVAR,LIX)*LNA(IVAR)
Calculate Concentrations:
C   components and soluble complexes
100   DO 300 LIA=1,NION
      C(LIA)=0.D0
      IF(NOLL(LIA)) GO TO 300
      LNC=LNA(LIA)-LNG(LIA)
      IF(LNC.GT.80.D0) LNC=80.D0
      IF(LNC.GT.-80.D0) C(LIA)=DEXP(LNC)
300   CONTINUE
      RETURN
      END
C9 -----------------------------------------------------------------------
      SUBROUTINE HALTA9
C  Calculation of activity coefficients
C  If the activity coeficients (which depend on the composition of the
C  fluid) have varied since they were last calculated, then recalculate
C  the equilibrium composition.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONT,NOLL
      INTEGER HUR
      DOUBLE PRECISION LBETA,P
      COMMON /HLTF0/ CONT,NA,NX,NF,LBETA(190),P(12,190),HUR(12),
     + NOLL(202)
      DOUBLE PRECISION LOGA,TOT,TOL,C,LOGF
      COMMON /HLTF1/ IERR,LOGA(202),TOT(12),TOL(12),C(202),LOGF(162)
      DOUBLE PRECISION LN10,LNBETA,LNKF,LNKMI,LNA,LNBA,LNG
      LOGICAL SINFAL,MONO,BER,FALLA,FALL
      COMMON /HLTF3/ LN10,INDIK,NVA,NFALL,NVAF,NUT,SINFAL,LNBETA(150),
     + LNKF(40),LNKMI(40),LNA(202),MONO(12),TOLY(12),LNG(162),FALL(40),
     + TOTVA(12),PVA(12,150),RUTA(40,40),RUT1(40,40),LNBA(150),IBE(12),
     + IFALL(12),IVAF(12),BER(12),FALLA(12),IVA(12),IVABRA(12),
     + IVANOV(12),ITERS,ITERC,IUT,NION,FSCAL(40)
      DOUBLE PRECISION F,TEMP,ATEMP,FL,CL,TOLLNG,OLDLNG(164)
      DATA FL/0.5D0/,CL/0.25D0/,TOLLNG/6.0D-4/
C  This version is specially adapted to avoid some oscillatoty behaviour
C  found in some aqueous solutions, like AlCl3, etc, where the activity
C  coefficients and the aqueous composition are affecting each other in
C  a strange way. The subroutine does give slightly slower convergence for
C  solutions of NaCl for example, but as it is now it seems more reliable
C  to converge.
C
      CALL FACTOR (NION,C,LNG)
C
      INDIK=2
      DO 4 LIA=1,NION
      TEMP=LNG(LIA)-OLDLNG(LIA)
      ATEMP=DABS(TEMP)
      IF(ATEMP.GT.TOLLNG) INDIK=1
      IF(ITERC.EQ.0) GO TO 3
      IF(C(LIA).LE.CL) GO TO 3
      F=FL/(DSQRT(C(LIA)))
      TEMP= TEMP*F
      LNG(LIA)=OLDLNG(LIA)+TEMP
3     OLDLNG(LIA)=LNG(LIA)
4     CONTINUE
      ITERC=ITERC+1
      IF(ITERC.LE.50) RETURN
      INDIK=2
      IERR=-5
      RETURN
      END
C4 ----------------------------------------------------------------------
C          INVERT
      SUBROUTINE HALTA4 (N,A,INDIC,PIVOT,IPIVOT,INDEX)
C Matrix inversion with accompanying solution of linear equations
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(40,N),PIVOT(N),SWAP,AMAX,T
      INTEGER IPIVOT(N),INDEX(N,2),INDIC
C
C     Initialization
C
      DO 20 J=1,N
20    IPIVOT(J)=0
      DO 550 I=1,N
C
C     Search for pivot element
C
      AMAX=0.D0
      DO 105 J=1,N
      IF(IPIVOT(J).EQ.1)GO TO 105
      DO 100 K=1,N
      IF(IPIVOT(K)-1)80,100,320
80    IF(DABS(AMAX)-DABS(A(J,K)))85,85,100
85    IROW=J
      ICOLUM=K
      AMAX=A(J,K)
100   CONTINUE
105   CONTINUE
      IPIVOT(ICOLUM)=IPIVOT(ICOLUM)+1
C
C     Interchange rows to put pivot element on diagonal
C
      IF(IROW.EQ.ICOLUM)GO TO 260
      DO 200 L=1,N
      SWAP=A(IROW,L)
      A(IROW,L)=A(ICOLUM,L)
200   A(ICOLUM,L)=SWAP
260   INDEX(I,1)=IROW
      INDEX(I,2)=ICOLUM
      PIVOT(I)=A(ICOLUM,ICOLUM)
C
C     Divide pivot row by pivot element
C
      IF(DABS(AMAX).GT.0.D0)GO TO 330
320   INDIC=1
      RETURN
330   A(ICOLUM,ICOLUM)=1.0D0
      DO 350 L=1,N
350   A(ICOLUM,L)=A(ICOLUM,L)/PIVOT(I)
C
C     Reduce non-pivot rows
C
      DO 550 L1=1,N
      IF(L1.EQ.ICOLUM)GO TO 550
      T=A(L1,ICOLUM)
      A(L1,ICOLUM)=0.D0
      DO 450 L=1,N
450   A(L1,L)=A(L1,L)-A(ICOLUM,L)*T
550   CONTINUE
C
C     Interchange columns
C
      DO 710 I=1,N
      L=N+1-I
      IF(INDEX(L,1).EQ.INDEX(L,2))GO TO 710
      JROW=INDEX(L,1)
      JCOLUM=INDEX(L,2)
      DO 705 K=1,N
      SWAP=A(K,JROW)
      A(K,JROW)=A(K,JCOLUM)
      A(K,JCOLUM)=SWAP
705   CONTINUE
710   CONTINUE
      RETURN
      END
C5 --------------------------------------------------------------------------
C               LNABER - Part 1
      SUBROUTINE HALTA5
Calculates Lnkmi and Lna(Ibe()) when some solid phase is assumed
C   to be present.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONT,NOLL
      INTEGER HUR
      DOUBLE PRECISION LBETA,P
      COMMON /HLTF0/ CONT,NA,NX,NF,LBETA(190),P(12,190),HUR(12),
     + NOLL(202)
      DOUBLE PRECISION LOGA,TOT,TOL,C,LOGF
      COMMON /HLTF1/ IERR,LOGA(202),TOT(12),TOL(12),C(202),LOGF(162)
      DOUBLE PRECISION LN10,LNBETA,LNKF,LNKMI,LNA,LNBA,LNG
      LOGICAL SINFAL,MONO,BER,FALLA,FALL
      COMMON /HLTF3/ LN10,INDIK,NVA,NFALL,NVAF,NUT,SINFAL,LNBETA(150),
     + LNKF(40),LNKMI(40),LNA(202),MONO(12),TOLY(12),LNG(162),FALL(40),
     + TOTVA(12),PVA(12,150),RUTA(40,40),RUT1(40,40),LNBA(150),IBE(12),
     + IFALL(12),IVAF(12),BER(12),FALLA(12),IVA(12),IVABRA(12),
     + IVANOV(12),ITERS,ITERC,IUT,NION,FSCAL(40)
Calculate LNKMI
      DO 2 LI=1,NFALL
      IF=IFALL(LI)
      LIAF=NX+IF
      W=LNKF(IF)
      DO 1 LIA=1,NA
      IF(.NOT.BER(LIA))  W=W-P(LIA,LIAF)*LNA(LIA)
1     CONTINUE
2     LNKMI(LI)=W
Calculate LNA(IBE())
      DO 5 LI=1,NFALL
      W=0.D0
      DO 3 LJ=1,NFALL
3     W=W+LNKMI(LJ)*RUTA(LJ,LI)
      IA=IBE(LI)
      LNA(IA)=W
5     CONTINUE
      RETURN
      END
C8 -----------------------------------------------------------------
C           LNABER - Part 2
      SUBROUTINE HALTA8
Calculate Totmi and Totva when a solid phase is assumed
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONT,NOLL
      INTEGER HUR
      DOUBLE PRECISION LBETA,P
      COMMON /HLTF0/ CONT,NA,NX,NF,LBETA(190),P(12,190),HUR(12),
     + NOLL(202)
      DOUBLE PRECISION LOGA,TOT,TOL,C,LOGF
      COMMON /HLTF1/ IERR,LOGA(202),TOT(12),TOL(12),C(202),LOGF(162)
      DOUBLE PRECISION LN10,LNBETA,LNKF,LNKMI,LNA,LNBA,LNG
      LOGICAL SINFAL,MONO,BER,FALLA,FALL
      COMMON /HLTF3/ LN10,INDIK,NVA,NFALL,NVAF,NUT,SINFAL,LNBETA(150),
     + LNKF(40),LNKMI(40),LNA(202),MONO(12),TOLY(12),LNG(162),FALL(40),
     + TOTVA(12),PVA(12,150),RUTA(40,40),RUT1(40,40),LNBA(150),IBE(12),
     + IFALL(12),IVAF(12),BER(12),FALLA(12),IVA(12),IVABRA(12),
     + IVANOV(12),ITERS,ITERC,IUT,NION,FSCAL(40)
      DOUBLE PRECISION TOTMI(12)
Calculate Totmi
      DO 5 LI=1,NFALL
      IA=IBE(LI)
5     TOTMI(LI)=TOT(IA)-C(IA)
Calculate TOTVA
      DO 7 LI=1,NVAF
      IA=IVAF(LI)
      W=TOT(IA)
      DO 6 LJ=1,NFALL
6     W=W-RUT1(LI,LJ)*TOTMI(LJ)
7     TOTVA(IA)=W
      RETURN
      END
C6 ----------------------------------------------------------------------------
C          LNBAS
      SUBROUTINE HALTA6(IVAR,X)
Calculates those parts of the activity of the soluble complexes that are
C independent of Lna(Ivar), (the Lna varied) and stores them in Lnba(ix).
C This subroutine is called once every time the value of Ivar is changed.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONT,NOLL
      INTEGER HUR
      DOUBLE PRECISION LBETA,P
      COMMON /HLTF0/ CONT,NA,NX,NF,LBETA(190),P(12,190),HUR(12),
     + NOLL(202)
      DOUBLE PRECISION LOGA,TOT,TOL,C,LOGF
      COMMON /HLTF1/ IERR,LOGA(202),TOT(12),TOL(12),C(202),LOGF(162)
      DOUBLE PRECISION LN10,LNBETA,LNKF,LNKMI,LNA,LNBA,LNG
      LOGICAL SINFAL,MONO,BER,FALLA,FALL
      COMMON /HLTF3/ LN10,INDIK,NVA,NFALL,NVAF,NUT,SINFAL,LNBETA(150),
     + LNKF(40),LNKMI(40),LNA(202),MONO(12),TOLY(12),LNG(162),FALL(40),
     + TOTVA(12),PVA(12,150),RUTA(40,40),RUT1(40,40),LNBA(150),IBE(12),
     + IFALL(12),IVAF(12),BER(12),FALLA(12),IVA(12),IVABRA(12),
     + IVANOV(12),ITERS,ITERC,IUT,NION,FSCAL(40)
      DOUBLE PRECISION X
      X=LNA(IVAR)
      DO 1 LIX=1,NX
      LNBA(LIX)=LNBETA(LIX)
      DO 1 LI=1,NA
      IF(LI.NE.IVAR)  LNBA(LIX)=LNBA(LIX)+P(LI,LIX)*LNA(LI)
1     CONTINUE
      RETURN
      END
C7 --------------------------------------------------------------------------
C             TOTBER
      SUBROUTINE HALTA7
Calculation of Y and comparison with Y0
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONT,NOLL
      INTEGER HUR
      DOUBLE PRECISION LBETA,P
      COMMON /HLTF0/ CONT,NA,NX,NF,LBETA(190),P(12,190),HUR(12),
     + NOLL(202)
      DOUBLE PRECISION LOGA,TOT,TOL,C,LOGF
      COMMON /HLTF1/ IERR,LOGA(202),TOT(12),TOL(12),C(202),LOGF(162)
      DOUBLE PRECISION LN10,LNBETA,LNKF,LNKMI,LNA,LNBA,LNG
      LOGICAL SINFAL,MONO,BER,FALLA,FALL
      COMMON /HLTF3/ LN10,INDIK,NVA,NFALL,NVAF,NUT,SINFAL,LNBETA(150),
     + LNKF(40),LNKMI(40),LNA(202),MONO(12),TOLY(12),LNG(162),FALL(40),
     + TOTVA(12),PVA(12,150),RUTA(40,40),RUT1(40,40),LNBA(150),IBE(12),
     + IFALL(12),IVAF(12),BER(12),FALLA(12),IVA(12),IVABRA(12),
     + IVANOV(12),ITERS,ITERC,IUT,NION,FSCAL(40)
      COMMON /HLTF2/ IVAR,X,Y,Y0,X1(12),X2(12),Y1(12),Y2(12),KARL(12)
     + ,STEG(12),STEP0,ITER(12)
      INTEGER ITERMA
      DATA ITERMA/50/
      IF(INDIK.EQ.10) GO TO 11
      Y=C(IVAR)
      IF(NFALL.EQ.0)  GO TO 4
C Some solid is assumed to be present
      DO 3 LIX=1,NX
      LIAX=NA+LIX
3     Y=Y+PVA(IVAR,LIX)*C(LIAX)
      Y=Y-TOTVA(IVAR)
      Y0=0.D0
      W=DABS(Y-Y0)
      IF(TOL(IVAR).GT.0.D0.AND.DABS(TOTVA(IVAR))
     + .GT.1.D0 .AND.W.LT.DABS(TOL(IVAR)*TOTVA(IVAR)))
     +      W=0.D0
      GO TO 6
C No solid phase assumed to be present
4     DO 5 LIX=1,NX
      LIAX=NA+LIX
5     Y=Y+P(IVAR,LIX)*C(LIAX)
      Y0=TOT(IVAR)
      W=DABS(Y-Y0)
Compare Y with Y0
6     IF(TOLY(IVAR).GT.W)  GO TO 11
C It was not OK
      IF(MONO(IVAR).AND.NFALL.EQ.0)  GO TO 10
7     INDIK=3
      ITER(IVAR)=ITER(IVAR)+1
      IF(ITER(IVAR).GE.ITERMA)  GO TO 111
      RETURN
C  Mononuclear component
10    IF(Y0.LE.0.D0  .OR.  Y.LE.0.D0) GO TO 7
      LNA(IVAR)=LNA(IVAR)+DLOG(Y0)-DLOG(Y)
      X=LNA(IVAR)
      INDIK=1
      ITER(IVAR)=ITER(IVAR)+1
      IF(ITER(IVAR).GE.ITERMA)  GO TO 111
      RETURN
C  It is OK
11    KARL(IVAR)=1
      INDIK=2
      ITER(IVAR)=0
      STEG(IVAR)=STEP0
      RETURN
C Too many iterations done with component Ivar
111   CONTINUE
      IF(IVABRA(IVAR).NE.0) GO TO 11
      IERR=-1
      IF(IUT.LE.0) GO TO 11
      IF(ITER(IVAR).NE.ITERMA)  GO TO 113
      WRITE(IUT,28) IVAR,(IVA(LI),LI=1,NVA+1)
      WRITE(IUT,30) (TOT(LI),LI=1,NA)
      WRITE(IUT,26) (C(LI),LI=1,NION+NF)
113   WRITE(IUT,25) IVAR,ITER(IVAR),(LNA(LI),LI=1,NA)
      WRITE(IUT,31) (LNG(LI),LI=1,NA)
      WRITE(IUT,26) (C(LI),LI=1,NION+NF)
      IF(NFALL.NE.0) GO TO 112
      WRITE(IUT,27) Y, Y0
      IF(.NOT.MONO(IVAR))
     + WRITE(IUT,29) X1(IVAR),X2(IVAR),Y1(IVAR),Y2(IVAR)
      GO TO 114
112   WRITE(IUT,41) NFALL,W
      WRITE(IUT,42) X1(IVAR),X2(IVAR),Y1(IVAR),Y2(IVAR)
114   CONTINUE
      IF(ITER(IVAR).GT.ITERMA+1)  GO TO 11
      RETURN
25    FORMAT(/,' Component:',I3,'   Iteration:',I3,/,
     +'  LnA()=',21(/5X,1P5G15.7))
26    FORMAT('  C()=',21(/5X,1P5G15.7))
27    FORMAT('   Tot(Calc)=',1PG17.9,'  Tot(Input)=',G17.9)
28    FORMAT(/' ?? Error: Too Many Iterations With Component:'I3,/,
     +'    Component numbers in the order they are tested:',/5X,12I4)
29    FORMAT('         Low LnA=',1PG23.16,'  High LnA=',G23.16
     + ,/,'   Low Tot(Calc)=',G17.9,'  High Tot(Calc)=',G17.9)
30    FORMAT('  TOT()=',21(/5X,1P5G15.7))
31    FORMAT('  ln f()=',21(/5X,1P5G15.7))
41    FORMAT('  No.Solids=',I5,'      Error in Tot.Conc.=',1PG17.9)
42    FORMAT('         Low LnA=',1PG23.16,'  High LnA=',G23.16
     + ,/,'   Low Err.Tot(Calc)=',G17.9,'  High Err.Tot(Calc)=',G17.9)
      END
C      SUBROUTINE FACTOR(NION,C,LNF)
C      LOGICAL GOOD
C      INTEGER NION
C      DOUBLE PRECISION C(NION),LNF(NION)
C      RETURN
C      END
CC Test of HALTA routine
CC Chemical Components are: H+, e-, Fe(s)
CC Soluble Complexes are: Fe2+, FeOH+, FeOH3 -, FeOH4 2-,
CC       Fe3+, FeOH 2+, FeOH2 +, FeOH3, Fe2OH2 4+, Fe3OH4 5+,
CC       OH-, O2(g), H2(g),
CC Solids are: Fe(s), FeOH2(s), Fe3O4(s), FeOOH(s)
CC - - - - - - - - - - -
C      LOGICAL CONT
C      INTEGER HUR
C      LOGICAL NOLL
C      DOUBLE PRECISION LBETA,P
C      COMMON /HLTF0/ CONT,NA,NX,NF,LBETA(50),P(12,50),HUR(12),NOLL(62)
C      DOUBLE PRECISION LOGA,TOT,TOL,C,LOGF
C      COMMON /HLTF1/ IERR,LOGA(62),TOT(12),TOL(12),C(62),LOGF(52)
CC - - - - - - - - - - -
C      DATA NA,NX,NF/3,13,4/, LBETA/14.9,5.1,-16.88,-30.8,2.43,-0.62,
C     1 -3.88,-10.49,1.9,1.52,-14.2,-84.49,-1.39,0.,-1.52,-8.6,1.12,
C     1 33*0./
C      DATA P/ 0.,-2.,1.,9*0., -1.,-2.,1.,9*0.,
C     1       -3.,-2.,1.,9*0., -4.,-2.,1.,9*0.,
C     2        0.,-3.,1.,9*0., -1.,-3.,1.,9*0.,
C     3       -2.,-3.,1.,9*0., -3.,-3.,1.,9*0.,
C     4       -2.,-6.,2.,9*0., -4.,-9.,3.,9*0.,
C     5       -1., 0.,0.,9*0., -4.,-4.,0.,9*0.,
C     6        2., 2.,0.,9*0.,
C     7        0., 0.,1.,9*0., -2.,-2.,1.,9*0.,
C     8       -8.,-8.,3.,9*0., -3.,-3.,1.,  405*0./
C      DATA TOL/12*1.D-5/,CONT/.FALSE./,IERR/5/,IUT/5/,
C     1 HUR/12*1/,TOT/1.D-5, -0.15, 0.1, 9*0./,NOLL/62*.FALSE./
CC  Fe(s) is a solid component
C      NOLL(3)=.TRUE.
CC  e- are not stable in aqueous solutions
C      NOLL(2)=.TRUE.
C      CALL HALTA
C      WRITE(IUT,28) IERR,(TOT(LI),LI=1,NA)
C      WRITE(IUT,25) (LOGA(LI),LI=1,NA+NX+NF)
C      WRITE(IUT,26) (C(LI),LI=1,NA+NX+NF)
C25    FORMAT('    LOGA()=',21(/8X,4G14.4))
C26    FORMAT('    C()=',21(/8X,1P4G14.4))
C28    FORMAT(/,' RESULTS:',/,' IERR = ',I3,/,
C     1 '    TOT()=',21(/8X,1P4G14.4))
C      STOP
C      END
C      SUBROUTINE FACTOR(GOOD,NION,LOGA,C,LOGF)
C      LOGICAL GOOD
C      INTEGER NION
C      DOUBLE PRECISION LOGA(NION),C(NION),LOGF(NION)
C      RETURN
C      END
CC Example of use of Subroutine HALTA
CC   Chemical Components: H2O, H+, Na+, Cl-, Ca2+, SO4 2-.
CC   Soluble Complexes: OH-, CaSO4, HSO4 -, CaOH+, NaSO4 -.
CC   Solid Complexes: Gypsum (CaSO4.2H2O), Ca(OH)2(s)
CC - - - - - - - - - - -
C      LOGICAL CONT
C      INTEGER HUR
C      LOGICAL NOLL
C      DOUBLE PRECISION LBETA,P
C      COMMON /HLTF0/ CONT,NA,NX,NF,LBETA(50),P(12,50),HUR(12),NOLL(62)
C      DOUBLE PRECISION LOGA,TOT,TOL,C,LOGF
C      COMMON /HLTF1/ IERR,LOGA(62),TOT(12),TOL(12),C(62),LOGF(52)
CC - - - - - - - - - - -
C      DOUBLE PRECISION Z,IONSTR,OSMCF,ELBLNC
C      COMMON /USER/ Z(52),IONSTR,OSMCF,ELBLNC
CC - - - - - - - - - - -
C      DATA NA,NX,NF/6,5,2/, HUR/2,1,1,1,1,1,6*0/,
C     1 LBETA/-13.998,2.309,1.987,-12.598,0.7,-4.602,22.80, 43*0./
CC            H2O  H+ Na+ Cl- Ca2+ SO4 2-
C      DATA P/ 0.,-1., 0., 0.,  0.,  0.,   6*0.,
C     1        0., 0., 0., 0.,  1.,  1.,   6*0.,
C     1        0., 1., 0., 0.,  0.,  1.,   6*0.,
C     1        1.,-1., 0., 0.,  1.,  0.,   6*0.,
C     1        0., 0., 1., 0.,  0.,  1.,   6*0.,
C     1        2., 0., 0., 0.,  1.,  1.,   6*0.,
C     1        2.,-2., 0., 0.,  1.,  0.,   522*0./
C      DATA NOLL/62*.FALSE./,LOGA/0.,-7.,0.,0.,-5.,-5., 56*0./,
C     1 Z/0.,1.,1.,-1.,2.,-2.,  -1.,0.,-1.,1.,-1.,  41*0./,
C     1 C/62*0./,LOGF/52*0./,TOL/12*1.D-5/,IUT/5/
CC Total concentrations:
C      DATA TOT/0., 0., 1., 1., 0.05, 0.05,  6*0./
C      NOLL(1)=.TRUE.
CC The ion pair NaSO4 - is not considered in the mass balance.
C      NOLL(11)=.TRUE.
C      CONT=.FALSE.
C      IERR=5
C      CALL HALTA
C      WRITE(IUT,28) IERR,IONSTR,(TOT(LI),LI=1,NA)
C      WRITE(IUT,25) (LOGA(LI),LI=1,NA+NX+NF)
C      WRITE(IUT,26) (C(LI),LI=1,NA+NX+NF)
C      WRITE(IUT,27) (LOGF(LI),LI=1,NA+NX)
C25    FORMAT('    LOGA()=',21(/8X,4G14.4))
C26    FORMAT('    C()=',21(/8X,1P4G14.4))
C28    FORMAT(/,' RESULTS:',/,' IERR = ',I3,/,
C     1 ' IONIC STRENGTH=',1PG14.4,/,'    TOT()=',21(/8X,1P4G14.4))
C27    FORMAT('    LOGF()=',21(/8X,1P4G14.4))
C      STOP
C      END
CC
C      SUBROUTINE FACTOR(GOOD,NION,LOGA,C,LOGF)
CCalculation of Activity Factors
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      LOGICAL GOOD
C      INTEGER NION
C      DOUBLE PRECISION ELBLNC,IONSTR,OSMCF,Z,LOGA(NION),C(NION),
C     1 LOGF(NION),OLDI,W,LN10
C      COMMON /USER/ Z(52),IONSTR,OSMCF,ELBLNC
C      DATA OLDI/0.0D0/,LN10/2.3025850929940D0/
CC -------------------------------------------------------------------------
CC   
CCalculate: Electrical balance of the solution and its Ionic strength.
CCheck if the Ionic Strenth has varied since it was last calculated.
CC  If it has varied, calculate activity coeff. and osmotic coeff.
CC  otherwise, set variable GOOD=.TRUE.
CC
CCalculate Electrical balance
C      ELBLNC=0.D0
C      DO 1 I=NION
C1     ELBLNC=ELBLNC+Z(I)*C(I)
CCalculate Ionic Strength.
C      IONSTR=DABS(ELBLNC)
C      DO 2 I=1,NION
C2     IONSTR=IONSTR + Z(I)*Z(I)*C(I)
C      IONSTR=0.5D0 * IONSTR
C      IF(IONSTR.GT.500.D0) IONSTR=500.D0
CCheck if the Ionic Strength is Constant.
C      DIFF=IONSTR-OLDI
C      OLDI=IONSTR
C      IF(IONSTR.GT.0.D0) DIFF=DIFF/IONSTR
C      IF(DIFF.LE.0.001) GOOD=.TRUE.
C      IF(GOOD) RETURN
CCalculate Activity coeff. and osmotic coeff.
CC
CCalculate Debye-Huckel Term
C      ROOTI=DSQRT(IONSTR)
C      DOWN= 1.0D0 + 1.5D0 * ROOTI
C      DH= -0.5107D0 * ROOTI / DOWN
CCalculate Logf (log of activity coeff.)
C      DO 10 I=1,NION
C      IF(Z(I).EQ.0.D0) GO TO 12
C      SUM=0.D0
C      IF((Z(I)*ELBLNC).GT.0.D0) SUM= 0.1D0 * Z(I) * ELBLNC
C         DO 11 J=1,NION
C         ZP=Z(I)*Z(J)
C         IF(ZP.GE.0.D0) GO TO 11
C         SUM=SUM - 0.1D0 * ZP * C(J)
C11       CONTINUE
C      LOGF(I)= Z(I)*Z(I)*DH + SUM
C      GO TO 14
C12    LOGF(I)= 0.1D0*IONSTR
C14    IF(LOGF(I).GT.15.D0) LOGF(I)=15.D0
C10    CONTINUE
CCalculate Osmotic Coeff.
CC   Debye-Huckel Term.
C      PHIDH= -0.69685D0 * (DOWN-2.D0*DLOG(DOWN)-(1.D0/DOWN))
CCalculte sum of ions and sum of products of conc. times specific inter. term
C      SUM=DABS(ELBLNC)
C      SUMPRD=0.D0
C      IF(ELBLNC.GE.0.D0) GO TO 18
C        DO 17 J=1,NION
C        IF(Z(J).GE.0.D0) GO TO 17
C        SUMPRD=SUMPRD + 0.23D0 *ELBLNC*Z(J)*C(J)
C17      CONTINUE
C18    DO 20 I=1,NION
C      IF(Z(I).EQ.0.D0) GO TO 20
C      SUM=SUM+C(I)
C      IF(Z(I).LT.0.D0) GO TO 20
C      IF(ELBLNC.GT.0.D0) SUMPRD=SUMPRD + 0.23D0 * Z(I)*C(I)*ELBLNC
C         DO 19 J=1,NION
C         IF(Z(J).GE.0.D0) GO TO 19
C         ZP=Z(I)*Z(J)
C         SUMPRD=SUMPRD - 0.23D0 * ZP * C(I) * C(J)
C19       CONTINUE
C20    CONTINUE
CCalculate Osmotic Coeff.
C      OSMCF= 1.D0 + (( PHIDH + SUMPRD) / SUM)
C      IF(OSMCF.GT.5.D0) OSMCF=5.D0
CCalculate water activity.
C      W = -OSMCF * SUM / 55.51
C      LOGA(1)=W/LN10
C      LOGF(1)= OSMCF
C      RETURN
C      END

