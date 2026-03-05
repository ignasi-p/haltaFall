@echo off
echo *************************************************************************
echo ****  Testing ECf (HaltaFall)  ****
echo *************************************************************************
echo.
echo *************************************************************************
echo ****  Tests using Davies eqn.
echo *************************************************************************
set _td_=.\tests-ECf\Davies_eqn
set _ts_=.\tests-ECf\SIT
echo on

EC "%_td_%\Fe_ECf.dat"  -m=0 -i=-1

EC %_td_%\CO3-Ca-NaCl_ECf.dat -m=0 -i=-1 -tol=1e-9 -tola=1e-8 -t=25

@echo off
echo *************************************************************************
echo **** Tests using the SIT model
echo *************************************************************************
echo.
echo using different command-line options...

set _name_=CaSO4-Na-Cl_ECf
echo on
EC %_ts_%\%_name_%.dat -m=1 -i=-1 -path="%_ts_%"
@echo off
copy /y "%_ts_%\%_name_%.csv" "%_ts_%\%_name_%(1).csv" >nul
del "%_ts_%\%_name_%.csv"

echo on
EC %_ts_%\%_name_%.dat -m=1 -i=-1 -tbls=; -path="%_ts_%"
@echo off
copy /y "%_ts_%\%_name_%.csv" "%_ts_%\%_name_%(2).csv" >nul
del "%_ts_%\%_name_%.csv"

echo on
EC %_ts_%\%_name_%.dat -m=1 -i=-1 -tbls=\t -path="%_ts_%"
@echo off
copy /y "%_ts_%\%_name_%.csv" "%_ts_%\%_name_%(3).csv" >nul
del "%_ts_%\%_name_%.csv"

echo on
EC %_ts_%\%_name_%.dat -m=1 -i=-1 -dbg -tbls= -path="%_ts_%"
@echo off
copy /y "%_ts_%\%_name_%.csv" "%_ts_%\%_name_%(4).csv" >nul
del "%_ts_%\%_name_%.csv"

:xit
@echo off
set _td_=
set _ts_=
set _name_=
echo.
echo *************************************************************************
echo **** End of testing ECf (HaltaFall)
pause
