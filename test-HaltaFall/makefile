# compiler
FC = gfortran

# compile flags
FCFLAGS = -I ./obj -J ./obj -fbounds-check -fbacktrace -W -Wunused -ffpe-trap=invalid,zero,overflow -ffpe-summary=none
# Using optimization (flag "-O") may give quite different results...
# link flags
FLFLAGS = -I ./obj -static

# "-I dir" Affects interpretation of the INCLUDE directive. This path is also used to
#   search for .mod files when previously compiled modules are required by a USE statement
# "-J dir" This option specifies where to put .mod files for compiled modules. It is
#   also added to the list of directories to searched by an USE statement.
#
# Note: hyphens not allowed in the "exe" name

go: go2 Test_HaltaFall

go2:  # write a heading to the user
	@echo Usage is:    mingw32-make  [Test_HaltaFall]  [clean]

clean:
	rmDir /q /s .\obj
	del /q Test_HaltaFall.exe
	del /q phreeqc.log

# ---- "Test_HaltaFall" does not need SIT nor Read_IR, but Factor_Module
#      needs Read_IR to be able to read the SIT-coefficients, etc
# These are the dendencies:
#   Read_IR+IO = READ  --> CHEM, SIT
#   CHEM --> SIT, FACTOR, HALTAFALL, Test_HaltaFall
#   SIT --> FACTOR
#   FACTOR --> HALTAFALL, Test_HaltaFall
#   HALTAFALL --> Test_HaltaFall


Test_HaltaFall: Test_HaltaFall.exe Modules 

Modules:  ./obj/Read_IR.o  ./obj/Chem_Module.o  ./obj/SIT.o  ./obj/Factor_Module.o  ./obj/HaltaFall.o

Test_HaltaFall.exe : ./obj/Test_HaltaFall.o ./obj/HaltaFall.o ./obj/haltaf_module.mod \
		./obj/Factor_Module.o ./obj/factor_module.mod ./obj/Chem_Module.o ./obj/chem.mod \
		./obj/SIT.o ./obj/sit.mod ./obj/Read_IR.o ./obj/readir.mod ./obj/io.mod
	$(FC) $(FLFLAGS) -o Test_HaltaFall.exe -I ./obj ./obj/Test_HaltaFall.o ./obj/HaltaFall.o ./obj/Read_IR.o ./obj/Chem_Module.o \
		./obj/Factor_Module.o ./obj/SIT.o


./obj/Test_HaltaFall.o : ./Test_HaltaFall.f90 ./obj/HaltaFall.o ./obj/haltaf_module.mod \
		./obj/Factor_Module.o ./obj/factor_module.mod ./obj/Chem_Module.o ./obj/chem.mod \
		./obj/SIT.o ./obj/sit.mod ./obj/Read_IR.o ./obj/readir.mod ./obj/io.mod
	$(FC) $(FCFLAGS) -o ./obj/Test_HaltaFall.o -c ./Test_HaltaFall.f90

./obj/Chem_Module.o : ../src/Chem_Module.f90 ./obj/Read_IR.o ./obj/readir.mod  ./obj/io.mod
	$(FC) $(FCFLAGS) -o ./obj/Chem_Module.o -c ../src/Chem_Module.f90

./obj/HaltaFall.o : ../src/HaltaFall.f90 ./obj/Chem_Module.o ./obj/chem.mod ./obj/Factor_Module.o ./obj/factor_module.mod
	$(FC) $(FCFLAGS) -o ./obj/HaltaFall.o -c ../src/HaltaFall.f90

./obj/SIT.o : ../src/SIT.f90 ./obj/Read_IR.o ./obj/io.mod ./obj/readir.mod ./obj/Chem_Module.o ./obj/chem.mod
	$(FC) $(FCFLAGS) -o ./obj/SIT.o -c ../src/SIT.f90

./obj/Factor_Module.o : ../src/Factor_Module.f90 ./obj/Chem_Module.o ./obj/chem.mod ./obj/SIT.o ./obj/sit.mod
	$(FC) $(FCFLAGS) -o ./obj/Factor_Module.o -c ../src/Factor_Module.f90

./obj/Read_IR.o : ../src/Read_IR.f90
	$(FC) $(FCFLAGS) -o ./obj/Read_IR.o -c ../src/Read_IR.f90
