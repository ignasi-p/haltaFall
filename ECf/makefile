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

go: go2 EC

go2:  # write a heading to the user
	@echo Usage is:    nmake  [EC]  [clean]

clean:
	del /q .\obj\*.*

EC: Modules EC.exe

Modules:  ./obj/Read_IR.o  ./obj/Chem_Module.o  ./obj/SIT.o  ./obj/Factor_Module.o  ./obj/HaltaFall.o

EC.exe : ./obj/EC.o ./obj/Chem_Module.o ./obj/HaltaFall.o ./obj/Factor_Module.o ./obj/SIT.o ./obj/Read_IR.o
	$(FC) $(FLFLAGS) -o EC.exe -I ./obj ./obj/EC.o ./obj/Read_IR.o ./obj/Chem_Module.o ./obj/HaltaFall.o \
    ./obj/Factor_Module.o ./obj/SIT.o


# Read_IR+IO = READ  --> CHEM, SIT, EC
# CHEM --> SIT, FACTOR, HALTAFALL, EC
# SIT --> FACTOR
# FACTOR --> HALTAFALL, EC
# HALTAFALL --> EC

./obj/Chem_Module.o : ../src/Chem_Module.f90 ./obj/Read_IR.o
	$(FC) $(FCFLAGS) -o ./obj/Chem_Module.o -c ../src/Chem_Module.f90

./obj/HaltaFall.o : ../src/HaltaFall.f90 ./obj/Chem_Module.o ./obj/Factor_Module.o
	$(FC) $(FCFLAGS) -o ./obj/HaltaFall.o -c ../src/HaltaFall.f90

./obj/SIT.o : ../src/SIT.f90 ./obj/Read_IR.o ./obj/Chem_Module.o
	$(FC) $(FCFLAGS) -o ./obj/SIT.o -c ../src/SIT.f90

./obj/Factor_Module.o : ../src/Factor_Module.f90 ./obj/Chem_Module.o ./obj/SIT.o
	$(FC) $(FCFLAGS) -o ./obj/Factor_Module.o -c ../src/Factor_Module.f90

./obj/Read_IR.o : ../src/Read_IR.f90
	$(FC) $(FCFLAGS) -o ./obj/Read_IR.o -c ../src/Read_IR.f90

./obj/EC.o : ./src/EC.f90 ./obj/Read_IR.o ./obj/Chem_Module.o ./obj/HaltaFall.o ./obj/Factor_Module.o
	$(FC) $(FCFLAGS) -o ./obj/EC.o -c ./src/EC.f90
