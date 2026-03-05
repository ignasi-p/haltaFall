# EC - Compile and Test

The file `Fortran_Windows.md` describes how to download and install
the [`MinGW64`](https://www.mingw-w64.org/) Fortran compiler for
64 bit Windows. It is expected that the installation folder
for the [GNU compiler collection](https://en.wikipedia.org/wiki/GNU_Compiler_Collection)
is `C:\bin\mingw64`. If not, you will need to
modify the `*.cmd` files with a text editor.

The source Fortran code and test files are found in sub-folders
`src` and `tests-ECf`, respectively.

This folder contains two files to create the Windows binary
(`exe`) file:  
  `1-compile.cmd` and `makefile`  
Running `1-compile.cmd` will create the `exe`-file based
on the instructions in `makefile`.

This folder also contains a file named `2-run_tests.cmd` that
will run several tests of EC.
