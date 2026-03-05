# Running EC

EC is a [command-line](https://en.wikipedia.org/wiki/Command-line_interface)
Fortran program that uses HaltaFall
[(Ingri et al., 1967)](https://doi.org/10.1016/0039-9140(67)80203-0)
as the calculation algorithm.

The file `Fortran_Windows.md` describes how to download and install
the [`MinGW64`](https://www.mingw-w64.org/) Fortran compiler for
64 bit Windows.

The file `compile_and_test.md` describes how to build and test EC.

The user must provide the input- and output-file names to EC
on the command-line.  Usage:  
```
EC  Input_File_Name  [-command:value]
```
Enclose `Input_File_Name` in double quotes ("file name")
if it contains blank space.  For a list of possible commands
enter:  `EC  -?`

Two output text files are created: a verbose output file
(with extension "out") and a table output file (with extension
"csv") where the calculation results are reported column-wise
with data separated by semicolons (";").
The table-output file may be imported into Excel or any other
plotting program.

## The input file

The input disk file may be a modified Medusa/Spana "dat"-file.
The format of the input file:

| Input:                                   | Comments:   |
| -------                                  | ----------- |
| `NA,NX,NSOL,SOLIDC,`                     |             |
| `IDENT(ia),`                             |(given for `ia`=1,`NA`)|
| `IDENT(ix),LOGBTA(ix),(A(ix,j),j=1,NA),` |(for `ix`=1,`NX`)      |
| `IDENT(if),LOGBTA(if),(A(if,j),j=1,NA),` |(for `if`=1,`NSOL`)    |
| `NPKT,`                                  |                       |
| `TEXT`                                   | any number of lines of text which will be printed as a heading |

Where:  
 - `NA`=number of components
 - `NX`=number of soluble reaction products
 - `NSOL`=number of solid reaction products
 - `SOLIDC`=number of components which are solids
 - `IDENT(i)`=identification for species `i`.  
    (the last `SOLIDC` components are assumed to be solids).
 - `LOGBTA(i)`=log(Beta) for reaction `i` (Beta=global equilibrium
            constant of formation)
 - `A(i,j)`=formula units for species `i` and component `j`.  
   That is, the stoichiometric coefficient for component `j`
   in reaction product `i`.
 - `NPKT`=number of points to calculate (equilibrium compositions).  
    May be negative (see text below).
 - `HOW`=one of: `T`, `TV`, `LTV`, `LA`, `LAV`.  
    - if `HOW`=`T` the total concentration for that component
      follows: `CONC(ia)`.  
    - if `HOW`=`TV` the total concentration for that component
      is varied within the limits that follow: `LOW(ia)`
      and `HIGH(ia)`.
    - if `HOW`=`LTV` the Log(Total Conc.) for that component
      is varied within the limits that follow: `LOW(ia)`
      and `HIGH(ia)`.
    - if `HOW`=`LA` the Log(activity) for that component
      follows: `CONC(ia)`.
    - if `HOW`=`LAV` the Log(activity) for that component is
      varied within the limits that follow: `LOW(ia)`
      and `HIGH(ia)`.

 _Notes:_  
 1. If `IDENT(i)` is either of: `E-`, `e-`, `H2O` or `h2o`,  
    or if `IDENT(i)` begins with `*` (for example:
    `*Cl2(g)`), then the concentration of that species is
    **not** taken into account (the conc. of aqueous
    electrons etc, is set equal to zero), although the activity
    for such species is calculated.
 2. If `IDENT(i)` is equal to either `H2O` or `h2o`, then the
    activity for that component is calculated using the osmotic
    coefficient, and its concentration is set equal to zero.
 3. If `NPKT` is less than zero, then `abs(NPKT)` points will be
    calculated, and the input should be:

    | Input:                 | Comments:   |
    | -------                | ----------- |
    | `-NPKT,`               | negative value of the number of <br> chemical compositions to calculate |
    | `HOW, CONC(ia),`       |(for `ia`=1,`NA`) for point 1 |
    | `HOW, CONC(ia),`       |(for `ia`=1,`NA`) for point 2 |
    |  &emsp; :  &emsp; &emsp; : |                          |
    | `HOW, CONC(ia),`       |(for `ia`=1,`NA`) for point `NPKT`|
    | `TEXT`                 | any number of lines of text which will be <br> printed as a heading |

### Input Example 1

A chemical system with 3 components and 8 reactions, of which
2 solids.

11 calculation points will be made at between pH = 4 and 8,
that is, at pH = 4.0, 4.4, 4.8, ... 7.2, 7.6, and 8.0.
The total concentrations = 0.01 for both carbonate and calcium.

Input file example:
<pre>3, 6, 2, 0,
H+
CO3 2-
Ca+2
OH-    , -14.0   -1  0  0
HCO3-  ,  10.327  1  1  0
CO2    ,  16.68   2  1  0
CaCO3  ,  3.224   0  1  1
CaHCO3+,  11.435  1  1  1
CaOH+  , -12.57  -1  0  1
CaCO3(cr)  ,  8.48    0  1  1
Ca(OH)2(cr), -22.75  -2  0  1
11
LAV,-4,-8,
T,0.01,
T,0.01,
Example nr 1
</pre>

### Input Example 2

The same chemical system as above, but including water (H2O)
as a component.  
Three calculations are requested, at pH = 4, 6 and 8.  
The total concentrations of carbonate and calcium are both 0.01.  

Note that the activity of water (H2O) is given (`LA,0`),
but it will be nevertheless calculated using the activity
coefficient model selected when running the EC program.  

Input file example:
<pre> 4, 6, 2, 0,
H+
CO3 2-
Ca+2
H2O
OH-         , -14.0     -1  0  0  1
HCO3-       ,  10.327    1  1  0  0
CO2         ,  16.68     2  1  0 -1
CaCO3       ,  3.224     0  1  1  0
CaHCO3+     ,  11.435    1  1  1  0
CaOH+       , -12.57    -1  0  1  1
Ca(OH)2(cr) , -22.75    -2  0  1  2
CaCO3(cr)   ,  8.48      0  1  1  0
-3
LA,-4, T,0.01, T,0.01, LA,0
LA,-6, T,0.01, T,0.01, LA,0
LA,-8, T,0.01, T,0.01, LA,0
Example nr 2
</pre>