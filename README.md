# HaltaFall
**HaltaFall** (Swedish: Halt = concentration, Falla= to precipitate)
is an old computer program which calculates the equilibrium concentrations of the species in mixtures of any number of components, which can form any number of complexes and solid phases, provided concentrations (or partial pressures) can be used in the equilibrium calculations. The equilibrium constants must be known and enough data be given about the gross composition. From the input information, the program itself devises an efficient plan for solving the simultaneous equations. The program also provides an efficient means of finding out which of many possible solid phases can actually appear in a certain equilibrium mixture

HaltaFall (and [LetaGrop](https://github.com/ignasi-p/LetaGrop))
were developed by Lars Gunnar Sillén and his co-workers at
the Department of Inorganic Chemistry, Royal Institute of
Technology ([KTH](https://www.kth.se/en)), Stockholm, Sweden.

The initial Fortran version by
[Ekelund et al. (1970)](#ekelund-r-sillén-l-g-and-wahlberg-o-1970)
appears to have been lost.

The algorithm, in the form of a [Fortran](https://en.wikipedia.org/wiki/Fortran)
subroutine from 1984 may be found in this repository.  A later
version written in Fortran 90 is provided here, together with
a program **"EC"** that calls the HaltaFall subroutine and
performs chemical equilibrium calculations.

A list of publications describing HaltaFall is given at the
end of this document. Instructions are also given here on how to download the Fortran compiler from
[`MinGW64`](https://www.mingw-w64.org/),
how to create the **"EC"** Windows binary (`exe`-file),
and how to run it and test it.

### Alternative software

See for example:
- https://en.wikipedia.org/wiki/Geochemical_modeling#Software_programs_in_common_use
- https://gemshub.github.io/site/
- https://www.eawag.ch/en/department/surf/projects/chemeql/
- http://www.hyperquad.co.uk/

## Download
All files are available from the [releases section](https://github.com/ignasi-p/letaGrop/releases/latest).

## References

#### [Ingri N, Kakołowicz W, Sillén L G and Warnqvist B (1967)](https://doi.org/10.1016/0039-9140(67)80203-0)
High-speed computers as a supplement to graphical
methods - V. HALTAFALL, a general program for calculating
the composition of equilibrium mixtures.  
_Talanta_ **14,** 1261–1286.

#### [Ingri N, Kakołowicz W, Sillén L G and Warnqvist B (1968)](https://doi.org/10.1016/0039-9140(68)80071-2)
Errata. _Talanta_ **15,** xi–xii.

#### [Elgquist B (1969)](https://doi.org/10.1016/0039-9140(69)80208-0)
A Fortran version of Haltafall for computing ionic equilibria.  
_Talanta_ **16,** 1502–1503.


#### [Anfält T and Jagner  D (1969)](https://doi.org/10.1016/S0003-2670(01)81805-5)
Interpretation of titration curves by means of the computer
program HaltaFall,  
_Anal. Chim. Acta_ **47,** 57-69.

#### [Ekelund R, Sillén L G and Wahlberg O (1970)](https://doi.org/10.3891/acta.chem.scand.24-3073)
Fortran editions of Haltafall and Letagrop.  
_Acta Chem. Scand._ **24,** 3073.

#### [Warnqvist B and Ingri N (1971)](https://doi.org/10.1016/0039-9140(71)80069-3)
The HALTAFALL program - some corrections, and comments on
recent experience.  
_Talanta_ **18,** 457–458.

------------
### See also
Ignasi's page on [water chemistry](https://sites.google.com/view/groundwatergeochemistry).
