# Internals — Structure of the code

This section gives a brief overview of the general structure of the code. It is
intended as a guide through the source code for deeper inspection of the
numerical internals or for debugging. The overview explains the function of the
different routines and how they interact within the program. It cannot substitute for a study of the sources, which are commented in relative detail.

## Code design

The general design of the program is directed towards high accuracy of the
computed line profiles. All errors in the different steps of the program are
explicitly user-controlled by setting thresholds. All discretizations necessary
to treat the problem numerically are performed in an adaptive way — there is no
predefined grid at all, and all grid parameters may change during the iteration
procedure.

Furthermore, the code was designed for high flexibility, i.e., the ability to
treat a very broad range of physical parameters with the same accuracy and
without numerical limitations. For example, the systematic velocities may range
from 0 to several times the turbulent velocity depending on the compiled
frequency field size. There is no inherent restriction to a particular range of
optical depths — the code has been tested for depths at line centre between
about -5 (moderate masering) and about 5000. However, convergence speed depends
dramatically on the optical depth, the level structure, and the local velocity
gradients.

The program is not optimized for speed. Although it runs about 60 times faster
than its predecessor by E. Krügel, other codes with lower inherent accuracy may
run another factor of 10 faster. Nevertheless, the code is suitable for
interactive work even on a modest PC.

## Main program

`LINETRANS` is the main program. It provides the loops for calling `MICRO` and
`LINE`, the routines for the computation of the molecular level populations and
the observed line profiles.

Both `MICRO` and `LINE` start by reading the parameters specified by the user.
`MICRO` then initializes the grids and starts the loop consisting of the
radiative transfer (`STPMIC`), the lambda iteration (`LAMBDA`, `LEVMIC`) and
the grid adjustment (`ADAPTGRID`). Upon convergence the results may be stored
(`WRITEDENS`).

`LINE` continues by setting up the grids for beam offset and frequency, calls
the radiative transfer for one line (`STPLINE`), the beam integration
(`INTBEAM`) and finally writes the results (`WRITELINE`).

`PRESENT` shows the opening and closing screen of the program.

## Input / Output

The routines for parameter input — `READPHYS`, `READNUM`, `READOBS` — have the
same general structure. One part handles file input and the other interactive
input. For interactive input, an inverse version of the file reader can be used
to produce a parameter file. For numerical parameters (`READNUM`), there is an
additional part allowing use of the built-in default values or the parameters
from a previous computation when called more than once.

The routine for the molecular data file input, `READSPECIES`, does not allow
direct terminal input and contains a translation between the large input fields
(sized to accommodate all reasonable molecular data files) and the internal
fields used in radiative transfer.

`WRITEDENS` stores the level populations in a file (see the
[Usage](usage.md#output-files) section for the file structure). `READDENS` reads
the populations from such a file. `WRITELINE` stores the integrated beam
temperatures in a file with variable structure. For FITS-format line output it
calls `SAVEFITS`, which requires linking with the cfitsio library and the
auxiliary routines `PRINTERROR` and `DELETEFILE` from `fitsadd.f`. The default
parameters for the file structure are set in `PREDEFINE`, which reads the
command line. `IOERROR` is called on file I/O exceptions; it produces the error
messages and controls the program flow.

## Radiative transfer for all lines

These routines perform the radiative transfer computation needed to obtain the
radiation energy density at each point and within each molecular transition.

`STPMIC` contains the loops over all displacement parameters and all radial
points. For each "ray", `LEFTKPL` is called to initialize the intensity fields
(with the background radiation via `UBACK`) and the absorption coefficient at
the edge of the cloud.

To carry out the transfer between two radial grid points, `FIRSTSTEP` is called
first. This subroutine determines the number of radiative transfer steps needed
between those points to keep the shifts in the line profiles sufficiently small.
If more than one step is required, `NEXTSTEP` is called in a loop to compute
the size of further steps.

`TRANSFER` is the core radiative transfer routine. It constrains the frequency
range where the line profile is non-vanishing and computes the new intensity in
that range. For small optical depths a linear approximation of the exponential
functions is used; otherwise the full expression is applied, assuming that
$\kappa$ and $\epsilon$ vary linearly within the integration range.
`CTRANSFER` computes the continuum radiative transfer within a central H II
region, assuming constant opacity and source function in that core.

At each radial point, `STPMIC` calls the integration routines `UNYINT` and
`UZINT` to compute the energy density.

## Radiative transfer for a single line

The subroutines `STPLINE`, `FIRSTLSTP`, `NEXTLSTP`, `LEFTKPL`, `TRANSFL`, and
`CTRANSFL` are reduced versions of `STPMIC`, `FIRSTSTEP`, `NEXTSTEP`, `LEFTKP`,
`TRANSFER`, and `CTRANSFER`, respectively. They compute the radiative transfer
only in a single line instead of all transitions of a molecule, and are used to
compute the observed emergent line profiles. No intensity integration is called.
The auxiliary routine `TAUCENTRAL` provides a control output of the central
opacity at line centre.

## Initialization and adjustment of the grids

`INITRAD` builds the initial radial grid. It checks the gradients in density and
velocity and inserts as many grid points as necessary to meet the conditions for
maximum changes specified by the numerical parameters. Additional points are
included at discontinuities detected in the density structure. To assign physical
parameters (density, temperature, systematic velocity, …) to a radial point,
`VALUES` is called, which computes these quantities from the physical input
parameters.

`INITN` initializes the level populations at the different radial grid points
according to the chosen first-guess approximation. For a central H II region, it
calls `UFROMHII`, which computes the integrated intensity at a given point due to
continuum radiation from the H II region.

`FILEDENS` takes both the initial radial grid and the level populations from a
precomputed file using `READDENS`. It thus replaces both `INITRAD` and `INITN`
when file input is used for the initial guess.

`INITFREQ` scans through all points and determines the local line width, from
which the maximum frequency/velocity step size at each point is computed.
`OVERLAP` checks for line overlap in the given velocity regime.

`ADAPTGRID` is responsible for updating the radial grid so that the differences
between points meet the conditions for maximum level population change and
maximum velocity difference. It calls `ADDPOINTS` to check for and insert
additional points, and `REMOVEPOINTS` to identify and remove superfluous points.
When new grid points are added, `VALUES` is called to compute the gas density
and velocities at the new radii.

`REMOVELEVELS` checks whether any levels can be dropped due to negligible
population, and calls `LINFROMLEV` to rearrange the field of transition
coefficients accordingly.

`MAKEPS` and `PSFORLINE` build the grid of displacement parameters providing
the "rays" on which radiative transfer is calculated. `MAKEPS` is responsible
for the grid used in `STPMIC`; `PSFORLINE` builds the grid for `STPLINE`. The
two grids differ: `MAKEPS` ensures that for a given radial point the relative
distance of the displacement points on the $z$ scale stays below the numerical
parameter set in `READNUM` (see [Numerical parameters](usage.md#numerical-parameters)).
`PSFORLINE` starts with the radial grid and makes it denser for the displacement
grid to achieve the same accuracy in covering the beam. In contrast to `MAKEPS`,
all radial points are used.

## Solution of the balance equations

`MATRIXUP` computes the coefficients and sets up the transition matrix. Together
with the normalization equation, the transition matrix forms the linear system of
balance equations. `LUSOLVER` solves this system by calling `LUDCMP` for LU
decomposition, `LUBKSB` for back substitution, and `MPROVE` for iterative
improvement of the solution.

`LEVMIC` calls either `ULEVMIC` or `DLEVMIC`, which store the energy densities
or level populations at all radial points for three cycles to be used in the
lambda acceleration. After three ordinary lambda iteration cycles, `LAMBDA` or
`DLAMBDA` is called to compute a second-order lambda acceleration step for the
energy densities or level populations. When the acceleration step is sufficiently
small (convergence criterion), `LEVMIC` terminates the iteration. `HIISMTH`
corrects the obtained level populations within the HII region to enable uniform
treatment across all grid points.

`KAPSRC` computes the new absorption coefficients and source functions from the
level populations using the Einstein coefficients from `AIJ`, `BIJ`, and `BJI`.
`EFFKAP` obtains the effective absorption coefficients for the turbulent medium
from the correlation length and clump absorptivities. `ATAU` is the corresponding
reduction function. For a central HII region, `HIIKAP` corrects the opacities at
the innermost point to obtain continuum values.

## The Sobolev approximation

`SOBOLEV` contains the loop over all radial points and all iterations to
determine the local energy densities and level populations in the Sobolev
approximation. It computes the velocity gradients used by `ESCAPE` to determine
the escape probability in a given direction and its first derivative with respect
to $n_j$. `GAUSSOB` integrates the escape probabilities over all directions using
a Gaussian quadrature whose grid points are initialized by `GAULEG`. From the
integrated escape probabilities, `SOBOLEV` computes the local energy densities
used by `SOBMATRIX` to obtain the level populations.

`SOBMATRIX` sets up the differential form of the balance equation system used
for the Newton-Raphson approach to the local level populations (see the
[Theory](theory.md#the-sobolev-approximation) section). It calls `AMATRIX` for
the original form of the matrix (right-hand side) and `DMATRIX` for the
differential form (left-hand side). The system is solved by `LUSOLVER`.
`SOBKAPSRC` obtains the absorptivities and source function and their derivatives
from the level populations. `SOBEFFKAP` computes the effective absorption
coefficient and its derivatives for a clumpy turbulent medium.

## Intensity integration

`UNYINT` sets the limits for the frequency integration at each radial point and
multiplies the intensities with the local line profile (see the energy density
integral in [Theory](theory.md#the-solution-of-the-radiative-transfer-problem)).
It then calls `QUINTGAUSS`, a simple equidistant integration routine, to carry
out the frequency integration. The line profile is computed by `PHI`.

`UZINT` prepares the fields for the $z$ integration (see the same integral).
For ranges with very small steps it calls the linear integration routine
`DINTLIN`; for larger steps it calls the cubic spline integration routine
`DINTCUB`. `DINTCUB` itself uses the cubic spline interpolation routine `SPLINE`.

`INTBEAM` is responsible for integrating the intensity over the telescope beam.
It computes the range of the cloud covered by a Gaussian beam, calls `SPECERR`
for the $\phi$ integration, and finally `DINTCUB` for the $p$ integration.
`SPECERR` computes a numerical approximation to the $\phi$ integral over the
beam (see [Theory](theory.md#computation-of-beam-temperatures)) that is
independent of the intensity.

`TRANSLAT` translates the angular beam width of the telescope into the projected
beam width on the cloud.

## Molecular constants

The following routines contain the full physics of the problem. All physical
quantities are specified within them; all other routines contain only numerics.
The fields used here are read in advance by `READSPECIES`.

`INITEINSTEIN` determines which molecule shall be considered and initializes the
fields containing molecule-specific quantities, especially $a_{j,j-1}$,
$b_{j,j-1}$, and $\nu_0$. The functions `AIJ`, `BIJ`, `BJI` and `CFREQ` simply
read the tables of $A_{j,l}$, $B_{j,l}$, and $\nu_0$ and return them as
function values.

`PROBTHERM`, `VTHERM`, and `UBACK` also use the tables produced by
`INITEINSTEIN`. These functions provide the population of a given level for a
thermal distribution, compute the mean thermal velocity of the molecules at a
given temperature, and determine the intensity of the cosmic background radiation
at the frequencies of the molecular transitions. `HIIINIT` computes the continuum
opacity and source function for all transitions of the considered molecule from
the table of frequencies and the electron density and temperature within a given
HII region.

`CIJ` obtains the collision rate coefficients at a given temperature for the
transition of interest. It calls `ILOCAT` to find the appropriate table entries
and interpolates the values for the given temperature.

## Source files

| File | Routines |
|------|---------|
| `ltr.f` | `LINETRANS`, `MICRO`, `LINE`, `PHI`, `TRANSLAT`, `PRESENT` |
| `ltrio.f` | `READPHYS`, `READNUM`, `READOBS`, `WRITEDENS`, `READDENS`, `WRITELINE`, `SAVEFITS`, `IOERROR`, `PREDEFINE` |
| `initial.f` | `INITRAD`, `INITN`, `FILEDENS`, `INITFREQ`, `ADAPTGRID`, `ADDPOINTS`, `REMOVEPOINTS`, `VALUES`, `UFROMHII`, `OVERLAP` |
| `stpdiff.f` | `STPMIC`, `FIRSTSTEP`, `NEXTSTEP`, `LEFTKP`, `TRANSFER`, `MAKEPS` |
| `stpline.f` | `STPLINE`, `FIRSTLSTP`, `NEXTLSTP`, `LEFTKPL`, `TRANSFL`, `PSFORLINE`, `TAUCENTRAL` |
| `matrix.f` | `LEVMIC`, `ULEVMIC`, `DLEVMIC`, `LAMBDA`, `DLAMBDA`, `REMOVELEVELS`, `LINFROMLEV`, `KAPSRC`, `EFFKAP`, `ATAU`, `HIIKAP`, `HIISMTH`, `MATRIXUP`, `LUSOLVER`, `LUDCMP`, `LUBKSB`, `MPROVE` |
| `sobolev.f` | `SOBOLEV`, `ESCAPE`, `GAUSSOB`, `GAULEG`, `SOBMATRIX`, `SOBKAPSRC`, `SOBEFFKAP` |
| `einstein.f` | `AIJ`, `BIJ`, `BJI`, `CFREQ`, `PROBTHERM`, `VTHERM`, `PLANCK`, `UBACK`, `HIIINIT` |
| `collrate.f` | `READSPECIES`, `CIJ`, `ILOCAT` |
| `ltrinteg.f` | `INTBEAM`, `SPECERR`, `UNYINT`, `QUINTGAUSS`, `UZINT`, `DINTCUB`, `DINTLIN`, `SPLINE` |
| `fitsadd.f` | `PRINTERROR`, `DELETEFILE` |

## Auxiliary programs

Several small auxiliary programs, primarily for testing, accompany the SimLine
code. `TESTCLOUD` (in `analyze/`) reads the level populations and energy densities
from a data file produced by the first part of SimLine and computes excitation
temperatures, radiation temperatures and a few additional quantities. The input
file must have the structure produced by `WRITEDENS`; the output file is
self-explanatory.

`CLINE` (in `src/`) reduces the two-dimensional data set produced by `WRITELINE`
and extracts the profile for zero beam offset, yielding a file containing only
one line profile. This is useful for users without a sophisticated graphics
front end such as GREG or IDL.

## Acknowledgements

The local statistical treatment of clumping and turbulence roughly follows the
formalism derived by H.M. Martin, D.B. Sanders & R.E. Hills (1985). The
accelerated lambda method was published by Lawrence Auer in 1987. The LU
decomposition method, the Gaussian quadrature, and the spline interpolation
algorithm were taken from "Numerical Recipes" (Press, Flannery, Teukolsky &
Vetterling, 1986).
