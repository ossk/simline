# User interface

## Program invocation

The radiative transfer code is started from the command line by:

```
simline [options]
```

where options can be:

| Option | Description |
|--------|-------------|
| `-fileinput` | Read all physical, observational and numerical parameters from files; disable interactive input |
| `-numdefaults` | Use the built-in default numerical parameters without prompting |
| `-populateonly` | Skip ray tracing; compute only the level populations |
| `-stoponerror` | Terminate on errors instead of prompting for corrections (for non-interactive use) |
| `-overwrite` | Overwrite existing output files without confirmation |
| `-autonames` | Automatically generate output file names from the input model name and transition |
| `-fits` | Write line profiles in FITS format |
| `-block` | Write line profiles in ASCII block format |
| `-lines` | Write line profiles in ASCII line-oriented format |
| `-nocomments` | Suppress comment lines in ASCII output |
| `-unixcomments` | Use `#` as comment character in ASCII output |
| `-vmscomments` | Use `!` as comment character in ASCII output |
| `-writetau` | Write additional optical depth profile files alongside line profile files |

All options are intended to suppress interactive questions during program
execution. If omitted, the corresponding parameters are determined interactively,
except for `-writetau`, which must be specified on the command line if optical
depth profiles are needed.

All output format options (`-fits`, `-block`, `-lines`, `-nocomments`,
`-unixcomments`, `-vmscomments`, `-overwrite`) can also be set at compile time
via preprocessor defines (e.g., `-Dnocomments -Doverwrite`).

## Terminal interaction

The program is designed for interactive computation of line profiles, with the
flow controlled by terminal input. However, all parameters necessary for the
computation may be given either interactively or read from data files, except
for the molecular constants which must be provided in a file.

There is a strict separation between parameters describing the molecule,
parameters for the cloud geometry and physics, parameters of the observation,
and numerical parameters controlling the accuracy of the different parts of the
computation.

Parameters entered interactively may be stored immediately in a file with the
appropriate structure. This file can be read at the next run, and individual
values can be modified with a text editor. If an error is made during interactive
input, the input procedure can be restarted by entering a non-numeric character
at the next prompt.

The two parts of the program may be run in arbitrary order and repeatedly. One
can first compute a new cloud model (the level populations at all radial points
throughout the cloud), or read in the results of a previous computation. Then
the observed line profiles can be computed for an arbitrary number of lines or
telescope parameters. For manual line fitting one can change cloud parameters
and repeat the whole procedure.

Using a graphics program in a parallel session/window to visualize the line
profiles immediately is recommended.

## Data files

### Input files

All input files are plain ASCII files that can be modified with any text editor.
Output is either plain ASCII or FITS to enable easy examination, modification,
and use in further reduction pipelines.

The input files for molecular, physical, observational and numerical parameters
consist of pairs of lines. The first line of each pair is treated as a comment
(typically describing the parameter); the second line contains the parameter
value in ASCII form. Exceptions are: the molecular data file (where levels,
transitions and collisional rates are given in tables with only one comment line
per table) and the cloud parameter file when a data table is used (which also has
only one comment line per table).

It is strongly recommended to generate a first set of input files using SimLine
itself through interactive input, which guarantees consistency of the file
structure. Subsequently, individual values can be modified with a text editor.

### Molecular data files

Molecular data are read from ASCII files with a simple table structure. Data
files for several molecules accompany the SimLine executables and serve as
examples.

The file starts with a comment line describing the molecule and ends with an
arbitrary number of comment lines for references. Entries within the file are
given as single values or tables, each preceded by a comment line explaining it.
The following entries occur:

**Relative molecule mass:**
Relative mass of the isotope given in atomic mass units.

**Number of levels:**
At most 50 levels can be read.

**Number of transitions:**
At most 200 radiative transitions can be read.

**Levels: index, statistical weight, energy [cm⁻¹], name**
A table with one line per level. Index numbering starts at 1. Each level is
denoted by a string of up to 10 characters.

**Transitions: from index, to index, Einstein-A [s⁻¹]**
A table with one line per radiative transition. Indices must correspond to those
in the level table.

**Collision rates / Number of temperatures:**
Two comment lines separate the collision rates from the radiation data. Collision
rate coefficients can be tabulated for at most 10 different temperatures.

**Table of temperature values:**
One temperature per line.

**Number of collisional transitions:**
At most 1000 collisional transitions can be read.

**Transitions from – to index, rate coefficients:**
A table with one line per transition, containing the indices of the two levels
and the rate coefficients for all tabulated temperatures.

**Comments and references:**
An arbitrary number of comment lines at the end of the file, providing
references to the source data.

As an alternative, SimLine can read LAMDA files in RADEX format. However, it is
currently limited to one collision partner; LAMDA files with more than one
partner are rejected to avoid ambiguities.

### Output files

**Level population files:**
The first output file type contains the level populations resulting from the
lambda iteration. It is intended only for further processing by SimLine. There
are no comment lines. The first line contains the molecular data file name. The
second line contains the number of radial points, the number of rotational
levels treated, and the background radiation temperature. The following block
stores the physical parameters at each radial point: the point number, its
radius, molecule density, systematic velocity, thermal+turbulent velocity
$\sigma$, and the correlation length normalized by $\sigma$. A second block
contains the level populations: for each radial point (identified by its
number) the populations at all computed levels are listed. The same file
structure must be used when reading a precomputed cloud model to calculate new
line profiles.

**Line profile files:**
The second file type contains the line profiles, designed to be readable by
programs such as IDL, GREG or GNUPLOT. Two general formats are available:

- **FITS:** Lines are stored in a 2D array with the first axis being the angular
  offset and the second axis the velocity. Beam offsets are in degrees.
- **ASCII block mode:** The 2D field of beam temperatures is stored as a block
  with rows representing different beam offsets and columns representing
  frequency/velocity points. Two preceding comment lines give the offset and
  frequency values. Beam offsets are in arcseconds.
- **ASCII line mode:** Each line contains three numbers — beam offset, frequency,
  and beam temperature. The file contains (number of offsets) × (number of
  frequency points) lines, with preceding comment lines specifying these numbers.
  This format requires more disk space but can be read by more programs. Beam
  offsets are in arcseconds.

The comment line format in ASCII files can be set to VMS style (`!`) or UNIX
style (`#`), or comment lines can be suppressed entirely.

**Optical depth files:**
When `-writetau` is specified, an optical depth profile file is written alongside
every line profile file. These have the same format as the line profile files.

Units in all output files: radii in pc, velocities/frequencies in km/s. Blue-
shifted frequencies (above line centre) are counted as negative velocities.

## Input parameters

### Physical cloud parameters

**Background radiation temperature [K]:**
For normal Galactic sources, the cosmic microwave background (2.73 K) dominates
the external radiation at the wavelengths of rotational transitions. For
early-universe objects, higher values should be used.

**Number of shells for the cloud representation:**
The cloud may be subdivided into shells with different laws for the physical
parameters. These input shells have no relation to the radial grid used in the
radiative transfer computation. Discontinuities at shell edges can be specified
for power law regions and table input (experienced users only). The maximum
number of shells is about 180.

**Shell type [0=power law / 1=data table]:**
Within each shell, the radial dependence of the physical parameters can be
specified either by power law exponents or as a data table (with interpolation).

**Radius of the inner core [pc]:**
To avoid the singularity of power laws at zero, an inner core must be defined.
Within the core all parameters are assumed constant. The core may contain an HII
region. For data table input, the core radius equals the innermost table radius.

**Electron density [e/cm³] and electron temperature [K]:**
Required for a central HII region in the core. Set electron density to 0 to
assume no central HII region.

For each shell, the following parameters must be specified:

**Outer radius of the shell [pc]**

**Hydrogen density [H₂/cm³]**

**Temperature [K]:**
The kinetic temperature determines the collision rates and contributes to the
local line width.

**Relative molecular abundance [X/H₂]:**
Possible molecular depletions should be taken into account here.

**FWHM of turbulent velocity [km/s]:**
Thermal velocities are added automatically to obtain the total kinetic line
width.

**Turbulence correlation length [pc]:**
The scale at which the autocorrelation function of the turbulence drops to 1/e.

**Systematic radial velocity [km/s]:**
Infall is counted with negative numbers.

For power law input, all parameters except the radius are characterized by the
value at the inner shell boundary and a radial exponent. For data table input,
only the values at the corresponding radial points are required.

### Observational parameters

**Upper level of the transition to be observed:**
The index of the upper level in the molecular data file. By default, only one
line is computed at a time.

**Lower level of the transition to be observed:**
The index of the lower level in the molecular data file. If both upper and lower
levels are set to zero, SimLine computes line maps for all excited transitions of
the molecule using the same observational parameters.

**Frequency resolution [km/s]:**
May be limited by the detector or the observational goal.

**Beam width [FWHM in arcsec]:**
The beam is assumed to be Gaussian.

**Step size for mapping [arcsec]:**
Typically chosen to be about half the beam width.

**Central offset for the first map point [arcsec]:**
For observational offsets or composite maps, the radial map can start at a
finite offset from the cloud centre.

**Number of map points [0 = map whole cloud]:**
1 means a single observation (no map). 0 causes the program to automatically
compute the number of points out to the outer cloud radius.

**Distance of the cloud [pc]:**
Used to convert the physical cloud size into angular units.

### Numerical parameters

The same set of numerical parameters controls the accuracy for both the level
population computation and the emergent line profile computation. In a few cases,
a parameter is interpreted slightly differently in the two stages.

**Cut-off for Gaussians in the line emission [sigma], default: 3.2:**
Determines at which argument a Gaussian is treated as zero. The absorptivity and
emissivity are set to zero for frequencies beyond this cut-off, leaving the
background field unchanged.

**Cut-off for Gaussians in the frequency integration [sigma], default: 4.2:**
Sets the limits for integration of intensity convolved with a Gaussian profile.
Applied in the frequency integration for the local energy density and in the
spatial beam integration.

**Resolution in scanning the Gaussians [sigma], default: 0.55:**
Defines the resolution of the frequency grid used in the radiative energy density
integration. The grid is locally adapted to be at least as fine as this factor
times the local $\sigma$.

**Maximum change factor for level densities between grid points, default: 1.3:**
Level populations $n_j$ are linearly interpolated between radial grid points.
This parameter limits the maximum fractional change between adjacent points;
additional grid points are inserted when this limit is exceeded. The factor must
be greater than 1.

Reducing this parameter often provides the greatest accuracy improvement, at
the cost of substantially increased computation time. It is therefore the most
important numerical parameter to tune.

**Hysteresis in the radial grid adjustment [0..1], default: 0.2:**
The radial grid is dynamically adjusted to the gradient of the level populations.
The hysteresis sets the relative difference between the thresholds for adding and
removing grid points, suppressing oscillations during iteration.

**Maximum change of radial velocity between grid points [sigma], default: 2.5:**
As with density, the velocity structure may require additional radial grid
points. This parameter limits the change in systematic velocity between two
radial grid points.

**Maximum relative distance of points on the spatial integration scale, default: 0.2:**
Constrains the maximum relative distance between two points on the $z$ axis
in the energy density integration. Many points along the $z$ axis are needed
for accurate integration.

**Maximum velocity shift between radiative transfer points [sigma], default: 0.55:**
The emissivity and absorptivity are interpolated linearly for each frequency
within each step. For large velocity gradients, this requires additional steps.
The parameter gives the maximum step length in units of the $z$-direction velocity
$v_z$. Setting this equal to the frequency resolution is usually a good choice.

**Initial populations [1=thermal, 2=no radiation, 3=Sobolev, 4=file input], default: 2:**
Determines the first guess for level populations before the lambda iteration.
Thermal distribution is appropriate for optically thick clouds; radiation-free
for optically thin clouds; Sobolev approximation for large velocity gradients.
Due to the high efficiency of the accelerated lambda iteration, the choice of
initial guess has relatively little influence on the total number of iterations.
A level population file from a previous run can also be used, which is especially
useful when scanning parameter space for line fitting.

**Neglection threshold for the level populations, default: $10^{-8}$:**
Sets the truncation of the equation system. When the excitation of a level and
all higher levels falls below this threshold throughout the entire cloud, the
level is dropped from the computation. Values below $10^{-14}$ are below machine
precision.

**Relative accuracy as convergence criterion, default: $10^{-6}$:**
Convergence is tested after each full lambda acceleration cycle. When the
acceleration step is smaller than this value, convergence is declared. To achieve
accuracies better than $2 \times 10^{-7}$, the program must be compiled with
`real*8` instead of `real*4`.

Due to the strong nonlinearity of the line radiative transfer problem, no general
error estimate can be given for a specific set of numerical parameters. Tests
over a broad range of cloud parameters show that with the default numerical
parameters, the error in line intensities is always below 5–6%.

When `-numdefaults` is specified on the command line, all computations use the
built-in default numerical parameters without prompting. This avoids confirmation
requests during automated fitting, but makes it impossible to change numerical
parameters within the same program run.

## Error messages and warnings

Due to the static field declarations in FORTRAN 77, the number of grid points
on each scale is bounded by the compiled field sizes. In the normal version:
180 radial points, 270 impact parameters, 223 frequency points, 20 molecular
levels, and 19 radiative transitions. This is suitable for most linear molecules.
The large memory version (`make mclarge`) provides larger field sizes for
molecules with up to 100 levels.

When a required field size exceeds the compilation values, the program either
returns to the main input loop or stops with exit code 1 (with `-stoponerror`).
Possible messages are:

**`Number of radial points insufficient to treat the problem with the required accuracy!`**
The gradient in the level densities requires more than 140 radial points, due to
a very steep density gradient or strong changes in the level populations (e.g., at
the edge of very thick clouds). Either increase the "maximum change factor for
level densities" parameter, or replace steep density gradients with discontinuities
(possible only with file input).

**`Number of ray points insufficient to treat the problem with the required accuracy!`**
This rarely occurs when the $z$ integration grid spacing is very small combined
with many radial grid points. It can also appear when computing emergent line
profiles with very large systematic velocities compared to local random
velocities, or when mapping a large cloud with a very fine beam. Increase the
"maximum relative distance of points on the spatial integration scale" parameter,
or for large velocities, increase the "resolution in scanning the Gaussians"
parameter.

**`Number of frequency points insufficient to treat the problem with the required accuracy!`**
More than 223 frequency points would be required given the maximum cloud
velocities and the required frequency resolution. In the level population
computation (first stage): increase the Gaussian scanning resolution or reduce
the frequency integration cut-off. In the emergent profile computation (second
stage): reduce the observational frequency resolution or reduce the frequency
integration cut-off. There is no solution if systematic velocities greatly
exceed the local line width.

**`Numerical singularity occurred for the given parameters!`**
The general worst-case check failed; the code can no longer compute non-singular
level populations. This typically results from unphysical input parameters such
as negative temperatures, temperatures outside the range of the collision data,
or negative radii. Carefully inspect all input parameters.

**`I cannot open the input file`**
The specified input file does not exist or is not readable. With `-stoponerror`,
exits with code 2.

**`Format error in input file`**
The line structure of the input file is unrecognizable: parameters are in the
wrong order or comment lines are missing. With `-stoponerror`, exits with code 3.

**`I cannot open the output file`**
The output file exists and `-overwrite` is not set, or permissions do not allow
writing, or the specified directory does not exist. With `-stoponerror`, exits
with code 4.

If the molecular data file contains more levels or transitions than the code can
handle, the molecule system is silently truncated at high energies until it fits
the internal fields. There is no explicit check whether the level count is
sufficient for a given problem; users should inspect the level density output
file and check the excitation of the last level. This can be a limiting factor
for warm, optically thick clouds.

Consistency check messages that may appear during input:

**`Line overlap detected: transition - from - to - at - Hz`**
The program detects overlapping lines but does not yet treat them correctly in
radiative transfer. Results will not be reliable for the affected lines.

**`There is no excited transition from level - to level -!`**
The selected transition does not exist as a radiative transition in the molecular
data file, or includes a level whose excitation fell below the neglection
threshold. Select a different transition.

**`Only finite beam widths are possible!`**
The beam FWHM must be greater than 0.

Note that most input parameters are not range-checked. Unphysical or absurd
parameters (e.g., cloud temperatures below the cosmic background temperature)
may lead to strange results or abnormal program termination.

## Hints for line fitting

- Start with the default numerical parameters. Only at a late stage should the
  accuracy be checked by tightening the numerical limits.

- When fitting multiple lines, always start with a simple one-shell model (unless
  a detailed physical model suggests otherwise). A single shell already allows a
  very wide variety of line shapes and intensity ratios. In complex models it is
  easy to get lost in parameter space. *Remember: in 1-D, there is never a unique
  solution to a line fitting problem!*

- To get a reasonable initial guess for the cloud parameters, use traditional
  methods such as comparing line strengths from different isotopologues of the
  considered molecule.

- When the cloud model spans a density range greater than $10^3$ or has a complex
  density structure, the radial grid may run out of points. For smooth density
  structures, relax the "maximum change factor for level densities" parameter.
  For strongly varying density profiles, it may be more appropriate to replace
  steep gradients with artificial density jumps (discontinuities can only be
  specified using file input).

- At the end of a density fit, check the consistency of the temperature estimate
  with the computed line temperature of the thick line used for temperature
  derivation. Even quite thick lines ($\tau > 10$) are often not completely
  thermalized. Using $^{12}$CO can provide a much more accurate temperature
  determination this way.

- As a first orientation for the turbulence parameters, Miesch & Bally (1994)
  found correlation lengths of about 0.1 pc for well-resolved clouds. Their
  values for the Kolmogorov exponent $\gamma$ (giving the radial dependence of
  the turbulent velocity dispersion) range from 0.74 for L1228 to around 0.5 for
  many other objects. For homogeneous isotropic turbulence, there is a theoretical
  anticorrelation between $\gamma$ and the correlation length:
  $\gamma \propto 1/\ln(r_{\rm c})$.

- The effective absorption coefficient approximation for a clumpy medium is valid
  for clumping in both real space and velocity space. The correlation length
  description corresponds to a homogeneous medium with pure velocity-space
  turbulence. For clumping in real space, there is no independent filling factor
  parameter; instead, use the clump hydrogen density for the density parameter,
  reduce the molecular abundance by the filling factor, and increase the
  correlation length by $1/f_{\rm fill}$.
