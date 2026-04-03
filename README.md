# SimLine

SimLine is a FORTRAN code for computing molecular transition line profiles in
spherically symmetric clouds with arbitrary density, temperature and velocity
structure. Turbulence and clumping are treated in a local statistical
approximation with a radial dependence of the correlation parameters.

The code solves two coupled problems:

1. Self-consistent level populations via accelerated lambda iteration
2. Emergent line profiles convolved with a finite Gaussian telescope beam

All numerical grids are fully adaptive; the user has direct control over all
accuracy thresholds. Optical depths at line centre may range from approximately
−5 (moderate masering) to +5000.

## Features

- Spherically symmetric radiative transfer in an arbitrary number of molecular transitions
- Adaptive radial, frequency and displacement grids
- Accelerated lambda iteration for fast convergence
- Local statistical treatment of turbulence and clumping (effective absorption coefficients)
- Optional central H II region (free-free continuum source)
- Sobolev approximation available as initial guess
- FITS and ASCII output formats
- Reads molecular data in native format or LAMDA/RADEX format
- Interactive or fully file-driven operation

## Building

The code requires the [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/) library.

```sh
cd src
make linux        # default Linux build (gfortran)
make mclarge      # large memory model for molecules with >100 levels
make debug        # debug build with runtime checks
```

Other targets: `aix`, `osf`, `hp-ux`, `absoft`, `sun-os`. Always run
`make clean` before switching targets.

## Quick start

```sh
cd src
make linux
./simline
```

The program runs interactively by default. Pass `-fileinput` to read all
parameters from files, and `-numdefaults` to skip numerical parameter prompts.
Example input files are in `examples/`.

## Documentation

| Document | Description |
|----------|-------------|
| [Introduction](docs/introduction.md) | Overview of the code |
| [Installation](docs/installation.md) | System requirements and compilation |
| [Usage](docs/usage.md) | Command-line options, input parameters, file formats, error messages |
| [Theory](docs/theory.md) | Radiative transfer equations, turbulence model, Sobolev approximation |
| [Implementation](docs/implementation.md) | Source code structure and routine descriptions |
| [Limitations](docs/limitations.md) | Known limitations and caveats |
| [Future plans](docs/plans.md) | Planned features |
| [References](docs/references.md) | Literature references and acknowledgements |
| [Changelog](docs/CHANGELOG.md) | Version history |

The original LaTeX manual is preserved as `docs/simline.tex` (and compiled PDF
`docs/simline.pdf`).

## Repository layout

```
simline/
├── README.md
├── docs/               # Documentation
├── src/                # FORTRAN source code and Makefile
│   ├── ltr.f           # Main program and core routines
│   ├── ltrio.f         # I/O routines
│   ├── initial.f       # Grid initialization and adjustment
│   ├── stpdiff.f       # Radiative transfer (all lines)
│   ├── stpline.f       # Radiative transfer (single line)
│   ├── matrix.f        # Balance equations and lambda iteration
│   ├── sobolev.f       # Sobolev approximation
│   ├── einstein.f      # Molecular constants
│   ├── collrate.f      # Collision rate input
│   ├── ltrinteg.f      # Intensity integration routines
│   ├── fitsadd.f       # FITS I/O helpers
│   ├── cline.f         # Auxiliary: extract zero-offset line profile
│   ├── trapfpe.c       # Optional: floating-point exception trap (Linux/glibc, debug use)
│   ├── fsizes_normal.inc
│   └── fsizes_large.inc
├── analyze/
│   └── testcloud.f     # Auxiliary: analyse level population files
├── special_cases/      # Special-case reference programs
└── examples/           # Example input files
```

## Author

Dr. V. Ossenkopf-Okada  
I. Physikalisches Institut, Universität zu Köln  
ossk@ph1.uni-koeln.de

## License

SimLine is released under the [GNU General Public License v3.0](LICENSE)
(GPL-3.0-or-later).

Copyright © Dr. V. Ossenkopf-Okada, I. Physikalisches Institut, Universität zu Köln.
