# Installation

## System requirements

The program is written very close to standard FORTRAN 77 so that most f77
compilers should be able to compile it. It has been tested on DEC Alpha running
OSF/1, DEC VAX running Ultrix, SUN Sparc running SunOS, IBM RISC running AIX,
HP with HP-UX, and Intel machines running Linux (with f77/f2c, g77, and
Absoft f77). On DEC VAX machines running VMS, four lines in `ltrio.f` have to
be changed due to a different record treatment compared to Unix systems.

The system requires the [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/) library
to be installed so that the program can be linked against it. No other special
libraries are needed.

## Compilation

Compilation is done through `make` by specifying the target platform, e.g.:

```sh
cd src
make linux
```

Other available targets: `aix`, `osf`, `hp-ux`, `absoft`, `sun-os`, `unknown`.

For Linux there is a second option: `make mclarge`, which uses the large memory
model needed for molecules with more than 100 levels. In this case the field
sizes from `fsizes_large.inc` are used instead of `fsizes_normal.inc`. The
large memory target may run somewhat slower due to the increased memory
management overhead. When switching between targets, always run `make clean`
first.

## Memory usage

Memory consumption is dominated by the large field of energy densities at all
radial points. It consists of $2 \times n_s \times n_{\rm stretch} \times n_{\rm tra}$
`real*4` numbers, which amounts to approximately 1.9 GB in the normal version
and 4.8 GB in the large memory version. Adding roughly 5% for all other fields
and the dynamically linked libraries gives a reasonable estimate of total memory
consumption.

To achieve accuracies better than $2 \times 10^{-7}$, the program must be
compiled with `real*8` instead of the default `real*4`.
