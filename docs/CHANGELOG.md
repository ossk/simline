# Changelog

## 2.17 → 2.18

- Allowed for longer level names (up to 10 characters) and molecules with up to 100 levels
- Fixed long filename handling for molecule file in level population file
- Fixed problem of artificial steps for table input leading to requirement of too many radial grid points
- Fixed possible NaN creation in rare cases of steep temperature gradients in the core
- Introduced large field version to be compiled with `make mclarge`
- Switched from included cfitsio library to linking against the system library

## 2.16 → 2.17

- Unified field size parameter handling for all fields into one `fsizes.inc`
- Increased length of input strings for file names to 250 characters
- Added `-writetau` option for additional writing of optical depth files

## 2.15 → 2.16

- Updated manual to contain all possible command line options
- Added rest frequency in all output formats with comments to allow easy conversion from temperature to energy scales

## 2.14 → 2.15

- Increased the field sizes in the routines from Numerical Recipes to allow for grids of a few thousand points defined in `fsizes.inc`

## 2.13 → 2.14

- Added exception handling forbidding electron temperatures below cosmic background in HII region

## 2.12 → 2.13

- Fixed condition for HII region that either electron density or temperature can be zero to specify a constant continuation

## 2.11 → 2.12

- Added exit codes to all stop commands for automatic fitting in MAGIX
- Truncation of Gaunt factor approximation for high frequencies to prevent negative optical depths
- Switched to composition of spline and trapezium integration for beam integration (splines created overshooting in an extreme case)
- Test for non-monotonous radial grid

## 2.10 → 2.11

- Adapted Makefile to new cfitsio library
- Switched from g77 to gfortran for Linux
- Changed indentation to be gfortran compatible

## 2.09 → 2.10

- Added exception handling for zero-width input shells
- Added reader for LAMDA files
- Added Windows libraries
- Upgraded Absoft compiler support to version 10.0

## 2.08 → 2.09

- Improved handling of grid of rays in radiative transfer for observed lines
- Small adjustment of radial grid treatment
- Reduced default value of epsz from 0.25 to 0.2
- Added error detection for temperatures below cosmic background in input

## 2.07 → 2.08

- Removed inconsistent accuracy handling in Sobolev introduced in 2.07
- Added options `-populateonly` and `-stoponerror`
- Fixed a bug in the initial grid adjustment that could lead to a numerical singularity in case of strong negative systematic velocity gradients
- Added exception handling for floating point overflows in irrelevant collision rates

## 2.06 → 2.07

- Added new atoms/molecules: OI, CII, H2CO, CH3OH
- Updated constants for CI
- Fixed bug in string handling of level names including spaces
- Introduced exception handling for levels with the same energy
- Restricted grid adjustment to significantly excited levels
- Increased accuracy in maser handling that could prevent convergence

## 2.05 → 2.06

- Fixed wrong velocity order in FITS files
- Fixed a bug in error handling when rewriting FITS files (some compilers)
- Reduced the range for the linear interpolation of the exponentials
- Added several molecules

## 2.02 → 2.05

- Added option `-autonames` for automatic generation of output file names based on input model name
- Optical depth output added as a comment in the line files
- Transition from 0 to 0 generates a loop over all excited lines
- Fixed small flaw in the call for the lambda iteration for level populations

## 2.01 → 2.02

- Improved accuracy of the small-tau approximation, extending its range — small speedup
- Fixed bug in interpolation of the source function that could result in strong localized errors in case of considerable masering
- Substituted linear extrapolation of collision rates at small temperatures by square root extrapolation
- Added more molecule files

## 2.0 → 2.01

- Fixed bug that produced a numeric overflow at very high effective opacities in the macroturbulent description
- Updated molecule files from JPL database
- Clean exception handling for zero collision rates introduced

## 1.51 → 2.0

- Name changed from LTR to SimLine since the original abbreviation became meaningless
- Molecular data moved out of the code into external files for greater flexibility
- Extension to arbitrary nonlinear molecules; changes in observational parameter input file format
- Introduction of data table treatment for physical cloud input parameters; changes in that input file format
- Introduction of FITS files as the default line output file format
- Input of precomputed level populations as initial guess
- All field sizes defined in a single include file
- Iterative improvement for large matrices; accelerated convergence for molecules with many transitions
- Introduction of a smart maser treatment for weakly saturated masers

## 1.5 → 1.51

- Removal of a format error writing a wrong comment line to some parameter files
- Introduction of the background temperature as a free parameter; this makes old physical parameter files and population files incompatible

## 1.4 → 1.5

- Introduction of new collision rates for HCO+ from T. Monteiro (1993) with up to 11 levels; all molecules can now be treated up to 300 K
- Removal of rarely used input parameters (lower and upper cut); cuts must now be introduced via additional shells
- Increase of the maximum shell number
- Acceleration of convergence in the Sobolev part
- Removal of two bugs that could prevent convergence in certain special cases on some machines
- Modularized error treatment

## 1.31 → 1.4

- Improvement of the accelerated lambda iteration; convergence speed-up for large optical depths, overshooting prevented
- Better extrapolation of CO–ortho-H₂ collision rates beyond 100 K, providing a smooth transition to the CO–para-H₂ values
- Addition of two new observational parameters for better map control (central offset and number of map points)
- Increase of all field sizes in frequency points

## 1.3 → 1.31

- Removal of a bug in the transfer on the central ray which wasted computing time
- Accuracy improvement for the routine computing the central optical depth (control output)

## 1.21 → 1.3

- Introduction of the treatment of a central HII region in the cloud core; new physical parameters for electron density and temperature within the HII region
- New molecule added (SiO in the ground vibrational state)
- New collision rates for CS from Turner et al. 1992, ApJ 399, 114; all molecules except HCO+ can now be treated at kinetic temperatures up to 300 K
- Improvement of ray spacing in the final computation of line profiles for large velocity gradients
- Removal of a bug in the spatial intensity integration appearing at density edges
- Further exception handling added for very high intensities (e.g., local masers)

## 1.2 → 1.21

- Removal of a severe bug in the central transfer code; addition of a second-order term important for large velocity gradient regimes

## 1.1 → 1.2

- Introduction of a first approximation for the treatment of turbulence and clumping in space or velocity space (effective optical depths) — Martin, Sanders & Hills 1984, MNRAS 208, 35
- Addition of the turbulence correlation length as an additional physical input parameter
- Possibility to eliminate most interactive questions during program execution via command line options or compiler directives

## 1.0 → 1.1

- Restructuring of field handling in the radiative transfer part; elimination of the emissivities → 20% speed-up
- Enabling compile-time and command line options giving default values for the output file structure and setting overwrite mode
- Removal of a bug in the Sobolev part that could slow down convergence in certain cases
