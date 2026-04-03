# Limitations and known issues

One should keep in mind that a main limitation in all radiative transfer
computations comes from the accuracy of available collision rates and the
temperature range where they are tabulated. At temperatures outside the interval
given in the molecular data file, the collision coefficients will be linearly
extrapolated. Moreover, one should also consider that the maximum number of
levels is restricted in the code, so that high temperatures combined with high
densities might produce non-vanishing excitation at levels that are neglected in
the computation. Users should verify that the maximum level number is sufficient
for their problem.

The program is not designed for strongly varying turbulence. The gradient in
both the absolute value of the turbulent velocity dispersion and in the
turbulence correlation length should not exceed the steepest gradient of either
the systematic velocity or the density.

Due to the limited number of frequency points, systematic velocities may be at
most one order of magnitude above the random velocities.
