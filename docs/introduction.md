# Introduction

SimLine is a FORTRAN code to compute the profiles of molecular transition lines
in spherically symmetric clouds with arbitrary density, temperature and velocity
structure. Turbulence and clumping effects are treated in a local statistical
approximation combined with a radial dependence of the correlation parameters.

The code consists of two parts: the self-consistent solution of the balance
equations for all level populations and energy densities at all radial points,
and the computation of the emergent line profiles observed from a telescope with
finite beam width and arbitrary offset. All numerical discretizations are done
in an adaptive way, the user has full error control and the line profiles may be
computed in almost arbitrary accuracy. An accelerated lambda iteration is used
to solve the equation system. The optical depths in the lines may vary between
about -5 and +5000. The input of physical cloud parameters, telescope
parameters, and numerical control parameters may be done either interactively
or as file input. Molecular data are read from table files in LAMDA or a similar
compressed format.
