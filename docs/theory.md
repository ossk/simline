# The basic theory

## The solution of the radiative transfer problem

The general formulation of the frequency-dependent radiative transfer equation

$$\frac{dI(\nu,\vec{r},\vec{n})}{ds} = -\kappa(\nu,\vec{r},\vec{n}) I(\nu,\vec{r},\vec{n}) + \epsilon(\nu,\vec{r},\vec{n})$$

may be reduced in the spherically symmetric case to three independent
coordinates. It is sufficient to consider radiation in the direction of the
$z$ axis and to neglect the $\phi$ dependence of the radius vector. Then the
radiative transfer equation can be written as

$$\frac{dI_z(\nu,p,z)}{dz} = -\kappa(\nu,p,z) I_z(\nu,p,z) + \epsilon(\nu,p,z)$$

where $p$ denotes the displacement variable ($p=\sqrt{x^2+y^2}$). The $z$ axis
is chosen in the direction towards the observer. Since the spatial grid is
adjusted at distinct radial points, the program does not use the form of the
radiative transfer equation with source function and derivative on the optical
depth, but with emission and absorption coefficients.

For the numerical solution of the radiative transfer equation it is convenient
to use stepwise the solution of the integral equation, allowing reduction of the
error to third order:

$$I_z(\nu,p,z_i) = \exp \left( -\int_{z_{i-1}}^{z_i}
\kappa(\nu,p,z) dz \right) \left[ I_z(\nu,p,z_{i-1}) +
\int_{z_{i-1}}^{z_i} \epsilon(\nu,p,z) \exp \left( \int_{z_{i-1}}^z
\kappa(\nu,p,z') dz' \right) dz \right]$$

The incident radiation at the outer boundary of the cloud is assumed to follow
a black-body spectrum with temperature $T_{\rm bg}$:

$$I_{\rm bg}(\nu) = \frac{2 h \nu^3}{c^2}\left[\exp\left(\frac{h\nu}{k T_{\rm bg}}\right) - 1\right]^{-1}$$

For the cosmic background radiation, $T_{\rm bg} = 2.73 \mathrm{K}$.

The program assumes complete redistribution of energy between the molecules at
a given level, so that one can define a unique set of level populations $n_j$
at any point $\vec{r}$ which is independent of frequency or the actual velocity
of a molecule. This also means that the line profiles for emission and
absorption are the same. In molecular clouds it is sufficient to consider pure
Doppler broadening so that the profiles are given by

$$\Phi(\nu) = \frac{1}{\sqrt{\pi}\sigma} \exp\left(-\frac{(\nu - \nu_0(1+v_z/c))^2}{\sigma^2}\right)$$

where $v_z$ is the systematic velocity in the $z$ direction, $\sigma$ is the
local line width, and $\nu$ is the frequency in the observer's frame. (Without
loss of generality, the program assumes that the centre of the cloud is at rest
in the observer's frame.)

With the assumption of complete redistribution, the absorption and emission
coefficients for a given transition are determined by

$$\begin{aligned}
\kappa_{j,l}(\nu) &= \frac{h\nu_0}{c}(n_l B_{l,j} - n_j B_{j,l})\Phi(\nu) \\
\epsilon_{j,l}(\nu) &= \frac{h\nu_0}{4\pi} n_j A_{j,l}\Phi(\nu)
\end{aligned}$$

where stimulated emission is treated as negative absorption.

The level populations at a given radius $r$ may be computed from the linear
system of balance equations

$$n_j \sum_{l\ne j}\left( A_{j,l}+B_{j,l}u_{j,l}+c_{j,l}\right) =
\sum_{l \ne j} n_l \left(A_{l,j}+B_{l,j}u_{j,l}+c_{l,j}\right)$$

where the radiative energy density $u_{j,l}$ for the transition between levels
$j$ and $l$ at a given point is computed from the integral

$$u_{j,l}(p,z) = \frac{1}{c} \int_0^{2\pi}d\phi \int_0^\pi d\theta \sin\theta
\int_{-\infty}^{\infty} d\nu \Phi_{\phi,\theta}(\nu,\vec{v})
I(\nu,p,z,\phi,\theta)$$

Exploiting spherical symmetry, this may be transformed into

$$u_{j,l}(r) = \frac{2\pi}{rc} \int_{-r}^r dz' \int_{-\infty}^\infty
\Phi(\nu) I_z\left(\nu,p'=\sqrt{r^2-z'^2},z'\right) d\nu$$

where the radius $r$ is given as $r=\sqrt{p^2+z^2}$. Due to the symmetry, only
radiation parallel to the $z$ axis is encountered. For the central line
frequency $\nu_0$ within the expression for $\Phi$, the frequency of the
transition $j \rightarrow l$ has to be taken when computing $u_{j,l}$.

The equation system is truncated whenever the excitation of a level falls below
the chosen accuracy limit or at the maximum level for which the collision rates
are known.

## The local turbulence approximation

The turbulence description uses two parameters for each point $\vec{r}$: the
total width of the velocity distribution $\sigma$ giving the local profile for
optically thin lines, and the correlation length $r_{\rm corr}$.

The width of the velocity distribution $\sigma$ is composed of a turbulent and
a thermal contribution:

$$\sigma = \frac{\nu_0}{c}\sqrt{\frac{2kT_{\rm kin}}{m} + \frac{2}{3}\langle v_{\rm turb}^2 \rangle}$$

where a Maxwellian distribution of turbulent velocities is assumed. In the data
input, the turbulent velocities are given as FWHM of the velocity distribution
$(\mathrm{FWHM}(v_{\rm turb})=\sqrt{8/3\times\ln 2 \langle v_{\rm turb}^2\rangle})$.
The thermal velocities are automatically computed from the kinetic gas temperature.

The long-range correlation of the turbulence spectrum described by a Kolmogorov
exponent can be simulated in 1D by a radial dependence of the turbulent velocity
dispersion $\langle v_{\rm turb}^2\rangle \propto r^\gamma$. Observations of
Fuller & Myers (1992) suggest exponents around $\gamma \approx 0.7$; several
other authors proposed somewhat smaller values around $\gamma \approx 0.5$.

For the local treatment of clumping in real space or velocity space, the
considered volume element is subdivided into numerous clumps with a thermal
internal velocity dispersion. From a Gaussian distribution of the density of
molecules with a certain velocity within the single fragments,

$$n(r) = n_0 \times \exp(-r^2/r_{\rm cl}^2)$$

we get an effective absorption coefficient for the whole medium given by

$$\kappa_{\rm eff} = n_{\rm cl} \times \pi r_{\rm cl}^2
\int_0^{\tau_{\rm cl}} \frac{1-\exp(-\tau)}{\tau} d\tau$$

(Martin et al. 1984), where $n_{\rm cl}$ is the number density of clumps
contributing at the considered frequency and
$\tau_{\rm cl}=\sqrt{\pi} \kappa r_{\rm cl}$ is their central opacity.

Assuming a Maxwellian turbulent velocity distribution of clumps and a thermal
velocity dispersion $\sigma_{\rm th}$ giving at most 1/3 of the total velocity
dispersion $\sigma$, we obtain an effective absorption coefficient

$$\kappa_{\rm eff}(\nu) = n_{\rm ges} \pi r_{\rm cl}^2 \times A(\tau_{\rm cl})
\times \frac{\sigma_{\rm th}}{\sigma}
\exp\left(-\frac{(\nu-\nu_0)^2}{\sigma^2}\right)$$

with

$$A(\tau) = \frac{1}{\sqrt{\pi}}\int_{-\infty}^{\infty}dv
\int_0^{\tau\exp(-v^2)} \frac{1-\exp(-\tau')}{\tau'} d\tau'$$

Here, $n_{\rm ges}$ is the total number of clumps. In case of pure turbulence
in velocity space, $n_{\rm ges}$ is given by the volume of the single clumps
and we obtain

$$\kappa_{\rm eff}(\nu) = \kappa(\nu) \times A(\tau_{\rm cl})/\tau_{\rm cl}$$

Consequently, $A(\tau_{\rm cl})/\tau_{\rm cl}$ is a measure of the reduction of
the opacity due to turbulence. For small clump sizes,
$A(\tau_{\rm cl}) = \tau_{\rm cl}$ so that we are in the microturbulent limit.
For $\tau_{\rm cl}\gg 1$ the macroturbulent limit is reached. In case of real
clumping in space, $\kappa_{\rm eff}(\nu)$ is further reduced by the filling
factor. In the program, this is simulated by a corresponding artificial
reduction of the molecular abundance.

The clump size $r_{\rm cl}$ is determined by the correlation length of the
velocity or density structure: $r_{\rm cl}$ is the length on which the abundance
of molecules within the same thermal velocity profile is reduced by the factor
$1/e$. Consequently, it can be computed from the turbulence correlation length
by $r_{\rm cl} = r_{\rm corr} \times \sigma_{\rm th}/\sigma$.

One should keep in mind that this is a statistical approximation, always
assuming a complete Maxwellian velocity distribution within each volume element
even if the single clumps become relatively large. Furthermore, the source
function $S=\epsilon(\nu)/\kappa(\nu)$ is not influenced in this local approach.

## The central H II region

For the simulation of internal heating produced by a central continuum source in
the cloud, it is possible to assume an H II region in the cloud core. The H II
region is characterized by two parameters: the electron density $N_{\rm e}$ and
the kinetic electron temperature $T_{\rm e}$.

The absorption coefficient for electron-ion bremsstrahlung in the Rayleigh-Jeans
approximation is given by:

$$\kappa(\nu) = \frac{8}{3\sqrt{2\pi}} \frac{e^6}{(4\pi\epsilon_0 m_{\rm e})^3 c}
\left(\frac{N_{\rm e}}{\nu}\right)^2 \left(\frac{m_{\rm e}}{k T_{\rm e}}\right)^{3/2}
\ln\Lambda$$

where the gas is assumed to be singly ionized and $\Lambda$ is given by

$$\Lambda = \left(\frac{2k T_{\rm e}}{\delta m_{\rm e}}\right)^{3/2}
\frac{4\pi\epsilon_0 m}{\pi\delta e^2 \nu} \approx
4.9573 \times 10^7 \left(\frac{T}{\mathrm{K}}\right)^{3/2} \frac{\mathrm{Hz}}{\nu}$$

for $T_{\rm e} < 3.2 \times 10^5 \mathrm{K}$. The quantities $e$ and
$m_{\rm e}$ denote the electron charge and mass, and $c$ is the speed of light.

For a thermal plasma, the emission coefficient follows from the Planck function:

$$\epsilon(\nu) = \kappa(\nu) \times B_\nu(T_{\rm e})$$

In the radiative transfer computations within the H II region, the small
frequency dependence of these continuum coefficients within the molecular
rotational lines is neglected.

## The Sobolev approximation

As an initial guess for the level populations, the Sobolev approximation may be
used. This approximation is strictly valid for clouds with large velocity
gradients and velocities increasing outwards, so that all radiative coupling is
restricted to a small spatial region. The level populations can then be computed
from local quantities only.

The radiative energy density at radius $r$ is determined by the local source
function $S_{j,l}(r)=\epsilon_{j,l}(r)/\kappa_{j,l}(r)$:

$$u_{j,l}(r) = \left(1-\beta_{j,l}(r)\right)S_{j,l}(r) + \beta_{j,l}(r)I_{\rm bg}(\nu_0)$$

where $\beta_{j,l}(r)$ denotes the photon escape probability for line $j$
emitted at $r$. It can be computed from an integral over all spatial directions
from a given point:

$$\beta_{j,l}(r) = \frac{1}{2}\int_{-1}^1
\frac{1-\exp(-\kappa_{j,l}(r)/Q(r,\mu))}{\kappa_{j,l}(r)/Q(r,\mu)} d\mu$$

with the line-integrated effective absorption coefficient $\kappa_{j,l}$ and the
velocity gradient in spherical symmetry:

$$Q(r,\mu) = \mu^2\frac{dv_r}{dr} + (1-\mu^2)\frac{v_r}{r}$$

Hence, the radial velocity gradient is the most important quantity for the local
energy density.

In contrast to the iterative solution of the balance equations, a more efficient
approach is possible here, since the energy densities are analytically coupled to
the level populations via the absorption and emission coefficients and the Sobolev
equations above. Using the derivatives of the matrix coefficients, the balance
equations can be solved by a Newton-Raphson approach. The corrections to the level
populations within each step are computed from the matrix equation

$$\sum_i\left(\sum_k \frac{\partial A_{jk}}{\partial n_i}n_k + A_{ji}\right)
\Delta n_i = \sum_i A_{ji} n_i$$

where the $A_{ij}$ are the coefficients of the matrix of balance equations.

## Computation of beam temperatures

When the level populations are known, the beam temperature relative to the
background is computed from the emergent intensity by integration over the
projection of the telescope beam on the cloud:

$$T_{\rm beam}(\nu) = \frac{c^2}{2k\nu_0^2}
 \frac{\displaystyle\int_0^{2\pi} d\phi\int_0^\infty dp p
(I_{\rm surf}(\nu,p)-I_{\rm bg}(\nu)) f_{\rm beam}(p,\phi)}
{\displaystyle\int_0^{2\pi} d\phi\int_0^\infty dp p f_{\rm beam}(p,\phi)}$$

The emergent intensity is the value on the cloud surface
$I_{\rm surf}(\nu,p)=I_z(\nu,p,\sqrt{R_{\rm cloud}^2-p^2})$. A Gaussian
profile for the beam is assumed:

$$f_{\rm beam}(p,\phi) = \exp\left(\frac{-(p-p_{\rm offset})^2(1+\phi^2)^2}{\sigma_{\rm beam}^2}\right)$$

The projected beam width is computed from the angular width by
$\sigma_{\rm beam} = \pi D^2/648000 \sigma_{\rm beam}['']$ where $D$ is the
distance of the cloud. The standard deviation $\sigma_{\rm beam}$ is coupled to
the full width at half maximum by
$\mathrm{FWHM} = 2\sqrt{\ln 2} \sigma_{\rm beam}$. The program computes a full
map, i.e., the line profiles at a given number of positions on a linear radial
scan through the cloud.

Within the program, all frequencies are treated in units of $(\nu-\nu_0)/\nu_0$,
velocities in units of $c$, and intensities in units of $\nu_0^3/c^2$, i.e. by
a factor $4\pi/h$ larger than the SI unit.

