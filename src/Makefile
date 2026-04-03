# Build SimLine

all:
	@echo ""
	@echo "Please run make with one of the following targets:"
	@echo "  linux    - default Linux build (gfortran)"
	@echo "  mclarge  - large memory model (molecules with >100 levels)"
	@echo "  debug    - debug build with runtime checks"
	@echo "  osf      - DEC Alpha / OSF1"
	@echo "  aix      - IBM RISC / AIX"
	@echo "  hp-ux    - HP with HP-UX"
	@echo "  absoft   - Absoft f77"
	@echo "  sun-os   - SUN Sparc / SunOS or Solaris"
	@echo "  unknown  - alias for osf"
	@echo ""

unknown: osf

solaris: sun-os

# Default Linux build
linux:
	@$(MAKE) simline F77="gfortran" FFLAGS="-O" LIBS="-lcfitsio -lm -lz" LDFLAGS="-Wl,--no-relax"

# Large memory model for molecules with more than 100 levels
mclarge:
	@ln -sf fsizes_large.inc fsizes.inc
	@$(MAKE) simline F77="gfortran" FFLAGS="-O -mcmodel=large" LIBS="-lcfitsio -lm -lz" LDFLAGS="-Wl,--no-relax"

# Debug build with runtime checks
debug:
	@$(MAKE) simline F77="gfortran" FFLAGS="-mcmodel=medium -g -fcheck=all -pedantic -ffpe-trap=overflow" LIBS="-lcfitsio -lm -lz" LDFLAGS="-Wl,--no-relax"

# Optional: floating-point exception trap library (Linux/glibc only).
# Causes the program to abort with a traceback on invalid operations,
# divide-by-zero, or overflow — useful alongside the debug build.
# Build with:  gcc -shared -fPIC -o libtrapfpe.so trapfpe.c
# Then run as: LD_PRELOAD=./libtrapfpe.so ./simline

absoft:
	@$(MAKE) simline F77="af77" FFLAGS="-O -U -N109 -N15" LIBS="-lcfitsio -lm -lz -lU77"

osf:
	@$(MAKE) simline F77="f77" FFLAGS="-O -u -fpe1" LIBS="-lcfitsio -lm -lz"

aix:
	@$(MAKE) simline F77="f77" FFLAGS="-O -u" LIBS="-lcfitsio -lm"

hp-ux:
	@$(MAKE) simline F77="f77" FFLAGS="-O +U77 -u" LIBS="-lcfitsio -lm -lz"

sun-os:
	@$(MAKE) simline F77="f77" FFLAGS="-O -u -Nl32" LIBS="-lcfitsio -lm -lz"

windows:
	@$(MAKE) simline F77="f77" FFLAGS="-O" LIBS="cfitsio.lib"

SRCS = ltr.f ltrio.f initial.f sobolev.f stpdiff.f stpline.f matrix.f \
       collrate.f einstein.f ltrinteg.f fitsadd.f

simline: $(SRCS) fsizes.inc
	$(F77) $(FFLAGS) $(SRCS) $(LIBS) $(LDFLAGS) -o simline

fsizes.inc:
	ln -s fsizes_normal.inc fsizes.inc

clean:
	rm -f simline *.o fsizes.inc
