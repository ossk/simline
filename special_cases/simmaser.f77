#g77 -O -Wimplicit ltr.f ltrio.f initial.f sobolev.f stpmaser.f stpline.f matrix.f collrate.f einstein.f ltrinteg.f fitsadd.f libfitsio.a -o simmaser
af77 -O -U -B101 -N109 -N1 ltr.f ltrio.f initial.f sobolev.f stpmaser.f stpline.f matrix.f collrate.f einstein.f ltrinteg.f fitsadd.f libfitsio_af77.a -lU77 -o simmaser
