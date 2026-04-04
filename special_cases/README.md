# Special cases

This directory contains modified versions of selected SimLine source files for
specific physical regimes or testing purposes. Each file is a drop-in replacement
for its counterpart in `src/` and must be compiled in place of it (see the build
commands at the top of `simmaser.f77`).

| File | Replaces | Description |
|------|----------|-------------|
| `stpmaser.f` | `src/stpdiff.f` | Normal simline linearizes maser emission, this routine allows for the exponential emplification, potentially leading to numerical instabilities |
| `stpline_exact.f` | `src/stpline.f` | Normal simline uses a stepwise approximation for the source function in linear gradients, here the correct analytic description is used |
| `stpdiff_exact.f` | `src/stpdiff.f` |  Normal simline uses a stepwise approximation for the source function in linear gradients, here the correct analytic description is used |
| `ltrinteg_simple.f` | `src/ltrinteg.f` | Coarser line integration routine, probably good enough for most cases |
| `simmaser.f77` | — | Build script and entry point for the maser variant of SimLine (uses `stpmaser.f`); contains compiler command lines for g77 and Absoft f77. |
