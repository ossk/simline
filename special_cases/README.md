# Special cases

This directory contains modified versions of selected SimLine source files for
specific physical regimes or testing purposes. Each file is a drop-in replacement
for its counterpart in `src/` and must be compiled in place of it (see the build
commands at the top of `simmaser.f77`).

| File | Replaces | Description |
|------|----------|-------------|
| `stpmaser.f` | `src/stpdiff.f` | _[one-sentence description]_ |
| `stpline_exact.f` | `src/stpline.f` | _[one-sentence description]_ |
| `stpdiff_exact.f` | `src/stpdiff.f` | _[one-sentence description]_ |
| `ltrinteg_simple.f` | `src/ltrinteg.f` | _[one-sentence description]_ |
| `simmaser.f77` | — | Build script and entry point for the maser variant of SimLine (uses `stpmaser.f`); contains compiler command lines for g77 and Absoft f77. |
