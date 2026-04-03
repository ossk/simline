# Examples

This directory contains example input files for the numerical parameters of
SimLine. They can be passed to the program with the `-fileinput` option.

| File | Description |
|------|-------------|
| `default.inp` | Default numerical parameters as built into the code |
| `noit.inp` | Numerical parameters with initial populations set to `-3` (no iteration — useful for quick tests using precomputed level populations) |

The recommended way to generate input files for the physical, observational and
molecular parameters is to run SimLine interactively and let it write the files
itself, which guarantees a consistent file structure.
