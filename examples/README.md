# Examples

This directory contains example input files for the numerical parameters of SimLine. They can be passed to the program with the `-fileinput` option or when interactively asked for the input file.

| File | Description |
|------|-------------|
| `pdr100_40_00_40_00_C34S.model` | Physical model from tabulated data from a PDR model |
| `pdr100_40_00_40_00_C34S.obs` | Observational data for observing that model with a fine pencil beam, all transitions |
| `c34s_and_h2.mol` | Molecular parameter file for C34S, native structure |
| `MonR2_test_HCN.mode`l | Physical model for a region with bright central HII region and steep gradients, numerically challenging |
| `MonR2_test_HCN.obs` |  Observational data for observing the  MonR2_test_HCN.model in the lowest transition only |
| `hcn.lamda` | Molecular parameter file for HCN, lamda structure |
| `default.inp` | Default numerical parameters as built into the code, for inspection |
| `noit.inp` | Numerical parameters with negative number for initial populations - switches off the lamda iteration - useful for quick tests |

The recommended way to generate input files for the physical, observational and molecular parameters is to run SimLine interactively and let it write the files itself, which guarantees a consistent file structure.
