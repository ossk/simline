# Contributing

Contributions to SimLine are welcome. Please follow these guidelines:

- For **bug reports** and **feature requests**, open a GitHub issue describing
  the problem or suggestion in as much detail as possible, including the version
  of SimLine used, the operating system and compiler, and a minimal example that
  reproduces the issue.

- For **code contributions**, please contact the author
  (ossk@ph1.uni-koeln.de) before starting work on substantial changes, to
  avoid duplicated effort and ensure the changes fit the overall design. Small
  bug fixes may be submitted directly as pull requests.

- Please do **not** change the numerical algorithms or physical approximations
  without discussion — correctness of the radiative transfer is the primary
  design goal.

- If you add support for a new platform or compiler, please update the
  `src/Makefile` and `docs/installation.md` accordingly.

- New molecular data files contributed by users are very welcome.
