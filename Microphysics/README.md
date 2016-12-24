# Microphysics

## N. B. This README is taken from the main BoxLib-Codes/Microphysics repository at

`https://github.com/BoxLib-Codes/Microphysics/`

*A collection of astrophysical microphysics routines with interfaces to
 the BoxLib codes*

To use this repository with BoxLib codes, set `MICROPHYSICS_HOME` to point
to the `Microphysics/` directory.

There are several core types of microphysics routines hosted here:

* `EOS/`: these are the equations of state.  All of them use a Fortran
  derived type `eos_t` to pass the thermodynamic state information in
  and out.

* `unit_test/`: code specific to unit tests within this repo.  In
  particular,

  - `test_eos` will test an equation of state by first calling
    it will (rho, T), and then calling it with other inputs
	to recover rho and/or T.  A cube of data, with rho, T, and
	X is tested.


These routines are written to be compatible with:

* Castro: http://boxlib-codes.github.io/Castro/

* Maestro: http://boxlib-codes.github.io/MAESTRO/


A user's guide for Microphysics can be found in `Docs/`.  Type `make`
to build it from its LaTeX source.

A PDF of the user's guide is available here:
http://bender.astro.sunysb.edu/Castro/staging/Microphysics/Docs/MicrophysicsUsersGuide.pdf
