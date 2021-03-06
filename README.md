##ECP CLASH Milestone #1 - Q1, FY 2017

This repository contains Boxlib-Codes and FLASH versions of an OpenACC accelerated Helmholtz equation of state. 
The accelerated equation of state routine itself is found in: 

`Microphysics/EOS/helmholtz/actual_eos.F90`

and 

`Helmholtz_starkiller/SpeciesBased/actual_eos_module.f90`

Copying these directories into a current CASTRO or MAESTRO (Microphyics) or FLASH (Helmholtz_starkiller) tree
will enable the use of the accelerated equation of state in either of these codes. Boxlib code users have the usual
freedom to place the Microphysics directory anywhere, as long as the build system is pointed to the correct path. 
FLASH users will want to copy Helmholtz_starkiller/ into:

`source/physics/Eos/EosMain` 

(i.e. as a peer to `source/physics/Eos/EosMain/Helmholtz`) 

Though the routines used in each case are essentially identical, there are some small differences in implementation in the current
versions of CASTRO and FLASH that lead to small differences. Primary among these are:

1) The Microphysics (Boxlib) version is set up to use the version of helm_table.dat that is packaged for Boxlib use, while the FLASH version
uses the helm_table.dat distributed with the public version of FLASH. These tables, though similar, have slightly different extents. (imax and jmax
control these extents.)

2) The FLASH version of the code has split up the actual_eos code and the actual_eos_init calls into separate files, in keeping with FLASH conventions.

3) CASTRO(Boxlib) and FLASH use different infrastructures to define physical constants and units. Though many constants used in the EoS routines
are declared locally, they are not necessarily consistent with values used in a build of the rest of each code. 

