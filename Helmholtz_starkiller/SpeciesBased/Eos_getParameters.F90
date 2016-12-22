!!****if* source/physics/Eos/EosMain/Helmholtz/SpeciesBased/Eos_getParameters
!!
!! NAME
!!
!!  Eos_getParameters
!!
!! SYNOPSIS
!!
!!  call Eos_getParameters(OPTIONAL,real(OUT)    :: eintSwitch,
!!                         OPTIONAL,logical(OUT) :: inputsAreUnchanged,
!!                         OPTIONAL,logical(OUT) :: inputTempIsGuess,
!!                         OPTIONAL,logical(OUT) :: constantGammaC,
!!                         OPTIONAL,logical(OUT) :: inputMassFracNeeded,
!!                         OPTIONAL,real(OUT) :: smalle)
!!
!! DESCRIPTION
!!
!!  Return information about the Eos implementation that may be of interest
!!  to other units.
!!
!! ARGUMENTS
!!
!!   eintSwitch -          value of the Eos unit's eintSwitch runtime parameter.
!!   inputsAreUnchanged -  Indicates whether calls to Eos (or Eos_wrapped, etc.) can result
!!                         in modification of some state variables that should be, strictly
!!                         speaking, input only to the EOS in the given mode.
!!                         If this is true, then calls to Eos with MODE_DENS_PRES can modify
!!                         the pressure, and calls with MODE_DENS_EI can modify enery variables.
!!   inputTempIsGuess -    Indicates whether the Eos implementation uses the temperature
!!                         provided on entry to a call as in initial gues in an iterative scheme.
!!   constantGammaC -      Indicates whether the gamc returned by Eos will always be constant.
!!   inputMassFracNeeded - Indicates whether the Eos implementation makes use of mass fractions.
!!   smalle -              value of the Eos unit's smallE runtime parameter.
!!
!! EXAMPLE
!!
!!     use Eos_interface, ONLY:  Eos_getParameters
!!     logical :: inputUnchanged, needTempInput
!!     ...     
!!     call Eos_getParameters(inputsAreUnchanged=inputUnchanged,inputTempIsGuess=needTempInput)
!!
!!
!!  NOTES
!!
!!  This interface assumes that Eos initialization has already taken place when
!!  Eos_getParameters is called.
!!
!!  Since this interface uses optional arguments, all routines calling this routine must include
!!  a use Eos_interface statement, possibly with "ONLY" attribute listing Eos_getParameters
!!  explicitly, e.g.
!!      use Eos_interface, ONLY:  Eos_getParameters
!!
!!***



subroutine Eos_getParameters(eintSwitch,inputsAreUnchanged,inputTempIsGuess,constantGammaC,&
     inputMassFracNeeded,smalle)

  use Eos_data, ONLY : eos_smalle, eos_eintSwitch
  use eos_helmData, ONLY : eos_forceConstantInput
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

#include "Flash.h"

  real,OPTIONAL,intent(OUT) :: eintSwitch
  logical,OPTIONAL,intent(OUT) :: inputsAreUnchanged
  logical,OPTIONAL,intent(OUT) :: inputTempIsGuess
  logical,OPTIONAL,intent(OUT) :: constantGammaC
  logical,OPTIONAL,intent(OUT) :: inputMassFracNeeded
  real,OPTIONAL,intent(OUT) :: smalle

  ! This assumes that runtime parameters have already been gotten.
  if (present(eintSwitch)) eintSwitch = eos_eintSwitch
  if (present(inputsAreUnchanged)) inputsAreUnchanged = eos_forceConstantInput
  if (present(inputTempIsGuess)) inputTempIsGuess = .TRUE.
  if (present(constantGammaC)) constantGammaC = .FALSE.

  if (present(inputMassFracNeeded)) inputMassFracNeeded = .TRUE.
  if (present(smalle)) smalle = eos_smalle

  return
end subroutine Eos_getParameters
