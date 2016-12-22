!!****if* source/physics/Eos/EosMain/Helmholtz/eos_vecAlloc
!!
!! NAME
!!
!!  eos_vecAlloc
!! 
!! SYNOPSIS
!!
!!  eos_vecAlloc(integer(IN) :: vecLen)
!!
!! DESCRIPTION
!!
!!  Allocate vectors for Helmholtz
!!
!! ARGUMENTS
!!
!!  vecLen -- size of vectors to be allocated
!!
!!
!!*** 

subroutine eos_vecAlloc(vecLen)

#include "constants.h"
#include "Flash.h"

  use eos_vecData, ONLY: tempRow, denRow, abarRow, zbarRow, &
       ptotRow, dptRow, dpdRow, etotRow,&
       detRow, dedRow,  deaRow, dezRow, &
       dstRow, dsdRow, stotRow, &
       pelRow, neRow, etaRow, cpRow, cvRow, gamcRow

  implicit none
  integer, intent(IN) :: vecLen

#ifndef FIXEDBLOCKSIZE

  allocate(tempRow(vecLen))
  allocate(denRow(vecLen))
  allocate(abarRow(vecLen))
  allocate(zbarRow(vecLen))

  !..totals and their derivatives
  allocate(ptotRow(vecLen))
  allocate(dptRow(vecLen))
  allocate(dpdRow(vecLen))
! UNUSED  allocate(dpaRow(vecLen))
! UNUSED  allocate(dpzRow(vecLen))
  allocate(etotRow(vecLen))
  allocate(detRow(vecLen))
  allocate(dedRow(vecLen))
  allocate(deaRow(vecLen))
  allocate(dezRow(vecLen))
! UNUSED  allocate(deaRow(vecLen))
! UNUSED  allocate(dezRow(vecLen))
  allocate(stotRow(vecLen))
  allocate(dstRow(vecLen))
  allocate(dsdRow(vecLen))
! UNUSED  allocate(dsaRow(vecLen))
! UNUSED  allocate(dszRow(vecLen))


  !..radiation contributions 
!UNUSED  allocate(pradRow(vecLen))
!UNUSED  allocate(eradRow(vecLen))
!UNUSED  allocate(sradRow(vecLen))

  !..ion contributions 

!UNUSED  allocate(pionRow(vecLen))
!UNUSED  allocate(eionRow(vecLen))
!UNUSED  allocate(sionRow(vecLen))
!UNUSED  allocate(xniRow(vecLen))

  !..electron-positron contributions -- most UNUSED and REMOVED
  !  See eos_fixedVecData if you want to see the original names
  allocate(pelRow(vecLen))
  allocate(neRow(vecLen))
  allocate(etaRow(vecLen))


  !..coulomb contributions -- all UNUSED and REMOVED

  !..thermodynamic consistency checks; maxwell relations -- all UNUSED and REMOVED

  !..derivative based quantities
  allocate(cpRow(vecLen))
  allocate(cvRow(vecLen))
  allocate(gamcRow(vecLen))

#endif

end subroutine eos_vecAlloc
