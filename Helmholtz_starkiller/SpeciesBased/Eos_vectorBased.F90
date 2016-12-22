!!****if* source/physics/Eos/EosMain/Helmholtz/SpeciesBased/Eos_vectorBased
!!
!! NAME
!!
!! Eos
!!
!! SYNOPSIS
!!
!! subroutine Eos(integer(IN) :: mode,
!!                integer(IN) :: vecLen,
!!                real(INOUT) :: eosData(vecLen*EOS_NUM),
!!            optional, integer(IN) :: vecBegin,
!!            optional, integer(IN) :: vecEnd,
!!      optional, real(IN)    :: massFrac(vecLen*NSPECIES),
!!      optional, logical(IN) :: mask(EOS_VARS+1:EOS_NUM)    )
!!
!! DESCRIPTION
!!
!!   Driver for the Helmholtz and Nadyozhin equations of state.
!!   See the NOTES section for important information about this implementation.
!!
!!  This routine applies the equation of state to thermodynamic 
!!  quantities at one or more grid cells.  The number of cells is 
!!  determined by the argument veclen.  Data is packaged for this 
!!  routine in the 1d array, eosData.  The data in eosData is organized: 
!!  1:vecLen points contain the first variable, vecLen+1:2*vecLen points 
!!  contain the second variable, and so on. The number and order of
!!  variables in the array is determined by the constants defined in Eos.h.
!!  
!!  The routine takes different quantities as givens depending on the
!!  value of the mode variable: if mode=MODE_DENS_TEMP, density and
!!  temperature are taken as given, and pressure and energy are generated
!!  as output; if mode=MODE_DENS_EI, density and internal energy are taken as
!!  givens, and pressure and temperature are generated as output.  If
!!  mode=MODE_DENS_PRES, density and pressure are taken as givens, and
!!  energy and temperature are generated as output.
!!  
!!  In addition to pressure, temperature, and internal energy, which are
!!  always thermodynamically consistent after this call, other quantities
!!  such as the various thermodynamic partial derivatives can be
!!  calculated based on the values in the argument, mask.  mask is a
!!  logical array with one entry per quantity, with the order determined
!!  by constants defined in Eos.h (the same as those for the eosData
!!  argument); .true. means return the quantity, .false. means don't.
!!
!!  ARGUMENTS
!! 
!!  mode :    Selects the mode of operation of the Eos unit.
!!             The valid values are MODE_DENS_EI, MODE_DENS_PRES and  
!!             MODE_DENS_TEMP as decribed above.
!!
!!  vecLen   : number of points for each input variable
!!
!!  eosData  : This array is the data structure through which variable values are 
!!             passed in and out of the Eos routine. The arrays is sized as 
!!             EOS_NUM*vecLen. EOS_NUM, and individual input and output
!!             Eos variables are defined in Eos.h. The array is organizes such that
!!             the first 1:vecLen entries represent the first Eos variable, vecLen+1:
!!             2*vecLen represent the second Eos variable and so on. 
!!
!!  vecBegin : Index of first cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             If not present, the default is 1.
!!             Values other than 1 are untested.
!!  vecEnd   : Index of last cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             If not present, the default is vecLen.
!!             Values other than vecLen are untested.
!!
!!  massFrac : Contains the mass fractions of the species included in
!!             the simulation. The array is sized as NSPECIES*vecLen.
!!
!!  mask     : Mask is a logical array the size of EOS_DERIVS (number
!!              of partial derivatives that can be computed, defined in
!!              Eos.h), where each index represents a specific partial derivative
!!              that can be calculated by the Eos unit. A .true. value in mask 
!!              results in the corresponding derivative being calculated and 
!!              returned. It should preferably be dimensioned as
!!              mask(EOS_VARS+1:EOS_NUM) in the calling routine 
!!              to exactly match the arguments declaration in Eos Unit.
!!             Note that the indexing of mask does not begin at 1, but rather at one past
!!             the number of variables.
!!
!!             The Helmholtz EOS kernel calculation ignores the mask setting and calculates
!!             all derivatives, whether needed or not.  This routine does not return
!!             derivatives if the mask is requested, but the calculation is not speeded up
!!             by setting the mask.
!!
!!
!! PARAMETERS
!!
!!  eos_tol    Controls the accuracy of the Newton Rhapson iterations for MODE_DENS_EI and 
!!             MODE_DENS_PRES.
!!
!! EXAMPLE
!!
!! --- A single-point at a time example, does not calculate derivatives (based on Cellular Simulation)---
!!
!!  #include "constants.h"   ! for MODE_DENS_TEMP
!!  #include "Flash.h"       ! for NSPECIES
!!  #include "Eos.h"         ! for EOS_VAR order
!!
!!  real  :: temp_zone, rho_zone, ptot, eint, gamma
!!  real, dimension(EOS_NUM)  :: eosData
!!  real, dimension(SPECIES_BEGIN:SPECIES_END) ::  massFraction  
!!  integer, dimension(2,MDIM)                 :: blockRange,blockExtent
!!
!!
!!  massFraction(:) = 1.0e-12        
!!  massFraction(C12_SPEC) = 1.0
!!
!!  .... initiale temp_zone, rho_zone
!!
!!  call Grid_getBlkIndexLimits(blockId,blockRange,blockExtent)
!!  do k = blockRange(LOW,KAXIS), blockRange(HIGH,KAXIS)
!!     do j = blockRange(LOW,JAXIS),blockRange(HIGH,JAXIS)
!!        do i = blockRange(LOW,IAXIS),blockRange(HIGH,IAXIS)
!!
!!           eosData(EOS_TEMP) = temp_zone
!!           eosData(EOS_DENS) = rho_zone
!!
!!           call Eos(MODE_DENS_TEMP,1,eosData,massFraction)
!!
!!           ptot = eosData(EOS_PRES)
!!           eint = eosData(EOS_EINT)
!!           gamma = eosData(EOS_GAMC)
!!           
!!           call Grid_putPointData(blockId,CENTER,TEMP_VAR,EXTERIOR,iPosition,temp_zone)
!!           call Grid_putPointData(blockId,CENTER,DENS_VAR,EXTERIOR,iPosition,rho_zone)
!!           call Grid_putPointData(blockId,CENTER,PRES_VAR,EXTERIOR,iPosition,ptot)
!!           call Grid_putPointData(blockId,CENTER,EINT_VAR,EXTERIOR,iPosition,eint)
!!               if you want ENER_VAR, calculate it from eint and kinetic energy
!!           call Grid_putPointData(blockId,CENTER,GAMC_VAR,EXTERIOR,iPosition,gamma)
!!           call Grid_putPointData(blockId,CENTER,GAME_VAR,EXTERIOR,iPosition,(ptot/(etot*sim_rhoAmbient) + 1.0))
!!
!!         enddo  ! end of k loop
!!     enddo     ! end of j loop
!!  enddo        ! end of i loop
!!
!! ------------------ Row at a time example, with derivates (based on Eos_unitTest) --------
!!
!!  #include "constants.h"   ! for MODE_DENS_TEMP
!!  #include "Flash.h"       ! for NSPECIES, EOS_NUM
!!  #include "Eos.h"         ! for EOS_VAR order
!!  integer veclen, isize, jsize, ksize, i,j,k, e
!!  real, dimension(:), allocatable :: eosData
!!  real, dimension(:), allocatable :: massFrac
!!  logical, dimension (EOS_VARS+1:EOS_NUM) :: mask
!!  real, allocatable, dimension(:,:,:,:) :: derivedVariables
!!  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
!!
!!   ! in the Eos_unitTest, this loops over all blocks.... here is a snippet from inside
!!     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
!!
!!    !  Allocate the necessary arrays for an entire block of data
!!    isize = (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1)
!!    jsize = (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1)
!!    ksize = (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1)
!!    vecLen=isize
!!    allocate(derivedVariables(isize,jsize,ksize,EOS_NUM))
!!    allocate(eosData(vecLen*EOS_NUM))
!!    allocate(massFrac(vecLen*NSPECIES))
!!    mask = .true.
!!
!!    ! indices into the first location for these variables
!!    pres = (EOS_PRES-1)*vecLen
!!    dens = (EOS_DENS-1)*vecLen
!!    temp = (EOS_TEMP-1)*vecLen
!!
!!
!!    call Grid_getBlkPtr(blockID,solnData)
!!    do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
!!        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH, JAXIS)
!!           do i = 1,vecLen
!!              massFrac((i-1)*NSPECIES+1:i*NSPECIES) = &
!!                   solnData(SPECIES_BEGIN:SPECIES_END,ib+i-1,j,k)
!!           end do
!!
!!           eosData(pres+1:pres+vecLen) =  solnData(PRES_VAR,ib:ie,j,k)
!!           eosData(dens+1:dens+vecLen) =  solnData(DENS_VAR,ib:ie,j,k)
!!           ! Eos Helmholtz needs a good initial estimate of temperature no matter what the mode
!!           eosData(temp+1:temp+vecLen) =  solnData(TEMP_VAR,ib:ie,j,k)
!!
!!           call Eos(MODE_DENS_PRES,vecLen,eosData,massFrac,mask)
!!
!!           do e=EOS_VARS+1,EOS_NUM
!!              m = (e-1)*vecLen
!!              derivedVariables(1:vecLen,j-NGUARD,k-NGUARD,e) =  eosData(m+1:m+vecLen)
!!           end do
!!        end do
!!     end do
!!
!! NOTES
!!
!!  NSPECIES is defined in Flash.h.
!!
!!  EOS_VARS and EOS_NUM  are defined in Eos.h.
!!  Calling funtions should included Eos.h, in order to get the definitions of
!!  Eos-specific constants to be able to populate the eosData and mask arrays.
!!  
!!  MODE_DENS_TEMP, MODE_DENS_EI, and MODE_DENS_PRES are defined in constants.h.
!!
!!  All routines calling this routine should include a 
!!  use Eos_interface 
!!  statement, preferable with "ONLY" attribute e.g.
!!  use Eos_interface, ONLY:  Eos
!!
!!  The Helmholtz equation of state calculations are iterative for any mode other
!!  than MODE_DENS_TEMP.  Therefore, the intial estimates for temperature and density
!!  must be pretty good upon entering Eos with any other MODE_....or the calculations will
!!  not converge.
!!
!!  This algorithm uses a data table helm_table.dat which contains the coefficients for
!!  one of the interpolating algorithms.  Upon first entry to the Eos, a binary version of this
!!  table (helm_table.bdat) is created for speed of access.  This binary file should NOT be
!!  carried across machine architectures or even compilers.
!!
!!  When USE_EOS_YE is defined, this routine is replaced by the one in 
!!  physics/Eos/EosMain/Helmholtz/Ye
!!
!!  When operating in MODE_DENS_EI, the INPUT energy is updated.  This change of an input parameter
!!     can be overridden by setting the runtime parameter eos_forceConstantInput to true.
!!     Noted below, see comments prefaced with ConstantInput.
!!  Similarly, when operating in MODE_DENS_PRES, the INPUT pressure is updated.  Physicists need
!!     to be aware of this.  Similarly can be overridden with the runtime parameter/
!!
!!  The accuracy can be adjusted with the parameter eos_tol.
!!
!! SEE ALSO
!! 
!!  Eos.h    defines the variables used.
!!  Eos_wrapped  sets up the data structure.
!!
!!
!!*** 


subroutine Eos(mode, vecLen, eosData, massFrac, mask, vecBegin,vecEnd)

  use Driver_interface, ONLY : Driver_abortFlash
  use Multispecies_interface, ONLY : Multispecies_getSumInv, &
    Multispecies_getSumFrac
  use Logfile_interface, ONLY:  Logfile_stampMessage
  use eos_helmInterface, ONLY : eos_helm

  use Eos_data, ONLY: eos_smallT, eos_tol, eos_maxNewton, &
       eos_forceConstantInput, eos_singleSpeciesA, eos_singleSpeciesZ
  use eos_vecData, ONLY:  tempRow, denRow, etotRow, abarRow, zbarRow, &
       gamcRow, ptotRow, deaRow, dezRow, &
       detRow, dptRow, dpdRow, dedRow, pelRow, neRow, etaRow, cvRow, cpRow

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#ifdef FLASH_MULTISPECIES
#include "Multispecies.h"
#endif

  !     Arguments
  integer, INTENT(in) :: mode, vecLen
  real, INTENT(inout), dimension(vecLen*EOS_NUM) :: eosData
  real, optional,INTENT(in), dimension(vecLen*NSPECIES) :: massFrac
  ! must correspond to dimensions of Eos_wrapped
  logical,optional, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::mask
  integer,optional,INTENT(in) :: vecBegin,vecEnd

  integer :: i, k
  integer :: ilo, ihi
  integer :: pres, temp, dens, gamc, eint
  integer :: abar, zbar
  integer :: dpt, dpd, det, ded, dea, dez, pel, ne, eta, c_v, c_p
  real    :: abarInv, zbarFrac

  ! declare some local storage for the results of the Newton iteration
  real,dimension(vecLen)::  ewantRow, tnew, error, runAll, pwantRow
  !  local storage for forcedConstantInput -- could be allocatable, but might be so slow
  !  that it's not worth the small storage save.
  real,dimension(vecLen)::  psaveRow, esaveRow
  logical  :: firstCall = .true.

  !      Fill the pipe with the initial temperature, density, and composition.
  !      The composition is parametrized by abar and zbar, which are computed
  !      from the mass fractions xn().

  !  if you are working with electron abundance mass scalars, then you don't
  !  necessarily have to have mass fractions.
  if(.not.present(massFrac)) then
     call Driver_abortFlash("[Eos] Helmholtz needs mass fractions")
  end if

   ! warn user if mask is used. Helmholtz EOS doesn't use it
  if (present(mask)) then
     if (firstCall) then
        
        call Logfile_stampMessage("[EOS Helmholtz] WARNING!  Mask setting does not speed up Eos Helmholtz calls")
        write(*,*)"[EOS Helmholtz] WARNING!  Mask setting does not speed up Eos Helmholtz calls"
        firstCall = .false.
     end if
  end if

  if (present(vecBegin)) then
     ilo = vecBegin
  else
     ilo = 1
  end if
  if (present(vecEnd)) then
     ihi = vecEnd
  else
     ihi = vecLen
  end if
#ifdef DEBUG_EOS
  if (ilo < 1 .OR. ilo > vecLen) then
     print*,'[Eos vectorBased] ilo is',ilo
     call Driver_abortFlash("[Eos vectorBased] invalid ilo")
  end if
  if (ihi < 1 .OR. ihi > vecLen) then
     print*,'[Eos vectorBased] ihi is',ihi
     call Driver_abortFlash("[Eos vectorBased] invalid ihi")
  end if
#endif

  runAll = 1.0  ! This array forces eos_helm to be completely run

  ! These integers are indexes into the lowest location in UNK that contain the appropriate variable
  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  eint = (EOS_EINT-1)*vecLen   ! in flash2 eos_helm, this is etot
  gamc = (EOS_GAMC-1)*vecLen   ! in flash2 eos_helm, this is gamc
  abar = (EOS_ABAR-1)*vecLen   
  zbar = (EOS_ZBAR-1)*vecLen   

  !! For allocatable arrays, set them up now.
#ifndef FIXEDBLOCKSIZE
  call eos_vecAlloc(vecLen)
#endif

  ! retrieve data from the eosData array; local variables of already of size vecLen
  tempRow(1:vecLen)    = eosData(temp+1:temp+vecLen)
  denRow(1:vecLen)     = eosData(dens+1:dens+vecLen)
  ewantRow   = eosData(eint+1:eint+vecLen)   ! store desired internal energy for mode=2 case
  pwantRow   = eosData(pres+1:pres+vecLen)   ! store desired pressure for mode=3 case

  ! DEV should we even be filling in this data instead of pulling it out of the eosData array?
  ! DEV in Eos.F90, we assume the user knows what he's doing.  Eos_wrapped does not.
#ifdef FLASH_MULTISPECIES
  do k = 1, vecLen
     !Calculate the inverse in a way that allows for zero mass fractions
     call Multispecies_getSumInv(A, abarInv,massFrac((k-1)*NSPECIES+1:k*NSPECIES))
     abarRow(k) = 1.e0 / abarInv

     call Multispecies_getSumFrac(Z, zbarFrac, massFrac((k-1)*NSPECIES+1:k*NSPECIES))
     zbarRow(k) = abarRow(k) * zbarFrac
  end do
#else
  ! No multispecies defined, use default values (same as Gamma formulation)
  abarRow(1:vecLen) = eos_singleSpeciesA
  zbarRow(1:vecLen) = eos_singleSpeciesZ
#endif
  ! fill up eosData with appropriate z,a
  eosData(abar+1:abar+vecLen) = abarRow(1:vecLen) 
  eosData(zbar+1:zbar+vecLen) = zbarRow(1:vecLen)


  !! Save the input parameters if eos_forceConstantInput is defined
  if (eos_forceConstantInput) then
     esaveRow = ewantRow
     psaveRow = pwantRow
  end if



  !==============================================================================

  !      MODE_DENS_TEMP  temperature and density given

  !      Crank the EOS on the pipes filled above, then fill the FLASH arrays
  !      with the thermodynamic quantities returned by the EOS.

  if (mode==MODE_DENS_TEMP) then
     ! No changes for MODE_DENS_TEMP to vectorize
     call eos_helm(ilo,ihi,runAll,mask)

     eosData(pres+1:pres+vecLen)=ptotRow(1:vecLen)
     eosData(eint+1:eint+vecLen)=etotRow(1:vecLen)
     eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)

     !==============================================================================
     !      MODE_DENS_EI  internal energy and density given

  else if (mode==MODE_DENS_EI) then
     ! Initialize the errors
     error = 0.0e0

     ! Do the first eos call with all the zones in the pipe
     !  NOTE that eos_helm can ONLY operate in the equivalent of
     !  MODE_DENS_TEMP, as it returns pressure, energy and derivatives only
     !  So if you send in a crappy temperature here, you'll get a crappy starting
     !  position and the iteration won't converge.
     !  Initial temperature here is what is stored in the grid, even though we 
     !    SUSPECT this is not in equilibrium (or we wouldn't be calling Eos if it was fine)

     call eos_helm(ilo,ihi,runAll,mask)
     !  Now eos_helm has returned ptotRow, etotRow, detRow, and gamcRow


     !  Create initial condition with vectorized version
     !  ewantRow is our desired EI input
     
     tnew = tempRow(1:vecLen) - ( (etotRow(1:vecLen)-ewantRow) / detRow(1:vecLen))

     ! Don't allow the temperature to change by more than an order of magnitude 
     ! in a single iteration
     tnew = min(tnew,tempRow(1:vecLen)*10.0)  
     tnew = max(tnew,tempRow(1:vecLen)*0.01)

     ! Compute the error.  Note error is already of size vecLen
     error = abs( (tnew - tempRow(1:vecLen)) / tempRow(1:vecLen) )

     ! Store the new temperature
     tempRow(1:vecLen) = tnew

     ! Check if we are freezing, if so set the temperature to smallt, and adjust 
     ! the error so we don't wait for this one
     where (tempRow(1:vecLen) .LT. eos_smallt)
        tempRow(1:vecLen) = eos_smallt
        error = 0.1*eos_tol
     end where


     ! Loop over newton iterations now
        do i = 2, eos_maxNewton
           if (maxval(error) .LT. eos_tol) goto 70

            ! At this point is would be ideal if we could run eos_helm ONLY for 
            !  locations where error > eos_tol
           call eos_helm(ilo,ihi,error,mask)

           tnew = tempRow(1:vecLen) - ((etotRow(1:vecLen) - ewantRow) &
                / detRow(1:vecLen))

           ! Don't allow the temperature to change by more than an order of magnitude 
           ! in a single iteration
            tnew = min(tnew,tempRow(1:vecLen)*10.0)  
            tnew = max(tnew,tempRow(1:vecLen)*0.01)

           ! Compute the error
           error = abs((tnew - tempRow(1:vecLen)) / tempRow(1:vecLen))

           ! Store the new temperature
           tempRow(1:vecLen) = tnew

           ! Check if we are freezing, if so set the temperature to eos_smallt, and adjust 
           ! the error so we don't wait for this one
           where (tempRow(1:vecLen) .LT. eos_smallt) 
              tempRow(1:vecLen) = eos_smallt
              error   = 0.1*eos_tol
           end where

        end do  ! end of Newton iterations loop.  Failure drops below, success goes to 70

        ! Land here if too many iterations are needed -- failure

        write(*,*) ' '
        write(*,*) 'Newton-Raphson failed in routine Helmholtz Eos'
        write(*,*) '(e and rho as input):'
        write(*,*) ' '
        write(*,*) 'too many iterations'
        write(*,*) ' '
        write(*,*) ' temp = ', tempRow(k)
        write(*,*) ' dens = ', denRow(k)
        write(*,*) ' pres = ', ptotRow(k)
        call Driver_abortFlash('[Eos] Error: too many iterations in Helmholtz Eos')


        ! Land here if the Newton iteration converged
        !  jumps out of the Newton iterations, 
        !  but then continues to the next vector location k

70      continue
!opt     end do  ! end of vector 1:vecLen

     ! Crank through the entire eos one last time
     
     call eos_helm(ilo,ihi,runAll,mask)

     ! Fill the FLASH arrays with the results.  
     !  In MODE_DENS_EI, we should be generating temperature and pressure (plus gamma)
     eosData(temp+1:temp+vecLen)=tempRow(1:vecLen)
     eosData(pres+1:pres+vecLen)=ptotRow(1:vecLen)
     eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
     !  Update the energy to be the true energy, instead of the energy we were trying to meet
     !  ConstantInput LBR and KW believe this is WRONG -- the input arrays should not be changed
     if (eos_forceConstantInput)  then
        eosData(eint+1:eint+vecLen) = esaveRow(1:vecLen)
     else
        eosData(eint+1:eint+vecLen) = etotRow(1:vecLen)
     end if


     !==============================================================================

     !      MODE_DENS_PRES  pressure and density given

  else if (mode==MODE_DENS_PRES) then

     ! Initialize the errors
     error = 0.0e0

     ! Do the first eos call with all the zones in the pipe
     call eos_helm(ilo,ihi,runAll,mask)

     ! begin initialization with vectorization
     tnew = tempRow(1:vecLen) - ( (ptotRow(1:vecLen)-pwantRow) / dptRow(1:vecLen))

     ! Don't allow the temperature to change by more than an order of magnitude 
     ! in a single iteration
     tnew = min(tnew,tempRow(1:vecLen)*10.0)  
     tnew = max(tnew,tempRow(1:vecLen)*0.01)

     ! Compute the error.  Note error is already of size vecLen
     error = abs( (tnew - tempRow(1:vecLen)) / tempRow(1:vecLen) )

     ! Store the new temperature
     tempRow(1:vecLen) = tnew


     ! Check if we are freezing, if so set the temperature to smallt, and adjust 
     ! the error so we don't wait for this one
     where (tempRow(1:vecLen) .LT. eos_smallt)
        tempRow(1:vecLen) = eos_smallt
        error = 0.1*eos_tol
     end where

     ! Loop over the zones individually now
     
     do i = 2, eos_maxNewton

        if (maxval(error) .LT. eos_tol) goto 170

        ! At this point is would be ideal if we could run eos_helm ONLY for 
        !  locations where error > eos_tol
        call eos_helm(ilo,ihi,error,mask)

        tnew = tempRow(1:vecLen) - ((ptotRow(1:vecLen) - pwantRow) &
             / dptRow(1:vecLen))


        ! Don't allow the temperature to change by more than an order of magnitude 
        ! in a single iteration
        tnew = min(tnew,tempRow(1:vecLen)*10.0)  
        tnew = max(tnew,tempRow(1:vecLen)*0.01)

        ! Compute the error
        error = abs((tnew - tempRow(1:vecLen)) / tempRow(1:vecLen))

        ! Store the new temperature
        tempRow(1:vecLen) = tnew

        ! Check if we are freezing, if so set the temperature to eos_smallt, and adjust 
        ! the error so we don't wait for this one
        where (tempRow(1:vecLen) .LT. eos_smallt) 
           tempRow(1:vecLen) = eos_smallt
           error   = 0.1*eos_tol
        end where

     enddo

     ! Land here if too many iterations are needed, otherwise success happens at 170

     write(*,*) ' '
     write(*,*) 'Newton-Raphson failed in Helmholtz Eos:'
     write(*,*) '(p and rho as input)'
     write(*,*) ' '
     write(*,*) 'too many iterations'
     write(*,*) ' '
     write(*,*) ' k    = ', k,ilo,ihi
     write(*,*) ' temp = ', tempRow(k)
     write(*,*) ' dens = ', denRow(k)
     write(*,*) ' pres = ', ptotRow(k)

     call Driver_abortFlash('[Eos] Error: too many iteration in Eos')

     ! Land here if the Newton iteration converged

170  continue

     ! Crank through the entire eos one last time

     call eos_helm(ilo,ihi,runAll,mask)

     ! Fill the FLASH arrays with the results.  
     eosData(temp+1:temp+vecLen)=tempRow(1:vecLen)
     eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
     eosData(eint+1:eint+vecLen)=etotRow(1:vecLen)
     ! Update the pressure to be the equilibrium pressure, instead of the pressure we were trying to meet
     !  ConstantInput LBR and KW believe this is wrong.  See notes at the top of the routine
     if (eos_forceConstantInput) then
        eosData(pres+1:pres+vecLen) = psaveRow(1:vecLen)
     else
        eosData(pres+1:pres+vecLen) = ptotRow(1:vecLen)
     end if


     !==============================================================================

     ! Unknown EOS mode selected

  else
     call Driver_abortFlash('[Eos] Error: unknown input in routine eos')
  end if

  ! Get the optional values
  if(present(mask)) then
     if(mask(EOS_DPT)) then
        dpt = (EOS_DPT-1)*vecLen
         eosData(dpt+1:dpt+vecLen) = dptRow(1:vecLen)
     end if
     if(mask(EOS_DPD)) then
        dpd = (EOS_DPD-1)*vecLen
        eosData(dpd+1:dpd+vecLen) = dpdRow(1:vecLen)
     end if
     if(mask(EOS_DET))then
        det = (EOS_DET-1)*vecLen
         eosData(det+1:det+vecLen) = detRow(1:vecLen)
     end if
     if(mask(EOS_DED))then 
        ded = (EOS_DED-1)*vecLen
        eosData(ded+1:ded+vecLen) = dedRow(1:vecLen)
     end if
     if(mask(EOS_DEA))then 
        dea = (EOS_DEA-1)*vecLen
        eosData(dea+1:dea+vecLen) = deaRow(1:vecLen)
     end if
     if(mask(EOS_DEZ))then 
        dez = (EOS_DEZ-1)*vecLen
        eosData(dez+1:dez+vecLen) = dezRow(1:vecLen)
     end if
     if(mask(EOS_PEL))then 
        pel = (EOS_PEL-1)*vecLen
        eosData(pel+1:pel+vecLen) = pelRow(1:vecLen)
     end if
     if(mask(EOS_NE))then 
        ne = (EOS_NE-1)*vecLen
        eosData(ne+1:ne+vecLen) = neRow(1:vecLen)
     end if
     if(mask(EOS_ETA))then 
        eta = (EOS_ETA-1)*vecLen 
        eosData(eta+1:eta+vecLen) = etaRow(1:vecLen)
     end if
     
     if(mask(EOS_CV))then
        if(mask(EOS_DET)) then
           c_v = (EOS_CV-1)*vecLen
           eosData(c_v+1:c_v+vecLen) = cvRow(1:vecLen)
        else
           call Driver_abortFlash("[Eos] cannot calculate C_V without DET.  Set mask appropriately.")
        end if
     end if
     
     if(mask(EOS_CP))then
        if(mask(EOS_CV).and.mask(EOS_DET)) then
           c_p = (EOS_CP-1)*vecLen
           eosData(c_p+1:c_p+vecLen) = cpRow(1:vecLen)
        else
           call Driver_abortFlash("[Eos] cannot calculate C_P without C_V and DET.  Set mask appropriately.")
        end if
     end if
  end if




  !! Close up arrays if previously allocated
#ifndef FIXEDBLOCKSIZE  
  call eos_vecDealloc()
#endif

  return

end subroutine Eos


