!!****if* source/physics/Eos/EosMain/Helmholtz/eos_helm
!!
!! NAME
!!
!!  eos_helm
!!
!!
!! SYNOPSIS
!!
!!  call eos_helm(integer(IN) :: eos_jlo,
!!                integer(IN) :: eos_jhi,
!!      optional, logical(IN) :: mask(EOS_VARS+1:EOS_NUM) )
!!
!!
!! DESCRIPTION
!!
!!  given a temperature temp [K], density den [g/cm**3], and a composition 
!!  characterized by abar and zbar, this routine returns some of the other 
!!  thermodynamic quantities. of prime interest is the pressure [erg/cm**3], 
!!  specific thermal energy [erg/gr], and their derivatives with respect
!!  to temperature and density.
!!
!!  NOTE that this routine does not have any Flash3 eosMode other than 
!!    MODE_DENS_TEMP.  So the initial temperature and density input to Eos() 
!!    for the Helmholtz unit must be a good guess, or the algorithms will never
!!    converge.
!!  
!!  for opacity purposes, the number density of electron-positrons,
!!  the chemical potential of the electrons, and the pressure due
!!  to electrons and positrons is also returned
!!  
!!  this routine uses planckian photons, an ideal gas of ions,
!!  and a fully ionized electron-positron gas with an arbitrary degree 
!!  of relativity and degeneracy. interpolation in a table of the 
!!  helmholtz free energy is used to return the electron-positron 
!!  thermodynamic quantities.
!!  
!!  references: cox & giuli chapter 24 ; timmes & swesty apj 1999
!!
!! ARGUMENTS
!!
!!  eos_jlo  -- low index of calculation range
!!  eos_jhi  -- high index of calculation range
!!  mask     --  Mask is a logical array the size of EOS_DERIVS (number
!!              of partial derivatives that can be computed, defined in
!!              Eos.h), where each index represents a specific partial derivative
!!              that can be calculated by the Eos unit. A .true. value in mask 
!!              results in the corresponding derivative being calculated and 
!!              returned. It should preferably be dimensioned as
!!              mask(EOS_VARS+1:EOS_NUM) in the calling routine 
!!              to exactly match the arguments declaration in Eos Unit.
!!             Note that the indexing of mask does not begin at 1, but rather at one past
!!             the number of basic variables.
!!
!!             The Helmholtz EOS kernel calculation ignores the mask setting and calculates
!!             all derivatives, whether needed or not.  This routine does not return
!!             derivatives if the mask is requested, but the calculation is not speeded up
!!             by setting the mask.
!!
!!
!! PARAMETERS
!!
!! NOTES
!!
!!  equation of state communication through eos_vecData
!!
!!  btemp    = temperature
!!  den      = density
!!  abar     = average number of nucleons per nuclei
!!  zbar     = average number of protons per nuclei
!!  z2bar    = square of zbar
!!  ytot1    = total number of moles per gram
!!  ye       = electron mole number
!!
!!  Since this subroutine has an optional argument, an explicit interface is needed.
!!  Calling program units should therefore have a line like
!!    use eos_helmInterface, ONLY : eos_helm
!!  no matter whether they call eos_helm with or without the optional mask argument.
!!***
#ifdef DEBUG_ALL
#define DEBUG_EOS
#endif

subroutine eos_helm(eos_jlo,eos_jhi,mask)

  use eos_helmData, ONLY: eos_f, eos_ft, eos_ftt, eos_fd, eos_fdd, eos_fdt, &
       eos_fddt, eos_fdtt, eos_fddtt, &
       eos_dpdf, eos_dpdft, eos_dpdfd, eos_dpdfdt, &
       eos_coulombMult, eos_dLo, eos_tLo, eos_t, eos_dt,  eos_dtInv, &
       eos_d, eos_dd, eos_ddInv, eos_tstpi, eos_dstpi, &
       eos_dtSqr, eos_dtSqrInv, eos_ddSqr, &
       eos_ef, eos_eft, eos_efd, eos_efdt, &
       eos_xf, eos_xft, eos_xfd, eos_xfdt, &
       EOSJMAX, EOSIMAX, eos_coulombAbort
  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Timers_interface, ONLY: Timers_start, Timers_stop
  use eos_vecData, ONLY: tempRow , denRow, abarRow, zbarRow, &
       ptotRow, etotRow, stotRow, &
       dpdRow, dptRow, dstRow, dedRow, detRow, dsdRow, &
       deaRow, dezRow, & !Calhoun
       pelRow, neRow, etaRow, gamcRow, cvRow, cpRow

  !! physical constants to high precision
  use eos_helmConstData, ONLY: kerg, kergavo, asoli3, avo, avoInv, sioncon, pi
  use Eos_data, ONLY : eos_meshMe

  implicit none


#include "constants.h" 
#include "Flash.h"  
#include "Eos.h"

!! Arguments
  integer, intent(IN) :: eos_jlo, eos_jhi
  logical,optional, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::mask
!! from eos_helmInterface.F90       logical,optional, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::mask

!! Local variables

  character*90 ::  internalFile
  integer      ::  i,j

  real         ::  btemp,den,abar,zbar,z2bar,ytot1,ye
  real         ::  local_coulombMult

  real         ::  x1,x2,x3,x4,x5,x6,x7, &
                   y0,y1,y2,y3,y4, &  ! don't reuse variables -- it's confusing!
       deni,tempi,kt, & 
       prad,dpraddd,dpraddt,erad,deraddd,deraddt, &
       srad, dsraddd, dsraddt, & 
       xni,pion,dpiondd,dpiondt,eion,deiondd, deiondt,& 
       sion,dsiondd, dsiondt, & 
       pele,dpepdd,dpepdt,eele,deepdd,deepdt, &
       sele, dsepdd,dsepdt, & 
       pres,dpresdd,dpresdt,ener,denerdt, &
       entr, & 
       presi,chit,chid,gamc,kavoy

 real :: cv, cp, etaele, xnefer,          &
         denerdd,dentrdd,dentrdt


!!  for the interpolations
  integer          iat,jat
  real             free,df_d,df_t,df_dd,df_tt,df_dt
  real             xt,xd,mxt,mxd, & 
       si0t,si1t,si2t,si0mt,si1mt,si2mt, & 
       si0d,si1d,si2d,si0md,si1md,si2md, & 
       dsi0t,dsi1t,dsi2t,dsi0mt,dsi1mt,dsi2mt, & 
       dsi0d,dsi1d,dsi2d,dsi0md,dsi1md,dsi2md, & 
       ddsi0t,ddsi1t,ddsi2t,ddsi0mt,ddsi1mt,ddsi2mt, & 
       zFunc,z0,z1,z2,z3,z4,z5,z6, & ! Split up confusing calculations
       psi0,dpsi0,ddpsi0,psi1,dpsi1,ddpsi1,psi2, & 
       dpsi2,ddpsi2,din,h5, & 
       xpsi0,xpsi1,h3dpd, & 
       w0t,w1t,w2t,w0mt,w1mt,w2mt, & 
       w0d,w1d,w2d,w0md,w1md,w2md
  real :: h3e, h3x
  real :: fi(36)


!!  for the coulomb corrections
  real ktinv,dxnidd,dsdd,lami,inv_lami,lamidd, & 
       s0,s1,s2,s3,s4, &   ! all temporary variables, reuse is confusing
       plasg,plasg_inv,plasgdd,plasgdt,a1,b1,c1,d1cc,e1cc,a2,b2,c2, & 
       ecoul,decouldd,decouldt,pcoul,dpcouldd,dpcouldt, & 
       scoul,dscouldd,dscouldt

!Added by Calhoun for calculations for the Aprox13t network
   real    :: deradda,dxnida,dpionda,deionda,dsepda,deepda,decoulda,&
     &     dsda,dsdda,lamida,plasgda,denerda,deraddz,deiondz,deepdz, &
     &     decouldz,dsepdz,plasgdz,denerdz
   logical :: bAprox13t  ! becomes true if variables for Aprox13t network are set



  real,   parameter :: third = 1.0e0/3.0e0, & 
       forth = 4.0e0/3.0e0, & 
       qe    = 4.8032068e-10,   & 
       esqu  = qe * qe


!!  for the uniform background coulomb correction
  data a1,b1,c1,d1cc,e1cc,a2,b2,c2 & 
       /-0.898004e0, 0.96786e0, 0.220703e0, -0.86097e0, & 
       2.5269e0  , 0.29561e0, 1.9885e0,    0.288675e0/


! --------------------------------------------------------------------------

  !! ***********Beginning of statement function declarations **********

  !! The next few statements are effectively function declarations.
  !! Normally they should have been declared as, for instance 
  !! real function psi0(z)
  !!   real,intent(in) :: z
  !!   psi0=z**3 * ( z * (-6.0e0*z + 15.0e0) -10.0e0) + 1.0e0
  !! end function psi0

!!  quintic hermite polynomial statement functions

!!  psi0 and its derivatives
  psi0(zFunc)   = zFunc**3 * ( zFunc * (-6.0e0*zFunc + 15.0e0) -10.0e0) + 1.0e0
  dpsi0(zFunc)  = zFunc**2 * ( zFunc * (-30.0e0*zFunc + 60.0e0) - 30.0e0)
  ddpsi0(zFunc) = zFunc* ( zFunc*( -120.0e0*zFunc + 180.0e0) -60.0e0)


!!  psi1 and its derivatives
  psi1(zFunc)   = zFunc* ( zFunc**2 * ( zFunc * (-3.0e0*zFunc + 8.0e0) - 6.0e0) + 1.0e0)
  dpsi1(zFunc)  = zFunc*zFunc * ( zFunc * (-15.0e0*zFunc + 32.0e0) - 18.0e0) +1.0e0
  ddpsi1(zFunc) = zFunc * (zFunc * (-60.0e0*zFunc + 96.0e0) -36.0e0)


!!  psi2  and its derivatives
  psi2(zFunc)   = 0.5e0*zFunc*zFunc*( zFunc* ( zFunc * (-zFunc + 3.0e0) - 3.0e0) + 1.0e0)
  dpsi2(zFunc)  = 0.5e0*zFunc*( zFunc*(zFunc*(-5.0e0*zFunc + 12.0e0) - 9.0e0) + 2.0e0)
  ddpsi2(zFunc) = 0.5e0*(zFunc*( zFunc * (-20.0e0*zFunc + 36.0e0) - 18.0e0) + 2.0e0)


! Define this choice if you want to decrease the operations.  It
!  does change results slightly.
#define EOS_LESSOPERATIONS
!!#undef EOS_LESSOPERATIONS
#ifndef EOS_LESSOPERATIONS
  h5(w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)= & 
       fi(1)  *w0d*w0t   + fi(2)  *w0md*w0t & 
       + fi(3)  *w0d*w0mt  + fi(4)  *w0md*w0mt & 
       + fi(5)  *w0d*w1t   + fi(6)  *w0md*w1t & 
       + fi(7)  *w0d*w1mt  + fi(8)  *w0md*w1mt & 
       + fi(9)  *w0d*w2t   + fi(10) *w0md*w2t & 
       + fi(11) *w0d*w2mt  + fi(12) *w0md*w2mt & 
       + fi(13) *w1d*w0t   + fi(14) *w1md*w0t & 
       + fi(15) *w1d*w0mt  + fi(16) *w1md*w0mt & 
       + fi(17) *w2d*w0t   + fi(18) *w2md*w0t & 
       + fi(19) *w2d*w0mt  + fi(20) *w2md*w0mt & 
       + fi(21) *w1d*w1t   + fi(22) *w1md*w1t & 
       + fi(23) *w1d*w1mt  + fi(24) *w1md*w1mt & 
       + fi(25) *w2d*w1t   + fi(26) *w2md*w1t & 
       + fi(27) *w2d*w1mt  + fi(28) *w2md*w1mt & 
       + fi(29) *w1d*w2t   + fi(30) *w1md*w2t & 
       + fi(31) *w1d*w2mt  + fi(32) *w1md*w2mt & 
       + fi(33) *w2d*w2t   + fi(34) *w2md*w2t & 
       + fi(35) *w2d*w2mt  + fi(36) *w2md*w2mt
#else
  !! Alternate attempt at reducing multiplications
  !! use define EOS_LESSOPERATIONS at the top of this file to make it work...
  h5(w0t,w1t,w2t,w0mt,w1mt,w2mt,w0d,w1d,w2d,w0md,w1md,w2md)= & 
         (fi(1)  *w0d  + fi(2)  *w0md   & 
       +  fi(17) *w2d  + fi(18) *w2md   & 
       +  fi(13) *w1d  + fi(14) *w1md) *w0t & 
       + (fi(3)  *w0d  + fi(4)  *w0md   & 
       +  fi(15) *w1d  + fi(16) *w1md   & 
       +  fi(19) *w2d  + fi(20) *w2md) *w0mt & 
       + (fi(5)  *w0d  + fi(6)  *w0md   & 
       +  fi(21) *w1d  + fi(22) *w1md   & 
       +  fi(25) *w2d  + fi(26) *w2md) *w1t & 
       + (fi(7)  *w0d  + fi(8)  *w0md   & 
       +  fi(23) *w1d  + fi(24) *w1md   & 
       +  fi(27) *w2d  + fi(28) *w2md) *w1mt & 
       + (fi(9)  *w0d  + fi(10) *w0md   & 
       +  fi(29) *w1d  + fi(30) *w1md   & 
       +  fi(33) *w2d  + fi(34) *w2md) *w2t & 
       + (fi(11) *w0d  + fi(12) *w0md   & 
       +  fi(31) *w1d  + fi(32) *w1md   & 
       +  fi(35) *w2d  + fi(36) *w2md) *w2mt
#endif

!!  cubic hermite polynomial statement functions
!!  psi0 & derivatives
  xpsi0(zFunc)  = zFunc * zFunc * (2.0e0*zFunc - 3.0e0) + 1.0
!UNUSED  xdpsi0(zFunc) = zFunc * (6.0e0*zFunc - 6.0e0)


!!  psi1 & derivatives
  xpsi1(zFunc)  = zFunc * ( zFunc * (zFunc - 2.0e0) + 1.0e0)
!UNUSED  xdpsi1(zFunc) = zFunc * (3.0e0*zFunc - 4.0e0) + 1.0e0



!!  bicubic hermite for the pressure derivative
#ifndef EOS_LESSOPERATIONS
  h3dpd(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) =  & 
       eos_dpdf(i,j)        *w0d*w0t   +   eos_dpdf(i+1,j)    *w0md*w0t  & 
       +   eos_dpdf(i,j+1)  *w0d*w0mt  +   eos_dpdf(i+1,j+1)  *w0md*w0mt & 
       +  eos_dpdft(i,j)    *w0d*w1t   +  eos_dpdft(i+1,j)    *w0md*w1t  & 
       +  eos_dpdft(i,j+1)  *w0d*w1mt  +  eos_dpdft(i+1,j+1)  *w0md*w1mt & 
       +  eos_dpdfd(i,j)    *w1d*w0t   +  eos_dpdfd(i+1,j)    *w1md*w0t  & 
       +  eos_dpdfd(i,j+1)  *w1d*w0mt  +  eos_dpdfd(i+1,j+1)  *w1md*w0mt & 
       + eos_dpdfdt(i,j)    *w1d*w1t   + eos_dpdfdt(i+1,j)    *w1md*w1t  & 
       + eos_dpdfdt(i,j+1)  *w1d*w1mt  + eos_dpdfdt(i+1,j+1)  *w1md*w1mt
#else
  h3dpd(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) =  & 
         (eos_dpdf(i,j)     *w0d  +    eos_dpdf(i+1,j)   *w0md        & 
       +  eos_dpdfd(i,j)    *w1d  +  eos_dpdfd(i+1,j)    *w1md) *w0t  & 
       + (eos_dpdf(i,j+1)   *w0d  +   eos_dpdf(i+1,j+1)  *w0md        & 
       +  eos_dpdfd(i,j+1)  *w1d  +  eos_dpdfd(i+1,j+1)  *w1md) *w0mt & 
       + (eos_dpdft(i,j)    *w0d  +  eos_dpdft(i+1,j)    *w0md        & 
       +  eos_dpdfdt(i,j)   *w1d  + eos_dpdfdt(i+1,j)    *w1md) *w1t  & 
       + (eos_dpdft(i,j+1)  *w0d  +  eos_dpdft(i+1,j+1)  *w0md        & 
       +  eos_dpdfdt(i,j+1) *w1d  + eos_dpdfdt(i+1,j+1)  *w1md) *w1mt
#endif


!!  bicubic hermite polynomial for the chemical potential
#ifndef EOS_LESSOPERATIONS
  h3e(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) =  & 
           eos_ef(i,j)    *w0d*w0t   +   eos_ef(i+1,j)    *w0md*w0t  & 
       +   eos_ef(i,j+1)  *w0d*w0mt  +   eos_ef(i+1,j+1)  *w0md*w0mt & 
       +  eos_eft(i,j)    *w0d*w1t   +  eos_eft(i+1,j)    *w0md*w1t  & 
       +  eos_eft(i,j+1)  *w0d*w1mt  +  eos_eft(i+1,j+1)  *w0md*w1mt & 
       +  eos_efd(i,j)    *w1d*w0t   +  eos_efd(i+1,j)    *w1md*w0t  & 
       +  eos_efd(i,j+1)  *w1d*w0mt  +  eos_efd(i+1,j+1)  *w1md*w0mt & 
       + eos_efdt(i,j)    *w1d*w1t   + eos_efdt(i+1,j)    *w1md*w1t  & 
       + eos_efdt(i,j+1)  *w1d*w1mt  + eos_efdt(i+1,j+1)  *w1md*w1mt
#else
  h3e(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) =  & 
          (eos_ef(i,j)    *w0d   +   eos_ef(i+1,j)    *w0md        & 
       +   eos_efd(i,j)   *w1d   +  eos_efd(i+1,j)    *w1md) *w0t  & 
       +  (eos_ef(i,j+1)  *w0d   +   eos_ef(i+1,j+1)  *w0md        & 
       +   eos_efd(i,j+1) *w1d   +  eos_efd(i+1,j+1)  *w1md) *w0mt & 
       +  (eos_eft(i,j)   *w0d   +  eos_eft(i+1,j)    *w0md        & 
       +   eos_efdt(i,j)  *w1d   + eos_efdt(i+1,j)    *w1md) *w1t  & 
       +  (eos_eft(i,j+1) *w0d   +  eos_eft(i+1,j+1)  *w0md        & 
       +   eos_efdt(i,j+1)*w1d   + eos_efdt(i+1,j+1)  *w1md) *w1mt
#endif



!!  bicubic hermite polynomial for electron positron number densities
#ifndef EOS_LESSOPERATIONS
  h3x(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) =  & 
           eos_xf(i,j)    *w0d*w0t   +   eos_xf(i+1,j)    *w0md*w0t  & 
       +   eos_xf(i,j+1)  *w0d*w0mt  +   eos_xf(i+1,j+1)  *w0md*w0mt & 
       +  eos_xft(i,j)    *w0d*w1t   +  eos_xft(i+1,j)    *w0md*w1t  & 
       +  eos_xft(i,j+1)  *w0d*w1mt  +  eos_xft(i+1,j+1)  *w0md*w1mt & 
       +  eos_xfd(i,j)    *w1d*w0t   +  eos_xfd(i+1,j)    *w1md*w0t  & 
       +  eos_xfd(i,j+1)  *w1d*w0mt  +  eos_xfd(i+1,j+1)  *w1md*w0mt & 
       + eos_xfdt(i,j)    *w1d*w1t   + eos_xfdt(i+1,j)    *w1md*w1t  & 
       + eos_xfdt(i,j+1)  *w1d*w1mt  + eos_xfdt(i+1,j+1)  *w1md*w1mt
#else
  h3x(i,j,w0t,w1t,w0mt,w1mt,w0d,w1d,w0md,w1md) =  & 
         (eos_xf(i,j)    *w0d   +   eos_xf(i+1,j)    *w0md        & 
       +  eos_xfd(i,j)   *w1d   +  eos_xfd(i+1,j)    *w1md) *w0t  & 
       + (eos_xf(i,j+1)  *w0d   +   eos_xf(i+1,j+1)  *w0md        & 
       +  eos_xfd(i,j+1) *w1d   +  eos_xfd(i+1,j+1)  *w1md) *w0mt & 
       + (eos_xft(i,j)   *w0d   +  eos_xft(i+1,j)    *w0md        & 
       +  eos_xfdt(i,j)  *w1d   + eos_xfdt(i+1,j)    *w1md) *w1t  & 
       + (eos_xft(i,j+1) *w0d   +  eos_xft(i+1,j+1)  *w0md        & 
       +  eos_xfdt(i,j+1)*w1d   + eos_xfdt(i+1,j+1)  *w1md) *w1mt
#endif

  !!************ This is the end statement function definitions *************
! -----------------------------------------------------------------------------

!!  popular format statements
03 format(1x,4(a,1pe11.3))
04 format(1x,4(a,i4))

! ------------------------------------------------------------------------------

! Initial testing of masks
!  Note that there should be things added here for entropy eventually
  bAprox13t = .false.  
  if (present(mask)) then
     if (mask(EOS_DEA).or.mask(EOS_DEZ)) &
          bAprox13t = .true.    ! we will do calculations for Aprox13t network
  end if


  !!  normal execution starts here, start of pipeline loop, no input checking
!!  biquintic hermite polynomial statement function
!  call Timers_start("eos_helm")

  do j=eos_jlo,eos_jhi

     btemp  = tempRow(j)
     den    = denRow(j)
     abar   = abarRow(j)
     zbar   = zbarRow(j)
     ytot1  = 1.0e0/abar
     ye     = ytot1 * zbar


     !!  frequent combinations
     deni    = 1.0e0/den
     tempi   = 1.0e0/btemp 
     kt      = kerg * btemp
     ktinv   = 1.0e0/kt
     kavoy   = kergavo * ytot1


     !!  radiation section:
     prad    = asoli3 * btemp * btemp * btemp * btemp
     dpraddt = 4.0e0 * prad * tempi
     dpraddd = 0.0e0

     x1      = prad * deni 
     erad    = 3.0e0 * x1
     deraddd = -erad*deni
     deraddt = 4.0e0 * erad * tempi
     ! Calhoun next two lines
     deradda = 0.0e0
     deraddz = 0.0e0

     srad    = (x1 + erad)*tempi
     dsraddd = (dpraddd*deni - x1*deni + deraddd)*tempi
     dsraddt = (dpraddt*deni + deraddt - srad)*tempi


     !!  ion section:
     dxnidd  = avo * ytot1
     xni     = dxnidd * den

     pion    = xni * kt
     dpiondd = avo * ytot1 * kt
     dpiondt = xni * kerg
     
   

     eion    = 1.5e0 * pion * deni
     deiondd = (1.5e0 * dpiondd - eion)*deni
     deiondt = 1.5e0 * xni * kerg *deni

     if (bAprox13t) then
        dxnida  = -xni*ytot1 !Calhoun
        dpionda = dxnida * kt !Calhoun
        deionda = 1.5e0 * dpionda*deni  !Calhoun
        deiondz = 0.0e0 !Calhoun
     end if

     !!  sackur-tetrode equation for the ion entropy of 
     !!  a single ideal gas characterized by abar
     x2      = abar*abar*sqrt(abar) * deni*avoinv
     y0      = sioncon * btemp
     z0      = x2 * y0 * sqrt(y0)
     sion    = (pion*deni + eion)*tempi + kavoy*log(z0)
     dsiondd = (dpiondd*deni - pion*deni*deni + deiondd)*tempi    &
          - kavoy * deni
     dsiondt = (dpiondt*deni + deiondt)*tempi   &
          - (pion*deni + eion) * tempi*tempi  &
          + 1.5e0 * kavoy * tempi




     !!  electron-positron section:
     !!  enter the table with ye*den, no checks of the input
     din = ye*den

     !!  hash locate this temperature and density
     jat = int((log10(btemp) - eos_tlo)*eos_tstpi) + 1
     jat = max(1,min(jat,EOSJMAX-1))
     iat = int((log10(din) - eos_dlo)*eos_dstpi) + 1
     iat = max(1,min(iat,EOSIMAX-1))
!     print *, 'jat = ',jat, ' iat= ', iat

     !!  access the table locations only once
     fi(1)  = eos_f(iat,jat)
     fi(2)  = eos_f(iat+1,jat)
     fi(3)  = eos_f(iat,jat+1)
     fi(4)  = eos_f(iat+1,jat+1)
     fi(5)  = eos_ft(iat,jat)
     fi(6)  = eos_ft(iat+1,jat)
     fi(7)  = eos_ft(iat,jat+1)
     fi(8)  = eos_ft(iat+1,jat+1)
     fi(9)  = eos_ftt(iat,jat)
     fi(10) = eos_ftt(iat+1,jat)
     fi(11) = eos_ftt(iat,jat+1)
     fi(12) = eos_ftt(iat+1,jat+1)
     fi(13) = eos_fd(iat,jat)
     fi(14) = eos_fd(iat+1,jat)
     fi(15) = eos_fd(iat,jat+1)
     fi(16) = eos_fd(iat+1,jat+1)
     fi(17) = eos_fdd(iat,jat)
     fi(18) = eos_fdd(iat+1,jat)
     fi(19) = eos_fdd(iat,jat+1)
     fi(20) = eos_fdd(iat+1,jat+1)
     fi(21) = eos_fdt(iat,jat)
     fi(22) = eos_fdt(iat+1,jat)
     fi(23) = eos_fdt(iat,jat+1)
     fi(24) = eos_fdt(iat+1,jat+1)
     fi(25) = eos_fddt(iat,jat)
     fi(26) = eos_fddt(iat+1,jat)
     fi(27) = eos_fddt(iat,jat+1)
     fi(28) = eos_fddt(iat+1,jat+1)
     fi(29) = eos_fdtt(iat,jat)
     fi(30) = eos_fdtt(iat+1,jat)
     fi(31) = eos_fdtt(iat,jat+1)
     fi(32) = eos_fdtt(iat+1,jat+1)
     fi(33) = eos_fddtt(iat,jat)
     fi(34) = eos_fddtt(iat+1,jat)
     fi(35) = eos_fddtt(iat,jat+1)
     fi(36) = eos_fddtt(iat+1,jat+1)


     !!  various differences
     xt  = max( (btemp - eos_t(jat))*eos_dtInv(jat), 0.0e0)
     xd  = max( (din - eos_d(iat))*eos_ddInv(iat), 0.0e0)
     mxt = 1.0e0 - xt
     mxd = 1.0e0 - xd


     !!  the density and temperature basis functions
     si0t =   psi0(xt)
     si1t =   psi1(xt)*eos_dt(jat)
     si2t =   psi2(xt)*eos_dtSqr(jat)

     si0mt =  psi0(mxt)
     si1mt = -psi1(mxt)*eos_dt(jat)
     si2mt =  psi2(mxt)*eos_dtSqr(jat)

     si0d =   psi0(xd)
     si1d =   psi1(xd)*eos_dd(iat)
     si2d =   psi2(xd)*eos_ddSqr(iat)

     si0md =  psi0(mxd)
     si1md = -psi1(mxd)*eos_dd(iat)
     si2md =  psi2(mxd)*eos_ddSqr(iat)


     !!  the first derivatives of the basis functions
     dsi0t =   dpsi0(xt)*eos_dtInv(jat)
     dsi1t =   dpsi1(xt)
     dsi2t =   dpsi2(xt)*eos_dt(jat)

     dsi0mt = -dpsi0(mxt)*eos_dtInv(jat)
     dsi1mt =  dpsi1(mxt)
     dsi2mt = -dpsi2(mxt)*eos_dt(jat)

     dsi0d =   dpsi0(xd)*eos_ddInv(iat)
     dsi1d =   dpsi1(xd)
     dsi2d =   dpsi2(xd)*eos_dd(iat)

     dsi0md = -dpsi0(mxd)*eos_ddInv(iat)
     dsi1md =  dpsi1(mxd)
     dsi2md = -dpsi2(mxd)*eos_dd(iat)


     !!  the second derivatives of the basis functions
     ddsi0t =   ddpsi0(xt)*eos_dtSqrInv(jat)
     ddsi1t =   ddpsi1(xt)*eos_dtInv(jat)
     ddsi2t =   ddpsi2(xt)

     ddsi0mt =  ddpsi0(mxt)*eos_dtSqrInv(jat)
     ddsi1mt = -ddpsi1(mxt)*eos_dtInv(jat)
     ddsi2mt =  ddpsi2(mxt)


     !!  the free energy
     free  = h5( & 
          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, & 
          si0d,   si1d,   si2d,   si0md,   si1md,   si2md)


     !!  derivative with respect to density
     df_d  = h5( & 
          si0t,   si1t,   si2t,   si0mt,   si1mt,   si2mt, & 
          dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)


     !!  derivative with respect to temperature
     df_t = h5( & 
          dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, & 
          si0d,   si1d,   si2d,   si0md,   si1md,   si2md)


     !!  second derivative with respect to temperature
     df_tt = h5( & 
          ddsi0t, ddsi1t, ddsi2t, ddsi0mt, ddsi1mt, ddsi2mt, & 
          si0d,   si1d,   si2d,   si0md,   si1md,   si2md)


     !!  second derivative with respect to temperature and density
     df_dt = h5( & 
          dsi0t,  dsi1t,  dsi2t,  dsi0mt,  dsi1mt,  dsi2mt, & 
          dsi0d,  dsi1d,  dsi2d,  dsi0md,  dsi1md,  dsi2md)




     !!  now get the pressure derivative with density, chemical potential, and 
     !!  electron positron number densities
     !!  get the interpolation weight functions
     si0t   =  xpsi0(xt)
     si1t   =  xpsi1(xt)*eos_dt(jat)

     si0mt  =  xpsi0(mxt)
     si1mt  =  -xpsi1(mxt)*eos_dt(jat)

     si0d   =  xpsi0(xd)
     si1d   =  xpsi1(xd)*eos_dd(iat)

     si0md  =  xpsi0(mxd)
     si1md  =  -xpsi1(mxd)*eos_dd(iat)


     !!  pressure derivative with density
     dpepdd  = h3dpd(iat,jat, & 
          si0t,   si1t,   si0mt,   si1mt, & 
          si0d,   si1d,   si0md,   si1md)
     dpepdd  = max(ye * dpepdd,0.0e0)



     !!  electron chemical potential etaele
     etaele  = h3e(iat,jat, & 
                    si0t,   si1t,   si0mt,   si1mt, & 
                    si0d,   si1d,   si0md,   si1md)


     !!  electron + positron number densities
     xnefer   = h3x(iat,jat, & 
                  si0t,   si1t,   si0mt,   si1mt, & 
                  si0d,   si1d,   si0md,   si1md)


     !!  the desired electron-positron thermodynamic quantities
     x3      = din * din
     pele    = x3 * df_d
     dpepdt  = x3 * df_dt

     sele    = -df_t * ye
     dsepdt  = -df_tt * ye
     dsepdd  = -df_dt * ye * ye

 
     eele    = ye * free + btemp * sele
     deepdt  = btemp * dsepdt
     deepdd  = ye*ye*df_d + btemp*dsepdd

     if (bAprox13t) then
        dsepda  = ytot1 * (ye * df_dt * din - sele)    !Calhoun
        dsepdz  = -ytot1 * (ye * df_dt * den  + df_t)  !Calhoun
        deepda  = -ye * ytot1 * (free +  df_d * din) + btemp * dsepda  !Calhoun
        deepdz  = ytot1* (free + ye * df_d * den) + btemp * dsepdz !Calhoun
     end if


     !!  coulomb section:
     !!  initialize
     pcoul    = 0.0e0
     dpcouldd = 0.0e0
     dpcouldt = 0.0e0
     ecoul    = 0.0e0
     decouldd = 0.0e0
     decouldt = 0.0e0
     decoulda = 0.0e0  !Calhoun
     decouldz = 0.0e0  !Calhoun
     scoul    = 0.0e0
     dscouldd = 0.0e0
     dscouldt = 0.0e0

     !! Set the coulomb multiplier to a local value -- we might change it only within this call
     local_coulombMult = eos_coulombMult



     !!  uniform background corrections & only the needed parts for speed
     !!  plasg is the plasma coupling parameter
     !!  split up calculations below -- they all used to depend upon a redefined z
     z1        = forth * pi
     s1        = z1 * xni
     dsdd      = z1 * dxnidd
     lami      = 1.0e0/s1**third
     inv_lami  = 1.0e0/lami
     z2        = -third * lami/s1
     lamidd    = z2 * dsdd


     plasg     = zbar*zbar*esqu*ktinv*inv_lami
     plasg_inv = 1.0e0/plasg  
     z3        = -plasg * inv_lami
     plasgdd   = z3 * lamidd
     plasgdt   = -plasg * ktinv * kerg

     if (bAprox13t) then
        dsda     = z1 * dxnida !Calhoun
        lamida   = z2 * dsda / s1  ! Calhoun
        plasgda  = z3 * lamida          !Calhoun
        plasgdz  = 2.0d0 * plasg/zbar  !Calhoun
     end if


     !!  yakovlev & shalybkov 1989 equations 82, 85, 86, 87
     if (plasg .ge. 1.0) then
        x4       = plasg**(0.25e0)
        z4       = c1/x4
        ecoul    = dxnidd * kt * (a1*plasg + b1*x4 + z4 + d1cc)
        pcoul    = third * den * ecoul
        scoul    = -kavoy*(3.0e0*b1*x4 - 5.0e0*z4 &
             + d1cc*(log(plasg) - 1.0e0) - e1cc)

        y1       = dxnidd * kt * (a1 + 0.25e0*plasg_inv*(b1*x4 - z4))
        decouldd = y1 * plasgdd 
        decouldt = y1 * plasgdt + ecoul * tempi
        dpcouldd = third * (ecoul + den * decouldd)
        dpcouldt = third * den  * decouldt

        if (bAprox13t) then
           decoulda = y1 * plasgda - ecoul/abar   !Calhoun
           decouldz = y1 * plasgdz                !Calhoun
        end if

        y2       = -kavoy*plasg_inv*(0.75e0*b1*x4 + 1.25e0*z4 + d1cc)
        dscouldd = y2 * plasgdd
        dscouldt = y2 * plasgdt



        !!  yakovlev & shalybkov 1989 equations 102, 103, 104
     else if (plasg .lt. 1.0) then
        x5       = plasg * sqrt(plasg)
        y3       = plasg**b2
        z5       = c2 * x5 - third * a2 * y3
        pcoul    = -pion * z5
        ecoul    = 3.0e0 * pcoul * deni
        scoul    = -kavoy*(c2*x5 - a2*(b2 - 1.0e0)/b2*y3)

        s2       = (1.5e0*c2*x5 - third*a2*b2*y3)*plasg_inv
        dpcouldd = -dpiondd*z5 - pion*s2*plasgdd
        dpcouldt = -dpiondt*z5 - pion*s2*plasgdt
        decouldd = 3.0e0*dpcouldd*deni - ecoul*deni
        decouldt = 3.0e0*dpcouldt*deni

        !! Equations for decoulda and decouldz can be added here....

        if (bAprox13t) then
           call Driver_abortFlash('[eos_helm] TRAGEDY decoulda and decouldz not defined!')
        end if

        s3       = -kavoy*plasg_inv*(1.5e0*c2*x5 - a2*(b2 - 1.0e0)*y3)
        dscouldd = s3 * plasgdd
        dscouldt = s3 * plasgdt
     end if

     s4 = prad + pion + pele
     x6 = s4 + pcoul*eos_coulombMult

     ! assume that NaN always compares as false in an inequality
     if ( .not. (x6 > 0.e0) ) then

        call Logfile_stampMessage('[eos_helm] Negative total pressure.')
        print *,'[eos_helm] Negative total pressure.'
        ! print some info to help operator figure out what is wrong
        ! this goes to stdout as well because DMT is unsure of the stampMessage infrastructure
        write(internalFile,*)' values: dens,temp: ',den,btemp
        call Logfile_stampMessage(internalFile)
        print *, internalFile
        write(internalFile,*)' values: abar,zbar: ',abar,zbar
        call Logfile_stampMessage(internalFile)
        print *, internalFile
        write(internalFile,*)' coulomb coupling parameter Gamma: ',plasg
        call Logfile_stampMessage(internalFile)
        print *, internalFile

        if ( .not. (abar > 0.e0) ) then
           write(internalFile,*) '  However, abar is negative, abar=',abar
           call Logfile_stampMessage(internalFile)
           print *, internalFile
           call Logfile_stampMessage('  It is possible that the mesh is of low quality.')
           print *, '      It is possible that the mesh is of low quality.'
           ! always abort if abar is negative as that is physically invalid
           call Driver_abortFlash('[eos_helm] ERROR: abar is negative.')
        endif

        if ( s4 > 0.e0 ) then
           write(internalFile,'(a,3es12.3)')  &
     &           ' nonpositive P caused by coulomb correction: Pnocoul,Pwithcoul: ',s4,x6
           call Logfile_stampMessage(internalFile)
           print *, internalFile

           if ( eos_coulombMult > 0.e0 ) then
              call Logfile_stampMessage(&
                       '  set runtime parameter eos_coulombMult to zero if plasma Coulomb corrections not important')
              print *, '  set runtime parameter eos_coulombMult to zero if plasma Coulomb corrections not important'
           end if
           if ( eos_coulombAbort) then
              call Driver_abortFlash('[eos_helm] ERROR: coulomb correction causing negative total pressure.')
           else
              call Logfile_stampMessage('Setting coulombMult to zero for this call, eos_coulombAbort=false.')
              print *, '  Setting coulombMult to zero for this call, eos_coulombAbort=false'
              local_coulombMult = 0.e0
           end if
        else
           write(internalFile,'(1p,4(a,e12.5))') &
                ' Prad  ',prad,   &
                ' Pion ',pion,    &
                ' Pele  ',pele,                  &
                ' Pcoul ',pcoul*eos_coulombMult
           call Logfile_stampMessage(internalFile)
           print *, internalFile
           write(internalFile,'(1p,2(a,e12.5))') &
                ' Ptot   ',x6,                     &
                ' df_d   ',df_d
           call Logfile_stampMessage(internalFile)
           print *, internalFile

           call Driver_abortFlash('[eos_helm] ERROR: negative total pressure.')
        end if

     end if

     pcoul    = pcoul * local_coulombMult
     dpcouldd = dpcouldd * local_coulombMult
     dpcouldt = dpcouldt * local_coulombMult

     ecoul    = ecoul * local_coulombMult
     decouldd = decouldd * local_coulombMult
     decouldt = decouldt * local_coulombMult

     scoul    = scoul * local_coulombMult
     dscouldd = dscouldd * local_coulombMult 
     dscouldt = dscouldt * local_coulombMult

     !!  sum all the components
     pres    = prad    + pion    + pele   + pcoul
     ener    = erad    + eion    + eele   + ecoul
     entr    = srad    + sion    + sele   + scoul

#ifdef DEBUG_EOS
999  format(1x,1P,'Eos(',G12.3,',',G12.3,') ',A1,':','I:',G12.3,1x,'E:',G12.3,1x,'R:',G12.3,1x,'c:',G12.3)
     print 999,btemp,den,'E', eion,eele,erad,ecoul
     print 999,btemp,den,'P', pion,pele,prad,pcoul
#endif

     dpresdd = dpraddd + dpiondd + dpepdd + dpcouldd
     dpresdt = dpraddt + dpiondt + dpepdt + dpcouldt

     denerdd = deraddd + deiondd + deepdd + decouldd
     denerdt = deraddt + deiondt + deepdt + decouldt

     dentrdd = dsraddd + dsiondd + dsepdd + dscouldd
     dentrdt = dsraddt + dsiondt + dsepdt + dscouldt

     !!  form gamma_1
     presi = 1.0e0/pres
     chit  = btemp*presi * dpresdt
     chid  = dpresdd * den*presi
     x7     = pres * deni * chit/(btemp * denerdt)
     gamc  = chit*x7 + chid
     cv    = denerdt
     cp    = cv*gamc/chid



     !!  store the output -- note that many of these are not used by the calling program!
     ptotRow(j)   = pres   !used by Eos as EOS_PRES = PRES_VAR
     etotRow(j)   = ener   !used by Eos as EOS_EINT = EINT_VAR
     stotRow(j)   = entr   !this is entropy, used by Eos as EOS_ENTR (masked)

     dpdRow(j)    = dpresdd  ! used as EOS_DPD
     dptRow(j)    = dpresdt  ! used as EOS_DPT ALWAYS used by MODE_DENS_PRES in Eos.F90

     dedRow(j)    = denerdd  ! used as EOS_DED
     detRow(j)    = denerdt  ! used as EOS_DET  ALWAYS used by MODE_DENS_EI in Eos.F90

     if (bAprox13t) then
        denerda = deradda + deionda + deepda + decoulda  !Calhoun
        denerdz = deraddz + deiondz + deepdz + decouldz  !Calhoun
        deaRow(j)    = denerda  !Calhoun EOS_DEA
        dezRow(j)    = denerdz  !Calhoun EOS_DEZ
     end if


     dsdRow(j)    = dentrdd  ! used as EOS_DSD  
     dstRow(j)    = dentrdt  ! used as EOS_DST  

     !UNUSED     pradRow(j)   = prad
     !UNUSED     eradRow(j)   = erad
     !UNUSED     sradRow(j)   = srad

     !UNUSED     pionRow(j)   = pion
     !UNUSED     eionRow(j)   = eion
     !UNUSED     sionRow(j)   = sion
 
     pelRow(j)   = pele     ! used as EOS_PEL
     !UNUSED     eeleRow(j)   = eele
     !UNUSED     seleRow(j)   = sele

     neRow(j)    = xnefer   ! used as EOS_NE  
     etaRow(j) = etaele     ! used as EOS_ETA 

     gamcRow(j)   = gamc    !used as EOS_GAMC = GAMC_VAR

     cvRow(j)     = cv      ! EOS_CV
     cpRow(j)     = cp      ! EOS_CP

     !!  end of vectorization loop
  enddo

!   call Timers_stop("eos_helm")

  return

end subroutine eos_helm
