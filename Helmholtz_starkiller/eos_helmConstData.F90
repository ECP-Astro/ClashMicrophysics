!!****if* source/physics/Eos/EosMain/Helmholtz/eos_helmConstData
!!
!! NAME
!!
!!  eos_helmConstData
!!
!!
!! SYNOPSIS
!!
!!  use eos_helmConstData
!!
!!
!! DESCRIPTION
!!
!! mathematical and physical constants (in cgs,except e0 which is in ev)
!!
!! the 1986 codta recommended valeus of the physical constants
!! by coehn & taylor 
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!
!!  NOTES
!!    See Timmes and Swesty, 2000, AJSS, "The Accuracy, Consistency, and Speed of an Electron-Positron
!!    Equation of State Based on Table Interpolation of the Helmholtz Free Energy"
!!
!!***


Module eos_helmConstData

  implicit none
  !
  !..math constants   
  real,parameter ::     pi       = 3.1415926535897932384e0, &
       eulercon = 0.577215664901532861e0
  real,parameter ::     a2rad    = pi/180.0e0,  rad2a = 180.0e0/pi
  !
  !..physical constants   
  real,parameter ::     g       = 6.67259e-8,                           &
       h       = 6.6260755e-27,                        &
       hbar    = 0.5 * h/pi,                           &
       e       = 4.8032068e-10,                        &
       avo     = 6.0221367e23,                         &
       c       = 2.99792458e10,                        &
       kerg    = 1.380658e-16,                         &
       kev     = 8.617385e-5,                          &
       amu     = 1.6605402e-24,                        &
       mn      = 1.6749286e-24,                        &
       mp      = 1.6726231e-24,                        &
       me      = 9.1093897e-28,                        &
       rbohr   = hbar*hbar/(me * e * e),               &
       fine    = e*e/(hbar*c),                         &
       hion    = 13.605698140e0,                       &
       ev2erg  = 1.602e-12,                             &
       ssol     = 5.67051e-5,                          &
       asol    = 4.0e0 * ssol / c,                     &
       weinlam = h*c/(kerg * 4.965114232e0),           &
       weinfre = 2.821439372e0*kerg/h,                 &
       rhonuc  = 2.342e14


  !..astronomical constants   
  real, parameter ::  msol    = 1.9892e33,                            &
       rsol    = 6.95997e10,                           &
       lsol    = 3.8268e33,                            &
       mearth  = 5.9764e27,                            &
       rearth  = 6.37e8,                               &
       ly      = 9.460528e17,                          &
       pc      = 3.261633e0 * ly,                      &
       au      = 1.495978921e13,                       &
       secyer  = 3.1558149984e7

  !.. derived constants
  real, parameter :: avoinv  = 1.0e0/avo, & 
       kergavo = kerg * avo, &
       asoli3  = asol/3.0e0, &
       sioncon = (2.0e0 * pi * amu * kerg)/(h*h)

end Module eos_helmConstData
