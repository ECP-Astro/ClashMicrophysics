!!****if* source/physics/Eos/EosMain/Helmholtz/SpeciesBased/eos_initHelmholtz
!!
!! NAME
!!
!!  eos_initHelmholtz
!!
!! SYNOPSIS
!!
!!  call eos_initHelmholtz()
!!
!! DESCRIPTION
!!
!!  Initialize the Helmholtz EOS.  The table data is read in on processor 0
!!  and broadcast to all other processors at the start of FLASH execution.
!!  This routine first checks for a binary copy of the table (helm_table.bdat), 
!!  and then for the ASCII version (helm_table.dat).  If the binary table 
!!  file was not found, it is created by this routine for subsequent use.
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!
!!   eos_coulombMult[Real, default 1.0] -- Coulomb correction multiplier. Set
!!               to zero to ignore the Coloumb correction.
!!   eos_tolerance[Real, 1.0e-8] -- Convergence tolerance for the Newton-Rhapson
!!               iterations
!!   eos_maxNewton[Integer, 50] -- Maximum number of Newton-Raphson iterations
!!   eos_forceConstantInput     -- This switch forces the Eos implementation
!!                                 to never modify EINT or PRES in MODE_DENS_EI
!!                                 and MODE_DENS_PRES. If this is .false. (the
!!                                 default), calls to Eos may slightly modify
!!                                 these input variables in order to preserve
!!                                 thermodynamic equilibrium.
!!
!!  NOTES
!!
!!  Helmholtz law Eos defines two mesh-based parameters GAMC_VAR and GAME_VAR in Flash.h
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_EOS
#endif

subroutine eos_initHelmholtz()

  use Eos_data, ONLY : eos_type, eos_meshMe, &
       eos_eintSwitch, eos_smallt
  use eos_helmData 
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use actual_eos_module
  use actual_network
  implicit none

  ! vector_eos.fh computes the vector length from nxb, nyb, nzb, so 
  ! this information must be provided

#include "Flash.h"
#include "constants.h"
#include "Eos.h"
  include 'Flash_mpi.h'
  integer:: unitEos =2
  integer :: i, j
  real(8) :: tstp, dstp
  integer :: istat, ierr

  double precision :: dth, dt2, dti, dt2i
  double precision :: dd, dd2, ddi, dd2i
  double precision :: tsav, dsav

  !get the runtime parameters

  call RuntimeParameters_get('smallt', eos_smallt)
  call RuntimeParameters_get('eos_tolerance', eos_tol)
  call RuntimeParameters_get('eos_maxNewton', eos_maxNewton)
  call RuntimeParameters_get('eos_coulombMult', eos_coulombMult)

#ifdef DEBUG_EOS
  print *, 'in eos_initHelmholtz'
#endif
  call RuntimeParameters_get('eos_coulombAbort', eos_coulombAbort)
#ifdef DEBUG_EOS
  print *, 'done with RuntimeParameters_get (eos_coulombAbort)'
#endif

#ifndef EINT_VAR
  if (eos_eintSwitch > 0.0) then
     call Driver_abortFlash("[Eos_init] eintSwitch is nonzero, but EINT_VAR not defined!")
  end if
#endif

#ifdef USE_EOS_YE
  write(*,*)"USE_EOS_YE should not be defined with Helmholtz/SpeciesBased EOS.  Use Helmholtz/Ye instead!"
  call Driver_abortFlash("[Eos_init] Use Helmholtz/Ye with USE_EOS_YE mode")
#endif

  call RuntimeParameters_get("eos_forceConstantInput",eos_forceConstantInput)




  





  if (eos_meshMe==MASTER_PE) then


     !! We will need to change the values of EOSIJMAX in order to use OpenACC version of helm_table.dat

     open(unit=unitEos,file='helm_table.bdat',status='old',iostat=istat)
     close(unit=unitEos)
     print*,'about to open file'
     istat=1
     if (istat.ne.0) then
        write(*,*) '[Eos_init] Cannot open helm_table.bdat!'
        write(*,*) '[Eos_init] Trying old helm_table.dat!'
        
        open(unit=unitEos,file='helm_table.dat',status='old',iostat=istat)

        if (istat .ne. 0) then
           write(*,*) '[Eos_init] ERROR: opening helm_table.dat!'
           call Driver_abortFlash("[Eos_init] ERROR: opening helm_table.dat")
        endif

        !..read the helmholtz free energy table
        do j=1,EOSJMAX
           do i=1,EOSIMAX
              read(unitEos,*) eos_f(i,j),eos_fd(i,j),eos_ft(i,j),&
                   eos_fdd(i,j),eos_ftt(i,j), & 
                   eos_fdt(i,j), eos_fddt(i,j),eos_fdtt(i,j),eos_fddtt(i,j)
           enddo
        enddo

        !..read the pressure derivative with density table
        do j=1,EOSJMAX
           do i=1,EOSIMAX
              read(unitEos,*) eos_dpdf(i,j),eos_dpdfd(i,j),&
                   eos_dpdft(i,j),eos_dpdfdt(i,j)
           enddo
        enddo

        !..read the electron chemical potential table
        do j=1,EOSJMAX
           do i=1,EOSIMAX
              read(unitEos,*) eos_ef(i,j),eos_efd(i,j),&
                   eos_eft(i,j),eos_efdt(i,j)
           enddo
        enddo

        !..read the number density table
        do j=1,EOSJMAX
           do i=1,EOSIMAX
              read(unitEos,*) eos_xf(i,j),eos_xfd(i,j),&
                   eos_xft(i,j),eos_xfdt(i,j)
           enddo
        enddo

        !..close up the data file and write a message
        close(unitEos)

        !..dump binary version of table for later use
        istat = EOSIMAX*EOSJMAX
        call eos_writeHfet(istat, & 
             eos_f,eos_fd,eos_ft,eos_fdd,&
             eos_ftt,eos_fdt,eos_fddt,eos_fdtt,eos_fddtt, & 
             eos_dpdf,eos_dpdfd,eos_dpdft,eos_dpdfdt, & 
             eos_ef,eos_efd,eos_eft,eos_efdt, & 
             eos_xf,eos_xfd,eos_xft,eos_xfdt)

        !..read binary version of table
     else
        istat = EOSIMAX*EOSJMAX
        call eos_readHfet(istat, & 
             eos_f,eos_fd,eos_ft,eos_fdd,&
             eos_ftt,eos_fdt,eos_fddt,eos_fdtt,eos_fddtt, & 
             eos_dpdf,eos_dpdfd,eos_dpdft,eos_dpdfdt, & 
             eos_ef,eos_efd,eos_eft,eos_efdt, & 
             eos_xf,eos_xfd,eos_xft,eos_xfdt)
     endif
  endif

  !..broadcast to rest of processors
  istat = EOSIMAX*EOSJMAX
  call MPI_BCAST(eos_f,      istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_fd,     istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_ft,     istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_fdd,    istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_ftt,    istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_fdt,    istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_fddt,   istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_fdtt,   istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_fddtt,  istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_dpdf,   istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_dpdfd,  istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_dpdft,  istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_dpdfdt, istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_ef,     istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_efd,    istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_eft,    istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_efdt,   istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_xf,     istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)

  call MPI_BCAST(eos_xfd,    istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_xft,    istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(eos_xfdt,   istat, FLASH_REAL, MASTER_PE, MPI_COMM_WORLD, ierr)


  !! After broadcast to all MPI ranks, we should NOW make copies into variables
  !! that will be used by "shared code".

  do j = 1, EOSJMAX
    do i = 1, EOSIMAX

      !! Make copies of helmholtz free energy table
      f(i,j)     = eos_f(i,j)
      fd(i,j)    = eos_fd(i,j)
      ft(i,j)    = eos_ft(i,j)
      fdd(i,j)   = eos_fdd(i,j)
      ftt(i,j)   = eos_ftt(i,j)
      fdt(i,j)   = eos_fdt(i,j)
      fddt(i,j)  = eos_fddt(i,j)
      fdtt(i,j)  = eos_fdtt(i,j)
      fddtt(i,j) = eos_fddtt(i,j)

      !! Make copies of pressure derivative with density table
      dpdf(i,j)   = eos_dpdf(i,j)
      dpdfd(i,j)  = eos_dpdfd(i,j)
      dpdft(i,j)  = eos_dpdft(i,j)
      dpdfdt(i,j) = eos_dpdfdt(i,j)

      !! Make copies of the electron chemical potential table
      ef(i,j)   = eos_ef(i,j)
      efd(i,j)  = eos_efd(i,j)
      eft(i,j)  = eos_eft(i,j)
      efdt(i,j) = eos_efdt(i,j)

      !! Make copies of the number density table
      xf(i,j)   = eos_xf(i,j)
      xfd(i,j)  = eos_xfd(i,j)
      xft(i,j)  = eos_xft(i,j)
      xfdt(i,j) = eos_xfdt(i,j)

    end do
  end do


  !! Make copies of helmholtz free energy table
!  f     = eos_f
!  fd    = eos_fd
!  ft    = eos_ft
!  fdd   = eos_fdd
!  ftt   = eos_ftt
!  fdt   = eos_fdt
!  fddt  = eos_fddt
!  fdtt  = eos_fdtt
!  fddtt = eos_fddtt 

  !! Make copies of pressure derivative with density table
!  dpdf   = eos_dpdf
!  dpdfd  = eos_dpdfd
!  dpdft  = eos_dpdft
!  dpdfdt = eos_dpdfdt

  !! Make copies of the electron chemical potential table
!  ef   = eos_ef
!  efd  = eos_efd
!  eft  = eos_eft
!  efdt = eos_efdt

  !! Make copies of the number density table
!  xf   = eos_xf
!  xfd  = eos_xfd
!  xft  = eos_xft
!  xfdt = eos_xfdt


  !! Most of this is probably not needed, but we will be using our
  !! own version regardless.
  eos_tlo   = 4.0e0
  tstp  = (11.0e0 - eos_tlo)/float(EOSJMAX-1)
  eos_tstpi = 1.0e0/tstp
  eos_dlo   = -10.0e0
  dstp  = (11.0e0 - eos_dlo)/float(EOSIMAX-1)
  eos_dstpi = 1.0e0/dstp
  do j=1,EOSJMAX
     eos_t(j) = 10.0e0**(eos_tlo + (j-1)*tstp)
     do i=1,EOSIMAX
        eos_d(i) = 10.0e0**(eos_dlo + (i-1)*dstp)
     enddo
  enddo


  !! This is new version used in OpenACC. We need to make sure to 
  !! alter some variable names so they are distinct from FLASH ones.

  itmax = imax
  jtmax = jmax

!  tlo   = 3.0d0
  tlo   = 4.0e0
!  thi   = 13.0d0
  thi   = 11.0e0
  tstpp  = (thi - tlo)/float(jmax-1)
  tstpi = 1.0d0/tstpp
!  dlo   = -12.0d0
  dlo   = -10.0e0
!  dhi   = 15.0d0
  dhi   = 11.0e0
  dstpp  = (dhi - dlo)/float(imax-1)
  dstpi = 1.0d0/dstpp
  do j=1,jmax
     tsav = tlo + (j-1)*tstpp
     t(j) = 10.0d0**(tsav)
     do i=1,imax
        dsav = dlo + (i-1)*dstpp
        d(i) = 10.0d0**(dsav)
     end do
  end do


  !! If we want to run some parts of the code with the original Helm EOS,
  !! we might need to leave some of this intact.

  !..store the temperature and density differences and their inverses 
  do j=1,EOSJMAX-1
     eos_dt(j)   = eos_t(j+1) - eos_t(j)
     eos_dtSqr(j)  = eos_dt(j)*eos_dt(j)
     eos_dtInv(j)  = 1.0e0/eos_dt(j)
     eos_dtSqrInv(j) = 1.0e0/eos_dtSqr(j)
  enddo
  do i=1,EOSIMAX-1
     eos_dd(i)   = eos_d(i+1) - eos_d(i)
     eos_ddSqr(i)  = eos_dd(i)*eos_dd(i)
     eos_ddInv(i)  = 1.0e0/eos_dd(i)
     eos_ddSqrInv(i) = 1.0e0/eos_ddSqr(i)
  enddo
  eos_type=EOS_HLM

  !..   construct the temperature and density deltas and their inverses
  do j = 1, jmax-1
     dth         = t(j+1) - t(j)
     dt2         = dth * dth
     dti         = 1.0d0/dth
     dt2i        = 1.0d0/dt2
     dt_sav(j)   = dth
     dt2_sav(j)  = dt2
     dti_sav(j)  = dti
     dt2i_sav(j) = dt2i
  end do
  do i = 1, imax-1
     dd          = d(i+1) - d(i)
     dd2         = dd * dd
     ddi         = 1.0d0/dd
     dd2i        = 1.0d0/dd2
     dd_sav(i)   = dd
     dd2_sav(i)  = dd2
     ddi_sav(i)  = ddi
     dd2i_sav(i) = dd2i
  end do


  ! Some initialization of constants
  esqu = qe * qe

  a2rad   = pi/180.0d0
  rad2a   = 180.0d0/pi

  hbar    = 0.5d0 * h/pi
  kev     = kerg/ev2erg_eos
  rbohr   = hbar*hbar/(me_eos * qe * qe)
  fine    = qe*qe/(hbar*clight)

  asol    = 4.0d0 * ssol / clight
  weinlam = h*clight/(kerg * 4.965114232d0)
  weinfre = 2.821439372d0*kerg/h
  pc      = 3.261633d0 * ly

  sioncon = (2.0d0 * pi * amu * kerg)/(h*h)
  forth   = 4.0d0/3.0d0
  forpi   = 4.0d0 * pi
  kergavo = kerg * avo_eos
  ikavo   = 1.0d0/kergavo
  asoli3  = asol/3.0d0
  light2  = clight * clight


  ! Set up the minimum and maximum possible densities.

  mintemp = 10.d0**tlo
  maxtemp = 10.d0**thi
  mindens = 10.d0**dlo
  maxdens = 10.d0**dhi

  input_is_constant = .true.
  do_coulomb        = .true.

  !$acc update &
  !$acc device(pi, a2rad, rad2a) &
  !$acc device(h, hbar, qe, avo_eos, clight, kerg) &
  !$acc device(ev2erg_eos, kev, amu, me_eos) &
  !$acc device(rbohr, fine) &
  !$acc device(ssol, asol, weinlam, weinfre) &
  !$acc device(ly, pc) &
  !$acc device(sioncon, forth, forpi, kergavo, ikavo, asoli3, light2) &
  !$acc device(a1, b1, c1, d1, e1, a2, b2, c2, onethird, esqu) &
  !$acc device(ttol, dtol, tlo, thi, dlo, dhi) &
  !$acc device(tstpp, tstpi, dstpp, dstpi) &
  !$acc device(itmax, jtmax, d, t) &
  !$acc device(f, fd, ft, fdd, ftt, fdt, fddt, fdtt, fddtt) &
  !$acc device(dpdf, dpdfd, dpdft, dpdfdt) &
  !$acc device(ef, efd, eft, efdt, xf, xfd, xft, xfdt) &
  !$acc device(dt_sav, dt2_sav, dti_sav, dt2i_sav) &
  !$acc device(dd_sav, dd2_sav, ddi_sav, dd2i_sav), &
  !$acc device(do_coulomb, input_is_constant) &
  !$acc device(max_newton)

  call actual_network_init()

!  allocate(eos_state(16))

end subroutine eos_initHelmholtz
