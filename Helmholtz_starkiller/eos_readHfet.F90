!!****if* source/physics/Eos/EosMain/Helmholtz/eos_readHfet
!!
!! NAME
!!
!!  eos_readHfet
!!
!! SYNOPSIS
!!
!!  eos_readHfet(integer(IN) :: n, 
!!            real(OUT) :: f(n), 
!!            real(OUT) :: fd(n), 
!!            real(OUT) :: ft(n), 
!!            real(OUT) :: fdd(n), 
!!            real(OUT) :: ftt(n),
!!            real(OUT) :: fdt(n), 
!!            real(OUT) :: fddt(n), 
!!            real(OUT) :: fdtt(n), 
!!            real(OUT) :: fddtt(n), 
!!            real(OUT) :: dpdf(n), 
!!            real(OUT) :: dpdfd(n),
!!            real(OUT) :: dpdft(n), 
!!            real(OUT) :: dpdfdt(n), 
!!            real(OUT) :: ef(n), 
!!            real(OUT) :: efd(n), 
!!            real(OUT) :: eft(n), 
!!            real(OUT) :: efdt(n), 
!!            real(OUT) :: xf(n), 
!!            real(OUT) :: xfd(n), 
!!            real(OUT) :: xft(n), 
!!            real(OUT) :: xfdt(n)  )
!!
!! DESCRIPTION
!!
!!  Read binary data file containing coefficients for tabular helmholtz
!!
!! ARGUMENTS
!!
!!     n -- number of variables in the arrays
!!     f --  Helmholtz free energy
!!     fd --  derivative of f wrt density
!!     ft --  derivative of f wrt temperature
!!     fdd --  second derivative of f wrt density
!!     ftt --  second derivative of f wrt temperature
!!     fdt --  second derivative of f wrt density and temperature
!!     fddt --  third derivative of f wrt density^2 and temperature
!!     fdtt --  third derivative of f wrt density and temperature^2 e.g. dF/(dd)(dt^2)
!!     fddtt --  fourth derivative of f wrt density^2 and temperature^2
!!     dpdf --  pressure derivative 
!!     dpdfd -- 
!!     dpdft --  
!!     dpdfdt --  
!!     ef --  electron chemical potential
!!     efd --  
!!     eft --  
!!     efdt --  
!!     xf --  number density
!!     xfd --  
!!     xft --  
!!     xfdt --
!!
!!  NOTE
!! 
!!    See Timmes and Swesty, 2000, AJSS, "The Accuracy, Consistency, and Speed of an Electron-Positron
!!    Equation of State Based on Table Interpolation of the Helmholtz Free Energy"
!!
!!***

subroutine eos_readHFet(n, f, fd, ft, fdd, ftt,  &
     &     fdt, fddt, fdtt, fddtt, dpdf, dpdfd,  &
     &   dpdft, dpdfdt, ef, efd, eft, efdt,      &
     &   xf, xfd, xft, xfdt)
   use Driver_interface, ONLY : Driver_abortFlash
   
   implicit none

#include "Flash.h"

  integer, intent(in) :: n
  real, intent(out) :: f(n),fd(n),ft(n),fdd(n),ftt(n),fdt(n),fddt(n),fdtt(n),fddtt(n), &
       &         dpdf(n),dpdfd(n),dpdft(n),dpdfdt(n), &
       &         ef(n),efd(n),eft(n),efdt(n),xf(n),xfd(n),xft(n),xfdt(n)

  !! Local variables
  !!  fileUnit is a file variable
  integer, parameter ::  fileUnit = 36
  integer          :: numRead, ioStat



#ifdef DEBUG
  if (n .lt. 0) then
     call Driver_abortFlash("[eos_readHfet]  n must be positive")
  endif
#endif

  !! Open the file
  open (fileUnit,FILE='helm_table.bdat',ACTION='READ',STATUS='OLD',FORM='UNFORMATTED',IOSTAT=ioStat)
  if (ioStat .NE. 0)  call Driver_abortFlash("[eos_readHfet]  file open failure!")

  read(fileUnit,END=101) f
  read(fileUnit,END=102) fd
  read(fileUnit,END=103) ft
  read(fileUnit,END=104) fdd
  read(fileUnit,END=105) ftt
  read(fileUnit,END=106) fdt
  read(fileUnit,END=107) fddt
  read(fileUnit,END=108) fdtt
  read(fileUnit,END=109) fddtt
  read(fileUnit,END=110) dpdf
  read(fileUnit,END=111) dpdfd
  read(fileUnit,END=112) dpdft
  read(fileUnit,END=113) dpdfdt
  read(fileUnit,END=114) ef
  read(fileUnit,END=115) efd
  read(fileUnit,END=116) eft
  read(fileUnit,END=117) efdt
  read(fileUnit,END=118) xf
  read(fileUnit,END=119) xfd
  read(fileUnit,END=120) xft
  read(fileUnit,END=121) xfdt


  !! close up and return
  close(fileUnit,IOSTAT=ioStat)
  if (ioStat .NE. 0)  call Driver_abortFlash("[eos_readHfet]  couldn't close file!")

  return

  !! Error messages on insufficient data

101 call Driver_abortFlash("[eos_readHfet]  failed read on f!")
102 call Driver_abortFlash("[eos_readHfet]  failed read on fd!")
103 call Driver_abortFlash("[eos_readHfet]  failed read on ft!")
104 call Driver_abortFlash("[eos_readHfet]  failed read on fdd!")
105 call Driver_abortFlash("[eos_readHfet]  failed read on ftt!")
106 call Driver_abortFlash("[eos_readHfet]  failed read on fdt!")
107 call Driver_abortFlash("[eos_readHfet]  failed read on fddt!")
108 call Driver_abortFlash("[eos_readHfet]  failed read on fdtt!")
109 call Driver_abortFlash("[eos_readHfet]  failed read on fddtt!")
110 call Driver_abortFlash("[eos_readHfet]  failed read on dpdf!")
111 call Driver_abortFlash("[eos_readHfet]  failed read on dpdfd!")
112 call Driver_abortFlash("[eos_readHfet]  failed read on dpdft!")
113 call Driver_abortFlash("[eos_readHfet]  failed read on dpdfdt!")
114 call Driver_abortFlash("[eos_readHfet]  failed read on ef!")
115 call Driver_abortFlash("[eos_readHfet]  failed read on efd!")
116 call Driver_abortFlash("[eos_readHfet]  failed read on eft!")
117 call Driver_abortFlash("[eos_readHfet]  failed read on efdt!")
118 call Driver_abortFlash("[eos_readHfet]  failed read on xf!")
119 call Driver_abortFlash("[eos_readHfet]  failed read on xfd!")
120 call Driver_abortFlash("[eos_readHfet]  failed read on xft!")
121 call Driver_abortFlash("[eos_readHfet]  failed read on xfdt!")



end subroutine eos_readHfet


!!#define SWAP(a,b)       temp=(a); (a)=(b); (b)=temp;


!!void byte_reverse (float *buf, int *nn){
!  int n;
!  char temp, *ptr;!


!#ifdef DEBUG
!  if (*nn<0)
!    call Driver_abortFlash("[temp] byte_reverse() :: n must be positive")
!#endif
  
!  for (ptr=(char *)buf,n=*nn; n--; ptr+=4)
!    {
!      SWAP(ptr[0],ptr[3])
!      SWAP(ptr[1],ptr[2])
!     }
!}
