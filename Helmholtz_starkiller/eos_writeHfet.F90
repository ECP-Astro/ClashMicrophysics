!!****if* source/physics/Eos/EosMain/Helmholtz/eos_writeHfet
!!
!! NAME
!!
!!  eos_writeHfet
!!
!! SYNOPSIS
!!
!!  eos_writeHfet(integer(IN) :: n, 
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

subroutine eos_writeHfet(n, f, fd, ft, fdd, ftt,  &
     &   fdt, fddt, fdtt, fddtt, dpdf, dpdfd,      &
     &   dpdft, dpdfdt, ef, efd, eft, efdt,        &
     &   xf, xfd, xft, xfdt)
   use Driver_interface, ONLY : Driver_abortFlash

   implicit none

#include "Flash.h"

   integer, intent(in) :: n
   real, intent(in) :: f(n),fd(n),ft(n),fdd(n),ftt(n),fdt(n),fddt(n),fdtt(n),fddtt(n), &
     &         dpdf(n),dpdfd(n),dpdft(n),dpdfdt(n), &
     &         ef(n),efd(n),eft(n),efdt(n),xf(n),xfd(n),xft(n),xfdt(n)

!! Local variables
!!  fileUnit is a file variable
   integer, parameter ::  fileUnit = 36
   integer          :: numRead, ioStat



#ifdef DEBUG
  if (n<0) then
    call Driver_abortFlash("[eos_writeHfet]  n must be positive")
  endif 
#endif

!! Open the file
  open (fileUnit,FILE='helm_table.bdat',ACTION='WRITE',STATUS='REPLACE',FORM='UNFORMATTED',IOSTAT=ioStat)
    if (ioStat .NE. 0)     call Driver_abortFlash("[eos_writeHfet]  file open failure!")

!! Start writing

  write(fileUnit,ERR=101) f
  write(fileUnit,ERR=102) fd  
  write(fileUnit,ERR=103) ft  
  write(fileUnit,ERR=104) fdd  
  write(fileUnit,ERR=105) ftt  
  write(fileUnit,ERR=106) fdt  
  write(fileUnit,ERR=107) fddt  
  write(fileUnit,ERR=108) fdtt  
  write(fileUnit,ERR=109) fddtt
  write(fileUnit,ERR=110) dpdf  
  write(fileUnit,ERR=111) dpdfd  
  write(fileUnit,ERR=112) dpdft
  write(fileUnit,ERR=113) dpdfdt
  write(fileUnit,ERR=114) ef
  write(fileUnit,ERR=115) efd 
  write(fileUnit,ERR=116) eft  
  write(fileUnit,ERR=117) efdt
  write(fileUnit,ERR=118) xf
  write(fileUnit,ERR=119) xfd
  write(fileUnit,ERR=120) xft
  write(fileUnit,ERR=121) xfdt
  
!! close up and return
  close (fileUnit,IOSTAT=ioStat)
   if (ioStat .NE. 0)  call Driver_abortFlash("[eos_writeHfet]  couldn't close file!")
  
  return

!!  Abort statements

101    call Driver_abortFlash("[eos_writeHfet]  failed write on f!")
102    call Driver_abortFlash("[eos_writeHfet]  failed write on fd!")
103    call Driver_abortFlash("[eos_writeHfet]  failed write on ft!")
104    call Driver_abortFlash("[eos_writeHfet]  failed write on fdd!")
105    call Driver_abortFlash("[eos_writeHfet]  failed write on ftt!")
106    call Driver_abortFlash("[eos_writeHfet]  failed write on fdt!")
107    call Driver_abortFlash("[eos_writeHfet]  failed write on fddt!")
108    call Driver_abortFlash("[eos_writeHfet]  failed write on fdtt!")
109    call Driver_abortFlash("[eos_writeHfet]  failed write on fddtt!")
110    call Driver_abortFlash("[eos_writeHfet]  failed write on dpdf!")
111    call Driver_abortFlash("[eos_writeHfet]  failed write on dpdfd!")
112    call Driver_abortFlash("[eos_writeHfet]  failed write on dpdft!")
113    call Driver_abortFlash("[eos_writeHfet]  failed write on dpdfdt!")
114    call Driver_abortFlash("[eos_writeHfet]  failed write on ef!")
115    call Driver_abortFlash("[eos_writeHfet]  failed write on efd!")
116    call Driver_abortFlash("[eos_writeHfet]  failed write on eft!")
117    call Driver_abortFlash("[eos_writeHfet]  failed write on efdt!") 
118    call Driver_abortFlash("[eos_writeHfet]  failed write on xf!") 
119    call Driver_abortFlash("[eos_writeHfet]  failed write on xfd!") 
120    call Driver_abortFlash("[eos_writeHfet]  failed write on xft!") 
121    call Driver_abortFlash("[eos_writeHfet]  failed write on xfdt!") 

end subroutine eos_writeHfet


