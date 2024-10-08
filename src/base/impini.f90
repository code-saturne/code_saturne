!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

subroutine impini () &
  bind(C, name='cs_f_impini')
!================


!===============================================================================
! Purpose:
!  ---------

! Print computation parameters after user changes in cs_user_parameters.cpp

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use radiat
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments


! Local variables

integer          igg, ige

!===============================================================================
! 1. Introduction
!===============================================================================

write(nfecra,1000)

if (ippmod(icod3p).ne.-1) then
  write(nfecra,1010)
  write(nfecra,1020) ippmod(icod3p)
  write(nfecra,1070) namgas, pcigas
  write(nfecra,1080) trim(nomcog(igfuel(1))), &
  -nreact(igoxy(1)), trim(nomcog(igoxy(1))), trim(nomcog(igprod(1)))
  write(nfecra,'(a20,10(1x,a14))') "Mass composition", "Fuel", "Oxydizer", "Products"
  write(nfecra,'(a20,10(1x,a14))') "----------------", "----", "--------", "--------"
  do ige = 1, ngaze
    write(nfecra,'(a15,10(1x,f14.5))') trim(nomcoe(ige)), (coefeg(ige,igg), igg=1, ngazg)
  enddo
  write(nfecra,1000)
  write(nfecra,'(a20,10(1x,a14))') "Molar composition", "Fuel", "Oxydizer", "Products"
  write(nfecra,'(a20,10(1x,a14))') "-----------------", "----", "--------", "--------"
  do ige = 1, ngaze
    write(nfecra,'(a15,10(1x,f14.5))') trim(nomcoe(ige)), (compog(ige,igg), igg=1, ngazg)
  enddo
else if (ippmod(islfm).ne.-1) then
  write(nfecra,1010)
  write(nfecra,1040) ippmod(islfm)
else if (ippmod(icoebu).ne.-1) then
  write(nfecra,1010)
  write(nfecra,1030) ippmod(icoebu), cebu
endif


 1000 format(                                                     &
                                                                /,&
' ===========================================================', /,&
                                                                /,&
'               CALCULATION PARAMETERS SUMMARY',                /,&
'               ==============================',                /,&
                                                                /,&
' -----------------------------------------------------------', /)

 9900 format(                                                     &
                                                                /,&
' -----------------------------------------------------------', /)
 1010 format(                                                     &
                                                                /,&
' ** SPECIFIC PHYSICS:',                                        /,&
'    ----------------',                                         /)
 1020 format(                                                     &
' --- Diffusion Flame: 3 Point Chemistry',                      /,&
'       OPTION = ',4x,i10                                       /)
 1030 format(                                                     &
' --- Premixed Flame: EBU Model',                               /,&
'       OPTION = ',4x,i10,                                      /,&
'       CEBU   = ',e14.5                                        /)
 1040 format(                                                     &
' --- Diffusion Flame: Steady laminar flamelet model',          /,&
'       OPTION = ',4x,i10,                                      /)
 1070 format(                                                     &
' --- Combustible characteristics',                             /,&
'       Combustible : ',4x,a,                                   /,&
'       PCI = ',4x,e14.5,  ' J/kg',                             /)
 1080 format(                                                     &
" --- Chemical reaction: ",                                     /,&
"       ", a," + ",f6.3," (",a,") --> ",a,                      /)

!===============================================================================
! DEFINITION GENERALE DU CAS
!===============================================================================

! --- Dimensions

write(nfecra,1500)
write(nfecra,1520) nvar,nscal,nscaus,nscapp

write(nfecra,9900)

 1500 format(                                                     &
                                                                /,&
' ** DIMENSIONS',                                               /,&
'    ----------',                                               /)
 1520 format(                                                     &
' --- Physics',                                                 /,&
'       NVAR   = ',4x,i10,    ' (Nb variables                )',/,&
'       NSCAL  = ',4x,i10,    ' (Nb scalars                  )',/,&
'       NSCAUS = ',4x,i10,    ' (Nb user scalars             )',/,&
'       NSCAPP = ',4x,i10,    ' (Nb specific physics scalars )',/)

!===============================================================================
! DISCRETISATION DES EQUATIONS
!===============================================================================

! --- Marche en temps

if (idtvar.ge.0) then

  ! Coefficient de relaxation de la masse volumique

  if (ippmod(icod3p).ge.0 .or. ippmod(islfm).ge.0 .or. &
      ippmod(icoebu).ge.0 .or. ippmod(icolwc).ge.0 .or. &
      ippmod(iccoal).ge.0) then
    write(nfecra,3000)
    write(nfecra,3050) srrom
    write(nfecra,9900)

  endif

endif

 3000 format(                                                     &
                                                                /,&
' ** TIME STEPPING',                                            /,&
'    -------------',                                            /)
 3050 format(                                                     &
'--- Relaxation coefficient',                                   /,&
'    RHO(n+1)=SRROM*RHO(n)+(1-SRROM)*RHO(n+1)',                 /,&
'       SRROM  = ',e14.5,                                       /)

return
end subroutine
