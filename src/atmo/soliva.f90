!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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
!> \file soliva.f90
!> \brief Atmospheric soil module - soil variables initialisation

!> \brief  Initialize soil model variables
!>  NB : data structures are define in module atsoil.f90
!-------------------------------------------------------------------------------

subroutine soliva

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use atincl
use atsoil
use ctincl, only: cp_a, cp_v
use mesh

implicit none

!===============================================================================

! Local variables

integer          ifac,iphysi
double precision rscp
double precision esaini,qsaini,huini,psini
double precision cpvcpa

!===============================================================================

!     1 - initialisation du tableau solva
!     -----------------------------------

cpvcpa = cp_v / cp_a

!  initialisation de t et qv en z0

if(qvsini.gt.1.d0) then
!  si qvsini>1 qvsini represente l'humidite relative en %
     esaini = 610.78d0*exp(17.2694d0*tsini/                           &
                      (tsini + tkelvi - 35.86d0))
     qsaini = esaini/(rvsra*p0 + esaini*(1.d0 - rvsra))
     qvsini = qvsini*qsaini/100.d0
endif

!==============================================================

if ( ippmod(iatmos).eq.1 ) iphysi = 0
if ( ippmod(iatmos).eq.2 ) iphysi = 3

do ifac = 1, nfmodsol

  psini = p0
  if(iphysi.gt.0)then
    rscp = (rair/cp0)*(1.d0 + (rvsra - cpvcpa)*qvsini)
  else
    rscp = (rair/cp0)
  endif

  solution_sol(ifac)%temp_sol = tsini

  solution_sol(ifac)%tempp = (tsini+tkelvi)*((ps/psini)**rscp)

  solution_sol(ifac)%total_water = 0.d0
  if(iphysi.gt.0) solution_sol(ifac)%total_water = qvsini
  solution_sol(ifac)%w1 = 0.d0
  solution_sol(ifac)%w2 = 0.d0

  if(iphysi.eq.3) then

   if(w1ini.lt.1.d-20) then
     !  si w1ini = 0 on fait un calcul approximatif de w1ini
     esaini = 610.78d0*exp(17.2694d0*tsini/                           &
                      (tsini + tkelvi - 35.86d0))
     qsaini = esaini/(rvsra*psini + esaini*(1.d0 - rvsra))
     huini  = qvsini/qsaini
     huini  = min(huini,1.d0)
     solution_sol(ifac)%w1 = acos(1.d0-2.d0*huini)/acos(-1.d0)
   else
     solution_sol(ifac)%w1  = w1ini
   endif

   if(w2ini.lt.1.d-20) then
!  si w2ini=0 on fixe le rapport w1ini/w2ini a 1
!  (equilibre entre les couches)
     solution_sol(ifac)%w2 = solution_sol(ifac)%w1
   else
     solution_sol(ifac)%w2 = w2ini
   endif

  endif

enddo

return
end subroutine
