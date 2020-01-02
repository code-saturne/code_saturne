!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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
!> \file raysze.f90
!> \brief 1D Radiative scheme - Solar data + zenithal angle)

!> \brief   Compute :
!>-     - zenithal angle
!>-     - solar contant (with correction for earth - solar length)
!>-     - albedo if above the sea
!>   ( utilisation des formules analytiques de paltrige et platt
!>                 dev.in atm. science no 5)
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]   xlat        latitude
!> \param[in]   xlong       longitude
!> \param[in]   jour        day in the year
!> \param[in]   heurtu      Universal time (hour)
!> \param[in]   imer        sea index
!> \param[in]   albe        albedo
!> \param[out]  muzero      cosin of zenithal angle
!> \param[out]  fo          solar constant
!-------------------------------------------------------------------------------
subroutine raysze  (xlat,xlong,jour,heurtu,imer,albe,muzero,fo)

implicit none

! Arguments

double precision xlat,xlong,jour,heurtu,muzero,fo,albe
integer imer

! Local variables

double precision pi,t00,decl,eqt,hr,ho,corfo,flat,flong,heure

!===============================================================================

!     1 - initialisations locales
!     ---------------------------
fo = 1370.d0
pi = 4.d0*atan(1.d0)

!       conversions sexagesimal-decimal

flat = xlat*pi/180.d0
flong = xlong*4.d0/60.d0

t00 = 2.d0*pi*jour/365.d0

!     2 - calcul de la declinaison (erreur maxi <3 mn)
!     ------------------------------------------------

decl = 0.006918d0 - 0.399912d0*dcos(t00) + 0.070257d0*dsin(t00)                   &
     - 0.006758d0*dcos(2.d0*t00) + 0.000907d0*dsin(2.d0*t00) - 0.002697d0*dcos(3.d0*t00) &
     + 0.001480d0*dsin(3.d0*t00)

!     3 - calcul de l'heure solaire locale
!     ------------------------------------
!   equation du temps     erreur maxi    35 secondes

eqt = (0.000075d0 + 0.001868d0*dcos(t00) - 0.032077d0*dsin(t00)                   &
    - 0.014615d0*dcos(2.d0*t00) - 0.040849d0*dsin(2.d0*t00))*12.d0/pi

heure = heurtu + flong + eqt

!   transfo    heure-radians

! On retire pi et on prend le modulo 2pi du resultat
if(heure.ge.12.d0) then
  hr = (heure - 12.d0)*pi/12.d0
else
  hr = (heure + 12.d0)*pi/12.d0
endif

!     4 - calcul du cosinus de l'angle zenithal
!     -----------------------------------------

muzero = dsin(decl)*dsin(flat) + dcos(decl)*dcos(flat)*dcos(hr)

!     5 - calcul de l'albedo sur mer qui depend de l'angle zenithal
!     -----------------------------------------

if(imer.eq.1) then
  ho = acos(muzero)
  ho = 180.d0*(pi/2.d0 - ho)/pi
  if(ho.lt.8.5d0) ho = 8.5d0
  if(ho.gt.60.d0) ho = 60.d0
  albe = 3.d0/ho
endif

!     6 - calcul de la constante solaire
!     ----------------------------------
!  correction distance terre-soleil
!                corfo=(r0/r)**2
!                          precision  mieux que e-04

corfo = 1.000110d0 + 0.034221d0*cos(t00) + 0.001280d0*sin(t00)                  &
      + 0.000719d0*cos(2.d0*t00) + 0.000077d0*sin(2.d0*t00)
fo = fo*corfo

return
end subroutine raysze
