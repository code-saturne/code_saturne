!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

subroutine laggeo
!================

!===============================================================================
! Purpose:
! ----------
!
!   Subroutine of the Lagrangian particle-tracking module:
!   ------------------------------------------------------
!
!
!   Deposition sub-model:
!
!   Construction of the geometric data needed by the model
!
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use entsor
use cstphy
use cstnum
use pointe
use period
use parall
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use cpincl
use radiat
use mesh

!===============================================================================

implicit none

! Arguments

! Local variables

integer          ifac , inoeud
double precision xs1,ys1,zs1,xs2,ys2,zs2,xs3,ys3,zs3
double precision xnor
double precision xn,yn,zn,xt,yt,zt,xtt,ytt,ztt

!===============================================================================
! 0.  Memory management and crossing counter
!===============================================================================


!===============================================================================
! 1.  Boundary faces planes computation
!     A X + B Y + C Z + D = 0
!===============================================================================

do ifac=1,nfabor

! Recover the 3 first face nodes

  inoeud=nodfbr(ipnfbr(ifac))
  xs1 = xyznod(1,inoeud)
  ys1 = xyznod(2,inoeud)
  zs1 = xyznod(3,inoeud)

  inoeud=nodfbr(ipnfbr(ifac)+1)
  xs2 = xyznod(1,inoeud)
  ys2 = xyznod(2,inoeud)
  zs2 = xyznod(3,inoeud)

  inoeud=nodfbr(ipnfbr(ifac)+2)
  xs3 = xyznod(1,inoeud)
  ys3 = xyznod(2,inoeud)
  zs3 = xyznod(3,inoeud)

! Face plane equation

  xnor = sqrt( surfbo(1,ifac)*surfbo(1,ifac)                      &
              +surfbo(2,ifac)*surfbo(2,ifac)                      &
              +surfbo(3,ifac)*surfbo(3,ifac) )
  dlgeo(ifac,1)= surfbo(1,ifac)/xnor
  dlgeo(ifac,2)= surfbo(2,ifac)/xnor
  dlgeo(ifac,3)= surfbo(3,ifac)/xnor
  dlgeo(ifac,4)=-( dlgeo(ifac,1)*xs1                              &
                  +dlgeo(ifac,2)*ys1                              &
                  +dlgeo(ifac,3)*zs1 )

! Matrix of Reference frame change

  xn = dlgeo(ifac,1)
  yn = dlgeo(ifac,2)
  zn = dlgeo(ifac,3)

  xnor = sqrt( (xs2-xs1)*(xs2-xs1)                                &
                +(ys2-ys1)*(ys2-ys1)                              &
                +(zs2-zs1)*(zs2-zs1) )
  xt = (xs2-xs1)/xnor
  yt = (ys2-ys1)/xnor
  zt = (zs2-zs1)/xnor

  xtt = yn*zt-zn*yt
  ytt = zn*xt-xn*zt
  ztt = xn*yt-yn*xt
  xnor=sqrt(xtt*xtt+ytt*ytt+ztt*ztt)
  xtt = xtt/xnor
  ytt = ytt/xnor
  ztt = ztt/xnor

  dlgeo(ifac, 5) = xn
  dlgeo(ifac, 6) = yn
  dlgeo(ifac, 7) = zn
  dlgeo(ifac, 8) = xt
  dlgeo(ifac, 9) = yt
  dlgeo(ifac,10) = zt
  dlgeo(ifac,11) = xtt
  dlgeo(ifac,12) = ytt
  dlgeo(ifac,13) = ztt

enddo

!===============================================================================

!--------
! FORMATS
!--------
!

!----
! FIN
!----

end subroutine
