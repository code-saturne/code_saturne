!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine etheq &
!================

 (pphy,thetal,qw,qldia,xnebdia,xnn,etheta,eq)

!===============================================================================
! FONCTION :
! ----------
! Calcul des coefficients etheta et eq dependant du degre de saturation

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!  pphy            !  tr! d  ! pression reelle en pascals                      !
! thetal           !  tr! d  ! temperature potentielle liquide transportee     !
!   qw             !  tr! d  ! quantite d'eau totale transportee               !
!  qldia           !  tr! d  ! quantite moyenne d'eau liquide                  !
! xnebdia          !  tr! d  ! nebulosite issue du diagnostique                !
!   xnn            !    !    ! moment d'ordre 2 "n" <s'ql'>/(2*sigma_s**2)     !
! etheta           !  tr! r  ! coefficient etheta                              !
!   eq             !  tr! r  ! coefficient eq                                  !
!------------------------------------------------------------------------------!
!    e             !  tr! d  ! energie cinetique turbulente                    !
!   eps            !  tr! d  ! dissipation                                     !
!------------------------------------------------------------------------------!
! model            !  e ! d  ! indice du modelisation des coef. etheta, eq     !
!__________________!____!____!_________________________________________________!
!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstphy
use ppthch
use entsor
use ppppar ! defines nozppm which is used in atincl
use atincl

!===============================================================================

implicit none

! Arguments

double precision pphy,xnn
double precision qldia,xnebdia
double precision thetal,qw
double precision etheta,eq

! Variables locales

double precision cp,lscp,beta1,a,alpha1
double precision tl,qsl,de,alpha2,beta2,qs
double precision rn,dd,aa,theta,t,q

! Declaration des fonctions

external         qsatliq
double precision qsatliq

!===============================================================================
! 0. The most probable and simplest case
!===============================================================================

if (qldia.le.0.d0.or.modsub.eq.0)then
  etheta = 1.d0
  eq = (rvsra-1.d0)*thetal
  return
endif

! ===============================================================================
! 1. Initialisation
! ===============================================================================

cp   = cp0
lscp = clatev/cp

etheta = 1.d0
eq = (rvsra - 1.d0)*thetal

! ===============================================================================
! 2. Modele de sous-maille
! ===============================================================================

! k-eps standard

! Calcul de la temperature "liquide"
tl = (thetal*(p0/pphy)**(-rair/cp))

! Constantes
beta1 = (clatev*clatev)/(cp*(rair*rvsra)*tl*tl)

! Calcul de ql de saturation
qsl = qsatliq(tl,pphy)

a = 1.d0/(1.d0 + beta1*qsl)

alpha1 = qsl*beta1*((pphy/p0)**(rair/cp))/lscp

! Temperature potentielle (seche)
theta = thetal+(clatev/cp)*((p0/pphy)**(rair/cp))*qldia

if (modsub.eq.0) then
  etheta = 1.d0
  eq     = (rvsra - 1.d0)*theta
  return
endif

! Calcul de la temperature "thermodynamique"
t = theta*(p0/pphy)**(-rair/cp)

! Calcul de q
q = qw - qldia

! Calcul de q saturation
qs = qsatliq(t,pphy)

! Calcul de de
de = (clatev/cp)*((p0/pphy)**(rair/cp)) - rvsra*theta

! Constantes pour le moddis = 3
beta2  = (clatev*clatev)/(cp*(rair*rvsra)*t*t)
aa     = 1.d0/(1.d0 + beta2*qs)
alpha2 = qs*(beta2*cp/clatev)*((pphy/p0)**(rair/cp))

! Nouveau calcul de d pour le moddis=3
dd = (clatev/cp)*((p0/pphy)**(rair/cp))                                         &
   * (1.d0 + (rvsra - 1.d0 )*q - qldia) - rvsra*theta

! Nebulosite
rn = xnebdia

if (modsub.eq.1) then
  ! Bechtold et al. 1995
  etheta = 1.d0 - a*alpha1*de*xnn
  eq     = (rvsra - 1.d0)*theta + a*de*xnn
elseif (modsub.eq.2) then
  ! Bouzereau et al. 2004
  etheta = 1.d0 + (rvsra-1.d0)*q - qldia - a*alpha1*dd*xnn
  eq     = (rvsra - 1.d0)*theta + a*dd*xnn
elseif (modsub.eq.3) then
  ! Cuijpers et Duynkerke 1993, etc.
  ! Interpolation lineaire entre le cas sature et non-sature
  ! (coefficient de nebulosite partielle r)
  etheta = 1.d0 + (rvsra - 1.d0)*q - rn*(qldia + aa*alpha2*dd)
  eq     = (rvsra - 1.d0)*theta + rn*aa*dd
endif

! ----
!  fin
! ----

end subroutine etheq

