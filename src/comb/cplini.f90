!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

subroutine cplini

!===============================================================================
! FONCTION :
! --------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
!   --------------------------------------------------------------

!    ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE

!      COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
!      TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON

!      INITIALISATION DES VARIABLES DE CALCUL

!      PENDANT DE USINIV.F

! Cette routine est appelee en debut de calcul (suite ou non)
!     avant le debut de la boucle en temps

! Elle permet d'INITIALISER ou de MODIFIER (pour les calculs suite)
!     les variables de calcul,
!     les valeurs du pas de temps


! On dispose ici de ROM et VISCL initialises par RO0 et VISCL0
!     ou relues d'un fichier suite
! On ne dispose des variables VISCLS, CP (quand elles sont
!     definies) que si elles ont pu etre relues dans un fichier
!     suite de calcul

! LA MODIFICATION DES PROPRIETES PHYSIQUES (ROM, VISCL, VISCLS, CP)
!     SE FERA EN STANDARD DANS LE SOUS PROGRAMME PPPHYV
!     ET PAS ICI

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
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
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use mesh
use field

!===============================================================================

implicit none

! Local variables

integer          iel, ige, mode, icha

double precision t1init, h1init, coefe(ngazem)
double precision f1mc(ncharm), f2mc(ncharm)
double precision xkent, xeent, d2s3

double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_phi
double precision, dimension(:), pointer :: cvar_fb, cvar_omg, cvar_nusa
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: cvar_r12, cvar_r13, cvar_r23
double precision, dimension(:), pointer :: cvar_scalt
double precision, dimension(:), pointer :: cvar_f1m, cvar_f2m
double precision, dimension(:), pointer :: cvar_f3m, cvar_f4p2m

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================

ipass = ipass + 1


d2s3 = 2.d0/3.d0

if (itytur.eq.2) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)
elseif (itytur.eq.3) then
  call field_get_val_s(ivarfl(ir11), cvar_r11)
  call field_get_val_s(ivarfl(ir22), cvar_r22)
  call field_get_val_s(ivarfl(ir33), cvar_r33)
  call field_get_val_s(ivarfl(ir12), cvar_r12)
  call field_get_val_s(ivarfl(ir13), cvar_r13)
  call field_get_val_s(ivarfl(ir23), cvar_r23)
  call field_get_val_s(ivarfl(iep), cvar_ep)
elseif (iturb.eq.50) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)
  call field_get_val_s(ivarfl(iphi), cvar_phi)
  call field_get_val_s(ivarfl(ifb), cvar_fb)
elseif (iturb.eq.60) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iomg), cvar_omg)
elseif (iturb.eq.70) then
  call field_get_val_s(ivarfl(inusa), cvar_nusa)
endif

call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)

!===============================================================================
! 2. INITIALISATION DES INCONNUES :
!      UNIQUEMENT SI ON NE FAIT PAS UNE SUITE
!===============================================================================

! RQ IMPORTANTE : pour la combustion CP, 1 seul passage suffit

if ( isuite.eq.0 .and. ipass.eq.1 ) then

! --> Initialisation de k et epsilon

  xkent = 1.d-10
  xeent = 1.d-10

  if (itytur.eq.2) then

    do iel = 1, ncel
      cvar_k(iel)  = xkent
      cvar_ep(iel) = xeent
    enddo

  elseif (itytur.eq.3) then

    do iel = 1, ncel
      cvar_r11(iel) = d2s3*xkent
      cvar_r22(iel) = d2s3*xkent
      cvar_r33(iel) = d2s3*xkent
      cvar_r12(iel) = 0.d0
      cvar_r13(iel) = 0.d0
      cvar_r23(iel) = 0.d0
      cvar_ep(iel)  = xeent
    enddo

  elseif (iturb.eq.50) then

    do iel = 1, ncel
      cvar_k(iel)   = xkent
      cvar_ep(iel)  = xeent
      cvar_phi(iel) = d2s3
      cvar_fb(iel)  = 0.d0
    enddo

  elseif (iturb.eq.60) then

    do iel = 1, ncel
      cvar_k(iel)   = xkent
      cvar_omg(iel) = xeent/cmu/xkent
    enddo

  elseif (iturb.eq.70) then

    do iel = 1, ncel
      cvar_nusa(iel) = cmu*xkent**2/xeent
    enddo

  endif

! --> On initialise tout le domaine de calcul avec de l'air a TINITK
!                   ================================================

!   Enthalpie

  t1init = t0

! ------ Variables de transport relatives au melange

  do ige = 1, ngazem
    coefe(ige) = zero
  enddo
  coefe(io2) = wmole(io2) / (wmole(io2)+xsi*wmole(in2))
  coefe(in2) = 1.d0 - coefe(io2)
  do icha = 1, ncharm
    f1mc(icha) = zero
    f2mc(icha) = zero
  enddo
  mode = -1
  call cpthp1                                                     &
  !==========
 ( mode   , h1init , coefe  , f1mc   , f2mc   ,                   &
   t1init )

  do iel = 1, ncel
    cvar_scalt(iel) = h1init
  enddo

! Fraction massique et variance

  do icha = 1, ncharb
    call field_get_val_s(ivarfl(isca(if1m(icha))), cvar_f1m)
    call field_get_val_s(ivarfl(isca(if2m(icha))), cvar_f2m)
    do iel = 1, ncel
      cvar_f1m(iel) = zero
      cvar_f2m(iel) = zero
    enddo
  enddo
  call field_get_val_s(ivarfl(isca(if3m)), cvar_f3m)
  call field_get_val_s(ivarfl(isca(if4p2m)), cvar_f4p2m)
  do iel = 1, ncel
    cvar_f3m(iel) = zero
    cvar_f4p2m(iel) = zero
  enddo

endif

!----
! FORMATS
!----

!----
! FIN
!----

return
end subroutine
