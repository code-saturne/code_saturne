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

!> \file paramx.f90
!> Module for definition of general parameters

module paramx

  !=============================================================================

  implicit none

  !=============================================================================

  ! nscamx : nombre max de scalaires
  ! nvarmx : nombre max de variables =
  !          nombre max de scalaires + 12 (u,v,w,P,Rij,e,alp)
  ! nprcmx : nombre max de proprietes physiques aux cellules (total) =
  !          nscamx (Lambda) + 7 (rho,Cp,viscl,visct,cou,fou,iprtot)
  !                          + 4 (estim)
  ! nprfmx : nombre max de proprietes physiques aux faces internes =
  !          nscamx (flumas) + 2*(flumas,alp)
  ! nprbmx : nombre max de proprietes physiques aux faces de bord =
  !          nscamx (flumab) + 3*(flumab,alp, romb)
  ! npromx : nombre max de proprietes physiques tout confondu
  !          majore par nprcmx+nprfmx+nprbmx
  ! ngrdmx : nombre max de grandeurs =
  !          nvarmx + npromx
  ! nsmamx : nombre max de cases pour les tableaux termes source de masse
  !          nvarmx + 1 pour smacel
  ! nvppmx : nombre de variables pour affichages
  !          ngrdmx + 20 (20 couvre dt, tpucou, et une marge de 16 ...)

  integer   nscamx, nvarmx, nprcmx, nprfmx, nprbmx, npromx
  integer   ngrdmx, nsmamx, nvppmx
  parameter(nscamx=200)
  parameter(nvarmx=nscamx+12)
  parameter(nprcmx=nscamx+11)
  parameter(nprfmx=nscamx+ 2)
  parameter(nprbmx=nscamx+ 3)
  parameter(npromx=nprcmx+ nprfmx+nprbmx)
  parameter(ngrdmx=nvarmx+ npromx)
  parameter(nsmamx=nvarmx+ 1)
  parameter(nvppmx=ngrdmx+20)

  ! ntypmx nombre de types de conditions aux limites possibles

  integer    ntypmx
  parameter(ntypmx=200)

  integer    iindef, ientre, isolib, isymet, iparoi,                &
             iparug, iesicf, isspcf, isopcf, ierucf,                &
             ieqhcf, icscpl

  parameter(iindef=1, ientre=2, isolib=3, isymet=4, iparoi=5,       &
            iparug=6, iesicf=7, isspcf=8, isopcf=9, ierucf=10,      &
            ieqhcf=11, icscpl=12)

  ! nestmx : nombre max d'estimateurs
  ! iespre, iesder, iescor, iestot : numeros
  integer    nestmx
  parameter (nestmx=4)
  integer    iespre  , iesder  , iescor  , iestot
  parameter (iespre=1, iesder=2, iescor=3, iestot=4)

  ! nbmomx : nombre max de moyennes (moments) calcule
  ! ndgmox : degre max des moments
  integer    nbmomx, ndgmox
  parameter (nbmomx = 50, ndgmox = 5)

  ! conditions aux limites possibles pour la vitesse de maillage en ale

  integer    ibfixe, igliss, ivimpo, ifresf
  parameter(ibfixe=1, igliss=2, ivimpo=3, ifresf=4)

  ! nstrmx : nombre de structures max en ale

  integer nstrmx
  parameter (nstrmx=200)

  !=============================================================================

end module paramx
