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


!> \file alstru.f90
!> Module for ALE structure movement with internal coupling

module alstru

  !=============================================================================

  use paramx

  implicit none

  !=============================================================================

  !  Methode ale - mouvement de structures en couplage interne

  ! nbstru : nombre de structures mobiles

  ! xmstru : matrice de masse des structures (kgl)
  ! xcstru : matrice de friction des structures (kg/s)
  ! xkstru : matrice de raideur des structures (kg/s2 = N/m)
  ! xstreq : vecteur ecart de la position des structure dans le maillage
  !          initial par rapport a leur position d'equilibre (m)
  ! xstr   : vecteur deplacement des structures par rapport a leur position
  !          dans le maillage initial (m)
  ! xpstr  : vecteur vitesse des structures (m/s)
  ! xppstr : vecteur acceleration des structures (m/s2)
  ! forstr : vecteur force exerce sur les structures (N)
  ! xsta   : valeur de xstr au pas de temps precedent
  ! xpsta  : valeur de xpstr au pas de temps precedent
  ! xppsta : valeur de xppstr au pas de temps precedent
  ! forsta : valeur de forstr au pas de temps precedent
  ! xstp   : valeur predite de xstr
  ! forstp : valeur predite de forstr
  ! dtstr  : pas de temps associe au mouvement des structures
  ! aexxst : coefficient de prediction du deplacement (sur xpstr)
  ! bexxst : coefficient de prediction du deplacement (sur xpstr-xpsta)
  ! cfopre : coefficient de prediction des efforts

  integer, save :: nbstru

  double precision, save :: xmstru(3,3,nstrmx)
  double precision, save :: xcstru(3,3,nstrmx)
  double precision, save :: xkstru(3,3,nstrmx)
  double precision, save :: xstr(3,nstrmx)  ,xsta(3,nstrmx)
  double precision, save :: xstp(3,nstrmx)  ,xstreq(3,nstrmx)
  double precision, save :: xpstr(3,nstrmx) ,xpsta(3,nstrmx)
  double precision, save :: xppstr(3,nstrmx),xppsta(3,nstrmx)
  double precision, save :: forstr(3,nstrmx),forsta(3,nstrmx)
  double precision, save :: forstp(3,nstrmx)
  double precision, save :: dtstr(nstrmx)
  double precision, save :: aexxst, bexxst, cfopre

  ! methode de Newmark hht

  !  alpnmk : coefficient alpha de la methode de Newmark hht
  !  betnmk : coefficient beta  de la methode de Newmark hht
  !  gamnmk : coefficient gamma de la methode de Newmark hht

  double precision, save :: alpnmk, gamnmk, betnmk

  !=============================================================================

end module alstru

