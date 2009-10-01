!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine inipst &
!=================

 ( ipstvl , ipstbo , ipstsy , ipstze,                             &
   ipstmd , ntpst  , fmtpst , optpst )

!===============================================================================
! FONCTION :
! --------

! RENVOIE LES PARAMETRES CONTENUS DANS LES COMMONS ET UTILES
! A L'INITIALISATION DU POST-TRAITEMENT

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ipstvl           ! e  ! --> ! indicateur pour le maillage volumique          !
! ipstbo           ! e  ! --> ! indicateur pour le maillage de peau            !
! ipstsy           ! e  ! --> ! indicateur pour les maillages de               !
!                  !    !     !  peau couples avec syrthes                     !
! ipstze           ! e  ! --> ! indicateur pour les maillages de               !
!                  !    !     !  zone d'echange aero                           !
! ipstmd           ! e  ! --> ! indicateur de maillage deformable :            !
!                  !    !     !  0 : pas de deformation ;                      !
!                  !    !     !  1 : deformation a topologie cste              !
!                  !    !     !  2 : topologie variable                        !
!                  !    !     ! 10 : idem que 0 avec ecriture d'un             !
!                  !    !     !                      champ de deplact          !
!                  !    !     ! 11 : idem que 1 avec ecriture d'un             !
!                  !    !     !                      champ de deplact          !
!                  !    !     ! 12 : idem que 2 avec ecriture d'un             !
!                  !    !     !                      champ de deplact          !
! ntpst            ! e  ! --> ! frequence des sorties post-traitement          !
! fmtpst           ! a  ! --> ! nom du format de post-traitement               !
! fmtpst           ! a  ! --> ! options du format de post-traitement           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "entsor.h"

!===============================================================================

! Arguments

integer          ipstvl , ipstbo , ipstsy , ipstze
integer          ipstmd , ntpst

character        fmtpst(32)
character        optpst(96)

! VARIABLES LOCALES

integer          ii

!===============================================================================

! Transmission des paramètres du common Fortran à la partie C.

ipstvl = ichrvl
ipstbo = ichrbo
ipstsy = ichrsy
ipstze = ichrze

ipstmd = ichrmd
ntpst = ntchr

do ii = 1, len(fmtchr)
  fmtpst(ii) = fmtchr(ii:ii)
enddo
do ii = 1, len(optchr)
  optpst(ii) = optchr(ii:ii)
enddo

!===============================================================================

return

end subroutine
