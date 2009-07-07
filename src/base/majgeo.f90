!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

subroutine majgeo &
!================

 ( ncel2  , ncele2 , nfac2  , nfabo2 , nsom2 ,           &
   lndfa2 , lndfb2 , ncelg2 , nfacg2 , nfbrg2 , nsomg2 , &
   nthdi2 , nthdb2 , ngrpi2 , ngrpb2 , idxfi  , idxfb )

!===============================================================================
! FONCTION :
! ---------

! PASSAGE DES DIMENSIONS DU MAILLAGE DU C AU FORTRAN.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ncel2            ! e  ! <-- ! nombre de cellules                             !
! ncele2           ! e  ! <-- ! nombre d'elements halo compris                 !
! nfac2            ! e  ! <-- ! nombre de faces internes                       !
! nfabo2           ! e  ! <-- ! nombre de faces de bord                        !
! nsom2            ! e  ! <-- ! nombre de sommets                              !
! lndfa2           ! e  ! <-- ! taille de lndfac                               !
! lndfb2           ! e  ! <-- ! taille de lndfbr                               !
! ncelg2           ! e  ! <-- ! nombre global de cellules                      !
! nfacg2           ! e  ! <-- ! nombre global de faces internes                !
! nfbrg2           ! e  ! <-- ! nombre global de faces de bord                 !
! nsomg2           ! e  ! <-- ! nombre global de sommets                       !
! nthdi2           ! e  ! <-- ! nb. max de threads par groupe de faces inter   !
! nthdb2           ! e  ! <-- ! nb. max de threads par groupe de faces de bord !
! ngrpi2           ! e  ! <-- ! nb. groupes de faces interieures               !
! ngrpb2           ! e  ! <-- ! nb. groupes de faces de bord                   !
! idxfi            ! e  ! <-- ! index pour faces internes                      !
! idxfb            ! e  ! <-- ! index pour faces de bord                       !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "dimens.h"
include "dimfbr.h"
include "paramx.h"
include "entsor.h"
include "parall.h"

! Arguments

integer          ncel2, ncele2, nfac2, nfabo2, nsom2
integer          lndfa2, lndfb2
integer          ncelg2, nfacg2 , nfbrg2, nsomg2
integer          nthdi2, nthdb2
integer          ngrpi2, ngrpb2
integer          idxfi(*), idxfb(*)

! VARIABLES LOCALES

integer          ii, jj

!===============================================================================

!===============================================================================
! 1. MISE A JOUR DU NOMBRE DE CELLULES
!===============================================================================

ncel = ncel2
ncelet = ncele2

!===============================================================================
! 2. MISE A JOUR DU NOMBRE DES FACES
!===============================================================================

nfac = nfac2
nfabor = nfabo2

lndfac = lndfa2
lndfbr = lndfb2

!     On remplit maintenant NDIMFB
if (nfabor.eq.0) then
  ndimfb = 1
else
  ndimfb = nfabor
endif

!===============================================================================
! 3. MISE A JOUR DU NOMBRE DES SOMMETS
!===============================================================================

nnod = nsom2

!===============================================================================
! 4. MISE A JOUR DES TAILLES GLOBALES
!===============================================================================

ncelgb = ncelg2
nfacgb = nfacg2
nfbrgb = nfbrg2
nsomgb = nsomg2

!===============================================================================
! 5. INITIALISATION DES INFORMATIONS SUR LES THREADS
!===============================================================================

do ii = 1, nthrd1
  do jj = 1, nthrd2
    iompli(1, ii, jj) = 0
    iompli(2, ii, jj) = 0
    iomplb(1, ii, jj) = 0
    iomplb(2, ii, jj) = 0
  enddo
enddo

! Pour le groupe j et le thread i, boucles sur les faces
! de iompl.(1, j, i) à iompl.(2, j, i).

! Par defaut (i.e. sans Open MP), 1 thread et un groupe

iompli(1, 1, 1) = 1
iompli(2, 1, 1) = nfac
iomplb(1, 1, 1) = 1
iomplb(2, 1, 1) = nfabor

! Numerotations pour boucles OpenMP sur les faces interieures

ngrpi = ngrpi2
nthrdi = nthdi2

if (nthrdi.gt.1 .and. ngrpi.gt.1) then

  do ii = 1, nthrdi
    do jj = 1, ngrpi

      iompli(1, jj, ii) = idxfi((ii-1)*ngrpi*2 + 2*jj-1) + 1
      iompli(2, jj, ii) = idxfi((ii-1)*ngrpi*2 + 2*jj)

    enddo
  enddo

endif

! Numerotations pour boucles OpenMP sur les faces de bord

ngrpb = ngrpb2
nthrdb = nthdb2

if (nthrdb.gt.1 .and. ngrpb.gt.1) then

  do ii = 1, nthrdb
    do jj = 1, ngrpb

      iomplb(1, jj, ii) = idxfb((ii-1)*ngrpb*2 + 2*jj-1) + 1
      iomplb(2, jj, ii) = idxfb((ii-1)*ngrpb*2 + 2*jj)

    enddo
  enddo

endif
return
end
