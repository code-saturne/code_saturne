!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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
   nthdi2 , nthdb2 , ngrpi2 , ngrpb2 , idxfi  , idxfb  , &
   iface2 , ifabo2 , ifmfb2 , ifmce2 , iprfm2,           &
   ipnfa2 , nodfa2 , ipnfb2 , nodfb2 ,                   &
   xyzce2 , surfa2 , surfb2 , cdgfa2 , cdgfb2 , xyzno2 , &
   volum2                                                &
)

!===============================================================================
! FONCTION :
! ---------

! PASSAGE DES DIMENSIONS DU MAILLAGE DU C AU FORTRAN.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
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

!===============================================================================
! Module files
!===============================================================================

use dimens
use paramx
use entsor
use parall
use mesh

!===============================================================================

implicit none

! Arguments

integer, intent(in) :: ncel2, ncele2, nfac2, nfabo2, nsom2
integer, intent(in) :: lndfa2, lndfb2
integer, intent(in) :: ncelg2, nfacg2 , nfbrg2, nsomg2
integer, intent(in) :: nthdi2, nthdb2
integer, intent(in) :: ngrpi2, ngrpb2

integer, dimension(*), intent(in) :: idxfi, idxfb

integer, dimension(2,nfac2), target :: iface2
integer, dimension(ncele2), target :: ifmce2
integer, dimension(nfabo2), target :: ifabo2, ifmfb2
integer, dimension(nfml,nprfml), target :: iprfm2
integer, dimension(nfac2+1), target :: ipnfa2
integer, dimension(lndfa2), target :: nodfa2
integer, dimension(nfabo2+1), target :: ipnfb2
integer, dimension(lndfb2), target :: nodfb2

double precision, dimension(3,ncele2), target :: xyzce2
double precision, dimension(3,nfac2), target :: surfa2, cdgfa2
double precision, dimension(3,nfabo2), target :: surfb2, cdgfb2
double precision, dimension(3,nsom2), target :: xyzno2
double precision, dimension(ncele2), target :: volum2

! Local variables

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

call init_fortran_omp(nfac, nfabor, &
                      nthdi2, nthdb2, ngrpi2, ngrpb2, idxfi, idxfb)

!===============================================================================
! 6. DEFINITION DES POINTEURS SUR LA STRUCTURE MAILLAGE
!===============================================================================

ifacel => iface2(1:2,1:nfac)
ifabor => ifabo2(1:nfabor)

ifmfbr => ifmfb2(1:nfabor)
ifmcel => ifmce2(1:ncelet)
iprfml => iprfm2(1:nfml,1:nprfml)

ipnfac => ipnfa2(1:nfac+1)
nodfac => nodfa2(1:lndfac)
ipnfbr => ipnfb2(1:nfabor+1)
nodfbr => nodfb2(1:lndfbr)

xyzcen => xyzce2(1:3,1:ncelet)

surfac => surfa2(1:3,1:nfac)
surfbo => surfb2(1:3,1:nfabor)
cdgfac => cdgfa2(1:3,1:nfac)
cdgfbo => cdgfb2(1:3,1:nfabor)

xyznod => xyzno2(1:3,1:nnod)

volume => volum2(1:ncelet)

return
end subroutine
