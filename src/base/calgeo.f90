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

subroutine calgeo &
!================

 ( idbia0 , idbra0 ,                                              &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idevel , ituser , ia     ,                                     &
   volmin , volmax , voltot ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
!  FONCTION  :
!  ---------

!  CALCUL DES ENTITES GEOMETRIQUES DEDUITES
!     DU JEU DE DONNEES MINIMAL

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! volmin           ! r  ! --> ! volume de controle minimal                     !
! volmax           ! r  ! --> ! volume de controle maximal                     !
! voltot           ! r  ! --> ! volume total du domaine                        !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use optcal
use pointe
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nideve , nrdeve , nituse , nrtuse

integer          idevel(nideve), ituser(nituse), ia(*)

double precision volmin, volmax, voltot
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer idebia, idebra

!===============================================================================
! 1. ON SAUVEGARDE LA POSITION DE LA MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 2. ON CALCULE LE VOLUME MIN et TOTAL DES ELEMENTS
!===============================================================================

call clvolc                                                       &
!==========
     ( ncelet , ncel   ,                                          &
       volmin , volmax , voltot , volume )

!===============================================================================
! 3. ON CALCULE LES SURFACES DES FACES
!===============================================================================

call clsurn                                                       &
!==========
     ( idebia , idebra ,                                          &
       nfac   , nfabor ,                                          &
       surfac , surfbo ,                                          &
       ra(isrfan) , ra(isrfbn) ,                                  &
       ia     , ra     )


!===============================================================================
! 4. ON CALCULE LE PRODUIT SCALAIRE DE LA NORMALE NORMEE A UNE FACE ET
!        DU VECTEUR DEFINI PAR LES VOISINS (VOISIN 1 : ORIGINE,
!        VOISIN 2 : EXTREMITE)
!               LA PONDERATION RESULTANTE   POND  = D2/(D1+D2)
!        OU D1 ET D2 SONT LES PROJETES SUR LA NORMALE A LA FACE DES
!        VECTEURS DEFINIS RESPECTIVEMENT PAR
!    D1 : (ORIGINE : VOISIN 1, EXTREMITE : CENTRE DE GRAVITE DE LA FACE)
!    D2 : (ORIGINE : CENTRE DE GRAVITE DE LA FACE, EXTREMITE : VOISIN 2)
!===============================================================================

call cldipo                                                       &
!==========
 ( idebia , idebra ,                                              &
   nfac   , nfabor , ncelet , ncel   ,                            &
   ifacel , ifabor ,                                              &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   ra(isrfan) , ra(isrfbn) ,                                      &
   ra(idist)  , ra(idistb) , ra(ipond) ,                          &
   ia     , ra     )

!===============================================================================
! 5. ON CALCULE LES VECTEURS IIP ET JJP POUR LES RECONSTRUCTIONS
!===============================================================================

call cldijp                                                       &
!==========
 ( idebia , idebra ,                                              &
   nfac   , nfabor , ncelet , ncel   ,                            &
   ifacel , ifabor ,                                              &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   ra(isrfan) , ra(isrfbn) ,                                      &
   ra(ipond)  ,                                                   &
   ra(idijpf) , ra(idiipb)  , ra(idofij) ,                        &
   ia     , ra     )


!===============================================================================
! 6. FILTRAGE DU VOISINAGE ETENDU POUR LE GRADIENT PAR MOINDRES CARRES
!===============================================================================

if (imrgra.eq.3) then

  call redvse (anomax)
  !==========

endif


!===============================================================================
! 8. FIN
!===============================================================================


return
end subroutine
