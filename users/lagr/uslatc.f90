!-------------------------------------------------------------------------------

!VERS


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

subroutine uslatc &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   ,          &
   nprfml , nnod   , lndfac , lndfbr , ncelbr ,                   &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   numpt  , itepa  , idevel , ituser , ia     ,                   &
   rep    , uvwr   , romf   , romp   , xnul   ,                   &
   xcp    , xrkl   , tauc   ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , propce , propfa , propfb ,                   &
   ettp   , ettpa  , tepa   ,                                     &
   rdevel , rtuser , ra        )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!    SOUS-PROGRAMME UTILISATEUR (INTERVENTION NON OBLIGATOIRE)

!    MODIFICATION DU CALCUL DU TEMPS DE RELAXATION THERMIQUES
!      DES PARTICULES EN FONCTION DE LA FORMULATION CHOISIE
!      LE CALCUL DU NUSSELT


!               m   Cp
!                p    p
!      Tau = ---------------
!         c          2
!               PI d    h
!                   p    e

!     Tau  : TEMPS DE RELAXATION THERMIQUE (VALEUR A CALCULER)
!        c

!     m    : MASSE DE LA PARTICULE
!      p

!     Cp   : CHALEUR SPECIFIQUE DE LA PARTICULE
!       p

!     d    : DIAMETRE DE LA PARTICULE
!      p

!     h    : COEFFICIENT D'ECHANGE THERMIQUE
!      e

!    LE COEFFICIENT D'ECHANGE THERMIQUE EST CALCULE A PARTIR
!      D'UN NOMBRE DE NUSSELT, LUI MEME EVALUE PAR UNE
!      CORRELATION (RANZ-MARSHALL EN STANDARD).

!            h  d
!             e  p
!     Nu = --------  = 2 + 0.55 Re **(0.5) Prt**(0.33)
!           Lambda                p

!     Lambda : CONDUCTIVITE THERMIQUE DE LA PHASE PORTEUSE

!     Re     : NOMBRE DE REYNOLDS PARTICULAIRES
!       p

!     Prt    : NOMBRE DE PRANDTL


!    CE SOUS PROGRAMME EST APPELE DANS UNE BOUCLE SUR
!      LES PARTICULES : ATTENTION DE NE PAS TROP LE CHARGER

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ndim             ! e  ! <-- ! dimension de l'espace                          !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nfml             ! e  ! <-- ! nombre de familles d entites                   !
! nprfml           ! e  ! <-- ! nombre de proprietese des familles             !
! nnod             ! e  ! <-- ! nombre de sommets                              !
! lndfac           ! e  ! <-- ! longueur du tableau nodfac                     !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr                     !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! numpt            ! e  ! <-- ! numero de la particule courante                !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! te ! --- ! macro tableau entier                           !
! rep              ! r  ! <-- ! nombre de reynolds particulaire                !
!                  !    !     ! rep = uvwr * ettp(numpt,jdp) / xnul            !
! uvwr             ! r  ! <-- ! vitesse relative de la particule               !
!                  !    !     ! uvwr= |vit fluide vu - vit particule|          !
! romf             ! r  ! <-- ! masse volumique du fluide a la                 !
!                  !    !     ! position de la particule                       !
! romp             ! r  ! <-- ! masse volumique de la particule                !
! xnul             ! r  ! <-- ! viscosite cinematique du fluide a la           !
!                  !    !     ! position de la particule                       !
! xcp              ! r  ! <-- ! chaleur specifique du fluide                   !
!                  !    !     ! a la position de la particule                  !
! xrkl             ! r  ! <-- ! coefficient de diffusion du fluide             !
!                  !    !     ! a la position de la particule                  !
! tauc             ! r  ! --> ! temps de relaxation thermique                  !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
! (ndim,ncelet     !    !     !                                                !
! surfac           ! tr ! <-- ! vecteur surface des faces internes             !
! (ndim,nfac)      !    !     !                                                !
! surfbo           ! tr ! <-- ! vecteur surface des faces de bord              !
! (ndim,nfabor)    !    !     !                                                !
! cdgfac           ! tr ! <-- ! centre de gravite des faces internes           !
! (ndim,nfac)      !    !     !                                                !
! cdgfbo           ! tr ! <-- ! centre de gravite des faces de bord            !
! (ndim,nfabor)    !    !     !                                                !
! xyznod           ! tr ! <-- ! coordonnes des noeuds                          !
! (ndim,nnod)      !    !     !                                                !
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet          !    !     !                                                !
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
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

include "paramx.h"
include "numvar.h"
include "cstnum.h"
include "cstphy.h"
include "optcal.h"
include "entsor.h"
include "lagpar.h"
include "lagran.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "cpincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          nideve , nrdeve , nituse , nrtuse
integer          numpt
integer          itepa(nbpmax,nivep)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision rep    , uvwr   , romf   , romp   , xnul
double precision xcp    , xrkl   , tauc
double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision dt(ncelet) , rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision rdevel(nrdeve), rtuser(nrtuse)
double precision ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ip

! VARIABLES LOCALES UTILISATEUR

double precision prt, fnus

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

ip = numpt

!===============================================================================
! 2. TEMPS DE RELAXATION THERMIQUE STANDARD
!===============================================================================

!   Cet exemple est desactive, il donne le temps de relaxation thermique
!   standard a titre indicatif.

if (1.eq.0) then

  prt  = xnul / xrkl

  fnus = 2.d0 + 0.55d0 * rep**0.5d0 * prt**(1.d0/3.d0)

  tauc = ettp(ip,jdp) *ettp(ip,jdp) * romp * ettp(ip,jcp)         &
           / ( fnus * 6.d0 * romf * xcp * xrkl )

endif

!==============================================================================

!--------
! FORMATS
!--------


!----
! FIN
!----

end
