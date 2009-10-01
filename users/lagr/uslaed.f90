!-------------------------------------------------------------------------------

!VERS


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

subroutine uslaed &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   itepa  , ibord  ,                                              &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , volume ,          &
   dt     , rtp    , propce , propfa , propfb ,                   &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct , tsvar  , &
   auxl1  , auxl2  , auxl3  , w1     , w2     , w3     ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     SOUS-PROGRAMME UTILISATEUR (INTERVENTION NON OBLIGATOIRE)

!     INTEGRATION DES EDS POUR LES VARIABLES UTILISATEUR
!       PAR DEFAUT LES VARIABLES SONT CONSTANTES


!                                         d T       T - PIP
!     LES EDS DOIVENT ETRE DE LA FORME : ----- = - ---------
!                                         d t         Tca


!     T : IIIIeme VARIABLE UTILISATEUR QUI EST REPEREE POUR LA
!         PARTICULE IP PAR :
!            T = ETTP(IP,JVLS(IIII))
!            T = ETTPA(IP,JVLS(IIII))

!     Tca : TEMPS CARACTERISTIQUE DE L'EDS
!              A DONNER DANS LE TABLEAU    AUXL1

!     PIP : COEFFICIENT DE L'EDS ("PSEUDO SECOND MEMBRE")
!              A DONNER DANS LE TABLEAU    AUXL2

!              SI LE SCHEMA CHOISI EST D'ORDRE 1 (NORDRE=1)
!              ALORS AU 1ER ET SEUL PASSAGE PIP S'EXPRIME
!              EN FONCTION DE DES GRANDEURS DU PAS DE TEMPS
!              PRECEDENT ETTPA

!              SI LE SCHEMA CHOISI EST D'ORDRE 2 (NORDRE=2)
!              ALORS AU 1ER PASSAGE (NOR=1) PIP S'EXPRIME
!              EN FONCTION DE DES GRANDEURS DU PAS DE TEMPS
!              PRECEDENT ETTPA, ET AU 2E PASSAGE (NOR=2)
!              PIP S'EXPRIME EN FONCTION DE DES GRANDEURS
!              DU PAS DE TEMPS COURANT ETTP


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
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! ibord            ! te ! <-- ! contient le numero de la                       !
!   (nbpmax)       !    !     !   face d'interaction part/frontiere            !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
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
! volume(ncelet    ! tr ! <-- ! volume d'un des ncelet elements                !
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
! taup(nbpmax)     ! tr ! <-- ! temps caracteristique dynamique                !
! tlag(nbpmax)     ! tr ! <-- ! temps caracteristique fluide                   !
! tempct           ! tr ! <-- ! temps caracteristique thermique et             !
!  (nbpmax,2)      !    !     !   ts implicite de couplage retour              !
! tsvar            ! tr ! <-- ! prediction 1er sous-pas pour la                !
! (nbpmax,nvp1)    !    !     !   variable ivar, utilise pour la               !
!                  !    !     !   correction au 2eme sous-pas                  !
! auxl1(nbpmax)    ! tr ! --- ! tableau de travail                             !
! auxl2(nbpmax)    ! tr ! --- ! tableau de travail                             !
! auxl3(nbpmax)    ! tr ! --- ! tableau de travail                             !
! w1..w3(ncelet    ! tr ! --- ! tableaux de travail                            !
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
include "cstphy.h"
include "cstnum.h"
include "optcal.h"
include "entsor.h"
include "lagpar.h"
include "lagran.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nideve , nrdeve , nituse , nrtuse
integer          itepa(nbpmax,nivep)  , ibord(nbpmax)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision volume(ncelet)
double precision dt(ncelet) , rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision taup(nbpmax) , tlag(nbpmax,3) , tempct(nbpmax,2)
double precision tsvar(nbpmax,nvp1)
double precision auxl1(nbpmax), auxl2(nbpmax), auxl3(nbpmax)
double precision w1(ncelet) ,  w2(ncelet) ,  w3(ncelet)
double precision rdevel(nrdeve) , rtuser(nrtuse)
double precision ra(*)

! VARIABLES LOCALES

integer          idebia , idebra
integer          npt , iel , iiii , ipl , iphas

!===============================================================================


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!===============================================================================

! On ne rentre en effet dans ce sous programme que si on a defini des
!   variables supplementaires dans uslag1 ; il faut donc necessairement
!   dire comment on les resout (non?).


if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS LE MODULE LAGRANGIEN             ',/,&
'@    =========                                               ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR uslaed DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Des variables supplementaires ont ete declarees dans      ',/,&
'@    uslag1 (NVLS=)                                          ',/,&
'@  Le sous-programme uslaed doit etre complete pour preciser ',/,&
'@    l''equation differentielle stochastique a resoudre.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


! FIXME : TODO ecrire un exemple utilisateur


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

idebia = idbia0
idebra = idbra0

iphas = ilphas

!===============================================================================
! 2. TEMPS CARACTERISTIQUE DE L'EDS COURANTE
!===============================================================================

! Boucle sur les variables supplementaires

do iiii = 1,nvls

!      Numero dans ETTP de la variable traitee

  ipl = jvls(iiii)

  do npt = 1,nbpart

    if ( itepa(npt,jisor).gt.0 ) then

      iel = itepa(npt,jisor)

!     Temps caracteristique Tca de l'equation differentielle
!     Cet exemple doit etre adapte au cas traite

      auxl1(npt) = 1.d0

!     Prediction au premier sous pas
!     Cet exemple doit etre adapte au cas traite

      if (nor.eq.1) then
        auxl2(npt) = ettpa(npt,ipl)
      else

!     Correction au deuxieme sous pas
!     Cet exemple doit etre adapte au cas traite

        auxl2(npt) = ettp(npt,ipl)
      endif

    endif
  enddo

!===============================================================================
! 3. INTEGRATION DE LA VARIABLE IPL
!===============================================================================

  call lagitg                                                     &
  !==========
   ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                 &
     ipl    ,                                                     &
     itepa(1,jisor)  , ibord  ,                                   &
     ettp   , ettpa  , auxl1  , auxl2  , tsvar  )

enddo

!===============================================================================

!----
! FIN
!----

end subroutine
