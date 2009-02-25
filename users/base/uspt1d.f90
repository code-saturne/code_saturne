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

                  subroutine uspt1d                               &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , nfpt1d , iphas  , iappel ,          &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ifpt1d , nppt1d , iclt1d ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   tppt1d , rgpt1d , eppt1d ,                                     &
   tept1d , hept1d , fept1d ,                                     &
   xlmt1d , rcpt1d , dtpt1d ,                                     &
   dt     , rtpa   ,                                              &
   propce , propfa , propfb ,                                     &
   coefa  , coefb  ,                                              &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

!     ENTREE DES DONNEES DU MODULE THERMIQUE EN 1D PAROI

! IAPPEL = 1 (un seul appel a l'initialisation) :
!             CALCUL DU NOMBRE DE CELLULES OU L'ON IMPOSE UNE PAROI

! IAPPEL = 2 (un seul appel a l'initialisation) :
!             REPERAGE DES CELLULES OU L'ON IMPOSE UNE PAROI
!             DONNEES RELATIVE AU MAILLAGE

! IAPPEL = 3 (appel a chaque pas de temps) :
!             VALEUR DES COEFFICIENTS PHYSIQUE DU CALCUL
!             TYPE DE CONDITION LIMITE EN PAROI EXTERIEURE
!             ICLT1D = 1 -> TEMPERATURE
!             ICLT1D = 3 -> FLUX
!             INITIALISATION DE LA TEMPERATURE


! IDENTIFICATION DES FACES DE BORD
! ================================
! L'identification des faces de bord concernees se fait grace
! a la commande GETFBR.

!  GETFBR(CHAINE,NLELT,LSTELT) :
!  - CHAINE est une chaine de caractere fournie par l'utilisateur
!    qui donne les criteres de selection
!  - NLTELT est renvoye par la commande. C'est un entier qui
!    correspond au nombre de faces de bord trouveees repondant au
!    critere
!  - LSTELT est renvoye par la commande. C'est un tableau d'entiers
!    de taille NLTELT donnant la liste des faces de bord trouvees
!    repondant au critere.

!  CHAINE peut etre constitue de :
!  - references de couleurs (ex. : 1, 8, 26, ...
!  - references de groupes (ex. : entrees, groupe1, ...)
!  - criteres geometriques (ex. X<0.1, Y>=0.25, ...)
!  Ces criteres peuvent etre combines par des operateurs logiques
!  (AND et OR) et des parentheses
!  ex. : '1 AND (groupe2 OR groupe3) AND Y<1' permettra de recuperer
!  les faces de bord de couleur 1, appartenant aux groupes 'groupe2'
!  ou 'groupe3' et de coordonnee Y inferieure a 1.

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nfml             ! e  ! <-- ! nombre de familles d entites                   !
! nprfml           ! e  ! <-- ! nombre de proprietese des familles             !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nfpt1d           ! e  ! <-- ! nombre de faces avec module therm 1d           !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! iappel           ! e  ! <-- ! indique les donnes a renvoyer                  !
! ifabor           ! te ! <-- ! element  voisin  d'une face de bord            !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! te ! <-- ! numero de famille d'une face de bord           !
! (nfabor)         !    !     !                                                !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
! nfml  ,nprfml    !    !     !                                                !
! maxelt           !  e ! <-- ! nb max d'elements (cell,fac,fbr)               !
! lstelt(maxelt) te ! --- ! tableau de travail                             !
! ifpt1d           ! te ! <-- ! numero de la face en traitement                !
!                  !    !     ! thermique en paroi                             !
! nppt1d           ! te ! <-- ! nombre de points de discretisation             !
!                  !    !     ! dans la paroi                                  !
! iclt1d           ! te ! <-- ! type de condition limite                       !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! eppt1d           ! tr ! <-- ! epaisseur de la paroi                          !
! rgpt1d           ! tr ! <-- ! raison du maillage                             !
! tppt1d           ! tr ! <-- ! temperature de paroi                           !
! tept1d           ! tr ! <-- ! temperature exterieure                         !
! hept1d           ! tr ! <-- ! coefficient d'echange exterieur                !
! fept1d           ! tr ! <-- ! flux exterieur                                 !
! xlmt1d           ! tr ! <-- ! conductivite thermique de la paroi             !
! rcpt1d           ! tr ! <-- ! rocp de la paroi                               !
! dtpt1d           ! tr ! <-- ! pas de temps de la paroi                       !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
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
include "entsor.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "parall.h"
include "period.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas  , nfpt1d
integer          nideve , nrdeve , nituse , nrtuse
integer          iphas  , iappel

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          ifpt1d(nfpt1d), nppt1d(nfpt1d), iclt1d(nfpt1d)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision eppt1d(nfpt1d) , rgpt1d(nfpt1d) , tppt1d(nfpt1d)
double precision tept1d(nfpt1d) , hept1d(nfpt1d) , fept1d(nfpt1d)
double precision xlmt1d(nfpt1d) , rcpt1d(nfpt1d) , dtpt1d(nfpt1d)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Variables locales

integer          idebia , idebra
integer          ifbt1d , ii , ifac
integer          ilelt, nlelt

!===============================================================================

idebia = idbia0
idebra = idbra0


! --- Relecture d'un fichier suite :
!     ISUIT1 = 0 ------> pas de relecture (reinitialisation du maillage et
!                                           de la temperature dans la paroi)
!     ISUIT1 = 1 ------> relecture du fichier suite de module thermique 1D
!     ISUIT1 = ISUITE -> relecture si le calcul fluide est une suite
!     L'initialisation de ISUIT1 dans uspt1d est obligatoire.

isuit1 = isuite

ifbt1d = 0

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

if(iappel.eq.1.or.iappel.eq.2) then


! --- Determination des faces avec module thermique 1D
!     NFPT1D     : nb total de faces avec module thermique 1D
!     IFPT1D(II) : numero de la IIeme face avec module thermque 1D

!     Remarque : lors de la relecture d'un fichier suite, NFPT1D et IFPT1D
!                sont compares aux valeurs issues du fichier suite. Une
!                concordance totale est necessaire pour continuer le calcul.
!                En ce qui concerne le test sur IFPT1D, il necessite qui
!                le tableau soit range dans un ordre croissant
!                   ( IFPT1D(JJ) > IFPT1D(II) si JJ > II ).
!                Si ce n'est pas possible, contacter l'equipe de developpement
!                pour desactiver le test.

  CALL GETFBR('3',NLELT,LSTELT)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)

    ifbt1d =ifbt1d + 1
    if (iappel.eq.2) ifpt1d(ifbt1d) = ifac

  enddo

endif

if (iappel.eq.1) then
   nfpt1d = ifbt1d
endif

! --- Remplissage des parametres du maillage et d'initialisation
!     (un seul passage en debut de calcul)
!     NPPT1D(II) : nombre de points de discretisation associes a la
!                                   IIeme face avec module thermique 1D
!     EPPT1D(II) : epaisseur de la paroi associee a la IIeme face avec
!                                                   module thermique 1D
!     RGPT1D(II) : raison geometrique du raffinement du maillage associe
!                              a la IIeme face avec module thermique 1D
!                       ( RGPT1D(II) > 1 => petites mailles cote fluide )
!     TPPT1D(II) : temperature d'initialisation de la paroi associee a la
!                                   IIeme face avec module thermique 1D

!     Remarque : lors de la relecture d'un fichier suite de module thermique
!                1D, TPPT1D n'est pas utilise. NFPT1D, EPPT1D et RGPT1D sont
!                compares aux valeurs issues du fichier suite. Une
!                correspondance exacte est necessaire pour continuer le calcul.
if (iappel.eq.2) then
   if(iphas.eq.1) then
      do ii = 1, nfpt1d
        ifac = ifpt1d(ii)
        nppt1d(ii) = 8
        eppt1d(ii) = 0.01144d0
        rgpt1d(ii) = 1.d0
        tppt1d(ii) = 25.d0
      enddo
   endif
endif


! --- Remplissage des conditions aux limites en paroi externe
!     ICLT1D(II) : type de condition a la limite
!                  = 1 : condition de dirichlet, avec coefficient d'echange
!                  = 3 : condition de flux
!     TEPT1D(II) : temperature exterieure
!     HEPT1D(II) : coefficient d'echange exterieur
!     FEPT1D(II) : flux applique a l'exterieur ( flux<0 = flux entrant)
!     XLMT1D(II) : coefficient de conductivite lambda de la paroi (W/m/°C)
!     RCPT1D(II) : coefficient rho*Cp de la paroi (J/m3/°C)
!     DTPT1D(II) : pas de temps de resolution de l'equation thermique dans
!                  la IIeme face de bord avec module thermique 1D (s)
if (iappel.eq.3) then
   if(iphas.eq.1) then
      do ii = 1, nfpt1d
         iclt1d(ii) = 1
!     parametres physiques
         ifac = ifpt1d(ii)
         if (cdgfbo(2,ifac).le.0.025d0) then
           iclt1d(ii) = 3
           fept1d(ii) = -1.d4
         else
           iclt1d(ii) = 3
           fept1d(ii) =  1.d4
         endif
         xlmt1d(ii) = 31.5d0
         rcpt1d(ii) = 3.5d6
         dtpt1d(ii) = 0.3d0
      enddo
   endif
endif

return

end

