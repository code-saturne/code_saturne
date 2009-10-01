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

subroutine usproj &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvep   , nivep  , ntersl , nvlsta , nvisbr , &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr , itepa  ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , ettpa  , tepa   , statis , stativ , tslagr , parbor , &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! MODIFICATION UTILISATEUR EN FIN DE PAS DE TEMPS
!   TOUT EST POSSIBLE,


! ON DONNE ICI PLUSIEURS EXEMPLES :

!  - CALCUL DE BILAN THERMIQUE
!    (au besoin, voir "ADAPTATION A UN SCALAIRE QUELCONQUE")

!  - CALCUL DES EFFORTS GLOBAUX SUR UN SOUS-ENSEMBLE DE FACES

!  - MODIFICATION ARBITRAIRE D'UNE VARIABLE DE CALCUL

!  - EXTRACTION D'UN PROFIL 1D

!  - IMPRESSION D'UN MOMENT

!  - EXEMPLES D'UTILISATION DES ROUTINES DE PARALLELISME

! CES EXEMPLES SONT DONNES EN SUPPOSANT UN CAS AVEC PERIODICITE
!  (IPERIO    .GT.0) ET PARALLELISME (IRANGP.GE.0).


! LE CALCUL DE BILAN THERMIQUE FOURNIT EN OUTRE UNE TRAME POUR
!  PLUSIEURS CHOSES
!  - CALCUL DE GRADIENT (AVEC LES PRECAUTIONS UTILES EN PARALLELE ET
!    PERIODIQUE)
!  - CALCUL DE GRANDEUR DEPENDANT DES VALEURS AUX CELLULES VOISINES
!    D'UNE FACE (AVEC LES PRECAUTIONS A PRENDRE EN PARALLELE ET
!    PERIODIQUE : VOIR L'ECHANGE DE DT ET DE CP)
!  - CALCUL D'UNE SOMME SUR LES PROCESSEURS LORS D'UN CALCUL
!    PARALLELE (PARSOM)


! Cells, boundary faces and interior faces identification
! =======================================================

! Cells, boundary faces and interior faces may be identified using
! the subroutines 'getcel', 'getfbr' and 'getfac' (respectively).

!  getfbr(string, nelts, eltlst):
!  - string is a user-supplied character string containing selection criteria;
!  - nelts is set by the subroutine. It is an integer value corresponding to
!    the number of boundary faces verifying the selection criteria;
!  - lstelt is set by the subroutine. It is an integer array of size nelts
!    containing the list of boundary faces verifying the selection criteria.

!  string may contain:
!  - references to colors (ex.: 1, 8, 26, ...)
!  - references to groups (ex.: inlet, group1, ...)
!  - geometric criteria (ex. x < 0.1, y >= 0.25, ...)
!  These criteria may be combined using logical operators ('and', 'or') and
!  parentheses.
!  Example: '1 and (group2 or group3) and y < 1' will select boundary faces
!  of color 1, belonging to groups 'group2' or 'group3' and with face center
!  coordinate y less than 1.

! Similarly, interior faces and cells can be identified using the 'getfac'
! and 'getcel' subroutines (respectively). Their syntax are identical to
! 'getfbr' syntax.

! For a more thorough description of the criteria syntax, it can be referred
! to the user guide.


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
! lndfac           ! e  ! <-- ! longueur du tableau nodfac (optionnel          !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr (optionnel          !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nbpmax           ! e  ! <-- ! nombre max de particules autorise              !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! ifacel           ! te ! <-- ! elements voisins d'une face interne            !
! (2, nfac)        !    !     !                                                !
! ifabor           ! te ! <-- ! element  voisin  d'une face de bord            !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! te ! <-- ! numero de famille d'une face de bord           !
! (nfabor)         !    !     !                                                !
! ifmcel           ! te ! <-- ! numero de famille d'une cellule                !
! (ncelet)         !    !     !                                                !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
! nfml  ,nprfml    !    !     !                                                !
! maxelt           !  e ! <-- ! nb max d'elements (cell,fac,fbr)               !
! lstelt(maxelt) te ! --- ! tableau de travail                             !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (nfac+1)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (lndfac)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!  (nfabor+1)      !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (lndfbr  )     !    !     !  (optionnel)                                   !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
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
! xyznod           ! tr ! <-- ! coordonnes des noeuds (optionnel)              !
! (ndim,nnod)      !    !     !                                                !
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet          !    !     !                                                !
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! rtp, rtpa        ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! statis           ! tr ! <-- ! moyennes statistiques                          !
!(ncelet,nvlsta    !    !     !                                                !
! stativ           ! tr ! <-- ! cumul pour les variances des                   !
!(ncelet,          !    !     !    statistiques volumiques                     !
!   nvlsta-1)      !    !     !                                                !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,ntersl    !    !     !   lagrangien sur la phase porteuse             !
! parbor           ! tr ! <-- ! infos sur interaction des particules           !
!(nfabor,nvisbr    !    !     !   aux faces de bord                            !
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

include "dimfbr.h"
include "paramx.h"
include "pointe.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "lagpar.h"
include "lagran.h"
include "parall.h"
include "period.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          itepa(nbpmax,nivep)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision statis(ncelet,nvlsta), stativ(ncelet,nvlsta-1)
double precision tslagr(ncelet,ntersl)
double precision parbor(nfabor,nvisbr)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          iel    , ielg   , ifac   , ifacg  , ivar
integer          iel1   , iel2   , ieltsm
integer          iortho , impout
integer          ifinia , ifinra
integer          igradx , igrady , igradz
integer          itravx , itravy , itravz , itreco
integer          inc    , iccocg
integer          nswrgp , imligp , iphydp , iwarnp
integer          iutile , iphas  , iclvar , iii
integer          ipcrom , ipcvst , iflmas , iflmab , ipccp, ipcvsl
integer          idimte , itenso , iscal
integer          ii     , nbr    , irangv , irang1 , npoint
integer          imom   , ipcmom , idtcm
integer          itab(3), iun
integer          ncesmp , icesmp , ismacp , itpsmp
integer          ilelt  , nlelt

double precision xrtpa  , xrtp
double precision xbilan , xbilvl , xbilpa , xbilpt
double precision xbilsy , xbilen , xbilso , xbildv
double precision xbilmi , xbilma
double precision epsrgp , climgp , extrap
double precision xfluxf , xgamma
double precision diipbx, diipby, diipbz, surfbn, distbr
double precision visct, flumab , xcp , xvsl, cp0iph, rrr
double precision xfor(3), xyz(3), xabs, xu, xv, xw, xk, xeps

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. INITIALISATION
!===============================================================================

! ---> Gestion memoire

idebia = idbia0
idebra = idbra0


!===============================================================================
! 2. EXEMPLE : CALCUL DE BILAN D'ENERGIE RELATIF A LA TEMPERATURE
!    ------------------------------------------------------------

!   On suppose que l'on souhaite faire des bilans (convectifs et
!     diffusifs) aux frontieres du domaine de calcul represente
!     ci-dessous (frontieres reperees par leurs couleurs).

!   Le scalaire considere est la temperature. On fera egalement
!     intervenir la chaleur massique (pour obtenir des bilans en Joules).



!   Domaine et couleurs
!   -------------------
!                  6
!      --------------------------
!      |                        |
!      |                        |
!   7  |           1            | 5
!      |     ^                  |
!      |     |                  |
!      --------------------------

!         2  3             4


!    2, 4, 7 : parois adiabatiques
!    6       : paroi  a temperature imposee
!    3       : entree
!    5       : sortie
!    1       : symetrie


!  ---------------------------------------------------------------------

!     POUR LA SIGNIFICATION PHYSIQUE DE CALCULS, IL EST UTILE D'ADOPTER
!       UN PAS DE TEMPS UNIFORME EN ESPACE (IDTVAR = 0 OU 1)
!     EN OUTRE, MEME EN SUITE DE CALCUL, LE BILAN EST FAUX SI INPDT0=1
!       (VISCT NON INITIALISE ET T(n-1) NON CONNU)

!  ---------------------------------------------------------------------

!     VARIABLE TEMPERATURE : IVAR = ISCA(ISCALT) (utiliser RTP(IEL,IVAR))

!  ---------------------------------------------------------------------

!     LE BILAN A L'INSTANT N VAUT :


!       n   iel=ncelet                n-1
!  BILAN  = SOMME { VOLUME(iel)*CP*ROM(iel)
!           iel=1
!                                      *(RTPA(iel,ivar)-RTP(iel,ivar)) }

!           ifac=nfabor
!         + SOMME { SURFBN(ifac)*DT(IFABOR(ifac))*CP
!           ifac=1

!              * [ VISLS0(ISCALT) + VISCT(IFABOR(ifac))/SIGMAS(ISCALT) ]
!              / DISTBR(ifac)
!              * [ COEFA(ifac,ICLVAR)
!                   + (COEFB(ifac,ICLVAR)-1.D0)*RTP(IFABOR(ifac,ivar)) ] }

!           ifac=nfabor
!         + SOMME { DT(IFABOR(ifac))*CP
!           ifac=1
!                              *RTP(IFABOR(ifac,ivar))*(-FLUMAB(ifac)) }


!  Le premier terme (nul en stationnaire) est negatif si la quantite
!    d'energie a decru dans le volume.
!  Les autres termes (convection, diffusion) sont positifs si la
!    quantite d'energie a augmente dans le volume par les
!    apports des conditions aux limites.

!  En regime stationnaire, un bilan positif indique donc
!   un gain d'energie.

!  ---------------------------------------------------------------------


!   AVEC ROM CALCULE PAR LA LOI REGISSANT LA MASSE VOLUMIQUE DANS USPHYV,
!     SOIT, PAR EXEMPLE :
!           n-1
!        ROM (iel) = P0 / [ RR * ( RTPA(iel,IVAR) + TKELV ) ]


!  ---------------------------------------------------------------------

!    CP ET LAMBDA/CP PEUVENT ETRE VARIABLES

!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------




!    ADAPTATION A UN SCALAIRE QUELCONQUE
!    -----------------------------------

!    L'approche peut s'utiliser pour realiser le bilan d'un scalaire
!      quelconque (mais les bilans ne sont plus attendus en Joules et
!      la chaleur massique n'intervient donc plus)

!    Pour cela :

!      - remplacer ISCALT(IPHAS) par le numero ISCAL du scalaire souhaite
!        ISCAL pouvant varier de 1 a NSCAL

!      - positionner IPCCP a 0 independamment de la valeur de ICP(IPHAS)
!        et affecter la valeur 1 a  CP0IPH (au lieu de CP0(IPHAS)).

!===============================================================================


!  Le bilan n'est pas valable si INPDT0=1
if (inpdt0.eq.0) then


! 2.1 INITIALISATION
! ==================

! --> Variables locales
!     -----------------

!       XBILVL : bilan volumique des termes instationnaires
!       XBILDV : bilan volumique du au terme en div(rho u)
!       XBILPA : bilan en paroi adiabatique
!       XBILPT : bilan en paroi a temperature imposee
!       XBILSY : bilan en symetrie
!       XBILEN : bilan en entree
!       XBILSO : bilan en sortie
!       XBILMI : bilan lie aux injections de masse
!       XBILMA : bilan lie aux aspirations de masse
!       XBILAN : bilan total

xbilvl = 0.d0
xbildv = 0.d0
xbilpa = 0.d0
xbilpt = 0.d0
xbilsy = 0.d0
xbilen = 0.d0
xbilso = 0.d0
xbilmi = 0.d0
xbilma = 0.d0
xbilan = 0.d0

! --- On travaillera sur la phase 1 uniquement
iphas = 1

! --- Le numero du scalaire temperature est ISCAL
iscal = iscalt(iphas)

! --- Le numero de variable temperature est IVAR
ivar = isca(iscal)

! --- Le numero pour les conditions aux limites est
iclvar = iclrtp(ivar,icoef)

! --- Numero des grandeurs physiques
ipcrom = ipproc(irom(iphas))
ipcvst = ipproc(ivisct(iphas))
iflmas = ipprof(ifluma(ivar))
iflmab = ipprob(ifluma(ivar))

! --- On stocke dans IPCCP un indicateur permettant de determiner si
!       la Chaleur massique est constante (=CP0) ou variable. Elle sera
!       utilisee pour realiser les bilans (XBILVL est en Joules).
if(icp(iphas).gt.0) then
  ipccp  = ipproc(icp   (iphas))
else
  ipccp  = 0
  cp0iph = cp0(iphas)
endif

! --- On stocke dans IPCVSL un indicateur permettant de determiner si
!       la diffusivite est constante (=VISLS0) ou variable. Elle sera
!       utilisee pour les termes diffusifs.
if(ivisls(iscal).gt.0) then
  ipcvsl = ipproc(ivisls(iscal))
else
  ipcvsl = 0
endif


! --> Echange de Cp et de Dt
!     ----------------------

!      Pour calculer les valeurs de flux aux faces internes, il est
!        necessaire d'avoir acces aux variables dans les cellules voisines.
!        En particulier, il faut connaitre la chaleur massique et la valeur
!        du pas de temps. Pour cela,
!        - dans les calculs paralleles, il est necessaire que
!          les cellules situees sur un bord de sous-domaine connaissent
!          la valeur de ces variables dans les cellules situees en
!          vis-a-vis sur le sous-domaine voisin.
!        - dans les calculs periodiques, il est necessaire que
!          les cellules periodiques aient acces a la valeur de ces0
!          variables dans les cellules periodiques correspondantes

!      Pour cela, il est necessaire d'appeler les routines de
!        communication PARCOM (parallelisme) et PERCOM (periodicite)
!        pour echanger les valeurs de Cp et de Dt avant de calculer le
!        gradient. L'appel de ces routines doit etre dait dans cet ordre
!        PARCOM puis PERCOM (si parallelisme et periodicite coexistent).

!      Si le calcul n'est ni periodique, ni parallele, on peut conserver
!        appels (les tests sur IPERIO et IRANGP assurent la generalite)



!    - Echange pour le parallelisme

  if(irangp.ge.0) then

!       Echange de Dt
    call parcom (dt)
    !==========

!       Echange de Cp si variable (sinon CP0(IPHAS) est utilise)
    if(ipccp.gt.0) then
      call parcom (propce(1,ipccp))
      !==========
    endif

  endif


!    - Echange pour la periodicite

  if(iperio.eq.1) then

    idimte = 0
    itenso = 0

!       Echange de Dt
    call percom                                                   &
    !==========
      ( idimte , itenso ,                                         &
        dt     , dt     , dt     ,                                &
        dt     , dt     , dt     ,                                &
        dt     , dt     , dt     )

!       Echange de Cp si variable (sinon CP0(IPHAS) est utilise)
    if(ipccp.gt.0) then
      call percom                                                 &
      !==========
      ( idimte , itenso ,                                         &
        propce(1,ipccp) , propce(1,ipccp) , propce(1,ipccp) ,     &
        propce(1,ipccp) , propce(1,ipccp) , propce(1,ipccp) ,     &
        propce(1,ipccp) , propce(1,ipccp) , propce(1,ipccp) )
    endif

  endif





! --> Calcul de la valeur reconstruite en I' pour les mailles de bord


!     Pour les maillages orthogonaux, elle doit etre egale
!        a la valeur au centre de la cellule
!     Cette valeur est calculee dans RA(ITRECO+IFAC-1)
!                                                   (avec IFAC=1,NFABOR)

!     Dans le cas de maillages orthogonaux, on peut simplifier :
!       il suffit d'affecter RTP(IEL,IVAR) a RA(ITRECO+IFAC-1),
!                                                  avec IEL=IFABOR(IFAC)
!       (cette option correspond a la deuxieme branche du if ci dessous,
!        avec IORTHO different de 0)


iortho = 0


! --> Cas des maillages non orthogonaux

if(iortho.eq.0) then

! --- Reservation de la memoire

  ifinia = idebia

  igradx = idebra
  igrady = igradx+ncelet
  igradz = igrady+ncelet
  itravx = igradz+ncelet
  itravy = itravx+ncelet
  itravz = itravy+ncelet
  itreco = itravz+ncelet
  ifinra = itreco+nfabor

! --- Verification de la disponibilite de la memoire

  CALL IASIZE('USPROJ',IFINIA)
  CALL RASIZE('USPROJ',IFINRA)


! --- Calcul du gradient de la temperature


!      Pour calculer le gradient de Temperature
!        - dans les calculs paralleles, il est necessaire que
!          les cellules situees sur un bord de sous-domaine connaissent
!          la valeur de temperature dans les cellules situees en
!          vis-a-vis sur le sous-domaine voisin.
!        - dans les calculs periodiques, il est necessaire que
!          les cellules periodiques aient acces a la valeur de la
!          temperature des cellules periodiques correspondantes

!      Pour cela, il est necessaire d'appeler les routines de
!        communication PARCOM (parallelisme) et PERCOM (periodicite)
!        pour echanger les valeurs de temperature avant de calculer le
!        gradient. L'appel a ces routines doit etre fait dans cet ordre
!        PARCOM puis PERCOM (pour les cas ou parallelisme et periodicite
!        coexistent).
!      En effet, on se situe ici a la fin du pas de temps n. Or,
!        les variables RTP ne seront echangees qu'en debut du pas de
!        temps n+1. Ici, seules les variables RTPA (obtenues a la fin
!        du pas de temps n-1) ont deja ete echangees.

!      Si le calcul n'est ni periodique, ni parallele, on peut conserver
!        appels (les tests sur IPERIO et IRANGP assurent la generalite)



!    - Echange pour le parallelisme

  if(irangp.ge.0) then

    call parcom (rtp(1,ivar))
    !==========

  endif

!    - Echange pour la periodicite

  if(iperio.eq.1) then

    idimte = 0
    itenso = 0
    call percom                                                   &
    !==========
      ( idimte , itenso ,                                         &
        rtp(1,ivar), rtp(1,ivar), rtp(1,ivar),                    &
        rtp(1,ivar), rtp(1,ivar), rtp(1,ivar),                    &
        rtp(1,ivar), rtp(1,ivar), rtp(1,ivar))

  endif


!    - Calcul du gradient

  inc = 1
  iccocg = 1
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  iphydp = 0

  call grdcel                                                     &
  !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   ra(itravx) , ra(itravx) , ra(itravx) ,                         &
   rtp(1,ivar) , coefa(1,iclvar) , coefb(1,iclvar) ,              &
   ra(igradx) , ra(igrady) , ra(igradz) ,                         &
!        ----------   ----------   ----------
   ra(itravx) , ra(itravy) , ra(itravz) ,                         &
   rdevel , rtuser , ra     )


!    - Calcul de la valeur reconstruite dans les cellules de bord

  do ifac = 1, nfabor
    iel = ifabor(ifac)
    iii = idiipb-1+3*(ifac-1)
    diipbx = ra(iii+1)
    diipby = ra(iii+2)
    diipbz = ra(iii+3)
    ra(itreco+ifac-1) = rtp(iel,ivar)                             &
         + diipbx*ra(igradx+iel-1)                                &
         + diipby*ra(igrady+iel-1)                                &
         + diipbz*ra(igradz+iel-1)
  enddo




! --> Cas des maillages orthogonaux

else

! --- Reservation de la memoire

  ifinia = idebia

  itreco = idebra
  ifinra = itreco+nfabor

! --- Verification de la disponibilite de la memoire

  CALL IASIZE('USPROJ',IFINIA)
  CALL RASIZE('USPROJ',IFINRA)

! ---  Calcul de la valeur reconstruite (en fait, ici, affectation
!     de la valeur non reconstruite)

  do ifac = 1, nfabor
    iel = ifabor(ifac)
    ra(itreco+ifac-1) = rtp(iel,ivar)
  enddo

endif





! 2.2 CALCUL DU BILAN A L'INSTANT N
! =================================

! --> Bilan sur les volumes internes
!     ------------------------------

!       La masse volumique ROM est celle qui a ete calculee en debut de
!       pas de temps a partir de la temperature au pas de temps
!       precedent (dans usphyv ou, si elle est constante RO0)


if(ipccp.gt.0) then
  do iel = 1, ncel
    xrtpa = rtpa(iel,ivar)
    xrtp  = rtp (iel,ivar)
    xbilvl = xbilvl                                               &
        + volume(iel) * propce(iel,ipccp) * propce(iel,ipcrom)    &
                                               * ( xrtpa - xrtp )
  enddo
else
  do iel = 1, ncel
    xrtpa = rtpa(iel,ivar)
    xrtp  = rtp (iel,ivar)
    xbilvl = xbilvl                                               &
        + volume(iel) * cp0iph * propce(iel,ipcrom)               &
                                               * ( xrtpa - xrtp )
  enddo
endif


! --> Bilan sur toutes les faces, internes et de bord, pour prendre
!     en compte le terme en div(rho u)-----------------------------
!     --------------------------------

!     Attention, on fait intervenir les valeurs de Cp et de Dt dans
!     les cellules voisines des faces internes, ce qui necessite d'avoir
!     pris ses precautions en periodicite et parallelisme.

! --- Si Cp est variable
!     Noter que si Cp est variable, ecrire ici un bilan sur l'equation
!     de la temerature n'est pas absolument correct

if(ipccp.gt.0) then
  do ifac = 1, nfac
    iel1 = ifacel(1,ifac)
    iel2 = ifacel(2,ifac)
    xbildv = xbildv + propfa(ifac,iflmas)                         &
          *(dt(iel1)*propce(iel1,ipccp)*rtp(iel1,ivar)            &
           -dt(iel2)*propce(iel2,ipccp)*rtp(iel2,ivar))
  enddo

  do ifac = 1, nfabor
    iel = ifabor(ifac)
    xbildv = xbildv + dt(iel)                                     &
         * propce(iel,ipccp)                                      &
         * propfb(ifac,iflmab)                                    &
         * rtp(iel,ivar)
  enddo

! --- Si Cp est constant

else
  do ifac = 1, nfac
    iel1 = ifacel(1,ifac)
    iel2 = ifacel(2,ifac)
    xbildv = xbildv + (dt(iel1)+ dt(iel2))*0.5d0                  &
         * cp0iph                                                 &
         * propfa(ifac,iflmas)                                    &
         * (rtp(iel1,ivar) - rtp(iel2,ivar))
  enddo

  do ifac = 1, nfabor
    iel = ifabor(ifac)
    xbildv = xbildv + dt(iel)                                     &
         * cp0iph                                                 &
         * propfb(ifac,iflmab)                                    &
         * rtp(iel,ivar)
  enddo
endif

!     En cas de terme source de masse, on ajoute la contribution en Gamma*Tn+1

ncesmp=ncetsm(iphas)
if (ncesmp.gt.0) then
  icesmp = iicesm(iphas)
  ismacp = ismace(iphas)
  itpsmp = iitpsm(iphas)
  do ieltsm = 1, ncesmp
    iel = ia(icesmp+ieltsm-1)
    xrtp  = rtp (iel,ivar)
    xgamma = ra( ismacp+ieltsm+ncesmp*(ipr(iphas)-1)-1)
    if(ipccp.gt.0) then
      xbildv = xbildv                                             &
           - volume(iel) * propce(iel,ipccp) * dt(iel)            &
           * xgamma * xrtp
    else
      xbildv = xbildv                                             &
           - volume(iel) * cp0iph * dt(iel)                       &
           * xgamma * xrtp
    endif
  enddo
endif



! --> Bilan sur les faces de bord
!     ---------------------------

!      On distingue ici les differents types de faces de bord
!        pour mieux analyser l'information, mais ce n'est pas oblige.


!     - Calcul de la contribution des parois de couleur 2, 4, 7
!         (ici adiabatiques, donc flux nul a priori)

CALL GETFBR('2 or 4 or 7',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

! ---   Element de bord

  iel    = ifabor(ifac)

! ---   Variables geometriques

  surfbn = ra(isrfbn-1+ifac)
  distbr = ra(idistb-1+ifac)

! ---   Variables physiques

  visct  = propce(iel,ipcvst)
  flumab = propfb(ifac,iflmab)

  if(ipccp.gt.0) then
    xcp = propce(iel,ipccp)
  else
    xcp    = cp0iph
  endif

  if(ipcvsl.gt.0) then
    xvsl = propce(iel,ipcvsl)
  else
    xvsl = visls0(iscal)
  endif

! ---   Calcul de la contribution au flux sur la facette courante
!         (flux de diffusion et de convection, negatif si entrant)

    xfluxf =          surfbn * dt(iel) * xcp *                    &
     (xvsl+visct/sigmas(iscal))/distbr *                          &
     (coefa(ifac,iclvar)+(coefb(ifac,iclvar)-1.d0)                &
                                               *ra(itreco+ifac-1))&
                    - flumab * dt(iel) * xcp *                    &
     (coefa(ifac,iclvar)+ coefb(ifac,iclvar)                      &
                                               *ra(itreco+ifac-1))

    xbilpa = xbilpa + xfluxf

enddo


!     - Calcul de la contribution des parois de couleur 6
!         (ici a temperature imposee ;
!                                   le flux convectif est nul a priori)

CALL GETFBR('6',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

! ---   Element de bord

  iel    = ifabor(ifac)

! ---   Variables geometriques

  surfbn = ra(isrfbn-1+ifac)
  distbr = ra(idistb-1+ifac)

! ---   Variables physiques

  visct  = propce(iel,ipcvst)
  flumab = propfb(ifac,iflmab)

  if(ipccp.gt.0) then
    xcp = propce(iel,ipccp)
  else
    xcp    = cp0iph
  endif

  if(ipcvsl.gt.0) then
    xvsl = propce(iel,ipcvsl)
  else
    xvsl = visls0(iscal)
  endif

! ---   Calcul de la contribution au flux sur la facette courante
!         (flux de diffusion et de convection, negatif si entrant)

    xfluxf =          surfbn * dt(iel) * xcp *                    &
     (xvsl+visct/sigmas(iscal))/distbr *                          &
     (coefa(ifac,iclvar)+(coefb(ifac,iclvar)-1.d0)                &
                                               *ra(itreco+ifac-1))&
                    - flumab * dt(iel) * xcp *                    &
     (coefa(ifac,iclvar)+ coefb(ifac,iclvar)                      &
                                               *ra(itreco+ifac-1))

    xbilpt = xbilpt + xfluxf

enddo


!     - Calcul de la contribution des symetries (couleur 1)
!         (a priori nul).

CALL GETFBR('1',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

! ---   Element de bord

  iel    = ifabor(ifac)

! ---   Variables geometriques

  surfbn = ra(isrfbn-1+ifac)
  distbr = ra(idistb-1+ifac)

! ---   Variables physiques

  visct  = propce(iel,ipcvst)
  flumab = propfb(ifac,iflmab)

  if(ipccp.gt.0) then
    xcp = propce(iel,ipccp)
  else
    xcp    = cp0iph
  endif

  if(ipcvsl.gt.0) then
    xvsl = propce(iel,ipcvsl)
  else
    xvsl = visls0(iscal)
  endif

! ---   Calcul de la contribution au flux sur la facette courante
!         (flux de diffusion et de convection, negatif si entrant)

    xfluxf =          surfbn * dt(iel) * xcp *                    &
     (xvsl+visct/sigmas(iscal))/distbr *                          &
     (coefa(ifac,iclvar)+(coefb(ifac,iclvar)-1.d0)                &
                                               *ra(itreco+ifac-1))&
                    - flumab * dt(iel) * xcp *                    &
     (coefa(ifac,iclvar)+ coefb(ifac,iclvar)                      &
                                               *ra(itreco+ifac-1))

    xbilsy = xbilsy + xfluxf

enddo


!     - Calcul de la contribution en entree (couleur 3)
!         (flux de diffusion et de convection)

CALL GETFBR('3',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

! ---   Element de bord

  iel    = ifabor(ifac)

! ---   Variables geometriques

  surfbn = ra(isrfbn-1+ifac)
  distbr = ra(idistb-1+ifac)

! ---   Variables physiques

  visct  = propce(iel,ipcvst)
  flumab = propfb(ifac,iflmab)

  if(ipccp.gt.0) then
    xcp = propce(iel,ipccp)
  else
    xcp    = cp0iph
  endif

  if(ipcvsl.gt.0) then
    xvsl = propce(iel,ipcvsl)
  else
    xvsl = visls0(iscal)
  endif

! ---   Calcul de la contribution au flux sur la facette courante
!         (flux de diffusion et de convection, negatif si entrant)

    xfluxf =          surfbn * dt(iel) * xcp *                    &
     (xvsl+visct/sigmas(iscal))/distbr *                          &
     (coefa(ifac,iclvar)+(coefb(ifac,iclvar)-1.d0)                &
                                               *ra(itreco+ifac-1))&
                    - flumab * dt(iel) * xcp *                    &
     (coefa(ifac,iclvar)+ coefb(ifac,iclvar)                      &
                                               *ra(itreco+ifac-1))

    xbilen = xbilen + xfluxf

enddo

!     - Calcul de la contribution en sortie (couleur 5)
!         (flux de diffusion et de convection)

CALL GETFBR('5',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

! ---   Element de bord

  iel    = ifabor(ifac)

! ---   Variables geometriques

  surfbn = ra(isrfbn-1+ifac)
  distbr = ra(idistb-1+ifac)

! ---   Variables physiques

  visct  = propce(iel,ipcvst)
  flumab = propfb(ifac,iflmab)

  if(ipccp.gt.0) then
    xcp = propce(iel,ipccp)
  else
    xcp    = cp0iph
  endif

  if(ipcvsl.gt.0) then
    xvsl = propce(iel,ipcvsl)
  else
    xvsl = visls0(iscal)
  endif

! ---   Calcul de la contribution au flux sur la facette courante
!         (flux de diffusion et de convection, negatif si entrant)

    xfluxf =          surfbn * dt(iel) * xcp *                    &
     (xvsl+visct/sigmas(iscal))/distbr *                          &
     (coefa(ifac,iclvar)+(coefb(ifac,iclvar)-1.d0)                &
                                               *ra(itreco+ifac-1))&
                    - flumab * dt(iel) * xcp *                    &
     (coefa(ifac,iclvar)+ coefb(ifac,iclvar)                      &
                                               *ra(itreco+ifac-1))

    xbilso = xbilso + xfluxf

enddo


! --> Bilan sur les termes sources de masse
!     -------------------------------------
!     On va séparer les injections de masse des aspirations
!     pour plus de generalite

ncesmp=ncetsm(iphas)
if (ncesmp.gt.0) then
  icesmp = iicesm(iphas)
  ismacp = ismace(iphas)
  itpsmp = iitpsm(iphas)
  do ieltsm = 1, ncesmp
!     suivant le type d'injection, on utilise la valeur SMACEL ou la valeur ambiante
!     de la température
    iel = ia(icesmp+ieltsm-1)
    xgamma = ra( ismacp+ieltsm+ncesmp*(ipr(iphas)-1)-1)
    if ( ia( itpsmp+ieltsm+ncesmp*(ivar-1)-1).eq.0                &
         .or. xgamma.lt.0.d0 ) then
      xrtp  = rtp (iel,ivar)
    else
      xrtp  = ra( ismacp+ieltsm+ncesmp*(ivar-1)-1)
    endif
    if(ipccp.gt.0) then
      if (xgamma.lt.0.d0) then
        xbilma = xbilma                                           &
             + volume(iel) * propce(iel,ipccp) * dt(iel)          &
             * xgamma * xrtp
      else
        xbilmi = xbilmi                                           &
             + volume(iel) * propce(iel,ipccp) * dt(iel)          &
             * xgamma * xrtp
      endif
    else
      if (xgamma.lt.0.d0) then
        xbilma = xbilma                                           &
             + volume(iel) * cp0iph * dt(iel)                     &
             * xgamma * xrtp
      else
        xbilmi = xbilmi                                           &
             + volume(iel) * cp0iph * dt(iel)                     &
             * xgamma * xrtp
      endif
    endif
  enddo
endif


!     - Somme des grandeurs sur tous les processeurs (calculs paralleles)

if (irangp.ge.0) then
  call parsom (xbilvl)
  call parsom (xbildv)
  call parsom (xbilpa)
  call parsom (xbilpt)
  call parsom (xbilsy)
  call parsom (xbilen)
  call parsom (xbilso)
  call parsom (xbilmi)
  call parsom (xbilma)
endif




! --> Bilan total
!     -----------

!      On ajoute les differentes contributions calculees plus haut.

xbilan = xbilvl + xbildv + xbilpa + xbilpt + xbilsy + xbilen      &
     + xbilso+ xbilmi + xbilma




! 2.3 ECRITURE DU BILAN A L'INSTANT N
! ===================================

write(nfecra,2000)                                                &
 ntcabs, xbilvl, xbildv, xbilpa, xbilpt, xbilsy, xbilen, xbilso,  &
     xbilmi, xbilma, xbilan


 2000 format(/,                                                   &
 3X,'** BILAN THERMIQUE **',/,                              &
 3X,'   ---------------',/,                                 &
 '---','------',                                                  &
 '------------------------------------------------------------',/,&
 'bt ','  ITER',                                                  &
 '   Volumique  Divergence  Paroi Adia  Paroi Timp    Symetrie',  &
 '      Entree      Sortie  Masse inj.  Masse asp.  Bil. total',/,&
 'bt ',I6,10E12.4,/,                                        &
 '---','------',                                                  &
 '------------------------------------------------------------')

!- Fin du test sur INPDT0
endif

!===============================================================================
! 3. EXEMPLE : CALCUL DES EFFORTS GLOBAUX SUR UN SOUS-ENSEMBLE DE FACES

!           A FAIRE AVEC PRECAUTIONS ...
!           L'UTILISATEUR PREND SES RESPONSABILITES.
!===============================================================================


! ----------------------------------------------

! Il est assez courant que l'on oublie d'eliminer cet exemple
!   de ce sous-programme.
! On a donc prevu le test suivant pour eviter les mauvaises surprises

iutile = 0

if(iutile.eq.0) return

! ----------------------------------------------

!     Si les efforts ont bien ete calcules :
if (ineedf.eq.1) then

  do ii = 1, ndim
    xfor(ii) = 0.d0
  enddo

  CALL GETFBR('2 or 3',NLELT,LSTELT)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)

    do ii = 1, ndim
      xfor(ii) = xfor(ii) + ra(iforbr + (ifac-1)*ndim + ii-1)
    enddo

  enddo

  if (irangp.ge.0) then
    call parrsm(ndim,xfor)
  endif

endif

!===============================================================================
! 4. EXEMPLE : MISE A 20 DE LA TEMPERATURE DANS UNE ZONE DONNEE
!              A PARTIR DU TEMPS 12s

!           A FAIRE AVEC PRECAUTIONS ...
!           L'UTILISATEUR PREND SES RESPONSABILITES.
!===============================================================================


! ----------------------------------------------

! Il est assez courant que l'on oublie d'eliminer cet exemple
!   de ce sous-programme.
! On a donc prevu le test suivant pour eviter les mauvaises surprises

iutile = 0

if(iutile.eq.0) return

! ----------------------------------------------


iphas = 1
iscal = iscalt(iphas)

if (ttcabs.ge.12.d0) then

  if (iscal.gt.0.and.iscal.le.nscal) then
    do iel = 1, ncel
      rtp(iel,isca(iscal)) = 20.d0
    enddo
  endif

  write(nfecra,3000)

endif

 3000 format(/,                                                   &
 ' MODIFICATION UTILISATEUR DES VARIABLES EN FIN DE PAS DE TEMPS',&
 /)

!===============================================================================
! 5. EXEMPLE : EXTRACTION D'UN PROFIL 1D
!    -----------------------------------

!    On cherche ici a extraire le profil de U, V, W, k et eps sur une
!     courbe 1D quelconque en fonction de l'abscisse curviligne.
!    Le profil est ecrit dans le fichier "profil.dat" (ne pas oublier de
!     prevoir son rapatriement par le script de lancement).

!    - la courbe utilisee ici est le segment [(0;0;0),(0;0.1;0)], mais la
!      generalisation a une courbe quelconque est simple.
!    - la routine est adaptee au parallelisme et a la periodicite, ainsi
!      qu'au differents modeles de turbulence.
!    - la courbe 1D est discretisee en NPOINT points. Pour chacun de ces
!      points, on calcule le centre de cellule le plus proche et on
!      imprime la valeur des variables au centre de cette cellule. Pour plus
!      de coherence, l'ordonnee imprimee est celle du centre de la cellule.
!    - on evite d'imprimer deux fois la meme cellule (si plusieurs points de
!      la courbe sont associes a la meme cellule).
!===============================================================================


! ----------------------------------------------

! Il est assez courant que l'on oublie d'eliminer cet exemple
!   de ce sous-programme.
! On a donc prevu le test suivant pour eviter les mauvaises surprises

iutile = 0

if(iutile.eq.0) return

! ----------------------------------------------


if (ntcabs.eq.ntmabs) then

!     Seul le processeur de rang 0 (parallele) ou -1 (scalaire) ecrit dans le
!     fichier. On utilise les unites "utilisateur".
  impout = impusr(1)
  if (irangp.le.0) then
    open(impout,file="profil.dat")
    write(impout,*)                                               &
         '# z(m) U(m/s) V(m/s) W(m/s) k(m2/s2) eps(m2/s3)'
  endif

  iphas  = 1
  npoint = 200
  iel1   = -999
  irang1 = -999
  do ii = 1, npoint

    xyz(1) = 0.d0
    xyz(2) = float(ii-1)/float(npoint-1)*0.1d0
    xyz(3) = 0.d0

    call findpt                                                   &
    !==========
      (ncelet, ncel, xyzcen,                                      &
       xyz(1), xyz(2), xyz(3), iel, irangv)

    if ((iel.ne.iel1).or.(irangv.ne.irang1)) then
      iel1   = iel
      irang1 = irangv

!     On remplit les variables temporaires XU, XV, ... pour le processeur
!      qui contient le point et on envoie ensuite l'info aux autres
!      processeurs.
      if (irangp.eq.irangv) then
        xabs = xyzcen(2,iel)
        xu   = rtp(iel,iu(iphas))
        xv   = rtp(iel,iv(iphas))
        xw   = rtp(iel,iw(iphas))
        xk   = 0.d0
        xeps = 0.d0
        if (itytur(iphas).eq.2 .or. iturb(iphas).eq.50            &
             .or. iturb(iphas).eq.60) then
          xk = rtp(iel,ik(iphas))
        elseif (itytur(iphas).eq.3) then
          xk = ( rtp(iel,ir11(iphas))+rtp(iel,ir22(iphas))+       &
               rtp(iel,ir33(iphas)) )/2.d0
        endif
        if (itytur(iphas).eq.2 .or. itytur(iphas).eq.3            &
             .or. iturb(iphas).eq.50) then
          xeps = rtp(iel,iep(iphas))
        elseif (iturb(iphas).eq.60) then
          xeps = cmu*rtp(iel,ik(iphas))*rtp(iel,iomg(iphas))
        endif
      else
        xabs = 0.d0
        xu   = 0.d0
        xv   = 0.d0
        xw   = 0.d0
        xk   = 0.d0
        xeps = 0.d0
      endif

!           Envoi aux autres processeurs si parallele
      if (irangp.ge.0) then
        iun = 1
        call parbcr(irangv,iun,xabs)
        call parbcr(irangv,iun,xu  )
        call parbcr(irangv,iun,xv  )
        call parbcr(irangv,iun,xw  )
        call parbcr(irangv,iun,xk  )
        call parbcr(irangv,iun,xeps)
      endif

      if (irangp.le.0)                                            &
           write(impout,99) xabs,xu,xv,xw,xk,xeps

 99         format(6g17.9)

    endif

  enddo

  if (irangp.le.0) close(impout)

endif


!===============================================================================
! 6. EXEMPLE : IMPRESSION DU PREMIER MOMENT CALCULE
!===============================================================================


! ----------------------------------------------

! Il est assez courant que l'on oublie d'eliminer cet exemple
!   de ce sous-programme.
! On a donc prevu le test suivant pour eviter les mauvaises surprises

iutile = 0

if(iutile.eq.0) return

! ----------------------------------------------


if(nbmomt.gt.0) then

!     Numero du moment : IMOM
  imom = 1

!     Position dans PROPCE du tableau de cumul temporel des moments
!       PROPCE(IEL,IPCMOM)
  ipcmom = ipproc(icmome(imom))

!     Le cumul temporel des moments doit etre divise par la variable
!       de cumul du temps qui est un tableau NCEL ou un reel :
!             un tableau NCEL   si IDTMOM(IMOM) > 0 : PROPCE(IEL,IDTCM)
!             ou un simple reel si IDTMOM(IMOM) < 0 : DTCMOM(IDTCM)

  if(idtmom(imom).gt.0) then
    idtcm = ipproc(icdtmo(idtmom(imom)))
    do iel = 1, ncel
      write(nfecra,4000) iel,propce(iel,ipcmom)/                  &
           max(propce(iel,idtcm),epzero)
    enddo
  elseif(idtmom(imom).lt.0) then
    idtcm = -idtmom(imom)
    do iel = 1, ncel
      write(nfecra,4000) iel,propce(iel,ipcmom)/                  &
           max(dtcmom(idtcm),epzero)
    enddo
  endif

endif

 4000 format(' Cellule ',I10,'   Premier moment ',E14.5)


!===============================================================================
! 7. EXEMPLE : UTILISATION DES ROUTINES DE CALCUL PARALLELE
!              POUR LES OPERATIONS SUIVANTES~:
!===============================================================================

!   Cet exemple n'a pas d'autre utilite que de fournir la liste des
!     routines utilisables pour simplifier certaines operations
!     globales en parallele.

!   ATTENTION, ces routines modifient leur argument


! ----------------------------------------------

! Il est assez courant que l'on oublie d'eliminer cet exemple
!   de ce sous-programme.
! On a donc prevu le test suivant pour eviter les mauvaises surprises

iutile = 0

if(iutile.eq.0) return

! ----------------------------------------------


! Maximum d'un compteur entier II, ici le nombre de cellules par processeur

!     Valeur locale
ii = ncel
!     Calcul du maximum sur les processeurs
if (irangp.ge.0) then
  call parcmx(ii)
endif
!     Ecriture de la valeur maximale renvoyee
write(nfecra,5010)ii
 5010 format(' USPROJ: Nombre de cellules max par processeur = ', I10)


! Somme d'un compteur entier II, ici le nombre de cellules

!     Valeur locale
ii = ncel
!     Calcul de la somme sur les processeurs
if (irangp.ge.0) then
  call parcpt(ii)
endif
!     Ecriture de la valeur somme renvoyee
write(nfecra,5020)ii
 5020 format(' USPROJ: Nombre de cellules total = ', I10)

! Somme d'un reel RRR, ici le volume

!     Valeur locale
rrr = 0.d0
do iel = 1, ncel
  rrr = rrr + volume(iel)
enddo
!     Calcul de la somme sur les processeurs
if (irangp.ge.0) then
  call parsom(rrr)
endif
!     Ecriture de la valeur somme renvoyee
write(nfecra,5030)rrr
 5030 format(' USPROJ: Volume du domaine total = ', E14.5)

! Maximum d'un reel RRR, ici le volume par processeur

!     Valeur locale
rrr = 0.d0
do iel = 1, ncel
  rrr = rrr + volume(iel)
enddo
!     Calcul du maximum sur les processeurs
if (irangp.ge.0) then
  call parmax(rrr)
endif
!     Ecriture de la valeur maximale renvoyee
write(nfecra,5040)rrr
 5040 format(' USPROJ: Volume max par processeur = ', E14.5)

! Minimum d'un reel RRR, ici le volume par processeur

!     Valeur locale
rrr = 0.d0
do iel = 1, ncel
  rrr = rrr + volume(iel)
enddo
!     Calcul du minimum sur les processeurs
if (irangp.ge.0) then
  call parmin(rrr)
endif
!     Ecriture de la valeur minimale renvoyee
write(nfecra,5050)rrr
 5050 format(' USPROJ: Volume min par processeur = ', E14.5)

!  Maximum d'un réel et valeurs reelles associées
!       ici le volume et sa localisation

!     Maximum local et sa localisation (NBR=3 coordonnees)
nbr = 3
rrr  = -1.d0
xyz(1) = 0.d0
xyz(2) = 0.d0
xyz(3) = 0.d0
do iel = 1, ncel
  if (rrr.lt.volume(iel)) then
    rrr = volume(iel)
    xyz(1) = xyzcen(1,iel)
    xyz(2) = xyzcen(2,iel)
    xyz(3) = xyzcen(3,iel)
  endif
enddo
!     Calcul du maximum et localisation sur les processeurs
if (irangp.ge.0) then
  call parmxl(nbr,rrr,xyz)
endif
!     Ecriture de la valeur maximale et localisation renvoyees
write(nfecra,5060)rrr,xyz(1),xyz(2),xyz(3)
 5060 format(' USPROJ: Volume max = ', E14.5,/,                   &
       '         Localisation x,y,z = ',3E14.5)

!  Minimum d'un réel et valeurs reelles associées
!       ici le volume et sa localisation

!     Minimum local et sa localisation (NBR=3 coordonnees)
nbr = 3
rrr  = 1.d+30
xyz(1) = 0.d0
xyz(2) = 0.d0
xyz(3) = 0.d0
do iel = 1, ncel
  if (rrr.gt.volume(iel)) then
    rrr = volume(iel)
    xyz(1) = xyzcen(1,iel)
    xyz(2) = xyzcen(2,iel)
    xyz(3) = xyzcen(3,iel)
  endif
enddo
!     Calcul du minimum et localisation sur les processeurs
if (irangp.ge.0) then
  call parmnl(nbr,rrr,xyz)
endif
!     Ecriture de la valeur minimale et localisation renvoyees
write(nfecra,5070)rrr,xyz(1),xyz(2),xyz(3)
 5070 format(' USPROJ: Volume min = ', E14.5,/,                   &
       '         Localisation x,y,z = ',3E14.5)

!  Somme d'un tableau d'entiers ici le nombre
!         de cellules, de faces et de faces de bord

!     Valeurs locales
nbr = 3
itab(1) = ncel
itab(2) = nfac
itab(3) = nfabor
!     Calcul de la somme sur les processeurs
if (irangp.ge.0) then
  call parism(nbr,itab)
endif
!     Ecriture de la valeur somme renvoyee
write(nfecra,5080)itab(1),itab(2),itab(3)
 5080 format(' USPROJ: Nombre de cellules       = ',I10,/,        &
       '         Nombre de faces internes = ',I10,/,        &
       '         Nombre de faces de bord  = ',I10)

!  Maximum d'un tableau d'entiers ici le nombre
!         de cellules, de faces et de faces de bord

!     Valeurs locales
nbr = 3
itab(1) = ncel
itab(2) = nfac
itab(3) = nfabor
!     Calcul du maximum sur les processeurs
if (irangp.ge.0) then
  call parimx(nbr,itab)
endif
!     Ecriture de la valeur maximale renvoyee
write(nfecra,5090)itab(1),itab(2),itab(3)
 5090 format(' USPROJ: Nombre de cellules max       par proc = ',I10,/, &
       '         Nombre de faces internes max par proc = ',I10,/, &
       '         Nombre de faces de bord  max par proc = ',I10)

!  Minimum d'un tableau d'entiers ici le nombre
!         de cellules, et de faces de bord
!     attention, une somme similaire pour compter les faces internes
!         compte 2 fois les faces internes situes sur les bords des
!         processeurs

!     Valeurs locales
nbr = 2
itab(1) = ncel
itab(2) = nfabor
!     Calcul du minimum sur les processeurs
if (irangp.ge.0) then
  call parimn(nbr,itab)
endif
!     Ecriture de la valeur minimale renvoyee
write(nfecra,5100)itab(1),itab(2)
 5100 format(' USPROJ: Nombre de cellules       min par proc = ',I10,/, &
       '         Nombre de faces de bord  min par proc = ',I10)

!  Somme d'un tableau de reels ici les trois composantes de la vitesse
!         (dans le but de faire une moyenne ensuite par exemple)

!     Valeurs locales
nbr = 3
xyz(1) = 0.d0
xyz(2) = 0.d0
xyz(3) = 0.d0
do iel = 1, ncel
  xyz(1) = xyz(1)+rtp(iel,iu(1))
  xyz(2) = xyz(2)+rtp(iel,iv(1))
  xyz(3) = xyz(3)+rtp(iel,iw(1))
enddo
!     Calcul de la somme sur les processeurs
if (irangp.ge.0) then
  call parrsm(nbr,xyz)
endif
!     Ecriture de la valeur somme renvoyee
write(nfecra,5110)xyz(1),xyz(2),xyz(3)
 5110 format(' USPROJ: Somme de U sur le domaine = ',E14.5,/,     &
       '         Somme de V sur le domaine = ',E14.5,/,     &
       '         Somme de W sur le domaine = ',E14.5)

!  Maximum d'un tableau de reels ici les trois composantes de la vitesse

!     Valeurs locales
nbr = 3
xyz(1) = rtp(1,iu(1))
xyz(2) = rtp(1,iv(1))
xyz(3) = rtp(1,iw(1))
do iel = 1, ncel
  xyz(1) = max(xyz(1),rtp(iel,iu(1)))
  xyz(2) = max(xyz(2),rtp(iel,iv(1)))
  xyz(3) = max(xyz(3),rtp(iel,iw(1)))
enddo
!     Calcul du maximum sur les processeurs
if (irangp.ge.0) then
  call parrmx(nbr,xyz)
endif
!     Ecriture de la valeur somme renvoyee
write(nfecra,5120)xyz(1),xyz(2),xyz(3)
 5120 format(' USPROJ: Maximum de U sur le domaine = ',E14.5,/,   &
       '         Maximum de V sur le domaine = ',E14.5,/,   &
       '         Maximum de W sur le domaine = ',E14.5)

!  Minimum d'un tableau de reels ici les trois composantes de la vitesse

!     Valeurs locales
nbr = 3
xyz(1) = rtp(1,iu(1))
xyz(2) = rtp(1,iv(1))
xyz(3) = rtp(1,iw(1))
do iel = 1, ncel
  xyz(1) = min(xyz(1),rtp(iel,iu(1)))
  xyz(2) = min(xyz(2),rtp(iel,iv(1)))
  xyz(3) = min(xyz(3),rtp(iel,iw(1)))
enddo
!     Calcul du minimum sur les processeurs
if (irangp.ge.0) then
  call parrmn(nbr,xyz)
endif
!     Ecriture de la valeur somme renvoyee
write(nfecra,5130)xyz(1),xyz(2),xyz(3)
 5130 format(' USPROJ: Minimum de U sur le domaine = ',E14.5,/,   &
       '         Minimum de V sur le domaine = ',E14.5,/,   &
       '         Minimum de W sur le domaine = ',E14.5)

!  Envoi d'un tableau de valeurs entieres locales aux autres processeurs
!     par exemple le nombre de cellules, de faces internes et de faces
!     de bord du processeur numero IRANGV=0.

!     Valeurs locales
irangv = 0
nbr = 3
itab(1) = ncel
itab(2) = nfac
itab(3) = nfabor
!     Envoi aux autres
if (irangp.ge.0) then
  call parbci(irangv,nbr,itab)
endif
!     Ecriture par chaque processeur de la valeur envoyee par le proc 0
write(nfecra,5140)irangv,itab(1),itab(2),itab(3)
 5140 format(' USPROJ: Sur le processeur ', I10 ,/,               &
       '         Nombre de cellules       = ',I10,/,        &
       '         Nombre de faces internes = ',I10,/,        &
       '         Nombre de faces de bord  = ',I10)

!  Envoi d'un tableau de valeurs reelles locales aux autres processeurs
!     par exemple trois valeurs de la vitesse
!     du processeur numero IRANGV=0.

!     Valeurs locales
irangv = 0
nbr = 3
xyz(1) = rtp(1,iu(1))
xyz(2) = rtp(1,iv(1))
xyz(3) = rtp(1,iw(1))
!     Envoi aux autres
if (irangp.ge.0) then
  call parbcr(irangv,nbr,xyz)
endif
!     Ecriture par chaque processeur de la valeur envoyee par le proc 0
write(nfecra,5150)irangv,xyz(1),xyz(2),xyz(3)
 5150 format(' USPROJ: Sur le processeur ', I10 ,/,               &
       '         Vitesse U dans la premiere cellule = ',E14.5,/,  &
       '         Vitesse V dans la premiere cellule = ',E14.5,/,  &
       '         Vitesse W dans la premiere cellule = ',E14.5)

!  Recuperation par tous les processeurs du numero global d'une cellule
!     appartenant a un processeur donne (exemple de la cellule IEL
!     du processeur IRANGV)
iel = 1
irangv = 0
if (irangp.ge.0) then
  call parcel(iel,irangv,ielg)
else
  ielg = -1
endif
!     Ecriture par chaque processeur du numero global de la cellule IEL
!       du processeur IRANGV
write(nfecra,5160)iel,irangv,ielg
 5160 format(' USPROJ: La cellule de numero local IEL    = ',I10,/,     &
       '            sur le processeur       IRANGV = ',I10,/,     &
       '            a pour numero global    IELG   = ',I10)

!  Recuperation par le processeur courant du numero global d'une cellule
!     IEL lui appartenant (exemple avec IEL=1, tous les processeurs ayant
!     au moins une cellule)
iel = 1
call parclg(iel,irangp,ielg)
!     Ecriture par chaque processeur du numero global de sa cellule IEL
!       (si le processeur a moins de IEL cellules, affiche 0)
!       (en sequentiel, affiche IEL)
write(nfecra,5170)iel,irangp,ielg
 5170 format(' USPROJ: La cellule de numero local IEL    = ',I10,/,     &
       '            sur le processeur       IRANGP = ',I10,/,     &
       '            a pour numero global    IELG   = ',I10)

!  Recuperation par le processeur courant du numero global d'une face
!     interne IFAC lui appartenant (exemple avec IFAC=1)
ifac = 1
call parfig(ifac,irangp,ifacg)
!     Ecriture par chaque processeur du numero global de sa face IFAC
!       (si le processeur a moins de IFAC faces internes, affiche 0)
!       (en sequentiel, affiche IFAC)
write(nfecra,5180)ifac,irangp,ifacg
 5180 format(' USPROJ: La face interne de numero local IFAC = ',I10,/,  &
       '            sur le processeur          IRANGP = ',I10,/,  &
       '            a pour numero global        IFACG = ',I10)

!  Recuperation par le processeur courant du numero global d'une face
!     de bord IFAC lui appartenant (exemple avec IFAC=1)
ifac = 1
call parfbg(ifac,irangp,ifacg)
!     Ecriture par chaque processeur du numero global de sa face IFAC
!       (si le processeur a moins de IFAC faces de bord , affiche 0)
!       (en sequentiel, affiche IFAC)
write(nfecra,5190)ifac,irangp,ifacg
 5190 format(' USPROJ: La face de bord de numero local IFAC = ',I10,/,  &
       '            sur le processeur          IRANGP = ',I10,/,  &
       '            a pour numero global        IFACG = ',I10)

return
end subroutine
