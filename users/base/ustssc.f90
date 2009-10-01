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

subroutine ustssc &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iscal  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr , icepdc , icetsm , itypsm , &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   crvexp , crvimp ,                                              &
   viscf  , viscb  , xam    ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     , w10    , w11    ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE UTILISATEUR
!   ON PRECISE LES TERMES SOURCES UTILISATEURS
!   POUR UN SCALAIRE SUR UN PAS DE TEMPS

! ON RESOUT RHO*VOLUME*D(VAR)/DT = CRVIMP*VAR + CRVEXP

! ON FOURNIT ICI CRVIMP ET CRVEXP (ILS CONTIENNENT RHO*VOLUME)
!    CRVEXP en kg variable/s :
!     ex : pour les temperatures      kg degres/s
!          pour les enthalpies        Joules/s
!    CRVIMP en kg /s :

! VEILLER A UTILISER UN CRVIMP NEGATIF
! (ON IMPLICITERA CRVIMP
!  IE SUR LA DIAGONALE DE LA MATRICE, LE CODE AJOUTERA :
!   MAX(-CRVIMP,0) EN SCHEMA STANDARD EN TEMPS
!       -CRVIMP    SI LES TERMES SOURCES SONT A L'ORDRE 2

! CES TABLEAUX SONT INITIALISES A ZERO AVANT APPEL A CE SOUS
!   PROGRAMME ET AJOUTES ENSUITE AUX TABLEAUX PRIS EN COMPTE
!   POUR LA RESOLUTION

! EN CAS D'ORDRE 2 DEMANDE SUR LES TERMES SOURCES, ON DOIT
!   FOURNIR CRVEXP A L'INSTANT N     (IL SERA EXTRAPOLE) ET
!           CRVIMP A L'INSTANT N+1/2 (IL EST  DANS LA MATRICE,
!                                     ON LE SUPPOSE NEGATIF)

! L'identification des faces de bord concernees se fait grace
! a la commande GETFBR.


! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
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
! ncepdp           ! e  ! <-- ! nombre de cellules avec pdc                    !
! ncesmp           ! e  ! <-- ! nombre de cellules a source de masse           !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! iscal            ! e  ! <-- ! numero du scalaire                             !
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
!   (lndfac)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !  (optionnel)                                   !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
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
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! crvexp(ncelet    ! tr ! --> ! tableau de travail pour part explicit          !
! crvimp(ncelet    ! tr ! --> ! tableau de travail pour terme instat           !
! viscf(nfac)      ! tr ! --- ! tableau de travail    faces internes           !
! viscb(nfabor     ! tr ! --- ! tableau de travail    faces de bord            !
! xam(nfac,2)      ! tr ! --- ! tableau de travail    faces de bord            !
! w1..11(ncelet    ! tr ! --- ! tableau de travail    cellules                 !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "pointe.h"
include "numvar.h"
include "entsor.h"
include "optcal.h"
include "cstphy.h"
include "parall.h"
include "period.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          nideve , nrdeve , nituse , nrtuse
integer          iscal

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision xam(nfac,2)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision w10(ncelet), w11(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

character*80     chaine
integer          idebia, idebra
integer          ivar, iiscvr, ipcrom, iel, iphas, iutile
integer          ilelt, nlelt

double precision tauf, prodf, volf, pwatt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

! --- Numero du scalaire a traiter : ISCAL

!     ISCAL est une donnee d'entree de ce sous programme
!       qui indique quel scalaire on est en train de traiter
!       en effet ce sous programme est appele successivement pour
!       tous les scalaires du calcul.
!     ISCAL ne doit pas etre modifie par l'utilisateur

!     Si le calcul comporte plusieurs scalaires (de 1 a n par exemple)
!       et que l'on souhaite imposer des termes sources pour un
!       scalaire p donne, on doit utiliser ISCAL ici :
!       CRVIMP et CRVEXP ne seront renseignes que "IF(ISCAL.EQ.p)"


! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

! --- Nom de la variable associee au scalaire a traiter ISCAL
chaine = nomvar(ipprtp(ivar))

! --- Indicateur de variance
!         Si ISCAVR = 0 :
!           le scalaire ISCAL n'est pas une variance
!         Si ISCAVR > 0 et ISCAVR < NSCAL + 1 :
!           le scalaire ISCAL est une variance associee
!           au scalaire ISCAVR
iiscvr = iscavr(iscal)

! --- Numero de phase associee au scalaire ISCAL
iphas = iphsca(iscal)

! --- Numero des grandeurs physiques (voir usclim)
ipcrom = ipproc(irom(iphas))

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif


!===============================================================================
! 2. EXEMPLE DE TERME SOURCE ARBITRAIRE, POUR LA VARIABLE F :

!                             S = A * F + B

!            APPARAISSANT DANS LES EQUATIONS SOUS LA FORME :

!                       RHO VOLUME D(F)/Dt = VOLUME*S


!   CE TERME A UNE PARTIE QU'ON VEUT IMPLICITER         : A
!           ET UNE PARTIE QU'ON VA TRAITER EN EXPLICITE : B


!   ICI PAR EXEMPLE ARBITRAIREMENT :

!     A = - RHO / TAUF
!     B =   RHO * PRODF
!        AVEC
!     TAUF   = 10.D0  [secondes  ] (TEMPS DE DISSIPATION DE F)
!     PRODF  = 100.D0 [variable/s] (PRODUCTION DE F PAR UNITE DE TEMPS)

!   ON A ALORS
!     CRVIMP(IEL) = VOLUME(IEL)* A = - VOLUME(IEL) (RHO / TAUF )
!     CRVEXP(IEL) = VOLUME(IEL)* B =   VOLUME(IEL) (RHO * PRODF)

!   ON REMPLIT CI-DESSOUS CRVIMP ET CRVEXP CORRESPONDANTS.

!===============================================================================



!     ATTENTION, L'EXEMPLE EST COMPLETEMENT ARBITRAIRE
!     =========
!         ET DOIT ETRE REMPLACE PAR LES TERMES UTILISATEURS ADEQUATS


! ----------------------------------------------

! Il est assez courant que l'on oublie d'eliminer cet exemple
!   de ce sous-programme.
! On a donc prevu le test suivant pour eviter les mauvaises surprises

iutile = 0

if(iutile.eq.0) return

! ----------------------------------------------

tauf  = 10.d0
prodf = 100.d0

do iel = 1, ncel
  crvimp(iel) = - volume(iel)*propce(iel,ipcrom)/tauf
enddo

do iel = 1, ncel
  crvexp(iel) =   volume(iel)*propce(iel,ipcrom)*prodf
enddo

!===============================================================================
! 3. EXEMPLE DE TERME SOURCE ARBITRAIRE, POUR LA VARIABLE F = ENTHALPIE :

!     IL S'AGIT D'UNE PARTICULARISATION DE L'EXEMPLE PRECEDENT AVEC

!                             S = B

!            APPARAISSANT DANS LES EQUATIONS SOUS LA FORME :

!                       RHO VOLUME D(F)/Dt = VOLUME*S


!   IL N'Y A RIEN A IMPLICITER, ON PEUT DONC IMPOSER

!     CRVIMP(IEL) = 0.D0


!   CE TERME A UNE PARTIE QU'ON VA TRAITER EN EXPLICITE : B

!     SI F EST UNE ENTHALPIE (J/kg), ON AURA ALORS B EN Watt/m3

!     SI ON CONNAIT LA PUISSANCE PWATT A INTRODUIRE UNIFORMEMENT DANS UN
!       VOLUME DONNE VOLF TEL QUE 0<X<1.2, 3.1<Y<4, ON POURRA ALORS
!       ECRIRE, POUR LES CELLULES DU MAILLAGE CONCERNEES :

!       CRVEXP(IEL) = VOLUME(IEL)* B =   VOLUME(IEL)*(PWATT/VOLF)


!===============================================================================



!     ATTENTION, L'EXEMPLE EST COMPLETEMENT ARBITRAIRE
!     =========
!         ET DOIT ETRE REMPLACE PAR LES TERMES UTILISATEURS ADEQUATS


! ----------------------------------------------

! Il est assez courant que l'on oublie d'eliminer cet exemple
!   de ce sous-programme.
! On a donc prevu le test suivant pour eviter les mauvaises surprises

iutile = 0

if(iutile.eq.0) return

! ----------------------------------------------

! PUISSANCE A INTRODUIRE DANS VOLF
!   ATTENTION on suppose qu'on travaille en enthalpie,
!             si on travaille en temperature, il ne faut pas oublier de
!               diviser PWATT par Cp.

pwatt = 100.d0

! CALCUL DE VOLF

volf  = 0.d0
CALL GETCEL('X > 0.0 and X < 1.2 and Y > 3.1 and'//               &
            'Y < 4.0',NLELT,LSTELT)

do ilelt = 1, nlelt
  iel = lstelt(ilelt)
  volf = volf + volume(iel)
enddo

do ilelt = 1, nlelt
  iel = lstelt(ilelt)
! PAS DE TERME IMPLICITE
  crvimp(iel) = 0.d0
! PUISSANCE EN EXPLICITE
  crvexp(iel) = volume(iel)*pwatt/volf
enddo

!--------
! FORMATS
!--------

 1000 format(' TERMES SOURCES UTILISATEURS POUR LA VARIABLE ',A8,/)

!----
! FIN
!----

return

end subroutine
