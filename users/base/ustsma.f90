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

subroutine ustsma &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ncesmp , iphas  , iappel ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr , icepdc , icetsm , itypsm , &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

!                    TERMES SOURCES DE MASSE
!                     POUR LA PHASE IPHAS

! IAPPEL = 1 :
!             CALCUL DU NOMBRE DE CELLULES AVEC SOURCE DE MASSE
!              NCESMP
! IAPPEL = 2 :
!             REPERAGE DES CELLULES AVEC SOURCE DE MASSE
!              ICETSM(NCESMP)
! IAPPEL = 3 :
!             CALCUL DES VALEURS DES COEFS DE SOURCE DE MASSE

!             L'equation de conservation de la masse devient :

!             d(rho)/dt + div(rho u) = GAMMA
!                                 -

!             L'equation d'une variable f devient :

!             d(f)/dt = ..... + GAMMA*(f_i - f)

!             discretisee en :

!             RHO*(f^(n+1) - f^(n))/DT = .....
!                                        + GAMMA*(f_i - f^(n+1))


!             Deux possibilites pour chaque variable f :
!                - Flux de masse a la valeur de f ambiante
!                           --> f_i = f^(n+1)
!                   (l'equation de f n'est alors pas modifiee)
!                - Flux de masse avec une valeur donnee pour f
!                           --> f_i specifie par l'utilisateur




!            VARIABLES A REMPLIR PAR L'UTILISATEUR :
!            =======================================

!             NCESMP : Nombre de cellules a source de masse

!             ICETSM(IELTSM) : Numero de la IELTSMieme cellule a
!                                source de masse (IELTSM<=NCESMP)

!             SMACEL(IEL,IPR(IPHAS)) : Valeur du flux de masse
!                                 GAMMA (en kg/m^3/s)
!                                 dans la IELieme cellule a source
!                                 de masse

!             ITYPSM(IEL,IVAR) : type de flux associe a la variable
!                                  IVAR dans la IELeme cellule a
!                                  source de masse (pour toutes les
!                                  variables sauf IVAR=IPR(IPHAS))
!                 + ITYPSM = 0 --> injection a la valeur ambiante
!                                  de IVAR
!                 + ITYPSM = 1 --> injection a une valeur donnee
!                                  de IVAR

!             SMACEL(IEL,IVAR) : Valeur de f_i pour la variable
!                                IVAR (pour toutes les variables
!                                sauf IVAR=IPR(IPHAS))


!            REMARQUES
!            =========
!            * Si ITYPSM(IEL,IVAR)=0, SMACEL(IEL,IVAR)
!               n'est pas utilise
!            * Si SMACEL(IEL,IPR(IPHAS))<0,
!               on enleve de la masse au systeme, Code_Saturne
!               utilise donc automatiquement f_i=f^(n+1)
!               QUELLES QUE SOIENT LES VALEURS DE ITYPSM(IEL)
!               et SMACEL(IEL,IVAR)

!            * Si une variable n'est pas reliee a la phase iphas
!               pour laquelle les informations precedentes ont ete
!               completees, aucun terme source ne sera impose.

!            * Pour un scalaire qui ne repondrait pas a l'equation
!               d(rho f)/dt + d(rho U f)/dx = ...
!               (champ convecteur different par exemple)
!               il est incorrect d'imposer le terme source de masse
!               comme fait ici (sauf en cas d'injection a la valeur
!               ambiante). Il faut introduire directement le terme
!               source comme un terme source scalaire dans ustssc.


! IDENTIFICATION DES CELLULES/FACES DE BORD/FACES INTERNES
! ========================================================

!  Les commandes GETCEL, GETFBR et GETFAC permettent d'identifier
!  respectivement les cellules, faces ou faces de bord en fonction
!  de differents criteres.

!  GETCEL(CHAINE,NLELT,LSTELT) :
!  - CHAINE est une chaine de caractere fournie par l'utilisateur
!    qui donne les criteres de selection
!  - NLTELT est renvoye par la commande. C'est un entier qui
!    correspond au nombre de cellules trouveees repondant au
!    critere
!  - LSTELT est renvoye par la commande. C'est un tableau d'entiers
!    de taille NLTELT donnant la liste des cellules trouvees
!    repondant au critere.

!  CHAINE peut etre constitue de :
!  - references de couleurs (ex. : 1, 8, 26, ...
!  - references de groupes (ex. : entrees, groupe1, ...)
!  - criteres geometriques (ex. X<0.1, Y>=0.25, ...)
!  Ces criteres peuvent etre combines par des operateurs logiques
!  (AND et OR) et des parentheses
!  ex. : '1 AND (groupe2 OR groupe3) AND Y<1' permettra de recuperer
!  les cellules de couleur 1, appartenant aux groupes 'groupe2'
!  ou 'groupe3' et de coordonnee Y inferieure a 1.

!  La syntaxe des commandes GETFBR et GETFAC est identique.

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
! ncesmp           ! e  ! <-- ! nombre de cellules avec tsm                    !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! iphas            ! e  ! <-- ! numero de phase                                !
! iappel           ! e  ! <-- ! indique les donnes a renvoyer                  !
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
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant            prec)          !
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
! (ncesmp,nvar)    !    !     !  source de masse                               !
!                  !    !     ! pour ivar=ipr, smacel=flux de masse            !
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
include "pointe.h"
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
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          nideve , nrdeve , nituse , nrtuse
integer          iphas  , iappel

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
double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6)
double precision smacel(ncesmp,nvar)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ieltsm
integer          ifac, iutile, ii
integer          ilelt, nlelt

double precision vent, vent2
double precision dh, ustar2
double precision xkent, xeent
double precision flucel
double precision vtot  , gamma

!===============================================================================

idebia = idbia0
idebra = idbra0

if(iappel.eq.1.or.iappel.eq.2) then

!===============================================================================
! 1. POUR CHAQUE PHASE : UN OU DEUX APPELS

!    PREMIER APPEL :

!        IAPPEL = 1 : NCESMP : CALCUL DU NOMBRE DE CELLULES
!                                AVEC TERME SOURCE DE MASSE


!    DEUXIEME APPEL (POUR LES PHASES AVEC NCESMP > 0) :

!        IAPPEL = 2 : ICETSM : REPERAGE DU NUMERO DES CELLULES
!                                AVEC TERME SOURCE DE MASSE


! REMARQUES :

!        Ne pas utiliser SMACEL dans cette section
!          (il est rempli au troisieme appel, IAPPEL = 3)

!        Ne pas utiliser ICETSM dans cette section
!           au premier appel (IAPPEL = 1)

!        On passe ici a chaque pas de temps
!           (ATTENTION au cout calcul de vos developpements)

!===============================================================================


!  1.1 A completer par l'utilisateur : selection des cellules
!  -----------------------------------------------------------

! Exemple 1 : Aucune source de masse (defaut)
  ieltsm = 0


! Exemple 2 : Sources de masse pour la phase 1
!               dans les cellules ayant une face de couleur 6
!            et dans les cellules dont le centre est a x < 2

!     Dans ce test en deux parties, il faut faire attention a
!       ne pas compter les cellules deux fois (une cellule dont une
!       face est couleur 6 peut aussi etre a x < 2, mais il ne
!       faut la compter qu'une seule fois). Il faut egalement noter
!       que lors du premier passage, le tableau ICETSM n'existe pas
!       encore (ne pas l'utiliser en dehors des tests IAPPEL.EQ.2).

  iutile = 0
  if(iutile.eq.1) then

    if(iphas.eq.1) then

      ieltsm = 0

!     Cellules dont le centre est a 250 < x < 500

      CALL GETCEL('X < 500.0 and X > 250.0',NLELT,LSTELT)
      do ilelt = 1, nlelt
        ii = lstelt(ilelt)
        ieltsm = ieltsm + 1
        if (iappel.eq.2) icetsm(ieltsm) = ii
      enddo


!     Cellules de bord dont une face est de couleur 3

      CALL GETFBR('3',NLELT,LSTELT)
      do ilelt = 1, nlelt
        ifac = lstelt(ilelt)
        ii   = ifabor(ifac)
!       On ne comptabilise cette cellule que si on ne l'a pas deja vue
!         pour cela on prend la negation du test precedent
        if (.not.(xyzcen(1,ii).lt.500.d0.and.                     &
                  xyzcen(1,ii).gt.250.d0)    )then
          ieltsm = ieltsm + 1
          if (iappel.eq.2) icetsm(ieltsm) = ii
        endif
      enddo

    else
      ieltsm = 0
    endif

  endif


!  1.2 Sous section generique a ne pas modifier
!  ---------------------------------------------

! --- Pour IAPPEL = 1,
!      Renseigner NCESMP, nombre de cellules avec terme source de masse
!      Le bloc ci dessous est valable pourles 2 exemples ci dessus

  if (iappel.eq.1) then
    ncesmp = ieltsm
  endif

!-------------------------------------------------------------------------------

elseif(iappel.eq.3) then

!===============================================================================

! 2. POUR CHAQUE PHASE AVEC NCESMP > 0 , TROISIEME APPEL

!      TROISIEME APPEL (POUR LES PHASES AVEC NCESMP > 0) :

!       IAPPEL = 3 : ITYPSM : TYPE DE SOURCE DE MASSE
!                    SMACEL : SOURCE DE MASSE


!   REMARQUE :

!      ATTENTION, Si on positionne ITYPSM(IEL,IVAR) A 1, IL FAUT
!      egalement renseigner SMACEL(IEL,IVAR)

!===============================================================================



!  2.1 A completer par l'utilisateur ITYPSM et SMACEL
!  -----------------------------------------------------

! Exemple 1 : simulation d'une entree par des termes de source de masse
!           et affichage du flux de masse total pour la phase 1

  if(iphas.eq.1) then

    vent = 0.1d0
    vent2 = vent**2
    dh     = 0.5d0
!         Calcul de la vitesse de frottement au carre (USTAR2)
!           et de k et epsilon en entree (XKENT et XEENT) a partir
!           de lois standards en conduite circulaire
!           (leur initialisation est inutile mais plus propre)
    ustar2 = 0.d0
    xkent  = epzero
    xeent  = epzero

    call keendb                                                   &
    !==========
      ( vent2, dh, ro0(iphas), viscl0(iphas), cmu, xkappa,        &
        ustar2, xkent, xeent )

    flucel = 0.d0
    do ieltsm = 1, ncesmp
      smacel(ieltsm,ipr(iphas)) = 30000.d0
      itypsm(ieltsm,iv(iphas)) = 1
      smacel(ieltsm,iv(iphas)) = vent
      if (itytur(iphas).eq.2) then
        itypsm(ieltsm,ik(iphas)) = 1
        smacel(ieltsm,ik(iphas)) = xkent
        itypsm(ieltsm,iep(iphas)) = 1
        smacel(ieltsm,iep(iphas)) = xeent
      else if (itytur(iphas).eq.3) then
        itypsm(ieltsm,ir11(iphas)) = 1
        itypsm(ieltsm,ir12(iphas)) = 1
        itypsm(ieltsm,ir13(iphas)) = 1
        itypsm(ieltsm,ir22(iphas)) = 1
        itypsm(ieltsm,ir23(iphas)) = 1
        itypsm(ieltsm,ir33(iphas)) = 1
        smacel(ieltsm,ir11(iphas)) = 2.d0/3.d0*xkent
        smacel(ieltsm,ir12(iphas)) = 0.d0
        smacel(ieltsm,ir13(iphas)) = 0.d0
        smacel(ieltsm,ir22(iphas)) = 2.d0/3.d0*xkent
        smacel(ieltsm,ir23(iphas)) = 0.d0
        smacel(ieltsm,ir33(iphas)) = 2.d0/3.d0*xkent
        itypsm(ieltsm,iep(iphas)) = 1
        smacel(ieltsm,iep(iphas)) = xeent
      else if (iturb(iphas).eq.50) then
        itypsm(ieltsm,ik(iphas)) = 1
        smacel(ieltsm,ik(iphas)) = xkent
        itypsm(ieltsm,iep(iphas)) = 1
        smacel(ieltsm,iep(iphas)) = xeent
        itypsm(ieltsm,iphi(iphas)) = 1
        smacel(ieltsm,iphi(iphas)) = 2.d0/3.d0
!     Il n'y a pas de terme source de masse dans l'equation de f_barre
      else if (iturb(iphas).eq.60) then
        itypsm(ieltsm,ik(iphas)) = 1
        smacel(ieltsm,ik(iphas)) = xkent
        itypsm(ieltsm,iomg(iphas))= 1
        smacel(ieltsm,iomg(iphas))= xeent/cmu/xkent
      endif
      if(nscal.gt.0) then
        do ii = 1, nscal
          if(iphsca(ii).eq.iphas) then
            itypsm(ieltsm,isca(ii)) = 1
            smacel(ieltsm,isca(ii)) = 1.d0
          endif
        enddo
      endif
      flucel = flucel+                                            &
                volume(icetsm(ieltsm))*smacel(ieltsm,ipr(iphas))
    enddo

    if (irangp.ge.0) then
      call parsom (flucel)
    endif

    if (iwarni(ipr(iphas)).ge.1) then
      write(nfecra,1000) iphas, flucel
    endif

  endif

!-------------------------------------------------------------------------------

! Exemple 2 : simulation d'un soutirage (par une pompe par exemple)
!               total de 80 000 kg/s
!             On suppose que l'on souhaite repartir uniformement ce
!               puits de masse sur les NCESMP cellules selectionnees
!               plus haut.


!     Le test sur IUTILE permet de ne pas passer dans l'exemple si
!       on oublie de l'eliminer, lors de la mise en place d'un calcul
!       reel.

  iutile = 0
  if(iutile.eq.1) then

    if(iphas.eq.1) then

!     Calcul du volume total de la zone ou est impose le terme source
!       (le cas des calculs paralleles est prevu avec parsom)
      vtot = 0.d0
      do ieltsm = 1, ncesmp
        vtot = vtot + volume(icetsm(ieltsm))
      enddo
      if (irangp.ge.0) then
        call parsom (vtot)
      endif

!     Le puits de masse est GAMMA = -80000/VTOT en kg/(m3 s)
!       (quel que soit le nombre de cellules NCESMP)
!     On l'impose ci-dessous (avec un test au cas ou VTOT=0)
!     On calcule au passage le terme puits total pour verification.

      if (vtot.gt.0.d0) then
        gamma = -80000.d0/vtot
      else
        write(nfecra,9000) iphas, vtot
        call csexit (1)
      endif

      flucel = 0.d0
      do ieltsm = 1, ncesmp
        smacel(ieltsm,ipr(iphas)) = gamma
        flucel = flucel+                                          &
                volume(icetsm(ieltsm))*smacel(ieltsm,ipr(iphas))
      enddo

      if (irangp.ge.0) then
        call parsom (flucel)
      endif

      if (iwarni(ipr(iphas)).ge.1) then
        write(nfecra,2000) iphas, flucel, vtot
      endif

    endif

  endif

!-------------------------------------------------------------------------------

endif

!--------
! FORMATS
!--------

 1000 format(/,'PHASE ',I3,                                             &
         ' : FLUX DE MASSE GENERE DANS LE DOMAINE : ',E14.5,/)

 2000 format(/,'PHASE ',I3,                                             &
         ' : FLUX DE MASSE GENERE DANS LE DOMAINE : ',E14.5,/,    &
'                           REPARTI SUR UN VOLUME : ',E14.5)

 9000 format(/,'PHASE ',I3,                                             &
         ' : ERREUR DANS USTSMA ',/,                        &
'   LE VOLUME DE LA ZONE AVEC PUITS DE MASSE VAUT = ',E14.5,/)

return

end
