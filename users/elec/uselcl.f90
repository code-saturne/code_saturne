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

subroutine uselcl &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , itrifb , itypfb , izfppp ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!    ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE
!                MODULE ELECTRIQUE
!   (Effet Joule, Arc Electrique, Conduction ionique)
!    REMPLISSAGE DU TABLEAU DE CONDITIONS AUX LIMITES
!    (ICODCL,RCODCL) POUR LES VARIABLES INCONNUES
!    PENDANT DE USCLIM.F



!    CE SOUS PROGRAMME UTILISATEUR EST OBLIGATOIRE
!    =============================================


! Introduction
! ============

! Here one defines boundary conditions on a per-face basis.

! Boundary faces may be identified using the 'getfbr' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


! Boundary condition types
! ========================

! Boundary conditions setup for standard variables (pressure, velocity,
! turbulence, scalars) is described precisely in the 'usclim' subroutine.

! Detailed explanation will be found in the theory guide.


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
!   (lndfac)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !  (optionnel)                                   !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! itrifb(nfabor    ! te ! <-- ! indirection pour tri des faces de brd          !
!  nphas      )    !    !     !                                                !
! itypfb(nfabor    ! te ! --> ! type des faces de bord                         !
!  nphas      )    !    !     !                                                !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
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
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
! w1,2,3,4,5,6     ! tr ! --- ! tableaux de travail                            !
!  (ncelet         !    !     !  (calcul du gradient de pression)              !
! coefu            ! tr ! --- ! tab de trav                                    !
!  (nfabor,3)      !    !     !  (calcul du gradient de pression)              !
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
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "elincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor,nphas), itypfb(nfabor,nphas)
integer          izfppp(nfabor)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision coefu(nfabor,ndim)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ifac, ii, iphas , iel
integer          idim
integer          izone,iesp
integer          ilelt, nlelt

double precision uref2, d2s3
double precision rhomoy, dhy, ustar2
double precision xkent, xeent
double precision z1   , z2

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!===============================================================================

if(1.eq.1) then
  write(nfecra,9001)
  call csexit (1)
endif

 9001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES COND. LIM.      ',/,&
'@    =========                                               ',/,&
'@                      MODULE ELECTRIQUE                     ',/,&
'@                                                            ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR uselcl DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@     Ce sous-programme utilisateur permet de definir les    ',/,&
'@       conditions aux limites. Il est indispensable.        ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  INITIALISATIONS

!===============================================================================

idebia = idbia0
idebra = idbra0

d2s3 = 2.d0/3.d0

!===============================================================================
! 2.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR LES FACES DE BORD
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!           ON IMPOSE LA CONDITION LIMITE

!          IMPOSER ICI LES CONDITIONS LIMITES SUR LES FACES DE BORD

!===============================================================================

iphas = 1


! --- On impose en couleur 1 une entree ; exemple de Cathode
!     ======================================================

CALL GETFBR('1',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  itypfb(ifac,iphas) = ientre

!      - Numero de zone (on numerote de 1 a n)
  izone = 1

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

  rcodcl(ifac,iu(iphas),1) = 0.d0
  rcodcl(ifac,iv(iphas),1) = 0.d0
  rcodcl(ifac,iw(iphas),1) = 0.d0

!         Turbulence

!     (ITYTUR est un indicateur qui vaut ITURB/10)
  if (itytur(iphas).eq.2 .or. itytur(iphas).eq.3                  &
       .or. iturb(iphas).eq.50 .or. iturb(iphas).eq.60) then

    uref2 = rcodcl(ifac,iu(iphas),1)**2                           &
           +rcodcl(ifac,iv(iphas),1)**2                           &
           +rcodcl(ifac,iw(iphas),1)**2
    uref2 = max(uref2,1.d-12)


!       Exemple de turbulence calculee a partir
!         de formules valables pour une conduite

!       On veillera a specifier le diametre hydraulique
!         adapte a l'entree courante.

!       On s'attachera egalement a utiliser si besoin une formule
!         plus precise pour la viscosite dynamique utilisee dans le
!         calcul du nombre de Reynolds (en particulier, lorsqu'elle
!         est variable, il peut etre utile de reprendre ici la loi
!         imposee dans USPHYV. On utilise ici par defaut la valeur
!         VISCL0 donnee dans USINI1
!       En ce qui concerne la masse volumique, on dispose directement
!         de sa valeur aux faces de bord (ROMB) et c'est celle que
!         utilise donc ici (elle est en particulier coherente avec
!         le traitement implante dans USPHYV, en cas de masse
!         volumique variable)

!         Diametre hydraulique
    dhy     = 0.075d0

!         Calcul de la vitesse de frottement au carre (USTAR2)
!           et de k et epsilon en entree (XKENT et XEENT) a partir
!           de lois standards en conduite circulaire
!           (leur initialisation est inutile mais plus propre)
    rhomoy = propfb(ifac,ipprob(irom(iphas)))
    ustar2 = 0.d0
    xkent  = epzero
    xeent  = epzero

    call keendb                                                   &
    !==========
     ( uref2, dhy, rhomoy, viscl0(iphas), cmu, xkappa,            &
       ustar2, xkent, xeent )

    if (itytur(iphas).eq.2) then

      rcodcl(ifac,ik(iphas),1)  = xkent
      rcodcl(ifac,iep(iphas),1) = xeent

    elseif(itytur(iphas).eq.3) then

      rcodcl(ifac,ir11(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir22(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir33(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir12(iphas),1) = 0.d0
      rcodcl(ifac,ir13(iphas),1) = 0.d0
      rcodcl(ifac,ir23(iphas),1) = 0.d0
      rcodcl(ifac,iep(iphas),1)  = xeent

    elseif (iturb(iphas).eq.50) then

      rcodcl(ifac,ik(iphas),1)   = xkent
      rcodcl(ifac,iep(iphas),1)  = xeent
      rcodcl(ifac,iphi(iphas),1) = d2s3
      rcodcl(ifac,ifb(iphas),1)  = 0.d0

    elseif (iturb(iphas).eq.60) then

      rcodcl(ifac,ik(iphas),1)   = xkent
      rcodcl(ifac,iomg(iphas),1) = xeent/cmu/xkent

    endif

  endif

! --- On traite les scalaires

!      Enthalpie en J/kg

  ii = ihm
  icodcl(ifac,isca(ii))   = 1
  rcodcl(ifac,isca(ii),1) = 1.d6

!  Potentiel electrique reel impose a 0. volts (exemple de Cathode en arc)

  ii = ipotr
  icodcl(ifac,isca(ii))   = 1
  rcodcl(ifac,isca(ii),1) = 0.d0

!  Fraction massique des (N-1) constituants

  if ( ngazg .gt. 1 ) then
    do iesp=1,ngazg-1
      ii = iycoel(iesp)
      icodcl(ifac,isca(ii))   = 1
      rcodcl(ifac,isca(ii),1) = 0.d0
    enddo
  endif

!  Specifique Version Effet Joule :

!       Potentiel Imaginaire impose a 0

  if ( ippmod(ieljou).ge. 2 ) then
    ii = ipoti
    icodcl(ifac,isca(ii))   = 1
    rcodcl(ifac,isca(ii),1) = 0.d0
  endif

!  Specifique Version Arc Electrique :

!       Potentiel vecteur : Flux nul

  if ( ippmod(ielarc).ge.2 ) then
    do idim= 1,ndimve
      ii = ipotva(idim)
      icodcl(ifac,isca(ii))   = 3
      rcodcl(ifac,isca(ii),3) = 0.d0
    enddo
  endif

enddo

! --- On impose en couleur 5 une entree/sortie ;
!     ====================================== exemple d'Electrode en Joule
!                                            ============================

CALL GETFBR('5',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          SORTIE : FLUX NUL VITESSE ET TEMPERATURE, PRESSION IMPOSEE
!            Noter que la pression sera recalee a P0
!                sur la premiere face de sortie libre (ISOLIB)

  itypfb(ifac,iphas)   = isolib

!      - Numero de zone (on numerote de 1 a n)
  izone = 2

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

! --- On traite les scalaires rattaches a la phase courante

!  Enthalpie en J/kg  (Par defaut flux nul avec ISOLIB)
!     Rien a faire

!  Fraction massique des (N-1) constituants (Par defaut flux nul avec ISOLIB)
!     Rien a faire

!  Specifique Version Effet Joule :

!     En effet Joule,
!       si l'on souhaite faire un calcul en recalant les conditions
!         aux limites (utiliser IELCOR=1 dans useli1)
!       pour atteindre la valeur de la puissance PUISIM
!         (a imposer dans useli1 en Ampere.Volt)
!       on multiplie la condition limite initiale sur le potentiel
!          reel (et sur le potentiel imaginaire s'il est pris en
!          compte) par le coefficient COEJOU.
!       COEJOU est determine automatiquement pour que la puissance
!          dissipee par effet Joule (partie reelle et partie
!          imaginaire si besoin) soit PUISIM
!       au debut du calcul, COEJOU vaut 1 ; COEJOU est transmis dans
!          les fichiers suites.

!     Si on ne souhaite pas faire un calcul avec recalage, on impose
!       directement une valeur adaptee.

  if ( ippmod(ieljou).ge. 1 ) then
    ii = ipotr
    icodcl(ifac,isca(ii))   = 1
    if(ielcor.eq.1) then
      rcodcl(ifac,isca(ii),1) = 500.d0*coejou
    else
      rcodcl(ifac,isca(ii),1) = 500.d0
    endif
  endif

  if ( ippmod(ieljou).ge. 2 ) then
    ii = ipoti
    icodcl(ifac,isca(ii))   = 1
    if(ielcor.eq.1) then
      rcodcl(ifac,isca(ii),1) = sqrt(3.d0)*500.d0*coejou
    else
      rcodcl(ifac,isca(ii),1) = sqrt(3.d0)*500.d0
    endif
  endif

enddo

! --- On impose en couleur 2 une entree/sortie ;
!     ============================== exemple d'Anode en arc electrique
!                                    =================================

CALL GETFBR('2',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          SORTIE : FLUX NUL VITESSE ET TEMPERATURE, PRESSION IMPOSEE
!            Noter que la pression sera recalee a P0
!                sur la premiere face de sortie libre (ISOLIB)

  itypfb(ifac,iphas)   = isolib

!      - Numero de zone (on numerote de 1 a n)
  izone = 3

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

! --- On traite les scalaires rattaches a la phase courante

!  Enthalpie en J/kg  (Par defaut flux nul avec ISOLIB)
!     Rien a faire

!  Potentiel electrique reel

!     En arc electrique,
!       si l'on souhaite faire un calcul en recalant le potentiel
!          de l'anode (utiliser IELCOR=1 dans useli1)
!       pour atteindre la valeur du courant COUIMP
!         (a imposer dans useli1 en Amperes)
!       on utilise alors la valeur DPOT comme condition limite
!       DPOT est en effet automatiquement adaptee par le calcul
!          pour que (j.E Volume/DPOT) = COUIMP
!          (initialiser DPOT dans useli1 en Volts avec une valeur
!           representative de la difference de potentiel imposee)

!     Si on ne souhaite pas faire un calcul avec recalage,  on impose
!       directement une valeur adaptee au cas
!       (par exemple, ici 1000 Volts ).

  ii = ipotr
  icodcl(ifac,isca(ii))   = 1

  if ( ippmod(ielarc).ge.1 .and. ielcor .eq.1) then
    rcodcl(ifac,isca(ii),1) = dpot
  else
    rcodcl(ifac,isca(ii),1) = 1000.d0
  endif


!  Fraction massique des (N-1) constituants (Par defaut flux nul avec ISOLIB)

!  Specifique Version Arc Electrique :
!      Potentiel vecteur : flux nul (par defaut)


enddo

! --- On impose en couleur 3 une paroi
!     ================================

CALL GETFBR('3',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          PAROI : DEBIT NUL (FLUX NUL POUR LA PRESSION)
!                  FROTTEMENT POUR LES VITESSES (+GRANDEURS TURB)
!                  FLUX NUL SUR LES SCALAIRES

!          Pour un calcul arc electrique 3D, on cale le potentiel
!            vecteur avec une condition de Dirichlet issue des valeurs
!            du potentiel vecteur au pas de temps precedent
!            dans une zone de paroi choisie
!            Par defaut, ailleurs, un flux nul s'applique (paroi isolee).

  itypfb(ifac,iphas)   = iparoi

!      - Numero de zone (on numerote de 1 a n)
  izone = 4

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

  if ( ippmod(ielarc).ge.2 ) then
    if ( cdgfbo(1,ifac) .le.  2.249d-2  .or.                      &
         cdgfbo(1,ifac) .ge.  2.249d-2  .or.                      &
         cdgfbo(3,ifac) .le. -2.249d-2  .or.                      &
         cdgfbo(3,ifac) .ge.  2.249d-2       ) then
      iel = ifabor(ifac)
      do idim = 1, ndimve
        ii = ipotva(idim)
        icodcl(ifac,isca(ii))   = 1
        rcodcl(ifac,isca(ii),1) = rtpa(iel,isca(ii))
      enddo
    endif
  endif

enddo

! --- On impose en couleur 51 : anode avec claquage
!     =============================================

CALL GETFBR('51',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  itypfb(ifac,iphas)   = iparoi

!      - Numero de zone (on numerote de 1 a n)
  izone = 5

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

! ---- Enthalpie (J/kg ) : coef echange impose

  ii=ihm
  icodcl(ifac,isca(ii))   = 1
  rcodcl(ifac,isca(ii),1) = 2.d4
  rcodcl(ifac,isca(ii),2) = 1.d5

!  Potentiel electrique reel

  ii = ipotr
  icodcl(ifac,isca(ii))   = 1

  if ( ippmod(ielarc).ge.1 .and. ielcor .eq.1) then
    rcodcl(ifac,isca(ii),1) = dpot
  else
    rcodcl(ifac,isca(ii),1) = 100.d0
  endif

!       Si CLAQUAGE : a adapter en fonction du cas et du
!                     sous-programme USELRC

  if ( ippmod(ielarc).ge.1 .and. ielcor .eq.1) then
    if(iclaq.eq.1 .and. ntcabs.le.ntdcla+30) then

      z1 = zclaq - 2.d-4
      if(z1.le.0.d0) z1 = 0.d0
      z2 = zclaq + 2.d-4
      if(z2.ge.2.d-2) z2 = 2.d-2

      if( cdgfbo(3,ifac).ge.z1 .and.                              &
           cdgfbo(3,ifac).le.z2       ) then
        icodcl(ifac,isca(ii))   = 1
        rcodcl(ifac,isca(ii),1) = dpot
      else
        icodcl(ifac,isca(ii))   = 3
        rcodcl(ifac,isca(ii),3) = 0.d0
      endif
    endif
  endif

!       Potentiel vecteur : Flux nul

  if ( ippmod(ielarc).ge.2 ) then
    do idim= 1,ndimve
      ii = ipotva(idim)
      icodcl(ifac,isca(ii))   = 3
      rcodcl(ifac,isca(ii),3) = 0.d0
    enddo
  endif

enddo

! --- On impose en couleur 4 une symetrie
!     ===================================

CALL GETFBR('4',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          SYMETRIES

  itypfb(ifac,iphas)   = isymet

!      - Numero de zone (on numerote de 1 a n)
  izone = 6

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

!     Par defaut tous les scalaires (potentiels en particulier)
!       recoivent une condition de flux nul.
!     En effet Joule, on peut souhaiter imposer une condition
!       d'antisymetrie sur le potentiel imaginaire selon la
!       configuration des electrodes :
  if ( ippmod(ieljou).ge. 2 ) then
    ii = ipoti
    icodcl(ifac,isca(ii))   = 1
    rcodcl(ifac,isca(ii),1) = 0.d0
  endif

enddo

!----
! FORMATS
!----

!----
! FIN
!----

return
end
