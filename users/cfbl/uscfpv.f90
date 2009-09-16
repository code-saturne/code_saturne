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

subroutine uscfpv &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , nphmx  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! ROUTINE UTILISATEUR : REMPLISSAGE DES VARIABLES PHYSIQUES
!    POUR LA PHYSIQUE PARTICULIERE : COMPRESSIBLE SANS CHOC
!    PENDANT DE USPHYV.F



! ATTENTION :
! =========


! Il est INTERDIT de modifier la viscosite turbulente VISCT ici
!        ========
!  (une routine specifique est dediee a cela : usvist)


!  Il FAUT AVOIR PRECISE ICP(IPHAS) = 1
!     ==================
!     automatiquement fait dans uscfth en fonction de la thermo
!     choisie


!  Il FAUT AVOIR PRECISE IVISLS(Numero de scalaire) = 1
!     ==================
!     dans usini1 si on souhaite une diffusivite VISCLS variable
!     pour le scalaire considere (sinon: ecrasement memoire).




! Remarques :
! ---------

! Cette routine est appelee au debut de chaque pas de temps

!    Ainsi, AU PREMIER PAS DE TEMPS (calcul non suite), les seules
!    grandeurs initialisees avant appel sont celles donnees
!      - dans usini1 :
!             . la masse volumique (initialisee a RO0(IPHAS))
!             . la viscosite       (initialisee a VISCL0(IPHAS))
!      - dans usiniv :
!             . les variables de calcul  (initialisees a 0 par defaut
!             ou a la valeur donnee dans usiniv)

! On peut donner ici les lois de variation aux cellules
!     - de la viscosite moleculaire                VISCL  kg/(m s)
!     - de la chaleur specifique associee          CP     J/(kg degres)
!     - de la conductivite thermique associee      LAMBDA W/(m degres)
!     - des "diffusivites" associees aux scalaires VISCLS kg/(m s)

! La masse volumique ne doit pas être renseignée : en compressible,
!   c'est une variable résolue, que l'on peut initialiser si
!   nécessaire dans USCFXI (RTP).

! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)


! Il est conseille de ne garder dans ce sous programme que
!    le strict necessaire.


! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


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
! nphmx            ! e  ! <-- ! nphsmx                                         !
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
! w1...3(ncelet    ! tr ! --- ! tableau de travail                             !
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
include "optcal.h"
include "cstphy.h"
include "entsor.h"
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
integer          nideve , nrdeve , nituse , nrtuse , nphmx

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml), maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ivart, iel, iphas
integer          ipcvis, ipcvsv, ipccp
integer          ipcvsl, ith, iscal, ii, iccfth, imodif
double precision varam, varbm, varcm, vardm
double precision varal, varbl, varcl, vardl
double precision varac, varbc
double precision xrtp

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) then
  iuscfp = 0
  return
endif

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 0. INITIALISATIONS A CONSERVER
!===============================================================================

! --- Initialisation memoire

idebia = idbia0
idebra = idbra0

!===============================================================================


!   LES EXEMPLES FANTAISISTES SUIVANTS SONT A ADAPTER PAR L'UTILISATEUR
!   ====================================================================

!                ILS SONT TOUS SUIVIS DE "CALL CSEXIT (1)" PAR SECURITE
!                                (stop a enlever a l'utilisation)


!   Il est conseille de ne garder dans ce sous programme que
!     le strict necessaire.


! EXEMPLE 1 : VISCOSITE VARIABLE EN FONCTION DE LA TEMPERATURE

! EXEMPLE 2 : VISCOSITE EN VOLUME VARIABLE EN FONCTION DE LA TEMPERATURE

! EXEMPLE 3 : CHALEUR SPECIFIQUE VARIABLE EN FONCTION DE LA TEMPERATURE

! EXEMPLE 4 : CONDUCTIVITE THERMIQUE VARIABLE
!                  EN FONCTION DE LA TEMPERATURE

! EXEMPLE 5 : DIFFUSIVITE VARIABLE EN FONCTION DE LA TEMPERATURE
!                  POUR LES SCALAIRES
!===============================================================================








!===============================================================================
!  EXEMPLE 1 : VISCOSITE VARIABLE EN FONCTION DE LA TEMPERATURE
! ===========
!    Ci dessous on donne pour toutes les phases la meme loi pour
!       la viscosite
!    Les valeurs de cette propriete doivent etre fournies au centre des
!       cellules.
!  ===================================================================

! --- Boucle sur les phases : debut
do iphas = 1, nphas


!   Positions des variables, coefficients
!   -------------------------------------

! --- Numero de variable température pour la phase courante iphas
!       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
!          IVART = ISCA(2)

  ivart = isca(itempk(iphas))

! --- Rang de la viscosite dynamique moleculaire de la phase IPHAS
!     dans PROPCE, prop. physiques au centre des elements       : IPCVIS

  ipcvis = ipproc(iviscl(iphas))

! --- Coefficients des lois choisis et imposes par l'utilisateur
!       Les valeurs donnees ici sont fictives

  varam = -3.4016d-9
  varbm =  6.2332d-7
  varcm = -4.5577d-5
  vardm =  1.6935d-3

!   Viscosite moleculaire dynamique en kg/(m s) au centre des cellules
!   ------------------------------------------------------------------
!       loi              MU        =
!                              T  *( T  *( AM  * T +  BM  )+ CM  )+ DM
!       soit    PROPCE(IEL,IPCVIS) =
!     &                       XRTP*(XRTP*(VARAM*XRTP+VARBM)+VARCM)+VARDM

  do iel = 1, ncel
    xrtp = rtp(iel,ivart)
    propce(iel,ipcvis) =                                          &
                       xrtp*(xrtp*(varam*xrtp+varbm)+varcm)+vardm
  enddo


enddo
! --- Boucle sur les phases : fin



! --- A enlever a l'utilisation
if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif





!===============================================================================
!  EXEMPLE 2: VISCOSITE EN VOLUME VARIABLE EN FONCTION DE LA TEMPERATURE
! ==========
!    Ci dessous on donne pour toutes les phases la meme loi pour
!       la viscosite
!    Les valeurs de cette propriete doivent etre fournies au centre des
!       cellules.
!  ===================================================================

! --- Boucle sur les phases : debut
do iphas = 1, nphas


!   Positions des variables, coefficients
!   -------------------------------------

! --- Numero de variable température pour la phase courante iphas
!       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
!          IVART = ISCA(2)

  ivart = isca(itempk(iphas))

! --- Rang de la viscosite dynamique moleculaire en volume de la phase IPHAS
!     dans PROPCE, prop. physiques au centre des elements       : IPCVSV

  if(iviscv(iphas).gt.0) then
    ipcvsv = ipproc(iviscv(iphas))
  else
    ipcvsv = 0
  endif

! --- Stop si non variable

  if(ipcvsv.le.0) then
    write(nfecra,2000) iphas, iphas, iviscv(iphas)
    call csexit (1)
  endif


! --- Coefficients des lois choisis et imposes par l'utilisateur
!       Les valeurs donnees ici sont fictives

  varam = -3.4016d-9
  varbm =  6.2332d-7
  varcm = -4.5577d-5
  vardm =  1.6935d-3

!   Viscosite en volume en kg/(m s) au centre des cellules
!   ------------------------------------------------------
!       loi            KAPPA       =
!                              T  *( T  *( AM  * T +  BM  )+ CM  )+ DM
!       soit    PROPCE(IEL,IPCVSV) =
!     &                       XRTP*(XRTP*(VARAM*XRTP+VARBM)+VARCM)+VARDM

  do iel = 1, ncel
    xrtp = rtp(iel,ivart)
    propce(iel,ipcvsv) =                                          &
                       xrtp*(xrtp*(varam*xrtp+varbm)+varcm)+vardm
  enddo

enddo
! --- Boucle sur les phases : fin



! --- A enlever a l'utilisation
if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif





!===============================================================================
!  EXEMPLE 3 : CHALEUR SPECIFIQUE VARIABLE EN FONCTION DE LA TEMPERATURE
! ===========

!    Ci dessous on donne pour toutes les phases la meme loi pour
!       la chaleur specifique
!    Les valeurs de cette propriete doivent etre fournies au centre des
!       cellules.

!===============================================================================

!                             ATTENTION           !

!     NE PAS RETIRER LA MISE A JOUR DE Cv A LA FIN DE L'EXEMPLE

!===============================================================================

!===============================================================================

! --- Boucle sur les phases : debut
do iphas = 1, nphas


!   Positions des variables, coefficients
!   -------------------------------------

! --- Numero de variable température pour la phase courante iphas
!       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
!          IVART = ISCA(2)

  ivart = isca(itempk(iphas))

! --- Rang de la chaleur specifique de la phase courante IPHAS
!     dans PROPCE, prop. physiques au centre des elements       : IPCCP

  if(icp(iphas).gt.0) then
    ipccp  = ipproc(icp   (iphas))
  else
    ipccp  = 0
  endif

! --- Stop si CP ou CV n'est pas variable

  if(ipccp.le.0) then
    write(nfecra,1000) iphas, iphas, icp(iphas)
    call csexit (1)
  endif
  if(icv(iphas).le.0) then
    write(nfecra,1001) iphas, iphas, icv(iphas)
    call csexit (1)
  endif


! --- Coefficients des lois choisis et imposes par l'utilisateur
!       Les valeurs donnees ici sont fictives

  varac = 0.00001d0
  varbc = 1000.0d0



!   Chaleur specifique J/(kg degres) au centre des cellules
!   --------------------------------------------------------
!       loi              CP        =  AC  * T   +  BM
!       soit    PROPCE(IEL,IPCCP ) = VARAC*XRTP + VARBC

  do iel = 1, ncel
    xrtp = rtp(iel,ivart)
    propce(iel,ipccp ) = varac*xrtp + varbc
  enddo


! --- Mise a jour de Cv

  iccfth = 432
  imodif = 0

  call uscfth                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   iccfth , imodif , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   propce(1, ipproc(icv(iphas)) )    , w1     , w2     , w3     , &
!        ---------------------------------
   rdevel , rtuser , ra     )

enddo
! --- Boucle sur les phases : fin



! --- A enlever a l'utilisation
if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif




!===============================================================================
!  EXEMPLE 4 : CONDUCTIVITE THERMIQUE VARIABLE EN FONCTION
! ===========    DE LA TEMPERATURE

!    Ci dessous on donne pour toutes les phases la meme loi pour lambda
!    Les valeurs de cette propriete doivent etre fournies au centre des
!       cellules.
!  ===================================================================

! --- Boucle sur les phases : debut
do iphas = 1, nphas


!   Positions des variables, coefficients
!   -------------------------------------

! --- Numero de variable température pour la phase courante iphas
!       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
!          IVART = ISCA(2)

  ivart = isca(itempk(iphas))

! --- Rang de Lambda de la temperature de phase courante IPHAS
!     dans PROPCE, prop. physiques au centre des elements       : IPCVSL

  if(ivisls(itempk(iphas)).gt.0) then
    ipcvsl = ipproc(ivisls(itempk(iphas)))
  else
    ipcvsl = 0
  endif

! --- Stop si Lambda n'est pas variable

  if(ipcvsl.le.0) then
    write(nfecra,1010)                                            &
      itempk(iphas), itempk(iphas), ivisls(itempk(iphas))
    call csexit (1)
  endif

! --- Coefficients des lois choisis et imposes par l'utilisateur
!       Les valeurs donnees ici sont fictives

  varal = -3.3283d-7
  varbl =  3.6021d-5
  varcl =  1.2527d-4
  vardl =  0.58923d0



!   Lambda en W/(m K) au centre des cellules
!   ----------------------------------------
!       loi       Lambda        =
!               {  T  *( T  *( AL  * T +  BL  )+ CL  )+ DL }
!       soit    PROPCE(IEL,IPCVSL) =
!     &         (XRTP*(XRTP*(VARAL*XRTP+VARBL)+VARCL)+VARDL)

  do iel = 1, ncel
    xrtp = rtp(iel,ivart)
    propce(iel,ipcvsl) =                                          &
         (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)
  enddo


enddo
! --- Boucle sur les phases : fin


! --- A enlever a l'utilisation
if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif




!===============================================================================
!  EXEMPLE 5 : DIFFUSIVITE VARIABLE EN FONCTION DE LA TEMPERATURE
! ===========      POUR LES SCALAIRES UTILISATEURS
!     A l'exclusion de
!        temperature, enthalpie (traites plus haut)
!        variances de fluctuations (propriete egale a celle du
!                                                      scalaire associe)

!    Ci dessous on donne pour tous les scalaires (aux exclusions
!      ci-dessus pres) la meme loi pour la diffusivite
!    Les valeurs de cette propriete doivent etre fournies au centre des
!       cellules.
!  ===================================================================

! --- Boucle sur les phases : debut
do ii = 1, nscaus

! --- Numero du scalaire utilisateur II dans la liste de tous les scalaires
  iscal = ii


! --- S'il s'agit de la temperature,
!                                   son cas a deja ete traite plus haut
  ith = 0
  do iphas = 1, nphas
    if (iscal.eq.itempk(iphas)) ith = 1
  enddo

! --- Si la variable est une fluctuation, sa diffusivite est
!       la meme que celle du scalaire auquel elle est rattachee :
!       il n'y a donc rien a faire ici : on passe directement
!       a la variable suivante sans renseigner PROPCE(IEL,IPCVSL).

  if (ith.eq.0.and.iscavr(iscal).le.0) then
! --- On ne traite ici que les variables non thermiques
!                                   et qui ne sont pas des fluctuations


!   Positions des variables, coefficients
!   -------------------------------------

! --- Numero de variable température pour la phase courante iphas
!       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
!          IVART = ISCA(2)

  ivart = isca(itempk(iphas))


! --- Rang de Lambda du scalaire
!     dans PROPCE, prop. physiques au centre des elements       : IPCVSL

    if(ivisls(iscal).gt.0) then
      ipcvsl = ipproc(ivisls(iscal))
    else
      ipcvsl = 0
    endif

! --- Stop si Lambda n'est pas variable

    if(ipcvsl.le.0) then
      write(nfecra,1010) iscal, iscal, ivisls(iscal)
      call csexit (1)
    endif

! --- Coefficients des lois choisis et imposes par l'utilisateur
!       Les valeurs donnees ici sont fictives

    varal = -3.3283d-7
    varbl =  3.6021d-5
    varcl =  1.2527d-4
    vardl =  0.58923d0


!   Lambda en kg/(m s) au centre des cellules
!   ------------------------------------------
!       loi       Lambda           =
!                  T  *( T  *( AL  * T +  BL  )+ CL  )+ DL
!       soit    PROPCE(IEL,IPCVSL) =
!     &          XRTP*(XRTP*(VARAL*XRTP+VARBL)+VARCL)+VARDL


    do iel = 1, ncel
      xrtp = rtp(iel,ivart)
      propce(iel,ipcvsl) =                                        &
           (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)
    enddo


  endif
! --- Tests sur ITH et ISCAVR : fin

enddo
! --- Boucle sur les scalaires : fin


! --- A enlever a l'utilisation
if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif





!===============================================================================

!===============================================================================
! FORMATS
!----

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DONNEES DE CALCUL INCOHERENTES                          ',/,&
'@                                                            ',/,&
'@    Pour la phase ',I10                                      ,/,&
'@      usini1 indique que la capacite calorifique a pression ',/,&
'@      constante est uniforme ICP(',I10   ,') = ',I10   ,'   ',/,&
'@      alors que l''on cherche a definir une capacite        ',/,&
'@      calorifique variable dans uscfpv.                     ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usini1 ou uscfpv.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DONNEES DE CALCUL INCOHERENTES                          ',/,&
'@                                                            ',/,&
'@    Pour la phase ',I10                                      ,/,&
'@      uscfth indique que la capacite calorifique a volume   ',/,&
'@      constant  est uniforme ICV(',I10   ,') = ',I10   ,'   ',/,&
'@      alors que l''on cherche a definir une capacite        ',/,&
'@      calorifique variable dans uscfpv.                     ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usini1 ou uscfpv.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DONNEES DE CALCUL INCOHERENTES                          ',/,&
'@                                                            ',/,&
'@    Pour le scalaire ',I10                                   ,/,&
'@      usini1 indique que la diffusivite est uniforme        ',/,&
'@        IVISLS(',I10   ,') = ',I10   ,' alors que           ',/,&
'@      uscfpv impose une diffusivite variable.               ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usini1 ou uscfpv.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DONNEES DE CALCUL INCOHERENTES                          ',/,&
'@                                                            ',/,&
'@    Pour la phase ',I10                                      ,/,&
'@      uscfx2 indique que la viscosite en volume est uniforme',/,&
'@        IVISCV(',I10   ,') = ',I10   ,' alors que           ',/,&
'@      uscfpv impose une viscosite en volume variable.       ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier uscfx2 ou uscfpv.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    APPEL A csexit DANS LE SOUS PROGRAMME uscfpv            ',/,&
'@                                                            ',/,&
'@    Un appel a csexit (arret) a ete rencontre dans le sous  ',/,&
'@      programme uscfpv. L''utilisateur est invite a verifier',/,&
'@      que les exemples standard fournis par defaut ont bien ',/,&
'@      ete elimines si besoin.                               ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier uscfpv.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

return
end
