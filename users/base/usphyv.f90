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

subroutine usphyv &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , nphmx  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , ibrom  ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   ,                                     &
   propce , propfa , propfb ,                                     &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! ROUTINE UTILISATEUR : REMPLISSAGE DES VARIABLES PHYSIQUES



! ATTENTION :
! =========


! Il est INTERDIT de modifier la viscosite turbulente VISCT ici
!        ========
!  (une routine specifique est dediee a cela : usvist)


!  Il FAUT AVOIR PRECISE ICP(IPHAS) = 1
!     ==================
!    dans usini1 si on souhaite imposer une chaleur specifique
!    CP variable pour la phase IPHAS (sinon: ecrasement memoire).


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
!     - de la masse volumique                      ROM    kg/m3
!         (et eventuellememt aux faces de bord     ROMB   kg/m3)
!     - de la viscosite moleculaire                VISCL  kg/(m s)
!     - de la chaleur specifique associee          CP     J/(kg degres)
!     - des "diffusivites" associees aux scalaires VISCLS kg/(m s)


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
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! nphmx            ! e  ! <-- ! nphsmx                                         !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! ibrom            ! te ! <-- ! indicateur de remplissage de romb              !
!   (nphmx   )     !    !     !                                                !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1...8(ncelet    ! tr ! --- ! tableau de travail                             !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "pointe.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "entsor.h"
include "parall.h"
include "period.h"

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
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr), ibrom(nphmx)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision w5(ncelet),w6(ncelet),w7(ncelet),w8(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          ivart, iclvar, iel, iphas
integer          ipcrom, ipbrom, ipcvis, ipccp
integer          ipcvsl, ith, iscal, ii
integer          iutile
double precision vara, varb, varc, varam, varbm, varcm, vardm
double precision                   varal, varbl, varcl, vardl
double precision                   varac, varbc
double precision xrtp

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

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

!   Chaque exemple est encadre par un test sur IUTILE, par securite.
!   Mettre IUTILE a 1 pour activer l'exemple.

!   Il est conseille de ne garder dans ce sous programme que
!     le strict necessaire.



!  EXEMPLE 1 : MASSE VOLUMIQUE VARIABLE EN FONCTION DE LA TEMPERATURE
!  EXEMPLE 2 : VISCOSITE       VARIABLE EN FONCTION DE LA TEMPERATURE
!  EXEMPLE 3 : CHALEUR SPECIFIQUE VARIABLE EN FONCTION DE LA TEMPERATURE
!  EXEMPLE 4 : Lambda/CP  VARIABLE EN FONCTION DE LA TEMPERATURE
!                  POUR LA TEMPERATURE OU L'ENTHALPIE
!  EXEMPLE 5 : DIFFUSIVITE VARIABLE EN FONCTION DE LA TEMPERATURE
!                  POUR LES SCALAIRES
!===============================================================================










!===============================================================================
!  EXEMPLE 1 : MASSE VOLUMIQUE VARIABLE EN FONCTION DE LA TEMPERATURE
! ===========
!    Ci dessous on donne pour toutes les phases la meme loi pour
!       la masse volumique
!    Les valeurs de cette propriete doivent etre fournies au centre des
!       cellules (et, de facon optionnelle, aux faces de bord).
!  ===================================================================

!     Le test sur IUTILE permet de desactiver les instructions (qui
!       ne sont fournies qu'a titre d'exemple a adapter)

iutile = 0
if(iutile.eq.1) then

! --- Boucle sur les phases : debut
  do iphas = 1, nphas


!   Positions des variables, coefficients
!   -------------------------------------

! --- Numero de variable thermique pour la phase courante iphas
!       (et de ses conditions limites)
!       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
!          IVART = ISCA(2)

    if (iscalt(iphas).gt.0) then
      ivart = isca(iscalt(iphas))
    else
      write(nfecra,9010) iscalt(iphas)
      call csexit (1)
    endif

! --- Position des conditions limites de la variable IVART

    iclvar = iclrtp(ivart,icoef)

! --- Rang de la masse volumique de la phase courante IPHAS
!     dans PROPCE, prop. physiques au centre des elements       : IPCROM
!     dans PROPFB, prop. physiques au centre des faces de bord  : IPBROM

    ipcrom = ipproc(irom(iphas))
    ipbrom = ipprob(irom(iphas))

! --- Coefficients des lois choisis et imposes par l'utilisateur
!       Les valeurs donnees ici sont fictives

    vara  = -4.0668d-3
    varb  = -5.0754d-2
    varc  =  1000.9d0



!   Masse volumique au centre des cellules
!   ---------------------------------------
!       loi              RHO       =   T  * (  A *  T +  B ) +   C
!       soit    PROPCE(IEL,IPCROM) = XRTP * (VARA*XRTP+VARB) + VARC

    do iel = 1, ncel
      xrtp = rtp(iel,ivart)
      propce(iel,ipcrom) = xrtp * (vara*xrtp+varb) + varc
    enddo



!   Masse volumique aux faces de bord
!   ----------------------------------

!       Par defaut, la valeur de rho au bord est la valeur prise
!         au centre des elements voisins. C'est l'approche conseillee.
!       Pour etre dans ce cas il suffit de ne rien faire :
!         ne pas prescrire de valeur pour PROPFB(IFAC,IPBROM) et
!         ne pas modifier IBROM(IPHAS)
!         ---------------

!       Pour les utilisateurs qui ne souhaiteraient pas suivre ce
!         conseil, on precise que la temperature au bord peut etre
!         fictive, simplement destinee a conserver un flux (c'est
!         en particulier le cas en paroi). La valeur de rho calculee
!         au bord avec en introduisant cette temperature fictive
!         dans une loi physique peut donc etre totalement fausse
!         (negative par exemple).

!       Si malgre tout on souhaite imposer la loi :
!                        RHO       =   T  * (  A *  T +  B ) +   C
!       soit   PROPFB(IFAC,IPBROM) = XRTP * (VARA*XRTP+VARB) + VARC
!       T etant la temperature prise au centre des faces de bord,

!         on peut utiliser les lignes de code suivantes (volontairement
!         desactivees, car a manier avec precaution) :

!         Noter bien que dans le cas ou l'on impose la masse volumique
!         au bord, il faut le faire sur TOUTES les faces de bord.
!                                       ======

!          IBROM(IPHAS) = 1
!          DO IFAC = 1, NFABOR
!            IEL = IFABOR(IFAC)
!            XRTP = COEFA(IFAC,ICLVAR)+RTP(IEL,IVART)*COEFB(IFAC,ICLVAR)
!            PROPFB(IFAC,IPBROM) = XRTP * (VARA*XRTP+VARB) + VARC
!          ENDDO

!         IFABOR(IFAC) est l'element en regard de la face de bord

!         Attention IBROM(IPHAS) = 1 est indispensable pour que la loi
!           soit  prise en compte.       -------------



  enddo
! --- Boucle sur les phases : fin
endif
! --- Test sur IUTILE : fin






!===============================================================================
!  EXEMPLE 2 : VISCOSITE       VARIABLE EN FONCTION DE LA TEMPERATURE
! ===========
!    Ci dessous on donne pour toutes les phases la meme loi pour
!       la viscosite
!    Les valeurs de cette propriete doivent etre fournies au centre des
!       cellules.
!  ===================================================================

!     Le test sur IUTILE permet de desactiver les instructions (qui
!       ne sont fournies qu'a titre d'exemple a adapter)

iutile = 0
if(iutile.eq.1) then

! --- Boucle sur les phases : debut
  do iphas = 1, nphas


!   Positions des variables, coefficients
!   -------------------------------------

! --- Numero de variable thermique pour la phase courante iphas
!       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
!          IVART = ISCA(2)

    if (iscalt(iphas).gt.0) then
      ivart = isca(iscalt(iphas))
    else
      write(nfecra,9010) iscalt(iphas)
      call csexit (1)
    endif

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
      propce(iel,ipcvis) =                                        &
           xrtp*(xrtp*(varam*xrtp+varbm)+varcm)+vardm
    enddo


  enddo
! --- Boucle sur les phases : fin
endif
! --- Test sur IUTILE : fin





!===============================================================================
!  EXEMPLE 3 : CHALEUR SPECIFIQUE VARIABLE EN FONCTION DE LA TEMPERATURE
! ===========

!    Ci dessous on donne pour toutes les phases la meme loi pour
!       la chaleur specifique
!    Les valeurs de cette propriete doivent etre fournies au centre des
!       cellules.
!  ===================================================================

!     Le test sur IUTILE permet de desactiver les instructions (qui
!       ne sont fournies qu'a titre d'exemple a adapter)

iutile = 0
if(iutile.eq.1) then

! --- Boucle sur les phases : debut
  do iphas = 1, nphas


!   Positions des variables, coefficients
!   -------------------------------------

! --- Numero de variable thermique pour la phase courante iphas
!       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
!          IVART = ISCA(2)

    if (iscalt(iphas).gt.0) then
      ivart = isca(iscalt(iphas))
    else
      write(nfecra,9010) iscalt(iphas)
      call csexit (1)
    endif

! --- Rang de la chaleur specifique de la phase courante IPHAS
!     dans PROPCE, prop. physiques au centre des elements       : IPCCP

    if(icp(iphas).gt.0) then
      ipccp  = ipproc(icp   (iphas))
    else
      ipccp  = 0
    endif

! --- Stop si CP n'est pas variable

    if(ipccp.le.0) then
      write(nfecra,1000) iphas, iphas, icp(iphas)
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


  enddo
! --- Boucle sur les phases : fin
endif
! --- Test sur IUTILE : fin






!===============================================================================
!  EXEMPLE 4 : Lambda/CP  VARIABLE EN FONCTION DE LA TEMPERATURE
! ===========      POUR LA TEMPERATURE OU L'ENTHALPIE

!    Ci dessous on donne pour toutes les phases la meme loi pour
!       le rapport lambda/Cp
!    Les valeurs de cette propriete doivent etre fournies au centre des
!       cellules.
!  ===================================================================

!     Le test sur IUTILE permet de desactiver les instructions (qui
!       ne sont fournies qu'a titre d'exemple a adapter)

iutile = 0
if(iutile.eq.1) then

! --- Boucle sur les phases : debut
  do iphas = 1, nphas


!   Positions des variables, coefficients
!   -------------------------------------

! --- Numero de variable thermique pour la phase courante iphas
!       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
!          IVART = ISCA(2)

    if (iscalt(iphas).gt.0) then
      ivart = isca(iscalt(iphas))
    else
      write(nfecra,9010) iscalt(iphas)
      call csexit (1)
    endif

! --- Rang de Lambda/CP de la variable thermique de phase courante IPHAS
!     dans PROPCE, prop. physiques au centre des elements       : IPCVSL

    if(ivisls(iscalt(iphas)).gt.0) then
      ipcvsl = ipproc(ivisls(iscalt(iphas)))
    else
      ipcvsl = 0
    endif

! --- Stop si Lambda/CP n'est pas variable

    if(ipcvsl.le.0) then
      write(nfecra,1010)                                          &
           iscalt(iphas), iscalt(iphas), ivisls(iscalt(iphas))
      call csexit (1)
    endif

! --- Rang de la chaleur specifique de la phase courante IPHAS
!     dans PROPCE, prop. physiques au centre des elements       : IPCCP

    if(icp(iphas).gt.0) then
      ipccp  = ipproc(icp   (iphas))
    else
      ipccp  = 0
    endif

! --- Coefficients des lois choisis et imposes par l'utilisateur
!       Les valeurs donnees ici sont fictives

    varal = -3.3283d-7
    varbl =  3.6021d-5
    varcl =  1.2527d-4
    vardl =  0.58923d0



!   Lambda/Cp en kg/(m s) au centre des cellules
!   ------------------------------------------
!       loi       Lambda/CP        =
!               {  T  *( T  *( AL  * T +  BL  )+ CL  )+ DL } / Cp
!       soit    PROPCE(IEL,IPCVSL) =
!     &         (XRTP*(XRTP*(VARAL*XRTP+VARBL)+VARCL)+VARDL)/CP0(IPHAS)

!       On suppose Cp renseigne au prealable.

    if(ipccp.le.0) then

! --- Si CP est uniforme, on utilise CP0(IPHAS)
      do iel = 1, ncel
        xrtp = rtp(iel,ivart)
        propce(iel,ipcvsl) =                                      &
             (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)         &
             /cp0(iphas)
      enddo

    else

! --- Si CP est non uniforme, on utilise PROPCE ci dessus
      do iel = 1, ncel
        xrtp = rtp(iel,ivart)
        propce(iel,ipcvsl) =                                      &
             (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)         &
             /propce(iel,ipccp)
      enddo

    endif


  enddo
! --- Boucle sur les phases : fin
endif
! --- Test sur IUTILE : fin





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

!     Le test sur IUTILE permet de desactiver les instructions (qui
!       ne sont fournies qu'a titre d'exemple a adapter)

iutile = 0
if(iutile.eq.1) then

! --- Boucle sur les scalaires : debut
  do ii = 1, nscaus

! --- Numero du scalaire utilisateur II dans la liste de tous les scalaires
    iscal = ii


! --- S'il s'agit d'une variable thermique,
!                                   son cas a deja ete traite plus haut
    ith = 0
    do iphas = 1, nphas
      if (iscal.eq.iscalt(iphas)) ith = 1
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

! --- Numero de variable thermique pour la phase courante iphas
!       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
!          IVART = ISCA(2)

      iphas = iphsca(iscal)
      if (iscalt(iphas).gt.0) then
        ivart = isca(iscalt(iphas))
      else
        write(nfecra,9010) iscalt(iphas)
        call csexit (1)
      endif

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
        propce(iel,ipcvsl) =                                      &
             (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)
      enddo


    endif
! --- Tests sur ITH et ISCAVR : fin

  enddo
! --- Boucle sur les scalaires : fin
endif
! --- Test sur IUTILE : fin






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
'@      usini1 indique que la chaleur specifique est uniforme ',/,&
'@        ICP(',I10   ,') = ',I10   ,' alors que              ',/,&
'@      usphyv impose une chaleur specifique variable.        ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usini1 ou usphyv.                              ',/,&
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
'@      usphyv impose une diffusivite variable.               ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usini1 ou usphyv.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    APPEL A csexit DANS LE SOUS PROGRAMME usphyv            ',/,&
'@                                                            ',/,&
'@    La variable dont dependent les proprietes physiques ne  ',/,&
'@      semble pas etre une variable de calcul.               ',/,&
'@    En effet, on cherche a utiliser la temperature alors que',/,&
'@      ISCALT(IPHAS) = ',I10                                  ,/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage de usphyv (et le test lors de la     ',/,&
'@      definition de IVART).                                 ',/,&
'@    Verifier la definition des variables de calcul dans     ',/,&
'@      usini1. Si un scalaire doit jouer le role de la       ',/,&
'@      temperature, verifier que ISCALT a ete renseigne.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

return
end subroutine
