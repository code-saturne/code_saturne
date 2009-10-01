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

subroutine ctphyv &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , nphmx  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , ibrom  , izfppp ,          &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!   REMPLISSAGE DES VARIABLES PHYSIQUES : Version Aerorefrigerants


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
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !  (optionnel)                                   !
! ibrom            ! te ! <-- ! indicateur de remplissage de romb              !
!   (nphmx   )     !    !     !                                                !
! izfppp           ! te ! <-- ! numero de zone de la face de bord              !
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
! w1...8(ncelet    ! tr ! --- ! tableau de travail                             !
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
include "cstnum.h"
include "cstphy.h"
include "entsor.h"
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
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr), ibrom(nphmx)
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
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision w5(ncelet),w6(ncelet),w7(ncelet),w8(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra, ifinia
integer          iel   , iphas
integer          ipcrom, ipcvis, ipccp , ipcray
integer          ipcvsl, ith   , iscal , ivart, ii
integer          iclvar
integer          iiii  , ipcsig, it
integer          iesp  , iesp1 , iesp2 , mode , isrrom
integer          maxelt, ils

double precision tp    , delt  , somphi, val
double precision rho   , r     , cpa   , cpe , cpv , del
double precision hv0 , hvti , rhoj , tti , xxi,  xsati , dxsati
double precision rho0 , t00 , p00 , t1

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
!===============================================================================
! 0 - INITIALISATIONS A CONSERVER
!===============================================================================

! --- Initialisation memoire

idebia = idbia0
idebra = idbra0

ipass = ipass + 1

!===============================================================================
! 1 - MASSE VOLUMIQUE
!===============================================================================

! --- Boucle sur les phases : debut
do iphas = 1, nphas


!   Positions des variables, coefficients
!   -------------------------------------

! --- Numero de variable thermique pour la phase courante iphas
!       (et de ses conditions limites)

  ivart = isca(itemp4)

! --- Position des conditions limites de la variable IVART

  iclvar = iclrtp(ivart,icoef)

! --- Rang de la masse volumique de la phase courante IPHAS
!     dans PROPCE, prop. physiques au centre des elements       : IPCROM

  ipcrom = ipproc(irom(iphas))

! --- Coefficients des lois choisis et imposes par l'utilisateur
!       Les valeurs donnees ici sont fictives

  rho0 = 1.293d0
  t00  = 273.15d0
  p00  = 101325.d0
  del  = 0.622d0
  t1   = 273.15d0
  r    = 8.3143d0


!   Masse volumique au centre des cellules
!   ---------------------------------------

  do iel = 1, ncel

    tti = rtp(iel,isca(itemp4))
    xxi = rtp(iel,isca(ihumid))

    call xsath(tti,xsati)
!         ==========

    if (xxi .le. xsati) then

      rho = rho0*t00/(tti+t1)*del/(del+xxi)

    else

      if (tti.le.0.d0) then
        rhoj = 917.0d0
      else
        rhoj = 998.36d0 - 0.4116d0*(tti-20.d0)                    &
               - 2.24d0*(tti-20.d0)*(tti-70.d0)/625.d0
      endif

      rho = 1.0d0/((tti+t1)*p00/(t00*p00*rho0*del)                &
            *(del+xsati)+(xxi-xsati)/rhoj)

    endif

      if (rho .lt. 0.1d0) rho = 0.1d0
      propce(iel,ipcrom) = rho

  enddo


enddo
! --- Boucle sur les phases : fin



!===============================================================================
! 2 - CHALEUR SPECIFIQUE
!===============================================================================

! --- Boucle sur les phases : debut
do iphas = 1, nphas


!   Positions des variables, coefficients
!   -------------------------------------

! --- Numero de variable thermique pour la phase courante iphas

  ivart = isca(itemp4)

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

  cpa    = 1006.0d0
  cpv    = 1831.0d0
  cpe    = 4179.0d0
  hv0    = 2501600.0d0


!   Chaleur specifique J/(kg degres) au centre des cellules
!   --------------------------------------------------------

  if (ippmod(iaeros).eq.1) then

    do iel = 1, ncel

      tti = rtp(iel,isca(itemp4))
      xxi = rtp(iel,isca(ihumid))

      call xsath(tti,xsati)
!           ==========

      if (xxi .le. xsati) then
        propce(iel,ipccp) = cpa + xxi*cpv
      else
        hvti = (cpv-cpe)*tti + hv0
        call dxsath(tti,dxsati)
        !==========
        propce(iel,ipccp) = cpa + xsati*cpv +                     &
                            (xxi-xsati)*cpe + dxsati*hvti
      endif

    enddo

  elseif (ippmod(iaeros).eq.2) then

    do iel = 1, ncel

      tti = rtp(iel,isca(itemp4))

      call xsath(tti,xsati)
!           ==========

      hvti = cpv*tti + hv0

      call dxsath(tti,dxsati)
      !==========

      propce(iel,ipccp) = cpa + xsati*cpv + dxsati*hvti

    enddo

  endif


enddo
! --- Boucle sur les phases : fin


!===============================================================================
! 3 - ON PASSE LA MAIN A L'UTILISATEUR
!===============================================================================



! La masse volumique au bord est traitee dans phyvar (recopie de la valeur
!     de la cellule de bord).

!--------
! FORMATS
!--------

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

!----
! FIN
!----

return
end subroutine
