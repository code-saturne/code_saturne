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

subroutine csc2cl &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nvcp   , nvcpto , nfbcpl , nfbncp ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , itrifb , itypfb ,                                     &
   lfbcpl , lfbncp ,                                              &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   rvcpfb , pndcpl , dofcpl ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!         TRADUCTION DE LA CONDITION ITYPFB(*,*) = ICSCPL

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
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
!                  !    !     ! = 10  -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle non bloquee :             !
!                  !    !     !  prescrire une valeur de dirichlet en          !
!                  !    !     !  prevision pour les scalaires k, eps,          !
!                  !    !     !  scal en plus du neumann usuel                 !
! itrifb(nfabor    ! te ! <-- ! indirection pour tri des faces de brd          !
!  nphas      )    !    !     !                                                !
! itypfb(nfabor    ! te ! --> ! type des faces de bord                         !
!  nphas      )    !    !     !                                                !
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
!                  !    !     !  flux (negatif si gain) w/m2                   !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
! w1,2,3,4,5,6     ! tr ! --- ! tableaux de travail                            !
!  (ncelet         !    !     !  (calcul du gradient de pression)              !
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
include "cstnum.h"
include "entsor.h"
include "parall.h"
include "period.h"
include "cplsat.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse
integer          nvcp   , nvcpto
integer          nfbcpl , nfbncp

integer          ifacel(2,nfac)  , ifabor(nfabor)
integer          ifmfbr(nfabor)  , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1)  , nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icodcl(nfabor,nvar)
integer          lfbcpl(nfbcpl)  , lfbncp(nfbncp)
integer          itrifb(nfabor,nphas), itypfb(nfabor,nphas)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

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
double precision rvcpfb(nfbcpl,nvcpto), pndcpl(nfbcpl)
double precision dofcpl(3,nfbcpl)
double precision rdevel(nrdeve), rtuser(nrtuse)
double precision ra(*)

! VARIABLES LOCALES


integer          idebia, idebra
integer          ifac, iel,isou, iphas
integer          inc, iccocg, iphydp, iclvar, nswrgp, imligp
integer          iwarnp, ivar
integer          ipt
integer          iii

double precision epsrgp, climgp, extrap
double precision xp
double precision xip   , xiip  , yiip  , ziip
double precision xjp
double precision xipf, yipf, zipf, ipf

double precision xif, yif, zif, xopf, yopf, zopf
double precision gradi

!===============================================================================


idebia = idbia0
idebra = idbra0


!================================================================================
! 1.  TRADUCTION DU COUPLAGE EN TERMES DE CONDITIONS AUX LIMITES
!================================================================================

! On rappelle que les variables sont reçues dans l'ordre de VARPOS ;
! il suffit dont de boucler sur les variables.

do ivar = 1, nvcp

!   --- Calcul du gradient de la variable si celle-ci est interpolée
!         Les échanges pour le parallélisme et la pédiocité (PARCOM
!         et PERCOM) ont déjà été fait dans CSCPFB.
!         Inutile de les refaire.

  inc    = 1
  iccocg = 1
  iphydp = 0

  iclvar = iclrtp(ivar,icoef)
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml,  &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp,  &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume,  &
   w4     , w4     , w4     ,                                     &
   rtp(1,ivar) , coefa(1,iclvar) , coefb(1,iclvar) ,              &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )

!   --- Traduction en termes de condition limite pour les faces de bord localisées
!         --> CL type Dirichlet

  do ipt = 1, nfbcpl

    ifac = lfbcpl(ipt)
    iel  = ifabor(ifac)

!         Information de l'instance en cours interpolée en I'
    iii = idiipb-1+3*(ifac-1)
    xiip = ra(iii+1)
    yiip = ra(iii+2)
    ziip = ra(iii+3)

    xif = cdgfbo(1,ifac) -xyzcen(1,iel)
    yif = cdgfbo(2,ifac) -xyzcen(2,iel)
    zif = cdgfbo(3,ifac) -xyzcen(3,iel)

    xipf = cdgfbo(1,ifac)-xiip - xyzcen(1,iel)
    yipf = cdgfbo(2,ifac)-yiip - xyzcen(2,iel)
    zipf = cdgfbo(3,ifac)-ziip - xyzcen(3,iel)

    ipf = sqrt(xipf**2+yipf**2+zipf**2)


    iii = idiipb-1+3*(ifac-1)
    xiip = ra(iii+1)
    yiip = ra(iii+2)
    ziip = ra(iii+3)

    xopf = dofcpl(1,ipt)
    yopf = dofcpl(2,ipt)
    zopf = dofcpl(3,ipt)

!         Informations locales interpolees en I'/O'

    xip =  rtp(iel,ivar)                                          &
      + (w1(iel)*(xiip+xopf) + w2(iel)*(yiip+yopf) +              &
         w3(iel)*(ziip+zopf))

!         Informations recues de l'instance distante en J'/O'
    xjp = rvcpfb(ipt,ivar)


    gradi = (w1(iel)*xipf+w2(iel)*yipf+w3(iel)*zipf)/ipf

    do iphas = 1, nphas
      itypfb(ifac,iphas)  = icscpl
    enddo

    if(ivar.ne.ipr(1)) then
      icodcl(ifac,ivar  ) = 1
      rcodcl(ifac,ivar,1) = 0.5d0*(xip+xjp)
    else
      icodcl(ifac,ivar  ) = 3
      rcodcl(ifac,ivar,3) = -0.5d0*dt(iel)*(gradi+xjp)
    endif


  enddo

! --- Faces de bord non localisées
!       --> CL type Neumann homogène

  do ipt = 1, nfbncp

    ifac = lfbncp(ipt)

    do iphas = 1, nphas
      itypfb(ifac,iphas)  = icscpl
    enddo

    icodcl(ifac,ivar  ) = 3
    rcodcl(ifac,ivar,3) = 0.d0

  enddo

enddo


!----
! FORMAT
!----

!----
! FIN
!----

return
end subroutine
