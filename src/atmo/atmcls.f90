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

subroutine atmcls &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , iphas  , ifac   , iel    , &
   uk     , utau   , yplus  ,                                     &
   uet    ,                                                       &
   gredu  , q0     , e0     , rib    , lmo    ,                   &
   cfnnu ,  cfnns  , cfnnk  , cfnne  ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , icodcl ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    ,          propce , propfa , propfb , rcodcl , &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FUNCTION :
! --------
! Compute u*, q0, e0, (momentum, sensible heat and latent heat fluxes)
!   for a non neutral atmospheric surface layer using the explicit formula
!   developped for the ECMWF by Louis (1982)

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
! iphas            ! e  ! <-- ! numero de phase                                !
! ifac             ! e  ! <-- ! face de bord traitee                           !
! iel              ! e  ! <-- ! cellule de bord en regard de la face           !
!                  !    !     !  traitee                                       !
! uk               ! r  ! <-- ! vitesse de frottement cf entete                !
! utau             ! r  ! <-- ! vitesse moyenne tangentielle                   !
! yplus            ! r  ! <-- ! distance adim a la paroi                       !
!                  !    !     !  calculee au moyen de uk                       !
! uet              ! r  ! <-- ! vitesse de frottement cf entete                !
! gredu            ! r  ! --> ! reduced gravity for non horizontal wall        !
! q0, e0           ! r  ! <-- ! sensible and latent heat flux                  !
! rib, lmo         ! r  ! <-- ! Richardson number and Monin-Obukhov length     !
! coeffu,s,k,e     ! r  ! <-- ! non neutral correction coefficients for        !
!                  !    !     !   profiles of momentum scalar turbulence       !
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
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
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
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
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
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "atincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse
integer          iphas  , ifac   , iel

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icodcl(nfabor,nvar)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision uk, utau, yplus, uet
double precision gredu, rib, lmo, q0, e0
double precision cfnnu, cfnns, cfnnk,cfnne
double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra

double precision tpot1,tpot2,tpotv1,tpotv2
double precision rscp1,rscp2
double precision actu,actt,b,c,d
double precision fm,fh,fmden1,fmden2,fhden
double precision rugd,rugt,distbf

!===============================================================================



!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

idebia = idbia0
idebra = idbra0

b = 5.d0
c = 5.d0
d = 5.d0
rib = 0.d0
lmo = 999.d0
q0 = 0.d0
e0 = 0.d0

rugd = rcodcl(ifac,iu(iphas),3)
distbf = yplus*rugd
rugt = rcodcl(ifac,iv(iphas),3)
actu = xkappa/log((distbf+rugd)/rugd)
actt = xkappa/log((distbf+rugt)/rugt)

! prise en compte de l'humidite dans le rapport r/cp
!if ( ippmod(iatmos).eq.2 ) then
!   rscp1=(rair/cp0(iphas))*(1.+(rvsra-cpvcpa)*qvs(ifac))
!   rscp2=(rair/cp0(iphas))*(1.+(rvsra-cpvcpa)*qv(iel))
!else
    rscp1 = rair/cp0(iphas)
    rscp2 = rair/cp0(iphas)
!endif

  tpot1 = rcodcl(ifac,isca(iscalt(iphas)),1)
  tpot2 = rtp(iel,isca(iscalt(iphas)))

!     .........    ............................................
!     3.2 - compute virtual potential temperature at two levels
!     .........    ............................................
!if ( ippmod(iatmos).eq.2 ) then
!    tpotv1=tpot1*(1.+(rvsra-1.)*qvs(ifac))
!    tpotv2=tpot2*(1.+(rvsra-1.)*qv(iel))
!else
  tpotv1 = tpot1
  tpotv2 = tpot2
!endif

!     .........    .....................................
!     3.3 - compute layer average Richardson number
!     .........    .....................................

! NB: rib =0 if thermal flux conditions are imposed and tpot1 not defined
if (abs(utau).le.epzero.or.icodcl(ifac,isca(iscalt(iphas))).eq.3) then
 rib = 0.d0
else
 rib = 2.d0*gredu*distbf*(tpotv2-tpotv1)/(tpotv1+tpotv2)/utau/utau
endif

!     .........    ..................................................
!     3.4 - compute corrction factors based on ECMWF parametrisation
!         Louis (1982)
!     ...............................................................

  if (rib.ge.0.d0) then
     fm = 1./(1.+2.*b*rib/sqrt(1.+d*rib))
     fh = 1/(1.+3.*b*rib*sqrt(1.+d*rib))
  else
     fmden1 = (distbf+rugt)*abs(rib)/rugt
     fmden2 = 1.+3.d0*b*c*actu*actt*sqrt(fmden1)
     fm = 1.d0-2.*b*rib/fmden2
     fhden = 3.d0*b*c*actu*actt*sqrt((distbf+rugt)/rugt)
     fh = 1.d0-(3.d0*b*rib)/(1.d0+fhden*sqrt(abs(rib)))
  endif

  cfnnu = 1.d0/sqrt(fm)
  cfnns = sqrt(fm)/fh
  cfnnk = sqrt(1.d0-rib)  ! +correction with turbulent Prandtl
  cfnne = (1.d0-rib)/sqrt(fm)

!     ------------------------------------
!     4 - compute friction velocity  uet
!     ------------------------------------

  uet = actu * utau * sqrt(fm)

!     -----------------------------------------
!     5 - compute surface sensible heat flux q0
!     -----------------------------------------

  q0 = (tpot1-tpot2) * uet * actt * fh / sqrt(fm)

!     -----------------------------------------------
!     6 - compute Monin-Obukhov lenght (Garratt p 38)
!     -----------------------------------------------

if (abs(gredu*q0).le.epzero) then
  lmo = -99999.d0
else
  lmo = -uet**3*(t0(iphas)+tkelvi)/(xkappa*abs(gredu)*q0)
endif

!     ---------------------------------------
!     7 - compute surface latent heat flux e0
!     ---------------------------------------

!  if ( ippmod(iatmos).eq.2 ) then
!    e0 = (qvs(ifac)-qv(iel)) * uet * actt * fh / sqrt(fm)
!  endif

!----
! fin
!----

return
end subroutine
