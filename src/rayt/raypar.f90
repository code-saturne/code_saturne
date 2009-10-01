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

subroutine raypar &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb ,          &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , isothp , izfrap ,                                     &
   idevel , ituser , ia     ,                                     &
   tmin   , tmax   , tx     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefa  , coefb  ,                                              &
   rdevel , rtuser ,                                              &
   tparop , qincip , textp  , tintp  ,                            &
   xlamp  , epap   , epsp   ,                                     &
   hfconp , flconp , tempkp ,                                     &

   ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
!   -----------------------------------------

!   CALCUL DES TEMPERATURES DE PAROIS AVEC UN BILAN DE FLUX

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
! iphas            ! e  ! <-- ! numero de la phase                             !
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
! itypfb(nfabor    ! te ! <-- ! type des faces de bord                         !
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
! isothp(nfabor    ! te ! <-- ! liste des frontieres isothermes                !
! izfrap(nfabor    ! te ! <-- ! numero de zone des faces de bord               !
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
! surfbo           ! tr ! <-- ! vecteur surface des faces de bord              !
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
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! tparop(nfabor    ! tr ! --> ! temperature de paroi en kelvin                 !
! qincip(nfabor    ! tr ! <-- ! densite de flux radiatif aux bords             !
! textp(nfabor)    ! tr ! <-- ! temperature de bord externe                    !
!                  !    !     ! en degres celcius                              !
! tintp(nfabor)    ! tr ! <-- ! temperature de bord interne                    !
!                  !    !     ! en degres celcius                              !
! xlamp(nfabor)    ! tr ! <-- ! coefficient de conductivite thermique          !
!                  !    !     ! des facettes de paroi (w/m/k)                  !
! epap(nfabor)     ! tr ! <-- ! epaisseur des facettes de paroi (m)            !
! epsp(nfabor)     ! tr ! <-- ! emissivite des facettes de bord                !
! hfconp(nfabor    ! tr ! <-- ! coefficient d'echange fluide aux               !
!                  !    !     ! faces de bord                                  !
! flconp(nfabor    ! tr ! <-- ! densite de flux convectif aux faces            !
! tempkp(ncelet    ! tr ! <-- ! temperature en kelvin                          !
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
include "numvar.h"
include "entsor.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "pointe.h"
include "parall.h"
include "ppppar.h"
include "radiat.h"


!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , iphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml) , itypfb(nfabor)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse), ia(*)

integer          isothp(nfabor), izfrap(nfabor)
integer          icodcl(nfabor,nvar)

double precision tmin , tmax , tx
double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision coefa(nfabor,*), coefb(nfabor,*)

double precision tparop(nfabor), qincip(nfabor)
double precision textp(nfabor), tintp(nfabor)
double precision xlamp(nfabor), epap(nfabor), epsp(nfabor)
double precision hfconp(nfabor) , flconp(nfabor)
double precision tempkp(ncelet)

double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)


! VARIABLES LOCALES

integer          idebia , idebra

integer          ivart, ifac, iel, izone
integer          ifacmx, ifacmn
integer          n1min,n1max,nrelax,nmoins,nplus
integer          iitpim,iipgrn,iipref,iifgrn,iifref
integer          indtp(nozrdm)
integer          nbrval, indtpm
integer          nb1int  ,nb2int  ,nbrrdp
parameter       (nb1int=5,nb2int=5,nbrrdp=5)
integer          inttm1(nb1int),inttm2(nb2int)

double precision esl,epp,sigt3,sigt4
double precision tpmin,tpmax,rapp,rapmax,abrapp
double precision qcmax,qrmax
double precision qcmin,qrmin
double precision qrayt,qconv,qinci,detep
double precision surfbn,und0,tp4
double precision xtpmax,ytpmax,ztpmax,xtpmin,ytpmin,ztpmin
double precision tzomax(nozrdm),tzomin(nozrdm),tzomoy(nozrdm)
double precision flunet(nozrdm),radios(nozrdm),surft(nozrdm)
double precision rdptmp(nbrrdp)


!===============================================================================

!===============================================================================
! 0. GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. INITIALISATION
!===============================================================================


und0 = 1.d0

ivart = isca(iscalt(iphas))


tpmax = -grand
tpmin =  grand

rapmax = 0.d0
n1min  = 0
n1max  = 0
nrelax = 0
nmoins = 0
nplus  = 0
iitpim = 0
iipgrn = 0
iipref = 0
iifgrn = 0
iifref = 0

ifacmx = 0
ifacmn = 0

!===============================================================================
! 2. Decompte du nombre de couleurs differentes de facettes de bord
!    On suppose qu'on en a au moins une
!===============================================================================


do izone = 1, nozrdm
  indtp (izone) =      0
  tzomax(izone) = -grand
  tzomin(izone) =  grand
enddo

!===============================================================================
! 3. Calcul des temperatures de paroi
!===============================================================================

!     3.1) parois isothermes
!     ----------------------


do ifac = 1, nfabor
  if (isothp(ifac).eq.itpimp) then
    izone = izfrap(ifac)

!      Calcul
    tparop(ifac)  = tintp(ifac)

!      Max-Min
    qconv  = flconp(ifac)
    qinci  = qincip(ifac)
    sigt4  = stephn*tparop(ifac)**4
    epp    = epsp(ifac)
    qrayt  = epp *( qinci - sigt4 )

    if (tpmax.le.tparop(ifac)) then
      ifacmx = ifac
      tpmax = tparop(ifac)
      qcmax = qconv
      qrmax = qrayt
    endif

    if (tpmin.ge.tparop(ifac)) then
      ifacmn = ifac
      tpmin = tparop(ifac)
      qcmin = qconv
      qrmin = qrayt
    endif

    tzomax(izone) = max(tzomax(izone),tparop(ifac))
    tzomin(izone) = min(tzomin(izone),tparop(ifac))

!      Reperage pour impression
    iitpim = 1
    indtp(izone) = itpimp


  endif
enddo


!     3.2) parois grises ou noires (non reflechissantes)
!     --------------------------------------------------

do ifac = 1, nfabor
  if (isothp(ifac).eq.ipgrno) then
    izone = izfrap(ifac)

!       Calcul
    esl    = epap(ifac) /xlamp(ifac)
    qconv  = flconp(ifac)
    qinci  = qincip(ifac)
    epp    = epsp(ifac)

    sigt3  = stephn*tparop(ifac)**3
    qrayt  = epp *( qinci -sigt3*tparop(ifac))

    detep  = ( esl*(qconv+qrayt) - (tparop(ifac)-textp(ifac)) )   &
           / ( 1.d0+ 4.d0*esl*epp*sigt3 + esl*hfconp(ifac) )

    rapp   = detep /tparop(ifac)
    abrapp = abs(rapp)

!       Relaxation
    if (abrapp.ge.tx) then
      nrelax = nrelax +1
      tparop(ifac) = tparop(ifac) * (1.d0+ tx*sign(und0,rapp))
    else
      tparop(ifac) = tparop(ifac) + detep
    endif

    rapmax = max(rapmax,abrapp)

    if (rapp.le.0.d0) then
      nmoins = nmoins+1
    else
      nplus  = nplus +1
    endif

!       Clipping
    if (tparop(ifac).lt.tmin) then
      n1min = n1min +1
      tparop(ifac) = tmin
    endif
    if (tparop(ifac).gt.tmax) then
      n1max = n1max +1
      tparop(ifac) = tmax
    endif

!       Max-Min
    if (tpmax.le.tparop(ifac)) then
      ifacmx = ifac
      tpmax = tparop(ifac)
      qcmax = qconv
      qrmax = qrayt
    endif

    if (tpmin.ge.tparop(ifac)) then
      ifacmn = ifac
      tpmin = tparop(ifac)
      qcmin = qconv
      qrmin = qrayt
    endif

    tzomax(izone) = max(tzomax(izone),tparop(ifac))
    tzomin(izone) = min(tzomin(izone),tparop(ifac))

!       Reperage pour impression
    iipgrn = 1
    indtp(izone) = ipgrno

  endif
enddo


!     3.3) parois reflechissantes
!     ---------------------------


do ifac = 1, nfabor
  if (isothp(ifac).eq.iprefl) then
    izone = izfrap(ifac)

!       Calcul
    esl    = epap(ifac) /xlamp(ifac)
    qconv  = flconp(ifac)

    detep  = ( esl*qconv - (tparop(ifac)-textp(ifac)) )           &
           / ( 1.d0 + esl*hfconp(ifac) )

    rapp = detep /tparop (ifac)
    abrapp = abs(rapp)

!       Relaxation
    if (abrapp.ge.tx) then
      nrelax = nrelax +1
      tparop(ifac) = tparop(ifac) * (1.d0+ tx*sign(und0,rapp))
    else
      tparop(ifac) = tparop(ifac) + detep
    endif

    rapmax = max(rapmax,abrapp)

    if (rapp.le.0.d0) then
      nmoins = nmoins+1
    else
      nplus  = nplus +1
    endif

!       Clipping
    if (tparop(ifac).lt.tmin) then
      n1min = n1min +1
      tparop(ifac) = tmin
    endif
    if (tparop(ifac).gt.tmax) then
      n1max = n1max +1
      tparop(ifac) = tmax
    endif

!       Max-Min
    if (tpmax.le.tparop(ifac)) then
      ifacmx = ifac
      tpmax = tparop(ifac)
      qcmax = qconv
      qrmax = 0.d0
    endif

    if (tpmin.ge.tparop(ifac)) then
      ifacmn = ifac
      tpmin = tparop(ifac)
      qcmin = qconv
      qrmin = 0.d0
    endif

    tzomax(izone) = max(tzomax(izone),tparop(ifac))
    tzomin(izone) = min(tzomin(izone),tparop(ifac))

!       Reperage pour impression
    iipref = 1
    indtp(izone) = iprefl

  endif
enddo



!     3.4) parois a flux de conduction impose non reflechissante
!          si le flux impose est nul, la paroi est adiabatique
!          (la transmition de chaleur par rayonnement est
!          equilibree avec celle de la convection)


do ifac = 1, nfabor
  if (isothp(ifac).eq.ifgrno) then
    izone = izfrap(ifac)

!       Calcul
    qconv  = flconp(ifac)
    qinci  = qincip(ifac)
    epp    = epsp(ifac)

    sigt3  = stephn*tparop(ifac)**3
    qrayt  = epp *( qinci -sigt3*tparop(ifac))

    detep  = ( qconv + qrayt - rcodcl(ifac,ivart,3) )             &
           / ( 4.d0*epp*sigt3 + hfconp(ifac) )

    rapp   = detep /tparop (ifac)
    abrapp = abs(rapp)

!       Relaxation
    if (abrapp.ge.tx) then
      nrelax = nrelax +1
      tparop(ifac) = tparop(ifac) * (1.d0+ tx*sign(und0,rapp))
    else
      tparop(ifac) = tparop(ifac) + detep
    endif

    rapmax = max(rapmax,abrapp)

    if (rapp.le.0.d0) then
      nmoins = nmoins+1
    else
      nplus  = nplus +1
    endif

!       Clipping
    if (tparop(ifac).lt.tmin) then
      n1min = n1min +1
      tparop(ifac) = tmin
    endif
    if (tparop(ifac).gt.tmax) then
      n1max = n1max +1
      tparop(ifac) = tmax
    endif

!       Max-Min
    if (tpmax.le.tparop(ifac)) then
      ifacmx = ifac
      tpmax = tparop(ifac)
      qcmax = qconv
      qrmax = qrayt
    endif

    if (tpmin.ge.tparop(ifac)) then
      ifacmn = ifac
      tpmin = tparop(ifac)
      qcmin = qconv
      qrmin = qrayt
    endif

    tzomax(izone) = max(tzomax(izone),tparop(ifac))
    tzomin(izone) = min(tzomin(izone),tparop(ifac))

!       Reperage pour impression
    iifgrn = 1
    indtp(izone) = ifgrno

  endif
enddo


!     3.5) parois a flux de conduction imposee et reflechissante
!          equivalent a impose un flux total au fluide



do ifac = 1, nfabor
  if (isothp(ifac).eq.ifrefl) then
    izone = izfrap(ifac)

!       Calcul
    iel = ifabor(ifac)

    tparop(ifac)= hfconp(ifac)*tempkp(iel)-rcodcl(ifac,ivart,3)
    tparop(ifac)= tparop(ifac)/max(hfconp(ifac),epzero)

!       Clipping
    if (tparop(ifac).lt.tmin) then
      n1min = n1min +1
      tparop(ifac) = tmin
    endif
    if (tparop(ifac).gt.tmax) then
      n1max = n1max +1
      tparop(ifac) = tmax
    endif

!       Max-Min
    qconv  = flconp(ifac)
    qrayt  = 0.d0

    if (tpmax.le.tparop(ifac)) then
      ifacmx = ifac
      tpmax = tparop(ifac)
      qcmax = qconv
      qrmax = qrayt
    endif

    if (tpmin.ge.tparop(ifac)) then
      ifacmn = ifac
      tpmin = tparop(ifac)
      qcmin = qconv
      qrmin = qrayt
    endif

    tzomax(izone) = max(tzomax(izone),tparop(ifac))
    tzomin(izone) = min(tzomin(izone),tparop(ifac))

!       Reperage pour impression
    iifref = 1
    indtp(izone) = ifrefl

  endif
enddo


!===============================================================================
! 4. IMPRESSIONS
!===============================================================================



! Test pour savoir s'il y a des zones de paroi
!  (donc si des impressions sont necessaires)

if (irangp.ge.0) then
  call parimx(nozarm,indtp )
  !==========
endif

indtpm = 0
do izone = 1, nozrdm
  if (indtp(izone).ne.0) then
    indtpm = 1
  endif
enddo

! Si il y a des parois

if(indtpm.gt.0) then

!   Si on veut imprimer en niveau 1 au moins

  if (iimpar.ge.1) then

!      Calcul de TZOMOY FLUNET RADIOS SURFT

    do izone = 1, nozrdm
      tzomoy(izone) = zero
      flunet(izone) = zero
      radios(izone) = zero
      surft (izone) = zero
    enddo
    do ifac = 1, nfabor
      surfbn = ra(isrfbn-1+ifac)
      izone = izfrap(ifac)
      if (indtp(izone).ne.0) then
        tp4 = tparop(ifac)**4
        tzomoy(izone)= tzomoy(izone) + tparop(ifac)*surfbn
        flunet(izone)= flunet(izone)                              &
             + epsp(ifac) *(qincip(ifac) -stephn*tp4 )*surfbn
        radios(izone)= radios(izone)                              &
             - ( epsp(ifac)       *stephn*tp4                     &
             +(1.d0-epsp(ifac))*qincip(ifac)      )*surfbn
        surft (izone) = surft(izone) + surfbn
      endif
    enddo

    if(irangp.ge.0) then
      call parrsm(nozarm,tzomoy)
      !==========
      call parrsm(nozarm,flunet)
      !==========
      call parrsm(nozarm,radios)
      !==========
      call parrsm(nozarm,surft )
      !==========
    endif

    do izone = 1, nozrdm
      if (indtp(izone).ne.0) then
        tzomoy(izone) = tzomoy(izone)/surft(izone)
        radios(izone) = radios(izone)/surft(izone)
      endif
    enddo


!      Determination de la TPMAX TPMIN et des grandeurs associees

    if(ifacmx.gt.0) then
      xtpmax = xyzcen(1,ifabor(ifacmx))
      ytpmax = xyzcen(2,ifabor(ifacmx))
      ztpmax = xyzcen(3,ifabor(ifacmx))
    else
      xtpmax = 0.d0
      ytpmax = 0.d0
      ztpmax = 0.d0
    endif
    if(ifacmn.gt.0) then
      xtpmin = xyzcen(1,ifabor(ifacmn))
      ytpmin = xyzcen(2,ifabor(ifacmn))
      ztpmin = xyzcen(3,ifabor(ifacmn))
    else
      xtpmin = 0.d0
      ytpmin = 0.d0
      ztpmin = 0.d0
    endif

    if(irangp.ge.0) then

      rdptmp(1) = xtpmax
      rdptmp(2) = ytpmax
      rdptmp(3) = ztpmax
      rdptmp(4) = qcmax
      rdptmp(5) = qrmax
      nbrval = nbrrdp
      call parmxl(nbrval,tpmax,rdptmp)
      !==========
      xtpmax =  rdptmp(1)
      ytpmax =  rdptmp(2)
      ztpmax =  rdptmp(3)
      qcmax  =  rdptmp(4)
      qrmax  =  rdptmp(5)

      rdptmp(1) = xtpmin
      rdptmp(2) = ytpmin
      rdptmp(3) = ztpmin
      rdptmp(4) = qcmin
      rdptmp(5) = qrmin
      nbrval = nbrrdp
      call parmnl(nbrval,tpmin,rdptmp)
      !==========
      xtpmin =  rdptmp(1)
      ytpmin =  rdptmp(2)
      ztpmin =  rdptmp(3)
      qcmin  =  rdptmp(4)
      qrmin  =  rdptmp(5)

    endif


!      Determination des compteurs et autres

    if(irangp.ge.0) then

      call parmax(rapmax)
      !==========

      inttm2(1) = nmoins
      inttm2(2) = nplus
      inttm2(3) = n1min
      inttm2(4) = n1max
      inttm2(5) = nrelax
      nbrval = nb2int
      call parism(nbrval,inttm2)
      !==========
      nmoins = inttm2(1)
      nplus  = inttm2(2)
      n1min  = inttm2(3)
      n1max  = inttm2(4)
      nrelax = inttm2(5)

      call parrmx(nozarm,tzomax)
      !==========
      call parrmn(nozarm,tzomin)
      !==========

      inttm1(1) = iitpim
      inttm1(2) = iipgrn
      inttm1(3) = iipref
      inttm1(4) = iifgrn
      inttm1(5) = iifref
      nbrval    = nb1int
      call parimx(nbrval,inttm1)
      !==========
      iitpim = inttm1(1)
      iipgrn = inttm1(2)
      iipref = inttm1(3)
      iifgrn = inttm1(4)
      iifref = inttm1(5)

    endif


!      Impressions

    write(nfecra,1000)
    write(nfecra,1010)

    if (nrelax.gt.0) then
      write(nfecra,2010) tx*100.d0 ,nrelax
      write(nfecra,1010)
    endif

    if ( n1min.gt.0 .or. n1max.gt.0 ) then
      write(nfecra,2020)
      write(nfecra,2030) n1min
      write(nfecra,2040) n1max
      write(nfecra,1010)
    endif

    if (rapmax.gt.0 .or. nmoins.gt.0 .or. nplus.gt.0 ) then
      write(nfecra,2050) rapmax * 100.d0
      write(nfecra,2060) nmoins
      write(nfecra,2070) nplus
      write(nfecra,1010)
    endif

    if (iitpim.eq.1) then
      write(nfecra,3020)
      do izone = 1, nozrdm
        if (indtp(izone).eq.itpimp) then
          write(nfecra,3000) izone,                               &
               tzomax(izone)-tkelvi,                              &
               tzomin(izone)-tkelvi,                              &
               tzomoy(izone)-tkelvi,                              &
               flunet(izone)
        endif
      enddo
      write(nfecra,1010)
    endif

    if (iipgrn.eq.1) then
      write(nfecra,3010)
      do izone = 1, nozrdm
        if (indtp(izone).eq.ipgrno) then
          write(nfecra,3000) izone,                               &
               tzomax(izone)-tkelvi,                              &
               tzomin(izone)-tkelvi,                              &
               tzomoy(izone)-tkelvi,                              &
               flunet(izone)
        endif
      enddo
      write(nfecra,1010)
    endif

    if (iipref.eq.1) then
      write(nfecra,3030)
      do izone = 1, nozrdm
        if (indtp(izone).eq.iprefl) then
          write(nfecra,3000) izone,                               &
               tzomax(izone)-tkelvi,                              &
               tzomin(izone)-tkelvi,                              &
               tzomoy(izone)-tkelvi,                              &
               flunet(izone)
        endif
      enddo
      write(nfecra,1010)
    endif

    if (iifgrn.eq.1) then
      write(nfecra,3040)
      do izone = 1, nozrdm
        if (indtp(izone).eq.ifgrno) then
          write(nfecra,3000) izone,                               &
               tzomax(izone)-tkelvi,                              &
               tzomin(izone)-tkelvi,                              &
               tzomoy(izone)-tkelvi,                              &
               flunet(izone)
        endif
      enddo
      write(nfecra,1010)
    endif

    if (iifref.eq.1) then
      write(nfecra,3050)
      do izone = 1, nozrdm
        if (indtp(izone).eq.ifrefl) then
          write(nfecra,3000) izone,                               &
               tzomax(izone)-tkelvi,                              &
               tzomin(izone)-tkelvi,                              &
               tzomoy(izone)-tkelvi,                              &
               flunet(izone)
        endif
      enddo
      write(nfecra,1010)
    endif

  endif


!   Si on veut imprimer EN PLUS en niveau 2
!                       =======

  if (iimpar.ge.2) then

    write(nfecra,5010) tpmax-tkelvi
    write(nfecra,5020) xtpmax, ytpmax, ztpmax
    write(nfecra,5030) qcmax
    write(nfecra,5040) qrmax
    write(nfecra,5050) tpmin-tkelvi
    write(nfecra,5060) xtpmin, ytpmin, ztpmin
    write(nfecra,5070) qcmin
    write(nfecra,5080) qrmin

  endif

endif

! -------
! FORMATS
! -------

 1000 FORMAT (/, 3X,'** INFORMATIONS SUR LA TEMPERATURE DES PAROIS',/,  &
           3X,'   ------------------------------------------')

 1010 format('------------------------------------'                     &
         ,'-----------------------------------')

 2010 format('ATTENTION, temperature de paroi relaxee a ',G7.2,'%'      &
      ,' en (',I8,' points)')

 2020 format('ATTENTION, temperature de paroi CLIPPE AU MIN-MAX :')

 2030 format('NOMBRE DE POINTS CLIPPE AU MINIMUM : ',I8)

 2040 format('NOMBRE DE POINTS CLIPPE AU MAXIMUM : ',I8)

 2050 format('Variation maximale : ',G9.4,'%')

 2060 format('Temperature de paroi diminuant  : ',I8,' faces de bord')

 2070 format('Temperature de paroi augmentant : ',I8,' faces de bord')

 3000 format(i10,8x,e10.4,4x,e10.4,4x,e10.4,4x,e10.4)

 3010 format('Grises ou noires Temp max (C)  '                          &
      ,'Temp min (C)  Temp moy (C)  Flux Net (W)')

 3020 format('Profils imposes  Temp max (C)  '                          &
      ,'Temp min (C)  Temp moy (C)  Flux Net (W)')

 3030 format('Parois a EPS=0   Temp max (C)  '                          &
      ,'Temp min (C)  Temp moy (C)  Flux Net (W)')

 3040 format('Flux imp EPS!=0  Temp max (C)  '                          &
      ,'Temp min (C)  Temp moy (C)  Flux Net (W)')

 3050 format('Flux imp EPS=0   Temp max (C)  '                          &
      ,'Temp min (C)  Temp moy (C)  Flux Net (W)')

 5010 format(/,12X,'TEMPERATURE PAROI MAXIMALE (Degre celsius) = ',     &
         g15.7)
 5020 format(15X,' AU POINT X Y Z =  ',E10.4,2X,E10.4,2X,E10.4)

 5030 format(15X,' FLUX CONVECTIF  = ',G15.7)

 5040 format(15X,' FLUX RADIATIF  = ',G15.7)

 5050 format(/,12X,'TEMPERATURE PAROI MINIMALE (Degre celsius) = ',     &
         g15.7)

 5060 format(15X,' AU POINT X Y Z =  ',E10.4,2X,E10.4,2X,E10.4)

 5070 format(15X,' FLUX CONVECTIF  = ',G15.7)

 5080 format(15X,' FLUX RADIATIF  = ',G15.7,/)


! ---
! FIN
! ---

return

end subroutine

!   ---------------------------------------------------------------------------------
!   Variation maximale : .1484E-13%
!   Temperature de paroi diminuant  :       33
!   Temperature de paroi augmentant :       27
!   ---------------------------------------------------------------------------------
!   Zones de parois  Temp max (C)    Temp min (C)    Temp moy (C)    Flux Net (W)
!           15        0.7941E+03      0.7410E+03      0.7688E+03      0.7688E+03
!   ---------------------------------------------------------------------------------
!   Profils imposes  Temp max (C)    Temp min (C)    Temp moy (C)    Flux Net (W)
!           15        0.7941E+03      0.7410E+03      0.7688E+03      0.7688E+03
!   ---------------------------------------------------------------------------------

