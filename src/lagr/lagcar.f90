!-------------------------------------------------------------------------------

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

subroutine lagcar &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   ,          &
   nprfml , nnod   , lndfac , lndfbr , ncelbr ,                   &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   itepa  , idevel , ituser , ia     ,                            &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod ,          &
   volume , dt     , rtp    , propce , propfa , propfb ,          &
   ettp   , ettpa  , tepa   , taup   , tlag   ,                   &
   piil   , bx     , tempct , statis ,                            &
   gradpr , gradvf , energi , dissip , romp   ,                   &
   rdevel , rtuser , ra        )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!    CALCUL DES CARACTERISTIQUES DES PARTICULES : Tp, TL et PI

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
! lndfac           ! e  ! <-- ! longueur du tableau nodfac                     !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr                     !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! te ! --- ! macro tableau entier                           !
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
! xyznod           ! tr ! <-- ! coordonnes des noeuds                          !
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
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! taup(nbpmax)     ! tr ! --> ! temps caracteristiques dynamique               !
! tlag(nbpmax)     ! tr ! --> ! temps caracteristiques fluide                  !
! piil(nbpmax,3    ! tr ! --> ! terme dans l'integration des eds up            !
! bx(nbpmax,3,2    ! tr ! --> ! caracteristiques de la turbulence              !
! tempct           ! tr ! --> ! temps caracteristique thermique                !
!   (nbpmax,2)     !    !     !                                                !
! statis(ncelet    ! tr ! <-- ! cumul des statistiques volumiques              !
!    nvlsta)       !    !     !                                                !
! gradpr           ! tr ! <-- ! gradient de pression                           !
! (ncelet,3)       !    !     !                                                !
! gradvf           ! tr ! <-- ! gradient de la vitesse du fluide               !
! (ncelet,3)       !    !     !                                                !
! energi(ncelet    ! tr ! --- ! tableau de travail                             !
! dissip(ncelet    ! tr ! --- ! tableau de travail                             !
! romp(nbpmax)     ! tr ! --- ! tableau de travail                             !
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
include "cstnum.h"
include "cstphy.h"
include "optcal.h"
include "entsor.h"
include "lagpar.h"
include "lagran.h"
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
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nideve , nrdeve , nituse , nrtuse
integer          itepa(nbpmax,nivep)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision dt(ncelet) , rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision taup(nbpmax) , tlag(nbpmax,3)
double precision piil(nbpmax,3) , bx(nbpmax,3,2)
double precision tempct(nbpmax,2)
double precision statis(ncelet,nvlsta)
double precision gradpr(ncelet,3) , gradvf(ncelet,9)
double precision energi(ncelet) , dissip(ncelet), romp(nbpmax)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          iel , ip , id , igvx , igvy , igvz , ivt
integer          iromf , iphas

double precision cd1 , cd2 , rec , cl , c0 , cb , cbcb
double precision upart , vpart , wpart
double precision uflui , vflui , wflui
double precision uvwdif , tl , uvwr
double precision rep , d2 , d3 , fdr , d1s3 , d3s444 , d6spi
double precision bb1 , bb2 , bb3 , ktil , bx1 , bx2 , bx3
double precision vpmx , vpmy , vpmz
double precision r11 , r22 , r33
double precision xnul , rom , prt , fnus , xrkl , xcp

!===============================================================================

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

iphas = ilphas

cd1  = 0.15d0
cd2  = 0.687d0
rec  = 1000.d0
c0   = 2.1d0
cl   = 1.d0 / (0.5d0 + (3.d0/4.d0)*c0 )
cb   = 0.8d0
cbcb = 0.64d0

d6spi = 6.d0 / pi
d1s3 = 1.d0 / 3.d0
d3s444 = 0.44d0 * 3.d0 / 4.d0

! Pointeur sur la masse volumique en fonction de l'ecoulement

if ( ippmod(icp3pl).ge.0 .or. ippmod(icfuel).ge.0 ) then
  iromf = ipproc(irom1)
else
  iromf = ipproc(irom(iphas))
endif

! Calcul de la masse volumique

do ip = 1,nbpart
  if ( itepa(ip,jisor).gt.0 ) then
    d3 = ettp(ip,jdp) * ettp(ip,jdp) * ettp(ip,jdp)
    romp(ip) = ettp(ip,jmp) * d6spi / d3
  endif
enddo

!===============================================================================
! 2. CALCUL DE Tp ET DE Tc SI THERMIQUE
!===============================================================================

do ip = 1,nbpart

  if ( itepa(ip,jisor) .gt.0 ) then

    iel = itepa(ip,jisor)

    rom  = propce(iel,iromf)
    xnul = propce(iel,ipproc(iviscl(iphas))) / rom

    uvwr = sqrt( ( ettp(ip,juf) -ettp(ip,jup) )*                  &
                 ( ettp(ip,juf) -ettp(ip,jup) )                   &
               + ( ettp(ip,jvf) -ettp(ip,jvp) )*                  &
                 ( ettp(ip,jvf) -ettp(ip,jvp) )                   &
               + ( ettp(ip,jwf) -ettp(ip,jwp) )*                  &
                 ( ettp(ip,jwf) -ettp(ip,jwp) )  )

!--->  CALCUL DU REYNOLDS LOCAL

    rep  = uvwr * ettp(ip,jdp) / xnul

!--->  CALCUL DU COEFFICIENT DE TRAINEE

    d2 = ettp(ip,jdp) * ettp(ip,jdp)

    if (rep.le.rec) then
      fdr = 18.d0 * xnul * (1.d0 + cd1 * rep**cd2) / d2
    else
      fdr = d3s444 * uvwr / ettp(ip,jdp)
    endif

!--->  CALCUL DE Tp

    taup(ip) = romp(ip) / rom / fdr

!--->  CALCUL UTILISATEUR DE Tp

    call uslatp                                                   &
    !==========
     ( idebia , idebra ,                                          &
       ndim   , ncelet , ncel   , nfac   , nfabor , nfml   ,      &
       nprfml , nnod   , lndfac , lndfbr , ncelbr ,               &
       nvar   , nscal  , nphas  ,                                 &
       nbpmax , nvp    , nvp1   , nvep   , nivep  ,               &
       nideve , nrdeve , nituse , nrtuse ,                        &
       ip     , itepa  , idevel , ituser , ia     ,               &
       rep    , uvwr   , rom    , romp(ip) , xnul , taup(ip) ,    &
       xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod ,      &
       volume , dt     , rtp    , propce , propfa , propfb ,      &
       ettp   , ettpa  , tepa   ,                                 &
       rdevel , rtuser , ra     )

!--->  CALCUL DE Tc

    if ( (iphyla.eq.1 .and. itpvar.eq.1) .or.                     &
         (iphyla.eq.2)                        ) then

!     CP fluide

      if (icp(iphas).gt.0) then
        xcp = propce( 1,ipproc(icp(iphas)) )
      else
        xcp = cp0(iphas)
      endif

!     CALCUL DU NUSSELT LOCAL

      if ( ippmod(icp3pl).ge.0 .or.                               &
           ippmod(icpl3c).ge.0 .or.                               &
           ippmod(icod3p).ge.1 .or.                               &
           ippmod(icoebu).eq.1 .or.                               &
           ippmod(icoebu).eq.3 .or.                               &
           ippmod(icfuel).ge.0 .or.                               &
           ippmod(ielarc).ge.0 .or.                               &
           ippmod(ieljou).ge.0      ) then
        ivt = ihm
      else
        ivt = iscalt(iphas)
      endif

! a priori en combustion gaz ou CP, la diffusvite est toujours constante

      if (ippmod(icoebu).eq.0 .or. ippmod(icoebu).eq.2) then
        xrkl = diftl0 / rom
      else if (ivisls(ivt).ge.1) then
        xrkl = propce(iel,ipproc(ivisls(ivt))) / rom
      else
        xrkl = visls0(ivt) / rom
      endif

      prt  = xnul / xrkl
      fnus = 2.d0 + 0.55d0 * rep**0.5d0 * prt**(d1s3)

! Calcul du temps caracteristique thermique Tc

      tempct(ip,1) = d2 * romp(ip) * ettp(ip,jcp)                 &
                   / ( fnus * 6.d0 * rom * xcp * xrkl )

!--->  CALCUL UTILISATEUR DE Tc

    call uslatc                                                   &
    !==========
     ( idebia , idebra ,                                          &
       ndim   , ncelet , ncel   , nfac   , nfabor , nfml   ,      &
       nprfml , nnod   , lndfac , lndfbr , ncelbr ,               &
       nvar   , nscal  , nphas  ,                                 &
       nbpmax , nvp    , nvp1   , nvep   , nivep  ,               &
       nideve , nrdeve , nituse , nrtuse ,                        &
       ip     , itepa  , idevel , ituser , ia     ,               &
       rep    , uvwr   , rom    , romp(ip) , xnul ,               &
       xcp    , xrkl   , tempct(ip,1) ,                           &
       xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod ,      &
       volume , dt     , rtp    , propce , propfa , propfb ,      &
       ettp   , ettpa  , tepa   ,                                 &
       rdevel , rtuser , ra     )

! Terme source implicite pour le couplage retour thermique

      tempct(ip,2) = fnus * pi * ettp(ip,jdp) * xrkl * rom

    endif

  endif

enddo


!===============================================================================
! 3. CALCUL DE TL
!===============================================================================

!--> Calcul de l'energie turbulente et de la dissipation
!      en fonction du modele de turbulence

if (idistu.eq.1) then

  if (itytur(iphas).eq.2 .or. iturb(iphas).eq.50) then
    do iel = 1,ncel
      energi(iel) = rtp(iel,ik(iphas))
      dissip(iel) = rtp(iel,iep(iphas))
    enddo
  else if (itytur(iphas).eq.3) then
    do iel = 1,ncel
      energi(iel) = 0.5d0*( rtp(iel,ir11(iphas))                  &
                           +rtp(iel,ir22(iphas))                  &
                           +rtp(iel,ir33(iphas)) )
      dissip(iel) = rtp(iel,iep(iphas))
    enddo
  else if (iturb(iphas).eq.60) then
    do iel = 1,ncel
      energi(iel) = rtp(iel,ik(iphas))
      dissip(iel) = cmu*energi(iel)*rtp(iel,iomg(iphas))
    enddo
  else
    write(nfecra,2000) iilagr, idistu, iphas, iturb(iphas)
    call csexit (1)
!              ======
  endif

!--> Calcul de TL et BX

  do ip = 1,nbpart

    if (itepa(ip,jisor).gt.0) then

      iel = itepa(ip,jisor)


      if (dissip(iel).gt.0.d0 .and.                               &
          energi(iel).gt.0.d0              ) then

      tl = cl * energi(iel) / dissip(iel)
      tl = max(tl,epzero)

      upart = ettp(ip,jup)
      vpart = ettp(ip,jvp)
      wpart = ettp(ip,jwp)
      uflui = ettp(ip,juf)
      vflui = ettp(ip,jvf)
      wflui = ettp(ip,jwf)

      if (modcpl.gt.0 .and. iplas.gt.modcpl) then
        if (statis(iel,ilpd).gt.seuil) then
          upart = statis(iel,ilvx) / statis(iel,ilpd)
          vpart = statis(iel,ilvy) / statis(iel,ilpd)
          wpart = statis(iel,ilvz) / statis(iel,ilpd)
          uflui = rtp(iel,iu(iphas))
          vflui = rtp(iel,iv(iphas))
          wflui = rtp(iel,iw(iphas))
        endif
      endif

      uvwdif = (uflui-upart) * (uflui-upart)                      &
             + (vflui-vpart) * (vflui-vpart)                      &
             + (wflui-wpart) * (wflui-wpart)

      uvwdif = (3.d0 * uvwdif) / (2.d0 *energi(iel))

      if (modcpl.gt.0 .and. iplas.gt.modcpl) then

        if (idirla.eq.1) then
          bb1 = sqrt( 1.d0 + cbcb *uvwdif )
          tlag(ip,1)= tl / bb1
          bb2  = sqrt( 1.d0 + 4.d0 *cbcb *uvwdif )
          tlag(ip,2)= tl / bb2
          bb3  = sqrt( 1.d0 + 4.d0 *cbcb *uvwdif )
          tlag(ip,3)= tl / bb3

        else if (idirla.eq.2) then
          bb1  = sqrt( 1.d0 + 4.d0 *cbcb *uvwdif )
          tlag(ip,1)= tl / bb1
          bb2  = sqrt( 1.d0 + cbcb *uvwdif )
          tlag(ip,2)= tl / bb2
          bb3  = sqrt( 1.d0 + 4.d0 *cbcb *uvwdif )
          tlag(ip,3)= tl / bb3

        else if (idirla.eq.3) then
          bb1  = sqrt( 1.d0 + 4.d0 *cbcb *uvwdif )
          tlag(ip,1)= tl / bb1
          bb2  = sqrt( 1.d0 + 4.d0 *cbcb *uvwdif )
          tlag(ip,2)= tl / bb2
          bb3  = sqrt( 1.d0 + cbcb *uvwdif )
          tlag(ip,3)= tl / bb3
        else
          write(nfecra,1000) idirla
          call csexit(1)
!                    ======
        endif

        if (itytur(iphas).eq.3) then
          r11 = rtp(iel,ir11(iphas))
          r22 = rtp(iel,ir22(iphas))
          r33 = rtp(iel,ir33(iphas))
          ktil = 3.d0 * ( r11*bb1 + r22*bb2 + r33*bb3  )          &
               / (2.d0 * (bb1+bb2+bb3) )
        else if (itytur(iphas).eq.2 .or. iturb(iphas).eq.50       &
                                    .or. iturb(iphas).eq.60) then
          ktil = energi(iel)
        endif

        bx1 = dissip(iel) * ( (c0*bb1*ktil/energi(iel))           &
                 +(2.d0 *(bb1*ktil/energi(iel) -1.d0)/3.d0) )
        bx2 = dissip(iel) * ( (c0*bb2*ktil/energi(iel))           &
                 +(2.d0 *(bb2*ktil/energi(iel) -1.d0)/3.d0) )
        bx3 = dissip(iel) * ( (c0*bb3*ktil/energi(iel))           &
                 +(2.d0 *(bb3*ktil/energi(iel) -1.d0)/3.d0) )

        if (bx1.gt.0.d0) then
          bx(ip,1,nor) = sqrt(bx1)
        else
          bx(ip,1,nor) = 0.d0
        endif

        if (bx2.gt.0.d0) then
          bx(ip,2,nor) = sqrt(bx2)
        else
          bx(ip,2,nor) = 0.d0
        endif

        if (bx3.gt.0.d0) then
          bx(ip,3,nor) = sqrt(bx3)
        else
          bx(ip,3,nor) = 0.d0
        endif

      else

        tlag(ip,1) = tl
        tlag(ip,2) = tl
        tlag(ip,3) = tl

        if (idiffl.eq.0) then
          uvwdif = sqrt(uvwdif)
          tlag(ip,1) = tl/(1.d0 + cb*uvwdif)
          tlag(ip,2) = tlag(ip,1)
          tlag(ip,3) = tlag(ip,1)
        endif

        bx(ip,1,nor) = sqrt(c0*dissip(iel))
        bx(ip,2,nor) = bx(ip,1,nor)
        bx(ip,3,nor) = bx(ip,1,nor)
      endif

      else

        tlag(ip,1) = epzero
        tlag(ip,2) = epzero
        tlag(ip,3) = epzero
        bx(ip,1,nor) = zero
        bx(ip,2,nor) = zero
        bx(ip,3,nor) = zero

      endif

    endif

  enddo

else

  do ip = 1,nbpart

    if ( itepa(ip,jisor) .gt.0 ) then
      tlag(ip,1) = epzero
      tlag(ip,2) = epzero
      tlag(ip,3) = epzero
      bx(ip,1,nor) = zero
      bx(ip,2,nor) = zero
      bx(ip,3,nor) = zero
    endif

  enddo

endif

!===============================================================================
! 4. CALCUL DE PII
!===============================================================================

do id = 1,3

  igvx = 3*(id-1)+1
  igvy = 3*(id-1)+2
  igvz = 3*(id-1)+3

  do ip = 1,nbpart

    if ( itepa(ip,jisor).gt.0 ) then

!--->   Calcul de II = ( -grad(P)/Rom(f)+grad(<Vf>)*(<Up>-<Uf>) )
!       ou
!       Calcul de II = ( -grad(P)/Rom(f) )

      iel = itepa(ip,jisor)

      piil(ip,id) = gradpr(iel,id)

      if (modcpl.gt.0 .and. iplas.gt.modcpl) then
        if ( statis(iel,ilpd) .gt. seuil ) then
          vpmx = statis(iel,ilvx) / statis(iel,ilpd)
          vpmy = statis(iel,ilvy) / statis(iel,ilpd)
          vpmz = statis(iel,ilvz) / statis(iel,ilpd)

          uflui = rtp(iel,iu(iphas))
          vflui = rtp(iel,iv(iphas))
          wflui = rtp(iel,iw(iphas))

          piil(ip,id) = gradpr(iel,id)                            &
                       +gradvf(iel,igvx) * (vpmx-uflui)           &
                       +gradvf(iel,igvy) * (vpmy-vflui)           &
                       +gradvf(iel,igvz) * (vpmz-wflui)

        endif
      endif

!--->  Terme purement explicite : probleme avec petit diametre
!      NE PAS EFFACER SVP           !

!            IF (IILAGR.EQ.2 .AND. ISTALA.GE.1) THEN

!              IF (STATIS(IEL,ILPD).GT.SEUIL) THEN

!                ROM   = PROPCE(IEL,IROMF)

!                FF = ROMP(IP) / ROM
!     &             *( STATIS(IEL,ILFV) / (DBLE(NPST)*VOLUME(IEL)) )
!     &             *( ETTP(IP,JUF+(ID-1)) - ETTP(IP,JUP+(ID-1)))
!     &                   /TAUP(IP)

!                PIIL(IP,ID) = PIIL(IP,ID) - FF

!              ENDIF

!            ENDIF

    endif

  enddo

enddo

!==============================================================================

!--------
! FORMATS
!--------

 1000   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE CHOIX DE LA DIRECTION DU MODELE COMPLET              ',/,&
'@       A UNE VALEUR NON PERMISE (LAGCAR).                   ',/,&
'@                                                            ',/,&
'@    IDIRLA DEVRAIT ETRE UN ENTIER EGAL A 1 2 OU 3           ',/,&
'@       (LA VALEUR 1 POUR UN ECOULEMENT SELON L''AXE X,      ',/,&
'@        LA VALEUR 2 POUR UN ECOULEMENT SELON L''AXE Y,      ',/,&
'@        LA VALEUR 3 POUR UN ECOULEMENT SELON L''AXE Z)      ',/,&
'@       IL VAUT ICI IDIRLA = ', I10                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IDIRLA dans la subroutine USLAG1.   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2000   format(                                                         &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    LE MODULE LAGRANGIEN EST INCOMPATIBLE AVEC LE MODELE    ',/,&
'@    DE TURBULENCE SELECTIONNE.                              ',/,&
'@                                                            ',/,&
'@   Le module Lagrangien a ete active avec IILAGR = ',I10     ,/,&
'@     et la dispersion turbulente est prise en compte        ',/,&
'@                                     avec IDISTU = ',I10     ,/,&
'@   Le modele de turbulence active pour la phase ',I6         ,/,&
'@     correspond a ITURB  = ',I10                             ,/,&
'@   Or, les seuls traitements de la turbulence compatibles   ',/,&
'@     avec le module Lagrangien et la dispersion turbulente  ',/,&
'@     sont k-epsilon et Rij-epsilon, v2f et k-omega.         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IILAGR et IDISTU dans la subroutine ',/,&
'@  USLAG1 et verifier la valeur de ITURB  dans la subroutine ',/,&
'@  USINI1.                                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

end
