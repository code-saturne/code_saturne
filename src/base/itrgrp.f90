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

subroutine itrgrp &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   init   , inc    , imrgra , iccocg , nswrgp , imligp ,          &
   iphydp , iwarnp , nfecra ,                                     &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   fextx  , fexty  , fextz  ,                                     &
   pvar   , coefap , coefbp , viscf  , viscb  ,                   &
   viselx , visely , viselz ,                                     &
   diverg ,                                                       &
   dpdx   , dpdy   , dpdz   , dpdxa  , dpdya  , dpdza  ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! INCREMENTATION DE LA DIVERGENCE DU FLUX DE MASSE
! A PARTIR DE GRAD(P)
! grad(P) = GRADIENT FACETTE

!  .          .        --         --->     -->
!  m   =  m     - \    Visc   grad(P) . n
!   ij     ij     /__      ij       ij   ij
!                  j

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
! init             ! e  ! <-- ! > 0 : initialisation du flux de masse          !
! inc              ! e  ! <-- ! indicateur = 0 resol sur increment             !
!                  !    !     !              1 sinon                           !
! imrgra           ! e  ! <-- ! indicateur = 0 gradrc 97                       !
!                  ! e  ! <-- !            = 1 gradmc 99                       !
! iccocg           ! e  ! <-- ! indicateur = 1 pour recalcul de cocg           !
!                  !    !     !              0 sinon                           !
! nswrgp           ! e  ! <-- ! nombre de sweep pour reconstruction            !
!                  !    !     !             des gradients                      !
! imligp           ! e  ! <-- ! methode de limitation du gradient              !
!                  !    !     !  < 0 pas de limitation                         !
!                  !    !     !  = 0 a partir des gradients voisins            !
!                  !    !     !  = 1 a partir du gradient moyen                !
! iwarnp           ! e  ! <-- ! niveau d'impression                            !
! iphydp           ! e  ! <-- ! indicateur de prise en compte de la            !
!                  !    !     ! pression hydrostatique                         !
! nfecra           ! e  ! <-- ! unite du fichier sortie std                    !
! epsrgp           ! r  ! <-- ! precision relative pour la                     !
!                  !    !     !  reconstruction des gradients 97               !
! climgp           ! r  ! <-- ! coef gradient*distance/ecart                   !
! extrap           ! r  ! <-- ! coef extrap gradient                           !
! pvar  (ncelet    ! tr ! <-- ! variable (pression)                            !
! coefap, b        ! tr ! <-- ! tableaux des cond lim pour pvar                !
!   (nfabor)       !    !     !  sur la normale a la face de bord              !
! viscf (nfac)     ! tr ! <-- ! "viscosite" face interne(dt*surf/dist          !
! viscb (nfabor    ! tr ! <-- ! "viscosite" face de bord(dt*surf/dist          !
! viselx(ncelet    ! tr ! <-- ! "viscosite" par cellule  dir x                 !
! visely(ncelet    ! tr ! <-- ! "viscosite" par cellule  dir y                 !
! viselz(ncelet    ! tr ! <-- ! "viscosite" par cellule  dir z                 !
! diverg(ncelet    ! tr ! <-- ! divergence du flux de masse                    !
! fextx,y,z        ! tr ! <-- ! force exterieure generant la pression          !
!   (ncelet)       !    !     !  hydrostatique                                 !
! dpd.(ncelet      ! tr ! --- ! tableau de travail pour le grad de p           !
! dpdxa,y,z        ! tr ! --- ! tableau de travail pour le grad de p           !
!    (ncelet       !    !     !  avec decentrement amont                       !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
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
include "period.h"
include "parall.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse
integer          init   , inc    , imrgra , iccocg
integer          nswrgp , imligp
integer          iwarnp , iphydp , nfecra
double precision epsrgp , climgp , extrap

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision pvar(ncelet), coefap(nfabor), coefbp(nfabor)
double precision viscf(nfac), viscb(nfabor)
double precision viselx(ncelet), visely(ncelet), viselz(ncelet)
double precision diverg(ncelet)
double precision fextx(ncelet),fexty(ncelet),fextz(ncelet)
double precision dpdx (ncelet),dpdy (ncelet),dpdz (ncelet)
double precision dpdxa(ncelet),dpdya(ncelet),dpdza(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ifac, ii, jj, iij, iii, ivar
integer          idimte, itenso
double precision pfac,pip
double precision dpxf  , dpyf  , dpzf  , flumas, flumab
double precision dijpfx, dijpfy, dijpfz
double precision diipbx, diipby, diipbz
double precision dist  , surfan
double precision dijx  , dijy  , dijz

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

if( init.ge.1 ) then
  do ii = 1, ncelet
    diverg(ii) = 0.d0
  enddo
elseif(init.eq.0.and.ncelet.gt.ncel) then
  do ii = ncel+1, ncelet
    diverg(ii) = 0.d0
  enddo
elseif(init.ne.0) then
  write(nfecra,1000) init
  call csexit (1)
endif

! ---> TRAITEMENT DU PARALLELISME

if(irangp.ge.0) call parcom (pvar)
                !==========

! ---> TRAITEMENT DE LA PERIODICITE

if(iperio.eq.1) then
  idimte = 0
  itenso = 0
  call percom                                                     &
  !==========
  ( idimte , itenso ,                                             &
    pvar   , pvar   , pvar  ,                                     &
    pvar   , pvar   , pvar  ,                                     &
    pvar   , pvar   , pvar  )
endif

!===============================================================================
! 2.  INCREMENT DU FLUX DE MASSE SS TECHNIQUE DE RECONSTRUCTION
!===============================================================================

if( nswrgp.le.1 ) then

!     FLUX DE MASSE SUR LES FACETTES FLUIDES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    flumas = viscf(ifac)*(pvar(ii) -pvar(jj))
    diverg(ii) = diverg(ii) + flumas
    diverg(jj) = diverg(jj) - flumas

  enddo

!     FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    pfac = inc*coefap(ifac) +coefbp(ifac)*pvar(ii)

    flumab = viscb(ifac)*( pvar(ii) -pfac )
    diverg(ii) = diverg(ii) + flumab

  enddo

endif


!===============================================================================
! 3.  INCREMENTATION DU FLUX DE MASSE AVEC TECHNIQUE DE
!         RECONSTRUCTION SI LE MAILLAGE EST NON ORTHOGONAL
!===============================================================================

if( nswrgp.gt.1 ) then

!     CALCUL DU GRADIENT

!     IVAR ne sert a GRDCEL que si la variable est une composante de la vitesse
!     ou de Rij pour la periodicite. Ici la variable est soit la pression
!     soit phi, donc on peut mettre IVAR=0
  ivar = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &

   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   fextx  , fexty  , fextz  ,                                     &
   pvar   , coefap , coefbp ,                                     &
   dpdx   , dpdy   , dpdz   ,                                     &
!        ------   ------   ------
   dpdxa  , dpdya  , dpdza  ,                                     &
   rdevel , rtuser , ra     )

! ---> TRAITEMENT DU PARALLELISME

if(irangp.ge.0) then
  call parcom (viselx)
  !==========
  call parcom (visely)
  !==========
  call parcom (viselz)
  !==========
endif

! ---> TRAITEMENT DE LA PERIODICITE

if(iperio.eq.1) then
  idimte = 1
  itenso = 0
  call percom                                                     &
  !==========
  ( idimte , itenso ,                                             &
    viselx , viselx , viselx ,                                    &
    visely , visely , visely ,                                    &
    viselz , viselz , viselz )
endif

!     FLUX DE MASSE SUR LES FACETTES FLUIDES

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    iij = idijpf-1+3*(ifac-1)

    dijpfx = ra(iij+1)
    dijpfy = ra(iij+2)
    dijpfz = ra(iij+3)
    dist   = ra(idist-1+ifac)
    surfan = ra(isrfan-1+ifac)

!---> DIJ = IJ - (IJ.N) N
    dijx = (xyzcen(1,jj)-xyzcen(1,ii))-dijpfx
    dijy = (xyzcen(2,jj)-xyzcen(2,ii))-dijpfy
    dijz = (xyzcen(3,jj)-xyzcen(3,ii))-dijpfz

    dpxf = 0.5d0*(viselx(ii)*dpdx(ii) + viselx(jj)*dpdx(jj))
    dpyf = 0.5d0*(visely(ii)*dpdy(ii) + visely(jj)*dpdy(jj))
    dpzf = 0.5d0*(viselz(ii)*dpdz(ii) + viselz(jj)*dpdz(jj))

    flumas = viscf(ifac)*( pvar(ii) -pvar(jj) )                   &
     + ( dpxf * dijx                                              &
     +   dpyf * dijy                                              &
     +   dpzf * dijz )*surfan/dist
    diverg(ii) = diverg(ii) + flumas
    diverg(jj) = diverg(jj) - flumas

  enddo

!     FLUX DE MASSE SUR LES FACETTES DE BORD

  do ifac = 1, nfabor

    ii = ifabor(ifac)

    iii = idiipb-1+3*(ifac-1)
    diipbx = ra(iii+1)
    diipby = ra(iii+2)
    diipbz = ra(iii+3)

    pip = pvar(ii) +                                              &
      dpdx(ii)*diipbx+dpdy(ii)*diipby+dpdz(ii)*diipbz
    pfac = inc*coefap(ifac) +coefbp(ifac)*pip

    flumab = viscb(ifac)*( pip -pfac )
    diverg(ii) = diverg(ii) + flumab

  enddo

endif

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format('ITRGRP APPELE AVEC INIT = ',I10)

#else

 1000 format('ITRGRP CALLED WITH INIT = ',I10)

#endif

!----
! FIN
!----

return

end subroutine
