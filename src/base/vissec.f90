!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

subroutine vissec &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iphas  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , icepdc , icetsm , itypsm , &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , smacel ,                            &
   trav   ,                                                       &
   viscf  , viscb  , vistot ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! AJOUT AU SECOND MEMBRE DES TERMES

!         GRAD( (K -2/3 MU) DIV(U) ) + DIV( MU (GRAD_TRANSPOSE(U)) )

! AVEC MU = MU_LAMINAIRE + MU_TURBULENT
!  ET  K = VISCOSITE EN VOLUME (NULLE EN GENERAL)

! DIV(U)            EST CALCULE A PARTIR DE FLUMAS/RHO
! GRAD_TRANSPOSE(U) EST UN GRADIENT CELLULE

! REMARQUES :
!  - Theoriquement le terme en div(u) devrait plutot etre calcule
!      par un gradient cellule, pour correspondre exactement au
!      terme en dUj/dxi. Mais comme la partie en dUi/dxj est
!      calculee completement autrement (gradient facette et implicitation)
!      de toute facon on n'aura jamais Trace(tau_ij)=0 exactement.
!  - Pour la meme raison, comme le terme en dUi/dxj est calcule sur les
!      elements de bord et pas celui en dUj/dxi, il est difficile de
!      traiter le terme en div(u) de maniere rigoureuse. Il est donc
!      conserve sur les elements de bord.
!  - En LES, le tenseur <(u-<u>)(u-<u>)> est modelise par mut <S>
!      et non pas par mut <S> - 2/3 mut Tr(<S>) Id + 2/3 k Id
!      de sorte que il n'apparait pas ici de mut div<u>
!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
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
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! iphas            ! i  ! <-- ! phase number                                   !
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
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
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
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant prec)                     !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! trav(ncelet,3    ! tr ! <-- ! tableau de travail pour sec mem                !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! vistot(ncelet    ! tr ! --- ! tableau de travail pour mu                     !
! w1...6(ncelet    ! tr ! --- ! tableau de travail                             !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "dimfbr.h"
include "paramx.h"
include "cstphy.h"
include "entsor.h"
include "numvar.h"
include "optcal.h"
include "vector.h"
include "pointe.h"
include "period.h"
include "parall.h"
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
integer          ncepdp , ncesmp
integer          nideve , nrdeve , nituse , nrtuse , iphas

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision trav(ncelet,3)
double precision viscf(nfac), viscb(nfabor)
double precision vistot(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          iccocg, inc, iel, ifac, ivar, isou, ii, jj, init
integer          idim
integer          iuiph, iviph, iwiph
integer          iclvar
integer          nswrgp, imligp, iwarnp
integer          ipcrom, ipbrom, ipcvis, ipcvst, iflmas, iflmab
integer          iphydp
integer          ipcvsv

double precision epsrgp, climgp, extrap
double precision romf, d2s3m, vecfac

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

iuiph  = iu(iphas)
iviph  = iv(iphas)
iwiph  = iw(iphas)

ipcrom = ipproc(irom  (iphas))
ipcvis = ipproc(iviscl(iphas))
ipcvst = ipproc(ivisct(iphas))


if(ippmod(icompf).ge.0) then
  if(iviscv(iphas).gt.0) then
    ipcvsv = ipproc(iviscv(iphas))
  else
    ipcvsv = 0
  endif
else
  ipcvsv = -1
endif


iflmas = ipprof(ifluma(iuiph))

ipbrom = ipprob(irom  (iphas))
iflmab = ipprob(ifluma(iuiph))


!     Si on extrapole les termes sources, on prend les prop a l'instant n
if(isno2t(iphas).gt.0) then
  if(iroext(iphas).gt.0) then
    ipcrom = ipproc(iroma (iphas))
    ipbrom = ipprob(iroma (iphas))
  endif
  if(iviext(iphas).gt.0) then
    ipcvis = ipproc(ivisla(iphas))
    ipcvst = ipproc(ivista(iphas))
  endif
!     Il faudrait aussi faire quelque chose pour le flux de masse, non ?
endif

! --- Calcul de la viscosite totale

if (itytur(iphas).eq.3) then
  do iel = 1, ncel
    vistot(iel) = propce(iel,ipcvis)
  enddo
else
  do iel = 1, ncel
    vistot(iel) = propce(iel,ipcvis) + propce(iel,ipcvst)
  enddo
endif

!    Pour la periodicite de rotation, il faut avoir calcule
!      le gradient avec grdcel. La seule solution consiste donc a
!      echanger VISTOT puis a faire le produit, y compris sur les
!      cellules halo (calcul sur le halo, exceptionnellement).
!    Pour le parallelisme, on s'aligne sur la sequence ainsi definie.

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(vistot)
  !==========
endif


!===============================================================================
! 2.  CALCUL DES TERMES EN GRAD_TRANSPOSE
!===============================================================================

do isou = 1, 3

  if (isou.eq.1) ivar = iuiph
  if (isou.eq.2) ivar = iviph
  if (isou.eq.3) ivar = iwiph

! Ceci pointe eventuellement sur ICLRTP(IVAR,ICOEF)
  iclvar = iclrtp(ivar,icoeff)

! --- Calcul du gradient de la vitesse

  iccocg = 1
  inc    = 1
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  iphydp = 0

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
   w4     , w4     , w4     ,                                     &
   rtpa(1,ivar)    , coefa(1,iclvar) , coefb(1,iclvar) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )


  do iel = 1, ncelet
    w6(iel) = 1.d0
  enddo
  do ifac = 1, nfabor
    w6(ifabor(ifac)) = 0.d0
  enddo

! --- Assemblage sur les faces internes

  do idim = 1, 3

! On a echange le gradient dans grdcel et vistot plus haut

    if(idim.eq.1) then
      do iel = 1, ncelet
        w4(iel) = vistot(iel)*w1(iel)
      enddo
    elseif(idim.eq.2) then
      do iel = 1, ncelet
        w4(iel) = vistot(iel)*w2(iel)
      enddo
    elseif(idim.eq.3) then
      do iel = 1, ncelet
        w4(iel) = vistot(iel)*w3(iel)
      enddo
    endif




! On initialise TRAV(NCEL+1, NCELET)
!     (valeur bidon, mais pas NaN : les calculs sur le halo sont
!      par principe denue de sens, sauf exception)
    if(ncelet.gt.ncel) then
      do iel = ncel+1, ncelet
        trav(iel,idim) = 0.d0
      enddo
    endif



    if(ivecti.eq.1) then

!CDIR NODEP
      do ifac = 1, nfac
        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)
!MO             VECFAC = SURFAC(ISOU,IFAC)
!MO     &                  *(POND(IFAC)*W4(II)+(1.D0-POND(IFAC))*W4(JJ))
        vecfac = surfac(isou,ifac)*(w4(ii)+w4(jj))*0.5d0
        trav(ii,idim) = trav(ii,idim) + vecfac*w6(ii)
        trav(jj,idim) = trav(jj,idim) - vecfac*w6(jj)
      enddo

    else

! VECTORISATION NON FORCEE
      do ifac = 1, nfac
        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)
!MO             VECFAC = SURFAC(ISOU,IFAC)
!MO     &                  *(POND(IFAC)*W4(II)+(1.D0-POND(IFAC))*W4(JJ))
        vecfac = surfac(isou,ifac)*(w4(ii)+w4(jj))*0.5d0
        trav(ii,idim) = trav(ii,idim) + vecfac*w6(ii)
        trav(jj,idim) = trav(jj,idim) - vecfac*w6(jj)
      enddo

    endif


! --- Assemblage sur les faces de bord

!MO          IF(IVECTB.EQ.1) THEN
!MOC
!MO!CDIR NODEP
!MO            DO IFAC = 1, NFABOR
!MO             II = IFABOR(IFAC)
!MO             TRAV(II,IDIM) = TRAV(II,IDIM) + SURFBO(ISOU,IFAC)*W4(II)
!MO            ENDDO
!MOC
!MO          ELSE
!MOC
!MOC VECTORISATION NON FORCEE
!MO            DO IFAC = 1, NFABOR
!MO             II = IFABOR(IFAC)
!MO             TRAV(II,IDIM) = TRAV(II,IDIM) + SURFBO(ISOU,IFAC)*W4(II)
!MO            ENDDO
!MOC
!MO          ENDIF

  enddo

enddo

!===============================================================================
! 3.  CALCUL DES TERMES EN DIV
!===============================================================================
!  Pour periodicite et parallelisme, ROM est echange dans phyvar.
!     ou apres avoir ete calcule dans cfmsvl en compressible

!  Ici pour l'ordre 2 en temps, il faudrait tout prendre en n...


do ifac = 1, nfac
  ii = ifacel(1,ifac)
  jj = ifacel(2,ifac)
!MO        ROMF = POND(IFAC)*PROPCE(II,IPCROM)
!MO     &      + (1.D0-POND(IFAC))*PROPCE(JJ,IPCROM)
  romf = (propce(ii,ipcrom)+propce(jj,ipcrom))*0.5d0
  viscf(ifac) = propfa(ifac,iflmas)/romf
enddo
do ifac = 1, nfabor
  viscb(ifac) = propfb(ifac,iflmab)/propfb(ifac,ipbrom)
enddo

init = 1
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
                                   ifacel,ifabor,viscf ,viscb ,w1)

d2s3m = -2.d0/3.d0


if(ipcvsv.gt.0) then
  do iel = 1, ncel
    w4(iel) = ( propce(iel,ipcvsv) + d2s3m*vistot(iel) )          &
            * w1(iel)/volume(iel)
  enddo
elseif(ipcvsv.eq.0) then
  do iel = 1, ncel
    w4(iel) = ( viscv0(iphas)      + d2s3m*vistot(iel) )          &
            * w1(iel)/volume(iel)
  enddo
else

  if( itytur(iphas).eq.4) then
    do iel = 1, ncel
      w4(iel) = d2s3m*propce(iel,ipcvis)*w1(iel)/volume(iel)
    enddo
  else
    do iel = 1, ncel
      w4(iel) = d2s3m*vistot(iel)*w1(iel)/volume(iel)
    enddo
  endif
endif

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(w4)
  !==========
endif


do ifac = 1, nfac
  ii = ifacel(1,ifac)
  jj = ifacel(2,ifac)
!MO        VISCF (IFAC) = (POND(IFAC)*W4(II)+(1.D0-POND(IFAC))*W4(JJ))
  viscf (ifac) = (w4(ii)+w4(jj))*0.5d0
enddo

do isou = 1, 3

  idim = isou


! --- Assemblage sur les faces internes

  if(ivecti.eq.1) then

!CDIR NODEP
    do ifac = 1, nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      vecfac = surfac(isou,ifac)*viscf(ifac)
      trav(ii,idim) = trav(ii,idim) + vecfac
      trav(jj,idim) = trav(jj,idim) - vecfac
    enddo

  else

! VECTORISATION NON FORCEE
    do ifac = 1, nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
      vecfac = surfac(isou,ifac)*viscf(ifac)
      trav(ii,idim) = trav(ii,idim) + vecfac
      trav(jj,idim) = trav(jj,idim) - vecfac
    enddo

  endif


! --- Assemblage sur les faces de bord

  if(ivectb.eq.1) then

!CDIR NODEP
    do ifac = 1, nfabor
     ii = ifabor(ifac)
     trav(ii,idim) = trav(ii,idim) + surfbo(isou,ifac)*w4(ii)
    enddo

  else

! VECTORISATION NON FORCEE
    do ifac = 1, nfabor
     ii = ifabor(ifac)
     trav(ii,idim) = trav(ii,idim) + surfbo(isou,ifac)*w4(ii)
    enddo

  endif

! --- Calcul des efforts aux bords (partie 4/5)

  if (ineedf.eq.1) then
    do ifac = 1, nfabor
      ii = ifabor(ifac)
      ra(iforbr+(ifac-1)*ndim+isou-1) =                           &
             ra(iforbr+(ifac-1)*ndim+isou-1)                      &
           + surfbo(isou,ifac)*w4(ii)
    enddo
  endif

enddo


return

end subroutine
