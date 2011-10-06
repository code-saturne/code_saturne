!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

subroutine visecv &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   vela   , propce , propfa , propfb ,                            &
   coefa  , coefb  , cofafu , cofbfu , ckupdc , smacel ,          &
   trav   ,                                                       &
   viscf  , viscb  )

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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! vela             ! tr ! <-- ! variables de calcul au centre des              !
! (3,ncelet)       !    !     !    cellules (instant prec)                     !
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
! trav(3,ncelet)   ! tr ! <-- ! tableau de travail pour sec mem                !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use cstphy
use entsor
use numvar
use optcal
use pointe, only: forbr
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision vela(3,ncelet)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision trav(3,ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision cofafu(3,ndimfb)
double precision cofbfu(3,3,ndimfb)

! Local variables

integer          iccocg, inc, iel, ifac, ivar, isou, jsou, ii, jj, init
integer          nswrgp, imligp, iwarnp
integer          ipcrom, ipbrom, ipcvis, ipcvst, iflmas, iflmab
integer          ipcvsv
logical          ilved

double precision epsrgp, climgp, extrap
double precision romf, d2s3m, vecfac

double precision, allocatable, dimension(:) :: vistot
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: w4, w6
double precision, dimension(:,:,:), allocatable :: gradv, gradva

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! Allocate temporary arrays
allocate(vistot(ncelet))

ipcrom = ipproc(irom  )
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)

if(ippmod(icompf).ge.0) then
  if(iviscv.gt.0) then
    ipcvsv = ipproc(iviscv)
  else
    ipcvsv = 0
  endif
else
  ipcvsv = -1
endif


iflmas = ipprof(ifluma(iu))

ipbrom = ipprob(irom  )
iflmab = ipprob(ifluma(iu))


!     Si on extrapole les termes sources, on prend les prop a l'instant n
if(isno2t.gt.0) then
  if(iroext.gt.0) then
    ipcrom = ipproc(iroma )
    ipbrom = ipprob(iroma )
  endif
  if(iviext.gt.0) then
    ipcvis = ipproc(ivisla)
    ipcvst = ipproc(ivista)
  endif
!     Il faudrait aussi faire quelque chose pour le flux de masse, non ?
endif

! --- Calcul de la viscosite totale

if (itytur.eq.3) then
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

! Allocate a temporary array for the gradient calculation
allocate(w4(ncelet), w6(ncelet))
allocate(gradv(3,3,ncelet), gradva(3,3,ncelet))

! Les coefficients pris sont ceux sur le flux de vitesse, qui sont
! eventuellement egaux a ces sur les gradients

! --- Calcul du gradient de la vitesse

inc    = 1
iccocg = 1
nswrgp = nswrgr(iu)
imligp = imligr(iu)
iwarnp = iwarni(iu)
epsrgp = epsrgr(iu)
climgp = climgr(iu)
extrap = extrag(iu)

ilved = .true.

! Vectorial gradient
call grdvec                                                       &
!==========
 ( iu     , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ilved  ,                                                       &
   vela   , cofafu , cofbfu ,                                     &
   gradv  )


! Ceci parait VRAIMENT faux et ne plus conserver la quantite de mouvement
do iel = 1, ncelet
  w6(iel) = 1.d0
enddo
do ifac = 1, nfabor
  w6(ifabor(ifac)) = 0.d0
enddo

do iel = 1, ncelet
! --- Assemblage sur les faces internes
  do isou = 1, 3

! On stoque vistot*gradv dans gradva
    do jsou = 1, 3
      gradva(isou,jsou,iel) = vistot(iel)*gradv(isou,jsou,iel)
    enddo
  enddo
enddo

! FIXME: Plutot que de faire (mu tGRADU)_ij=0.5*(mu tGRADU)_i+mu tGRADU)_j)
!        on pourrait faire (mu tGRADU)_ij=mu_ij*0.5*( tGRADU_i+ tGRADU_j)
!        pour etre consistant avec le terme grad U

! On initialise TRAV(NCEL+1, NCELET)
!     (valeur bidon, mais pas NaN : les calculs sur le halo sont
!      par principe denue de sens, sauf exception)

if(ncelet.gt.ncel) then
  do iel = ncel+1, ncelet
    do isou = 1, 3
      trav(isou,iel) = 0.d0
    enddo
  enddo
endif

do ifac = 1, nfac

  ii = ifacel(1,ifac)
  jj = ifacel(2,ifac)

  do isou = 1, 3
    do jsou = 1, 3
      vecfac = 0.5d0*surfac(jsou,ifac)*                      &
      (gradva(jsou,isou,ii)+gradva(jsou,isou,jj))

      trav(isou,ii) = trav(isou,ii) + vecfac*w6(ii)
      trav(isou,jj) = trav(isou,jj) - vecfac*w6(jj)
    enddo
  enddo
enddo

! Free memory
deallocate(w6)

!===============================================================================
! 3.  CALCUL DES TERMES EN DIV
!===============================================================================
!  Pour periodicite et parallelisme, ROM est echange dans phyvar.
!     ou apres avoir ete calcule dans cfmsvl en compressible

!  Ici pour l'ordre 2 en temps, il faudrait tout prendre en n...

! Allocate a temporary array
allocate(w1(ncelet))

! Il serait sans doute mieux ici de prendre la trace de (-2/3mu GRADU )

do ifac = 1, nfac
  ii = ifacel(1,ifac)
  jj = ifacel(2,ifac)
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
    w4(iel) = (viscv0 + d2s3m*vistot(iel)) * w1(iel)/volume(iel)
  enddo
else

  if( itytur.eq.4) then
    do iel = 1, ncel
      w4(iel) = d2s3m*propce(iel,ipcvis)*w1(iel)/volume(iel)
    enddo
  else
    do iel = 1, ncel
      w4(iel) = d2s3m*vistot(iel)*w1(iel)/volume(iel)
    enddo
  endif
endif

! Free memory
deallocate(vistot)
deallocate(w1)

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(w4)
  !==========
endif


do ifac = 1, nfac
  ii = ifacel(1,ifac)
  jj = ifacel(2,ifac)
  viscf (ifac) = (w4(ii)+w4(jj))*0.5d0
enddo



! --- Assemblage sur les faces internes

do ifac = 1, nfac
  ii = ifacel(1,ifac)
  jj = ifacel(2,ifac)
  do isou = 1, 3
    vecfac = surfac(isou,ifac)*viscf(ifac)
    trav(isou,ii) = trav(isou,ii) + vecfac
    trav(isou,jj) = trav(isou,jj) - vecfac
  enddo
enddo

! --- Assemblage sur les faces de bord

do ifac = 1, nfabor
  ii = ifabor(ifac)
  do isou = 1, 3
    trav(isou,ii) = trav(isou,ii) + surfbo(isou,ifac)*w4(ii)
  enddo
enddo

! --- Calcul des efforts aux bords (partie 4/5)

if (ineedf.eq.1) then
  do ifac = 1, nfabor
    ii = ifabor(ifac)
    do isou = 1, 3
      forbr(isou,ifac) = forbr(isou,ifac) + surfbo(isou,ifac)*w4(ii)
    enddo
  enddo
endif

! Free memory
deallocate(w4)
deallocate(gradv,gradva)

return
end subroutine
