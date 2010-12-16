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

subroutine cfdivs &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   iphas  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   rtp    , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , smacel ,                            &
   diverg , ux     , uy     , uz     ,                            &
   vistot ,                                                       &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   ra     )

!===============================================================================
! FONCTION :
! ----------
!                                   v
! CALCULE DIVERG = DIVERG + DIV(SIGMA .U)

!          v               t
! AVEC SIGMA = MU (GRAD(U) + GRAD(U)) + (KAPPA - 2/3 MU) DIV(U) Id

! ET MU = MU_LAMINAIRE + MU_TURBULENT

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
!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iphas            ! i  ! <-- ! phase number                                   !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! ia(*)            ! ia ! --- ! main integer work array                        !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant)                  !
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
! diverg(ncelet    ! tr ! --> ! div(sigma.u)                                   !
! ux,y,z(ncelet    ! tr ! <-- ! composantes du vecteur u                       !
! vistot(ncelet    ! tr ! --- ! tableau de travail pour mu                     !
! w1...6(ncelet    ! tr ! --- ! tableau de travail                             !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use cstphy
use entsor
use numvar
use optcal
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          iphas

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ia(*)

double precision rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision diverg(ncelet)
double precision ux(ncelet), uy(ncelet), uz(ncelet)
double precision vistot(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          iccocg, inc, iel, ifac, ivar, isou, ii, jj
integer          iuiph, iviph, iwiph
integer          iclvar
integer          nswrgp, imligp, iwarnp
integer          ipcvis, ipcvst, ipcvsv
integer          iphydp
double precision epsrgp, climgp, extrap
double precision vecfac, visttt

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

iuiph  = iu(iphas)
iviph  = iv(iphas)
iwiph  = iw(iphas)

ipcvis = ipproc(iviscl(iphas))
ipcvst = ipproc(ivisct(iphas))
if(iviscv(iphas).gt.0) then
  ipcvsv = ipproc(iviscv(iphas))
else
  ipcvsv = 0
endif


! --- Calcul de la viscosite totale

if (itytur(iphas).eq.3 ) then
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

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(vistot)
  !==========
  if (ipcvsv.gt.0) then
    call synsca(propce(1,ipcvsv))
    !==========
  endif
endif


!===============================================================================
! 2.  CALCUL DES TERMES DE LA DIVERGENCE
!===============================================================================

! --- Boucle sur les composantes de vitesse Ui

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
   nphas  ,                                                       &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   w6     , w6     , w6     ,                                     &
   rtp(1,ivar)     , coefa(1,iclvar) , coefb(1,iclvar) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   ra     )



! --- Assemblage sur les faces internes

! On a echange le gradient dans grdcel et vistot plus haut

  if(ipcvsv.gt.0) then
    if    (isou.eq.1) then
      do iel = 1, ncelet
        visttt = propce(iel,ipcvsv) - 2.d0/3.d0*vistot(iel)
        w4(iel) = vistot(iel)*( 2.d0*w1(iel)*ux(iel)              &
                                   + w2(iel)*uy(iel)              &
                                   + w3(iel)*uz(iel) )            &
                     + visttt*w1(iel)*ux(iel)
        w5(iel) = vistot(iel)*w2(iel)*ux(iel)                     &
                     + visttt*w1(iel)*uy(iel)
        w6(iel) = vistot(iel)*w3(iel)*ux(iel)                     &
                     + visttt*w1(iel)*uz(iel)
      enddo

    elseif(isou.eq.2) then
      do iel = 1, ncelet
        visttt = propce(iel,ipcvsv) - 2.d0/3.d0*vistot(iel)
        w4(iel) = vistot(iel)*w1(iel)*uy(iel)                     &
                     + visttt*w2(iel)*ux(iel)
        w5(iel) = vistot(iel)*(      w1(iel)*ux(iel)              &
                              + 2.d0*w2(iel)*uy(iel)              &
                                   + w3(iel)*uz(iel) )            &
                     + visttt*w2(iel)*uy(iel)
        w6(iel) = vistot(iel)*w3(iel)*uy(iel)                     &
                     + visttt*w2(iel)*uz(iel)
      enddo

    elseif(isou.eq.3) then
      do iel = 1, ncelet
        visttt = propce(iel,ipcvsv) - 2.d0/3.d0*vistot(iel)
        w4(iel) = vistot(iel)*w1(iel)*uz(iel)                     &
                     + visttt*w3(iel)*ux(iel)
        w5(iel) = vistot(iel)*w2(iel)*uz(iel)                     &
                     + visttt*w3(iel)*uy(iel)
        w6(iel) = vistot(iel)*(      w1(iel)*ux(iel)              &
                                   + w2(iel)*uy(iel)              &
                              + 2.d0*w3(iel)*uz(iel) )            &
                     + visttt*w3(iel)*uz(iel)
      enddo

    endif

  else

    if    (isou.eq.1) then
      do iel = 1, ncelet
        visttt = viscv0(iphas) - 2.d0/3.d0*vistot(iel)
        w4(iel) = vistot(iel)*( 2.d0*w1(iel)*ux(iel)              &
                                   + w2(iel)*uy(iel)              &
                                   + w3(iel)*uz(iel) )            &
                     + visttt*w1(iel)*ux(iel)
        w5(iel) = vistot(iel)*w2(iel)*ux(iel)                     &
                     + visttt*w1(iel)*uy(iel)
        w6(iel) = vistot(iel)*w3(iel)*ux(iel)                     &
                     + visttt*w1(iel)*uz(iel)
      enddo

    elseif(isou.eq.2) then
      do iel = 1, ncelet
        visttt = viscv0(iphas) - 2.d0/3.d0*vistot(iel)
        w4(iel) = vistot(iel)*w1(iel)*uy(iel)                     &
                     + visttt*w2(iel)*ux(iel)
        w5(iel) = vistot(iel)*(      w1(iel)*ux(iel)              &
                              + 2.d0*w2(iel)*uy(iel)              &
                                   + w3(iel)*uz(iel) )            &
                     + visttt*w2(iel)*uy(iel)
        w6(iel) = vistot(iel)*w3(iel)*uy(iel)                     &
                     + visttt*w2(iel)*uz(iel)
      enddo

    elseif(isou.eq.3) then
      do iel = 1, ncelet
        visttt = viscv0(iphas) - 2.d0/3.d0*vistot(iel)
        w4(iel) = vistot(iel)*w1(iel)*uz(iel)                     &
                     + visttt*w3(iel)*ux(iel)
        w5(iel) = vistot(iel)*w2(iel)*uz(iel)                     &
                     + visttt*w3(iel)*uy(iel)
        w6(iel) = vistot(iel)*(      w1(iel)*ux(iel)              &
                                   + w2(iel)*uy(iel)              &
                              + 2.d0*w3(iel)*uz(iel) )            &
                     + visttt*w3(iel)*uz(iel)
      enddo

    endif

  endif



! On initialise DIVERG(NCEL+1, NCELET)
!     (valeur bidon, mais pas NaN : les calculs sur le halo sont
!      par principe denue de sens, sauf exception)
  if(ncelet.gt.ncel) then
    do iel = ncel+1, ncelet
      diverg(iel) = 0.d0
    enddo
  endif



  if(ivecti.eq.1) then

!CDIR NODEP
    do ifac = 1, nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
!MO             VECFAC = SURFAC(ISOU,IFAC)
!MO     &                  *(POND(IFAC)*W4(II)+(1.D0-POND(IFAC))*W4(JJ))
      vecfac = surfac(1,ifac)*(w4(ii)+w4(jj))*0.5d0               &
             + surfac(2,ifac)*(w5(ii)+w5(jj))*0.5d0               &
             + surfac(3,ifac)*(w6(ii)+w6(jj))*0.5d0
      diverg(ii) = diverg(ii) + vecfac
      diverg(jj) = diverg(jj) - vecfac
    enddo

  else

! VECTORISATION NON FORCEE
    do ifac = 1, nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
!MO             VECFAC = SURFAC(ISOU,IFAC)
!MO     &                  *(POND(IFAC)*W4(II)+(1.D0-POND(IFAC))*W4(JJ))
      vecfac = surfac(1,ifac)*(w4(ii)+w4(jj))*0.5d0               &
             + surfac(2,ifac)*(w5(ii)+w5(jj))*0.5d0               &
             + surfac(3,ifac)*(w6(ii)+w6(jj))*0.5d0
      diverg(ii) = diverg(ii) + vecfac
      diverg(jj) = diverg(jj) - vecfac
    enddo

  endif


! --- Assemblage sur les faces de bord

  if(ivectb.eq.1) then

!CDIR NODEP
    do ifac = 1, nfabor
      ii = ifabor(ifac)
      vecfac = surfbo(1,ifac)*w4(ii)                              &
             + surfbo(2,ifac)*w5(ii)                              &
             + surfbo(3,ifac)*w6(ii)
      diverg(ii) = diverg(ii) + vecfac
    enddo

  else

! VECTORISATION NON FORCEE
    do ifac = 1, nfabor
      ii = ifabor(ifac)
      vecfac = surfbo(1,ifac)*w4(ii)                              &
             + surfbo(2,ifac)*w5(ii)                              &
             + surfbo(3,ifac)*w6(ii)
      diverg(ii) = diverg(ii) + vecfac
    enddo

  endif

enddo


return

end subroutine
