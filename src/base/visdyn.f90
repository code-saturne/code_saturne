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

subroutine visdyn &
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
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smagor ,                                                       &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     , w9     , w10    , xmij   , &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! CALCUL DE LA VISCOSITE "TURBULENTE" POUR
! UN MODELE LES SMAGORINSKI DYNAMIQUE

! SMAGO = LijMij/MijMij

! PROPCE(1,IVISCT(IPHAS)) = ROM * SMAGO  * L**2 * SQRT ( 2 * Sij.Sij )
!       Sij = (DUi/Dxj + DUj/Dxi)/2

! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)

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
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
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
! smagor(ncelet    ! tr ! <-- ! constante de smagorinsky dans le cas           !
! , nphas)         !    !     ! d'un modlele dynamique                         !
! w1..10(ncelet    ! tr ! --- ! tableau de travail                             !
! xmij(ncelet,6    ! tr ! --- ! tableau de travail                             !
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
include "numvar.h"
include "cstnum.h"
include "optcal.h"
include "cstphy.h"
include "entsor.h"
include "period.h"
include "parall.h"

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
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision smagor(ncelet)
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision w5(ncelet),w6(ncelet),w7(ncelet),w8(ncelet)
double precision w9(ncelet),w10(ncelet),xmij(ncelet,6)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          ii, iel, iccocg, inc
integer          iuiph, iviph, iwiph
integer          ipcliu, ipcliv, ipcliw
integer          ipcrom, ipcvst, iphydp
integer          iclipc
double precision coef, radeux, deux, delta, deltaf
double precision s11, s22, s33, s11f, s22f, s33f
double precision dudy, dudz, dvdx, dvdz, dwdx, dwdy
double precision dudyf, dudzf, dvdxf, dvdzf, dwdxf, dwdyf
double precision xfil, xa, xb, xfil2, xsmgmx
double precision aij, bij
double precision xl11, xl22, xl33, xl12, xl13, xl23
double precision xm11, xm22, xm33, xm12, xm13, xm23
double precision smagma, smagmn, smagmy

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! --- Memoire
idebia = idbia0
idebra = idbra0

! --- Numero des variables (dans RTP)
iuiph = iu(iphas)
iviph = iv(iphas)
iwiph = iw(iphas)

! --- Rang des variables dans PROPCE (prop. physiques au centre)
ipcvst = ipproc(ivisct(iphas))
ipcrom = ipproc(irom  (iphas))

! --- Rang des c.l. des variables dans COEFA COEFB
!        (c.l. std, i.e. non flux)
ipcliu = iclrtp(iuiph,icoef)
ipcliv = iclrtp(iviph,icoef)
ipcliw = iclrtp(iwiph,icoef)

! --- Pour le calcul de la viscosite de sous-maille
xfil   = xlesfl(iphas)
xfil2  = xlesfd(iphas)
xa     = ales(iphas)
xb     = bles(iphas)
deux   = 2.d0
radeux = sqrt(deux)
xsmgmx = smagmx(iphas)

!===============================================================================
! 2.  CALCUL DES GRADIENTS DE VITESSE ET DE
!       S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2)
!===============================================================================

!     Les RTPA ont ete echange pour les calculs en parallele,
!       au debut du pas de temps (donc pas utile de le refaire ici)

iccocg = 1
inc = 1
iphydp = 0

! W1 = DUDX, W2 = DUDY, W3=DUDZ

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iuiph  , imrgra , inc    , iccocg ,                            &
   nswrgr(iuiph) , imligr(iuiph) , iphydp , iwarni(iuiph) ,       &
   nfecra , epsrgr(iuiph) , climgr(iuiph) , extrag(iuiph) ,       &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w6     , w6     , w6     ,                                     &
   rtpa(1,iuiph) , coefa(1,ipcliu) , coefb(1,ipcliu) ,            &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w6     , w7     , w8     ,                                     &
   rdevel , rtuser , ra     )

! Filtrage de W1, LE RESULTAT EST DANS W6

call cfiltr                                                       &
!==========
 ( w1     , w6     , w7     , w8     )

do iel = 1, ncel
  s11   = w1(iel)
  s11f  = w6(iel)
  xmij(iel,1)        = s11
  propce(iel,ipcvst) = s11**2
  w9(iel) = s11f**2
enddo


!            W2 = DUDY, W3=DUDZ
! W4 = DVDX, W1 = DVDY, W5=DVDZ

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iviph  , imrgra , inc    , iccocg ,                            &
   nswrgr(iviph) , imligr(iviph) , iphydp , iwarni(iviph) ,       &
   nfecra , epsrgr(iviph) , climgr(iviph) , extrag(iviph) ,       &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w6     , w6     , w6     ,                                     &
   rtpa(1,iviph) , coefa(1,ipcliv) , coefb(1,ipcliv) ,            &
   w4     , w1     , w5     ,                                     &
!        ------   ------   ------
   w6     , w7     , w8     ,                                     &
   rdevel , rtuser , ra     )

call cfiltr                                                       &
!==========
 ( w1     , w6     , w7     , w8    )

do iel = 1, ncel
  s22  = w1(iel)
  s22f = w6(iel)
  xmij(iel,2) = s22
  propce(iel,ipcvst) = propce(iel,ipcvst) + s22**2
  w9(iel) = w9(iel) + s22f**2
enddo

call cfiltr                                                       &
!==========
 ( w2     , w6     , w8     , w1     )

call cfiltr                                                       &
!==========
 ( w4     , w7     , w8     , w1     )

do iel = 1, ncel
  dudy  = w2(iel)
  dvdx  = w4(iel)
  dudyf = w6(iel)
  dvdxf = w7(iel)
  xmij(iel,4) = 0.5d0*(dudy+dvdx)
  propce(iel,ipcvst) = propce(iel,ipcvst) + 0.5d0*(dudy+dvdx)**2
  w9(iel) = w9(iel) + 0.5d0*(dudyf+dvdxf)**2
enddo

!                       W3=DUDZ
!                       W5=DVDZ
! W2 = DWDX, W4 = DWDY, W1=DWDZ

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iwiph  , imrgra , inc    , iccocg ,                            &
   nswrgr(iwiph) , imligr(iwiph) , iphydp , iwarni(iwiph) ,       &
   nfecra , epsrgr(iwiph) , climgr(iwiph) , extrag(iwiph) ,       &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w6     , w6     , w6     ,                                     &
   rtpa(1,iwiph) , coefa(1,ipcliw) , coefb(1,ipcliw) ,            &
   w2     , w4     , w1     ,                                     &
!        ------   ------   ------
   w6     , w7     , w8     ,                                     &
   rdevel , rtuser , ra     )

call cfiltr                                                       &
!==========
 ( w1     , w6     , w7     , w8     )

do iel = 1, ncel
  s33  = w1(iel)
  s33f = w6(iel)
  xmij(iel,3) = s33
  propce(iel,ipcvst) = propce(iel,ipcvst) + s33**2
  w9(iel) = w9(iel) + s33f**2
enddo

call cfiltr                                                       &
!==========
 ( w2     , w1     , w7     , w8     )

call cfiltr                                                       &
!==========
 ( w3     , w6     , w7     , w8     )

do iel = 1, ncel
  dudz  = w3(iel)
  dwdx  = w2(iel)
  dudzf = w6(iel)
  dwdxf = w1(iel)
  xmij(iel,5) = 0.5d0*(dudz+dwdx)
  propce(iel,ipcvst) =                                            &
    propce(iel,ipcvst) + 0.5d0*(dudz+dwdx)**2
  w9(iel) = w9(iel) + 0.5d0*(dudzf+dwdxf)**2
enddo

call cfiltr                                                       &
!==========
 ( w4     , w1     , w7     , w8     )

call cfiltr                                                       &
!==========
 ( w5     , w6     , w7     , w8     )

do iel = 1, ncel
  dvdz  = w5(iel)
  dwdy  = w4(iel)
  dvdzf = w6(iel)
  dwdyf = w1(iel)
  xmij(iel,6) = 0.5d0*(dvdz+dwdy)
  propce(iel,ipcvst) =                                            &
    propce(iel,ipcvst) + 0.5d0*(dvdz+dwdy)**2
  w9(iel) = w9(iel) + 0.5d0*(dvdzf+dwdyf)**2
enddo

do iel = 1, ncel
  propce(iel,ipcvst) = radeux*sqrt(propce(iel,ipcvst))
  w9(iel)            = radeux*sqrt(w9(iel)           )
enddo

!     Ici XMIJ contient Sij
!         PROPCE(IEL,IPCVST) contient ||S||
!            SQRT(2)*SQRT(S11^2+S22^2+S33^2+2(S12^2+S13^2+S23^2))
!         W9                 contient ||SF||
!            SQRT(2)*SQRT(S11F^2+S22F^2+S33F^2+2(S12F^2+S13F^2+S23F^2))

!===============================================================================
! 3.  CALCUL DE Mij
!===============================================================================

do iel = 1, ncel
  w7(iel) = xfil *(xa*volume(iel))**xb
enddo

do ii = 1, 6

  call cfiltr                                                     &
  !==========
 ( xmij(1,ii) , w1     , w2     , w3     )

  do iel = 1, ncel
    delta = w7(iel)
    w2(iel) = -deux*delta**2*propce(iel,ipcvst)*xmij(iel,ii)
  enddo

  call cfiltr                                                     &
  !==========
 ( w2     , w3     , w4     , w5     )

  do iel = 1, ncel
    delta = w7(iel)
    deltaf = xfil2*delta
    aij    = -deux*deltaf**2*w9(iel)*w1(iel)
    bij    =  w3(iel)
    xmij(iel,ii) = aij - bij
  enddo

enddo

!     Ici Aij contient alpha_ij, Bij contient beta_ij tilde
!        et XMIJ contient M_ij

!===============================================================================
! 4.  CALCUL DE LA CONSTANTE DE SMAGORINSKY DYNAMIQUE
!===============================================================================

! FILTRAGE DE LA VITESSE ET DE SON CARRE


! U**2
do iel = 1,ncel
  w9(iel) = rtp(iel,iuiph)*rtp(iel,iuiph)
enddo
call cfiltr                                                       &
!==========
 ( w9     , w1     , w7     , w8     )

! V**2
do iel = 1,ncel
  w9(iel) = rtp(iel,iviph)*rtp(iel,iviph)
enddo
call cfiltr                                                       &
!==========
 ( w9     , w2     , w7     , w8     )

! W**2
do iel = 1,ncel
  w9(iel) = rtp(iel,iwiph)*rtp(iel,iwiph)
enddo
call cfiltr                                                       &
!==========
 ( w9     , w3     , w7     , w8     )

! UV
do iel = 1,ncel
  w9(iel) = rtp(iel,iuiph)*rtp(iel,iviph)
enddo
call cfiltr                                                       &
!==========
 ( w9     , w4     , w7     , w8     )

! UW
do iel = 1,ncel
  w9(iel) = rtp(iel,iuiph)*rtp(iel,iwiph)
enddo
call cfiltr                                                       &
!==========
 ( w9     , w5     , w7     , w8     )

! VW
do iel = 1,ncel
  w9(iel) = rtp(iel,iviph)*rtp(iel,iwiph)
enddo
call cfiltr                                                       &
!==========
 ( w9     , w6     , w7     , w8     )

! U
call cfiltr                                                       &
!==========
 ( rtp(1,iuiph)    , w7     , w8     , w9     )

! V
call cfiltr                                                       &
!==========
 ( rtp(1,iviph)    , w8     , w9     , smagor )

! W
call cfiltr                                                       &
!==========
 ( rtp(1,iwiph)    , w9     , smagor , w10    )

do iel = 1, ncel

! --- Calcul de Lij
  xl11 = w1(iel) - w7(iel) * w7(iel)
  xl22 = w2(iel) - w8(iel) * w8(iel)
  xl33 = w3(iel) - w9(iel) * w9(iel)
  xl12 = w4(iel) - w7(iel) * w8(iel)
  xl13 = w5(iel) - w7(iel) * w9(iel)
  xl23 = w6(iel) - w8(iel) * w9(iel)

  xm11 = xmij(iel,1)
  xm22 = xmij(iel,2)
  xm33 = xmij(iel,3)
  xm12 = xmij(iel,4)
  xm13 = xmij(iel,5)
  xm23 = xmij(iel,6)
! ---Calcul de Mij :: Lij
  w1(iel) = xm11 * xl11 + 2.d0* xm12 * xl12 + 2.d0* xm13 * xl13 + &
                                xm22 * xl22 + 2.d0* xm23 * xl23 + &
                                                    xm33 * xl33
! ---Calcul de Mij :: Mij
  w2(iel) = xm11 * xm11 + 2.d0* xm12 * xm12 + 2.d0* xm13 * xm13 + &
                                xm22 * xm22 + 2.d0* xm23 * xm23 + &
                                                    xm33 * xm33

enddo

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(w1)
  !==========
  call synsca(w2)
  !==========
endif

!     Par defaut on fait une moyenne locale du numerateur et du
!     denominateur, puis seulement on fait le rapport.
!     L'utilisateur peut faire autrement dans USSMAG

call cfiltr                                                       &
!==========
 ( w1     , w3     , w5     , w6     )

call cfiltr                                                       &
!==========
 ( w2     , w4     , w5     , w6     )

do iel = 1, ncel
  if(abs(w4(iel)).le.epzero) then
    smagor(iel) = xsmgmx**2
  else
    smagor(iel) = w3(iel)/w4(iel)
  endif
enddo

call ussmag                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iphas  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , icepdc , icetsm , itypsm , &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smagor , w1     , w2     ,                                     &
   w3     , w4     , w5     , w6     ,                            &
   rdevel , rtuser , ra     )

iclipc = 0
do iel = 1, ncel
  if(smagor(iel).ge.xsmgmx**2) then
    smagor(iel) = xsmgmx**2
    iclipc = iclipc + 1
  elseif(smagor(iel).le.-xsmgmx**2) then
    smagor(iel) = -xsmgmx**2
    iclipc = iclipc + 1
  endif
enddo

!===============================================================================
! 3.  CALCUL DE LA VISCOSITE (DYNAMIQUE)
!===============================================================================

! On clippe en (mu + mu_t)>0 dans phyvar

do iel = 1, ncel
  coef = smagor(iel)
  delta  = xfil * (xa*volume(iel))**xb
  propce(iel,ipcvst) = propce(iel,ipcrom)                         &
       * coef * delta**2 * propce(iel,ipcvst)
enddo

!     Quelques impressions
if(iwarni(iuiph).ge.1) then

  smagma = -1.0d12
  smagmn =  1.0d12
  smagmy =  0.d0
  do iel = 1, ncel
    smagma = max(smagma,smagor(iel))
    smagmn = min(smagmn,smagor(iel))
    smagmy = smagmy + smagor(iel)*volume(iel)
  enddo
  if(irangp.ge.0) then
    call parmax(smagma)
    !==========
    call parmin(smagmn)
    !==========
    call parsom(smagmy)
    !==========
    call parcpt(iclipc)
    !==========
  endif
  smagmy = smagmy / voltot
  write(nfecra,1000) iclipc
  write(nfecra,2001) iphas
  write(nfecra,2002) smagma, smagmn, smagmy
  write(nfecra,2003)

endif

!----
! FORMAT
!----

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
' Nb Clipping Constante Smagorinsky par valeurs maximales ',I10,/)
 2001 format(                                                           &
' --- Phase : ',I10                                            ,/,&
' --- Informations sur la constante de Smagorinsky^2          ',/,&
' ----------------------------------                          ',/,&
' Valeur moy  Valeur min  Valeur max                          ',/,&
' ----------------------------------                          '  )
 2002 format(                                                           &
 e12.4    ,      e12.4,      e12.4                               )
 2003 format(                                                           &
' ----------------------------------                          ',/)

#else

 1000 format(                                                           &
' Nb of clipping of the Smagorinsky constant by max values',I10,/)
 2001 format(                                                           &
' --- Phase: ',I10                                             ,/,&
' --- Informations on the squared Smagorinsky constant'        ,/,&
' --------------------------------'                            ,/,&
' Mean value  Min value  Max value'                            ,/,&
' --------------------------------'                              )
 2002 format(                                                           &
 e12.4    ,      e12.4,      e12.4                               )
 2003 format(                                                           &
' --------------------------------'                            ,/)

#endif

!----
! FIN
!----

return
end subroutine
