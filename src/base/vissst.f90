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

subroutine vissst &
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
   coefa  , coefb  , ckupdc , smacel , s2kw   , divukw ,          &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! CALCUL DE LA VISCOSITE TURBULENTE POUR
!          LE MODELE K-OMEGA SST

! VISCT = ROM * A1 * K /MAX(A1*W ; SQRT(S2KW)*F2)
! AVEC S2KW =  2 * Sij.Sij
!       Sij = (DUi/Dxj + DUj/Dxi)/2

! ET F2 = TANH(ARG2**2)
! ARG2**2 = MAX(2*SQRT(K)/CMU/W/Y ; 500*NU/W/Y**2)

! DIVU EST CALCULE EN MEME TEMPS QUE S2KW POUR ETRE REUTILISE
! DANS TURBKW

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
! ncepdp           ! e  ! <-- ! nombre de cellules avec pdc                    !
! ncesmp           ! e  ! <-- ! nombre de cellules a source de masse           !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! iphas            ! e  ! <-- ! numero de phase                                !
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
! s2kw(ncelet)     ! tr ! --> ! 2 sij.sij                                      !
! divukw(ncelet    ! tr ! --> ! divergence du u pour utilisation dans          !
!                  !    !     ! turbkw                                         !
! w1...8(ncelet    ! tr ! --- ! tableau de travail                             !
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
include "cstnum.h"
include "pointe.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "entsor.h"

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
double precision s2kw(ncelet), divukw(ncelet)
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision w5(ncelet),w6(ncelet),w7(ncelet),w8(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          iel, iccocg, inc
integer          iuiph, iviph, iwiph, ikiph, iomgip
integer          ipcliu, ipcliv, ipcliw
integer          ipcrom, ipcvis, ipcvst
integer          nswrgp, imligp, iwarnp, iphydp
integer          ifacpt
double precision epsrgp, climgp, extrap
double precision xk, xw, rom, xmu, xdist, xarg2, xf2

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! --- Memoire
idebia = idbia0
idebra = idbra0

! --- Numero des variables (dans RTP)
iuiph  = iu  (iphas)
iviph  = iv  (iphas)
iwiph  = iw  (iphas)
ikiph  = ik  (iphas)
iomgip = iomg(iphas)

! --- Rang des variables dans PROPCE (prop. physiques au centre)
ipcvis = ipproc(iviscl(iphas))
ipcvst = ipproc(ivisct(iphas))
ipcrom = ipproc(irom  (iphas))

! --- Rang des c.l. des variables dans COEFA COEFB
!        (c.l. std, i.e. non flux)
ipcliu = iclrtp(iuiph,icoef)
ipcliv = iclrtp(iviph,icoef)
ipcliw = iclrtp(iwiph,icoef)

!===============================================================================
! 2.  CALCUL DES GRADIENTS DE VITESSE ET DE
!       S2KW = 2* (S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2)
!===============================================================================

iccocg = 1
inc = 1


! SMBRK  = DUDX ,W4 = DUDY ,W5 = DUDZ

nswrgp = nswrgr(iuiph)
imligp = imligr(iuiph)
iwarnp = iwarni(iuiph)
epsrgp = epsrgr(iuiph)
climgp = climgr(iuiph)
extrap = extrag(iuiph)
iphydp = 0

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iuiph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w6     , w6     , w6     ,                                     &
   rtpa(1,iuiph)   , coefa(1,ipcliu) , coefb(1,ipcliu) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w6     , w7     , w8     ,                                     &
   rdevel , rtuser , ra     )


! S2KW    = (S11)**2
! DIVUKW = S11

do iel = 1, ncel
  s2kw   (iel)   = w1(iel)**2
  divukw(iel)   = w1(iel)
enddo


nswrgp = nswrgr(iviph)
imligp = imligr(iviph)
iwarnp = iwarni(iviph)
epsrgp = epsrgr(iviph)
climgp = climgr(iviph)
extrap = extrag(iviph)
iphydp = 0

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iviph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w6     , w6     , w6     ,                                     &
   rtpa(1,iviph)   , coefa(1,ipcliv) , coefb(1,ipcliv) ,          &
   w1     , w4     , w5     ,                                     &
!        ------   ------   ------
   w6     , w7     , w8     ,                                     &
   rdevel , rtuser , ra     )


! S2KW    = 2 (S11)**2 + 2 (S22)**2 + (2 S12)**2
! DIVUKW = S11 + S22

do iel = 1, ncel
  s2kw  (iel)   = 2.d0*(s2kw(iel) + w4(iel)**2)                   &
           +  (w1(iel)+w2(iel))**2
  divukw(iel)   = divukw(iel) + w4(iel)
enddo

nswrgp = nswrgr(iwiph)
imligp = imligr(iwiph)
iwarnp = iwarni(iwiph)
epsrgp = epsrgr(iwiph)
climgp = climgr(iwiph)
extrap = extrag(iwiph)
iphydp = 0

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iwiph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w6     , w6     , w6     ,                                     &
   rtpa(1,iwiph)   , coefa(1,ipcliw) , coefb(1,ipcliw) ,          &
   w1     , w2     , w4     ,                                     &
!        ------   ------   ------
   w6     , w7     , w8     ,                                     &
   rdevel , rtuser , ra     )

! S2KW    =  2 (S11)**2 + 2 (S22)**2 + 2 (S33)**2
!        + (2 S12)**2 + (2 S13)**2 + (2 S23)**2 )
! DIVUKW = S11 + S22 + S33

do iel = 1, ncel
  s2kw   (iel)   = s2kw(iel) + 2.d0*w4(iel)**2                    &
           + (w1(iel)+w3(iel))**2                                 &
           + (w2(iel)+w5(iel))**2
  divukw(iel)   = divukw(iel) + w4(iel)
enddo

!===============================================================================
! 3.  CALCUL DE LA DISTANCE A LA PAROI
!===============================================================================

if(abs(icdpar).eq.2) then
  do iel = 1 , ncel
    ifacpt = ia(iifapa(iphas)-1+iel)
    if (ifacpt.gt.0) then
      w1(iel) =                                                   &
            (cdgfbo(1,ifacpt)-xyzcen(1,iel))**2                   &
           +(cdgfbo(2,ifacpt)-xyzcen(2,iel))**2                   &
           +(cdgfbo(3,ifacpt)-xyzcen(3,iel))**2
      w1(iel) = sqrt(w1(iel))
    else
      w1(iel) = grand
    endif
  enddo
else
  do iel = 1 , ncel
    w1(iel) =  max(ra(idipar+iel-1),epzero)
  enddo
endif

!===============================================================================
! 4.  CALCUL DE LA VISCOSITE
!===============================================================================

do iel = 1, ncel

  xk = rtpa(iel,ikiph)
  xw = rtpa(iel,iomgip)
  rom = propce(iel,ipcrom)
  xmu = propce(iel,ipcvis)
  xdist = w1(iel)
  xarg2 = max (                                                   &
       2.d0*sqrt(xk)/cmu/xw/xdist,                                &
       500.d0*xmu/rom/xw/xdist**2 )
  xf2 = tanh(xarg2**2)

  propce(iel,ipcvst) = rom*ckwa1*xk                               &
       /max( ckwa1*xw , sqrt(s2kw(iel))*xf2 )

enddo


!----
! FORMAT
!----


!----
! FIN
!----

return
end subroutine
