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

subroutine viswal &
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
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   w9     , w10    , w11    , w12    ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! CALCUL DE LA VISCOSITE "TURBULENTE" POUR
!          UN MODELE LES WALE

! PROPCE(1,IVISCT(IPHAS)) = ROM * (CWALE*L)**2 *
!   [(Sijd.Sijd)**(3/2)] / [(Sij.Sij)**(5/2) + (Sijd.Sijd)**(5/4)]
!
! avec
!       Sij = 0.5*[DUi/Dxj + DUj/Dxi]
! et
!       Sijd = 0.5*[DUi/Dxk.DUk/Dxj + DUj/Dxk.DUk/Dxi] - 1/3*Delta_ij.Gkk**2

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
! ckupdc(ncepdp    ! tr ! <-- ! tableau de travail pour pdc                    !
!     , 6)!        !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! w1..12(ncelet    ! tr ! --- ! tableau de travail                             !
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
include "optcal.h"
include "cstphy.h"
include "entsor.h"
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
double precision w1(ncelet)   , w2(ncelet)   , w3(ncelet)
double precision w4(ncelet)   , w5(ncelet)   , w6(ncelet)
double precision w7(ncelet)   , w8(ncelet)   , w9(ncelet)
double precision w10(ncelet)  , w11(ncelet)  , w12(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra, ifinra
integer          iel, iccocg, inc
integer          iuiph, iviph, iwiph
integer          ipcliu, ipcliv, ipcliw
integer          ipcrom, ipcvst, iphydp, ipcvis
integer          i, j, k

double precision coef, deux, delta, tiers
double precision sij, sijd, s, sd, sinv
double precision xfil, xa  , xb  , radeux, con
double precision dudx(ndim,ndim), kdelta(ndim,ndim)

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
ipcvis = ipproc(iviscl(iphas))
ipcvst = ipproc(ivisct(iphas))
ipcrom = ipproc(irom  (iphas))

! --- Rang des c.l. des variables dans COEFA COEFB
!        (c.l. std, i.e. non flux)
ipcliu = iclrtp(iuiph,icoef)
ipcliv = iclrtp(iviph,icoef)
ipcliw = iclrtp(iwiph,icoef)
! --- Pour le calcul de la viscosite de sous-maille
xfil   = xlesfl(iphas)
xa     = ales(iphas)
xb     = bles(iphas)
deux   = 2.d0
radeux = sqrt(deux)
tiers  = 1.d0/3.d0


!===============================================================================
! 2.  CALCUL DU GRADIENT DE VITESSE
!       W1 = DU/DX, W2 = DU/DY, W3 = DU/DZ
!       W4 = DV/DX, W5 = DV/DY, W6 = DV/DZ
!       W7 = DW/DX, W8 = DW/DY, W9 = DW/DZ
!===============================================================================

iccocg = 1
inc = 1
iphydp = 0

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
   w12    , w12    , w12    ,                                     &
   rtpa(1,iuiph) , coefa(1,ipcliu) , coefb(1,ipcliu) ,            &
   w1            , w2              , w3              ,            &
!        ------   ------   ------
   w10           , w11             , w12             ,            &
   rdevel , rtuser , ra     )

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
   w12    , w12    , w12    ,                                     &
   rtpa(1,iviph) , coefa(1,ipcliv) , coefb(1,ipcliv) ,            &
   w4            , w5              , w6              ,            &
!        ------   ------   ------
   w10           , w11             , w12             ,            &
   rdevel , rtuser , ra     )

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
   w12    , w12    , w12    ,                                     &
   rtpa(1,iwiph) , coefa(1,ipcliw) , coefb(1,ipcliw) ,            &
   w7            , w8              , w9            ,              &
!        ------   ------   ------
   w10           , w11             , w12             ,            &
   rdevel , rtuser , ra     )

! Kronecker delta Dij

kdelta(1,1) = 1
kdelta(1,2) = 0
kdelta(1,3) = 0
kdelta(2,1) = 0
kdelta(2,2) = 1
kdelta(2,3) = 0
kdelta(3,1) = 0
kdelta(3,2) = 0
kdelta(3,3) = 1

coef = cwale(iphas)**2 * radeux

do iel = 1, ncel

  dudx(1,1) = w1(iel)
  dudx(1,2) = w2(iel)
  dudx(1,3) = w3(iel)
  dudx(2,1) = w4(iel)
  dudx(2,2) = w5(iel)
  dudx(2,3) = w6(iel)
  dudx(3,1) = w7(iel)
  dudx(3,2) = w8(iel)
  dudx(3,3) = w9(iel)

  s  = 0.d0
  sd = 0.d0

  do i = 1, ndim
    do j = 1, ndim

      ! Sij = 0.5 * (dUi/dXj + dUj/dXi)

      sij = 0.5d0*(dudx(i,j)+dudx(j,i))

      s = s + sij**2

      do k = 1, ndim

!  traceless symmetric part of the square of the velocity gradient tensor
!    Sijd = 0.5 * ( dUi/dXk dUk/dXj + dUj/dXk dUk/dXi) - 1/3 Dij dUk/dXk dUk/dXk

        sijd = 0.5d0*(dudx(i,k)*dudx(k,j)+ dudx(j,k)*dudx(k,i)) &
              -tiers*kdelta(i,j)*dudx(k,k)**2

        sd = sd + sijd**2

      enddo
    enddo
  enddo

!===============================================================================
! 3.  CALCUL DE LA VISCOSITE TURBULENTE
!===============================================================================

  ! Turbulent inverse time scale =
  !   (Sijd Sijd)^3/2 / [ (Sij Sij)^5/2 + (Sijd Sijd)^5/4 ]

  sinv = (s**2.5d0 + sd**1.25d0)
  if (sinv.gt.0.d0) then
    con = sd**1.5d0 / sinv
  else
    con = 0.d0
  endif

  delta = xfil* (xa*volume(iel))**xb
  delta = coef * delta**2

  propce(iel,ipcvst) = propce(iel,ipcrom) * delta * con

enddo

!----
! FORMAT
!----

!----
! FIN
!----

return

end subroutine
