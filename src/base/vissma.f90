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

subroutine vissma &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   ra     )

!===============================================================================
! FONCTION :
! --------

! CALCUL DE LA VISCOSITE "TURBULENTE" POUR
!          UN MODELE LES SMAGORINSKI

! PROPCE(1,IVISCT) = ROM * (SMAGO  * L) **2 * SQRT ( 2 * Sij.Sij )
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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! ia(*)            ! ia ! --- ! main integer work array                        !
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
use dimens, only: ndimfb
use numvar
use optcal
use cstphy
use entsor
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          iel, iccocg, inc
integer          ipcliu, ipcliv, ipcliw
integer          ipcrom, ipcvst

double precision coef, deux, delta
double precision s11, s22, s33
double precision dudy, dudz, dvdx, dvdz, dwdx, dwdy
double precision xfil, xa  , xb  , radeux

double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w4, w5

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(w4(ncelet), w5(ncelet))

! --- Memoire
idebia = idbia0
idebra = idbra0

! --- Rang des variables dans PROPCE (prop. physiques au centre)
ipcvst = ipproc(ivisct)
ipcrom = ipproc(irom  )

! --- Rang des c.l. des variables dans COEFA COEFB
!        (c.l. std, i.e. non flux)
ipcliu = iclrtp(iu,icoef)
ipcliv = iclrtp(iv,icoef)
ipcliw = iclrtp(iw,icoef)
! --- Pour le calcul de la viscosite de sous-maille
xfil   = xlesfl
xa     = ales
xb     = bles
deux   = 2.d0
radeux = sqrt(deux)

!===============================================================================
! 2.  CALCUL DES GRADIENTS DE VITESSE ET DE
!       S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2)
!===============================================================================

iccocg = 1
inc = 1

! W1 = DUDX, W2 = DUDY, W3=DUDZ

call grdcel &
!==========
 ( iu  , imrgra , inc    , iccocg ,                      &
   nswrgr(iu) , imligr(iu) , iwarni(iu) ,                &
   nfecra , epsrgr(iu) , climgr(iu) , extrag(iu) ,       &
   ia     ,                                              &
   rtpa(1,iu) , coefa(1,ipcliu) , coefb(1,ipcliu) ,      &
   w1     , w2     , w3     ,                            &
!        ------   ------   ------
   ra     )

do iel = 1, ncel
  s11  = w1(iel)
  propce(iel,ipcvst) = s11**2
enddo


!            W2 = DUDY, W3=DUDZ
! W4 = DVDX, W1 = DVDY, W5=DVDZ

call grdcel &
!==========
 ( iv  , imrgra , inc    , iccocg ,                      &
   nswrgr(iv) , imligr(iv) , iwarni(iv) ,                &
   nfecra , epsrgr(iv) , climgr(iv) , extrag(iv) ,       &
   ia     ,                                              &
   rtpa(1,iv) , coefa(1,ipcliv) , coefb(1,ipcliv) ,      &
   w4     , w1     , w5     ,                            &
!        ------   ------   ------
   ra     )

do iel = 1, ncel
  s22 = w1(iel)
  propce(iel,ipcvst) = propce(iel,ipcvst) + s22**2
enddo
do iel = 1, ncel
  dudy = w2(iel)
  dvdx = w4(iel)
  propce(iel,ipcvst) = propce(iel,ipcvst) + 0.5d0*(dudy+dvdx)**2
enddo


!                       W3=DUDZ
!            W1 = DVDY, W5=DVDZ
! W2 = DWDX, W4 = DWDY, W1=DWDZ

call grdcel &
!==========
 ( iw  , imrgra , inc    , iccocg ,                      &
   nswrgr(iw) , imligr(iw) , iwarni(iw) ,                &
   nfecra , epsrgr(iw) , climgr(iw) , extrag(iw) ,       &
   ia     ,                                              &
   rtpa(1,iw) , coefa(1,ipcliw) , coefb(1,ipcliw) ,      &
   w2     , w4     , w1     ,                            &
!        ------   ------   ------
   ra     )

do iel = 1, ncel
  s33 = w1(iel)
  propce(iel,ipcvst) = propce(iel,ipcvst) + s33**2
enddo
do iel = 1, ncel
  dudz = w3(iel)
  dwdx = w2(iel)
  dvdz = w5(iel)
  dwdy = w4(iel)
  propce(iel,ipcvst) =                                            &
    propce(iel,ipcvst) + 0.5d0*((dudz+dwdx)**2+(dvdz+dwdy)**2)
enddo


!===============================================================================
! 3.  CALCUL DE LA VISCOSITE (DYNAMIQUE)
!===============================================================================

coef = csmago**2 * radeux

do iel = 1, ncel
  delta  = xfil* (xa*volume(iel))**xb
  delta  = coef * delta**2
  propce(iel,ipcvst) =                                            &
    propce(iel,ipcrom) * delta * sqrt(propce(iel,ipcvst))
enddo

! Free memory
deallocate(w1, w2, w3)
deallocate(w4, w5)

!----
! FORMAT
!----


!----
! FIN
!----

return
end subroutine
