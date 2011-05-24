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
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel , s2kw   , divukw ,          &
   ra     )

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
! s2kw(ncelet)     ! tr ! --> ! 2 sij.sij                                      !
! divukw(ncelet    ! tr ! --> ! divergence du u pour utilisation dans          !
!                  !    !     ! turbkw                                         !
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
use cstnum
use pointe
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
double precision s2kw(ncelet), divukw(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          iel, iccocg, inc
integer          iuiph, iviph, iwiph, ikiph, iomgip
integer          ipcliu, ipcliv, ipcliw
integer          ipcrom, ipcvis, ipcvst
integer          nswrgp, imligp, iwarnp
integer          ifacpt

double precision epsrgp, climgp, extrap
double precision xk, xw, rom, xmu, xdist, xarg2, xf2

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

! --- Numero des variables (dans RTP)
iuiph  = iu
iviph  = iv
iwiph  = iw
ikiph  = ik
iomgip = iomg

! --- Rang des variables dans PROPCE (prop. physiques au centre)
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)
ipcrom = ipproc(irom  )

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

call grdcel                                                       &
!==========
 ( iuiph  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   rtpa(1,iuiph)   , coefa(1,ipcliu) , coefb(1,ipcliu) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   ra     )


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

call grdcel                                                       &
!==========
 ( iviph  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   rtpa(1,iviph)   , coefa(1,ipcliv) , coefb(1,ipcliv) ,          &
   w1     , w4     , w5     ,                                     &
!        ------   ------   ------
   ra     )


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

call grdcel                                                       &
!==========
 ( iwiph  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   rtpa(1,iwiph)   , coefa(1,ipcliw) , coefb(1,ipcliw) ,          &
   w1     , w2     , w4     ,                                     &
!        ------   ------   ------
   ra     )

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
    ifacpt = ia(iifapa-1+iel)
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
