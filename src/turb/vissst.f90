!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine vissst &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel )

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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
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
use pointe, only: s2kw, divukw, ifapat, dispar, coefau, coefbu
use numvar
use optcal
use cstphy
use entsor
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)

! Local variables

integer          iel, iccocg, inc
integer          ipcliu
integer          ipcrom, ipcvis, ipcvst
integer          nswrgp, imligp, iwarnp
integer          ifacpt

double precision d1s3, d2s3
double precision epsrgp, climgp, extrap
double precision xk, xw, rom, xmu, xdist, xarg2, xf2

logical          ilved

double precision, allocatable, dimension(:) :: w1
double precision, dimension(:,:,:), allocatable :: gradv

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! --- Rang des variables dans PROPCE (prop. physiques au centre)
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)
ipcrom = ipproc(irom  )

! --- Rang des c.l. des variables dans COEFA COEFB
!        (c.l. std, i.e. non flux)
ipcliu = iclrtp(iu,icoef)

d1s3 = 1.d0/3.d0
d2s3 = 2.d0/3.d0

!===============================================================================
! 2. Compute the scalar s2kw rate SijSij and the trace of the velocity
!    gradient

!      (Sij^D) (Sij^D)  is stored in    s2kw (deviatoric s2kw tensor rate)
!      tr(Grad u)       is stored in    divukw
!===============================================================================


! Allocate temporary arrays for gradients calculation
allocate(gradv(ncelet,3,3))

iccocg = 1
inc = 1

nswrgp = nswrgr(iu)
imligp = imligr(iu)
iwarnp = iwarni(iu)
epsrgp = epsrgr(iu)
climgp = climgr(iu)
extrap = extrag(iu)

ilved = .false.

call grdvec &
!==========
( iu  , imrgra , inc    ,                               &
  nswrgr(iu) , imligr(iu) , iwarni(iu) ,                &
  nfecra , epsrgr(iu) , climgr(iu) , extrag(iu) ,       &
  ilved  ,                                              &
  rtpa(1,iu) ,  coefau , coefbu,                        &
  gradv  )

! s2kw = Stain rate of the deviatoric part of the s2kw tensor
!      = 2 (Sij^D).(Sij^D)
! divukw   = trace of the velocity gradient
!          = dudx + dvdy + dwdz

do iel = 1, ncel

  s2kw(iel) = 2.d0                                                           &
    *( ( d2s3*gradv(iel,1,1) - d1s3*gradv(iel,2,2) - d1s3*gradv(iel,3,3))**2   &
     + (-d1s3*gradv(iel,1,1) + d2s3*gradv(iel,2,2) - d1s3*gradv(iel,3,3))**2   &
     + (-d1s3*gradv(iel,1,1) - d1s3*gradv(iel,2,2) + d2s3*gradv(iel,3,3))**2   &
     )                                                                         &
    + (gradv(iel,2,1) + gradv(iel,1,2))**2                                     &
    + (gradv(iel,3,1) + gradv(iel,1,3))**2                                     &
    + (gradv(iel,3,2) + gradv(iel,2,3))**2

  divukw(iel) = gradv(iel,1,1) + gradv(iel,2,2) + gradv(iel,3,3)

enddo

! Free memory
deallocate(gradv)

!===============================================================================
! 3.  CALCUL DE LA DISTANCE A LA PAROI
!===============================================================================

! Allocate a work array
allocate(w1(ncelet))

if(abs(icdpar).eq.2) then
  do iel = 1 , ncel
    ifacpt = ifapat(iel)
    if (ifacpt.gt.0) then
      w1(iel) =  (cdgfbo(1,ifacpt)-xyzcen(1,iel))**2           &
               + (cdgfbo(2,ifacpt)-xyzcen(2,iel))**2           &
               + (cdgfbo(3,ifacpt)-xyzcen(3,iel))**2
      w1(iel) = sqrt(w1(iel))
    else
      w1(iel) = grand
    endif
  enddo
else
  do iel = 1 , ncel
    w1(iel) =  max(dispar(iel),epzero)
  enddo
endif

!===============================================================================
! 4.  CALCUL DE LA VISCOSITE
!===============================================================================

do iel = 1, ncel

  xk = rtpa(iel,ik)
  xw = rtpa(iel,iomg)
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
deallocate(w1)

!-------
! Format
!-------

!----
! End
!----

return
end subroutine
