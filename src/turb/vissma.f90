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

subroutine vissma &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel )

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
use numvar
use optcal
use cstphy
use entsor
use pointe, only: coefau, coefbu
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
integer          ipcliu, ipcliv, ipcliw
integer          ipcrom, ipcvst

double precision coef, deux, delta
double precision s11, s22, s33
double precision dudy, dudz, dvdx, dvdz, dwdx, dwdy
double precision xfil, xa  , xb  , radeux

logical          ilved

double precision, dimension(:,:,:), allocatable :: gradv

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! Allocate temporary arrays for gradients calculation
allocate(gradv(ncelet,3,3))

! --- Memoire

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

ilved = .false.

call grdvec &
!==========
( iu  , imrgra , inc    ,                               &
  nswrgr(iu) , imligr(iu) , iwarni(iu) ,                &
  nfecra , epsrgr(iu) , climgr(iu) , extrag(iu) ,       &
  ilved  ,                                              &
  rtpa(1,iu) ,  coefau , coefbu,                        &
  gradv  )

do iel = 1, ncel

  ! gradv(iel, xyz, uvw)
  s11  = gradv(iel,1,1)
  s22  = gradv(iel,2,2)
  s33  = gradv(iel,3,3)
  dudy = gradv(iel,2,1)
  dvdx = gradv(iel,1,2)
  dudz = gradv(iel,3,1)
  dwdx = gradv(iel,1,3)
  dvdz = gradv(iel,3,2)
  dwdy = gradv(iel,2,3)

  propce(iel,ipcvst) = s11**2 + s22**2 + s33**2       &
                     + 0.5d0*((dudy+dvdx)**2          &
                     +        (dudz+dwdx)**2          &
                     +        (dvdz+dwdy)**2)
enddo

! Free memory
deallocate(gradv)

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

!----
! FORMAT
!----


!----
! FIN
!----

return
end subroutine
