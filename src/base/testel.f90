!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine testel &
!================

 ( nvar   ,                 &
   rtp    , coefa  , coefb  )

!===============================================================================
! FONCTION :
! --------

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1,2,3,4,5,6     ! ra ! --- ! work arrays                                    !
!  (ncelet)        !    !     !  (computation of pressure gradient)            !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use optcal
use cstphy
use cstnum
use pointe
use entsor
use albase
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar


double precision rtp(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)

! Local variables

integer          ntpost
integer          ifac  , iel   , ivar
integer          inc   , iccocg
integer          nswrgp, imligp, iwarnp
integer          ipclip, ialold
integer          indwri, indact, ipart, idimt, ientla, ivarpr

double precision ttpost
double precision epsrgp, climgp, extrap
double precision xx, yy, zz
double precision rbid(1)

character*32     namevr

double precision, allocatable, dimension(:,:) :: grad

!===============================================================================

!===============================================================================
! 0.  INITIALISATIONS
!===============================================================================

! Allocate temporary arrays
allocate(grad(ncelet,3))

! On positionne l'indicateur ALE a 1 de maniere a forcer le recalcul
! de la contribution des cellules de bord a chaque appel de GRDCEL
ialold = iale
iale = 1

! Postprocessing should be time-independent
ntpost = -1
ttpost = 0.d0

! Symmetry type:
! value 0 avoids extrapolating the gradient on boundary faces.
do ifac = 1, nfabor
   isympa(ifac) = 0
enddo

!===============================================================================
! 1. FONCTION ANALYTIQUE SIN(X+2Y+3Z)
!===============================================================================

ivar   = ipr
ipclip = iclrtp(ivar,icoef)

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  rtp(iel,ivar) = sin(xx+2.d0*yy+3.d0*zz)
enddo

do ifac = 1, nfabor
  xx = cdgfbo(1,ifac)
  yy = cdgfbo(2,ifac)
  zz = cdgfbo(3,ifac)
  coefa(ifac,ipclip) = sin(xx+2.d0*yy+3.d0*zz)
enddo

do ifac = 1, nfabor
  coefb(ifac,ipclip) = 0.d0
enddo

! On active le writer standard

indwri = -1
indact = 1
call pstact(indwri, indact)
!==========

! Options de sorties des variables (gradient non entrelaces)

ipart = -1
idimt = 3
ientla = 0
ivarpr = 1

!===============================================================================
! 2. CALCUL DU GRADIENT DE LA FONCTION ANALYTIQUE

!    NE PAS CHANGER L'ORDRE DE CALCUL DES GRADIENTS:
!      * IMRGRA = 0
!      * IMRGRA = 1 (voisinage standard)
!      * IMRGRA = 2 (voisinage etendu)
!      * IMRGRA = 4 (voisinage etendu)
!      * IMRGRA = 3 (reduction du voisinage etendu)
!===============================================================================

inc = 1
iccocg = 1
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
iwarnp = iwarni(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)

!  2.1 APPEL A GRDCEL AVEC IMRGRA = 0
!  ==================================

imrgra = 0
imligp = -1

call grdcel &
!==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   grad   )

! On sort le gradient

namevr = 'Grad_RC'
call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
!==========
            ntpost, ttpost, grad, rbid, rbid)

! Calcul de l'erreur absolue

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  grad(iel,1) = grad(iel,1)-     cos(xx+2.d0*yy+3.d0*zz)
  grad(iel,2) = grad(iel,2)-2.d0*cos(xx+2.d0*yy+3.d0*zz)
  grad(iel,3) = grad(iel,3)-3.d0*cos(xx+2.d0*yy+3.d0*zz)
enddo

! On sort l'erreur

namevr = 'Err_Grad_RC'
call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
!==========
            ntpost, ttpost, grad, rbid, rbid)


!  2.2 APPEL A GRDCEL AVEC IMRGRA = 1
!  ==================================

imrgra = 1
imligp = 1

call grdcel &
!==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   grad   )

! On sort le gradient

namevr = 'Grad_LSQ'
call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
!==========
            ntpost, ttpost, grad, rbid, rbid)

! Calcul de l'erreur absolue

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  grad(iel,1) = grad(iel,1)-     cos(xx+2.d0*yy+3.d0*zz)
  grad(iel,2) = grad(iel,2)-2.d0*cos(xx+2.d0*yy+3.d0*zz)
  grad(iel,3) = grad(iel,3)-3.d0*cos(xx+2.d0*yy+3.d0*zz)
enddo

! On sort l'erreur

namevr = 'Err_Grad_LSQ'
call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
!==========
            ntpost, ttpost, grad, rbid, rbid)


!  2.3 APPEL A GRDCEL AVEC IMRGRA = 2
!  ==================================

imrgra = 2
imligp = 1

call grdcel &
!==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   grad   )

! On sort le gradient

namevr = 'Grad_LSQ_Ext'
call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
!==========
            ntpost, ttpost, grad, rbid, rbid)

! Calcul de l'erreur absolue

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  grad(iel,1) = grad(iel,1)-     cos(xx+2.d0*yy+3.d0*zz)
  grad(iel,2) = grad(iel,2)-2.d0*cos(xx+2.d0*yy+3.d0*zz)
  grad(iel,3) = grad(iel,3)-3.d0*cos(xx+2.d0*yy+3.d0*zz)
enddo

! On sort l'erreur

namevr = 'Err_Grad_LSQ_Ext'
call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
!==========
              ntpost, ttpost, grad, rbid, rbid)


!  2.4 APPEL A GRDCEL AVEC IMRGRA = 4
!  ==================================

imrgra = 4
imligp = -1

call grdcel &
!==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   grad   )

! On sort le gradient

namevr = 'Grad_LSQ_RC'
call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
!==========
            ntpost, ttpost, grad, rbid, rbid)

! Calcul de l'erreur absolue

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  grad(iel,1) = grad(iel,1)-     cos(xx+2.d0*yy+3.d0*zz)
  grad(iel,2) = grad(iel,2)-2.d0*cos(xx+2.d0*yy+3.d0*zz)
  grad(iel,3) = grad(iel,3)-3.d0*cos(xx+2.d0*yy+3.d0*zz)
enddo

! On sort l'erreur

namevr = 'Err_Grad_LSQ_RC'
call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
!==========
            ntpost, ttpost, grad, rbid, rbid)


!  2.5 APPEL A GRDCEL AVEC IMRGRA = 3
!  ==================================

! Reduction du voisinage etendu

call redvse(anomax)
!==========

imrgra = 3
imligp = 1

call grdcel &
!==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   grad   )

! On sort le gradient

namevr = 'Grad_LSQ_ExtRed'
call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
!==========
            ntpost, ttpost, grad, rbid, rbid)

! Calcul de l'erreur absolue

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  grad(iel,1) = grad(iel,1)-     cos(xx+2.d0*yy+3.d0*zz)
  grad(iel,2) = grad(iel,2)-2.d0*cos(xx+2.d0*yy+3.d0*zz)
  grad(iel,3) = grad(iel,3)-3.d0*cos(xx+2.d0*yy+3.d0*zz)
enddo

! On sort l'erreur

namevr = 'Err_Grad_LSQ_ExtRed'
call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
!==========
            ntpost, ttpost, grad, rbid, rbid)

! Reset ALE flag to old value
! de la contribution des cellules de bord a chaque appel de GRDCEL
iale = ialold

! Free memory
deallocate(grad)

!----
! FIN
!----

return
end subroutine
