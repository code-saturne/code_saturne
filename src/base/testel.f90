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

subroutine testel &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   ,                                                       &
   ia     ,                                                       &
   rtp    , coefa  , coefb  ,                                     &
   ra     )

!===============================================================================
! FONCTION :
! --------

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ia(*)            ! ia ! --- ! main integer work array                        !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1,2,3,4,5,6     ! ra ! --- ! work arrays                                    !
!  (ncelet)        !    !     !  (computation of pressure gradient)            !
! ra(*)            ! ra ! --- ! main real work array                           !
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

integer          idbia0 , idbra0
integer          nvar

integer          ia(*)

double precision rtp(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ra(*)

! Local variables

integer          idebia, idebra, ifinia, ifinra
integer          ifac  , iel   , ivar
integer          inc   , iccocg
integer          iuiph , iviph , iwiph
integer          nswrgp, imligp, iwarnp
integer          ipclip
integer          indwri, indact, ipart, idimt, ientla, ivarpr

double precision epsrgp, climgp, extrap
double precision xx, yy, zz
double precision rbid(1)

character*32     namevr

!===============================================================================

!===============================================================================
! 0.  INITIALISATIONS
!===============================================================================

ifinia = idbia0

! On positionne l'indicateur ALE a 1 de maniere a forcer le recalcul
! de la contribution des cellules de bord a chaque appel de GRDCEL
iale = 1

! Symmetry type:
! value 0 avoids extrapolating the gradient on boundary faces.
do ifac = 1, nfabor
   ia(iisymp-1+ifac) = 0
enddo

!===============================================================================
! 1. FONCTION ANALYTIQUE SIN(X+2Y+3Z)
!===============================================================================

iuiph = iu
iviph = iv
iwiph = iw

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

call grdcel                                                       &
!==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   rtp(1,iuiph)    , rtp(1,iviph)    , rtp(1,iwiph)    ,          &
   ra     )

! On sort le gradient

namevr = 'Grad_RC'
if (ichrvl.eq.1) then
  call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
  !==========
              ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)
endif

! Calcul de l'erreur absolue

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  rtp(iel,iuiph) = rtp(iel,iuiph)-     cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iviph) = rtp(iel,iviph)-2.d0*cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iwiph) = rtp(iel,iwiph)-3.d0*cos(xx+2.d0*yy+3.d0*zz)
enddo

! On sort l'erreur

namevr = 'Err_Grad_RC'
if (ichrvl.eq.1) then
  call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
  !==========
              ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)
endif


!  2.2 APPEL A GRDCEL AVEC IMRGRA = 1
!  ==================================

imrgra = 1
imligp = 1

call grdcel                                                       &
!==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   rtp(1,iuiph)    , rtp(1,iviph)    , rtp(1,iwiph)    ,          &
   ra     )


! On sort le gradient

namevr = 'Grad_LSQ'
if (ichrvl.eq.1) then
  call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
  !==========
              ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)
endif

! Calcul de l'erreur absolue

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  rtp(iel,iuiph) = rtp(iel,iuiph)-     cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iviph) = rtp(iel,iviph)-2.d0*cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iwiph) = rtp(iel,iwiph)-3.d0*cos(xx+2.d0*yy+3.d0*zz)
enddo

! On sort l'erreur

namevr = 'Err_Grad_LSQ'
if (ichrvl.eq.1) then
  call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
  !==========
              ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)
endif


!  2.3 APPEL A GRDCEL AVEC IMRGRA = 2
!  ==================================

imrgra = 2
imligp = 1

call grdcel                                                       &
!==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   rtp(1,iuiph)    , rtp(1,iviph)    , rtp(1,iwiph)    ,          &
   ra     )

! On sort le gradient

namevr = 'Grad_LSQ_Ext'
if (ichrvl.eq.1) then
  call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
  !==========
              ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)
endif

! Calcul de l'erreur absolue

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  rtp(iel,iuiph) = rtp(iel,iuiph)-     cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iviph) = rtp(iel,iviph)-2.d0*cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iwiph) = rtp(iel,iwiph)-3.d0*cos(xx+2.d0*yy+3.d0*zz)
enddo

! On sort l'erreur

namevr = 'Err_Grad_LSQ_Ext'
if (ichrvl.eq.1) then
  call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
  !==========
              ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)
endif


!  2.4 APPEL A GRDCEL AVEC IMRGRA = 4
!  ==================================

imrgra = 4
imligp = -1

call grdcel                                                       &
!==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   rtp(1,iuiph)    , rtp(1,iviph)    , rtp(1,iwiph)    ,          &
   ra     )

! On sort le gradient

namevr = 'Grad_LSQ_RC'
if (ichrvl.eq.1) then
  call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
  !==========
              ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)
endif

! Calcul de l'erreur absolue

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  rtp(iel,iuiph) = rtp(iel,iuiph)-     cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iviph) = rtp(iel,iviph)-2.d0*cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iwiph) = rtp(iel,iwiph)-3.d0*cos(xx+2.d0*yy+3.d0*zz)
enddo

! On sort l'erreur

namevr = 'Err_Grad_LSQ_RC'
if (ichrvl.eq.1) then
  call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
  !==========
              ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)
endif


!  2.5 APPEL A GRDCEL AVEC IMRGRA = 3
!  ==================================

! Reduction du voisinage etendu

call redvse(anomax)
!==========

imrgra = 3
imligp = 1

call grdcel                                                       &
!==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   rtp(1,iuiph)    , rtp(1,iviph)    , rtp(1,iwiph)    ,          &
   ra     )

! On sort le gradient

namevr = 'Grad_LSQ_ExtRed'
if (ichrvl.eq.1) then
  call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
  !==========
              ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)
endif

! Calcul de l'erreur absolue

do iel = 1, ncelet
  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)
  zz = xyzcen(3,iel)
  rtp(iel,iuiph) = rtp(iel,iuiph)-     cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iviph) = rtp(iel,iviph)-2.d0*cos(xx+2.d0*yy+3.d0*zz)
  rtp(iel,iwiph) = rtp(iel,iwiph)-3.d0*cos(xx+2.d0*yy+3.d0*zz)
enddo

! On sort l'erreur

namevr = 'Err_Grad_LSQ_ExtRed'
if (ichrvl.eq.1) then
  call psteva(ipart , namevr, idimt, ientla, ivarpr,    &
  !==========
              ntcabs, ttcabs, rtp(1,iuiph), rbid, rbid)
endif

!----
! FIN
!----

return
end subroutine
