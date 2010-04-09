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
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  , nvar   ,          &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtp    ,                                                       &
   coefa  , coefb  ,                                              &
   rdevel , rtuser , ra     )

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
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
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
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1,2,3,4,5,6     ! ra ! --- ! work arrays                                    !
!  (ncelet)        !    !     !  (computation of pressure gradient)            !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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
include "cstnum.h"
include "pointe.h"
include "entsor.h"
include "albase.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr , nphas , nvar
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision rtp(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra, ifinia, ifinra
integer          ifac  , iel   , ivar  , iphas
integer          inc   , iccocg, iphydp
integer          iuiph , iviph , iwiph
integer          nswrgp, imligp, iwarnp
integer          ipclip
integer          iw1   , iw2   , iw3
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

iw1    = idbra0
iw2    = iw1    + ncelet
iw3    = iw2    + ncelet
ifinra = iw3    + ncelet

call rasize('testel', ifinra)
!==========

! Symmetry type:
! value 0 avoids extrapolating the gradient on boundary faces.
do ifac = 1, nfabor
   ia(iisymp-1+ifac) = 0
enddo

!===============================================================================
! 1. FONCTION ANALYTIQUE SIN(X+2Y+3Z)
!===============================================================================

iphas = 1
iuiph = iu(iphas)
iviph = iv(iphas)
iwiph = iw(iphas)

ivar   = ipr(iphas)
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
iphydp = 0

!  2.1 APPEL A GRDCEL AVEC IMRGRA = 0
!  ==================================

imrgra = 0
imligp = -1

call grdcel                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,  iphydp ,&
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   ra(iw1), ra(iw1), ra(iw1),                                     &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   rtp(1,iuiph)    , rtp(1,iviph)    , rtp(1,iwiph)    ,          &
   ra(iw1), ra(iw2), ra(iw3),                                     &
   rdevel , rtuser , ra     )

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
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,  iphydp ,&
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   ra(iw1), ra(iw1), ra(iw1),                                     &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   rtp(1,iuiph)    , rtp(1,iviph)    , rtp(1,iwiph)    ,          &
   ra(iw1), ra(iw2), ra(iw3),                                     &
   rdevel , rtuser , ra     )


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
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,  iphydp ,&
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   ra(iw1), ra(iw1), ra(iw1),                                     &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   rtp(1,iuiph)    , rtp(1,iviph)    , rtp(1,iwiph)    ,          &
   ra(iw1), ra(iw2), ra(iw3),                                     &
   rdevel , rtuser , ra     )

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
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,  iphydp ,&
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   ra(iw1), ra(iw1), ra(iw1),                                     &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   rtp(1,iuiph)    , rtp(1,iviph)    , rtp(1,iwiph)    ,          &
   ra(iw1), ra(iw2), ra(iw3),                                     &
   rdevel , rtuser , ra     )

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
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,  iphydp ,&
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   ra(iw1), ra(iw1), ra(iw1),                                     &
   rtp(1,ivar)     , coefa(1,ipclip) , coefb(1,ipclip) ,          &
   rtp(1,iuiph)    , rtp(1,iviph)    , rtp(1,iwiph)    ,          &
   ra(iw1), ra(iw2), ra(iw3),                                     &
   rdevel , rtuser , ra     )

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
