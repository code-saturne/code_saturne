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

!===============================================================================
! Function:
! ---------

!> \file alelav.f90
!>
!> \brief This subroutine performs the solving of a Poisson equation
!> on the mesh velocity for ALE module.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in,out] rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in,out] propfb        physical properties at boundary face centers
!_______________________________________________________________________________

subroutine alelav &
 ( rtp    , rtpa   , propce , propfb )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use entsor
use optcal
use cstnum
use cstphy
use pointe
use albase, only: ialtyb
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

double precision rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(ndimfb,*)

! Local variables

integer          iel   , isou  , jsou  , ifac
integer          ipcvmx, ipcvmy, ipcvmz, ippu  , ippv  , ippw
integer          iflmas, iflmab, ipbrom
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp, ivisep
integer          iswdyp, idftnp

double precision blencp, epsilp, epsrgp, climgp, extrap, thetv
double precision epsrsp, prosrf
double precision relaxp
double precision hint, distbf, srfbn2
double precision rinfiv(3), pimpv(3)

double precision rvoid(1)

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:,:) :: smbr, meshv, meshva
double precision, allocatable, dimension(:,:,:) :: fimp
double precision, dimension(:), pointer :: imasfl, bmasfl

!===============================================================================

!===============================================================================
! 1. INITIALIZATION
!===============================================================================

! Allocate temporary arrays for the radiative equations resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(smbr(3,ncelet))
allocate(meshv(3,ncelet), meshva(3,ncelet))
allocate(fimp(3,3,ncelet))

rinfiv(1) = rinfin
rinfiv(2) = rinfin
rinfiv(3) = rinfin
ipcvmx = ipproc(ivisma(1))
ipcvmy = ipproc(ivisma(2))
ipcvmz = ipproc(ivisma(3))
! The mass flux is necessary to call coditv but not used (ICONV=0)
! Except for the free surface, where it is used as a Boundary condition
call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

if(iwarni(iuma).ge.1) then
  write(nfecra,1000)
endif

! We compute the boundary condition on the mesh velocity at the free surface
! from the new mass flux.

! Density at the boundary
ipbrom = ipprob(irom)

! The mesh move in the direction of the gravity in case of free-surface
do ifac = 1, nfabor
  if (ialtyb(ifac) .eq. ifresf) then

    iel = ifabor(ifac)
    distbf = distb(ifac)
    srfbn2 = surfbn(ifac)**2
    if (ipcvmx.eq.ipcvmy) then
      hint = propce(iel,ipproc(ivisma(1)))/distbf
    else !FIXME
      hint = ( propce(iel,ipproc(ivisma(1)))*surfbo(1,ifac)**2    &
             + propce(iel,ipproc(ivisma(2)))*surfbo(2,ifac)**2    &
             + propce(iel,ipproc(ivisma(3)))*surfbo(3,ifac)**2 )  &
           /distbf/srfbn2
    endif

    prosrf = gx*surfbo(1,ifac) + gy*surfbo(2,ifac) + gz*surfbo(3,ifac)
    pimpv(1) = gx*propfb(ifac,iflmab)/(propfb(ifac,ipbrom)*prosrf)
    pimpv(2) = gy*propfb(ifac,iflmab)/(propfb(ifac,ipbrom)*prosrf)
    pimpv(3) = gz*propfb(ifac,iflmab)/(propfb(ifac,ipbrom)*prosrf)

    call set_dirichlet_vector &
         !====================
       ( claale(1,ifac)  , cfaale(1,ifac)  ,             &
         clbale(1,1,ifac), cfbale(1,1,ifac),             &
         pimpv           , hint            , rinfiv )
  endif
enddo

!===============================================================================
! 2. SOLVING OF THE MESH VELOCITY EQUATION
!===============================================================================


ippu = ipprtp(iuma)
ippv = ipprtp(ivma)
ippw = ipprtp(iwma)

if(iwarni(iuma).ge.1) then
  write(nfecra,1100) nomvar(ippu)
  write(nfecra,1100) nomvar(ippv)
  write(nfecra,1100) nomvar(ippw)
endif

do iel = 1, ncelet
  meshva(1,iel) = rtpa(iel,iuma)
  meshva(2,iel) = rtpa(iel,ivma)
  meshva(3,iel) = rtpa(iel,iwma)
  do isou = 1, 3
    smbr(isou,iel)   = 0.d0
    do jsou = 1, 3
      fimp(isou,jsou,iel) = 0.d0
    enddo
  enddo
enddo

if (ipcvmx.eq.ipcvmy) then
  call viscfa &
  !==========
( imvisf ,                                                       &
  propce(1,ipcvmx),                                              &
  viscf  , viscb  )
else
!TODO viortv
  call visort &
  !==========
( imvisf ,                                                       &
  propce(1,ipcvmx), propce(1,ipcvmy), propce(1,ipcvmz),          &
  viscf  , viscb  )
endif

iconvp = iconv (iuma)
idiffp = idiff (iuma)
ireslp = iresol(iuma)
ndircp = ndircl(iuma)
nitmap = nitmax(iuma)
nswrsp = nswrsm(iuma)
nswrgp = nswrgr(iuma)
imligp = imligr(iuma)
ircflp = ircflu(iuma)
ischcp = ischcv(iuma)
isstpp = isstpc(iuma)
iescap = 0
idftnp = idften(iuma)
iswdyp = iswdyn(iuma)
imgrp  = imgr  (iuma)
ncymxp = ncymax(iuma)
nitmfp = nitmgf(iuma)
iwarnp = iwarni(iuma)
blencp = blencv(iuma)
epsilp = epsilo(iuma)
epsrsp = epsrsm(iuma)
epsrgp = epsrgr(iuma)
climgp = climgr(iuma)
extrap = extrag(iuma)
relaxp = 1.d0
thetv  = 1.d0

! we do not take into account the transpose of grad
ivisep = 0

call coditv &
!==========
 ( idtvar , iuma   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp , ivisep ,          &
   ischcp , isstpp , iescap , idftnp , iswdyp ,                   &
   imgrp  , ncymxp , nitmfp , ippu   , ippv   , ippw   , iwarnp , &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   meshva , meshva ,                                              &
   claale , clbale , cfaale , cfbale ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  , viscf  , viscb  ,          &
   fimp   ,                                                       &
   smbr   ,                                                       &
   meshv  ,                                                       &
   rvoid  )

do iel = 1, ncelet
   rtp(iel,iuma) = meshv(1,iel)
   rtp(iel,ivma) = meshv(2,iel)
   rtp(iel,iwma) = meshv(3,iel)
enddo

! Free memory
deallocate(viscf, viscb)
deallocate(smbr, fimp)
deallocate(meshv, meshva)

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000    format(/,                                                &
'   ** RESOLUTION DE LA VITESSE DE MAILLAGE                   ',/,&
'      ------------------------------------                   ',/)
 1100    format(/,'           RESOLUTION POUR LA VARIABLE ',A8,/)

#else

 1000    format(/,                                                &
'   ** SOLVING MESH VELOCITY'                                  ,/,&
'      ---------------------'                                  ,/)
 1100    format(/,'           SOLVING VARIABLE ',A8           ,/)

#endif

!----
! END
!----

return

end subroutine
