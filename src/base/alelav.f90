!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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
!> \param[in]     propce        physical properties at cell centers
!_______________________________________________________________________________

subroutine alelav &
 ( propce )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
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
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

double precision propce(ncelet,*)

! Local variables

character(len=80) :: chaine
integer          iel   , isou  , jsou  , ifac
integer          ipcvmx, ipcvmy, ipcvmz
integer          iflmas, iflmab
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          ivisep
integer          iswdyp, idftnp, icvflb
integer          ivoid(1)

double precision blencp, epsilp, epsrgp, climgp, thetv
double precision epsrsp, prosrf
double precision relaxp
double precision hint, distbf, srfbn2
double precision rinfiv(3), pimpv(3)

double precision rvoid(1)

double precision, dimension(:,:), pointer :: cfaale, claale
double precision, dimension(:,:,:), pointer :: cfbale, clbale

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:,:) :: smbr
double precision, allocatable, dimension(:,:,:) :: fimp
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: brom
double precision, dimension(:,:), pointer :: mshvel, mshvela

!===============================================================================

!===============================================================================
! 1. INITIALIZATION
!===============================================================================

! Allocate temporary arrays for the radiative equations resolution
allocate(viscf(nfac), viscb(nfabor))
allocate(smbr(3,ncelet))
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

call field_get_val_v(ivarfl(iuma), mshvel)
call field_get_val_prev_v(ivarfl(iuma), mshvela)

if(iwarni(iuma).ge.1) then
  write(nfecra,1000)
endif

call field_get_coefa_v(ivarfl(iuma), claale)
call field_get_coefb_v(ivarfl(iuma), clbale)
call field_get_coefaf_v(ivarfl(iuma), cfaale)
call field_get_coefbf_v(ivarfl(iuma), cfbale)

! We compute the boundary condition on the mesh velocity at the free surface
! from the new mass flux.

! Density at the boundary
call field_get_val_s(ibrom, brom)

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
    pimpv(1) = gx*bmasfl(ifac)/(brom(ifac)*prosrf)
    pimpv(2) = gy*bmasfl(ifac)/(brom(ifac)*prosrf)
    pimpv(3) = gz*bmasfl(ifac)/(brom(ifac)*prosrf)

    call set_dirichlet_vector &
         !====================
       ( claale(:,ifac)  , cfaale(:,ifac)  ,             &
         clbale(:,:,ifac), cfbale(:,:,ifac),             &
         pimpv           , hint            , rinfiv )
  endif
enddo

!===============================================================================
! 2. SOLVING OF THE MESH VELOCITY EQUATION
!===============================================================================

if (iwarni(iuma).ge.1) then
  call field_get_name(ivarfl(iuma), chaine)
  write(nfecra,1100) chaine(1:16)
endif

do iel = 1, ncelet
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
  call visort &
  !==========
( imvisf ,                                                       &
  propce(1,ipcvmx), propce(1,ipcvmy), propce(1,ipcvmz),          &
  viscf  , viscb  )
endif

iconvp = iconv (iuma)
idiffp = idiff (iuma)
ndircp = ndircl(iuma)
nswrsp = nswrsm(iuma)
nswrgp = nswrgr(iuma)
imligp = imligr(iuma)
ircflp = ircflu(iuma)
ischcp = ischcv(iuma)
isstpp = isstpc(iuma)
iescap = 0
idftnp = idften(iuma)
iswdyp = iswdyn(iuma)
iwarnp = iwarni(iuma)
blencp = blencv(iuma)
epsilp = epsilo(iuma)
epsrsp = epsrsm(iuma)
epsrgp = epsrgr(iuma)
climgp = climgr(iuma)
relaxp = 1.d0
thetv  = 1.d0
! all boundary convective flux with upwind
icvflb = 0

! we do not take into account the transpose of grad
ivisep = 0

call coditv &
!==========
 ( idtvar , iuma   , iconvp , idiffp , ndircp ,                   &
   imrgra , nswrsp , nswrgp , imligp , ircflp , ivisep ,          &
   ischcp , isstpp , iescap , idftnp , iswdyp ,                   &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp ,                   &
   relaxp , thetv  ,                                              &
   mshvela         , mshvela         ,                            &
   claale , clbale , cfaale , cfbale ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  , viscf  , viscb  ,          &
   icvflb , ivoid  ,                                              &
   fimp   ,                                                       &
   smbr   ,                                                       &
   mshvel ,                                                       &
   rvoid  )

! Free memory
deallocate(viscf, viscb)
deallocate(smbr, fimp)

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000    format(/,                                                &
'   ** RESOLUTION DE LA VITESSE DE MAILLAGE'                   ,/,&
'      ------------------------------------'                   ,/)
 1100    format(/,'           RESOLUTION POUR LA VARIABLE ', a16,/)

#else

 1000    format(/,                                                &
'   ** SOLVING MESH VELOCITY'                                  ,/,&
'      ---------------------'                                  ,/)
 1100    format(/,'           SOLVING VARIABLE ', a16          ,/)

#endif

!----
! END
!----

return

end subroutine
