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

!> \file radf08.f90
!>
!> \brief Subroutine of the ADF radiation model
!>
!>  Determination of the radiation coeffcients of the ADF model as well as the
!>  corresponding weights.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     pco2          CO2 volume fraction
!> \param[in]     ph2o          H2O volume fraction
!> \param[in]     fv            Soot volume fraction
!> \param[in]     teloc         Gas temperature
!> \param[in]     kloc          Radiation coeffcient of the i different gases
!> \param[in]     aloc          Weights of the i different gases in cells
!> \param[in]     alocbo        Weights of the i different gases at boundaries
!_______________________________________________________________________________


subroutine radf08 &
 ( pco2 , ph2o , fv , teloc, kloc, aloc, alocbo)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use radiat
use mesh
use field
!===============================================================================

implicit none

! Arguments
double precision pco2(ncelet),ph2o(ncelet),fv(ncelet),teloc(ncelet)
double precision kloc(ncelet,nwsgg),aloc(ncelet,nwsgg),alocbo(nfabor,nwsgg)

! Local variables

integer          iel,ifac,i,j,k,l,it,ix,ipass,ntsto,nysto
double precision tref,xh2oref,rt,rx,kmloc

double precision, allocatable, dimension(:) ::     y,tsto,ysto
double precision, allocatable, dimension(:,:,:) :: asto,ksto2

double precision, dimension(:), pointer :: tpaadf ! Points to the walls
! temperature table

data ipass /0/

save ipass,ksto2,asto,ntsto,nysto,tsto,ysto ! These local values must be stored.
! Otherwise they have to be read again for every iteration cycle.
!===============================================================================
! 0 - GESTION MEMOIRE
!===============================================================================

allocate(y(ncelet))

call field_get_val_s(itparo,tpaadf)

!===============================================================================
!  1 - COEFFICIENT D'ABSORPTION DU MELANGE GAZEUX (m-1)
!===============================================================================
ipass = ipass + 1
! The ADF data base is read only once during the very first iteration
if(ipass.eq.1) then
  open(unit=10,file='dp_radiat_ADF8')
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*) ntsto  ! Number of tabulated gas phase temperatures
  !!
  !!     allocate parameters
  !!
  allocate   (tsto(ntsto))
  read(10,*)
  read(10,*) (tsto(i),i=1,ntsto) ! Tabulated gas phase temperatures
  read(10,*)
  read(10,*)  nysto ! Number of tabulated ratios ph2o/pco2
  allocate   (ysto(nysto))
  read(10,*)
  read(10,*) (ysto(i),i=1,nysto) ! Tabulated ratios ph2o/pco2
  read(10,*)
  read(10,*) tref,xh2oref ! Reference temperature and h2o volume fraction
  read(10,*)
  allocate   (asto(nwsgg,nysto,ntsto))
  allocate   (ksto2(nwsgg,nysto,ntsto))
  !!
  !!   READING THE PARAMETERS
  !!
  do i=1,nwsgg
     read(10,*)
     do j=1,ntsto
       read(10,*) (ksto2(i,k,j),k=1,nysto),(asto(i,k,j),k=1,nysto)
                  !ksto2(i,k,j): Radiation coeffcient of the i-th grey gas,
                  !              the k-th ratio of ph2o/pco2, and
                  !              the j-th tabulated temperature
                  !asto2(i,k,j): Weight of the i-th grey gas,
                  !              the k-th ratio of ph2o/pco2, and
                  !              the j-th tabulated temperature

     enddo
  enddo
endif

do iel = 1,ncel
  if(pco2(iel).gt.0.d0) then
    y(iel) = ph2o(iel)/pco2(iel)
  else
    y(iel) = ysto(nysto)
  endif

! Interpolation temperature
  if(teloc(iel).le.tsto(1)) then
    rt=0.d0
    it=1
  else if(teloc(iel).ge.tsto(ntsto)) then
    rt=1.d0
    it=ntsto-1
  else
     l=1
     do while(teloc(iel).gt.tsto(l))
       l=l+1
     enddo
     it=l-1
     rt=(teloc(iel)-tsto(it))/(tsto(it+1)-tsto(it))
  endif
! Interpolation H2O-molefraction
  if(y(iel).le.ysto(1)) then
    rx=0.d0
    ix=1
  else if(y(iel).ge.ysto(nysto)) then
    rx=1.d0
    ix=nysto-1
  else
     l=1
     do while(y(iel).gt.ysto(l))
       l=l+1
     enddo
     ix=l-1
     rx=(y(iel)-ysto(ix))/(ysto(ix+1)-ysto(ix))
  endif
  ! Absortion Coeffcient
  do i = 1, nwsgg
    kmloc = (1._8-rt)*(1._8-rx)*ksto2(i,ix,it)+(1._8-rt)*(rx)*ksto2(i,ix+1,it)+&
                 rt*(1._8-rx)*ksto2(i,ix,it+1)+rt*rx*ksto2(i,ix+1,it+1)
    kloc(iel,i) = pco2(iel)*kmloc*100.d0*(p0/1.d5) ! Local radiation coeffcient of the i-th grey gas
    aloc(iel,i) = (1._8-rt)*(1._8-rx)*asto(i,ix,it)+(1._8-rt)*(rx)*asto(i,ix+1,it)+&
                 rt*(1._8-rx)*asto(i,ix,it+1)+rt*rx*asto(i,ix+1,it+1)
               ! Local weight of the i-th grey gas
  enddo
enddo
do ifac = 1, nfabor
  iel = ifabor(ifac)
  if (pco2(iel).gt.0.d0) then
    y(iel) = ph2o(iel)/pco2(iel)
  else
    y(iel) = ysto(nysto)
  endif
! Interpolation temperature
  if (tpaadf(ifac).le.tsto(1)) then
    rt=0.d0
    it=1
  else if(tpaadf(ifac).ge.tsto(ntsto)) then
    rt=1.d0
    it=ntsto-1
  else
     l=1
     do while(tpaadf(ifac).gt.tsto(l))
       l=l+1
     enddo
     it=l-1
     rt=(tpaadf(ifac)-tsto(it))/(tsto(it+1)-tsto(it))
  endif
! Interpolation H2O-molefraction
  if(y(iel).le.ysto(1)) then
    rx=0.d0
    ix=1
  else if(y(iel).ge.ysto(nysto)) then
    rx=1.d0
    ix=nysto-1
  else
     l=1
     do while(y(iel).gt.ysto(l))
       l=l+1
     enddo
     ix=l-1
     rx=(y(iel)-ysto(ix))/(ysto(ix+1)-ysto(ix))
  endif
! Absortion Coeffcient
  do i=1,nwsgg
    alocbo(ifac,i)=(1._8-rt)*(1._8-rx)*asto(i,ix,it)+(1._8-rt)*(rx)*           &
                   asto(i,ix+1,it)+rt*(1._8-rx)*asto(i,ix,it+1)+               &
                   rt*rx*asto(i,ix+1,it+1)
                   ! Local weight of the i-th grey gas
  enddo
enddo

return

end subroutine
