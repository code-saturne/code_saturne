!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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
! --------
!> \file cs_fuel_physprop.f90
!>
!> \brief Specific physic routine: pulverized coal combustion.
!>        Calculation of \f$\rho\f$ of the mixture
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     mbrom         filling indicator of romb
!> \param[in]     izfppp        zone number of the edge face for the
!>                              specific physic modul
!> \param[in]     rtp           calculated variables at cell centers
!>                              (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!______________________________________________________________________________!

subroutine cs_fuel_physprop &
 ( mbrom  , izfppp ,                                              &
   rtp    , propce )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use dimens, only: nvar
use cstphy
use cstnum
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh
use field
!===============================================================================

implicit none

! Arguments

integer          mbrom
integer          izfppp(nfabor)

double precision rtp(ncelet,nflown:nvar)
double precision propce(ncelet,*)

! Local variables

integer          iel, icla, ipcro2
integer          izone, ifac
integer          iromf , ioxy , nbclip1,nbclip2

double precision x1sro1, x2sro2, srrom1, uns1pw
double precision x2tot, wmolme, unsro1
double precision ff3min,ff3max,valmin,valmax

integer          ipass
data             ipass /0/
save             ipass

integer          iok1,iok2,iok3
double precision , dimension ( : )     , allocatable :: f1m,f2m,f3m,f4m,f5m
double precision , dimension ( : )     , allocatable :: f6m,f7m,f8m,f9m
double precision , dimension ( : )     , allocatable :: enth1 , fvp2m
double precision , dimension ( : )     , allocatable :: xoxyd,enthox
double precision, dimension(:), pointer ::  brom, crom
double precision, dimension(:), pointer :: cpro_x1
!===============================================================================
!
!===============================================================================
! 0. We count the passages
!===============================================================================

ipass = ipass + 1

!===============================================================================
! 1. Initializations to be kept
!===============================================================================

! Massic fraction of gas
call field_get_val_s_by_name("x_c", cpro_x1)

!===============================================================================
! Deallocation dynamic arrays
!----
allocate(f1m(1:ncelet),f2m(1:ncelet),f3m(1:ncelet),              STAT=iok1)
allocate(f4m(1:ncelet),f5m(1:ncelet),                            STAT=iok1)
allocate(f6m(1:ncelet),f7m(1:ncelet),f8m(1:ncelet),f9m(1:ncelet),STAT=iok2)
allocate(enth1(1:ncel),fvp2m(1:ncel),              STAT=iok3)
!----
if ( iok1 > 0 .or. iok2 > 0 .or. iok3 > 0) then
  write(nfecra,*) ' Memory allocation error inside: '
  write(nfecra,*) '     cs_fuel_physprop            '
  call csexit(1)
endif
if ( ieqnox .eq. 1 ) then
!----
  allocate(xoxyd(1:ncelet),enthox(1:ncelet),STAT=iok1)
  !----
  if ( iok1 > 0 ) then
    write(nfecra,*) ' Memory allocation error inside:         '
    write(nfecra,*) '   cs_fuel_physprop for xoxyd and enthox '
  endif
endif
!===============================================================================

!     Pointer on mass density of the cells gas
iromf = ipproc(irom1)
!
!===============================================================================
! 2. Calculation of physic properties of the dispersed phase
!                    cell values
!                    -----------
!    Mass fraction of the solid
!    Diameter
!    Mass density
!===============================================================================
!
call cs_fuel_physprop2 ( ncelet , ncel , rtp , propce )
!=====================

!===============================================================================
! 3. Calculation of physic properties of the gaseous phase
!                    Cell values
!                    -----------
!    Temperature
!    Mass density
!    Concentraions of the gaseous species
!===============================================================================

! --- Calculation of the gas enthalpy  enth1
!        of F1M
!        of F2M
!        of F3M                      in W3=1-F1M-F2M-F4M-F5M-F6M-F7M-F8M-F9M
!        of F4M
!        of F5M
!        of F6M
!        of F7M
!        of F8M
!        of F9M
!        of FVP2M
!
! Initialization of fm
! f1m is always equal to zero. In the context of the fuel oil combustion there
! is only one gaseous fuel.
f1m( : ) = 0.d0
! f2m = ifvap
f2m( : ) = 0.d0
! f3m, f4m, f5m = Oxidizers
f3m( : ) = 0.d0
f4m( : ) = 0.d0
f5m( : ) = 0.d0
! Water vapour
f6m( : ) = 0.d0
! Heterogeneous combustion
f7m( : ) = 0.d0
! f8m, f9m is always equal to zero
f8m( : ) = 0.d0
f9m( : ) = 0.d0

! Mass fraction of gas phase
! - nclacp = nclafu
! - cpro_x1 is directly calculated in function with iyfol
do iel = 1, ncel
  cpro_x1(iel) = 1.d0
enddo

do icla = 1, nclafu
  do iel = 1, ncel
    cpro_x1(iel) = cpro_x1(iel) - rtp(iel,isca(iyfol(icla)))
  enddo
enddo

do iel = 1, ncel
  f2m(iel) =  f2m(iel) + rtp(iel,isca(ifvap))
enddo

if ( ieqnox .eq. 1 ) then
  do iel = 1, ncel
    xoxyd(iel)= (cpro_x1(iel))-f1m(iel)-f2m(iel)
  enddo
endif

ff3min = 1.d+20
ff3max =-1.d+20
nbclip1= 0
nbclip2= 0
valmin = 1.d+20
valmax =-1.d+20

do iel = 1, ncel
  uns1pw = 1.d0/cpro_x1(iel)

   ! -At the moment, the variable <noxyd> does not exist in the fuel oil version
   !  because it is not declared in fulecd.f90. We have to add it in this
   !  subroutine (keyword: oxycombustion).
  if ( noxyd .ge. 2 ) then
    f4m(iel) = rtp(iel,isca(if4m))
    if ( noxyd .eq. 3 ) then
      f5m(iel) = rtp(iel,isca(if5m))
    endif
  endif

  ! - Relative scalar for the heterogeneous combustion.
  f7m(iel) =  rtp(iel,isca(if7m))

  ! The transported variance.
  fvp2m(iel) = rtp(iel,isca(ifvp2m))

  ! Units: [kg scalars / kg gas]
  f1m(iel)  = f1m(iel)    *uns1pw
  f2m(iel)  = f2m(iel)    *uns1pw
  f4m(iel)  = f4m(iel)    *uns1pw
  f5m(iel)  = f5m(iel)    *uns1pw
  f6m(iel)  = f6m(iel)    *uns1pw
  f7m(iel)  = f7m(iel)    *uns1pw
  f8m(iel)  = f8m(iel)    *uns1pw
  f9m(iel)  = f9m(iel)    *uns1pw

  fvp2m(iel)= fvp2m(iel)*uns1pw

  f3m(iel) = 1.d0                                        &
           -( f1m(iel)+f2m(iel)+f4m(iel)+f5m(iel)        &
             +f6m(iel)+f7m(iel)+f8m(iel)+f9m(iel) )

  ff3max = max(ff3max,f3m(iel))
  ff3min = min(ff3min,f3m(iel))

  if ( ieqnox .eq. 1 ) then
    enthox(iel) = rtp(iel,isca(ihox))/xoxyd(iel)
  endif

enddo

if ( irangp .ge. 0 ) then
  call parmin(ff3min)
  call parmax(ff3max)
  call parcpt(nbclip1)
  call parcpt(nbclip2)
  call parmin(valmin)
  call parmax(valmax)
endif
write(nfecra,*) ' Values of f3 min and max: ',FF3MIN,FF3MAX
if ( nbclip1 .gt. 0 ) then
  write(nfecra,*) ' Clipping phase gas variance in min:',nbclip1,valmin
endif
if ( nbclip2 .gt. 0 ) then
  write(nfecra,*) ' Clipping phase gas variance in max:',nbclip1,valmin
endif

! ---- Gas enthalpy H1
enth1( : ) =0.D0
do icla = 1, nclafu
  do iel = 1, ncel
    enth1(iel) =  enth1(iel) + rtp(iel,isca(ih2(icla)))
  enddo
enddo
do iel = 1, ncel
  enth1(iel) = (rtp(iel,isca(iscalt))-enth1(iel))/ cpro_x1(iel)
enddo

call cs_fuel_physprop1 &
!=====================
 ( ncelet , ncel   ,                                      &
   f1m    , f2m    , f3m    , f4m    , f5m    ,           &
   f6m    , f7m    , f8m    , f9m    , fvp2m  ,           &
   enth1  , enthox ,                                      &
   rtp    , propce , propce(1,iromf)   )

!===============================================================================
! 4. Calculation of physics properties of the dispersed phase
!                    Cell values
!                    -----------
!    Temperature
!===============================================================================

! --- Transport of H2

call  cs_fuel_thfieldconv2 ( ncelet , ncel , rtp , propce )
!=========================

!===============================================================================
! 5. Calculation of the mixture physic properties
!                    Cell values
!                    -----------
!    Mass density
!===============================================================================
! --- Calculation of Rho of the mixture: 1/Rho = X1/Rho1 + Sum(X2/Rho2)
!     We under discharge when we have a rho n at our disposal, ie
!       from the second passage or
!       from the first one if we are in continuation of the calculation and
!         we have reread the mass density in the file continuation.

call field_get_val_s(icrom, crom)

if (ipass.gt.1.or.(isuite.eq.1.and.initro.eq.1)) then
  srrom1 = srrom
else
  srrom1 = 0.d0
endif

do iel = 1, ncel
  x2sro2 = zero

  do icla = 1, nclafu
    ipcro2 = ipproc(irom2(icla))
    x2sro2 = x2sro2 + rtp(iel,isca(iyfol(icla))) / propce(iel,ipcro2)
  enddo
  x1sro1 = cpro_x1(iel) / propce(iel,iromf)
  ! ---- Under eventual relaxation to give in ppini1.F
  crom(iel) = srrom1*crom(iel)                  &
                     + (1.d0-srrom1)/(x1sro1+x2sro2)
enddo


!===============================================================================
! 6. Calculation of the rho of the mixture
!                      Facet value
!                      -----------
!===============================================================================

mbrom = 1
call field_get_val_s(ibrom, brom)
call field_get_val_s(icrom, crom)
! ---> Mass density at the edge for all facets
!      The input facets will be recalculated.
!
do ifac = 1, nfabor
  iel = ifabor(ifac)
  brom(ifac) = crom(iel)
enddo

! ---> Mass density at the edge for ONLY all input facets
!     The test on izone serves with the calculation resumption

if ( ipass.gt.1 .or. isuite.eq.1 ) then
  do ifac = 1, nfabor

    izone = izfppp(ifac)
    if(izone.gt.0) then
      if ( ientat(izone).eq.1 .or. ientfl(izone).eq.1 ) then
        x2sro2 = zero
        x2tot  = zero
        do icla = 1, nclafu
          x2sro2 = x2sro2 + x20(izone,icla)/rho0fl
          x2tot  = x2tot  + x20(izone,icla)
        enddo

        ioxy = inmoxy(izone)
        wmolme =( oxyo2(ioxy)+oxyn2(ioxy)                         &
                 +oxyh2o(ioxy)+oxyco2(ioxy))                      &
               /( wmole(io2) *oxyo2(ioxy)                         &
                 +wmole(in2) *oxyn2(ioxy)                         &
                 +wmole(ih2o)*oxyh2o(ioxy)                        &
                 +wmole(ico2)*oxyco2(ioxy) )

        unsro1 = (wmolme*rr*timpat(izone)) / p0
        x1sro1 = (1.d0-x2tot) * unsro1
        brom(ifac) = 1.d0 / (x1sro1+x2sro2)
      endif
    endif

  enddo
endif
!--------
! Formats
!--------

!===============================================================================
! Deallocation dynamic arrays
deallocate(f1m,f2m,f3m,f4m,f5m,STAT=iok1)
deallocate(f6m,f7m,f8m,f9m,    STAT=iok2)
deallocate(enth1,fvp2m,     STAT=iok3)

if (iok1 > 0 .or. iok2 > 0 .or. iok3 > 0) then
  write(nfecra,*) ' Memory deallocation error inside: '
  write(nfecra,*) '     cs_fuel_physprop              '
  call csexit(1)
endif
if (ieqnox .eq. 1) then

  deallocate(xoxyd,enthox)

  if (iok1 > 0) then
    write(nfecra,*) ' Memory deallocation error inside:      '
    write(nfecra,*) '   cs_fuel_physprop for xoxyd and enthox'
    call csexit(1)
  endif
endif
!===============================================================================

!----
! End
!----
return
end subroutine
