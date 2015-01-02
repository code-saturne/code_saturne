!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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
!> \file vandri.f90
!> \brief Imposition of an amortization of Van Driest type for the LES.
!>        \f$ \nu_T \f$ is aborsorbed by \f$ (1-\exp(\dfrac{-y^+}{d^+}))^2 \f$
!>        where \f$ d^+ \f$ is set at 26.

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
! modename        name          role
!______________________________________________________________________________!
!> \param[in]     itypfb        boundary face types
!> \param[in]     ifapat        number of the edge face code 5 the nearst
!>                                 (R_{ij} and wall echo)
!> \param[in]     visvdr        dynamic viscosity in edge cells after
!>                               driest velocity amortization
!> \param[in]     yplusc        \f$ y^+\f$ value in cells in the
!>                                 case \c abs(icdpar).eq.1
!______________________________________________________________________________!

subroutine vandri &
 (  itypfb , ifapat , visvdr , yplusc )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use entsor
use cstphy
use parall
use pointe, only: uetbor
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          itypfb(nfabor)
integer          ifapat(ncelet)

double precision visvdr(ncelet)
double precision yplusc(ncelet)

! Local variables

integer          iel   , ifac
double precision yplus , yminpa, viscos
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: viscl, visct

!===============================================================================

call field_get_val_s(icrom, crom)
call field_get_val_s(iprpfl(iviscl), viscl)
call field_get_val_s(iprpfl(ivisct), visct)

!     Direst calculation of the distance to the wall
!     (non compatible parall/perio)
if(abs(icdpar).eq.2) then

  !     In sequentiel, nothing to report
  if(irangp.lt.0) then
    do iel = 1, ncel
      ifac = ifapat(iel)
      viscos = viscl(iel)/crom(iel)
      yminpa = sqrt((cdgfbo(1,ifac)-xyzcen(1,iel))**2             &
           +        (cdgfbo(2,ifac)-xyzcen(2,iel))**2             &
           +        (cdgfbo(3,ifac)-xyzcen(3,iel))**2)
      yplus = uetbor(ifac) * yminpa/ viscos
      visct(iel) = visct(iel)*                                    &
           (1.0d0-exp(-yplus/cdries))**2
    enddo
  !     In parallel, we absorb only the first mesh of the wall:
  !     dangereous but a priori useless (because the use of
  !     icdpar=+/-2 in parallel is bloqued in verini)
  else
    write(nfecra,1000)
    do ifac = 1, nfabor
      if(itypfb(ifac).eq.iparoi .or.                              &
         itypfb(ifac).eq.iparug ) then
        iel = ifabor(ifac)
        viscos = viscl(iel)/crom(iel)
        yminpa = sqrt((cdgfbo(1,ifac)-xyzcen(1,iel))**2           &
             +        (cdgfbo(2,ifac)-xyzcen(2,iel))**2           &
             +        (cdgfbo(3,ifac)-xyzcen(3,iel))**2)
        yplus = uetbor(ifac) * yminpa/ viscos
        visct(iel) = visct(iel)*                                  &
             (1.0d0-exp(-yplus/cdries))**2
      endif
    enddo
endif

!     Nex calculation mode: it is simpler
elseif(abs(icdpar).eq.1) then
  do iel = 1, ncel
    yplus = yplusc(iel)
    visct(iel) = visct(iel)*                                      &
           (1.0d0-exp(-yplus/cdries))**2
  enddo
endif

!     For the wall cells we add the turbulent viscosity which was absorbed
!     in clptur and which has served to calculate the boundary conditions
do iel = 1, ncel
  if (visvdr(iel).gt.-900.d0)                                     &
       visct(iel) = visvdr(iel)
enddo

!--------
! Formats
!--------
#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ Attention: Dans le cas de la LES avec amortissement     ',/,&
'@    =========                                               ',/,&
'@  L''amortissement de Van Driest n''est fait que            ',/,&
'@  sur la premiere cellule a la paroi en cas de parallelisme ',/,&
'@                                                            ',/,&
'@  Le calcul se poursuit.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                           &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: In case of LES with damping'                    ,/,&
'@    ========'                                                ,/,&
'@    Van Driest damping is only effective on the first cell'  ,/,&
'@    off-wall in case of parallelism'                         ,/,&
'@'                                                            ,/,&
'@  The calculation will be run.'                              ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#endif
!----
! End
!----


return
end subroutine
