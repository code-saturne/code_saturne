!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

!> \file cs_tagmri.f90
!>
!> \brief The 1D thermal model to compute the temperature to impose
!> at the cold wall. This one is used by the COPAIN model to estimate
!> the heat flux at the wall where the condensation occurs.
!>
!> This subroutine is used to compute at each face the
!> \f$ \mbox{tmur} \f$ at cold wall.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nfabor        number of boundary faces
!> \param[in]     isvtb         scalar number coupled
!> \param[in,out] icodcl        boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradt \, \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right) \f$
!_______________________________________________________________________________

subroutine cs_tagmri &
 ( nfabor ,                                                       &
   isvtb  , icodcl ,                                              &
   rcodcl )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstnum
use cstphy
use dimens, only: nvar
use entsor
use pointe
use field
use mesh, only:ifabor
use cs_nz_condensation, only: izzftcd, iztag1d, ztpar
use cs_nz_tagmr, only: ztmur

!===============================================================================

implicit none

! Arguments

integer          nfabor
integer          isvtb  , icodcl(nfabor,nvar)
double precision rcodcl(nfabor,nvar,3)

! Local variables


integer          ii, iz, ivar
integer          ifac, iel
integer          icldef
double precision temper, enthal

double precision, dimension(:), pointer :: cpro_cp

!===============================================================================


! Sans specification, une face couplee est une face de type paroi

icldef = 5

ivar = isca(isvtb)

do ii = 1, nfbpcd

  ifac = ifbpcd(ii)
  iz = izzftcd(ii)
  if(iztag1d(iz).eq.1) then
    icodcl(ifac,ivar)   = 1
    rcodcl(ifac,ivar,1) = ztmur(ii,1)
    rcodcl(ifac,ivar,2) = rinfin
    rcodcl(ifac,ivar,3) = 0.d0
  else
    icodcl(ifac,ivar)   = 1
    rcodcl(ifac,ivar,1) = ztpar(iz)
    rcodcl(ifac,ivar,2) = rinfin
    rcodcl(ifac,ivar,3) = 0.d0
  endif
enddo

! Conversion eventuelle temperature -> enthalpie

if (isvtb.eq.iscalt .and. itherm.eq.2) then

! --- Specific heat
if (icp.ge.0) then
  call field_get_val_s(icp, cpro_cp)
! --- Stop if Cp is not variable
else
  write(nfecra,1000) icp
  call csexit (1)
endif

  do ii = 1, nfbpcd

    ifac = ifbpcd(ii)
    iel = ifabor(ifac)

    temper = rcodcl(ifac,ivar,1)

    enthal = (temper + tkelvi) * cpro_cp(iel)

    rcodcl(ifac,ivar,1) = enthal

  enddo

endif

!--------
! Formats
!--------

 1000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@    Inconsistent calculation data',/,                           &
'@',/,                                                            &
'@      usipsu specifies that the specific heat is uniform',/,    &
'@        icp = ',i10   ,' while',/,                              &
'@      cs_user_physical_properties prescribes a variable specific heat.',/,&
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipsu or cs_user_physical_properties.',/,           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

!----
! End
!----

return
end subroutine
