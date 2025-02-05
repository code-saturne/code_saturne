!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

subroutine lwcphy

!===============================================================================
! FONCTION :
! --------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME DE PREMELANGE MODELE LWC
! Calcul de RHO adiabatique ou permeatique (transport de H)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use coincl
use ppincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

! Local variables

integer          igg, iel
integer          ifac
double precision coefg(ngazgm)
double precision, dimension(:), pointer :: brom,  crom
double precision, dimension(:), pointer :: bsval
double precision, dimension(:), pointer :: cvar_yfm, cvar_yfp2m
double precision, dimension(:), pointer :: cvar_fm, cvar_fp2m
double precision, dimension(:), pointer :: cvar_coyfp, cpro_ymgg

!===============================================================================

interface

  subroutine cs_combustion_boundary_conditions_density_ebu_lw()  &
    bind(C, name='cs_combustion_boundary_conditions_density_ebu_lw')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_combustion_boundary_conditions_density_ebu_lw

end interface

!===============================================================================
! 1. INITIALISATIONS A CONSERVER
!===============================================================================

! --- Initialisation memoire


! ---> Initialisation

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

! ---> Positions des variables, coefficients

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

call field_get_val_s(ifm, cvar_fm)
call field_get_val_s(ifp2m, cvar_fp2m)
call field_get_val_s(iyfm, cvar_yfm)
call field_get_val_s(iyfp2m, cvar_yfp2m)
if (ippmod(icolwc).ge.2) call field_get_val_s(icoyfp, cvar_coyfp)

!===============================================================================
! 2. DETERMINATION DES GRANDEURS THERMOCHIMIQUES MOYENNES
!===============================================================================


if ((ippmod(icolwc).eq.0) .or. (ippmod(icolwc).eq.1)) then

  call pdflwc(ncelet, ncel, cvar_fm, cvar_fp2m, cvar_yfm, cvar_yfp2m)

endif

if ((ippmod(icolwc).eq.2) .or. (ippmod(icolwc).eq.3)) then

  call pdfpp3(ncelet, ncel, cvar_fm, cvar_fp2m, cvar_yfm, cvar_yfp2m, cvar_coyfp)

endif

if ((ippmod(icolwc).eq.4).or.(ippmod(icolwc).eq.5)) then

  call pdfpp4(ncelet, ncel, cvar_fm, cvar_fp2m, cvar_yfm, cvar_yfp2m, cvar_coyfp)

endif

!===============================================================================
! 3. CALCUL DE RHO ET DES FRACTIONS MASSIQUES DES ESPECES GLOBALES
!    SUR LES BORDS
!===============================================================================

! --> Masse volumique au bord

call cs_combustion_boundary_conditions_density_ebu_lw()

! --> Fractions massiques des especes globales au bord
do igg = 1, ngazg
  call field_get_val_s(ibym(igg), bsval)
  call field_get_val_s(iym(igg), cpro_ymgg)
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    bsval(ifac) = cpro_ymgg(iel)
  enddo
enddo

!----
! FIN
!----

return
end subroutine
