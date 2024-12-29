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

!===============================================================================
! Function:
! ---------

!> \file d3pphy.f90
!>
!> \brief Specific physic subroutine: diffusion flame.
!>
!> Calculation of mean density
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine d3pphy ()

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use coincl
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use coincl
use ppincl
use radiat
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

integer iel
double precision had

double precision :: phi_t(nvar_turb)
double precision, dimension(:), pointer :: kir
double precision, dimension(:), pointer :: cvar_fm, cvar_fp2m, cvar_scalt
double precision, dimension(:), pointer :: cpro_temp, cpro_rho
double precision, dimension(:), pointer :: cpro_t4m, cpro_tem2
double precision, dimension(:), pointer :: cpro_fuel, cpro_oxyd, cpro_prod

integer       ipass
data          ipass /0/
save          ipass

interface

  subroutine cs_combustion_boundary_conditions_inlet_density()  &
    bind(C, name='cs_combustion_boundary_conditions_inlet_density')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_combustion_boundary_conditions_inlet_density

end interface

!===============================================================================
! 0. ON COMPTE LES PASSAGES
!===============================================================================

ipass = ipass + 1

!===============================================================================
! 1. INITIALISATIONS A CONSERVER
!===============================================================================

!     On reserve la memoire pour le tabeau INDPDF : passage ou non
!     par les PDF (on le garde pendant tout le sous-programme

call field_get_val_s(ivarfl(isca(ifm)), cvar_fm)
call field_get_val_s(ivarfl(isca(ifp2m)), cvar_fp2m)

call field_get_val_s(itemp, cpro_temp)
call field_get_val_s(icrom, cpro_rho)

if (ippmod(icod3p).eq.1) then
  call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)
endif

if ( iirayo.gt.0 ) then
  call field_get_val_s(it4m, cpro_t4m)
endif

call field_get_val_s(it2m, cpro_tem2)
call field_get_val_s_by_name('ym_fuel',  cpro_fuel)
call field_get_val_s_by_name('ym_oxyd',  cpro_oxyd)
call field_get_val_s_by_name('ym_prod',  cpro_prod)

! Calcul de l'enthalpie defect, Kir

allocate (kir(ncelet)); kir(:)=0.d0

if ( iirayo.gt.0 ) then
  do iel =1,ncel
    had = cvar_fm(iel) * hinfue + (1.d0-cvar_fm(iel))*hinoxy
    kir(iel) = max(-(cvar_scalt(iel) - had), 0.d0)
  enddo
endif

do iel = 1, ncel

  call cs_compute_burke_schumann_properties(cvar_fm(iel), cvar_fp2m(iel), kir(iel), phi_t)

  cpro_rho(iel)  = phi_t(4)
  cpro_temp(iel) = phi_t(5)

  cpro_fuel(iel) = phi_t(6)
  cpro_oxyd(iel) = phi_t(7)
  cpro_prod(iel) = phi_t(8)

  cpro_tem2(iel) = phi_t(9)
  cpro_t4m(iel)  = phi_t(10)

enddo

! Free memory
deallocate(kir)

call cs_combustion_boundary_conditions_inlet_density

!----
! END
!----

return
end subroutine
