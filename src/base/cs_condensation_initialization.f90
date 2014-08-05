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

!> \file cs_condensation_initialization.f90
!> \brief Initialization of calculation variables for condensation modelling
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step
!______________________________________________________________________________

subroutine cs_condensation_initialization &
 ( nvar   , nscal  , dt )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use pointe
use entsor
use parall
use period
use ihmpre
use ppppar
use ppthch
use ppincl
use mesh
use field
use cfpoin, only:ithvar
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

integer          iel, iesp
integer          iok
integer          f_id

double precision valmax, valmin, vfmin , vfmax
double precision volgas, volvap

type(severe_acc_species_prop) s_h2o_g, s_k

double precision, allocatable, dimension(:) :: mix_mol_mas
double precision, dimension(:), pointer :: cvar_enth, cvar_yk
double precision, dimension(:), pointer :: y_h2o_g
double precision, dimension(:), pointer :: cpro_rho, cpro_viscl, cpro_cp, cpro_venth

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================


iok = 0

!-- Specific heat
if (icp.gt.0) then
  call field_get_val_s(iprpfl(icp), cpro_cp)

!-- Stop if Cp is not variable
else
  call csexit (1)
endif

! Enthaltpy
call field_get_val_s(ivarfl(isca(iscalt)), cvar_enth)

! Condensable gas H20 (vapor)
call field_get_val_s_by_name("y_h2o_g", y_h2o_g)
call field_get_id("y_h2o_g", f_id)
call field_get_key_struct_severe_acc_species_prop(f_id, s_h2o_g)

!===============================================================================
! 2. User initialization
!===============================================================================

call cs_user_initialization &
( nvar   , nscal  ,                                            &
  dt     )

!===============================================================================
! 3. Deduce the mass fraction of vapor H2O from the mass fraction of
!    non-condensable gases
!===============================================================================

allocate(mix_mol_mas(ncelet))

! Initialization
volgas = 0.d0
volvap = 0.d0

do iel = 1, ncel
  y_h2o_g(iel) = 1.d0

  ! Specific Heat of the mix
  cpro_cp(iel) = 0.d0

  ! Mixing molar mass
  mix_mol_mas(iel) = 0.d0
enddo

do iesp = 1, nscasp
  ! mass fraction array of the different
  ! species (O2, N2)
  !-------------------------------------------------
  call field_get_val_s(ivarfl(isca(iscasp(iesp))), cvar_yk)

  call field_get_key_struct_severe_acc_species_prop( &
       ivarfl(isca(iscasp(iesp))), s_k)

  do iel = 1, ncel
    if (cvar_yk(iel).gt.1.d0.or.cvar_yk(iel).lt.0.d0) then
      iok = iok + 1
    endif
    y_h2o_g(iel) = y_h2o_g(iel)-cvar_yk(iel)

    ! Mixing specific heat (Cp_m0)
    cpro_cp(iel) = cpro_cp(iel) + cvar_yk(iel)*s_k%cp

    mix_mol_mas(iel) = mix_mol_mas(iel) + cvar_yk(iel)/s_k%mol_mas
  enddo
enddo

! Finalization and check
do iel = 1, ncel
  if (y_h2o_g(iel).gt.1.d0.or.y_h2o_g(iel).lt.0.d0) then
    iok = iok + 1
  endif
  y_h2o_g(iel) = min(max(y_h2o_g(iel), 0.d0), 1.d0)

  ! Mixing specific heat (Cp_m0)
  cpro_cp(iel) = cpro_cp(iel) + y_h2o_g(iel)*s_h2o_g%cp

  !-- enthalpy initialization
  cvar_enth(iel) = cpro_cp(iel)*t0

  mix_mol_mas(iel) = mix_mol_mas(iel) + y_h2o_g(iel)/s_h2o_g%mol_mas
  mix_mol_mas(iel) = 1.d0/mix_mol_mas(iel)

  !-- Helium volume and Total gas volume injected
  volvap = volvap + volume(iel)*    &
          (y_h2o_g(iel)/18.d0) * (mix_mol_mas(iel)*1000.d0)
  volgas = volgas +volume(iel)

enddo

if (irangp.ge.0) then
 call parsom(volgas)
 call parsom(volvap)
endif


! Deallocate the temporary array
deallocate(mix_mol_mas)

!===============================================================================
! 4. Print to the listing to Check the variables intialization
!===============================================================================
write(nfecra, 200)
write(nfecra, 203) volgas , volvap

!===============================================================================
! 5. Stop if problem
!===============================================================================

if (iok.gt.0) then
  write(nfecra,3090) iok
  call csexit (1)
endif

write(nfecra,3000)

!--------
! Formats
!--------

 200 format                                                   &
 (/,                                                          &
 5x,'-------------------------------------------------' ,/,&
 5x,'** Condensation: Check variables initialization**' ,/,&
 5x,'-------------------------------------------------' ,/)

 203 format( &
 3x, '   Total gas Volume:', 3x, g17.9                     ,/,&
 3x, '   Steam gas Volume:', 3x, g17.9                     ,/,&
 3x,                                                          &
 '========================================================' )

 3000 format(/,/,                                                 &
'-------------------------------------------------------------',/)
 3090 format(                                                     &
'@',                                                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@',                                                            /,&
'@    THE VARIABLES INITIALIZATION IS INCOMPLETE OR',           /,&
'@    INCOHERENT WITH THE PARAMETERS VALUE OF THE CALCULATION', /,&
'@',                                                            /,&
'@  The calculation will not be run (',i10,' errors).',         /,&
'@',                                                            /,&
'@  Refer to the previous warnings for further information.',   /,&
'@  Pay attention to the initialization of',                    /,&
'@                                the time-step',               /,&
'@                                the turbulence',              /,&
'@                                the scalars and variances',   /,&
'@                                the time averages',           /,&
'@',                                                            /,&
'@  Verify the initialization and the restart file.',           /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

!----
! End
!----

return
end subroutine
