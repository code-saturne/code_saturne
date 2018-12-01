!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!> \file cs_gas_mix_initialization.f90
!> \brief Initialization of calculation variables for gas mixture modelling
!> in presence of the steam gas or another gas used as variable deduced
!> and not solved.
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

subroutine cs_gas_mix_initialization &
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

character(len=80) :: name_d

double precision volgas, vol_d

type(gas_mix_species_prop) s_d, s_k

double precision, dimension(:), pointer :: cvar_enth, cvar_yk
double precision, dimension(:), pointer :: y_d
double precision, dimension(:), pointer :: cpro_cp
double precision, dimension(:), pointer :: mix_mol_mas

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

iok = 0

!-- Specific heat
if (icp.ge.0) then
  call field_get_val_s(icp, cpro_cp)

!-- Stop if Cp is not variable
else
  call csexit (1)
endif

if (ippmod(icompf).lt.0) then
  ! Enthalpy
  call field_get_val_s(ivarfl(isca(iscalt)), cvar_enth)
endif

!Deduced species (h2o_g) with steam gas
! or Helium or Hydrogen  with noncondensable gases
if (ippmod(igmix).eq.0) then
  name_d = "y_he"
elseif (ippmod(igmix).eq.1) then
  name_d = "y_h2"
elseif (ippmod(igmix).ge.2.and.ippmod(igmix).lt.5) then
  name_d = "y_h2o_g"
else ! ippmod(igmix).eq.5
  name_d = "y_o2"
endif
call field_get_val_s_by_name(name_d, y_d)
call field_get_id(name_d, f_id)
call field_get_key_struct_gas_mix_species_prop(f_id, s_d)

call field_get_val_s(igmxml, mix_mol_mas)

!===============================================================================
! 2. User initialization
!===============================================================================

call cs_user_f_initialization &
( nvar   , nscal  ,                                            &
  dt     )

!===============================================================================
! 3. Deduce the mass fraction (y_d) from the mass fractions (yk) of
!    the noncondensable gases transported
!===============================================================================

if (isuite.eq.0) then
  ! Initialization
  volgas = 0.d0
  vol_d  = 0.d0

  do iel = 1, ncel
    y_d(iel) = 1.d0

    ! Mixture specific Heat
    cpro_cp(iel) = 0.d0

    ! Mixture molar mass
    mix_mol_mas(iel) = 0.d0
  enddo

  do iesp = 1, nscasp
    ! Mass fraction array of the different species
    call field_get_val_s(ivarfl(isca(iscasp(iesp))), cvar_yk)

    call field_get_key_struct_gas_mix_species_prop( &
         ivarfl(isca(iscasp(iesp))), s_k)

    do iel = 1, ncel
      if (cvar_yk(iel).gt.1.d0.or.cvar_yk(iel).lt.0.d0) then
        iok = iok + 1
      endif
      y_d(iel) = y_d(iel)-cvar_yk(iel)

      ! Mixture specific heat (Cp_m0)
      cpro_cp(iel) = cpro_cp(iel) + cvar_yk(iel)*s_k%cp

      mix_mol_mas(iel) = mix_mol_mas(iel) + cvar_yk(iel)/s_k%mol_mas
    enddo
  enddo

  ! Finalization and check
  do iel = 1, ncel
    if (y_d(iel).gt.1.d0.or.y_d(iel).lt.0.d0) then
      iok = iok + 1
    endif
    y_d(iel) = min(max(y_d(iel), 0.d0), 1.d0)

    ! specific heat (Cp_m0) of the gas mixture
    cpro_cp(iel) = cpro_cp(iel) + y_d(iel)*s_d%cp

    if (ippmod(icompf).lt.0) then
      ! Enthalpy initialization
      cvar_enth(iel) = cpro_cp(iel)*t0
    endif

    mix_mol_mas(iel) = mix_mol_mas(iel) + y_d(iel)/s_d%mol_mas
    mix_mol_mas(iel) = 1.d0/mix_mol_mas(iel)

    ! Gas deduced and Total gas volumes injected
    vol_d = vol_d + volume(iel)*    &
            (y_d(iel)/s_d%mol_mas)*mix_mol_mas(iel)
    volgas = volgas +volume(iel)

  enddo

  if (irangp.ge.0) then
   call parsom(volgas)
   call parsom(vol_d)
  endif

  !===============================================================================
  ! 4. Print to the log to check the variables intialization
  !===============================================================================

  write(nfecra, 200)
  write(nfecra, 203) volgas , vol_d

endif

!===============================================================================
! 5. Stop if problem
!===============================================================================

if (iok.gt.0) then
  write(nfecra,3090) iok
  call csexit (1)
endif

!--------
! Formats
!--------

 200 format                                                         &
 (/,                                                                &
 5x,'----------------------------------------------------------' ,/,&
 5x,'**     Gas mixture : Check variables initialization     **' ,/,&
 5x,'----------------------------------------------------------' ,/)

 203 format( &
 3x, '   Total   gas Volume:', 3x, g17.9                         ,/,&
 3x, '   Deduced gas Volume:', 3x, g17.9                         ,/,&
 3x,                                                                &
 3x,'----------------------------------------------------------' )

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
end subroutine cs_gas_mix_initialization
