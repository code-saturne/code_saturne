!-------------------------------------------------------------------------------

!VERS

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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
! Purpose:
! -------

!> \file cs_user_physical_properties-compressible_flow.f90
!>
!> \brief Set (variable) physical properties for the compressible flow scheme.
!>
!> This subroutine is called at the beginning of each time step.
!>
!> At the very first time step (not at restart), the only variables that
!> have been initialized are those provided:
!>   - in the GUI and in the user subroutines \ref usipsu and \ref uscfx2; ex.:
!>     - the density             (set to ro0)
!>     - the molecular viscosity (set to viscl0)
!>     - the volumetric molecular viscosity (set to viscv0)
!>     - the molecular thermal conductivity (set to diffusivity_ref of itempk)
!>   - in the user subroutine \ref cs_user_initialization; ex.:
!>     - the unknown variables (null by default)
!>
!> This subroutine allows the user to set the cell values for:
!>   - the molecular viscosity:                            cpro_viscl  kg/(m s)
!>   - the isobaric specific heat
!>   (\f$ C_p = \left. \dfrac{\dd h}{\dd T}\right|_P \f$): cpro_cp     J/(kg degree)
!>   - the molecular thermal conductivity:                 lambda W/(m degree)
!>   - the molecular diffusivity for user-defined scalars: viscls kg/(m s)
!>
!> \section Warnings
!>
!> The density <b> must not </b> be set here: for the compressible scheme,
!> it is one of the unknowns, and it can be initialized as such in the user
!> subroutine \ref cs_user_initialization.
!>
!> The turbulent viscosity <b> must not </b> be modified here (to modify this
!> variable, use the user subroutine \ref usvist)
!>
!> To set a variable isobaric specific heat, the integer \c icp must
!> have been set to 0: the value for \c icp is set automatically in the
!> subroutine \ref cf_set_thermo_options, depending on the thermodynamics laws
!> selected by the user.
!>
!> To set a variable diffusivity for a given user-defined scalar, the
!> variable integer key kivisl (diffusivity_id) must have been set to 0
!> or a field id for the matching field in the user subroutine \ref usipsu or
!> in the GUI (otherwise, a memory problem is expected).
!>
!> Examples are provided in the present subroutine (but they do not have
!> any physical signification).
!>
!> \section cell_id Cells identification
!>
!> Cells may be identified using the \ref getcel subroutine.
!> The syntax of this subroutine is described in the
!> \ref cs_user_boundary_conditions subroutine,
!> but a more thorough description can be found in the user guide.
!>
!> The type of the boundary faces at the previous time step is available
!> (except at the first time step, since the arrays \c itypfb and \c itrifb have
!> not yet been set);
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     mbrom         indicator of filling of romb array
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine usphyv &
 ( nvar   , nscal  ,                                              &
   mbrom  ,                                                       &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use field
use mesh
use cs_cf_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          mbrom

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          ivart, iel
integer          ith, iscal, ifcvsl
double precision varam, varbm, varcm, vardm
double precision varal, varbl, varcl, vardl
double precision varac, varbc
double precision xvart

double precision, dimension(:), pointer :: cpro_viscl, cpro_viscv
double precision, dimension(:), pointer :: cpro_vtmpk, cpro_vscal
double precision, dimension(:), pointer :: cpro_cp, cpro_cv, mix_mol_mas
double precision, dimension(:), pointer :: cvar_scalt
!< [loc_var_dec]

!===============================================================================
! 1. Mandatory initializations
!===============================================================================

!===============================================================================

! Warning: the examples provided below are physically meaningless.
! =======

! These examples must be adapted by the user

! It is adviced to discard all the examples that are not necessary, so
! as to minimize the risk of error.

! List of examples
! ================

! Ex. 1: molecular viscosity varying with temperature
! Ex. 2: molecular volumetric viscosity varying with temperature
! Ex. 3: isobaric specific heat varying with temperature
! Ex. 4: molecular thermal conductivity varying with temperature
! Ex. 5: molecular diffusivity of user-defined scalars varying with temperature

!===============================================================================


!===============================================================================
! Ex. 1: molecular viscosity varying with temperature
! =====
!    The values of the molecular viscosity are provided as a function of
!    the temperature. All variables are evaluated at the cell centers.
!===============================================================================

!< [example_1]
!     To refer to the user-defined scalar number 2 instead, for example, use
!     ivart = isca(2)

ivart = isca(itempk)
call field_get_val_s(ivarfl(ivart), cvar_scalt)

! --- Molecular dynamic viscosity 'cpro_viscl'
!     (physical properties at the cell centers)

call field_get_val_s(iviscl, cpro_viscl)

! --- User-defined coefficients for the selected law.
!     The values hereafter are provided as a mere example. They
!     are physically meaningless.

varam = -3.4016d-9
varbm =  6.2332d-7
varcm = -4.5577d-5
vardm =  1.6935d-3

! --- Molecular dynamic viscosity mu at the cell centers, kg/(m s)
!     In this example, mu is provided as a function of the temperature T:
!       mu(T)              =    T  *( T  *( am  * T +  bm  )+ cm  )+ dm
!     that is:
!       cpro_viscl(iel) =   xvart*(xvart*(varam*xvart+varbm)+varcm)+vardm

do iel = 1, ncel
  xvart = cvar_scalt(iel)
  cpro_viscl(iel) = xvart*(xvart*(varam*xvart+varbm)+varcm)+vardm
enddo
!< [example_1]

!===============================================================================
! Ex. 2: molecular volumetric viscosity varying with temperature
! =====
!    The values of the molecular volumetric viscosity are provided as a function
!    of the temperature. All variables are evaluated at the cell centers.
!===============================================================================

!     To refer to the user-defined scalar number 2 instead, for example, use
!     ivart = isca(2)

!< [example_2]
ivart = isca(itempk)
call field_get_val_s(ivarfl(ivart), cvar_scalt)

! --- Molecular volumetric viscosity

if (iviscv.ge.0) then
  call field_get_val_s(iviscv, cpro_viscv)
else
  cpro_viscv => NULL()
endif

! --- Stop if the volumetric viscosity has not been defined as variable

if (iviscv.lt.0) then
  write(nfecra,2000) iviscv
  call csexit (1)
endif

! --- User-defined coefficients for the selected law.
!     The values provided hereafter are provided as a mere example. They
!     are physically meaningless.

varam = -3.4016d-9
varbm =  6.2332d-7
varcm = -4.5577d-5
vardm =  1.6935d-3

! --- Molecular dynamic volumetric viscosity kappa at the cell centers, kg/(m s)
!     In this example, kappa is provided as a function of the temperature T:
!       kappa(T)           =    T  *( T  *( am  * T +  bm  )+ cm  )+ dm
!     that is:
!       cpro_viscv(iel) =   xvart*(xvart*(varam*xvart+varbm)+varcm)+vardm

do iel = 1, ncel
  xvart = cvar_scalt(iel)
  cpro_viscv(iel) = xvart*(xvart*(varam*xvart+varbm)+varcm)+vardm
enddo
!< [example_2]

!===============================================================================
! Ex. 3: isobaric specific heat varying with temperature
! =====
!    The values of the isobaric specific heat values are provided as a function
!    of the temperature. All variables are evaluated at the cell centers.
!===============================================================================

!     To refer to the user-defined scalar number 2 instead, for example, use
!     ivart = isca(2)

!< [example_3]
ivart = isca(itempk)
call field_get_val_s(ivarfl(ivart), cvar_scalt)

! --- Isobaric specific heat

if (icp.ge.0) call field_get_val_s(icp, cpro_cp)

! --- Stop if the isobaric or isochoric specific heat (cpro_cp or cpro_cv)
!     has not been defined as variable

if (icp.lt.0) then
  write(nfecra,1000) icp
  call csexit (1)
endif
if (icv.lt.0) then
  write(nfecra,1001) icv
  call csexit (1)
endif

! --- User-defined coefficients for the selected law.
!     The values provided hereafter are provided as a mere example. They
!     are physically meaningless.

varac = 0.00001d0
varbc = 1000.0d0

! --- Isobaric specific heat cpro_cp at the cell centers, J/(kg degree)
!     In this example, cpro_cp is provided as a function of the temperature T:
!       cpro_cp(T)              =      ac * T  + ab
!     that is:
!       cpro_cp(iel)            =    varac*xvart+varbc

do iel = 1, ncel
  xvart = cvar_scalt(iel)
  cpro_cp(iel) = varac*xvart + varbc
enddo

! --- The isochoric specific heat is deduced from the isobaric specific heat

call field_get_val_s(icv, cpro_cv)
call field_get_val_s(igmxml, mix_mol_mas)
call cs_cf_thermo_cv(cpro_cp, mix_mol_mas, cpro_cv, ncel)

!< [example_3]

!===============================================================================
! Ex. 4: molecular thermal conductivity varying with temperature
! =====
!    The values of the molecular thermal conductivity are provided as a function
!    of the temperature. All variables are evaluated at the cell centers.
!===============================================================================

!     To refer to the user-defined scalar number 2 instead, for example, use
!     ivart = isca(2)

!< [example_4]
ivart = isca(itempk)
call field_get_val_s(ivarfl(ivart), cvar_scalt)

! --- Molecular thermal conductivity

call field_get_key_int(ivarfl(isca(itempk)), kivisl, ifcvsl)
if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, cpro_vtmpk)
else
  cpro_vtmpk => NULL()
endif

! --- Stop if the molecular thermal conductivity has not
!     been defined as variable

if (ifcvsl.lt.0) then
  write(nfecra,1010) itempk
  call csexit (1)
endif

! --- User-defined coefficients for the selected law.
!     The values provided hereafter are provided as a mere example. They
!     are physically meaningless.

varal = -3.3283d-7
varbl =  3.6021d-5
varcl =  1.2527d-4
vardl =  0.58923d0

! --- Molecular thermal conductivity lambda at the cell centers, W/(m degree)
!     In this example, lambda is provided as a function of the temperature T:
!       lambda(T)          =    T  *( T  *( al  * T +  bl  )+ cl  )+ dl
!     that is:
!       cpro_vtmpk(iel) =   xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl

do iel = 1, ncel
  xvart = cvar_scalt(iel)
  cpro_vtmpk(iel) = (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)
enddo
!< [example_4]

!===============================================================================
! Ex. 5: molecular diffusivity of user-defined scalars varying with temperature
! =====
!    The molecular diffusivity can be set for all the user-defined scalars
!    ** except **:
!      - temperature and enthalpy (already dealt with above: for these
!        variables, the 'diffusivity' is the thermal conductivity)
!      - variances of the fluctuations of another scalar variable (the
!        diffusivity is assumed to be equal to that of the associated
!        scalar)
!    The values of the molecular diffusivity are provided as a function
!    of the temperature. All variables are evaluated at the cell centers.
!===============================================================================

!< [example_5]
! --- Loop on the scalars
do iscal = 1, nscaus

  ! --- If the scalar is the temperature, it is marked by ith = 1
  !     so that it will be skipped.

  ith = 0
  if (iscal.eq.itempk) ith = 1

  ! --- If the variable represents the variance of the fluctuations of
  !     another scalar variable (iscavr <= 0), it is simply skipped.

  if (ith.eq.0.and.iscavr(iscal).le.0) then

    ! --- Here, iscal points to any scalar variable except the temperature,
    !     the enthalpy and the variance of the fluctuations of another
    !     scalar variable.

    ivart = isca(itempk)
    call field_get_val_s(ivarfl(ivart), cvar_scalt)

    ! --- Molecular diffusivity of the current scalar iscal

    call field_get_key_int(ivarfl(isca(iscal)), kivisl, ifcvsl)
    if (ifcvsl.ge.0) then
      call field_get_val_s(ifcvsl, cpro_vscal)
    else
      cpro_vscal => NULL()
    endif

    ! --- Stop if the molecular diffusivity has not been defined as variable

    if (ifcvsl.lt.0) then
      write(nfecra,1010) iscal
      call csexit (1)
    endif

    ! --- User-defined coefficients for the selected law.
    !     The values provided hereafter are provided as a mere example. They
    !     are physically meaningless.

    varal = -3.3283d-7
    varbl =  3.6021d-5
    varcl =  1.2527d-4
    vardl =  0.58923d0

    ! --- Molecular diffusivity lambda at the cell centers, kg/(m s)
    !     In this example, lambda is provided as a function of the temperature T:
    !       lambda(T)          =    T  *( T  *( al  * T +  bl  )+ cl  )+ dl
    !     that is:
    !       cpro_vscal(iel) =   xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl

    do iel = 1, ncel
      xvart = cvar_scalt(iel)
      cpro_vscal(iel) = (xvart*(xvart*(varal*xvart+varbl)+varcl)+vardl)
    enddo


  endif
  ! --- End of the tests on ith and iscavr

enddo
! --- End of the loop on the scalars
!< [example_5]

!--------
! Formats
!--------

 1000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@       in the GUI or in the user subroutine ''usipsu'', the',/, &
'@         isobaric specific heat is declared as a property',/,   &
'@         uniform in space: icp = ',i10   ,/,                    &
'@       in the user subroutine ''usphyv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the GUI input data or the',/, &
'@    user subroutines ''usipsu'' or ''usphyv''.',/,              &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 1001 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@       in the GUI or in the user subroutine ''usipsu'', the',/, &
'@         isochoric specific heat is declared as a property',/,  &
'@         uniform in space: icv = ',i10   ,/,                    &
'@       in the user subroutine ''usphyv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@@',/,                                                           &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the GUI input data or the',/, &
'@    user subroutines ''usipsu'' or ''usphyv''.',/,              &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 1010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@     For the scalar ',i10,/,                                    &
'@       in the GUI or in the user subroutine ''usipsu'', the',/, &
'@         molecular diffusivity is declared as a property',/,    &
'@         uniform in space.',/,                                  &
'@       in the user subroutine ''usphyv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@@',/,                                                           &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the GUI input data or the',/, &
'@    user subroutines ''usipsu'' or ''usphyv''.',/,              &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 2000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@       in the user subroutine ''uscfx2'', the molecular',/,     &
'@         volumetric viscosity is declared as a property',/,     &
'@         uniform in space: iviscv = ',i10,/,                    &
'@       in the user subroutine ''usphyv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@@',/,                                                           &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the user subroutines',/,      &
'@    ''uscfx2'' or ''usphyv''.',/,                               &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

!----
! End
!----

return
end subroutine usphyv
