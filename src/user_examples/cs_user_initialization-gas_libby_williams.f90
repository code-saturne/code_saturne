!-------------------------------------------------------------------------------

!VERS

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

!> \file cs_user_initialization-gas_libby_williams.f90
!>
!> \brief Libby-Williams gas example
!>
!> See \ref cs_user_initialization for examples.
!>
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________


subroutine cs_user_f_initialization &
 ( nvar   , nscal  ,                                              &
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
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          iel, igg, izone
double precision hinit, coefg(ngazgm)
double precision sommqf, sommqt, sommq, tentm, fmelm


character(len=80) :: chaine
integer           :: iscal, ivar, ii
double precision  :: valmax, valmin

double precision, dimension(:), pointer :: cvar_yfm, cvar_fm, cvar_cyfp2m
double precision, dimension(:), pointer :: cvar_fp2m, cvar_coyfp
double precision, dimension(:), pointer :: cvar_scalt, cvar_scal
!< [loc_var_dec]

!===============================================================================

!< [init]
! Variables initialization:
!   ONLY done if there is no restart computation

if (isuite.gt.0) return

! Control output
write(nfecra,9001)

call field_get_val_s(ivarfl(isca(iyfm)), cvar_yfm)
call field_get_val_s(ivarfl(isca(ifm)), cvar_fm)
call field_get_val_s(ivarfl(isca(iyfp2m)), cvar_cyfp2m)
call field_get_val_s(ivarfl(isca(ifp2m)), cvar_fp2m)
call field_get_val_s(ivarfl(isca(icoyfp)), cvar_coyfp)
call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

! a. Preliminary calculations

sommqf = zero
sommq  = zero
sommqt = zero

! Deal with multiple inlets
do izone = 1, nozapm
  sommqf = sommqf + qimp(izone)*fment(izone)
  sommqt = sommqt + qimp(izone)*tkent(izone)
  sommq  = sommq  + qimp(izone)
enddo

if (abs(sommq).gt.epzero) then
  fmelm = sommqf / sommq
  tentm = sommqt / sommq
else
  fmelm = zero
  tentm = t0
endif

! Calculation of the Enthalpy of the gas mixture
! (unburned mean gas)

if (ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3               &
                        .or. ippmod(icolwc).eq.5) then
  coefg(1) = fmelm
  coefg(2) = (1.d0-fmelm)
  coefg(3) = zero

  !   Converting the mean temperatur boundary conditions into
  !   enthalpy values
  hinit = cs_gas_combustion_t_to_h(coefg, tentm)
endif

do iel = 1, ncel

  ! b. Initialisation

  ! Mass fraction of Unburned (fresh) Gas
  cvar_yfm(iel)  = 0.0d0*fmelm
  ! Mean Mixture Fraction
  cvar_fm(iel)   = 0.d0*fmelm
  ! Variance of fuel Mass fraction
  cvar_cyfp2m(iel) = zero
  ! Variance of Mixture Fraction
  cvar_fp2m(iel)  = zero

  ! Covariance for NDIRAC >= 3

  if (ippmod(icolwc).ge. 2) then
    cvar_coyfp(iel)   = zero
  endif

  ! Enthalpy

  if (ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3             &
                          .or. ippmod(icolwc).eq.5) then
    cvar_scalt(iel) = hinit
  endif

enddo

! Control Output of the user defined initialization values

write(nfecra,2000)

do ii  = 1, nscapp
  iscal = iscapp(ii)
  ivar  = isca(iscal)
  call field_get_val_s(ivarfl(isca(ivar)), cvar_scal)
  valmax = -grand
  valmin =  grand
  do iel = 1, ncel
    valmax = max(valmax,cvar_scal(iel))
    valmin = min(valmin,cvar_scal(iel))
  enddo
  if ( irangp.ge.0 ) then
    call parmax(valmax)
    call parmin(valmin)
  endif

  call field_get_label(ivarfl(ivar), chaine)
  write(nfecra,2010)chaine(1:8),valmin,valmax
enddo

write(nfecra,2020)

!--------
! Formats
!--------

 9001 format(                                                   /,&
'  cs_user_initialization: variables initialization by user'   ,/,&
                                                                /)

 2000 format(                                                   /,&
                                                                /,&
' -----------------------------------------------------------' ,/,&
                                                                /,&
' ** INITIALIZATION OF VARIABLES FOR Libby-Williams model'     ,/,&
'    --------------------------------------------------------' ,/,&
'           ONLY ONE PASS'                                     ,/,&
' ---------------------------------'                           ,/,&
'  Variable  Valeur min  Valeur max'                           ,/,&
' ---------------------------------'                             )

 2010 format(                                                     &
 2x,     a8,      e12.4,      e12.4                              )

 2020 format(                                                     &
' ---------------------------------'                           ,/)

!< [init]

!----
! End
!----

return
end subroutine cs_user_f_initialization
