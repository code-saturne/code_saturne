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

subroutine cothht &
 ( mode   , nespec , nespem , xespec ,                            &
   npo    , npot   , th     , eh     ,                            &
   enthal , temper )

!===============================================================================
!  FONCTION  :
!  --------

! CETTE FONCTION CALCULE L'ENTHALPIE A PARTIR DE LA
!  COMPOSITION ET DE LA VALEUR DE LA TEMPERATURE
!  SPECIFIQUE A LA COMBUSTION GAZ

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! mode             ! e  ! <-- !  -1 : t -> h  ;   1 : h -> t                   !
! nespec           ! e  ! <-- ! nb de constituants                             !
! nespem           ! e  ! <-- ! nb maximal de constituants                     !
! xespec           ! tr ! <-- ! fraction massique des constituants             !
! npo              ! e  ! <-- ! nombre de ponits de tabulation                 !
! npot             ! e  ! <-- ! nombre maximal de ponits                       !
!                  !    !     !   de tabulation                                !
! th               ! tr ! <-- ! tabulation temperature en kelvin               !
! eh               ! tr ! <-- ! tabulation enthalpie - temperature             !
! enthal           ! r  ! <-- ! enthalpie massique j/kg                        !
! temper           ! r  ! <-- ! temperature en kelvin                          !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
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
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          mode , npo , npot , nespec , nespem

double precision xespec(nespem)
double precision th(npot) , eh(nespem,npot)
double precision temper , enthal

!===============================================================================
! 1. CALCUL DE L'ENTHALPIE A PARTIR DE LA TEMPERATURE
!===============================================================================

if (mode.eq.-1) then

  enthal = cs_gas_combustion_t_to_h(xespec, temper)

!===============================================================================
! 2. CALCUL DE LA TEMPERATURE A PARTIR DE l'ENTHALPIE
!===============================================================================

else if (mode.eq.1) then

  temper = cs_gas_combustion_h_to_t(xespec, enthal)

else

  write(nfecra,1000) mode
  call csexit (1)

endif

!--------
! Formats
!--------

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR DANS COTHHT                          ',/,&
'@    =========                                               ',/,&
'@    VALEUR INCORRECTE DE L''ARGUMENT MODE                   ',/,&
'@    CE DOIT ETRE UN ENTIER EGAL A 1 OU -1                   ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

return
end subroutine

!===============================================================================
! Function :
! --------

!> \brief Convert enthalpy to temperature at cells for gas combustion.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     xespec       mass fraction of constituents
!> \param[in]     enthal       enthalpy at cells
!> \param[out]    temper       temperature at cells
!_______________________________________________________________________________


function cs_gas_combustion_h_to_t(xespec, enthal) result(temper)  &
   bind(C, name='cs_gas_combustion_h_to_t')

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding
use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppthch
use coincl

!===============================================================================

implicit none

! Arguments

real(kind=c_double), dimension(*) :: xespec
real(c_double), value :: enthal
real(c_double) :: temper

! Local variables

integer          it , iesp

double precision eh1 , eh0

it  = npo-1
eh1 = zero
do iesp = 1, ngazg
  eh1 = eh1 + xespec(iesp)*ehgazg(iesp,it+1)
enddo
if (enthal.ge.eh1) temper = th(it+1)

it  = 1
eh0 = zero
do iesp = 1, ngazg
  eh0 = eh0 + xespec(iesp)*ehgazg(iesp,it  )
enddo
if (enthal.le.eh0) temper = th(it)

do it = 1, npo-1
  eh0 = zero
  eh1 = zero
  do iesp = 1, ngazg
    eh0 = eh0 + xespec(iesp)*ehgazg(iesp,it  )
    eh1 = eh1 + xespec(iesp)*ehgazg(iesp,it+1)
  enddo
  if (enthal.ge.eh0 .and. enthal.le.eh1)                        &
    temper = th(it) + (enthal-eh0)*(th(it+1)-th(it))/(eh1-eh0)
enddo

!----
! End
!----

return
end function cs_gas_combustion_h_to_t

!===============================================================================
! Function :
! --------

!> \brief Convert a temperature value to enthalpy for gas combustion.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     xespec       mass fraction of constituants
!> \param[in]     temper       temperature at cells
!> \param[out]    enthal       enthalpy at cells
!_______________________________________________________________________________


function cs_gas_combustion_t_to_h(xespec, temper) result(enthal)  &
   bind(C, name='cs_gas_combustion_t_to_h')

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding
use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppthch
use coincl

!===============================================================================

implicit none

! Arguments

real(kind=c_double), dimension(*) :: xespec
real(c_double), value :: temper
real(c_double) :: enthal

! Local variables

integer          it , iesp

double precision eh1 , eh0

!===============================================================================
! 1. CALCUL DE L'ENTHALPIE A PARTIR DE LA TEMPERATURE
!===============================================================================

it = npo
if (temper.ge.th(it)) then
  enthal = zero
  do iesp = 1, ngazg
    enthal = enthal + xespec(iesp)*ehgazg(iesp,it)
  enddo
  go to 11
endif

it = 1
if (temper.le.th(it)) then
  enthal = zero
  do iesp = 1, ngazg
    enthal = enthal + xespec(iesp)*ehgazg(iesp,it)
  enddo
  go to 11
endif
10 continue

it = it + 1
if (temper.le.th(it)) then
  eh0 = zero
  eh1 = zero
  do iesp = 1, ngazg
    eh0 = eh0 + xespec(iesp)*ehgazg(iesp,it-1)
    eh1 = eh1 + xespec(iesp)*ehgazg(iesp,it  )
  enddo
  enthal = eh0 + (eh1-eh0)*(temper-th(it-1))/(th(it)-th(it-1))
  goto 11
endif
goto 10
11 continue

!----
! End
!----

return
end function cs_gas_combustion_t_to_h

!===============================================================================
! Function :
! --------

!> \brief Convert enthalpy to temperature at boundary for gas combustion

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     h_b           enthalpy at boundary
!> \param[in,out] t_b           temperature at boundary
!_______________________________________________________________________________


subroutine coh2tb(h_b, t_b)

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use pointe
use ppppar
use ppthch
use coincl
use ppincl
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

double precision, dimension(nfabor), intent(in) :: h_b
double precision, dimension(nfabor), intent(out), target :: t_b

! Local variables

integer          iel , ifac, izone
integer          igg

double precision coefg(ngazgm)
double precision hbl

double precision, dimension(:), pointer :: bym1, bym2, bym3, cpro_temp

!===============================================================================

! Non-specific physics

if (ippmod(icoebu).ge.0 .or. ippmod(icod3p).ge.0) then

  call field_get_val_s(ibym(1), bym1)
  call field_get_val_s(ibym(2), bym2)
  call field_get_val_s(ibym(3), bym3)

  do ifac = 1, nfabor
    iel = ifabor(ifac)
    hbl = h_b(ifac)

    do igg = 1, ngazgm
      coefg(igg) = zero
    enddo
    coefg(1) = bym1(ifac)
    coefg(2) = bym2(ifac)
    coefg(3) = bym3(ifac)
    t_b(ifac) = cs_gas_combustion_h_to_t(coefg, hbl)
  enddo

elseif (ippmod(islfm).ge.0) then

  call field_get_val_s(itemp, cpro_temp)

  do ifac = 1, nfabor
    izone = izfppp(ifac)

    if (ientfu(izone).eq.1) then
      t_b(ifac) = tinfue
    elseif (ientox(izone).eq.1) then
      t_b(ifac) = tinoxy
    else
      iel = ifabor(ifac)
      t_b(ifac) = cpro_temp(iel)
    endif
  enddo

endif

return
end subroutine coh2tb

!===============================================================================
! Function :
! --------

!> \brief Convert temperature to enthalpy at boundary for gas combustion.

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     n_faces       number of faces in list
!> \param[in]     face_ids      list of boundary faces at which conversion
!>                              is requested (0-based numbering)
!> \param[in]     t_b           temperature at boundary
!> \param[out]    h_b           enthalpy at boundary
!_______________________________________________________________________________


subroutine cot2hb(n_faces, face_ids, t_b, h_b)

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use pointe
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

! Arguments

integer          n_faces
integer          face_ids(n_faces)

double precision, dimension(nfabor),intent(in) :: t_b
double precision, dimension(nfabor), intent(out), target :: h_b

! Local variables

integer          iel , ilst, ifac, igg, izone

double precision coefg(ngazgm)
double precision tbl

double precision, dimension(:), pointer :: bym1, bym2, bym3, cvar_h

!===============================================================================

if (ippmod(icoebu).ge.0 .or. ippmod(icod3p).ge.0) then

  ! Mappings for gas combustion
  call field_get_val_s(ibym(1), bym1)
  call field_get_val_s(ibym(2), bym2)
  call field_get_val_s(ibym(3), bym3)

  ! Now loop on faces

  do ilst = 1, n_faces

    ifac = face_ids(ilst) + 1
    iel = ifabor(ifac)

    tbl = t_b(ifac)

    do igg = 1, ngazgm
      coefg(igg) = zero
    enddo
    coefg(1) = bym1(ifac)
    coefg(2) = bym2(ifac)
    coefg(3) = bym3(ifac)
    h_b(ifac) = cs_gas_combustion_t_to_h(coefg, tbl)

  enddo

elseif (ippmod(islfm).ge.0) then

  ! Mappings for gas combustion
  call field_get_val_s(ivarfl(isca(iscalt)), cvar_h)

  ! Now loop on faces

  do ilst = 1, n_faces

    ifac = face_ids(ilst) + 1
    izone = izfppp(ifac)

    if (ientfu(izone).eq.1) then
      h_b(ifac) = hinfue
    elseif (ientox(izone).eq.1) then
      h_b(ifac) = hinoxy
    else
      iel = ifabor(ifac)
      h_b(ifac) = cvar_h(iel)
    endif
  enddo

endif

return
end subroutine cot2hb
