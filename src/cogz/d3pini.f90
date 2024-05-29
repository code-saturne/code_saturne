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

subroutine d3pini

!===============================================================================
! FONCTION :
! --------

! INITIALISATION DES VARIABLES DE CALCUL
!    POUR LA PHYSIQUE PARTICULIERE :
!        COMBUSTION GAZ - FLAMME DE DIFFUSION CHIMIE 3 POINTS
!    PENDANT DE USINIV.F

! Cette routine est appelee en debut de calcul (suite ou non)
!     avant le debut de la boucle en temps

! Elle permet d'INITIALISER ou de MODIFIER (pour les calculs suite)
!     les variables de calcul,
!     les valeurs du pas de temps


! On dispose ici de ROM et VISCL initialises par RO0 et VISCL0
!     ou relues d'un fichier suite
! On ne dispose des variables VISCLS, CP (quand elles sont
!     definies) que si elles ont pu etre relues dans un fichier
!     suite de calcul

! LA MODIFICATION DES PROPRIETES PHYSIQUES (ROM, VISCL, VISCLS, CP)
!     SE FERA EN STANDARD DANS LE SOUS PROGRAMME PPPHYV
!     ET PAS ICI

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

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
use parall
use period
use ppppar
use ppthch
use coincl
use ppincl
use mesh
use field
use cs_c_bindings
use iso_c_binding

!===============================================================================

implicit none

! Local variables

integer          ii, iel, igg
double precision coefg(ngazgm), hair
double precision ye(ngaze), yg(ngazg), ytot

double precision, dimension(:), pointer :: cvar_scalt
double precision, dimension(:), pointer :: cvar_fm, cvar_fp2m
double precision, dimension(:), pointer :: cvar_npm, cvar_fsm

!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================

call field_get_val_s(ivarfl(isca(ifm)), cvar_fm)
call field_get_val_s(ivarfl(isca(ifp2m)), cvar_fp2m)

if (ippmod(icod3p).eq.1) then
  call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)
endif

if (isoot.ge.1) then
  call field_get_val_s(ivarfl(isca(inpm)), cvar_npm)
  call field_get_val_s(ivarfl(isca(ifsm)), cvar_fsm)
endif

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

!===============================================================================
! 2. INITIALISATION DES INCONNUES :
!      UNIQUEMENT SI ON NE FAIT PAS UNE SUITE
!===============================================================================

if ( isuite.eq.0 ) then

  ! ---> Initialization a with air at T0

  yg(1) = 0.d0
  yg(2) = 1.d0
  yg(3) = 0.d0

  do ii = 1, ngaze-1
    ye(ii) = 0.d0
    do igg = 1 ,ngazg
        ye(ii) = ye(ii) + coefeg(ii,igg) * yg(igg)
    enddo
  enddo

  ytot = 0.d0
  do ii = 1, ngaze-1
    ytot = ytot + ye(ii)
  enddo
  ye(ngaze) = 1.d0 - ytot

  hair = cs_compute_burke_schumann_enthalpy(t0, ye)

  do iel = 1, ncel

    ! Mean mixture fraction and its variance
    cvar_fm(iel)   = zero
    cvar_fp2m(iel) = zero

    ! Enthalpy
    if ( ippmod(icod3p).eq.1 ) then
      cvar_scalt(iel) = hair
    endif

    ! Soot
    if (isoot.ge.1) then
      cvar_npm(iel) = 0.d0
      cvar_fsm(iel) = 0.d0
    endif

  enddo

  ! ---> Parallelism and periodic exchange
  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(cvar_fm)
    call synsca(cvar_fp2m)
    if ( ippmod(icod3p).eq.1 ) then
      call synsca(cvar_scalt)
    endif
    if (isoot.ge.1) then
      call synsca(cvar_npm)
      call synsca(cvar_fsm)
    endif
  endif

endif

!----
! FORMATS
!----

!----
! End
!----

return
end subroutine
