!-------------------------------------------------------------------------------

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

subroutine cs_steady_laminar_flamelet_init

!===============================================================================
! FONCTION :
! --------

! INITIALISATION DES VARIABLES DE CALCUL
!    POUR LA PHYSIQUE PARTICULIERE :
!        COMBUSTION GAZ - FLAMME DE DIFFUSION Steady laminar flamelet
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
use cpincl
use ppincl
use mesh
use field

!===============================================================================

implicit none

! Local variables

integer          iel

double precision, dimension(:), pointer :: cvar_scalt
double precision, dimension(:), pointer :: cvar_fm, cvar_fp2m, cvar_fsqm
double precision, dimension(:), pointer :: cvar_prog, cpro_prog
double precision, dimension(:), pointer :: cvar_npm, cvar_fsm

!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================

call field_get_val_s(ivarfl(isca(ifm)), cvar_fm)

if (mode_fp2m.eq.0) then
  call field_get_val_s(ivarfl(isca(ifp2m)), cvar_fp2m)
else if(mode_fp2m.eq.1) then
  call field_get_val_s(ivarfl(isca(ifsqm)), cvar_fsqm)
endif

if ( ippmod(islfm).ge.2 ) then
  call field_get_val_s(ivarfl(isca(ipvm)), cvar_prog)
  call field_get_val_s(iym(ngazgm), cpro_prog)
endif

if (ippmod(islfm).eq.1 .or. ippmod(islfm).eq.3) then
  call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)
endif

if (isoot.eq.1) then
  call field_get_val_s(ivarfl(isca(inpm)), cvar_npm)
  call field_get_val_s(ivarfl(isca(ifsm)), cvar_fsm)
endif

!===============================================================================
! 2. INITIALISATION DES INCONNUES :
!      UNIQUEMENT SI ON NE FAIT PAS UNE SUITE
!===============================================================================

if ( isuite.eq.0 ) then

  do iel = 1, ncel

    ! Mean mixture fraction and its variance
    cvar_fm(iel)   = zero

    if (mode_fp2m.eq.0) then
      cvar_fp2m(iel) = zero
    else if(mode_fp2m.eq.1) then
      cvar_fsqm(iel) = zero
    endif

    ! Progress variable
    if ( ippmod(islfm).ge.2 ) then
      cvar_prog(iel) = cpro_prog(iel)
    endif

    ! Enthalpy
    if ( ippmod(islfm).eq.1 .or. ippmod(islfm).eq.3 ) then
      cvar_scalt(iel) = hinoxy
    endif

    ! Soot
    if (isoot.eq.1) then
      cvar_npm(iel) = 0.d0
      cvar_fsm(iel) = 0.d0
    endif

  enddo

  ! ---> Parallelism and periodic exchange
  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(cvar_fm)
    if (mode_fp2m.eq.0) call synsca(cvar_fp2m)
    if (mode_fp2m.eq.1) call synsca(cvar_fsqm)
    if ( ippmod(islfm).ge.2 ) then
      call synsca(cvar_prog)
    endif
    if ( ippmod(islfm).eq.1 .or. ippmod(islfm).eq.3 ) then
      call synsca(cvar_scalt)
    endif
    if (isoot.eq.1) then
      call synsca(cvar_npm)
      call synsca(cvar_fsm)
    endif
  endif

endif


!----
! FORMATS
!----

!----
! Fin
!----

return
end subroutine
