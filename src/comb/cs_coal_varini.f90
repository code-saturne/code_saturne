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

subroutine cs_coal_varini

!===============================================================================
! FONCTION :
! --------

! INITIALISATION DES VARIABLES DE CALCUL
!    POUR LA PHYSIQUE PARTICULIERE : COMBUSTION CP
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
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use mesh
use field

!===============================================================================

implicit none

! Local variables

integer          iel, ige, icha, ifac

double precision t1init, h1init, coefe(ngazem)
double precision t2init
double precision f1mc(ncharm), f2mc(ncharm)

integer ioxy
double precision wmh2o,wmco2,wmn2,wmo2,dmas

double precision, dimension(:), pointer :: cvar_scalt, cvar_xch
double precision, dimension(:), pointer :: cvar_yco2, cvar_yhcn, cvar_ynh3
double precision, dimension(:), pointer :: cvar_yno, cvar_hox
double precision, dimension(:), pointer :: cpro_x1, bpro_x1

!===============================================================================
! Interfaces
!===============================================================================

interface

  function cs_coal_ht_convert_t_to_h_gas_by_yi   &
    (tp, xesp, f1mc, f2mc)  result(eh) &
    bind(C, name='cs_coal_ht_convert_t_to_h_gas_by_yi')
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double), value :: tp
    real(c_double), dimension(*) :: xesp
    real(c_double), dimension(*) :: f1mc, f2mc
    real(c_double) :: eh
  end function cs_coal_ht_convert_t_to_h_gas_by_yi

end interface

!===============================================================================
! 1. Initializations
!===============================================================================

! Massic fraction of gas
call field_get_val_s_by_name("x_c", cpro_x1)
call field_get_val_s_by_name("b_x_c", bpro_x1)

call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)
call field_get_val_s_by_name("x_c_h", cvar_xch)

if (ieqco2.ge.1) then
  call field_get_val_s(iyco2, cvar_yco2)
endif
if (ieqnox.eq.1) then
  call field_get_val_s(iyhcn, cvar_yhcn)
  call field_get_val_s(iynh3, cvar_ynh3)
  call field_get_val_s(iyno, cvar_yno)
  call field_get_val_s(ihox, cvar_hox)
endif

!===============================================================================
! 2. Variable initialization
!===============================================================================

if (isuite.eq.0) then

  ! --> All the domain is filled with the first oxidizer at TINITK
  !                   ============================================

  ! ---- Computation of H1INIT and H2INIT

  t1init = t0
  t2init = t0

  ! ------ Transported variables for the mix (solid+carrying gas)^2

  do ige = 1, ngazem
    coefe(ige) = zero
  enddo

  ! Oxidizer are mix of O2, N2 (air), CO2 and H2O (recycled exhaust)
  ! the composition of the fisrt oxidiser is taken in account
  coefe(io2) = wmole(io2)*oxyo2(1)                                &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(ih2o) = wmole(ih2o)*oxyh2o(1)                             &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(ico2) = wmole(ico2)*oxyco2(1)                             &
              /( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)       &
                +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))
  coefe(in2) = 1.d0-coefe(io2)-coefe(ih2o)-coefe(ico2)

  do icha = 1, ncharm
    f1mc(icha) = zero
    f2mc(icha) = zero
  enddo

  h1init = cs_coal_ht_convert_t_to_h_gas_by_yi(t1init, coefe, f1mc, f2mc)

  do iel = 1, ncel
    cvar_scalt(iel) = h1init
    cvar_xch(iel)   = h1init
  enddo

  ! Transported variables for the mix (passive scalars, variance)

  if (ieqco2.ge.1) then
    ioxy   = 1
    wmo2   = wmole(io2)
    wmco2  = wmole(ico2)
    wmh2o  = wmole(ih2o)
    wmn2   = wmole(in2)
    dmas = (  oxyo2 (ioxy)*wmo2 +oxyn2 (ioxy)*wmn2              &
            + oxyh2o(ioxy)*wmh2o+oxyco2(ioxy)*wmco2)
    xco2 = oxyco2(ioxy)*wmco2/dmas

    do iel = 1, ncel
      cvar_yco2(iel) = oxyco2(ioxy)*wmco2/dmas
    enddo
  endif

  if (ieqnox.eq.1) then
    do iel = 1, ncel
      cvar_hox(iel) = h1init
    enddo
  endif

endif

do iel = 1, ncel
  ! Initialization of the continuous mass fraction
  cpro_x1(iel) = 1.d0
enddo

! Initialization of the continuous mass fraction AT the BCs
do ifac = 1, nfabor
  bpro_x1(ifac) = 1.d0
enddo

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
