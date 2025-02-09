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

subroutine ebuphy

!===============================================================================
! FONCTION :
! --------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME DE PREMELAMGE MODELE EBU
! Calcul de RHO adiabatique ou permeatique (transport de H)


! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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

! Local variables

integer          igg, iel
integer          ifac
double precision coefg(ngazgm), ygfm, ygbm, epsi
double precision nbmol , temsmm , fmel , ckabgf, ckabgb
double precision masmgb, hgb, tgb, masmgf

! Enthalpie massique gaz frais identique pour premelange frais et dilution
double precision hgf

double precision, allocatable, dimension(:) :: yfuegf, yoxygf, yprogf
double precision, allocatable, dimension(:) :: yfuegb, yoxygb, yprogb
double precision, allocatable, dimension(:) :: temp, masmel
double precision, dimension(:), pointer :: brom, crom
double precision, dimension(:), pointer :: bsval
double precision, dimension(:), pointer :: cvar_ygfm, cvar_fm
double precision, dimension(:), pointer :: cvar_scalt
double precision, dimension(:), pointer :: cpro_temp,cpro_ymgg
double precision, dimension(:), pointer :: cpro_ym1, cpro_ym2, cpro_ym3
double precision, dimension(:), pointer :: cpro_ckabs, cpro_t4m, cpro_t3m

integer       ipass
data          ipass /0/
save          ipass

!===============================================================================

interface

  subroutine cs_combustion_boundary_conditions_density_ebu_lw()  &
    bind(C, name='cs_combustion_boundary_conditions_density_ebu_lw')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_combustion_boundary_conditions_density_ebu_lw

end interface

!===============================================================================
! 0. ON COMPTE LES PASSAGES
!===============================================================================

ipass = ipass + 1

!===============================================================================
! 1. INITIALISATIONS A CONSERVER
!===============================================================================

! Allocate temporary arrays
allocate(yfuegf(ncelet), yoxygf(ncelet), yprogf(ncelet))
allocate(yfuegb(ncelet), yoxygb(ncelet), yprogb(ncelet))
allocate(temp(ncelet), masmel(ncelet))

! Initialize variables to avoid compiler warnings

ckabgb = 0.d0
ckabgf = 0.d0

! --- Initialisation memoire


! ---> Initialisation

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

! ---> Positions des variables, coefficients

call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)
call field_get_val_s(itemp, cpro_temp)
call field_get_val_s(iym(1), cpro_ym1)
call field_get_val_s(iym(2), cpro_ym2)
call field_get_val_s(iym(3), cpro_ym3)

if ( iirayo.gt.0 ) then
  call field_get_val_s(ickabs, cpro_ckabs)
  call field_get_val_s(it4m, cpro_t4m)
  call field_get_val_s(it3m, cpro_t3m)
endif

call field_get_val_s(iygfm, cvar_ygfm)
if ( ippmod(icoebu).ne.0 .and. ippmod(icoebu).ne.1 ) then
  call field_get_val_s(ifm, cvar_fm)
endif
if ( ippmod(icoebu).eq.1 .or. ippmod(icoebu).eq.3 ) then
  call field_get_val_s(ihm, cvar_scalt)
endif

!===============================================================================
! 2. DETERMINATION DES GRANDEURS THERMOCHIMIQUES
!===============================================================================

! ---> Grandeurs GAZ FRAIS

! ----- Fournies par l'utilisateur
!       FMEL         --> Taux de melange
!                        constant pour les options 0 et 1
!                        variable sinon
!       TGF          --> Temperature gaz frais en K identique
!                        pour premelange frais et dilution
! ----- Deduites
!       YFUEGF(    .)    --> Fraction massique fuel gaz frais
!       YOXYGF(    .)    --> Fraction massique oxydant gaz frais
!       YPROGF(    .)    --> Fraction massique produits gaz frais
!       HGF          --> Enthalpie massique gaz frais identique
!                        pour premelange frais et dilution
!       MASMGF       --> Masse molaire gaz frais
!       CKABGF       --> Coefficient d'absorption

! ---> Grandeurs GAZ BRULES

! ----- Deduites
!       TGB          --> Temperature gaz brules en K
!       YFUEGB(    .)    --> Fraction massique fuel gaz brules
!       YOXYGB(    .)    --> Fraction massique oxydant gaz brules
!       YPROGB(    .)    --> Fraction massique produits gaz brules
!       MASMGB       --> Masse molaire gaz brules
!       CKABGB       --> Coefficient d'absorption

! ---> Grandeurs MELANGE

!       masmel           --> Masse molaire du melange
!       cpro_temp        --> Temperature du melange
!       crom             --> Masse volumique du melange
!       cpro_.f,o,p      --> Fractions massiques en F, O, P
!       cpro_ckabs       --> Coefficient d'absorption
!       cpro_t4m         --> terme T^4
!       cpro_t3m         --> terme T^3


! ---> Fractions massiques des gaz frais et brules en F, O, P

do iel = 1, ncel

  if ( ippmod(icoebu).eq.0 .or. ippmod(icoebu).eq.1 ) then
    fmel = frmel
  else
    fmel = cvar_fm(iel)
  endif

  yfuegf(iel) = fmel
  yoxygf(iel) = 1.d0-fmel
  yprogf(iel) = 0.d0

  yfuegb(iel) = max(zero,(fmel-fs(1))/(1.d0-fs(1)))
  yprogb(iel) = (fmel-yfuegb(iel))/fs(1)
  yoxygb(iel) = 1.d0 - yfuegb(iel) - yprogb(iel)

enddo

epsi = 1.d-06

do iel = 1, ncel

! ---> Coefficients d'absorption des gaz frais et brules

  if ( iirayo.gt.0 ) then
     ckabgf = yfuegf(iel)*ckabsg(1) + yoxygf(iel)*ckabsg(2)       &
            + yprogf(iel)*ckabsg(3)
     ckabgb = yfuegb(iel)*ckabsg(1) + yoxygb(iel)*ckabsg(2)       &
            + yprogb(iel)*ckabsg(3)
  endif

! ---> Masse molaire des gaz frais

  coefg(1) = yfuegf(iel)
  coefg(2) = yoxygf(iel)
  coefg(3) = yprogf(iel)
  nbmol = 0.d0
  do igg = 1, ngazg
    nbmol = nbmol + coefg(igg)/wmolg(igg)
  enddo
  masmgf = 1.d0/nbmol

! ---> Calcul de l'enthalpie des gaz frais

  hgf = cs_gas_combustion_t_to_h(coefg, tgf)

! ---> Masse molaire des gaz brules

  coefg(1) = yfuegb(iel)
  coefg(2) = yoxygb(iel)
  coefg(3) = yprogb(iel)
  nbmol = 0.d0
  do igg = 1, ngazg
    nbmol = nbmol + coefg(igg)/wmolg(igg)
  enddo
  masmgb = 1.d0/nbmol

  ygfm = cvar_ygfm(iel)
  ygbm = 1.d0 - ygfm

! ---> Masse molaire du melange

  masmel(iel) = 1.d0 / ( ygfm/masmgf + ygbm/masmgb )

! ---> Calcul Temperature des gaz brules

  if ( ippmod(icoebu).eq.0 .or. ippmod(icoebu).eq.2 ) then
! ---- EBU Standard et modifie en conditions adiabatiques (sans H)
    hgb = hgf
  else if ( ippmod(icoebu).eq.1 .or. ippmod(icoebu).eq.3 ) then
! ---- EBU Standard et modifie en conditions permeatiques (avec H)
    hgb = hgf
    if ( ygbm.gt.epsi ) then
      hgb = ( cvar_scalt(iel)-hgf*ygfm ) / ygbm
    endif
  endif

  tgb = cs_gas_combustion_h_to_t(coefg, hgb)

! ---> Temperature du melange
!      Rq PPl : Il serait plus judicieux de ponderer par les CP (GF et GB)
  cpro_temp(iel) = ygfm*tgf + ygbm*tgb

! ---> Temperature / Masse molaire

  temsmm = ygfm*tgf/masmgf + ygbm*tgb/masmgb

! ---> Masse volumique du melange

  if (ipass.gt.1.or.(isuite.eq.1.and.initro.eq.1)) then
    crom(iel) = srrom*crom(iel)                          &
              + (1.d0-srrom)*                            &
                ( p0/(cs_physical_constants_r*temsmm) )
  endif

! ---> Fractions massiques des especes globales

  cpro_ym1(iel) = yfuegf(iel)*ygfm + yfuegb(iel)*ygbm
  cpro_ym2(iel) = yoxygf(iel)*ygfm + yoxygb(iel)*ygbm
  cpro_ym3(iel) = yprogf(iel)*ygfm + yprogb(iel)*ygbm

! ---> Grandeurs relatives au rayonnement

  if ( iirayo.gt.0 ) then
    cpro_ckabs(iel) = ygfm*ckabgf + ygbm*ckabgb
    cpro_t4m(iel)  = ygfm*tgf**4 + ygbm*tgb**4
    cpro_t3m(iel)  = ygfm*tgf**3 + ygbm*tgb**3
  endif

enddo

!===============================================================================
! 3. CALCUL DE RHO ET DES FRACTIONS MASSIQUES DES ESPECES GLOBALES
!    SUR LES BORDS
!===============================================================================

! --> Masse volumique au bord

call cs_combustion_boundary_conditions_density_ebu_lw()

! --> Fractions massiques des especes globales au bord
!     Uniquement si transport de H

do igg = 1, ngazg
  call field_get_val_s(ibym(1), bsval)
  call field_get_val_s(iym(igg),cpro_ymgg)
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    bsval(ifac) = cpro_ymgg(iel)
  enddo
enddo

! Free memory
deallocate(yfuegf, yoxygf, yprogf)
deallocate(yfuegb, yoxygb, yprogb)
deallocate(temp, masmel)

!===============================================================================
! FORMATS
!----


!----
! FIN
!----

return
end subroutine
