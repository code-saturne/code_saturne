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
! Function :
! --------

!> \file pptycl.f90
!>
!> \brief Boundary conditions for specific physics modules.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     init          partial treatment (before time loop) if true
!> \param[in,out] itypfb        boundary face types
!> \param[in,out] izfppp        index of the zone for the boundary faces
!>                               (for the specific physics)
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine pptycl                 &
 ( init , itypfb , izfppp , dt )  &
 bind(C, name='cs_f_pptycl')

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

! Arguments

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use mesh
use field
use cs_c_bindings
use cs_c_bindings
use dimens, only: nvar

!===============================================================================

implicit none

! Arguments

logical(c_bool), value :: init

integer(c_int) ::  itypfb(nfabor)
integer(c_int) ::  izfppp(nfabor)

real(c_double) :: dt(ncelet)

! Local variables

integer          ifac, iok, ifvu, ii, izone, izonem
integer, pointer, dimension(:,:) :: icodcl
double precision, pointer, dimension(:,:,:) :: rcodcl

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_coal_boundary_conditions(bc_type)  &
    bind(C, name='cs_coal_boundary_conditions')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), dimension(*) :: bc_type
  end subroutine cs_coal_boundary_conditions

  subroutine cs_ctwr_bcond()  &
    bind(C, name='cs_ctwr_bcond')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_ctwr_bcond

  subroutine cs_atmo_bcond()  &
    bind(C, name='cs_atmo_bcond')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_atmo_bcond

end interface

!===============================================================================
! 1. Zones list (for some models)
!===============================================================================

if (ippmod(icompf).lt.0.and.ippmod(iatmos).eq.-1.and.ippmod(iaeros).eq.-1) then
  ! --> faces all belong to a boundary zone
  iok = 0

  do ifac = 1, nfabor
    if(izfppp(ifac).le.0.or.izfppp(ifac).gt.nozppm) then
      iok = iok + 1
      izfppp(ifac) = min(0, izfppp(ifac))
    endif
  enddo

  if (irangp.ge.0) call parcmx(iok)

  if (iok.gt.0) then
    write(nfecra,1000) nozppm
    call boundary_conditions_error(izfppp)
  endif

  ! a list gathering numbers of boundary zones is built
  ! (list is local to a sub-domain in parallel)
  nzfppp = 0
  do ifac = 1, nfabor
    ifvu = 0
    do ii = 1, nzfppp
      if (ilzppp(ii).eq.izfppp(ifac)) then
        ifvu = 1
      endif
    enddo
    if (ifvu.eq.0) then
      nzfppp = nzfppp + 1
      if (nzfppp.le.nozppm) then
        ilzppp(nzfppp) = izfppp(ifac)
      else
        write(nfecra,1001) nozppm
        write(nfecra,1002)(ilzppp(ii),ii=1,nozppm)
        call csexit (1)
      endif
    endif
  enddo

  ! maximum zone number

  izonem = 0
  do ii = 1, nzfppp
    izone = ilzppp(ii)
    izonem = max(izonem,izone)
  enddo
  if (irangp.ge.0) then
    call parcmx(izonem)
  endif
  nozapm = izonem
endif

!===============================================================================
! 2. Call to boundary conditions computations, model by model.
!    Computations should not be called under initialization for most
!    models; those for which it should be called are tested first.
!===============================================================================

call field_build_bc_codes_all(icodcl, rcodcl) ! Get map

! Atmospheric flows
if (ippmod(iatmos).ge.0) then
  call attycl(itypfb, icodcl, rcodcl)
  call cs_atmo_bcond()
endif

! Cooling towers
if (ippmod(iaeros).ge.0) then
  call cs_ctwr_bcond()
endif

if (init .eqv. .true.) return

! ---> Chimie 3 points : USD3PC

if (ippmod(icod3p).ge.0) then
  call d3ptcl(itypfb, izfppp, icodcl, rcodcl)

! ---> Steady laminar flamelet

elseif (ippmod(islfm).ge.0) then
  call cs_steady_laminar_flamelet_bcond(itypfb, izfppp, icodcl, rcodcl)

! ---> Combustion gaz USEBUC
!      Flamme de premelange modele EBU

elseif (ippmod(icoebu).ge.0) then
  call ebutcl(itypfb, izfppp, rcodcl)

! ---> Combustion gaz USLWCC
!      Flamme de premelange modele LWC

elseif (ippmod(icolwc).ge.0) then
  call lwctcl(itypfb, izfppp, rcodcl)

! ---> Combustion charbon pulverise

elseif (ippmod(iccoal).ge.0) then
  call cs_coal_boundary_conditions(itypfb)

endif

!--------
! Formats
!--------

 1000 format(                                                           &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : PHYSIQUE PARTICULIERE'                       ,/,&
'@    ========='                                               ,/,&
'@    LES CONDITIONS AUX LIMITES SONT INCOMPLETES OU ERRONEES' ,/,&
'@'                                                            ,/,&
'@  Le numero de zone associee a certaines faces doit etre'    ,/,&
'@    un entier strictement positif et inferieur ou egal a'    ,/,&
'@    NOZPPM = ',I10                                           ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@'                                                            ,/,&
'@  Verifier les conditions aux limites.'                      ,/,&
'@'                                                            ,/,&
'@  Vous pouvez visualiser les faces de bord sorties en'       ,/,&
'@  erreur.'                                                   ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PHYSIQUE PARTICULIERE                       ',/,&
'@    =========                                               ',/,&
'@    PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/,&
'@                                                            ',/,&
'@  Le nombre maximal de zones frontieres qui peuvent etre    ',/,&
'@    definies par l''utilisateur est NBZPPM = ',I10           ,/,&
'@    Il a ete depasse.                                       ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites.                      ',/,&
'@                                                            ',/,&
'@  Les NBZPPM premieres zones frontieres                     ',/,&
'@    portent ici les numeros suivants :                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1002 format(i10)

!----
! End
!----

return
end subroutine
