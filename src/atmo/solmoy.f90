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
!> \file solmoy.f90
!> \brief Atmospheric soil module - Initialize ground level parameters from land use
!
!> \brief
!>   Compute soil-atmosphere coefficients
!>   we knowx :
!>    - the percentage of each category per face (boundary zone field)
!>    - coefficients values for each category
!>
!>   So we compute the average value for all faces of the zone
!>    e.g. for 3 categories of soil (water, forest, building)
!>     albedo(face) = albedo(water)    * % surface_water   (face)
!>                  + albedo(forest)   * % surface_forest  (face)
!>                  + albedo(building) * % surface_building(face)
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[out]   ierreu         code error
!-------------------------------------------------------------------------------

subroutine solmoy ( ierreu )

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
use ppincl
use atincl
use atsoil
use field
use mesh

implicit none

!===============================================================================

integer          ierreu

! Local variables

integer iirugdy,iirugth,iialbed,iiemiss,iicsol,iiveget
integer iic1w,iic2w,iir1,iir2
integer isol, n, ifac
double precision codinv
double precision rugdij,rugtij,albeij,emisij
double precision vegeij,c1wij,c2wij,csolij
double precision r1ij,r2ij
double precision surf_zone
double precision solmax(10),solmea(10),solmin(10)
character(len=12) :: solnom(10)

integer, dimension(:), pointer :: elt_ids

double precision, pointer, dimension(:)   :: bpro_rugdyn ! all boundary faces
double precision, pointer, dimension(:)   :: bpro_rugthe ! all boundary faces
double precision, pointer, dimension(:)   :: bpro_albedo ! all boundary faces
double precision, pointer, dimension(:)   :: bpro_emissi ! all boundary faces
double precision, pointer, dimension(:)   :: bpro_vegeta ! only zone faces
double precision, pointer, dimension(:)   :: bpro_c1w    ! only zone faces
double precision, pointer, dimension(:)   :: bpro_c2w    ! only zone faces
double precision, pointer, dimension(:)   :: bpro_csol   ! only zone faces
double precision, pointer, dimension(:)   :: bpro_r1     ! only zone faces
double precision, pointer, dimension(:)   :: bpro_r2     ! only zone faces
double precision, pointer, dimension(:)   :: bpro_tprof  ! only zone faces
double precision, pointer, dimension(:,:) :: bpro_pourcent_sol ! only zone faces

!  ================================================================
!     1) initialisations
!  ================================================================

call field_get_val_s_by_name("boundary_roughness", bpro_rugdyn)
call field_get_val_s_by_name("boundary_thermal_roughness", bpro_rugthe)
call field_get_val_s_by_name("boundary_albedo", bpro_albedo)
call field_get_val_s_by_name("boundary_emissivity", bpro_emissi)
call field_get_val_s_by_name("boundary_vegetation", bpro_vegeta)
call field_get_val_s_by_name("soil_water_capacity", bpro_c1w   )
call field_get_val_s_by_name("soil_water_ratio", bpro_c2w   )
call field_get_val_s_by_name("soil_thermal_capacity", bpro_csol  )
call field_get_val_s_by_name("soil_r1", bpro_r1    )
call field_get_val_s_by_name("soil_r2", bpro_r2    )
call field_get_val_s_by_name("soil_temperature_deep", bpro_tprof )

call field_get_val_v_by_name("atmo_soil_percentages", bpro_pourcent_sol )

call atmo_get_soil_zone(nfmodsol, nbrsol, elt_ids)

codinv = -999.d0
do isol = 1, nfmodsol

  ifac = elt_ids(isol) + 1 ! C > Fortran
  bpro_rugdyn(ifac) = codinv
  bpro_rugthe(ifac) = codinv
  bpro_albedo(ifac) = codinv
  bpro_emissi(ifac) = codinv

  bpro_vegeta(isol) = codinv
  bpro_c1w   (isol) = codinv
  bpro_c2w   (isol) = codinv
  bpro_csol  (isol) = codinv
  bpro_r1    (isol) = codinv
  bpro_r2    (isol) = codinv
enddo

iirugdy = 1
iirugth = 2
iialbed = 3
iiemiss = 4
iicsol  = 5
iiveget = 6
iic1w   = 7
iic2w   = 8
iir1    = 9
iir2    = 10

solnom(iirugdy)='z0 dynamique'
solnom(iirugth)='z0 thermique'
solnom(iialbed)='albedo      '
solnom(iiemiss)='emissivite  '
solnom(iicsol )='csol (x1e+6)'
solnom(iiveget)='vegetation  '
solnom(iic1w  )='c1w         '
solnom(iic2w  )='c2w         '
solnom(iir1   )='r1          '
solnom(iir2   )='r2          '

!  ================================================================
!     2) calcul des coefficients pour chaque maille
!  ================================================================

do isol = 1, nfmodsol
  rugdij = 0.d0
  rugtij = 0.d0
  albeij = 0.d0
  emisij = 0.d0
  csolij = 0.d0
  vegeij = 0.d0
  c1wij  = 0.d0
  c2wij  = 0.d0
  r1ij   = 0.d0
  r2ij   = 0.d0

  ! zrrel = z(2)*(z(km)-CDGFBO(3,IFAC))/z(km)

  do n = 1, nbrsol
    ! shifting of 1 for the percentage because starting by default value
    rugdij = rugdij + tab_sol(n)%rugdyn*bpro_pourcent_sol(n+1,isol)/100.d0
    rugtij = rugtij + tab_sol(n)%rugthe*bpro_pourcent_sol(n+1,isol)/100.d0
    albeij = albeij + tab_sol(n)%albedo*bpro_pourcent_sol(n+1,isol)/100.d0
    emisij = emisij + tab_sol(n)%emissi*bpro_pourcent_sol(n+1,isol)/100.d0
    csolij = csolij + tab_sol(n)%csol  *bpro_pourcent_sol(n+1,isol)/100.d0
    vegeij = vegeij + tab_sol(n)%vegeta*bpro_pourcent_sol(n+1,isol)/100.d0
    c1wij  = c1wij  + tab_sol(n)%c1w   *bpro_pourcent_sol(n+1,isol)/100.d0
    c2wij  = c2wij  + tab_sol(n)%c2w   *bpro_pourcent_sol(n+1,isol)/100.d0
    r1ij   = r1ij   + tab_sol(n)%r1    *bpro_pourcent_sol(n+1,isol)/100.d0
    r2ij   = r2ij   + tab_sol(n)%r2    *bpro_pourcent_sol(n+1,isol)/100.d0
  enddo

  ifac = elt_ids(isol) + 1 ! C > Fortran
  bpro_rugdyn(ifac) = rugdij
  bpro_rugthe(ifac) = rugtij
  bpro_albedo(ifac) = albeij
  bpro_emissi(ifac) = emisij

  bpro_csol  (isol) = csolij
  bpro_vegeta(isol) = vegeij
  bpro_c1w   (isol) = c1wij
  bpro_c2w   (isol) = c2wij
  bpro_r1    (isol) = r1ij
  bpro_r2    (isol) = r2ij

  ! Pour temperatures profondes dans un premier temps on initialise a tprini
  bpro_tprof(isol) = tprini
enddo

!  ================================================================
! 3) controle
!  ================================================================

ierreu = 0
do isol = 1, nfmodsol
  ifac = elt_ids(isol) + 1 ! C > Fortran
  if(bpro_rugdyn(ifac) .eq. codinv) ierreu = ierreu + 1
  if(bpro_rugthe(ifac) .eq. codinv) ierreu = ierreu + 1
  if(bpro_albedo(ifac) .eq. codinv) ierreu = ierreu + 1
  if(bpro_emissi(ifac) .eq. codinv) ierreu = ierreu + 1

  if(bpro_csol  (isol) .eq. codinv) ierreu = ierreu + 1
  if(bpro_vegeta(isol) .eq. codinv) ierreu = ierreu + 1
  if(bpro_c1w   (isol) .eq. codinv) ierreu = ierreu + 1
  if(bpro_c2w   (isol) .eq. codinv) ierreu = ierreu + 1
  if(bpro_r1    (isol) .eq. codinv) ierreu = ierreu + 1
  if(bpro_r2    (isol) .eq. codinv) ierreu = ierreu + 1
enddo

if (irangp.ge.0) then
  call parcpt(ierreu)
endif

! impression eventuelle d'un message d'erreur

if (ierreu.ne.0) then
  write(nfecra,9999)
  write(nfecra,9991) ierreu

  ! ou impression de controle

else
  do n = 1, 10
    solmin(n) = +999999.d0
    solmea(n) = 0.d0
    solmax(n) = -999999.d0
  enddo
  surf_zone = 0.d0

  do isol = 1, nfmodsol

    ifac = elt_ids(isol) + 1 ! C > Fortran
    if (bpro_rugdyn(ifac) .gt. solmax(1)) solmax(1) &
         = bpro_rugdyn(ifac)
    if (bpro_rugthe(ifac) .gt. solmax(2)) solmax(2) &
         = bpro_rugthe(ifac)
    if (bpro_albedo(ifac) .gt. solmax(3)) solmax(3) &
         = bpro_albedo(ifac)
    if (bpro_emissi(ifac) .gt. solmax(4)) solmax(4) &
         = bpro_emissi(ifac)

    if (bpro_csol(isol)   .gt. solmax(5)) solmax(5) &
         = bpro_csol(isol)
    if (bpro_vegeta(isol) .gt. solmax(6)) solmax(6) &
         = bpro_vegeta(isol)
    if (bpro_c1w(isol)    .gt. solmax(7)) solmax(7) &
         = bpro_c1w(isol)
    if (bpro_c2w(isol)    .gt. solmax(8)) solmax(8) &
         = bpro_c2w(isol)
    if (bpro_r1(isol)     .gt. solmax(9)) solmax(9) &
         = bpro_r1(isol)
    if (bpro_r2(isol)     .gt. solmax(10))solmax(10)&
         = bpro_r2(isol)

    ifac = elt_ids(isol) + 1 ! C > Fortran
    if (bpro_rugdyn(ifac) .lt. solmin(1)) solmin(1) &
         = bpro_rugdyn(ifac)
    if (bpro_rugthe(ifac) .lt. solmin(2)) solmin(2) &
         = bpro_rugthe(ifac)
    if (bpro_albedo(ifac) .lt. solmin(3)) solmin(3) &
         = bpro_albedo(ifac)
    if (bpro_emissi(ifac) .lt. solmin(4)) solmin(4) &
         = bpro_emissi(ifac)

    if (bpro_csol(isol)   .lt. solmin(5)) solmin(5) &
         = bpro_csol(isol)
    if (bpro_vegeta(isol) .lt. solmin(6)) solmin(6) &
         = bpro_vegeta(isol)
    if (bpro_c1w(isol)    .lt. solmin(7)) solmin(7) &
         = bpro_c1w(isol)
    if (bpro_c2w(isol)    .lt. solmin(8)) solmin(8) &
         = bpro_c2w(isol)
    if (bpro_r1(isol)     .lt. solmin(9)) solmin(9) &
         = bpro_r1(isol)
    if (bpro_r2(isol)     .lt. solmin(10))solmin(10)&
         = bpro_r2(isol)

    ifac = elt_ids(isol) + 1 ! C > Fortran

    solmea(1) = solmea(1) + surfbn(ifac) * bpro_rugdyn(ifac)
    solmea(2) = solmea(2) + surfbn(ifac) * bpro_rugthe(ifac)
    solmea(3) = solmea(3) + surfbn(ifac) * bpro_albedo(ifac)
    solmea(4) = solmea(4) + surfbn(ifac) * bpro_emissi(ifac)

    solmea(5) = solmea(5) + surfbn(ifac) * bpro_csol  (isol)
    solmea(6) = solmea(6) + surfbn(ifac) * bpro_vegeta(isol)
    solmea(7) = solmea(7) + surfbn(ifac) * bpro_c1w   (isol)
    solmea(8) = solmea(8) + surfbn(ifac) * bpro_c2w   (isol)
    solmea(9) = solmea(9) + surfbn(ifac) * bpro_r1    (isol)
    solmea(10)= solmea(10)+ surfbn(ifac) * bpro_r2    (isol)

    ! Surface of the zone, could use it directly in C
    surf_zone = surf_zone + surfbn(ifac)
  enddo

  if (irangp.ge.0) then
    do n = 1, 10
      call parsom(solmea(n))
      call parmin(solmin(n))
      call parmax(solmax(n))
    enddo
    call parsom(surf_zone)
  endif

  do n = 1, 10
    solmea(n)= solmea(n) / surf_zone
  enddo

  write(nfecra,3001)
  write(nfecra,3002)
  n = iirugdy
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iirugth
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iialbed
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iiemiss
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iicsol
  write(nfecra,3003)solnom(n),solmin(n)*1.d6,solmea(n)*1.d6                     &
       ,solmax(n)*1.d6
  n = iiveget
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iic1w
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iic2w
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iir1
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  n = iir2
  write(nfecra,3003)solnom(n),solmin(n),solmea(n),solmax(n)
  write(nfecra,3004)
endif

! sortie

return

!--------
! formats
!--------

9999 format(//,5x,'%% erreur solmoy: erreur numero  1')
9991 format( 23x,'initialisation incorrecte des coefficients de ',&
     'l''interface sol atmosphere',                         &
             /,23x,'sur',i9,' valeurs il y en a',i9,' non initialisees')
3001 format(//,8x,' ** ========================================= **',&
             /,8x,' ** interface sol atmosphere                  **',      &
             /,8x,' ** valeurs des constantes calculees          **',      &
             /,8x,' ** ========================================= **',/)
3002 format(   ' *            * minimum* moyenne* maximum*')
3003 format(   ' *',a12,'*',3(f8.4,'*'))
3004 format(/,8x,' ** ========================================= **',//)

end subroutine solmoy
