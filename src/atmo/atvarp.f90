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

!===============================================================================
! Function:
! ---------

!> \file atvarp.f90
!> \brief Declare additional transported variables for atmospheric module.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine atvarp

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use ppincl
use ihmpre
use atincl
use atchem
use field
use siream

!===============================================================================

implicit none

! Local variables

integer        ii, jj, isc, f_id
integer        kscmin, kscmax
character(len=80) :: name, label
character(len=2) :: nbin
character(len=10), dimension(nesp_aer) :: f_esp_siream

integer(c_int) :: n_chem_species
integer(c_int), dimension(65) :: species_f_id

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_field_pointer_map_atmospheric(n_chem_species, species_f_id)  &
    bind(C, name='cs_field_pointer_map_atmospheric')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value        :: n_chem_species
    integer(c_int), dimension(*) :: species_f_id
  end subroutine cs_field_pointer_map_atmospheric

end interface

!===============================================================================
! 0. Definitions for fields
!===============================================================================

! Key ids for clipping
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)

!===============================================================================
! 1. GUI and model information
!===============================================================================

if (iihmpr.eq.1) then
  call uiati1 (imeteo, ficmet, len(ficmet))
endif

! Set base default model parameters and call usati1

call atini0

!===============================================================================
! 1. Add variables
!===============================================================================

! 1.1  Dry atmosphere
! =====================

if (ippmod(iatmos).eq.1) then

  ! Potential temperature, in Kelvin
  itherm = 1
  itpscl = 1
  call add_model_scalar_field('temperature', 'PotTemp', iscalt)
  f_id = ivarfl(isca(iscalt))
  call field_set_key_double(f_id, kscmin, 0.d0)

endif

! 1.2  Humid atmosphere
! =====================

if (ippmod(iatmos).eq.2) then

  ! Potential temperature, in Kelvin
  itherm = 1
  itpscl = 1
  call add_model_scalar_field('temperature', 'LqPotTmp', iscalt)
  f_id = ivarfl(isca(iscalt))
  call field_set_key_double(f_id, kscmin, 200.d0)

  ! total water content
  call add_model_scalar_field('total_water', 'TotWater', itotwt)
  f_id = ivarfl(isca(itotwt))
  call field_set_key_double(f_id, kscmin, 0.d0)

  ! total number of droplets
  call add_model_scalar_field('number_of_droplets', 'TotDrop', intdrp)
  f_id = ivarfl(isca(intdrp))
  call field_set_key_double(f_id, kscmin, 0.d0)

endif

!===============================================================================
! 3. Chemistry variables
!===============================================================================

! Atmospheric gaseous chemistry
! Do not change this order
if (iaerosol.eq.1) ichemistry = 3
! if a chemical scheme is solved, a concentration profiles
! file must be used
if (ichemistry.ge.1) ifilechemistry = ichemistry
if (inogaseouschemistry.eq.1) ichemistry = 0

if (ifilechemistry.ge.1) then

  ! logical unit and name of the chemical profiles file
  impmec = 27
  ficmec = 'chemistry'

  ! Initialization of the chemical scheme
  ! quasi steady equilibrium NOx scheme with 4 species and 5 reactions
  if (ifilechemistry.eq.1) then
    nrg = 5
    nespg = 4
    allocate(isca_chem(nespg))
    allocate(dmmk(nespg))
    allocate(chempoint(nespg))
    call add_model_scalar_field('chemistry_no',  'NO',  isca_chem(1))
    call add_model_scalar_field('chemistry_no2', 'NO2', isca_chem(2))
    call add_model_scalar_field('chemistry_o3',  'O3',  isca_chem(3))
    call add_model_scalar_field('chemistry_o3p', 'O3P', isca_chem(4))
    dmmk(1)=30.d-3 ! Molar mass NO
    dmmk(2)=46.d-3 ! Molar mass NO2
    dmmk(3)=48.d-3 ! Molar mass O3
    dmmk(4)=16.d-3 ! Molar mass O3P
    chempoint = (/ 4, 3, 2, 1 /)
  ! scheme with 20 species and 34 reactions
  else if (ifilechemistry.eq.2) then
    nrg = 34
    nespg = 20
    allocate(isca_chem(nespg))
    allocate(dmmk(nespg))
    allocate(chempoint(nespg))
    call add_model_scalar_field('species_no',    'NO',    isca_chem(1))
    call add_model_scalar_field('species_no2',   'NO2',   isca_chem(2))
    call add_model_scalar_field('species_o3',    'O3',    isca_chem(3))
    call add_model_scalar_field('species_o3p',   'O3P',   isca_chem(4))
    call add_model_scalar_field('species_o1d',   'O1D',   isca_chem(5))
    call add_model_scalar_field('species_oh',    'OH',    isca_chem(6))
    call add_model_scalar_field('species_ho2',   'HO2',   isca_chem(7))
    call add_model_scalar_field('species_h2o2',  'H2O2',  isca_chem(8))
    call add_model_scalar_field('species_no3',   'NO3',   isca_chem(9))
    call add_model_scalar_field('species_n2o5',  'N2O5',  isca_chem(10))
    call add_model_scalar_field('species_hono',  'HONO',  isca_chem(11))
    call add_model_scalar_field('species_hno3',  'HNO3',  isca_chem(12))
    call add_model_scalar_field('species_co',    'CO',    isca_chem(13))
    call add_model_scalar_field('species_hcho',  'HCHO',  isca_chem(14))
    call add_model_scalar_field('species_ald2',  'ALD2',  isca_chem(15))
    call add_model_scalar_field('species_c2o3',  'C2O3',  isca_chem(16))
    call add_model_scalar_field('species_pan',   'PAN',   isca_chem(17))
    call add_model_scalar_field('species_xo2',   'XO2',   isca_chem(18))
    call add_model_scalar_field('species_so2',   'SO2',   isca_chem(19))
    call add_model_scalar_field('species_h2so4', 'H2SO4', isca_chem(20))
    dmmk(1)=30.d-3     ! Molar mass NO
    dmmk(2)=46.d-3     ! Molar mass NO2
    dmmk(3)=48.d-3     ! Molar mass O3
    dmmk(4)=16.d-3     ! Molar mass O3P
    dmmk(5)=16.d-3     ! Molar mass O1D
    dmmk(6)=17.01d-3   ! Molar mass OH
    dmmk(7)=33.01d-3   ! Molar mass HO2
    dmmk(8)=34.01d-3   ! Molar mass H2O2
    dmmk(9)=62.01d-3   ! Molar mass NO3
    dmmk(10)=108.01d-3 ! Molar mass N2O5
    dmmk(11)=47.01d-3  ! Molar mass HONO
    dmmk(12)=63.01d-3  ! Molar mass HNO3
    dmmk(13)=28.01d-3  ! Molar mass CO
    dmmk(14)=30.03d-3  ! Molar mass HCHO
    dmmk(15)=44.05d-3  ! Molar mass ALD2
    dmmk(16)=75.04d-3  ! Molar mass C2O3
    dmmk(17)=121.05d-3 ! Molar mass PAN
    dmmk(18)=47.03d-3  ! Molar mass XO2
    dmmk(19)=64.06d-3  ! Molar mass SO2
    dmmk(20)=98.08d-3  ! Molar mass H2SO4
    chempoint = (/ 20, 19, 16, 17, 2, 15, 14, 3, 18, 7, 8, 9, 4, &
                   10, 1, 12, 11, 13, 5, 6 /)
  ! scheme CB05 with 52 species and 155 reactions
  else if (ifilechemistry.eq.3) then
    if (iaerosol.eq.1) then
      nrg = 162
      nespg = 65
    else
      nrg = 155
      nespg = 52
    endif
    allocate(isca_chem(nespg))
    allocate(dmmk(nespg))
    allocate(chempoint(nespg))

    call add_model_scalar_field('species_no',    'NO',    isca_chem(1))
    call add_model_scalar_field('species_no2',   'NO2',   isca_chem(2))
    call add_model_scalar_field('species_o3',    'O3',    isca_chem(3))
    call add_model_scalar_field('species_o3p',   'O3P',   isca_chem(4))
    call add_model_scalar_field('species_o1d',   'O1D',   isca_chem(5))
    call add_model_scalar_field('species_oh',    'OH',    isca_chem(6))
    call add_model_scalar_field('species_ho2',   'HO2',   isca_chem(7))
    call add_model_scalar_field('species_h2o2',  'H2O2',  isca_chem(8))
    call add_model_scalar_field('species_no3',   'NO3',   isca_chem(9))
    call add_model_scalar_field('species_n2o5',  'N2O5',  isca_chem(10))
    call add_model_scalar_field('species_hono',  'HONO',  isca_chem(11))
    call add_model_scalar_field('species_hno3',  'HNO3',  isca_chem(12))
    call add_model_scalar_field('species_hno4',  'HNO4',  isca_chem(13))
    call add_model_scalar_field('species_co',    'CO',    isca_chem(14))
    call add_model_scalar_field('species_hcho',  'HCHO',  isca_chem(15))
    call add_model_scalar_field('species_ald2',  'ALD2',  isca_chem(16))
    call add_model_scalar_field('species_c2o3',  'C2O3',  isca_chem(17))
    call add_model_scalar_field('species_pan',   'PAN',   isca_chem(18))
    call add_model_scalar_field('species_aldx',  'ALDX',  isca_chem(19))
    call add_model_scalar_field('species_cxo3',  'CXO3',  isca_chem(20))
    call add_model_scalar_field('species_panx',  'PANX',  isca_chem(21))
    call add_model_scalar_field('species_xo2',   'XO2',   isca_chem(22))
    call add_model_scalar_field('species_xo2n',  'XO2N',  isca_chem(23))
    call add_model_scalar_field('species_ntr',   'NTR',   isca_chem(24))
    call add_model_scalar_field('species_etoh',  'ETOH',  isca_chem(25))
    call add_model_scalar_field('species_ch4',   'CH4',   isca_chem(26))
    call add_model_scalar_field('species_meo2',  'MEO2',  isca_chem(27))
    call add_model_scalar_field('species_meoh',  'MEOH',  isca_chem(28))
    call add_model_scalar_field('species_mepx',  'MEPX',  isca_chem(29))
    call add_model_scalar_field('species_facd',  'FACD',  isca_chem(30))
    call add_model_scalar_field('species_etha',  'ETHA',  isca_chem(31))
    call add_model_scalar_field('species_rooh',  'ROOH',  isca_chem(32))
    call add_model_scalar_field('species_aacd',  'AACD',  isca_chem(33))
    call add_model_scalar_field('species_pacd',  'PACD',  isca_chem(34))
    call add_model_scalar_field('species_par',   'PAR',   isca_chem(35))
    call add_model_scalar_field('species_ror',   'ROR',   isca_chem(36))
    call add_model_scalar_field('species_eth',   'ETH',   isca_chem(37))
    call add_model_scalar_field('species_ole',   'OLE',   isca_chem(38))
    call add_model_scalar_field('species_iole',  'IOLE',  isca_chem(39))
    call add_model_scalar_field('species_isop',  'ISOP',  isca_chem(40))
    call add_model_scalar_field('species_ispd',  'ISPD',  isca_chem(41))
    call add_model_scalar_field('species_terp',  'TERP',  isca_chem(42))
    call add_model_scalar_field('species_tol',   'TOL',   isca_chem(43))
    call add_model_scalar_field('species_xyl',   'XYL',   isca_chem(44))
    call add_model_scalar_field('species_cres',  'CRES',  isca_chem(45))
    call add_model_scalar_field('species_to2',   'TO2',   isca_chem(46))
    call add_model_scalar_field('species_open',  'OPEN',  isca_chem(47))
    call add_model_scalar_field('species_cro',   'CRO',   isca_chem(48))
    call add_model_scalar_field('species_mgly',  'MGLY',  isca_chem(49))
    call add_model_scalar_field('species_so2',   'SO2',   isca_chem(50))
    call add_model_scalar_field('species_h2so4', 'H2SO4', isca_chem(51))
    call add_model_scalar_field('species_hco3',  'HCO3',  isca_chem(52))
    dmmk(1)=30.d-3     ! Molar mass NO
    dmmk(2)=46.d-3     ! Molar mass NO2
    dmmk(3)=48.d-3     ! Molar mass O3
    dmmk(4)=16.d-3     ! Molar mass O3P
    dmmk(5)=16.d-3     ! Molar mass O1D
    dmmk(6)=17.01d-3   ! Molar mass OH
    dmmk(7)=33.01d-3   ! Molar mass HO2
    dmmk(8)=34.01d-3   ! Molar mass H2O2
    dmmk(9)=62.01d-3   ! Molar mass NO3
    dmmk(10)=108.01d-3 ! Molar mass N2O5
    dmmk(11)=47.01d-3  ! Molar mass HONO
    dmmk(12)=63.01d-3  ! Molar mass HNO3
    dmmk(13)=79.01d-3  ! Molar mass HNO4
    dmmk(14)=28.01d-3  ! Molar mass CO
    dmmk(15)=30.03d-3  ! Molar mass HCHO
    dmmk(16)=44.05d-3  ! Molar mass ALD2
    dmmk(17)=75.04d-3  ! Molar mass C2O3
    dmmk(18)=121.05d-3 ! Molar mass PAN
    dmmk(19)=43.04d-3  ! Molar mass ALDX
    dmmk(20)=74.04d-3  ! Molar mass CXO3
    dmmk(21)=120.04d-3 ! Molar mass PANX
    dmmk(22)=47.03d-3  ! Molar mass XO2
    dmmk(23)=47.03d-3  ! Molar mass XO2N
    dmmk(24)=77.04d-3  ! Molar mass NTR
    dmmk(25)=46.07d-3  ! Molar mass ETOH
    dmmk(26)=16.04d-3  ! Molar mass CH4
    dmmk(27)=47.03d-3  ! Molar mass MEO2
    dmmk(28)=32.04d-3  ! Molar mass MEOH
    dmmk(29)=48.04d-3  ! Molar mass MEPX
    dmmk(30)=46.03d-3  ! Molar mass FACD
    dmmk(31)=30.07d-3  ! Molar mass ETHA
    dmmk(32)=47.03d-3  ! Molar mass ROOH
    dmmk(33)=60.05d-3  ! Molar mass AACD
    dmmk(34)=76.05d-3  ! Molar mass PACD
    dmmk(35)=15.03d-3  ! Molar mass PAR
    dmmk(36)=16.d-3    ! Molar mass ROR
    dmmk(37)=28.05d-3  ! Molar mass ETH
    dmmk(38)=27.05d-3  ! Molar mass OLE
    dmmk(39)=56.11d-3  ! Molar mass IOLE
    dmmk(40)=68.18d-3  ! Molar mass ISOP
    dmmk(41)=70.09d-3  ! Molar mass ISPD
    dmmk(42)=136.24d-3 ! Molar mass TERP
    dmmk(43)=92.14d-3  ! Molar mass TOL
    dmmk(44)=106.16d-3 ! Molar mass XYL
    dmmk(45)=108.14d-3 ! Molar mass CRES
    dmmk(46)=141.15d-3 ! Molar mass TO2
    dmmk(47)=48.04d-3  ! Molar mass OPEN
    dmmk(48)=108.14d-3 ! Molar mass CRO
    dmmk(49)=72.06d-3  ! Molar mass MGLY
    dmmk(50)=64.06d-3  ! Molar mass SO2
    dmmk(51)=98.08d-3  ! Molar mass H2SO4
    dmmk(52)=63.03d-3  ! Molar mass HCO3
    if (iaerosol.eq.1) then
      call add_model_scalar_field('species_hc8',      'HC8',      isca_chem(53))
      call add_model_scalar_field('species_api',      'API',      isca_chem(54))
      call add_model_scalar_field('species_lim',      'LIM',      isca_chem(55))
      call add_model_scalar_field('species_cvar01',   'CVARO1',   isca_chem(56))
      call add_model_scalar_field('species_cvar02',   'CVARO2',   isca_chem(57))
      call add_model_scalar_field('species_cvalk1',   'CVALK1',   isca_chem(58))
      call add_model_scalar_field('species_cvole1',   'CVOLE1',   isca_chem(59))
      call add_model_scalar_field('species_cvapi1',   'CVAPI1',   isca_chem(60))
      call add_model_scalar_field('species_cvapi2',   'CVAPI2',   isca_chem(61))
      call add_model_scalar_field('species_cvlim1',   'CVLIM1',   isca_chem(62))
      call add_model_scalar_field('species_cvlim2',   'CVLIM2',   isca_chem(63))
      call add_model_scalar_field('species_nh3',      'NH3',      isca_chem(64))
      call add_model_scalar_field('species_hcl',      'HCL',      isca_chem(65))
      call add_model_scalar_field('species_cvbibmp',  'CVBIBMP',  isc)
      call add_model_scalar_field('species_cvanclp',  'CVANCLP',  isc)
      call add_model_scalar_field('species_cvbiiso1', 'CVBIISO1', isc)
      call add_model_scalar_field('species_cvbiiso2', 'CVBIISO2', isc)
      dmmk(53)=114.0d0   ! Molar mass HC8
      dmmk(54)=136.0d0   ! Molar mass API
      dmmk(55)=136.0d0   ! Molar mass LIM
      dmmk(56)=150.0d0   ! Molar mass CVARO1
      dmmk(57)=150.0d0   ! Molar mass CVARO2
      dmmk(58)=140.0d0   ! Molar mass CVALK1
      dmmk(59)=140.0d0   ! Molar mass CVOLE1
      dmmk(60)=184.0d0   ! Molar mass CVAPI1
      dmmk(61)=184.0d0   ! Molar mass CVAPI2
      dmmk(62)=200.0d0   ! Molar mass CVLIM1
      dmmk(63)=200.0d0   ! Molar mass CVLIM2
      dmmk(64)=17.0d0    ! Molar mass NH3
      dmmk(65)=36.5d0    ! Molar mass HCL
    endif
    if (iaerosol.eq.1) then
      chempoint = (/ 64, 65, 59, 57, 3, 55, 61, 20, 56, 16, 23, 50,&
                     17, 51, 54, 60, 62, 13, 48, 58, 18, 63, 52, 46,&
                     4, 5, 53, 14, 22, 35, 6, 32, 49, 33, 47, 19, 34,&
                     36, 37, 44, 45, 39, 7, 8, 40, 15, 41, 43, 42, 9,&
                     10, 21, 11, 24, 27, 30, 31, 12, 38, 25, 26, 28, 29,&
                     1, 2 /)
    else
      chempoint = (/ 48, 52, 47, 43, 1, 42, 50, 17, 44, 9, 15, 38, 13, 37,&
                     41, 45, 51, 10, 35, 46, 14, 49, 39, 33, 2, 3, 40, 11,&
                     19, 20, 4, 21, 36, 22, 34, 16, 23, 24, 25, 31, 32, 26,&
                     5, 6, 27, 12, 28, 30, 29, 7, 8, 18 /)
    endif
 endif

endif

! Atmospheric aerosol chemistry
if (iaerosol.eq.1) then

  write (nfecra,*) "Atmospheric aerosol chemistry activated"

  ! logical unit and name of the chemical profiles file
  impmea=28
  ficmea='aerosols'

  ! 4 species not followed in gaseous scheme 3 but necessary to siream
  nespg_siream=nespg+4

  ! Verification
  if (ifilechemistry.ne.3) then
    write(nfecra,1004)
    call csexit (1)
  endif

  ! Filling the names array
  esp_siream=(/'MD    ', 'BC    ', 'NA    ', 'H2SO4 ', 'NH3   ', 'HNO3  ',&
 'HCL   ','ARO1  ', 'ARO2  ', 'ALK1  ', 'OLE1  ', 'API1  ', 'API2  ',&
 'LIM1  ', 'LIM2  ', 'ANCLP ', 'BIISO1', 'BIISO2',&
 'BIBMP ', 'POA   ', 'H2O   '/)

  f_esp_siream=(/'md    ', 'bc    ', 'na    ', 'h2so4 ', 'nh3   ', 'hno3  ',&
 'hcl   ','aro1  ', 'aro2  ', 'alk1  ', 'ole1  ', 'api1  ', 'api2  ',&
 'lim1  ', 'lim2  ', 'anclp ', 'biiso1', 'biiso2',&
 'bibmp ', 'poa   ', 'h2o   '/)

  do jj = 1, nesp_aer
    do ii = 1, nbin_aer
      write(nbin,"(i2)") ii
      name = 'species_'//trim(f_esp_siream(jj))//'_bin'//trim(adjustl(nbin))
      label = trim(esp_siream(jj))//'_bin'//trim(adjustl(nbin))
      call add_model_scalar_field(name, label,  isc)
    enddo
  enddo

  do ii = 1, nbin_aer
    write(nbin,"(i2)") ii
    name = 'species_'//'naero_bin'//trim(adjustl(nbin))
    label = 'Naero_bin'//trim(adjustl(nbin))
    call add_model_scalar_field(name, label,  isc)
  enddo

endif

!===============================================================================
! 4. Map to fields and GUI
!===============================================================================

n_chem_species = nespg
do isc = 1, nespg
  species_f_id(isc) = ivarfl(isca(isca_chem(isc)))
enddo

call cs_field_pointer_map_atmospheric(n_chem_species, species_f_id)

!===============================================================================
! 5. General field and physical properties
!===============================================================================

! Cp is constant
icp = -1

!----
! Formats
!----

#if defined(_CS_LANG_FR)

 1004 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    CHIMIE ATMOSPHERIQUE PARTICULAIRE DEMANDEE              ',/,&
'@                                                            ',/,&
'@  Lorsque le modele de chimie particulaire siream est utilise',/,&
'@   un schema gazeux complet(CB05) est automatiquement utilise',/,&
'@  L''utilisateur ne doit pas en specifier d''autre (ichemistry)',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usati1 (cs_user_parameters.f90)                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1004 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ATMOSPHERIC AEROSOL CHEMISTRY MODULE                    ',/,&
'@                                                            ',/,&
'@  When aerosol chemistry model siream is used               ',/,&
'@   a full gaseous scheme (CB05) is automatically used       ',/,&
'@  The user cannot specify any other scheme (ichemistry)     ',/,&
'@  Computation CAN NOT run.                                  ',/,&
'@                                                            ',/,&
'@  Check the input data given through the User Interface     ',/,&
'@   or in cs_user_parameters.f90.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! End
!----

return
end subroutine atvarp
