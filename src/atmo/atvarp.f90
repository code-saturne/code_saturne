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
use atincl
use atchem
use field
use sshaerosol
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

integer        ii, jj, isc, f_id
integer        kscmin, kscmax

integer(c_int) :: n_chem_species
integer(c_int), dimension(:), allocatable :: species_f_id

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
! 1. Model information
!===============================================================================


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
  call add_model_scalar_field('ym_water', 'Ym water', iymw)
  f_id = ivarfl(isca(iymw))
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
if (iaerosol.eq.CS_ATMO_AEROSOL_SSH) ichemistry = 4
! if a chemical scheme is solved, a concentration profiles
! file  named 'chemistry' should be used
if (ichemistry.ge.1) ifilechemistry = ichemistry
if (nogaseouschemistry .and. iaerosol.eq.CS_ATMO_AEROSOL_OFF) ichemistry = 0

if (ifilechemistry.ge.1) then

  ! Set the name of the chemical profiles file
  call cs_atmo_set_chem_conc_file_name('chemistry'//c_null_char)

  ! Initialization of the chemical scheme
  ! quasi steady equilibrium NOx scheme with 4 species and 5 reactions
  if (ifilechemistry.eq.1) then
    nrg = 5
    nespg = 4

    ! Map isca_chem, dmmk, chempoint and allocate it if needed
    call init_chemistry_pointers()

    call add_model_scalar_field('species_no',  'NO',  isca_chem(1))
    call add_model_scalar_field('species_no2', 'NO2', isca_chem(2))
    call add_model_scalar_field('species_o3',  'O3',  isca_chem(3))
    call add_model_scalar_field('species_o3p', 'O3P', isca_chem(4))
    dmmk(1)=30.d0  ! Molar mass (g/mol) NO
    dmmk(2)=46.d0  ! Molar mass (g/mol) NO2
    dmmk(3)=48.d0  ! Molar mass (g/mol) O3
    dmmk(4)=16.d0  ! Molar mass (g/mol) O3P
    chempoint = (/ 4, 3, 2, 1 /)

  ! scheme with 20 species and 34 reactions! Note pas de COV
  else if (ifilechemistry.eq.2) then
    nrg = 34
    nespg = 20

    ! Map isca_chem, dmmk, chempoint and allocate it if needed
    call init_chemistry_pointers()

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
    dmmk(1)=30.d0      ! Molar mass (g/mol) NO
    dmmk(2)=46.d0      ! Molar mass (g/mol) NO2
    dmmk(3)=48.d0      ! Molar mass (g/mol) O3
    dmmk(4)=16.d0      ! Molar mass (g/mol) O3P
    dmmk(5)=16.d0      ! Molar mass (g/mol) O1D
    dmmk(6)=17.01d0    ! Molar mass (g/mol) OH
    dmmk(7)=33.01d0    ! Molar mass (g/mol) HO2
    dmmk(8)=34.01d0    ! Molar mass (g/mol) H2O2
    dmmk(9)=62.01d0    ! Molar mass (g/mol) NO3
    dmmk(10)=108.01d0  ! Molar mass (g/mol) N2O5
    dmmk(11)=47.01d0   ! Molar mass (g/mol) HONO
    dmmk(12)=63.01d0   ! Molar mass (g/mol) HNO3
    dmmk(13)=28.01d0   ! Molar mass (g/mol) CO
    dmmk(14)=30.03d0   ! Molar mass (g/mol) HCHO
    dmmk(15)=44.05d0   ! Molar mass (g/mol) ALD2
    dmmk(16)=75.04d0   ! Molar mass (g/mol) C2O3
    dmmk(17)=121.05d0  ! Molar mass (g/mol) PAN
    dmmk(18)=47.03d0   ! Molar mass (g/mol) XO2
    dmmk(19)=64.06d0   ! Molar mass (g/mol) SO2
    dmmk(20)=98.08d0   ! Molar mass (g/mol) H2SO4
    chempoint = (/ 20, 19, 16, 17, 2, 15, 14, 3, 18, 7, 8, 9, 4, &
                   10, 1, 12, 11, 13, 5, 6 /)

  ! scheme CB05 with 52 species and 155 reactions
  else if (ifilechemistry.eq.3) then
    nrg = 155
    nespg = 52

    ! Map isca_chem, dmmk, chempoint and allocate it if needed
    call init_chemistry_pointers()

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
    dmmk(1)=30.d0      ! Molar mass (g/mol) NO
    dmmk(2)=46.d0      ! Molar mass (g/mol) NO2
    dmmk(3)=48.d0      ! Molar mass (g/mol) O3
    dmmk(4)=16.d0      ! Molar mass (g/mol) O3P
    dmmk(5)=16.d0      ! Molar mass (g/mol) O1D
    dmmk(6)=17.01d0    ! Molar mass (g/mol) OH
    dmmk(7)=33.01d0    ! Molar mass (g/mol) HO2
    dmmk(8)=34.01d0    ! Molar mass (g/mol) H2O2
    dmmk(9)=62.01d0    ! Molar mass (g/mol) NO3
    dmmk(10)=108.01d0  ! Molar mass (g/mol) N2O5
    dmmk(11)=47.01d0   ! Molar mass (g/mol) HONO
    dmmk(12)=63.01d0   ! Molar mass (g/mol) HNO3
    dmmk(13)=79.01d0   ! Molar mass (g/mol) HNO4
    dmmk(14)=28.01d0   ! Molar mass (g/mol) CO
    dmmk(15)=30.03d0   ! Molar mass (g/mol) HCHO
    dmmk(16)=44.05d0   ! Molar mass (g/mol) ALD2
    dmmk(17)=75.04d0   ! Molar mass (g/mol) C2O3
    dmmk(18)=121.05d0  ! Molar mass (g/mol) PAN
    dmmk(19)=43.04d0   ! Molar mass (g/mol) ALDX
    dmmk(20)=74.04d0   ! Molar mass (g/mol) CXO3
    dmmk(21)=120.04d0  ! Molar mass (g/mol) PANX
    dmmk(22)=47.03d0   ! Molar mass (g/mol) XO2
    dmmk(23)=47.03d0   ! Molar mass (g/mol) XO2N
    dmmk(24)=77.04d0   ! Molar mass (g/mol) NTR
    dmmk(25)=46.07d0   ! Molar mass (g/mol) ETOH
    dmmk(26)=16.04d0   ! Molar mass (g/mol) CH4
    dmmk(27)=47.03d0   ! Molar mass (g/mol) MEO2
    dmmk(28)=32.04d0   ! Molar mass (g/mol) MEOH
    dmmk(29)=48.04d0   ! Molar mass (g/mol) MEPX
    dmmk(30)=46.03d0   ! Molar mass (g/mol) FACD
    dmmk(31)=30.07d0   ! Molar mass (g/mol) ETHA
    dmmk(32)=47.03d0   ! Molar mass (g/mol) ROOH
    dmmk(33)=60.05d0   ! Molar mass (g/mol) AACD
    dmmk(34)=76.05d0   ! Molar mass (g/mol) PACD
    dmmk(35)=15.03d0   ! Molar mass (g/mol) PAR
    dmmk(36)=16.d0     ! Molar mass (g/mol) ROR
    dmmk(37)=28.05d0   ! Molar mass (g/mol) ETH
    dmmk(38)=27.05d0   ! Molar mass (g/mol) OLE
    dmmk(39)=56.11d0   ! Molar mass (g/mol) IOLE
    dmmk(40)=68.18d0   ! Molar mass (g/mol) ISOP
    dmmk(41)=70.09d0   ! Molar mass (g/mol) ISPD
    dmmk(42)=136.24d0  ! Molar mass (g/mol) TERP
    dmmk(43)=92.14d0   ! Molar mass (g/mol) TOL
    dmmk(44)=106.16d0  ! Molar mass (g/mol) XYL
    dmmk(45)=108.14d0  ! Molar mass (g/mol) CRES
    dmmk(46)=141.15d0  ! Molar mass (g/mol) TO2
    dmmk(47)=48.04d0   ! Molar mass (g/mol) OPEN
    dmmk(48)=108.14d0  ! Molar mass (g/mol) CRO
    dmmk(49)=72.06d0   ! Molar mass (g/mol) MGLY
    dmmk(50)=64.06d0   ! Molar mass (g/mol) SO2
    dmmk(51)=98.08d0   ! Molar mass (g/mol) H2SO4
    dmmk(52)=63.03d0   ! Molar mass (g/mol) HCO3
    chempoint = (/ 48, 52, 47, 43, 1, 42, 50, 17, 44, 9, 15, 38, 13, 37,&
                   41, 45, 51, 10, 35, 46, 14, 49, 39, 33, 2, 3, 40, 11,&
                   19, 20, 4, 21, 36, 22, 34, 16, 23, 24, 25, 31, 32, 26,&
                   5, 6, 27, 12, 28, 30, 29, 7, 8, 18 /)

  ! User defined chemistry using SPACK file and routines
  else if (ifilechemistry.eq.4) then

    ! This function read the number of species, their molar mass
    ! and creates variables
    call cs_atmo_declare_chem_from_spack()

    ! Read the number of reactions
    call ssh_dimensions(ii, nrg, jj)

    ! Verification
    if (ii.ne.nespg) then
      write(nfecra,1003)
      call csexit (1)
    endif

    ! Map isca_chem, dmmk, chempoint and allocate it if needed
    call init_chemistry_pointers()

  endif

  ! Finish initialization of C chemistry when SPACK was not used
  if (ifilechemistry.ne.4) then
    call cs_atmo_chem_init_c_chemistry()
  endif

endif

! Atmospheric aerosol chemistry
if (iaerosol.ne.CS_ATMO_AEROSOL_OFF) then

  ! Set the name of the chemical profiles file
  call cs_atmo_set_aero_conc_file_name('aerosols')

  ! Verification
  if (ifilechemistry.ne.4) then
    write(nfecra,1004)
    call csexit (1)
  endif

  ! Load shared library
  ! Initialise external aerosol code
  ! Create variables
  call cs_atmo_aerosol_initialize()

  ! Remap pointers following bft_realloc when initializing aerosols
  call init_aerosol_pointers()

endif

! Set clippings for gas and aerosol species
if (ichemistry.ge.1) then
  do ii = 1, nespg
    f_id = ivarfl(isca(isca_chem(ii)))
    call field_set_key_double(f_id, kscmin, 0.d0)
  enddo
endif
if (iaerosol.ne.CS_ATMO_AEROSOL_OFF) then
  do ii = nespg + 1, nespg + n_aer * (nlayer_aer + 1)
    f_id = ivarfl(isca(isca_chem(ii)))
    call field_set_key_double(f_id, kscmin, 0.d0)
  enddo
  ! Allow large aerosol numbers
  do ii = nespg + n_aer * nlayer_aer + 1, nespg + n_aer * (nlayer_aer + 1)
    f_id = ivarfl(isca(isca_chem(ii)))
    call field_set_key_double(f_id, kscmax, 1.d40)
  enddo
endif

!===============================================================================
! 4. Map to fields and GUI
!===============================================================================

allocate(species_f_id(nespg))

n_chem_species = nespg
do isc = 1, nespg
  species_f_id(isc) = ivarfl(isca(isca_chem(isc)))
enddo

call cs_field_pointer_map_atmospheric(n_chem_species, species_f_id)

deallocate(species_f_id)

!===============================================================================
! 5. General field and physical properties
!===============================================================================

! Cp is constant
icp = -1

!----
! Formats
!----

 1003 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ATMOSPHERIC CHEMISTRY FROM SPACK                        ',/,&
'@                                                            ',/,&
'@  The number of gaseous species read from the SPACK file    ',/,&
'@  is not equal to the one read in the SPACK source file     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1004 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ATMOSPHERIC AEROSOL CHEMISTRY MODULE                    ',/,&
'@                                                            ',/,&
'@  When aerosol chemistry model is used                      ',/,&
'@   a full gaseous scheme (CB05) is automatically used       ',/,&
'@  The user cannot specify any other scheme (ichemistry)     ',/,&
'@  Computation CAN NOT run.                                  ',/,&
'@                                                            ',/,&
'@  Check the input data given through the User Interface     ',/,&
'@   or in cs_user_parameters.f90.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

return
end subroutine atvarp
