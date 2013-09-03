!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine cs_coal_varpos
!========================

!===============================================================================
!  FONCTION  :
!  ---------
!       INIT DES POSITIONS DES VARIABLES TRANSPORTEES POUR
!                COMBUSTION CHARBON PULVERISE
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use ihmpre
use cs_coal_incl
use field

!===============================================================================

implicit none

integer          icla,  is, icha, isc , is1
integer          f_id, itycat, ityloc, idim1, idim3
integer          keyccl, keydri, keyvis, keylbl, kscmin, kscmax
integer          iscdri, iopchr

logical          ilved, iprev, inoprv

character*80     f_name, name

!===============================================================================
! 0. Definitions for fields
!===============================================================================

! The itycat variable is used to define field categories. It is used in Fortran
! code with hard-coded values, but in the C API, those values are based on
! (much clearer) category mask definitions in cs_field.h.

itycat = FIELD_INTENSIVE + FIELD_VARIABLE  ! for most variables
ityloc = 1 ! variables defined on cells
idim1  = 1
idim3  = 3
ilved  = .false.   ! not interleaved by default
iprev  = .true.    ! variables have previous value
inoprv = .false.   ! variables have no previous value
iopchr = 1         ! Postprocessing level for variables

name = 'post_vis'
call field_get_key_id(name, keyvis)

name = 'label'
call field_get_key_id(name, keylbl)

! Key id of the coal scalar class
call field_get_key_id("scalar_class", keyccl)

! Key id for drift scalar
call field_get_key_id("drift_scalar_model", keydri)

! Key id for scamin and scamax
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)

!===============================================================================
! 1. DEFINITION DES POINTEURS
!===============================================================================

! ---> Variables propres a la suspension gaz - particules
is = 1
ihm   = iscapp(is)

! Activate the drift: 0 (no activation), 1 (activation)
iscdri = i_coal_drift

! Dispersed phase variables
!--------------------------

do icla = 1, nclacp
  is = 1+icla
  ! Field number of particle
  inp(icla) = iscapp(is)
  write(f_name,'(a5,i2.2)')'Np_CP' ,icla
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, rinfin)

  ! Scalar with drift: DO create additional mass flux
  if (i_coal_drift.eq.1) then
    iscdri = ibset(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)

  is = 1+1*nclacp+icla
  ! Field Xch
  ixch(icla)= iscapp(is)
  write(f_name,'(a6,i2.2)')'Xch_CP' ,icla
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_coal_drift.eq.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)

  is = 1+2*nclacp+icla
  ! Field xck
  ixck(icla) = iscapp(is)
  write(f_name,'(a6,i2.2)')'xck_cp' ,icla
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_coal_drift.eq.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)

  ! With drying
  if (ippmod(iccoal).eq.1) then
    is = 1+3*nclacp+icla
    ! Field Xwt
    ixwt(icla) = iscapp(is)
    write(f_name,'(a6,i2.2)')'Xwt_CP', icla
    call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
    call field_set_key_str(f_id, keylbl, f_name)
    ! Set the index of the scalar class in the field structure
    call field_set_key_int(f_id, keyccl, icla)
    ! Set min and max clipping
    call field_set_key_double(f_id, kscmin, 0.d0)
    call field_set_key_double(f_id, kscmax, 1.d0)

    ! Scalar with drift: BUT Do NOT create additional mass flux
    if (i_coal_drift.eq.1) then
      iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
      call field_set_key_int(f_id, keydri, iscdri)
    endif

    ! For post-processing
    call field_set_key_int(f_id, keyvis, iopchr)

    is = 1+4*nclacp+icla
    ! Field h2
    ih2(icla) = iscapp(is)
    write(f_name,'(a6,i2.2)')'Ent_CP' ,icla
    call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
    call field_set_key_str(f_id, keylbl, f_name)
    ! Set the index of the scalar class in the field structure
    call field_set_key_int(f_id, keyccl, icla)
    ! Set min and max clipping
    call field_set_key_double(f_id, kscmin, -grand)
    call field_set_key_double(f_id, kscmax, grand)

    ! Scalar with drift: BUT Do NOT create additional mass flux
    if (i_coal_drift.eq.1) then
      iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
      call field_set_key_int(f_id, keydri, iscdri)
    endif

    ! For post-processing
    call field_set_key_int(f_id, keyvis, iopchr)

    if (i_coal_drift.eq.1) then
      is = 1+5*nclacp+icla
      ! Field X_Age
      iagecp_temp(icla) = iscapp(is)
      write(f_name,'(a8,i2.2)')'X_Age_CP' ,icla
      call field_create(f_name, itycat, ityloc, idim1, ilved, iprev, f_id)
      call field_set_key_str(f_id, keylbl, f_name)
      ! Set the index of the scalar class in the field structure
      call field_set_key_int(f_id, keyccl, icla)
      ! Set min and max clipping
      call field_set_key_double(f_id, kscmin, 0.d0 )
      call field_set_key_double(f_id, kscmax, grand)

      ! Scalar with drift: BUT Do NOT create additional mass flux
      iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
      call field_set_key_int(f_id, keydri, iscdri)

      ! For post-processing
      call field_set_key_int(f_id, keyvis, iopchr)
    endif
  else
    is = 1+3*nclacp+icla
    ! Field h2
    ih2(icla) = iscapp(is)
    write(f_name,'(a6,i2.2)')'Ent_CP' ,icla
    call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
    call field_set_key_str(f_id, keylbl, f_name)
    ! Set the index of the scalar class in the field structure
    call field_set_key_int(f_id, keyccl, icla)
    ! Set min and max clipping
    call field_set_key_double(f_id, kscmin, -grand)
    call field_set_key_double(f_id, kscmax, grand)

    ! Scalar with drift: BUT Do NOT create additional mass flux
    if (i_coal_drift.eq.1) then
      iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
      call field_set_key_int(f_id, keydri, iscdri)
    endif

    ! For post-processing
    call field_set_key_int(f_id, keyvis, iopchr)

    if (i_coal_drift.eq.1) then
      is = 1+4*nclacp+icla
      ! Field X_Age
      iagecp_temp(icla) = iscapp(is)
      write(f_name,'(a8,i2.2)')'X_Age_CP' ,icla
      call field_create(f_name, itycat, ityloc, idim1, ilved, iprev, f_id)
      call field_set_key_str(f_id, keylbl, f_name)
      ! Set the index of the scalar class in the field structure
      call field_set_key_int(f_id, keyccl, icla)

      ! Scalar with drift: BUT Do NOT create additional mass flux
      iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
      call field_set_key_int(f_id, keydri, iscdri)

      ! For post-processing
      call field_set_key_int(f_id, keyvis, iopchr)
    endif
  endif
enddo

! Continuous phase variables
!---------------------------

! Class index for the gas scalars
icla = -1

! Matiere volatiles legeres (F8) et lourdes (F9)
is1 = is
do icha = 1, ncharb
  is = is1+icha
  ! Field Fr_mv1
  if1m(icha)  = iscapp(is)
  write(f_name,'(a6,i2.2)')'Fr_mv1' ,icha
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! The first gas salar contains the drift flux, the others
  ! Scalar with drift: DO create additional mass flux
  if (i_coal_drift.eq.1) then
    if (icha.eq.1) then
      iscdri = ibset(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    ! Scalar with drift: BUT Do NOT create additional mass flux
    else
      iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    endif
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)

  is          = is1+ncharb+icha
  ! Field Fr_mv2
  if2m(icha)  = iscapp(is)
  write(f_name,'(a6,i2.2)')'Fr_mv2' ,icha
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_coal_drift.eq.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)
enddo

! Oxydant 2
if (noxyd .ge. 2) then
  is = is+1
  ! Field FR_OXYD2
  if4m  = iscapp(is)
  f_name = 'FR_OXYD2'
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_coal_drift.eq.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)
endif

! Oxydant 3
if (noxyd .ge. 3) then
  is = is+1
  ! Field FR_OXYD3
  if5m  = iscapp(is)
  f_name = 'FR_OXYD3'
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_coal_drift.eq.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)
endif

! Humidite
if (ippmod(iccoal).eq.1) then
  is = is+1
  ! Field FR_H20
  if6m  = iscapp(is)
  f_name = 'FR_H2O'
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_coal_drift.eq.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)
endif

! Produits de la combustion du coke par O2
is = is+1
! Field FR_HET_O2
if7m  = iscapp(is)
f_name = 'FR_HET_O2'
call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
call field_set_key_str(f_id, keylbl, f_name)
! Set the index of the scalar class in the field structure
call field_set_key_int(f_id, keyccl, icla)
! Set min and max clipping
call field_set_key_double(f_id, kscmin, 0.d0)
call field_set_key_double(f_id, kscmax, 1.d0)

! Scalar with drift: BUT Do NOT create additional mass flux
if (i_coal_drift.eq.1) then
  iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
  call field_set_key_int(f_id, keydri, iscdri)
endif

! For post-processing
call field_set_key_int(f_id, keyvis, iopchr)

! Produits de la combustion du coke par CO2
if (ihtco2.eq.1) then
  is = is+1
  ! Field FR_HET_CO2
  if8m  = iscapp(is)
  f_name = 'FR_HET_CO2'
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_coal_drift.eq.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)
endif

! Produits de la combustion du coke par H2O
if (ihth2o.eq.1) then
  is = is+1
  ! Field FR_HET_CO2
  if9m  = iscapp(is)
  f_name = 'FR_HET_H2O'
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_coal_drift.eq.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)
endif

! variance
is = is+1
! Field Var_F1F2
ifvp2m = iscapp(is)
f_name = 'Var_F1F2'
call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
call field_set_key_str(f_id, keylbl, f_name)
! Set the index of the scalar class in the field structure
call field_set_key_int(f_id, keyccl, icla)
! Set min and max clipping
call field_set_key_double(f_id, kscmin, 0.d0)
call field_set_key_double(f_id, kscmax, 0.25d0)

! Scalar with drift: BUT Do NOT create additional mass flux
if (i_coal_drift.eq.1) then
  iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
  call field_set_key_int(f_id, keydri, iscdri)
endif

! For post-processing
call field_set_key_int(f_id, keyvis, iopchr)

! Transport du CO ou du CO2
if (ieqco2.ge.1) then
  is    = is+1
  ! Field
  iyco2 = iscapp(is)
  f_name = 'FR_CO2'
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_coal_drift.eq.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)
endif

! Transport du NOx : HCN, NOx et Tair
if (ieqnox.eq.1) then
  is     = is+1
  ! Field
  iyhcn  = iscapp(is)
  f_name = 'FR_HCN'
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_coal_drift.eq.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)

  ! On ajoute le NH3 comme variable transportee
  is     = is+1
  ! Field
  iynh3  = iscapp(is)
  f_name =  'FR_NH3'
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_coal_drift.eq.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)

  is     = is+1
  ! Field
  iyno   = iscapp(is)
  f_name =  'FR_NO'
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_coal_drift.eq.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)

  is     = is+1
  ! Field
  ihox   = iscapp(is)
  f_name =  'Enth_Ox'
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin,-grand)
  call field_set_key_double(f_id, kscmax, grand)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  if (i_coal_drift.eq.1) then
    iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
    call field_set_key_int(f_id, keydri, iscdri)
  endif

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)
endif

if (i_coal_drift.eq.1) then
  is     = is+1
  ! Field
  iaggas_temp = iscapp(is)
  f_name = 'X_Age_Gas'
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0 )
  call field_set_key_double(f_id, kscmax, grand)

  ! Scalar with drift: BUT Do NOT create additional mass flux
  iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
  call field_set_key_int(f_id, keydri, iscdri)

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)
endif
!   - Interface Code_Saturne
!     ======================
!     Construction de l'indirection entre la numerotation du noyau et XML

if (iihmpr.eq.1) then
!
  call uicpsc (ncharb, nclacp, noxyd, ippmod,               &
               iccoal, ieqnox, ieqco2, ihtco2,              &
               ihth2o, ihm, inp, ixch, ixck, ixwt, ih2,     &
               if1m, if2m, if4m, if5m, if6m,                &
               if7m, if8m, ifvp2m, iyco2, if9m,     &
               iyhcn, iyno, ihox, iynh3)
endif

!===============================================================================
! 2. PROPRIETES PHYSIQUES
!    A RENSEIGNER OBLIGATOIREMENT (sinon pb dans varpos)
!    - PROPRES AUX SCALAIRES   : IVISLS, ISCAVR
!      Rq : pas de variance associee a un scalaire dans notre cas
!    - PROPRES A LA SUSPENSION : ICP
!===============================================================================

do isc = 1, nscapp

  if ( iscavr(iscapp(isc)).le.0 ) then

! ---- Viscosite dynamique de reference relative au scalaire
!      ISCAPP(ISC)
    ivisls(iscapp(isc)) = 0

  endif

enddo

! ---- Bien que l'on soit en enthalpie on conserve un CP constant

icp    = 0

return
end subroutine
