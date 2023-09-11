!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

subroutine ecrava

!===============================================================================
! Purpose:
! --------

! Write main and auxiliary restart files

!-------------------------------------------------------------------------------
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

use, intrinsic :: iso_c_binding

use paramx
use dimens, only: nvar
use numvar
use cstphy
use entsor
use pointe
use optcal
use albase
use alaste
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use ppcpfu
use cplsat
use field
use atincl, only: init_at_chem
use atchem, only: ichemistry
use sshaerosol, only: iaerosol, CS_ATMO_AEROSOL_OFF
use mesh
use cs_c_bindings
use cs_nz_condensation, only:nfbpcd, izzftcd, nztag1d, ztpar, ifbpcd
use cs_nz_tagmr, only:znmurx, znmur, ztmur

!===============================================================================

implicit none

! Arguments

! Local variables

character        rubriq*64,car2*2,car54*54
character        cindfc*2,cindfl*4
character        ficsui*32
integer          f_id, t_id, ivar
integer          icha
integer          ii    , ivers
integer          itysup, nbval
integer          nfmtsc, nfmtfl, nfmtch, nfmtcl
integer          ilecec, iecr
integer          ifac
integer          iz, kk
integer          ival(1)
integer          key_t_ext_id
double precision rval(1)

type(c_ptr) :: rp

double precision, allocatable, dimension(:,:) :: tmurbf
double precision, allocatable, dimension(:) :: tparbf
double precision, pointer, dimension(:) :: dt_s

type(var_cal_opt) :: vcopt

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_ale_restart_write(r)  &
    bind(C, name='cs_ale_restart_write')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), value :: r
  end subroutine cs_ale_restart_write

  subroutine cs_mobile_structures_restart_write(r)  &
    bind(C, name='cs_mobile_structures_restart_write')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), value :: r
  end subroutine cs_mobile_structures_restart_write

  ! Interface to C function writing notebook variables

  subroutine cs_restart_write_notebook_variables(r)  &
    bind(C, name='cs_restart_write_notebook_variables')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), value :: r
  end subroutine cs_restart_write_notebook_variables

end interface

!===============================================================================
!     A noter :
!        Lorsque qu'il est necessaire d'utiliser un ordre implicite
!        de rangement des variables, on a choisi :
!          U, V, W,
!          P,
!          turbulence
!          scalaires

!          avec turbulence = k, epsilon
!                       ou   R11, R22, R33, R12, R13, R23, epsilon
!                       ou   k, epsilon, phi, f_barre
!                       ou   k, omega

!        Ceci est par exemple utilise pour relier les flux de masse aux
!        variables

!===============================================================================
! 0. Initialisation
!===============================================================================

! Time extrapolation?
call field_get_key_id("time_extrapolated", key_t_ext_id)

!===============================================================================
! 1. VERIFICATIONS DE BASE ET CODAGE DES CHAINES DE CARACTERES
!===============================================================================

!  --->  On code en chaine le numero des phases et scalaires
!        ----------------------------------------------------

!     Nombre de scalaires, de flux, et de charbons
!       max pour les formats choisis
nfmtsc = 9999
nfmtfl = 9999
nfmtch = 99
nfmtcl = 9999

!     Indefini (on met qqch de different de lecamo (pour generer une
!       erreur a la lecture)
cindfc = 'XX'
cindfl = 'XXXX'

!     Verifications pour les formats et les numeros
!       de scalaire en chaine.
!     Avertissement (certaines infos passent a la trappe)
if (nscamx.gt.nfmtsc) then
  write(nfecra,7001)nfmtsc,nscamx
endif
if (nvarmx.gt.nfmtfl) then
  write(nfecra,7002)nfmtfl,nvarmx
endif
if (ncharm.gt.nfmtch) then
  write(nfecra,7004)nfmtch,ncharm
endif
if (ncpcmx.gt.nfmtcl) then
  write(nfecra,7005)nfmtcl,ncpcmx
endif

!===============================================================================
! 2. OUVERTURE FICHIER SUITE DE BASE
!===============================================================================
! ILECEC = 2 : ecriture

write(nfecra,1000)

ilecec = 2

ficsui = 'main.csc'
call restart_create(ficsui, '', 1, rp)

!===============================================================================
! 3. ECRITURE FICHIER SUITE DE BASE
!===============================================================================

write(nfecra,1100)

! Restart version (for version x.y.z, xxyyzz)

ivers  = 400000
itysup = 0
nbval  = 1
ival(1) = ivers
rubriq = 'code_saturne:checkpoint:main:version'
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

! Main field metadata

call restart_write_field_info(rp)

! 3.1 OPTIONS (Celles servant a donner le nombre de tableaux a lire)
!============================================================================
! Remarque : ces variables ne sont pas toutes utilies pour la relecture
!            en revanche, elles peuvent completer les infos au sein
!            du fichier suite

!  ---> Nombre de pas de temps, instant precedent
rubriq = 'nbre_pas_de_temps'
itysup = 0
nbval  = 1
ival(1) = ntcabs
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

rubriq = 'instant_precedent'
itysup = 0
nbval  = 1
rval(1) = ttcabs
call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

!  ---> Modeles de turbulence
rubriq = 'turbulence_model'
itysup = 0
nbval  = 1
ival(1) = iturb
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

!  ---> Methode ALE
rubriq = 'methode_ALE'
itysup = 0
nbval  = 1
ival(1) = iale
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

!  ---> VOF
rubriq = 'vof'
itysup = 0
nbval  = 1
ival(1) = ivofmt
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

call turbomachinery_restart_write(rp)

if (ichemistry.gt.0.or.iaerosol.ne.CS_ATMO_AEROSOL_OFF) then
  rubriq = 'atmospheric_chem'
  itysup = 0
  nbval  = 1
  ival(1) = init_at_chem
  call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)
endif

car54 =' End writing the options                              '
write(nfecra,1110) car54

! 3.2 VARIABLES "PRINCIPALES"
!============================

call restart_write_variables(rp, 0)

f_id = -1
do ivar = 1, nvar
  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
  if (vcopt%ibdtso.gt.1) then
    if (f_id.ne.ivarfl(ivar)) then
      f_id = ivarfl(ivar)
      do t_id = 1, vcopt%ibdtso - 1
        call restart_write_field_vals(rp, f_id, t_id)
      enddo
    endif
  endif
enddo

call restart_write_fields(rp, RESTART_MAIN)

! 3.3 Notebook variables
!================================

call cs_restart_write_notebook_variables(rp)

!===============================================================================
! 4. FERMETURE FICHIER SUITE DE BASE
!===============================================================================

! Fermeture du fichier suite principal
call restart_destroy(rp)

write(nfecra,1200)

!===============================================================================
! 5. ECRITURE FICHIER SUITE AUXILIAIRE
!===============================================================================

!     Si  l'ecriture du fichier suite auxiliaire est demandee
if (iecaux.eq.1) then

! 5.0. OUVERTURE FICHIER SUITE AUXILIAIRE
!================================================

  write(nfecra,2000)

  ilecec = 2
  ficsui = 'auxiliary.csc'
  call restart_create(ficsui, '', 1, rp)

  write(nfecra,1100)

! Restart version (for version x.y.z, xxyyzz)

  ivers  = 400000
  itysup = 0
  nbval  = 1
  ival(1) = ivers
  rubriq = 'code_saturne:checkpoint:auxiliary:version'
  call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

! 5.1 DIMENSIONS : les dimensions geometriques sont ecrites
!===============   automatiquement lors de l'ouverture du fichier

!  ---> Indicateur de pas de temps variable
  rubriq = 'indic_dt_variable'
  itysup = 0
  nbval  = 1
  ival(1) = idtvar
  call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

  rubriq = 'methode_ALE'
  itysup = 0
  nbval  = 1
  ival(1) = iale
  call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

  !  ---> VOF
  rubriq = 'vof'
  itysup = 0
  nbval  = 1
  ival(1) = ivofmt
  call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

  car54 =' End writing the dimensions and options               '
  write(nfecra,1110)car54

! 5.3 ECRITURE DES VARIABLES
!===================================

! --->  Proprietes physiques

  !     Point de reference pour la pression totale
  !     On n'ecrit que si XYZP0 a ete specifie par l'utilisateur ou
  !       calcule a partir de faces de sorties ou de Dirichlet
  if (ixyzp0.eq.1) then
    rubriq = 'ref_presstot01'
    itysup = 0
    nbval  = 3
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,xyzp0)
  endif

  ! The physical variables here below are required for the low-Mach algorithm

  if (idilat.eq.3.or.ipthrm.eq.1) then

    !the reference density updated with the low-Mach algorithm
    rubriq = 'ro001'
    itysup = 0
    nbval  = 1
    rval(1) = ro0
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

    ! the thermodynamic pressure for the previous time step
    rubriq = 'pther01'
    itysup = 0
    nbval  = 1
    rval(1) = pther
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)
  endif

  call restart_write_linked_fields(rp, "diffusivity_id", iecr)

  car54 =' End writing the physical properties                  '
  write(nfecra,1110)car54

! ---> Pas de temps

  call field_get_id('dt', f_id)
  if (idtvar.eq.1) then
    call field_get_val_s(f_id, dt_s)
    rubriq = 'dt_variable_temps'
    itysup = 0
    nbval  = 1
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,dt_s)
  endif

  car54 =' End writing the time step                            '
  write(nfecra,1110)car54

  ! Mass fluxes

  call restart_write_linked_fields(rp, "inner_mass_flux_id", iecr)
  call restart_write_linked_fields(rp, "boundary_mass_flux_id", iecr)

  ! Boundary condition coefficients

  call restart_write_bc_coeffs(rp)

! ---> Termes sources
!      Lorsqu'ils sont extrapoles (pour les versions elec, voir plus bas)

  call restart_write_linked_fields(rp, "source_term_prev_id", iecr)

  if (iecr.ne.0) then
    car54=' End writing the source terms                         '
    write(nfecra,1110)car54
  endif

! ---> Moyennes (cumuls)

  call time_moment_restart_write(rp)

  ! if (iilagr .gt. 0 .and. istala .gt. 0) then
  !    call lagr_moment_restart_write(rp)
  ! endif

! ---> Distance a la paroi
!      On pourra ecrire ici la distance a la paroi

  iecr = 0

  if (iecr.ne.0) then
    car54=' End writing the wall distance                        '
    write(nfecra,1110)car54
  endif

!----------------------------------------------------------------
! ---> Wall temperature associated to the condensation model
!      with or without the 1D thermal model tag1D
!----------------------------------------------------------------

  if (icondb.eq.0) then

    if (nztag1d.eq.1) then

      ! Wall temperature array in the thickness for the 1D thermal model
      allocate(tmurbf(nfabor,znmurx))

      tmurbf(:,:) = 0.d0

      rubriq = 'tmur_bf_prev'
      itysup = 3
      nbval  = znmurx

      do ii = 1, nfbpcd
        ifac= ifbpcd(ii) + 1 ! C numbering
        iz  = izzftcd(ii) + 1 ! C numbering
        do kk = 1, znmur(iz)
          tmurbf(ifac, kk) = ztmur(ii,kk)
        enddo
      enddo

      call restart_write_section_real_t(rp,rubriq,itysup,nbval,tmurbf)

      ! Free memory
      deallocate(tmurbf)
    else
      ! Wall temperature array in the thickness for the 1D thermal model
      allocate(tparbf(nfabor))

      tparbf(:) = 0.d0

      rubriq = 'tpar_bf_prev'
      itysup = 3
      nbval  = 1

      do ii = 1, nfbpcd
        ifac= ifbpcd(ii) + 1 ! C numbering
        iz  = izzftcd(ii) + 1 ! C numbering
        tparbf(ifac) = ztpar(iz)
      enddo
      call restart_write_section_real_t(rp,rubriq,itysup,nbval,tparbf)

      ! Free memory
      deallocate(tparbf)
    endif

  endif

! ---> Methode ALE

  if (iale.ge.1) then

    call cs_ale_restart_write(rp)
    call cs_mobile_structures_restart_write(rp)

    car54=' End writing the ALE data              '
    write(nfecra,1110)car54

  endif

! ---> Grandeurs complementaires pour la combustion gaz

!     Modele COD3P :
!     ============

  if ( ippmod(icod3p).ge.0 ) then

    rubriq = 'hinfue_cod3p'
    itysup = 0
    nbval  = 1
    rval(1) = hinfue
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

    rubriq = 'hinoxy_cod3p'
    itysup = 0
    nbval  = 1
    rval(1) = hinoxy
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

    rubriq = 'tinfue_cod3p'
    itysup = 0
    nbval  = 1
    rval(1) = tinfue
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

    rubriq = 'tinoxy_cod3p'
    itysup = 0
    nbval  = 1
    rval(1) = tinoxy
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

!       Numero des zones
    itysup = 3
    nbval  = 1
    rubriq = 'num_zone_fb_cod3p'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,izfppp)

!       Entree Fuel (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientfu_zone_bord_cod3p'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ientfu)

!       Entree oxydant (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientox_zone_bord_cod3p'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ientox)

    car54=' End writing combustion information (COD3P)         '
    write(nfecra,1110)car54

  endif

!     Modele SLFM :
!     ============

  if ( ippmod(islfm).ge.0 ) then

    rubriq = 'hinfue_slfm'
    itysup = 0
    nbval  = 1
    rval(1) = hinfue
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

    rubriq = 'hinoxy_slfm'
    itysup = 0
    nbval  = 1
    rval(1) = hinoxy
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

    rubriq = 'tinfue_slfm'
    itysup = 0
    nbval  = 1
    rval(1) = tinfue
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

    rubriq = 'tinoxy_slfm'
    itysup = 0
    nbval  = 1
    rval(1) = tinoxy
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

!       Numero des zones
    itysup = 3
    nbval  = 1
    rubriq = 'num_zone_fb_slfm'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,izfppp)

!       Entree Fuel (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientfu_zone_bord_slfm'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ientfu)

!       Entree oxydant (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientox_zone_bord_slfm'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ientox)

    car54=' End writing combustion information (SLFM)         '
    write(nfecra,1110)car54

  endif

!      Modele EBU :
!      ==========

  if ( ippmod(icoebu).ge.0 ) then

    rubriq = 'temperature_gaz_frais_ebu'
    itysup = 0
    nbval  = 1
    rval(1) = tgf
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

    rubriq = 'frmel_ebu'
    itysup = 0
    nbval  = 1
    rval(1) = frmel
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

    ! Numero des zones
    itysup = 3
    nbval  = 1
    rubriq = 'num_zone_fb_ebu'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,izfppp)

!       Entree Gaz brule(si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientgb_zone_bord_ebu'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ientgb)

!       Entree gaz frais (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientgf_zone_bord_ebu'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ientgf)

!       FMENT (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    rubriq = 'fment_zone_bord_ebu'
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,fment)

!       TKENT (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    rubriq = 'tkent_zone_bord_ebu'
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,tkent)

    car54=' End writing the combustion information (EBU)      '
    write(nfecra,1110)car54

  endif

!      Modele LWC :
!      ==========

  if ( ippmod(icolwc).ge.0 ) then

    rubriq = 'fmin_lwc'
    itysup = 0
    nbval  = 1
    rval(1) = fmin
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

    rubriq = 'fmax_lwc'
    itysup = 0
    nbval  = 1
    rval(1) = fmax
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

    rubriq = 'hmin_lwc'
    itysup = 0
    nbval  = 1
    rval(1) = hmin
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

    rubriq = 'hmax_lwc'
    itysup = 0
    nbval  = 1
    rval(1) = hmax
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

!       Numero des zones
    itysup = 3
    nbval  = 1
    rubriq = 'num_zone_fb_lwc'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,izfppp)

!       Entree Gaz brule(si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientgb_zone_bord_lwc'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ientgb)

!       Entree gaz frais (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientgf_zone_bord_lwc'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ientgf)

!       FMENT (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    rubriq = 'fment_zone_bord_lwc'
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,fment)

!       TKENT (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    rubriq = 'tkent_zone_bord_lwc'
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,tkent)

    car54=' End writing combustion information (LWC)          '
    write(nfecra,1110)car54

  endif

! ---> Grandeurs complementaires pour la combustion CP

  if (ippmod(iccoal).ge.0) then

!     Charbon PuLVerise : masse vol des charbons

    itysup = 0
    nbval  = 1
    do icha = 1, ncharb
      if (icha.le.nfmtch) then
        write(car2,'(I2.2)') icha
      else
        car2 = cindfc
      endif
      rubriq = 'masse_volumique_charbon'//car2
      rval(1) = rhock(icha)
      call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)
    enddo

    car54=' End writing combustion information (CP)            '
    write(nfecra,1110)car54

  endif

! ---> Grandeurs complementaires pour les versions electriques

  iecr = 0

!     Recalage des CL pot des versions electriques

  if ( ippmod(ieljou).ge.1       ) then
    if (ielcor.eq.1) then

      iecr   = 1
      rubriq = 'coeff_recalage_joule'
      itysup = 0
      nbval  = 1
      rval(1) = coejou

      call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

    endif
  endif

  if ( ippmod(ielarc).ge.1 .or. ippmod(ieljou).ge.1 ) then
    if (ielcor.eq.1) then

      iecr   = 1
      rubriq = 'ddpot_recalage_arc_elec'
      itysup = 0
      nbval  = 1
      rval(1) = pot_diff

      call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

      rubriq = 'elcou_recalage_arc_elec'
      rval(1) = elcou
      call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)
    endif
  endif

  if (iecr.ne.0) then
    car54=' End writing the electric information                '
    write(nfecra,1110)car54
  endif

  call restart_write_fields(rp, RESTART_AUXILIARY)

  ! Fermeture du fichiers suite auxiliaire
  call restart_destroy(rp)

  write(nfecra,1200)

endif
!     Fin de l'ecriture du fichier suite auxiliaire


return

!===============================================================================
! 6. FORMATS
!===============================================================================

 1000 format(3x,'** Writing the main restart file',/,             &
             3x,'   -----------------------------',/)
 1100 format(' Start writing'                                      )
 1110 format('  ',a54                                              )
 1200 format(' End writing'                                        )
 2000 format(/,3x,'** Writing the auxiliary restart file',/,      &
               3x,'   ----------------------------------',/)

 7001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHILE WRITING THE RESTART FILE                 ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@      The maximum number of scalars NSCAMX handled by the   ',/,&
'@        restart file writing format is                      ',/,&
'@        NFMTSC = ',i10                                       ,/,&
'@      The current maximum number of scalars is greater.     ',/,&
'@        NSCAMX = ',i10                                       ,/,&
'@      The scalars with a larger number will not be read.    ',/,&
'@                                                            ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@    Refer to the subroutine ecrava.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7002 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHILE WRITING THE RESTART FILE                 ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@      The maximum number of mass flux NVARMX handled by the ',/,&
'@        restart file writing format is                      ',/,&
'@        NFMTFL = ',i10                                       ,/,&
'@      The current maximum number of mass fluxes is greater. ',/,&
'@        NVARMX = ',i10                                       ,/,&
'@      The fluxes with a larger number will not be read.     ',/,&
'@                                                            ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@    Refer to the subroutine ecrava.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7004 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHILE WRITING THE RESTART FILE                 ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@      The maximum number of coals NCHARM handled by the     ',/,&
'@        restart file writing format is                      ',/,&
'@        NFMTCH = ',i10                                       ,/,&
'@      The current maximum number of coals is greater.       ',/,&
'@        NCHARM = ',i10                                       ,/,&
'@      Some information relative to coals with a greater     ',/,&
'@        number will not be read.                            ',/,&
'@                                                            ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@    Refer to the subroutine ecrava.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7005 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHILE WRITING THE RESTART FILE                 ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@      The number of coal classes NCPCMX handled by the      ',/,&
'@        restart file writing format is                      ',/,&
'@        NFMTCL = ',i10                                       ,/,&
'@      The current number of coal classes is greater.        ',/,&
'@        NCPCMX = ',i10                                       ,/,&
'@      Some information relative to classes with a greater   ',/,&
'@        number will not be read.                            ',/,&
'@                                                            ',/,&
'@    The calculation will be run.                            ',/,&
'@                                                            ',/,&
'@    Refer to the subroutine ecrava.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

end subroutine
