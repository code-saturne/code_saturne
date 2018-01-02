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
use dimens, only: nvar, nscal
use numvar
use cstphy
use entsor
use pointe
use optcal
use albase
use alstru
use alaste
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use cs_fuel_incl
use ppcpfu
use cplsat
use field
use atincl, only: init_at_chem
use atchem, only: ichemistry
use siream, only: iaerosol
use mesh
use cs_c_bindings
use cs_nz_condensation, only:izzftcd, nztag1d, ztpar
use cs_nz_tagmr, only:znmurx, znmur, ztmur

!===============================================================================

implicit none

! Arguments

! Local variables

character        rubriq*64,car2*2,car4*4,car54*54
character        cindfc*2,cindfl*4
character        cstruc(nstrmx)*2, cindst*2
character        ficsui*32
integer          f_id, t_id, ivar, iscal
integer          idecal, iclapc, icha  , icla
integer          ii    , ivers
integer          itysup, nbval
integer          ipcefj, ipcla
integer          nfmtsc, nfmtfl, nfmtch, nfmtcl
integer          nfmtst
integer          ilecec, iecr
integer          ngbstr(2)
integer          ifac, istr
integer          iz, kk
integer          ival(1)
double precision rval(1), tmpstr(27)

type(c_ptr) :: rp

double precision, allocatable, dimension(:,:) :: tmurbf
double precision, allocatable, dimension(:) :: tparbf
double precision, pointer, dimension(:) :: dt_s
double precision, pointer, dimension(:,:) :: disale

type(var_cal_opt) :: vcopt

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

ficsui = 'main'
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

!  ---> Cavitation
rubriq = 'cavitation'
itysup = 0
nbval  = 1
ival(1) = icavit
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

!  ---> VOF
rubriq = 'vof'
itysup = 0
nbval  = 1
ival(1) = ivofmt
call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

rubriq = 'instant_mobile_precedent'
itysup = 0
nbval  = 1
rval(1) = ttcmob
call restart_write_section_real_t(rp,rubriq,itysup,nbval,rval)

if (ichemistry.gt.0.or.iaerosol.gt.0) then
  rubriq = 'atmospheric_chem'
  itysup = 0
  nbval  = 1
  ival(1) = init_at_chem
  call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)
endif

#if defined(_CS_LANG_FR)
car54 =' Fin de l''ecriture des options                       '
#else
car54 =' End writing the options                              '
#endif
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
  ficsui = 'auxiliary'
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

  !  ---> Cavitation
  rubriq = 'cavitation'
  itysup = 0
  nbval  = 1
  ival(1) = icavit
  call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

  !  ---> VOF
  rubriq = 'vof'
  itysup = 0
  nbval  = 1
  ival(1) = ivofmt
  call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

#if defined(_CS_LANG_FR)
  car54 =' Fin de l''ecriture des dimensions et des options     '
#else
  car54 =' End writing the dimensions and options               '
#endif
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

  !     Masse volumique si elle est variable uniquement
  !     La masse volumique est egalement ecrite pour l'algo. VOF
  if (irovar.eq.1.or.ivofmt.ge.0) then
    ! Masse volumique - cellules
    call restart_write_field_vals(rp, icrom, 0)

    ! Masse volumique du pdt precedent - cellules
    ! only for VOF algo. and dilatable models (idilat = 4, 5)
    if (ivofmt.ge.0.or.idilat.ge.4) then
      call restart_write_field_vals(rp, icrom, 1)
    endif

    ! Masse volumique - faces de bord
    call restart_write_field_vals(rp, ibrom, 0)

    ! Scalar source terms for dilatable model (idilat = 4, 5)
    if (idilat.ge.4) then
      do iscal = 1, nscal
        f_id = iustdy(iscal)
        call restart_write_field_vals(rp, f_id, 0)
      enddo
    endif

  endif

  !     On n'ecrit les proprietes physiques que si on les extrapole ou
  !     pour le modele de cavitation
  !       On pourrait les ecrire a tous les coups en prevision d'une
  !       suite avec extrapolation, mais
  !          - c'est rare
  !          - si on demarre un calcul a l'ordre deux a partir d'un calcul
  !            a l'ordre 1, on peut estimer que les premiers pas de temps
  !            sont a jeter de toute facon.
  !       Une exception : on ecrit egalement Cp en effet joule pour
  !         pouvoir calculer la temperature H/Cp en debut de calcul

  if (iviext.gt.0.or.ivofmt.ge.0) then
    !  Viscosite moleculaire - cellules (si variable ou cavitation)
    if (ivivar.eq.1.or.ivofmt.ge.0) then
      call restart_write_field_vals(rp, iviscl, 0)
    endif

    if (iviext.gt.0) then
      ! Viscosite turbulente ou de sous-maille - cellules
      call restart_write_field_vals(rp, ivisct, 0)
    endif
  endif

  if ((icpext.gt.0.and.icp.ge.0).or.              &
       (ippmod(ieljou).ge.1.and.icp.ge.0))  then
    !  Chaleur massique - cellules
    call restart_write_field_vals(rp, icp, 0)
  endif

  call restart_write_linked_fields(rp, "scalar_diffusivity_id", iecr)

#if defined(_CS_LANG_FR)
  car54 =' Fin de l''ecriture des proprietes physiques          '
#else
  car54 =' End writing the physical properties                  '
#endif
  write(nfecra,1110)car54

! ---> Pas de temps

  call field_get_id('dt', f_id)
  if (idtvar.eq.2) then
    call restart_write_field_vals(rp, f_id, 0)
  elseif (idtvar.eq.1) then
    call field_get_val_s(f_id, dt_s)
    rubriq = 'dt_variable_temps'
    itysup = 0
    nbval  = 1
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,dt_s)
  endif

#if defined(_CS_LANG_FR)
  car54 =' Fin de l''ecriture du pas de temps                   '
#else
  car54 =' End writing the time step                            '
#endif
  write(nfecra,1110)car54

  ! Mass fluxes

  call restart_write_linked_fields(rp, "inner_mass_flux_id", iecr)
  call restart_write_linked_fields(rp, "boundary_mass_flux_id", iecr)

  ! Symmetry flag (used for least-squares gradients,
  ! with extrapolation at boundary).

  rubriq = 'isympa_fb_phase01'
  itysup = 3
  nbval  = 1
  call restart_write_section_int_t(rp,rubriq,itysup,nbval,isympa)

  ! Boundary condition coefficients

  call restart_write_bc_coeffs(rp)

! ---> Termes sources
!      Lorsqu'ils sont extrapoles (pour les versions elec, voir plus bas)

  call restart_write_linked_fields(rp, "source_term_prev_id", iecr)

  if (iecr.ne.0) then
#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des termes sources                '
#else
    car54=' End writing the source terms                         '
#endif
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

  if (ineedy.eq.1) then
    iecr   = 1
    itysup = 1
    nbval  = 1
    rubriq = 'dist_fac_par_ce_phase01'
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,dispar)
  endif

  if (iecr.ne.0) then
#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture de la distance a la paroi         '
#else
    car54=' End writing the wall distance                        '
#endif
    write(nfecra,1110)car54
  endif

! ---> Force exterieure

  if (iphydr.eq.1) then
    call field_get_id('volume_forces', f_id)
    call restart_write_field_vals(rp, f_id, 0)
  endif

! ---> Pression hydrostatique predite

  if (iphydr.eq.2) then
    call field_get_id('hydrostatic_pressure_prd', f_id)
    call restart_write_field_vals(rp, f_id, 0)
  endif

!----------------------------------------------------------------
! ---> Wall temperature associated to the condensation model
!      with or wihtout the 1D thermal model tag1D
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
        ifac= ifbpcd(ii)
        iz  = izzftcd(ii)
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
        ifac= ifbpcd(ii)
        iz  = izzftcd(ii)
        tparbf(ifac) = ztpar(iz)
      enddo
      call restart_write_section_real_t(rp,rubriq,itysup,nbval,tparbf)

      ! Free memory
      deallocate(tparbf)
    endif

  endif

! ---> Methode ALE

  if (iale.eq.1) then

    itysup = 4
    nbval  = 3

    call field_get_val_v(fdiale, disale)

    rubriq = 'vertex_displacement'
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,disale)

!     Viscosite de maillage (elle est souvent definie geometriquement sur le
!       maillage initial ... il est donc plus facile de la relire ensuite)

    rubriq = 'type_visc_mail'
    itysup = 0
    nbval  = 1
    ival(1) = iortvm
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ival)

    call restart_write_field_vals(rp, ivisma, 0)

#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des donnees ALE    '
#else
    car54=' End writing the ALE data              '
#endif
    write(nfecra,1110)car54

    ngbstr(1) = nbstru
    ngbstr(2) = nbaste

    rubriq = 'nombre_structures'
    itysup = 0
    nbval  = 2
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ngbstr)

    if (nbstru.gt.0) then

      nfmtst = 99
      cindst = 'XX'
!     Codage en chaine de caracteres du numero de la structure
      do istr = 1, min(nbstru ,nfmtst)
        write(cstruc(istr),'(I2.2)') istr
      enddo
      do istr = min(nbstru ,nfmtst)+1,nbstru
        cstruc(istr) = cindst
      enddo

      do istr = 1, nbstru
        rubriq = 'donnees_structure_'//cstruc(istr)
        itysup = 0
        nbval  = 27

        do ii = 1, 3
          tmpstr(   ii) = xstr  (ii,istr)
          tmpstr(3 +ii) = xpstr (ii,istr)
          tmpstr(6 +ii) = xppstr(ii,istr)
          tmpstr(9 +ii) = xsta  (ii,istr)
          tmpstr(12+ii) = xpsta (ii,istr)
          tmpstr(15+ii) = xppsta(ii,istr)
          tmpstr(18+ii) = xstp  (ii,istr)
          tmpstr(21+ii) = forstr(ii,istr)
          tmpstr(24+ii) = forsta(ii,istr)
        enddo

        call restart_write_section_real_t(rp,rubriq,itysup,nbval,tmpstr)
      enddo

#if defined(_CS_LANG_FR)
      car54=' Fin de l''ecriture des donnees des structures (ALE)'
#else
      car54=' End writing the structures data (ALE)              '
#endif
      write(nfecra,1110)car54

    endif
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

#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des informations combustion COD3P'
#else
    car54=' End writing combustion information (COD3P)         '
#endif
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

#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des informations combustion EBU '
#else
    car54=' End writing the combustion information (EBU)      '
#endif
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

#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des informations combustion LWC '
#else
    car54=' End writing combustion information (LWC)          '
#endif
    write(nfecra,1110)car54

  endif

! ---> Grandeurs complementaires pour la combustion CP

  if (ippmod(icpl3c).ge.0 .or. ippmod(iccoal).ge.0) then

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


!     Charbon PuLVerise : type de zones de bord, ientat, inmoxy, ientcp, timpat
!       x20, pour le calcul de rho au bord en entree

!       Numero des zones
    itysup = 3
    nbval  = 1
    rubriq = 'num_zone_fb_charbon_pulverise'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,izfppp)

!       Type entree air ou cp (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientat_zone_bord_charbon_pulverise'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ientat)

!       ientat, inmoxy et x20 ne servent pas pour le CP couple Lagrangien (cplphy)
    if (ippmod(iccoal).ge.0) then

      itysup = 0
      nbval  = nozppm
      rubriq = 'ientcp_zone_bord_charbon_pulverise'
      call restart_write_section_int_t(rp,rubriq,itysup,nbval,ientcp)

      itysup = 0
      nbval  = nozppm
      rubriq = 'inmoxy_zone_bord_charbon_pulverise'
      call restart_write_section_int_t(rp,rubriq,itysup,nbval,inmoxy)

      itysup = 0
      nbval  = nozppm

      idecal = 0
      do icha = 1, ncharb
        do iclapc = 1, nclpch(icha)
          icla = iclapc + idecal
          if (icha.le.nfmtch.and.iclapc.le.nfmtcl) then
            write(car2,'(I2.2)') icha
            write(car4,'(I4.4)') iclapc
          else
            car2 = cindfc
            car4 = cindfl
          endif
          rubriq = 'x20_zone_bord_charbon'//car2//'_classe'//car4
          call restart_write_section_real_t(rp,rubriq,itysup,nbval,x20(:,icla))
        enddo
      enddo

    endif

!       Temperature
    itysup = 0
    nbval  = nozppm
    rubriq = 'timpat_zone_bord_charbon_pulverise'
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,timpat)

#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des informations combustion CP    '
#else
    car54=' End writing combustion information (CP)            '
#endif
    write(nfecra,1110)car54

  endif

! ---> Grandeurs complementaires pour la FUEL

  if ( ippmod(icfuel).ge.0 ) then

!     Fioul : type de zones de bord, ientat, ientfl, timpat
!       qimpat et qimpfl  pour le calcul de rho au bord en entree

!       Numero des zones
    itysup = 3
    nbval  = 1
    rubriq = 'num_zone_fb_fuel'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,izfppp)

!       Type entree air ou fuel (si ce n'est pas NOZPPM, erreur)
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientat_zone_bord_fuel'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ientat)

    itysup = 0
    nbval  = nozppm
    rubriq = 'ientfl_zone_bord_fuel'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,ientfl)

    itysup = 0
    nbval  = nozppm
    rubriq = 'inmoxy_zone_bord_fuel'
    call restart_write_section_int_t(rp,rubriq,itysup,nbval,inmoxy)

!       Timpat
    itysup = 0
    nbval  = nozppm
    rubriq = 'timpat_zone_bord_fuel'
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,timpat)

!       Qimpat
    itysup = 0
    nbval  = nozppm
    rubriq = 'qimpat_zone_bord_fuel'
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,qimpat)

!       Qimpfl
    itysup = 0
    nbval  = nozppm
    rubriq = 'qimpfl_zone_bord_fuel'
    call restart_write_section_real_t(rp,rubriq,itysup,nbval,qimpfl)

#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des informations combustion FUEL  '
#else
    car54=' End writing combustion information (FUEL)           '
#endif
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

!     Termes sources des versions electriques

  if ( ippmod(ieljou).ge.1 .or.                                   &
       ippmod(ielarc).ge.1       ) then

    call field_get_id('joule_power', ipcefj)
    iecr   = 1

    call restart_write_field_vals(rp, ipcefj, 0)

  endif

  if (ippmod(ielarc).ge.1) then

    iecr   = 1
    call field_get_id('laplace_force', ipcla)

    call restart_write_field_vals(rp, ipcla, 0)

  endif

  if (iecr.ne.0) then
#if defined(_CS_LANG_FR)
    car54=' Fin de l''ecriture des informations electriques      '
#else
    car54=' End writing the electric information                '
#endif
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

#if defined(_CS_LANG_FR)

 1000 format(3x,'** Ecriture du fichier suite principal',/,       &
             3x,'   ----------------------------------- ',/)
 1100 format(' Debut de l''ecriture')
 1110 format('  ',A54)
 1200 format(' Fin de l''ecriture')
 2000 format(/,3x,'** Ecriture du fichier suite auxiliaire',/,    &
               3x,'   ------------------------------------ ',/)

 7001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ECRITURE DU FICHIER SUITE        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      Le nombre de scalaires maximal NSCAMX supporte par le ',/,&
'@        format d''ecriture du fichier suite est             ',/,&
'@        NFMTSC = ',i10                                       ,/,&
'@      On a ici un nombre de scalaires maximal superieur     ',/,&
'@        NSCAMX = ',i10                                       ,/,&
'@      On ne pourra pas relire les scalaires dont le numero  ',/,&
'@        est superieur                                       ',/,&
'@                                                            ',/,&
'@    Le calcul sera execute.                                 ',/,&
'@                                                            ',/,&
'@    Voir le sous-programme ecrava.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7002 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ECRITURE DU FICHIER SUITE        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      Le nombre de flux de masse max NVARMX supporte par le ',/,&
'@        format d''ecriture du fichier suite est             ',/,&
'@        NFMTFL = ',i10                                       ,/,&
'@      On a ici un nombre de flux      maximal superieur     ',/,&
'@        NVARMX = ',i10                                       ,/,&
'@      On ne pourra pas relire les flux      dont le numero  ',/,&
'@        est superieur                                       ',/,&
'@                                                            ',/,&
'@    Le calcul sera execute.                                 ',/,&
'@                                                            ',/,&
'@    Voir le sous-programme ecrava.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7004 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ECRITURE DU FICHIER SUITE        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      Le nombre de charbons      max NCHARM supporte par le ',/,&
'@        format d''ecriture du fichier suite est             ',/,&
'@        NFMTCH = ',i10                                       ,/,&
'@      On a ici un nombre de charbons  maximal superieur     ',/,&
'@        NCHARM = ',i10                                       ,/,&
'@      On ne pourra pas relire certaines informations        ',/,&
'@        relatives aux charbons dont le numero               ',/,&
'@        est superieur                                       ',/,&
'@                                                            ',/,&
'@    Le calcul sera execute.                                 ',/,&
'@                                                            ',/,&
'@    Voir le sous-programme ecrava.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7005 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ECRITURE DU FICHIER SUITE        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      Le nombre de classes par charbon max NCPCMX supporte  ',/,&
'@        par le format d''ecriture du fichier suite est      ',/,&
'@        NFMTCL = ',i10                                       ,/,&
'@      On a ici un nombre de classes par charbon superieur   ',/,&
'@        NCPCMX = ',i10                                       ,/,&
'@      On ne pourra pas relire certaines informations        ',/,&
'@        relatives aux classes  dont le numero               ',/,&
'@        est superieur                                       ',/,&
'@                                                            ',/,&
'@    Le calcul sera execute.                                 ',/,&
'@                                                            ',/,&
'@    Voir le sous-programme ecrava.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

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

#endif


end subroutine
