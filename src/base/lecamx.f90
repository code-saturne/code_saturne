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

!===============================================================================
! Function :
! --------

!> \file lecamx.f90
!>
!> \brief Reading of auxiliary restart file.
!>
!> Here ileaux = 1.
!> We stop if:
!>  - the file cannot be opened
!>  - the file is not an auxiliary restart file
!>  - ncel is not correct
!>  - jdtvar is not readable
!>  - a required temporal moment is not readable
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     oflmap        pointer to old field map
!_______________________________________________________________________________


subroutine lecamx &
 ( oflmap )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

use paramx
use cstphy
use cstnum
use entsor
use optcal
use pointe
use numvar
use albase
use alaste
use parall
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use ppcpfu
use mesh
use field
use cavitation
use vof
use cs_c_bindings
use cs_nz_condensation, only: nfbpcd, izzftcd, nztag1d, ztpar, ifbpcd
use cs_nz_tagmr, only: znmurx, znmur, ztmur

!===============================================================================

implicit none

! Arguments

type(c_ptr)      oflmap

! Local variables

character        rubriq*64,car2*2
character        car54*54
character        cindfp*2,cindfs*4,cindff*4,cindfm*4
character        cindfc*2,cindfl*4
character        ficsui*32
logical          lprev
integer          iel   , ifac, ii, nlfld
integer          iz, kk
integer          icha
integer          jdtvar
integer          ierror, itysup, nbval
integer          nberro, inierr, ivers(1)
integer          ilu   , ierrch
integer          nfmtsc, nfmtfl, nfmtch, nfmtcl
integer          jale, jvolfl
integer          f_id, iflmas, iflmab, iflvoi, iflvob
integer          key_t_ext_id, icpext
integer          iviext
integer          ival(1)
integer          iread_status_loc, iread_status_glob
double precision rval(1)

logical(kind=c_bool) :: ncelok, nfaiok, nfabok, nsomok

type(c_ptr) :: rp

double precision, dimension(:), pointer :: sval
double precision, dimension(:), pointer :: voidfl

double precision, allocatable, dimension(:,:) :: tmurbf
double precision, allocatable, dimension(:) :: tparbf

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_ale_restart_read(r)  &
    bind(C, name='cs_ale_restart_read')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), value :: r
  end subroutine cs_ale_restart_read

  subroutine cs_mobile_structures_restart_read(r)  &
    bind(C, name='cs_mobile_structures_restart_read')
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr), value :: r
  end subroutine cs_mobile_structures_restart_read

end interface

!===============================================================================
! 0. Initialisation
!===============================================================================

!  ---> Banniere
write(nfecra,1000)

!  --->  On code en chaine le numero des phases et scalaires

!     Nombre max pour les formats choisis
nfmtsc = 9999
nfmtfl = 9999
nfmtch = 99
nfmtcl = 9999

!     Indefini a 2 et 4 caracteres
cindfp='YY'
cindfs='YYYY'
cindff='YYYY'
cindfm='YYYY'
cindfc='YY'
cindfl='YYYY'

! Time extrapolation?
call field_get_key_id("time_extrapolated", key_t_ext_id)

!===============================================================================
! 1. OUVERTURE DU FICHIER SUITE AUXILIAIRE
!===============================================================================

ficsui = 'auxiliary.csc'
call restart_create(ficsui, '', 0, rp)

! ---> Debut de la lecture
write(nfecra,1100)

!===============================================================================
! 2. ENTETES DU FICHIER SUITE OU STOP
!===============================================================================

!  Rubrique "fichier suite aux"
!        Pourrait porter le numero de version si besoin.
!        On ne se sert pas de IVERS pour le moment.

itysup = 0
nbval  = 1

call restart_read_int_t_compat(rp,                                           &
                               'code_saturne:checkpoint:auxiliary:version',  &
                               'version_fichier_suite_auxiliaire',           &
                               itysup, nbval, ivers, ierror)


if (ierror.ne.0) then
  write(nfecra,9100)ficsui
  call csexit (1)
endif

!     Supports

call restart_check_base_location(rp,ncelok,nfaiok,nfabok,nsomok)

if (ncelok.eqv..false.) then
  write(nfecra,9101)
  call csexit (1)
endif

if (nfaiok.eqv..false.) write(nfecra,8200)'internes','internes'

if (nfabok.eqv..false.) write(nfecra,8200)'de bord ','de bord '

!     Methode ALE

nberro = 0

rubriq = 'methode_ALE'
itysup = 0
nbval  = 1
call restart_read_section_int_t(rp,rubriq,itysup,nbval,ival,ierror)
jale = ival(1)
nberro=nberro+ierror

! --->  Message si erreur (pas de stop pour compatibilite avec les fichiers anterieurs)
!       -> on n'affiche le message que si IALE>=1 (sinon RAS)
if (nberro.ne.0) then
  if (iale.ge.1) write(nfecra,9210)
  jale = 0
endif

! ---> Pas d'iteration d'initialisation si suite de calcul ALE
if (italin.eq.-999) then
  if (iale.ge.1 .and. jale.ge.1) then
    italin = 0
  else if (iale.ge.1) then
    italin = 1
  else
    italin = 0
  endif
endif

!     VOF

nberro = 0

rubriq = 'vof'
itysup = 0
nbval  = 1
call restart_read_section_int_t(rp,rubriq,itysup,nbval,ival,ierror)
jvolfl = ival(1)
nberro=nberro+ierror

! warning if error found (backward compatibility)
! and only if VoF model is enabled
if (nberro.ne.0) then
  if (ivofmt.gt.0) write(nfecra,9221)
  jvolfl = 0
endif

car54 =' Finished reading options.                            '
write(nfecra,1110)car54

!===============================================================================
! 3. PROPRIETES PHYSIQUES
!===============================================================================

nberro = 0

!     On lit les infos des phases communes

! ---> Point de reference de pression
!     On lit les coordonnees si l'utilisateur n'a rien specifie, i.e.
!       si IXYZP0=-1, et on met IXYZP0 a 1 pour eviter de le changer ensuite.
if (ixyzp0.eq.-1) then
  rubriq = 'ref_presstot01'
  itysup = 0
  nbval  = 3
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,      &
                                   xyzp0,ierror)
  nberro = nberro+ierror
  if (ierror.eq.0) then
    write(nfecra,7000) (xyzp0(ii),ii=1,3)
    ixyzp0 = 1
  endif
endif

! Here the physical variables below are required for the low-Mach algorithm
if (idilat.eq.3.or.ipthrm.eq.1) then

  !the reference density updated with the low-Mach algorithm
  rubriq = 'ro001'
  itysup = 0
  nbval  = 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  ro0 = rval(1)
  nberro=nberro+ierror

  ! the thermodynamic pressure for the previous time step
  rubriq = 'pther01'
  itysup = 0
  nbval  = 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  pther = rval(1)
  nberro=nberro+ierror
endif

! ---> Masse volumique
!     On la lit, qu'elle soit extrapolee ou pas,
!       pour permettre les sous-relaxations
!     Pour les suites a rho constant, cependant, on ne la lit pas,
!       afin de ne pas ecraser la valeur RO0 fixee par l'utilisateur
!       et eventuellement modifiee.
!     La masse volumique est egalement lue pour le modele de cavitation
!     et les modeles dilatables (idilat = 4, 5)

inierr = 0

if (    irovar.eq.1                   &
    .or.(ivofmt.gt.0.and.jvolfl.ge.0)) then

  ! Si on a reussi a initialiser la masse volumique aux cellules ET
  ! aux faces de bord, on l'indique (pour schtmp)
  iread_status_loc  = 1
  iread_status_glob = 1

  call restart_get_field_read_status(icrom, iread_status_loc)
  if (iread_status_loc.eq.0) then
    iread_status_glob = 0
  endif
  if (nfabok.eqv..true.) then
    call restart_get_field_read_status(ibrom, iread_status_loc)
    if (iread_status_loc.eq.0) then
      iread_status_glob = 0
    endif
  endif

  if (nfabok.eqv..true..and.iread_status_glob.eq.1) then
    initro = 1
  endif

else
  !     Si la masse volumique est constante, elle est initialisee
  !       correctement
  initro = 1
endif

! ---> Viscosite moleculaire et "turbulente" ou de "sous-maille"
!     Si elle est extrapolee en temps, on la lit
!     La viscosite moleculaire est egalement lue pour le modele de cavitation
!     Si on reussit, on l'indique


call field_get_key_int(ivisct, key_t_ext_id, iviext)
if (iviext.gt.0.or.(ivofmt.gt.0.and.jvolfl.ge.0)) then
  !     Si on a initialise les viscosites, on l'indique (pour schtmp)
  iread_status_glob = 1
  call restart_get_field_read_status(iviscl, iread_status_loc)
  if (iread_status_loc.eq.0) then
    iread_status_glob = 0
  endif
  call restart_get_field_read_status(ivisct, iread_status_loc)
  if (iread_status_loc.eq.0) then
    iread_status_glob = 0
  endif
  if (iread_status_loc.eq.1) then
    initvi = 1
  endif
endif


! ---> Chaleur massique
!     On cherche a la lire si elle est variable
!       et qu'on l'extrapole ou qu'on est en Joule
!     Si on reussit, on l'indique
!       (Ca sert quand elle est extrapolee en temps
!        et quand l'utilisateur peut s'en servir pour passer
!        de H a T, comme en effet Joule par exemple).

if (icp.ge.0) then

  call field_get_key_int(icp, key_t_ext_id, icpext)

  if (icpext.gt.0.or.ippmod(ieljou).ge.1) then

    inierr = 0
    ! Si on a initialise Cp, on l'indique (pour schtmp)
    call restart_get_field_read_status(icp, iread_status_loc)
    if (iread_status_loc.eq.1) then
      initcp = 1
    endif

  endif
endif

!     Si on a des scalaires, on lit a diffusivite
!       si le scalaire a un correspondant et
!       si on doit extrapoler la diffusivite
!       (et qu'elle est variable, et que le scalaire n'est pas une variance)

nlfld = 0
call restart_read_linked_fields(rp, oflmap, "diffusivity_id", nlfld)

!     Si erreur, on previent mais pas stop :
!       auparavant on n'avait pas stocke les prop phy
!         si on n'etait pas a l'ordre 2
!       c'est discutable pour rho

if (nberro.ne.0) then
  car54 = 'Reading physical properties                         '
  write(nfecra,8300)car54
endif

car54 = ' Finished reading physical properties.                '
write(nfecra,1110)car54

!===============================================================================
! 4. PAS DE TEMPS
!===============================================================================

!  ---> Indicateur de pas de temps variable
rubriq = 'indic_dt_variable'
itysup = 0
nbval  = 1
call restart_read_section_int_t(rp,rubriq,itysup,nbval,ival,ierror)
jdtvar = ival(1)

!  ---> On s'arrete si erreur
!     Si on ne peut pas relire un entier, c'est que le fichier n'est pas bon
if (ierror.ne.0) then
  car54 ='Error reading the time stepping mode                  '
  write(nfecra,9200)car54
  call csexit(1)
endif

!  ---> Pas de temps
!     Rq : si marche en temps differente, on conserve la valeur par defaut
!          DTREF imposee dans INIVA0

nberro = 0

call field_get_id('dt', f_id)

if (idtvar.ne.jdtvar) then
  write(nfecra,8400)idtvar,jdtvar,dtref

elseif (idtvar.eq.1) then
  rubriq = 'dt_variable_temps'
  itysup = 0
  nbval  = 1
  call field_get_val_s(f_id, sval)
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,sval,ierror)
  nberro=nberro+ierror
  do iel = 1, ncel
    sval(iel) = sval(1)
  enddo

elseif (idtvar.eq.2) then
  call restart_read_field_vals(rp, f_id, 0, ierror)
  nberro=nberro+ierror
endif

!     Si erreur, on previent mais pas stop :
!       le pas de temps n'est pas une donnee majeure
!       c'est discutable

if (nberro.ne.0) then
  car54 = 'Reading time step                                     '
  write(nfecra,8300)car54
endif

car54 = ' Finished reading time step.                          '
write(nfecra,1110)car54

!===============================================================================
! 5. FLUX DE MASSE
!===============================================================================
!     Pour retrouver la correspondance entre les variables
!     et les flux de masse, on utilise le nom de chaque variable
!     (nomflu(I)= nom de la ieme variable)
!     Ensuite, pour chaque variable, si on a deja le flux, on ne fait
!       rien, sinon on lit quel est le numero local du
!       flux qui lui est associe (en pratique 1 ou 2) et le flux lui meme.

!     Les flux de masse ne sont a lire que si les supports des faces
!     de bord ou des faces internes coincident

!     On lit d'abord le flux de masse (a l'instant n) et
!       ensuite a l'instant precedent si on est en schema en temps
!       particulier (ISTMPF NE 1)

if (nfaiok.eqv..true. .or. nfabok.eqv..true.) then

  ! restart_read_linked_fields reads all time steps present and required

  nlfld = 0
  call restart_read_linked_fields(rp, oflmap, "inner_mass_flux_id", nlfld)
  call restart_read_linked_fields(rp, oflmap, "boundary_mass_flux_id", nlfld)

  nberro=0

  ! Initialization of the void fraction convective flux, if required
  if (ivofmt.gt.0.and.jvolfl.lt.0) then

    ! Interior faces

    call field_get_key_int(ivarfl(iu), kimasf, iflmas)
    call field_get_key_int(ivarfl(ivolf2), kimasf, iflvoi)

    call field_get_val_s(iflmas, sval)
    call field_get_val_s(iflvoi, voidfl)
    do ifac = 1, nfac
      voidfl(ifac) = sval(ifac)/rho1
    enddo

    call field_have_previous(iflmas, lprev)
    if (lprev .neqv. .false.) then
      call field_get_val_prev_s(iflmas, sval)
      call field_get_val_prev_s(iflvoi, voidfl)
      do ifac = 1, nfac
        voidfl(ifac) = sval(ifac)/rho1
      enddo
    endif

    call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
    call field_get_val_s(iflmab, sval)
    call field_get_key_int(ivarfl(ivolf2), kbmasf, iflvob)
    call field_get_val_s(iflvob, voidfl)
    do ifac = 1, nfabor
      voidfl(ifac) = sval(ifac)/rho1
    enddo

    ! Boundary faces

    call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
    call field_get_key_int(ivarfl(ivolf2), kbmasf, iflvob)

    call field_get_val_s(iflmab, sval)
    call field_get_val_s(iflvob, voidfl)
    do ifac = 1, nfabor
      voidfl(ifac) = sval(ifac)/rho1
    enddo

    call field_have_previous(iflmas, lprev)
    if (lprev .neqv. .false.) then
      call field_get_val_prev_s(iflmab, sval)
      call field_get_val_prev_s(iflvob, voidfl)
      do ifac = 1, nfabor
        voidfl(ifac) = sval(ifac)/rho1
      enddo
    endif

  endif

  ! In case of error, warn but do not stop
  if (nberro.ne.0) then
    car54 = 'Reading mass flux                                     '
    write(nfecra,8300) car54
  endif

  car54 = ' Finished reading mass flux                           '
  write(nfecra,1110) car54

endif

! fin de "s'il faut lire les flux de masse (ie. supports coincidents)"

!===============================================================================
! 6. CONDITIONS AUX LIMITES
!===============================================================================

! A ne relire que si les supports sont identiques (faces de bord)

ilu = 0

if (nfabok.eqv..true.) then

  ilu = 1

  nberro=0

  ! Variable BC coefficients

  call restart_read_bc_coeffs(rp)

endif
!     fin du test "si les supports des faces de bord sont identiques"

if (ilu.eq.1) then

!     Si erreur, on previent mais pas stop :
!       (on n'a pas forcement les coefs 2
!        si on prend les fichiers d'une version anterieure)
  if (nberro.ne.0) then
    car54 = 'Reading boundary conditions                           '
    write(nfecra,8300)car54
  endif

  car54 = ' Finished reading boundary conditions.                '
  write(nfecra,1110)car54

endif

!===============================================================================
! 7. TERMES SOURCES EXTRAPOLES et TERMES SOURCES SCHEMA EN TEMPS
!===============================================================================

! Variables at previous time steps

nberro = 0
ilu = 0

! ---> Termes sources

! Do not use iscold for scalars here, as we use field names and not
! scalar numbers here, and names are assumed stable

call restart_read_linked_fields(rp, oflmap, "source_term_prev_id", ilu)

if (ilu.ne.0) then
  car54 =' Finished reading source terms.                       '
  write(nfecra,1110)car54
endif

!===============================================================================
! 8. MOYENNES
!===============================================================================

call time_moment_restart_read(rp)

!===============================================================================
! 9. DISTANCE A LA PAROI
!===============================================================================

ilu = 0
nberro = 0

! On la lit si on en a besoin uniquement.

! Si l'utilisateur a force le recalcul, on ne la lit pas
!   il faudra la mettre a jour (sauf si zero pas de temps).

! Sinon, on cherche a la lire.
!   On pourrait la relire aussi quand le nombre de faces a
!   change, mais il vaut mieux la recalculer au cas ou des faces de
!   paroi auraient disparu
!   Si on arrive a la lire, on note qu'elle est a jour (sauf si ALE).

if (ineedy.eq.1) then
  if (icdpar.gt.0) then
    if (nfabok.eqv..true.) then
      if (ierror.eq.0 .and. iale.eq.0 ) then
        imajdy = 1
      endif
      ilu   = ilu + 1
    endif
  endif
endif

if (nberro.ne.0) then
  car54 = 'Reading wall distance                                 '
  write(nfecra,8300)car54
endif

if (ilu.ne.0) then
  car54 = ' Finished reading wall distance.                      '
  write(nfecra,1110)car54
endif

!===============================================================================
! 11. Predicted hydrostatic pressure
!===============================================================================

!===============================================================================
! 12.  Wall temperature associated to the condensation model
!      with or wihtout the 1D thermal model tag1D
!===============================================================================

if (icondb.eq.0) then

  if (nztag1d.eq.1) then

    !tmur array associated at each boundary face in the wall tickness
    nberro = 0
    itysup = 3
    nbval  = znmurx
    rubriq = 'tmur_bf_prev'

    allocate(tmurbf(nfabor,znmurx))

    call restart_read_section_real_t(rp,rubriq,itysup,nbval,tmurbf,ierror)
    nberro=nberro+ierror

    do ii = 1, nfbpcd
      ifac = ifbpcd(ii) + 1 ! C numbering
      iz  = izzftcd(ii) + 1 ! C numbering
      do kk = 1, znmur(iz)
        ztmur(ii,kk) =  tmurbf(ifac,kk)
      enddo
    enddo
    ! Free memory
    deallocate(tmurbf)
  else
    !tmur array associated at each boundary face without 1D thermal model
    nberro = 0
    itysup = 3
    nbval  = 1
    rubriq = 'tpar_bf_prev'

    allocate(tparbf(nfabor))

    call restart_read_section_real_t(rp,rubriq,itysup,nbval,tparbf,ierror)
    nberro=nberro+ierror

    do ii = 1, nfbpcd
      ifac= ifbpcd(ii) + 1 ! C numbering
      iz  = izzftcd(ii)+ 1 ! C numbering)
      do kk = 1, znmur(iz)
        ztpar(iz) =  tparbf(ifac)
      enddo
    enddo

    ! Free memory
    deallocate(tparbf)

  endif

endif

!===============================================================================
! 13.  DEPLACEMENT AUX NOEUDS EN ALE
!===============================================================================

if (iale.ge.1 .and. jale.ge.1) then

  call cs_ale_restart_read(rp)
  call cs_mobile_structures_restart_read(rp)

  car54 =' Finished reading ALE information.                    '
  write(nfecra,1110) car54

endif

!===============================================================================
! 14. LECTURE DES INFORMATIONS COMPLEMENTAIRES COMBUSTION GAZ, CP
!===============================================================================

nberro = 0
ilu = 0

!     Modele COD3P :
!     ============

if ( ippmod(icod3p).ge.0 ) then

  rubriq = 'hinfue_cod3p'
  itysup = 0
  nbval  = 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  hinfue = rval(1)
  nberro=nberro+ierror
  ilu = ilu + 1
  if (ierror.ne.0) then
    write(nfecra,9400)
  endif

  rubriq = 'hinoxy_cod3p'
  itysup = 0
  nbval  = 1
  ilu    = ilu + 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  hinoxy = rval(1)
  nberro=nberro+ierror
  ilu = ilu + 1
  if (ierror.ne.0) then
    write(nfecra,9400)
  endif

  rubriq = 'tinfue_cod3p'
  itysup = 0
  nbval  = 1
  ilu    = ilu + 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  tinfue = rval(1)
  nberro=nberro+ierror
  ilu = ilu + 1
  if (ierror.ne.0) then
    write(nfecra,9400)
  endif

  rubriq = 'tinoxy_cod3p'
  itysup = 0
  nbval  = 1
  ilu    = ilu + 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  tinoxy = rval(1)
  nberro=nberro+ierror
  ilu = ilu + 1
  if (ierror.ne.0) then
    write(nfecra,9400)
  endif

!       Il faut le meme nbr de faces de bord, sinon on ne lit pas
  if (nfabok.eqv..true.) then

    ilu = ilu + 1

    ierrch = 0

!       Numero des zones
    itysup = 3
    nbval  = 1
    rubriq = 'num_zone_fb_cod3p'
    call restart_read_section_int_t(rp,rubriq,itysup,nbval,izfppp,ierror)
    nberro=nberro+ierror

!       Type entree Fuel
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientfu_zone_bord_cod3p'
    call restart_read_section_int_t(rp,rubriq,itysup,nbval,ientfu,ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       Type entree Oxydant
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientox_zone_bord_cod3p'
    call restart_read_section_int_t(rp,rubriq,itysup,nbval,ientox,ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!     Par securite, si on ne parvient pas a lire la
!       IENTCPFU ou IENTOX, on remet a zero le numero des zones IZFPPP
!       car il a peut etre ete lu.
!       Ceci permettra d'eviter de se servir des valeurs par defaut

    if (ierrch.ne.0) then
      do ifac = 1, nfabor
        izfppp(ifac) = 0
      enddo
    endif

  endif

endif

!     Modele SLFM :
!     ===========

if ( ippmod(islfm).ge.0 ) then

  rubriq = 'hinfue_slfm'
  itysup = 0
  nbval  = 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  hinfue = rval(1)
  nberro=nberro+ierror
  ilu = ilu + 1
  if (ierror.ne.0) then
    write(nfecra,9400)
  endif

  rubriq = 'hinoxy_slfm'
  itysup = 0
  nbval  = 1
  ilu    = ilu + 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  hinoxy = rval(1)
  nberro=nberro+ierror
  ilu = ilu + 1
  if (ierror.ne.0) then
    write(nfecra,9400)
  endif

  rubriq = 'tinfue_slfm'
  itysup = 0
  nbval  = 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  tinfue = rval(1)
  nberro=nberro+ierror
  ilu = ilu + 1
  if (ierror.ne.0) then
    write(nfecra,9400)
  endif

  rubriq = 'tinoxy_slfm'
  itysup = 0
  nbval  = 1
  ilu    = ilu + 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  tinoxy = rval(1)
  nberro=nberro+ierror
  ilu = ilu + 1
  if (ierror.ne.0) then
    write(nfecra,9400)
  endif

!       Il faut le meme nbr de faces de bord, sinon on ne lit pas
  if (nfabok.eqv..true.) then

    ilu = ilu + 1

    ierrch = 0

!       Numero des zones
    itysup = 3
    nbval  = 1
    rubriq = 'num_zone_fb_slfm'
    call restart_read_section_int_t(rp,rubriq,itysup,nbval,izfppp,ierror)
    nberro=nberro+ierror

!       Type entree Fuel
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientfu_zone_bord_slfm'
    call restart_read_section_int_t(rp,rubriq,itysup,nbval,ientfu,ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       Type entree Oxydant
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientox_zone_bord_slfm'
    call restart_read_section_int_t(rp,rubriq,itysup,nbval,ientox,ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!     Par securite, si on ne parvient pas a lire la
!       IENTCPFU ou IENTOX, on remet a zero le numero des zones IZFPPP
!       car il a peut etre ete lu.
!       Ceci permettra d'eviter de se servir des valeurs par defaut

    if (ierrch.ne.0) then
      do ifac = 1, nfabor
        izfppp(ifac) = 0
      enddo
    endif

  endif

endif

!     Modele EBU :
!     ==========

if ( ippmod(icoebu).ge.0 ) then

  rubriq = 'temperature_gaz_frais_ebu'
  itysup = 0
  nbval  = 1
  ilu    = 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  tgf = rval(1)
  nberro=nberro+ierror
  if (ierror.ne.0) then
    write(nfecra,9500)
  endif

  rubriq = 'frmel_ebu'
  itysup = 0
  nbval  = 1
  ilu    = 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  frmel = rval(1)
  nberro=nberro+ierror
  if (ierror.ne.0) then
    write(nfecra,9500)
  endif

!       Il faut le meme nbr de faces de bord, sinon on ne lit pas
  if (nfabok.eqv..true.) then

    ilu = ilu + 1

    ierrch = 0

!       Numero des zones
    itysup = 3
    nbval  = 1
    rubriq = 'num_zone_fb_ebu'
    call restart_read_section_int_t(rp,rubriq,itysup,nbval,izfppp,ierror)
    nberro=nberro+ierror

!       Type entree Gaz brulee
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientgb_zone_bord_ebu'
    call restart_read_section_int_t(rp,rubriq,itysup,nbval,ientgb,ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       Type entree gaz frais
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientgf_zone_bord_ebu'
    call restart_read_section_int_t(rp,rubriq,itysup,nbval,ientgf,ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       FMENT
    itysup = 0
    nbval  = nozppm
    rubriq = 'fment_zone_bord_ebu'
    call restart_read_section_real_t(rp,rubriq,itysup,nbval,fment,ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       TKENT
    itysup = 0
    nbval  = nozppm
    rubriq = 'tkent_zone_bord_ebu'
    call restart_read_section_real_t(rp,rubriq,itysup,nbval,tkent,ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!     Par securite, si on ne parvient pas a lire la
!       IENTCPFU ou IENTOX, on remet a zero le numero des zones IZFPPP
!       car il a peut etre ete lu.
!       Ceci permettra d'eviter de se servir des valeurs par defaut

    if (ierrch.ne.0) then
      do ifac = 1, nfabor
        izfppp(ifac) = 0
      enddo
    endif

  endif

endif

!     Modele LWC :
!     ==========

if ( ippmod(icolwc).ge.0 ) then

  rubriq = 'fmin_lwc'
  itysup = 0
  nbval  = 1
  ilu    = 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  fmin = rval(1)
  nberro=nberro+ierror
  if (ierror.ne.0) then
    write(nfecra,9600)
  endif

  rubriq = 'fmax_lwc'
  itysup = 0
  nbval  = 1
  ilu    = 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  fmax = rval(1)
  nberro=nberro+ierror
  if (ierror.ne.0) then
    write(nfecra,9600)
  endif

  rubriq = 'hmin_lwc'
  itysup = 0
  nbval  = 1
  ilu    = 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  hmin = rval(1)
  nberro=nberro+ierror
  if (ierror.ne.0) then
    write(nfecra,9600)
  endif

  rubriq = 'hmax_lwc'
  itysup = 0
  nbval  = 1
  ilu    = 1
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  hmax = rval(1)
  nberro=nberro+ierror
  if (ierror.ne.0) then
    write(nfecra,9600)
  endif

!       Il faut le meme nbr de faces de bord, sinon on ne lit pas
  if (nfabok.eqv..true.) then

    ilu = ilu + 1

    ierrch = 0

!       Numero des zones
    itysup = 3
    nbval  = 1
    rubriq = 'num_zone_fb_lwc'
    call restart_read_section_int_t(rp,rubriq,itysup,nbval,izfppp,ierror)
    nberro=nberro+ierror

!       Type entree Gaz brulee
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientgb_zone_bord_lwc'
    call restart_read_section_int_t(rp,rubriq,itysup,nbval,ientgb,ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       Type entree gaz frais
    itysup = 0
    nbval  = nozppm
    rubriq = 'ientgf_zone_bord_lwc'
    call restart_read_section_int_t(rp,rubriq,itysup,nbval,ientgf,ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       FMENT
    itysup = 0
    nbval  = nozppm
    rubriq = 'fment_zone_bord_lwc'
    call restart_read_section_real_t(rp,rubriq,itysup,nbval,fment,ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!       TKENT
    itysup = 0
    nbval  = nozppm
    rubriq = 'tkent_zone_bord_lwc'
    call restart_read_section_real_t(rp,rubriq,itysup,nbval,tkent,ierror)
    ierrch=ierrch+ierror
    nberro=nberro+ierror

!     Par securite, si on ne parvient pas a lire la
!       IENTCPFU ou IENTOX, on remet a zero le numero des zones IZFPPP
!       car il a peut etre ete lu.
!       Ceci permettra d'eviter de se servir des valeurs par defaut

    if (ierrch.ne.0) then
      do ifac = 1, nfabor
        izfppp(ifac) = 0
      enddo
    endif

  endif

endif

!     Charbon PuLVerise : masse vol des charbons
if (ippmod(iccoal).ge.0) then
  itysup = 0
  nbval  = 1
  ierrch = 0
  do icha = 1, ncharb
    if (icha.le.nfmtch) then
      write(car2,'(I2.2)')icha
    else
      car2 = cindfc
    endif
    rubriq = 'masse_volumique_charbon'//car2
    call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
    rhock(icha) = rval(1)
    ierrch = ierrch + ierror
    nberro = nberro + ierror
    ilu = ilu + 1
  enddo
  if (ierrch.ne.0) then
    write(nfecra,8611)
    do icha = 1, ncharb
      write(nfecra,8612)icha,rhock(icha)
    enddo
    write(nfecra,8613)
  endif

endif

if (nberro.ne.0) then
  car54 = 'Reading combustion information                        '
  write(nfecra,8300)car54
endif

if (ilu.ne.0) then
  car54=' Finished reading combustion information.             '
  write(nfecra,1110)car54
endif

!===============================================================================
! 15. LECTURE DES INFORMATIONS COMPLEMENTAIRES ELECTRIQUES
!===============================================================================

nberro=0
ilu  = 0

!     Recalage des CL pot des versions electriques

if ( ippmod(ieljou).ge.1       ) then
  if (ielcor.eq.1) then
    ilu = ilu + 1
    rubriq = 'coeff_recalage_joule'
    itysup = 0
    nbval  = 1
    call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
    coejou = rval(1)
    nberro=nberro+ierror
  endif
endif
if ( ippmod(ielarc).ge.1  .or. ippmod(ieljou).ge.1 ) then
  if (ielcor.eq.1) then
    ilu = 1
    rubriq = 'ddpot_recalage_arc_elec'
    itysup = 0
    nbval  = 1
    rval(1) = pot_diff
    call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
    pot_diff = rval(1)
    nberro=nberro+ierror

    rubriq = 'elcou_recalage_arc_elec'
    rval(1) = elcou
    call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
    elcou = rval(1)
    nberro=nberro+ierror
  endif
endif

! ---> Termes sources des versions electriques

if (nberro.ne.0) then
  car54 = 'Reading electric information                          '
  write(nfecra,8300)car54
endif

if (ilu.ne.0) then
  car54=' Finished reading electric information.               '
  write(nfecra,1110)car54
endif

!===============================================================================
! 16. LECTURE DES CHAMPS PAR CLEF RESTART
!===============================================================================

call restart_read_fields(rp, RESTART_AUXILIARY)

!===============================================================================
! 17.  FERMETURE DU FICHIER SUITE AUXILAIRE
!===============================================================================

!     Fermeture du fichier suite auxilaire
call restart_destroy(rp)

write(nfecra,1200)

!===============================================================================
! 18. SORTIE
!===============================================================================

return

!===============================================================================
! 19. FORMATS
!===============================================================================

! --- ETAPES

 1000 format(/,                                                   &
     3X,'      READING THE AUXILIARY RESTART FILE             ',/)
 1100 format(' Start reading                                          ')
 1110 format('  ', a54                                                 )
 1200 format(' End reading                                            ')

! --- INFORMATIONS

 7000 format(/,                                                   &
'   Apdatation of the reference point for the total pressure  ',/,&
'       by reading the restart file                           ',/,&
'    XYZP0 = ',       E14.5,        E14.5,        E14.5        ,/)

! --- MISES EN GARDE

 8200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHEN READING THE AUXILIARY RESTART FILE        ',/,&
'@    =======                                                 ',/,&
'@      PREVIOUS and PRESENT INPUT DATA ARE DIFFERENT         ',/,&
'@                                                            ',/,&
'@    The number of the faces ',A8  ,' has been modified      ',/,&
'@                                                            ',/,&
'@    The run can continue but the data on the                ',/,&
'@      faces ',A8  ,' will not be reread                     ',/,&
'@      in the suite file.                                    ',/,&
'@    They will be initialized by the default values.         ',/,&
'@                                                            ',/,&
'@     This situation can occur when the restart file         ',/,&
'@      originates from a run using different options         ',/,&
'@      to join the grids or when the periodicity boundary    ',/,&
'@      conditions have been modified.                        ',/,&
'@     This situation can also be generated when the          ',/,&
'@      run is conducted on a different machine               ',/,&
'@      in which case the precision of the machine modifies   ',/,&
'@      the number of faces generated when joinning the grids.',/,&
'@                                                            ',/,&
'@     Finally, this situation can be due to the fact that    ',/,&
'@      the auxiliary restart file does not correspond to     ',/,&
'@      the present case.                                     ',/,&
'@                                                            ',/,&
'@    Verify that the auxiliary restart file being used       ',/,&
'@      corresponds to the present case.                      ',/,&
'@                                                            ',/,&
'@     The run will continue...                               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8300 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHEN READING THE AUXILIARY RESTART FILE        ',/,&
'@    =======                                                 ',/,&
'@      ', a54                                                 ,/,&
'@                                                            ',/,&
'@    It was not possible to read some values from the        ',/,&
'@      auxiliary restart file.                               ',/,&
'@    They will be initialized by the default values.         ',/,&
'@                                                            ',/,&
'@     This situation can occur when the restart file         ',/,&
'@      originates from a run realised with different         ',/,&
'@      options or when the file is damaged.                  ',/,&
'@                                                            ',/,&
'@    The run will continue...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8400 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHEN READING THE AUXILIARY RESTART FILE        ',/,&
'@    =======                                                 ',/,&
'@      THE RUN RESTARTED            WITH IDTVAR = ',I10       ,/,&
'@      FROM RUN CONDUCTED WITH            IDTVAR = ',I10      ,/,&
'@                                                            ',/,&
'@    The variable time step method has been modified.        ',/,&
'@    The (uniform) value of the time step is                 ',/,&
'@      DTREF = ',E12.4                                        ,/,&
'@                                                            ',/,&
'@    It is advised however in this case to                   ',/,&
'@      verify the value of IDTVAR.                           ',/,&
'@                                                            ',/,&
'@    Verify that the auxiliary restart file being used       ',/,&
'@      corresponds  to the present case.                     ',/,&
'@                                                            ',/,&
'@    The run will continue...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8611 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: WHEN READING THE AUXILIARY RESTART FILE        ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@      Combustion model for pulverised coal                  ',/,&
'@      The densities of the coals can not be found           ',/,&
'@        in the restart file. This is normal if the          ',/,&
'@        previous run did not include pulverised coal.       ',/,&
'@        The default value is used as an initial             ',/,&
'@        value :                                             ',/,&
'@         Coal           rho                                 '  )
 8612 format(                                                     &
'@        ',I10   ,'  ',E14.5                                    )
 8613 format(                                                     &
'@                                                            ',/,&
'@    The run will continue...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

! --- ERREURS

 9100 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@      FILE TYPE IS INCORRECT                                ',/,&
'@                                                            ',/,&
'@    The file ',A13      ,' does not appear to be an         ',/,&
'@      auxiliary file.                                       ',/,&
'@                                                            ',/,&
'@    The run can not be executed.                            ',/,&
'@                                                            ',/,&
'@    Verify that the restart file used corresponds to        ',/,&
'@        an auxiliary restart file.                          ',/,&
'@    If necessary, it is possible to deactivate the reading  ',/,&
'@      of the auxiliary restart file by ILEAUX = 0.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9101 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@      INCOHERENT PREVIOUS NAD ACTUAL DATA                   ',/,&
'@                                                            ',/,&
'@    The number of cells was modified                        ',/,&
'@                                                            ',/,&
'@    The run can not be executed.                            ',/,&
'@                                                            ',/,&
'@    Verify that the restart file used corresponds to        ',/,&
'@        the present case.                                   ',/,&
'@    If necessary, it is possible to deactivate the reading  ',/,&
'@      of the auxiliary restart file by ILEAUX = 0.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@      ', a54                                                 ,/,&
'@                                                            ',/,&
'@    The run can not be executed.                            ',/,&
'@                                                            ',/,&
'@    Verify that the restart file used has not been damaged. ',/,&
'@                                                            ',/,&
'@    If necessary, it is possible to deactivate the reading  ',/,&
'@      of the auxiliary restart file by ILEAUX = 0.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9210 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ERROR WHILE READING THE AUXILIARY RESTART FILE ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@      ERROR WHEN READING THE INDICATOR OF THE ALE METHOD    ',/,&
'@                                                            ',/,&
'@    It is possible that the file read corresponds to an old ',/,&
'@      version of code_saturne, without the ALE method.      ',/,&
'@    The run will be executed with reinitialising all        ',/,&
'@      ALE data.                                             ',/,&
'@                                                            ',/,&
'@    Verify that the restart file used has not been damaged. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9221 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ERROR WHILE READING THE AUXILIARY RESTART FILE ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@      ERROR WHEN READING THE INDICATOR OF THE VOF MODEL     ',/,&
'@                                                            ',/,&
'@    It is possible that the file read corresponds to an old ',/,&
'@      version of code_saturne, without the VOF model.       ',/,&
'@    The run will be executed with reinitialising all        ',/,&
'@      VOF model data.                                       ',/,&
'@                                                            ',/,&
'@    Verify that the restart file used has not been damaged. ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9400 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@    Error reading the CO3DP variables                       ',/,&
'@                                                            ',/,&
'@    The run can not be executed.                            ',/,&
'@                                                            ',/,&
'@    Verify the auxiliary restart file                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9500 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@    Error reading the EBU variables                         ',/,&
'@                                                            ',/,&
'@    The run will not be executed.                           ',/,&
'@                                                            ',/,&
'@    Verify the auxiliary restart file                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9600 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: STOP WHILE READING THE AUXILIARY RESTART FILE  ',/,&
'@    =======                                                 ',/,&
'@                                                            ',/,&
'@    Error reading the LWC variables                         ',/,&
'@                                                            ',/,&
'@    The run will not be executed.                           ',/,&
'@                                                            ',/,&
'@    Verify the auxiliary restart file                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

end subroutine
