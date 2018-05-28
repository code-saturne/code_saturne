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

!> \file varpos.f90
!> \brief Variables location initialization, according to calculation type
!> selected by the user.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!______________________________________________________________________________

subroutine varpos

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use cstnum
use entsor
use albase
use lagran
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use radiat
use ihmpre
use mesh
use post
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

integer          iscal , id, ityloc, itycat, ifcvsl, pflag
integer          ii
integer          iok
integer          f_id
integer          ivisph, iest
integer          key_buoyant_id, is_buoyant_fld

double precision gravn2

character(len=80) :: f_name, f_label, s_label, s_name

type(var_cal_opt) :: vcopt

!===============================================================================
! Verification
!===============================================================================

! ALE viscosity
if (iale.eq.1) then
  if (iortvm.ne.0 .and. iortvm.ne.1) then
    write(nfecra, 8022) 'iortvm', iortvm
    call csexit (1)
  endif
endif

!===============================================================================
! Initialization
!===============================================================================

! Key id for buoyant field (inside the Navier Stokes loop)
call field_get_key_id("is_buoyant", key_buoyant_id)

! Determine itycor now that irccor is known (iturb/itytur known much earlier)
! type of rotation/curvature correction for turbulent viscosity models
if (irccor.eq.1.and.(itytur.eq.2.or.itytur.eq.5)) then
  itycor = 1
else if (irccor.eq.1.and.(iturb.eq.60.or.iturb.eq.70)) then
  itycor = 2
endif

pflag = POST_ON_LOCATION + POST_MONITOR

!===============================================================================
! Additional physical properties
!===============================================================================

! CP when variable
if (icp.ge.0) then
  call add_property_field_1d('specific_heat', 'Specific Heat', icp)

  call field_set_key_int(icp, keyvis, 1)
  call field_set_key_int(icp, keylog, 1)
endif

! ALE mesh viscosity
if (iale.eq.1) then
  if (iortvm.eq.0) then
    call add_property_field('mesh_viscosity', 'Mesh Visc', 1, .false., ivisma)
  else
    call add_property_field('mesh_viscosity', 'Mesh Visc', 3, .false., ivisma)
  endif
endif

if (irccor.eq.1) then
  if (idtvar.ge.0) then
     call add_property_field('strain_rate_tensor', 'Strain Rate Tensor', 6, &
                             .false., istraio)
     call hide_property(istraio)
  endif
endif

!===============================================================================
! Time-scheme related properties
!===============================================================================

! Dans le cas ou on a un schema en temps d'ordre 2, il faut aussi
!   prevoir les proprietes au temps n-1. Ce sera fait au dernier appel

! Dans le cas ou on calcule des moments, il faut en prevoir le nombre
!   et prevoir le nombre de tableaux necessaires pour stocker le
!   temps cumule de moyenne. On suppose que l'on peut s'aider
!   d'infos situees en tete de fichier suite (si on fait une
!   suite avec des moments non reinitialises).

! 1.1 PROPRIETES ADDITIONNELLES POUR LES ET SCHEMA EN TEMPS
! ---------------------------------------------------------

! Initialisations par defaut eventuelles et verifications
!   des options utilisees ci-dessous pour decider si l'on
!   reserve des tableaux supplementaires pour des grandeurs
!   au pas de temps precedent

iok = 0

!  Pression hydrostatique
if (iphydr.eq.0.or.iphydr.eq.2) then
  icalhy = 0
else if (iphydr.eq.1.and.icalhy.eq.-1) then
  gravn2 = gx**2+gy**2+gz**2
  if (gravn2.lt.epzero**2) then
    icalhy = 0
  else
    icalhy = 1
  endif
endif

!   Schemas en temps
!       en LES : Ordre 2 ; sinon Ordre 1
!       (en particulier, ordre 2 impossible en k-eps couple)
if (ischtp.eq.-999) then
  if (itytur.eq.4) then
    ischtp = 2
  else
    ischtp = 1
  endif
endif

!   Schemas en temps : variables deduites
!   Schema pour le Flux de masse
if (istmpf.eq.-999) then
  if (ischtp.eq.1) then
    istmpf = 1
  else if (ischtp.eq.2) then
    istmpf = 2
  endif
endif
!     Masse volumique
if (iroext.eq.-999) then
  if (ischtp.eq.1) then
    iroext = 0
  else if (ischtp.eq.2) then
    !       Pour le moment par defaut on ne prend pas l'ordre 2
    !              IROEXT = 1
    iroext = 0
  endif
endif
!     Viscosite
if (iviext.eq.-999) then
  if (ischtp.eq.1) then
    iviext = 0
  else if (ischtp.eq.2) then
    !       Pour le moment par defaut on ne prend pas l'ordre 2
    !              IVIEXT = 1
    iviext = 0
  endif
endif
!     Chaleur massique
if (icpext.eq.-999) then
  if (ischtp.eq.1) then
    icpext = 0
  else if (ischtp.eq.2) then
    !       Pour le moment par defaut on ne prend pas l'ordre 2
    !              ICPEXT = 1
    icpext = 0
  endif
endif
!     Termes sources NS,
if (isno2t.eq.-999) then
  if (ischtp.eq.1) then
    isno2t = 0
    !            ELSE IF (ISCHTP.EQ.2.AND.IVISSE.EQ.1) THEN
  else if (ischtp.eq.2) then
    !       Pour le moment par defaut on prend l'ordre 2
    isno2t = 1
    !              ISNO2T = 0
  endif
endif
!     Termes sources turbulence (k-eps, Rij, v2f ou k-omega)
!     On n'autorise de changer ISTO2T qu'en Rij (sinon avec
!       le couplage k-eps/omega il y a pb)
if (isto2t.eq.-999) then
  if (ischtp.eq.1) then
    isto2t = 0
  else if (ischtp.eq.2) then
    !       Pour le moment par defaut on ne prend pas l'ordre 2
    !              ISTO2T = 1
    isto2t = 0
  endif
else if (itytur.eq.2.or.iturb.eq.50.or.iturb.ne.60) then
  write(nfecra,8132) iturb,isto2t
  iok = iok + 1
endif

do iscal = 1, nscal
  ! Termes sources Scalaires,
  if (isso2t(iscal).eq.-999) then
    if (ischtp.eq.1) then
      isso2t(iscal) = 0
    else if (ischtp.eq.2) then
      ! Pour coherence avec Navier Stokes on prend l'ordre 2
      ! mais de toute facon qui dit ordre 2 dit LES et donc
      ! generalement pas de TS scalaire a interpoler.
      isso2t(iscal) = 1
    endif
  endif
  ! Diffusivite scalaires
  if (ivsext(iscal).eq.-999) then
    if (ischtp.eq.1) then
      ivsext(iscal) = 0
    else if (ischtp.eq.2) then
      ! Pour le moment par defaut on ne prend pas l'ordre 2
      ivsext(iscal) = 0
    endif
  endif

  ! Model for turbulent fluxes u'T' (SGDH, GGDH, AFM, DFM)
  ityturt(iscal) = iturt(iscal)/10

  if (iscal.eq.iscalt) then
    if (iturt(iscalt).gt.0.and.irovar.eq.1) then
      call add_property_field_1d('thermal_expansion', 'Beta', ibeta)
    endif
  endif

enddo

! Pression hydrostatique
if (iphydr.ne.0.and.iphydr.ne.1.and.iphydr.ne.2) then
  write(nfecra,8121) 'IPHYDR ',iphydr
  iok = iok + 1
endif

! Viscosite secondaire
ivisph = ivisse
if (ivisph.ne.0.and.ivisph.ne.1) then
  write(nfecra,8022) 'IVISSE ',ivisph
  iok = iok + 1
endif

! Schemas en temps

!     Schema en temps global.
if (ischtp.ne. 1.and.ischtp.ne.2) then
  write(nfecra,8101) 'ISCHTP',ischtp
  iok = iok + 1
endif
if (ischtp.eq. 2.and.idtvar.ne.0) then
  write(nfecra,8111) ischtp,idtvar
  iok = iok + 1
endif
if (ischtp.eq. 2.and.itytur.eq.2) then
  write(nfecra,8112) ischtp,iturb
  iok = iok + 1
endif
if (ischtp.eq.1.and.itytur.eq.4) then
  write(nfecra,8113) ischtp,iturb
endif
if (ischtp.eq. 2.and.iturb.eq.50) then
  write(nfecra,8114) ischtp,iturb
  iok = iok + 1
endif
if (ischtp.eq. 2.and.iturb.eq.51) then
  write(nfecra,8117) ischtp,iturb
  iok = iok + 1
endif
if (ischtp.eq. 2.and.iturb.eq.60) then
  write(nfecra,8115) ischtp,iturb
  iok = iok + 1
endif
if (ischtp.eq. 2.and.iturb.eq.70) then
  write(nfecra,8116) ischtp,iturb
  iok = iok + 1
endif

! Schema en temps pour le flux de masse
if (istmpf.ne. 2.and.istmpf.ne.0.and.istmpf.ne. 1) then
  write(nfecra,8121) 'ISTMPF',istmpf
  iok = iok + 1
endif

! Schema en temps pour les termes sources de NS
if (isno2t.ne.0.and.isno2t.ne. 1.and.isno2t.ne.2) then
  write(nfecra,8131) 'ISNO2T',isno2t
  iok = iok + 1
endif
! Schema en temps pour les termes sources des grandeurs turbulentes
if (isto2t.ne.0.and.isto2t.ne. 1.and.isto2t.ne.2) then
  write(nfecra,8131) 'ISTO2T',isto2t
  iok = iok + 1
endif

! Schema en temps pour la masse volumique
if (iroext.ne.0.and.iroext.ne. 1.and.iroext.ne.2) then
  write(nfecra,8131) 'IROEXT',iroext
  iok = iok + 1
endif

! Schema en temps pour la viscosite
if (iviext.ne.0.and.iviext.ne. 1.and.iviext.ne.2) then
  write(nfecra,8131) 'IVIEXT',iviext
  iok = iok + 1
endif

! Schema en temps pour la chaleur specifique
if (icpext.ne.0 .and. icpext.ne.1 .and.icpext.ne.2) then
  write(nfecra,8131) 'ICPEXT',icpext
  iok = iok + 1
endif

do iscal = 1, nscal
  ! Schema en temps pour les termes sources des scalaires
  if (isso2t(iscal).ne.0.and.                                    &
      isso2t(iscal).ne. 1.and.isso2t(iscal).ne.2) then
    write(nfecra,8141) iscal,'ISSO2T',isso2t(iscal)
    iok = iok + 1
  endif
  ! Schema en temps pour la viscosite
  if (ivsext(iscal).ne.0.and.                                    &
      ivsext(iscal).ne. 1.and.ivsext(iscal).ne.2) then
    write(nfecra,8141) iscal,'IVSEXT',ivsext(iscal)
    iok = iok + 1
  endif
enddo

! Stop si probleme
if (iok.gt.0) then
  call csexit(1)
endif

! Source term for weakly compressible algorithm (semi analytic scheme)
if (idilat.ge.4) then
  do iscal = 1, nscal
    id = ivarfl(isca(iscal))
    call field_get_name(id, s_name)
    call field_get_label(id, s_label)
    f_name  = trim(s_name) // '_dila_st'
    call add_property_field_1d(f_name, '', iustdy(iscal))
    id = iustdy(iscal)
    call field_set_key_int(id, keyvis, 0)
  enddo
  itsrho = nscal + 1
  call add_property_field_1d('dila_st', '', iustdy(itsrho))
  id = iustdy(iscal)
  call field_set_key_int(id, keyvis, 0)
endif

! Density at the second previous time step for VOF algorithm
! or dilatable algorithm
if (ivofmt.ge.0.or.idilat.gt.1) then
  call field_set_n_previous(icrom, 2)
  call field_set_n_previous(ibrom, 2)
  ! The density at the previous time step is required if
  ! we perform a hydrostatic pressure correction (icalhy=1)
else if (iroext.gt.0.or.icalhy.eq.1.or.ipthrm.eq.1.or.ippmod(icompf).ge.0) then
  call field_set_n_previous(icrom, 1)
  call field_set_n_previous(ibrom, 1)
endif
! Dans le cas d'une extrapolation de la viscosite totale
if (iviext.gt.0) then
  call field_set_n_previous(iviscl, 1)
  call field_set_n_previous(ivisct, 1)
endif

! CP s'il est variable
if (icp.ge.0) then
  if (icpext.gt.0) then
    call field_set_n_previous(icp, 1)
  endif
endif

! On a besoin d'un tableau pour les termes sources de Navier Stokes
!  a extrapoler. Ce tableau est NDIM dans le cas general et NDIM+1
!  si on extrapole aussi les termes sources de l equation sur le taux
!  de vide pour l'algo. VOF.
if (isno2t.gt.0) then
  call add_source_term_prev_field(ivarfl(iu))
  if (ivofmt.ge.0) then
    call add_source_term_prev_field(ivarfl(ivolf2))
  endif
endif

if (isto2t.gt.0) then
  ! The dimension of this array depends on turbulence model:
  if (itytur.eq.2) then
    call add_source_term_prev_field(ivarfl(ik))
    call add_source_term_prev_field(ivarfl(iep))
  else if (itytur.eq.3) then
    call add_source_term_prev_field(ivarfl(ir11))
    call add_source_term_prev_field(ivarfl(ir22))
    call add_source_term_prev_field(ivarfl(ir33))
    call add_source_term_prev_field(ivarfl(ir12))
    call add_source_term_prev_field(ivarfl(ir13))
    call add_source_term_prev_field(ivarfl(ir23))
    call add_source_term_prev_field(ivarfl(iep))
    if (iturb.eq.32) then
      call add_source_term_prev_field(ivarfl(ial))
    endif
  else if (itytur.eq.5) then
    call add_source_term_prev_field(ivarfl(ik))
    call add_source_term_prev_field(ivarfl(iep))
    call add_source_term_prev_field(ivarfl(iphi))
    if (iturb.eq.50) then
      call add_source_term_prev_field(ivarfl(ifb))
    else if (iturb.eq.51) then
      call add_source_term_prev_field(ivarfl(ial))
    endif
  else if (iturb.eq.60) then
    call add_source_term_prev_field(ivarfl(ik))
    call add_source_term_prev_field(ivarfl(iomg))
  else if (iturb.eq.70) then
    call add_source_term_prev_field(ivarfl(inusa))
  endif
endif

! Proprietes des scalaires : termes sources pour theta schema
!   et VISCLS si elle est variable
if (nscal.ge.1) then
  do ii = 1, nscal
    if (isso2t(ii).gt.0) then
      ! For buoyant scalars, save the current user source term
      call field_get_key_int(ivarfl(isca(ii)), key_buoyant_id, is_buoyant_fld)
      if (is_buoyant_fld.eq.1) then
        call add_source_term_field(ivarfl(isca(ii)))
      endif
      call add_source_term_prev_field(ivarfl(isca(ii)))
    endif
    ! Only usefull for Min/Max limiter
    call field_get_key_struct_var_cal_opt(ivarfl(isca(ii)), vcopt)
    if (vcopt%isstpc.eq.2) then
      call add_source_term_field(ivarfl(isca(ii)))
    endif
    ! Diffusivity
    call field_get_key_int (ivarfl(isca(ii)), kivisl, ifcvsl)
    if (ifcvsl.ge.0.and.iscavr(ii).le.0) then
      if (ivsext(ii).gt.0) then
        call field_set_n_previous(ifcvsl, 1)
      endif
    endif
    ! Density
    call field_get_key_int (ivarfl(isca(ii)), kromsl, ifcvsl)
    if (ifcvsl.ge.0.and.iscavr(ii).le.0) then
      if (ivsext(ii).gt.0) then
        call field_set_n_previous(ifcvsl, 1)
      endif
    endif
  enddo
endif

! Porosity
ityloc = 1 ! cells
itycat = FIELD_INTENSIVE + FIELD_PROPERTY

if (iporos.ge.1) then
  f_name = 'porosity'
  call field_create(f_name, itycat, ityloc, 1, .false., ipori)
  call field_set_key_int(ipori, keylog, 1)
  call field_set_key_int(ipori, keyvis, pflag)
  if (iporos.eq.2) then
    f_name = 'tensorial_porosity'
    call field_create(f_name, itycat, ityloc, 6, .false., iporf)
  endif
  if (iporos.eq.3) then
    f_name = 'poro_div_duq'
    call field_create(f_name,&
                      itycat,&
                      1,& ! location: cell
                      3,& ! dimension
                      .false.,&
                      f_id)

    f_name = 'i_poro_duq_0'
    call field_create(f_name,&
                      itycat,&
                      2,& ! location: inner faces
                      1,& ! dimension
                      .false.,&
                      f_id)

    f_name = 'i_poro_duq_1'
    call field_create(f_name,&
                      itycat,&
                      2,& ! location: inner faces
                      1,& ! dimension
                      .false.,&
                      f_id)

    f_name = 'b_poro_duq'
    call field_create(f_name,&
                      itycat,&
                      3,& ! location: boundary faces
                      1,& ! dimension
                      .false.,&
                      f_id)

  endif
endif

!===============================================================================
! Local time step and postprocessing fields
!===============================================================================

! Local time step

ityloc = 1 ! cells
itycat = FIELD_INTENSIVE

call field_create('dt', itycat, ityloc, 1, .false., id)
call field_set_key_str(id, keylbl, 'Local Time Step')
if (idtvar.gt.0) then
  if (idtvar.eq.2) then
    call field_set_key_int(id, keylog, 1)
    call field_set_key_int(id, keyvis, pflag)
  endif
endif

itycat = FIELD_INTENSIVE

! Transient velocity/pressure coupling, postprocessing field
! (variant used for computation is a tensorial field, not this one)

if (ipucou.ne.0 .or. ncpdct.gt.0) then
  call field_create('dttens', itycat, ityloc, 6, .false., idtten)
  call field_set_key_int(idtten, keyvis, POST_ON_LOCATION)
  call field_set_key_int(idtten, keylog, 1)
endif

! Error estimators

do iest = 1, nestmx
  iestim(iest) = -1
enddo

if (iescal(iespre).gt.0) then
  write(f_name,  '(a14,i1)') 'est_error_pre_', iescal(iespre)
  write(f_label, '(a5,i1)') 'EsPre', iescal(iespre)
  call add_property_field(f_name, f_label, 1, .false., iestim(iespre))
endif
if (iescal(iesder).gt.0) then
  write(f_name,  '(a14,i1)') 'est_error_der_', iescal(iesder)
  write(f_label, '(a5,i1)') 'EsDer', iescal(iesder)
  call add_property_field(f_name, f_label, 1, .false., iestim(iesder))
endif
if (iescal(iescor).gt.0) then
  write(f_name,  '(a14,i1)') 'est_error_cor_', iescal(iescor)
  write(f_label, '(a5,i1)') 'EsCor', iescal(iescor)
  call add_property_field(f_name, f_label, 1, .false., iestim(iescor))
endif
if (iescal(iestot).gt.0) then
  write(f_name,  '(a14,i1)') 'est_error_tot_', iescal(iestot)
  write(f_label, '(a5,i1)') 'EsTot', iescal(iestot)
  call add_property_field(f_name, f_label, 1, .false., iestim(iestot))
endif

!===============================================================================
! Map to field pointers
!===============================================================================

call cs_field_pointer_map_base
call cs_field_pointer_map_boundary

return

!===============================================================================
! Formats
!===============================================================================

#if defined(_CS_LANG_FR)

 8022 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0 OU 1                ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8101 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 1 ou 2                ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8111 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    AVEC UN SCHEMA EN TEMPS D ORDRE 2 : ISCHTP = ', I10      ,/,&
'@    IL FAUT UTILISER UN PAS DE TEMPS CONSTANT ET UNIFORME   ',/,&
'@    OR IDTVAR = ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8112 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,&
'@    EN K-EPSILON (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   La version courante ne supporte pas l ordre 2 avec le    ',/,&
'@   couplage des termes sources du k-epsilon.                ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8113 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 1 (ISCHTP = ',I10   ,/,&
'@    EN LES (ITURB = ',I10,' )'                               ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute.                                   ',/,&
'@                                                            ',/,&
'@  Il est conseille de verifier les parametres.              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8114 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,&
'@    EN PHI_FBAR (ITURB = ',I10,' )'                          ,/,&
'@                                                            ',/,&
'@   La version courante ne supporte pas l''ordre 2 avec le   ',/,&
'@   couplage des termes sources du k-epsilon.                ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8117 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,&
'@    EN BL-V2/K  (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   La version courante ne supporte pas l''ordre 2 avec le   ',/,&
'@   couplage des termes sources du k-epsilon.                ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8115 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,&
'@    EN K-OMEGA   (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   La version courante ne supporte pas l''ordre 2 avec le   ',/,&
'@   couplage des termes sources du k-omega.                  ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8116 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,&
'@    EN SPALART   (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   La version courante ne supporte pas l''ordre 2 avec le   ',/,&
'@   couplage des termes sources de Spalart-Allmaras.         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8121 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A  0, 1 OU 2            ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8131 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0, 1 OU 2             ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8132 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  Avec le modele de turbulence choisi, ITURB = ',I10         ,/,&
'@    la valeur de ISTO2T (extrapolation des termes sources   ',/,&
'@    pour les variables turbulentes) ne doit pas etre modifie',/,&
'@    or ISTO2T a ete force a ',I10                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8141 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES SCALAIRE ',I10 ,/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0, 1 OU 2             ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 8022 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0 OR 1               ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8101 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 1 OR 2               ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verifier parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8111 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    WITH A SECOND ORDER SCHEME IN TIME: ISCHTP = ', I10      ,/,&
'@    IT IS NECESSARY TO USE A CONSTANT AND UNIFORM TIME STEP ',/,&
'@    BUT IDTVAR = ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verifier parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8112 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2ND ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    WITH K-EPSILON (ITURB = ',I10,' )'                       ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-epsilon.               ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8113 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   :      AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 1st ORDER SCHEME HAS BEEN IMPOSSED   (ISCHTP = ',I10   ,/,&
'@    FOR LES (ITURB = ',I10,' )'                              ,/,&
'@                                                            ',/,&
'@  The calculation will   be executed                        ',/,&
'@                                                            ',/,&
'@  It is recommended to verify  parameters.                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8114 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR PHI_FBAR (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-epsilon.               ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8117 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA                    ',/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR BL-V2/K  (ITURB = ',I10,' )'                        ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-epsilon.               ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8115 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR K-OMEGA   (ITURB = ',I10,' )'                        ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-omega.                 ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8116 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR SPALART   (ITURB = ',I10,' )'                        ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of Spalart-Allmaras.        ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8121 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0, 1 OR 2            ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8131 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0, 1 OR 2            ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8132 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  With the chosen turbulence model   , ITURB = ',I10         ,/,&
'@    the value of ISTO2T (extrapolation of the source terms  ',/,&
'@    for the turbulent variables) cannot be modified         ',/,&
'@    yet ISTO2T has been forced to ',I10                      ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8141 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITAL DATA FOR SCALARS    ',I10 ,/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0, 1 OR 2            ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!===============================================================================
! 5. End
!===============================================================================

return
end subroutine varpos

!===============================================================================

!> \brief add field defining previous source term values for a given field
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     f_id          base field id
!_______________________________________________________________________________

subroutine add_source_term_prev_field &
 ( f_id )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field

!===============================================================================

implicit none

! Arguments

integer, intent(in) :: f_id

! Local variables

character(len=64) :: f_name

integer :: type_flag, location_id, st_id, f_dim
logical :: has_previous

!===============================================================================

type_flag = FIELD_EXTENSIVE + FIELD_PROPERTY
location_id = 1 ! variables defined on cells
has_previous = .false.

! Define asscociated field

call field_get_dim(f_id, f_dim)
call field_get_name (f_id, f_name)

call field_create(trim(f_name)//'_prev_st', type_flag,               &
                  location_id, f_dim, has_previous, st_id)

call field_set_key_int(f_id, kstprv, st_id)

return

end subroutine add_source_term_prev_field


!===============================================================================

!> \brief add field defining current source term values for a given field
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     f_id          base field id
!_______________________________________________________________________________

subroutine add_source_term_field &
 ( f_id )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field

!===============================================================================

implicit none

! Arguments

integer, intent(in) :: f_id

! Local variables

character(len=64) :: f_name

integer :: type_flag, location_id, st_id, f_dim
logical :: has_previous

!===============================================================================

type_flag = FIELD_EXTENSIVE + FIELD_PROPERTY
location_id = 1 ! variables defined on cells
has_previous = .false.

! Define asscociated field

call field_get_dim(f_id, f_dim)
call field_get_name (f_id, f_name)

call field_create(trim(f_name)//'_st', type_flag,               &
                  location_id, f_dim, has_previous, st_id)

call field_set_key_int(f_id, kst, st_id)

return

end subroutine add_source_term_field
