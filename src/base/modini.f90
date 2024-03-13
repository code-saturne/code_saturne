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

!> \file modini.f90
!> \brief Modify calculation parameters after user changes (module variables)
!>
!------------------------------------------------------------------------------

subroutine modini

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use field
use numvar
use optcal
use cstphy
use entsor
use albase
use cplsat
use post
use ppincl
use rotation
use turbomachinery
use vof
use cs_c_bindings

!===============================================================================

implicit none

! Arguments


! Local variables

integer          f_id, f_id_d
integer          ii, jj, iok
integer          nbccou
integer          nscacp, iscal
integer          iclvfl, kclvfl
integer          iscacp, kcpsyr, icpsyr
integer          nfld, f_type
integer          key_t_ext_id, icpext, kscmin, kscmax
integer          iviext, isso2t, kisso2t, kthetss, kthetvs, kcdtvar
integer          kturt

double precision relxsp, clvfmn, clvfmx, visls_0, visls_cmp
double precision scminp, thetss, thetvs, cdtvar

character(len=80) :: name

type(var_cal_opt) :: vcopt , vcopt1

procedure() :: hide_property, indsui, nbccpl

!===============================================================================

! Indicateur erreur (0 : pas d'erreur)
iok = 0

call field_get_key_id("syrthes_coupling", kcpsyr)

call field_get_key_id("time_extrapolated", key_t_ext_id)

call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)
call field_get_key_id('turbulent_flux_model', kturt)
call field_get_key_id("st_exp_extrapolated", kthetss)
call field_get_key_id("diffusivity_extrapolated", kthetvs)
call field_get_key_id("scalar_time_scheme", kisso2t)

call field_get_key_id("variance_clipping", kclvfl)
call field_get_key_id("time_step_factor", kcdtvar)

!===============================================================================
! 1. ENTREES SORTIES entsor
!===============================================================================

! Logging and postprocessing output

if (irovar.eq.0) then
  call hide_property(icrom)
  call hide_property(ibrom)
endif

if (ivivar.eq.0) then
  call hide_property(iviscl)
endif

!===============================================================================
! 2. OPTIONS DU CALCUL : TABLEAUX DE optcal
!===============================================================================

! time scheme

if (ntmabs.eq.-1 .and. ttmabs.lt.-0.5) then
  ntmabs = 10
endif

! restart

call indsui(isuite)

if (isuit1.eq.-1) isuit1 = isuite

!    -- Proprietes physiques
call field_get_key_int(iviscl, key_t_ext_id, iviext)
if (abs(thetvi+999.d0).gt.epzero) then
  write(nfecra,1011) 'IVIEXT', iviext, 'THETVI'
  iok = iok + 1
elseif (iviext.eq.0) then
  thetvi = 0.0d0
elseif (iviext.eq.1) then
  thetvi = 0.5d0
elseif (iviext.eq.2) then
  thetvi = 1.d0
endif

if (icp.ge.0) then
  call field_get_key_int(icp, key_t_ext_id, icpext)
  if (abs(thetcp+999.d0).gt.epzero) then
    write(nfecra,1011) 'ICPEXT', icpext, 'THETCP'
    iok = iok + 1
  elseif (icpext.eq.0) then
    thetcp = 0.0d0
  elseif (icpext.eq.1) then
    thetcp = 0.5d0
  elseif (icpext.eq.2) then
    thetcp = 1.d0
  endif
endif

!    -- Termes sources NS
if (abs(thetsn+999.d0).gt.epzero) then
  write(nfecra,1011) 'ISNO2T', isno2t, 'THETSN'
  iok = iok + 1
elseif (isno2t.eq.1) then
  thetsn = 0.5d0
elseif (isno2t.eq.2) then
  thetsn = 1.d0
elseif (isno2t.eq.0) then
  thetsn = 0.d0
endif

!    -- Termes sources grandeurs turbulentes
if (abs(thetst+999.d0).gt.epzero) then
  write(nfecra,1011) 'ISTO2T', isto2t, 'THETST'
  iok = iok + 1
elseif (isto2t.eq.1) then
  thetst = 0.5d0
elseif (isto2t.eq.2) then
  thetst = 1.d0
elseif (isto2t.eq.0) then
  thetst = 0.d0
endif

do iscal = 1, nscal
!    -- Termes sources des scalaires
  f_id = ivarfl(isca(iscal))
  call field_get_key_int(f_id, kisso2t, isso2t)
  call field_get_key_double(f_id, kthetss, thetss)
  if (abs(thetss+1.d0).gt.epzero) then
    write(nfecra,1021) f_id, 'ISSO2T', isso2t, 'THETSS'
    iok = iok + 1
  elseif (isso2t.eq.1) then
    thetss = 0.5d0
    call field_set_key_double(f_id, kthetss, thetss)
  elseif (isso2t.eq.2) then
    thetss = 1.d0
    call field_set_key_double(f_id, kthetss, thetss)
  elseif (isso2t.eq.0) then
    thetss = 0.d0
    call field_set_key_double(f_id, kthetss, thetss)
  endif
  ! Scalars diffusivity
  call field_get_key_int(f_id, kivisl, f_id_d)
  if (f_id_d.ge.0) then
    call field_get_key_int(f_id_d, key_t_ext_id, iviext)
  else
    iviext = 0
  endif
  call field_get_key_double(f_id, kthetvs, thetvs)
  if (abs(thetvs+1.d0).gt.epzero) then
    write(nfecra,1021) iscal, 'IVSEXT', iviext, 'THETVS'
    iok = iok + 1
  elseif (iviext.eq.0) then
    thetvs = 0.d0
    call field_set_key_double(f_id, kthetvs, thetvs)
  elseif (iviext.eq.1) then
    thetvs = 0.5d0
    call field_set_key_double(f_id, kthetvs, thetvs)
  elseif (iviext.eq.2) then
    thetvs = 1.d0
    call field_set_key_double(f_id, kthetvs, thetvs)
  endif
enddo

! Loop on on field variables
call field_get_n_fields(nfld)

do f_id = 0, nfld - 1
  call field_get_type(f_id, f_type)
  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    call field_get_key_struct_var_cal_opt(f_id, vcopt)
    if (abs(vcopt%thetav+1.d0).gt.epzero) then
      call field_get_name(f_id, name)
      write(nfecra,1131) trim(name),'THETAV'
    else
      ! For the pressure, no theta-scheme
      if (f_id.eq.ivarfl(ipr)) then
        vcopt%thetav = 1.d0
      else if (vcopt%istat.eq.0) then
        vcopt%thetav = 1.d0
      else if (ischtp.eq.1) then
        vcopt%thetav = 1.d0
      else if (ischtp.eq.2) then
        vcopt%thetav = 0.5d0
      endif
    endif
    call field_set_key_struct_var_cal_opt(f_id, vcopt)
  endif
enddo

! Diffusivity model:
! Daly Harlow (GGDH) on Rij and epsilon by default
if (itytur.eq.3) then

  call field_get_key_struct_var_cal_opt(ivarfl(irij), vcopt1)
  call field_get_key_struct_var_cal_opt(ivarfl(iep), vcopt)

  ! Diffusivity model:
  ! Daly Harlow (GGDH) on Rij and epsilon by default
  if (idirsm.ne.0) then
     vcopt1%idften = ANISOTROPIC_RIGHT_DIFFUSION
     vcopt%idften  = ANISOTROPIC_RIGHT_DIFFUSION
     ! Scalar diffusivity (Shir model) elswhere (idirsm = 0)
  else
     vcopt1%idften = ISOTROPIC_DIFFUSION
     vcopt%idften  = ISOTROPIC_DIFFUSION
  endif

  call field_set_key_struct_var_cal_opt(ivarfl(irij), vcopt1)
  call field_set_key_struct_var_cal_opt(ivarfl(iep), vcopt)
endif

! ---> ISSTPC
!        Si l'utilisateur n'a rien specifie pour le test de pente (=-1),
!        On impose 1 (ie sans) pour la vitesse en LES
!                  0 (ie avec) sinon

if (itytur.eq.4) then
  ii = iu
  call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  if (vcopt%isstpc.eq.-999) then
    vcopt%isstpc = 1
    call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  endif
  do jj = 1, nscal
    ii = isca(jj)
    call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)
    if (vcopt%isstpc.eq.-999) then
      vcopt%isstpc = 0
      call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)
    endif
 enddo
endif

do f_id = 0, nfld - 1
  call field_get_type(f_id, f_type)
  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    call field_get_key_struct_var_cal_opt(f_id, vcopt)
    if (vcopt%isstpc.eq.-999) then
      vcopt%isstpc = 0
      call field_set_key_struct_var_cal_opt(f_id, vcopt)
    endif
  endif
enddo

! ---> BLENCV
!        Si l'utilisateur n'a rien specifie pour le schema convectif
!                  1 (ie centre) pour les vitesses
!                                     les scalaires utilisateurs
!                                     le scalaire thermique
!                  0 (ie upwind pur) pour le reste
!   (en particulier, en L.E.S. toutes les variables sont donc en centre)

!  Pour le modele de cavitation on force dans tous les cas le taux de vide en
!  upwind et on affiche un message si l'utilisateur avait specifie autre chose

ii = iu
call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)
if (abs(vcopt%blencv+1.d0).lt.epzero) then
  vcopt%blencv = 1.d0
  call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)
endif

if (iand(ivofmt,VOF_MERKLE_MASS_TRANSFER).ne.0) then
  ii = ivolf2
  call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  if (abs(vcopt%blencv+1.d0).lt.epzero) then
    if (abs(vcopt%blencv+1.d0).gt.epzero) &
         write(nfecra,3000) 0.d0, vcopt%blencv
    vcopt%blencv = 0.d0
    call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  endif
else if (ivofmt.gt.0) then
  ii = ivolf2
  call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  if (abs(vcopt%blencv+1.d0).lt.epzero) then
    vcopt%blencv = 1.d0
    call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  endif
endif

do jj = 1, nscaus
  ii = isca(jj)
  call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  if (abs(vcopt%blencv+1.d0).lt.epzero) then
    vcopt%blencv = 1.d0
    call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  endif
enddo

if (iscalt.gt.0) then
  ii = isca(iscalt)
  call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  if (abs(vcopt%blencv+1.d0).lt.epzero) then
    vcopt%blencv = 1.d0
    call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  endif
endif

do f_id = 0, nfld - 1
  call field_get_type(f_id, f_type)
  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    call field_get_key_struct_var_cal_opt(f_id, vcopt)
    if (abs(vcopt%blencv+1.d0).lt.epzero) then
      vcopt%blencv = 0.d0
      call field_set_key_struct_var_cal_opt(f_id, vcopt)
    endif
  endif
enddo

! ---> NSWRSM, EPSRSM ET EPSILO
!        Si l'utilisateur n'a rien specifie  (NSWRSM=-1),
!        On impose
!           a l'ordre 1 :
!                  2  pour la pression
!                  1  pour les autres variables
!                  on initialise EPSILO a 1.d-8
!                     pour la pression
!                  on initialise EPSILO a 1.d-5
!                     pour les autres variables
!                  on initialise EPSRSM a 10*EPSILO
!           a l'ordre 2 :
!                  5  pour la pression
!                  10 pour les autres variables
!                  on initialise EPSILO a 1.D-5
!                  on initialise EPSRSM a 10*EPSILO
!     Attention aux tests dans verini

if (ischtp.eq.2) then
  ii = ipr
  call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  if (vcopt%nswrsm.eq.-1) vcopt%nswrsm = 5
  if (abs(vcopt%epsilo+1.d0).lt.epzero) vcopt%epsilo = 1.d-5
  if (abs(vcopt%epsrsm+1.d0).lt.epzero) vcopt%epsrsm = 10.d0*vcopt%epsilo
  call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  ii = iu
  call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  if (vcopt%nswrsm.eq.-1) vcopt%nswrsm = 10
  if (abs(vcopt%epsilo+1.d0).lt.epzero) vcopt%epsilo = 1.d-5
  if (abs(vcopt%epsrsm+1.d0).lt.epzero) vcopt%epsrsm = 10.d0*vcopt%epsilo
  call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  do jj = 1, nscal
    ii = isca(jj)
    call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)
    if (vcopt%nswrsm.eq.-1) vcopt%nswrsm = 10
    if (abs(vcopt%epsilo+1.d0).lt.epzero) vcopt%epsilo = 1.d-5
    if (abs(vcopt%epsrsm+1.d0).lt.epzero) vcopt%epsrsm = 10.d0*vcopt%epsilo
    call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  enddo
endif

! For the pressure, default solver precision 1e-8
! because the mass conservation is up to this precision
ii = ipr
call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)
if (vcopt%nswrsm.eq.-1) vcopt%nswrsm = 2
if (abs(vcopt%epsilo+1.d0).lt.epzero) vcopt%epsilo = 1.d-8
call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)

do f_id = 0, nfld - 1
  call field_get_type(f_id, f_type)
  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    call field_get_key_struct_var_cal_opt(f_id, vcopt)
    if (vcopt%nswrsm.eq.-1) vcopt%nswrsm = 1
    if (abs(vcopt%epsilo+1.d0).lt.epzero) vcopt%epsilo = 1.d-5
    if (abs(vcopt%epsrsm+1.d0).lt.epzero) vcopt%epsrsm = 10.d0*vcopt%epsilo
    call field_set_key_struct_var_cal_opt(f_id, vcopt)
  endif
enddo

! ---> dtmin dtmax cdtvar

if (dtmin.le.-grand) then
  dtmin = 0.1d0*dtref
endif
if (dtmax.le.-grand) then
  dtmax = 1.0d3*dtref
endif

! Init. of time step factor for velocity, pressure and turbulent variables
! FIXME time step factor is used ONLY for additional variables (user or model)

call field_get_key_double(ivarfl(iu), kcdtvar, cdtvar)
call field_set_key_double(ivarfl(iv), kcdtvar, cdtvar)
call field_set_key_double(ivarfl(iw), kcdtvar, cdtvar)
call field_set_key_double(ivarfl(ipr), kcdtvar, cdtvar)

if (itytur.eq.2) then
  call field_get_key_double(ivarfl(ik), kcdtvar, cdtvar)
  call field_set_key_double(ivarfl(iep), kcdtvar, cdtvar)

elseif (itytur.eq.3) then
  call field_get_key_double(ivarfl(ir11), kcdtvar, cdtvar)
  call field_set_key_double(ivarfl(ir22), kcdtvar, cdtvar)
  call field_set_key_double(ivarfl(ir33), kcdtvar, cdtvar)
  call field_set_key_double(ivarfl(ir12), kcdtvar, cdtvar)
  call field_set_key_double(ivarfl(ir13), kcdtvar, cdtvar)
  call field_set_key_double(ivarfl(ir23), kcdtvar, cdtvar)
  call field_set_key_double(ivarfl(iep), kcdtvar, cdtvar)

  ! cdtvar for ial is useless because no time dependance in the equation of alpha.
  if (iturb.eq.32) then
    call field_get_key_double(ivarfl(ir11), kcdtvar, cdtvar)
    call field_set_key_double(ivarfl(ial), kcdtvar, cdtvar)
  endif
elseif (itytur.eq.5) then
  call field_get_key_double(ivarfl(ik), kcdtvar, cdtvar)
  call field_set_key_double(ivarfl(iep), kcdtvar, cdtvar)
  call field_set_key_double(ivarfl(iphi), kcdtvar, cdtvar)

!     CDTVAR pour IFB/IAL est en fait inutile car pas de temps dans
!     l'eq de f_barre/alpha
  if (iturb.eq.50) then
    call field_get_key_double(ivarfl(ik), kcdtvar, cdtvar)
    call field_set_key_double(ivarfl(ifb), kcdtvar, cdtvar)
  elseif (iturb.eq.51) then
    call field_get_key_double(ivarfl(ik), kcdtvar, cdtvar)
    call field_set_key_double(ivarfl(ial), kcdtvar, cdtvar)
  endif
elseif (iturb.eq.60) then
  call field_get_key_double(ivarfl(ik), kcdtvar, cdtvar)
  call field_set_key_double(ivarfl(iomg), kcdtvar, cdtvar)
elseif (iturb.eq.70) then
  ! cdtvar est a 1.0 par defaut dans cs_parameters.c
  call field_get_key_double(ivarfl(inusa), kcdtvar, cdtvar)
  call field_set_key_double(ivarfl(inusa), kcdtvar, cdtvar)
endif

! ---> IWALLF
! For laminar cases or when using low Reynolds model: no wall function.
! When using mixing length, Spalart-Allmaras or LES: one scale log law.
! When using EB-RSM : all y+ wall functions
! In all other cases: 2 scales log law.
! Here iwallf is set automatically only if it wasn't set in the gui or
! in a user subroutine.

if (iwallf.eq.-999) then
  if (    iturb.eq.10.or.iturb.eq.70 &
      .or.itytur.eq.4) then
    iwallf = 2
  elseif (iturb.eq.0.or.itytur.eq.5) then
    iwallf = 0
  elseif (iturb.eq.32) then
    iwallf = 7
  else
    iwallf = 3
  endif
endif

! ---> IWALFS
! If the wall function for the velocity is the two scales wall function using
! Van Driest mixing length (iwallf=5), then the corresponding wall function for
! scalar should be used (iwalfs=1).
! For atmospheric Flows, it is by default Louis, or Monin-Obukhov
! Here iwalfs is set automatically only if it wasn't set in a user subroutine.

if (iwalfs.eq.-999) then
  if (ippmod(iatmos).ge.0) then
    iwalfs = 2
  else if (iwallf.eq.5) then
    iwalfs = 1
  else
    iwalfs = 0
  endif
endif

! ---> YPLULI
! 1/XKAPPA est la valeur qui assure la continuite de la derivee
! entre la zone lineaire et la zone logarithmique.

! Dans le cas des lois de paroi invariantes, on utilise la valeur de
! continuite du profil de vitesse, 10.88.

! Pour la LES, on remet 10.88, afin d'eviter des clic/clac quand on est a
! la limite (en modele a une echelle en effet, YPLULI=1/XKAPPA ne permet pas
! forcement de calculer u* de maniere totalement satisfaisante).
! Idem en Spalart-Allmaras.

if (ypluli.lt.-grand) then
  if (iwallf.eq.4 .or. itytur.eq.4 .or. iturb.eq.70.or.iwallf.eq.6.or.iturb.eq.60 &
      .or. iturb.eq.22 ) then
    ypluli = 10.88d0
  else
    ypluli = 1.d0/xkappa
  endif
endif

! ---> ICPSYR
!      Si l'utilisateur n'a pas modifie ICPSYR, on prend par defaut :
!        s'il n y a pas de couplage
!          0 pour tous les scalaires
!        sinon
!          1 pour le scalaire ISCALT s'il existe
!          0 pour les autres
!      Les modifs adequates devront etre ajoutees pour les physiques
!        particulieres
!      Les tests de coherence seront faits dans verini.

if (nscal.gt.0) then

!     On regarde s'il y a du couplage

  nbccou = cs_syr_coupling_n_couplings()

!     S'il y a du couplage
  if (nbccou .ne. 0) then

!       On compte le nombre de scalaires couples
    nscacp = 0
    do iscal = 1, nscal
      call field_get_key_int(ivarfl(isca(iscal)), kcpsyr, icpsyr)
      if (icpsyr.eq.1) then
        nscacp = nscacp + 1
      endif
    enddo

!       Si l'utilisateur n'a pas couple de scalaire,
    if (nscacp.eq.0) then

!         On couple le scalaire temperature de la phase
      if (iscalt.gt.0.and.iscalt.le.nscal) then
        icpsyr = 1
        call field_set_key_int(ivarfl(isca(iscalt)), kcpsyr, icpsyr)
        goto 100
      endif
 100        continue

    endif

  endif

endif

! ---> "is_temperature"
!      If the user has not modified "is_temperature", we take by default:
!        passive scalar of scalars other than iscalt
!         = 0 : passive, enthalpy, or energy
!         = 1 : temperature

if (nscal.gt.0) then
  do ii = 1, nscal
    call field_get_key_int(ivarfl(isca(ii)), kscacp, iscacp)
    if (iscacp.eq.-1) then
      if (ii.eq.iscalt .and. itherm.eq.1) then
        iscacp = 1
      else
        iscacp = 0
      endif
      call field_set_key_int(ivarfl(isca(ii)), kscacp, iscacp)
    endif
  enddo
endif

! ---> ICALHY
!      Calcul de la pression hydrostatique en sortie pour les conditions de
!        Dirichlet sur la pression. Se deduit de IPHYDR et de la valeur de
!        la gravite (test assez arbitraire sur la norme).
!      ICALHY est initialise a -1 (l'utilisateur peut avoir force
!        0 ou 1 et dans ce cas, on ne retouche pas)

if (icalhy.ne.-1.and.icalhy.ne.0.and.icalhy.ne.1) then
  write(nfecra,1061) icalhy
  iok = iok + 1
endif

! ---> IKECOU
!      If the fluid_solid option is enabled, we force ikecou to 0.
if (fluid_solid) then
  if(ikecou .eq. 1) then
    ikecou = 0
    write(nfecra,5000)
  endif
endif

! ---> RELAXV
if (idtvar.lt.0) then
  relxsp = 1.d0-relxst
  if (relxsp.le.epzero) relxsp = relxst
  call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
  if (abs(vcopt%relaxv+1.d0).le.epzero) then
    vcopt%relaxv = relxsp
    call field_set_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
  endif
  do f_id = 0, nfld - 1
    call field_get_type(f_id, f_type)
    ! Is the field of type FIELD_VARIABLE?
    if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
      call field_get_key_struct_var_cal_opt(f_id, vcopt)
      if (abs(vcopt%relaxv+1.d0).le.epzero) then
        vcopt%relaxv = relxst
        call field_set_key_struct_var_cal_opt(f_id, vcopt)
      endif
    endif
  enddo
else
  do f_id = 0, nfld - 1
    call field_get_type(f_id, f_type)
    ! Is the field of type FIELD_VARIABLE?
    if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
      call field_get_key_struct_var_cal_opt(f_id, vcopt)
      if (abs(vcopt%relaxv+1.d0).le.epzero) then
        vcopt%relaxv = 1.d0
        call field_set_key_struct_var_cal_opt(f_id, vcopt)
      endif
    endif
  enddo
endif

! Options specific to steady case
if (idtvar.lt.0) then
  ipucou = 0
  dtref = 1.d0
  dtmin = 1.d0
  dtmax = 1.d0
  do f_id = 0, nfld - 1
    call field_get_type(f_id, f_type)
    ! Is the field of type FIELD_VARIABLE?
    if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
      call field_get_key_struct_var_cal_opt(f_id, vcopt)
      vcopt%istat = 0
      call field_set_key_struct_var_cal_opt(f_id, vcopt)
    endif
  enddo
  call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)
  arak = arak/max(vcopt%relaxv,epzero)
endif

! With a staggered approach no Rhie and Chow correction is needed
if (staggered.eq.1) then
  arak = 0.d0
endif

!===============================================================================
! 3. TABLEAUX DE cstphy
!===============================================================================

! ---> ICLVFL
!      Si l'utilisateur n'a pas modifie ICLVFL, on prend par defaut :
!        0 pour les variances
!      Les modifs adequates devront etre ajoutees pour les physiques
!        particulieres
!      If the user gives a value we put iclcfl to 2.

do iscal = 1, nscal
  if (iscavr(iscal).gt.0) then
    call field_get_key_int(ivarfl(isca(iscal)), kclvfl, iclvfl)
    ! Get the min clipping
    call field_get_key_double(ivarfl(isca(iscal)), kscmin, scminp)
    ! If modified put 2
    if (iclvfl.eq.-1 .and. abs(scminp+grand).ge.epzero) then
      call field_set_key_int(ivarfl(isca(iscal)), kclvfl, 2)

    else if (iclvfl.eq.-1) then
      call field_set_key_int(ivarfl(isca(iscal)), kclvfl, 0)
    endif

    ! Min for variances is 0 or greater
    call field_get_key_double(ivarfl(isca(iscal)), kscmin, scminp)
    ! set min clipping to 0
    scminp = max(0.d0, scminp)
    call field_set_key_double(ivarfl(isca(iscal)), kscmin, scminp)
  endif
enddo

do ii = 1, nscal
  f_id = ivarfl(isca(ii))
  call field_get_key_double(f_id, kvisl0, visls_0)

  ! For scalars which are not variances, define the reference diffusivity
  if (iscavr(ii).le.0 .and. visls_0.lt.-grand) then
    call field_get_key_int(f_id, kscacp, iscacp)
    if (iscacp.gt.0) then
      ! For temperature, the diffusivity factor is directly the thermal conductivity
      ! lambda = Cp * mu / Pr
      ! where Pr is the (molecular) Prandtl number
      visls_0 = viscl0 * cp0
    else
      visls_0 = viscl0
    endif
    call field_set_key_double(f_id, kvisl0, visls_0)
  endif

  ! For fluctuation variances, the diffusivity is that of the associated scalar.
  iscal = iscavr(ii)
  if (iscal.gt.0.and.iscal.le.nscal)then
    call field_get_key_double(ivarfl(isca(iscal)), kvisl0, visls_0)
    call field_get_key_double(f_id, kvisl0, visls_cmp)
    call field_set_key_double(f_id, kvisl0, visls_0)
    if (visls_cmp.gt.-grand) then
      write(nfecra,1071) ii, iscal, ii, iscal, visls_0
    endif
  endif
enddo

! xyzp0 : reference point for hydrostatic pressure
! The user should specify the 3 coordinates, otherwise
! it is set to (0.,0.,0.).

if (xyzp0(1).gt.-0.5d0*rinfin.and. &
    xyzp0(2).gt.-0.5d0*rinfin.and. &
    xyzp0(3).gt.-0.5d0*rinfin       ) then
  ixyzp0 = 1
else
  do ii = 1, 3
    xyzp0(ii) = 0.d0
  enddo
endif

! VoF model enabled
if (ivofmt.gt.0) then
  if (rho2.gt.rho1) then
    ro0    = rho2
    viscl0 = mu2
  else
    ro0    = rho1
    viscl0 = mu1
  endif

  ! VOF algorithm: continuity of the flux across internal faces
  call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
  vcopt%imvisf = 1
  call field_set_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
endif

!===============================================================================
! 4. ELEMENTS DE albase
!===============================================================================

if (iale.ge.1) then
  if (isuite.eq.0 .and. italin.eq.-999 ) italin = 1
else
  italin = 0
endif

!===============================================================================
! 4. PARAMETRES DE cplsat
!===============================================================================

! Get number of couplings

call nbccpl(nbrcpl)

if (nbrcpl.ge.1.and.iturbo.ne.0) then
  ifaccp = 1
endif

!===============================================================================
! 5. Define Min/Max clipping values of void fraction of VOF model
!===============================================================================

if (ivofmt.gt.0) then
  call field_get_key_double(ivarfl(ivolf2), kscmin, clvfmn)
  call field_get_key_double(ivarfl(ivolf2), kscmax, clvfmx)

  if (clvfmn.lt.-0.5d0*grand) then
    clvfmn = 0.d0
    if (iand(ivofmt,VOF_MERKLE_MASS_TRANSFER).ne.0) clvfmn = epzero
  endif
  if (clvfmx.gt.0.5d0*grand) then
    clvfmx = 1.d0
    if (iand(ivofmt,VOF_MERKLE_MASS_TRANSFER).ne.0) clvfmx = 1.d0-epzero
  endif

  call field_set_key_double(ivarfl(ivolf2), kscmin, clvfmn)
  call field_set_key_double(ivarfl(ivolf2), kscmax, clvfmx)
endif

!===============================================================================
! 6. STOP SI PB
!===============================================================================

if (iok.ne.0) then
  call csexit(1)
endif

 1011 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION',                /,&
'@    ========',                                                /,&
'@    ',a6,' = ',   i10,                                        /,&
'@    ',a6,' WILL BE INITIALIZED AUTOMATICALLY',                /,&
'@    DO NOT MODIFY IT.,'                                       /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  Check cs_user_parameters.f90',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 1021 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION',                /,&
'@    ========',                                                /,&
'@    SCALAR ',   i10,' ',a6,' = ',   i10,                      /,&
'@    ',a6,' WILL BE INITIALIZED AUTOMATICALLY',                /,&
'@    DO NOT MODIFY IT.,'                                       /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  Check cs_user_parameters.f90',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 1131 format( &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ADVANCED MODIFICATION FOR',                      /,&
'@    ========,'                                                /,&
'@    ',a17,' OF THE VARIABLE'                                  /,&
'@    ',a6,'.'                                                  /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 1061 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION',                /,&
'@    ========',                                                /,&
'@    ICALHY must be an integer equal to 0 or 1',               /,&
'@',                                                            /,&
'@  Its value is ',i10,                                         /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  Check cs_user_parameters.f90',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 1071 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING:       IN THE DATA SPECIFICATION',                /,&
'@    ========',                                                /,&
'@',                                                            /,&
'@  The scalar ',i10,   ' is the fluctuations variance',        /,&
'@    of the scalar ',i10,                                      /,&
'@',                                                            /,&
'@  The diffusivity_ref value of the scalar ', i10,             /,&
'@    must not be set:',                                        /,&
'@    it is automatically set equal to the scalar',             /,&
'@    diffusivity ', i10,   ' i.e. ',e14.5,                     /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3000 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING:       IN THE DATA SPECIFICATION',                /,&
'@    ========',                                                /,&
'@',                                                            /,&
'@  The cavitation model requires an upwind convection scheme' ,/,&
'@    for the void fraction (BLENCV(IVOLF2)=',e14.5,').',       /,&
'@  The user has set BLENCV(IVOLF2)=',e14.5,                    /,&
'@',                                                            /,&
'@  The upwind scheme for the void fraction is forced.',        /,&
'@',                                                            /,&
'@  The calculation will be run.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 5000 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING:       IN THE DATA SPECIFICATION',                /,&
'@    ========',                                                /,&
'@',                                                            /,&
'@  The pseudo coupling of turbulent dissipation and turbulent',/,&
'@  kinetic energy (ikecou = 1) is not compatible with the use',/,&
'@  of fluid/solid option to disable the dynamic in the solid ',/,&
'@  cells (fluid_solid =1). ',                                  /,&
'@',                                                            /,&
'@  The parameter ikecou is forced to 0 (no coupling)',         /,&
'@',                                                            /,&
'@  The calculation will be run.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
return
end subroutine
