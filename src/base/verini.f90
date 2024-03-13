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

subroutine verini &
 ( iok    )

!===============================================================================
! Purpose:
! --------

! Check computation parameters after user modification (modules)
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use albase
use parall
use period
use ppthch
use ppppar
use ppincl
use lagran
use radiat
use cplsat
use mesh
use field
use turbomachinery
use cs_c_bindings
use cfpoin

!===============================================================================

implicit none

! Arguments

integer          iok

! Local variables

character        chaine*80
integer          ii    , iis   , jj    , iisct
integer          iscal , ivar
integer          f_id, n_fields
integer          iiidef, istop, iclvfl
integer          kscmin, kscmax, kclvfl
integer          keyvar, keysca
integer          key_t_ext_id, icpext
integer          iviext, iscacp, kcdtvar
integer          iroext, isso2t
integer          ivisext, krvarfl, kisso2t
integer          kturt, turb_flux_model, kthetss, kthetvs
double precision scmaxp, rvarfl, thetss, thetvs, cdtvar
double precision turb_schmidt, visls_0

character(len=3), dimension(3) :: nomext3
character(len=4), dimension(3) :: nomext63

type(var_cal_opt) :: vcopt, vcopt1, vcopt2, vcopt3

!===============================================================================

! Initialize variables to avoid warnings

jj = 0

nomext3 = (/'[X]', '[Y]', '[Z]'/)
nomext63 = (/'[11]', '[22]', '[33]'/)

call field_get_n_fields(n_fields)

call field_get_key_id("scalar_id", keysca)
call field_get_key_id("variable_id", keyvar)

! Key ids for clippings
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)
call field_get_key_id("variance_clipping", kclvfl)

call field_get_key_id("st_exp_extrapolated", kthetss)
call field_get_key_id("diffusivity_extrapolated", kthetvs)

! Time extrapolation?
call field_get_key_id("time_extrapolated", key_t_ext_id)

call field_get_key_id("scalar_time_scheme", kisso2t)

call field_get_key_id("variance_dissipation", krvarfl)
call field_get_key_id("time_step_factor", kcdtvar)

if (icp.ge.0) then
  call field_get_key_int(icp, key_t_ext_id, icpext)
else
  icpext = 0
endif
call field_get_key_int(iviscl, key_t_ext_id, iviext)
call field_get_key_int(icrom, key_t_ext_id, iroext)

!===============================================================================
! 1. OPTIONS DU CALCUL : TABLEAUX DE optcal : formats 2000
!===============================================================================

! --- Dimensions

if (nscal.lt.0.or.nscal.gt.nscamx) then
  write(nfecra,2000)'NSCAL ',nscamx,nscal
  iok = iok + 1
endif
if (nscaus.lt.0.or.nscaus.gt.nscamx) then
  write(nfecra,2000)'NSCAUS',nscamx,nscaus
  iok = iok + 1
endif
if (nscapp.lt.0.or.nscapp.gt.nscamx) then
  write(nfecra,2000)'NSCAPP',nscamx,nscapp
  iok = iok + 1
endif
if (nvar.lt.0.or.nvar.gt.nvarmx) then
  write(nfecra,2000)'NVAR  ',nvarmx,nvar
  iok = iok + 1
endif

if (itytur.eq.4.or.ischtp.eq.2) then
  call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)
  call field_get_label(ivarfl(iu), chaine)
  chaine = trim(chaine)
  iiidef = 10
  if (vcopt%nswrsm.lt.iiidef) then
    write(nfecra,2125) chaine(1:16),iiidef,vcopt%nswrsm
  endif
  jj = ipr
  call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
  call field_get_label(ivarfl(jj), chaine)
  iiidef = 5
  if (vcopt%nswrsm.lt.iiidef) then
    write(nfecra,2125) chaine(1:16),iiidef,vcopt%nswrsm
  endif

  ! Scalars
  do ii = 1, nscal
    iiidef = 10
    call field_get_label(ivarfl(isca(ii)), chaine)
    call field_get_key_struct_var_cal_opt(ivarfl(isca(ii)), vcopt)
    chaine = trim(chaine)
    if (vcopt%nswrsm.lt.iiidef) then
      write(nfecra,2125) chaine(1:16),iiidef,vcopt%nswrsm
    endif
  enddo
endif

!     Test du theta de la viscosite secondaire, du flux de masse et
!     de la viscosite par rapport a celui de la vitesse
jj = iu
call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)
if (abs(vcopt%thetav-1.d0).lt.epzero.and.               &
     (istmpf.eq.2.or.                                   &
     isno2t.ne.0.or.                                    &
     isto2t.ne.0.or.                                    &
     iroext.ne.0.or.                                    &
     iviext.ne.0.or.                                    &
     icpext.ne.0   ) ) then
  write(nfecra,2131) vcopt%thetav,                      &
       istmpf,isno2t,isto2t,                            &
       iroext,iviext,icpext
endif
if (abs(vcopt%thetav-0.5d0).lt.epzero.and.             &
     (istmpf.ne.2.or.                                   &
     isno2t.ne.1.or.                                    &
     isto2t.ne.1.or.                                    &
     iroext.ne.1.or.                                    &
     iviext.ne.1.or.                                    &
     icpext.ne.1   ) ) then
  write(nfecra,2132) vcopt%thetav,                      &
       istmpf,isno2t,isto2t,                            &
       iroext,iviext,icpext
endif
do iscal = 1, nscal
  call field_get_key_int(ivarfl(isca(iscal)), kisso2t, isso2t)
  if (isso2t.ne.isno2t) then
    write(nfecra,2133) isno2t
  endif
  call field_get_key_int (ivarfl(isca(iscal)), kivisl, f_id)
  if (f_id.ge.0.and.iscavr(iscal).le.0) then
    call field_get_key_int(f_id, key_t_ext_id, ivisext)
    if (ivisext.ne.iviext) then
      write(nfecra,2134) iscal,ivisext,iviext
    endif
  endif
enddo

if (ischtp.eq.2.and.vcopt%ibdtso.gt.1) then
  ! NB: this test does not prevent from incompatible user modifications of
  ! isno2t, thetav, etc.
  write(nfecra,1135)
  iok = iok + 1
endif

!     A priori, pour le moment, l'ordre 2 en temps
!       n'est pas pris en compte en k-eps, v2f ou k-omega couple : on s'arrete
if (itytur.eq.2 .and.ikecou.eq.1) then
  call field_get_key_struct_var_cal_opt(ivarfl(ik), vcopt)
  call field_get_key_struct_var_cal_opt(ivarfl(iep), vcopt1)
  if ((     thetst       .gt.0.d0 ).or.                    &
       (    isto2t       .gt.0    ).or.                    &
       (abs(vcopt%thetav-1.0d0).gt.epzero).or.             &
       (abs(vcopt1%thetav-1.0d0).gt.epzero) ) then
    write(nfecra,2142)iturb,ikecou,                        &
         thetst,isto2t,                                    &
         vcopt%thetav,vcopt1%thetav
    iok = iok + 1
  endif
endif
if (iturb.eq.50.and.ikecou.eq.1) then
  call field_get_key_struct_var_cal_opt(ivarfl(ik), vcopt)
  call field_get_key_struct_var_cal_opt(ivarfl(iep), vcopt1)
  call field_get_key_struct_var_cal_opt(ivarfl(iphi), vcopt2)
  call field_get_key_struct_var_cal_opt(ivarfl(ifb), vcopt3)
  if ((     thetst       .gt.0.d0 ).or.                    &
      (    isto2t       .gt.0    ).or.                     &
      (abs(vcopt%thetav-1.0d0).gt.epzero).or.              &
      (abs(vcopt1%thetav-1.0d0).gt.epzero).or.             &
      (abs(vcopt2%thetav-1.0d0).gt.epzero).or.             &
      (abs(vcopt3%thetav-1.0d0).gt.epzero) ) then
    write(nfecra,2143)iturb,ikecou,                        &
         thetst,isto2t,                                    &
         vcopt%thetav,vcopt1%thetav,                       &
         vcopt2%thetav,vcopt3%thetav
    iok = iok + 1
  endif
endif
if (iturb.eq.51.and.ikecou.eq.1) then
  call field_get_key_struct_var_cal_opt(ivarfl(ik), vcopt)
  call field_get_key_struct_var_cal_opt(ivarfl(iep), vcopt1)
  call field_get_key_struct_var_cal_opt(ivarfl(iphi), vcopt2)
  call field_get_key_struct_var_cal_opt(ivarfl(ial), vcopt3)
  if ((    thetst       .gt.0.d0 ).or.                     &
     (    isto2t       .gt.0    ).or.                      &
      (abs(vcopt%thetav-1.0d0).gt.epzero).or.              &
      (abs(vcopt1%thetav-1.0d0).gt.epzero).or.             &
      (abs(vcopt2%thetav-1.0d0).gt.epzero).or.             &
      (abs(vcopt3%thetav-1.0d0).gt.epzero) ) then
    write(nfecra,2143)iturb,ikecou,                        &
         thetst,isto2t,                                    &
         vcopt%thetav,vcopt1%thetav,                       &
         vcopt2%thetav,vcopt3%thetav
    iok = iok + 1
  endif
endif
if (iturb.eq.60.and.ikecou.eq.1) then
  call field_get_key_struct_var_cal_opt(ivarfl(ik), vcopt)
  call field_get_key_struct_var_cal_opt(ivarfl(iomg), vcopt1)
  if ((    thetst       .gt.0.d0 ).or.                     &
       (    isto2t       .gt.0    ).or.                    &
       (abs(vcopt%thetav-1.0d0).gt.epzero).or.             &
       (abs(vcopt1%thetav-1.0d0).gt.epzero) ) then
    write(nfecra,2144)iturb,ikecou,                        &
         thetst,isto2t,                                    &
         vcopt%thetav,vcopt1%thetav
    iok = iok + 1
  endif
endif
if (iturb.eq.70) then
  call field_get_key_struct_var_cal_opt(ivarfl(inusa), vcopt)
  if ((thetst .gt.0.d0).or.                                &
      (isto2t .gt. 0).or.                                  &
      (abs(vcopt%thetav-1.0d0).gt.epzero)) then
    write(nfecra,2145)iturb,thetst,isto2t,vcopt%thetav
    iok = iok + 1
  endif
endif

!     A priori, pour le moment, l'ordre 2 en temps
!       (rho, visc, cp, termes sources N.S., Turb., Scal., theta)
!     est incompatible avec les physiques particulieres
!     Ici on s'arrete si on n'est pas dans le cas du schema std

if (ippmod(iphpar).ge.1) then
  istop = 0
  do ivar = 1, nvar
    call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
    if ( (abs(vcopt%thetav-1.0d0).gt.1.d-3) ) istop = 1
  enddo
  if ((thetsn .gt.0.d0).or.                                &
      (isno2t .gt.0   ).or.                                &
      (iroext .gt.0   ).or.                                &
      (thetvi .gt.0.d0).or.                                &
      (iviext .gt.0   ).or.                                &
      (thetcp .gt.0.d0).or.                                &
      (icpext .gt.0   )) istop = 1
  do iscal = 1, nscal
    call field_get_key_int(ivarfl(isca(iscal)), kisso2t, isso2t)
    call field_get_key_double(ivarfl(isca(iscal)), kthetss, thetss)
    call field_get_key_double(ivarfl(isca(iscal)), kthetvs, thetvs)
    if ((thetss.gt.0.d0).or.            &
        (isso2t.gt.0   ).or.            &
        (thetvs.gt.0.d0)) istop = 1
  enddo

  if (istop.ne.0) then
    write(nfecra,2146)
  endif
endif

!     A priori, pour le moment, l'ordre 2 en temps
!       n'est pas pris en compte pour les termes issus du Lagrangien.
!       On pourrait le signaler et continuer : on s'arrete.
if (iilagr .eq. 2) then
  if ((    thetsn       .gt.0.d0 ).or.                    &
       (    isno2t       .gt.0    ).or.                   &
       (    thetst       .gt.0.d0 ).or.                   &
       (    isto2t       .gt.0    ) ) then
    write(nfecra,2147)thetsn,isno2t,thetst,isto2t
    iok = iok + 1
  endif
  if ((itherm.eq.1 .and. itpscl.eq.1) .or. itherm.eq.2) then
    call field_get_key_int(ivarfl(isca(iscalt)), kisso2t, isso2t)
    call field_get_key_double(ivarfl(isca(iscalt)), kthetss, thetss)
    if (thetss.gt.0.d0 .or. isso2t.gt.0) then
      write(nfecra,2148)                                           &
        'lagrangian ',ivarfl(isca(iscal)),thetss,isso2t, 'cs_user_lagr_model'
      iok = iok + 1
    endif
  endif
endif

!     A priori, pour le moment, l'ordre 2 en temps
!       n'est pas pris en compte pour les termes issus du rayonnement.
!       On pourrait le signaler et continuer : on s'arrete.
if (iirayo.gt.0) then
  if (iscalt.gt.0) then
    call field_get_key_int(ivarfl(isca(iscalt)), kisso2t, isso2t)
    call field_get_key_double(ivarfl(isca(iscalt)), kthetss, thetss)
    if (thetss.gt.0.d0 .or. isso2t.gt.0) then
      write(nfecra,2148) 'rayonnement',iscal,thetss,'usray1'
      iok = iok + 1
    endif
  endif
endif

! Time step multiplier

do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyvar, ii)
  if (ii.ge.1) then
    call field_get_key_double(ivarfl(ii), kcdtvar, cdtvar)
    if (cdtvar.le.0.d0) then
      call field_get_label(f_id, chaine)
      write(nfecra,2530) chaine(1:16),ivarfl(ii),cdtvar
      iok = iok + 1
    endif
  endif
enddo

! --- Turbulence

!    Nb de variables

if (nscal.ge.1) then
  if (iscalt.gt.nscal) then
    write(nfecra,2610)                                   &
         'NUMERO DU SCALAIRE TEMPERATURE ',iscalt,       &
         'NOMBRE DE SCALAIRES  ',          nscal
    iok = iok + 1
  endif
  if ( (nvar.lt. 4+nscal               ) .or.            &
       (nvar.lt. 6+nscal.and.itytur.eq.2).or.            &
       (nvar.lt.11+nscal.and.iturb.eq.30).or.            &
       (nvar.lt.11+nscal.and.iturb.eq.31).or.            &
       (nvar.lt.12+nscal.and.iturb.eq.32).or.            &
       (nvar.lt. 8+nscal.and.itytur.eq.5).or.            &
       (nvar.lt. 6+nscal.and.iturb.eq.60).or.            &
       (nvar.lt. 5+nscal.and.iturb.eq.70)      ) then
    write(nfecra,2610)                                   &
         'NOMBRE DE VARIABLES  ',          nvar,         &
         'NOMBRE DE SCALAIRES  ',          nscal
    iok = iok + 1
  endif


  do iscal = 1, nscal
    call field_get_key_id('turbulent_flux_model', kturt)
    call field_get_key_int(ivarfl(isca(iscal)), kturt, turb_flux_model)

    ! Turbulent flux model for scalar
    if (turb_flux_model.ne. 0.and.turb_flux_model.ne.10 .and. &
        turb_flux_model.ne.20.and.turb_flux_model.ne.30 .and. &
        turb_flux_model.ne.11.and.turb_flux_model.ne.21 .and. &
        turb_flux_model.ne.31                              &
                     ) then
      write(nfecra,2604) 'turbulent_flux_model  ', turb_flux_model
      write(nfecra,2610)                         &
           'Index of the scalar: ', iscal,       &
           'Number of scalars: ', nscal

      iok = iok + 1
    endif
  enddo

endif

!      Specifique k-epsilon, v2f et k-omega

if (itytur.eq.2 .or. iturb.eq.50                          &
     .or. iturb.eq.60 ) then
  if ( (nvar.le.5.and.itytur.eq.2) .or.                   &
       (nvar.le.7.and.itytur.eq.5) .or.                   &
       (nvar.le.5.and.iturb.eq.60)     ) then
    write(nfecra,2610)                                    &
         'NOMBRE DE VARIABLES  ',          nvar,          &
         'OPTION POUR LA TURBULENCE      ',iturb
    iok = iok + 1
  endif

  !        IF ( IGRAKE.EQ.1.AND.(GX**2+GY**2+GZ**2).LE.EPZERO**2 ) THEN
  !          WRITE(NFECRA,2620)'IGRAKE',IGRAKE,GX,GY,GZ
  !          IOK = IOK + 1
  !        ENDIF
  if (nscal.gt.0) then
    if (iscalt.le.0.and. (gx**2+gy**2+gz**2).ge.epzero**2) then
      write(nfecra,2621) gx,gy,gz,iscalt
      if (igrake.eq.1) then
        write(nfecra,2622)'IGRAKE',igrake
      endif
    endif
  endif
endif

!     Specifique Rij-epsilon

if (itytur.eq.3) then
  if (nvar.le.10) then
    write(nfecra,2610)                                    &
         'NOMBRE DE VARIABLES            ', nvar,         &
         'OPTION POUR LA TURBULENCE      ', iturb
    iok = iok + 1
  endif
  !        IF ( IGRARI.EQ.1.AND.(GX**2+GY**2+GZ**2).LE.EPZERO**2 ) THEN
  !          WRITE(NFECRA,2620)'IGRARI',IGRARI,GX,GY,GZ
  !          IOK = IOK + 1
  !        ENDIF
  if (nscal.gt.0) then
    if (iscalt.le.0.and.(gx**2+gy**2+gz**2).ge.epzero**2) then
      write(nfecra,2621)gx,gy,gz,iscalt
      if (igrari.eq.1) then
        write(nfecra,2622)'IGRARI',igrari
      endif
    endif
  endif
endif

!===============================================================================
! 2. TABLEAUX DE cstphy : formats 4000
!===============================================================================

! --- Scalaires

if (nscal.gt.0) then

!     Scalaire passif, temperature, enthalpie, energie
  do ii = 1, nscal
    call field_get_key_int(ivarfl(isca(ii)), kscacp, iscacp)
    if (iscacp.lt.0.or.iscacp.gt.2) then
      call field_get_label(ivarfl(isca(ii)), chaine)
      write(nfecra,4300)chaine(1:16),iscacp
      iok = iok + 1
    endif
  enddo

!     Valeur du Schmidt turbulent positif
  do ii = 1, nscal
    call field_get_key_double(ivarfl(isca(ii)), ksigmas, turb_schmidt)
    if (turb_schmidt.le.0d0) then
      call field_get_label(ivarfl(isca(ii)), chaine)
      write(nfecra,4350)chaine(1:16),ii,turb_schmidt
      iok = iok + 1
    endif
  enddo

!     Valeur de la borne sup de clipping si on l'utilise
  do ii = 1, nscal
    ! Get the max clipping
    call field_get_key_double(ivarfl(isca(ii)), kscmax, scmaxp)
    call field_get_key_int(isca(ii), kclvfl, iclvfl)
    if (iscavr(ii).gt.0.and.iscavr(ii).le.nscal.and.               &
       iclvfl.eq.2.and.scmaxp.le.0.d0) then
      call field_get_label(ivarfl(isca(ii)), chaine)
      write(nfecra,4370)chaine(1:16),ii,scmaxp
      iok = iok + 1
    endif
  enddo

!     Warning sur Rvarfl < 0
  do ii = 1, nscal
    call field_get_key_double(ivarfl(isca(ii)), krvarfl, rvarfl)

    if (iscavr(ii).gt.0.and.iscavr(ii).le.nscal.and.rvarfl.le.0.d0) then
      call field_get_label(ivarfl(isca(ii)), chaine)
      write(nfecra,4380)chaine(1:16),ii,rvarfl
      iok = iok + 1
    endif
  enddo


!  Rien a verifier sur SCAMIN SCAMAX

!     Si CP0 est utilise (resolution d'un scalaire en temperature
!       et CP constant), il doit etre positif
  if (icp.eq.-1) then
    if (cp0.lt.0.d0) then
      iisct = 0
      do iis = 1, nscal
        call field_get_key_int(ivarfl(isca(iis)), kscacp, iscacp)
        if (iscacp.gt.0) then
          iisct = 1
        endif
      enddo
      if (iisct.eq.1) then
        write(nfecra,2511)'CP0   ',cp0
        iok = iok + 1
      endif
    endif
  endif

endif

!===============================================================================
! 3. TABLEAUX DE period : formats 5000
!===============================================================================

! --- periodicite de rotation incompatible avec couplage
!       renforce vitesse pression et l'ALE
if (iperot.gt.0.and. (ipucou.ne.0.or.iale.ne.0)) then
  write(nfecra,5003)iperio,ipucou,iale
  iok = iok + 1
endif

! --- periodicite incompatible avec le mode de calcul
!       direct de la distance a la paroi
if (iperio.eq.1.and.ineedy.eq.1.and.icdpar.eq.2) then
  write(nfecra,5005)iperio, icdpar
  iok = iok + 1
endif

! --- periodicite de rotation incompatible avec le rayonnement DOM
if (iperot.gt.0.and.iirayo.gt.0) then
  if (iirayo.eq.1) then
    write(nfecra,5008) iperio,  iirayo
    iok = iok + 1
  endif
endif

!===============================================================================
! 4. TABLEAUX DE parall : formats 6000 (limitations)
!===============================================================================

! --- parallelisme incompatible avec le mode de calcul
!       direct de la distance a la paroi
if (irangp.ge.0.and.ineedy.eq.1.and.icdpar.eq.2) then
  write(nfecra,6005)irangp, icdpar
  iok = iok + 1
endif

!===============================================================================
! 5. METHODE ALE (albase)
!===============================================================================

if (iale.ge.1) then

  if (nalinf.lt.0) then
    write(nfecra,7010)nalinf
    iok = iok + 1
  endif

  if (nalimx.le.0) then
    write(nfecra,7030)nalimx
    iok = iok + 1
  endif
  if (epalim.le.0.d0) then
    write(nfecra,7040)epalim
    iok = iok + 1
  endif

  if (italin.ne.-999 .and. italin.ne.0 .and. italin.ne.1) then
    write(nfecra,7050)italin
    iok = iok + 1
  endif

endif

!===============================================================================
! 6. COMPRESSIBLE : formats 8000
!===============================================================================

if (ippmod(icompf).ge.0) then
  if (t0.le.0.d0.or.p0.le.0.d0) then
    write(nfecra,8000)t0,p0
    iok = iok + 1
  endif
  call field_get_key_double(ivarfl(isca(itempk)), kvisl0, visls_0)
  if (visls_0.le.0.d0) then
    write(nfecra,8010) visls_0
    iok = iok + 1
  endif
  if (viscv0.lt.0.d0) then
    write(nfecra,8020) viscv0
    iok = iok + 1
  endif
  if (ieos.lt.1.or.ieos.gt.4) then
    write(nfecra,8030) 'IEOS (Equation of state. )',ieos
    iok = iok + 1
  endif
  if (ieos.eq.2.and.gammasg.lt.1.d0) then
    write(nfecra,8040) gammasg
    iok = iok + 1
  endif
  if (ieos.eq.1.and.cp0.lt.cv0) then
    write(nfecra,8050) cp0, cv0
    iok = iok + 1
  endif
  if ((ieos.eq.1.or.ieos.eq.3).and.psginf.ne.0.d0) then
    write(nfecra,8060) psginf
    iok = iok + 1
  endif
endif

!===============================================================================
! 7. Unsteady rotor/stator coupling: 9000 formats
!===============================================================================

if (iturbo.eq.2) then
  ! Unsteady rotor/stator coupling is not compatible with the
  !   steady algorithm...
  if (idtvar.lt.0) then
    write(nfecra,9010) idtvar
    iok = iok + 1
  endif
  ! ... nor with the time/space variable time steps
  ! TODO-> make it ok for idtvar=1
  if (idtvar.eq.1.or.idtvar.eq.2) then
    write(nfecra,9011) idtvar
    iok = iok + 1
  endif
endif

!===============================================================================
! 8. Cavitation modelling: 9100 formats
!===============================================================================

if (ivofmt.gt.0) then
  ! VOF method is not compatible with dilatable or low-mach algorithms
  if (idilat.gt.1) then
    write(nfecra,9120) idilat
    iok = iok + 1
  endif
endif

!===============================================================================
! 9. FORMATS VERIFICATION
!===============================================================================

 2000 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  ERROR:   STOP WHILE READING INPUT DATA',                 /,&
'@     =====',                                                  /,&
'@',  a6, ' must be an integer',                                /,&
'@    larger than 0 or less than or equal to ', i10,            /,&
'@    it has value', i10,                                       /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2125 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING:      WHEN READING INPUT DATA',                   /,&
'@    =======',                                                 /,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME',                    /,&
'@',                                                            /,&
'@   ORDER 2 IN TIME or LES',                                   /,&
'@   THE VALUE RECOMMENDED FOR THE PARAMETER NSWRSM FOR',       /,&
'@        VARIABLE ', a16, ' IS',    i10,                       /,&
'@     NSWRSM IS NOW IMPOSED AS',  i10,                         /,&
'@',                                                            /,&
'@  computation will go on',                                    /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 1135 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:       WHEN READING INPUT DATA',                    /,&
'@    =====',                                                   /,&
'@   CHOICE OF TIME-SCHEME ISCHTP = 2 IS NOT COMPATIBLE WITH',  /,&
'@   IBDTSO > 1',                                               /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2131 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHEN READING INPUT DATA',                /,&
'@    =========',                                               /,&
'@   CHOICE OF TIME-SCHEME',                                    /,&
'@',                                                            /,&
'@     TIME-SCHEME FOR VELOCITY IS FIRST ORDER',                /,&
'@       (THETAV = ', e10.2, ')',                               /,&
'@     SOME TERMS ARE HOWEVER SECOND ORDER IN TIME WITH',       /,&
'@       THE FOLLOWING SETTINGS:',                              /,&
'@',                                                            /,&
'@ parameters       ISTMPF ISNO2T ISTO2T IROEXT IVIEXT ICPEXT', /,&
'@ Values  entered', 6I7,                                       /,&
'@',                                                            /,&
'@  computation will go on.',                                   /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2132 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHEN READING INPUT DATA',                /,&
'@    =========',                                               /,&
'@   CHOICE OF TIME-SCHEME',                                    /,&
'@',                                                            /,&
'@     TIME-SCHEME FOR VELOCITY IS SECOND ORDER',               /,&
'@       (THETAV = ', e10.2, ')',                               /,&
'@     SOME TERMS ARE HOWEVER FIRST ORDER IN TIME  WITH',       /,&
'@       THE FOLLOWING SETTINGS:',                              /,&
'@',                                                            /,&
'@ parameters       ISTMPF ISNO2T ISTO2T IROEXT IVIEXT ICPEXT', /,&
'@ Values  entered', 6I7,                                       /,&
'@',                                                            /,&
'@  computation will go on.',                                   /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2133 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHEN READING INPUT DATA',                /,&
'@    =========',                                               /,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME',                    /,&
'@',                                                            /,&
'@   SCALAR ',   i10,' ISSO2T = ', i10,                         /,&
'@       IS DIFFERENT FROM ISNO2T',                             /,&
'@     ISNO2T = ', i10,                                         /,&
'@',                                                            /,&
'@  computation will go on',                                    /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2134 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHEN READING INPUT DATA',                /,&
'@    =========',                                               /,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME',                    /,&
'@',                                                            /,&
'@   SCALAIRE ', i10,' IVSEXT = ', i10,                         /,&
'@       IS DIFFERENT FROM IVIEXT',                             /,&
'@     IVIEXT = ', i10,                                         /,&
'@',                                                            /,&
'@  computation will go on',                                    /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2142 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME',   /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  With the k-epsilon turbulence',                             /,&
'@  model (ITURB = ', i10, ') solved coupled(IKECOU = ', i10,'):',/,&
'@    the current version does not allow second order',         /,&
'@    in time resolution of k-epsilon equations in a coupled',  /,&
'@    manner.',                                                 /,&
'@',                                                            /,&
'@   Thus one or more of the values below are not permited',    /,&
'@',                                                            /,&
'@       THETST    ISTO2T     THETA K   THETA EPS',             /,&
'@',      e12.4,      i10,      e12.4,      e12.4,              /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2143 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME',   /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  With the   V2F     turbulence',                             /,&
'@  model (ITURB = ', i10, ') solved coupled(IKECOU = ', i10,'):',/,&
'@    the current version does not allow second order',         /,&
'@    in time resolution of V2F equations in a coupled',        /,&
'@    manner.',                                                 /,&
'@',                                                            /,&
'@   Thus one or more of the values below are not permited',    /,&
'@',                                                            /,&
'@       THETST    ISTO2T     THETA K   THETA EPS',             /,&
'@',      e12.4,      i10,      e12.4,      e12.4,              /,&
'@     THETA PHI    THETA FB',                                  /,&
'@',       e12.4,      e12.4,                                   /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2144 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME',   /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  With the k-omega   turbulence',                             /,&
'@  model (ITURB = ', i10, ') solved coupled(IKECOU = ', i10,'):',/,&
'@    the current version does not allow second order',         /,&
'@    in time resolution of k-omega equations in a coupled',    /,&
'@    manner.',                                                 /,&
'@',                                                            /,&
'@   Thus one or more of the values below are not permited',    /,&
'@',                                                            /,&
'@       THETST    ISTO2T     THETA K   THETA OMEGA',           /,&
'@',      e12.4,      i10,      e12.4,      e12.4,              /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2145 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME',   /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  With the Spalart-Allmaras',                                 /,&
'@  turbulence model (ITURB = ', i10, ')',                      /,&
'@    the current version does not allow second order',         /,&
'@    in time resolution.',                                     /,&
'@',                                                            /,&
'@   Thus one or more of the values below are not permited',    /,&
'@',                                                            /,&
'@       THETST    ISTO2T     THETA K   THETA OMEGA',           /,&
'@',      e12.4,      i10,      e12.4,      e12.4,              /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2146 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:',                                               /,&
'@    =========',                                               /,&
'@    OPTIONS NOT VALIDATED WITH TIME DICRETISATION SCHEME',    /,&
'@',                                                            /,&
'@  The current version is not validated with this time',       /,&
'@  discretisation scheme when a specific physics is active',   /,&
'@    (combustion, coal, electrical, ...).',                    /,&
'@',                                                            /,&
'@  Verify   the parameters.',                                  /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2147 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME',   /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  The source terms stemming from the',                        /,&
'@    Lagrangien module will not be computed as second order',  /,&
'@    in this version, despite user settings chosen below',     /,&
'@',                                                            /,&
'@',                                                            /,&
'@     For Navier-Stokes       For   turbulence',               /,&
'@       THETSN    ISNO2T      THETST    ISTO2T',               /,&
'@',      e12.4,      i10,      e12.4,      i10,                /,&
'@',                                                            /,&
'@  (other source terms could be second order in time)',        /,&
'@',                                                            /,&
'@  Verify   the parameters.',                                  /,&
'@  and cs_user_lagr_model.',                                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2148 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME',   /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  source terms coming from module', A11,                      /,&
'@    will not be computed as second order',                    /,&
'@    in this version, despite user settings chosen below',     /,&
'@',                                                            /,&
'@       For scalar ',       i10,                               /,&
'@       THETSS    ISSO2T',                                     /,&
'@',      e12.4,      i10,                                      /,&
'@',                                                            /,&
'@  (other source terms could be second order in time)',        /,&
'@',                                                            /,&
'@  Verify   the parameters given by the interface,',           /,&
'@  cs_user_parameters.f90, and', a6,                           /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2511 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE A POSITIVE REAL',                        /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2530 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    WITH ID i10 CDTVAR MUST BE A  STRICTLY POSITIVE REAL',    /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  CDTVAR multiplier coefficient applied to the',              /,&
'@  timestep for the  resolution of variable I.',               /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
2604 format( &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA'               ,/,&
'@    ========='                                               ,/,&
'@    ',a21,' MUST BE AN INTEGER EQUAL TO 0, 10, 11, 20, 21,'  ,/,&
'@    30 OR 31'                                                ,/,&
'@   IT HAS VALUE ',i10                                        ,/,&
'@'                                                            ,/,&
'@   The calculation could NOT run.'                           ,/,&
'@'                                                            ,/,&
'@ Check the input data.',                                      /,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 2610 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    DONNEES INCOHERENTES',                                    /,&
'@',    a31,i10,                                                /,&
'@',    a31,i10,                                                /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2621 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@  Gravity is taken into account',                             /,&
'@',             3e14.5,                                        /,&
'@    without solving for temperature or energy',               /,&
'@    (ISCALT = ', i10,   ')',                                  /,&
'@',                                                            /,&
'@  The calculation will run.',                                 /,&
'@',                                                            /,&
'@  The above options are not incompatible',                    /,&
'@    but gravity is more often activated when density is',     /,&
'@    variable with temperature (natural convection)',          /,&
'@   this could be an error',                                   /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2622 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@     Gravity is taken into account',                          /,&
'@   in the turbulence source terms  (',a6,' = ', i10,   ')',   /,&
'@    without solving for temperature or energy',               /,&
'@',                                                            /,&
'@  The calculation will run.',                                  /&
'@',                                                            /,&
'@  gravity usualy afects turbulence only via density effects', /,&
'@    Check that density is variable.',                         /,&
'@     other than via temperature',                             /,&
'@  this could be by user defined density.',                    /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 4300 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@ Scalar ', a16,                                               /,&
'@  "is_temperature" must be an integer equal to 1 or 0',       /,&
'@   it has value', i1,                                         /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4350 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    SCALAR ', a16,                                            /,&
'@    SIGMAS(',i10,   ') MUST BE   A  POSITIVE REAL',           /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  SIFMAS(I) is the turbulent Prandtl turbulent',              /,&
'@    associated to scalar I.',                                 /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4370 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    SCALAR ', a16,                                            /,&
'@    SCAMAX(',i10,   ') MUST BE A STRICTLY POSITIVE REAL',     /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  SCAMAX(I) is the maximum acceptable value for',             /,&
'@    scalar   I, which is a variance',                         /,&
'@  with ICLVFL(I) = 2, value SCAMAX must be',                  /,&
'@   strictly  positive.',                                      /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4380 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    SCALAR ', a16,                                            /,&
'@    RVARFL MUST BE A STRICTLY POSITIVE REAL',     /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  RVARFL is the coefficient R for the scalar I (which is', /,&
'@ a variance) related to the dissipation equation source term',/,&
'@    - (1/R) rho scalaire epsilon/k',                          /,&
'@',                                                            /,&
'@ Check the input data.',                                      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 5003 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    ROTATIONAL PERIODICITY IS NOT COMPATIBLE WITH THE',       /,&
'@    ENHANCED PRESSSURE-VELOCITY COUPLING  or ALE METHOD',     /,&
'@      IN THE CURRENT VERSION',                                /,&
'@',                                                            /,&
'@  The calculation CAN NOT run.',                              /,&
'@',                                                            /,&
'@  At least one rotational periodicity has been defined',      /,&
'@  The flag IPUCOU is defined as',  i10,                       /,&
'@    (enhanced coupling for IPUCOU=1).',                       /,&
'@  The ALE fag  IALE is defined as',  i10,                     /,&
'@    (method activated if IALE=1)',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 5005 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    PERIODICITY IS INCOMPATIBLE WITH THIS METHOD FOR',        /,&
'@    COMPUTING DISTANCE TO WALL IN THE CURRENT VERSION',       /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  At least one periodicity has been defined',                 /,&
'@  the parameters specified need the calculation of the',      /,&
'@  distance to the wall (Rij-eps LRR with wall echo term,   ', /,&
'@     van Driest damping   or k-omega SST).',                  /,&
'@  The method for computing wall distance  :',                 /,&
'@    ICDPAR = ', i10,   ' is not compatible with',             /,&
'@     periodicity',                                            /,&
'@',                                                            /,&
'@  Recommendation: use ICDPAR = 1 or -1.',                     /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 5008 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@ ROTATION PERIODICITY IS NOT COMPATIBLE WITH RADIATIVE HEAT', /,&
'@      TRANSFER IN SEMI TRANSPARENT MEDIA',                    /,&
'@      (in the current version)',                              /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  At least one rotational periodicity has been defined',      /,&
'@  Flag IIRAYO is equal to', i10, '.',                         /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
6005 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    THIS METHOD FOR COMPUTING WALL DISTANCES',                /,&
'@    IS NOT COMPATIBLE WITH PARALLEL COMPUTING',               /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  The present CPU has rank', i10,                             /,&
'@',                                                            /,&
'@  Wall distance is necessary for RSTM wall echo terms, or,',  /,&
'@     van Driest damping or k-omega SST).',                    /,&
'@  The method for computing wall distance defined by',         /,&
'@    ICDPAR = ', i10, ' does not allow for parallel computing', /,&
'@',                                                            /,&
'@  Use ICDPAR = 1 or -1.',                                     /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 7010 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    NUMBER OF ITERATIONS FOR FLUID INITIALZATION WITH ALE',   /,&
'@',                                                            /,&
'@  NALINF MUST BE A POSITIVE INTEGER',                         /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  Verify the parameters.',                                    /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 7030 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    MAX number of iterations for implicit ALE method',        /,&
'@',                                                            /,&
'@  NALIMX MUST BE A POSITIVE INTEGER',                         /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  Verify the parameters given in the interface or usipsu.',   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 7040 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    COUPLING PRECISION FOR ALE',                              /,&
'@',                                                            /,&
'@  EPALIM MUST BE A REAL NUMBER,  STRICTLY  POSITIVE',         /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 7050 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@   INITIALISATION  ITERATION  FOR ALE',                       /,&
'@',                                                            /,&
'@  ITALIN must be =   0 or 1',                                 /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  Verify   the parameters  given  in usipsu.',                /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 8000 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@               COMPRESSIBLE FLOW MODULE',                     /,&
'@    T0 AND  P0 MUST BE STRICTLY POSITIVE REAL NUMBERS',       /,&
'@    Here they have values:',                                  /,&
'@                   T0 = ', e14.5,                             /,&
'@                   P0 = ', e14.5,                             /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 8010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@               COMPRESSIBLE FLOW MODULE',                     /,&
'@                                                            ',/,&
'@    THE THERMAL CONDUCTIVITY MUST BE                        ',/,&
'@    A STRICTLY POSITIVE REAL NUMBER                         ',/,&
'@    IT HAS VALUE ',E12.4                                     ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  Check uscfx2.                                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8020 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@               COMPRESSIBLE FLOW MODULE',                     /,&
'@                                                            ',/,&
'@    THE VOLUMIC VISCOSITY MUST BE                           ',/,&
'@    A STRICTLY POSITIVE REAL NUMBER                         ',/,&
'@    IT HAS VALUE ',E12.4                                     ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  Check uscfx2.                                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
8030 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a33,                          ' MUST BE AN INTEGER',   /, &
'@    BETWEEN 1 and 3',                                         /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 8040 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@               COMPRESSIBLE FLOW MODULE',                     /,&
'@                                                            ',/,&
'@    THE POLYTROPIC COEFFICIENT FOR THE STIFFENED GAS LAW    ',/,&
'@    MUST BE A REAL NUMBER SUPERIOR TO 1.                    ',/,&
'@    IT HAS VALUE ',E12.4                                     ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  Check uscfx2.                                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8050 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@               COMPRESSIBLE FLOW MODULE',                     /,&
'@                                                            ',/,&
'@    THE SPECIFIC HEAT RATIO (CP0 / CV0)                     ',/,&
'@    MUST BE A REAL NUMBER STRICTLY SUPERIOR TO 1.           ',/,&
'@    CP0 HAS VALUE ',E12.4                                    ,/,&
'@    CV0 HAS VALUE ',E12.4                                    ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  Check uscfx2.                                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8060 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@               COMPRESSIBLE FLOW MODULE',                     /,&
'@                                                            ',/,&
'@    THE LIMIT PRESSSURE OF THE STIFFENED GAS LAW MUST BE    ',/,&
'@    ZERO IN IDEAL GAS OR IDEAL GAS MIX.                     ',/,&
'@    PSGINF HAS VALUE ',E12.4                                 ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  Check uscfx2.                                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9010  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@   UNSTEADY ROTOR/STATOR COUPLING IS NOT COMPATIBLE',         /,&
'@     WITH THE STEADY ALGORITHM',                              /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  Integer parameter IDTVAR was set to', i10,                  /,&
'@    through the User Interface or in cs_user_parameters.f90.',/,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9011  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@   UNSTEADY ROTOR/STATOR COUPLING IS NOT COMPATIBLE',         /,&
'@     WITH THE SPACE OR TIME VARIABLE TIME STEPS',             /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  Integer parameter IDTVAR was set to', i10,                  /,&
'@    through the User Interface or in cs_user_parameters.f90.',/,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9120  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@   THE VOF METHOD IS NOT COMPATIBLE WITH THE',                /,&
'@     DILATABLE OR LOW-MACH FLOWS ALGORITHMS',                 /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  Integer parameter IDILAT was set to', i10,                  /,&
'@    through the User Interface or in cs_user_parameters.f90.',/,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

!----
! End
!----

return
end subroutine
