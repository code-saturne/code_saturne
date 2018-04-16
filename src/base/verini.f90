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
use alstru
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

!===============================================================================

implicit none

! Arguments

integer          iok

! Local variables

character        chaine*80
integer          ii    , iis   , jj    , iisct
integer          iscal , iest  , iiesca, ivar
integer          f_id, n_fields
integer          indest, iiidef, istop
integer          kscmin, kscmax, ifcvsl
integer          keyvar, keysca
double precision scmaxp, scminp
double precision turb_schmidt

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

!===============================================================================
! 1. ENTREES SORTIES entsor : formats 1000
!===============================================================================

! --- Suite, Historiques, Listing

if (nthist.le.0.and.nthist.ne.-1) then
  write(nfecra,1210) 'NTHIST (Periode   Sortie Histo. )',nthist
  iok = iok + 1
endif

!===============================================================================
! 2. OPTIONS DU CALCUL : TABLEAUX DE optcal : formats 2000
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

!     Extrap de rho : necessairement rho variable
if (irovar.eq.0.and.iroext.gt.0) then
  write(nfecra,2005)iroext,irovar
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
endif

do ii = 1, nscal
  if (itytur.eq.4.or.ischtp.eq.2) then
    iiidef = 10
    if (vcopt%nswrsm.lt.iiidef) then
      write(nfecra,2125) chaine(1:16),iiidef,vcopt%nswrsm
    endif
  endif
enddo

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
  if (isso2t(iscal).ne.isno2t) then
    write(nfecra,2133) iscal,isso2t(iscal),isno2t
  endif
  if (ivsext(iscal).ne.iviext) then
    write(nfecra,2134) iscal,ivsext(iscal),iviext
  endif
enddo

if (ischtp.eq.2.and.vcopt%ibdtso.gt.1) then
  ! NB: this test does not prevent from incompatible user modifications of
  ! isno2t, thetav, etc.
  write(nfecra,1135)
  iok = iok + 1
endif

!     Test du theta de la diffusivite des scalaires et de Cp : ils doivent etre
!       variables en (en espace) si on les extrapole (en temps) (...)
if ( icpext.gt.0 .and. icp.lt.0 ) then
  write(nfecra,2135) icpext, icp
  iok = iok + 1
endif
do iscal = 1, nscal
  if (ivsext(iscal).gt.0) then
    call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
    if (ifcvsl.lt.0) then
      write(nfecra,2136) iscal, ivsext(iscal), ifcvsl
      iok = iok + 1
    endif
  endif
enddo

!     Pour les tests suivants : Utilise-t-on un estimateur d'erreur ?
indest = 0
do iest = 1, nestmx
  iiesca = iescal(iest)
  if (iiesca.gt.0) then
    indest = 1
  endif
enddo

!     Estimateurs incompatibles avec calcul a champ de vitesse
!       fige (on ne fait rien, sauf ecrire des betises dans le listing)
if (indest.eq.1.and.iccvfg.eq.1) then
  write(nfecra,2137)
  iok = iok + 1
endif

!     A priori, pour le moment, l'ordre 2 en temps
!       (rho, visc, termes sources N.S, theta vitesse)
!     est incompatible avec
!       - estimateurs
!       - ipucou
!       - iphydr et icalhy
!       - dt variable en espace ou en temps et stationnaire
!     Ici on s'arrete si on n'est pas dans le cas du schema std

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)

if ( (abs(vcopt%thetav-1.0d0).gt.1.d-3).or.               &
     (    thetsn       .gt.0.d0 ).or.                     &
     (    isno2t       .gt.0    ).or.                     &
     (    iroext       .gt.0    ).or.                     &
     (    thetvi       .gt.0.d0 ).or.                     &
     (    iviext       .gt.0    )    ) then
  if (indest.eq.1.or.ipucou.eq.1.or.                      &
                     iphydr.eq.2.or.                      &
      idtvar.eq.1.or.idtvar.eq.2.or.idtvar.lt.0) then
    write(nfecra,2140)                                    &
         vcopt%thetav,                                    &
         isno2t,thetsn,                                   &
         iroext,                                          &
         iviext,thetvi
    iok = iok + 1
  endif
endif

if (nterup.gt.1) then

  if (ipucou.eq.1.or.indest.eq.1.or.                              &
      ippmod(icompf).ge.0.or.iccvfg.eq.1.or.                     &
      idtvar.eq.-1) then
    write(nfecra,2141) nterup
    iok = iok + 1
  endif

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
    if ((    thetss(iscal)       .gt.0.d0 ).or.            &
        (    isso2t(iscal)       .gt.0    ).or.            &
        (    thetvs(iscal).gt.0.d0 ).or.                   &
        (    ivsext(iscal).gt.0    )    ) istop = 1
  enddo

  if (istop.ne.0) then
    write(nfecra,2146)
    iok = iok + 1
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
    if (thetss(iscalt).gt.0.d0 .or. isso2t(iscalt).gt.0) then
      write(nfecra,2148)                                           &
        'lagrangian ',iscal,thetss(iscalt),isso2t(iscalt), 'cs_user_lagr_model'
      iok = iok + 1
    endif
  endif
endif

!     A priori, pour le moment, l'ordre 2 en temps
!       n'est pas pris en compte pour les termes issus du rayonnement.
!       On pourrait le signaler et continuer : on s'arrete.
if (iirayo.gt.0) then
  if (iscalt.gt.0) then
    if (thetss(iscalt).gt.0.d0 .or. isso2t(iscalt).gt.0) then
      write(nfecra,2148)                                           &
        'rayonnement',iscal,thetss(iscal),isso2t(iscal),'usray1'
      iok = iok + 1
    endif
  endif
endif

! --- Solveurs iteratifs

! Il n'y a pas besoin de test sur les epsilons
!   Ce sont simplement des reels
!   Une valeur negative indique qu'on veut atteindre
!   le nombre d'iterations maximal

do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyvar, ii)
  if (ii.ge.1) then
    if (idircl(ii).ne.0.and.idircl(ii).ne.1) then
      call field_get_label(f_id, chaine)
      write(nfecra,2401) chaine(1:16),ii,idircl(ii)
      iok = iok + 1
    endif
  endif
enddo

! --- Suite de calcul

if (ileaux.ne.0.and.ileaux.ne.1) then
  write(nfecra,2200) 'ILEAUX',ileaux
  iok = iok + 1
endif
if (iecaux.ne.0.and.iecaux.ne.1) then
  write(nfecra,2200) 'IECAUX',iecaux
  iok = iok + 1
endif
! En LES, on previent que ce n'est pas malin de ne pas relire le fichier
!   auxiliaire
if (iecaux.eq.0.or.ileaux.eq.0) then
  if (itytur.eq.4) then
    write(nfecra,2420) iturb,ileaux,iecaux
  endif
endif


do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyvar, ii)
  if (ii.ge.1) then
    if (cdtvar(ii).le.0.d0) then
      call field_get_label(f_id, chaine)
      write(nfecra,2530) chaine(1:16),ii,cdtvar(ii)
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

  ! Turbulent flux model for scalar
  if (iturt(iscal).ne. 0.and.iturt(iscal).ne.10 .and. &
      iturt(iscal).ne.20.and.iturt(iscal).ne.30 .and. &
      iturt(iscal).ne.11.and.iturt(iscal).ne.21 .and. &
      iturt(iscal).ne.31                              &
                   ) then
    write(nfecra,2604) 'iturt  ',iturt(iscal)
    write(nfecra,2610)                         &
         'Index of the scalar: ', iscal,       &
         'Number of scalars: ', nscal

    iok = iok + 1
  endif

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

! --- Estimateurs  d'erreur pour Navier-Stokes

do iest = 1, nestmx
  iiesca = iescal(iest)
  if (iiesca.ne.0.and.iiesca.ne.1.and.iiesca.ne.2) then
    write(nfecra,2664) iest,iest,iiesca,            &
         iespre,iesder,iescor,iestot
    iok = iok + 1
  endif
enddo


! --- Distance a la paroi

if (ineedy.eq.1) then

  if (abs(icdpar).ne.1) then
    write(nfecra,2700) icdpar
    iok = iok + 1
  endif
  if (ntcmxy.lt.1) then
    write(nfecra,3100) 'NTCMXY',ntcmxy
    iok = iok + 1
  endif

  if (coumxy.le.0.d0) then
    write(nfecra,2740) 'COUMXY',coumxy
    iok = iok + 1
  endif
  if (yplmxy.le.0.d0) then
    write(nfecra,2740) 'YPLMXY',yplmxy
    iok = iok + 1
  endif

endif

!===============================================================================
! 3. TABLEAUX DE cstphy : formats 4000
!===============================================================================

!     LES
if (itytur.eq.4) then
  if (xlesfl.lt.0.d0) then
    write(nfecra,2511) 'XLESFL', xlesfl
    iok = iok + 1
  endif
  if (ales  .lt.0.d0) then
    write(nfecra,2511) 'ALES  ', ales
    iok = iok + 1
  endif
  if (bles  .lt.0.d0) then
    write(nfecra,2511) 'BLES  ', bles
    iok = iok + 1
  endif
  if (csmago.lt.0.d0) then
    write(nfecra,2511) 'CSMAGO', csmago
    iok = iok + 1
  endif
  if (cwale.lt.0.d0) then
    write(nfecra,2511) 'CWALE', cwale
    iok = iok + 1
  endif
  if (idries.eq.1.and.cdries.lt.0) then
    write(nfecra,2511) 'CDRIES', cdries
    iok = iok + 1
  endif
  if (iturb.eq.41) then
    if (xlesfd.lt.0.d0) then
      write(nfecra,2511) 'XLESFD', xlesfd
      iok = iok + 1
    endif
    if (smagmx.lt.0.d0) then
      write(nfecra,2511) 'SMAGMX', smagmx
      iok = iok + 1
    endif
  endif
endif

! --- Scalaires

if (nscal.gt.0) then

!     Scalaire passif, temperature, enthalpie, energie
  do ii = 1, nscal
    if (iscacp(ii).lt.0.or.iscacp(ii).gt.1) then
      call field_get_label(ivarfl(isca(ii)), chaine)
      write(nfecra,4300)chaine(1:16),ii,iscacp(ii)
      iok = iok + 1
    endif
  enddo

!     Mode de clipping dans le cas des variances
!       pour les scalaires non variance, on n'utilise pas ICLVFL.
!       On demande donc a l'utilisateur de ne pas y toucher
!       (ca permet d'etre sur qu'il sait ce qu'il fait)
  do ii = 1, nscal
    if (iscavr(ii).le.nscal.and.iscavr(ii).gt.0) then
      if (iclvfl(ii).ne.0.and.                                     &
         iclvfl(ii).ne.1.and.iclvfl(ii).ne.2) then
        call field_get_label(ivarfl(isca(ii)), chaine)
        write(nfecra,4330)chaine(1:16),ii,iclvfl(ii)
        iok = iok + 1
      endif
    elseif (iscavr(ii).eq.0) then
      if (iclvfl(ii).ne.-1) then
        call field_get_label(ivarfl(isca(ii)), chaine)
        write(nfecra,4331)chaine(1:16),ii,iclvfl(ii)
        iok = iok + 1
      endif
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

!     Si on n'utilise pas la borne inf de clipping, on demande
!      a l'utilisateur de ne pas y toucher (ca permet d'etre sur
!      qu'il sait ce qu'il fait)
  do ii = 1, nscal
    ! Get the min clipping
    call field_get_key_double(ivarfl(isca(ii)), kscmin, scminp)

    if (iscavr(ii).gt.0.and.iscavr(ii).le.nscal.and.               &
       iclvfl(ii).ne.2.and.abs(scminp+grand).ge.epzero) then
      call field_get_label(ivarfl(isca(ii)), chaine)
      write(nfecra,4360)chaine(1:16),ii,scminp,ii,iclvfl(ii)
      iok = iok + 1
    endif
  enddo

!     Si on n'utilise pas la borne sup de clipping, on demande
!      a l'utilisateur de ne pas y toucher (ca permet d'etre sur
!      qu'il sait ce qu'il fait)
  do ii = 1, nscal
    ! Get the max clipping
    call field_get_key_double(ivarfl(isca(ii)), kscmax, scmaxp)

    if (iscavr(ii).gt.0.and.iscavr(ii).le.nscal.and.               &
       iclvfl(ii).ne.2.and.abs(scmaxp-grand).ge.epzero) then
      call field_get_label(ivarfl(isca(ii)), chaine)
      write(nfecra,4361)chaine(1:16),ii,scmaxp,ii,iclvfl(ii)
      iok = iok + 1
    endif
  enddo

!     Valeur de la borne sup de clipping si on l'utilise
  do ii = 1, nscal
    ! Get the max clipping
    call field_get_key_double(ivarfl(isca(ii)), kscmax, scmaxp)

    if (iscavr(ii).gt.0.and.iscavr(ii).le.nscal.and.               &
       iclvfl(ii).eq.2.and.scmaxp.le.0.d0) then
      call field_get_label(ivarfl(isca(ii)), chaine)
      write(nfecra,4370)chaine(1:16),ii,scmaxp
      iok = iok + 1
    endif
  enddo

!     Warning sur Rvarfl < 0
  do ii = 1, nscal
    if (iscavr(ii).gt.0.and.iscavr(ii).le.nscal.and.               &
                           rvarfl(ii).le.0.d0) then
      call field_get_label(ivarfl(isca(ii)), chaine)
      write(nfecra,4380)chaine(1:16),ii,rvarfl(ii)
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
        if (iscacp(iis).eq.1) then
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
! 4. TABLEAUX DE period : formats 5000
!===============================================================================

! --- periodicite de rotation incompatible avec couplage
!       renforce vitesse pression et l'ALE
if (iperot.gt.0.and. (ipucou.ne.0.or.iale.ne.0)) then
  write(nfecra,5003)iperio,ipucou,iale
  iok = iok + 1
endif

! --- periodicite incompatible avec le mode de calcul
!       direct de la distance a la paroi
if (iperio.eq.1.and.ineedy.eq.1.and.abs(icdpar).eq.2) then
  write(nfecra,5005)iperio, icdpar
  iok = iok + 1
endif


! --- periodicite de rotation incompatible avec le rayonnement DOM
!if (iperio.gt.0.and.iirayo.gt.0) then ! de translation aussi ?
if (iperot.gt.0.and.iirayo.gt.0) then
  if (iirayo.eq.1) then
    write(nfecra,5008) iperio,  iirayo
    iok = iok + 1
  endif
endif

! --- periodicite de rotation douteuse avec rij
!      (et donc a fortiori avec ordre 2 sur Rij)
if (iperot.gt.0) then
  if (itytur.eq.3) then
    write(nfecra,5009)iperio,iturb
  endif
endif

!===============================================================================
! 5. TABLEAUX DE parall : formats 6000 (limitations)
!===============================================================================

! --- parallelisme incompatible avec le mode de calcul
!       direct de la distance a la paroi
if (irangp.ge.0.and.ineedy.eq.1.and.abs(icdpar).eq.2) then
  write(nfecra,6005)irangp, icdpar
  iok = iok + 1
endif

!===============================================================================
! 6. METHODE ALE (albase, alstru) : formats 7000
!===============================================================================

if (iale.ne.0 .and. iale.ne.1) then
  write(nfecra,7000)iale
  iok = iok + 1
endif

if (iale.eq.1) then

  if (nalinf.lt.0) then
    write(nfecra,7010)nalinf
    iok = iok + 1
  endif

  if (alpnmk.lt.0.d0 .or. alpnmk.gt.1.d0  .or.                   &
      gamnmk.lt.0.d0 .or. gamnmk.gt.1.d0  .or.                   &
      betnmk.lt.0.d0 .or. betnmk.gt.0.5d0) then
    write(nfecra,7020)alpnmk,betnmk,gamnmk
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
! 7. COMPRESSIBLE : formats 8000
!===============================================================================

if (ippmod(icompf).ge.0) then
  if (t0.le.0.d0.or.p0.le.0.d0) then
    write(nfecra,8000)t0,p0
    iok = iok + 1
  endif
  if (visls0(itempk).le.0.d0) then
    write(nfecra,8010) visls0(itempk)
    iok = iok + 1
  endif
  if (viscv0.lt.0.d0) then
    write(nfecra,8020) viscv0
    iok = iok + 1
  endif
  if (ieos.lt.1.or.ieos.gt.3) then
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
! 8. Unsteady rotor/stator coupling: 9000 formats
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
! 10. Cavitation modelling: 9100 formats
!===============================================================================

if (ivofmt.ge.0) then
  ! VOF method is not compatible with dilatable or low-mach algorithms
  if (idilat.gt.1) then
    write(nfecra,9120) idilat
    iok = iok + 1
  endif
endif

!===============================================================================
! 9. FORMATS VERIFICATION
!===============================================================================

#if defined(_CS_LANG_FR)

 1210 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a33,                          ' DOIT ETRE UN ENTIER',   /,&
'@    STRICTEMENT POSITIF OU EGAL A -1',                        /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 1135 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERREUR :     ARRET A L''ENTREE DES DONNEES',              /,&
'@    ========',                                                /,&
'@   LE CHOIX DU SCHEMA EN TEMPS ISCHTP = 2 N''EST PAS',        /,&
'@   COMPATIBLE AVEC IBDTSO > 1',                               /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2000 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE UN ENTIER',                              /,&
'@      STRICTEMENT POSITIF ET INFERIEUR OU EGAL A', i10,       /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute',                            /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2005 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    ON DEMANDE UNE EXTRAPOLATION TEMPORELLE DE RHO AVEC',     /,&
'@      IROEXT = ', i10,                                        /,&
'@    CECI EST INCOMPATIBLE AVEC RHO CONSTANT',                 /,&
'@      IROVAR = ', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute',                            /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2125 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES',                /,&
'@    =========',                                               /,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS',                    /,&
'@',                                                            /,&
'@   ORDRE 2 EN TEMPS OU LES',                                  /,&
'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE NSWRSM POUR',      /,&
'@     LA VARIABLE', a16, ' EST',  i10,                        /, &
'@     NSWRSM A ETE IMPOSE ICI A', i10,                         /,&
'@',                                                            /,&
'@  Le calcul sera execute',                                    /,&
'@',                                                            /,&
'@  Il est conseille de verifier les parametres donnes via',    /,&
'@  l''interface ou cs_user_parameters.f90.',                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2131 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES',                /,&
'@    =========',                                               /,&
'@   CHOIX DU SCHEMA EN TEMPS',                                 /,&
'@',                                                            /,&
'@     LE SCHEMA EN TEMPS POUR LA VITESSE EST D ORDRE 1',       /,&
'@       (THETAV = ', e10.2, ')',                               /,&
'@     CERTAINS TERMES SONT CEPENDANT PRIS A L''ORDRE 2 AVEC',  /,&
'@       LES CHOIX SUIVANTS :',                                 /,&
'@',                                                            /,&
'@ Parametres       ISTMPF ISNO2T ISTO2T IROEXT IVIEXT ICPEXT', /,&
'@ Valeurs entrees', 6I7,                                       /,&
'@',                                                            /,&
'@  Le calcul sera execute.',                                   /,&
'@',                                                            /,&
'@  Il est conseille de verifier les parametres donnes via',    /,&
'@  l''interface ou cs_user_parameters.f90.',                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2132 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES',                /,&
'@    =========',                                               /,&
'@   CHOIX DU SCHEMA EN TEMPS',                                 /,&
'@',                                                            /,&
'@     LE SCHEMA EN TEMPS POUR LA VITESSE EST D ORDRE 2',       /,&
'@       (THETAV = ', e10.2, ')',                               /,&
'@     CERTAINS TERMES SONT CEPENDANT PRIS A L''ORDRE 1 AVEC',  /,&
'@       LES CHOIX SUIVANTS :',                                 /,&
'@',                                                            /,&
'@ Parametres       ISTMPF ISNO2T ISTO2T IROEXT IVIEXT ICPEXT', /,&
'@ Valeurs entrees', 6I7,                                       /,&
'@',                                                            /,&
'@  Le calcul sera execute.',                                   /,&
'@',                                                            /,&
'@  Il est conseille de verifier les parametres donnes via',    /,&
'@  l''interface ou cs_user_parameters.f90.',                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2133 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES',                /,&
'@    =========',                                               /,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS',                    /,&
'@',                                                            /,&
'@   SCALAIRE ', i10,' ISSO2T = ', i10,                         /,&
'@     EST DIFFERENT DE ISNO2T',                                /,&
'@     ISNO2T = ', i10,                                         /,&
'@',                                                            /,&
'@  Le calcul sera execute',                                    /,&
'@',                                                            /,&
'@  Il est conseille de verifier les parametres donnes via',    /,&
'@  l''interface ou cs_user_parameters.f90.',                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2134 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES',                /,&
'@    =========',                                               /,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS',                    /,&
'@',                                                            /,&
'@   SCALAIRE ', i10,' IVSEXT = ', i10,                         /,&
'@     EST DIFFERENT DE IVIEXT',                                /,&
'@     IVIEXT = ', i10,                                         /,&
'@',                                                            /,&
'@  Le calcul sera execute',                                    /,&
'@',                                                            /,&
'@  Il est conseille de verifier les parametres donnes via',    /,&
'@  l''interface ou cs_user_parameters.f90.',                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2135 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES',                /,&
'@    =========',                                               /,&
'@   CHOIX INCOMPATIBLE POUR LE SCHEMA EN TEMPS',               /,&
'@',                                                            /,&
'@     La  chaleur massique est extrapolee en temps avec',      /,&
'@       ICPEXT = ', i10,                                       /,&
'@     Pour cela, elle doit etre variable, or',                 /,&
'@       ICP    = ', i10,                                       /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute',                             /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@    - desactiver le choix d''extrapolation de Cp en temps',   /,&
'@      ou',                                                    /,&
'@    - imposer Cp variable',                                   /,&
'@         (et le renseigner alors via l''interface ou usphyv)',/,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2136 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES',                /,&
'@    =========',                                               /,&
'@   CHOIX INCOMPATIBLE POUR LE SCHEMA EN TEMPS',               /,&
'@',                                                            /,&
'@   Scalaire ISCAL = ', i10,                                   /,&
'@     La  diffusivite      est extrapolee en temps avec',      /,&
'@       IVSEXT(ISCAL) = ', i10,                                /,&
'@     Pour cela, elle doit etre variable, or',                 /,&
'@       scalar_diffusivity_id = ', i10, ' pour ce champ.'      /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute',                             /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@    - desactiver le choix d''extrapolation en temps',         /,&
'@                                     de la diffusivite',      /,&
'@      ou',                                                    /,&
'@    - imposer la diffusivite variable',                       /,&
'@         (et le renseigner alors via l''interface ou usphyv)',/,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2137 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES',                /,&
'@    =========',                                               /,&
'@   CHOIX INCOMPATIBLE POUR LES ESTIMATEURS D''ERREUR',        /,&
'@',                                                            /,&
'@  On a active un ou plusieurs estimateurs d''erreur pour',    /,&
'@    Navier-Stokes dans un calcul a champ de vitesse fige.',   /,&
'@    Le ou les estimateurs ne seront pas calcules.',           /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute',                             /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@      desactiver les estimateurs d erreur ou',                /,&
'@               le calcul a champ de vitesse fige (ICCVFG)',   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2140 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS', /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  On souhaite utiliser un schema en temps d''ordre 2 :',      /,&
'@      U : THETA = ', e12.4,                                   /,&
'@      Termes sources Navier-Stokes: ISNO2T = ', i10,          /,&
'@                                    THETSN = ', e12.4,        /,&
'@      Masse volumique             : IROEXT = ', i10,          /,&
'@                                    THETRO = ', e12.4,        /,&
'@      Viscosite                   : IVIEXT = ', i10,          /,&
'@                                    THETVI = ', e12.4,        /,&
'@  La version actuelle ne le permet pas lorsque l''une des',   /,&
'@    options suivantes a ete activee (c''est le cas ici) :',   /,&
'@    - utilisation d''un estimateur d''erreur (IESCAL)',       /,&
'@    - couplage instationnaire (IPUCOU)',                      /,&
'@    - prise en compte specifique de la pression',             /,&
'@      hydrostatique (IPHYDR et ICALHY)',                      /,&
'@    - pas de temps variable en temps ou en espace  ou',       /,&
'@      algorithme stationnaire (IDTVAR)',                      /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2141 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS', /,&
'@',                                                            /,&
'@  On souhaite utiliser un couplage',                          /,&
'@    vitesse-pression par point fixe (NTERUP = ', i10,         /,&
'@  La version actuelle ne le permet pas lorsque l''une des',   /,&
'@    options suivantes a ete activee (c''est le cas ici) :',   /,&
'@    - utilisation d''un estimateur d''erreur (IESCAL)',       /,&
'@    - couplage instationnaire (IPUCOU)',                      /,&
'@    - algorithme stationnaire (IDTVAR=-1)',                   /,&
'@    - module compressible (IPPMOD(ICOMPF)>=0)',               /,&
'@    - champ de vitesse fige (ICCVFG=1)',                      /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2142 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS', /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Avec le modele de turbulence',                              /,&
'@    k-epsilon (ITURB = ', i10, ') couple (IKECOU = ', i10,') :',/,&
'@    la version courante ne permet pas de traiter les',        /,&
'@    equations du modele k-epsilon a l''ordre 2 en temps avec',/,&
'@    couplage.',                                               /,&
'@    Une ou plusieurs valeurs parmi les suivantes ne sont',    /,&
'@    donc pas permises :',                                     /,&
'@',                                                            /,&
'@       THETST    ISTO2T     THETA K   THETA EPS',             /,&
'@',      e12.4,      i10,      e12.4,      e12.4,              /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2143 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS', /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Avec le modele de turbulence',                              /,&
'@    v2f (ITURB = ', i10,   ') couple (IKECOU = ', i10,    ') :',/,&
'@    la version courante ne permet pas de traiter les',        /,&
'@    equations du modele v2f a l''ordre 2 en temps avec',      /,&
'@    couplage.',                                               /,&
'@    Une ou plusieurs valeurs parmi les suivantes ne sont',    /,&
'@    donc pas permises :',                                     /,&
'@',                                                            /,&
'@       THETST    ISTO2T     THETA K   THETA EPS',             /,&
'@',      e12.4,      i10,      e12.4,      e12.4,              /,&
'@     THETA PHI    THETA FB',                                  /,&
'@',       e12.4,      e12.4,                                   /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2144 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS', /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Avec le modele de turbulence',                              /,&
'@    k-omega (ITURB = ', i10,   ') couple (IKECOU = ', i10,') :',/,&
'@    la version courante ne permet pas de traiter les',        /,&
'@    equations du modele k-omega a l''ordre 2 en temps avec',  /,&
'@    couplage.',                                               /,&
'@    Une ou plusieurs valeurs parmi les suivantes ne sont',    /,&
'@    donc pas permises :',                                     /,&
'@',                                                            /,&
'@       THETST    ISTO2T     THETA K   THETA OMEGA',           /,&
'@',      e12.4,      i10,      e12.4,      e12.4,              /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2145 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS', /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Avec le modele de turbulence',                              /,&
'@    Spallart-Allmaras (ITURB = ', i10,   ')',                 /,&
'@    la version courante ne permet pas de traiter les',        /,&
'@    l''ordre 2 en temps.',                                    /,&
'@    couplage.',                                               /,&
'@    Une ou plusieurs valeurs parmi les suivantes ne sont',    /,&
'@    donc pas permises :',                                     /,&
'@',                                                            /,&
'@       THETST    ISTO2T     THETA K   THETA OMEGA',           /,&
'@',      e12.4,      i10,      e12.4,      e12.4,              /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2146 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS', /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  La version actuelle ne permet pas de modifier le schema en',/,&
'@    temps lorsqu''une physique particuliere est activee',     /,&
'@    (combustion, charbon, electrique).',                      /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2147 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS', /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Les termes sources provenant du module',                    /,&
'@    Lagrangien ne sont pas traites a l''ordre 2 en temps',    /,&
'@    dans la version courante malgre le choix utilisateur',    /,&
'@    suivant :',                                               /,&
'@',                                                            /,&
'@     Pour Navier-Stokes    Pour la turbulence',               /,&
'@       THETSN    ISNO2T      THETST    ISTO2T',               /,&
'@',      e12.4,      i10,      e12.4,      i10,                /,&
'@',                                                            /,&
'@  (Les autres termes sources pourraient etre traites a',      /,&
'@   l''ordre 2)',                                              /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface,',          /,&
'@    cs_user_parameters.f90, et cs_user_lagr_model.',          /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2148 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS', /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Les termes sources provenant du module', A11,               /,&
'@    ne sont pas traites a l''ordre 2 en temps',               /,&
'@    dans la version courante malgre le choix utilisateur',    /,&
'@    suivant :',                                               /,&
'@',                                                            /,&
'@       Pour le scalaire ', i10,                               /,&
'@       THETSS    ISSO2T',                                     /,&
'@',      e12.4,      i10,                                      /,&
'@',                                                            /,&
'@  (Les autres termes sources pourraient etre traites a',      /,&
'@   l''ordre 2)',                                              /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface,',          /,&
'@    cs_user_parameters.f90, et', a6,                          /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2200 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE UN ENTIER EGAL A 0 OU 1',                /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2401 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    IDIRCL(',i10,   ') DOIT ETRE UN ENTIER EGAL A 0 OU 1',    /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  IDIRCL(I) indique si le code doit decaler la diagonale de', /,&
'@    la matrice de la variable I en l''absence de Dirichlet',  /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2420 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    RISQUE DE PERTE D''INFORMATION EN CALCUL SUITE',          /,&
'@',                                                            /,&
'@  Le calcul sera engage.',                                    /,&
'@',                                                            /,&
'@  Un modele de LES a ete active par ITURB = ', i10,           /,&
'@    mais on a desactive l''ecriture ou la lecture du fichier',/,&
'@    suite auxiliaire :',                                      /,&
'@    ILEAUX = ', i10,   '    IECAUX = ', i10,                  /,&
'@  Bien que ce fichier ne soit pas necessaire a la poursuite', /,&
'@    d''un calcul, il contient neanmoins des informations',    /,&
'@    qui permettent d''eviter les perturbations numeriques',   /,&
'@    au moment des suites.',                                   /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2511 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE UN REEL POSITIF',                        /,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2530 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    CDTVAR(',i10,   ') DOIT ETRE UN REEL STRICTEMENT POSITIF',/,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  CDTVAR(I) est le coefficient multiplicatif applique au pas',/,&
'@    de temps pour la resolution de la variable I.',           /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2604 format(&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    ========='                                               ,/,&
'@    ',a7,' DOIT ETRE UN ENTIER EGAL A 0, 10,         20,'    ,/,&
'@       OU 30'                                                ,/,&
'@    IL VAUT ICI ',i10                                        ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@'                                                            ,/,&
'@  Verifier les parametres donnes via l''interface'           ,/,&
'@    ou cs_user_parameters.f90.'                              ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 2610 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    DONNEES INCOHERENTES',                                    /,&
'@',    a31,i10,                                                /,&
'@',    a31,i10,                                                /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2621 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@  ON DEMANDE LA PRISE EN COMPTE DE L''ACCELERATION DE LA',    /,&
'@    PESANTEUR', 3e14.5,                                       /,&
'@    SANS RESOLUTION D''UNE VARIABLE TEMPERATURE OU ENERGIE',  /,&
'@    (ISCALT = ', i10,   ')',                                  /,&
'@',                                                            /,&
'@  Le calcul sera engage.',                                    /,&
'@',                                                            /,&
'@  Il n''y a pas d''incompatibilite a prendre en compte',      /,&
'@    l''acceleration de la pesanteur sans effets thermiques,', /,&
'@    mais, souvent, la masse volumique depend essentiellement',/,&
'@    de la temperature et la combinaison des options du',      /,&
'@    present calcul est caracteristique d''un oubli.',         /,&
'@',                                                            /,&
'@  Il est conseille de verifier les parametres donnes via',    /,&
'@  l''interface ou cs_user_parameters.f90.',                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2622 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@  ON DEMANDE LA PRISE EN COMPTE DES TERMES DE GRAVITE DANS',  /,&
'@    LES EQUATIONS DE LA TURBULENCE (',a6,' = ', i10,   ')',   /,&
'@    SANS RESOLUTION D''UNE VARIABLE TEMPERATURE OU ENERGIE',  /,&
'@',                                                            /,&
'@  Le calcul sera engage.',                                    /,&
'@',                                                            /,&
'@  Si des effets de gravite sont recherches, il convient de',  /,&
'@    s''assurer que la masse volumique est variable.',         /,&
'@  Le nombre de Prandtl turbulent sera pris egal a 1.',        /,&
'@  Elle peut varier en fonction d''autres grandeurs que',      /,&
'@    la temperature ou l''enthalpie ; si c''est le cas, ce',   /,&
'@    message pourra etre ignore ; sinon, verifier usipsu',     /,&
'@    ou imposer une variation de masse volumique.',            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2664 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    CHOIX DU CALCUL DES ESTIMATEURS D''ERREUR',               /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  L''indicateur IESCAL relatif a',                            /,&
'@    l''estimateur d''erreur numero IEST = ', i10,   ' pour',  /,&
'@    Navier-Stokes doit etre un entier egal a 0, 1 ou 2.',     /,&
'@  Il vaut ici : IESCAL(',i10,  ') = ', i10,                   /,&
'@',                                                            /,&
'@  Rq. : les valeurs possibles de IEST sont :',                /,&
'@        IESPRE = ', i10,                                      /,&
'@        IESDER = ', i10,                                      /,&
'@        IESCOR = ', i10,                                      /,&
'@        IESTOT = ', i10,                                      /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2700 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    CHOIX DU MODE DE CALCUL DE LA DISTANCE A LA PAROI',       /,&
'@',                                                            /,&
'@  ICDPAR DOIT ETRE UN ENTIER EGAL A -2, -1, 1 ou 2',          /,&
'@  IL VAUT ICI',  i10,                                         /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2740 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE UN REEL STRICTEMENT POSITIF.',           /,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 3100 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE UN ENTIER STRICTEMENT POSITIF',          /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 4300 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    SCALAIRE ', a16,                                          /,&
'@    ISCACP(',i10,   ') DOIT ETRE UN ENTIER EGAL A 0 ou 1',    /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4330 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    SCALAIRE ', a16,                                          /,&
'@    ICLVFL(',i10,   ') DOIT ETRE UN ENTIER EGAL A 0, 1 OU 2', /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  ICLVFL(I) indique le mode de clipping du scalaire I',       /,&
'@    lorsqu il s agit d une variance de fluctuations.',        /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4331 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    SCALAIRE ', a16,                                          /,&
'@    ICLVFL(',i10,   ') N EST UTILISE QUE POUR LES VARIANCES', /,&
'@    IL VAUT ICI', i10,                                        /,&
'@    ALORS QUE LE SCALAIRE N EST PAS UNE VARIANCE.',           /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  ICLVFL(I) indique le mode de clipping du scalaire I',       /,&
'@    lorsqu il s agit d une variance de fluctuations.',        /,&
'@    Il n est pas utilise pour les autres scalaires.',         /,&
'@  L utilisateur est invite a ne pas modifier ICLVFL pour',    /,&
'@    les scalaires qui ne sont pas des variances.',            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4350 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    SCALAIRE ', a16,                                          /,&
'@    SIGMAS(',i10,   ') DOIT ETRE UN REEL POSITIF',            /,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  SIFMAS(I) est le nombre de Prandtl turbulent associe',      /,&
'@    au scalaire I.',                                          /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4360 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    SCALAIRE ', a16,                                          /,&
'@    SCAMIN(',i10,   ') VAUT ICI', e14.5,                      /,&
'@      AVEC ICLVFL(',i10,   ') = ', i10,                       /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  SCAMIN(I) est la valeur minimale acceptee pour le',         /,&
'@    scalaire I. Lorsque le scalaire est une variance',        /,&
'@    (iscavr(I) > 0) la valeur de SCAMIN n est prise en',      /,&
'@    compte que si ICLVFL(I) = 2',                             /,&
'@  Si l utilisateur souhaite effectivement que le',            /,&
'@    scalaire I (en fait, une variance) soit limite a SCAMIN', /,&
'@    (positif) il faut imposer ICLVFL = 2 dans usipsu.',       /,&
'@  Si l utilisateur souhaite utiliser l option ICLVFL = 1',    /,&
'@    il est invite a ne pas modifier SCAMIN dans usipsu.',     /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4361 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    SCALAIRE ', a16,                                          /,&
'@    SCAMAX(',i10,   ') VAUT ICI', e14.5,                      /,&
'@      AVEC ICLVFL(',i10,   ') = ', i10,                       /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  SCAMAX(I) est la valeur maximale acceptee pour le',         /,&
'@    scalaire I. Lorsque le scalaire est une variance',        /,&
'@    (iscavr(I) > 0) la valeur de SCAMAX n est prise en',      /,&
'@    compte que si ICLVFL(I) = 2',                             /,&
'@  Si l utilisateur souhaite effectivement que le',            /,&
'@    scalaire I (en fait, une variance) soit limite a SCAMAX', /,&
'@    (positif) il faut imposer ICLVFL = 2 dans usipsu.',       /,&
'@  Si l utilisateur souhaite utiliser l option ICLVFL = 1',    /,&
'@    il est invite a ne pas modifier SCAMAX dans usipsu.',     /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4370 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    SCALAIRE ', a16,                                          /,&
'@    SCAMAX(',i10,   ') DOIT ETRE UN REEL STRICTEMENT POSITIF',/,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  SCAMAX(I) est la valeur maximale acceptee pour le',         /,&
'@    scalaire I, ici une variance',                            /,&
'@  Avec ICLVFL(I) = 2, la valeur de SCAMAX doit donc etre',    /,&
'@   strictment positive.',                                     /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4380 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    SCALAIRE ', a16,                                          /,&
'@    RVARFL(',i10,   ') DOIT ETRE UN REEL STRICTEMENT POSITIF',/,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  RVARFL(I) est le coefficient R pour le scalaire I (qui est',/,&
'@    une variance) intervenant dans le terme de dissipation :',/,&
'@    - (1/R) rho scalaire epsilon/k',                          /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 5003 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    LA PERIODICITE DE ROTATION N''EST PAS COMPATIBLE AVEC LE',/,&
'@      COUPLAGE VITESSE PRESSION RENFORCE  OU LA METHODE',     /,&
'@      ALE DANS LA VERSION COURANTE',                          /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Au moins une periodicite de rotation a ete definie.',       /,&
'@  L''indicateur IPUCOU a ete positionne a', i10,              /,&
'@    dans l''interface ou usipsu (couplage renforce pour',     /,&
'@    IPUCOU=1).',                                              /,&
'@  L''indicateur IALE a ete positionne a', i10,                /,&
'@    dans l''interface ou usipsu (methode activee si IALE=1)', /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 5005 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    MODE DE CALCUL DE LA DISTANCE A LA PAROI INCOMPATIBLE',   /,&
'@      AVEC LA PERIODICITE',                                   /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Au moins une periodicite a ete definie.',                   /,&
'@  Les parametres de calcul specifies necessitent le calcul',  /,&
'@    la distance a la paroi (Rij-eps LRR avec echo de paroi,', /,&
'@    LES avec amortissement de van Driest ou k-omega SST).',   /,&
'@  Le mode de calcul de la distance a la paroi defini par',    /,&
'@    ICDPAR = ', i10,   ' ne permet pas de prendre en compte', /,&
'@    la periodicite.',                                         /,&
'@',                                                            /,&
'@  Utiliser ICDPAR = 1 ou -1.',                                /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 5008 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    LA PERIODICITE DE ROTATION N''EST PAS COMPATIBLE AVEC LE',/,&
'@      RAYONNEMENT SEMI-TRANSPARENT  (ORDONNEES DISCRETES)',   /,&
'@      DANS LA VERSION COURANTE',                              /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Au moins une periodicite de rotation a ete definie.',       /,&
'@  L''indicateur IIRAYO a ete positionne a', i10,              /,&
'@    dans l''interface ou usray1.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 5009 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    DES DEFAUTS PEUVENT SE RENCONTRER LORS DE L''UTILISATION',/,&
'@      DE LA PERIODICITE DE ROTATION EN RIJ-EPSILON.',         /,&
'@',                                                            /,&
'@  Au moins une periodicite de rotation a ete definie.',       /,&
'@  L''indicateur ITURB a ete positionne a', i10,               /,&
'@',                                                            /,&
'@  Le calcul peut etre execute.',                              /,&
'@    Les defauts eventuels evoques proviennent de la prise en',/,&
'@    compte de la rotation du tenseur de viscosite orthotrope',/,&
'@    Il a cependant en general une influence faible de sorte', /,&
'@    que les tests menes jusqu''a present n''ont pas fait',    /,&
'@    apparaitre de probleme.',                                 /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 6005 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    MODE DE CALCUL DE LA DISTANCE A LA PAROI INCOMPATIBLE',   /,&
'@      AVEC LE PARALLELISME',                                  /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Le processeur courant est de rang', i10,                    /,&
'@  Les parametres de calcul specifies necessitent le calcul',  /,&
'@    la distance a la paroi (Rij-eps LRR avec echo de paroi,', /,&
'@    LES avec amortissement de van Driest ou k-omega SST).',   /,&
'@  Le mode de calcul de la distance a la paroi defini par',    /,&
'@    ICDPAR = ', i10, ' ne permet pas de prendre en compte',   /,&
'@    le parallelisme.',                                        /,&
'@',                                                            /,&
'@  Utiliser ICDPAR = 1 ou -1.',                                /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 7000 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    INDICATEUR DE METHODE ALE',                               /,&
'@',                                                            /,&
'@  IALE DOIT VALOIR 0 OU 1',                                   /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface ou usipph.',/,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 7010 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    NOMBRE D''ITERATIONS D''INITIALISATION DU FLUIDE EN ALE', /,&
'@',                                                            /,&
'@  NALINF DOIT ETRE UN ENTIER POSITIF',                        /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface ou usipsu.',/,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 7020 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    COEFFICIENTS DE LA METHODE DE NEWMARK NON VALIDES',       /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  ALPNMK doit etre compris entre 0 et 1',                     /,&
'@  BETNMK doit etre compris entre 0 et 1/2',                   /,&
'@  GAMNMK doit etre compris entre 0 et 1',                     /,&
'@  On a ici :',                                                /,&
'@',                                                            /,&
'@       ALPNMK      BETNMK      GAMNMK',                       /,&
'@',      e12.4,      e12.4,      e12.4,                        /,&
'@',                                                            /,&
'@  Verifier les parametres.',                                  /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 7030 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    NOMBRE D''ITERATIONS MAX DE COUPLAGE IMPLICITE EN ALE',   /,&
'@',                                                            /,&
'@  NALIMX DOIT ETRE UN ENTIER POSITIF',                        /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres.',                                  /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 7040 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    PRECISION DU COUPLAGE IMPLICITE EN ALE',                  /,&
'@',                                                            /,&
'@  EPALIM DOIT ETRE UN REEL STRICTEMENT POSITIF',              /,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface ou usipsu.',/,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 7050 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    ITERATION D''INITIALISATION DE L''ALE',                   /,&
'@',                                                            /,&
'@  ITALIN DOIT VALOIR 0 OU 1',                                 /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes dans usipsu.',               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 8000 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========   MODULE COMPRESSIBLE',                         /,&
'@',                                                            /,&
'@    T0 ET P0 DOIVENT ETRE DES REELS STRICTEMENT POSITIFS',    /,&
'@    ILS VALENT ICI :',                                        /,&
'@                   T0 = ', e14.5,                             /,&
'@                   P0 = ', e14.5,                             /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 8010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========   MODULE COMPRESSIBLE                         ',/,&
'@                                                            ',/,&
'@    LA CONDUCTIVITE THERMIQUE DOIT ETRE                     ',/,&
'@    UN REEL POSITIF STRICTEMENT                             ',/,&
'@    ELLE A POUR VALEUR ',E12.4                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier uscfx2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8020 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========   MODULE COMPRESSIBLE                         ',/,&
'@                                                            ',/,&
'@    LA VISCOSITE EN VOLUME DOIT ETRE                        ',/,&
'@    UN REEL POSITIF                                         ',/,&
'@    ELLE A POUR VALEUR ',E12.4                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier uscfx2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8030 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a33,                          ' DOIT ETRE UN ENTIER',   /,&
'@    ENTRE 1 ET 3',                                            /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 8040 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========   MODULE COMPRESSIBLE                         ',/,&
'@                                                            ',/,&
'@    LE COEFFICIENT POLYTROPIQUE DE LA LOI ''STIFFENED GAS'' ',/,&
'@    DOIT ETRE UN REEL SUPERIEUR A 1.                        ',/,&
'@    IL A POUR VALEUR ',E12.4                                 ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier uscfx2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8050 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========   MODULE COMPRESSIBLE                         ',/,&
'@                                                            ',/,&
'@    LE RAPPORT DES CHALEURS SPECIFIQUES (CP0 / CV0)         ',/,&
'@    DOIT ETRE UN REEL STRICTEMENT SUPERIEUR A 1.            ',/,&
'@    CP0 A POUR VALEUR ',E12.4                                ,/,&
'@    CV0 A POUR VALEUR ',E12.4                                ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier uscfx2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8060 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========   MODULE COMPRESSIBLE                         ',/,&
'@                                                            ',/,&
'@    LA PRESSION LIMITE DE LA LOI ''STIFFENED GAS'' DOIT ETRE',/,&
'@    NULLE EN GAZ PARFAIT OU MELANGE DE GAZ PARFAIT          ',/,&
'@    PSGINF A POUR VALEUR ',E12.4                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier uscfx2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9010  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@   LE COUPLAGE ROTOR/STATOR INSTATIONNAIRE N''EST PAS',       /,&
'@     COMPATIBLE AVEC L''ALGORITHME STATIONNAIRE',             /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  L''indicateur IDTVAR a ete positionne a', i10,              /,&
'@    par l''interface ou dans cs_user_parameter.f90',          /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9011  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@   LE COUPLAGE ROTOR/STATOR INSTATIONNAIRE N''EST PAS',       /,&
'@     COMPATIBLE AVEC LES PAS DE TEMPS VARIABLES EN ESPACE',   /,&
'@     OU EN TEMPS',                                            /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  L''indicateur IDTVAR a ete positionne a', i10,              /,&
'@    par l''interface ou dans cs_user_parameter.f90',          /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9120  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@   LA METHODE VOF N''EST PAS COMPATIBLE AVEC LES',            /,&
'@     ALGORITHMES D''ECOULEMENTS DILATABLES OU BAS-MACH',      /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  L''indicateur IDILAT a ete positionne a', i10,              /,&
'@    par l''interface ou dans cs_user_parameter.f90',          /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

#else

 1210 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a33,                          ' MUST BE AN INTEGER',    /,&
'@    LARGER OR EQUAL TO  1 or EQUAL TO  -1',                   /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2000 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' MUST BE AN INTEGER',                               /,&
'@      STRICTLY POSITIVE AND LESS THAN or EGAL TO',i10,        /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@  The calculation could NOT run.',                            /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2005 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    TEMPORAL EXTRAPOLATION OF DENSITY RHO REQUESTED,',        /,&
'@      BUT IROEXT = ', i10,                                    /,&
'@    THIS IS INCOMPATIBLE WITH RHO = CONSTANT',                /,&
'@      IROVAR = ', i10,                                        /,&
'@',                                                            /,&
'@  The calculation could NOT run.',                            /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2125 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHEN READING INPUT DATA',                /,&
'@    =========',                                               /,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME',                    /,&
'@',                                                            /,&
'@   ORDRE 2 EN TEMPS or the',                                  /,&
'@   THE VALUE RECOMMANDED FOR THE PARAMETER NSWRSM FOR',       /,&
'@        VARIABLE', a16, ' IS',    i10,                        /,&
'@     NSWRSM IS NOW IMPOSED AS',  i10,                         /,&
'@',                                                            /,&
'@  computation will go on',                                    /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 1135 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR :      WHEN READING INPUT DATA',                    /,&
'@    ========',                                                /,&
'@   CHOICE OF TIME-SCHEME ISCHTP = 2 IS NOT COMPATIBLE WITH',  /,&
'@   IBDTSO > 1',                                               /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
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
'@     CERTAIN TERMS ARE HOWEVER SECOND ORDER IN TIME WITH',    /,&
'@       THE FOLLOWING SETTINGS:',                              /,&
'@',                                                            /,&
'@ parameters       ISTMPF ISNO2T ISTO2T IROEXT IVIEXT ICPEXT', /,&
'@ Values  entered', 6I7,                                       /,&
'@',                                                            /,&
'@  computation will go on.',                                   /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
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
'@     CERTAIN TERMS ARE HOWEVER FIRST ORDER IN TIME  WITH',    /,&
'@       THE FOLLOWING SETTINGS:',                              /,&
'@',                                                            /,&
'@ parameters       ISTMPF ISNO2T ISTO2T IROEXT IVIEXT ICPEXT', /,&
'@ Values  entered', 6I7,                                       /,&
'@',                                                            /,&
'@  computation will go on.',                                   /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
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
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
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
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2135 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHEN READING INPUT DATA',                /,&
'@    =========',                                               /,&
'@  INCOMPATIBILITY FOR TIME DISCRETISATION SCHEME',            /,&
'@',                                                            /,&
'@     Specific heat is extrapolated in time with',             /,&
'@       ICPEXT = ', i10,                                       /,&
'@    in which case it should be variable, or',                 /,&
'@       ICP    = ', i10,                                       /,&
'@',                                                            /,&
'@  Computation will NOT go on',                                /,&
'@',                                                            /,&
'@  Verify   the parameters',                                   /,&
'@    - deactivate xtrapolation of Cp in time',                 /,&
'@      or',                                                    /,&
'@    - define Cp as variable',                                 /,&
'@         (define its variation law in the GUI ou usphyv)',    /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2136 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHEN READING INPUT DATA',                /,&
'@    =========',                                               /,&
'@  INCOMPATIBILITY FOR TIME DISCRETISATION SCHEME',            /,&
'@',                                                            /,&
'@   Scalar   ISCAL = ', i10,                                   /,&
'@    Diffusivity   is extrapolated in time wih',               /,&
'@       IVSEXT(ISCAL) = ', i10,                                /,&
'@     it should thus be a variable, or',                       /,&
'@       scalar_diffusivity_id = ', i10, ' for this field.'     /,&
'@',                                                            /,&
'@ Computation will  NOT  proceed',                             /,&
'@',                                                            /,&
'@  Verify the parameters',                                     /,&
'@    - deactivate intepolation in time',                       /,&
'@                                     for diffusivity',        /,&
'@      or',                                                    /,&
'@    - impose diffusivite variable',                           /,&
'@         (define its variation law in the GUI ou usphyv)',    /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2137 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHEN READING INPUT DATA',                /,&
'@    =========',                                               /,&
'@   CHOICE INCOMPATIBLE FOR ERROR ESTIMATES',                  /,&
'@',                                                            /,&
'@  One or several error estimates are activated for',          /,&
'@    Navier-Stokes in simulation with frozen velocity field.', /,&
'@    estimates will not be computed',                          /,&
'@',                                                            /,&
'@ Computation will  NOT  proceed',                             /,&
'@',                                                            /,&
'@  Verify   the parameters :',                                 /,&
'@      desactivate  ERROR ESTIMATES  (ICCVFG)',                /,&
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2140 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME',   /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  A second ordre time-scheme was requested:',                 /,&
'@      U,V,W : THETA = ', e12.4,                               /,&
'@     Source terme in Navier-Stokes: ISNO2T = ', i10,          /,&
'@                                    THETSN = ', e12.4,        /,&
'@      Density                     : IROEXT = ', i10,          /,&
'@                                    THETRO = ', e12.4,        /,&
'@      Viscosity                   : IVIEXT = ', i10,          /,&
'@                                    THETVI = ', e12.4,        /,&
'@  Current version does not allow this in combination with',   /,&
'@  one of the following option (which has been activated ):',  /,&
'@    - Error estimation (IESCAL)',                             /,&
'@    - reinforced U-P coupling (IPUCOU)',                      /,&
'@    - specific treatment of hydrostatic pressure',            /,&
'@      contribution  (IPHYDR et ICALHY)',                      /,&
'@    - time-step variable with space or iteration or',         /,&
'@      steady-state   algorithm(IDTVAR)',                      /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2141 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME',   /,&
'@',                                                            /,&
'@  Pressure-Velocity coupling by fixed-point method',          /,&
'@    was selected  by setting   NTERUP = ', i10,               /,&
'@  Current version does not allow this in combination with',   /,&
'@  one of the following options (which has been activated ):', /,&
'@    - Error estimation (IESCAL)',                             /,&
'@    - reinforced U-P coupling (IPUCOU)',                      /,&
'@    - time-step variable with space or iteration or',         /,&
'@      steady-state   algorithm(IDTVAR=-1)',                   /,&
'@    - compressible module (IPPMOD(ICOMPF)>=0)',               /,&
'@    - frozen velocity field (ICCVFG=1)',                      /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
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
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
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
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
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
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
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
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2146 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME',   /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  The current version does not allow changing the time',      /,&
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
'@  (other termes sources could be second order in time',       /,&
'@             )',                                              /,&
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
'@  (Les autres termes sources pourraient etre traites a',      /,&
'@   l''ordre 2)',                                              /,&
'@',                                                            /,&
'@  Verify   the parameters given by the interface,',           /,&
'@  cs_user_parameters.f90, and', a6,                           /,&
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2200 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' MUST BE AN INTEGER EQUAL  0  OR 1',                /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2401 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    IDIRCL(',i10,   ') MUST BE AN INTEGER EQUAL  0  OR 1',    /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  IDIRCL(I) tells if the diagonal of the matrix for variable',/,&
'@  I should be shifted in the absence of Dirichlet condition', /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2420 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    RISQUE DE PERTE D''INFORMATION EN CALCUL SUITE',          /,&
'@',                                                            /,&
'@  The calculation will run.',                                  /&
'@',                                                            /,&
'@ A turbulence model was activated by ITURB = ', i10,          /,&
'@    but writing to auxiliary restart file was de-activated',  /,&
'@',                                                            /,&
'@    ILEAUX = ', i10,   '    IECAUX = ', i10,                  /,&
'@  Although this file is not necessary to restart',            /,&
'@  a computation, it does contain information that avoid',     /,&
'@   numerical perturbations when restarting a computation',    /,&
'@',                                                            /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
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
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
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
'@    CDTVAR(',i10,   ') MUST BE A  STRICTLY POSITIVE REAL',    /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  CDTVAR(I) multiplier coefficient applied to the',           /,&
'@  timestep for the  resolution of variable I.',               /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
2604 format( &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA'               ,/,&
'@    ========='                                               ,/,&
'@    ',a7,' MUST BE AN INTEGER EQUAL TO 0, 10, 11, 20, 21,'   ,/,&
'@    30 OR 31'                                                ,/,&
'@   IT HAS VALUE ',i10                                        ,/,&
'@'                                                            ,/,&
'@   The calculation could NOT run.'                           ,/,&
'@'                                                            ,/,&
'@ Check the input data given through the User Interface'      ,/,&
'@   or in cs_user_parameters.f90.'                            ,/,&
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
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
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
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
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
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2664 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    CHOIX DU CALCUL DES ESTIMATEURS D''ERREUR',               /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@   Flag IESCAL related to',                                   /,&
'@   error estimate number IEST = ', i10,   ' for',             /,&
'@    Navier-Stokes MUST BE AN INTEGER egal to 0, 1 or 2.',     /,&
'@  It IS EQUAL : IESCAL(',i10,  ') = ', i10,                   /,&
'@',                                                            /,&
'@  Rq. : the possible values of IEST are    :',                /,&
'@        IESPRE = ', i10,                                      /,&
'@        IESDER = ', i10,                                      /,&
'@        IESCOR = ', i10,                                      /,&
'@        IESTOT = ', i10,                                      /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2700 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    Choice of method for computing distance to the wall',     /,&
'@',                                                            /,&
'@  ICDPAR MUST BE AN INTEGER  EQUAL TO -2, -1, 1 or 2',        /,&
'@  IL IS EQUAL',  i10,                                         /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2740 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' MUST BE   A  STRICTLY POSITIVE REAL.',             /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 3100 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' MUST BE AN INTEGER  STRICTLY  POSITIVE',           /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 4300 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    SCALAR ', a16,                                            /,&
'@    ISCACP(',i10,   ') MUST BE AN INTEGER EQUAL TO 1 OR 0',   /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4330 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    SCALAR ', a16,                                            /,&
'@    ICLVFL(',i10,   ') MUST BE AN INTEGER  EGAL A 0, 1 or 2', /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  ICLVFL(I) defines the type of clipping of scalar I',        /,&
'@    when it is a variance of  fluctuations.',                 /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4331 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    SCALAR ', a16,                                            /,&
'@    ICLVFL(',i10,   ') is only used for  VARIANCES',          /,&
'@   IT HAS VALUE', i10,                                        /,&
'@    BUT THE SCALAR IS NOT A VARIANCE',                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  ICLVFL(I) flags the type of clipping for scalar I',         /,&
'@',                                                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
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
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4360 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    SCALAR ', a16,                                            /,&
'@    SCAMIN(',i10,   ') IS EQUAL', e14.5,                      /,&
'@      AVEC ICLVFL(',i10,   ') = ', i10,                       /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  SCAMIN(I) is the  minimale acceptable value for',           /,&
'@    scalaire I. When this scalar is a variance',              /,&
'@    (iscavr(I) > 0) value of SCAMIN is only used if',         /,&
'@                  ICLVFL(I) = 2',                             /,&
'@',                                                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4361 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    SCALAR ', a16,                                            /,&
'@    SCAMAX(',i10,   ') IS EQUAL', e14.5,                      /,&
'@      AVEC ICLVFL(',i10,   ') = ', i10,                       /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  SCAMAX(I) is the maximum acceptable value for',             /,&
'@    scalar  I. When this is a variance',                      /,&
'@    (iscavr(I) > 0) the value of SCAMAX is only used',        /,&
'@               if ICLVFL(I) = 2',                             /,&
'@',                                                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
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
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
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
'@    RVARFL(', i10,  ') MUST BE A STRICTLY POSITIVE REAL',     /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  RVARFL(I) is the coefficient R for the scalar I (which is', /,&
'@ a variance) related to the dissipation equation sourceterme',/,&
'@    - (1/R) rho scalaire epsilon/k',                          /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
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
'@  Flag IIRAYO is equal to', i10,                              /,&
'@    in usray1.',                                              /,&
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 5009 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    DEFECTS CAN APPEAR WHEN USING A COMBINATION OF',/,          &
'@     ANGULAR PERIODICITY (ROTATION) AND RSTM RIJ-EPSILON.',   /,&
'@',                                                            /,&
'@  At least one rotational periodicity has been defined',      /,&
'@  Flag for turb ITURB is = ', i10,                            /,&
'@',                                                            /,&
'@  Job can run.',                                              /,&
'@',                                                            /,&
'@  The defects are related to the turbulent transport terms',  /,&
'@  in the Re stress equations (equivalent to an anisotropic',  /,&
'@  diffusion tensor), but these terms are generaly small',     /,&
'@',                                                            /,&
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
 7000 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    FLAG FOR ALE  METHOD',                                    /,&
'@',                                                            /,&
'@  IALE should be = 0 or 1',                                   /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  Verify   the parameters given  in   interface or usipph.',  /,&
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
 7020 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    NON VALID COEFFICIENTS IN  NEWMARK METHOD',               /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  ALPNMK MUST BE   between  0 and  1',                        /,&
'@  BETNMK MUST BE   between  0 and 1/2',                       /,&
'@  GAMNMK MUST BE   between  0 and 1',                         /,&
'@  We have here:',                                             /,&
'@',                                                            /,&
'@       ALPNMK      BETNMK      GAMNMK',                       /,&
'@',      e12.4,      e12.4,      e12.4,                        /,&
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
'@  Verify the parameters given in the interface or usipsu.',   /,&
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

#endif

!----
! End
!----

return
end subroutine
