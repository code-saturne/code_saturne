!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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

!===============================================================================

implicit none

! Arguments

integer          iok

! Local variables

character        chaine*80, chain2*80
integer          ii    , iis   , jj    , iisct, kval
integer          iscal , iest  , iiesca, ivar
integer          nbsccp
integer          c_id, f_id, f_dim, n_fields, ippf
integer          ipp   , nbccou, imrgrl
integer          iokpre, indest, iiidef, istop
integer          kscmin, kscmax, ifcvsl
integer          keyvar, keysca
double precision arakfr, scmaxp, scminp

character(len=3), dimension(3) :: nomext3
character(len=4), dimension(3) :: nomext63

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

if (ncapt.lt.0.or.ncapt.gt.ncaptm) then
  write(nfecra,1230)ncaptm,ncapt
  iok = iok + 1
endif

if (nthist.le.0.and.nthist.ne.-1) then
  write(nfecra,1210) 'NTHIST (Periode   Sortie Histo. )',nthist
  iok = iok + 1
endif

do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyipp, ippf)
  if (ippf.le.1) cycle
  call field_get_dim (f_id, f_dim)
  do c_id = 1, min(f_dim, 3)
    ipp = ippf + c_id - 1
    if (ihisvr(ipp,1).gt.ncapt.or.                                  &
      (ihisvr(ipp,1).lt.0.and.ihisvr(ipp,1).ne.-1) ) then
      call field_get_label(f_id, chaine)
      if (f_dim .eq. 3) then
        chaine = trim(chaine) // nomext3(c_id)
      else if (f_dim .eq. 6) then
        chaine = trim(chaine) // nomext63(c_id)
      endif
      write(nfecra,1240)chaine(1:16),ipp,ncapt,ihisvr(ipp,1)
      iok = iok + 1
    endif
  enddo
enddo

do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyipp, ippf)
  if (ippf.le.1) cycle
  call field_get_dim (f_id, f_dim)
  do c_id = 1, min(f_dim, 3)
    ipp = ippf + c_id - 1
    if ((ihisvr(ipp,1).gt.0.and.ihisvr(ipp,1).lt.ncapt)) then
      do jj = 1, ihisvr(ipp,1)
        if (ihisvr(ipp,jj+1).le.0.or.ihisvr(ipp,jj+1).gt.ncapt) then
          call field_get_label(f_id, chaine)
          if (f_dim .eq. 3) then
            chaine = trim(chaine) // nomext3(c_id)
          else if (f_dim .eq. 6) then
            chaine = trim(chaine) // nomext63(c_id)
          endif
          write(nfecra,1250) chaine(1:16),ipp,jj+1,ncapt,ihisvr(ipp,jj+1)
          iok = iok + 1
        endif
      enddo
    endif
  enddo
enddo

do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keylog, kval)
  if (kval.ne.0 .and. kval.ne.1) then
    call field_get_label(f_id, chaine)
    write(nfecra,1260)chaine(1:16),kval
    iok = iok + 1
  endif
enddo

if (ntlist.ne.-1.and.ntlist.le.0) then
  write(nfecra,1210) 'NTLIST (Periode   Sortie Listing)',ntlist
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

! --- Thermal model

if (itherm.lt.0 .or. itherm.gt.3) then
  write(nfecra,2050) itherm
  iok = iok + 1
endif

! --- Rho et visc constants ou variables

if (irovar.ne.0.and.irovar.ne.1) then
  write(nfecra,2201)'IROVAR',irovar
  iok = iok + 1
endif
if (ivivar.ne.0.and.ivivar.ne.1) then
  write(nfecra,2201)'IVIVAR',ivivar
  iok = iok + 1
endif

! --- Definition des equations, schema en temps, schema convectif

do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyvar, ii)
  if (ii.ge.1) then
    if ((iconv (ii).ne.0.and.iconv (ii).ne.1).or.                            &
        (istat (ii).ne.0.and.istat (ii).ne.1).or.                            &
        (idiff (ii).ne.0.and.idiff (ii).ne.1).or.                            &
        (idifft(ii).ne.0.and.idifft(ii).ne.1).or.                            &
        (thetav(ii).gt.1.d0.or.thetav(ii).lt.0.d0).or.                       &
        (blencv(ii).gt.1.d0.or.blencv(ii).lt.0.d0).or.                       &
        (ischcv(ii).lt.0.and.ischcv(ii).gt.3).or.                            &
        (isstpc(ii).lt.0.and.isstpc(ii).gt.3)   ) then
      call field_get_label(f_id, chaine)
      write(nfecra,2100) chaine(1:16),                            &
                         istat(ii),iconv(ii),                     &
                         idiff(ii),idifft(ii),                    &
                         thetav(ii),blencv(ii),ischcv(ii),        &
                         isstpc(ii)
      iok = iok + 1
    endif
  endif
enddo


!     Extrap de rho : necessairement rho variable
if (irovar.eq.0.and.iroext.gt.0) then
  write(nfecra,2005)iroext,irovar
  iok = iok + 1
endif

!     Coherence des vitesses
!       Pour le moment, theta est fixe automatiquement dans modini
!       On conserve quand meme le test pour plus tard puisqu'il est ecrit.
if (abs(thetav(iv)-thetav(iu)).gt.epzero.or.       &
     abs(thetav(iw)-thetav(iu)).gt.epzero) then
  write(nfecra,2111) thetav(iu),thetav(iv), &
       thetav(iw)
  iok = iok + 1
endif


!     Theta pression : vaut 1
!       Pour le moment, theta est fixe automatiquement dans modini
!         (on ne devrait donc jamais voir cet affichage)
!       On conserve quand meme le test pour plus tard puisqu'il est ecrit.
jj = ipr
if (abs(thetav(jj)-1.0d0).gt.epzero) then
  call field_get_label(ivarfl(jj), chaine)
  write(nfecra,2112) thetav(jj)
  iok = iok + 1
endif

!     En LES il y a des verification de coherence supplementaires.
!       (simple avertissement si on s'ecarte des choix std)
!        mais stop si on fait plus de 5% d'upwind
!     Schema centre sans/avec test de pente, nwsrsm
if (itytur.eq.4) then
  do ii = 1,3
    if (ii.eq.1) jj = iu
    if (ii.eq.2) jj = iv
    if (ii.eq.3) jj = iw
    call field_get_label(ivarfl(iu), chaine)
    chaine = trim(chaine) // nomext3(ii)
    if (abs(thetav(jj)-0.5d0).gt.epzero) then
      write(nfecra,2121) chaine(1:16),thetav(jj)
    endif
    if (blencv(jj).lt.0.95d0) then
      write(nfecra,2127) chaine(1:16),blencv(jj)
      iok = iok + 1
    elseif (abs(blencv(jj)-1.d0).gt.epzero) then
      write(nfecra,2122) chaine(1:16),blencv(jj)
    endif
    if (isstpc(jj).eq.0) then
      write(nfecra,2123) chaine(1:16),isstpc(jj)
    endif
  enddo
endif
if (itytur.eq.4.or.ischtp.eq.2) then
  do ii = 1,3
    if (ii.eq.1) jj = iu
    if (ii.eq.2) jj = iv
    if (ii.eq.3) jj = iw
    call field_get_label(ivarfl(iu), chaine)
    chaine = trim(chaine) // nomext3(ii)
    iiidef = 10
    if (nswrsm(jj).ne.iiidef) then
      write(nfecra,2125) chaine(1:16),iiidef,nswrsm(jj)
    endif
  enddo
  jj = ipr
  call field_get_label(ivarfl(jj), chaine)
  iiidef = 5
  if (nswrsm(jj).ne.iiidef) then
    write(nfecra,2125) chaine(1:16),iiidef,nswrsm(jj)
  endif
endif
do ii = 1, nscal
  if (itytur.eq.4) then
    jj    = isca(ii)
    call field_get_label(ivarfl(jj), chaine)
    if (abs(thetav(jj)-0.5d0).gt.epzero) then
      write(nfecra,2121) chaine(1:16),thetav(jj)
    endif
    if (blencv(jj).lt.0.95d0) then
      write(nfecra,2127) chaine(1:16),blencv(jj)
      iok = iok + 1
    elseif (abs(blencv(jj)-1.d0).gt.epzero) then
      write(nfecra,2122) chaine(1:16),blencv(jj)
    endif
    if (isstpc(jj).eq.1) then
      write(nfecra,2124) chaine(1:16),isstpc(jj)
    endif
  endif
  if (itytur.eq.4.or.ischtp.eq.2) then
    iiidef = 10
    if (nswrsm(jj).ne.iiidef) then
      write(nfecra,2125) chaine(1:16),iiidef,nswrsm(jj)
    endif
  endif
enddo

!     Test du theta de la viscosite secondaire, du flux de masse et
!     de la viscosite par rapport a celui de la vitesse
jj = iu
if (abs(thetav(jj)-1.d0).lt.epzero.and.                 &
     (istmpf.eq.2.or.                                   &
     isno2t.ne.0.or.                                    &
     isto2t.ne.0.or.                                    &
     iroext.ne.0.or.                                    &
     iviext.ne.0.or.                                    &
     icpext.ne.0   ) ) then
  write(nfecra,2131) thetav(jj),                        &
       istmpf,isno2t,isto2t,                            &
       iroext,iviext,icpext
endif
if (abs(thetav(jj)-0.5d0).lt.epzero.and.                &
     (istmpf.ne.2.or.                                   &
     isno2t.ne.1.or.                                    &
     isto2t.ne.1.or.                                    &
     iroext.ne.1.or.                                    &
     iviext.ne.1.or.                                    &
     icpext.ne.1   ) ) then
  write(nfecra,2132) thetav(jj),                        &
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

if (ischtp.eq.2.and.ibdtso(iu).gt.1) then
  ! NB: this test does not prevent from incompatible user modifications of
  ! isno2t, thetav, etc.
  write(nfecra,1135)
  iok = iok + 1
endif

!     Test du theta de la diffusivite des scalaires et de Cp : ils doivent etre
!       variables en (en espace) si on les extrapole (en temps) (...)
if ( icpext.gt.0 .and. icp.le.0 ) then
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
if ( (abs(thetav(iu)-1.0d0).gt.1.d-3).or.                 &
     (abs(thetav(iv)-1.0d0).gt.1.d-3).or.                 &
     (abs(thetav(iw)-1.0d0).gt.1.d-3).or.                 &
     (    thetsn       .gt.0.d0 ).or.                     &
     (    isno2t       .gt.0    ).or.                     &
     (    thetro       .gt.0.d0 ).or.                     &
     (    iroext       .gt.0    ).or.                     &
     (    thetvi       .gt.0.d0 ).or.                     &
     (    iviext       .gt.0    )    ) then
  if (indest.eq.1.or.ipucou.eq.1.or.                      &
      iphydr.eq.1.or.iphydr.eq.2.or.icalhy.eq.1.or.       &
      idtvar.eq.1.or.idtvar.eq.2.or.idtvar.lt.0) then
    write(nfecra,2140)                                    &
         thetav(iu),thetav(iv),thetav(iw),                &
         isno2t,thetsn,                                   &
         iroext,thetro,                                   &
         iviext,thetvi
    iok = iok + 1
  endif
endif

!     Iterations sur navsto
!     Doit etre un entier superieur ou egal a 1
!     Pour le moment, on interdit NTERUP > 1 et
!       estimateurs, matrices poids, pression hydrostatique
!       et algo stationnaire

if (nterup.le.0) then
  WRITE(NFECRA,3100) 'NTERUP', NTERUP
  iok = iok + 1
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
  if ((     thetst       .gt.0.d0 ).or.                    &
       (    isto2t       .gt.0    ).or.                    &
       (abs(thetav(ik )-1.0d0).gt.epzero).or.              &
       (abs(thetav(iep)-1.0d0).gt.epzero) ) then
    write(nfecra,2142)iturb,ikecou,                        &
         thetst,isto2t,                                    &
         thetav(ik ),thetav(iep)
    iok = iok + 1
  endif
endif
if (iturb.eq.50.and.ikecou.eq.1) then
  if ((     thetst       .gt.0.d0 ).or.                    &
      (    isto2t       .gt.0    ).or.                     &
      (abs(thetav(ik  )-1.0d0).gt.epzero).or.              &
      (abs(thetav(iep )-1.0d0).gt.epzero).or.              &
      (abs(thetav(iphi)-1.0d0).gt.epzero).or.              &
      (abs(thetav(ifb )-1.0d0).gt.epzero) ) then
    write(nfecra,2143)iturb,ikecou,                        &
         thetst,isto2t,                                    &
         thetav(ik  ),thetav(iep ),                        &
         thetav(iphi),thetav(ifb )
    iok = iok + 1
  endif
endif
if (iturb.eq.51.and.ikecou.eq.1) then
  if ((    thetst       .gt.0.d0 ).or.                     &
     (    isto2t       .gt.0    ).or.                      &
     (abs(thetav(ik  )-1.0d0).gt.epzero).or.               &
     (abs(thetav(iep )-1.0d0).gt.epzero).or.               &
     (abs(thetav(iphi)-1.0d0).gt.epzero).or.               &
     (abs(thetav(ial )-1.0d0).gt.epzero) ) then
    write(nfecra,2143)iturb,ikecou,                        &
         thetst,isto2t,                                    &
         thetav(ik  ),thetav(iep ),                        &
         thetav(iphi),thetav(ial )
    iok = iok + 1
  endif
endif
if (iturb.eq.60.and.ikecou.eq.1) then
  if ((    thetst       .gt.0.d0 ).or.                     &
       (    isto2t       .gt.0    ).or.                    &
       (abs(thetav(ik  )-1.0d0).gt.epzero).or.             &
       (abs(thetav(iomg)-1.0d0).gt.epzero) ) then
    write(nfecra,2144)iturb,ikecou,                        &
         thetst,isto2t,                                    &
         thetav(ik  ),thetav(iomg)
    iok = iok + 1
  endif
endif
if (iturb.eq.70) then
  if ((thetst .gt.0.d0).or.                                &
      (isto2t .gt. 0).or.                                  &
      (abs(thetav(inusa)-1.0d0).gt.epzero)) then
    write(nfecra,2145)iturb,thetst,isto2t,thetav(inusa)
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
    if ( (abs(thetav(ivar)-1.0d0).gt.1.d-3) ) istop = 1
  enddo
  if ((thetsn .gt.0.d0).or.                                &
      (isno2t .gt.0   ).or.                                &
      (thetro .gt.0.d0).or.                                &
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
        'lagrangian ',iscal,thetss(iscalt),isso2t(iscalt), 'uslag1'
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

! --- Algorithme stationnaire
if (idtvar.lt.0) then
  do f_id = 0, n_fields-1
    call field_get_key_int(f_id, keyvar, ii)
    if (ii.ge.1) then
      if (relaxv(ii).gt.1d0.or.relaxv(ii).lt.0d0) then
        call field_get_label(f_id, chaine)
        write(nfecra,2149) chaine(1:16),relaxv(ii)
        iok = iok + 1
      endif
    endif
  enddo
  if ((relaxv(iv).ne.relaxv(iu))                  &
       .or.(relaxv(iw).ne.relaxv(iu)) ) then
    write(nfecra,2150) relaxv(iu),relaxv(iv),     &
         relaxv(iw)
    iok = iok + 1
  endif
!       L'algorithme stationnaire n'est pas compatible avec le module Lagrangien
  if (iilagr.ne.0) then
    write(nfecra,2151) iilagr
    iok = iok + 1
  endif
!       L'algorithme stationnaire n'est pas compatible avec la LES
  if (itytur.eq.4) then
    write(nfecra,2152) iturb
    iok = iok + 1
  endif
endif

! --- Reconstruction des gradients
if (imrgra.gt.16 .or. imrgra.lt.-16) then
  write(nfecra,2205) 'IMRGRA', imrgra
  iok = iok + 1
endif

imrgrl = abs(imrgra)
imrgrl = modulo(imrgrl,10)

! On verifie l'angle de non orthogonalite de selection du
!   voisinage etendu dans le cas du moindre carre qui l'utilise

if (imrgrl.eq.3.or.imrgrl.eq.6) then
  if (anomax.gt.pi*0.5d0.or.anomax.lt.0.d0) then
    write(nfecra,2206) anomax, imrgra
  endif
endif

! Extrapolation : indetermination possible par mc,
!     necessitant un traitement particulier dans gradmc,
!     pour lequel on fait certaines hypotheses
if (imrgrl.ne.0.and.imrgrl.ne.4) then
  if (      (abs(extrag(ipr)-1.d0).gt.epzero)             &
      .and. (abs(extrag(ipr)     ).gt.epzero)) then
    write(nfecra,2207) imrgra, extrag(ipr)
    iok = iok + 1
  endif
endif

! Les nombres de sweeps n'ont pas a etre verifies :
!  ce sont simplement des entiers (negatifs si on veut etre sur de ne
!  *jamais* entrer dans les boucles)

do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyvar, ii)
  if (ii.ge.1) then
    if (imligr(ii).gt.1) then
      call field_get_label(f_id, chaine)
      write(nfecra,2300) chaine(1:16),ii,imligr(ii)
      iok = iok + 1
    endif
  endif
enddo

do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyvar, ii)
  if (ii.ge.1) then
    if (ircflu(ii).ne.1.and.ircflu(ii).ne.0) then
      call field_get_label(f_id, chaine)
      write(nfecra,2310) chaine(1:16),ii,ircflu(ii)
      iok = iok + 1
    endif
  endif
enddo
! Non reconstruction des flux en SOLU n'a pas de sens pour la convection
do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyvar, ii)
  if (ii.ge.1) then
    if (ircflu(ii).eq.0.and.ischcv(ii).eq.0.and.blencv(ii).ne.0.d0) then
      call field_get_label(f_id, chaine)
      write(nfecra,2311) chaine(1:16),ii,ircflu(ii),ii,ischcv(ii)
      iok = iok + 1
    endif
  endif
enddo

! Il n'y a pas besoin de test sur les epsilons
!   Ce sont simplement des reels
!   Une valeur negative indique qu'on veut atteindre
!   le nombre d'iterations maximal

do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyvar, ii)
  if (ii.ge.1) then
    if (climgr(ii).lt.1.d0) then
      call field_get_label(f_id, chaine)
      write(nfecra,2320) chaine(1:16),ii,climgr(ii)
      iok = iok + 1
    endif
  endif
enddo


! EXTRAG non nul permis uniquement pour la pression.
!        et dans ce cas egal a 1
do f_id = 0, n_fields-1
  call field_get_key_int(f_id, keyvar, ii)
  if (ii.ge.1) then
    if (abs(extrag(ii)     ).ge.epzero) then
      iokpre = 0
      if (ii.eq.ipr) then
        iokpre = 1
        if (abs(extrag(ii)-1.d0).ge.epzero) then
          call field_get_label(f_id, chaine)
          write(nfecra,2330) chaine(1:16),ii,extrag(ii)
          iok = iok + 1
        endif
      endif
      if (iokpre.eq.0) then
        call field_get_label(f_id, chaine)
        write(nfecra,2331) chaine(1:16),ii,extrag(ii)
        iok = iok + 1
      endif
    endif
  endif
enddo


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

! --- Reperage du temps et marche en temps

! Le nombre de pas de temps pourrait bien etre negatif : pas de test.

if (inpdt0.ne.0.and.inpdt0.ne.1) then
  write(nfecra,2200) 'INPDT0',inpdt0
  iok = iok + 1
endif

if (idtvar.lt.-1.or.idtvar.gt.2) then
  write(nfecra,2500) 'IDTVAR',idtvar
  iok = iok + 1
endif

if (idtvar.gt.0.and.varrdt.lt.0.d0) then
  write(nfecra,2510)'VARRDT', varrdt
  iok = iok + 1
endif

if (dtref .lt.0.d0) then
  write(nfecra,2510)'DTREF ', dtref
  iok = iok + 1
endif

if (dtmin.le.0.d0 .or. dtmax.le.0.d0 .or. dtmin.gt.dtmax) then
  write(nfecra,2520) dtmin, dtmax
  if (idtvar.gt.0) then
    iok = iok + 1
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

if (iptlro.ne.0.and.iptlro.ne.1) then
  write(nfecra,2200) 'IPTLRO',iptlro
  iok = iok + 1
endif

!     Si on est en pas de temps constant, on ne touche pas le pas de temps,
!     mais on indique juste les depassements de critere local.
if (idtvar.eq.0.and.iptlro.eq.1) then
  write(nfecra,2540) iptlro,idtvar
endif
if (idtvar.eq.-1.and.iptlro.eq.1) then
  write(nfecra,2541) iptlro,idtvar
endif

! --- Turbulence

!    Modele

if ( iturb.ne. 0.and.iturb.ne.10.and.iturb.ne.20.and.        &
     iturb.ne.21.and.iturb.ne.30.and.iturb.ne.31.and.        &
     iturb.ne.32.and.                                        &
     iturb.ne.40.and.iturb.ne.41.and.iturb.ne.42.and.        &
     iturb.ne.50.and.iturb.ne.51.and.iturb.ne.60.and.iturb.ne.70 ) then
  write(nfecra,2600) 'ITURB  ',iturb
  iok = iok + 1
endif

! Rotation curvature correction for eddy-viscosity models
if ( irccor.ne.0.and.irccor.ne.1 ) then
  write(nfecra,2601) 'IRCCOR ',irccor
  iok = iok + 1
endif

! Rotation curvature correction compatible only with RANS eddy-viscosity models
if ( irccor.eq.1.and.(itytur.ne.2 .and.itytur.ne.5 .and.      &
                      iturb .ne.60.and.iturb .ne.70) ) then
  write(nfecra,2602) iturb
  iok = iok + 1
endif

! In lagrangian with two-way coupling, k-omega SST is forbidden (not
! properly implemented)
if (iturb.eq.60 .and. iilagr.eq.2) then
  write(nfecra,2603) iilagr
  iok = iok + 1
endif

!     Methode des vortex pour la LES

if (ivrtex.ne.0 .and.ivrtex.ne.1) then
  write(nfecra,2200) 'IVRTEX ',ivrtex
  iok = iok + 1
endif
if (ivrtex.eq.1.and.itytur.ne.4) then
  write(nfecra,2606)itytur,ivrtex
  ivrtex = 0
endif

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
      iturt(iscal).ne.20.and.iturt(iscal).ne.30       &
                   ) then
    write(nfecra,2604) 'iturt  ',iturt(iscal)
    write(nfecra,2610)                         &
         'Index of the scalar: ', iscal,       &
         'Number of scalars: ', nscal

    iok = iok + 1
  endif

endif

if (iwallf.lt.0.or.iwallf.gt.6) then
  write(nfecra,2001)'IWALLF',6,iwallf
  iok = iok + 1
endif
if (iwallf.gt.2 .and.                                    &
    (iturb.eq.0 .or. iturb.eq.10 .or.                    &
    itytur.eq.4 .or. iturb.eq.70)) then
  write(nfecra,2209)iturb,iwallf
  iok = iok + 1
endif
if (iwallt.lt.0.or.iwallt.gt.1) then
  write(nfecra,2201)'IWALLT',iwallt
  iok = iok + 1
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
  !     Le choix de ICLKEP n'est possible qu'en k-eps ou v2f
  if (iturb.ne.60) then
    if (iclkep.ne.0.and.iclkep.ne.1) then
      WRITE(NFECRA,2201)'ICLKEP',ICLKEP
      iok = iok + 1
    endif
  endif
  if (ikecou.ne.0.and.ikecou.ne.1) then
    WRITE(NFECRA,2201)'IKECOU',IKECOU
    iok = iok + 1
  endif
  !     En k-eps a prod lin et en v2f on force IKECOU a 0
  if (ikecou.eq.1 .and. (iturb.eq.21 .or. itytur.eq.5)) then
    write(nfecra,2208)iturb,ikecou
    iok = iok + 1
  endif
  !     En stationnaire on force IKECOU a 0
  if (ikecou.ne.0.and.idtvar.lt.0) then
    write(nfecra,2210)ikecou
    iok = iok + 1
  endif

  if (igrhok.ne.0.and.igrhok.ne.1) then
    write(nfecra,2201)'IGRHOK',igrhok
    iok = iok + 1
  endif
  if (igrake.ne.0.and.igrake.ne.1) then
    write(nfecra,2201)'IGRAKE',igrake
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

  !     Si RELAXV(IK) a ete modifie par l'utilisateur mais que IKECOU n'est
  !     pas egal a 0, on previent que ce sera sans effet
  !     Sinon on verifie que RELAXV(IK) est bien compris entre 0 et 1
  !     (en stationnaire cela a deja ete fait plus haut)
  if (itytur.eq.6) then
    if ( (abs(relaxv(ik)+999.d0).gt.epzero .or.                 &
          abs(relaxv(iomg)+999.d0).gt.epzero ) .and.            &
         ikecou.ne.0) write(nfecra,2623)                        &
         relaxv(ik),relaxv(iomg)
    if (ikecou.eq.0 .and. idtvar.ge.0) then
      if (relaxv(ik).gt.1.d0.or.relaxv(ik).lt.0.d0 .or.         &
           relaxv(iomg).gt.1.d0.or.relaxv(iomg).lt.0.d0) then
        write(nfecra,2624) relaxv(ik),relaxv(iomg)
        iok = iok + 1
      endif
    endif
  else
    if ( (abs(relaxv(ik)+999.d0).gt.epzero .or.                 &
          abs(relaxv(iep)+999.d0).gt.epzero ) .and.             &
         ikecou.ne.0) write(nfecra,2623)                        &
         relaxv(ik),relaxv(iep)
    if (ikecou.eq.0 .and. idtvar.ge.0) then
      if (relaxv(ik).gt.1.d0.or.relaxv(ik).lt.0.d0 .or.         &
           relaxv(iep).gt.1.d0.or.relaxv(iep).lt.0.d0) then
        write(nfecra,2624) relaxv(ik),relaxv(iep)
        iok = iok + 1
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
  if (irijnu.ne.0.and.irijnu.ne.1) then
    write(nfecra,2201)'IRIJNU',irijnu
    iok = iok + 1
  endif
  if (irijrb.ne.0.and.irijrb.ne.1) then
    write(nfecra,2201)'IRIJRB',irijrb
    iok = iok + 1
  endif
  if (iturb.eq.30) then
    !     echo de paroi et implicitation speciale de la diffusion de epsilon
    !     seulement en Rij standard
    if (irijec.ne.0.and.irijec.ne.1) then
      write(nfecra,2201)'IRIJEC',irijec
      iok = iok + 1
    endif
    if (idifre.ne.0.and.idifre.ne.1) then
      write(nfecra,2201)'IDIFRE',idifre
      iok = iok + 1
    endif
  endif
  if (igrari.ne.0.and.igrari.ne.1) then
    write(nfecra,2201)'IGRARI',igrari
    iok = iok + 1
  endif
  !        IF ( IGRARI.EQ.1.AND.(GX**2+GY**2+GZ**2).LE.EPZERO**2 ) THEN
  !          WRITE(NFECRA,2620)'IGRARI',IGRARI,GX,GY,GZ
  !          IOK = IOK + 1
  !        ENDIF
  if (iclsyr.ne.0.and.iclsyr.ne.1) then
    write(nfecra,2201)'ICLSYR',iclsyr
    iok = iok + 1
  endif
  if (iclptr.ne.0.and.iclptr.ne.1) then
    write(nfecra,2201)'ICLPTR',iclptr
    iok = iok + 1
  endif
  if (nscal.gt.0) then
    if (iscalt.le.0.and.(gx**2+gy**2+gz**2).ge.epzero**2) then
      write(nfecra,2621)gx,gy,gz,iscalt
      if (igrari.eq.1) then
        write(nfecra,2622)'IGRARI',igrari
      endif
    endif
  endif
endif

!     Specifique LES

if (itytur.eq.4) then
  if (idries.ne.1.and.idries.ne.0) then
    write(nfecra,2201)'IDRIES',idries
    iok = iok + 1
  endif
  if (idries.ne.0.and.(iturb.eq.41.or.iturb.eq.42)) then
    write(nfecra,2630) idries,iturb
    iok = iok + 1
  endif
  !         La reduction du voisinage etendu peut degrader
  !         les resultats du modele dynamique en LES
  if (     iturb.eq.41                                                       &
      .and.(imrgrl.eq.3.or.imrgrl.eq.6)) then
    write(nfecra,2607) iturb, imrgra
  endif
endif

! --- Stokes


if (iprco .ne.0.and.iprco .ne.1) then
  write(nfecra,2200) 'IPRCO ',iprco
  iok = iok + 1
endif
if (iprco.eq.1) then
  if (irevmc.ne.0.and.irevmc.ne.1.and.irevmc.ne.2) then
    write(nfecra,2211) 'IREVMC',irevmc
    iok = iok + 1
  endif
  arakfr = arak
  if (idtvar.lt.0) arakfr=arakfr*relaxv(iu)
  if (arakfr.gt.1.d0 .or. arakfr.lt.0.d0) then
    write(nfecra,2640) 'ARAK  ',arakfr
    iok = iok + 1
  endif
  if (relaxv(ipr).gt.1d0 .or. relaxv(ipr).lt.0d0) then
    write(nfecra,2625) relaxv(ipr)
    iok = iok + 1
  endif
endif

! --- Couplage U-P

if (ipucou.ne.0.and.ipucou.ne.1) then
  write(nfecra,2200) 'IPUCOU ',ipucou
  iok = iok + 1
endif

! Incompatibilite pour le moment des matrices poids avec un theta schema
! Si theta n'est pas egal a 1 pour la vitesse (incompatibilite du
! pas de temps variable aussi)

jj = iu
if ((abs(thetav(jj)-1.0d0).gt.epzero).and.                       &
     ((idtvar.ne.0).or.(ipucou.eq.1))) then
  write(nfecra,2204) thetav(jj),idtvar,ipucou
endif

! --- Champ de vitesse fige

if (iccvfg.ne.0.and.iccvfg.ne.1) then
  write(nfecra,2200) 'ICCVFG ',iccvfg
  iok = iok + 1
endif

! --- Interpolation face des viscosites

if (imvisf.ne.0.and.imvisf.ne.1) then
  write(nfecra,2200) 'IMVISF',imvisf
  iok = iok + 1
endif

! --- Traitement de la temperature pour couplage SYRTHES
!     Verification du nombre de couplages

if (itbrrb.ne.0 .and. itbrrb.ne.1       ) then
  write(nfecra,2200) 'ITBRRB',itbrrb
  iok = iok + 1
endif

!     On regarde si ICPSYR a des valeurs realistes
if (nscal.gt.0) then
  do iscal = 1, nscal
    if (icpsyr(iscal).ne.0.and.icpsyr(iscal).ne.1) then
      call field_get_label(ivarfl(isca(iscal)), chaine)
      write(nfecra,2650)chaine(1:16),'ICPSYR',iscal,icpsyr(iscal)
      iok = iok + 1
    endif
  enddo
endif

!     On compte le nombre de scalaires couples
nbsccp = 0
if (nscal.gt.0) then
  do iscal = 1, nscal
    nbsccp = nbsccp+icpsyr(iscal)
  enddo
endif

! On regarde s'il y a du couplage

call nbcsyr (nbccou)

! S'il n'y a pas de couplage
if (nbccou.eq.0) then

  !  et qu'il n'y a pas zero scalaire couple, on s'arrete
  if (nbsccp.ne.0) then
    write(nfecra,2660)nbsccp,nscal
    iok = iok + 1
  endif

  ! Sinon, s'il y a du couplage
else

  ! et qu'il n'y a pas un et un seul scalaire couple, on s'arrete
  if (nbsccp.ne.1) then
    write(nfecra,2661)nscal,nbsccp
    iok = iok + 1
  endif

  ! que le scalaire couple n'est pas la temperature, on s'arrete
  ! attention : tout est pret, mais il n'y a pas eu de valid.
  ! en outre, ca permet de bien verifier que l'utilisateur
  ! ne s'est pas trompe dans 99% des cas d'utilisation
  ! (en compressible, on couple l'energie)
  do iscal = 1, nscal
    if (icpsyr(iscal).eq.1) then
      if (ippmod(icompf).lt.0) then
        if (abs(iscacp(iscal)).ne.1) then
          write(nfecra,2662)iscal,iscal,iscacp(iscal)
          iok = iok + 1
        endif
      else
        if (iscal.eq.iscalt .and. itherm.ne.3) then
          write(nfecra,2663)iscal,iscalt
          iok = iok + 1
        endif
      endif
    endif
  enddo

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
  if (nitmay.lt.1) then
    write(nfecra,3100) 'NITMAY',nitmay
    iok = iok + 1
  endif
  if (imligy.gt.1) then
    write(nfecra,2750) 'IMLIGY',imligy
    iok = iok + 1
  endif
  if (ircfly.ne.1.and.ircfly.ne.0) then
    write(nfecra,2200) 'IRCFLY',ircfly
    iok = iok + 1
  endif
  if (ischcy.ne.1.and.ischcy.ne.0) then
    write(nfecra,2200) 'ISCHCY',ischcy
    iok = iok + 1
  endif
  if (isstpy.ne.1.and.isstpy.ne.0) then
    write(nfecra,2200) 'ISSTPY',isstpy
    iok = iok + 1
  endif
  if (ntcmxy.lt.1) then
    write(nfecra,3100) 'NTCMXY',ntcmxy
    iok = iok + 1
  endif

  if (blency.gt.1.d0.or.blency.lt.0.d0) then
    write(nfecra,2710) 'BLENCY',blency
    iok = iok + 1
  endif
  if (climgy.lt.1.d0) then
    write(nfecra,2720) 'CLIMGY',climgy
    iok = iok + 1
  endif
  if (abs(extray-1.d0).gt.epzero.and.abs(extray).gt.epzero) then
    write(nfecra,2730) 'EXTRAY',extray
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

! --- Dynamic relaxation option

do ivar = 1, nvar
  if (iswdyn(ivar).ge.1) then
    ! The number of reconstruction sweeps is set to 20 at least
    if (nswrsm(ivar).lt.20) then
      nswrsm(ivar) = 20
      write(nfecra,2742) nswrsm(ivar)
    endif
    if (ivar.eq.iu.or.ivar.eq.iv.or.ivar.eq.iw) then
      if (isstpc(iu).eq.0.or.isstpc(iv).eq.0.or.isstpc(iw).eq.0) then
        isstpc(iu) = 1
        isstpc(iv) = 1
        isstpc(iw) = 1

        write(nfecra, 2743)
      endif
    endif
  endif
enddo

do ivar = 1, nvar
  if (nswrsm(ivar).le.0) then
    call field_get_label(ivarfl(ivar), chaine)
    write(nfecra,2747) chaine(1:16), nswrsm(ivar), 1
    nswrsm(ivar) = 1
  endif
enddo

!===============================================================================
! 3. TABLEAUX DE cstphy : formats 4000
!===============================================================================

! --- Constantes physiques de chaque phase

if (ro0.lt.0d0) then
  write(nfecra,2511) 'RO0   ', ro0
  iok = iok + 1
endif
if (viscl0.lt.0d0) then
  write(nfecra,2511) 'VISCL0', viscl0
  iok = iok + 1
endif

! --- Turbulence

!    On a besoin de UREF si on initialise la turbulence (debut de calcul
!      en turbulence ou suite laminaire->turbulent)
!    Ici on met juste un avertissement, car sans UREF l'utilisateur peut
!      initialiser la turbulence a la main. Un test complementaire sera fait
!      dans inivar.
if (itytur.eq.2.or.itytur.eq.3                    &
     .or.itytur.eq.5.or.iturb.eq.60               &
     .or.iturb.eq.70) then
  if (uref  .lt.0.d0) then
    write(nfecra,4100) uref
  endif
endif

if (iturb.eq.10) then
  if (xlomlg.le.0.d0) then
    write(nfecra,2511)'XLOMLG', xlomlg
    iok = iok + 1
  endif
endif

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

!     Scalaire associe dans le cas des variances
  do ii = 1, nscal
    if (iscavr(ii).gt.nscal.or.iscavr(ii).lt.0) then
      call field_get_label(ivarfl(isca(ii)), chaine)
      write(nfecra,4320)chaine(1:16),ii,nscal,iscavr(ii)
      iok = iok + 1
    endif
  enddo

!     On verifie que le scalaire associe a une fluctuation n'est
!     pas une fluctuation
  do ii = 1, nscal
    if (iscavr(ii).gt.0) then
      if (iscavr(iscavr(ii)).gt.0) then
        call field_get_label(ivarfl(isca(ii)), chaine)
        call field_get_label(ivarfl(isca(iscavr(ii))), chain2)
        write(nfecra,4321)chaine(1:16),chain2(1:16),ii,iscavr(ii),  &
             iscavr(ii),iscavr(iscavr(ii))
        iok = iok + 1
      endif
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

!     Valeur de la diffusivite positive (si cste)
!        si pas cste, on verifiera apres usphyv
  do ii = 1, nscal
    call field_get_key_int (ivarfl(isca(ii)), kivisl, ifcvsl)
    if (ifcvsl.lt.0.and.visls0(ii).lt.0d0) then
      call field_get_label(ivarfl(isca(ii)), chaine)
      write(nfecra,4340)chaine(1:16),ii,ifcvsl,visls0(ii)
      iok = iok + 1
    endif
  enddo

!     Valeur du sigma positif
  do ii = 1, nscal
    if (sigmas(ii).le.0d0) then
      call field_get_label(ivarfl(isca(ii)), chaine)
      write(nfecra,4350)chaine(1:16),ii,sigmas(ii)
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
  if (icp.eq.0) then
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
! 8. Rotating frame and unsteady rotor/stator coupling: 9000 formats
!===============================================================================

if (icorio.ne.0 .and. icorio.ne.1) then
  write(nfecra,9000) icorio
  iok = iok + 1
endif

if (imobil.eq.1 .or. iturbo.eq.2) then
  ! Unsteady rotor/stator coupling is not compatible with the
  !   steady algorithm...
  if (idtvar.lt.0) then
    write(nfecra,9010) idtvar
    iok = iok + 1
  endif
  ! ... nor with the time/space variable time steps
  if (idtvar.eq.1.or.idtvar.eq.2) then
    write(nfecra,9011) idtvar
    iok = iok + 1
  endif
endif

! Unsteady rotor/stator coupling is not compatible with the Lagrangian module
if (iturbo.eq.2.and.iilagr.ne.0) then
    write(nfecra,9012)
    iok = iok + 1
endif

if (icorio.ne.0 .and. iturbo.ne.0) then
  write(nfecra,9020) icorio
  iok = iok + 1
endif

!===============================================================================
! 10. Cavitation modelling: 9100 formats
!===============================================================================

if (icavit.ge.0) then
  ! For now, cavitation model is not compatible with the handling of
  ! hydrostatic pressure
  if (iphydr.ne.0) then
    write(nfecra,9110) iphydr
    iok = iok + 1
  endif
  ! Cavitation model is not compatible with dilatable or low-mach algorithms
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
 1230 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    NCAPT  DOIT ETRE UN ENTIER INFERIEUR OU EGAL A', i10,     /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  NCAPT  est le nombre de sondes utilisees pour les',         /,&
'@    historiques.',                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 1240 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    IHISVR(',i10,   ',1) DOIT ETRE UN ENTIER EGAL A -1',      /,&
'@      POSITIF OU NUL ET INFERIEUR OU EGAL A NCAPT = ',i10,    /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  IHISVR(I,1) indique les sondes a utiliser pour la variable',/,&
'@    I (-1 signifiant que toutes sont utilisees)',             /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 1250 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    IHISVR(',i10,   ',',i10,   ') DOIT ETRE UN ENTIER',       /,&
'@      STRICTEMENT POSITIF ET',                                /,&
'@      INFERIEUR OU EGAL A NCAPT = ', i10,                     /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute',                            /,&
'@',                                                            /,&
'@  IHISVR(I,j+1) indique le numero de la jieme sonde a',       /,&
'@    utiliser pour la variable a post-traiter numero I',       /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 1260 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    CHAMP', a16,                                              /,&
'@      MOT CLE ''log'' DOIT ETRE UN ENTIER EGAL A 0 OU 1',    /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute',                            /,&
'@',                                                            /,&
'@  ILISVR(I) indique si la variable I sera suivie lors des',   /,&
'@    impressions dans le listing',                             /,&
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
 2001 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE UN ENTIER POSITIF',                      /,&
'@      ET INFERIEUR OU EGAL A', i10,                           /,&
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
 2050 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    ITHERM DOIT ETRE UN ENTIER EGAL A 0, 1, 2 or 3',      /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2100 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@     PARAMETRES DU SCHEMA NUMERIQUE POUR LA VARIABLE', a16,   /,&
'@',                                                            /,&
'@ Parametre               ISTAT     ICONV',                    /,&
'@ Valeurs acceptees     0 ou  1   0 ou  1',                    /,&
'@ Valeurs entrees ici',i10,     i10,                           /,&
'@',                                                            /,&
'@ Parametre                         IDIFF     IDIFFT',         /,&
'@ Valeurs acceptees               0 ou  1   0 ou  1',          /,&
'@ Valeurs entrees ici',10X,     i10,      i10,                 /,&
'@',                                                            /,&
'@ Parametre               THETAV   BLENCV    ISCHCV    ISSTPC',/,&
'@ Valeurs acceptees     [0.; 1.] [0.; 1.]   0,1,2,3   0,1,2,3',/,&
'@ Valeurs entrees ici',      e9.2,     e9.2,i10,  i10,         /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2111 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@   INCOMPATIBILITE POUR LE SCHEMA EN TEMPS',                  /,&
'@',                                                            /,&
'@   Schema en temps pour la vitesse phase',                    /,&
'@      THETA n''a pas la meme valeur pour les 3 composantes',  /,&
'@',                                                            /,&
'@ Parametre THETAV              U          V          W',      /,&
'@ Valeurs entrees ici', e10.2,e10.2,e10.2,                     /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2112 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@   DONNEES NON ADMISSIBLES POUR LE SCHEMA EN TEMPS',          /,&
'@',                                                            /,&
'@  LE PARAMETRE THETAV POUR LA PRESSION DOIT VALOIR 1',        /,&
'@',                                                            /,&
'@  Il vaut ici', e14.5,                                        /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2121 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES',                /,&
'@    =========',                                               /,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS',                    /,&
'@',                                                            /,&
'@   EN L.E.S.',                                                /,&
'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE THETAV DU SCHEMA', /,&
'@     EN TEMPS DE LA VARIABLE', a16, ' EST 0.5',               /,&
'@     THETAV A ETE IMPOSE ICI A', e14.5,                       /,&
'@',                                                            /,&
'@  Le calcul sera execute',                                    /,&
'@',                                                            /,&
'@  Il est conseille de verifier les parametres donnes via',    /,&
'@  l''interface ou cs_user_parameters.f90.',                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2122 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES',                /,&
'@    =========',                                               /,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS',                    /,&
'@',                                                            /,&
'@   EN L.E.S.',                                                /,&
'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE BLENCV DU SCHEMA', /,&
'@     CONVECTIF DE LA VARIABLE', a16, ' EST 1.0',              /,&
'@     BLENCV A ETE IMPOSE ICI A', e14.5,                       /,&
'@',                                                            /,&
'@  Le calcul sera execute',                                    /,&
'@',                                                            /,&
'@  Il est conseille de verifier les parametres donnes via',    /,&
'@  l''interface ou cs_user_parameters.f90.',                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2123 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES',                /,&
'@    =========',                                               /,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS',                    /,&
'@',                                                            /,&
'@   EN L.E.S.',                                                /,&
'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE ISSTPC DU SCHEMA', /,&
'@     CONVECTIF DE LA VARIABLE', a16, ' EST 1',                /,&
'@     ISSTPC A ETE IMPOSE ICI A', i10,                         /,&
'@',                                                            /,&
'@  Le calcul sera execute',                                    /,&
'@',                                                            /,&
'@  Il est conseille de verifier les parametres donnes via',    /,&
'@  l''interface ou cs_user_parameters.f90.',                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2124 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES',                /,&
'@    =========',                                               /,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS',                    /,&
'@',                                                            /,&
'@   EN L.E.S.',                                                /,&
'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE ISSTPC DU SCHEMA', /,&
'@     CONVECTIF DE LA VARIABLE', a16, ' EST 0',                /,&
'@     ISSTPC A ETE IMPOSE ICI A', i10,                         /,&
'@',                                                            /,&
'@  Le calcul sera execute',                                    /,&
'@',                                                            /,&
'@  Il est conseille de verifier les parametres donnes via',    /,&
'@  l''interface ou cs_user_parameters.f90.',                   /,&
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
 2127 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS',                    /,&
'@',                                                            /,&
'@   EN L.E.S.',                                                /,&
'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE BLENCV DU SCHEMA', /,&
'@     CONVECTIF DE LA VARIABLE', a16, ' EST 1.0',              /,&
'@     BLENCV A ETE IMPOSE ICI A', e14.5,                       /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute',                             /,&
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
'@      U,V,W : THETA = ', 3e12.4,                              /,&
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
'@    cs_user_parameters.f90, et uslag1.',                      /,&
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
 2149 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@     ALGORITHME STATIONNAIRE',                                /,&
'@     COEFFICIENT DE RELAXATION POUR LA VARIABLE', a16,        /,&
'@',                                                            /,&
'@ RELAXV doit etre un reel compris entre 0 et 1',              /,&
'@ Il vaut ici', e14.5,                                         /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2150  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@   INCOMPATIBILITE ALGORITHME STATIONNAIRE',                  /,&
'@',                                                            /,&
'@   Coefficient de relaxation vitesse phase', i10,             /,&
'@      RELAXV n''a pas la meme valeur pour les 3 composantes', /,&
'@',                                                            /,&
'@ Parametre RELAXV              U          V          W',      /,&
'@ Valeurs entrees ici', e10.2,e10.2,e10.2,                     /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2151  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@   L''ALGORITHME STATIONNAIRE N''EST PAS COMPATIBLE AVEC',    /,&
'@    LE MODULE LAGRANGIEN QUI EST UNE APPROCHE INSTATIONNAIRE',/,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  L''indicateur IILAGR a ete positionne a', i10,              /,&
'@    dans uslag1 (module lagrangien active pour IILAGR>0).',   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2152  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@   L''ALGORITHME STATIONNAIRE N''EST PAS COMPATIBLE AVEC',    /,&
'@    LA L.E.S. QUI EST UNE MODELISATION INSTATIONNAIRE DE',    /,&
'@    LA TURBULENCE',                                           /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  L''indicateur ITURB a ete positionne a', i10,               /,&
'@    dans l''interface ou usipph.',                            /,&
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
 2201 format(                                                     &
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
 2204 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@   DONNEES NON ADMISSIBLES POUR LE SCHEMA EN TEMPS',          /,&
'@',                                                            /,&
'@  ON DEMANDE LA PRISE UN SCHEMA EN TEMPS POUR LA VITESSE',    /,&
'@    D''ORDRE 2 AVEC UN PAS DE TEMPS NON CONSTANT, UN',        /,&
'@    ALGORITHME STATIONNAIRE OU LES MATRICES POIDS',           /,&
'@',                                                            /,&
'@  THETAV VAUT ICI', e14.5,' POUR LA VITESSE',                 /,&
'@  ALORS QUE IDTVAR VAUT', i10,                                /,&
'@  ET IPUCOU VAUT',        i10,                                /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2205 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE UN ENTIER COMPRIS ENTRE -6 et 6',        /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2206 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    ANOMAX DOIT ETRE UN REEL POSITIF OU NUL ET',              /,&
'@                             INFERIEUR OU EGAL A PI/2',       /,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  On demande la reconstruction des gradients par moindres',   /,&
'@    carres sur voisinage etendu reduit (IMRGRA = ', i10,  ').',/,&
'@    Le critere est base sur l''angle de non orthogonalite',   /,&
'@    des faces ANOMAX qui doit etre fourni en radians et',     /,&
'@    compris dans les bornes indiquees ci-dessus.',            /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2207 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    L''UTILISATION DE LA METHODE DE CALCUL DE GRADIENT PAR',  /,&
'@      MOINDRES CARRES EST IMPOSSIBLE AVEC',                   /,&
'@        EXTRAG(IPR) DIFFERENT DE 0 ET 1',                     /,&
'@',                                                            /,&
'@    ON A ICI',                                                /,&
'@        IMRGRA         = ', i10,                              /,&
'@        EXTRAG(IPR) = ', e14.5,                               /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2208 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    EN K-EPS PROD LIN (ITURB=21) ET EN V2F (ITURB=50/51)',    /,&
'@    IKECOU DOIT ETRE EGAL A 0',                               /,&
'@    ITURB  VAUT ICI', i10,                                    /,&
'@    IKECOU VAUT ICI', i10,                                    /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2209 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    LE MODELE DE PAROI A DEUX ECHELLES (IWALLF=3, 4 OU 5)',   /,&
'@    EST INCOMPATIBLE AVEC UN CALCUL EN LAMINAIRE, EN',        /,&
'@    LONGUEUR DE MELANGE, EN SPALART-ALLMARAS OU EN L.E.S.',   /,&
'@    ON A ICI ITURB=',i10,                                     /,&
'@         ET IWALLF=',i10,                                     /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2210 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    L''ALGORITHME STATIONNAIRE EST INCOMPATIBLE AVEC LE',     /,&
'@    COUPLAGE DES TERMES SOURCES EN K-EPS, V2F OU K-OMEGA',    /,&
'@',                                                            /,&
'@    ON A ICI IKECOU=',i10,                                    /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2211 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE UN ENTIER EGAL A 0, 1 OU 2',             /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2300 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    IMLIGR(',i10,   ') DOIT ETRE UN ENTIER',                  /,&
'@      INFERIEUR OU EGAL A 1',                                 /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  IMLIGR(I) indique le mode de limitation des gradients',     /,&
'@    pour la variable I',                                      /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2310 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    IRCFLU(',i10,   ') DOIT ETRE UN ENTIER EGAL A 0 OU 1',    /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  IRCFLU(I) indique si les flux sont reconstruits',           /,&
'@    pour la variable I',                                      /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2311 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    IRCFLU(',i10,   ') = ', i10,   ' EST INCOMPATIBLE AVEC',  /,&
'@    ISCHCV(',i10,   ') = ', i10,                              /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  IRCFLU(I) = 0 indique que les flux ne sont pas',            /,&
'@    reconstruits pour la variable I.',                        /,&
'@  ISCHCV(I) = 0 (schema SOLU) demande une reconstruction',    /,&
'@    pour les flux convectifs.',                               /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2320 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    CLIMGR(',i10,   ') DOIT ETRE UN REEL',                    /,&
'@      SUPERIEUR OU EGAL A 1',                                 /,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  CLIMGR(I) est le coefficient de limitation des gradients',  /,&
'@    pour la variable I',                                      /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2330 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    EXTRAG(',i10,   ') DOIT ETRE UN REEL EGAL A 0 OU 1',      /,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2331 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    EXTRAG(',i10,   ') DOIT ETRE NUL',                        /,&
'@      (Des valeurs non nulles sont autorisees pour la',       /,&
'@       pression uniquement)',                                 /,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
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
 2500 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE UN ENTIER EGAL A -1, 0, 1 OU 2',         /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2510 format(                                                     &
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
 2520 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@     DTMIN DOIT ETRE SUPERIEUR OU EGAL A 0. ET',              /,&
'@                     INFERIEUR OU EGAL A DTMAX',              /,&
'@      ICI DTMIN = ', e14.5,     ' ET DTMAX = ', e14.5,        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Quand le pas de temps n est pas uniforme et constant,',     /,&
'@    les reels DTMIN et DTMAX bornent ses variations.',        /,&
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
 2540 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@  ON DEMANDE UNE LIMITATION DU PAS DE TEMPS LIEE AUX EFFETS', /,&
'@    DE DENSITE (IPTLRO = ', i10,   ') AVEC UNE OPTION DE',    /,&
'@    PAS DE TEMPS FIXE (IDTVAR = ', i10,   ')',                /,&
'@',                                                            /,&
'@  Le calcul sera engage, mais le pas de temps ne sera pas',   /,&
'@    clippe. Le code indiquera juste les eventuels',           /,&
'@    depassements de critere local.',                          /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2541 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@  ON DEMANDE UNE LIMITATION DU PAS DE TEMPS LIEE AUX EFFETS', /,&
'@    DE DENSITE (IPTLRO = ', i10,   ') AVEC UN ALGORITHME',    /,&
'@    STATIONNAIRE (IDTVAR = ', i10,   ')',                     /,&
'@',                                                            /,&
'@  Le calcul sera engage, l''option IPTLRO ignoree.',          /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2600 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE UN ENTIER EGAL A 0, 10, 20, 21, 30, 31,',/,&
'@    40, 41, 42, 50, 51 OU 60',                                /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2601 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    ========='                                               ,/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0 OU 1'                ,/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@'                                                            ,/,&
'@  Verifier les parametres donnes via l''interface'           ,/,&
'@    ou cs_user_parameters.f90.'                              ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2602 format( &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    ========='                                               ,/,&
'@    L''OPTION IRCCOR = 1 N''EST COMPATIBLE QU''AVEC'         ,/,&
'@    L''OPTION ITURB = 20, 21, 50, 51, 60 ou 70'              ,/,&
'@    ITURB VAUT ICI ',I10                                     ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@'                                                            ,/,&
'@  Verifier les parametres donnes via l''interface'           ,/,&
'@    ou cs_user_parameters.f90.'                              ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2603 format( &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    ========='                                               ,/,&
'@    LE MODELE DE TURBULENCE K-OMEGA SST N''EST PAS'          ,/,&
'@     COMPATIBLE AVEC LE LAGRANGIEN EN COUPLAGE INVERSE'      ,/,&
'@'                                                            ,/,&
'@  Le calcul ne peut etre execute.'                           ,/,&
'@'                                                            ,/,&
'@  Verifier les parametres donnes via l''interface'           ,/,&
'@    ou cs_user_parameters.f90.'                              ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
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
 2606 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    LA METHODE DES VORTEX NE CONCERNE QUE LES CALCULS EN LES',/,&
'@    ON A ICI',                                                /,&
'@    ITURB   = ',i10,                                          /,&
'@    IVRETX  = ',i10,                                          /,&
'@',                                                            /,&
'@  Le calcul sera execute en ignorant le mot clef IVRTEX.',    /,&
'@    (il est repositione a 0)',                                /,&
'@',                                                            /,&
'@  Il est conseille de verifier les parametres donnes via',    /,&
'@  l''interface ou cs_user_parameters.f90.',                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2607 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    ON DEMANDE LA REDUCTION DU VOISINAGE ETENDU POUR LE',     /,&
'@    CALCUL DES GRADIENTS PAR MOINDRE CARRES, POUR UN',        /,&
'@    CALCUL EN L.E.S AVEC LE MODELE DYNAMIQUE.',               /,&
'@    MODELE DYNAMIQUE.',                                       /,&
'@      ITURB = ',i10,                                          /,&
'@      IMRGRA= ',i10,                                          /,&
'@',                                                            /,&
'@  Le calcul sera engage.',                                    /,&
'@',                                                            /,&
'@  Le calcul de la moyenne locale de la constante de',         /,&
'@    Smagorinsky dynamique peut etre degrade.',                /,&
'@',                                                            /,&
'@  Il est conseille :',                                        /,&
'@    - soit d''utiliser une methode de calculs des gradients', /,&
'@      par moindre carres sur voisinage etendu complet',       /,&
'@      (IMRGRA = 2)',                                          /,&
'@    - soit de calculer sa propre moyenne de la constante',    /,&
'@      dynamique dans la subroutine USSMAG.',                  /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
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
 2623 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',                                                            /,&
'@ Le coefficient RELAXV des variables de la turbulence a ete', /,&
'@ modifie alors que IKECOU ne vaut pas 0. Il vaut ici',        /,&
'@ - pour k                  :', e12.4,                         /,&
'@ - pour epsilon (ou omega) :', e12.4,                         /,&
'@',                                                            /,&
'@ La modification sera sans effet (RELAXV n''est utile que',   /,&
'@ si IKECOU=0)',                                               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2624 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',                                                            /,&
'@ Le coefficient RELAXV des variables de la turbulence doit',  /,&
'@ etre un reel compris entre 0 et 1. Il vaut ici :',           /,&
'@ - pour k                  :', e12.4,                         /,&
'@ - pour epsilon (ou omega) :', e12.4,                         /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2625 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',                                                            /,&
'@ Le coefficient RELAXV de la pression doit etre un reel',     /,&
'@ compris entre 0 et 1. Il vaut ici :', e12.4,                 /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2630 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',                                                            /,&
'@  On demande la prise en compte de l''amortissement de',      /,&
'@    Van Driest (IDRIES = ',  i10,')',                         /,&
'@    avec un modele LES incompatible (ITURB = ', i10,')',      /,&
'@    (modele dynamique et modele WALE)',                       /,&
'@',                                                            /,&
'@  Le calcul ne sera pas realise.',                            /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2640 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE UN REEL INCLUS DANS L''INTERVALLE [0;1]',/,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2650 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    SCALAIRE ', a16,                                          /,&
'@',    a6,'(',i10,   ') DOIT ETRE UN ENTIER EGAL A 0 OU 1',    /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  ICPSYR(I) est l indicateur de couplage du scalaire I avec', /,&
'@    SYRTHES.',                                                /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2660 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    INCOHERENCE POUR LE COUPLAGE SYRTHES',                    /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Aucun couplage avec SYRTHES n''a ete defini.',              /,&
'@  Le nombre de scalaires couples est cependant', i10,         /,&
'@    (Le nombre de scalaires total est ici', i10,   ')',       /,&
'@  Verifier le couplage SYRTHES-Noyau.',                       /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2661 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    INCOHERENCE POUR LE COUPLAGE SYRTHES',                    /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Un couplage avec SYRTHES a ete defini.',                    /,&
'@  Le couplage avec SYRTHES necessite de disposer d''un',      /,&
'@    scalaire couple (et d''un seul).',                        /,&
'@  Le nombre de scalaires total   est ici', i10,               /,&
'@  Le nombre de scalaires couples est ici', i10,               /,&
'@  Verifier le couplage SYRTHES-Noyau.',                       /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2662 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    INCOHERENCE POUR LE COUPLAGE SYRTHES',                    /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Un couplage avec SYRTHES a ete defini.',                    /,&
'@  Si le code est couple a SYRTHES, le scalaire couple doit',  /,&
'@    etre la temperature.',                                    /,&
'@  Le scalaire couple est ici le scalaire ', i10,              /,&
'@    Ce n''est pas une temperature car',                       /,&
'@                    ISCCACP(',i10,   ') = ', i10,             /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@  Pour coupler un scalaire qui n''est pas la temperature,',   /,&
'@    contacter l''equipe de developpement.',                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2663 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    INCOHERENCE POUR LE COUPLAGE SYRTHES',                    /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Un couplage avec SYRTHES a ete defini.',                    /,&
'@  En compressible, si le code est couple a SYRTHES, le',      /,&
'@    scalaire couple doit etre l''energie.',                   /,&
'@  Le scalaire couple est ici le scalaire ', i10,              /,&
'@    Ce n''est pas l''energie car',                            /,&
'@                    ISCALT = ', i10,                          /,&
'@',                                                            /,&
'@  Contacter l''equipe de developpement.',                     /,&
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
 2710 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE UN REEL INCLUS DANS LINTERVALLE [0.;1.]',/,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2720 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE UN REEL SUPERIEUR OU EGAL A 1.',         /,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2730 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE UN REEL EGAL A 0. OU A 1.',              /,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
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
 2742 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :',                                             /,&
'@    =========',                                               /,&
'@    OPTION iswdyn = 1  : ON SOUHAITE NSWRSM(IPR) > 19',       /,&
'@    NSWRSM(IPR) VAUT ICI', i10,                               /,&
'@',                                                            /,&
'@  Le calcul continue avec nswrsm(ipr) = 20.',                 /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2743 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :',                                             /,&
'@    =========',                                               /,&
'@    OPTION SWPDYN = 1  : LE TEST DE PENTE SUR U, V, W DOIT',  /,&
'@    ETRE DESACTIVE.',                                         /,&
'@',                                                            /,&
'@  Le calcul continue avec :',                                 /,&
'@  isstpc(iu) = isstpc(iu) = isstpc(iu) = 1.',                 /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2747 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES',                /,&
'@    =========',                                               /,&
'@   LA VALEUR POUR LE PARAMETRE NSWRSM POUR'                   /,&
'@     LA VARIABLE', a16, 'DOIT ETRE POSITIVE, VAUT ICI',  i10, /,&
'@     NSWRSM A ETE IMPOSE A', i10,                             /,&
'@',                                                            /,&
'@  Le calcul sera execute',                                    /,&
'@',                                                            /,&
'@  Il est conseille de verifier les parametres donnes via',    /,&
'@  l''interface ou cs_user_parameters.f90.',                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2750 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@',    a6,' DOIT ETRE UN ENTIER INFERIEUR OU EGAL A 1',        /,&
'@    IL VAUT ICI', i10,                                        /,&
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

 4100 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ENTREE DES DONNEES',                          /,&
'@    =========',                                               /,&
'@    LA VITESSE DE REFERENCE UREF N''A PAS ETE INITIALISEE',   /,&
'@    OU A ETE MAL INITIALISEE (VALEUR NEGATIVE).',             /,&
'@    ELLE VAUT ICI', e14.5,                                    /,&
'@',                                                            /,&
'@  Le calcul pourra etre execute si la turbulence est',        /,&
'@  initialisee a partir d''un fichier suite de calcul ou par', /,&
'@  l''interface ou par la routine cs_user_initialization',     /,&
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
 4320 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    SCALAIRE ', a16,                                          /,&
'@    iscavr(',i10,   ') DOIT ETRE UN ENTIER',                  /,&
'@      POSITIF OU NUL ET',                                     /,&
'@      INFERIEUR OU EGAL A NSCAL = ', i10,                     /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Si iscavr(I) est nul, le scalaire I n est pas une variance',/,&
'@  Si iscavr(I) est positif, le scalaire I est une variance :',/,&
'@    il s agit de la variance des fluctuations du scalaire J', /,&
'@    dont le numero est iscavr(I)',                            /,&
'@  Verifier les parametres donnes via l''interface',           /,&
'@    ou cs_user_parameters.f90.',                              /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4321 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    LE SCALAIRE ', a16, 'EST DEFINI COMME LA FLUCTUATION',    /,&
'@    DU SCALAIRE ', a16,                                       /,&
'@    (iscavr(',i10,   ') = ', i10,   '),',                     /,&
'@    QUI EST LUI-MEME DEFINI COMME UNE FLUCTUATION',           /,&
'@    (iscavr(',i10,   ') = ', i10,   ').',                     /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Si iscavr(I) est positif, le scalaire I est une variance :',/,&
'@    il s agit de la variance des fluctuations du scalaire J', /,&
'@    dont le numero est iscavr(I) et on a donc forcement',     /,&
'@    iscavr(J) = 0',                                           /,&
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
 4340 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    SCALAIRE ', a16,                                          /,&
'@    VISLS0(',i10,   ') DOIT ETRE UN REEL POSITIF',            /,&
'@      QUAND scalar_diffusivity_id = ', i10,                   /,&
'@    IL VAUT ICI', e14.5,                                      /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  VISLS0(I) est le coefficient de diffusion moleculaire du',  /,&
'@    scalaire et doit etre positif quand la valeur de son',    /,&
'@    mot cle scalar_diffusivity_id est inferieur a 0',         /,&
'@    (dans ce cas, un coefficient variable est donne',         /,&
'@    via l''interface ou usphyv).',                            /,&
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
 9000 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@    INDICATEUR DE PRISE EN COMPTE DE LA ROTATION'             /,&
'@',                                                            /,&
'@  ICORIO DOIT VALOIR 0 OU 1',                                 /,&
'@    IL VAUT ICI', i10,                                        /,&
'@',                                                            /,&
'@  Le calcul ne peut etre execute.',                           /,&
'@',                                                            /,&
'@  Verifier les parametres donnes via cs_user_parameter.f90'  ,/,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
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
 9012  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@   LE COUPLAGE ROTOR/STATOR INSTATIONNAIRE N''EST PAS',       /,&
'@     COMPATIBLE AVEC LE MODULE LAGRANGIEN',                   /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9020  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@   LE CALCUL EN REFERENTIEL RELATIF (ICORIO = 1) N''EST PAS', /,&
'@     COMPATIBLE AVEC LE MODULE TURBOMACHINE (ITURBO > 0)',    /,&
'@',                                                            /,&
'@  Si necessaire, utiliser uniquement le module turbomachine', /,&
'@    en definissant plusieurs rotors',                         /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9110  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES',               /,&
'@    =========',                                               /,&
'@   LE MODELE DE CAVITATION N''EST PAS COMPATIBLE AVEC LES',   /,&
'@     ALGORITHMES DE PRISE EN COMPTE DE LA PRESSION',          /,&
'@     HYDROSTATIQUE',                                          /,&
'@',                                                            /,&
'@  Le calcul ne sera pas execute.',                            /,&
'@',                                                            /,&
'@  L''indicateur IPHYDR a ete positionne a', i10,              /,&
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
'@   LE MODELE DE CAVITATION N''EST PAS COMPATIBLE AVEC LES',   /,&
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
 1230 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    NCAPT  MUST BE AN INTEGER LESS THAN or EGAL A', i10,     /, &
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The computation CANNOT  start.',                           /,&
'@',                                                            /,&
'@  NCAPT  is the number of probes for history/ time series',   /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 1240 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    IHISVR(',i10,   ',1) MUST BE AN INTEGER EGAL to -1',      /,&
'@      or Zero,    or positive but less than NCAPT = ',i10,    /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@    The calculation could NOT run.',                          /,&
'@',                                                            /,&
'@  IHISVR(I,1) is the number of probes for variable  I',       /,&
'@      (-1 means all probes are used)',                        /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 1250 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    IHISVR(',i10,   ',',i10,   ') MUST BE AN INTEGER',        /,&
'@    LARGER OR EQUAL TO 1 AND',                                /,&
'@      LESS OR EQUAL TO NCAPT = ', i10,                        /,&
'@   IT HAS VALUE = ',i10,                                      /,&
'@',                                                            /,&
'@  The calculation could NOT run.',                            /,&
'@',                                                            /,&
'@  IHISVR(I,j+1) gives the number of the j-ieth probe',        /,&
'@    to be used with variable number I to post-process',       /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 1260 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    FIELD', a16,                                              /,&
'@      KEY WORD ''log'' MUST BE AN INTEGER EQUAL  0 OR  1',    /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@  The calculation could NOT run.',                            /,&
'@',                                                            /,&
'@  ILISVR(I) tells if variable (I)should be included',         /,&
'@    in the printed listing',                                  /,&
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
 2001 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' MUST BE AN INTEGER POSITIVE ',                     /,&
'@      AND LESS THAN or EGAL TO',i10,                          /,&
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
 2050 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    ITHERM MUST BE AN INTEGER EGAL TO 0, 1, 2 or 3',      /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2100 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@     PARAMETERS OF NUMERICAL SCHEME FOR VARIABLE',     a16,   /,&
'@',                                                            /,&
'@ Parameter               ISTAT     ICONV',                    /,&
'@ Values   accepted     0 or  1   0 or  1',                    /,&
'@ Values entered here',i10,     i10,                           /,&
'@',                                                            /,&
'@ PARAMETER                         IDIFF     IDIFFT',         /,&
'@ Values   accepted               0 or  1   0 or  1',          /,&
'@ Values entered here',10X,     i10,      i10,                 /,&
'@',                                                            /,&
'@ PARAMETER               THETAV   BLENCV    ISCHCV    ISSTPC',/,&
'@ Values   accepted     [0.; 1.] [0.; 1.]   0,1,2,3   0,1,2,3',/,&
'@ Values entered here',      e9.2,     e9.2,i10,  i10,         /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2111 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@   INCOMPATIBILITY FOR TIME DISCRETISATION SCHEME',           /,&
'@',                                                            /,&
'@   TIME DISCRETISATION SCHEME for velocity',                  /,&
'@      THETA DOES NOT HAVE SAME VALUE FOR ALL 3 COMPONENTS',   /,&
'@',                                                            /,&
'@ Parameter THETAV              U          V          W',      /,&
'@ Values  entered here',e10.2,e10.2,e10.2,                     /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2112 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@  INCOMPATIBILITY FOR TIME DISCRETISATION SCHEME',            /,&
'@',                                                            /,&
'@  PARAMETER THETAV FOR PRESSURE MUST BE EQUAL TO    1',       /,&
'@',                                                            /,&
'@  It has value', e14.5,                                       /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2121 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING  :      WHEN READING INPUT DATA',                /&
'@    =========',                                               /,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME',                    /,&
'@',                                                            /,&
'@   IN L.E.S.',                                                /,&
'@   RECOMMENDED VALUE FOR PARAMETER THETAV IN TIME-SCHEME',    /,&
'@   FOR VARIABLE',               a16, ' IS 0.5',               /,&
'@     THETAV IS NOW IMPOSED AT',  e14.5,                       /,&
'@',                                                            /,&
'@  computation will go on',                                    /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2122 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHEN READING INPUT DATA',                /,&
'@    =========',                                               /,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME',                    /,&
'@',                                                            /,&
'@   WITH L.E.S.',                                              /,&
'@   THE VALUE RECOMMANDED FOR THE PARAMETER BLENCV',  /,   &
'@   FOR CONVECTION OF VARIABLE', a16, ' EST 1.0',              /,&
'@     BLENCV IS NOW IMPOSED AS', e14.5,                        /,&
'@',                                                            /,&
'@  computation will go on',                                    /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2123 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHEN READING INPUT DATA',                /,&
'@    =========',                                               /,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME',                    /,&
'@',                                                            /,&
'@   EN L.E.S.',                                                /,&
'@   THE VALUE RECOMMANDED FOR THE PARAMETER ISSTPC IN SCHEME', /,&
'@     CONVECTION OF VARIABLE', a16, ' IS    1',                /,&
'@     ISSTPC IS NOW IMPOSED AS',  i10,                         /,&
'@',                                                            /,&
'@  Computation will go on',                                    /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2124 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHEN READING INPUT DATA',                /,&
'@    =========',                                               /,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME',                    /,&
'@',                                                            /,&
'@   EN L.E.S.',                                                /,&
'@   THE VALUE RECOMMANDED FOR THE PARAMETER ISSTPC IN SCHEME', /,&
'@     CONVECTION OF VARIABLE',   a16, ' IS  0',                /,&
'@     ISSTPC IS NOW IMPOSED AS',  i10,                         /,&
'@',                                                            /,&
'@  Computation will go on',                                    /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
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
 2127 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME',                    /,&
'@',                                                            /,&
'@   EN L.E.S.',                                                /,&
'@   THE VALUE RECOMMANDED FOR THE PARAMETER BLENCV IN',        /,&
'@     CONVECTION OF VARIABLE',   a16, ' IS  1.0',              /,&
'@     BLENCV IS NOW IMPOSED AS',  e14.5,                       /,&
'@',                                                            /,&
'@  Computation will NOT proceed',                              /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
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
'@      U,V,W : THETA = ', 3e12.4,                              /,&
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
'@  and uslag1.',                                               /,&
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
 2149 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@     ALGORITHME STATIONNAIRE',                                /,&
'@     COEFFICIENT DE RELAXATION POUR LA VARIABLE', a16,        /,&
'@',                                                            /,&
'@ RELAXV MUST BE A  REAL comprised between    0 et 1',         /,&
'@ it has here a value of', e14.5,                              /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2150  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@   INCOMPATIBILITE ALGORITHME STATIONNAIRE',                  /,&
'@',                                                            /,&
'@   relaxation Coefficient of velocity  phase', i10,           /,&
'@      RELAXV should have the same value for all components',  /,&
'@',                                                            /,&
'@ PARAMETER RELAXV              U          V          W',      /,&
'@ Values entered here', e10.2,e10.2,e10.2,                     /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2151  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@   L''ALGORITHME STATIONNAIRE N''EST PAS COMPATIBLE AVEC',    /,&
'@    LE MODULE LAGRANGIEN QUI EST UNE APPROCHE INSTATIONNAIRE',/,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@   Integer parameter IILAGR was set to',    i10,              /,&
'@    in uslag1 (lagrangian module active for IILAGR>0).',      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2152  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@   STEADY STATE SOLVER SCHEME IS NOT COMPATIBLE WITH',        /,&
'@   L.E.S. TRUBULENCE MODELLING WHICH IS TIME-DEPENDENT',      /,&
'@    BY NATURE',                                               /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  Integer parameter ITURB was set to', i10,                   /,&
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
 2201 format(                                                     &
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
 2204 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@  INCOMPATIBILITY FOR TIME DISCRETISATION SCHEME',            /,&
'@',                                                            /,&
'@  ON DEMANDE LA PRISE UN SCHEMA EN TEMPS FOR    VELOCITY',    /,&
'@    D''ORDRE 2 AVEC UN PAS DE TEMPS NON CONSTANT, UN',        /,&
'@    ALGORITHME STATIONNAIRE or the MATRICES POIDS',           /,&
'@',                                                            /,&
'@  THETAV IS EQUAL', e14.5,' FOR    VELOCITY',                 /,&
'@  ALORS QUE IDTVAR VAUT', i10,                                /,&
'@  ET IPUCOU VAUT',        i10,                                /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2205 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' MUST BE AN INTEGER BETWEEN -6 and 6',              /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2206 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    ANOMAX MUST BE A REAL POSITIVE or NUL AND',               /,&
'@                             LESS THAN or EQUAL TO PI/2',     /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  On demande la reconstruction des gradients par moindres',   /,&
'@    carres sur voisinage etendu reduit (IMRGRA = ', i10,  ').',/,&
'@    Le critere est base sur l''angle de non orthogonalite',   /,&
'@    des faces ANOMAX qui doit etre fourni en radians et',     /,&
'@    compris dans the bornes indiquees ci-dessus.',            /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2207 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    L''UTILISATION DE LA METHODE DE CALCUL DE GRADIENT PAR',  /,&
'@      MOINDRES CARRES EST IMPOSSIBLE AVEC',                   /,&
'@        EXTRAG(IPR) DIFFERENT DE 0 ET 1',                     /,&
'@',                                                            /,&
'@    ON A ICI',                                                /,&
'@        IMRGRA      = ', i10,                                 /,&
'@        EXTRAG(IPR) = ', e14.5,                               /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2208 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    EN K-EPS PROD LIN (ITURB=21) ET EN V2F (ITURB=50/51)',    /,&
'@    IKECOU DOIT ETRE EGAL A 0',                               /,&
'@    ITURB  IS EQUAL', i10,                                    /,&
'@    IKECOU IS EQUAL', i10,                                    /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2209 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    2u*scales version of WALL FUNCTION (IWALLF=3, 4 or 5)',   /,&
'@    EST INCOMPATIBLE AVEC UN CALCUL EN LAMINAIRE, EN',        /,&
'@    LONGUEUR DE MELANGE or EN L.E.S.',                        /,&
'@    ON A ICI ITURB=',i10,                                     /,&
'@         ET IWALLF=',i10,                                     /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2210 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    SOLVE STEADY-STATE EQN. OPTION IS NOT COMPATIBLE WITH',   /,&
'@    COUPLING OF SOURCE TERMS IN K-EPS, V2F or K-OMEGA',       /, &
'@',                                                            /,&
'@    WE HAVE  IKECOU=',i10,                                    /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2211 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' MUST BE AN INTEGER EGAL A 0, 1 or 2',             /, &
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2300 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    IMLIGR(',i10,   ') MUST BE AN INTEGER',                  /, &
'@      LESS THAN or EGAL A 1',                                 /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  IMLIGR(I) indique le mode de limitation des gradients',     /,&
'@    pour la variable I',                                      /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2310 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    IRCFLU(',i10,   ') MUST BE AN INTEGER EQUAL  0  OR 1',    /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  IRCFLU(I) flags if fluxes are reconstructed for',            /&
'@             variable I',                                      /&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2311 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    IRCFLU(',i10,   ') = ', i10,  ' IS  INCOMPATIBLE WITH',   /,&
'@    ISCHCV(',i10,   ') = ', i10,                              /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  IRCFLU(I) = 0 fluxes are not reconstructed for variable I', /,&
'@',                            /,                          &
'@  ISCHCV(I) = 0 (schema SOLU) requests a  reconstruction',    /,&
'@   for convective fluxes',                                    /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2320 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    CLIMGR(',i10,   ') DOIT ETRE A  REAL',                    /,&
'@      SUPERIEUR or EGAL A 1',                                 /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  CLIMGR(I) is the coefficient limiting the gradients',       /,&
'@    for variable I',                                          /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2330 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    EXTRAG(',i10,   ')  MUST BE REAL and equal 0 or 1',       /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2331 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    VARIABLE', a16,                                           /,&
'@    EXTRAG(',i10,   ') MUST BE ZERO',                         /,&
'@      (non zero values are only possible for pressure',       /,&
'@',                                                            /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
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
 2500 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' MUST BE AN INTEGER equal to -1, 0, 1 or 2',        /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2510 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' must be a positive real number',                   /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
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
 2520 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@       DTMIN must be > or =  0.     and',                     /,&
'@              DTMIN must be  < or = DTMAX',                   /,&
'@      Here DTMIN = ', e14.5,     ' and DTMAX = ', e14.5,      /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  When the time-step is non uniforme and constant',           /,&
'@              DTMIN <= timestep <= DTMAX',                    /,&
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
 2540 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@  ON DEMANDE UNE LIMITATION DU PAS DE TEMPS LIEE AUX EFFETS', /,&
'@    DE DENSITE (IPTLRO = ', i10,   ') AVEC UNE OPTION DE',    /,&
'@    PAS DE TEMPS FIXE (IDTVAR = ', i10,   ')',                /,&
'@',                                                            /,&
'@  Le calcul sera engage, mais le pas de temps ne sera pas',   /,&
'@    clippe. Le code indiquera juste the eventuels',           /,&
'@    depassements de critere local.',                          /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2541 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@  A time-step reduction in liaison with variable density',    /,&
'@    effects (IPTLRO = ', i10,   ') was requested while',      /,&
'@   steady-state algorythm is selected (IDTVAR = ', i10,   ')',/,&
'@',                                                            /,&
'@  Computation will run, but option IPTLRO is ignored',        /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2600 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' MUST BE AN INTEGER EGAL A 0, 10, 20, 21, 30, 31,', /,&
'@    40, 41, 42, 50, 51 or 60',                                /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2601 format( &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    ',A6,' MUST BE AN INTEGER EGAL A 0 OR 1',                 /,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2602 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING:   STOP WHILE READING INPUT DATA',                /,&
'@    =========',                                               /,&
'@    IRCCOR = 1 OPTION IS ONLY COMPATIBLE WITH',               /,&
'@    ITURB = 20, 21, 50, 51, 60 OR 70',                        /,&
'@    ITURB HAS VALUE ',I10                                    ,/,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)
 2603 format( &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    THE K-OMEGA SST TURBULENCE MODEL IS NOT COMPATIBLE WITH', /,&
'@    TWO-WAY COUPLING IN LAGRANGIAN MODELLING',                /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
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
 2606 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    Synthetic Vortex method is only for for LES',             /,&
'@    here we have',                                            /,&
'@    ITURB   = ',i10,                                          /,&
'@    IVRETX  = ',i10,                                          /,&
'@',                                                            /,&
'@  computation will go on while ignoring keyword  IVRTEX.',    /,&
'@    (it is reset to       0)',                                /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2607 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@   A reduction of the extened neighborhood was selected',     /,&
'@   for the caculation of the gradients by least squares.',    /,&
'@   However this will also be applied to the averaging in',    /,&
'@    the LES Dynamic model (also selected)',                   /,&
'@      ITURB = ',i10,                                          /,&
'@      IMRGRA= ',i10,                                          /,&
'@',                                                            /,&
'@  Computation will run, but',                                 /,&
'@',                                                            /,&
'@  averaging of the Smagorinsky constant can be',              /,&
'@  degraded, as it uses the same reduced neighborhood.',       /,&
'@',                                                            /,&
'@  Recommendation',                                            /,&
'@    - use extended neighborhood',                             /,&
'@',                                                            /,&
'@      (IMRGRA = 2)',                                          /,&
'@    - user defines (yourself) the averaging of the dynamic',  /,&
'@   Smagorinsky constant via  subroutine USSMAG.',             /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
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
 2623 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',                                                            /,&
'@ Coefficient RELAXV for turbulence variables was modified',   /,&
'@ although IKECOU is not = 0.      It is in fact to',          /,&
'@ - for k                   :', e12.4,                         /,&
'@ - for  epsilon (or omega) :', e12.4,                         /,&
'@',                                                            /,&
'@ The  modification will be ignored (RELAXV is only useful',   /,&
'@ if IKECOU=0)',                                               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2624 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',                                                            /,&
'@ LCoefficient RELAXV for turbulence variables must',          /,&
'@ be a REAL comprised between 0 and 1. IT IS EQUAL :',         /,&
'@ - for  k                  :', e12.4,                         /,&
'@ - for  epsilon (ou omega) :', e12.4,                         /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2625 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',                                                            /,&
'@ Coefficient RELAXV for the pressure must be a REAL',         /,&
'@ between  0 et 1. It  IS EQUAL :', e12.4,                     /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2630 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',                                                            /,&
'@   Van Driest near wall damping was selected',                /,&
'@    Van Driest (IDRIES = ', i10,')',                          /,&
'@    with a LES model not compatible (ITURB = ', i10,')',      /,&
'@    (dynamic model or WALE model)',                           /,&
'@',                                                            /,&
'@  The calculation could NOT run.',                            /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2640 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' must be a REAL number in the range  [0;1]',        /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2650 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    SCALAR ',   a16,                                          /,&
'@',    a6,'(',i10,   ') MUST BE AN INTEGER EQUAL  0  OR 1',    /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  ICPSYR(I) is a flag for coupling scalar I with',            /,&
'@    SYRTHES  (solid walls modeller)',                         /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2660 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    INCOHERENCE POUR LE COUPLAGE SYRTHES',                    /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  No coupling with  SYRTHES has been defined',                /,&
'@  The number of coupled scalars is however',     i10,         /,&
'@    (the total number of scalars is',i10,   ')',              /,&
'@  Verify  the couplage SYRTHES-Noyau.',                       /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2661 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    INCOHERENCE POUR LE COUPLAGE SYRTHES',                    /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  A coupling with   SYRTHES was defined',                     /,&
'@  This requires a coupled scalar (and only one).',            /,&
'@  The total number of scalars   is', i10,                     /,&
'@  The number of coupled scalars is', i10,                     /,&
'@  Verify   le couplage SYRTHES-Noyau.',                       /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2662 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    INCOHERENCE for coupling with SYRTHES',                   /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@',                                                            /,&
'@  The coupled scalar must be the temperature',                /,&
'@',                                                            /,&
'@  Here it is scalar number :', i10,                           /,&
'@    which is not temperature because',                        /,&
'@                    ISCACP(',i10,   ') = ', i10,              /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2663 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    INCOHERENCE IN COUPLING WITH SYRTHES',                    /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@',                                                            /,&
'@  For compressible and SYRTHES coupling',                     /,&
'@    the coupled scalar must be energy.',                      /,&
'@  here it is scalar ',                      i10,              /,&
'@    which is not the energy  here since',                     /,&
'@                    ISCALT = ', i10,                          /,&
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
 2710 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' MUST BE   A  REAL INCLUD IN SEGMENT       [0.;1.]',/,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2720 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' MUST BE   A  REAL SUPERIEUR or EGAL A 1.',         /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2730 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' MUST BE   A  REAL EGAL A 0. or A 1.',              /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
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
 2742 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING:',                                                /,&
'@    =========',                                               /,&
'@    iswdyn = 1 OPTION: NSWRSM(IPR) > 19 IS REQUIRED',         /,&
'@    NSWRSM(IPR) HAS VALUE', i10,                              /,&
'@',                                                            /,&
'@  The calculation continue with nswrsm(ipr) = 20.',           /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2743 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING:',                                                /,&
'@    =========',                                               /,&
'@    SWPDYN = 1 OPTION: THE SLOPE TEST ON U, V, W MUST BE',    /,&
'@    DEACTIVATED.',                                            /,&
'@',                                                            /,&
'@  The calculation continue with:',                            /,&
'@  isstpc(iu) = isstpc(iu) = isstpc(iu) = 1.',                 /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 2747 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@   WARNING :      WHEN READING INPUT DATA',                /,&
'@    =========',                                               /,&
'@   THE PARAMETER NSWRSM FOR THE',                             /,&
'@        VARIABLE', a16, 'MUST BE POSITIVE, HAS VALUE',  i10,  /,&
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
 2750 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@',    a6,' MUST BE AN INTEGER  LESS THAN or EGAL A 1',        /,&
'@   IT HAS VALUE', i10,                                        /,&
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

 4100 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ATTENTION : ENTREE DES DONNEES',                          /,&
'@    =========',                                               /,&
'@    REFERENCE VELOCITY  UREF WAS NOT DEFINED',                /,&
'@    or is ill defined  (NEGATIVE value ? ).',                 /,&
'@   It IS EQUAL', e14.5,                                       /,&
'@',                                                            /,&
'@  calculation can only run if turbulence is defined',         /,&
'@  via a restart file or the interface or the',                /,&
'@  cs_user_initialization subroutine.',                        /,&
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
 4320 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    SCALAR ', a16,                                            /,&
'@    iscavr(',i10,   ') MUST BE AN INTEGER',                   /,&
'@      POSITIVE or NUL AND',                                   /,&
'@      LESS THAN or EGAL A NSCAL = ', i10,                     /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  If iscavr(I) =0,  le scalare  I is not a  variance',        /,&
'@  If iscavr(I) is POSITIVE, scalare I is a variance :',       /,&
'@    it is the variance of fluctuations of scalaire J',        /,&
'@    who''s number is   iscavr(I)',                            /,&
'@ Check the input data given through the User Interface',      /,&
'@   or in cs_user_parameters.f90.',                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 4321 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    THE SCALAR ', a16, 'IS DEFINED AS        FLUCTUATION',    /,&
'@    OF SCALAR ', a16,                                         /,&
'@    (iscavr(',i10,   ') = ', i10,   '),',                     /,&
'@    WHICH ITSELF IS DEFINED AS A FLUCTUATION',                /,&
'@    (iscavr(',i10,   ') = ', i10,   ').',                     /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  If iscavr(I) is POSITIVE, scalar  I is a  variance :',      /,&
'@    variance of fluctuations of scalar J',                    /,&
'@    who''s number is iscavr(I),  so we must have',            /,&
'@    iscavr(J) = 0',                                           /,&
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
 4340 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    SCALAR ', a16,                                            /,&
'@    VISLS0(',i10,   ') MUST BE   A  POSITIVE REAL',           /,&
'@      WHILE scalar_diffusivity_id = ', i10,                   /,&
'@   IT HAS VALUE', e14.5,                                      /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  VISLS0(I) is the molecular diffusion coefficient  of the',  /,&
'@    scalar and  MUST BE POSITIVE when that field''s',         /,&
'@    scalar_diffusivity_id key value is less than 0',          /,&
'@    (in which case a variable coefficient is given',          /,&
'@    using the User Interface or usphyv).',                    /,&
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
 9000 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@    FLAG FOR CONSIDERATION OF ROTATION',                      /,&
'@',                                                            /,&
'@  ICORIO should be = 0 or 1',                                 /,&
'@   IT HAS VALUE', i10,                                        /,&
'@',                                                            /,&
'@   The calculation could NOT run.',                           /,&
'@',                                                            /,&
'@  Verify the parameters given in cs_user_parameter.f90',      /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
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
 9012  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@   UNSTEADY ROTOR/STATOR COUPLING IS NOT COMPATIBLE',         /,&
'@     WITH THE LAGRANGIAN MODULE',                             /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9020  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@   THE COMPUTATION IN A RELATIVE FRAME OF REFERENCE',         /,&
'@     (ICORIO = 1) IS NOT COMPATIBLE WITH THE TURBOMACHINERY', /,&
'@     MODULE (ITURBO > 0)',                                    /,&
'@',                                                            /,&
'@  If necessary, use only the turbomachinery module defining', /,&
'@    several rotors',                                          /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9110  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA',               /,&
'@    =========',                                               /,&
'@   THE CAVITATION MODEL IS NOT COMPATIBLE WITH THE HANDLING', /,&
'@     OF HYDROSTATIC PRESSURE ALGORITHMS',                     /,&
'@',                                                            /,&
'@  Computation CAN NOT run',                                   /,&
'@',                                                            /,&
'@  Integer parameter IPHYDR was set to', i10,                  /,&
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
'@   THE CAVITATION MODEL IS NOT COMPATIBLE WITH THE',          /,&
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
