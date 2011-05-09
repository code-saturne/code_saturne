!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine modini
!================


!===============================================================================
!  FONCTION  :
!  ---------

! MODIFICATION DES PARAMETRES DE CALCUL
!    APRES INTERVENTION UTILISATEUR
!   (COMMONS)
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
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use albase
use alstru
use cplsat

!===============================================================================

implicit none

! Arguments


! Local variables

integer          ii, jj, ivar, iphas, iok, iest, imom, ikw
integer          icompt, ipp, nbccou, nn
integer          nscacp, iscal
double precision relxsp
double precision omgnrm, ctheta, stheta
double precision ux, uy, uz

!===============================================================================

! Indicateur erreur (0 : pas d'erreur)
iok = 0

!===============================================================================
! 1. ENTREES SORTIES entsor.h
!===============================================================================

! ---> Niveau d'impression listing
!       Non initialise -> standard
do ii = 1, nvarmx
  if(iwarni(ii).eq.-10000) then
    iwarni(ii) = 0
  endif
enddo

!---> Variables de calcul ITRSVR = ivar. Sinon 0 (valeur d'initialisation).

do ivar = 1, nvar
  itrsvr(ipprtp(ivar)) = ivar
enddo

!---> sorties chrono?
!     Sauf mention contraire de l'utilisateur, on sort a la fin les
!        variables de calcul, la viscosite, rho, le pas de temps s'il
!        est variable, les estimateurs s'ils sont actives, les moments
!        s'il y en a et la viscosite de maillage en ALE.

do ii = 2, nvppmx
  if(itrsvr(ii).ge.1.and.ichrvr(ii).eq.-999) then
    ichrvr(ii) = 1
  endif
enddo
do iphas = 1, nphas
  ipp = ipppro(ipproc(irom(iphas)))
  if(                       ichrvr(ipp).eq.-999) ichrvr(ipp) = 1
  ipp = ipppro(ipproc(ivisct(iphas)))
  if( (iturb(iphas).eq.10 .or. itytur(iphas).eq.2                 &
       .or. iturb(iphas).eq.50 .or. iturb(iphas).eq.60            &
       .or. iturb(iphas).eq.70 )                                  &
       .and.ichrvr(ipp).eq.-999) ichrvr(ipp) = 1
  if (idtvar.lt.0) then
    ichrvr(ipppro(ipproc(icour(iphas)))) = 0
    ichrvr(ipppro(ipproc(ifour(iphas)))) = 0
  endif
  do iest = 1, nestmx
    if(iescal(iest,iphas).gt.0) then
      ipp = ipppro(ipproc(iestim(iest,iphas)))
      if(                     ichrvr(ipp).eq.-999) ichrvr(ipp) = 1
    endif
  enddo
enddo
if(idtvar.eq.2.and.ichrvr(ippdt).eq.-999) ichrvr(ippdt) = 1
if(ipucou.ne.1) then
  ichrvr(ipptx) = 0
  ichrvr(ippty) = 0
  ichrvr(ipptz) = 0
endif

if(nbmomt.gt.0) then
  do imom = 1, nbmomt
    ipp = ipppro(ipproc(icmome(imom)))
    if(ichrvr(ipp).eq.-999) ichrvr(ipp) = 1
  enddo
endif
if (iale.eq.1) then
  nn = 1
  if (iortvm.eq.1) nn = 3
  do ii = 1, nn
    ipp = ipppro(ipproc(ivisma(ii)))
    if(ichrvr(ipp).eq.-999) ichrvr(ipp) = 1
  enddo
endif

do ii = 1, nvppmx
  if(ichrvr(ii).eq.-999) then
    ichrvr(ii) = 0
  endif
enddo

icompt = 0
do ii = 2, nvppmx
  if(ichrvr(ii).eq.1) icompt = icompt+1
enddo
if(icompt.eq.0) then
  ntchr = -1
  frchr = -1.d0
endif

! Adapt the output frequency parameters according to the time scheme.
if (idtvar.lt.0.or.idtvar.eq.2) then
  frchr = -1.d0
else
  if (frchr > 0.d0) then
    ntchr = -1
  endif
endif

!---> sorties historiques ?
!      Si une valeur non modifiee par l'utilisateur (=-999)
!        on la met a sa valeur par defaut
!      On sort toutes les variables a tous les pas de temps par defaut
!      IHISVR nb de sonde et numero par variable (-999 non initialise)
!             -1 : toutes les sondes
!      NTHIST = -1 : on ne sort pas d'historiques
!      NTHIST =  n : on sort des historiques tous les n pas de temps
!      NTHSAV = -1 : on sauvegarde a la fin uniquement
!      NTHSAV =  0 : periode par defaut (voir caltri)
!             > 0  : periode

do ii = 2, nvppmx
  if(itrsvr(ii).ge.1.and.ihisvr(ii,1).eq.-999) then
    ihisvr(ii,1) = -1
  endif
enddo
if(ihisvr(ippdt ,1).eq.-999) ihisvr(ippdt ,1) = -1
if(ipucou.ne.1) then
  ihisvr(ipptx,1) = 0
  ihisvr(ippty,1) = 0
  ihisvr(ipptz,1) = 0
endif
do iphas = 1, nphas
  ipp = ipppro(ipproc(ivisct(iphas)))
  if( (iturb(iphas).eq.10 .or. itytur(iphas).eq.2                 &
       .or. iturb(iphas).eq.50 .or. iturb(iphas).eq.60            &
       .or. iturb(iphas).eq.70 )                                  &
       .and.ihisvr(ipp,1).eq.-999) ihisvr(ipp,1) = -1
  if (idtvar.lt.0) then
    ihisvr(ipppro(ipproc(icour(iphas))),1) = 0
    ihisvr(ipppro(ipproc(ifour(iphas))),1) = 0
  endif
enddo
if(nbmomt.gt.0) then
  do imom = 1, nbmomt
    ipp = ipppro(ipproc(icmome(imom)))
    if(ihisvr(ipp,1).eq.-999) ihisvr(ipp,1) = -1
  enddo
endif

do ii = 1, nvppmx
  if(ihisvr(ii,1).eq.-999) then
    ihisvr(ii,1) = 0
  endif
enddo

!     Si on est en ALE, on a un test equivalent dans strini.F
if (iale.eq.0) then
  icompt = 0
  do ii = 2, nvppmx
    if(ihisvr(ii,1).ne.0) icompt = icompt+1
  enddo

  if(icompt.eq.0.or.ncapt.eq.0) then
    nthist = -1
    frhist = -1.d0
  endif
endif

! Adapt the output frequency parameters according to the time scheme.
if (idtvar.lt.0.or.idtvar.eq.2) then
  frhist = -1.d0
else
  if (frhist > 0.d0) then
    nthist = -1
  endif
endif

! ---> Nom des variables

do iphas = 1, nphas

  IF(NOMVAR(IPPRTP(IPR   (IPHAS))) .EQ.' ') THEN
    WRITE(NOMVAR(IPPRTP(IPR   (IPHAS))),'(A6,I2.2)')'PresPh',IPHAS
  endif
  IF(NOMVAR(IPPRTP(IU    (IPHAS))) .EQ.' ') THEN
    WRITE(NOMVAR(IPPRTP(IU    (IPHAS))),'(A6,I2.2)')'VitesX',IPHAS
  endif
  IF(NOMVAR(IPPRTP(IV    (IPHAS))) .EQ.' ') THEN
    WRITE(NOMVAR(IPPRTP(IV    (IPHAS))),'(A6,I2.2)')'VitesY',IPHAS
  endif
  IF(NOMVAR(IPPRTP(IW    (IPHAS))) .EQ.' ') THEN
    WRITE(NOMVAR(IPPRTP(IW    (IPHAS))),'(A6,I2.2)')'VitesZ',IPHAS
  endif
  if(itytur(iphas).eq.2) then
    IF(NOMVAR(IPPRTP(IK    (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(IK    (IPHAS))),'(A6,I2.2)')            &
                                                    'EnTurb',IPHAS
    endif
    IF(NOMVAR(IPPRTP(IEP   (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(IEP   (IPHAS))),'(A6,I2.2)')            &
                                                    'Dissip',IPHAS
    endif
  elseif(itytur(iphas).eq.3) then
    IF(NOMVAR(IPPRTP(IR11  (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(IR11  (IPHAS))),'(A6,I2.2)')            &
                                                    'R11pha',IPHAS
    endif
    IF(NOMVAR(IPPRTP(IR22  (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(IR22  (IPHAS))),'(A6,I2.2)')            &
                                                    'R22pha',IPHAS
    endif
    IF(NOMVAR(IPPRTP(IR33  (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(IR33  (IPHAS))),'(A6,I2.2)')            &
                                                    'R33pha',IPHAS
    endif
    IF(NOMVAR(IPPRTP(IR12  (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(IR12  (IPHAS))),'(A6,I2.2)')            &
                                                    'R12pha',IPHAS
    endif
    IF(NOMVAR(IPPRTP(IR13  (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(IR13  (IPHAS))),'(A6,I2.2)')            &
                                                    'R13pha',IPHAS
    endif
    IF(NOMVAR(IPPRTP(IR23  (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(IR23  (IPHAS))),'(A6,I2.2)')            &
                                                    'R23pha',IPHAS
    endif
    IF(NOMVAR(IPPRTP(IEP   (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(IEP   (IPHAS))),'(A6,I2.2)')            &
                                                    'Dissip',IPHAS
    endif
  elseif(iturb(iphas).eq.50) then
    IF(NOMVAR(IPPRTP(IK    (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(IK    (IPHAS))),'(A6,I2.2)')            &
                                                    'EnTurb',IPHAS
    endif
    IF(NOMVAR(IPPRTP(IEP   (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(IEP   (IPHAS))),'(A6,I2.2)')            &
                                                    'Dissip',IPHAS
    endif
    IF(NOMVAR(IPPRTP(IPHI  (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(IPHI  (IPHAS))),'(A6,I2.2)')            &
                                                    'phipha',IPHAS
    endif
    IF(NOMVAR(IPPRTP(IFB   (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(IFB   (IPHAS))),'(A6,I2.2)')            &
                                                    'fbarre',IPHAS
    endif
  elseif(iturb(iphas).eq.60) then
    IF(NOMVAR(IPPRTP(IK    (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(IK    (IPHAS))),'(A6,I2.2)')            &
                                                    'EnTurb',IPHAS
    endif
    IF(NOMVAR(IPPRTP(IOMG  (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(IOMG  (IPHAS))),'(A5,I2.2)')            &
                                                    'Omega',IPHAS
    endif
  elseif(iturb(iphas).eq.70) then
    IF(NOMVAR(IPPRTP(INUSA (IPHAS))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPRTP(INUSA (IPHAS))),'(A6,I2.2)')            &
                                                    'NuTild',IPHAS
    endif
  endif

  IF(NOMVAR(IPPPRO(IPPROC(IROM  (IPHAS)))) .EQ.' ') THEN
    WRITE(NOMVAR(IPPPRO(IPPROC(IROM  (IPHAS)))),'(A6,I2.2)')      &
                                                    'MasVol',IPHAS
  endif
  IF(NOMVAR(IPPPRO(IPPROC(IVISCT(IPHAS)))) .EQ.' ') THEN
    WRITE(NOMVAR(IPPPRO(IPPROC(IVISCT(IPHAS)))),'(A6,I2.2)')      &
                                                    'VisTur',IPHAS
  endif
  IF(NOMVAR(IPPPRO(IPPROC(IVISCL(IPHAS)))) .EQ.' ') THEN
    WRITE(NOMVAR(IPPPRO(IPPROC(IVISCL(IPHAS)))),'(A6,I2.2)')      &
                                                    'VisMol',IPHAS
  endif
  if (ismago(iphas).gt.0) then
    IF(NOMVAR(IPPPRO(IPPROC(ISMAGO(IPHAS)))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPPRO(IPPROC(ISMAGO(IPHAS)))),'(A6,I2.2)')    &
                                                    'Csdyn2',IPHAS
    endif
  endif
  if(icp   (iphas).gt.0) then
    IF(NOMVAR(IPPPRO(IPPROC(ICP   (IPHAS)))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPPRO(IPPROC(ICP   (IPHAS)))),'(A6,I2.2)')    &
                                                    'ChalSp',IPHAS
    endif
  endif
  if(iescal(iespre,iphas).gt.0) then
    ipp = ipppro(ipproc(iestim(iespre,iphas)))
    IF(NOMVAR(IPP) .EQ.' ') THEN
      WRITE(NOMVAR(IPP),'(A5,I1,I2.2)')                           &
                             'EsPre',IESCAL(IESPRE,IPHAS),IPHAS
    endif
  endif
  if(iescal(iesder,iphas).gt.0) then
    ipp = ipppro(ipproc(iestim(iesder,iphas)))
    IF(NOMVAR(IPP) .EQ.' ') THEN
      WRITE(NOMVAR(IPP),'(A5,I1,I2.2)')                           &
                             'EsDer',IESCAL(IESDER,IPHAS),IPHAS
    endif
  endif
  if(iescal(iescor,iphas).gt.0) then
    ipp = ipppro(ipproc(iestim(iescor,iphas)))
    IF(NOMVAR(IPP) .EQ.' ') THEN
      WRITE(NOMVAR(IPP),'(A5,I1,I2.2)')                           &
                             'EsCor',IESCAL(IESCOR,IPHAS),IPHAS
    endif
  endif
  if(iescal(iestot,iphas).gt.0) then
    ipp = ipppro(ipproc(iestim(iestot,iphas)))
    IF(NOMVAR(IPP) .EQ.' ') THEN
      WRITE(NOMVAR(IPP),'(A5,I1,I2.2)')                           &
                             'EsTot',IESCAL(IESTOT,IPHAS),IPHAS
    endif
  endif

enddo

do jj = 1, nscaus
  ii = jj
  IF(NOMVAR(IPPRTP(ISCA  (II   ))) .EQ.' ') THEN
    WRITE(NOMVAR(IPPRTP(ISCA  (II   ))),'(A5,I3.3)')'Scaus' ,II
  endif
enddo
do jj = 1, nscapp
  ii = iscapp(jj)
  IF(NOMVAR(IPPRTP(ISCA  (II   ))) .EQ.' ') THEN
    WRITE(NOMVAR(IPPRTP(ISCA  (II   ))),'(A5,I3.3)')'Scapp' ,II
  endif
enddo

if(nbmomt.gt.0) then
  do imom = 1, nbmomt
    ipp = ipppro(ipproc(icmome(imom)))
    IF(NOMVAR(IPP) .EQ.' ') THEN
      WRITE(NOMVAR(IPP),'(A6,I2.2)')'MoyTps',IMOM
    endif
  enddo
endif

IF(NOMVAR(IPPDT                ) .EQ.' ') THEN
  WRITE(NOMVAR(IPPDT                ),'(A8     )')'PasDeTmp'
endif

IF(NOMVAR(IPPTX                ) .EQ.' ') THEN
  WRITE(NOMVAR(IPPTX                ),'(A8     )')'Tx      '
endif
IF(NOMVAR(IPPTY                ) .EQ.' ') THEN
  WRITE(NOMVAR(IPPTY                ),'(A8     )')'Ty      '
endif
IF(NOMVAR(IPPTZ                ) .EQ.' ') THEN
  WRITE(NOMVAR(IPPTZ                ),'(A8     )')'Tz      '
endif

if (iale.eq.1) then
  IF(NOMVAR(IPPRTP(IUMA)) .EQ.' ') THEN
    WRITE(NOMVAR(IPPRTP(IUMA)),'(A8)')'VitmailX'
  endif
  IF(NOMVAR(IPPRTP(IVMA)) .EQ.' ') THEN
    WRITE(NOMVAR(IPPRTP(IVMA)),'(A8)')'VitmailY'
  endif
  IF(NOMVAR(IPPRTP(IWMA)) .EQ.' ') THEN
    WRITE(NOMVAR(IPPRTP(IWMA)),'(A8)')'VitmailZ'
  endif
  if (iortvm.eq.0) then
    IF(NOMVAR(IPPPRO(IPPROC(IVISMA(1)))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPPRO(IPPROC(IVISMA(1)))),'(A8)')'ViscMail'
    endif
  else
    IF(NOMVAR(IPPPRO(IPPROC(IVISMA(1)))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPPRO(IPPROC(IVISMA(1)))),'(A8)')'ViscMaiX'
    endif
    IF(NOMVAR(IPPPRO(IPPROC(IVISMA(2)))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPPRO(IPPROC(IVISMA(2)))),'(A8)')'ViscMaiY'
    endif
    IF(NOMVAR(IPPPRO(IPPROC(IVISMA(3)))) .EQ.' ') THEN
      WRITE(NOMVAR(IPPPRO(IPPROC(IVISMA(3)))),'(A8)')'ViscMaiZ'
    endif
  endif
endif

! ---> Sorties listing

do iphas = 1, nphas
  ipp = ipppro(ipproc(irom  (iphas)))
  if(irovar(iphas).eq.1.and.ilisvr(ipp).eq.-999) ilisvr(ipp) = 1
  ipp = ipppro(ipproc(ivisct(iphas)))
  if( (iturb(iphas).eq.10 .or. itytur(iphas).eq.2                 &
       .or. iturb(iphas).eq.50 .or. iturb(iphas).eq.60            &
       .or. iturb(iphas).eq.70 )                                  &
       .and.ilisvr(ipp).eq.-999) ilisvr(ipp) = 1
  if (inusa(iphas) .gt. 0) then
    ipp = ipppro(ipproc(inusa(iphas)))
    if(iturb(iphas).eq.70.and.ilisvr(ipp).eq.-999) ilisvr(ipp) = 1
  endif
  ipp = ipppro(ipproc(icour(iphas)))
  if (ilisvr(ipp).eq.-999 .or. idtvar.lt.0) ilisvr(ipp) = 0
  ipp = ipppro(ipproc(ifour(iphas)))
  if (ilisvr(ipp).eq.-999 .or. idtvar.lt.0) ilisvr(ipp) = 0
  if(iescal(iespre,iphas).gt.0) then
    ipp = ipppro(ipproc(iestim(iespre,iphas)))
    if(                     ilisvr(ipp).eq.-999) ilisvr(ipp) = 1
  endif
  if(iescal(iesder,iphas).gt.0) then
    ipp = ipppro(ipproc(iestim(iesder,iphas)))
    if(                     ilisvr(ipp).eq.-999) ilisvr(ipp) = 1
  endif
  if(iescal(iescor,iphas).gt.0) then
    ipp = ipppro(ipproc(iestim(iescor,iphas)))
    if(                     ilisvr(ipp).eq.-999) ilisvr(ipp) = 1
  endif
  if(iescal(iestot,iphas).gt.0) then
    ipp = ipppro(ipproc(iestim(iestot,iphas)))
    if(                     ilisvr(ipp).eq.-999) ilisvr(ipp) = 1
  endif
enddo

if(nbmomt.gt.0) then
  do imom = 1, nbmomt
    ipp = ipppro(ipproc(icmome(imom)))
    if(ilisvr(ipp).eq.-999) ilisvr(ipp) = 1
  enddo
endif

if(ilisvr(ippdt).eq.-999) ilisvr(ippdt)  = 1
if(ipucou.ne.1 .or. idtvar.lt.0) then
  ilisvr(ipptx) = 0
  ilisvr(ippty) = 0
  ilisvr(ipptz) = 0
endif

do ii = 2, nvppmx
  if(itrsvr(ii).ge.1.and.ilisvr(ii).eq.-999) then
    ilisvr(ii) = 1
  endif
enddo
do ii = 1, nvppmx
  if(ilisvr(ii).eq.-999) then
    ilisvr(ii) = 0
  endif
enddo



!===============================================================================
! 2. POSITION DES VARIABLES DE numvar.h
!===============================================================================

! ---> Reperage des variables qui disposeront de deux types de CL

!     Fait dans varpos.
!     Si l'utilisateur y a touche ensuite, on risque l'incident.

!===============================================================================
! 3. OPTIONS DU CALCUL : TABLEAUX DE optcal.h
!===============================================================================

! ---> restart

call indsui(isuite)
!==========

! ---> Schema en temps

do iphas = 1, nphas

!   -- Flux de masse
  if(abs(thetfl(iphas)+999.d0).gt.epzero) then
    write(nfecra,1001) iphas,istmpf(iphas)
    iok = iok + 1
  elseif(istmpf(iphas).eq.0) then
    thetfl(iphas) = 0.d0
  elseif(istmpf(iphas).eq.2) then
    thetfl(iphas) = 0.5d0
  endif

!    -- Proprietes physiques
  if(abs(thetro(iphas)+999.d0).gt.epzero) then
    WRITE(NFECRA,1011) IPHAS,'IROEXT',IROEXT(IPHAS),'THETRO'
    iok = iok + 1
  elseif(iroext(iphas).eq.0) then
    thetro(iphas) = 0.0d0
  elseif(iroext(iphas).eq.1) then
    thetro(iphas) = 0.5d0
  elseif(iroext(iphas).eq.2) then
    thetro(iphas) = 1.d0
  endif
  if(abs(thetvi(iphas)+999.d0).gt.epzero) then
    WRITE(NFECRA,1011) IPHAS,'IVIEXT',IVIEXT(IPHAS),'THETVI'
    iok = iok + 1
  elseif(iviext(iphas).eq.0) then
    thetvi(iphas) = 0.0d0
  elseif(iviext(iphas).eq.1) then
    thetvi(iphas) = 0.5d0
  elseif(iviext(iphas).eq.2) then
    thetvi(iphas) = 1.d0
  endif
  if(abs(thetcp(iphas)+999.d0).gt.epzero) then
    WRITE(NFECRA,1011) IPHAS,'ICPEXT',ICPEXT(IPHAS),'THETCP'
    iok = iok + 1
  elseif(icpext(iphas).eq.0) then
    thetcp(iphas) = 0.0d0
  elseif(icpext(iphas).eq.1) then
    thetcp(iphas) = 0.5d0
  elseif(icpext(iphas).eq.2) then
    thetcp(iphas) = 1.d0
  endif

!    -- Termes sources NS
  if(abs(thetsn(iphas)+999.d0).gt.epzero) then
    WRITE(NFECRA,1011) IPHAS,'ISNO2T',ISNO2T(IPHAS),'THETSN'
    iok = iok + 1
  elseif(isno2t(iphas).eq.1) then
    thetsn(iphas) = 0.5d0
  elseif(isno2t(iphas).eq.2) then
    thetsn(iphas) = 1.d0
  elseif(isno2t(iphas).eq.0) then
    thetsn(iphas) = 0.d0
  endif

!    -- Termes sources grandeurs turbulentes
  if(abs(thetst(iphas)+999.d0).gt.epzero) then
    WRITE(NFECRA,1011) IPHAS,'ISTO2T',ISTO2T(IPHAS),'THETST'
    iok = iok + 1
  elseif(isto2t(iphas).eq.1) then
    thetst(iphas) = 0.5d0
  elseif(isto2t(iphas).eq.2) then
    thetst(iphas) = 1.d0
  elseif(isto2t(iphas).eq.0) then
    thetst(iphas) = 0.d0
  endif

enddo

do iscal = 1, nscal
!    -- Termes sources des scalaires
  if(abs(thetss(iscal)+999.d0).gt.epzero) then
    WRITE(NFECRA,1021) ISCAL,'ISSO2T',ISSO2T(ISCAL),'THETSS'
    iok = iok + 1
  elseif(isso2t(iscal).eq.1) then
    thetss(iscal) = 0.5d0
  elseif(isso2t(iscal).eq.2) then
    thetss(iscal) = 1.d0
  elseif(isso2t(iscal).eq.0) then
    thetss(iscal) = 0.d0
  endif
!    -- Diffusivite des scalaires
  if(abs(thetvs(iscal)+999.d0).gt.epzero) then
    WRITE(NFECRA,1021) ISCAL,'IVSEXT',IVSEXT(ISCAL),'THETVS'
    iok = iok + 1
  elseif(ivsext(iscal).eq.0) then
    thetvs(iscal) = 0.d0
  elseif(ivsext(iscal).eq.1) then
    thetvs(iscal) = 0.5d0
  elseif(ivsext(iscal).eq.2) then
    thetvs(iscal) = 1.d0
  endif
enddo

!     Ici on interdit que l'utilisateur fixe lui meme THETAV, par securite
!       mais on pourrait le laisser faire
!       (enlever le IOK, modifier le message et les tests dans verini)
do iphas = 1, nphas

!     Vitesse pression (la pression est prise sans interp)
  if(abs(thetav(iu (iphas))+999.d0).gt.epzero.or.                 &
     abs(thetav(iv (iphas))+999.d0).gt.epzero.or.                 &
     abs(thetav(iw (iphas))+999.d0).gt.epzero.or.                 &
     abs(thetav(ipr(iphas))+999.d0).gt.epzero) then
    WRITE(NFECRA,1031) IPHAS,'VITESSE-PRESSION ','THETAV'
    iok = iok + 1
  elseif(ischtp(iphas).eq.1) then
    thetav(iu (iphas)) = 1.d0
    thetav(iv (iphas)) = 1.d0
    thetav(iw (iphas)) = 1.d0
    thetav(ipr(iphas)) = 1.d0
  elseif(ischtp(iphas).eq.2) then
    thetav(iu (iphas)) = 0.5d0
    thetav(iv (iphas)) = 0.5d0
    thetav(iw (iphas)) = 0.5d0
    thetav(ipr(iphas)) = 1.d0
  endif

!     Turbulence (en k-eps : ordre 1)
  if(itytur(iphas).eq.2) then
    if(abs(thetav(ik (iphas))+999.d0).gt.epzero.or.               &
       abs(thetav(iep(iphas))+999.d0).gt.epzero) then
      WRITE(NFECRA,1031) IPHAS,'VARIABLES   K-EPS','THETAV'
      iok = iok + 1
    elseif(ischtp(iphas).eq.1) then
      thetav(ik (iphas)) = 1.d0
      thetav(iep(iphas)) = 1.d0
    elseif(ischtp(iphas).eq.2) then
!     pour le moment, on ne peut pas passer par ici (cf varpos)
      thetav(ik (iphas)) = 0.5d0
      thetav(iep(iphas)) = 0.5d0
    endif
  elseif(itytur(iphas).eq.3) then
    if(abs(thetav(ir11(iphas))+999.d0).gt.epzero.or.              &
       abs(thetav(ir22(iphas))+999.d0).gt.epzero.or.              &
       abs(thetav(ir33(iphas))+999.d0).gt.epzero.or.              &
       abs(thetav(ir12(iphas))+999.d0).gt.epzero.or.              &
       abs(thetav(ir13(iphas))+999.d0).gt.epzero.or.              &
       abs(thetav(ir23(iphas))+999.d0).gt.epzero.or.              &
       abs(thetav(iep (iphas))+999.d0).gt.epzero) then
      WRITE(NFECRA,1031) IPHAS,'VARIABLES  RIJ-EP','THETAV'
      iok = iok + 1
    elseif(ischtp(iphas).eq.1) then
      thetav(ir11(iphas)) = 1.d0
      thetav(ir22(iphas)) = 1.d0
      thetav(ir33(iphas)) = 1.d0
      thetav(ir12(iphas)) = 1.d0
      thetav(ir13(iphas)) = 1.d0
      thetav(ir23(iphas)) = 1.d0
      thetav(iep (iphas)) = 1.d0
    elseif(ischtp(iphas).eq.2) then
      thetav(ir11(iphas)) = 0.5d0
      thetav(ir22(iphas)) = 0.5d0
      thetav(ir33(iphas)) = 0.5d0
      thetav(ir12(iphas)) = 0.5d0
      thetav(ir13(iphas)) = 0.5d0
      thetav(ir23(iphas)) = 0.5d0
      thetav(iep (iphas)) = 0.5d0
    endif
  elseif(iturb(iphas).eq.50) then
    if(abs(thetav(ik  (iphas))+999.d0).gt.epzero.or.              &
       abs(thetav(iep (iphas))+999.d0).gt.epzero.or.              &
       abs(thetav(iphi(iphas))+999.d0).gt.epzero.or.              &
       abs(thetav(ifb (iphas))+999.d0).gt.epzero) then
      WRITE(NFECRA,1031) IPHAS,'VARIABLES     V2F','THETAV'
      iok = iok + 1
    elseif(ischtp(iphas).eq.1) then
      thetav(ik  (iphas)) = 1.d0
      thetav(iep (iphas)) = 1.d0
      thetav(iphi(iphas)) = 1.d0
      thetav(ifb (iphas)) = 1.d0
    elseif(ischtp(iphas).eq.2) then
!     pour le moment, on ne peut pas passer par ici (cf varpos)
      thetav(ik  (iphas)) = 0.5d0
      thetav(iep (iphas)) = 0.5d0
      thetav(iphi(iphas)) = 0.5d0
      thetav(ifb (iphas)) = 0.5d0
    endif
  elseif(iturb(iphas).eq.60) then
    if(abs(thetav(ik  (iphas))+999.d0).gt.epzero.or.              &
       abs(thetav(iomg(iphas))+999.d0).gt.epzero ) then
      WRITE(NFECRA,1031) IPHAS,'VARIABLES K-OMEGA','THETAV'
      iok = iok + 1
    elseif(ischtp(iphas).eq.1) then
      thetav(ik  (iphas)) = 1.d0
      thetav(iomg(iphas)) = 1.d0
    elseif(ischtp(iphas).eq.2) then
!     pour le moment, on ne peut pas passer par ici (cf varpos)
      thetav(ik  (iphas)) = 0.5d0
      thetav(iomg(iphas)) = 0.5d0
    endif
  elseif(iturb(iphas).eq.70) then
    if(abs(thetav(inusa(iphas))+999.d0).gt.epzero) then
      WRITE(NFECRA,1031) IPHAS,'VARIABLE NU_tilde de SA','THETAV'
      iok = iok + 1
    elseif(ischtp(iphas).eq.1) then
      thetav(inusa(iphas)) = 1.d0
    elseif(ischtp(iphas).eq.2) then
!     pour le moment, on ne peut pas passer par ici (cf varpos)
      thetav(inusa(iphas)) = 0.5d0
    endif
  endif

enddo

!     Scalaires
do iscal = 1, nscal
  iphas = 1
  ivar  = isca(iscal)
  if(abs(thetav(ivar)+999.d0).gt.epzero) then
    WRITE(NFECRA,1041) IPHAS,'SCALAIRE',ISCAL,'THETAV'
    iok = iok + 1
  elseif(ischtp(iphas).eq.1) then
    thetav(ivar) = 1.d0
  elseif(ischtp(iphas).eq.2) then
    thetav(ivar) = 0.5d0
  endif
enddo

!     Vitesse de maillage en ALE
if (iale.eq.1) then
  iphas = 1
  if(abs(thetav(iuma)+999.d0).gt.epzero.or.                       &
     abs(thetav(ivma)+999.d0).gt.epzero.or.                       &
     abs(thetav(iwma)+999.d0).gt.epzero) then
    WRITE(NFECRA,1032) 'THETAV'
    iok = iok + 1
  elseif(ischtp(iphas).eq.1) then
    thetav(iuma) = 1.d0
    thetav(ivma) = 1.d0
    thetav(iwma) = 1.d0
  elseif(ischtp(iphas).eq.2) then
!     pour le moment, on ne peut pas passer par ici (cf varpos)
    thetav(iuma) = 0.5d0
    thetav(ivma) = 0.5d0
    thetav(iwma) = 0.5d0
  endif
endif

! ---> ISSTPC
!        Si l'utilisateur n'a rien specifie pour le test de pente (=-999),
!        On impose 1 (ie sans) pour la vitesse en LES
!                  0 (ie avec) sinon

do iphas = 1, nphas
  if(itytur(iphas).eq.4) then
    ii = iu(iphas)
    if(isstpc(ii).eq.-999) isstpc(ii) = 1
    ii = iv(iphas)
    if(isstpc(ii).eq.-999) isstpc(ii) = 1
    ii = iw(iphas)
    if(isstpc(ii).eq.-999) isstpc(ii) = 1
    do jj = 1, nscal
      ii = isca(jj)
      if(isstpc(ii).eq.-999) isstpc(ii) = 0
    enddo
  endif
enddo

do ii = 1, nvarmx
  if (isstpc(ii).eq.-999) then
    isstpc(ii) = 0
  endif
enddo

! ---> BLENCV
!        Si l'utilisateur n'a rien specifie pour le schema convectif
!                  1 (ie centre) pour les vitesses et
!                                      les scalaires utilisateurs
!                  0 (ie upwind pur) pour le reste
!   (en particulier, en L.E.S. toutes les variables sont donc en centre)

do iphas = 1, nphas
  ii = iu(iphas)
  if(abs(blencv(ii)+999.d0).lt.epzero) blencv(ii) = 1.d0
  ii = iv(iphas)
  if(abs(blencv(ii)+999.d0).lt.epzero) blencv(ii) = 1.d0
  ii = iw(iphas)
  if(abs(blencv(ii)+999.d0).lt.epzero) blencv(ii) = 1.d0
enddo
do jj = 1, nscaus
  ii = isca(jj)
  if(abs(blencv(ii)+999.d0).lt.epzero) blencv(ii) = 1.d0
enddo

do ii = 1, nvarmx
  if (abs(blencv(ii)+999.d0).lt.epzero) then
    blencv(ii) = 0.d0
  endif
enddo


! ---> NSWRSM, EPSRSM ET EPSILO
!        Si l'utilisateur n'a rien specifie  (NSWRSM=-999),
!        On impose
!           a l'ordre 1 :
!                  2  pour la pression
!                  1  pour les autres variables
!                  on initialise EPSRSM a 1.D-8
!                  on initialise EPSILO a 1.D-8
!           a l'ordre 2 :
!                  5  pour la pression
!                  10 pour les autres variables
!                  on initialise EPSRSM a 1.D-5
!                  on initialise EPSILO a 1.D-5
!     Attention aux tests dans verini

do iphas = 1, nphas
  if(ischtp(iphas).eq.2) then
    ii = ipr(iphas)
    if(nswrsm(ii).eq.-999) nswrsm(ii) = 5
    if(abs(epsrsm(ii)+999.d0).lt.epzero) epsrsm(ii) = 1.d-5
    if(abs(epsilo(ii)+999.d0).lt.epzero) epsilo(ii) = 1.d-5
    ii = iu(iphas)
    if(nswrsm(ii).eq.-999) nswrsm(ii) = 10
    if(abs(epsrsm(ii)+999.d0).lt.epzero) epsrsm(ii) = 1.d-5
    if(abs(epsilo(ii)+999.d0).lt.epzero) epsilo(ii) = 1.d-5
    ii = iv(iphas)
    if(nswrsm(ii).eq.-999) nswrsm(ii) = 10
    if(abs(epsrsm(ii)+999.d0).lt.epzero) epsrsm(ii) = 1.d-5
    if(abs(epsilo(ii)+999.d0).lt.epzero) epsilo(ii) = 1.d-5
    ii = iw(iphas)
    if(nswrsm(ii).eq.-999) nswrsm(ii) = 10
    if(abs(epsrsm(ii)+999.d0).lt.epzero) epsrsm(ii) = 1.d-5
    if(abs(epsilo(ii)+999.d0).lt.epzero) epsilo(ii) = 1.d-5
    do jj = 1, nscal
      ii = isca(jj)
      if(nswrsm(ii).eq.-999) nswrsm(ii) = 10
      if(abs(epsrsm(ii)+999.d0).lt.epzero) epsrsm(ii) = 1.d-5
      if(abs(epsilo(ii)+999.d0).lt.epzero) epsilo(ii) = 1.d-5
    enddo
  endif
  ii = ipr(iphas)
  if(nswrsm(ii).eq.-999) nswrsm(ii) = 2
enddo

do ii = 1, nvarmx
  if (nswrsm(ii).eq.-999) nswrsm(ii) = 1
  if(abs(epsrsm(ii)+999.d0).lt.epzero) epsrsm(ii) = 1.d-8
  if(abs(epsilo(ii)+999.d0).lt.epzero) epsilo(ii) = 1.d-8
enddo

! ---> ANOMAX
!        Si l'utilisateur n'a rien specifie pour l'angle de non
!          orthogonalite pour la selection du voisinage etendu,
!          on impose pi/4 (utile aussi en mode verifications)

if (anomax.le.-grand) then
  anomax = pi*0.25d0
endif

! ---> IMLIGR
!        Si l'utilisateur n'a rien specifie pour la limitation des
!          gradients (=-999),
!        On impose -1 avec gradrc (pas de limitation)
!               et  1 avec gradmc (limitation)

if (imrgra.eq.0.or.imrgra.eq.4) then
  do ii = 1, nvarmx
    if (imligr(ii).eq.-999) then
      imligr(ii) = -1
    endif
  enddo
elseif (imrgra.eq.1.or.imrgra.eq.2.or.imrgra.eq.3) then
  do ii = 1, nvarmx
    if (imligr(ii).eq.-999) then
      imligr(ii) = 1
    endif
  enddo
endif

! ---> DTMIN DTMAX CDTVAR


if(dtmin.le.-grand) then
  dtmin = 0.1d0*dtref
endif
if(dtmax.le.-grand) then
  dtmax = 1.0d3*dtref
endif

do iphas = 1, nphas

!     Ici, ce n'est pas grave pour le moment,
!      etant entendu que ces coefs ne servent pas
!      s'ils servaient, attention dans le cas a plusieurs phases avec
!      une seule pression : celle ci prend le coef de la derniere phase
  cdtvar(iv (iphas)) = cdtvar(iu(iphas))
  cdtvar(iw (iphas)) = cdtvar(iu(iphas))
  cdtvar(ipr(iphas)) = cdtvar(iu(iphas))

  if(itytur(iphas).eq.2) then
    cdtvar(iep (iphas)) = cdtvar(ik  (iphas))
  elseif(itytur(iphas).eq.3) then
    cdtvar(ir22(iphas)) = cdtvar(ir11(iphas))
    cdtvar(ir33(iphas)) = cdtvar(ir11(iphas))
    cdtvar(ir12(iphas)) = cdtvar(ir11(iphas))
    cdtvar(ir13(iphas)) = cdtvar(ir11(iphas))
    cdtvar(ir23(iphas)) = cdtvar(ir11(iphas))
    cdtvar(iep (iphas)) = cdtvar(ir11(iphas))
  elseif(iturb(iphas).eq.50) then
    cdtvar(iep (iphas)) = cdtvar(ik  (iphas))

    cdtvar(iphi(iphas)) = cdtvar(ik  (iphas))
!     CDTVAR(IFB) est en fait inutile car pas de temps dans l'eq de f_barre
    cdtvar(ifb (iphas)) = cdtvar(ik  (iphas))
  elseif(iturb(iphas).eq.60) then
    cdtvar(iomg(iphas)) = cdtvar(ik  (iphas))
  elseif(iturb(iphas).eq.70) then
! cdtvar est à 1.0 par defaut dans iniini.f90
    cdtvar(inusa(iphas))= cdtvar(inusa(iphas))
  endif

enddo

! ---> IDEUCH, YPLULI
!      En laminaire, longueur de melange, Spalar-Allmaras et LES,
!      une echelle de vitesse.
!      Sinon, 2 echelles, sauf si l'utilisateur choisit 1 echelle.
!      On a initialise IDEUCH a -999 pour voir si l'utilisateur essaye
!        de choisir deux echelles quand ce n'est pas possible et le
!        prevenir dans la section verification.

do iphas = 1, nphas
  if(ideuch(iphas).eq.-999) then
    if(iturb(iphas).eq. 0.or.                                     &
       iturb(iphas).eq.10.or.                                     &
       itytur(iphas).eq.4.or.                                     &
       iturb(iphas).eq.70) then
      ideuch(iphas) = 0
    else
      ideuch(iphas) = 1
    endif
  endif

  ! Pour YPLULI, 1/XKAPPA est la valeur qui assure la continuite de la derivee
  ! entre la zone lineaire et la zone logarithmique.

  ! Dans le cas des lois de paroi invariantes, on utilise la valeur de
  ! continuite du profil de vitesse, 10.88.

  ! Pour la LES, on remet 10.88, afin d'eviter des clic/clac quand on est a
  ! la limite (en modele a une echelle en effet, YPLULI=1/XKAPPA ne permet pas
  ! forcement de calculer u* de maniere totalement satisfaisante).
  ! Idem en Spalart-Allmaras.

  if (ypluli(iphas).lt.-grand) then
    if (ideuch(iphas).eq.2 .or. itytur(iphas).eq.4 .or.           &
        iturb(iphas).eq.70    ) then
      ypluli(iphas) = 10.88d0
    else
      ypluli(iphas) = 1.d0/xkappa
    endif
  endif
enddo


! ---> Van Driest
do iphas = 1, nphas
  if(idries(iphas).eq.-1) then
!   On met 1 en supposant qu'en periodicite ou parallele on utilise le
!     mode de calcul de la distance a la paroi qui les prend en charge
!     (ICDPAR=+/-1, valeur par defaut)
    if(iturb(iphas).eq.40) then
      idries(iphas) = 1
    elseif(iturb(iphas).eq.41) then
      idries(iphas) = 0
    elseif(iturb(iphas).eq.42) then
      idries(iphas) = 0
    endif
  endif
enddo


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

if(nscal.gt.0) then

!     On regarde s'il y a du couplage

  call nbcsyr (nbccou)
  !==========

!     S'il y a du couplage
  if (nbccou .ne. 0) then

!       On compte le nombre de scalaires couples
    nscacp = 0
    do iscal = 1, nscal
      if(icpsyr(iscal).eq.1) then
        nscacp = nscacp + 1
      endif
    enddo

!       Si l'utilisateur n'a pas couple de scalaire,
    if(nscacp.eq.0) then

!         On couple le scalaire temperature de la premiere phase
      do iphas = 1, nphas
        if(iscalt(iphas).gt.0.and.iscalt(iphas).le.nscal) then
          icpsyr(iscalt(iphas)) = 1
          goto 100
        endif
      enddo
 100        continue

    endif

  endif

!     Pour tous les autres scalaires, non renseignes pas l'utilisateur
!       on ne couple pas
  do iscal = 1, nscamx
    if(icpsyr(iscal).eq.-999) then
      icpsyr(iscal) = 0
    endif
  enddo

endif


! ---> ISCSTH
!      Si l'utilisateur n'a pas modifie ISCSTH, on prend par defaut :
!        scalaire passif  pour les scalaires autres que ISCALT
!      Les modifs adequates devront etre ajoutees pour les physiques
!        particulieres
!      Noter en outre que, par defaut, si on choisit temperature
!        elle est en K (ceci n'est utile que pour le rayonnement et les pp)

!         =-10: non renseigne
!         =-1 : temperature en C
!         = 0 : passif
!         = 1 : temperature en K
!         = 2 : enthalpie



if(nscal.gt.0) then
  do ii = 1, nscal
    if(iscsth(ii).eq.-10)then
      iphas = 1
      if(ii.ne.iscalt(iphas)) then
        iscsth(ii) = 0
      endif
    endif
  enddo
endif

! ---> ICALHY
!      Calcul de la pression hydrostatique en sortie pour les conditions de
!        Dirichlet sur la pression. Se deduit de IPHYDR et de la valeur de
!        la gravite (test assez arbitraire sur la norme).
!      ICALHY est initialise a -1 (l'utilisateur peut avoir force
!        0 ou 1 et dans ce cas, on ne retouche pas)

if(icalhy.ne.-1.and.icalhy.ne.0.and.icalhy.ne.1) then
  write(nfecra,1061)icalhy
  iok = iok + 1
endif


! ---> IDGMOM
!      Calcul du degre des moments

do imom = 1, nbmomx
  idgmom(imom) = 0
enddo
do imom = 1, nbmomt
  do ii = 1, ndgmox
    if(idfmom(ii,imom).ne.0) idgmom(imom) = idgmom(imom) + 1
  enddo
enddo

! ---> ICDPAR
!      Calcul de la distance a la paroi. En standard, on met ICDPAR a -1, au cas
!      ou les faces de bord auraient change de type d'un calcul a l'autre. En k-omega,
!      il faut la distance a la paroi pour une suite propre, donc on initialise a 1 et
!      on avertit (dans verini).
ikw = 0
do iphas = 1, nphas
  if (iturb(iphas).eq.60) ikw = 1
enddo
if (icdpar.eq.-999) then
  icdpar = -1
  if (ikw.eq.1) icdpar = 1
  if (isuite.eq.1 .and. ikw.eq.1) write(nfecra,2000)
endif
if (icdpar.eq.-1 .and. ikw.eq.1 .and. isuite.eq.1)                &
     write(nfecra,2001)

! ---> INEEDY, IMLIGY
!      Calcul de la distance a la paroi
!       (une seule phase ...)

ineedy = 0
do iphas = 1, nphas
  if((iturb(iphas).eq.30.and.irijec(iphas).eq.1).or.              &
     (itytur(iphas).eq.4.and.idries(iphas).eq.1).or.              &
      iturb(iphas).eq.60.or.iturb(iphas).eq.70      ) then
    ineedy = 1
  endif
enddo

if (imrgra.eq.0 .or. imrgra.eq.4) then
  if (imligy.eq.-999) then
    imligy = -1
  endif
elseif (imrgra.eq.1.or.imrgra.eq.2.or.imrgra.eq.3) then
  if (imligy.eq.-999) then
    imligy = 1
  endif
endif

!     Warning : non initialise => comme la vitesse
if(iwarny.eq.-999) then
  iwarny = iwarni(iu(1))
endif


! ---> IKECOU
!     En k-eps prod lin, v2f ou k-omega, on met IKECOU a 0 par defaut,
!     sinon on le laisse a 1
!     Dans verini on bloquera le v2f et le k-eps prod lin si IKECOU.NE.0
!     On bloquera aussi le stationnaire si IKECOU.NE.0
do iphas = 1, nphas
  if (ikecou(iphas).eq.-999) then
    if (idtvar.lt.0) then
      ikecou(iphas) = 0
    else if (iturb(iphas).eq.21 .or. iturb(iphas).eq.50           &
        .or. iturb(iphas).eq.60 ) then
      ikecou(iphas) = 0
    else
      ikecou(iphas) = 1
    endif
  endif
enddo

! ---> RELAXV
if (idtvar.lt.0) then
  relxsp = 1.d0-relxst
  if (relxsp.le.epzero) relxsp = relxst
  do iphas = 1, nphas
    if (abs(relaxv(ipr(iphas))+999.d0).le.epzero)                 &
         relaxv(ipr(iphas)) = relxsp
  enddo
  do ii = 1, nvarmx
    if (abs(relaxv(ii)+999.d0).le.epzero) relaxv(ii) = relxst
  enddo
else
  do iphas = 1, nphas
    if ( ikecou(iphas).eq.0) then
      if (itytur(iphas).eq.2 .or. itytur(iphas).eq.5) then
        if (abs(relaxv(ik(iphas))+999.d0).lt.epzero)              &
             relaxv(ik(iphas)) = 0.7d0
        if (abs(relaxv(iep(iphas))+999.d0).lt.epzero)             &
             relaxv(iep(iphas)) = 0.7d0
      else if (itytur(iphas).eq.6) then
        if (abs(relaxv(ik(iphas))+999.d0).lt.epzero)              &
             relaxv(ik(iphas)) = 0.7d0
        if (abs(relaxv(iomg(iphas))+999.d0).lt.epzero)            &
             relaxv(iomg(iphas)) = 0.7d0
      endif
    endif
    if(iturb(iphas).eq.70) then
      if(abs(relaxv(inusa(iphas))+999.d0).lt.epzero) then
        relaxv(inusa(iphas)) = 1.D0
      endif
    endif
    if (abs(relaxv(ipr(iphas))+999.d0).lt.epzero)                 &
             relaxv(ipr(iphas)) = 1.d0
  enddo
endif

! ---> SPECIFIQUE STATIONNAIRE
if (idtvar.lt.0) then
  dtref = 1.d0
  dtmin = 1.d0
  dtmax = 1.d0
  do ii = 1, nvarmx
    istat(ii) = 0
  enddo
  do iphas = 1, nphas
    arak(iphas) = arak(iphas)/max(relaxv(iu(iphas)),epzero)
  enddo
endif

! ---> INEEDF
!     Si on a demande un posttraitement des efforts aux bords, on
!     les calcule !
if (mod(ipstdv,ipstfo).eq.0) then
  ineedf = 1
endif
!     Si on est en ALE, par defaut on calcule les efforts aux bords
!     (un test eventuel sur la presence de structures viendrait trop
!     tard)
if (iale.eq.1) ineedf = 1

!===============================================================================
! 4. TABLEAUX DE cstphy.h
!===============================================================================

! ---> Constantes
!    Ca fait un calcul en double, mais si qqn a bouge cmu, apow, bpow,
!     ca servira.

cpow    = apow**(2.d0/(1.d0-bpow))
dpow    = 1.d0/(1.d0+bpow)
cmu025 = cmu**0.25d0

! ---> ICLVFL
!      Si l'utilisateur n'a pas modifie ICLVFL, on prend par defaut :
!        0 pour les variances
!      Les modifs adequates devront etre ajoutees pour les physiques
!        particulieres

do iscal = 1, nscal
  if(iscavr(iscal).gt.0) then
    if(iclvfl(iscal).eq.-1) then
      iclvfl(iscal) = 0
    endif
  endif
enddo


! ---> VISLS0 (IVISLS ont ete verifies dans varpos)
!      Pour les variances de fluctuations, les valeurs du tableau
!        precedent ne doivent pas avoir ete modifiees par l'utilisateur
!        Elles sont prises egales aux valeurs correspondantes pour le
!        scalaire associe.

if(nscal.gt.0) then
  do ii = 1, nscal
    iscal = iscavr(ii)
    if(iscal.gt.0.and.iscal.le.nscal)then
      if(visls0(ii).lt.-grand) then
        visls0(ii) = visls0(iscal)
      else
        write(nfecra,1071)ii,                                     &
          ii,iscal,ii,iscal,                                      &
          ii,ii,iscal,visls0(iscal)
        iok = iok + 1
      endif
    endif
  enddo
endif

! ---> XYZP0 : reference pour la pression hydrostatique
!      On considere que l'utilisateur a specifie la reference
!      a partir du moment ou il a specifie une coordonnee.
!      Pour les coordonnees non specifiees, on met 0.

do iphas = 1, nphas
  do ii = 1, 3
    if (xyzp0(ii,iphas).gt.-0.5d0*rinfin) then
      ixyzp0(iphas) = 1
    else
      xyzp0(ii,iphas) = 0.d0
    endif
  enddo
enddo

! Vecteur rotation et matrice(s) associees

omgnrm = sqrt(omegax**2 + omegay**2 + omegaz**2)

if (omgnrm.ge.epzero) then

  ! Normalized rotation vector

  ux = omegax / omgnrm
  uy = omegay / omgnrm
  uz = omegaz / omgnrm

  ! Matrice de projection sur l'axe de rotation

  prot(1,1) = ux**2
  prot(2,2) = uy**2
  prot(3,3) = uz**2

  prot(1,2) = ux*uy
  prot(2,1) = prot(1,2)

  prot(1,3) = ux*uz
  prot(3,1) = prot(1,3)

  prot(2,3) = uy*uz
  prot(3,2) = prot(2,3)

  ! Antisymetrc representation of Omega

  qrot(1,1) = 0.d0
  qrot(2,2) = 0.d0
  qrot(3,3) = 0.d0

  qrot(1,2) = -uz
  qrot(2,1) = -qrot(1,2)

  qrot(1,3) =  uy
  qrot(3,1) = -qrot(1,3)

  qrot(2,3) = -ux
  qrot(3,2) = -qrot(2,3)

  ! Matrice de rotation

  ctheta = cos(dtref*omgnrm)
  stheta = sin(dtref*omgnrm)

  do ii = 1, 3
    do jj = 1, 3
      irot(ii,jj) = 0.d0
    enddo
    irot(ii,ii) = 1.d0
  enddo

  do ii = 1, 3
    do jj = 1, 3
      rrot(ii,jj) = ctheta*irot(ii,jj) + (1.d0 - ctheta)*prot(ii,jj) &
                                       +         stheta *qrot(ii,jj)
    enddo
  enddo

else

  do ii = 1, 3
    do jj = 1, 3
      irot(ii,jj) = 0.d0
      prot(ii,jj) = 0.d0
      qrot(ii,jj) = 0.d0
      rrot(ii,jj) = 0.d0
    enddo
    irot(ii,ii) = 1.d0
    rrot(ii,ii) = 1.d0
  enddo

endif


!===============================================================================
! 5. ELEMENTS DE albase.h
!===============================================================================

if (iale.eq.1) then
  if (isuite.eq.0 .and. italin.eq.-999 ) italin = 1
else
  italin = 0
endif

!===============================================================================
! 6. COEFFICIENTS DE alstru.h
!===============================================================================

if (betnmk.lt.-0.5d0*grand)                                       &
     betnmk = (1.d0-alpnmk)**2/4.d0
if (gamnmk.lt.-0.5d0*grand)                                       &
     gamnmk = (1.d0-2.d0*alpnmk)/2.d0
if (aexxst.lt.-0.5d0*grand) aexxst = 0.5d0
if (bexxst.lt.-0.5d0*grand) bexxst = 0.0d0
if (cfopre.lt.-0.5d0*grand) cfopre = 2.0d0

!===============================================================================
! 7. PARAMETRES DE cplsat.h
!===============================================================================

! Get coupling number

call nbccpl(nbrcpl)
!==========

if (nbrcpl.ge.1) then
  ! Si on est en couplage rotor/stator avec resolution en repere absolu
  omgnrm = sqrt(omegax**2 + omegay**2 + omegaz**2)
  if (omgnrm.ge.epzero) then
    ! Couplage avec interpolation aux faces
    ifaccp = 1
    ! Maillage mobile
    if (icorio.eq.0) then
      imobil = 1
      ichrmd = 1
    endif
  endif
endif

!===============================================================================
! 8. STOP SI PB
!===============================================================================

if(iok.ne.0) then
  call csexit (1)
endif

#if defined(_CS_LANG_FR)

 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHASE ',   I10,' ISTMPF = ',   I10                       ,/,&
'@    THETFL SERA INITIALISE AUTOMATIQUEMENT.                 ',/,&
'@    NE PAS LE MODIFIER DANS usini1.                         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1011 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHASE ',   I10,' ',A6,' = ',   I10                       ,/,&
'@    ',A6,' SERA INITIALISE AUTOMATIQUEMENT.                 ',/,&
'@    NE PAS LE MODIFIER DANS usini1.                         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1021 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',   I10,' ',A6,' = ',   I10                    ,/,&
'@    ',A6,' SERA INITIALISE AUTOMATIQUEMENT.                 ',/,&
'@    NE PAS LE MODIFIER DANS usini1.                         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1031 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHASE ',   I10,' ',A17                                   ,/,&
'@    ',A6,' SERA INITIALISE AUTOMATIQUEMENT.                 ',/,&
'@    NE PAS LE MODIFIER DANS usini1.                         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1032 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VITESSE DE MAILLAGE EN ALE                              ',/,&
'@    ',A6,' SERA INITIALISE AUTOMATIQUEMENT.                 ',/,&
'@    NE PAS LE MODIFIER DANS usini1.                         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1041 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHASE ',   I10,' ',A8,' ',I10                            ,/,&
'@    ',A6,' SERA INITIALISE AUTOMATIQUEMENT.                 ',/,&
'@    NE PAS LE MODIFIER DANS usini1.                         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1061 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ICALHY doit etre un entier egal a 0 ou 1                ',/,&
'@                                                            ',/,&
'@  Il vaut ici ',I10                                          ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1071 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',I10   ,' NE PAS MODIFIER LA DIFFUSIVITE      ',/,&
'@                                                            ',/,&
'@  Le scalaire ',I10   ,' represente la variance des         ',/,&
'@    fluctuations du scalaire ',I10                           ,/,&
'@                             (ISCAVR(',I10   ,') = ',I10     ,/,&
'@  La diffusivite VISLS0(',I10   ,') du scalaire ',I10        ,/,&
'@    ne doit pas etre renseignee :                           ',/,&
'@    elle sera automatiquement prise egale a la diffusivite  ',/,&
'@    du scalaire ',I10   ,' soit ',E14.5                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  Le modele de turbulence k-omega a ete choisi. Pour gerer  ',/,&
'@    correctement la suite de calcul, l''indicateur ICDPAR a ',/,&
'@    ete mis a 1 (relecture de la distance a la paroi dans le',/,&
'@    fichier suite).                                         ',/,&
'@  Si cette initialisation pose probleme (modification du    ',/,&
'@    nombre et de la position des faces de paroi depuis le   ',/,&
'@    calcul precedent), forcer ICDPAR=-1 dans usini1 (il peut',/,&
'@    alors y avoir un leger decalage dans la viscosite       ',/,&
'@    turbulente au premier pas de temps).                    ',/,&
'@                                                            ',/,&
'@  Le calcul sera execute.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  Le modele de turbulence k-omega a ete choisi, avec        ',/,&
'@    l''option de recalcul de la distance a la paroi         ',/,&
'@    (ICDPAR=-1). Il se peut qu''il y ait un leger decalage  ',/,&
'@    dans la viscosite turbulente au premier pas de temps.   ',/,&
'@                                                            ',/,&
'@  Le calcul sera execute.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    PHASE ',   I10,' ISTMPF = ',   I10                       ,/,&
'@    THETFL WILL BE AUTOMATICALLY INITIALIZED.               ',/,&
'@    DO NOT MODIFY IT IN usini1.                             ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usini1.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1011 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    PHASE ',   I10,' ',A6,' = ',   I10                       ,/,&
'@    ',A6,' WILL BE INITIALIZED AUTOMATICALLY                ',/,&
'@    DO NOT MODIFY IT IN usini1.                             ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usini1.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1021 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    SCALAR ',   I10,' ',A6,' = ',   I10                      ,/,&
'@    ',A6,' WILL BE INITIALIZED AUTOMATICALLY                ',/,&
'@    DO NOT MODIFY IT IN usini1.                             ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usini1.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1031 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    PHASE ',   I10,' ',A17                                   ,/,&
'@    ',A6,' WILL BE INITIALIZED AUTOMATICALLY                ',/,&
'@    DO NOT MODIFY IT IN usini1.                             ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usini1.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1032 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    MESH VELOCITY IN ALE                                    ',/,&
'@    ',A6,' WILL BE INITIALIZED AUTOMATICALLY                ',/,&
'@    DO NOT MODIFY IT IN usini1.                             ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usini1.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1041 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    PHASE ',   I10,' ',A8,' ',I10                            ,/,&
'@    ',A6,' WILL BE INITIALIZED AUTOMATICALLY                ',/,&
'@    DO NOT MODIFY IT IN usini1.                             ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usini1.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1061 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    ICALHY must be an integer equal to 0 or 1               ',/,&
'@                                                            ',/,&
'@  Its value is ',I10                                         ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usini1.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1071 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    SCALAR ',I10   ,' DO NOT MODIFY THE DIFFUSIVITY         ',/,&
'@                                                            ',/,&
'@  The scalar ',I10   ,' is the fluctuations variance        ',/,&
'@    of the scalar ',I10                                      ,/,&
'@                             (ISCAVR(',I10   ,') = ',I10     ,/,&
'@  The diffusivity VISLS0(',I10   ,') of the scalar ',I10     ,/,&
'@    must not be set:                                        ',/,&
'@    it will be automatically set equal to the scalar        ',/,&
'@    diffusivity ',I10   ,' i.e. ',E14.5                      ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usini1.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING:       IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@  The k-omega turbulence model has been chosen. In order to ',/,&
'@    have a correct calculation restart, the ICDPAR indicator',/,&
'@    has been set to 1 (read the wall distance in the restart',/,&
'@    file).                                                  ',/,&
'@  If this initialization raises any issue (modification of  ',/,&
'@    the number and position of the wall faces since the     ',/,&
'@    previous calcuation), force ICDPAR=1 in usini1 (there   ',/,&
'@    maight be a small shift in the turbulent viscosity at   ',/,&
'@    the first time-step).                                   ',/,&
'@                                                            ',/,&
'@  The calculation will be run.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING:       IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@  The k-omega turbulence model has been chosen, with the    ',/,&
'@    option for a re-calculation of the wall distance        ',/,&
'@    (ICDPAR=-1). There might be a small shift in the        ',/,&
'@    turbulent viscosity at the first time-step.             ',/,&
'@                                                            ',/,&
'@  The calculation will be run.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

return
end subroutine
