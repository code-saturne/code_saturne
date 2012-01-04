!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

integer          ii, jj, ivar, iok, iest, imom, ikw
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
ipp = ipppro(ipproc(irom))
if(                       ichrvr(ipp).eq.-999) ichrvr(ipp) = 1
ipp = ipppro(ipproc(ivisct))
if( (iturb.eq.10 .or. itytur.eq.2                 &
     .or. itytur.eq.5 .or. iturb.eq.60            &
     .or. iturb.eq.70 )                                  &
     .and.ichrvr(ipp).eq.-999) ichrvr(ipp) = 1
if (idtvar.lt.0) then
  ichrvr(ipppro(ipproc(icour))) = 0
  ichrvr(ipppro(ipproc(ifour))) = 0
endif
do iest = 1, nestmx
  if(iescal(iest).gt.0) then
    ipp = ipppro(ipproc(iestim(iest)))
    if(                     ichrvr(ipp).eq.-999) ichrvr(ipp) = 1
  endif
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
ipp = ipppro(ipproc(ivisct))
if( (iturb.eq.10 .or. itytur.eq.2                 &
     .or. itytur.eq.5 .or. iturb.eq.60            &
     .or. iturb.eq.70 )                                  &
     .and.ihisvr(ipp,1).eq.-999) ihisvr(ipp,1) = -1
if (idtvar.lt.0) then
  ihisvr(ipppro(ipproc(icour)),1) = 0
  ihisvr(ipppro(ipproc(ifour)),1) = 0
endif
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

IF(NOMVAR(IPPRTP(IPR   )) .EQ.' ') THEN
  NOMVAR(IPPRTP(IPR   )) = 'Pres'
endif
IF(NOMVAR(IPPRTP(IU    )) .EQ.' ') THEN
  NOMVAR(IPPRTP(IU    )) = 'VitesX'
endif
IF(NOMVAR(IPPRTP(IV    )) .EQ.' ') THEN
  NOMVAR(IPPRTP(IV    )) = 'VitesY'
endif
IF(NOMVAR(IPPRTP(IW    )) .EQ.' ') THEN
  NOMVAR(IPPRTP(IW    )) = 'VitesZ'
endif
if(itytur.eq.2) then
  IF(NOMVAR(IPPRTP(IK    )) .EQ.' ') THEN
    NOMVAR(IPPRTP(IK    )) = 'EnTurb'
  endif
  IF(NOMVAR(IPPRTP(IEP   )) .EQ.' ') THEN
    NOMVAR(IPPRTP(IEP   )) = 'Dissip'
  endif
elseif(itytur.eq.3) then
  IF(NOMVAR(IPPRTP(IR11  )) .EQ.' ') THEN
    NOMVAR(IPPRTP(IR11  )) =  'R11'
  endif
  IF(NOMVAR(IPPRTP(IR22  )) .EQ.' ') THEN
    NOMVAR(IPPRTP(IR22  )) = 'R22'
  endif
  IF(NOMVAR(IPPRTP(IR33  )) .EQ.' ') THEN
    NOMVAR(IPPRTP(IR33  )) = 'R33'
  endif
  IF(NOMVAR(IPPRTP(IR12  )) .EQ.' ') THEN
    NOMVAR(IPPRTP(IR12  )) = 'R12'
  endif
  IF(NOMVAR(IPPRTP(IR13  )) .EQ.' ') THEN
    NOMVAR(IPPRTP(IR13  )) = 'R13'
  endif
  IF(NOMVAR(IPPRTP(IR23  )) .EQ.' ') THEN
    NOMVAR(IPPRTP(IR23  )) = 'R23'
  endif
  IF(NOMVAR(IPPRTP(IEP   )) .EQ.' ') THEN
    NOMVAR(IPPRTP(IEP   )) = 'Dissip'
  endif
elseif(itytur.eq.5) then
  IF(NOMVAR(IPPRTP(IK    )) .EQ.' ') THEN
    NOMVAR(IPPRTP(IK    )) = 'EnTurb'
  endif
  IF(NOMVAR(IPPRTP(IEP   )) .EQ.' ') THEN
    NOMVAR(IPPRTP(IEP   )) = 'Dissip'
  endif
  IF(NOMVAR(IPPRTP(IPHI  )) .EQ.' ') THEN
    NOMVAR(IPPRTP(IPHI  )) = 'phi'
  endif
  if(iturb.eq.50) then
    IF(NOMVAR(IPPRTP(IFB   )) .EQ.' ') THEN
      NOMVAR(IPPRTP(IFB   )) = 'fbarre'
    endif
  elseif(iturb.eq.51) then
    IF(NOMVAR(IPPRTP(IAL   )) .EQ.' ') THEN
      NOMVAR(IPPRTP(IAL   )) = 'Alpha'
    endif
  endif
elseif(iturb.eq.60) then
  IF(NOMVAR(IPPRTP(IK    )) .EQ.' ') THEN
    NOMVAR(IPPRTP(IK    )) = 'EnTurb'
  endif
  IF(NOMVAR(IPPRTP(IOMG  )) .EQ.' ') THEN
    NOMVAR(IPPRTP(IOMG  )) = 'Omega'
  endif
elseif(iturb.eq.70) then
  IF(NOMVAR(IPPRTP(INUSA )) .EQ.' ') THEN
    NOMVAR(IPPRTP(INUSA )) = 'NuTild'
  endif
endif

IF(NOMVAR(IPPPRO(IPPROC(IROM  ))) .EQ.' ') THEN
  NOMVAR(IPPPRO(IPPROC(IROM  ))) = 'MasVol'
endif
IF(NOMVAR(IPPPRO(IPPROC(IVISCT))) .EQ.' ') THEN
  NOMVAR(IPPPRO(IPPROC(IVISCT))) = 'VisTur'
endif
IF(NOMVAR(IPPPRO(IPPROC(IVISCL))) .EQ.' ') THEN
  NOMVAR(IPPPRO(IPPROC(IVISCL))) = 'VisMol'
endif
if (ismago.gt.0) then
  IF(NOMVAR(IPPPRO(IPPROC(ISMAGO))) .EQ.' ') THEN
    NOMVAR(IPPPRO(IPPROC(ISMAGO))) = 'Csdyn2'
  endif
endif
if(icp   .gt.0) then
  IF(NOMVAR(IPPPRO(IPPROC(ICP   ))) .EQ.' ') THEN
    NOMVAR(IPPPRO(IPPROC(ICP   ))) = 'ChalSp'
  endif
endif
if(iescal(iespre).gt.0) then
  ipp = ipppro(ipproc(iestim(iespre)))
  IF(NOMVAR(IPP) .EQ.' ') THEN
    WRITE(NOMVAR(IPP),'(A5,I1)') 'EsPre',IESCAL(IESPRE)
  endif
endif
if(iescal(iesder).gt.0) then
  ipp = ipppro(ipproc(iestim(iesder)))
  IF(NOMVAR(IPP) .EQ.' ') THEN
    WRITE(NOMVAR(IPP),'(A5,I1)') 'EsDer',IESCAL(IESDER)
  endif
endif
if(iescal(iescor).gt.0) then
  ipp = ipppro(ipproc(iestim(iescor)))
  IF(NOMVAR(IPP) .EQ.' ') THEN
    WRITE(NOMVAR(IPP),'(A5,I1)') 'EsCor',IESCAL(IESCOR)
  endif
endif
if(iescal(iestot).gt.0) then
  ipp = ipppro(ipproc(iestim(iestot)))
  IF(NOMVAR(IPP) .EQ.' ') THEN
    WRITE(NOMVAR(IPP),'(A5,I1)') 'EsTot',IESCAL(IESTOT)
  endif
endif

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

ipp = ipppro(ipproc(irom  ))
if(irovar.eq.1.and.ilisvr(ipp).eq.-999) ilisvr(ipp) = 1
ipp = ipppro(ipproc(ivisct))
if( (iturb.eq.10 .or. itytur.eq.2                 &
     .or. itytur.eq.5 .or. iturb.eq.60            &
     .or. iturb.eq.70 )                                  &
     .and.ilisvr(ipp).eq.-999) ilisvr(ipp) = 1
if (inusa .gt. 0) then
  ipp = ipppro(ipproc(inusa))
  if(iturb.eq.70.and.ilisvr(ipp).eq.-999) ilisvr(ipp) = 1
endif
ipp = ipppro(ipproc(icour))
if (ilisvr(ipp).eq.-999 .or. idtvar.lt.0) ilisvr(ipp) = 0
ipp = ipppro(ipproc(ifour))
if (ilisvr(ipp).eq.-999 .or. idtvar.lt.0) ilisvr(ipp) = 0
if(iescal(iespre).gt.0) then
  ipp = ipppro(ipproc(iestim(iespre)))
  if(                     ilisvr(ipp).eq.-999) ilisvr(ipp) = 1
endif
if(iescal(iesder).gt.0) then
  ipp = ipppro(ipproc(iestim(iesder)))
  if(                     ilisvr(ipp).eq.-999) ilisvr(ipp) = 1
endif
if(iescal(iescor).gt.0) then
  ipp = ipppro(ipproc(iestim(iescor)))
  if(                     ilisvr(ipp).eq.-999) ilisvr(ipp) = 1
endif
if(iescal(iestot).gt.0) then
  ipp = ipppro(ipproc(iestim(iestot)))
  if(                     ilisvr(ipp).eq.-999) ilisvr(ipp) = 1
endif

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

!   -- Flux de masse
if(abs(thetfl+999.d0).gt.epzero) then
  write(nfecra,1001) istmpf
  iok = iok + 1
elseif(istmpf.eq.0) then
  thetfl = 0.d0
elseif(istmpf.eq.2) then
  thetfl = 0.5d0
endif

!    -- Proprietes physiques
if(abs(thetro+999.d0).gt.epzero) then
  WRITE(NFECRA,1011) 'IROEXT',IROEXT,'THETRO'
  iok = iok + 1
elseif(iroext.eq.0) then
  thetro = 0.0d0
elseif(iroext.eq.1) then
  thetro = 0.5d0
elseif(iroext.eq.2) then
  thetro = 1.d0
endif
if(abs(thetvi+999.d0).gt.epzero) then
  WRITE(NFECRA,1011) 'IVIEXT',IVIEXT,'THETVI'
  iok = iok + 1
elseif(iviext.eq.0) then
  thetvi = 0.0d0
elseif(iviext.eq.1) then
  thetvi = 0.5d0
elseif(iviext.eq.2) then
  thetvi = 1.d0
endif
if(abs(thetcp+999.d0).gt.epzero) then
  WRITE(NFECRA,1011) 'ICPEXT',ICPEXT,'THETCP'
  iok = iok + 1
elseif(icpext.eq.0) then
  thetcp = 0.0d0
elseif(icpext.eq.1) then
  thetcp = 0.5d0
elseif(icpext.eq.2) then
  thetcp = 1.d0
endif

!    -- Termes sources NS
if(abs(thetsn+999.d0).gt.epzero) then
  WRITE(NFECRA,1011) 'ISNO2T',ISNO2T,'THETSN'
  iok = iok + 1
elseif(isno2t.eq.1) then
  thetsn = 0.5d0
elseif(isno2t.eq.2) then
  thetsn = 1.d0
elseif(isno2t.eq.0) then
  thetsn = 0.d0
endif

!    -- Termes sources grandeurs turbulentes
if(abs(thetst+999.d0).gt.epzero) then
  WRITE(NFECRA,1011) 'ISTO2T',ISTO2T,'THETST'
  iok = iok + 1
elseif(isto2t.eq.1) then
  thetst = 0.5d0
elseif(isto2t.eq.2) then
  thetst = 1.d0
elseif(isto2t.eq.0) then
  thetst = 0.d0
endif

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

!     Vitesse pression (la pression est prise sans interp)
if(abs(thetav(iu )+999.d0).gt.epzero.or.                 &
     abs(thetav(iv )+999.d0).gt.epzero.or.                 &
     abs(thetav(iw )+999.d0).gt.epzero.or.                 &
     abs(thetav(ipr)+999.d0).gt.epzero) then
  WRITE(NFECRA,1031) 'VITESSE-PRESSION ','THETAV'
  iok = iok + 1
elseif(ischtp.eq.1) then
  thetav(iu ) = 1.d0
  thetav(iv ) = 1.d0
  thetav(iw ) = 1.d0
  thetav(ipr) = 1.d0
elseif(ischtp.eq.2) then
  thetav(iu ) = 0.5d0
  thetav(iv ) = 0.5d0
  thetav(iw ) = 0.5d0
  thetav(ipr) = 1.d0
endif

!     Turbulence (en k-eps : ordre 1)
if(itytur.eq.2) then
  if(abs(thetav(ik )+999.d0).gt.epzero.or.               &
       abs(thetav(iep)+999.d0).gt.epzero) then
    WRITE(NFECRA,1031) 'VARIABLES   K-EPS','THETAV'
    iok = iok + 1
  elseif(ischtp.eq.1) then
    thetav(ik ) = 1.d0
    thetav(iep) = 1.d0
  elseif(ischtp.eq.2) then
    !     pour le moment, on ne peut pas passer par ici (cf varpos)
    thetav(ik ) = 0.5d0
    thetav(iep) = 0.5d0
  endif
elseif(itytur.eq.3) then
  if(abs(thetav(ir11)+999.d0).gt.epzero.or.              &
       abs(thetav(ir22)+999.d0).gt.epzero.or.              &
       abs(thetav(ir33)+999.d0).gt.epzero.or.              &
       abs(thetav(ir12)+999.d0).gt.epzero.or.              &
       abs(thetav(ir13)+999.d0).gt.epzero.or.              &
       abs(thetav(ir23)+999.d0).gt.epzero.or.              &
       abs(thetav(iep )+999.d0).gt.epzero) then
    WRITE(NFECRA,1031) 'VARIABLES  RIJ-EP','THETAV'
    iok = iok + 1
  elseif(ischtp.eq.1) then
    thetav(ir11) = 1.d0
    thetav(ir22) = 1.d0
    thetav(ir33) = 1.d0
    thetav(ir12) = 1.d0
    thetav(ir13) = 1.d0
    thetav(ir23) = 1.d0
    thetav(iep ) = 1.d0
  elseif(ischtp.eq.2) then
    thetav(ir11) = 0.5d0
    thetav(ir22) = 0.5d0
    thetav(ir33) = 0.5d0
    thetav(ir12) = 0.5d0
    thetav(ir13) = 0.5d0
    thetav(ir23) = 0.5d0
    thetav(iep ) = 0.5d0
  endif
elseif(iturb.eq.50) then
  if(abs(thetav(ik  )+999.d0).gt.epzero.or.              &
       abs(thetav(iep )+999.d0).gt.epzero.or.              &
       abs(thetav(iphi)+999.d0).gt.epzero.or.              &
       abs(thetav(ifb )+999.d0).gt.epzero) then
    WRITE(NFECRA,1031) 'VARIABLES     V2F','THETAV'
    iok = iok + 1
  elseif(ischtp.eq.1) then
    thetav(ik  ) = 1.d0
    thetav(iep ) = 1.d0
    thetav(iphi) = 1.d0
    thetav(ifb ) = 1.d0
  elseif(ischtp.eq.2) then
    !     pour le moment, on ne peut pas passer par ici (cf varpos)
    thetav(ik  ) = 0.5d0
    thetav(iep ) = 0.5d0
    thetav(iphi) = 0.5d0
    thetav(ifb ) = 0.5d0
  endif
elseif(iturb.eq.51) then
  if(abs(thetav(ik  )+999.d0).gt.epzero.or.              &
       abs(thetav(iep )+999.d0).gt.epzero.or.              &
       abs(thetav(iphi)+999.d0).gt.epzero.or.              &
       abs(thetav(ial )+999.d0).gt.epzero) then
    WRITE(NFECRA,1031) 'VARIABLES BL-V2/K','THETAV'
    iok = iok + 1
  elseif(ischtp.eq.1) then
    thetav(ik  ) = 1.d0
    thetav(iep ) = 1.d0
    thetav(iphi) = 1.d0
    thetav(ial ) = 1.d0
  elseif(ischtp.eq.2) then
    !     pour le moment, on ne peut pas passer par ici (cf varpos)
    thetav(ik  ) = 0.5d0
    thetav(iep ) = 0.5d0
    thetav(iphi) = 0.5d0
    thetav(ial ) = 0.5d0
  endif
elseif(iturb.eq.60) then
  if(abs(thetav(ik  )+999.d0).gt.epzero.or.              &
       abs(thetav(iomg)+999.d0).gt.epzero ) then
    WRITE(NFECRA,1031) 'VARIABLES K-OMEGA','THETAV'
    iok = iok + 1
  elseif(ischtp.eq.1) then
    thetav(ik  ) = 1.d0
    thetav(iomg) = 1.d0
  elseif(ischtp.eq.2) then
    !     pour le moment, on ne peut pas passer par ici (cf varpos)
    thetav(ik  ) = 0.5d0
    thetav(iomg) = 0.5d0
  endif
elseif(iturb.eq.70) then
  if(abs(thetav(inusa)+999.d0).gt.epzero) then
    WRITE(NFECRA,1031) 'VARIABLE NU_tilde de SA','THETAV'
    iok = iok + 1
  elseif(ischtp.eq.1) then
    thetav(inusa) = 1.d0
  elseif(ischtp.eq.2) then
    !     pour le moment, on ne peut pas passer par ici (cf varpos)
    thetav(inusa) = 0.5d0
  endif
endif

!     Scalaires
do iscal = 1, nscal
  ivar  = isca(iscal)
  if(abs(thetav(ivar)+999.d0).gt.epzero) then
    WRITE(NFECRA,1041) 'SCALAIRE',ISCAL,'THETAV'
    iok = iok + 1
  elseif(ischtp.eq.1) then
    thetav(ivar) = 1.d0
  elseif(ischtp.eq.2) then
    thetav(ivar) = 0.5d0
  endif
enddo

!     Vitesse de maillage en ALE
if (iale.eq.1) then
  if(abs(thetav(iuma)+999.d0).gt.epzero.or.                       &
     abs(thetav(ivma)+999.d0).gt.epzero.or.                       &
     abs(thetav(iwma)+999.d0).gt.epzero) then
    WRITE(NFECRA,1032) 'THETAV'
    iok = iok + 1
  elseif(ischtp.eq.1) then
    thetav(iuma) = 1.d0
    thetav(ivma) = 1.d0
    thetav(iwma) = 1.d0
  elseif(ischtp.eq.2) then
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

if(itytur.eq.4) then
  ii = iu
  if(isstpc(ii).eq.-999) isstpc(ii) = 1
  ii = iv
  if(isstpc(ii).eq.-999) isstpc(ii) = 1
  ii = iw
  if(isstpc(ii).eq.-999) isstpc(ii) = 1
  do jj = 1, nscal
    ii = isca(jj)
    if(isstpc(ii).eq.-999) isstpc(ii) = 0
  enddo
endif

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

ii = iu
if(abs(blencv(ii)+999.d0).lt.epzero) blencv(ii) = 1.d0
ii = iv
if(abs(blencv(ii)+999.d0).lt.epzero) blencv(ii) = 1.d0
ii = iw
if(abs(blencv(ii)+999.d0).lt.epzero) blencv(ii) = 1.d0
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

if(ischtp.eq.2) then
  ii = ipr
  if(nswrsm(ii).eq.-999) nswrsm(ii) = 5
  if(abs(epsrsm(ii)+999.d0).lt.epzero) epsrsm(ii) = 1.d-5
  if(abs(epsilo(ii)+999.d0).lt.epzero) epsilo(ii) = 1.d-5
  ii = iu
  if(nswrsm(ii).eq.-999) nswrsm(ii) = 10
  if(abs(epsrsm(ii)+999.d0).lt.epzero) epsrsm(ii) = 1.d-5
  if(abs(epsilo(ii)+999.d0).lt.epzero) epsilo(ii) = 1.d-5
  ii = iv
  if(nswrsm(ii).eq.-999) nswrsm(ii) = 10
  if(abs(epsrsm(ii)+999.d0).lt.epzero) epsrsm(ii) = 1.d-5
  if(abs(epsilo(ii)+999.d0).lt.epzero) epsilo(ii) = 1.d-5
  ii = iw
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
ii = ipr
if(nswrsm(ii).eq.-999) nswrsm(ii) = 2

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

!     Ici, ce n'est pas grave pour le moment,
!      etant entendu que ces coefs ne servent pas
!      s'ils servaient, attention dans le cas a plusieurs phases avec
!      une seule pression : celle ci prend le coef de la derniere phase
cdtvar(iv ) = cdtvar(iu)
cdtvar(iw ) = cdtvar(iu)
cdtvar(ipr) = cdtvar(iu)

if(itytur.eq.2) then
  cdtvar(iep ) = cdtvar(ik  )
elseif(itytur.eq.3) then
  cdtvar(ir22) = cdtvar(ir11)
  cdtvar(ir33) = cdtvar(ir11)
  cdtvar(ir12) = cdtvar(ir11)
  cdtvar(ir13) = cdtvar(ir11)
  cdtvar(ir23) = cdtvar(ir11)
  cdtvar(iep ) = cdtvar(ir11)
elseif(itytur.eq.5) then
  cdtvar(iep ) = cdtvar(ik  )
  cdtvar(iphi) = cdtvar(ik  )
!     CDTVAR(IFB/IAL) est en fait inutile car pas de temps dans
!     l'eq de f_barre/alpha
  if(iturb.eq.50) then
    cdtvar(ifb ) = cdtvar(ik  )
  elseif(iturb.eq.51) then
    cdtvar(ial ) = cdtvar(ik  )
  endif
elseif(iturb.eq.60) then
  cdtvar(iomg) = cdtvar(ik  )
elseif(iturb.eq.70) then
  ! cdtvar est à 1.0 par defaut dans iniini.f90
  cdtvar(inusa)= cdtvar(inusa)
endif

! ---> IDEUCH, YPLULI
!      En laminaire, longueur de melange, Spalar-Allmaras et LES,
!      une echelle de vitesse.
!      Sinon, 2 echelles, sauf si l'utilisateur choisit 1 echelle.
!      On a initialise IDEUCH a -999 pour voir si l'utilisateur essaye
!        de choisir deux echelles quand ce n'est pas possible et le
!        prevenir dans la section verification.

if(ideuch.eq.-999) then
  if(iturb.eq. 0.or.                                     &
       iturb.eq.10.or.                                     &
       itytur.eq.4.or.                                     &
       iturb.eq.70) then
    ideuch = 0
  else
    ideuch = 1
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

if (ypluli.lt.-grand) then
  if (ideuch.eq.2 .or. itytur.eq.4 .or.           &
       iturb.eq.70    ) then
    ypluli = 10.88d0
  else
    ypluli = 1.d0/xkappa
  endif
endif


! ---> Van Driest
if(idries.eq.-1) then
  !   On met 1 en supposant qu'en periodicite ou parallele on utilise le
  !     mode de calcul de la distance a la paroi qui les prend en charge
  !     (ICDPAR=+/-1, valeur par defaut)
  if(iturb.eq.40) then
    idries = 1
  elseif(iturb.eq.41) then
    idries = 0
  elseif(iturb.eq.42) then
    idries = 0
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

!         On couple le scalaire temperature de la phase
      if(iscalt.gt.0.and.iscalt.le.nscal) then
        icpsyr(iscalt) = 1
        goto 100
      endif
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
      if(ii.ne.iscalt) then
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
if (iturb.eq.60) ikw = 1
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
if((iturb.eq.30.and.irijec.eq.1).or.              &
     (itytur.eq.4.and.idries.eq.1).or.              &
     iturb.eq.60.or.iturb.eq.70      ) then
  ineedy = 1
endif

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
  iwarny = iwarni(iu)
endif


! ---> IKECOU
!     En k-eps prod lin, v2f ou k-omega, on met IKECOU a 0 par defaut,
!     sinon on le laisse a 1
!     Dans verini on bloquera le v2f et le k-eps prod lin si IKECOU.NE.0
!     On bloquera aussi le stationnaire si IKECOU.NE.0
if (ikecou.eq.-999) then
  if (idtvar.lt.0) then
    ikecou = 0
  else if (iturb.eq.21 .or. itytur.eq.5           &
       .or. iturb.eq.60 ) then
    ikecou = 0
  else
    ikecou = 1
  endif
endif

! ---> RELAXV
if (idtvar.lt.0) then
  relxsp = 1.d0-relxst
  if (relxsp.le.epzero) relxsp = relxst
  if (abs(relaxv(ipr)+999.d0).le.epzero)                 &
       relaxv(ipr) = relxsp
  do ii = 1, nvarmx
    if (abs(relaxv(ii)+999.d0).le.epzero) relaxv(ii) = relxst
  enddo
else
  if ( ikecou.eq.0) then
    if (itytur.eq.2 .or. itytur.eq.5) then
      if (abs(relaxv(ik)+999.d0).lt.epzero)              &
           relaxv(ik) = 0.7d0
      if (abs(relaxv(iep)+999.d0).lt.epzero)             &
           relaxv(iep) = 0.7d0
    else if (itytur.eq.6) then
      if (abs(relaxv(ik)+999.d0).lt.epzero)              &
           relaxv(ik) = 0.7d0
      if (abs(relaxv(iomg)+999.d0).lt.epzero)            &
           relaxv(iomg) = 0.7d0
    endif
  endif
  if(iturb.eq.70) then
    if(abs(relaxv(inusa)+999.d0).lt.epzero) then
      relaxv(inusa) = 1.D0
    endif
  endif
  if (abs(relaxv(ipr)+999.d0).lt.epzero)                 &
       relaxv(ipr) = 1.d0
endif

! ---> SPECIFIQUE STATIONNAIRE
if (idtvar.lt.0) then
  dtref = 1.d0
  dtmin = 1.d0
  dtmax = 1.d0
  do ii = 1, nvarmx
    istat(ii) = 0
  enddo
  arak = arak/max(relaxv(iu),epzero)
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

do ii = 1, 3
  if (xyzp0(ii).gt.-0.5d0*rinfin) then
    ixyzp0 = 1
  else
    xyzp0(ii) = 0.d0
  endif
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
      call pstdfm
      !==========
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
'@    ISTMPF = ',   I10                                        ,/,&
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
'@    ',A6,' = ',   I10                                        ,/,&
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
'@    ',A17                                                    ,/,&
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
'@    ',A8,' ',I10                                             ,/,&
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
'@    ISTMPF = ',   I10                                        ,/,&
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
'@    ',A6,' = ',   I10                                        ,/,&
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
'@    ',A17                                                    ,/,&
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
'@    ',A8,' ',I10                                             ,/,&
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
