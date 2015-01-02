!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

subroutine rayir &
!===============

 ( k1,kmray,ico2,emis,                            &
   qqv, qqqv, qqvinf, zqq,                        &
   acinfe, dacinfe, aco2, daco2, acsup, dacsup,   &
   zray,temray,qvray,qlray,fnerir,                &
   romray,preray,aeroso,foir,rayi )

!==============================================================================
!  Purpose:
!  --------

!    Atmospheric module subroutine.

!
!   calcul dans le domaine spectral infra-rouge:
!      - du profil de la divergence du flux,
!      - du flux descendant au sol.
!
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.__________________________________________________.
! !    nom    !type!mode!                   role                                 !
!_!___________!____!____!________________________________________________________!
! !  k1       ! e  ! d  ! indice du premier point du segment vertical            !
! !           !    !    ! considere                                              !
! !  kmray    ! e  ! d  ! nombre de niveaux verticaux pour les modules           !
! !           !    !    ! de rayonnement                                         !
! !  ico2     ! e  ! d  ! ico2=1 -> calcul fct d'absorption par co2              !
! !  emis     ! r  ! d  ! emissivite de la surface terrestre                     !
! !  qqv      ! tr ! d  ! absorption par la vapeur d'eau + dimeres               !
! !  qqqv     ! tr ! d  !   idem niveaux intermediaires                          !
! !  qqvinf   ! tr ! d  !   idem mais contribution > 11000m                      !
! !  zzq      ! tr ! d  ! vertical coordinate                                    !
! !  acinfe   ! tr ! d  ! absorption par CO2 + 03                                !
! !  dacinfe  ! tr ! d  ! absorption differentielle par CO2 + 03                 !
! !  aco2     ! tr ! d  ! idem pour CO2 seul                                     !
! !  daco2    ! tr ! d  ! idem pour CO2 seul                                     !
! !  acsup    ! tr ! d  ! idem acinfe, flux descendant                           !
! !  dacsup   ! tr ! d  ! idem acinfe, flux descendant                           !
! !  zray     ! tr ! d  ! altitude (maillage physique)                           !
! !  temray   ! tr ! d  ! temperature en Celsius                                 !
! !  qvray    ! tr ! d  ! humidite                                               !
! !  qlray    ! tr ! d  ! teneur en eau liquide                                  !
! !  fnerir   ! tr ! d  ! nebulosite                                             !
! !  romray   ! tr ! d  ! masse volumique                                        !
! !  preray   ! tr ! d  ! pression sur maillage vitesse                          !
! !  aeroso   ! tr ! d  ! contenu en aerosol                                     !
! !  foir     ! r  ! r  ! flux IR descendant au sol                              !
! !  rayi     ! tr ! r  ! divergence du flux IR                                  !
!_!___________!____!____!________________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use atincl, only: kmx, rvsra, cpvcpa

!===============================================================================

implicit none

! Arguments

integer k1,kmray,ico2
double precision emis
double precision qqv(kmx+1), qqqv(kmx+1), qqvinf, zqq(kmx)
double precision acinfe(kmx), dacinfe(kmx), aco2(kmx,kmx), daco2(kmx,kmx)
double precision acsup(kmx), dacsup(kmx)
double precision temray(kmx),romray(kmx),preray(kmx)
double precision qvray(kmx),qlray(kmx),fnerir(kmx)
double precision zray(kmx)
double precision aeroso(kmx)
double precision foir,rayi(kmx)

!Local variables

integer n,k,kk,kp1,ineb,inua,iaer
double precision sig
double precision qqcinf,cetyp,dzs8
double precision xqqvinf,xqqcinf,xqqlinf,fo,t4zt
double precision tauv,dtauv,taul,ctray,a1
double precision qqqqv,qqqqc,qqqql,a2,tvinfe
double precision dtvinfe,tvsup,dtvsup,dul
double precision t41,tlinfe,tlsup,qqlinf
double precision fn,fns,fni,u,tco2,zz,zzk
double precision pzp0,zqq0,corp
double precision beta,wh2ol
double precision zbas
double precision caero

double precision,allocatable :: rov(:),roc(:),rol(:),qv0(:),qc(:)
double precision,allocatable :: qqc(:),qql(:)
double precision,allocatable :: qqqc(:),qqql(:)
double precision,allocatable :: pspo(:),pspoqq(:),dt4dz(:),tqq(:)
double precision,allocatable :: dz0(:)
double precision,allocatable :: kliq(:)

!===============================================================================

allocate(rov(kmx),roc(kmx),rol(kmx),qv0(kmx),qc(kmx))
allocate(qqc(kmx),qql(kmx))
allocate(qqqc(kmx),qqql(kmx))
allocate(pspo(kmx),pspoqq(kmx),dt4dz(kmx),tqq(kmx))
allocate(dz0(kmx))
allocate(kliq(kmx))

!   0- initialisations locales
!   --------------------------
xqqlinf = 0.d0
taul = 0.d0
sig = stephn ! Bolztman constant 5.6703d-8

do n = 1, kmx
  rov(n) = 0.d0
  roc(n) = 0.d0
  rol(n) = 0.d0
  qv0(n) = 0.d0
  qc (n) = 0.d0
  qqv(n) = 0.d0
  qqc(n) = 0.d0
  qql(n) = 0.d0
  qqqv(n) = 0.d0
  qqqc(n) = 0.d0
  qqql(n) = 0.d0
  pspo(n) = 0.d0
  pspoqq(n) = 0.d0
  dt4dz(n) = 0.d0
  tqq(n) = 0.d0
  rayi(n) = 0.d0
enddo
qqv(kmx+1) = 0.d0
qqqv(kmx+1) = 0.d0

foir = 0.d0

! contribution des couches superieures de l atmosphere (11000 a 44000)
! aux quantites d'absorbant

qqvinf = 0.002515d0
qqcinf = 2.98d-7
qqlinf = 0.d0

cetyp = rvsra/1.134d0

! coefficient pour prendre en compte l'ensemble des angles de diffusion
! (deja pris en compte dans les fonctions d'emissivite)
beta = 1.66d0

! definition du coefficient d'absorption par l'eau nuageuse :
! (on prend le meme rayon equivalent que dans le rayonnement solaire)
!
do k = k1, kmray
! densite de l'eau liquide en g/m3 dans les couches considerees
  wh2ol = 1.d3*romray(k)*qlray(k)
! coefficient d'absorption suppose constant pour des rayons de
! particules allant de 0.01 a 10 microns
  kliq(k) = 120.d0
enddo

inua = 0
iaer = 0
caero = 1.d-9

do k = k1, kmray-1
 dz0(k) = zray(k+1) - zray(k)
enddo

dz0(kmray) = 0.d0
zbas = zray(k1)

!   1- calcul des quantites d'absorbant pour la vapeur d'eau et les dime
!   --------------------------------------------------------------------
do k = k1, kmray
  if(qlray(k).gt.1.d-8) inua = 1
  if(aeroso(k).gt.1.d-8) iaer = 1
  pspo(k) = preray(k)/preray(k1)
  corp = pspo(k)*preray(k1)/101300.d0
  qv0(k) = qvray(k)*corp*sqrt(tkelvi/(temray(k) + tkelvi))
  rov(k) = romray(k)*qv0(k)
  qc(k) = qvray(k)*qvray(k)*corp*exp(1800.d0*                         &
      (1.d0/(temray(k) + tkelvi)-1.d0/296.d0))*cetyp
  roc(k) = romray(k)*qc(k)
enddo

! en presence de couches nuageuse on calcule la quantite d'absorbant
! pour l'eau liquide
! on prend en compte les aerosols uniquement en ciel clair

if(inua.eq.1) then
  do k = k1, kmray
    rol(k) = romray(k)*(qlray(k) + caero*aeroso(k))
  enddo
elseif(iaer.eq.1) then
  do k = k1, kmray
    rol(k) = caero*romray(k)*aeroso(k)
    if(aeroso(k).gt.1.d-10) fnerir(k) = 1.d0
  enddo
endif

!  2 -calcul des chemins optiques pour la vapeur d'eau et les dimeres
!  ------------------------------------------------------------------
! les variables qqq correspondent aux niveaux standards, les variables
! qq aux niveaux intermediaires

do k = k1+1, kmray
  dzs8 = dz0(k-1)/8.d0
  tqq(k) = (temray(k)+temray(k-1))/2.d0 + tkelvi
  zqq(k) = (zray(k) + zray(k-1))/2.d0
  dt4dz(k) = 4.d0*(tqq(k)**3.d0)*(temray(k) - temray(k-1))/dz0(k-1)
  pspoqq(k) = (pspo(k) + pspo(k-1))/2.d0
  qqv(k) = qqqv(k-1) + dzs8*(rov(k) + 3.d0*rov(k-1))
  qqc(k) = qqqc(k-1) + dzs8*(roc(k) + 3.d0*roc(k-1))
  qqqv(k) = qqv(k) + dzs8*(3.d0*rov(k) + rov(k-1))
  qqqc(k) = qqc(k) + dzs8*(3.d0*roc(k) + roc(k-1))
enddo

zqq(k1) = zray(k1)
pspoqq(k1) = pspo(k1)
xqqvinf = qqqv(kmray) + qqvinf
xqqcinf = qqqc(kmray) + qqcinf

!   3 -en presence de nuages, calcul des chemins optiques pour l'eau
!      nuageuse en utilisant les memes notations
!   ----------------------------------------------------------------

if(inua.eq.1) then
  do k = k1+1,kmray
    dzs8 = dz0(k-1)/8.d0
    qql(k) = qqql(k-1) + dzs8*(rol(k) + 3.d0*rol(k-1))
    qqql(k) = qql(k) + dzs8*(3.d0*rol(k) + rol(k-1))
  enddo
  xqqlinf = qqql(kmray) + qqlinf
endif

!   4 -au premier pas de temps, on calcule les fonctions d'absorption
!      pour le CO2 et O3

if(ico2.eq.1) then

  do k = k1,kmray
    zz = zray(k)
    if(k.ne.k1) then
      do kk = k1+1, k
        u = qqqv(k) - qqv(kk)
        tco2 = (temray(k) + tkelvi + tqq(kk))/2.d0
        call rayigc(zbas,zz,pspo(k),zqq(kk),pspoqq(kk),aco2(k,kk),    &
               daco2(k,kk),qv0(k),u,tco2,romray(k))
      enddo
    endif
    kp1 = k+1
    do kk = kp1, kmray
      u = qqv(kk) - qqqv(k)
      tco2 = (temray(k) + tkelvi + tqq(kk))/2.d0
      call rayigc(zbas,zz,pspo(k),zqq(kk),pspoqq(kk),aco2(k,kk),      &
             daco2(k,kk),qv0(k),u,tco2,romray(k))
    enddo
    u = xqqvinf - qqqv(k)
    tco2 = 375.d0
    pzp0 = 0.d0
    zqq0 = 44000.d0
    call rayigc(zbas,zz,pspo(k),zqq0,pzp0,acsup(k),dacsup(k),         &
           qv0(k),u,tco2,romray(k))
    if(k.ne.k1) then
      u = qqqv(k)
      tco2 = (temray(k1)+temray(k))/2.d0 + tkelvi
      zz = zray(k1)
      zzk = zray(k)
      call rayigc(zbas,zz,pspo(k1),zzk,pspo(k),acinfe(k),dacinfe(k),  &
             qv0(k),u,tco2,romray(k))
    endif
  enddo

endif

!   5 -calcul du flux infrarouge descendant au sol
!   ----------------------------------------------

fo = 0.d0
t4zt = (temray(kmray) + tkelvi)**4

if(inua.ne.1) then

!  cas du ciel clair
!  =================
  do k = k1+1, kmray
! calcul de la transmission par la vapeur d'eau et ses dimeres
    call rayive(tauv,dtauv,qqv(k),qv0(k1),qqc(k),qc(k1),romray(k))
    fo = fo - (1.d0 - tauv + aco2(k1,k))*dt4dz(k)*dz0(k-1)
  enddo
  call rayive(tauv,dtauv,xqqvinf,qv0(k1),xqqcinf,qc(k1),            &
           romray(kmray))
  foir = sig*(fo + t4zt*(1.d0 - tauv + acsup(k1)))

else

!  cas du ciel nuageux
!  ===================
  fn = fnerir(k1+1)
  do k = k1+1, kmray
! calcul de la nebulosite pour la couche consideree. on prend le maximum
    fn = max(fn,fnerir(k))
    if(fn.lt.1.d-3) then
      fn = 0.d0
      taul = 0.d0
    else
      taul = exp(-kliq(k)*qql(k)/fn)
    endif
    call rayive(tauv,dtauv,qqv(k),qv0(k1),qqc(k),qc(k1),romray(k))
    fo = fo - (1.d0 - (1.d0 + fn*(taul - 1.d0))                                 &
       * (tauv - aco2(k1,k)))*dt4dz(k)*dz0(k-1)
  enddo
  call rayive(tauv,dtauv,xqqvinf,qv0(k1),xqqcinf,qc(k1),            &
           romray(kmray))
  foir = sig*(fo + t4zt*(1.d0 - (1.d0 + fn*(taul - 1.d0))*(tauv - acsup(k1))))
endif

!   6 -calcul du refroidissement
!   ----------------------------

t41 = (temray(k1) + tkelvi)**4
if(inua.ne.1) then

!  par ciel clair
!  ==============
  do k = k1+1, kmray

    ctray = sig/romray(k)/(cp0*(1.d0 + (cpvcpa - 1.d0)*qvray(k)))
! calcul de a1 contribution de 0 a z
    a1 = 0.d0
    do kk = k1+1, k
      qqqqv = qqqv(k) - qqv(kk)
      qqqqc = qqqc(k) - qqc(kk)
      call rayive(tauv,dtauv,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
      a1 = a1 + dt4dz(kk)*(daco2(k,kk) - dtauv)*dz0(kk-1)
    enddo
    kp1 = k+1
! calcul de a2 contribution de z a ztop
    a2 = 0.d0
    do kk = kp1, kmray
      qqqqv = qqv(kk) - qqqv(k)
      qqqqc = qqc(kk) - qqqc(k)
      call rayive(tauv,dtauv,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
      a2 = a2 + dt4dz(kk)*(daco2(k,kk) - dtauv)*dz0(kk-1)
    enddo
! contribution de z a l'infini
    qqqqv = xqqvinf - qqqv(k)
    qqqqc = xqqcinf - qqqc(k)
    call rayive(tvsup,dtvsup,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
! contribution du sol transmise par les couches inferieures (0-z)
    qqqqv = qqqv(k)
    qqqqc = qqqc(k)
    call rayive(tvinfe,dtvinfe,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
! calcul du refoidissement
    rayi(k) = ctray*(a1 - a2 + t4zt*(dacsup(k) - dtvsup) -(1.d0 - emis)         &
              *(t41 - foir/sig)*(dtvinfe - dacinfe(k)))
  enddo
else

!  par ciel nuageux
!  ================
  do k =k1+1, kmray

    ctray = sig/romray(k)/(cp0*(1.d0 + (cpvcpa - 1.d0)*qvray(k)))
! cacul de a1 contribution de 0 a z
    a1 = 0.d0
    dul = kliq(k)*rol(k)
    do kk = k1+1, k
      qqqqv = qqqv(k) - qqv(kk)
      qqqqc = qqqc(k) - qqc(kk)
      qqqql = qqql(k) - qql(kk)
      call rayive(tauv,dtauv,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
      fn = fnerir(kk)
      do ineb = kk, k
        fn = max(fn,fnerir(ineb))
      enddo
      if(fn.lt.1.d-3) then
        taul = 0.d0
        fn = 0.d0
      else
        taul = exp(-kliq(k)*qqqql/fn)
      endif
      a1 = a1 + dt4dz(kk)*((1.d0 + fn*(taul - 1.d0))*(daco2(k,kk) - dtauv)      &
         + taul*dul*(tauv - aco2(k,kk)))*dz0(kk-1)
    enddo
    kp1 = k+1
! calcul de a2 contribution de z a ztop
    a2 = 0.d0
    do kk = kp1, kmray
      qqqqv = qqv(kk) - qqqv(k)
      qqqqc = qqc(kk) - qqqc(k)
      qqqql = qql(kk) - qqql(k)
      call rayive(tauv,dtauv,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
      fn = fnerir(k)
      do ineb = k, kk
        fn = max(fn,fnerir(ineb))
      enddo
      if(fn.lt.1.d-3) then
        taul = 0.d0
        fn = 0.d0
      else
        taul = exp(-kliq(k)*qqqql/fn)
      endif
      a2 = a2 + dt4dz(kk)*((1.d0 + fn*(taul - 1.d0))*(daco2(k,kk)-dtauv)        &
         + taul*dul*(tauv - aco2(k,kk)))*dz0(kk-1)
    enddo
! contribution de z a l'infini
    qqqqv = xqqvinf - qqqv(k)
    qqqqc = xqqcinf - qqqc(k)
    qqqql = xqqlinf - qqql(k)
    call rayive(tvsup,dtvsup,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
    fns = fnerir(k)
    do ineb = k, kmray
      fns = max(fns,fnerir(ineb))
    enddo
    if(fns.lt.1.d-3) then
      tlsup = 0.d0
      fns = 0.d0
    else
      tlsup = exp(-kliq(k)*qqqql/fns)
    endif
! contribution du sol transmise par les couches inferieures (0-z)
    qqqqv = qqqv(k)
    qqqqc = qqqc(k)
    qqqql = qqql(k)
    call rayive(tvinfe,dtvinfe,qqqqv,qv0(k),qqqqc,qc(k),romray(k))
    fni = fnerir(k1+1)
    do ineb =2, k
      fni = max(fni,fnerir(ineb))
    enddo
    if(fni.lt.1.d-3) then
      tlinfe = 0.d0
      fni = 0.d0
    else
      tlinfe=exp(-kliq(k)*qqqql/fni)
    endif
! calcul du refoidissement
    rayi(k) = ctray*(a1 - a2 + t4zt*((1.d0 + fns*(tlsup - 1.d0))                &
            *(dacsup(k) - dtvsup) + dul*tlsup*(tvsup - acsup(k)))               &
            -(1.d0 - emis)*(t41 - foir/sig)*((1.d0 + fni*(tlinfe - 1.d0))       &
            *(dtvinfe - dacinfe(k)) - dul*tlinfe*(tvinfe - acinfe(k))))

  enddo

endif

deallocate(rov,roc,rol,qv0,qc)
deallocate(qqc,qql)
deallocate(qqqc,qqql)
deallocate(pspo,pspoqq,dt4dz,tqq)
deallocate(dz0)
deallocate(kliq)

return
end subroutine rayir
