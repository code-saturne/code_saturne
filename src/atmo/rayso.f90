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

subroutine rayso &
!===============

 (k1, kmray, heuray, imer1, albe,                &
  qqv, qqqv, qqvinf, zqq,                        &
  zray, qvray, qlray, fneray,                    &
  romray, aeroso, fos, rayst)

!==============================================================================
!  Purpose:
!  --------

!    Atmospheric module subroutine.

!
!    - ce programme calcule les flux solaires en atmosphere claire et
!    nuageuse suivant le schema de Lacis et Hansen 1974 (LH 74)
!    - la diffusion multiple est prise en compte par la methode
!    d'addition des couches adjacentes
!    par rapport au schema original de LH 74 ont ete ajoutes:
!    - la prise en compte d'un albedo de simple diffusion des gouttes
!      suivant Fouquart et Bonnel, 1980
!    - la prise en compte d'une nebulosite fractionnaire suivant la
!      regle de recouvrememt maximum, Geleyn et Hollingworth,1979
!
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.__________________________________________________.
! !    name   !type!mode!                   role                                 !
!_!___________!____!____!________________________________________________________!
! !  k1       ! e  ! d  ! indice du premier point du segment vertical            !
! !           !    !    ! considere                                              !
! !  kmray    ! e  ! d  ! nombre de niveaux verticaux pour les modules           !
! !           !    !    ! de rayonnement                                         !
! !  heuray   ! r  ! d  ! heure TU utilisee ds module de rayonnement             !
! !  imer1    ! e  ! d  ! indice de presence de mer                              !
! !  albe     ! r  ! m  ! albedo (modifie ds le cas ou imer=1)                   !
! !  qqv      ! tr ! d  ! epaisseur optique vapeur eau                           !
! !  qqqv     ! tr ! d  ! idem niveaux intermédiaires                            !
! !  qqvinf   ! tr ! d  ! idem mais contribution > 11000m                        !
! !  zqq      ! tr ! d  ! niveaux verticaux                                      !
! !  zray     ! tr ! d  ! altitude (maillage physique)                           !
! !  qvray    ! tr ! d  ! humidite specifique                                    !
! !  qlray    ! tr ! d  ! teneur en eau liquide                                  !
! !  fneray   ! tr ! d  ! nebulosite                                             !
! !  romray   ! tr ! d  ! masse volumique                                        !
! !  aeroso   ! tr ! d  ! contenu en aerosol                                     !
! !  fos      ! r  ! r  ! rayonnement solaire absorbe par le sol                 !
! !  rayst    ! tr ! r  ! divergence du flux de rayonnement solaire              !
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
use atincl, only: kmx, nbmett, squant, xlat, xlon, cpvcpa

!===============================================================================

implicit none

! Arguments

integer k1,kmray,imer1
double precision albe,heuray,fos
double precision qqv(kmx+1), qqqv(kmx+1), qqvinf, zqq(kmx+1)
double precision qlray(kmx),fneray(kmx),zray(kmx)
double precision qvray(kmx)
double precision aeroso(kmx)
double precision rayst(kmx),romray(kmx)

! Local variables

integer i,k,n,l,inua,k1p1,iaer
integer itop,ibase,itopp1,itopp2,ibasem1
double precision muzero,fo,rr1,m,mbar,rabar,rabar2,rbar
double precision rabarc,rabar2c,rbarc
double precision qqvtot,ym1,y,ystarm1,ystar
double precision zqm1,zq,xm1,x,xstar,xstarm1,fabs
double precision rrbar,rrbar2s,foo3,foo3c,foh2o
double precision tauctot,wh2ol,rm,req,deltaz
double precision extlnp,extlnm,pioc,zbas,dud1
double precision gasym,drt,tln,drtt1,dtrb1
double precision kn(8),pkn(8),dowtot1,dqqv
double precision zaero,reaero,piaero,gaero,caero
double precision cphum,qureel
double precision waero,taua(kmx+1),tauatot
double precision rabara,rabar2a,rbara,raero
double precision niaer,nraer
double precision s3,gama1,gama2,kt,gas,fas

double precision, allocatable:: fabsh2o(:),fabso3(:),tauc(:)
double precision, allocatable:: tau(:,:),pic(:,:),ref(:,:)
double precision, allocatable:: reft(:,:),trat(:,:),refts(:,:)
double precision, allocatable:: refs(:,:),tras(:,:),trats(:,:)
double precision, allocatable:: refb(:,:),trab(:,:),upw(:,:)
double precision, allocatable:: refbs(:,:),fabso3c(:,:),tra(:,:)
double precision, allocatable:: dow(:,:),atln(:,:),absn(:,:)
double precision, allocatable:: fnebmax(:),fneba(:)

! data pour la distribution pkn et kn

data kn/4.d-5,0.002,0.035,0.377,1.95,9.40,44.6,190./
data pkn/0.6470,0.0698,0.1443,0.0584,0.0335,0.0225,0.0158,0.0087/

!========================================================================

allocate(fabsh2o(kmx+1),fabso3(kmx+1),tauc(kmx+1))
allocate(tau(kmx+1,8),pic(kmx+1,8),ref(kmx+1,8))
allocate(reft(kmx+1,8),trat(kmx+1,8),refts(kmx+1,8))
allocate(refs(kmx+1,8),tras(kmx+1,8),trats(kmx+1,8))
allocate(refb(kmx+1,8),trab(kmx+1,8),upw(kmx+1,8))
allocate(refbs(kmx+1,8),fabso3c(kmx+1,2),tra(kmx+1,8))
allocate(dow(kmx+1,8),atln(kmx+1,8),absn(kmx+1,8))
allocate(fnebmax(kmx+1),fneba(kmx+1))

!     1 - initialisations locales
!     ===========================

inua = 0
iaer = 0
ibase = 0

do k = 1,kmray
  if(qlray(k).gt.1.d-8.or.aeroso(k).gt.1.d-8) inua = 1
  if(aeroso(k).gt.1.d-10) iaer = 1
enddo

do k = 1, kmx+1
  fabsh2o(k) = 0.d0
  fabso3(k) = 0.d0
  tauc(k) = 0.d0
  taua(k) = 0.d0
  fnebmax(k) = 0.d0
  if(iaer.ge.1) then
    fneba(k) = 1.d0
  endif
  do l = 1, 2
    fabso3c(k,l) = 0.d0
  enddo
  do l = 1, 8
    tau(k,l) = 0.d0
    pic(k,l) = 0.d0
    atln(k,l) = 0.d0
    absn(k,l) = 0.d0
    ref(k,l) = 0.d0
    reft(k,l) = 0.d0
    refts(k,l) = 0.d0
    refs(k,l) = 0.d0
    refb(k,l) = 0.d0
    refbs(k,l) = 0.d0
    trat(k,l) = 0.d0
    tras(k,l) = 0.d0
    trab(k,l) = 0.d0
    trats(k,l) = 0.d0
    tra(k,l) = 0.d0
    upw(k,l) = 0.d0
    dow(k,l) = 0.d0
  enddo
enddo

!  constantes caracteristiques des aerosols
!  hauteur couche d'aerosols
zaero  = 11000.d0
raero  = 0.1d0
! Leighton 1980
! (M.Tombette 2008   piaero = 0.92)
piaero = 0.84d0

gaero  = 0.66d0
nraer  = 1.55d0
niaer  = 0.01d0

! caero = 0. pour annuler les aerosols
caero = 1.d-9

k1p1 = k1+1

!  2 - calcul de muzero et de la constante solaire fo
!  ===================================================
!        calcul de muzero et de la constante solaire fo
!        muzero = cosinus de l*angle zenithal
!        fo = constante solaire en watt/m2

!
!  attention : 0. h < heuray < 24. h
!  ---------

qureel = float(squant)
call raysze(xlat,xlon,qureel,heuray,imer1,albe,muzero,fo)

! si muzero est negatif ou nul il fait nuit et les calculs
! de rayonnement solaire ne sont pas effectues

if(muzero.gt.0.d0) then

! correction pour les angles zenithal rasant

  rr1 = 0.1255d-2
  muzero = rr1/(sqrt(muzero**2 + rr1*(rr1 + 2.d0)) - muzero)
  m = 35.d0/sqrt(1224.d0*muzero*muzero + 1.d0)
  mbar = 1.9d0

!  3 - calcul des differents albedos intervenant dans le rayonnement
!      solaire pour l'ozone  et la diffusion rayleigh
!  ==================================================================
  rabar = 0.219d0/(1.d0 + 0.816d0*muzero)
  rabar2 = 0.144d0
  rbar = rabar + (1.d0 - rabar)*(1.d0 - rabar2)*albe/(1.d0 - rabar2*albe)
  rrbar = 0.28d0/(1.d0 + 6.43d0*muzero)
  rrbar2s = 0.0685d0

!  4 - ajout d'un niveau supplementaire pour le ryt solaire
!  ========================================================   d0
  zqq(kmray+1) = 16000.d0
  qqvtot = qqvinf + qqqv(kmray)
  qqv(kmray+1) = qqvtot - qqvinf/4.d0

! test sur la presence des couches nuageuses ou d'aerosols

  if((inua.eq.0).and.(iaer.eq.0)) then

!   5 -  calcul du rayonnement solaire par ciel clair
!       (ciel clair = sans particules (ni gouttes ni aero)
!============================================================================

! calcul du rechauffement solaire dans les couches par ciel clair

    rayst(k1) = 0.d0

! calcul pour la vapeur d'eau
! dans ce cas on utilise les epaisseurs optiques pour la vapeur d'eau
! qui ont ete calculees pour le rayonnement ir

    do i = k1p1, kmray
      ym1 = m*(qqvtot - qqv(i))
      if(i.eq.k1p1) ym1 = m*qqvtot
      y = m*(qqvtot - qqv(i+1))
      ystarm1 = m*qqvtot + 5.d0/3.d0*qqv(i)
      if(i.eq.k1p1) ystarm1 = m*qqvtot
      ystar = m*qqvtot + 5.d0/3.d0*qqv(i+1)
      fabsh2o(i) = muzero*fo*(raysve(ym1) - raysve(y) + albe*(raysve(ystar)   &
                 - raysve(ystarm1)))

! rechauffement par l ozone

      zqm1 = zqq(i)
      zq = zqq(i+1)
      if(i.eq.k1p1) zqm1 = zray(k1)
      xm1 = m*rayuoz(zqm1)
      x = m*rayuoz(zq)
      zbas = zray(k1)
      xstarm1 = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zqm1))
      xstar = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zq))
      fabso3(i) = muzero*fo*(raysoz(xm1) - raysoz(x) + rbar                  &
                 *(raysoz(xstar) - raysoz(xstarm1)))

! rechauffement total

      fabs = fabsh2o(i) + fabso3(i)
      cphum = cp0*(1.d0 + (cpvcpa - 1.d0)*qvray(i))
      rayst(i) = fabs/romray(i)/cphum/(zq-zqm1)

    enddo

! calcul du flux de rayonnement solaire au sol

    foo3 = muzero*fo*(0.647d0 - rrbar - raysoz(m*rayuoz(zbas)))*             &
         (1.d0 - albe)/(1.d0 - rrbar2s*albe)
    foh2o = muzero*fo*(0.353d0 - raysve(m*qqvtot))*(1.d0 - albe)

    fos = foh2o + foo3

  else

!     6 - calcul du rayonnement solaire par ciel nuageux
!     la prise en compte de la nebulosite fractionnaire impose de
!     realiser les calculs de diffusion multiple dans les cas de ciel
!     couvert (indice 1 des tableaux de transmission tra et de reflexion
!     ref) et de ciel clair (indice 2)

!============================================================================

!  6.1 recherche de la position du nuage, on determine le sommet du nuage
!      le plus haut et la base du plus bas
!      ........    .........................................................

    itop = 0

    do i = kmray, k1p1, -1
      if(qlray(i).gt.1.d-8.or.aeroso(i).gt.1.d-5) then
        if(itop.eq.0) then
          itop = i
          ibase = i
        else
          ibase = i
        endif
      endif
    enddo

! si itop = 0 il n'y a pas de nuage mais on peut tout de meme executer
! la methode d'addition des couches adjacentes

    if(itop.eq.0) then
      itop = k1
      ibase = k1
    endif
    itopp1 = itop +1
    itopp2 = itop +2
    ibasem1 = ibase -1

! 6.2 calcul des parametres optiques du nuage, albedo de simple
!     diffusion et epaisseur optique suivant Fouquart et Bonnel (1980)
!     .........    .....................................................

    fnebmax(kmray+1) = 0.d0
    tauctot = 0.d0
    tauatot = 0.d0
    do i = kmray, k1p1, -1
      if((i.ge.ibasem1).and.(i.le.itopp1)) then
! densite de l'eau liquide en g/m3 dans les couches considerees
        wh2ol = 1.d3*(romray(i)*qlray(i)) + aeroso(i)
! rayon moyen des gouttes en microns
        rm = 30.d0*wh2ol + 2.d0
! on fixe le rayon moyen max a 10 microns ds le domaine etudie
! et a 2 microns au-dessus
        if(i.le.nbmett) then
          rm = min(10.d0,rm)
        else
          rm = min(2.d0,rm)
        endif
        req = 3.d0*rm/2.d0
        deltaz = zqq(i+1) - zqq(i)
        if(i.eq.k1p1) deltaz = zqq(i+1) - zray(k1)
! le rayon req doit etre exprime en microns
        tauc(i) = 1.5d0*wh2ol*deltaz/req
        tauctot = tauctot + tauc(i)
      else
        tauc(i) = 0.d0
      endif

! calcul des parametres optiques pour les aerosols que l'on assimile
! a des couches de nuage de proprietes differentes

      fneba(i) = 0.d0
      if((iaer.eq.1).and.(zray(i).le.zaero)) then
        fneba(i) = 1.d0
        waero = 1.e3*romray(i)*caero*aeroso(i)
        deltaz = zqq(i+1) - zqq(i)
        if(i.eq.k1p1) deltaz=zqq(i+1) - zray(k1)
! le rayon raero doit etre exprime en microns
        reaero = 3.d0*raero/2.d0
        taua(i) = 1.5d0*waero*deltaz/reaero
        tauatot = tauatot + taua(i)
      endif

! calcul de la nebulosite maximum

      fnebmax(i) = max(fnebmax(i+1),fneray(i))
    enddo

    fnebmax(k1) = fnebmax(k1p1)

! albedo de simple diffusion de l'ensemble du nuage
!     par defaut (eau pure pioc=1)
    pioc = 0.9988d0
    tauc(kmray+1) = 0.d0

! 6.3 calcul de l'absorption par l'ozone en presence de nuages
!     .........    ...............................................

! calcul des differents albedos intervenant dans le rayonnement
! solaire pour l'ozone
! (Stephens, 74)

! facteur d'assymetrie
    gasym = 0.85d0

    s3 = sqrt(3.d0)
    rabarc = s3*(1.d0 - gasym)*tauctot/(2.d0 + s3*(1.d0 - gasym)*tauctot)
    rabar2c = rabarc
    rbarc = rabarc + (1.d0 - rabarc)*(1.d0 - rabar2c)*albe/(1.d0 - rabar2c*albe)
    rabara = s3*(1.d0 - gaero)*tauatot/(2.d0 + s3*(1.d0 - gaero)*tauatot)
    rabar2a = rabara
    rbara = rabara + (1.d0 - rabara)*(1.d0 - rabar2a)*albe/(1.d0 - rabar2a*albe)

! dans le cas d'un ciel totalement couvert,
! on ne calcule l'absorption que pour les couches au dessus de la
! surface reflective la plus haute, ici le sommet de la couche nuageuse

! l'absorption est calculee en ponderant la contribution des flux
! suivant la nebulosite fractionnaire maximum rencontree selon le
! chemin parcouru

    do i = itop, kmray

      zqm1 = zqq(i)
      if(i.eq.k1p1) zqm1 = zray(k1)
      zq = zqq(i+1)
      xm1 = m*rayuoz(zqm1)
      x = m*rayuoz(zq)
!  ciel couvert
      zbas = zray(itop)
      xstarm1 = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zqm1))
      xstar = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zq))
      fabso3c(i,1) = muzero*fo*(fnebmax(i-1)*                                &
                   (raysoz(xm1) - raysoz(x))                                 &
                   + fnebmax(k1p1)*rbarc*(raysoz(xstar) - raysoz(xstarm1)))
!  ciel clair
      zbas = zray(k1)
      xstarm1 = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zqm1))
      xstar = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zq))
      fabso3c(i,2) = muzero*fo*((1.d0 - fnebmax(i-1))*                        &
                     (raysoz(xm1) - raysoz(x))                                &
                   + (1.d0 - fnebmax(k1p1))*rbar                              &
                   * (raysoz(xstar) - raysoz(xstarm1)))

!   les aerosols sont pris en compte dans les calculs ciel clair
!   avec une nebulosite de 1
      if((iaer.eq.1).and.(zqm1.gt.zaero)) then
        zbas = zaero
        xstarm1 = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zqm1))
        xstar = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zq))
        fabso3c(i,2) = muzero*fo*((1.d0 - fnebmax(i-1))*                        &
                       (raysoz(xm1) - raysoz(x))                                &
                     + (1.d0 - fnebmax(k1p1))*rbara                             &
                     *(raysoz(xstar) - raysoz(xstarm1)))
      endif

      if((iaer.eq.1).and.(zqm1.le.zaero)) then
        x = m*rayuoz(zaero)
        xstar = m*rayuoz(zaero)
        fabso3c(i,2) = 0.d0
      endif

      fabso3(i) = fabso3c(i,1) + fabso3c(i,2)

    enddo

    do i = k1p1, itop-1

!  ciel couvert
      fabso3c(i,1) = 0.d0
!  ciel clair
      zqm1 = zqq(i)
      if(i.eq.k1p1) zqm1 = zray(k1)
      zq = zqq(i+1)
      xm1 = m*rayuoz(zqm1)
      x = m*rayuoz(zq)
      zbas = zray(k1)
      xstarm1 = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zqm1))
      xstar = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zq))
      fabso3c(i,2) = muzero*fo*((1.d0 - fnebmax(i-1))*                        &
                     (raysoz(xm1) - raysoz(x))                                &
                   + (1.d0-fnebmax(k1p1))*rbar                                &
                   * (raysoz(xstar)-raysoz(xstarm1)))

      if((iaer.eq.1).and.(zqm1.gt.zaero)) then
        zbas = zaero
        xstarm1 = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zqm1))
        xstar = m*rayuoz(zbas) + mbar*(rayuoz(zbas) - rayuoz(zq))
        fabso3c(i,2) = muzero*fo*((1.d0 - fnebmax(i-1))*                     &
                     (raysoz(xm1) - raysoz(x))                               &
                     + (1.d0 - fnebmax(k1p1))*rbara                          &
                     * (raysoz(xstar) - raysoz(xstarm1)))
      endif

      if((iaer.eq.1).and.(zqm1.le.zaero)) then
        x = m*rayuoz(zaero)
        xstar = m*rayuoz(zaero)
        fabso3c(i,2) = 0.d0
      endif

      fabso3(i) = fabso3c(i,1) + fabso3c(i,2)

    enddo

! 6.4 calcul de l'absorption par la vapeur d'eau et l'eau nuageuse
!     .........    ...................................................

! dans ce cas il faut resoudre le probleme de la diffusion multiple
! ceci est realise par une methode d'addition des couches adjacentes
! ou "adding method" suivant Lacis et Hansen, 1974

! calcul des fonctions de reflexion et de transmission de chacune
! des couches considerees

    do n = 1, 8
      do l = k1p1, kmray
        dqqv = kn(n)*(qqv(l+1) - qqv(l))/10.d0
        if(l.eq.k1p1) dqqv = kn(n)*qqv(l+1)/10.d0
! ciel couvert
        tau(l,n) = tauc(l) + dqqv + taua(l)
        if(tauc(l).ne.0.) then
          pioc = 0.9988d0
          pic(l,n) = (pioc*tauc(l) + piaero*taua(l))/tau(l,n)
          gas = (pioc*tauc(l)*gasym + piaero*taua(l)*gaero)               &
              /(pic(l,n)*tau(l,n))
! correction Joseph, 1976 pour angle rasant
          fas = gas*gas
          tau(l,n) = (1.d0 - pic(l,n)*fas)*tau(l,n)
          pic(l,n) = pic(l,n)*(1.d0 - fas)/(1.d0 - pic(l,n)*fas)
          gas = (gas - fas)/(1.d0 - fas)

          gama1 = (s3/2.d0)*(2.d0 - pic(l,n)*(1.d0 + gas))
          gama2 = (s3*pic(l,n)/2.d0)*(1.d0 - gas)
          kt = sqrt(gama1*gama1 - gama2*gama2)
          tln = kt*tau(l,n)
          extlnp = exp(tln)
          extlnm = exp(-tln)
          drt = (kt + gama1)*extlnp + (kt-gama1)*extlnm
          ref(l,n) = fneray(l)*gama2*(extlnp - extlnm)/drt
          tra(l,n) = fneray(l)*2.d0*kt/drt                                &
                    + (1.d0 - fneray(l))*exp(-5.d0*dqqv/3.d0)
          refs(l,n) = ref(l,n)
          tras(l,n) = tra(l,n)
        else
! dans les couches ciel clair
          ref(l,n) = 0.d0
          tra(l,n) = exp(-5.d0*tau(l,n)/3.d0)
          refs(l,n) = ref(l,n)
          tras(l,n) = tra(l,n)
          if(l.ge.itopp1) tra(l,n) = exp(-m*tau(l,n))
! dans les couches d'aerosols
          if((iaer.eq.1).and.(zray(l).le.zaero)) then
            tau(l,n) = taua(l) + dqqv
            pioc = piaero
            pic(l,n) = pioc*taua(l)/tau(l,n)
            gas = gaero
! correction Joseph, 1976 pour angle rasant
            fas = gas*gas
            tau(l,n) = (1.d0 - pic(l,n)*fas)*tau(l,n)
            pic(l,n) = pic(l,n)*(1.d0 - fas)/(1.d0 - pic(l,n)*fas)
            gas = (gas - fas)/(1.d0 - fas)
            gama1 = (s3/2.d0)*(2.d0 - pic(l,n)*(1.d0 + gas))
            gama2 = (s3*pic(l,n)/2.d0)*(1.d0 - gas)
            kt = sqrt(gama1*gama1 - gama2*gama2)
            tln = kt*tau(l,n)
            extlnp = exp(tln)
            extlnm = exp(-tln)
            drt = (kt+gama1)*extlnp + (kt - gama1)*extlnm
            ref(l,n) = fneba(l)*gama2*(extlnp - extlnm)/drt
            tra(l,n) = fneba(l)*2.d0*kt/drt                                &
                     + (1.d0 - fneba(l))*exp(-5.d0*dqqv/3.d0)
          endif
        endif

      enddo

! conditions aux limites

      tau(kmray+1,n) = kn(n)*qqvinf/40.d0
      tra(kmray+1,n) = exp(-m*tau(kmray+1,n))
      ref(kmray+1,n) = 0.d0
      tras(kmray+1,n) = tra(kmray+1,n)
      refs(kmray+1,n) = ref(kmray+1,n)
      tra(k1,n) = 0.d0
      ref(k1,n) = albe
      tras(k1,n) = 0.d0
      refs(k1,n) = 0.d0

! addition des couches en descendant

      trat(kmray+1,n) = tra(kmray+1,n)
      reft(kmray+1,n) = ref(kmray+1,n)
      refts(kmray+1,n) = refs(kmray+1,n)
      trats(kmray+1,n) = tras(kmray+1,n)
      fneray(k1) = 0.d0

      do l = kmray, k1, -1
        drtt1 = 1.d0 - refts(l+1,n)*ref(l,n)
        reft(l,n) = reft(l+1,n)                                          &
                  + trat(l+1,n)*ref(l,n)*trats(l+1,n)/drtt1
        trat(l,n) = trat(l+1,n)*tra(l,n)/drtt1
        if(l.gt.k1) then
          refts(l,n) = refs(l,n)                                         &
                     + tras(l,n)*refts(l+1,n)*tra(l,n)/drtt1
          trats(l,n) = trats(l+1,n)*tras(l,n)/drtt1
        endif
      enddo

! addition des couches en montant
      refb(k1,n) = ref(k1,n)
      refbs(k1,n) = refs(k1,n)

      do l = k1p1, kmray
        dtrb1 = 1.d0 - refb(l-1,n)*refs(l,n)
        refb(l,n) = ref(l,n) + tra(l,n)*refb(l-1,n)*tras(l,n)/dtrb1
      enddo

! calcul des flux montant et descendant

      do l = kmray+1, k1p1, -1
        dud1 = 1.d0 - refts(l,n)*refb(l-1,n)
        upw(l,n) = trat(l,n)*refb(l-1,n)/dud1
        dow(l,n) = trat(l,n)/dud1

! l'absorption est calculee en ponderant la contribution
! des flux montants et descendants par la nebulosite maximale
! rencontree entre l'infini et z

        atln(l,n) = pkn(n)*((1.d0 - reft(k1,n)) + upw(l,n)              &
                  - dow(l,n))
      enddo

! calcul de l'absorption dans les couches individuelles

      do l = kmray, k1p1, -1
        absn(l,n) = atln(l,n) - atln(l+1,n)
      enddo

    enddo

! sommation sur les frequences et calcul de l'absorption integree sur le
! spectre

    do l = kmray, k1p1, -1
      fabsh2o(l) = 0.d0
      do n = 1, 8
        fabsh2o(l) = fabsh2o(l)+absn(l,n)
      enddo
      fabsh2o(l) = fabsh2o(l)*fo*muzero
    enddo


! 6.5 calcul du rechauffement dans les couches
!     .........    ...............................

    rayst(k1) = 0.d0
    do i = k1p1, kmray
      deltaz = zqq(i+1) - zqq(i)
      if(i.eq.k1p1) deltaz = zqq(i+1) - zray(k1)
      cphum = cp0*(1.d0 + (cpvcpa - 1.d0)*qvray(i))
      rayst(i) = (fabsh2o(i) + fabso3(i))/deltaz/romray(i)/cphum
    enddo

! 6.6 calcul du flux solaire descendant au sol avec ponderation
!     suivant la nebulosite fractionnaire
!     .........    .................................................

    dowtot1 = 0.d0
    do n = 2, 8
      dowtot1 = dowtot1 + pkn(n)*dow(k1p1,n)
    enddo
    foh2o = fo*muzero*(1.d0 - albe)*dowtot1
    foo3c = fo*muzero*(0.647d0 - rrbar - raysoz(m*rayuoz(zray(itop))))      &
           *(1.d0 - rabarc)*(1.d0 - albe)/(1.d0 - rabarc*albe)
    foo3 = muzero*fo*(0.647d0 - rrbar - raysoz(m*rayuoz(zbas)))             &
           *(1.d0 - albe)/(1.d0 - rrbar2s*albe)

    if(iaer.ne.0.d0) then
      foo3 = fo*muzero*(0.647d0 - rrbar - raysoz(m*rayuoz(zaero)))          &
          *(1.d0 - rabara)*(1.d0 - albe)/(1.d0 - rabara*albe)
    endif

    foo3 = foo3c*fnebmax(k1p1) + foo3*(1.d0 - fnebmax(k1p1))
    fos = foh2o + foo3

  endif

! si muzero est inferieur a zero il fait nuit
else

  muzero = 0.d0
  do k = k1, kmray
    rayst(k) = 0.d0
  enddo

endif

deallocate(fabsh2o,fabso3,tauc)
deallocate(tau,pic,ref)
deallocate(reft,refts)
deallocate(refs,tras,trats)
deallocate(refb,trab,upw)
deallocate(refbs,fabso3c,tra)
deallocate(dow,atln,absn)
deallocate(fnebmax,fneba)

return

!========================
!7. additional functions:
!========================
contains


! 7.1 Computes ozone concentration for the altitude zh
!-------------------------------------------------------
function rayuoz(zh)

  implicit none
  double precision, intent(in) :: zh  ! absolute altitude
  double precision ::  rayuoz
  double precision ::  a, b, c

  a = 0.4d0
  b = 20000.d0
  c = 5000.d0

  rayuoz = a*(1.d0 + exp(-b/c))/(1.d0 + exp((zh-b)/c))

end function rayuoz

! 7.2 Aborption function of the solar radiation by water vapor
!-------------------------------------------------------------------------

function raysve(y)

  implicit none
  double precision, intent(in) :: y        ! specific humidity
  double precision :: raysve

  raysve = 0.29d0*y/((1.d0 + 14.15d0*y)**0.635d0 + 0.5925d0*y)

end function raysve

! 7.3 Aborption function of the solar radiation by ozone
!--------------------------------------------------------

function raysoz(x)
  implicit none
  double precision, intent(in) :: x
  double precision :: raysoz
  double precision :: ao3vis, ao3uv

  ao3vis = 0.02118d0*x/(1.d0 + (0.042d0 + 0.000323d0*x)*x)
  ao3uv = 1.082d0*x/(1.d0 + 138.6d0*x)**0.805d0                     &
        + 0.0658d0*x/(1.d0 + (103.6d0*x)**3)
  raysoz = ao3vis + ao3uv

end function raysoz

end subroutine rayso

