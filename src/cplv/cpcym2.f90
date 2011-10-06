!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
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

subroutine cpcym2 &
!================

 ( ncelet , ncel   , nrtbmc ,                                     &
   indpdf ,                                                       &
   rtp    ,                                                       &
   f1m    , f2m    , f3m   , f4m   , f5m  , f6m , f7m ,           &
   f3max  ,                                                       &
   f1cl   , f2cl   , f3cl  , f4cl  , f5cl , f6cl , f7cl ,         &
   f4m1   , f4m2   , d4cl  , d4f4  , hrec ,                       &
   rtbmc  , w1     ,                                              &
   fuel1  , fuel2  , fuel3 , oxyd  , prod1 , prod2  ,             &
   xiner  )

!===============================================================================
! FONCTION :
! --------

! CALCUL DES CONCENTRATIONS GAZEUSES MOYENNES
!  VALEURS CELLULES
!  ----------------

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant)                  !
! fvap             ! tr ! <-- ! moyenne du traceur 1 mvl [fovm+co]             !
! fhet             ! tr !  ->           ! moyenne du traceur 3 c héterogène    !
! f4p2m            ! tr ! <-- ! variance du traceur 4 (air)                    !
! indpdf           ! te ! <-- ! passage par les pdf                            !
! f1m              ! tr ! <-- ! moyenne du traceur 1 mvl [chx1m+co]            !
! f2m              ! tr ! <-- ! moyenne du traceur 2 mvl [chx2m+co]            !
! f3m              ! tr ! <-- ! moyenne du traceur 3 (co c.het)                !
! f4m              ! tr ! <-- ! moyenne du traceur 4 (air)                     !
! f5m              ! tr ! <-- ! moyenne du traceur 5 (h2o)                     !
! fuel1            ! tr !  <- ! fraction massique chx1m                        !
! fuel2            ! tr !  <- ! fraction massique chx2m                        !
! fuel3            ! tr !  <- ! fraction massique co                           !
! oxyd             ! tr !  <- ! fraction massique o2                           !
! prod1            ! tr !  <- ! fraction massique co2                          !
! prod2            ! tr !  <- ! fraction massique h2o                          !
! xiner            ! tr !  <- ! fraction massique n2                           !
! f4m1             ! tr ! <-- ! borne minimum                                  !
! f4m2             ! tr ! <-- ! borne max                                      !
! d4cl             ! tr ! <-- ! amplitude du pic de dirac en f4cl              !
! d4f4             ! tr ! <-- ! amplitude du pic de dirac en 1                 !
!                  !    !     ! (air pur)                                      !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHAMNUMERIQUE), T (TABLEAU)
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
use cstnum
use entsor
use parall
use pointe
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel , nrtbmc
integer          indpdf(ncelet)

double precision rtp(ncelet,*)
double precision f1m(ncelet)  , f2m(ncelet)
double precision f3m(ncelet)  , f4m(ncelet) , f5m(ncelet)
double precision f6m(ncelet)  , f7m(ncelet) , f3max(ncelet)

double precision f1cl(ncelet) , f2cl(ncelet)
double precision f3cl(ncelet) , f4cl(ncelet), f5cl(ncelet)
double precision f6cl(ncelet) , f7cl(ncelet)

double precision f4m1(ncelet) , f4m2(ncelet) , d4cl(ncelet)
double precision d4f4(ncelet) , hrec(ncelet)

double precision fuel1(ncelet), fuel2(ncelet) , fuel3(ncelet)
double precision oxyd(ncelet) , prod1(ncelet), prod2(ncelet)
double precision xiner(ncelet)

double precision rtbmc(ncelet,nrtbmc)
double precision w1(ncelet)

! Local variables

integer          iel , icha , ii

integer          n1 , n2 , n3 , n4 , n5 , n6 , n7 , n8 , n9
integer          n10
integer          itbr , icla
integer          nclo2,nclco,nclo21,nclco1
integer          nbrint
parameter       (nbrint = 10)
integer          inttmp(nbrint)

double precision anu1 , anu2 , anu3
double precision somm

double precision zchx10 , zchx20 , zco0   , zco20
double precision zn20   , zh2o0  , zo20
double precision ychx10 , ychx20
double precision zchx11 , zchx21 , zco1   , zco21
double precision zn21   , zh2o1  , zo21
double precision zchx12 , zchx22 , zco2   , zco22
double precision zn22   , zh2o2  , zo22
double precision zchx13 , zchx23 , zco3   , zco23
double precision zn23   , zh2o3  , zo23
double precision reac1  , reac2  , reac3
double precision x2     , xch    , xck    , xash , xwat

double precision zcof10 , zcof20 , zh2o10 , zh2o20
double precision cx1m   , cx2m
double precision wmchx1 , wmchx2 , wmco   , wmco2 , wmh2o
double precision wmo2   , wmn2   , wmc
double precision achx1f1, acof1  , achx2f2, acof2
double precision yo2ox
double precision ah2of1 , ah2of2

!     Variables pour le support de la Pdf et sa description

double precision f42m
double precision zzcl(7),zzs1(7),zzs2(7),zzs3(7)
double precision zz(7)  ,zzoxy(7)

double precision f4s1 , f4s2 , f4s3
double precision bb1  , bb2  , den1 , den2
double precision fp1 , fp2 , fp3 , fp4 , fp5 , zco2t , zcot

double precision f4loc,f6loc,f7loc,s4loc,s6loc,s7loc
double precision foxy ,foxycl

double precision sommin,sommax,ao2min,ao2max

double precision cmaxs1,cmins1,cmaxa1,cmina1
double precision cmaxs2,cmins2,cmaxa2,cmina2

!===============================================================================

!===============================================================================
! 0. INITIALISATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings

zchx13 = 0.d0
zchx23 = 0.d0
zco23 = 0.d0
zco2t = 0.d0
zco3 = 0.d0
zh2o3 = 0.d0
zn23 = 0.d0
zo23 = 0.d0

! --- Initialisation tableau de travail

do iel = 1, ncel
  w1(iel) = zero
  do itbr = 1, nrtbmc
    rtbmc(iel,itbr) = zero
  enddo
enddo

!===============================================================================
! 1. CALCULS PRELIMINAIRES
!===============================================================================

! --- Masses molaires

wmco   = wmole(ico  )
wmo2   = wmole(io2  )
wmco2  = wmole(ico2 )
wmh2o  = wmole(ih2o )
wmn2   = wmole(in2)
wmc    = wmolat(iatc)

! --- Calcul de F1M de chaque charbon dans RTBMC
!            de F2M de chaque charbon dans RTBMC

! ------ W1 = - Somme des X2(icla)

do icla = 1, nclacp
  do iel = 1, ncel
    xck  = rtp(iel,isca(ixck(icla)))
    xch  = rtp(iel,isca(ixch(icla)))
    xash = rtp(iel,isca(inp (icla)))*xmash(icla)
    if ( ippmod(icp3pl) .eq. 1 ) then
      xwat = rtp(iel,isca(ixwt(icla)))
    else
      xwat = 0.d0
    endif
    x2   = xch + xck + xash + xwat
    w1(iel) =  w1(iel) - x2
  enddo
enddo

! ------ RTBMC(IF1MC(ICHA)) = F1M(ICHA)
!        RTBMC(IF2MC(ICHA)) = F2M(ICHA)

do icha = 1, ncharb
  do iel = 1, ncel
    rtbmc(iel,if1mc(icha)) = rtp(iel,isca(if1m(icha)))            &
                           / (1.d0+w1(iel))
    rtbmc(iel,if2mc(icha)) = rtp(iel,isca(if2m(icha)))            &
                           / (1.d0+w1(iel))
  enddo
enddo


! --- Calcul de CX1M, CX2M
!     et des coefficient relatifs a l'expression de
!        ZCHx1m0, ZCHx2m0, ZCO0, ZO20 et ZN20
!        introduction de ZH2O0
!        fonctions de F1, F2, F3 et F4

do iel = 1, ncel

! ---- YCHX1,20 en kg/kg de melange et
!      ZCHX1,20 en nb de moles/ kg de melange

  ychx10 = zero
  ychx20 = zero
  zchx10 = zero
  zchx20 = zero
  zcof10 = zero
  zcof20 = zero
  zh2o10 = zero
  zh2o20 = zero
  do icha = 1, ncharb
    den1   = 1.d0                                                 &
           / ( a1(icha)*wmole(ichx1c(icha))                       &
              +b1(icha)*wmole(ico)                                &
              +c1(icha)*wmole(ih2o) )
    ychx10 = ychx10 + den1 *                                      &
           ( rtbmc(iel,if1mc(icha))*a1(icha)*wmole(ichx1c(icha)) )
    zchx10 = zchx10 + den1 *                                      &
           ( rtbmc(iel,if1mc(icha))*a1(icha) )
    zcof10 = zcof10 + den1 *                                      &
           ( rtbmc(iel,if1mc(icha))*b1(icha) )
    zh2o10 = zh2o10 + den1 *                                      &
           ( rtbmc(iel,if1mc(icha))*c1(icha) )
    den2   = 1.d0                                                 &
           / ( a2(icha)*wmole(ichx2c(icha))                       &
              +b2(icha)*wmole(ico)                                &
              +c2(icha)*wmole(ih2o) )
    ychx20 = ychx20 + den2 *                                      &
           ( rtbmc(iel,if2mc(icha))*a2(icha)*wmole(ichx2c(icha)) )
    zchx20 = zchx20 + den2 *                                      &
           ( rtbmc(iel,if2mc(icha))*a2(icha) )
    zcof20 = zcof20 + den2 *                                      &
           ( rtbmc(iel,if2mc(icha))*b2(icha) )
    zh2o20 = zh2o20 + den2 *                                      &
           ( rtbmc(iel,if2mc(icha))*c2(icha) )
  enddo
  rtbmc(iel,ix1mc) = 4.d0
  if ( zchx10.gt.epzero )                                         &
    rtbmc(iel,ix1mc) = ( ychx10/zchx10-wmolat(iatc) )             &
                     / wmolat(iath)
  rtbmc(iel,ix2mc) = 2.d0
  if ( zchx20.gt.epzero )                                         &
    rtbmc(iel,ix2mc) = ( ychx20/zchx20-wmolat(iatc) )             &
                     / wmolat(iath)

  rtbmc(iel,ichx1f1) = 0.d0
  rtbmc(iel,icof1  ) = 0.d0
  rtbmc(iel,ih2of1 ) = 0.d0
  if ( f1m(iel).gt.epzero ) then
    rtbmc(iel,ichx1f1) = zchx10 / f1m(iel)
    rtbmc(iel,icof1  ) = zcof10 / f1m(iel)
    rtbmc(iel,ih2of1 ) = zh2o10 / f1m(iel)
  endif
  rtbmc(iel,ichx2f2) = 0.d0
  rtbmc(iel,icof2  ) = 0.d0
  rtbmc(iel,ih2of2 ) = 0.d0
  if ( f2m(iel).gt.epzero ) then
    rtbmc(iel,ichx2f2) = zchx20 / f2m(iel)
    rtbmc(iel,icof2  ) = zcof20 / f2m(iel)
    rtbmc(iel,ih2of2 ) = zh2o20 / f2m(iel)
  endif

enddo

! ---- Initialisation des fractions massiques avec de l'air

do iel = 1, ncel
  fuel1(iel)  = zero
  fuel2(iel)  = zero
  fuel3(iel)  = zero
  oxyd(iel)   = wmo2*ao2f4
  prod1(iel)  = zero
  prod2(iel)  = zero
  xiner(iel)  = wmn2*an2f4
enddo

nclo2  = 0
nclco  = 0
nclo21 = 0
nclco1 = 0

cmina1 =  1.d+20
cmaxa1 = -1.d+20
cmina2 =  1.d+20
cmaxa2 = -1.d+20
cmins1 =  1.d+20
cmaxs1 = -1.d+20
cmins2 =  1.d+20
cmaxs2 = -1.d+20

!===============================================================================
! 3. CALCUL DE LA COMPOSITION DU MELANGE SANS LES PDF
!    SI LES FLUCTUATIONS DE S SONT TROP FAIBLES
!===============================================================================

do iel = 1, ncel

  if ( indpdf(iel).eq.0 ) then

! --> Calculs preliminaires

    if ( ieqco2 .eq. 1 ) then
      zco2t = (rtp(iel,isca(iyco2))/ (1.d0+w1(iel)))/wmco2
    else if ( ieqco2 .eq. 2 ) then
      zcot  = (rtp(iel,isca(iyco2))/ (1.d0+w1(iel)))/wmco
    endif

    cx1m    = rtbmc(iel,ix1mc)
    cx2m    = rtbmc(iel,ix2mc)
    wmchx1  = wmolat(iatc) + cx1m*wmolat(iath)
    wmchx2  = wmolat(iatc) + cx2m*wmolat(iath)
    achx1f1 = rtbmc(iel,ichx1f1)
    achx2f2 = rtbmc(iel,ichx2f2)
    acof1   = rtbmc(iel,icof1)
    acof2   = rtbmc(iel,icof2)
    ah2of1  = rtbmc(iel,ih2of1)
    ah2of2  = rtbmc(iel,ih2of2)

! --> Composition de la phase gazeuse avant combustion (indice 0)

!          ZCHX10 = ACHX1F1*F1
!          ZCHX20 = ACHX2F2*F2
!          ZCO0   = ACOF1*F1 + ACOF2*F2 + ACOF3*F3
!          ZCO20  = ZERO
!          ZH2O0  = ZERO
!          ZO20   = AO2F4*F4 + AO2F3*F3
!          ZN20   = AN2F4*F4

    zchx10  = achx1f1*f1m(iel)
    zchx20  = achx2f2*f2m(iel)
    zco0    = (acof1*f1m(iel)+acof2*f2m(iel) +                    &
               acof3*f3m(iel))
    zo20    = ao2f4*f4m(iel)+ao2f6*f6m(iel)+ao2f7*f7m(iel)        &
             +ao2f3*f3m(iel)
    zn20    = an2f4*f4m(iel)+an2f6*f6m(iel)+an2f7*f7m(iel)
    zco20   = aco2f4*f4m(iel)+aco2f6*f6m(iel)+aco2f7*f7m(iel)
    zh2o0   = ah2of4*f4m(iel)+ah2of6*f6m(iel)+ah2of7*f7m(iel)     &
             +ah2of5*f5m(iel)+ah2of1*f1m(iel)+ah2of2*f2m(iel)

! --> Calcul de la composition du melange

    anu1   = 0.25d0*(cx1m-cx2m)
    reac1  = min(zchx10, (zo20/anu1))
    zchx11 = zchx10 -           reac1
    zo21   = zo20   -      anu1*reac1
    zchx21 = zchx20 +           reac1
    zco1   = zco0
    zco21  = zco20
    zh2o1  = zh2o0  + 2.d0*anu1*reac1
    zn21   = zn20

    anu2   = 0.25d0*(2.d0 + cx2m)
    reac2  = min(zchx21, (zo21/anu2))
    zchx12 = zchx11
    zchx22 = zchx21 -           reac2
    zo22   = zo21   -      anu2*reac2
    zco2   = zco1   +           reac2
    zco22  = zco21
    zh2o2  = zh2o1  + .5d0*cx2m*reac2
    zn22   = zn21

    anu3   = 0.5d0
    if ( ieqco2 .eq. 0 ) then
      reac3  = min(zco2, (zo22/anu3))
      zchx13 = zchx12
      zchx23 = zchx22
      zo23   = zo22  - anu3*reac3
      zco3   = zco2  -      reac3
      zco23  = zco22 +      reac3
      zh2o3  = zh2o2
      zn23   = zn22

    else if ( ieqco2 .eq. 1 ) then

!      Transport de CO2

      zchx13 = zchx12
      zchx23 = zchx22

      if ( (zco2-(zco2t-zco22)) .lt. -epzero ) then
        nclco = nclco +1
        cmins1 = min((zco2-zco2t),cmins1)
        cmaxs1 = max((zco2-zco2t),cmins1)
      endif
      if ( (zo22-anu3*(zco2t-zco22)) .lt. -epzero ) then
        nclo2 = nclo2 +1
        cmins2 = min((zo22-anu3*zco2t),cmins2)
        cmaxs2 = max((zo22-anu3*zco2t),cmins2)
      endif
      reac3 = min(max(zco2t-zco22,0.d0),zco2,zo22/anu3)

!            RTP(IEL,ISCA(IYCO2))= ZCO2T*(1.D0+W1(IEL))*WMCO2

      zo23   = zo22  - anu3*reac3
      zco3   = zco2  -      reac3
      zco23  = zco22 +      reac3

      zh2o3  = zh2o2
      zn23   = zn22

    else

      write(nfecra,3000) ieqco2
      call csexit(1)

    endif

    fuel1(iel) = zchx13 * wmchx1
    fuel2(iel) = zchx23 * wmchx2
    fuel3(iel) = zco3   * wmco
    prod1(iel) = zco23  * wmco2
    prod2(iel) = zh2o3  * wmh2o
    oxyd(iel)  = zo23   * wmo2
    xiner(iel) = zn23   * wmn2

  endif

enddo

!===============================================================================
! 4. CALCUL DE LA COMPOSITION DU MELANGE AVEC LES PDF
!===============================================================================

do iel = 1, ncel

  if ( indpdf(iel).eq.1 ) then

! --> Calculs preliminaires

    if ( ieqco2 .eq. 1 ) then
      zco2t = (rtp(iel,isca(iyco2))/(1.d0+w1(iel)))/wmco2
    else if ( ieqco2 .eq. 2 ) then
      zcot  = (rtp(iel,isca(iyco2))/(1.d0+w1(iel)))/wmco
    endif

    cx1m    = rtbmc(iel,ix1mc)
    cx2m    = rtbmc(iel,ix2mc)
    wmchx1  = wmolat(iatc) + cx1m*wmolat(iath)
    wmchx2  = wmolat(iatc) + cx2m*wmolat(iath)
    achx1f1 = rtbmc(iel,ichx1f1)
    achx2f2 = rtbmc(iel,ichx2f2)
    acof1   = rtbmc(iel,icof1)
    acof2   = rtbmc(iel,icof2)
    ah2of1  = rtbmc(iel,ih2of1)
    ah2of2  = rtbmc(iel,ih2of2)

    f4loc = f4m(iel)/(f4m(iel)+f6m(iel)+f7m(iel))
    f6loc = f6m(iel)/(f4m(iel)+f6m(iel)+f7m(iel))
    f7loc = f7m(iel)/(f4m(iel)+f6m(iel)+f7m(iel))
    s4loc = oxyo2 (1)*wmo2 +oxyn2 (1)*wmn2                        &
           +oxyh2o(1)*wmh2o+oxyco2(1)*wmco2
    s6loc = oxyo2 (2)*wmo2 +oxyn2 (2)*wmn2                        &
           +oxyh2o(2)*wmh2o+oxyco2(2)*wmco2
    s7loc = oxyo2 (3)*wmo2 +oxyn2 (3)*wmn2                        &
           +oxyh2o(3)*wmh2o+oxyco2(3)*wmco2

    yo2ox  =  wmo2                                                &
            /(f4loc*s4loc+f6loc*s6loc+f7loc*s7loc)
    an2f3  = (1.d0-yo2ox)/wmole(in2) * (1.d0-f3max(iel))

    anu1   = 0.25d0*(cx1m-cx2m)
    anu2   = 0.25d0*(2.d0 + cx2m)
    anu3   = 0.5d0

    foxy   = f4m(iel)+f5m(iel)+f6m(iel)+f7m(iel)
    foxycl = f4cl(iel)+f5cl(iel)+f6cl(iel)+f7cl(iel)

!         Concentrations aux extrémités

    do ii = 1,7
      zzcl(ii)  = zero
      zzoxy(ii) = zero
    enddo

    zzcl(ichx1)= achx1f1*f1cl(iel)
    zzcl(ichx2)= achx2f2*f2cl(iel)
    zzcl(ico)  = acof1 *f1cl(iel)+acof2*f2cl(iel)                 &
                                 +acof3*f3cl(iel)
    zzcl(io2)  = ao2f4*f4cl(iel)+ao2f6*f6cl(iel)+ao2f7*f7cl(iel)  &
                +ao2f3*f3cl(iel)
    zzcl(in2)  = an2f4*f4cl(iel)+an2f6*f6cl(iel)+an2f7*f7cl(iel)
    zzcl(ih2o) = ah2of4*f4cl(iel)+ah2of6*f6cl(iel)                &
                +ah2of7*f7cl(iel)                                 &
                +ah2of1*f1cl(iel)+ah2of2*f2cl(iel)                &
                +ah2of5*f5cl(iel)
    zzcl(ico2) = aco2f4*f4cl(iel)+aco2f6*f6cl(iel)                &
                +aco2f7*f7cl(iel)

    zzoxy(io2)  = (ao2f4*f4m(iel)+ao2f6*f6m(iel)+ao2f7*f7m(iel))  &
                 /foxy
    zzoxy(in2)  = (an2f4*f4m(iel)+an2f6*f6m(iel)+an2f7*f7m(iel))  &
                 /foxy
    zzoxy(ih2o) = ( ah2of4*f4m(iel)+ah2of6*f6m(iel)               &
                   +ah2of7*f7m(iel)                               &
                   +ah2of5*f5m(iel) )                             &
                 /foxy
    zzoxy(ico2) = (aco2f4*f4m(iel)+aco2f6*f6m(iel)                &
                                  +aco2f7*f7m(iel))               &
                 /foxy

!         Calcul de la stoechiométrie de la première réaction

    if ( zzcl(ichx1)  .gt. 0.d0 .or.                              &
         zzoxy(io2)   .gt. 0.d0         ) then
      f4s1 = ( anu1 * zzcl(ichx1) + foxycl * zzoxy(io2) )         &
           / ( anu1 * zzcl(ichx1) +          zzoxy(io2) )
    else
      f4s1 = 1.d0
    endif

!         Calcul des concentrations au premier point stoechio
!         avant réaction

    do ii = 1,7
      zzs1(ii) = zzcl(ii)                                         &
                + (f4s1-foxycl)/(1.d0-foxycl)                     &
           *(zzoxy(ii)-zzcl(ii))
    enddo

    zzs1(ichx2) = zzs1(ichx2) + zzs1(ichx1)
    zzs1(io2)   = zero
    zzs1(ih2o)  = zzs1(ih2o) + 2.d0*anu1*zzs1(ichx1)
    zzs1(ichx1) = zero

!         Calcul de la stoechiometrie de la deuxième réaction

    if ( zzs1(ichx2) .gt. 0.d0 .or.                               &
         zzoxy(io2)  .gt. 0.d0         ) then
      f4s2 = ( anu2 * zzs1(ichx2) + f4s1 * zzoxy(io2) )           &
            /( anu2 * zzs1(ichx2) +        zzoxy(io2) )
    else
      f4s2 = 1.d0
    endif

!         Calcul des concentrations au deuxième point stoechio
!         avant réaction

    if ( abs(1.d0-f4s1).gt. 0.d0 ) then
      do ii = 1,7
         zzs2(ii) = zzs1(ii)                                      &
           + (f4s2-f4s1) / (1.d0-f4s1) * (zzoxy(ii) - zzs1(ii))
      enddo
    else
      do ii = 1,7
         zzs2(ii) = zzs1(ii)
      enddo
    endif

    zzs2(ico)   = zzs2(ico)  + zzs2(ichx2)
    zzs2(ih2o)  = zzs2(ih2o) + (2.d0*anu2-1.d0)*zzs2(ichx2)
    zzs2(io2)   = zero
    zzs2(ichx2) = zero

!         Calcul de la stoechiometrie de la troisième réaction
!         CO + O2 => CO2
!         Les concentrations sont désormais linéaires entre F4S2 et 1
!         ZZ(F4,II)  = ZZS2(II) + (F4-F4S2)/(1-F4S2)*(ZZF4(II)-ZZS2(II))
!         On cherche le point tel que
!         ZZ(f4s3,ICO) = ZZ(F4S3,IO2)/ANU3
!         La formule est semblable à la précedente en substituant
!         ANU2 par ANU3

    if ( zzs2(ico) .gt. 0.d0 .or.                                 &
         zzs2(io2) .gt. 0.d0         ) then
      f4s3 = ( anu3 * zzs2(ico ) + f4s2 * zzoxy(io2) )            &
           / ( anu3 * zzs2(ico ) +        zzoxy(io2) )
    else
      f4s3 = 1.d0
    endif

!         Calcul des concentrations au troisième point stoechio
!         avant réaction

    if ( abs(1.d0-f4s2).gt. 0.d0 ) then
      do ii = 1,7
         zzs3(ii) = zzs2(ii)                                      &
           + (f4s3-f4s2) / (1.d0-f4s2) * (zzoxy(ii) - zzs2(ii))
      enddo
    else
      do ii = 1,7
         zzs3(ii) = zzs2(ii)
      enddo
    endif

!         Le nombre de moles de réaction est le nombre de moles de CO

    if ( ieqco2 .eq. 0 ) then
      zzs3(ico2) = zzs3(ico2) + zzs3(ico)
      zzs3(io2)  = zero
      zzs3(ico)  = zero
    endif

!        Désormais on connait les concentrations en 5 points
!        cl,s1,s2,s3,f4
!        les concentrations intermédiaires sont linéaires par morceaux
!        et les paramêtres de la pdf D4cl,D4f4, F4m1,F4m2,HREC

!        On initialise par les concentrations aux extrémités
   do ii = 1,7
     zz(ii) = d4cl(iel)*zzcl(ii)+d4f4(iel)*zzoxy(ii)
   enddo

!        Intégration sur le premier intervalle de richesse (entre cl et s1)

   bb1 = max(f4m1(iel),foxycl)
   bb2 = min(f4m2(iel),f4s1)
   if( bb2.gt.bb1 ) then
     do ii = 1,7
       zz(ii) =  zz(ii) + hrec(iel)*(bb2-bb1)/(f4s1-foxycl)       &
               *( (zzcl(ii)*f4s1-zzs1(ii)*foxycl)                 &
              +   (zzs1(ii)-zzcl(ii))*(bb1+bb2)*0.5d0)
     enddo
   endif

!        Intégration sur le deuxième intervalle de (entre s1 et s2)

   bb1 = max(f4m1(iel),f4s1)
   bb2 = min(f4m2(iel),f4s2)
   if( bb2.gt.bb1 ) then
     do ii = 1,7
       zz(ii) = zz(ii) + hrec(iel)*(bb2-bb1)/(f4s2-f4s1)          &
             * ( (zzs1(ii)*f4s2-zzs2(ii)*f4s1)                    &
           +   (zzs2(ii)-zzs1(ii))*(bb1+bb2)*0.5d0)
     enddo
   endif

!        Intégration de s2 à s3

   bb1 = max(f4m1(iel),f4s2)
   bb2 = min(f4m2(iel),f4s3)
   if( bb2.gt.bb1 ) then
     do ii = 1,7
       zz(ii) = zz(ii) + hrec(iel)*(bb2-bb1)/(f4s3-f4s2)          &
           * ( (zzs2(ii)*f4s3-zzs3(ii)*f4s2)                      &
           +   (zzs3(ii)-zzs2(ii))*(bb1+bb2)*0.5d0)
     enddo
   endif

!       Intégration de s3 à f4

   bb1 = max(f4m1(iel),f4s3)
   bb2 = min(f4m2(iel),1.d0)
   if ( bb2.gt.bb1 ) then
     do ii = 1,7
      zz(ii) = zz(ii) + hrec(iel)*(bb2-bb1)/(1.d0-f4s3)           &
             * ( (zzs3(ii)*1.d0-zzoxy(ii)*f4s3)                   &
                +(zzoxy(ii)-zzs3(ii))*(bb1+bb2)*0.5d0)
     enddo
   endif

   if ( ieqco2 .eq. 1 ) then

!  Transport de CO2

     if ( (zz(ico) - (zco2t-zz(ico2))) .lt. -epzero ) then
       nclco1 = nclco1 +1
       cmina1 = min((zz(ico)- (zco2t-zz(ico2))),cmina1)
       cmaxa1 = max((zz(ico)- (zco2t-zz(ico2))),cmaxa1)
     endif
     if ( (zz(io2)-anu3*(zco2t-zz(ico2))) .lt. -epzero ) then
       nclo21 = nclo21 +1
       cmina2 = min((zz(io2) -anu3*(zco2t-zz(ico2))),cmina2)
       cmaxa2 = max((zz(io2) -anu3*(zco2t-zz(ico2))),cmaxa2)
     endif

     reac3 = min(max(zco2t-zz(ico2),0.d0),zz(ico),                &
                                          zz(io2)/anu3)

     zz(ico ) = zz(ico) -      reac3
     zz(io2 ) = zz(io2) - anu3*reac3
     zz(ico2) = zz(ico2)+      reac3

!           RTP(IEL,ISCA(IYCO2))= ZCO2T*(1.D0+W1(IEL))*WMCO2

   else if ( ieqco2 .ne. 0 ) then

     write(nfecra,3000) ieqco2
     call csexit(1)

   endif

   fuel1(iel) = zz(ichx1) * wmchx1
   fuel2(iel) = zz(ichx2) * wmchx2
   fuel3(iel) = zz(ico  ) * wmco
   prod1(iel) = zz(ico2 ) * wmco2
   prod2(iel) = zz(ih2o ) * wmh2o
   oxyd(iel)  = zz(io2  ) * wmo2
   xiner(iel) = zz(in2  ) * wmn2

  endif

enddo

n2 = 0
do iel = 1, ncel
  if ( fuel3(iel) .lt. (-epzero) ) then
    n2= n2 + 1
  endif
enddo

if ( ieqco2 .eq. 1 .or. ieqco2 .eq.2 ) then


  if (irangp.ge.0) then

    call parcpt(n2)
    call parcpt(nclco)
    call parcpt(nclco1)
    call parcpt(nclo2)
    call parcpt(nclo21)

    call parmin(cmins1)
    call parmin(cmina1)
    call parmin(cmins2)
    call parmin(cmina2)

    call parmax(cmaxs1)
    call parmax(cmaxa1)
    call parmax(cmaxs2)
    call parmax(cmaxa2)

  endif

  WRITE(NFECRA,*) ' Point a CO < 0 ',N2
  WRITE(NFECRA,*) ' Clipping CO2 - CO sans avec : ',NCLCO,NCLCO1
  WRITE(NFECRA,*) ' Min Max sans                : ',CMINS1,CMAXS1
  WRITE(NFECRA,*) ' Min Max avec                : ',CMINA1,CMAXA1
  WRITE(NFECRA,*) ' Clipping CO2 - O2 sans avec : ',NCLO2,NCLO21
  WRITE(NFECRA,*) ' Min Max sans                : ',CMINS2,CMAXS2
  WRITE(NFECRA,*) ' Min Max avec                : ',CMINA2,CMAXA2

endif

!===============================================================================
! 5.  IMPRESSION
!===============================================================================

n2  = 0
n3  = 0
n4  = 0
n5  = 0
n6  = 0
n7  = 0
n8  = 0
n9  = 0
n10 = 0

! --> Controle des parametres de la pdf

n1 = ncel

do iel = 1, ncel
  if ( indpdf(iel).ne.0 ) then
    n2 = n2 +1
  endif
enddo

! --> Controle des differentes valeurs des fractions massiques

sommin = 1.d+20
sommax =-1.d+20
do iel = 1, ncel

  somm = fuel1(iel) + fuel2(iel) + fuel3(iel)                     &
       + oxyd(iel)                                                &
       + prod1(iel) + prod2(iel) + xiner(iel)

  sommin = min(sommin,somm)
  sommax = max(sommax,somm)

  if ( abs(somm-1.d0).lt.epsicp )   then
    n3  = n3 +1
  endif

  if ( fuel1(iel).lt.(-epzero) .or.                               &
       fuel1(iel).gt.(1.d0+epzero) ) n4 = n4+1
  if ( fuel2(iel).lt.(-epzero) .or.                               &
       fuel2(iel).gt.(1.d0+epzero) ) n5 = n5+1
  if ( fuel3(iel).lt.(-epzero) .or.                               &
       fuel3(iel).gt.(1.d0+epzero) ) n6 = n6+1
  if ( oxyd(iel).lt.(-epzero) .or.                                &
       oxyd(iel).gt.(1.d0+epzero)  ) n7 = n7+1
  if ( xiner(iel).lt.(-epzero) .or.                               &
       xiner(iel).gt.(1.d0+epzero) ) n8 = n8+1
  if ( prod1(iel).lt.(-epzero) .or.                               &
       prod1(iel).gt.(1.d0+epzero) ) n9 = n9+1
  if ( prod2(iel).lt.(-epzero) .or.                               &
       prod2(iel).gt.(1.d0+epzero) ) n10 = n10+1

enddo

n1 = ncel

if (irangp.ge.0) then
  inttmp( 1) = n1
  inttmp( 2) = n2
  inttmp( 3) = n3
  inttmp( 4) = n4
  inttmp( 5) = n5
  inttmp( 6) = n6
  inttmp( 7) = n7
  inttmp( 8) = n8
  inttmp( 9) = n9
  inttmp(10) = n10
  call parism(nbrint,inttmp)
  !==========
endif

write(nfecra,1000) n1 , n2

write(nfecra,2200) n3 , n4, n5, n6, n7, n8,                       &
                        n9, n10

if ( irangp .ge. 0 ) then
  call parmin(sommin)
  call parmax(sommax)
endif

write(NFECRA,*) ' Somme Min MAX = ',SOMMIN,SOMMAX

!-------
! FORMAT
!-------

 1000 format (/,                                                  &
'MODELISATION DE LA COMBUSTION AVEC LE MODELE DE DIFFUSION ',     &
'TURBULENTE (CPCYM2)',/,                                    &
'CHIMIE RAPIDE A 3 CORPS - EXTENSION A 3 COMBUSTIBLES ',          &
'(Application au FUEL)',/,                                  &
'==========================================================',     &
'==================',/,                                     &
' Nb de points de calculs                                     = ',&
   i9,/,                                                    &
' Nb de points turbulents (passage par les PDF)               = ',&
   i9)

 2200 format(/,                                                   &
'CONTROLE DES VALEURS DES FRACTIONS MASSIQUES',/,           &
' Nb de points de calculs qui respectent somme des Yi = 1    = ', &
   i9,/,                                                    &
' Nb de points YCHX1 ,YCHX2    < 0 ou > 1 = ',I9,I9,/,      &
' Nb de points YC0             < 0 ou > 1 = ',I9,/,         &
' Nb de points YO2  ,YN2       < 0 ou > 1 = ',I9,I9,/,      &
' Nb de points YCO2 ,YH2O      < 0 ou > 1 = ',I9,I9,/)

 3000 format(/,                                                   &
'MODELE DE CO CHARBON ACTIVE : ',/,                         &
'       AVEC  IEQCO2 = ',I2,/,                              &
'HORS SEUL LES OPTIONS 0 , 1 et 2 SONT DISPONIBLES',/,      &
'    ARRET DU CALCUL : VERIFIER USPPMO')


return
end subroutine
