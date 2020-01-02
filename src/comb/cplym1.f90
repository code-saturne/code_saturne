!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

subroutine cplym1 &
!================

 ( ncelet , ncel   , nitbmc , nrtbmc ,                            &
   f1m    , f2m    , f3m    , f4m    ,                            &
   indpdf , si7    , si8    , sp2m   , f4i7   ,                   &
   dsi7   , dsi8   , sdeb   , sfin   , haut   ,                   &
   fuel1  , fuel2  , fuel3  , oxyd   , prod1  , prod2  ,          &
   xiner  ,                                                       &
   itbmc  , rtbmc  ,                                              &
   ipi7i8 ,                                                       &
   sp     , si     , sr     , w1     )

!===============================================================================
! FONCTION :
! --------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
!   --------------------------------------------------------------

!    ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE

!      COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
!      TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON

!      CALCUL DES CONCENTRATIONS GAZEUSES MOYENNES
!       VALEURS CELLULES
!       ----------------

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nitbmc           ! e  ! <-- ! taille du macro tableau mc entiers             !
! nrtbmc           ! e  ! <-- ! taille du macro tableau mc reels               !
! f1m              ! tr ! <-- ! moyenne du traceur 1 mvl [chx1m+co]            !
! f2m              ! tr ! <-- ! moyenne du traceur 2 mvl [chx2m+co]            !
! f3m              ! tr ! <-- ! moyenne du traceur 3 (co c.het)                !
! f4m              ! tr ! <-- ! moyenne du traceur 4 (air)                     !
! indpdf           ! te ! <-- ! passage par les pdf                            !
! si7              ! tr ! <-- ! abscisse curviligne au point i7                !
! si8              ! tr ! <-- ! abscisse curviligne au point i8                !
! sp2m             ! tr ! <-- ! variance droite support                        !
! f4i7             ! tr ! <-- ! f4 au point i7                                 !
! dsi7             ! tr ! <-- ! dirac au point i7                              !
! dsi8             ! tr ! <-- ! dirac au point i8                              !
! sdeb             ! tr ! <-- ! abscisse debut rectangle                       !
! sfin             ! tr ! <-- ! abscisse fin rectangle                         !
! hrec             ! tr ! <-- ! hauteur rectangle                              !
! fuel1            ! tr !  <- ! fraction massique chx1m                        !
! fuel2            ! tr !  <- ! fraction massique chx2m                        !
! fuel3            ! tr !  <- ! fraction massique co                           !
! oxyd             ! tr !  <- ! fraction massique o2                           !
! prod1            ! tr !  <- ! fraction massique co2                          !
! prod2            ! tr !  <- ! fraction massique h2o                          !
! xiner            ! tr !  <- ! fraction massique n2                           !
! itbmc            ! tr ! <-- ! macro tableau entier mc travail                !
! rtbmc            ! tr ! <-- ! macro tableau reel   mc travail                !
! ipi7i8           ! te ! <-- ! position des points i7 et i8 dans              !
!                  !    !     ! les differents domaines de richesse            !
! sp               ! tr ! <-- ! delimitation sur le support des                !
!                  !    !     ! domaines pauvre et intermediaire               !
! si               ! tr ! <-- ! delimitation sur le support des                !
!                  !    !     ! domaines intermediaire et riche                !
! sr               ! tr ! <-- ! delimitation sur le support des                !
!                  !    !     ! domaines riche et hyper-riche                  !
! w1               ! tr ! <-- ! tableau de travail                             !
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
use field

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel
integer          nitbmc , nrtbmc
integer          itbmc(ncelet,nitbmc)
integer          indpdf(ncelet) , ipi7i8(ncelet)

double precision f1m(ncelet)  , f2m(ncelet)
double precision f3m(ncelet)  , f4m(ncelet)
double precision dsi7(ncelet) , dsi8(ncelet) , sdeb(ncelet)
double precision sfin(ncelet) , haut(ncelet)
double precision si7(ncelet)  , si8(ncelet)  , sp2m(ncelet)
double precision f4i7(ncelet)
double precision fuel1(ncelet), fuel2(ncelet) , fuel3(ncelet)
double precision oxyd(ncelet) , prod1(ncelet) , prod2(ncelet)
double precision xiner(ncelet)
double precision rtbmc(ncelet,nrtbmc)
double precision sp(ncelet) , si(ncelet) , sr(ncelet)
double precision w1(ncelet)

! Local variables

integer          iel, icha, itbr, itbi, ige

integer               n1 , n2 , n3 , n4 , n5 , n6 , n7 , n8 , n9
integer          n10, n11, n12, n13, n14, n15, n16, n17, n18, n19
integer          n20, n21, n22, n23, n24, n25, n26, n27, n28, n29
integer          n30, n31, n32
integer          nbrint
parameter        (nbrint = 32)
integer          inttmp(nbrint)

double precision zcof10 , zcof20
double precision cx1m   , cx2m
double precision wmchx1 , wmchx2 , wmco   , wmco2 , wmh2o
double precision wmo2   , wmn2   , wmc
double precision achx1f1, acof1  , achx2f2, acof2
double precision cc1    , cc2    , cc3    , cc4
double precision f1i7   , f2i7   , f3i7
double precision vars
double precision zc(ngazem)
double precision ac0(ngazem) , bc0(ngazem)
double precision ac1(ngazem) , bc1(ngazem)
double precision ac2(ngazem) , bc2(ngazem)
double precision ac3(ngazem) , bc3(ngazem)
double precision ai(ngazem)  , bi(ngazem)
double precision anu1 , anu2 , anu3
double precision s1, s2, s3, sm1, sm2
double precision den1   , den2   , somm
#ifdef DEBUG
double precision sommf  , mom0   , mom1   , mom2
#endif
double precision zchx10 , zchx20 , zco0   , zco20
double precision zn20   , zh2o0  , zo20
double precision ychx10 , ychx20
#ifdef DEBUG
double precision yco0   , yco20
double precision yn20   , yh2o0  , yo20
#endif
double precision zchx11 , zchx21 , zco1   , zco21
double precision zn21   , zh2o1  , zo21
double precision zchx12 , zchx22 , zco2   , zco22
double precision zn22   , zh2o2  , zo22
double precision zchx13 , zchx23 , zco3   , zco23
double precision zn23   , zh2o3  , zo23
double precision reac1  , reac2  , reac3

double precision, dimension(:), pointer :: cvar_f1m, cvar_f2m

!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings

s1 = 0.d0

! --- Initialisation tableau de travail

do iel = 1, ncel
  w1(iel) = zero
  do itbr = 1, nrtbmc
    rtbmc(iel,itbr) = zero
  enddo
  do itbi = 1, nitbmc
    itbmc(iel,itbi) = 0
  enddo
enddo

!===============================================================================
! 2. CALCULS PRELIMINAIRES
!===============================================================================

! --- Calcul de F1M de chaque charbon dans RTBMC
!            de F2M de chaque charbon dans RTBMC

! ------ RTBMC(IF1MC(ICHA)) = F1M(ICHA)
!        RTBMC(IF2MC(ICHA)) = F2M(ICHA)

do icha = 1, ncharb
  call field_get_val_s(ivarfl(isca(if1m(icha))), cvar_f1m)
  call field_get_val_s(ivarfl(isca(if2m(icha))), cvar_f2m)
  do iel = 1, ncel
    rtbmc(iel,if1mc(icha)) = cvar_f1m(iel)
    rtbmc(iel,if2mc(icha)) = cvar_f2m(iel)
  enddo
enddo


! --- Masses molaires

wmco   = wmole(ico  )
wmo2   = wmole(io2  )
wmco2  = wmole(ico2 )
wmh2o  = wmole(ih2o )
wmn2   = wmole(in2  )
wmc    = wmolat(iatc)


! --- Calcul de CX1M, CX2M
!     et des coefficient relatifs a l'expression de
!        ZCHx1m0, ZCHx2m0, ZCO0, ZO20 et ZN20
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
  do icha = 1, ncharb
    den1   = 1.d0                                                 &
           / ( a1(icha)*wmole(ichx1c(icha))+b1(icha)*wmole(ico))
    ychx10 = ychx10 + den1 *                                      &
           ( rtbmc(iel,if1mc(icha))*a1(icha)*wmole(ichx1c(icha)) )
    zchx10 = zchx10 + den1 *                                      &
           ( rtbmc(iel,if1mc(icha))*a1(icha) )
    zcof10 = zcof10 + den1 *                                      &
           ( rtbmc(iel,if1mc(icha))*b1(icha) )
    den2   = 1.d0                                                 &
           / ( a2(icha)*wmole(ichx2c(icha))+b2(icha)*wmole(ico))
    ychx20 = ychx20 + den2 *                                      &
           ( rtbmc(iel,if2mc(icha))*a2(icha)*wmole(ichx2c(icha)) )
    zchx20 = zchx20 + den2 *                                      &
           ( rtbmc(iel,if2mc(icha))*a2(icha) )
    zcof20 = zcof20 + den2 *                                      &
           ( rtbmc(iel,if2mc(icha))*b2(icha) )
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
  if ( f1m(iel).gt.epzero ) then
    rtbmc(iel,ichx1f1) = zchx10 / f1m(iel)
    rtbmc(iel,icof1  ) = zcof10 / f1m(iel)
  endif
  rtbmc(iel,ichx2f2) = 0.d0
  rtbmc(iel,icof2  ) = 0.d0
  if ( f2m(iel).gt.epzero ) then
    rtbmc(iel,ichx2f2) = zchx20 / f2m(iel)
    rtbmc(iel,icof2  ) = zcof20 / f2m(iel)
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
  ipi7i8(iel) = 0
  sp(iel)     = zero
  si(iel)     = zero
  sr(iel)     = zero
enddo



!===============================================================================
! 3. CALCUL DE LA COMPOSITION DU MELANGE SANS LES PDF
!    SI LES FLUCTUATIONS DE S SONT TROP FAIBLES
!    SI SI7 ET SI8 SONT TROP PROCHE DE 0
!    ALORS F1 = F1M, F2 = F2M , F3 = F3M ET F4 = F4M
!===============================================================================

do iel = 1, ncel

  if ( indpdf(iel).eq.0 ) then

! --> Calculs preliminaires

    cx1m    = rtbmc(iel,ix1mc)
    cx2m    = rtbmc(iel,ix2mc)
    wmchx1  = wmolat(iatc) + cx1m*wmolat(iath)
    wmchx2  = wmolat(iatc) + cx2m*wmolat(iath)
    achx1f1 = rtbmc(iel,ichx1f1)
    achx2f2 = rtbmc(iel,ichx2f2)
    acof1   = rtbmc(iel,icof1)
    acof2   = rtbmc(iel,icof2)

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
    zo20    = (ao2f4*f4m(iel)+ao2f3*f3m(iel))
    zn20    = (an2f4*f4m(iel))
    zco20   = zero
    zh2o0   = zero
! ---- Test
#ifdef DEBUG
    ychx10  = zchx10*wmchx1
    ychx20  = zchx20*wmchx2
    yco0    = zco0 * wmco
    yo20    = zo20*wmo2
    yn20    = zn20*wmn2
    yco20   = zco20*wmco2
    yh2o0   = zh2o0*wmh2o
    somm    = ychx10 + ychx20 + yco0 + yo20 + yn20 +              &
              yco20  + yh2o0
    if (abs(somm-1.d0).gt.epsicp) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'SOMME0', SOMM
      write(NFECRA,*) 'F1M',F1M(IEL)
      write(NFECRA,*) 'F2M',F2M(IEL)
      write(NFECRA,*) 'F3M',F3M(IEL)
      write(NFECRA,*) 'F4M',F4M(IEL)
    endif
#endif
! ---- Fin test

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
    reac3  = min(zco2, (zo22/anu3))
    zchx13 = zchx12
    zchx23 = zchx22
    zo23   = zo22  - anu3*reac3
    zco3   = zco2  -      reac3
    zco23  = zco22 +      reac3
    zh2o3  = zh2o2
    zn23   = zn22

    fuel1(iel) = zchx13 * wmchx1
    fuel2(iel) = zchx23 * wmchx2
    fuel3(iel) = zco3   * wmco
    prod1(iel) = zco23  * wmco2
    prod2(iel) = zh2o3  * wmh2o
    oxyd(iel)  = zo23   * wmo2
    xiner(iel) = zn23   * wmn2

! ---- Test
#ifdef DEBUG
    somm = fuel1(iel) + fuel2(iel) + fuel3(iel)                   &
         + prod1(iel) + prod2(iel)                                &
         + oxyd(iel)  + xiner(iel)
!          IF (ABS(SOMM-1.D0).GT.EPSICP) THEN
    if ( (abs(somm-1.d0).gt.epsicp) .or.                          &
         (prod2(iel).lt.(-epzero) .or.                            &
          prod2(iel).gt.(1.d0+epzero)) ) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'INDPDF0',INDPDF(IEL)
      write(NFECRA,*) 'F1M',F1M(IEL)
      write(NFECRA,*) 'F2M',F2M(IEL)
      write(NFECRA,*) 'F3M',F3M(IEL)
      write(NFECRA,*) 'F4M',F4M(IEL)
      sommf = f1m(iel) + f2m(iel) + f3m(iel) + f4m(iel)
      write(NFECRA,*) 'SOMME F', SOMMF
      write(NFECRA,*) 'AO2F4, AO2F3 ', AO2F4, AO2F3
      write(NFECRA,*) 'ZCHX1, 0,1,2,3',                           &
                       zchx10,zchx11,zchx12,zchx13
      write(NFECRA,*) 'ZCHX2, 0,1,2,3',                           &
                       zchx20,zchx21,zchx22,zchx23
      write(NFECRA,*) 'ZCO, 0,1,2,3',                             &
                       zco0,zco1,zco2,zco3
      write(NFECRA,*) 'ZCO2, 0,1,2,3',                            &
                       zco20,zco21,zco22,zco23
      write(NFECRA,*) 'ZH2O, 0,1,2,3',                            &
                       zh2o0,zh2o1,zh2o2,zh2o3
      write(NFECRA,*) 'ZO2, 0,1,2,3',                             &
                       zo20,zo21,zo22,zo23
      write(NFECRA,*) 'ZN2, 0,1,2,3',                             &
                       zn20,zn21,zn22,zn23
      write(NFECRA,*) 'Fuel1',FUEL1(IEL)
      write(NFECRA,*) 'Fuel2',FUEL2(IEL)
      write(NFECRA,*) 'Fuel3',FUEL3(IEL)
      write(NFECRA,*) 'Oxyd',OXYD(IEL)
      write(NFECRA,*) 'Iner',XINER(IEL)
      write(NFECRA,*) 'Prod1',PROD1(IEL)
      write(NFECRA,*) 'Prod2',PROD2(IEL)
      write(NFECRA,*) 'SOMME Y', SOMM
    endif
#endif
! ---- Fin Test

  endif

enddo



!===============================================================================
! 4. CALCUL DE LA COMPOSITION DU MELANGE AVEC LES PDF
!    ATTENTION : ON SE LIMITE POUR L'INSTANT A INDPDF = 3 SUPPORT (I4,M)
!===============================================================================


! -->  I7 = I4 : F1I7 = F2I7 = F3I7 = 0
!                F4I7 = 1

! -->  RELATIONS SUR LE SUPPORT DE LA PDF
!        FJ = FJM + CJ*S avec CJ = (FJM-FJI7)/SI7 (J = 1 a 4)

! Devectorisation de la boucle pour VPP 5000

!CDIR NOVECTOR

do iel = 1, ncel

  if ( indpdf(iel).eq.3 ) then

! --> Calculs preliminaires

    cx1m    = rtbmc(iel,ix1mc)
    cx2m    = rtbmc(iel,ix2mc)
    wmchx1  = wmolat(iatc) + cx1m*wmolat(iath)
    wmchx2  = wmolat(iatc) + cx2m*wmolat(iath)
    achx1f1 = rtbmc(iel,ichx1f1)
    achx2f2 = rtbmc(iel,ichx2f2)
    acof1   = rtbmc(iel,icof1)
    acof2   = rtbmc(iel,icof2)

! ---- Calcul de F1I7, F2I7, F3I7

    f1i7   = zero
    f2i7   = zero
    f3i7   = 1.d0 - f4i7(iel)

! ----- Calcul de C1, C2, C3 et C4

    cc1 = (f1i7     -f1m(iel)) / si7(iel)
    cc2 = (f2i7     -f2m(iel)) / si7(iel)
    cc3 = (f3i7     -f3m(iel)) / si7(iel)
    cc4 = (f4i7(iel)-f4m(iel)) / si7(iel)

! ---- Rq : les autres pointeurs sont deja definis dans cplecd.F

! --> Initialisation

    do ige = 1, (ngaze-2*ncharb)
      zc(ige) = zero
    enddo

! --> Calculs des coefficients avant combustion phase gaz
!       ZC(IGE) = AC0(IGE) + BC0(IGE) * S

    ac0(ichx1) = achx1f1*f1m(iel)
    bc0(ichx1) = achx1f1*cc1
    ac0(ichx2) = achx2f2*f2m(iel)
    bc0(ichx2) = achx2f2*cc2
    ac0(ico  ) = acof1*f1m(iel)+acof2*f2m(iel) +                  &
                 acof3*f3m(iel)
    bc0(ico  ) = acof1*cc1     +acof2*cc2      +                  &
                 acof3*cc3
    ac0(ico2 ) = zero
    bc0(ico2 ) = zero
    ac0(ih2o ) = zero
    bc0(ih2o ) = zero
    ac0(io2  ) = ao2f3*f3m(iel)+ao2f4*f4m(iel)
    bc0(io2  ) = ao2f3*cc3     +ao2f4*cc4
    ac0(in2  ) = an2f4*f4m(iel)
    bc0(in2  ) = an2f4*cc4

! ---- Test
#ifdef DEBUG
    wmole(ichx1) = wmchx1
    wmole(ichx2) = wmchx2
    somm  = 0.d0
    sommf = 0.d0
    do ige = 1, (ngaze-2*ncharb)
      somm  = somm  + ac0(ige)*wmole(ige)
      sommf = sommf + bc0(ige)*wmole(ige)
    enddo
    if (abs(somm-1.d0).gt.epsicp) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'F4I7 ', F4I7(IEL)
      write(NFECRA,*) 'SOMM AC0', SOMM
    endif
    if (abs(sommf).gt.epsicp) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'SOMM BC0', SOMMF
    endif
#endif
! ---- Fin Test

! --> Les pics de Dirac sont situes dans des domaines
!     ou la reaction n'est pas possible (points froids)

    do ige = 1, (ngaze-2*ncharb)
      zc(ige) = dsi7(iel)*(ac0(ige)+bc0(ige)*si7(iel))            &
              + dsi8(iel)*(ac0(ige)+bc0(ige)*si8(iel))
    enddo

! --> Initialisation des coefficients avant integration

    do ige = 1, (ngaze-2*ncharb)
      ac1(ige) = ac0(ige)
      bc1(ige) = bc0(ige)
      ac2(ige) = ac0(ige)
      bc2(ige) = bc0(ige)
      ac3(ige) = ac0(ige)
      bc3(ige) = bc0(ige)
    enddo

! --> Integration dans le domaine hyper riche

! Rq : On rajoute le test suivant CAR
!      si F1M = 0 alors le domaine hyper riche n'existe pas

    if ( f1m(iel).gt.epsicp ) then

    anu1 = 0.25d0*(cx1m-cx2m)
    s1   = - (ac0(io2)-ac0(ichx1)*anu1)                           &
           / (bc0(io2)-bc0(ichx1)*anu1)
    sm1  = max(sdeb(iel),s1)
    sm2  = sfin(iel)

    do ige = 1, (ngaze-2*ncharb)
      ai(ige) = ac0(ige)
      bi(ige) = bc0(ige)
    enddo
    bi(ichx1) = (ac0(ichx1)+bc0(ichx1)*si8(iel))/(si8(iel)-s1)
    ai(ichx1) = -s1*bi(ichx1)
    ai(ichx2) = ac0(ichx2) +           (ac0(ichx1)-ai(ichx1))
    bi(ichx2) = bc0(ichx2) +           (bc0(ichx1)-bi(ichx1))
    ai(io2  ) = ac0(io2  ) -      anu1*(ac0(ichx1)-ai(ichx1))
    bi(io2  ) = bc0(io2  ) -      anu1*(bc0(ichx1)-bi(ichx1))
    ai(ih2o ) = ac0(ih2o ) + 2.d0*anu1*(ac0(ichx1)-ai(ichx1))
    bi(ih2o ) = bc0(ih2o ) + 2.d0*anu1*(bc0(ichx1)-bi(ichx1))

    if ( sm1.lt.sm2 ) then

      do ige = 1, (ngaze-2*ncharb)
        zc(ige) = zc(ige)                                         &
          + haut(iel)*(sm2-sm1)*(ai(ige)+.5d0*bi(ige)*(sm1+sm2))
      enddo

    endif

    do ige = 1, (ngaze-2*ncharb)
      ac1(ige) = ( ac0(ige)*s1-ai(ige)*si7(iel) +                 &
                     s1*si7(iel)*(bc0(ige)-bi(ige)) )             &
                 / ( s1-si7(iel) )
      bc1(ige) = (ai(ige)+bi(ige)*s1-ac1(ige)) / s1
    enddo

! ---- Apres la premiere reaction, on se trouve hors du domaine hyper riche
!      CHX1 ne doit plus exister AC1(ICHX1) = BC1(ICHX1) = 0

! ---- Test
#ifdef DEBUG
    somm  = 0.d0
    sommf = 0.d0
    do ige = 1, (ngaze-2*ncharb)
      somm  = somm  + ac1(ige)*wmole(ige)
      sommf = sommf + bc1(ige)*wmole(ige)
    enddo
      if ( abs(ac1(ichx1)).gt.epsicp .or.                         &
         abs(bc1(ichx1)).gt.epsicp ) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'AC1(ICHX1)', AC1(ICHX1)
      write(NFECRA,*) 'BC1(ICHX1)', BC1(ICHX1)
    endif
      if ( abs(somm-1.d0).gt.epsicp ) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'SOMM AC1', SOMM
    endif
      if ( abs(sommf).gt.epsicp ) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'SOMM BC1', SOMMF
    endif
      if ( abs(ac1(io2)+bc1(io2)*s1).gt.epsicp ) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'AC01(IO2)', AC0(IO2),AC1(IO2)
      write(NFECRA,*) 'BC01(IO2)', BC0(IO2),BC1(IO2)
    endif
#endif
! ---- Fin Test

    endif

! --> Integration dans le domaine riche

    anu2 = 0.25d0*(2.d0 + cx2m)
    s2   = - (ac1(io2)-ac1(ichx2)*anu2)                           &
           / (bc1(io2)-bc1(ichx2)*anu2)

    sm1  = max(sdeb(iel),s2)
    sm2  = min(sfin(iel),s1)

    do ige = 1, (ngaze-2*ncharb)
      ai(ige) = ac1(ige)
      bi(ige) = bc1(ige)
    enddo

    bi(ichx2) = (ac1(ichx2)+bc1(ichx2)*s1)/(s1-s2)
    ai(ichx2) = -s2*bi(ichx2)
    ai(io2  ) = ac1(io2 ) -      anu2*(ac1(ichx2)-ai(ichx2))
    bi(io2  ) = bc1(io2 ) -      anu2*(bc1(ichx2)-bi(ichx2))
    ai(ico  ) = ac1(ico ) +           (ac1(ichx2)-ai(ichx2))
    bi(ico  ) = bc1(ico ) +           (bc1(ichx2)-bi(ichx2))
    ai(ih2o ) = ac1(ih2o) + .5d0*cx2m*(ac1(ichx2)-ai(ichx2))
    bi(ih2o ) = bc1(ih2o) + .5d0*cx2m*(bc1(ichx2)-bi(ichx2))

    if ( sm1.lt.sm2 ) then

      do ige = 1, (ngaze-2*ncharb)
        zc(ige) = zc(ige)                                         &
          + haut(iel)*(sm2-sm1)*(ai(ige)+.5d0*bi(ige)*(sm1+sm2))
      enddo

    endif

    do ige = 1, (ngaze-2*ncharb)
!              AC2(IGE) = (AC1(IGE)*S2-AI(IGE)*SI7(IEL)) +
!     &                    S2*SI7(IEL)*(BC1(IGE)-BI(IGE))/(S2-SI7(IEL))
!              BC2(IGE) = (AI(IGE)+BI(IGE)*S2-AC2(IGE))/(S2-SI7(IEL))
      ac2(ige) = ( ac1(ige)*s2-ai(ige)*si7(iel) +                 &
                     s2*si7(iel)*(bc1(ige)-bi(ige)) )             &
                 / ( s2-si7(iel) )
      bc2(ige) = (ai(ige)+bi(ige)*s2-ac2(ige)) / s2
    enddo

! ---- Apres la seconde reaction, on se trouve hors du domaine riche
!      CHX1 et CHX2 ne doivent plus exister
!      AC2(ICHX1) = BC2(ICHX1) = AC2(ICHX2) = BC2(ICHX2) = 0

! ---- Test
#ifdef DEBUG
    somm  = 0.d0
    sommf = 0.d0
    do ige = 1, (ngaze-2*ncharb)
      somm  = somm  + ac2(ige)*wmole(ige)
      sommf = sommf + bc2(ige)*wmole(ige)
    enddo
    if (abs(ac2(ichx1)).gt.epsicp .or.                            &
        abs(bc2(ichx1)).gt.epsicp .or.                            &
        abs(ac2(ichx2)).gt.epsicp .or.                            &
        abs(bc2(ichx2)).gt.epsicp      ) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'AC2(ICHX1)', AC2(ICHX1)
      write(NFECRA,*) 'BC2(ICHX1)', BC2(ICHX1)
      write(NFECRA,*) 'AC2(ICHX2)', AC2(ICHX2)
      write(NFECRA,*) 'BC2(ICHX2)', BC2(ICHX2)
    endif
    if (abs(somm-1.d0).gt.epsicp) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'SOMM AC2', SOMM
    endif
    if (abs(sommf).gt.epsicp) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'SOMM BC2', SOMMF
    endif
    if (abs(ac2(io2)+bc2(io2)*s2).gt.epsicp) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'AC12(IO2)', AC1(IO2),AC2(IO2)
      write(NFECRA,*) 'BC12(IO2)', BC1(IO2),BC2(IO2)
    endif
#endif
! ---- Fin Test

! --> Integration dans le domaine intermediaire

    anu3 = 0.5d0
    s3   = - (ac2(io2)-ac2(ico)*anu3)                             &
         / (bc2(io2)-bc2(ico)*anu3)

    sm1  = max(sdeb(iel),s3)
    sm2  = min(sfin(iel),s2)

    do ige = 1, (ngaze-2*ncharb)
      ai(ige) = ac2(ige)
      bi(ige) = bc2(ige)
    enddo

    bi(ico)   = (ac2(ico)+bc2(ico)*s2)/(s2-s3)
    ai(ico)   = -s3*bi(ico)
    ai(io2)   = ac2(io2 ) - anu3*(ac2(ico)-ai(ico))
    bi(io2)   = bc2(io2 ) - anu3*(bc2(ico)-bi(ico))
    ai(ico2)  = ac2(ico2) +      (ac2(ico)-ai(ico))
    bi(ico2)  = bc2(ico2) +      (bc2(ico)-bi(ico))

    if ( sm1.lt.sm2 ) then

      do ige = 1, (ngaze-2*ncharb)
        zc(ige) = zc(ige)                                         &
          + haut(iel)*(sm2-sm1)*(ai(ige)+.5d0*bi(ige)*(sm1+sm2))
      enddo

    endif

    do ige = 1, (ngaze-2*ncharb)
      ac3(ige) = ( ac2(ige)*s3-ai(ige)*si7(iel) +                 &
                     s3*si7(iel)*(bc2(ige)-bi(ige)) )             &
                 / ( s3-si7(iel) )
      bc3(ige) = (ai(ige)+bi(ige)*s3-ac3(ige)) / s3
    enddo

! ---- Apres la troisieme reaction, om se trouve dans le domaine pauvre
!      CHX1,CHX2 et CO doivent etre nuls
!      AC3(ICHX1) = BC3(ICHX1) = 0
!      AC3(ICHX2) = BC3(ICHX2) = 0
!      AC3(ICO  ) = BC3(ICO  ) = 0

! ---- Test
#ifdef DEBUG
    somm  = 0.d0
    sommf = 0.d0
    do ige = 1, (ngaze-2*ncharb)
      somm  = somm  + ac3(ige)*wmole(ige)
      sommf = sommf + bc3(ige)*wmole(ige)
    enddo
    if (abs(ac3(ichx1)).gt.epsicp .or.                            &
        abs(bc3(ichx1)).gt.epsicp .or.                            &
        abs(ac3(ichx2)).gt.epsicp .or.                            &
        abs(bc3(ichx2)).gt.epsicp .or.                            &
        abs(ac3(ico)).gt.epsicp .or.                              &
        abs(bc3(ico)).gt.epsicp     ) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'AC3(ICHX1)', AC3(ICHX1)
      write(NFECRA,*) 'BC3(ICHX1)', BC3(ICHX1)
      write(NFECRA,*) 'AC3(ICHX2)', AC3(ICHX2)
      write(NFECRA,*) 'BC3(ICHX2)', BC3(ICHX2)
      write(NFECRA,*) 'AC3(ICO)', AC3(ICO)
      write(NFECRA,*) 'BC3(ICO)', BC3(ICO)
    endif
    if (abs(somm-1.d0).gt.epsicp) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'SOMM AC3', SOMM
    endif
    if (abs(sommf).gt.epsicp) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'SOMM BC3', SOMMF
    endif
    if (abs(ac3(io2)+bc3(io2)*s3).gt.epsicp) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'AC23(IO2)', AC2(IO2),AC3(IO2)
      write(NFECRA,*) 'BC23(IO2)', BC2(IO2),BC3(IO2)
    endif
#endif
! ---- Fin Test

! --> Integration dans le domaine pauvre

    sm1 = sdeb(iel)
    sm2  = min(sfin(iel),s3)

    do ige = 1, (ngaze-2*ncharb)
      ai(ige) = ac3(ige)
      bi(ige) = bc3(ige)
    enddo

    if ( sm1.lt.sm2 ) then

      do ige = 1, (ngaze-2*ncharb)
        zc(ige) = zc(ige)                                         &
          + haut(iel)*(sm2-sm1)*(ai(ige)+.5d0*bi(ige)*(sm1+sm2))
      enddo

    endif

    fuel1(iel) = zc(ichx1) * wmchx1
    fuel2(iel) = zc(ichx2) * wmchx2
    fuel3(iel) = zc(ico  ) * wmco
    prod1(iel) = zc(ico2 ) * wmco2
    prod2(iel) = zc(ih2o ) * wmh2o
    oxyd(iel)  = zc(io2  ) * wmo2
    xiner(iel) = zc(in2  ) * wmn2

! ---- Test
#ifdef DEBUG
    somm = fuel1(iel) + fuel2(iel) + fuel3(iel)                   &
         + prod1(iel) + prod2(iel)                                &
         + oxyd(iel)  + xiner(iel)
    if ( (abs(somm-1.d0).gt.epsicp) .or.                          &
         (prod2(iel).lt.(-epzero) .or.                            &
          prod2(iel).gt.(1.d0+epzero)) ) then
      write(NFECRA,*) 'PB CELLULE ',IEL
      write(NFECRA,*) 'INDPDF3',INDPDF(IEL)
      write(NFECRA,*) 'ACHX1F1',ACHX1F1
      write(NFECRA,*) 'ACHX2F2',ACHX2F2
      write(NFECRA,*) 'ACOF1',ACOF1
      write(NFECRA,*) 'ACOF2',ACOF2
      write(NFECRA,*) 'AO2F4',AO2F4
      write(NFECRA,*) 'AO2F3',AO2F3
      write(NFECRA,*) 'AN2F4',AN2F4
      write(NFECRA,*) 'F1M',F1M(IEL)
      write(NFECRA,*) 'F2M',F2M(IEL)
      write(NFECRA,*) 'F3M',F3M(IEL)
      write(NFECRA,*) 'F4M',F4M(IEL)
      sommf = f1m(iel) + f2m(iel) + f3m(iel) + f4m(iel)
      write(NFECRA,*) 'SOMME F', SOMMF
      write(NFECRA,*) 'SI7  ',SI7(IEL)
      write(NFECRA,*) 'SI8  ',SI8(IEL)
      write(NFECRA,*) 'SDEB  ',SDEB(IEL)
      write(NFECRA,*) 'SFIN  ',SFIN(IEL)
      write(NFECRA,*) 'HAUT  ',HAUT(IEL)
      write(NFECRA,*) 'DSI7  ',DSI7(IEL)
      write(NFECRA,*) 'DSI8  ',DSI8(IEL)
      anu1 = 0.25d0*(cx1m-cx2m)
      s1   = - (anu1*ac0(io2)-ac0(ichx1))                         &
           / (anu1*bc0(io2)-bc0(ichx1))
      write(NFECRA,*) 'S1 ', S1
      anu2 = 0.25d0*(2.d0 + cx2m)
      s2   = - (anu2*ac1(io2)-ac1(ichx2))                         &
           / (anu2*bc1(io2)-bc1(ichx2))
      write(NFECRA,*) 'S2 ', S2
      anu3 = 0.5d0
      s3   = - (anu3*ac2(io2)-ac2(ico))                           &
           / (anu3*bc2(io2)-bc2(ico))
      write(NFECRA,*) 'S3 ', S3
      mom0 = haut(iel)*(sfin(iel)-sdeb(iel))                      &
           + dsi8(iel)+dsi7(iel)
      write(NFECRA,*) 'MOM0 ', MOM0, 1.D0
      mom1 = haut(iel)*(sfin(iel)**2-sdeb(iel)**2)/2.d0           &
           + dsi7(iel)*si7(iel)+dsi8(iel)*si8(iel)
      write(NFECRA,*) 'MOM1 ', MOM1, 0.D0
      mom2 = haut(iel)*(sfin(iel)**3-sdeb(iel)**3)/3.d0           &
           + dsi7(iel)*si7(iel)**2+dsi8(iel)*si8(iel)**2
      write(NFECRA,*) 'MOM2 ', MOM2, SP2M(IEL)
      write(NFECRA,*) 'Fuel1',FUEL1(IEL)
      write(NFECRA,*) 'Fuel2',FUEL2(IEL)
      write(NFECRA,*) 'Fuel3',FUEL3(IEL)
      write(NFECRA,*) 'Oxyd',OXYD(IEL)
      write(NFECRA,*) 'Iner',XINER(IEL)
      write(NFECRA,*) 'Prod1',PROD1(IEL)
      write(NFECRA,*) 'Prod2',PROD2(IEL)
      write(NFECRA,*) 'SOMME', SOMM
    endif
#endif
! ---- Fin Test

  endif

enddo



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
n11 = 0
n12 = 0
n13 = 0
n14 = 0
n15 = 0
n16 = 0
n17 = 0
n18 = 0
n19 = 0
n20 = 0
n21 = 0
n22 = 0
n23 = 0
n24 = 0
n25 = 0
n26 = 0
n27 = 0
n28 = 0
n29 = 0
n30 = 0
n31 = 0
n32 = 0

! --> Controle des differentes valeurs des fractions massiques

do iel = 1, ncel
  somm = fuel1(iel) + fuel2(iel) + fuel3(iel) + oxyd(iel)         &
       + prod1(iel) + prod2(iel) + xiner(iel)
  if ( fuel1(iel).lt.(-epzero) .or.                               &
       fuel1(iel).gt.(1.d0+epzero) ) n14 = n14+1
  if ( fuel2(iel).lt.(-epzero) .or.                               &
       fuel2(iel).gt.(1.d0+epzero) ) n15 = n15+1
  if ( fuel3(iel).lt.(-epzero) .or.                               &
       fuel3(iel).gt.(1.d0+epzero) ) n16 = n16+1
  if ( oxyd(iel).lt.(-epzero) .or.                                &
       oxyd(iel).gt.(1.d0+epzero)  ) n17 = n17+1
  if ( xiner(iel).lt.(-epzero) .or.                               &
       xiner(iel).gt.(1.d0+epzero) ) n18 = n18+1
  if ( prod1(iel).lt.(-epzero) .or.                               &
       prod1(iel).gt.(1.d0+epzero) ) n19 = n19+1
  if ( prod2(iel).lt.(-epzero) .or.                               &
       prod2(iel).gt.(1.d0+epzero) ) n20 = n20+1

  if ( abs(somm-1.d0).lt.epsicp )    n3  = n3 +1
enddo

! --> Controle des parametres de la pdf

do iel = 1, ncel
  if ( indpdf(iel).ne.0 ) then
    vars = sp2m(iel) / (-si7(iel)*si8(iel))
    if ( vars.gt.epzero .and.                                     &
      vars.le.(1.-epzero) ) n2 = n2+1
    if ( dsi7(iel).lt.epzero .and. dsi8(iel).lt.epzero )          &
      n4 = n4+1
    if ( dsi7(iel).gt.epzero .and. dsi8(iel).lt.epzero )          &
      n5 = n5+1
    if ( dsi7(iel).lt.epzero .and. dsi8(iel).gt.epzero )          &
      n6 = n6+1
    if ( dsi7(iel).gt.epzero .and. dsi8(iel).gt.epzero )          &
      n7 = n7+1
    if ( dsi7(iel).lt.(-epzero) )   n8 = n8+1
    if ( dsi8(iel).lt.(-epzero) )   n9 = n9+1
    if ( haut(iel).lt.(-epzero) )  n10 = n10+1
    n27 = n27+1
    if ( sdeb(iel).lt.si7(iel) .or. sdeb(iel).gt.si8(iel) )       &
      n11 = n11+1
    if ( sfin(iel).lt.si7(iel) .or. sfin(iel).gt.si8(iel) )       &
      n12 = n12+1
    if ( (sfin(iel)-sdeb(iel)).lt.zero ) n13 = n13+1
    if ( indpdf(iel).eq.12 ) n28 = n28+1
    if ( indpdf(iel).eq.2  ) n29 = n29+1
    if ( indpdf(iel).eq.13 ) n31 = n31+1
    if ( indpdf(iel).eq.11 ) n30 = n30+1
    if ( indpdf(iel).eq.3  ) n32 = n32+1
    if ( ipi7i8(iel).gt.10 .and. ipi7i8(iel).lt.15 ) n21 = n21+1
    if ( ipi7i8(iel).gt.20 .and. ipi7i8(iel).lt.25 ) n22 = n22+1
    if ( ipi7i8(iel).eq.11 .or.  ipi7i8(iel).eq.21 ) n23 = n23+1
    if ( ipi7i8(iel).eq.12 .or.  ipi7i8(iel).eq.22 ) n24 = n24+1
    if ( ipi7i8(iel).eq.13 .or.  ipi7i8(iel).eq.23 ) n25 = n25+1
    if ( ipi7i8(iel).eq.14 .or.  ipi7i8(iel).eq.24 ) n26 = n26+1
  endif
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
  inttmp(11) = n11
  inttmp(12) = n12
  inttmp(13) = n13
  inttmp(14) = n14
  inttmp(15) = n15
  inttmp(16) = n16
  inttmp(17) = n17
  inttmp(18) = n18
  inttmp(19) = n19
  inttmp(20) = n20
  inttmp(21) = n21
  inttmp(22) = n22
  inttmp(23) = n23
  inttmp(24) = n24
  inttmp(25) = n25
  inttmp(26) = n26
  inttmp(27) = n27
  inttmp(28) = n28
  inttmp(29) = n29
  inttmp(30) = n30
  inttmp(31) = n31
  inttmp(32) = n32
  call parism(nbrint,inttmp)
  !==========
endif

write(nfecra,1000) n1 , n27
write(nfecra,1100) n28, n29, n31, n30, n32, n21, n22, n23,        &
                   n24, n25, n26
write(nfecra,2000) n4 , n5 , n6 , n7
write(nfecra,2100) n27, n2 , n8 , n9 , n10, n11, n12, n13
write(nfecra,2200) n1 , n3 , n14, n15, n16, n17, n18, n19, n20

!-------
! FORMAT
!-------

 1000 format (/,                                                  &
'MODELISATION DE LA COMBUSTION AVEC LE MODELE DE DIFFUSION ',     &
'TURBULENTE (CPLYM1)',/,                                    &
'CHIMIE RAPIDE A 3 CORPS - EXTENSION A 3 COMBUSTIBLES ',          &
'(Application au CP)',/,                                    &
'==========================================================',     &
'==================',/,                                     &
' Nb de points de calculs                                     = ',&
   i9,/,                                                    &
' Nb de points turbulents (passage par les PDF)               = ',&
   i9)

 1100 format(                                                           &
! ..v.7..1....v    ....2....v....3....v....4....v....5....v....6....v....7.I
' Nb de points turbulents pour lesquels I7 app. [I4,L3] T12   = ',&
   i9,/,                                                    &
' - - - - - - - - - - - - pour lesquels I7 app. [I4,L5] T2    = ',&
   i9,/,                                                    &
' - - - - - - - - - - - - pour lesquels I7 app. [I4,L5] T13   = ',&
   i9,/,                                                    &
' - - - - - - - - - - - - pour lesquels I7 app. [L5,I3max] T11= ',&
   i9,/,                                                    &
' - - - - - - - - - - - - pour lesquels I7 = I4  T3           = ',&
   i9,/,                                                    &
' - - - - - - - - - - - - pour lesquels I7 domaine P     = ',     &
   i9,/,                                                    &
' - - - - - - - - - - - - pour lesquels I7 domaine I     = ',     &
   i9,/,                                                    &
' - - - - - - - - - - - - pour lesquels I8 domaine P     = ',     &
   i9,/,                                                    &
' - - - - - - - - - - - - pour lesquels I8 domaine I     = ',     &
   i9,/,                                                    &
' - - - - - - - - - - - - pour lesquels I8 domaine R     = ',     &
   i9,/,                                                    &
' - - - - - - - - - - - - pour lesquels I8 domaine HR    = ',     &
   i9)

 2000 format(                                                           &
'PDF CONJOINTE DEGENEREE EN PDF MONODIMENSIONNELLE',/,      &
' Nb de points PDF rectangle sans Dirac            = ',I9,/,&
' - - - - - - - - - - - - -  et Dirac en I7        = ',I9,/,&
' - - - - - - - - - - - - -  et Dirac en I8        = ',I9,/,&
' - - - - - - - - - - - - -  et Diracs en I8 et I7 = ',I9)

 2100 format(                                                           &
'CONTROLE DES PARAMETRES DE LA PDF',/,                      &
' Nb de points turbulents (passage par les PDF)              = ', &
   i9,/,                                                    &
' Nb de points turbulents pour lesquels VARS est physique    = ', &
   i9,/                                                           &
' Nb de points Dirac_I7, Dirac_I8, HAUT < 0  = ',I9,I9,I9,/,&
' Nb de points Abs.Deb < SI7 ou > SI8        = ', I9,/,     &
' Nb de points Abs.Fin < SI7 ou > SI8        = ', I9,/,     &
' Nb de points (Abs.Fin-Abs.Deb) < 0         = ', I9)

 2200 format(                                                           &
'CONTROLE DES VALEURS DES FRACTIONS MASSIQUES',/,           &
' Nb de points de calculs                                    = ', &
   i9,/,                                                    &
' Nb de points de calculs qui respectent somme des Yi = 1    = ', &
   i9,/,                                                    &
' Nb de points YCHx1m,YCHx2m,YCO < 0 ou > 1 = ',I9,I9,I9,/, &
' Nb de points YO2,YN2        < 0 ou > 1    = ',I9,I9,/,    &
' Nb de points YCO2,YH2O      < 0 ou > 1    = ',I9,I9)


return
end subroutine
