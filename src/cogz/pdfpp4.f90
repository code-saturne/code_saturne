!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

subroutine pdfpp4 &
!================

 ( ncelet , ncel  ,                                               &
   fm     , fp2m  , yfm    , yfp2m , coyfp )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DES PARAMETRES DE LA PDF
! PDF LIBBY - WILLIAMS 4  POINTS AVEC HYPOTHESE DE CURL
!             OU CURL MODIFIE


! COMMENTAIRES :
! ------------

!    Dans un diagramme (F, Yf), on construit deux droites:
!         - La droite de combustion complete
!         - La droite de melange

!    Dans ce domaine, nous allons construire  deux pics sur F qui
!    seront dedoubler chaque un en deux avec un Curl sur Yf


! LE RESULTAT EST :
! ---------------

!    CALCUL DES PARAMETRES ASSOCIES AUX FONCTIONS DIRAC

!      Les positions des Pics sont :
!         [F(1),Y1(1)] et [F(1),Y1(2)]
!         [F(2),Y2(1)] et [F(2),Y2(2)]
!      Leurs amplitudes respectives sont :
!               D1(1) et D1(2)
!               D2(1) et D2(2)
!      Pour chaque dirac, on calcule :
!      la temperature Ti(j),
!      la masse volumique RHOi(j),
!      le terme source chimique Wi(j),
!           i etant la positiion sur F du pic de Dirac
!           j etant la positiion sur Yf du pic de Dirac

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! fm               ! tr ! <-- ! moyenne de la fraction de melange              !
! fp2m             ! tr ! <-- ! variance de la fraction de melange             !
! yfm              ! tr ! <-- ! moyenne de la fraction massique                !
! yfp2m            ! tr ! <-- ! variance de la fraction massique               !
! coyfp            ! tr ! <-- ! covariance                                     !
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
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use pointe
use ppppar
use ppthch
use ppincl
use coincl
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments
!===============================================================================

integer          ncelet, ncel

double precision fm(ncelet)   , fp2m(ncelet)
double precision yfm(ncelet)  , yfp2m(ncelet)
double precision coyfp(ncelet)

!===============================================================================
! Local variables
!===============================================================================

integer          iel, igg, idirac

!      INTEGER          IPCKAB, IPT4 , IPT3

double precision f(ndracm), y(ndracm), d(ndracm)
double precision h(ndracm), teml(ndracm)
double precision maml(ndracm), w(ndracm)
double precision rhol(ndracm), theta(ndracm)

double precision coefg(ngazgm), epsi
double precision yfuel
double precision yoxyd, yo2
double precision yprod,fmp,fp2mp,yfmp,yfp2mp,coyfpp
double precision yfp2max

double precision nbmol,  temsmm
double precision sum7, sum8, sum9, sum10, sum11, sum12, sum17
double precision sum1, sum2, sum3, sum4, sum5, sum6, sum16, sum15
double precision wmin, wmax, tetmin, tetmax, o2min, o2max

double precision y1, y2, f1, f2
double precision cstfa1, cstfa2, cst
double precision cstvar
!     DOUBLE PRECISION CSTVA2

double precision ymin(2), ymax(2)
double precision y2p(2)
double precision y2pmin(2), y2pmax(2)

! ---> variables pour clipping

integer          icpt1,icpt2
integer          cly2p1, cly2p2, cly2p3,cly2p4
integer          clfp2m, clyf21, clyf22
integer          clcyf1, clcyf2

double precision vymx, vfmx
!      DOUBLE PRECISION VYMN
double precision vcyfmx, vcyfmn
double precision mxcyfp,mxcfp, mncyfp
double precision mxccyf, mnccyf
double precision mxyfp2,maxfp2, mnyfp2
double precision mxcoyf, mncoyf
double precision mcy2p1, mcy2p3
double precision mcy2p2, mcy2p4
double precision my2p1 , my2p3, my2p2, my2p4
double precision, dimension(:), pointer ::  crom
double precision, dimension(:), pointer ::  cpro_temp, cpro_tsc, cpro_mam
double precision, dimension(:), pointer ::  cpro_ym1, cpro_ym2, cpro_ym3
type(pmapper_double_r1), dimension(:), pointer :: cpro_fmel, cpro_fmal, cpro_teml
type(pmapper_double_r1), dimension(:), pointer :: cpro_tscl, cpro_rhol, cpro_maml
type(pmapper_double_r1), dimension(:), pointer :: cpro_ampl

integer ipass
data    ipass /0/
save    ipass

!===============================================================================

!===============================================================================
!     0. POSITION ET INITIALISATION DES VARIABLES
!===============================================================================

call field_get_val_s(icrom, crom)
call field_get_val_s(itemp, cpro_temp)
call field_get_val_s(iym(1), cpro_ym1)
call field_get_val_s(iym(2), cpro_ym2)
call field_get_val_s(iym(3), cpro_ym3)
call field_get_val_s(itsc, cpro_tsc)
call field_get_val_s(imam, cpro_mam)

allocate(cpro_fmel(ndirac))
allocate(cpro_fmal(ndirac))
allocate(cpro_teml(ndirac))
allocate(cpro_tscl(ndirac))
allocate(cpro_rhol(ndirac))
allocate(cpro_maml(ndirac))
allocate(cpro_ampl(ndirac))

do idirac = 1, ndirac
  call field_get_val_s(iampl(idirac), cpro_ampl(idirac)%p)
  call field_get_val_s(ifmel(idirac), cpro_fmel(idirac)%p)
  call field_get_val_s(ifmal(idirac), cpro_fmal(idirac)%p)
  call field_get_val_s(iteml(idirac), cpro_teml(idirac)%p)
  call field_get_val_s(imaml(idirac), cpro_maml(idirac)%p)
  call field_get_val_s(irhol(idirac), cpro_rhol(idirac)%p)
  call field_get_val_s(itscl(idirac), cpro_tscl(idirac)%p)
enddo

!if ( iirayo.gt.0 ) then
!  call field_get_val_s(ickabs, cpro_ckabs)
!  call field_get_val_s(it4m, cpro_t4m)
!  call field_get_val_s(it3m, cpro_t3m)
!endif

! ---> Initialisation des variables et compteurs

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

epsi = 1.d-10

clfp2m  = 0
clyf21 = 0
clyf22 = 0
clcyf1 = 0
clcyf2 = 0

cly2p1 = 0
cly2p2 = 0
cly2p3 = 0
cly2p4 = 0

wmax = -1.d+10
wmin =  1.d+10
tetmin = 1.d+10
tetmax = -1.d+10
o2max = -1.d+10
o2min =  1.d+10

mxcfp = -1.d+10
mxcyfp = -1.d+10
mncyfp = -1.d+10
mxccyf = -1.d+10
mnccyf = -1.d+10

maxfp2= -1.d+10
mxyfp2= -1.d+10
mnyfp2= -1.d+10
mxcoyf= -1.d+10
mncoyf= -1.d+10

mcy2p1=-1.d+10
mcy2p3=-1.d+10
mcy2p2=-1.d+10
mcy2p4=-1.d+10

my2p1=-1.d+10
my2p3=-1.d+10
my2p2=-1.d+10
my2p4=-1.d+10

icpt1 = 0
icpt2 = 0

!===============================================================================
!     1. BOUCLE SUR LES CELLULES
!===============================================================================

ipass = ipass + 1

do iel = 1, ncel

! ---> position des variables

fmp=fm(iel)
yfmp=max(yfm(iel), zero)
fp2mp=fp2m(iel)
yfp2mp=yfp2m(iel)
coyfpp=coyfp(iel)

!===============================================================================
!    TEST DE PASSAGE OU NON PAR LA PDF
!    (non passage par pdf si aucune variance ni en f ni en y)
!    (on peut si on veut toujours passer par la pdf,
!     le cas des variances nulles y est traite)
!===============================================================================

 if (((fp2mp.lt.epsi).and.(yfp2mp.lt.epsi))) then

!===============================================================================
!    2.    NON PASSAGE PAR PDF
!===============================================================================

    icpt1 = icpt1 +1

    sum1 = zero
    sum2 = zero
    sum3 = zero
    sum4 = zero
    sum5 = zero
    sum6 = zero
    sum15 = zero
    sum16 = zero

! ---> boucle sur chaque DIRAC

    do idirac =1, ndirac

! ---> calcul de f, Y, Amplitude de chaque DIRAC == Val moy

      d(idirac) = 1.d0 / ndirac
      f(idirac) = fmp
      y(idirac) = yfmp

! ---> Calcul de l'enthalpie

      h(idirac) = ((hmax-hmin)*f(idirac) + hmin*fmax - hmax*fmin) &
                / (fmax-fmin)

! ---> Calcul de la fraction massique des gaz (F, O et P)

      yfuel = y(idirac)
      yoxyd = 1.d0 - (coeff3+1.0d0)*f(idirac) + coeff3*y(idirac)
      yprod = 1.d0 - yfuel - yoxyd
      yo2   = coeff1 - (coeff1 + coeff2) * f(idirac)              &
                      + coeff2 * y(idirac)

! ---> Calcul de la masse molaire et de la temperature

      coefg(1) = yfuel
      coefg(2) = yoxyd
      coefg(3) = yprod

! ---> Masse molaire

      nbmol = zero
      do igg = 1, ngazg
        nbmol = nbmol + coefg(igg)/wmolg(igg)
      enddo
      maml(idirac) = 1.d0/nbmol

! ---> Calcul de la temperature pour chaque pic

      teml(idirac) = cs_gas_combustion_h_to_t(coefg, h(idirac))

! ---> Calcul de la masse volumique pour chaque pic

      if ( ipass.gt.1 .or.                                        &
          (isuite.eq.1.and.initro.eq.1) ) then
        rhol(idirac) = p0 * maml(idirac)                   &
             / (cs_physical_constants_r*teml(idirac))
      else
        rhol(idirac) = ro0
      endif

! ---> Calcul du terme source du scalaire YFM pour chaque pic

      theta(idirac) = ta / teml(idirac)                           &
           * (1.d0 - teml(idirac) / tstar)

      w(idirac) = vref / lref * (- d(idirac)                      &
           * yfuel*yo2                                            &
           * exp( -theta(idirac) ))

! ---> Controle du signe de W

      w(idirac) = min( w(idirac), zero)

! ---> Masse molaire du melange

      sum1 = sum1 + d(idirac)*maml(idirac)

! ---> Temperature du melange

      sum2 = sum2 + d(idirac)*teml(idirac)

! ---> Temperature / Masse molaire

      sum3 = sum3 + d(idirac)*teml(idirac)/maml(idirac)

! ---> Fractions massiques des especes globales

      sum4 = sum4 + yfuel*d(idirac)
      sum5 = sum5 + yoxyd*d(idirac)
      sum6 = sum6 + yprod*d(idirac)
      sum15 = sum15 +rhol(idirac)*d(idirac)
      sum16 = sum16 +w(idirac)

! ---> Stockage des proprietes

      cpro_ampl(idirac)%p(iel) = d(idirac)
      cpro_fmel(idirac)%p(iel) = f(idirac)
      cpro_fmal(idirac)%p(iel) = y(idirac)
      cpro_teml(idirac)%p(iel) = teml(idirac)
      cpro_maml(idirac)%p(iel) = maml(idirac)
      cpro_rhol(idirac)%p(iel) = rhol(idirac)
      cpro_tscl(idirac)%p(iel) = w(idirac)

!fin de boucle sur chaque DIRAC
    enddo

    cpro_mam(iel)  = sum1
    cpro_temp(iel) = sum2
    temsmm         = sum3
    cpro_ym1(iel)  = sum4
    cpro_ym2(iel)  = sum5
    cpro_ym3(iel)  = sum6
    cpro_tsc(iel)  = sum16

! ---> Masse volumique du melange

    if ( ipass.gt.1.or.                                           &
        (isuite.eq.1.and.initro.eq.1) ) then
      crom(iel) = srrom*crom(iel)               &
           +(1.d0-srrom)*(p0/(cs_physical_constants_r*temsmm))
    endif

  else

!==================================================================================
!    3.    PASSAGE PAR PDF
!==================================================================================

    icpt2 = icpt2 + 1

! ---> Clipping sur la variance en f

    vfmx = (fmax-fmp)*(fmp-fmin)

    if (fp2mp.gt.vfmx) then
      if ((fp2mp-vfmx).gt.mxcfp) then
        mxcfp=(fp2mp-vfmx)
        maxfp2=vfmx
      endif
      fp2mp=vfmx
      clfp2m=clfp2m+1
    endif

! ---> Calcul des positions et amplitudes en F des Pics avec lWCURL en F

    cst = 1.d0

    call lwcurl                                                   &
!         =========
   ( cst  , fmp , fp2mp  ,                                        &
       fmin , fmax ,                                              &
       f1 , f2 , cstfa1 , cstfa2 )

    f(1) = f1
    f(2) = f1
    f(3) = f2
    f(4) = f2

! ---> Determination des valeur max et min de Yf en F1 et F2

    ymin(1) = max(zero , ((f1- fs(1))/(1d0-fs(1))))
    ymax(1) = f1
    ymin(2) = max(zero , ((f2- fs(1))/(1d0-fs(1))))
    ymax(2) = f2

! ---> clipping COVARIANCE en fonction des bornes des moyennes conditionnelles

    vcyfmx=min(ymax(2)*cstfa2*(f2-f1)-yfmp*(fmp-f1)               &
         ,yfmp*(f2-fmp)-ymin(1)*cstfa1*(f2-f1))

    vcyfmn=max(ymin(2)*cstfa2*(f2-f1)-yfmp*(fmp-f1)               &
         ,yfmp*(f2-fmp)-ymax(1)*cstfa1*(f2-f1))

    if (coyfpp.gt.vcyfmx) then
      if ((coyfpp-vcyfmx).gt.mxccyf) then
        mxccyf=(coyfpp-vcyfmx)
        mxcoyf=vcyfmx
      endif
      coyfpp=vcyfmx
      clcyf1=clcyf1+1
    elseif (coyfpp.lt.vcyfmn) then
      if ((vcyfmn-coyfpp).gt.mnccyf) then
        mnccyf=(vcyfmn-coyfpp)
        mncoyf=vcyfmn
      endif
      coyfpp=vcyfmn
      clcyf2=clcyf2+1
    endif

! --->  On calcul les Moyennes conditionnelles Y1, Y2

    if ((f2-f1).gt.epsi) then
      y2 = (yfmp*(fmp- f1) + coyfpp)                              &
           /(cstfa2*(f2 - f1))
      y1 = (yfmp*(f2-fmp) - coyfpp)                               &
           /(cstfa1*(f2 - f1))
    else
      y2 = yfmp
      y1 = yfmp
    endif
    if ((fmp-yfmp).lt.epsi) then
      y2=f2
      y1=f1
    endif
    if ((yfmp-max(zero,(fmp -fs(1))                               &
         /(1.d0 - fs(1)))).lt.epsi) then
      y2=max(zero,(f2 -fs(1))/(1.d0 - fs(1)))
      y1=max(zero,(f1 -fs(1))/(1.d0 - fs(1)))
    endif

! ---> Determination des valeurs max et min des variances de y en F1 et F2

    y2pmax(1) =(y1-ymin(1))*(ymax(1) - y1)
    y2pmin(1) = zero
    y2pmax(2) = (y2-ymin(2))*(ymax(2) - y2)
    y2pmin(2) = zero

! ---> calcul de la variance en Y max avec Y2PMAX 1 et 2

    yfp2max=cstfa1*(y1**2+y2pmax(1))                              &
         +cstfa2*(y2**2+y2pmax(2))-yfmp**2
    vymx=yfp2max

! ---> rapport des variances conditionnelles

    cstvar = 0.d0

    if (((ymax(2)-y2).gt.epsi).and.((y2-ymin(2)).gt.epsi)         &
         .and.((ymax(1)-y1).gt.epsi).and.((y1-ymin(1)).gt.epsi))  &
         then
      cstvar = ((ymax(2)-y2)*(y2-ymin(2)))                        &
           /((ymax(1)-y1)*(y1 -ymin(1)))
    endif


! ---> clipping VARIANCE Y
!      (on peut soit clipper la variance en y en fonction
!       des valeurs extremes des variances conditionnelles
!       ou clipper directement les variances conditionnelles)

!         CSTVA2 =  (CSTFA1*Y1**2+CSTFA2*Y2**2)-YFMP**2

!         IF (((YMAX(2)-Y2).GT.EPSI).AND.((Y2-YMIN(2)).GT.EPSI)
!     &  .AND.((YMAX(1)-Y1).GT.EPSI).AND.((Y1-YMIN(1)).GT.EPSI))
!     &       THEN

!      VYMX=MIN(Y2PMAX(1)*(CSTFA1+CSTFA2*CSTVAR)+CSTVA2
!     &         ,Y2PMAX(2)*(CSTFA1/CSTVAR+CSTFA2)+CSTVA2)
!       VYMN=MAX(Y2PMIN(1)*(CSTFA1+CSTFA2*CSTVAR)+CSTVA2
!     &           ,Y2PMIN(2)*(CSTFA1/CSTVAR+CSTFA2)+CSTVA2)

!         ELSEIF (((YMAX(2)-Y2).GT.EPSI)
!     &         .AND.((Y2-YMIN(2)).GT.EPSI)) THEN
!  VYMX=Y2PMAX(2)*CSTFA2+CSTVA2
!  VYMN=Y2PMIN(2)*CSTFA2+CSTVA2

!         ELSEIF (((YMAX(1)-Y1).GT.EPSI)
!     &         .AND.((Y1-YMIN(1)).GT.EPSI)) THEN
!  VYMX=Y2PMAX(1)*CSTFA1+CSTVA2
!  VYMN=Y2PMIN(1)*CSTFA1+CSTVA2
!  ENDIF

    if (yfp2mp.gt.vymx) then
      if ((yfp2mp-vymx).gt.mxcyfp) then
        mxcyfp=(yfp2mp-vymx)
        mxyfp2=vymx
      endif
      yfp2mp=vymx
      clyf21=clyf21+1
!   ELSEIF (YFP2M(IEL).LT.VYMN) THEN
!   IF ((VYMN-YFP2MP).GT.MNCYFP) THEN
!      MNCYFP=(VYMN-YFP2MP)
!      MNYFP2=VYMN
!   ENDIF
!   YFP2MP=VYMN
!   CLYF22=CLYF22+1
    endif

! ---> calcul des variances conditionnelles

    if (((ymax(2)-y2).gt.epsi).and.((y2-ymin(2)).gt.epsi)         &
         .and.((ymax(1)-y1).gt.epsi).and.((y1-ymin(1)).gt.epsi))  &
         then

      y2p(1) = ((yfmp**2)+yfp2mp                                  &
           -cstfa2*(y2**2)-cstfa1*(y1**2))                        &
           /(cstfa1 + cstfa2*cstvar)

      y2p(2) = ((yfmp**2)+yfp2mp                                  &
           -cstfa2*(y2**2)-cstfa1*(y1**2))                        &
           /(cstfa1/cstvar + cstfa2)

    elseif (((ymax(2)-y2).gt.epsi)                                &
           .and.((y2-ymin(2)).gt.epsi)) then

      y2p(1) = zero

      y2p(2) = ((yfmp**2)+yfp2mp                                  &
           -cstfa2*(y2**2)-cstfa1*(y1**2))                        &
           /cstfa2

    elseif (((ymax(1)-y1).gt.epsi)                                &
           .and.((y1-ymin(1)).gt.epsi)) then

      y2p(2) = zero

      y2p(1) = ((yfmp**2)+yfp2mp                                  &
           -cstfa2*(y2**2)-cstfa1*(y1**2))                        &
           /cstfa1

    else
      y2p(1) = zero
      y2p(2) = zero
    endif

! ---> clipping pour variances conditionnelles

    if (y2p(1).gt.y2pmax(1)) then
      if ((y2p(1)-y2pmax(1)).gt.mcy2p1) then
        mcy2p1=(y2p(1)-y2pmax(1))
        my2p1=y2pmax(1)
      endif
      y2p(1) = y2pmax(1)
      y2p(2) = (((yfmp**2)+yfp2mp-cstfa1*                         &
           ((y1**2)+y2p(1)))/cstfa2)-(y2**2)
      cly2p1 =  cly2p1 + 1
    elseif (y2p(1).lt.y2pmin(1)) then
      if ((y2pmin(1)-y2p(1)).gt.mcy2p3) then
        mcy2p3=(y2pmin(1)-y2p(1))
        my2p3=y2pmin(1)
      endif
      y2p(1) = y2pmin(1)
      y2p(2) = (((yfmp**2)+yfp2mp-cstfa1*                         &
           ((y1**2)+y2p(1)))/cstfa2)-(y2**2)
      cly2p3 = cly2p3 + 1
    endif
    if (y2p(2).gt.y2pmax(2)) then
      if ((y2p(2)-y2pmax(2)).gt.mcy2p2) then
        mcy2p2=(y2p(2)-y2pmax(2))
        my2p2=y2pmax(2)
      endif
      y2p(2) = y2pmax(2)
      y2p(1) = (((yfmp**2)+yfp2mp-cstfa2*                         &
           ((y2**2)+y2p(2)))/cstfa1)-(y1**2)
      cly2p2 =  cly2p2 + 1
    elseif (y2p(2).lt.y2pmin(2)) then
      if ((y2pmin(2)-y2p(2)).gt.mcy2p4) then
        mcy2p4=(y2pmin(2)-y2p(2))
        my2p4=y2pmin(2)
      endif
      y2p(2) = y2pmin(2)
      y2p(1) = (((yfmp**2)+yfp2mp-cstfa2*                         &
           ((y2**2)+y2p(2)))/cstfa1)-(y1**2)
      cly2p4 = cly2p4 + 1
    endif

! ---> calcul des positions et amplitudes des pics en Y sur F1

    call lwcurl                                                   &
!         =========
  ( cstfa1  , y1   , y2p(1)  ,                                    &
    ymin(1) , ymax(1) ,                                           &
    y(1) , y(2) , d(1) , d(2) )

! ---> calcul des positions et amplitudes des pics en Y sur F2

    call lwcurl                                                   &
!         =========
  ( cstfa2  , y2   , y2p(2)  ,                                    &
    ymin(2) , ymax(2) ,                                           &
    y(3) , y(4) , d(3) , d(4) )


!===============================================================================
! 2.  DETERMINATION DES GRANDEURS THERMOCHIMIQUES DES DEUX PICS
!===============================================================================

! ---> Calcul de l'enthalpies en 1 et 2

  sum7  = zero
  sum8  = zero
  sum9  = zero
  sum10 = zero
  sum11 = zero
  sum12 = zero
  sum17 = zero

! ---> boucle sur chaque DIRAC

  do idirac = 1, ndirac

! ---> Calcul de l'enthalpie

    h(idirac) = ( (hmax-hmin)*f(idirac)                           &
         + hmin*fmax - hmax*fmin) / (fmax-fmin)

! ---> Calcul de la fraction massique des gaz (F, O et P) pour chaque pic

    yfuel = y(idirac)
    yoxyd = 1.d0 - (coeff3+1.0d0)*f(idirac)                       &
         + coeff3*y(idirac)
    yprod = 1.d0 - yfuel - yoxyd
    yo2   = coeff1 - (coeff1 + coeff2) * f(idirac)                &
         + coeff2 * y(idirac)

! ---> Coefficients d'absorption pour chaque pic

!            IF ( IIRAYO .GT. 0  ) THEN
!              KABSGF = YFUEGF(IEL)*KABSG(1) + YOXYGF(IEL)*KABSG(2)
!     &                +YPROGF(IEL)*KABSG(3)
!              KABSGB = YFUEGB(IEL)*KABSG(1) + YOXYGB(IEL)*KABSG(2)
!     &                +YPROGB(IEL)*KABSG(3)
!            ENDIF


! ---> Calcul de la masse molaire et de la temperature pour chaque pic

    coefg(1) = yfuel
    coefg(2) = yoxyd
    coefg(3) = yprod

! ---> Masse molaire pour chaque pic

    nbmol = zero
    do igg = 1, ngazg
      nbmol = nbmol + coefg(igg)/wmolg(igg)
    enddo
    maml(idirac) = 1.d0/nbmol

! --->Calcul de la temperature pour chaque pic

    teml(idirac) = cs_gas_combustion_h_to_t(coefg, h(idirac))

! ---> Calcul de la masse volumique  pour chaque pic

    if ( ipass.gt.1.or.                                           &
        (isuite.eq.1.and.initro.eq.1)) then
      rhol(idirac) = p0 * maml(idirac)                     &
           /(cs_physical_constants_r*teml(idirac))
    else
      rhol(idirac) = ro0
    endif

! ---> Calcul du terme source du scalaire YFM  pour chaque pic

    theta(idirac) = ta / teml(idirac)                             &
         *(1.d0 - teml(idirac)/tstar)

    tetmax = max(theta(idirac),tetmax)
    tetmin = min(theta(idirac),tetmin)

    w(idirac) = vref / lref                                       &
         *(- d(idirac)                                            &
         * yfuel*yo2                                              &
         * exp( -theta(idirac) ))

    wmax = max(w(idirac),wmax)
    wmin = min(w(idirac),wmin)
    o2max = max(yo2,o2max)
    o2min = min(yo2,o2min)

! ---> Controle du signe de W

    w(idirac) = min( w(idirac), zero)

! ---> Masse molaire du melange

    sum7 = sum7 + d(idirac)*maml(idirac)

! ---> Temperature du melange

    sum8 = sum8 + d(idirac)*teml(idirac)

! ---> Temperature / Masse molaire

    sum9 = sum9 + d(idirac)*teml(idirac)                          &
         /maml(idirac)

! ---> Fractions massiques des especes globales

    sum10 = sum10 + yfuel*d(idirac)

    sum11 = sum11 + yoxyd*d(idirac)

    sum12 = sum12 + yprod*d(idirac)

    sum17 = sum17 + w(idirac)

! ---> Stockage des proprietes

    cpro_ampl(idirac)%p(iel) = d(idirac)
    cpro_fmel(idirac)%p(iel) = f(idirac)
    cpro_fmal(idirac)%p(iel) = y(idirac)
    cpro_maml(idirac)%p(iel) = maml(idirac)
    cpro_teml(idirac)%p(iel) = teml(idirac)
    cpro_rhol(idirac)%p(iel) = rhol(idirac)
    cpro_tscl(idirac)%p(iel) = w(idirac)

! ---> Grandeurs relatives au rayonnement

!           if ( iirayo.gt.0 ) then
!             cpro_ckab(iel) = ygfm*kabsgf + ygbm*kabsgb
!             cpro_t4m(iel)  = ygfm*tgf**4+ygbm*tgb**4
!             cpro_t3m(iel)  = ygfm*tgf**3+ygbm*tgb**3
!           endif

! fin de boucle sur chaque DIRAC
  enddo

  cpro_mam(iel)  = sum7
  cpro_temp(iel) = sum8
  temsmm         = sum9
  cpro_ym1(iel)  = sum10
  cpro_ym2(iel)  = sum11
  cpro_ym3(iel)  = sum12
  cpro_tsc(iel)  = sum17

! ---> Masse volumique du melange

  if ( ipass.gt.1 .or.                                            &
      (isuite.eq.1.and.initro.eq.1) ) then
    crom(iel) = srrom * crom(iel)               &
         + (1.d0-srrom) * (p0/(cs_physical_constants_r*temsmm))
  endif

! de passage ou non par la PDF
endif

! fin de boucle sur les cellules
enddo

deallocate(cpro_fmel, cpro_fmal, cpro_teml)
deallocate(cpro_tscl, cpro_rhol, cpro_maml)
deallocate(cpro_ampl)

! ---> impression clipping

if(irangp.ge.0) then
  call parcpt(clfp2m)
  call parmax(mxcfp)
  call parmax(maxfp2)
endif

WRITE(NFECRA,*)' nombre de clip haut sur la variance en f =',     &
               clfp2m
WRITE(NFECRA,*)' ecart maximum (valeur atteinte-valeur max) =',   &
               mxcfp
WRITE(NFECRA,*)' valeur max (pour l ecart max) =', MAXFP2
WRITE(NFECRA,*)'     '

  if(irangp.ge.0) then
  call parcpt(clyf21)
  call parmax(mxcyfp)
  call parmax(mxyfp2)
  endif

WRITE(NFECRA,*)' nombre de clip haut sur la variance en y =',     &
               clyf21
WRITE(NFECRA,*)' ecart maximum (valeur atteinte-valeur max)  =',  &
               mxcyfp
WRITE(NFECRA,*)' valeur max (pour l ecart max) =', MXYFP2
WRITE(NFECRA,*)'     '

!     WRITE(NFECRA,*)' nombre de clip bas sur la variance en y =',
!    &               CLYF22
!     WRITE(NFECRA,*)' ecart maximum (valeur min-valeur atteinte) =',
!    &               MNCYFP
!     WRITE(NFECRA,*)' valeur min (pour l ecart max) =', MNYFP2
!     WRITE(NFECRA,*)'     '

if(irangp.ge.0) then
  call parcpt(clcyf1)
  call parmax(mxccyf)
  call parmax(mxcoyf)
endif

WRITE(NFECRA,*)' nombre de clip haut sur la covariance =',        &
               clcyf1
WRITE(NFECRA,*)' ecart maximum (valeur atteinte-valeur max)F =',  &
               mxccyf
WRITE(NFECRA,*)' valeur max (pour l ecart max) =', MXCOYF
WRITE(NFECRA,*)'     '

if(irangp.ge.0) then
  call parcpt(clcyf2)
  call parmax(mnccyf)
  call parmin(mncoyf)
endif

WRITE(NFECRA,*)' nombre de clip bas sur la covariance =',         &
               clcyf2
WRITE(NFECRA,*)' ecart maximum (valeur min-valeur atteinte) =',   &
               mnccyf
WRITE(NFECRA,*)' valeur min (pour l ecart max) =', MNCOYF
WRITE(NFECRA,*)'     '

if(irangp.ge.0) then
  call parcpt(cly2p1)
  call parmax(mcy2p1)
  call parmax(my2p1)
endif

write(nfecra,*)                                                   &
' nombre de clip haut sur la variance conditionnelle 1 =', CLY2P1
WRITE(NFECRA,*)' ecart maximum (valeur atteinte-valeur max) =',   &
               mcy2p1
WRITE(NFECRA,*)' valeur max (pour l ecart max) =', MY2P1
WRITE(NFECRA,*)'     '

if(irangp.ge.0) then
  call parcpt(cly2p3)
  call parmax(mcy2p3)
  call parmin(my2p3)
endif

write(nfecra,*)                                                   &
' nombre de clip bas sur la variance conditionnelle 1 =', CLY2P3
WRITE(NFECRA,*)' ecart maximum (valeur min-valeur atteinte) =',   &
               mcy2p3
WRITE(NFECRA,*)' valeur min (pour l ecart max) =', MY2P3
WRITE(NFECRA,*)'     '

if(irangp.ge.0) then
  call parcpt(cly2p2)
  call parmax(mcy2p2)
  call parmax(my2p2)
endif

write(nfecra,*)                                                   &
' nombre de clip haut sur la variance conditionnelle 2 =', CLY2P2
WRITE(NFECRA,*)' ecart maximum (valeur atteinte-valeur max) =',   &
               mcy2p2
WRITE(NFECRA,*)' valeur max (pour l ecart max) =', MY2P2
WRITE(NFECRA,*)'     '

if(irangp.ge.0) then
  call parcpt(cly2p4)
  call parmax(mcy2p4)
  call parmin(my2p4)
endif

write(nfecra,*)                                                   &
' nombre de clip bas sur la variance conditionnelle 2=', CLY2P4
WRITE(NFECRA,*)' ecart maximum (valeur min-valeur atteinte) =',   &
               mcy2p4
WRITE(NFECRA,*)' valeur min (pour l ecart max) =', MY2P4
WRITE(NFECRA,*)'     '

if(irangp.ge.0) then
  call parcpt(icpt1)
  call parcpt(icpt2)
  call parmax(o2max)
  call parmin(o2min)
  call parmax(tetmax)
  call parmin(tetmin)
  call parmax(wmax)
  call parmin(wmin)
endif

WRITE(NFECRA,*) ' POINT NON  PDF = ',ICPT1
WRITE(NFECRA,*) ' POINT AVEC PDF = ',ICPT2
WRITE(NFECRA,*) ' MIN MAX O2     = ',O2MIN,O2MAX
WRITE(NFECRA,*) ' MIN MAX THETA  = ',TETMIN,TETMAX
WRITE(NFECRA,*) ' MIN MAX W      = ',WMIN,WMAX

 end
