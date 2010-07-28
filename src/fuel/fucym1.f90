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

subroutine fucym1 &
!================

 ( ncelet , ncel   ,                                              &
   indpdf ,                                                       &
   rtp    ,                                                       &
   f1m    , f3m    , f4m   , f1cl  , f3cl , f4cl ,                &
   f4m1   , f4m2   , d4cl  , d4f4  , hrec ,                       &
   fuel1  , fuel3  , oxyd  , prod1 , prod2  ,                     &
   xiner  , xh2s   , xso2  , f4s3no )

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
! fuel1            ! tr !  <- ! fraction massique fovm                         !
! fuel3            ! tr !  <- ! fraction massique co                           !
! oxyd             ! tr !  <- ! fraction massique o2                           !
! prod1            ! tr !  <- ! fraction massique co2                          !
! prod2            ! tr !  <- ! fraction massique h2o                          !
! xiner            ! tr !  <- ! fraction massique n2                           !
! xh2s             ! tr !  <- ! fraction massique h2s                          !
! xso2             ! tr !  <- ! fraction massique so2                          !
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

implicit none

!==============================================================================
! Common blocks
!==============================================================================

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "parall.h"
include "pointe.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "fuincl.h"
include "ppincl.h"
include "ppcpfu.h"

!===============================================================================

! Arguments

integer          ncelet , ncel
integer          indpdf(ncelet)

double precision rtp(ncelet,*)
double precision f1m(ncelet) , f3m(ncelet)  , f4m(ncelet)
double precision f1cl(ncelet), f3cl(ncelet) , f4cl(ncelet)

double precision f4m1(ncelet) , f4m2(ncelet) , d4cl(ncelet)
double precision d4f4(ncelet) , hrec(ncelet)

double precision fuel1(ncelet), fuel3(ncelet)
double precision oxyd(ncelet) , prod1(ncelet), prod2(ncelet)
double precision xiner(ncelet)
double precision xh2s(ncelet) , xso2(ncelet)
double precision f4s3no(ncelet)

! Local variables

integer          iel, icla

integer          n1 , n2 , n3 , n4 , n5 , n6 , n7 , n8 , n9
integer          n10, n11
integer          nbrint
parameter        (nbrint = 11)
integer          inttmp(nbrint)

double precision wmfov , wmco   , wmco2 , wmh2o
double precision wmo2   , wmn2   , wmc
double precision wmh2s , wmso2

double precision anu1 , anu2 , anu3
double precision somm

double precision zfov0 , zco0   , zco20
double precision zn20  , zh2o0  , zo20
double precision zfov1 , zco1   , zco21
double precision zn21  , zh2o1  , zo21
double precision zfov2 , zco2   , zco22
double precision zn22  , zh2o2  , zo22
double precision zfov3 , zco3   , zco23
double precision zn23  , zh2o3  , zo23
double precision reac1 , reac2  , reac3

!     Variables pour le support de la Pdf et sa description

double precision zzcl(8),zzs1(8),zzs2(8),zzs3(8),zzf4(8)
double precision zz(8)

integer ii
double precision bb1,bb2
double precision zh2s0,zh2s1,zh2s2,zh2s3
double precision zso20,zso21,zso22,zso23
double precision f4s1,f4s2,f4s3,zco2t
double precision x2t

!===============================================================================

!===============================================================================
! 1. CALCULS PRELIMINAIRES
!===============================================================================

! Initialize variables to avoid compiler warnings

zco23 = 0.d0
zco2t = 0.d0
zco3 = 0.d0
zh2o3 = 0.d0
zh2s3 = 0.d0
zn23 = 0.d0
zo23 = 0.d0
zfov3 = 0.d0
zso23 = 0.d0

! --- Masses molaires

wmco   = wmole(ico  )
wmo2   = wmole(io2  )
wmco2  = wmole(ico2 )
wmh2o  = wmole(ih2o )
wmn2   = wmole(in2)
wmfov  = wmole(ifov)
wmh2s  = wmole(ih2s)
wmso2  = wmole(iso2)

wmc    = wmolat(iatc)


! --- Calcul des coefficient relatifs a l'expression de
!        ZCHx1m0,  ZCO0, ZO20 et ZN20
!        fonctions de F1,  F3 et F4

!     Rappel : 1 = vapeur de fuel, 2 = non utilisé
!              3 = C sous forme CO provenant de l'oxydation du coke
!              4 = air

! ---- Initialisation des fractions massiques avec de l'air

do iel = 1, ncel
  fuel1(iel)  = zero
  fuel3(iel)  = zero
  oxyd(iel)   = wmo2*ao2f4
  prod1(iel)  = zero
  prod2(iel)  = zero
  xiner(iel)  = wmn2*an2f4
  xh2s (iel)  = zero
  xso2 (iel)  = zero
enddo

!  ANU1 : nb de moles d'oxygène pour convertir
!         une mole de FOV en CO et H2O

anu1   = 0.5d0+0.25d0*nhcfov

!  ANU2 nb de moles d'oxygènes pour la réaction H2S -> SO2 + H2O

anu2 = 1.5d0

!  ANU3 nb de mole d'oxygene pour la reaction CO -> CO2

anu3   = 0.5d0

if ( ieqnox .eq. 1 ) then
  do iel=1,ncel
    f4s3no(iel) = 0.d0
  enddo
endif

!===============================================================================
! 3. CALCUL DE LA COMPOSITION DU MELANGE SANS LES PDF
!    SI LES FLUCTUATIONS DE S SONT TROP FAIBLES
!===============================================================================


do iel = 1, ncel

  if ( indpdf(iel).eq.0 ) then

    if ( ieqco2 .eq. 1 ) then

      x2t = 0.d0
      do icla = 1, nclafu
        x2t =  x2t+rtp(iel,isca(iyfol(icla)))
      enddo

      zco2t = ( rtp(iel,isca(iyco2))/(1.d0-x2t) )/wmco2

    endif

! --> Composition de la phase gazeuse avant combustion (indice 0)

!          ZFOV0  = AFOVF1*F1
!          ZCO0   = ACOF1*F1 +  ACOF3*F3
!          ZCO20  = ZERO
!          ZH2S0  = AH2SF1*F1 + AH2SF3*F3
!          ZSO20  = ZERO
!          ZH2O0  = AH2OF3*F3
!          ZO20   = AO2F4*F4 + AO2F3*F3
!          ZN20   = AN2F4*F4

! ..v.7..x....v    ....x....v....x....v....x....v....L....v....x....v....x.I
    zfov0  = afovf1* f1m(iel)
    zco0   = acof1 * f1m(iel) + acof3 * f3m(iel)
    zh2s0  = ah2sf1* f1m(iel) + ah2sf3* f3m(iel)
    zo20   = ao2f4 * f4m(iel) + ao2f3 * f3m(iel)
    zn20   = an2f4 * f4m(iel) + an2f3 * f3m(iel)
    zco20  = zero
    zso20  = zero
    zh2o0  = ah2of3 * f3m(iel)

! --> Calcul de la composition du melange

    reac1  = min(zfov0, (zo20/anu1))
    zfov1 = zfov0  -           reac1
    zo21   = zo20  -      anu1*reac1
    zco1   = zco0  +           reac1
    zco21  = zco20
    zh2o1  = zh2o0 + 0.5d0*nhcfov*reac1
    zn21   = zn20
    zh2s1  = zh2s0
    zso21  = zso20

    reac2 = min(zh2s1,(zo21/anu2))
    zfov2 = zfov1
    zco2 = zco1
    zo22  = zo21  - anu2 * reac2
    zh2s2 = zh2s1 -        reac2
    zh2o2 = zh2o1 +        reac2
    zso22 = zso21 +        reac2
    zco22 = zco21
    zn22  = zn21

    if ( ieqco2 .eq. 0 ) then

      reac3 = min(zco2,(zo22/anu3))
      zfov3 = zfov2
      zo23  = zo22  - anu3*reac3
      zco3  = zco2  -      reac3
      zco23 = zco22 +      reac3
      zh2o3 = zh2o2
      zn23  = zn22
      zh2s3 = zh2s2
      zso23 = zso22

    else if ( ieqco2 .eq. 1 ) then

!           Transport de CO2

      zfov3 = zfov2
      zco2t = min(zco2t,zco2,zo22/anu3)

      zo23  = zo22  - anu3*zco2t
      zco3  = zco2  -      zco2t
      zco23 = zco2t
      zh2o3 = zh2o2
      zn23  = zn22
      zh2s3 = zh2s2
      zso23 = zso22

    else if ( ieqco2 .ne. 0 ) then

      write(nfecra,3000) ieqco2
      call csexit(1)

    endif

    fuel1(iel) = zfov3 * wmfov
    fuel3(iel) = zco3  * wmco
    prod1(iel) = zco23 * wmco2
    prod2(iel) = zh2o3 * wmh2o
    oxyd(iel)  = zo23  * wmo2
    xiner(iel) = zn23  * wmn2
    xh2s (iel) = zh2s3 * wmh2s
    xso2 (iel) = zso23 * wmso2

  endif

enddo



!===============================================================================
! 4. CALCUL DE LA COMPOSITION DU MELANGE AVEC LES PDF
!===============================================================================

!  14/08/2006
!  En toute rigueur, on pourrait calculer l'enthalpie de mélange aux 5 pts
!  en déduire la température et intégrer l'inverse de la masse volumique
!  mais ... les flammes fuligineuses sont rarement adiabatiques.
!  Il faudrait donc rechercher l'enthalpie aux points stoechio comme variable
!  auxilliaire pour retrouver l'enthalpie moyenne (qui tient compte des pertes)
!  En mettant les choses au mieux (i.e. on saurait formuler une hypothèse de
!  répartition des pertes entre les 3 points) on en déduirait une température
!  linéaire par morceaux, les fractions mmassiques étant linéaires,
!  1/Rho serait quadratique.
!  Déjà fait une fois (flamme de diff) mais pénible.

!  Pour l'instant on conserve le calcul de la température et de la masse
!  volumique d'aprés les concentrations moyennes dans FUPHY1
!  comme pour le cp


do iel = 1, ncel

  if ( indpdf(iel).eq.1 ) then

    if ( ieqco2 .eq. 1 ) then

      x2t = 0.d0
      do icla = 1, nclafu
        x2t =  x2t+rtp(iel,isca(iyfol(icla)))
      enddo

      zco2t = ( rtp(iel,isca(iyco2))/(1.d0-x2t) )/wmco2

    endif

!         Concentrations aux extrémités

    do ii = 1,8
      zzcl(ii) = zero
      zzf4(ii) = zero
    enddo

    zzcl(ifov) = afovf1*f1cl(iel)
    zzcl(ico)  = acof1 *f1cl(iel)+acof3 *f3cl(iel)
    zzcl(ih2s) = ah2sf1*f1cl(iel)+ah2sf3*f3cl(iel)
    zzcl(ih2o) =                  ah2of3*f3cl(iel)
    zzcl(in2)  =                  an2f3 *f3cl(iel)
    zzf4(io2)  = ao2f4
    zzf4(in2)  = an2f4

!         Calcul de la stoechiométrie de la première réaction
!         FOV + O2 => CO + H2O
!         FOV et O2 sont, avant la première réaction linéaire en F4
!         ZZ(F4,II) = ZZcl(II) + (F4-F4cl)/(1-F4cl) * (ZZF4(II)-ZZcl(II))
!         On cherche le point tel que
!         ZZ(F4s1,IFOV) = ZZ(F4s1,IO2) / (0.5+0.25*nHCFOV)
!         ANU1*ZZ(F4s1,IFOV) = ZZ(F4s1,IO2)
!         ANU1*(1-F4cl)*ZZcl(IFOV) + ANU1*(F4s1-F4cl)*(ZZF4(IFOV)-ZZcl(IFOV)) =
!              (1-F4cl)*ZZcl(IO2)  +      (F4s1-F4cl)*(ZZF4(IO2) -ZZcl(IO2))
!         on remarque que ZZcl(IO2) = ZZF4(IFOV) = 0

    f4s1 = ( anu1 * zzcl(ifov) + f4cl(iel) * zzf4(io2) )          &
         / ( anu1 * zzcl(ifov) +             zzf4(io2) )

!         Calcul des concentrations au premier point stoechio
!         avant réaction

    do ii = 1,8
      zzs1(ii) = zzcl(ii)                                         &
           + (f4s1-f4cl(iel))/(1.d0-f4cl(iel))                    &
           *(zzf4(ii)-zzcl(ii))
    enddo

!         Le nombre de moles de réaction est le nombre de moles de fov

    zzs1(ico) = zzs1(ico) + zzs1(ifov)
    zzs1(io2) = zero
    zzs1(ih2o)= zzs1(ih2o) + 0.5d0*nhcfov*zzs1(ifov)
    zzs1(ifov)= zero

!         Calcul de la stoechiometrie de la deuxième réaction
!         H2S + O2 => SO2 + H2O
!         Les concentrations sont désormais linéaires entre F4s1 et 1
!         ZZ(F4,II)  = ZZs1(II) + (F4-F4s1)/(1-F4s1)*(ZZF4(II)-ZZs1(II))
!         On cherche le point tel que
!         ZZ(f4s2,IH2S) = ZZ(F4S2,IO2)/ANU2
!         a nouveau ZZs1(IO2) = ZZF4(IH2S) = 0
!         La formule est semblable à la précedente en substituant
!         ANU1 par ANU2
!         IFOV par IH2S
!         cl par s1
!         s1 par s2

    f4s2 = ( anu2 * zzs1(ih2s) + f4s1 * zzf4(io2) )               &
          /( anu2 * zzs1(ih2s) +        zzf4(io2) )

!     Calcul des concentrations au deuxième point stoechio
!     avant réaction

    do ii = 1,8
      zzs2(ii) = zzs1(ii)                                         &
           + (f4s2-f4s1) / (1.d0-f4s1) * (zzf4(ii) - zzs1(ii))
    enddo

!         Le nombre de moles de réaction est le nombre de moles de H2S

    zzs2(iso2) = zzs2(iso2) + zzs2(ih2s)
    zzs2(ih2o) = zzs2(ih2o) + zzs2(ih2s)
    zzs2(io2)  = zero
    zzs2(ih2s) = zero

!         Calcul de la stoechiometrie de la troisième réaction
!         CO + O2 => CO2
!         Les concentrations sont désormais linéaires entre F4S2 et 1
!         ZZ(F4,II)  = ZZS2(II) + (F4-F4S2)/(1-F4S2)*(ZZF4(II)-ZZS2(II))
!         On cherche le point tel que
!         ZZ(f4s3,ICO) = ZZ(F4S3,IO2)/ANU3
!         La formule est semblable à la précedente en substituant
!         ANU2 par ANU3
!         IH2S par ICO
!         s1 par s2
!         s2 par s3

    f4s3 = ( anu3 * zzs2(ico ) + f4s2 * zzf4(io2) )               &
         / ( anu3 * zzs2(ico ) +        zzf4(io2) )

! si calcul de NOx, on stocke F4S3 pour l'utiliser dans fucyno

    if ( ieqnox .eq. 1 ) then
      f4s3no(iel) = f4s3
    endif

!         Calcul des concentrations au troisième point stoechio
!         avant réaction

    do ii = 1,8
      zzs3(ii) = zzs2(ii)                                         &
           + (f4s3-f4s2) / (1.d0-f4s2) * (zzf4(ii) - zzs2(ii))
    enddo

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
    do ii = 1,8
      zz(ii) = d4cl(iel)*zzcl(ii)+d4f4(iel)*zzf4(ii)
    enddo

!        Intégration sur le premier intervalle de richesse (entre cl et s1)

    bb1 = max(f4m1(iel),f4cl(iel))
    bb2 = min(f4m2(iel),f4s1)
    if( bb2.gt.bb1 ) then
      do ii = 1,8
        zz(ii) =  zz(ii) + hrec(iel)*(bb2-bb1)/(f4s1-f4cl(iel))   &
             *( (zzcl(ii)*f4s1-zzs1(ii)*f4cl(iel))                &
             +   (zzs1(ii)-zzcl(ii))*(bb1+bb2)*0.5d0)
      enddo
    endif

!        Intégration sur le deuxième intervalle de (entre s1 et s2)

    bb1 = max(f4m1(iel),f4s1)
    bb2 = min(f4m2(iel),f4s2)
    if( bb2.gt.bb1 ) then
      do ii = 1,8
        zz(ii) = zz(ii) + hrec(iel)*(bb2-bb1)/(f4s2-f4s1)         &
             * ( (zzs1(ii)*f4s2-zzs2(ii)*f4s1)                    &
             +   (zzs2(ii)-zzs1(ii))*(bb1+bb2)*0.5d0)
      enddo
    endif

!        Intégration de s2 à s3

    bb1 = max(f4m1(iel),f4s2)
    bb2 = min(f4m2(iel),f4s3)
    if( bb2.gt.bb1 ) then
      do ii = 1,8
        zz(ii) = zz(ii) + hrec(iel)*(bb2-bb1)/(f4s3-f4s2)         &
             * ( (zzs2(ii)*f4s3-zzs3(ii)*f4s2)                    &
             +   (zzs3(ii)-zzs2(ii))*(bb1+bb2)*0.5d0)
      enddo
    endif

!       Intégration de s3 à f4

    bb1 = max(f4m1(iel),f4s3)
    bb2 = min(f4m2(iel),1.d0)
    if ( bb2.gt.bb1 ) then
      do ii = 1,8
        zz(ii) = zz(ii) + hrec(iel)*(bb2-bb1)/(1.d0-f4s3)         &
             * ( (zzs3(ii)*1.d0-zzf4(ii)*f4s3)                    &
             +(zzf4(ii)-zzs3(ii))*(bb1+bb2)*0.5d0)
      enddo
    endif

   if ( ieqco2 .eq. 1 ) then

!  Transport de CO2

     zco2t = min(zco2t,zz(ico),2.d0*zz(io2))

     zz(ico ) = zz(ico) - zco2t
     zz(io2 ) = zz(io2) - anu3*zco2t
     zz(ico2) = zco2t

   else if ( ieqco2 .ne. 0 ) then

     write(nfecra,3000) ieqco2
     call csexit(1)

   endif

    fuel1(iel) = zz(ifov ) * wmfov
    fuel3(iel) = zz(ico  ) * wmco
    prod1(iel) = zz(ico2 ) * wmco2
    prod2(iel) = zz(ih2o ) * wmh2o
    oxyd(iel)  = zz(io2  ) * wmo2
    xiner(iel) = zz(in2  ) * wmn2
    xh2s(iel)  = zz(ih2s ) * wmh2s
    xso2(iel)  = zz(iso2 ) * wmso2

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

! --> Controle des parametres de la pdf

n1 = ncel

do iel = 1, ncel
  if ( indpdf(iel).ne.0 ) then
    n2 = n2 +1
  endif
enddo

! --> Controle des differentes valeurs des fractions massiques

do iel = 1, ncel

  somm = fuel1(iel) + fuel3(iel) + oxyd(iel)                      &
       + prod1(iel) + prod2(iel) + xiner(iel)                     &
       + xh2s(iel) + xso2(iel)

  if ( abs(somm-1.d0).lt.epsicp )    n3  = n3 +1
  if ( fuel1(iel).lt.(-epzero) .or.                               &
       fuel1(iel).gt.(1.d0+epzero) ) n4 = n4+1
  if ( fuel3(iel).lt.(-epzero) .or.                               &
       fuel3(iel).gt.(1.d0+epzero) ) n5 = n5+1
  if ( oxyd(iel).lt.(-epzero) .or.                                &
       oxyd(iel).gt.(1.d0+epzero)  ) n6 = n6+1
  if ( xiner(iel).lt.(-epzero) .or.                               &
       xiner(iel).gt.(1.d0+epzero) ) n7 = n7+1
  if ( prod1(iel).lt.(-epzero) .or.                               &
       prod1(iel).gt.(1.d0+epzero) ) n8 = n8+1
  if ( prod2(iel).lt.(-epzero) .or.                               &
       prod2(iel).gt.(1.d0+epzero) ) n9 = n9+1
  if ( xh2s(iel).lt.(-epzero) .or.                                &
       xh2s(iel).gt.(1.d0+epzero) ) n10 = n10+1
  if ( xso2(iel).lt.(-epzero) .or.                                &
       xso2(iel).gt.(1.d0+epzero) ) n11 = n11+1

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
  call parism(nbrint,inttmp)
  !==========
endif

write(nfecra,1000) n1 , n2

write(nfecra,2200) n3 , n4, n5, n6, n7, n8,                       &
                        n9, n10, n11


!-------
! FORMAT
!-------

 1000 format (/,                                                  &
'MODELISATION DE LA COMBUSTION AVEC LE MODELE DE DIFFUSION ',     &
'TURBULENTE (FUCYM1)',/,                                    &
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
' Nb de points YFOV ,YCO    < 0 ou > 1 = ',I9,I9,/,         &
' Nb de points YO2  ,YN2    < 0 ou > 1 = ',I9,I9,/,         &
' Nb de points YCO2 ,YH2O   < 0 ou > 1 = ',I9,I9,/,         &
' Nb de points YH2S ,YSO2   < 0 ou > 1 = ',I9,I9,/)

 3000 format(/,                                                   &
'MODELE DE CO FUEL ACTIVE : ',/,                            &
'       AVEC  IEQCO2 = ',I2,/,                              &
'HORS SEUL LES OPTIONS 0 et 1 SONT DISPONIBLES',/,          &
'    ARRET DU CALCUL : VERIFIER USPPMO')

return
end subroutine
