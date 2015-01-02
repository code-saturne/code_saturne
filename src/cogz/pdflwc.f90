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

subroutine pdflwc &
!================

 ( ncelet , ncel  ,                                               &
   fm     , fp2m  , yfm    , yfp2m ,                              &
   propce )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DES PARAMETRES DE LA PDF
! PDF LIBBY - WILLIAMS 2 POINTS AVEC HYPOTHESE DE CURL


! COMMENTAIRES :
! ------------

!    Dans un diagramme (F, Yf), on construit deux droites:
!         - La droite de combustion complete
!         - La droite de melange

!    Dans ce domaine, nous allons trouver deux pics qui
!    definiront une troisieme droite sur laquelle on definit
!    une abscisse curviligne G.


! LE RESULTAT EST :
! ---------------

!    CALCUL DES PARAMETRES ASSOCIES AUX FONCTIONS DIRAC

!      Les Diracs sont en position [F(.,1),Y(.,1)] et [F(.,2),Y(.,2)]
!      Leurs amplitudes respectives sont D(.,1) et D(.,2)
!      Pour chaque dirac,
!          on calcule la temperature [T(.,1), T(.,2)]
!                  la masse volumique [RHO(.,1), RHO(.,2)]
!                  le terme source chimique [W(.,1),W(.,2)]


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
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!==============================================================================
! Module files
!==============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use pointe
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use field
!===============================================================================

implicit none

! Arguments

integer          ncelet, ncel

double precision fm(ncelet)   , fp2m(ncelet)
double precision yfm(ncelet)  , yfp2m(ncelet)
double precision propce(ncelet,*)

! Local variables

integer          iel, igg, idirac
integer          mode

double precision coefg(ngazgm), epsi
double precision yfuel
double precision yoxyd, yo2
double precision yprod,fmp,fp2mp,yfmp,yfp2mp

! --- Expression des droites

double precision aa1, aa2
double precision aa3, aa6
double precision bb1, bb2
double precision bb3, bb6
double precision f12, f13, f14
double precision f15, f16

! --- Extrema

double precision gmin, gmax, gp2m
double precision g1max, g2max, g3max
double precision g1min, g2min
double precision g3min, g4min

! --- Position des Diracs

double precision f(ndracm)    , y(ndracm) , d(ndracm)
double precision g(ndracm)    , h(ndracm)
double precision teml(ndracm) , maml(ndracm)
double precision w(ndracm)    , rhol(ndracm)
double precision theta(ndracm)

! --- Pointeurs & autres

integer ipctem, ipcmam
integer ipampl(ndracm), ipfmel(ndracm)
integer ipfmal(ndracm), ipteml(ndracm)
integer ipmaml(ndracm)
integer iprhol(ndracm), iptscl(ndracm)
integer ipcfue, ipcoxy, ipcpro, ipctsc
!      INTEGER IPCKAB, IPT4 , IPT3
double precision nbmol,  temsmm
double precision sum1, sum2, sum3, sum4,  sum5,  sum6 , sum16
double precision sum7, sum8, sum9, sum10, sum11, sum12, sum15
double precision sum17
double precision, dimension(:), pointer ::  crom
integer ipass
data    ipass /0/
save    ipass

!===============================================================================

! Initialize variables to avoid compiler warnings

gmin =  grand
gmax = -grand

!===============================================================================
! 0.  POSITION DES VARIABLES
!===============================================================================

do idirac = 1, ndirac
  ipampl(idirac) = ipproc(iampl(idirac))
  ipfmel(idirac) = ipproc(ifmel(idirac))
  ipfmal(idirac) = ipproc(ifmal(idirac))
  ipteml(idirac) = ipproc(iteml(idirac))
  ipmaml(idirac) = ipproc(imaml(idirac))
  iprhol(idirac) = ipproc(irhol(idirac))
  iptscl(idirac) = ipproc(itscl(idirac))
enddo
ipcfue = ipproc(iym(1))
ipcoxy = ipproc(iym(2))
ipcpro = ipproc(iym(3))
ipctsc = ipproc(itsc)
ipctem = ipproc(itemp)
call field_get_val_s(icrom, crom)
ipcmam = ipproc(imam)

!      IF ( IIRAYO.GT.0 ) THEN
!        IPCKAB = IPPROC(ICKABS)
!        IPT4  = IPPROC(IT4M)
!        IPT3  = IPPROC(IT3M)
!      ENDIF

! ---> Initialisation

do igg = 1, ngazgm
  coefg(igg) = zero
enddo
epsi = 1.d-09

!===============================================================================
! 1.  CALCULS PRELIMINAIRES
!===============================================================================

ipass = ipass + 1

do iel = 1, ncel


  fmp    =  max(min(fmax,fm(iel)),fmin)
  yfmp   =  max(min(yfm(iel),fmp),                                &
                    zero,(fmp-fs(1))/(1.d0-fs(1)))
  fp2mp  =  max(min(fp2m(iel),                                    &
               (min(fmax,(1-fs(1))*yfmp+fs(1))-fmp)               &
             *      (fmp-max(fmin,yfmp))),zero)
  yfp2mp =  max(min(yfp2m(iel),(fmp-yfmp)                         &
             *      (yfmp-max(zero,                               &
                    (fmp-fs(1))/(1.d0-fs(1))))),zero)


!===============================================================================
! 2. NON PASSAGE PAR LA PDF
!===============================================================================

  if ( (fp2mp.le.epsi).and.(yfp2mp.le.epsi )) then

    sum1 = zero
    sum2 = zero
    sum3 = zero
    sum4 = zero
    sum5 = zero
    sum6 = zero
    sum15 = zero
    sum16 = zero

    do idirac =1, ndirac
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

! ------ Masse molaire

      nbmol = 0.d0
      do igg = 1, ngazg
        nbmol = nbmol + coefg(igg)/wmolg(igg)
      enddo
      maml(idirac) = 1.d0/nbmol

! ------ Calcul de la temperature pour le pic 1 et 2

      mode    = 1
      call cothht                                                 &
      !==========
      ( mode   , ngazg , ngazgm  , coefg  ,                       &
        npo    , npot  , th      , ehgazg ,                       &
        h(idirac)     , teml(idirac)    )

! ---> Calcul de la masse volumique en 1 et 2

      if ( ipass.gt.1.or.                                         &
          (isuite.eq.1.and.initro.eq.1) ) then
        rhol(idirac) = p0 * maml(idirac)                   &
                     / (rr*teml(idirac))
      else
        rhol(idirac) = ro0
      endif

! ---> Calcul du terme source en 1 et 2 du scalaire YFM

      theta(idirac) = ta / teml(idirac)                           &
                    * (1.d0 - teml(idirac) / tstar)

      w(idirac) = vref / lref * (- d(idirac)*rhol(idirac)         &
                  * yfuel*yo2                                     &
                  * exp( -theta(idirac) ))

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

! ---> Stockage des proprietes via PROPCE

      propce(iel,ipampl(idirac)) = d(idirac)
      propce(iel,ipfmel(idirac)) = f(idirac)
      propce(iel,ipfmal(idirac)) = y(idirac)
      propce(iel,ipteml(idirac)) = teml(idirac)
      propce(iel,ipmaml(idirac)) = maml(idirac)
      propce(iel,iprhol(idirac)) = rhol(idirac)
      propce(iel,iptscl(idirac)) = w(idirac)

    enddo

    propce(iel,ipcmam) = sum1
    propce(iel,ipctem) = sum2
    temsmm             = sum3
    propce(iel,ipcfue) = sum4
    propce(iel,ipcoxy) = sum5
    propce(iel,ipcpro) = sum6
    propce(iel,ipctsc) = sum16

! ---> Masse volumique du melange

    if ( ipass.gt.1 .or.                                          &
        (isuite.eq.1.and.initro.eq.1) ) then
      crom(iel) = srrom * crom(iel)             &
                         + (1.d0-srrom) * (p0/(rr*temsmm))
    endif

  else

!===============================================================================
! 3.  PASSAGE PAR LA PDF
!===============================================================================

    if  ( ( (fp2mp.gt.epsi) .and. (yfp2mp.lt.epsi) )              &
      .or.( (fp2mp.lt.epsi) .and. (yfp2mp.gt.epsi) ) )then

! -- Choix suivant les differents cas

! --- CAS 1 : Droite Verticale
!     =====
      if (fp2mp .lt. epsi) then

! ---> Choix des extrema

        gmin =   max ( -  yfmp,                                   &
                    -  (yfmp-(fmp-fs(1))/(1.d0-fs(1))))
        gmax = fmp - yfmp

! ---> Calcul des amplitudes des Diracs

        d(1) = gmax / ( gmax-gmin )
        d(2) = (1.d0 - d(1))

! ---> Calcul de la variance de l'abscisse curviligne (GP2M)

        gp2m = fp2mp + yfp2mp

! ---> Test sur GP2M

        g(1) = -sqrt( -gmin/gmax*gp2m )
        g(2) =  -d(1) * g(1) / d(2)

        do idirac = 1, ndirac
          f(idirac) =                                             &
                    fmp + g(idirac) * sqrt( fp2mp/gp2m)
          y(idirac) =                                             &
                    yfmp + g(idirac) * sqrt( yfp2mp/gp2m)
        enddo

! --- CAS 2 : Droite Horizontale
!     =====

      elseif    (yfp2mp .lt. epsi) then

! ---> Calcul des differentes intersections des droites

        f12 = yfmp
        f13 = fs(1) + yfmp * (1-fs(1))
        f14 = fmin
        f15 = fmax

! ---> Calcul des extrema de l'abscisse curviligne
!      Appel de la function G

        call lwcgfu(g2max, f15, fmp, yfp2mp, fp2mp)
        call lwcgfu(g3max, f13, fmp, yfp2mp, fp2mp)
        call lwcgfu(g1min, f12, fmp, yfp2mp, fp2mp)
        call lwcgfu(g2min, f14, fmp, yfp2mp, fp2mp)

! ---> Choix des extrema

        gmin =   max(  g1min,  g2min)
        gmax =   min(  g2max,  g3max)
! ---> Calcul des amplitudes des Diracs

        d(1) = gmax / ( gmax-gmin )
        d(2) = (1.d0 - d(1))

! ---> Calcul de la variance de l'abscisse curviligne (GP2M)

        gp2m = fp2mp + yfp2mp

        g(1) = -sqrt( -gmin/gmax*gp2m )
        g(2) =  -d(1) * g(1) / d(2)

        do idirac = 1, ndirac
          f(idirac) =                                             &
              fmp + g(idirac) * sqrt( fp2mp/gp2m)
          y(idirac) =                                             &
              yfmp + g(idirac) * sqrt( yfp2mp/gp2m)
        enddo
      endif

    elseif ( (fp2mp.gt.epsi) .and. (yfp2mp.gt.epsi) ) then

! --- CAS 3 : Parallelisme avec la Droite de Melange
!     =====

      if ( ((yfp2mp / fp2mp) .lt. 1.d0 + epsi)                    &
        .and. ((yfp2mp / fp2mp) .gt. 1.d0 - epsi) ) then

        aa1 = 1.d0
        bb1 = yfmp - fmp

        aa3 = 1.d0 / (1.d0-fs(1))
        bb3 = -fs(1) / (1.d0-fs(1))

        aa6 = zero
        bb6 = zero

! ---> Calcul des differentes intersections de ces droites

        f13 = (bb3-bb1) / (aa1-aa3)
        f14 = fmin
        f15 = fmax
        f16 = (bb6-bb1) / (aa1-aa6)

! ---> Calcul des extrema de l'abscisse curviligne
!      Appel de la function G

        call lwcgfu(g2max,f15, fmp, yfp2mp, fp2mp)
        call lwcgfu(g3max,f13, fmp, yfp2mp, fp2mp)
        call lwcgfu(g2min,f14, fmp, yfp2mp, fp2mp)
        call lwcgfu(g3min,f16, fmp, yfp2mp, fp2mp)

        gmin =   max(  g2min,  g3min)
        gmax =   min(  g2max,  g3max)

! ---> Calcul des amplitudes des Diracs

        d(1) = gmax / ( gmax-gmin )
        d(2) = (1.d0 - d(1))

! ---> Calcul de la variance de l'abscisse curviligne (GP2M)

        gp2m = fp2mp + yfp2mp

! ---> Test sur GP2M

        g(1) = -sqrt( -gmin/gmax*gp2m )
        g(2) =  -d(1) * g(1) / d(2)

        do idirac = 1, ndirac
          f(idirac) =                                             &
              fmp + g(idirac) * sqrt( fp2mp/gp2m)
          y(idirac) =                                             &
              yfmp + g(idirac) * sqrt( yfp2mp/gp2m)
        enddo

! --- CAS 4 : Parallelisme avec la Droite de Combustion Complete
!     =====
      elseif ( ((sqrt(yfp2mp/fp2mp) * (1.d0-fs(1)))               &
                                        .lt. 1.d0 + epsi)         &
        .and.  ((sqrt(yfp2mp/fp2mp) * (1.d0-fs(1)))               &
                                        .gt. 1.d0 - epsi) ) then

        aa1 = sqrt( yfp2mp/fp2mp )
        bb1 = yfmp - sqrt( yfp2mp/fp2mp ) * fmp

        aa2 = 1.d0
        bb2 = zero

        aa6 = zero
        bb6 = zero

! ---> Calcul des differentes intersections de ces droites

        f12 = (bb2-bb1) / (aa1-aa2)
        f14 = fmin
        f15 = fmax
        f16 = (bb6-bb1) / (aa1-aa6)

! ---> Calcul des extrema de l'abscisse curviligne
!      Appel de la function G

        call lwcgfu(g2max,f15, fmp, yfp2mp, fp2mp)
        call lwcgfu(g1max,f12, fmp, yfp2mp, fp2mp)
        call lwcgfu(g2min,f14, fmp, yfp2mp, fp2mp)
        call lwcgfu(g3min,f16, fmp, yfp2mp, fp2mp)

        gmin =   max(  g2min,  g3min)
        gmax =   min(  g1max,  g2max)

! ---> Calcul des amplitudes des Diracs

        d(1) = gmax / ( gmax-gmin )
        d(2) = (1.d0 - d(1))

! ---> Calcul de la variance de l'abscisse curviligne (GP2M)

        gp2m = fp2mp + yfp2mp

! ---> Test sur GP2M

        g(1) = -sqrt( -gmin/gmax*gp2m )
        g(2) =  -d(1) * g(1) / d(2)

        do idirac = 1, ndirac
          f(idirac) =                                             &
              fmp + g(idirac) * sqrt( fp2mp/gp2m)
          y(idirac) =                                             &
             yfmp + g(idirac) * sqrt( yfp2mp/gp2m)
        enddo

! --- CAS GENERAL
!==========
      else

! ---> Expression des differentes droites: Y=AX+B


        aa1 = -sqrt( yfp2mp/fp2mp )
        bb1 = yfmp - aa1 * fmp

        aa2 = 1.d0
        bb2 = zero

        aa3 = 1.d0 / (1.d0-fs(1))
        bb3 = -fs(1) / (1.d0-fs(1))

        aa6 = zero
        bb6 = zero

! ---> Calcul des differentes intersections de ces droites

        f12 = (bb2-bb1) / (aa1-aa2)
        f13 = (bb3-bb1) / (aa1-aa3)
        f14 = fmin
        f15 = fmax
        f16 = (bb6-bb1) / (aa1-aa6)

! ---> Calcul des extrema de l'abscisse curviligne
!      Appel de la function G

        call lwcgfu(g1max, f12, fmp, yfp2mp, fp2mp)
        call lwcgfu(g2max, f15, fmp, yfp2mp, fp2mp)
        call lwcgfu(g3max, f13, fmp, yfp2mp, fp2mp)
        call lwcgfu(g1min, f12, fmp, yfp2mp, fp2mp)
        call lwcgfu(g2min, f14, fmp, yfp2mp, fp2mp)
        call lwcgfu(g3min, f16, fmp, yfp2mp, fp2mp)
        call lwcgfu(g4min, f13, fmp, yfp2mp, fp2mp)

!===============================================================================
! 4.  CALCUL DES PARAMETRES DES DEUX PICS DE DIRAC
!===============================================================================

! ---> Choix des extrema suivant les pentes des droites

        if ( aa1 .gt. aa3 ) then
          gmin =   max (  g2min,  g3min,  g4min )
          gmax =   min ( g1max, g2max )

        elseif (( aa1 .gt. aa2 ).and. ( aa1 .le. aa3 )) then
          gmin =   max (  g2min,  g3min )
          gmax =   min ( g1max, g2max, g3max )

        elseif ( aa1 .le. aa2 ) then
          gmin =   max (  g1min,  g2min,  g3min )
          gmax =   min ( g2max, g3max )
        endif

! ---> Calcul des amplitudes des Diracs

        d(1) = gmax / ( gmax-gmin )
        d(2) = (1.d0 - d(1))

! ---> Calcul de la variance de l'abscisse curviligne (GP2M)

        gp2m = fp2mp + yfp2mp

        g(1) = -sqrt( gmin/gmax*gp2m )
        g(2) =  -d(1) * g(1) / d(2)

        do idirac = 1, ndirac
          f(idirac) =                                             &
              fmp + g(idirac) * sqrt( fp2mp/gp2m)
          y(idirac) =                                             &
              yfmp + g(idirac) * sqrt( yfp2mp/gp2m)
        enddo

      endif
    endif

!===============================================================================
! 5.  DETERMINATION DES GRANDEURS THERMOCHIMIQUES DES DEUX PICS
!===============================================================================

! ---> Calcul de l'enthalpies en 1 et 2

    sum7  = zero
    sum8  = zero
    sum9  = zero
    sum10 = zero
    sum11 = zero
    sum12 = zero
    sum17 = zero

    do idirac = 1, ndirac
      h(idirac) = ((hmax-hmin)*f(idirac) + hmin*fmax - hmax*fmin) &
                / (fmax-fmin)

! ---> Calcul de la fraction massique des gaz (F, O et P) en 1 et 2

      yfuel = y(idirac)
      yoxyd = 1.d0 - (coeff3+1.0d0)*f(idirac) + coeff3*y(idirac)
      yprod = 1.d0 - yfuel - yoxyd
      yo2   = coeff1 - (coeff1 + coeff2) * f(idirac)              &
                      + coeff2 * y(idirac)

! ---> Coefficients d'absorption pour les pics 1 et 2

!            IF ( IIRAYO .GT. 0  ) THEN
!              KABSGF = YFUEGF(IEL)*KABSG(1) + YOXYGF(IEL)*KABSG(2)
!     &                                      + YPROGF(IEL)*KABSG(3)
!              KABSGB = YFUEGB(IEL)*KABSG(1) + YOXYGB(IEL)*KABSG(2)
!     &                                      + YPROGB(IEL)*KABSG(3)
!            ENDIF


! ---> Calcul de la masse molaire et de la temperature en 1 et 2

      coefg(1) = yfuel
      coefg(2) = yoxyd
      coefg(3) = yprod

! ------ Masse molaire pour le pic 1 et 2

      nbmol = 0.d0
      do igg = 1, ngazg
        nbmol = nbmol + coefg(igg)/wmolg(igg)
      enddo
      maml(idirac) = 1.d0/nbmol

! ------ Calcul de la temperature pour le pic 1 et 2

      mode    = 1
      call cothht                                                 &
      !==========
      ( mode   , ngazg , ngazgm  , coefg  ,                       &
        npo    , npot  , th      , ehgazg ,                       &
        h(idirac)     , teml(idirac)    )

! ---> Calcul de la masse volumique en 1 et 2

      if ( ipass.gt.1 .or.                                        &
          (isuite.eq.1.and.initro.eq.1) ) then
        rhol(idirac) = p0 * maml(idirac)                   &
                     / (rr*teml(idirac))
      else
        rhol(idirac) = ro0
      endif

! ---> Calcul du terme source en 1 et 2 du scalaire YFM

      theta(idirac) = ta / teml(idirac)                           &
                    * (1.d0 - teml(idirac) / tstar)

      w(idirac) = vref / lref * (- d(idirac)*rhol(idirac)         &
                  * yfuel*yo2                                     &
                  * exp( -theta(idirac) ))

! ---> Masse molaire du melange

      sum7 = sum7 + d(idirac)*maml(idirac)

! ---> Temperature du melange

      sum8 = sum8 + d(idirac)*teml(idirac)

! ---> Temperature / Masse molaire

      sum9 = sum9 + d(idirac)*teml(idirac)/maml(idirac)

! ---> Fractions massiques des especes globales

      sum10 = sum10 + yfuel*d(idirac)

      sum11 = sum11 + yoxyd*d(idirac)

      sum12 = sum12 + yprod*d(idirac)

      sum17 = sum17 +w(idirac)

! ---> Stockage des proprietes via PROPCE

      propce(iel,ipampl(idirac)) = d(idirac)
      propce(iel,ipfmel(idirac)) = f(idirac)
      propce(iel,ipfmal(idirac)) = y(idirac)
      propce(iel,ipmaml(idirac)) = maml(idirac)
      propce(iel,ipteml(idirac)) = teml(idirac)
      propce(iel,iprhol(idirac)) = rhol(idirac)
      propce(iel,iptscl(idirac)) = w(idirac)

! ---> Grandeurs relatives au rayonnement

!            IF ( IIRAYO .GT. 0  ) THEN
!              PROPCE(IEL,IPCKAB) = YGFM*KABSGF + YGBM*KABSGB
!              PROPCE(IEL,IPT4)    = YGFM*TGF**4+YGBM*TGB**4
!              PROPCE(IEL,IPT3)    = YGFM*TGF**3+YGBM*TGB**3
!            ENDIF

    enddo

    propce(iel,ipcmam) = sum7
    propce(iel,ipctem) = sum8
    temsmm             = sum9
    propce(iel,ipcfue) = sum10
    propce(iel,ipcoxy) = sum11
    propce(iel,ipcpro) = sum12
    propce(iel,ipctsc) = sum17

! ---> Masse volumique du melange

  if (ipass.gt.1.or.(isuite.eq.1.and.initro.eq.1)) then
      crom(iel) = srrom * crom(iel)             &
                       + (1.d0-srrom) * (p0/(rr*temsmm))
    endif

  endif

enddo

!----
! FIN
!----

end subroutine
