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

subroutine pdflwc &
!================

 ( ncelet , ncel  ,                                               &
   fm     , fp2m  , yfm    , yfp2m )

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
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          ncelet, ncel

double precision fm(ncelet)   , fp2m(ncelet)
double precision yfm(ncelet)  , yfp2m(ncelet)

! Local variables

integer          iel, igg, idirac

double precision coefg(ngazgm), epsi
double precision yfuel
double precision yoxyd, yo2
double precision yprod,fmp,fp2mp,yfmp,yfp2mp

! --- Lines expression

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

! --- Diracs position

double precision f(ndracm)    , y(ndracm) , d(ndracm)
double precision g(ndracm)    , h(ndracm)
double precision teml(ndracm) , maml(ndracm)
double precision w(ndracm)    , rhol(ndracm)
double precision theta(ndracm)

double precision nbmol,  temsmm
double precision sum1, sum2, sum3, sum4,  sum5,  sum6 , sum16
double precision sum7, sum8, sum9, sum10, sum11, sum12, sum15
double precision sum17
double precision, dimension(:), pointer ::  crom
double precision, dimension(:), pointer ::  cpro_temp, cpro_tsc, cpro_mam
double precision, dimension(:), pointer ::  cpro_ym1, cpro_ym2, cpro_ym3
type(pmapper_double_r1), dimension(:), pointer :: cpro_fmel, cpro_fmal
type(pmapper_double_r1), dimension(:), pointer :: cpro_tscl, cpro_rhol
type(pmapper_double_r1), dimension(:), pointer :: cpro_ampl, cpro_teml
type(pmapper_double_r1), dimension(:), pointer :: cpro_maml
integer ipass
data    ipass /0/
save    ipass

!===============================================================================

! Initialize variables to avoid compiler warnings

gmin =  grand
gmax = -grand

!===============================================================================
! 0.  initialization
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

! ---> Initialisation

do igg = 1, ngazgm
  coefg(igg) = zero
enddo
epsi = 1.d-09

!===============================================================================
! 1.  preliminary calculations
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
! 2. do not pass through PDF
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

! ---> enthalpy computation

      h(idirac) = ((hmax-hmin)*f(idirac) + hmin*fmax - hmax*fmin) &
                / (fmax-fmin)

! ---> computation of mass fractions of gas (F, O et P)

      yfuel = y(idirac)
      yoxyd = 1.d0 - (coeff3+1.0d0)*f(idirac) + coeff3*y(idirac)
      yprod = 1.d0 - yfuel - yoxyd
      yo2   = coeff1 - (coeff1 + coeff2) * f(idirac)              &
                      + coeff2 * y(idirac)

! ---> molar mass and temperature computation

      coefg(1) = yfuel
      coefg(2) = yoxyd
      coefg(3) = yprod

! ------ molar mass

      nbmol = 0.d0
      do igg = 1, ngazg
        nbmol = nbmol + coefg(igg)/wmolg(igg)
      enddo
      maml(idirac) = 1.d0/nbmol

! ------ temperature computation for peak 1 and 2

      teml(idirac) = cs_gas_combustion_h_to_t(coefg, h(idirac))

! ---> density computation in 1 and 2

      if ( ipass.gt.1.or.                                         &
          (isuite.eq.1.and.initro.eq.1) ) then
        rhol(idirac) = p0 * maml(idirac)                   &
                     / (cs_physical_constants_r*teml(idirac))
      else
        rhol(idirac) = ro0
      endif

! ---> source term computation in 1 and 2 for scalar yfm

      theta(idirac) = ta / teml(idirac)                           &
                    * (1.d0 - teml(idirac) / tstar)

      w(idirac) = vref / lref * (- d(idirac)*rhol(idirac)         &
                  * yfuel*yo2                                     &
                  * exp( -theta(idirac) ))

! ---> mix molar mass

      sum1 = sum1 + d(idirac)*maml(idirac)

! ---> mix temperature

      sum2 = sum2 + d(idirac)*teml(idirac)

! ---> temperature / molar mass

      sum3 = sum3 + d(idirac)*teml(idirac)/maml(idirac)

! ---> global species mass fractions

      sum4 = sum4 + yfuel*d(idirac)

      sum5 = sum5 + yoxyd*d(idirac)

      sum6 = sum6 + yprod*d(idirac)

      sum15 = sum15 +rhol(idirac)*d(idirac)

      sum16 = sum16 +w(idirac)

! ---> property storage in fields

      cpro_ampl(idirac)%p(iel) = d(idirac)
      cpro_fmel(idirac)%p(iel) = f(idirac)
      cpro_fmal(idirac)%p(iel) = y(idirac)
      cpro_teml(idirac)%p(iel) = teml(idirac)
      cpro_maml(idirac)%p(iel) = maml(idirac)
      cpro_rhol(idirac)%p(iel) = rhol(idirac)
      cpro_tscl(idirac)%p(iel) = w(idirac)

    enddo

    cpro_mam(iel)  = sum1
    cpro_temp(iel) = sum2
    temsmm         = sum3
    cpro_ym1(iel)  = sum4
    cpro_ym2(iel)  = sum5
    cpro_ym3(iel)  = sum6
    cpro_tsc(iel)  = sum16

! ---> mix density

    if ( ipass.gt.1 .or.              &
        (isuite.eq.1.and.initro.eq.1) ) then
      crom(iel) = srrom * crom(iel)                                   &
                + (1.d0-srrom) * (p0/(cs_physical_constants_r*temsmm))
    endif

  else

!===============================================================================
! 3.  pass through PDF
!===============================================================================

    if  ( ( (fp2mp.gt.epsi) .and. (yfp2mp.lt.epsi) )              &
      .or.( (fp2mp.lt.epsi) .and. (yfp2mp.gt.epsi) ) )then

! --- case 1 : vertical straight line
      if (fp2mp .lt. epsi) then

! ---> extrema choices

        gmin =   max ( -  yfmp,                                   &
                    -  (yfmp-(fmp-fs(1))/(1.d0-fs(1))))
        gmax = fmp - yfmp

! ---> Diracs amplitudes computation

        d(1) = gmax / ( gmax-gmin )
        d(2) = (1.d0 - d(1))

! ---> Calcul de la variance de l'abscisse curviligne (GP2M)

        gp2m = fp2mp + yfp2mp

! ---> test on gp2m

        g(1) = -sqrt( -gmin/gmax*gp2m )
        g(2) =  -d(1) * g(1) / d(2)

        do idirac = 1, ndirac
          f(idirac) =                                             &
                    fmp + g(idirac) * sqrt( fp2mp/gp2m)
          y(idirac) =                                             &
                    yfmp + g(idirac) * sqrt( yfp2mp/gp2m)
        enddo

! --- case 2 : horizontal straight line

      elseif    (yfp2mp .lt. epsi) then

! ---> computation of the different intersections of the lines

        f12 = yfmp
        f13 = fs(1) + yfmp * (1-fs(1))
        f14 = fmin
        f15 = fmax

! ---> curviline coordinate computation
!      call to G function

        call lwcgfu(g2max, f15, fmp, yfp2mp, fp2mp)
        call lwcgfu(g3max, f13, fmp, yfp2mp, fp2mp)
        call lwcgfu(g1min, f12, fmp, yfp2mp, fp2mp)
        call lwcgfu(g2min, f14, fmp, yfp2mp, fp2mp)

! ---> extrema choices

        gmin =   max(  g1min,  g2min)
        gmax =   min(  g2max,  g3max)
! ---> computation of Diracs amplitudes

        d(1) = gmax / ( gmax-gmin )
        d(2) = (1.d0 - d(1))

! ---> computation of curviline coordinate variance (gp2m)

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

! --- case 3 : parallelism and mix line

      if ( ((yfp2mp / fp2mp) .lt. 1.d0 + epsi)                    &
        .and. ((yfp2mp / fp2mp) .gt. 1.d0 - epsi) ) then

        aa1 = 1.d0
        bb1 = yfmp - fmp

        aa3 = 1.d0 / (1.d0-fs(1))
        bb3 = -fs(1) / (1.d0-fs(1))

        aa6 = zero
        bb6 = zero

! ---> computation of the different intersections of the lines

        f13 = (bb3-bb1) / (aa1-aa3)
        f14 = fmin
        f15 = fmax
        f16 = (bb6-bb1) / (aa1-aa6)

! ---> curviline coordinate computation
!      call to G function

        call lwcgfu(g2max,f15, fmp, yfp2mp, fp2mp)
        call lwcgfu(g3max,f13, fmp, yfp2mp, fp2mp)
        call lwcgfu(g2min,f14, fmp, yfp2mp, fp2mp)
        call lwcgfu(g3min,f16, fmp, yfp2mp, fp2mp)

        gmin =   max(  g2min,  g3min)
        gmax =   min(  g2max,  g3max)

! ---> computation of Diracs amplitudes

        d(1) = gmax / ( gmax-gmin )
        d(2) = (1.d0 - d(1))

! ---> computation of curviline coordinate variance (gp2m)

        gp2m = fp2mp + yfp2mp

! ---> Test on gp2m

        g(1) = -sqrt( -gmin/gmax*gp2m )
        g(2) =  -d(1) * g(1) / d(2)

        do idirac = 1, ndirac
          f(idirac) =                                             &
              fmp + g(idirac) * sqrt( fp2mp/gp2m)
          y(idirac) =                                             &
              yfmp + g(idirac) * sqrt( yfp2mp/gp2m)
        enddo

! --- case 4 : parallelism with complete combustion line

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

! ---> computation of the different intersections of the lines

        f12 = (bb2-bb1) / (aa1-aa2)
        f14 = fmin
        f15 = fmax
        f16 = (bb6-bb1) / (aa1-aa6)

! ---> curviline coordinate extrema computation
!      call to G function

        call lwcgfu(g2max,f15, fmp, yfp2mp, fp2mp)
        call lwcgfu(g1max,f12, fmp, yfp2mp, fp2mp)
        call lwcgfu(g2min,f14, fmp, yfp2mp, fp2mp)
        call lwcgfu(g3min,f16, fmp, yfp2mp, fp2mp)

        gmin =   max(  g2min,  g3min)
        gmax =   min(  g1max,  g2max)

! ---> Diracs amplitudes computations

        d(1) = gmax / ( gmax-gmin )
        d(2) = (1.d0 - d(1))

! ---> computation of curviline coordinate variance (gp2m)

        gp2m = fp2mp + yfp2mp

! ---> Test on gp2m

        g(1) = -sqrt( -gmin/gmax*gp2m )
        g(2) =  -d(1) * g(1) / d(2)

        do idirac = 1, ndirac
          f(idirac) =                                             &
              fmp + g(idirac) * sqrt( fp2mp/gp2m)
          y(idirac) =                                             &
             yfmp + g(idirac) * sqrt( yfp2mp/gp2m)
        enddo

! --- general case
      else

! ---> expression for the different lines: Y=AX+B


        aa1 = -sqrt( yfp2mp/fp2mp )
        bb1 = yfmp - aa1 * fmp

        aa2 = 1.d0
        bb2 = zero

        aa3 = 1.d0 / (1.d0-fs(1))
        bb3 = -fs(1) / (1.d0-fs(1))

        aa6 = zero
        bb6 = zero

! ---> computation of the different intersections of the lines

        f12 = (bb2-bb1) / (aa1-aa2)
        f13 = (bb3-bb1) / (aa1-aa3)
        f14 = fmin
        f15 = fmax
        f16 = (bb6-bb1) / (aa1-aa6)

! ---> curviline coordinate extrema computation
!      call to G function

        call lwcgfu(g1max, f12, fmp, yfp2mp, fp2mp)
        call lwcgfu(g2max, f15, fmp, yfp2mp, fp2mp)
        call lwcgfu(g3max, f13, fmp, yfp2mp, fp2mp)
        call lwcgfu(g1min, f12, fmp, yfp2mp, fp2mp)
        call lwcgfu(g2min, f14, fmp, yfp2mp, fp2mp)
        call lwcgfu(g3min, f16, fmp, yfp2mp, fp2mp)
        call lwcgfu(g4min, f13, fmp, yfp2mp, fp2mp)

!===============================================================================
! 4.  computation of the parameters of the two Dirac peaks
!===============================================================================

! ---> extrema choice according to the lines slope

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

! ---> Dirac amplitudes computations

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
! 5.  computations of the thermochemical quantities at the two peaks
!===============================================================================

! ---> enthalpy computation in 1 and 2

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

! ---> computation of the mass fraction of the gas (F, O et P) in 1 and 2

      yfuel = y(idirac)
      yoxyd = 1.d0 - (coeff3+1.0d0)*f(idirac) + coeff3*y(idirac)
      yprod = 1.d0 - yfuel - yoxyd
      yo2   = coeff1 - (coeff1 + coeff2) * f(idirac)              &
                      + coeff2 * y(idirac)

! ---> absorption coefficients for peaks 1 and 2

!            IF ( IIRAYO .GT. 0  ) THEN
!              KABSGF = YFUEGF(IEL)*KABSG(1) + YOXYGF(IEL)*KABSG(2)
!     &                                      + YPROGF(IEL)*KABSG(3)
!              KABSGB = YFUEGB(IEL)*KABSG(1) + YOXYGB(IEL)*KABSG(2)
!     &                                      + YPROGB(IEL)*KABSG(3)
!            ENDIF


! ---> computation of molar mass and temperature in 1 and 2

      coefg(1) = yfuel
      coefg(2) = yoxyd
      coefg(3) = yprod

! ------ molar mass for peaks 1 and 2

      nbmol = 0.d0
      do igg = 1, ngazg
        nbmol = nbmol + coefg(igg)/wmolg(igg)
      enddo
      maml(idirac) = 1.d0/nbmol

! ------ computation of temperature for peaks 1 and 2

      teml(idirac) = cs_gas_combustion_h_to_t(coefg, h(idirac))

! ---> computation of density in 1 and 2

      if ( ipass.gt.1 .or.                                        &
          (isuite.eq.1.and.initro.eq.1) ) then
        rhol(idirac) = p0 * maml(idirac)                   &
                     / (cs_physical_constants_r*teml(idirac))
      else
        rhol(idirac) = ro0
      endif

! ---> term source computation in 1 and 2 for scalar yfm

      theta(idirac) = ta / teml(idirac)                           &
                    * (1.d0 - teml(idirac) / tstar)

      w(idirac) = vref / lref * (- d(idirac)*rhol(idirac)         &
                  * yfuel*yo2                                     &
                  * exp( -theta(idirac) ))

! ---> mix molar mass

      sum7 = sum7 + d(idirac)*maml(idirac)

! ---> mix temperature

      sum8 = sum8 + d(idirac)*teml(idirac)

! ---> temperature / molar mass

      sum9 = sum9 + d(idirac)*teml(idirac)/maml(idirac)

! ---> global species mass fractions

      sum10 = sum10 + yfuel*d(idirac)

      sum11 = sum11 + yoxyd*d(idirac)

      sum12 = sum12 + yprod*d(idirac)

      sum17 = sum17 +w(idirac)

! ---> property storage in fields

      cpro_ampl(idirac)%p(iel) = d(idirac)
      cpro_fmel(idirac)%p(iel) = f(idirac)
      cpro_fmal(idirac)%p(iel) = y(idirac)
      cpro_maml(idirac)%p(iel) = maml(idirac)
      cpro_teml(idirac)%p(iel) = teml(idirac)
      cpro_rhol(idirac)%p(iel) = rhol(idirac)
      cpro_tscl(idirac)%p(iel) = w(idirac)

! ---> radiative transfer quantities

!           if ( iirayo.gt.0 ) then
!             cpro_ckab(iel) = ygfm*kabsgf + ygbm*kabsgb
!             cpro_t4m(iel)  = ygfm*tgf**4+ygbm*tgb**4
!             cpro_t3m(iel)  = ygfm*tgf**3+ygbm*tgb**3
!           endif

    enddo

    cpro_mam(iel)  = sum7
    cpro_temp(iel) = sum8
    temsmm         = sum9
    cpro_ym1(iel)  = sum10
    cpro_ym2(iel)  = sum11
    cpro_ym3(iel)  = sum12
    cpro_tsc(iel)  = sum17

! ---> mix density

  if (ipass.gt.1.or.(isuite.eq.1.and.initro.eq.1)) then
      crom(iel) = srrom * crom(iel)             &
                       + (1.d0-srrom) * (p0/(cs_physical_constants_r*temsmm))
    endif

  endif

enddo

deallocate(cpro_fmel, cpro_fmal, cpro_teml)
deallocate(cpro_tscl, cpro_rhol, cpro_maml)
deallocate(cpro_ampl)

!----
! FIN
!----

end subroutine
