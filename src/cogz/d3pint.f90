!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

!===============================================================================
! Function:
! ---------

!> \file d3pint.f90
!>
!> \brief Specific physic subroutine: diffusion flame.
!>
!> Integration of thermodynamical variables function of mixture fraction

!  NB : Temperature integration could be ponderated by Cp
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     indpdf        indicator for pdf integration or mean value
!> \param[in]     dirmin        Dirac's peak value at \f$ f_{min} \f$
!> \param[in]     dirmax        Dirac's peak value at \f$ f_{max} \f$
!> \param[in]     fdeb          abscissa of rectangle low boundary
!> \param[in]     ffin          abscissa of rectangle high boundary
!> \param[in]     hrec          rectangle height
!> \param[out]    tpdf          indicator for pdf shape:
!                               - 0: Dirac at mean value
!                               - 1: rectangle
!                               - 2: Dirac's peak at \f$ f_{min} \f$
!                               - 3: Dirac's peak at \f$ f_{max} \f$
!                               - 4: rectangle and 2 Dirac's pics
!> \param[in]     w1            work array
!_______________________________________________________________________________

subroutine d3pint &
 ( indpdf ,                                                       &
   dirmin , dirmax , fdeb   , ffin   , hrec   , tpdf ,            &
   w1      )                                                      &
  bind(C, name='cs_f_d3pint')

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppthch
use coincl
use ppincl
use radiat
use mesh
use field
use pointe

!===============================================================================

implicit none

! Arguments

integer          indpdf(ncelet)
double precision dirmin(ncelet), dirmax(ncelet)
double precision fdeb(ncelet), ffin(ncelet), hrec(ncelet), tpdf(ncelet)
double precision w1(ncelet)


! Local variables

integer          iel, icg
integer          ih, if, jh, jf
double precision aa1, bb1, aa2, bb2, f1, f2, a, b, fmini, fmaxi
double precision u, v, c, d, temsmm, fsir
double precision fm, fp2m

double precision dtsmdf  , dd1df  , dd2df  , df1df  , df2df  , dhrecdf
double precision dtsmdfp2, dd1dfp2, dd2dfp2, df1dfp2, df2dfp2, dhrecdfp2
double precision dtsmdd1, dtsmdd2, dtsmdf1, dtsmdf2, dtsmdhrec, dtsmdhs
double precision dadhs, dbdhs, yprod
double precision, dimension(:), pointer :: cpro_rho, cpro_tsrho
double precision, dimension(:), pointer :: cproaa_rho
double precision, dimension(:), pointer :: cvar_scalt, cpro_tsscalt
double precision, dimension(:), pointer :: cvar_fm, cvar_fp2m
double precision, dimension(:), pointer :: cpro_tsfm, cpro_tsfp2m, cpro_ym3
double precision, dimension(:), pointer :: cpro_temp
double precision, dimension(:), pointer :: cpro_ckabs, cpro_t4m, cpro_t3m
type(pmapper_double_r1), dimension(:), pointer :: cpro_csca

character(len=80) :: th_name, th_st_name

!===============================================================================

integer       ipass
data          ipass /0/
save          ipass

!===============================================================================
! 0. ON COMPTE LES PASSAGES
!===============================================================================

allocate(cpro_csca(ngazg))

do icg = 1, ngazg
  call field_get_val_s(iym(icg), cpro_csca(icg)%p)
enddo

! Initialize variables to avoid compiler warnings

aa1 = 0.d0
aa2 = 0.d0
bb1 = 0.d0
bb2 = 0.d0

ipass = ipass + 1

if (ippmod(icod3p).eq.1) call field_get_val_s(ihm, cvar_scalt)
call field_get_val_s(ifm, cvar_fm)
call field_get_val_s(ifp2m, cvar_fp2m)

!===============================================================================
! 1. INTEGRATION DES NGAZG FRACTIONS MASSIQUES D'ESPECES GLOBALES
!===============================================================================

! ---> En flamme de diffusion chimie 3 points :
!      - il n'y a qu'une seule reaction globale (IR= )
!      - le taux de melange f varie entre 0 et 1
fsir = fs(1)
fmini = zero
fmaxi = 1.d0

do iel = 1, ncel

  fm   = cvar_fm(iel)
  fp2m = cvar_fp2m(iel)

  do icg = 1, ngazg

! --->  Determination des parametres des droites Variables(f)
!       Y = A + B F
!       En flamme de diffusion, il n'y a qu'une seule reaction globale
!         (IR=1)
!       Par definition les fractions massiques des especes globales
!         sont alors

!         F       0  FS  1
!         YFUEL   0  0   1
!         YOXYD   1  0   0
!         YPROD   0  1   0

    if ( icg.eq.1 ) then
!         Fuel
      aa1 =  zero
      bb1 =  zero
      aa2 = -fsir/(1.d0-fsir)
      bb2 =  1.d0/(1.d0-fsir)

    elseif ( icg.eq.2 ) then
!         Oxydant
      aa1 =  1.d0
      bb1 = -1.d0/fsir
      aa2 =  zero
      bb2 =  zero
    elseif ( icg.eq.3 ) then
!         Produits
      aa1 =  zero
      bb1 =  1.d0/fsir
      aa2 =  1.d0/(1.d0-fsir)
      bb2 = -1.d0/(1.d0-fsir)
    endif

    if (indpdf(iel) .eq. 1) then

! ---> Integration de la PDF

      cpro_csca(icg)%p(iel) = dirmin(iel) * ( aa1 + bb1 * fmini )  &
                            + dirmax(iel) * ( aa2 + bb2 * fmaxi )
      if (fdeb(iel).lt.fsir) then
        f1 = fdeb(iel)
        f2 = min( fsir,ffin(iel) )
        cpro_csca(icg)%p(iel) = cpro_csca(icg)%p(iel)                 &
                              + hrec(iel)*(f2-f1)*(aa1+bb1*5.d-1*(f2+f1))
      endif
      if (ffin(iel).gt.fsir) then
        f1 = max(fsir,fdeb(iel))
        f2 = ffin(iel)
        cpro_csca(icg)%p(iel) = cpro_csca(icg)%p(iel)                 &
                              + hrec(iel)*(f2-f1)*(aa2+bb2*5.d-1*(f2+f1))
      endif
    else

! ---> Degenerescence sur la valeur moyenne

      if (fm.le.fsir) then
        cpro_csca(icg)%p(iel) = aa1+bb1*fm
      else
        cpro_csca(icg)%p(iel) = aa2+bb2*fm
      endif

    endif

  enddo

enddo


!===============================================================================
! 2. DETERMINATION LOCALE DE L'ENTHALPIE DES GAZ STOECHIOMETRIQUES
!    BRULES EN PERMEATIQUE (non adiab)
!    (stocke dans le tableau W1)
!===============================================================================

! ---> Calcul de HSTOE dans W1

! ---- Initialisation

do iel = 1, ncel
  w1(iel) = hstoea
enddo

if (ippmod(icod3p).eq.1) then
  call d3phst(ncel, indpdf, dirmin, dirmax, fdeb, ffin, hrec,   &
              cvar_fm, cvar_scalt, w1)
endif

!===============================================================================
! 3. INTEGRATION a) DE LA TEMPERATURE
!                b) DU COEFFICIENT D'ABSORPTION si rayonnement
!                c) DES TERME T^4 et T^3 si rayonnement
!                d) DE LA MASSE VOLUMIQUE
!===============================================================================

! ---> Positions des variables, coefficients

call field_get_val_s(itemp, cpro_temp)
call field_get_val_s(icrom, cpro_rho)

if ( iirayo.gt.0 ) then
  call field_get_val_s(ickabs, cpro_ckabs)
  call field_get_val_s(it4m, cpro_t4m)
  call field_get_val_s(it3m, cpro_t3m)
endif

do iel = 1, ncel

  fm   = cvar_fm(iel)
  fp2m = cvar_fp2m(iel)

  if (indpdf(iel).eq.1) then

! ---> Integration de la PDF

    ih = 1
    do jh = 1,(nmaxh-1)
      if (w1(iel).gt.hh(jh+1) .and. w1(iel).le.hh(jh))        &
           ih = jh
    enddo
    if (w1(iel) .ge. hh(1)) ih = 1
    if (w1(iel) .le. hh(nmaxh)) ih = nmaxh-1
    cpro_temp(iel) = dirmin(iel)*tinoxy + dirmax(iel)*tinfue
    temsmm = dirmin(iel)/wmolg(2)*tinoxy + dirmax(iel)/wmolg(1)*tinfue

    if (iirayo.gt.0) then
      cpro_ckabs(iel) = dirmin(iel)*ckabsg(2)  + dirmax(iel)*ckabsg(1)
      cpro_t4m(iel) =   dirmin(iel)*tinoxy**4 + dirmax(iel)*tinfue**4
      cpro_t3m(iel) =   dirmin(iel)*tinoxy**3 + dirmax(iel)*tinfue**3
    endif
    if = 1
    do jf = 1, (nmaxf-1)
      if (fdeb(iel).ge.ff(jf) .and.                             &
          fdeb(iel).lt.ff(jf+1)) if = jf
    enddo
    if (fdeb(iel) .le. ff(1)) if = 1
    if (fdeb(iel) .ge. ff(nmaxf)) if = nmaxf-1
    f2 = zero
    f1 = fdeb(iel)

    do while ( (ffin(iel)-f2).gt.epzero )
      f2 = min(ff(if+1),ffin(iel))
!          Dans le tableau TFH,
!           on extrait sur chaque ligne i : T = Ai+Bi*F
!           et on construit pour la valeur courante de HSTOE (W1)
!                               T = A+B*F
      aa1 = tfh(if,ih)
      bb1 = (tfh(if+1,ih)-tfh(if,ih))/(ff(if+1)-ff(if))
      aa2 = tfh(if,ih+1)
      bb2 = (tfh(if+1,ih+1)-tfh(if,ih+1))/(ff(if+1)-ff(if))
      a = aa1 + (w1(iel)-hh(ih))*(aa2-aa1)/(hh(ih+1)-hh(ih))
      b = bb1 + (w1(iel)-hh(ih))*(bb2-bb1)/(hh(ih+1)-hh(ih))
      a = a - b*ff(if)

! ----- Calcul de la temperature par integration

      cpro_temp(iel) = cpro_temp(iel) + hrec(iel)*(f2-f1)*(a+b*(f1+f2)/2.d0)

! ----- Preparation aux calculs du coefficient d'absorption
!                               de T^4 et de T^3
!         Cote pauvre
!           UNSMM = (FS-F)/FS / WMOLG(2)+ F/FS / WMOLG(3)
!           CKABS = (FS-F)/FS * CKABSG(2) + F/FS * CKABSG(3)
!         Cote riche
!           UNSMM = (F-FS)/(1-FS)/WMOLG(1) + (1-F)/(1-FS)/WMOLG(3)
!           CKABS = (F-FS)/(1-FS)*CKABSG(1)  + (1-F)/(1-FS)*CKABSG(3)
!         Partout
!           UNSMM = c + df
!           CKABS = u + vF
!           TEMSMM = T*UNSMM = (c+df)*(a+bf) = ca +(cb+ad)f + bd f^2
!           T^4 = (a+bf)^4
!               = a4 + 4a3b f + 6a2b2 f^2 + 4ab3 f^3 + b4 f^4
!           T^3 = (a+bf)^3
!               = a3 + 3a2b f + 3ab2 f^2 + b3 f^3


      if ( f1.lt.fsir ) then
!         On a demarre cote pauvre
        c =   1.d0/wmolg(2)
        d = (-1.d0/wmolg(2)+1.d0/wmolg(3))/fsir
      else
!         On termine cote riche (en commencant avec f1=fs)
        c = (  -fsir/wmolg(1)+1.d0/wmolg(3))/(1.d0-fsir)
        d = (   1.d0/wmolg(1)-1.d0/wmolg(3))/(1.d0-fsir)
      endif

      if ( iirayo.gt.0 ) then
        if ( f1.lt.fsir ) then
!         On a demarre cote pauvre
          u =   ckabsg(2)
          v = (-ckabsg(2)+ ckabsg(3))/fsir
        else
!         On termine cote riche (en commencant avec f1=fs)
          u = (-fsir*ckabsg(1)+ ckabsg(3))/(1.d0-fsir)
          v = (      ckabsg(1)- ckabsg(3))/(1.d0-fsir)
        endif

! ----- Calcul du coefficient d'absorption
!           et des termes T^4 et de T^3 (si rayonnement)

        cpro_ckabs(iel) = cpro_ckabs(iel) +               &
          hrec(iel)*( u*(f2-f1) + v*(f2**2-f1**2)*0.5d0 )

        cpro_t4m(iel) = cpro_t4m(iel) +                             &
                        hrec(iel)*                                  &
                        (      a**4            * (f2-f1)            &
                      +(4.d0*a**3  *b      ) * (f2**2-f1**2)/2.d0   &
                      +(6.d0*(a**2)*(b**2) ) * (f2**3-f1**3)/3.d0   &
                      +(4.d0*a     *(b**3) ) * (f2**4-f1**4)/4.d0   &
                      +(            (b**4) ) * (f2**5-f1**5)/5.d0  )

        cpro_t3m(iel) = cpro_t3m(iel) +                              &
                        hrec(iel)*                                   &
                        (      (a**3)          * (f2-f1)             &
                      +   (3.d0*(a**2)*b      ) * (f2**2-f1**2)/2.d0 &
                      +   (3.d0*a     *(b**2) ) * (f2**3-f1**3)/3.d0 &
                      +   (            (b**3) ) * (f2**4-f1**4)/4.d0  )

      endif

! ----- Calcul du terme Temperature/masse molaire

      temsmm = temsmm + hrec(iel)*                               &
        ( a*c       * (f2-f1)                                     &
        + (c*b+a*d) * (f2**2-f1**2)/2.d0                          &
        +  b*d      * (f2**3-f1**3)/3.d0 )

      if = if+1
      f1 = f2
    enddo

  else

! ---> Degenerescence sur la valeur moyenne

    ih = 1
    do jh = 1, (nmaxh-1)
      if (w1(iel).gt.hh(jh+1) .and. w1(iel).le.hh(jh))        &
            ih = jh
    enddo
    if (w1(iel) .ge. hh(1)) ih =1
    if (w1(iel) .le. hh(nmaxh)) ih =nmaxh-1
    if = 1
    do jf = 1, (nmaxf-1)
      if (fm.ge.ff(jf) .and. fm.lt.ff(jf+1))        &
           if = jf
    enddo
    if (fm .le. ff(1)) if = 1
    if (fm .ge. ff(nmaxf)) if = nmaxf-1
    aa1 = tfh(if,ih)
    bb1 = (tfh(if+1,ih)-tfh(if,ih))/(ff(if+1)-ff(if))
    aa2 = tfh(if,ih+1)
    bb2 = (tfh(if+1,ih+1)-tfh(if,ih+1))/(ff(if+1)-ff(if))
    a  = aa1 + (w1(iel)-hh(ih))*(aa2-aa1)/(hh(ih+1)-hh(ih))
    b  = bb1 + (w1(iel)-hh(ih))*(bb2-bb1)/(hh(ih+1)-hh(ih))
    a  = a - b*ff(if)

! ----- Calcul de la temperature a partir de la valeur moyenne

    cpro_temp(iel) = a+b*fm

    if (fm.lt.fsir) then
!         On a demarre cote pauvre
      c =   1.d0/wmolg(2)
      d = (-1.d0/wmolg(2)+1.d0/wmolg(3))/fsir
    else
!         On termine cote riche (en commencant avec f1=fs)
      c = (  -fsir/wmolg(1)+1.d0/wmolg(3))/(1.d0-fsir)
      d = (   1.d0/wmolg(1)-1.d0/wmolg(3))/(1.d0-fsir)
    endif

    if (iirayo.gt.0) then
      if (fm.lt.fsir) then
!         On a demarre cote pauvre
        u =   ckabsg(2)
        v = (-ckabsg(2)+ ckabsg(3))/fsir
      else
!         On termine cote riche (en commencant avec f1=fs)
        u = (-fsir*ckabsg(1)+ ckabsg(3))/(1.d0-fsir)
        v = (      ckabsg(1)- ckabsg(3))/(1.d0-fsir)
      endif

! ----- Calcul du coefficient d'absorption
!         et des termes T^4 et de T^3
!         a partir de la valeur moyenne (si rayonnement)

      cpro_ckabs(iel) = u + v*fm
      cpro_t4m(iel) = a**4                                     &
                    + (4.d0*(a**3)*b      ) * fm               &
                    + (6.d0*(a**2)*(b**2) ) * fm**2            &
                    + (4.d0*a     *(b**3) ) * fm**3            &
                    + (            (b**4) ) * fm**4

      cpro_t3m(iel) = a**3                                     &
                    + ( 3.d0*(a**2)*b      ) * fm              &
                    + ( 3.d0*a     *(b**2) ) * fm**2           &
                    + (             (b**3) ) * fm**3

    endif

! ----- Calcul du terme Temperature/masse molaire

    temsmm = a*c +(c*b+a*d)*fm + b*d*fm**2

  endif

! ---> Calcul de la masse volumique

  if (ipass.ge.1.or.(isuite.eq.1.and.initro.eq.1)) then
    cpro_rho(iel) = srrom*cpro_rho(iel)               &
                  + (1.d0-srrom)*                         &
                  ( pther/(cs_physical_constants_r*temsmm) )
  endif

enddo

deallocate(cpro_csca)

return
end subroutine
