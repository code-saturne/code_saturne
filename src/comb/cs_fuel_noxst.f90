!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

subroutine cs_fuel_noxst &
!=======================

 ( ncelet , ncel   ,                                                      &
   indpdf ,                                                               &
   pdfm1  , pdfm2  , doxyd  , dfuel  , hrec ,                             &
   f3m    , f4m    , f5m    , f6m    , f7m  , f8m , f9m ,                 &
   fs3no  , fs4no  , yfs4no , enthox )

!===============================================================================
! FONCTION :
! --------
! CALCUL DE :
!
!           K1 exp(-E1/RT)   (conversion HCN en N2)
!           K2 exp(-E2/RT)   (conversion HCN en NO)
!           K3 exp(-E3/RT)   (NO thermique)
!
!----------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! indpdf           ! te ! <-- ! passage par les pdf                            !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHAMNUMERIQUE), T (TABLEAU)
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
use parall
use pointe
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use field
use cs_fuel_incl

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel
integer          indpdf(ncel)

double precision pdfm1(ncel) , pdfm2(ncel) , dfuel(ncel)
double precision f3m(ncel)   , f4m(ncel)   , f5m(ncel)
double precision f6m(ncel)   , f7m(ncel)   , f8m(ncel) ,f9m(ncel)
double precision doxyd(ncel) , hrec(ncel)
double precision fs3no(ncel) , fs4no(ncel) , yfs4no(ncel,ngazg)
double precision enthox(ncel)

!
! VARIABLES LOCALES
!
integer          iel , icla , i
integer          ipdf1 , ipdf2 , ipdf3
integer          npart
parameter        (npart = 200 )

double precision ee1,ee2,ee3,kk1,kk2,kk3
double precision wmo2,tg,xo2,bb,xmx2
double precision bb1 , bb2 , bb3 , bb4
double precision dirac , tfuel , tmpgaz , lrf ,lro
double precision qqq , rrr  , sss , ttt , uuu
double precision gt1,gt2,gt3,gt10,gt20,gt30
double precision ts4,ts4num,ts4den
double precision dgs
double precision val(npart+1),tt(npart+1) , gs(npart+1) , yyo2(npart+1)
!
integer          icha,ige,mode
double precision hhf , hfuel,tfs4ad,hfs4ad
double precision aa
double precision coefe(ngazem),f1mc(ncharm),f2mc(ncharm)
double precision yo2moy,yo2ox,yo2cb,yo24den,yo24num,yo2s4,yo2
double precision toxyd,hoxyd , som , react , deltamol
!
integer          i300 ,i000 ,imini,i2500,i2600,i2700,i2800,i3000,i3500
integer          i4000,i5000,imaxi,inok
integer          nbpt , nbclip1 , nbclip2
integer          nbclip30,nbclip31
double precision ts4min,ts4max
double precision somm  , sommax , sommin , ts4admin,ts4admax
double precision yo2min,yo2max,yo2min1,yo2max1
double precision yo2oxmin,yo2oxmax
double precision toxmin,toxmax

double precision, dimension(:), pointer :: cpro_exp1, cpro_exp2, cpro_exp3
double precision, dimension(:), pointer :: cpro_temp1, cpro_yo2, cpro_mmel
type(pmapper_double_r1), dimension(:), allocatable :: cvar_yfolcl, cvar_h2cl
type(pmapper_double_r1), dimension(:), allocatable :: cpro_temp2

!
!
!===============================================================================
! 1. CALCULS PRELIMINAIRES
!===============================================================================

! --- Masses molaires
!
wmo2   = wmole(io2  )
!
! --- pointeurs
!
call field_get_val_s(ighcn1,cpro_exp1)
call field_get_val_s(ighcn2,cpro_exp2)
call field_get_val_s(ignoth,cpro_exp3)
!
! Parametres des lois d'Arrhenius
!
kk1 = 3.0e12
ee1 = 3.0e4
kk2 = 1.2e10
ee2 = 3.35e4
kk3 = 3.4e12
ee3 = 6.69e4
!
! Pour les termes, indicateur de calcul par PDF ou non
!       = 1 --> passage par pdf
!       = 0 --> on ne passe pas par les pdf
!
ipdf1 = 0
ipdf2 = 0
ipdf3 = 1
!
! Initialisation
!
do iel = 1, ncel
  cpro_exp1(iel)  = zero
  cpro_exp2(iel)  = zero
  cpro_exp3(iel)  = zero
enddo

!===============================================================================
! 2. CALCUL SANS LES PDF
!===============================================================================

call field_get_val_s(itemp1,cpro_temp1)
call field_get_val_s(iym1(io2),cpro_yo2)
call field_get_val_s(immel,cpro_mmel)

do iel = 1, ncel
!
  tg  = cpro_temp1(iel)
  yo2 = cpro_yo2(iel)
  xo2 = yo2*cpro_mmel(iel)/wmo2
!
! Reaction HCN + NO + 1/4 O2 ---> N2 + 1/2 H2O + CO
!
  cpro_exp1(iel)  = kk1*exp(-ee1/tg)
!
! Reaction HCN + 5/4 O2 --> NO + 1/2 H2O  + CO
!
  if ( xo2 .gt. 0.018d0 ) then
    bb = 0.d0
  else if ( xo2 .lt. 0.0025d0 ) then
    bb= 1.d0
  else
    bb=(0.018d0-xo2)/(0.018d0-0.0025d0)
  endif
  cpro_exp2(iel)  = kk2*exp(-ee2/tg)*(xo2**bb)
!
! No thermique : Zeldovich
!
  cpro_exp3(iel)  = kk3*exp(-ee3/tg)*(xo2**0.5d0)
!
enddo
!
!===============================================================================
! 3. CALCUL AVEC LES PDF
!===============================================================================
!
if ( ipdf1 .eq. 1 .or. ipdf2 .eq. 1 .or. ipdf3 .eq. 1 ) then

  ! Arrays of pointers containing the fields values for each class
  ! (loop on cells outside loop on classes)
  allocate(cvar_yfolcl(nclafu), cvar_h2cl(nclafu))
  allocate(cpro_temp2(nclafu))
  do icla = 1, nclafu
    call field_get_val_s(ivarfl(isca(iyfol(icla))), cvar_yfolcl(icla)%p)
    call field_get_val_s(ivarfl(isca(ih2(icla))), cvar_h2cl(icla)%p)
    call field_get_val_s(itemp2(icla), cpro_temp2(icla)%p)
  enddo

!
  inok = 0
  i300 = 0
  i000 = 0
  imini= 0
  i2500= 0
  i2600= 0
  i2700= 0
  i2800= 0
  i3000= 0
  i3500= 0
  i4000= 0
  i5000= 0
  imaxi= 0
  nbpt = 0
  nbclip1 =0
  nbclip2 = 0
  ts4min= 1.d+20
  ts4max=-1.d+20
  sommin= 1.d+20
  sommax=-1.d+20
  ts4admin= 1.d+20
  ts4admax=-1.d+20
  toxmin = 1.d+20
  toxmax =-1.d+20
  yo2oxmin= 1.d+20
  yo2oxmax=-1.d20
!
  nbclip30 = 0
  nbclip31 = 0
  yo2min  = 1.d+20
  yo2max  =-1.d+20
  yo2min1 = 1.d+20
  yo2max1 =-1.d+20
!
  do iel=1,ncel
!
    if ( indpdf(iel) .eq. 1 ) then
!
!  Calcul de Yo2 dans l'oxydant
!            Yo2 en fs4
!
      bb1 = max(0.d0      ,pdfm1(iel))
      bb2 = min(fs3no(iel),pdfm2(iel))
      if ( bb2 .gt. bb1 ) then
        lro = hrec(iel)
      else
        lro = 0.d0
      endif
      qqq = bb2**2 - bb1**2
      rrr = bb2    - bb1
      gt1 = lro*qqq/(2.d0*fs3no(iel))
      gt2 = lro*rrr
      gt3 = doxyd(iel)
      yo2cb = 0.d0

      if (  cpro_yo2(iel) .gt. 0.d0 ) then
        yo2ox = cpro_yo2(iel)/(-gt1+gt2+gt3)
!
        yo2oxmin = min(yo2oxmin,yo2ox)
        yo2oxmax = max(yo2oxmax,yo2ox)
!
        yo2moy = cpro_yo2(iel)
        dirac  =  dfuel(iel)*yo2cb + doxyd(iel)*yo2ox
!
        bb1 = max(0.D0      ,pdfm1(iel))
        bb2 = min(fs4no(iel),pdfm2(iel))
        bb3 = max(fs4no(iel),pdfm1(iel))
        bb4 = min(fs3no(iel),pdfm2(iel))
        if ( bb2 .gt. bb1 ) then
          lro = hrec(iel)
        else
          lro = 0.d0
        endif
        if ( bb4 .gt. bb3 ) then
          lrf = hrec(iel)
        else
          lrf = 0.d0
        endif
!
        qqq = bb2**2 - bb1**2
        rrr = bb2    - bb1
        sss = bb4**2 - bb3**2
        ttt = bb4    - bb3
        uuu = fs4no(iel)-fs3no(iel)
!
        gt1 = lro*qqq/(2.d0*fs4no(iel))
        gt2 = lro*rrr
!
        gt10= lrf*sss/(2.d0*uuu)
        gt20= lrf*ttt*fs3no(iel)/uuu
!
        yo24num = yo2moy - dirac + yo2ox*(gt1 -gt2)
        yo24den = gt1+gt10-gt20
!
        yo2s4  = yo24num/yo24den
!
      else
        yo2ox = 0.d0
        yo2s4 = 0.d0
      endif
!
!     Affichage et clipping
!
      yo2min = min(yo2s4,yo2min)
      yo2max = max(yo2s4,yo2max)
      if ( yo2s4 .lt. 0.d0 ) then
         nbclip30 = nbclip30+1
         yo2s4 = 0.d0
      endif
      if ( yo2s4 .gt. yo2ox ) then
         nbclip31 = nbclip31+1
         yo2s4 = yo2ox
      endif
      yo2min1 = min(yo2s4,yo2min1)
      yo2max1 = max(yo2s4,yo2max1)

!
!     Calcul de Tfioul moyen
!     ----------------------
!
      tfuel = 0.d0
      xmx2  = 0.d0
!
      do icla = 1, nclafu
        xmx2 = xmx2 + cvar_yfolcl(icla)%p(iel)
      enddo
      if ( xmx2 .gt. 0.d0 ) then
        do icla=1,nclafu
          tfuel = tfuel +cvar_yfolcl(icla)%p(iel)*cpro_temp2(icla)%p(iel)
        enddo
!
        tfuel = tfuel/xmx2
!
      else
        tfuel = cpro_temp1(iel)
      endif
!
!    On recupere la valeur de Toxyd a partir de hoxyd
!
      hoxyd = enthox(iel)
!
      do ige = 1, ngazem
        coefe(ige) = zero
      enddo
      coefe(io2) =( af3(io2 )*f3m(iel)+af4(io2 )*f4m(iel)+af5(io2 )*f5m(iel) &
                   +af6(io2 )*f6m(iel)+af7(io2 )*f7m(iel)+af8(io2 )*f8m(iel) &
                   +af9(io2 )*f9m(iel) ) * wmole(io2)

      coefe(in2) =( af3(in2 )*f3m(iel)+af4(in2 )*f4m(iel)+af5(in2 )*f5m(iel) &
                   +af6(in2 )*f6m(iel)+af7(in2 )*f7m(iel)+af8(in2 )*f8m(iel) &
                   +af9(in2 )*f9m(iel) ) * wmole(in2)

      coefe(ico2)=( af3(ico2)*f3m(iel)+af4(ico2)*f4m(iel)+af5(ico2)*f5m(iel) &
                   +af6(ico2)*f6m(iel)+af7(ico2)*f7m(iel)+af8(ico2)*f8m(iel) &
                   +af9(ico2)*f9m(iel) ) * wmole(ico2)

      coefe(ih2o)=( af3(ih2o)*f3m(iel)+af4(ih2o)*f4m(iel)+af5(ih2o)*f5m(iel) &
                   +af6(ih2o)*f6m(iel)+af7(ih2o)*f7m(iel)+af8(ih2o)*f8m(iel) &
                   +af9(ih2o)*f9m(iel) ) * wmole(ih2o)

      coefe(ico)=(  af3(ico)*f3m(iel)+af4(ico)*f4m(iel)+af5(ico)*f5m(iel)    &
                   +af6(ico)*f6m(iel)+af7(ico)*f7m(iel)+af8(ico)*f8m(iel)    &
                   +af9(ico)*f9m(iel) ) * wmole(ico)

      coefe(ihy)=(  af3(ihy)*f3m(iel)+af4(ihy)*f4m(iel)+af5(ihy)*f5m(iel)    &
                   +af6(ihy)*f6m(iel)+af7(ihy)*f7m(iel)+af8(ihy)*f8m(iel)    &
                   +af9(ihy)*f9m(iel) ) * wmole(ihy)
!
!   Nbre de mole qui a reagis avec CO pour faire du CO2
!
      deltamol   = (coefe(io2)-yo2ox)/wmole(io2)
      react      = min(2.d0*deltamol,coefe(ico)/wmole(ico))
      coefe(io2 )= coefe(io2 )-react/2.d0*wmole(io2)
      coefe(ico )= coefe(ico )-react     *wmole(ico)
      coefe(ico2)= coefe(ico2)+react     *wmole(ico2)

!
      som = coefe(io2)+coefe(in2)+coefe(ico2)     &
           +coefe(ih2o)+coefe(ico)+coefe(ihy)
      coefe(io2 ) = coefe(io2 ) /som
      coefe(in2 ) = coefe(in2 ) /som
      coefe(ico2) = coefe(ico2) /som
      coefe(ih2o) = coefe(ih2o) /som
      coefe(ico ) = coefe(ico ) /som
      coefe(ihy ) = coefe(ihy ) /som
!
      do icha = 1, ncharm
        f1mc(icha) = zero
        f2mc(icha) = zero
      enddo
!
      mode = 1
      call cs_fuel_htconvers1 ( mode , hoxyd , coefe , toxyd )
      !======================
!
      toxmin = min(toxmin,toxyd)
      toxmax = max(toxmax,toxyd)
!
      if ( toxyd .gt. cpro_temp1(iel) ) then
        toxyd = cpro_temp1(iel)
      endif
!
!    On initialise par les temperatures Toxy et Tfuel aux extremites
!
      dirac  =  dfuel(iel)*tfuel + doxyd(iel)*toxyd
!
!    On recupere la valeur de la temperature moyenne
!
      tmpgaz = cpro_temp1(iel)
!
      bb1 = max(0.D0      ,pdfm1(iel))
      bb2 = min(fs4no(iel),pdfm2(iel))
      bb3 = max(fs4no(iel),pdfm1(iel))
      bb4 = min(1.d0      ,pdfm2(iel))
!
      if ( bb2 .gt. bb1 ) then
        lro = hrec(iel)
      else
        lro = 0.d0
      endif
      if ( bb4 .gt. bb3 ) then
        lrf = hrec(iel)
      else
        lrf = 0.d0
      endif
!
      qqq = bb2**2 - bb1**2
      rrr = bb2    - bb1
      sss = bb4**2 - bb3**2
      ttt = bb4    - bb3
      uuu = 1.d0   - fs4no(iel)
!
      gt1 = lro*qqq/(2.d0*fs4no(iel))
      gt2 = lro*rrr
!
      gt10= lrf*sss/(2.d0*uuu)
      gt20= lrf*ttt
      gt30= lrf*ttt/uuu
!
      ts4num = tmpgaz - dirac + toxyd*(gt1 -gt2      )    &
                              - tfuel*(gt10+gt20-gt30)
      ts4den = gt1-gt10+gt30
!
      ts4 = ts4num/ts4den
!
!   Calcul de l'enthalpie des vapeur a Tfuel
!
      hhf = 0.D0
      do icla=1,nclafu
        hhf = hhf + cvar_h2cl(icla)%p(iel)
      enddo
!
      if ( xmx2 .gt. epsifl ) then
!
        hfuel = hhf / xmx2
!
        hfs4ad = fs4no(iel)*hfuel+(1.d0-fs4no(iel))*hoxyd
!
        do ige = 1, ngazem
          coefe(ige) = zero
        enddo
        do ige = 1, ngazg
          coefe(ige) = yfs4no(iel,ige)
        enddo
!
        mode = 1
        call cs_fuel_htconvers1 ( mode , hfs4ad , coefe , tfs4ad )
        !======================

!
! Calcul pour affichage
!
        ts4min = min(ts4min,ts4)
        ts4max = max(ts4max,ts4)
!
        ts4admin= min(ts4admin,tfs4ad)
        ts4admax= max(ts4admax,tfs4ad)
!
        somm = 0.d0
        do ige = 1, ngazg
          somm = somm + yfs4no(iel,ige)
        enddo
        sommin=min(sommin,somm)
        sommax=max(sommax,somm)
!
        if ( (ts4 .gt. min(toxyd,tmpgaz,tfuel)) .and. &
             (ts4 .lt. 2400.d0) ) inok = inok + 1
        if ( ts4 .lt. min(toxyd,tmpgaz,tfuel) ) then
          if ( ts4 .ge. 300.d0 ) then
            i300  = i300  + 1
          else if ( ts4 .gt. 0 ) then
            i000  = i000  + 1
          else
            imini = imini + 1
          endif
        endif
!
        if ( ts4 .gt. 2400.d0 ) then
          if ( ts4 .lt. 2500 ) then
            i2500 = i2500 +1
          else if ( ts4 .lt. 2600 ) then
            i2600 = i2600 +1
          else if ( ts4 .lt. 2700 ) then
            i2700 = i2700 +1
          else if ( ts4 .lt. 2800 ) then
            i2800 = i2800 +1
          else if ( ts4 .lt. 3000 ) then
            i3000 = i3000 +1
          else if ( ts4 .lt. 3500 ) then
            i3500 = i3500 +1
          else if ( ts4 .lt. 4000 ) then
            i4000 = i4000 +1
          else if ( ts4 .lt. 5000 ) then
            i5000 = i5000 +1
          else
            imaxi = imaxi + 1
          endif
        endif
!
! Fin calcul pour affichage
!
!
! Clipping de Ts4 : a min(toxyd,tmpgaz,tfuel) en min
!                   a ts4ad                   en max
!
        nbpt = nbpt + 1
        if ( ts4 .lt. min(toxyd,tmpgaz,tfuel) ) then
           nbclip1 = nbclip1 + 1
           ts4 = min(toxyd,tmpgaz,tfuel)
        endif
!
        if ( ts4 .gt. tfs4ad ) then
           nbclip2 = nbclip2 + 1
           ts4 = tfs4ad
        endif
!
!   Concentration oxygene
!
        xo2 = cpro_yo2(iel)                 &
             *cpro_mmel(iel)/wmole(io2)
!
!  Integration
!
        do i = 1, npart+1
          gs(i) = pdfm1(iel)+dble(i-1)/dble(npart)*(pdfm2(iel)-pdfm1(iel))
!        calcul de T
          if( gs(i) .lt. fs4no(iel) ) then
            tt(i) = (ts4-toxyd)/fs4no(iel)* gs(i) + toxyd
          else
            tt(i) = (tfuel-ts4)/(1.d0-fs4no(iel))*gs(i)       &
                   + tfuel - (tfuel-ts4)/(1.d0-fs4no(iel))
          endif
!        calul de yo2
          if ( gs(i) .lt. fs4no(iel) )  then
            yyo2(i) = (yo2s4-yo2ox)/fs4no(iel) * gs(i) + yo2ox
          else if ( gs(i) .lt. fs3no(iel) ) then
            aa = yo2s4/(fs4no(iel)-fs3no(iel))
            yyo2(i) = aa * ( gs(i) -fs3no(iel) )
          else
            yyo2(i) = 0.d0
          endif
        enddo
!
!     pas d'integration
!
        dgs = ( pdfm2(iel) - pdfm1(iel) ) / dble(npart)
!
! Calcul de K1*EXP(-E1/T)
!
        if ( ipdf1 .eq. 1 ) then
!
          cpro_exp1(iel) = kk1*exp(-ee1/toxyd)*doxyd(iel)        &
                            +kk1*exp(-ee1/tfuel)*dfuel(iel)

          do i = 1, npart+1
            val(i) = kk1*exp(-ee1/tt(i))*hrec(iel)
          enddo

          do i = 1, npart
            cpro_exp1(iel) = cpro_exp1(iel)                 &
                               +0.5d0*dgs*(val(i)+val(i+1))
          enddo

        endif
!
!  Calcul de K2*EXP(-E2/T)
!
          if ( ipdf2 .eq. 1 ) then
!
            if ( xo2 .gt. 0.d0 ) then
!
              if(xo2.gt.0.018d0) then
                bb=0.d0
              else if(xo2 .lt. 0.0025d0) then
                bb=1.d0
              else
                bb=(0.018d0-xo2)/(0.018d0-0.0025d0)
              endif
!
              cpro_exp2(iel) = kk2*exp(-ee2/toxyd)*doxyd(iel)*(xo2**bb) &
                                 +kk2*exp(-ee2/tfuel)*dfuel(iel)*(xo2**bb)
!
              do i = 1, npart+1
                val(i) = kk2*exp(-ee2/tt(i))*hrec(iel)
              enddo

              do i = 1, npart
                cpro_exp2(iel) = cpro_exp2(iel)                      &
                                   +0.5d0*dgs*(val(i)+val(i+1))*(xo2**bb)
              enddo
            else
              cpro_exp2(iel) = 0.d0
            endif

          endif
!
!  Calcul de K3*EXP(-E3/T)
!
          if ( ipdf3 .eq. 1 ) then

            if ( xo2 .gt. 0.d0 ) then
!
              cpro_exp3(iel) = kk3*exp(-ee3/toxyd)                  &
                                         *doxyd(iel)*(yo2ox**0.5d0)    &
                                 +kk3*exp(-ee3/tfuel)                  &
                                         *dfuel(iel)*(yo2cb**0.5d0)
!
              do i = 1, npart+1
                if (yyo2(i).gt.0.d0) then
                  if (gs(i).le.fs3no(iel)) then
                    val(i) = kk3*exp(-ee3/tt(i))*hrec(iel)*(yyo2(i)**0.5d0)
                  else
                    val(i) = 0.d0
                  endif
                else
                  val(i) = 0.d0
                endif
              enddo

              do i = 1, npart
                cpro_exp3(iel) = cpro_exp3(iel)                  &
                                   +0.5d0*dgs*(val(i)+val(i+1))
              enddo
!
            else
              cpro_exp3(iel) = 0.d0
            endif
          endif
!
      endif
!
    endif
  enddo

  deallocate(cvar_yfolcl, cvar_h2cl)
  deallocate(cpro_temp2)

!
  if ( irangp .ge. 0 ) then
    call parcpt(inok)
    call parcpt(i300)
    call parcpt(i000)
    call parcpt(imini)
    call parcpt(i2500)
    call parcpt(i2600)
    call parcpt(i2700)
    call parcpt(i2800)
    call parcpt(i3000)
    call parcpt(i3500)
    call parcpt(i4000)
    call parcpt(i5000)
    call parcpt(imaxi)
    call parcpt(nbpt)
    call parcpt(nbclip1)
    call parcpt(nbclip2)
    call parmin(ts4min)
    call parmax(ts4max)
    call parmin(ts4admin)
    call parmax(ts4admax)
    call parmin(sommin)
    call parmax(sommax)
!
    call parcpt(nbclip30)
    call parcpt(nbclip31)
    call parmin(yo2min)
    call parmax(yo2max)
    call parmin(yo2min1)
    call parmax(yo2max1)
!
    call parmin(toxmin)
    call parmax(toxmax)
    call parmin(yo2oxmin)
    call parmax(yo2oxmax)
  endif
!===============================================================================
  write(nfecra,*) ' '
  write(nfecra,*) ' Min max de TSox ',toxmin,toxmax
  write(nfecra,*) ' Min max de TS4 ',ts4min,ts4max
  write(nfecra,*) '    nombre de pts a Tmini < ts4 < 2400  ',inok
  write(nfecra,*) '    nombre de pts a         ts4 < 0     ',imini
  write(nfecra,*) '    nombre de pts a 0     < ts4 < 300   ',i000
  write(nfecra,*) '    nombre de pts a 300   < ts4 < Tmini ',i300
  write(nfecra,*) '    nombre de pts a 2400  < ts4 < 2500  ',i2500
  write(nfecra,*) '    nombre de pts a 2500  < ts4 < 2600  ',i2600
  write(nfecra,*) '    nombre de pts a 2600  < ts4 < 2700  ',i2700
  write(nfecra,*) '    nombre de pts a 2700  < ts4 < 2800  ',i2800
  write(nfecra,*) '    nombre de pts a 2800  < ts4 < 3000  ',i3000
  write(nfecra,*) '    nombre de pts a 3000  < ts4 < 3500  ',i3500
  write(nfecra,*) '    nombre de pts a 3500  < ts4 < 4000  ',i4000
  write(nfecra,*) '    nombre de pts a 4000  < ts4 < 5000  ',i5000
  write(nfecra,*) '    nombre de pts a 5000  < ts4         ',imaxi
  write(nfecra,*) ' Min max de TS4ad ',ts4admin,ts4admax
  write(nfecra,*) ' Min max concentration en fs4 ',sommin,sommax
  write(nfecra,*) ' Clipping en min en ',nbclip1,' points sur ',nbpt,' points'
  write(nfecra,*) ' Clipping en max en ',nbclip2,' points sur ',nbpt,' points'
!
  write(nfecra,*) ' '
  write(nfecra,*) ' Min max de Yo2ox en 0 ',yo2oxmin,yo2oxmax
  write(nfecra,*) ' Min max de Yo2 en fs4 avant clipping ',yo2min,yo2max
  write(nfecra,*) ' Clipping en min sur Yo2 en fs4       ',nbclip30
  write(nfecra,*) ' Clipping en max sur Yo2 en fs4       ',nbclip31
  write(nfecra,*) ' Min max de Yo2 en fs4 apres clipping ',yo2min1,yo2max1
!===============================================================================
!
endif

return
end subroutine
