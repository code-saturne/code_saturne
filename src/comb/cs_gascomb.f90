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

subroutine cs_gascomb &
!====================
 ( ncelet , ncel   , icb1 , icb2 ,                                        &
   indpdf ,                                                               &
   f1m    , f2m    , f3m      , f4m    , f5m  , f6m , f7m , f8m , f9m ,   &
   pdfm1  , pdfm2  , doxyd    , dfuel  , hrec ,                           &
   af1    , af2    , cx1m     , cx2m   , wmf1   , wmf2 ,                  &
   fuel1  , fuel2  , fuel3 , fuel4 , fuel5 ,fuel6 , fuel7  ,              &
   oxyd   , prod1  , prod2  , prod3 , xiner ,                             &
   fs3no  , fs4no  , yfs4no    )

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
! indpdf           ! te ! <-- ! passage par les pdf                            !
! f1m              ! tr ! <-- ! moyenne du traceur 1 mvl [chx1m+co]            !
! f2m              ! tr ! <-- ! moyenne du traceur 2 mvl [chx2m+co]            !
! f3m              ! tr ! <-- ! moyenne du traceur 3 (oxydant 1)               !
! f4m              ! tr ! <-- ! moyenne du traceur 4 (oxydant 2)               !
! f5m              ! tr ! <-- ! moyenne du traceur 5 (oxydant 3)               !
! f6m              ! tr ! <-- ! moyenne du traceur 6 (humidite )               !
! f7m              ! tr ! <-- ! moyenne du traceur 7 (C + O2  )                !
! f8m              ! tr ! <-- ! moyenne du traceur 8 (C + CO2 )                !
! f9m              ! tr ! <-- ! moyenne du traceur 9 (C + H2O )                !
! pdfm1            ! tr ! <-- ! borne minimum de la pdf                        !
! pdfm2            ! tr ! <-- ! borne maximum de la pdf                        !
! dfuel            ! tr ! <-- ! amplitude du pic de dirac en 0                 !
! doxyd            ! tr ! <-- ! amplitude du pic de dirac en 1                 !
! hrec             ! tr ! <-- ! hauteur du rectangle de la pdf                 !
! fuel1            ! tr !  <- ! fraction massique chx1m                        !
! fuel2            ! tr !  <- ! fraction massique chx2m                        !
! fuel3            ! tr !  <- ! fraction massique co                           !
! fuel4            ! tr !  <- ! fraction massique h2s                          !
! fuel5            ! tr !  <- ! fraction massique h2                           !
! fuel6            ! tr !  <- ! fraction massique hcn                          !
! oxyd             ! tr !  <- ! fraction massique o2                           !
! prod1            ! tr !  <- ! fraction massique co2                          !
! prod2            ! tr !  <- ! fraction massique h2o                          !
! prod3            ! tr !  <- ! fraction massique so2                          !
! prod4            ! tr !  <- ! fraction massique nh3                          !
! xiner            ! tr !  <- ! fraction massique n2                           !
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

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel
integer          icb1   , icb2
integer          indpdf(ncel)

double precision f1m(ncel)  , f2m(ncel) , f3m(ncel)
double precision f4m(ncel)  , f5m(ncel) , f6m(ncel)
double precision f7m(ncel)  , f8m(ncel) , f9m(ncel)


double precision pdfm1(ncel), pdfm2(ncel) , dfuel(ncel)
double precision doxyd(ncel) , hrec(ncel)

double precision af1(ncel,ngazg), af2(ncel,ngazg)
double precision cx1m(ncel), cx2m(ncel) , wmf1(ncel) , wmf2(ncel)

double precision fuel1(ncel), fuel2(ncel) , fuel3(ncel)
double precision fuel4(ncel), fuel5(ncel) , fuel6(ncel), fuel7(ncel)
double precision oxyd(ncel)
double precision prod1(ncel), prod2(ncel), prod3(ncel)
double precision xiner(ncel)
!
double precision fs3no(ncel),fs4no(ncel)
double precision yfs4no(ncel,ngazg)

! Local variables
integer           NBPRINT
parameter        (NBPRINT=15)
integer          INTTMP(NBPRINT)
!
integer          ii,iel
integer          n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15

double precision bb1   , bb2  , fs1 , fs2 , fs3 , fs4
double precision anu1  , anu2  , anu3  , anu4  , anu5
double precision reac1 , reac2 , reac3 , reac4 , reac5
double precision scomb , soxy
double precision wmco  , wmh2s , wmh2  , wmhcn
double precision wmo2  , wmco2 , wmh2o , wmso2 , wmnh3 , wmn2
double precision zco2t
double precision somm  , sommax , sommin

double precision zz(ngazem) , zzox(ngazem) , zzcl(ngazem)
double precision zzs1(ngazem),zzs2(ngazem),zzs3(ngazem),zzs4(ngazem)

double precision, dimension(:), pointer :: x1
double precision, dimension(:), pointer :: cvar_yco2

!===============================================================================
! 1. CALCULS PRELIMINAIRES
!===============================================================================

if (ieqco2.eq.1) then
  call field_get_val_s(ivarfl(isca(iyco2)), cvar_yco2)
endif

! Massic fraction of gas
call field_get_val_s_by_name("x_c", x1)

! --- Masses molaires
wmco   = wmole(ico  )
wmh2s  = wmole(ih2s )
wmh2   = wmole(ihy  )
wmhcn  = wmole(ihcn )
wmnh3  = wmole(inh3 )
wmo2   = wmole(io2  )
wmco2  = wmole(ico2 )
wmh2o  = wmole(ih2o )
wmso2  = wmole(iso2 )
wmn2   = wmole(in2  )

! ---- Initialisation des fractions massiques avec l'oxydant 1

do iel = 1, ncel
  fuel1(iel)  = wmf1(iel)*af3(icb1)
  fuel2(iel)  = wmf2(iel)*af3(icb2)
  fuel3(iel)  = wmco     *af3(ico)
  fuel4(iel)  = wmh2s    *af3(ih2s)
  fuel5(iel)  = wmh2     *af3(ihy)
  fuel6(iel)  = wmhcn    *af3(ihcn)
  fuel7(iel)  = wmnh3    *af3(inh3)
  oxyd(iel)   = wmo2     *af3(io2)
  prod1(iel)  = wmco2    *af3(ico2)
  prod2(iel)  = wmh2o    *af3(ih2o)
  prod3(iel)  = wmso2    *af3(iso2)
  xiner(iel)  = wmn2     *af3(in2)
enddo

!===============================================================================
! 3. CALCUL DE LA COMPOSITION DU MELANGE SANS LES PDF
!    SI LES FLUCTUATIONS DE S SONT TROP FAIBLES
!===============================================================================

do iel = 1, ncel

  if ( indpdf(iel).eq.0 ) then

! --> Calculs preliminaires
    if ( ieqco2 .eq. 1 ) then
      zco2t = (cvar_yco2(iel)/ x1(iel))/wmco2
    else
      zco2t = 0.d0
    endif

! --> Composition de la phase gazeuse avant combustion
!
    do ii=1,ngazg
      zz(ii) = af1(iel,ii)*f1m(iel)+af2(iel,ii)*f2m(iel)                &
              +af3(ii)*f3m(iel)+af4(ii)*f4m(iel)+af5(ii)*f5m(iel)       &
              +af6(ii)*f6m(iel)+af7(ii)*f7m(iel)+af8(ii)*f8m(iel)       &
              +af9(ii)*f9m(iel)
    enddo
!
! --> Calcul de la composition du melange
!
!  1ere reaction :
!
    anu1   = 0.5d0
    reac1  = min(zz(ihy),(zz(io2)/anu1))
    zz(ihy) = zz(ihy)  -      reac1
    zz(io2) = zz(io2)  - anu1*reac1
    zz(ih2o)= zz(ih2o) +      reac1
!
!  2eme reaction :
!
    anu2   = 0.25d0*abs(cx1m(iel)-cx2m(iel))
    reac2  = min(zz(icb1), (zz(io2)/anu2))
    zz(icb1)= zz(icb1) -           reac2
    zz(icb2)= zz(icb2) +           reac2
    zz(io2)  = zz(io2)   -      anu2*reac2
    zz(ih2o) = zz(ih2o)  + 2.d0*anu2*reac2
!
!  3eme reaction :
!
    anu3  = 0.25d0*(2.d0 + cx2m(iel))
    reac3 = min(zz(icb2), (zz(io2)/anu3))
    zz(icb2) = zz(icb2) -                 reac3
    zz(ico)   = zz(ico)   +                 reac3
    zz(io2)   = zz(io2)   -            anu3*reac3
    zz(ih2o)  = zz(ih2o)  + 0.5d0*cx2m(iel)*reac3
!
!  4eme reaction :
!
    anu4  = 1.5D0
    reac4 = min(zz(ih2s), (zz(io2)/anu4))
    zz(ih2s) = zz(ih2s)   -      reac4
    zz(io2)  = zz(io2)    - anu4*reac4
    zz(ih2o) = zz(ih2o)   +      reac4
    zz(iso2) = zz(iso2)   +      reac4
!
!  5eme reaction :
!
    anu5   = 0.5d0
    if ( ieqco2 .eq. 0 ) then
      reac5  = min(zz(ico), (zz(io2)/anu5))
    else if ( ieqco2 .eq. 1 ) then
      reac5  = min(max(zco2t-zz(ico2),0.d0),zz(ico),zz(io2)/anu5)
    endif
    zz(ico) = zz(ico)    -      reac5
    zz(io2) = zz(io2)    - anu5*reac5
    zz(ico2)= zz(ico2)   +      reac5
!
    fuel1(iel) = zz(icb1) * wmf1(iel)
    fuel2(iel) = zz(icb2) * wmf2(iel)
    fuel3(iel) = zz(ico ) * wmco
    fuel4(iel) = zz(ih2s) * wmh2s
    fuel5(iel) = zz(ihy ) * wmh2
    fuel6(iel) = zz(ihcn) * wmhcn
    fuel7(iel) = zz(inh3) * wmnh3
    oxyd (iel) = zz(io2 ) * wmo2
    prod1(iel) = zz(ico2) * wmco2
    prod2(iel) = zz(ih2o) * wmh2o
    prod3(iel) = zz(iso2) * wmso2
    xiner(iel) = zz(in2 ) * wmn2
!
  endif

enddo

!===============================================================================
! 4. CALCUL DE LA COMPOSITION DU MELANGE AVEC LES PDF
!===============================================================================
!
do iel = 1, ncel

  if ( indpdf(iel).eq.1 ) then

! --> Calculs preliminaires

    if ( ieqco2 .eq. 1 ) then
      zco2t = (cvar_yco2(iel)/x1(iel))/wmco2
    else
      zco2t = 0.d0
    endif
! Oxydant local
    soxy = f3m(iel)+f4m(iel)+f5m(iel)+f6m(iel)+f7m(iel)+f8m(iel)+f9m(iel)
!
    zzox(icb1) = 0.d0
    zzox(icb2) = 0.d0
    do ii=3,ngazg
      zzox(ii)   = (  af3(ii)*f3m(iel)+af4(ii)*f4m(iel)     &
                     +af5(ii)*f5m(iel)+af6(ii)*f6m(iel)     &
                     +af7(ii)*f7m(iel)+af8(ii)*f8m(iel)     &
                     +af9(ii)*f9m(iel)                        ) / soxy
    enddo
! Combustible local
    scomb       = f1m(iel)+f2m(iel)
    do ii=1,ngazg
      zzcl(ii)  = ( af1(iel,ii)*f1m(iel)+af2(iel,ii)*f2m(iel) ) / scomb
    enddo
!
!  1ere reaction : recombinaison de l'hydrogene
!
    reac1 = min(zzox(ihy),(2.d0*zzox(io2)))
    zzox(ihy ) = zzox(ihy ) -       reac1
    zzox(ih2o) = zzox(ih2o) +       reac1
    zzox(io2 ) = zzox(io2 ) - 0.5d0*reac1
!
!  2eme reaction : CHx1 + (x1-x2)/4 O2  => CHx2 + (x1-x2)/2 H2O
!
    if ( zzcl(icb1) .gt. 0.d0  .and. zzox(io2) .gt. 0.d0 ) then
      fs1 = zzox(io2)/( abs(cx1m(iel)-cx2m(iel))/4.d0      &
           *zzcl(icb1)+zzox(io2) )
    else
      fs1 = 1.d0
    endif
!
    do ii=1,ngazg
      zzs1(ii) = fs1*zzcl(ii)+(1.d0-fs1)*zzox(ii)
    enddo
    zzs1(icb2) = zzs1(icb1)+zzs1(icb2)
    zzs1(ih2o) = zzs1(ih2o) + 0.5d0*abs(cx1m(iel)-cx2m(iel))*zzs1(icb1)
    zzs1(icb1) = 0.d0
    zzs1(io2 ) = 0.d0

!
!  3eme reaction : CHx2 + (2+x2)/4 O2 => CO + x2/2 H2O
!
    if ( zzs1(icb2).gt. 0.d0 .and. zzox(io2).gt. 0.d0 ) then
      fs2 = fs1 * zzox(io2) / ( (2.d0+cx2m(iel))/4.d0   &
           *zzs1(icb2)+zzox(io2) )
    else
      fs2 = fs1
    endif
!
    do ii=1,ngazg
      zzs2(ii) = (fs2/fs1)*zzs1(ii)+(1.d0-(fs2/fs1))*zzox(ii)
    enddo
    zzs2(ico )= zzs2(ico)+zzs2(icb2)
    zzs2(ih2o)= zzs2(ih2o)+cx2m(iel)/2.d0*zzs2(icb2)
    zzs2(icb2)= 0.d0
    zzs2(io2 )= 0.d0
!
!  4eme reaction :  H2S + 3/2 O2 => H2O + SO2
!
    if ( zzs2(ih2s).gt. 0.d0 .and. zzox(io2).gt. 0.d0 ) then
      fs3 = fs2 * zzox(io2) / ( 1.5d0*zzs2(ih2s)+zzox(io2) )
    else
      fs3 = fs2
    endif
!
    do ii=1,ngazg
      zzs3(ii) = (fs3/fs2)*zzs2(ii)+(1.d0-(fs3/fs2))*zzox(ii)
    enddo
    zzs3(iso2) = zzs3(iso2)+zzs3(ih2s)
    zzs3(ih2o) = zzs3(ih2o)+zzs3(ih2s)
    zzs3(ih2s) = 0.d0
    zzs3(io2 ) = 0.d0
!
!  5eme reaction CO+1/2 O2 => CO2
!
    if ( ( zzs3(ico).gt. 0.d0 .and. zzox(io2).gt. 0.d0 )   &
        .and. ieqco2 .eq. 0                              ) then
      fs4 = fs3 * zzox(io2) / ( zzs3(ico) + zzox(io2) )
    else
      fs4 = fs3
    endif
!
    do ii=1,ngazg
      zzs4(ii) = (fs4/fs3)*zzs3(ii)+(1.d0-(fs4/fs3))*zzox(ii)
    enddo
!
    if ( ieqco2 .eq. 0 ) then
      zzs4(ico2) = zzs4(ico2)+zzs4(ico)
      zzs4(io2 ) = zzs4(io2)/2.d0
      zzs4(ico ) = 0.d0
    endif
!
!  Stockage de fs3 , fs4 et des concentrations en fs4 pour le modele de NOx
!
    if ( ieqnox .eq. 1 ) then
!
      fs3no(iel) = fs3
!
      if ( zzs3(ico).gt. 0.d0 .and. zzox(io2).gt. 0.d0 ) then
        fs4no(iel) = fs3 * zzox(io2) / ( zzs3(ico) + zzox(io2) )
      else
        fs4no(iel) = fs3
      endif
!
      do ii=1,ngazg
        yfs4no(iel,ii) = (fs4no(iel)/fs3)*zzs3(ii)              &
                       +(1.d0-(fs4no(iel)/fs3))*zzox(ii)
      enddo
      yfs4no(iel,ico2) = yfs4no(iel,ico2)+yfs4no(iel,ico)
      yfs4no(iel,io2 ) = yfs4no(iel,io2)/2.d0
      yfs4no(iel,ico ) = 0.d0

      yfs4no(iel,icb1) = yfs4no(iel,icb1) * wmf1(iel)
      yfs4no(iel,icb2) = yfs4no(iel,icb2) * wmf2(iel)
      do ii=ico,ngazg
        yfs4no(iel,ii) = yfs4no(iel,ii)*wmole(ii)
      enddo
!
    endif
!
!        Désormais on connait les concentrations e
!        cl,s1,s2,s3,s4,Ox
!        les concentrations intermédiaires sont linéaires par morceaux
!        et les paramêtres de la pdf dfuel,doxyd, pdfm1,pdfm2,hrec
!
!        On initialise par les concentrations aux extrémités
    do ii = 1,ngazg
      zz(ii) = dfuel(iel)*zzcl(ii)+doxyd(iel)*zzox(ii)
    enddo
! Integration sur le premier intervalle de richesse (entre s1 et 1)
    bb1 = max(pdfm1(iel),fs1)
    bb2 = min(pdfm2(iel),1.d0)
    if( bb2.gt.bb1 ) then
      do ii = 1,ngazg
        zz(ii) =  zz(ii) + hrec(iel)*(bb2-bb1)/(1.d0-fs1)       &
                *(  zzs1(ii)-fs1*zzcl(ii)                       &
                  +(zzcl(ii)-zzs1(ii))*(bb1+bb2)*0.5d0 )
      enddo
    endif
! Intégration sur le deuxième intervalle de (entre s2 et s1)
    bb1 = max(pdfm1(iel),fs2)
    bb2 = min(pdfm2(iel),fs1)
    if( bb2.gt.bb1 ) then
      do ii = 1,ngazg
        zz(ii) =  zz(ii) + hrec(iel)*(bb2-bb1)/(fs1-fs2)         &
                *( fs1*zzs2(ii)-fs2*zzs1(ii)                     &
                  +(zzs1(ii)-zzs2(ii))*(bb1+bb2)*0.5d0 )
      enddo
    endif
! Integration sur le troisième intervalle de (entre s3 et s2)
    bb1 = max(pdfm1(iel),fs3)
    bb2 = min(pdfm2(iel),fs2)
    if( bb2.gt.bb1 ) then
      do ii = 1,ngazg
        zz(ii) =  zz(ii) + hrec(iel)*(bb2-bb1)/(fs2-fs3)         &
                *( fs2*zzs3(ii)-fs3*zzs2(ii)                     &
                  +(zzs2(ii)-zzs3(ii))*(bb1+bb2)*0.5d0 )
      enddo
    endif
!
! Integration sur le quatrieme intervalle de (entre s4 et s3)
    bb1 = max(pdfm1(iel),fs4)
    bb2 = min(pdfm2(iel),fs3)
    if ( bb2.gt.bb1 ) then
      do ii = 1,ngazg
        zz(ii) =  zz(ii) + hrec(iel)*(bb2-bb1)/(fs3-fs4)         &
                *( fs3*zzs4(ii)-fs4*zzs3(ii)                     &
                  +(zzs3(ii)-zzs4(ii))*(bb1+bb2)*0.5d0 )
      enddo
    endif
!  Integration sur le cinquieme intervalle de (entre 0 et s4)
    bb1 = max(pdfm1(iel),0.d0)
    bb2 = min(pdfm2(iel),fs4)
    if ( bb2.gt.bb1 ) then
      do ii = 1,ngazg
        zz(ii) =  zz(ii) + hrec(iel)*(bb2-bb1)/(fs4-0.d0)         &
                *( fs4*zzox(ii)-0.d0*zzs4(ii)                     &
                  +(zzs4(ii)-zzox(ii))*(bb1+bb2)*0.5d0 )
      enddo
    endif
!
    if ( ieqco2 .eq. 1 ) then
!
!   transport de co2
!
      anu5   = 0.5d0
      reac5  = min(max(zco2t-zz(ico2),0.d0),zz(ico),zz(io2)/anu5)
!
      zz(ico ) = zz(ico ) -      reac5
      zz(io2 ) = zz(io2 ) - anu5*reac5
      zz(ico2) = zz(ico2) +      reac5
!
    endif

    fuel1(iel) = zz(icb1) * wmf1(iel)
    fuel2(iel) = zz(icb2) * wmf2(iel)
    fuel3(iel) = zz(ico ) * wmco
    fuel4(iel) = zz(ih2s) * wmh2s
    fuel5(iel) = zz(ihy ) * wmh2
    fuel6(iel) = zz(ihcn) * wmhcn
    fuel7(iel) = zz(inh3) * wmnh3
    oxyd(iel)  = zz(io2 ) * wmo2
    prod1(iel) = zz(ico2) * wmco2
    prod2(iel) = zz(ih2o) * wmh2o
    prod3(iel) = zz(iso2) * wmso2
    xiner(iel) = zz(in2 ) * wmn2


  endif

enddo

!===============================================================================
! 5.  IMPRESSION
!===============================================================================
n1  = ncel
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

! --> Controle des parametres de la pdf

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
       + fuel4(iel) + fuel5(iel) + fuel6(iel) + fuel7(iel)        &
       + oxyd(iel)                                                &
       + prod1(iel) + prod2(iel) + prod3(iel)                     &
       + xiner(iel)

  sommin = min(sommin,somm)
  sommax = max(sommax,somm)

  if ( abs(somm-1.d0).lt.epsicp )   then
    n3  = n3 +1
  endif

  if ( fuel1(iel).lt.(-epzero) .or.                               &
       fuel1(iel).gt.(1.d0+epzero) ) n4  = n4 +1
  if ( fuel2(iel).lt.(-epzero) .or.                               &
       fuel2(iel).gt.(1.d0+epzero) ) n5  = n5 +1
  if ( fuel3(iel).lt.(-epzero) .or.                               &
       fuel3(iel).gt.(1.d0+epzero) ) n6  = n6 +1
  if ( fuel4(iel).lt.(-epzero) .or.                               &
       fuel4(iel).gt.(1.d0+epzero) ) n7  = n7 +1
  if ( fuel5(iel).lt.(-epzero) .or.                               &
       fuel5(iel).gt.(1.d0+epzero) ) n8  = n8 +1
  if ( fuel6(iel).lt.(-epzero) .or.                               &
       fuel6(iel).gt.(1.d0+epzero) ) n9  = n9 +1
  if ( fuel7(iel).lt.(-epzero) .or.                               &
       fuel7(iel).gt.(1.d0+epzero) ) n10 = n10+1
  if ( oxyd(iel).lt.(-epzero) .or.                                &
       oxyd(iel).gt.(1.d0+epzero)  ) n11 = n11+1
  if ( xiner(iel).lt.(-epzero) .or.                               &
       xiner(iel).gt.(1.d0+epzero) ) n12 = n12+1
  if ( prod1(iel).lt.(-epzero) .or.                               &
       prod1(iel).gt.(1.d0+epzero) ) n13 = n13+1
  if ( prod2(iel).lt.(-epzero) .or.                               &
       prod2(iel).gt.(1.d0+epzero) ) n14 = n14+1
  if ( prod3(iel).lt.(-epzero) .or.                               &
       prod3(iel).gt.(1.d0+epzero) ) n15 = n15+1

enddo

n1 = ncel

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
if (irangp.ge.0) then
  call parism(nbprint,inttmp)
  !==========
endif

write(nfecra,1000) inttmp(1),inttmp(2)
write(nfecra,2200) (inttmp(ii),ii=3,15)

if ( irangp .ge. 0 ) then
  call parmin(sommin)
  call parmax(sommax)
endif

write(nfecra,*) ' Somme Min MAX = ', sommin, sommax

!--------
! Formats
!--------

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
' Nb de points YCHX1 , YCHX2   < 0 ou > 1 = ',I9,I9,/,      &
' Nb de points YC0   , YH2S    < 0 ou > 1 = ',I9,I9,/,      &
' Nb de points YH2   , YHCN    < 0 ou > 1 = ',I9,I9,/,      &
' Nb de points YNH3            < 0 ou > 1 = ',I9,   /,      &
' Nb de points YO2   , YN2     < 0 ou > 1 = ',I9,I9,/,      &
' Nb de points YCO2  , YH2O    < 0 ou > 1 = ',I9,I9,/,      &
' Nb de points YSO2            < 0 ou > 1 = ',I9  )

!----
! End
!----

return
end subroutine
