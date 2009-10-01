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

subroutine fucyno &
!================

 ( ncelet , ncel   ,                                              &
   indpdf ,                                                       &
   rtp    , propce ,                                              &
   f1m    , f3m    , f4m   , f1cl  , f3cl , f4cl ,                &
   f4m1   , f4m2   , d4cl  , d4f4  , hrec , f4s3 ,                &
   fuel1  , fuel3  , oxyd  , prod1 , prod2  ,                     &
   xiner  , xh2s   , xso2   )

!===============================================================================
! FONCTION :
! --------

! CALCUL DE :

!          K1 exp(-E1/RT)   (conversion HCN en NO)
!          K2 exp(-E2/RT)   (conversion HCN en N2)
!          K3 exp(-E3/RT)   (NO thermique)

!  VALEURS CELLULES
!  ----------------

! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant)                  !
! fvap             ! tr ! <-- ! moyenne du traceur 1 mvl [fovm+co]             !
! fhet             ! tr !  ->           ! moyenne du traceur 3 c hÈterogËne    !
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
!     DONNEES EN COMMON
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

double precision rtp(ncelet,*), propce(ncelet,*)
double precision f1m(ncelet) , f3m(ncelet)  , f4m(ncelet)
double precision f1cl(ncelet), f3cl(ncelet) , f4cl(ncelet)

double precision f4m1(ncelet) , f4m2(ncelet) , d4cl(ncelet)
double precision d4f4(ncelet) , hrec(ncelet) , f4s3(ncel)

double precision fuel1(ncelet), fuel3(ncelet)
double precision oxyd(ncelet) , prod1(ncelet), prod2(ncelet)
double precision xiner(ncelet)
double precision xh2s(ncelet) , xso2(ncelet)

! VARIABLES LOCALES

integer          iel , icla , i

integer          n1 , n2 , n3 , n4 , n5 , n6 , n7 , n8 , n9
integer          n10, n11
integer          nbrint
integer          ipdf1 , ipdf2 , ipdf3
integer          iexp1 , iexp2 , iexp3
parameter        (nbrint = 11)
integer          inttmp(nbrint)
integer          npart
parameter        (npart = 200 )

double precision wmo2

double precision ee1,ee2,ee3,kk1,kk2,kk3
double precision ts3 , ts3den , ts3num , u , t , ts3s , dirac
double precision q , r , s , lpa , lri , tmpgaz , tfuel , xmx2
double precision bb1 , bb2 , bb3 , bb4 , wmco , xo2 , bb
double precision gh1, gh2 ,gh3 ,gh4 ,gh5,gh6, gh7 , tg
double precision gh8, gh9 ,gh10,gh11
double precision dgs , rsf4 , rscl , bmoy , p , spdf
double precision val(npart+1),tt(npart+1) , gs(npart+1)

!===============================================================================


!===============================================================================
! 0. INITIALISATIONS
!===============================================================================


!===============================================================================
! 1. CALCULS PRELIMINAIRES
!===============================================================================

! --- Masses molaires

wmo2   = wmole(io2)

! --- pointeurs

iexp1 = ipproc(ighcn1)
iexp2 = ipproc(ighcn2)
iexp3 = ipproc(ignoth)

! Parametres des lois d'Arrhenius

kk1 = 3.0e12
ee1 = 3.0e4
kk2 = 1.2e10
ee2 = 3.35e4
kk3 = 3.4e12
ee3 = 6.69e4

! Pour les termes, indicateur de calcul par PDF ou non
!       = 1 --> passage par pdf
!       = 0 --> on ne passe pas par les pdf

ipdf1 = 0
ipdf2 = 0
ipdf3 = 0

! ---- Initialisation

do iel = 1, ncel
!        PROPCE(IEL,IEXP1)  = ZERO
!        PROPCE(IEL,IEXP2)  = ZERO
!        PROPCE(IEL,IEXP3)  = ZERO
enddo

!===============================================================================
! 3. CALCUL SANS LES PDF
!===============================================================================

do iel = 1, ncel

  tg  = propce(iel,ipproc(itemp1))
  propce(iel,iexp1)  = kk1*exp(-ee1/tg)

enddo

do iel = 1, ncel

  tg  = propce(iel,ipproc(itemp1))
  xo2 = oxyd(iel)*propce(iel,ipproc(immel))                       &
       /wmo2
  if ( xo2 .gt. 0.018d0 ) then
    bb = 0.d0
  else if ( xo2 .lt. 0.0025d0 ) then
    bb= 1.d0
  else
    bb=(0.018d0-xo2)/(0.018d0-0.0025d0)
  endif

  propce(iel,iexp2)  = kk2*exp(-ee2/tg)*(xo2**bb)

enddo

do iel = 1, ncel

  tg  = propce(iel,ipproc(itemp1))
  xo2 = oxyd(iel)*propce(iel,ipproc(immel))                       &
       /wmo2

  propce(iel,iexp3)  = kk3*exp(-ee3/tg)*(xo2**0.5d0)

enddo

!===============================================================================
! 3. CALCUL AVEC LES PDF
!===============================================================================

if ( ipdf1 .eq. 1 .or. ipdf2 .eq. 1 .or. ipdf3 .eq. 1 ) then

  do iel=1,ncel

    if ( indpdf(iel).eq.1 ) then

!         Calcul de Tfioul moyen

      xmx2  = 0.d0
      tfuel = 0.d0

      do icla=1,nclafu
        xmx2 = xmx2 + rtp(iel,isca(iyfol(icla)))
      enddo

      if ( xmx2 .gt. 0.d0 ) then
        do icla=1,nclafu
          tfuel = tfuel + rtp(iel,isca(iyfol(icla)))              &
                         *propce(iel,ipproc(itemp3(icla)))
        enddo
        tfuel = tfuel/xmx2
      else
        tfuel = propce(iel,ipproc(itemp1))
      endif

!  On recupere la valeur de TAIR

      taire = rtp(iel,isca(itaire))

! On initialise par les temperatures Tair (en f4) et Tcl (TFUEL) aux extrÅÈmitÅÈs

      dirac  = d4cl(iel)*tfuel + d4f4(iel)*taire

!1 On recupere la valeur de la temperature moyenne

      tmpgaz = propce(iel,ipproc(itemp1))

!1 Integration premier intervalle, entre CL et S3 (stoechio derniÅËre rÅÈaction)

!1 Parametres de la droite entre CL et S3

!1 Detection des bornes d'integration (celles de la pdf) du domaine pauvre

      if(tfuel.gt.epsicp .and. tmpgaz.gt.epsicp) then

        if(f4s3(iel).le.f4m1(iel)) then
          p=((f4m1(iel)+f4m2(iel))/2.-f4s3(iel))/(1.-f4s3(iel))
          ts3s=((tmpgaz-dirac)/hrec(iel)                          &
                              /(f4m2(iel)-f4m1(iel))-taire*p)     &
                              /(1.d0-p)
        else if(f4s3(iel).ge.f4m2(iel)) then
          ts3s=((tmpgaz-dirac)/hrec(iel)                          &
               /(f4m2(iel)-f4m1(iel))-tfuel)                      &
               *(f4s3(iel)-f4cl(iel))                             &
               /((f4m1(iel)+f4m2(iel))/2.d0-f4cl(iel))+tfuel
        else
          ts3s=(2.*(tmpgaz-dirac)/hrec(iel)                       &
               - tfuel*(f4s3(iel)-f4m1(iel))**2                   &
                /(f4s3(iel)-f4cl(iel))                            &
               - taire*(f4m2(iel)-f4s3(iel))**2                   &
                /(1.d0-f4s3(iel)))                                &
      /((2.*f4m2(iel)-f4m1(iel)-f4s3(iel))                        &
               +(f4m1(iel)-f4cl(iel))*(f4s3(iel)-f4m1(iel))       &
                                     /(f4s3(iel)-f4cl(iel))       &
               -(f4m2(iel)-f4s3(iel))**2/(1.d0-f4s3(iel)))
        endif

        bb1 = max( f4m1(iel) , f4cl(iel) )
        bb2 = min( f4m2(iel) , f4s3(iel) )
        bb3 = max( f4m1(iel) , f4s3(iel) )
        bb4 = min( f4m2(iel) , 1.d0 )

! On definit quelques parametres intermediaires pour simplifier le code

        if (bb2 .gt. bb1 ) then
          lri = hrec(iel)
        else
          lri = 0.d0
        endif
        if (bb4 .gt. bb3 ) then
          lpa = hrec(iel)
        else
          lpa = 0.d0
        endif

        p = f4s3(iel) - f4cl(iel)
        q = bb2**2 - bb1**2
        r = bb2 - bb1
        s = 1.d0 - f4s3(iel)
        t = bb4**2 - bb3**2
        u = bb4 - bb3

! On calcul de TS3

!           TS3=(TMPGAZ-DIRAC+TFUEL*(((LRI*Q)/(2.D0*P))-((LRI*F4S3*R)/P))
!     &         - TAIRE*( ((LPA*T)/(2.D0*S)) +(LPA*U) -((LPA*U)/S)))
!     &         /
!     &         (((LRI*Q)/(2.D0*P))+(LRI*R)-((LRI*F4S3*R)/P)-((LPA*T)
!     &         /(2.D0*S)) +((LPA*U)/S) )

        gh1 = -lri*q/(p*2.d0)
        gh2 =  lri*r*f4s3(iel)/p
        gh3 =  lpa*t/(2.d0*s)
        gh4 =  lpa*u
        gh5 =  lpa*u/s
        gh6 =  lri*q/(p*2.d0)
        gh7 =  lri*r
        gh8 =  f4s3(iel)*lri*r/p
        gh9 = lpa*t/(s*2.d0)
        gh10= lpa*u/s
        gh11 = tmpgaz - dirac

        ts3num = gh11 - tfuel*(gh1+gh2) - taire*(gh3+gh4-gh5)
        ts3den = gh6 + gh7 -gh8 - gh9 + gh10


        ts3 = ts3num / ts3den
        if(abs(ts3-ts3s).gt.1.d0) then
          WRITE(NFECRA,*) 'TS3 paul-TS3 sandro ', IEL,TS3,TS3S
        endif
      endif

      spdf = hrec(iel) * (f4m2(iel) - f4m1(iel))

! Integration

      xo2 = oxyd(iel)*propce(iel,ipproc(immel))                   &
           /wmo2

      do i = 1, npart+1
        gs(i) = f4m1(iel) + dble(i-1)/dble(npart)                 &
                           *(f4m2(iel)-f4m1(iel))

        if( gs(i) .le. f4s3(iel) ) then
          tt(i) = (ts3-tfuel)/(f4s3(iel)-f4cl(iel))* gs(i)        &
               + ts3 - f4s3(iel)*(ts3 - tfuel)                    &
                                / ( f4s3(iel) - f4cl(iel) )

        else
          tt(i) = (taire -ts3)/(1.d0-f4s3(iel))*gs(i)             &
               + taire - (taire - ts3)/(1.d0-f4s3(iel))

        endif
      enddo

      if(xo2.gt.0.018d0) then
        bmoy=0.d0
      else if(xo2 .lt. 0.0025d0) then
        bmoy=1.d0
      else
        bmoy=(0.018d0-xo2)/(0.018d0-0.0025d0)
      endif

!      DGS est le pas d'integration

      dgs = ( f4m2(iel) - f4m1(iel) ) / dble(npart)

! Calcul de K1*EXP(-E1/T)

      if ( ipdf1 .eq. 1 ) then

        rsf4 = kk1*exp(-ee1/taire)*d4f4(iel)
        rscl = kk1*exp(-ee1/tfuel)*d4cl(iel)

        do i = 1, npart+1
          val(i) = kk1 * exp(-ee1/tt(i)) * hrec(iel)
        enddo

        propce(iel,iexp1)= rsf4 + rscl

        do i = 1, npart
          propce(iel,iexp1) = propce(iel,iexp1)                   &
                           +0.5d0*dgs*(val(i)+val(i+1))
        enddo

      endif

!  Calcul de K2*EXP(-E2/T)

      if ( ipdf2 .eq. 1 ) then

        if ( xo2 .gt. 0.d0 ) then
          propce(iel,iexp2) = kk2*exp(-ee2/taire)                 &
                                 *d4f4(iel)*xo2**bmoy             &
                             +kk2*exp(-ee2/tfuel)                 &
                                 *d4cl(iel)*xo2**bmoy

          do i = 1, npart+1
            val(i) = kk2*exp(-ee2/tt(i))*hrec(iel)
          enddo

          do i = 1, npart
            propce(iel,iexp2) = propce(iel,iexp2)                 &
                               +0.5d0*dgs*(val(i)+val(i+1))       &
                                     *xo2**bmoy
          enddo
        else
          propce(iel,iexp2) = 0.d0
        endif

      endif

!  Calcul de K3*EXP(-E3/T)

      if ( ipdf3 .eq. 1 ) then

        if ( xo2 .gt. 0.d0 ) then
          propce(iel,iexp3) = kk3*exp(-ee3/taire)                 &
                                 *d4f4(iel)*xo2**(0.5d0)          &
                             +kk3*exp(-ee3/tfuel)                 &
                                 *d4cl(iel)*xo2**(0.5d0)

          do i = 1, npart+1
            val(i) = kk3*exp(-ee3/tt(i))*hrec(iel)
          enddo

          do i = 1, npart
            propce(iel,iexp3)= propce(iel,iexp3)                  &
                     + 0.5d0*dgs*(val(i)+val(i+1))*xo2**(0.5d0)
          enddo
        else
          propce(iel,iexp3)= 0.d0
        endif
      endif

    endif

  enddo

endif

return
end subroutine
