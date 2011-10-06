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

subroutine d3pint &
!================

 ( ncelet , ncel   , indpdf ,                                     &
   dirmin , dirmax , fdeb   , ffin   , hrec   ,                   &
   fm     , hm     , p      ,                                     &
   propce ,                                                       &
   w1      )

!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE PHYSIQUE PARTICULIERE : FLAMME DE DIFFUSION
! Integration des variables thermodynamiques en fonction de
!  la fraction de melange
!  Rq : Il serait judicieux de ponderer l'integration de la
!           temperature par les CP

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! indpdf           ! ti ! <-- ! indicteur passage ou non par les pdf           !
! dirmin           ! tr ! <-- ! pdf : dirac en fmin                            !
! dirmax           ! tr ! <-- ! pdf : dirac en fmax                            !
! fdeb             ! tr ! <-- ! pdf : abscisse debut rectangle                 !
! ffin             ! tr ! <-- ! pdf : abscisse fin rectangle                   !
! hrec             ! tr ! <-- ! pdf : hauteur rectangle                        !
! fm               ! tr ! <-- ! fraction de melange moyenne                    !
! hm               ! tr ! <-- ! enthalpie massique moyenne                     !
!                  !    !     !  si ecoulement permeatique                     !
! p                ! tr ! <-- ! pression                                       !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     ! cellules ( concentrations, temp. )   !                                  !
! w1               ! tr ! --- ! tableau de tavail                              !
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
use entsor
use cstnum
use ppppar
use ppthch
use coincl
use ppincl
use radiat

!===============================================================================

implicit none

! Arguments

integer          ncelet, ncel
integer          indpdf(ncelet)
double precision dirmin(ncelet), dirmax(ncelet)
double precision fdeb(ncelet), ffin(ncelet), hrec(ncelet)
double precision fm(ncelet), hm(ncelet), p(ncelet)
double precision propce(ncelet,*), w1(ncelet)


! Local variables

integer          icel, icg
integer          ih, if, jh, jf, ipcrom
integer          ipcsca, ipctem, ipckab, ipct4, ipct3
double precision aa1, bb1, aa2, bb2, f1, f2, a, b, fmini, fmaxi
double precision u, v, c, d, temsmm, fsir


!===============================================================================

integer       ipass
data          ipass /0/
save          ipass


!===============================================================================

!===============================================================================
! 0. ON COMPTE LES PASSAGES
!===============================================================================

! Initialize variables to avoid compiler warnings

ipckab = 0
ipct3 = 0
ipct4 = 0

aa1 = 0.d0
aa2 = 0.d0
bb1 = 0.d0
bb2 = 0.d0


ipass = ipass + 1


!===============================================================================
! 1. INTEGRATION DES NGAZG FRACTIONS MASSIQUES D'ESPECES GLOBALES
!===============================================================================

! ---> En flamme de diffusion chimie 3 points :
!      - il n'y a qu'une seule reaction globale (IR= )
!      - le taux de melange f varie entre 0 et 1
fsir = fs(1)
fmini = zero
fmaxi = 1.d0

do icel = 1, ncel

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

    ipcsca = ipproc(iym(icg))

    if ( indpdf(icel) .eq. 1 ) then

! ---> Integration de la PDF

      propce(icel,ipcsca) = dirmin(icel) * ( aa1 + bb1 * fmini )  &
                          + dirmax(icel) * ( aa2 + bb2 * fmaxi )
      if ( fdeb(icel).lt.fsir ) then
        f1 = fdeb(icel)
        f2 = min( fsir,ffin(icel) )
        propce(icel,ipcsca) = propce(icel,ipcsca)                 &
             + hrec(icel)*(f2-f1)*(aa1+bb1*5.d-1*(f2+f1))
      endif
      if ( ffin(icel).gt.fsir ) then
        f1 = max(fsir,fdeb(icel))
        f2 = ffin(icel)
        propce(icel,ipcsca) = propce(icel,ipcsca)                 &
             + hrec(icel)*(f2-f1)*(aa2+bb2*5.d-1*(f2+f1))
      endif
    else

! ---> Degenerescence sur la valeur moyenne

      if ( fm(icel).le.fsir ) then
        propce(icel,ipcsca) = aa1+bb1*fm(icel)
      else
        propce(icel,ipcsca) = aa2+bb2*fm(icel)
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

do icel = 1, ncel
  w1(icel) = hstoea
enddo

if ( ippmod(icod3p).eq.1 ) then

  call d3phst                                                     &
  !==========
  ( ncelet , ncel    , indpdf ,                                   &
    dirmin , dirmax  , fdeb   , ffin   , hrec   ,                 &
    fm     , hm      ,                                            &
    w1      )

endif


!===============================================================================
! 3. INTEGRATION a) DE LA TEMPERATURE
!                b) DU COEFFICIENT D'ABSORPTION si rayonnement
!                c) DES TERME T^4 et T^3 si rayonnement
!                d) DE LA MASSE VOLUMIQUE
!===============================================================================


! ---> Positions des variables, coefficients

ipctem = ipproc(itemp)
if ( iirayo.gt.0 ) then
  ipckab = ipproc(ickabs)
  ipct4  = ipproc(it4m)
  ipct3  = ipproc(it3m)
endif

ipcrom = ipproc(irom)

do icel = 1, ncel

  if ( indpdf(icel) .eq. 1 ) then

! ---> Integration de la PDF

    ih = 1
    do jh = 1,(nmaxh-1)
      if ( w1(icel).gt.hh(jh+1) .and. w1(icel).le.hh(jh) )        &
           ih = jh
    enddo
    if ( w1(icel) .ge. hh(1)     ) ih = 1
    if ( w1(icel) .le. hh(nmaxh) ) ih = nmaxh-1
    propce(icel,ipctem) = dirmin(icel)*tinoxy +                   &
                          dirmax(icel)*tinfue
    temsmm = dirmin(icel)/wmolg(2)*tinoxy                         &
           + dirmax(icel)/wmolg(1)*tinfue
    if ( iirayo.gt.0 ) then
      propce(icel,ipckab) =                                       &
        dirmin(icel)*ckabsg(2)  + dirmax(icel)*ckabsg(1)
      propce(icel,ipct4) =                                        &
        dirmin(icel)*tinoxy**4 + dirmax(icel)*tinfue**4
      propce(icel,ipct3) =                                        &
        dirmin(icel)*tinoxy**3 + dirmax(icel)*tinfue**3
    endif
    if = 1
    do jf = 1, (nmaxf-1)
      if ( fdeb(icel).ge.ff(jf) .and.                             &
           fdeb(icel).lt.ff(jf+1) ) if = jf
    enddo
    if ( fdeb(icel) .le. ff(1)     ) if = 1
    if ( fdeb(icel) .ge. ff(nmaxf) ) if = nmaxf-1
    f2 = zero
    f1 = fdeb(icel)
    do while ( (ffin(icel)-f2).gt.epzero )
      f2 = min(ff(if+1),ffin(icel))
!          Dans le tableau TFH,
!           on extrait sur chaque ligne i : T = Ai+Bi*F
!           et on construit pour la valeur courante de HSTOE (W1)
!                               T = A+B*F
      aa1 = tfh(if,ih)
      bb1 = (tfh(if+1,ih)-tfh(if,ih))/(ff(if+1)-ff(if))
      aa2 = tfh(if,ih+1)
      bb2 = (tfh(if+1,ih+1)-tfh(if,ih+1))/(ff(if+1)-ff(if))
      a = aa1 + (w1(icel)-hh(ih))*(aa2-aa1)/(hh(ih+1)-hh(ih))
      b = bb1 + (w1(icel)-hh(ih))*(bb2-bb1)/(hh(ih+1)-hh(ih))
      a = a - b*ff(if)

! ----- Calcul de la temperature par integration

      propce(icel,ipctem) = propce(icel,ipctem)                   &
         + hrec(icel)*(f2-f1)*(a+b*(f1+f2)/2.d0)

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

        propce(icel,ipckab) = propce(icel,ipckab) +               &
          hrec(icel)*( u*(f2-f1) + v*(f2**2-f1**2)*0.5d0 )

        propce(icel,ipct4) = propce(icel,ipct4) +                 &
          hrec(icel)*                                             &
      (      a**4            * (f2-f1)                            &
   +   (4.d0*a**3  *b      ) * (f2**2-f1**2)/2.d0                 &
   +   (6.d0*(a**2)*(b**2) ) * (f2**3-f1**3)/3.d0                 &
   +   (4.d0*a     *(b**3) ) * (f2**4-f1**4)/4.d0                 &
   +   (            (b**4) ) * (f2**5-f1**5)/5.d0  )

        propce(icel,ipct3) = propce(icel,ipct3) +                 &
          hrec(icel)*                                             &
      (      (a**3)          * (f2-f1)                            &
   +   (3.d0*(a**2)*b      ) * (f2**2-f1**2)/2.d0                 &
   +   (3.d0*a     *(b**2) ) * (f2**3-f1**3)/3.d0                 &
   +   (            (b**3) ) * (f2**4-f1**4)/4.d0  )

      endif

! ----- Calcul du terme Temperature/masse molaire

      temsmm = temsmm + hrec(icel)*                               &
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
      if ( w1(icel).gt.hh(jh+1) .and. w1(icel).le.hh(jh) )        &
            ih = jh
    enddo
    if ( w1(icel) .ge. hh(1)     ) ih =1
    if ( w1(icel) .le. hh(nmaxh) ) ih =nmaxh-1
    if = 1
    do jf = 1, (nmaxf-1)
      if ( fm(icel).ge.ff(jf) .and. fm(icel).lt.ff(jf+1) )        &
           if = jf
    enddo
    if ( fm(icel) .le. ff(1)     ) if = 1
    if ( fm(icel) .ge. ff(nmaxf) ) if = nmaxf-1
    aa1 = tfh(if,ih)
    bb1 = (tfh(if+1,ih)-tfh(if,ih))/(ff(if+1)-ff(if))
    aa2 = tfh(if,ih+1)
    bb2 = (tfh(if+1,ih+1)-tfh(if,ih+1))/(ff(if+1)-ff(if))
    a  = aa1 + (w1(icel)-hh(ih))*(aa2-aa1)/(hh(ih+1)-hh(ih))
    b  = bb1 + (w1(icel)-hh(ih))*(bb2-bb1)/(hh(ih+1)-hh(ih))
    a  = a - b*ff(if)

! ----- Calcul de la temperature a partir de la valeur moyenne

    propce(icel,ipctem) = a+b*fm(icel)

    if ( fm(icel).lt.fsir ) then
!         On a demarre cote pauvre
      c =   1.d0/wmolg(2)
      d = (-1.d0/wmolg(2)+1.d0/wmolg(3))/fsir
    else
!         On termine cote riche (en commencant avec f1=fs)
      c = (  -fsir/wmolg(1)+1.d0/wmolg(3))/(1.d0-fsir)
      d = (   1.d0/wmolg(1)-1.d0/wmolg(3))/(1.d0-fsir)
     endif

    if ( iirayo.gt.0 ) then
      if ( fm(icel).lt.fsir ) then
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

      propce(icel,ipckab) = u + v*fm(icel)
      propce(icel,ipct4) = a**4                                   &
       + (4.d0*(a**3)*b      ) * fm(icel)                         &
       + (6.d0*(a**2)*(b**2) ) * fm(icel)**2                      &
       + (4.d0*a     *(b**3) ) * fm(icel)**3                      &
       + (            (b**4) ) * fm(icel)**4

      propce(icel,ipct3) = a**3                                   &
       + ( 3.d0*(a**2)*b      ) * fm(icel)                        &
       + ( 3.d0*a     *(b**2) ) * fm(icel)**2                     &
       + (             (b**3) ) * fm(icel)**3

    endif

! ----- Calcul du terme Temperature/masse molaire

    temsmm = a*c +(c*b+a*d)*fm(icel) + b*d*fm(icel)**2

  endif

! ---> Calcul de la masse volumique

  if (ipass.gt.1.or.(isuite.eq.1.and.initro.eq.1)) then
    propce(icel,ipcrom) = srrom*propce(icel,ipcrom)               &
                        + (1.d0-srrom)*                           &
                        ( p0/(rr*temsmm) )
  endif

enddo


return
end subroutine

