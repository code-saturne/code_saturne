!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine cs_fuel_physprop1 &
!===========================
 ( ncelet , ncel   ,                                      &
   f1m    , f2m    , f3m    , f4m    , f5m    ,           &
   f6m    , f7m    , f8m    , f9m    , fvp2m  ,           &
   enth   , enthox ,                                      &
   rtp    , propce , rom1   )

!===============================================================================
! FONCTION :
! --------

! CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE GAZEUSE
!  VALEURS CELLULES
!  ----------------
!  TEMPERATURE, MASSE VOLUMIQUE ET CONCENTRATIONS MOYENNES
!  (UTILISATION D'UNE PDF RECTANGLE-DIRAC)
! ==> CHIMIE RAPIDE MODELE EN 3 POINTS
!     EXTENSION A TROIS COMBUSTIBLES POUR LE CHARBON PULVERISE
!                                         --------------------
! REACTIONS HETEROGENES
!   - Pyrolyse
!     Composition elementaire de la mole de matieres volatiles
!     Le charbon reactif s'ecrit C(1)H(ALPHA)O(BETA)
!       -(k1)-> ALPHA/4 CH4  + BETA CO + (1-ALPHA/4-BETA)    Coke
!     Charbon reactif
!       -(k2)-> ALPHA/Y CXHY + BETA CO + (1-ALPHA/RYSX-BETA) Coke
!       Avec RYSX = Y/X
!   - Combustion heterogene
!     Coke + 1/2 (O2 + XSI N2) -> CO + XSI/2 N2
!   - Reactions en phase gaz
! (4/(4-RYSX)) CH4 + (O2 + XSI N2)   -(1)->  4/X/(4-RYSX)*CXHY + 2 H2O
!                                           + XSI N2
! CXHY + X/4*(2+RYSX) (O2 + XSI N2)  -(2)->  X CO + Y/2 H2O
!                                           + X/4*(2+RYSX)*XSI N2
!           CO + 1/2 (O2 + XSI N2)  -(3)->  CO2 + XSI/2 N2
! CHOIX DES VARIABLES
!  F1 est la fractions massique des matieres volatiles : CH4  + CO
!  F2 est la fractions massique des matieres volatiles : CXHY + CO
!  F3 est la fraction massique de carbone venant de la combustion
!    heterogene

!  Soit Y les fractions massiques et Z les concentrations (moles/kg)
!    indice f avant reaction, b final

! PDF CONJOINTE DEGENERE EN UNE PDF 1D DE TYPE RECTANGLE - DIRAC

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nitbcp           ! e  ! <-- ! taille du macro tableau cp entiers             !
! nrtbcp           ! e  ! <-- ! taille du macro tableau cp reels               !
! nitbmc           ! e  ! <-- ! taille du macro tableau mc entiers             !
! nrtbmc           ! e  ! <-- ! taille du macro tableau mc reels               !
! nitbwo           ! e  ! <-- ! taille du macro tableau work entiers           !
! nrtbwo           ! e  ! <-- ! taille du macro tableau work reels             !
! pa               ! tr ! <-- ! pression absolue en pascals                    !
! f1m              ! tr ! <-- ! moyenne du traceur 1 mvl [chx1+co]             !
! f2m              ! tr ! <-- ! moyenne du traceur 2 mvl [chx2+co]             !
! f3m              ! tr ! <-- ! moyenne du traceur 3 (co c.het)                !
! f4m              ! tr ! <-- ! moyenne du traceur 4 (air)                     !
! f4m              ! tr ! <-- ! moyenne du traceur 5 (h2o)                     !
! f4p2m            ! tr ! <-- ! variance du traceur 4 (air)                    !
! enth             ! tr ! <-- ! enthalpie en j/kg  soit du gaz                 !
!                  !    !     !                    soit du melange             !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant)                  !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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
use dimens, only: nvar
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_fuel_incl

!===============================================================================

implicit none

! Arguments
integer          ncelet , ncel
!
double precision f1m(ncelet), f2m(ncelet) , f3m(ncelet)
double precision f4m(ncelet), f5m(ncelet) , f6m(ncelet)
double precision f7m(ncelet), f8m(ncelet) , f9m(ncelet)
double precision fvp2m(ncelet)
double precision enth(ncelet),enthox(ncelet)
double precision rtp(ncelet,nflown:nvar), propce(ncelet,*)
double precision rom1(ncelet)
! Local variables
integer          iel , ii ,ice , icla
integer          ipcyf1,ipcyf2,ipcyf3,ipcyf4,ipcyf5,ipcyf6,ipcyf7
integer          ipcyox,ipcyp1,ipcyp2,ipcyp3,ipcyin
integer          ipcte1

double precision wmolme
double precision somch , cfolc , cfolh , cfolo

integer          iok1 , iok2 , iok3 , iok4 , iok5
integer          , dimension ( : )     , allocatable :: intpdf
double precision , dimension ( : )     , allocatable :: fmini,fmaxi,ffuel
double precision , dimension ( : )     , allocatable :: dfuel,doxyd,pdfm1,pdfm2,hrec
double precision , dimension ( : )     , allocatable :: x2,cx1m,cx2m,wmf1,wmf2
double precision , dimension ( : , : ) , allocatable :: af1    , af2
double precision , dimension ( : )     , allocatable :: fs3no  , fs4no
double precision , dimension ( : , : ) , allocatable :: yfs4no
double precision, allocatable, dimension(:) :: tpdf

integer          ipass
data ipass / 0 /

!===============================================================================
! 0. Memory allocation
!===============================================================================

allocate(intpdf(1:ncel)                                       ,stat=iok1)
allocate(fmini(1:ncel)      ,fmaxi(1:ncel)      ,ffuel(1:ncel),stat=iok2)
allocate(dfuel(1:ncel)      ,doxyd(1:ncel)      ,pdfm1(1:ncel),stat=iok3)
allocate(pdfm2(1:ncel)      ,hrec(1:ncel)                     ,stat=iok3)
allocate(x2(1:ncel)         ,cx1m(1:ncel)       ,cx2m(1:ncel) ,stat=iok4)
allocate(wmf1(1:ncel)       ,wmf2(1:ncel)                     ,stat=iok4)
allocate(af1(1:ncel,1:ngazg),af2(1:ncel,1:ngazg)              ,stat=iok5)
!----
if ( iok1 > 0 .or. iok2 > 0 .or. iok3 > 0 .or. iok4 > 0 .or. iok5 > 0 ) then
  write(nfecra,*) ' Memory allocation error inside: '
  write(nfecra,*) '     cs_fuel_physprop1           '
  call csexit(1)
endif
!
if ( ieqnox .eq. 1 ) then
  !----
  allocate(fs3no(1:ncel) , fs4no(1:ncel),stat=iok1)
  !----
  if ( iok1 > 0 ) then
    write(nfecra,*) ' Memory allocation error inside: '
    write(nfecra,*) '     cs_fuel_physprop1           '
    write(nfecra,*) ' for fs3no and fs4no arrays      '
    call csexit(1)
  endif
  !----
  allocate(yfs4no(1:ncel,1:ngazg),stat=iok2)
  !----
  if ( iok2 > 0 ) then
    write(nfecra,*) ' Memory allocation error inside: '
    write(nfecra,*) '     cs_fuel_physprop1           '
    write(nfecra,*) ' for fs4 concentrations array    '
    call csexit(1)
  endif
endif
!===============================================================================
!
!===============================================================================
! 1. INITIALISATION
!===============================================================================
!
! pointeur
ipcyf1 = ipproc(iym1(ifo0))
ipcyf2 = ipproc(iym1(ifov))
ipcyf3 = ipproc(iym1(ico  ))
ipcyf4 = ipproc(iym1(ih2s ))
ipcyf5 = ipproc(iym1(ihy  ))
ipcyf6 = ipproc(iym1(ihcn ))
ipcyf7 = ipproc(iym1(inh3 ))
ipcyox = ipproc(iym1(io2  ))
ipcyp1 = ipproc(iym1(ico2 ))
ipcyp2 = ipproc(iym1(ih2o ))
ipcyp3 = ipproc(iym1(iso2 ))
ipcyin = ipproc(iym1(in2  ))
!
ipass = ipass + 1

!===============================================================================
! 2. DETERMINATION DU TYPE DE PDF
!===============================================================================
!
do iel = 1, ncel
!  bornes min et max de la pdf
  fmini(iel) = 0
  fmaxi(iel) = 1.d0
!  Somme de F1+F2 (pour le fuel, f1=0)
  ffuel(iel)=f1m(iel)+f2m(iel)
enddo

allocate(tpdf(ncelet))

call pppdfr &
!==========
 ( ncelet ,ncel   , intpdf ,                                      &
   tpdf   ,                                                       &
   ffuel  , fvp2m ,                                               &
   fmini  , fmaxi ,                                               &
   doxyd  , dfuel , pdfm1 , pdfm2 , hrec )

! Free memory
deallocate(tpdf)

!===============================================================================
! 2.CALCUL DES CONCENTRATIONS MOYENNES
!===============================================================================
! - Les masses molaires des combustibles gazeuse sont calculées dans
!   cs_fuel_readata.f90.
!
do iel = 1, ncel
  wmf1(iel) = wmole(ifo0)
  wmf2(iel) = wmole(ifov)
enddo
!
x2( : ) = 0.d0
do icla = 1, nclafu
  do iel = 1, ncel
    x2(iel) = x2(iel) + rtp(iel,isca(iyfol(icla)))
  enddo
enddo
!
do iel=1,ncel
!
  do ii=1,ngazg
    af1(iel,ii) = zero
    af2(iel,ii) = zero
  enddo
!
  af2(iel,ifov) = chfov/wmole(ifov)
  af2(iel,ico ) = cofov/wmole(ico )
  af2(iel,ih2s) = hsfov/wmole(ih2s)
!
  cx1m(iel) = 0.d0
  cx2m(iel) = nhcfov
!
enddo
!
call cs_gascomb &
!==============
 ( ncelet , ncel   , ifo0 , ifov ,                                &
   intpdf ,                                                       &
   rtp    , x2  ,                                                 &
   f1m    , f2m , f3m , f4m , f5m , f6m , f7m , f8m , f9m ,       &
   pdfm1  , pdfm2  , doxyd    , dfuel  , hrec ,                   &
   af1    , af2    , cx1m     , cx2m   , wmf1   , wmf2 ,          &
   propce(1,ipcyf1) , propce(1,ipcyf2) , propce(1,ipcyf3) ,       &
   propce(1,ipcyf4) , propce(1,ipcyf5) , propce(1,ipcyf6) ,       &
   propce(1,ipcyf7) , propce(1,ipcyox) , propce(1,ipcyp1) ,       &
   propce(1,ipcyp2) , propce(1,ipcyp3) , propce(1,ipcyin) ,       &
   fs3no , fs4no , yfs4no )
!
! --> Clipping eventuel des fractions massiques
!
do iel = 1, ncel
  do ice = 1, ngazg
    if ( propce(iel,ipproc(iym1(ice))) .lt. zero )  then
       propce(iel,ipproc(iym1(ice))) = zero
    endif
!
    if ( propce(iel,ipproc(iym1(ice))) .gt. 1.d0 )  then
       propce(iel,ipproc(iym1(ice))) = 1.d0
    endif
!
    if ( abs(propce(iel,ipproc(iym1(ice)))) .lt. epsicp )  then
       propce(iel,ipproc(iym1(ice))) = zero
    endif
  enddo
enddo
!
!===============================================================================
! 4. CALCUL DE LA TEMPERATURE ET DE LA MASSE VOLUMIQUE
!===============================================================================
!
ipcte1 = ipproc(itemp1)
!
call cs_fuel_thfieldconv1 &
!========================
 ( ncelet , ncel   ,                                              &
   enth   ,                                                       &
   propce(1,ipcyf1), propce(1,ipcyf2), propce(1,ipcyf3),          &
   propce(1,ipcyf4), propce(1,ipcyf5), propce(1,ipcyf6),          &
   propce(1,ipcyf7), propce(1,ipcyox), propce(1,ipcyp1),          &
   propce(1,ipcyp2), propce(1,ipcyp3), propce(1,ipcyin),          &
   propce(1,ipcte1) )
!
do iel = 1, ncel
  wmolme = propce(iel,ipcyf1)/wmole(ifo0)                         &
         + propce(iel,ipcyf2)/wmole(ifov)                         &
         + propce(iel,ipcyf3)/wmole(ico )                         &
         + propce(iel,ipcyf4)/wmole(ih2s)                         &
         + propce(iel,ipcyf5)/wmole(ihy )                         &
         + propce(iel,ipcyf6)/wmole(ihcn)                         &
         + propce(iel,ipcyf7)/wmole(inh3)                         &
         + propce(iel,ipcyox)/wmole(io2 )                         &
         + propce(iel,ipcyp1)/wmole(ico2)                         &
         + propce(iel,ipcyp2)/wmole(ih2o)                         &
         + propce(iel,ipcyp3)/wmole(iso2)                         &
         + propce(iel,ipcyin)/wmole(in2 )

! stockage de la masse molaire du melange

  propce(iel,ipproc(immel)) = 1.d0 / wmolme

! ---- On ne met pas la pression mecanique RTP(IEL,IPR)
!      mais P0

  rom1(iel) = p0 / (wmolme*rr*propce(iel,ipcte1))
enddo
!
!===============================================================================
! 5. MODELE DE NOX
!===============================================================================
!
!
! MODEL NOx : on y passe pas a la 1ere iter relatif

if ( ieqnox .eq. 1 .and. ipass .gt. 1 ) then
!
  call cs_fuel_noxst &
 !==================
 ( ncelet , ncel   ,                                              &
   intpdf ,                                                       &
   pdfm1  , pdfm2  , doxyd  , dfuel  , hrec ,                     &
   f3m    , f4m    , f5m    , f6m    , f7m  , f8m , f9m ,         &
   fs3no  , fs4no  , yfs4no , enthox ,                            &
   rtp    , propce  )
!
else if ( ieqnox .eq. 1 ) then
!
  do iel = 1, ncel
    propce(iel,ipproc(ighcn1)) = 0.d0
    propce(iel,ipproc(ighcn2)) = 0.d0
    propce(iel,ipproc(ignoth)) = 0.d0
  enddo
!
endif

!===============================================================================
! 6. CALCUL DES BILANS en C , O et H
!===============================================================================
!
do iel=1,ncel
  propce(iel,ipproc(ibcarbone )) = (1.d0-x2(iel))                      &
            *( propce(iel,ipcyf1)*wmolat(iatc)/wmole(ifo0)             &
              +propce(iel,ipcyf2)*wmolat(iatc)/wmole(ifov)             &
              +propce(iel,ipcyf3)*wmolat(iatc)/wmole(ico )             &
              +propce(iel,ipcyf6)*wmolat(iatc)/wmole(ihcn)             &
              +propce(iel,ipcyp1)*wmolat(iatc)/wmole(ico2) )

  propce(iel,ipproc(iboxygen  )) = (1.d0-x2(iel))                      &
            *( propce(iel,ipcyf3)*     wmolat(iato)/wmole(ico )        &
              +propce(iel,ipcyox)*2.d0*wmolat(iato)/wmole(io2 )        &
              +propce(iel,ipcyp1)*2.d0*wmolat(iato)/wmole(ico2)        &
              +propce(iel,ipcyp2)*     wmolat(iato)/wmole(ih2o)        &
              +propce(iel,ipcyp3)*2.d0*wmolat(iato)/wmole(iso2) )

  propce(iel,ipproc(ibhydrogen)) = (1.d0-x2(iel))                      &
            *( propce(iel,ipcyf1)*nhcfov*wmolat(iath)/wmole(ifo0)      &
              +propce(iel,ipcyf2)*nhcfov*wmolat(iath)/wmole(ifov)      &
              +propce(iel,ipcyf4)*2.d0     *wmolat(iath)/wmole(ih2s)   &
              +propce(iel,ipcyf5)*2.d0     *wmolat(iath)/wmole(ihy )   &
              +propce(iel,ipcyf6)*          wmolat(iath)/wmole(ihcn)   &
              +propce(iel,ipcyf7)*3.d0     *wmolat(iath)/wmole(inh3)   &
              +propce(iel,ipcyp2)*2.d0     *wmolat(iath)/wmole(ih2o) )
!
enddo

do icla = 1, nclafu
  do iel = 1, ncel
!
    somch = cfol + hfol + ofol
    cfolc  = cfol/somch
    cfolh  = hfol/somch
    cfolo  = ofol/somch
!
    propce(iel,ipproc(ibcarbone )) = propce(iel,ipproc(ibcarbone ))    &
                                    +rtp(iel,isca(iyfol(icla)))*cfolc
!
    propce(iel,ipproc(iboxygen )) = propce(iel,ipproc(iboxygen ))    &
                                   +rtp(iel,isca(iyfol(icla)))*cfolo
!
    propce(iel,ipproc(ibhydrogen)) = propce(iel,ipproc(ibhydrogen))  &
                                    +rtp(iel,isca(iyfol(icla)))*cfolh
!
  enddo
enddo
!----
! Formats
!----

!===============================================================================
! Deallocation dynamic arrays
!----
deallocate(intpdf,                      stat=iok1)
deallocate(fmini,fmaxi,ffuel,           stat=iok2)
deallocate(dfuel,doxyd,pdfm1,pdfm2,hrec,stat=iok3)
deallocate(x2,cx1m,cx2m,wmf1,wmf2,      stat=iok4)
deallocate(af1, af2,                    stat=iok5)
!----
if ( iok1 > 0 .or. iok2 > 0 .or. iok3 > 0 .or. iok4 > 0 .or. iok5 > 0 ) then
  write(nfecra,*) ' Memory deallocation error inside: '
  write(nfecra,*) '     cs_fuel_physprop1             '
  call csexit(1)
endif
if ( ieqnox .eq. 1 ) then
  !----
  deallocate(fs3no,fs4no,stat=iok1)
  !----
  if ( iok1 > 0 ) then
    write(nfecra,*) ' Memory deallocation error inside: '
    write(nfecra,*) '     cs_fuel_physprop1             '
    write(nfecra,*) ' for fs3no and fs4no arrays        '
    call csexit(1)
  endif
  !----
  deallocate(yfs4no,stat=iok2)
  !----
  if ( iok2 > 0 ) then
    write(nfecra,*) ' Memory deallocation error inside: '
    write(nfecra,*) '     cs_fuel_physprop1             '
    write(nfecra,*) ' for fs4 concentrations array      '
    call csexit(1)
  endif
endif
!===============================================================================

!----
! End
!----

return
end subroutine
