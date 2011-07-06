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

subroutine cpphy1 &
!================

 ( ncelet , ncel   ,                                              &
   nitbcp , nrtbcp , nitbmc , nrtbmc , nitbwo , nrtbwo ,          &
   f1m    , f2m    , f3m    , f4m    , f5m    ,                   &
   f6m    , f7m    , f4p2m  , f3max  ,                            &
   enth   ,                                                       &
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

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel
integer          nitbcp , nrtbcp
integer          nitbmc , nrtbmc
integer          nitbwo , nrtbwo

double precision f1m(ncelet), f2m(ncelet)
double precision f3m(ncelet), f4m(ncelet) , f5m(ncelet)
double precision f6m(ncelet), f7m(ncelet)
double precision f4p2m(ncelet)
double precision f3max(ncelet)
double precision enth(ncelet)
double precision rtp(ncelet,*), propce(ncelet,*)
double precision rom1(ncelet)

! Local variables

integer          iel    , ice
integer          iitbcp , iitbmc , iitbwo
integer          ipcte1
integer          ipcyf1 , ipcyf2 , ipcyf3 , ipcyox
integer          ipcyp1 , ipcyp2 , ipcyin , ipcyce
integer          ipdf   , nbf3

double precision wmolme , wmchx1 , wmchx2
double precision f1cl , f2cl , f3cl , f4cl , f5cl
double precision f6cl , f7cl
double precision fp1  , fp2  , fp4  , fp5, fp6, fp7
double precision xalpha , xbeta , xgam
double precision f4cmin,f4cmax,xden

integer, allocatable, dimension(:,:) :: itbcp, itbmc, itbwo

double precision, allocatable, dimension(:,:) :: rtbcp, rtbmc, rtbwo

!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

! Allocate temporary arrays
allocate(itbcp(ncelet,nitbcp))
allocate(itbmc(ncelet,nitbmc))
allocate(itbwo(ncelet,nitbwo))

allocate(rtbcp(ncelet,nrtbcp))
allocate(rtbmc(ncelet,nrtbmc))
allocate(rtbwo(ncelet,nrtbwo))

! --- Initialisation memoire


! --- Initialisation des tableaux d'entiers de travail

do iel = 1, ncel
  do iitbcp = 1, nitbcp
    itbcp(iel,iitbcp) = 0
  enddo
  do iitbmc = 1, nitbmc
    itbmc(iel,iitbmc) = 0
  enddo
  do iitbwo = 1, nitbwo
    itbwo(iel,iitbwo) = 0
  enddo
enddo

! --- Initialisation des tableaux de reels de travail

do iel = 1, ncel
  do iitbcp = 1, nrtbcp
    rtbcp(iel,iitbcp) = zero
  enddo
  do iitbmc = 1, nrtbmc
    rtbmc(iel,iitbmc) = zero
  enddo
  do iitbwo = 1, nrtbwo
    rtbwo(iel,iitbwo) = zero
  enddo
enddo

!===============================================================================
! 2. DETERMINATION DU TYPE DE PDF
!===============================================================================

nbf3 = 0
do iel = 1, ncel

!       Calcul de Fcl

  if ( f1m(iel)+f2m(iel) .gt. 0.d0 ) then
    fp1 = f1m(iel)/(f1m(iel)+f2m(iel))
    fp2 = f2m(iel)/(f1m(iel)+f2m(iel))
  else
    fp1 = 0.d0
    fp2 = 0.d0
  endif
  if (f4m(iel)+f5m(iel)+f6m(iel)+f7m(iel) .gt. 0.d0 ) then
    fp4 = f4m(iel)/(f4m(iel)+f5m(iel)+f6m(iel)+f7m(iel))
    fp5 = f5m(iel)/(f4m(iel)+f5m(iel)+f6m(iel)+f7m(iel))
    fp6 = f6m(iel)/(f4m(iel)+f5m(iel)+f6m(iel)+f7m(iel))
    fp7 = f7m(iel)/(f4m(iel)+f5m(iel)+f6m(iel)+f7m(iel))
  else
    fp4 = 0.d0
    fp5 = 0.d0
    fp6 = 0.d0
    fp7 = 0.d0
  endif

!      Calcul de F3MAX

  if ( ( ao2f4*f4m(iel)+ao2f6*f6m(iel)                            &
        +ao2f7*f7m(iel) ) .gt. 0.d0     ) then
    f3max(iel) = 1.d0                                             &
     /(1.d0+0.5d0*(1.d0/wmolat(iatc))                             &
                 *(f4m(iel)+f5m(iel)+f6m(iel)+f7m(iel))           &
                 /(ao2f4*f4m(iel)+ao2f6*f6m(iel)+ao2f7*f7m(iel)) )
  else
    f3max(iel) = 1.d0
  endif

  if ( f3m(iel) .gt. f3max(iel) ) then
    f3m(iel) = f3max(iel)
    nbf3 = nbf3 + 1
  endif

  xden = f3m(iel)/f3max(iel)+f1m(iel)+f2m(iel)

  if ( xden .gt. 1.d-15 ) then

    xalpha = 1.d0-1.d0/xden

    f1cl = (1.d0-xalpha)*f1m(iel)
    f2cl = (1.d0-xalpha)*f2m(iel)
    f3cl = (1.d0-xalpha)*f3m(iel)
    f4cl = xalpha*fp4 + (1.d0-xalpha)*f4m(iel)
    f5cl = xalpha*fp5 + (1.d0-xalpha)*f5m(iel)
    f6cl = xalpha*fp6 + (1.d0-xalpha)*f6m(iel)
    f7cl = xalpha*fp7 + (1.d0-xalpha)*f7m(iel)

  else
    f1cl = fp1
    f2cl = fp2
    f3cl = 0.d0
    f4cl = 0.d0
    f5cl = 0.d0
    f6cl = 0.d0
    f7cl = 0.d0
  endif

  rtbcp(iel,1) = f1cl
  rtbcp(iel,2) = f2cl
  rtbcp(iel,3) = f3cl
  rtbcp(iel,4) = f4cl
  rtbcp(iel,5) = f5cl
  rtbcp(iel,6) = f6cl
  rtbcp(iel,7) = f7cl

!       bornes min et max de la pdf : F4CL a 1

  rtbcp(iel,8) = 1.d0

! Somme de F4+F5

  rtbcp(iel,14) = f4m(iel)+f5m(iel)+f6m(iel)+f7m(iel)

enddo

if ( irangp .ge. 0 ) then
  call parcpt(nbf3)
endif
WRITE(NFECRA,*) ' Nombre de clipping sur F3 : ',NBF3

call pppdfr                                                       &
!==========
 ( ncelet,ncel,                                                   &
   itbcp(1,1) ,                                                   &
   rtbcp(1,14), f4p2m ,                                           &
!          F4+F5
   rtbcp(1,4), rtbcp(1,8),                                        &
!           FMINI        FMAXI
   rtbcp(1,9) , rtbcp(1,10) , rtbcp(1,11) , rtbcp(1,12),          &
!           D4CL         D4F4          F4M1        F4M2
    rtbcp(1,13) )
!           HREC

!===============================================================================
! 2.CALCUL DES CONCENTRATIONS MOYENNES
!===============================================================================


ipcyf1 = ipproc(iym1(ichx1))
ipcyf2 = ipproc(iym1(ichx2))
ipcyf3 = ipproc(iym1(ico  ))
ipcyox = ipproc(iym1(io2  ))
ipcyp1 = ipproc(iym1(ico2 ))
ipcyp2 = ipproc(iym1(ih2o ))
ipcyin = ipproc(iym1(in2  ))

call cpcym2                                                       &
!==========
 ( ncelet , ncel   , nrtbmc ,                                     &
   itbcp(1,1) ,                                                   &
!         INTPDF
   rtp    ,                                                       &
   f1m    , f2m , f3m , f4m , f5m , f6m , f7m ,                   &
   f3max  ,                                                       &
   rtbcp(1,1),rtbcp(1,2),rtbcp(1,3),rtbcp(1,4),rtbcp(1,5),        &
!         F1CL        F2CL       F3CL      F4CL        F5CL
  rtbcp(1,6) ,  rtbcp(1,7) ,                                      &
!         F6CL         F7CL
   rtbcp(1,11) , rtbcp(1,12) , rtbcp(1,9) ,                       &
!           F4M1         F4M2        D4CL
   rtbcp(1,10) ,rtbcp(1,13) ,                                     &
!           D4F4         HREC
   rtbmc      , rtbwo(1,1)  ,                                     &
   propce(1,ipcyf1) , propce(1,ipcyf2) , propce(1,ipcyf3) ,       &
   propce(1,ipcyox) , propce(1,ipcyp1) , propce(1,ipcyp2) ,       &
   propce(1,ipcyin)  )

! --> Clipping eventuel des fractions massiques

do iel = 1, ncel
  do ice = 1, (ngaze-2*ncharb)
    ipcyce = ipproc(iym1(ice))
    if ( abs(propce(iel,ipcyce)) .lt. epsicp )                    &
         propce(iel,ipcyce) = zero
  enddo
enddo

!===============================================================================
! 4. CALCUL DE LA TEMPERATURE ET DE LA MASSE VOLUMIQUE
!===============================================================================

ipcte1 = ipproc(itemp1)

! --- Transport d'H2

call cpteh1                                                       &
!==========
 ( ncelet , ncel   , nitbmc , nrtbmc ,                            &
   enth,                                                          &
   propce(1,ipcyf1), propce(1,ipcyf2), propce(1,ipcyf3),          &
   propce(1,ipcyox), propce(1,ipcyp1), propce(1,ipcyp2),          &
   propce(1,ipcyin),                                              &
   propce(1,ipcte1),                                              &
   itbmc      , rtbmc      ,                                      &
!          MACRO TABLEAU MULTI CHARBONS ENTIERS REELS
   rtbwo(1,1) , rtbwo(1,2) )
!          TABLEAUX DE TRAVAIL

ipcte1 = ipproc(itemp1)
do iel = 1, ncel
  wmchx1 = wmolat(iatc)+rtbmc(iel,ix1mc)*wmolat(iath)
  wmchx2 = wmolat(iatc)+rtbmc(iel,ix2mc)*wmolat(iath)
  wmolme = propce(iel,ipcyf1)/wmchx1                              &
         + propce(iel,ipcyf2)/wmchx2                              &
         + propce(iel,ipcyf3)/wmole(ico )                         &
         + propce(iel,ipcyox)/wmole(io2 )                         &
         + propce(iel,ipcyp1)/wmole(ico2)                         &
         + propce(iel,ipcyp2)/wmole(ih2o)                         &
         + propce(iel,ipcyin)/wmole(in2 )

! stockage de la masse molaire du melange

  propce(iel,ipproc(immel)) = 1.d0 / wmolme

! ---- On ne met pas la pression mecanique RTP(IEL,IPR)
!      mais P0

  rom1(iel) = p0 / (wmolme*rr*propce(iel,ipcte1))
enddo

! Free memory
deallocate(itbcp, itbmc, itbwo)
deallocate(rtbcp, rtbmc, rtbwo)

!===============================================================================

!--------
! FORMATS
!--------


!----
! FIN
!----

return
end subroutine
