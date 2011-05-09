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

subroutine fuphy1 &
!================

 ( idbia0 , idbra0 ,                                              &
   ncelet , ncel   ,                                              &
   nitbfu , nrtbfu , nitbwo , nrtbwo ,                            &
   fvap   , fhtf   , f4p2m  ,                                     &
   enth   ,                                                       &
   rtp    , propce , rom1   ,                                     &
   itbfu  , rtbfu  ,                                              &
   itbwo  , rtbwo   )

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
!   - Evaporation
!     Composition de la vapeur (FOV pour Fuel Oil Vapor)

!        Le FOV est supposé être un mélange de H2S, CO, CHn
!          Les fractions massiques sont HSFOV pour H2S
!                                       COFOV      CO
!                                       CHFOV      CHn
!          l'hydrocarbure moyen est détermine par nHCFOV

!   - Combustion heterogene
!     La composition massque élémentairee du coke est donnée par
!          CKF, HKF, OKF, SKF
!           (et InKF inertes qui resteront dans l'inclusion)
!      lors de la réaction héterognène, on dégaze H2S, H2O, CO

!           Attention, ceci signifie qu'en presence de FHET il y a
!           eut prélevement d'O2 dans l'air environnant
!           (avant les réaction homogènes).

!   - Reactions en phase gaz

!     Avec l'O2 restant dans l'air (après dilution et oxydation hétérogène)
!     on considère des réaction s successives dans leur ordre
!          de priorité pour l'accès à l'O2
!           CHn + (1/2+n/4)O2 -(1)-> CO  + n/2 H2O
!                H2S + 3/2 O2      -(2)-> SO2 + H2O
!           CO + 1/2 O2       -(3)-> CO2

! CHOIX DES VARIABLES


!  Soit Y les fractions massiques et Z les concentrations (moles/kg)
!    indice f avant reaction, b final


! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nitbfu           ! e  ! <-- ! taille du macro tableau fuel entiers           !
! nrtbfu           ! e  ! <-- ! taille du macro tableau fuel reels             !
! nitbwo           ! e  ! <-- ! taille du macro tableau work entiers           !
! nrtbwo           ! e  ! <-- ! taille du macro tableau work reels             !
! pa               ! tr ! <-- ! pression absolue en pascals                    !
! fvap             ! tr ! <-- ! moyenne du traceur 1 fov [chn+co]              !
! fhtf             ! tr ! <-- ! moyenne du traceur 3 (co c.het)                !
! f4p2m            ! tr ! <-- ! variance du traceur 4 (air)                    !
! enth             ! tr ! <-- ! enthalpie en j/kg  soit du gaz                 !
!                  !    !     !                    soit du melange             !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant)                  !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! itbfu            ! tr ! <-- ! macro tableau entier fuel travail              !
! rtbfu            ! tr ! <-- ! macro tableau reel   fuel travail              !
! itbwo            ! tr ! <-- ! macro tableau entier travail                   !
! rtbwo            ! tr ! <-- ! macro tableau reel   travail                   !
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
use ppppar
use ppthch
use coincl
use cpincl
use fuincl
use ppincl
use ppcpfu

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          ncelet , ncel
integer          nitbfu , nrtbfu
integer          nitbwo , nrtbwo
integer          itbfu(ncelet,nitbfu)
integer          itbwo(ncelet,nitbwo)

double precision fvap(ncelet), fhtf(ncelet)
double precision f4p2m(ncelet), enth(ncelet)
double precision rtp(ncelet,*), propce(ncelet,*)
double precision rom1(ncelet)
double precision rtbfu(ncelet,nrtbfu)
double precision rtbwo(ncelet,nrtbwo)

! Local variables

integer          idebia , idebra
integer          iel    , iphas  , ice
integer          ipcte1
integer          ipcyf1 , ipcyf3 , ipcyox
integer          ipcyp1 , ipcyp2 , ipcyin , ipcyce
integer          ipcy2s , ipcyso
double precision wmolme
double precision f1m,f3m,f4m,f1cl,f3cl,f4cl

!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

! --- Initialisation memoire

idebia = idbia0
idebra = idbra0

!===============================================================================
! 2. DETERMINATION DU TYPE DE PDF
!===============================================================================

do iel = 1, ncel

!       Traceur virtuel au point moyen

  f1m =  fvap(iel)
  f3m =  fhtf(iel)/ff3max
  f4m = 1.d0 - f1m - f3m

!       Calcul des caractéristiques du point correspondant au
!       combustible moyen
!       F3cl : fraction de masse provenant de F3max et non de F3

  if ( (f3m+f1m*ff3max) .gt. 0.d0 ) then
    f1cl = f1m*ff3max/(f3m+f1m*ff3max)
  else
    f1cl = 0.d0
  endif

  f3cl = 1.d0-f1cl
  f4cl = (1.d0-ff3max)*f3cl

  rtbfu(iel,1) = f3m
  rtbfu(iel,2) = f4m
  rtbfu(iel,3) = f1cl
  rtbfu(iel,4) = f3cl
  rtbfu(iel,5) = f4cl

!       bornes min et max de la pdf : F4CL a 1

  rtbfu(iel,6) = 1.d0

enddo

call pppdfr                                                       &
!==========
 ( ncelet,ncel,                                                   &
   itbfu(1,1) ,                                                   &
   rtbfu(1,2), rtp(1,isca(if4p2m)),                               &
!           F4M
   rtbfu(1,5), rtbfu(1,6),                                        &
!           FMINI        FMAXI
   rtbfu(1,7) , rtbfu(1,8) , rtbfu(1,9) , rtbfu(1,10),            &
!           D4CL         D4F4          F4M1        F4M2
    rtbfu(1,11) )
!           HREC

!===============================================================================
! 2.CALCUL DES CONCENTRATIONS MOYENNES
!===============================================================================


ipcyf1 = ipproc(iym1(ifov))
ipcyf3 = ipproc(iym1(ico  ))
ipcyox = ipproc(iym1(io2  ))
ipcyp1 = ipproc(iym1(ico2 ))
ipcyp2 = ipproc(iym1(ih2o ))
ipcyin = ipproc(iym1(in2  ))
ipcy2s = ipproc(iym1(ih2s ))
ipcyso = ipproc(iym1(iso2 ))

 call fucym1                                                      &
!!==========
 ( ncelet , ncel   ,                                              &
   itbfu(1,1) ,                                                   &
!         INTPDF
   rtp    ,                                                       &
   fvap   ,   rtbfu(1,1) , rtbfu(1,2) ,                           &
!         F1M           F3M         F4M
  rtbfu(1,3) , rtbfu(1,4) ,rtbfu(1,5) ,                           &
!         F1CL         F3CL         F4CL

   rtbfu(1,9) , rtbfu(1,10) , rtbfu(1,7) ,                        &
!           F4M1         F4M2        D4CL
   rtbfu(1,8) ,rtbfu(1,11) ,                                      &
!           D4F4         HREC
   propce(1,ipcyf1) , propce(1,ipcyf3) ,                          &
   propce(1,ipcyox) , propce(1,ipcyp1) , propce(1,ipcyp2) ,       &
   propce(1,ipcyin) ,                                             &
   propce(1,ipcy2s) , propce(1,ipcyso) ,                          &
   rtbfu(1,12) )
!         F4S3 pour NOx

! --> Clipping eventuel des fractions massiques

do iel = 1, ncel
  do ice = 1, ngaze
    ipcyce = ipproc(iym1(ice))
    if ( abs(propce(iel,ipcyce)) .lt. epsifl )                    &
         propce(iel,ipcyce) = zero
  enddo
enddo

! MODEL NOx : on y passe pas a la 1ere iter

if ( ieqnox .eq. 1 .and. ntcabs .gt.1) then
  call fucyno                                                     &
  !==========
 ( ncelet , ncel   ,                                              &
   itbfu(1,1) ,                                                   &
!         INTPDF
   rtp    , propce ,                                              &
   fvap   ,   rtbfu(1,1) , rtbfu(1,2) ,                           &
!         F1M           F3M         F4M
  rtbfu(1,3) , rtbfu(1,4) ,rtbfu(1,5) ,                           &
!         F1CL         F3CL         F4CL

   rtbfu(1,9) , rtbfu(1,10) , rtbfu(1,7) ,                        &
!           F4M1         F4M2        D4CL
   rtbfu(1,8) ,rtbfu(1,11) , rtbfu(1,12) ,                        &
!           D4F4         HREC
   propce(1,ipcyf1) , propce(1,ipcyf3) ,                          &
   propce(1,ipcyox) , propce(1,ipcyp1) , propce(1,ipcyp2) ,       &
   propce(1,ipcyin) ,                                             &
   propce(1,ipcy2s) , propce(1,ipcyso) )

else if ( ieqnox .eq. 1 ) then

  write(*,*) ' passage init ',IGHCN1,IGHCN2,IGNOTH
  do iel = 1, ncel
    propce(iel,ipproc(ighcn1)) = 0.d0
    propce(iel,ipproc(ighcn2)) = 0.d0
    propce(iel,ipproc(ignoth)) = 0.d0
  enddo

endif

!===============================================================================
! 3. CALCUL DE LA TEMPERATURE ET DE LA MASSE VOLUMIQUE
!===============================================================================

ipcte1 = ipproc(itemp1)

call futeh1                                                       &
!==========
 ( ncelet , ncel   ,                                              &
   enth,                                                          &
   propce(1,ipcyf1), propce(1,ipcyf3),                            &
   propce(1,ipcyox), propce(1,ipcyp1), propce(1,ipcyp2),          &
   propce(1,ipcyin), propce(1,ipcy2s), propce(1,ipcyso),          &
   propce(1,ipcte1),                                              &
   rtbwo(1,1) , rtbwo(1,2) )

!          TABLEAUX DE TRAVAIL

iphas  = 1
ipcte1 = ipproc(itemp1)
do iel = 1, ncel
  wmolme = propce(iel,ipcyf1) / wmole(ifov)                       &
         + propce(iel,ipcyf3) / wmole(ico )                       &
         + propce(iel,ipcyox) / wmole(io2 )                       &
         + propce(iel,ipcyp1) / wmole(ico2)                       &
         + propce(iel,ipcyp2) / wmole(ih2o)                       &
         + propce(iel,ipcyin) / wmole(in2 )                       &
         + propce(iel,ipcy2s) / wmole(ih2s)                       &
         + propce(iel,ipcyso) / wmole(iso2)

! stockage de la masse molaire du melange

  propce(iel,ipproc(immel)) = 1.d0 / wmolme

! ---- On ne met pas la pression mecanique RTP(IEL,IPR)
!      mais P0

  rom1(iel) = p0 / (wmolme * rr * propce(iel,ipcte1) )
enddo

!===============================================================================
! FORMATS
!----



!----
! FIN
!----

return
end subroutine
