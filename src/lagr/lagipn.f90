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

subroutine lagipn &
!================

 ( idbia0 , idbra0 ,                                              &
   ncelet , ncel   ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   npar1  , npar2  ,                                              &
   itepa  ,                                                       &
   ia     ,                                                       &
   rtp    ,                                                       &
   ettp   , tepa   , vagaus ,                                     &
   w1     , w2     , w3     ,                                     &
   ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     Initialisation de la vitesse fluide vu pour les nouvelles
!     particules.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! npar1 ,npar2     ! e  ! <-- ! borne min et max des particules                !
!                  !    !     !    a initialiser                               !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! ia(*)            ! ia ! --- ! main integer work array                        !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! vagaus           ! tr ! --> ! variables aleatoires gaussiennes               !
!(nbpmax,nvgaus    !    !     !                                                !
! w1...w3(ncel)    ! tr ! --- ! tableau de travail                             !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use cstnum
use numvar
use optcal
use entsor
use cstphy
use pointe
use parall
use period
use lagpar
use lagran

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          ncelet , ncel
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          npar1 , npar2
integer          itepa(nbpmax,nivep)
integer          ia(*)

double precision rtp(ncelet,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision vagaus(nbpmax,*)
double precision w1(ncelet) ,  w2(ncelet) ,  w3(ncelet)
double precision ra(*)

! Local variables

integer          idebia , idebra
integer          iel , npt , nomb
double precision tu , d2s3

!===============================================================================

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. INITIALISATION
!===============================================================================

d2s3 = 2.d0 / 3.d0

!===============================================================================
! 2. SIMULATION DES VITESSES TURBULENTES FLUIDES INSTANTANNEES VUES
!    PAR LES PARTICULES SOLIDES LE LONG DE LEUR TRAJECTOIRE.
!===============================================================================

if (idistu.eq.1) then

  if (itytur.eq.2 .or. iturb.eq.50                  &
       .or. iturb.eq.60) then
    do iel = 1,ncel
      w1(iel) = rtp(iel,ik)
    enddo
  else if (itytur.eq.3) then
    do iel = 1,ncel
      w1(iel) = 0.5d0 * ( rtp(iel,ir11)                    &
                        + rtp(iel,ir22)                    &
                        + rtp(iel,ir33) )
    enddo
  else
    write(nfecra,9000) iilagr, idistu, iturb
    call csexit (1)
    !==========
  endif
else
  do iel = 1,ncel
    w1(iel) = 0.d0
  enddo
endif


!---> CALCUL DES TIRAGES ALEATOIRES
!     CALCUL DU TEMPS CARACTERISTIQUE DES PARTICULES
!     remarque : NORMALEN est dans le fichier ZUFALL.F
!     ^^^^^^^^

if (idistu.eq.1) then
  nomb = npar2-npar1+1
  call normalen (nomb,vagaus(npar1,1))
  call normalen (nomb,vagaus(npar1,2))
  call normalen (nomb,vagaus(npar1,3))
else
  do npt = npar1,npar2
    vagaus(npt,1) = 0.d0
    vagaus(npt,2) = 0.d0
    vagaus(npt,3) = 0.d0
  enddo
endif

do npt = npar1,npar2

  iel = itepa(npt,jisor)

  tu = sqrt( d2s3*w1(iel) )

  ettp(npt,juf) = rtp(iel,iu) + vagaus(npt,1)*tu
  ettp(npt,jvf) = rtp(iel,iv) + vagaus(npt,2)*tu
  ettp(npt,jwf) = rtp(iel,iw) + vagaus(npt,3)*tu

enddo

!--------
! FORMATS
!--------

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGIPN)                                    ',/,&
'@                                                            ',/,&
'@    LE MODULE LAGRANGIEN EST INCOMPATIBLE AVEC LE MODELE    ',/,&
'@    DE TURBULENCE SELECTIONNE.                              ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@   Le module Lagrangien a ete active avec IILAGR = ',I10     ,/,&
'@     et la dispersion turbulente est prise en compte        ',/,&
'@                                     avec IDISTU = ',I10     ,/,&
'@   Le modele de turbulence                                  ',/,&
'@     correspond a ITURB  = ',I10                             ,/,&
'@   Or, les seuls traitements de la turbulence compatibles   ',/,&
'@     avec le module Lagrangien et la dispersion turbulente  ',/,&
'@     sont k-epsilon et Rij-epsilon, v2f et k-omega.         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier la valeur de IILAGR et IDISTU dans la subroutine ',/,&
'@  USLAG1 et verifier la valeur de ITURB  dans la subroutine ',/,&
'@  USINI1.                                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

end subroutine
