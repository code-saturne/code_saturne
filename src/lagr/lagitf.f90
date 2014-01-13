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

subroutine lagitf &
!================

 ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   itepa  , ibord  ,                                              &
   rtp    , propce ,                                              &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct , tsvar  , &
   auxl1  )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     INTEGRATION DES EDS POUR LA TEMPERATURE FLUIDE
!      VU PAR LES PARTICULES.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! ibord            ! te ! <-- ! contient le numero de la                       !
!   (nbpmax)       !    !     !   face d'interaction part/frontiere            !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! ettp             ! tr ! --> ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! tlag(nbpmax)     ! tr ! <-- ! temps caracteristique fluide                   !
! tempct           ! tr ! <-- ! temps caracteristique thermique                !
!  (nbpmax,2)      !    !     !                                                !
! tsvar            ! tr ! <-- ! prediction 1er sous-pas pour la                !
! (nbpmax,nvp1)    !    !     !   variable ivar, utilise pour la               !
!                  !    !     !   correction au 2eme sous-pas                  !
! auxl1(nbpmax)    ! tr ! --- ! tableau de travail                             !
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
use cstphy
use cstnum
use optcal
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nbpmax , nvp , nvp1 , nvep , nivep

integer          itepa(nbpmax,nivep) , ibord(nbpmax)

double precision rtp(ncelet,*)
double precision propce(ncelet,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision taup(nbpmax), tlag(nbpmax,3) , tempct(nbpmax,2)
double precision tsvar(nbpmax,nvp1)
double precision auxl1(nbpmax)

! Local variables

integer          npt   , iel   , mode

double precision ct    , aux1  , aux2   , ter1   , ter2
double precision energ , dissip

double precision, allocatable, dimension(:) :: tempf

!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

! Allocate a temporary array
allocate(tempf(ncelet))

! Initialize variables to avoid compiler warnings

dissip = 0.d0
energ = 0.d0

ct = 1.d0

mode = 1

!===============================================================================
! 2. Temperature moyenne Fluide en degres Celsius
!===============================================================================

if ( ippmod(iccoal).ge.0 .or.                                     &
     ippmod(icpl3c).ge.0 .or.                                     &
     ippmod(icfuel).ge.0      ) then

  do iel = 1,ncel
    tempf(iel) = propce(iel,ipproc(itemp1)) - tkelvi
  enddo

else if ( ippmod(icod3p).ge.0 .or.                                &
          ippmod(icoebu).ge.0 .or.                                &
          ippmod(ielarc).ge.0 .or.                                &
          ippmod(ieljou).ge.0      ) then

  do iel = 1,ncel
    tempf(iel) = propce(iel,ipproc(itemp)) - tkelvi
  enddo

else if (itherm.eq.1 .and. itpscl.eq.2) then
  do iel = 1,ncel
    tempf(iel) = rtp(iel,isca(iscalt))
  enddo

else if (itherm.eq.1 .and. itpscl.eq.1) then
  do iel = 1,ncel
    tempf(iel) = rtp(iel,isca(iscalt)) - tkelvi
  enddo

else if (itherm.eq.2) then
  do iel = 1,ncel
    call usthht (mode, rtp(iel,isca(iscalt)), tempf(iel))
    !==========
  enddo
endif

!===============================================================================
! 3. INTEGRATION DE L'EDS SUR LES PARTICULES
!===============================================================================

do npt = 1,nbpart

  if ( itepa(npt,jisor).gt.0 ) then

    iel = itepa(npt,jisor)

    if (itytur.eq.2 .or. itytur.eq.3 .or.           &
        iturb.eq.50 .or. iturb.eq.60) then

      if ( itytur.eq.2 .or. iturb.eq.50 ) then

        energ  = rtp(iel,ik)
        dissip = rtp(iel,iep)

      else if ( itytur.eq.3 ) then

        energ  = 0.5d0 * ( rtp(iel,ir11)                   &
                         + rtp(iel,ir22)                   &
                         + rtp(iel,ir33) )
        dissip = rtp(iel,iep)

      else if (iturb.eq.60) then

        energ  = rtp(iel,ik)
        dissip = cmu*rtp(iel,ik)*rtp(iel,iomg)

      endif

      auxl1(npt) = energ / ( ct*dissip )
      auxl1(npt) = max( auxl1(npt),epzero )

    else

      auxl1(npt) = epzero

    endif

  endif
enddo

if (nor.eq.1) then

  do npt = 1,nbpart

    if (itepa(npt,jisor).gt.0) then

      iel = itepa(npt,jisor)

      aux1 = -dtp/auxl1(npt)
      aux2 = exp(aux1)

      ter1 = ettpa(npt,jtf) * aux2
      ter2 = tempf(iel) * ( 1.d0-aux2 )

      ettp(npt,jtf) = ter1 + ter2

!            Pour le cas NORDRE= 2, on calcule en plus TSVAR pour NOR= 2

      tsvar(npt,jtf) = 0.5d0 * ter1                               &
                     + tempf(iel) * ( -aux2 +(aux2-1.d0) / aux1 )
    endif
  enddo

else if (nor.eq.2) then

  do npt = 1,nbpart

    if ( itepa(npt,jisor).gt.0 .and. ibord(npt).eq.0 ) then

      iel = itepa(npt,jisor)

      aux1 = -dtp/auxl1(npt)
      aux2 = exp(aux1)

      ter1 = 0.5d0 * ettpa(npt,jtf) * aux2
      ter2 = tempf(iel) * (1.d0 - (aux2-1.d0) / aux1)

      ettp(npt,jtf) = tsvar(npt,jtf) + ter1 + ter2
    endif
  enddo
endif

! Free memory
deallocate(tempf)

!===============================================================================

!=======
! FORMAT
!=======

!----
! FIN
!----

end subroutine
