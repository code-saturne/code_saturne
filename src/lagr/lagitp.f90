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

subroutine lagitp &
!================

 ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   itepa  , ibord  ,                                              &
   propce ,                                                       &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct ,          &
   tsvar  , auxl1  , auxl2  )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     INTEGRATION DES EDS POUR LA TEMPERATURE DES PARTICULES

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
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! ettp             ! tr ! --> ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! ettpa            ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape precedente              !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! taup(nbpmax)     ! tr ! <-- ! temps caracteristique dynamique                !
! tlag(nbpmax)     ! tr ! <-- ! temps caracteristique fluide                   !
! tempct           ! tr ! <-- ! temps caracteristique thermique                !
!  (nbpmax,2)      !    !     !                                                !
! tsvar            ! tr ! <-- ! prediction 1er sous-pas pour la                !
! (nbpmax,nvp1)    !    !     !   variable ivar, utilise pour la               !
!                  !    !     !   correction au 2eme sous-pas                  !
! auxl1(nbpmax)    ! tr ! --- ! tableau de travail                             !
! auxl2(nbpmax)    ! tr ! --- ! tableau de travail                             !
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
use radiat
use mesh

!===============================================================================

implicit none

! Arguments

integer          nbpmax , nvp , nvp1 , nvep , nivep

integer          itepa(nbpmax,nivep) , ibord(nbpmax)

double precision propce(ncelet,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision taup(nbpmax) , tlag(nbpmax,3) , tempct(nbpmax,2)
double precision tsvar(nbpmax,nvp1)
double precision auxl1(nbpmax) , auxl2(nbpmax)

! Local variables

integer          npt , iel
double precision srad

!===============================================================================

!===============================================================================
!     REMPLISSAGE DU TEMPS CARACTERISTIQUE ET DU "PSEUDO SECOND MEMBRE"
!===============================================================================

do npt = 1,nbpart

  if (itepa(npt,jisor).gt.0) then

    auxl1(npt) = tempct(npt,1)

    if (nor.eq.1) then
      auxl2(npt) = ettpa(npt,jtf)
    else
      auxl2(npt) = ettp(npt,jtf)
    endif

  endif

enddo

!===============================================================================
!     PRISE EN COMPTE DU RAYONNEMENT S'IL Y A LIEU
!===============================================================================

if (iirayo.gt.0) then

  do npt = 1,nbpart

    iel = itepa(npt,jisor)

    if (iel.gt.0) then

      if (nor.eq.1) then

        srad = pi *ettpa(npt,jdp) *ettpa(npt,jdp)                 &
                  *tepa(npt,jreps) *(propce(iel,ipproc(ilumin))   &
                        -4.d0 *stephn *ettpa(npt,jtp)**4 )
        auxl2(npt) = ettpa(npt,jtf)                               &
               +auxl1(npt) *srad /ettpa(npt,jcp) /ettpa(npt,jmp)
      else

        srad = pi *ettp(npt,jdp) *ettp(npt,jdp) *tepa(npt,jreps)  &
                *(propce(iel,ipproc(ilumin))                      &
                -4.d0 *stephn *ettp(npt,jtp)**4 )
        auxl2(npt) = ettp(npt,jtf)                                &
                +auxl1(npt) *srad /ettp(npt,jcp) /ettp(npt,jmp)

      endif

    endif

  enddo

endif

!===============================================================================
!     INTEGRATION
!===============================================================================

call lagitg                                                       &
!==========
 ( nbpmax , nvp    , nvp1   ,                                     &
   jtp    ,                                                       &
   itepa(1,jisor)  , ibord  ,                                     &
   ettp   , ettpa  , auxl1  , auxl2  , tsvar  )

!===============================================================================

!----
! FIN
!----

end subroutine
