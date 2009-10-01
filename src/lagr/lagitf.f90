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

subroutine lagitf &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  , ibord  , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , volume ,          &
   dt     , rtp    , propce , propfa , propfb ,                   &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct , tsvar  , &
   auxl1  , auxl2  , tempf  ,                                     &
   ra     )

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
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ndim             ! e  ! <-- ! dimension de l'espace                          !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nfml             ! e  ! <-- ! nombre de familles d entites                   !
! nprfml           ! e  ! <-- ! nombre de proprietese des familles             !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! ibord            ! te ! <-- ! contient le numero de la                       !
!   (nbpmax)       !    !     !   face d'interaction part/frontiere            !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
! (ndim,ncelet     !    !     !                                                !
! surfac           ! tr ! <-- ! vecteur surface des faces internes             !
! (ndim,nfac)      !    !     !                                                !
! surfbo           ! tr ! <-- ! vecteur surface des faces de bord              !
! (ndim,nfabor)    !    !     !                                                !
! cdgfac           ! tr ! <-- ! centre de gravite des faces internes           !
! (ndim,nfac)      !    !     !                                                !
! cdgfbo           ! tr ! <-- ! centre de gravite des faces de bord            !
! (ndim,nfabor)    !    !     !                                                !
! xyznod           ! tr ! <-- ! coordonnes des noeuds                          !
! (ndim,nnod)      !    !     !                                                !
! volume(ncelet    ! tr ! <-- ! volume d'un des ncelet elements                !
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
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
! tempf(ncelet)    ! tr ! --- ! tableau de travail                             !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "numvar.h"
include "cstphy.h"
include "cstnum.h"
include "optcal.h"
include "entsor.h"
include "lagpar.h"
include "lagran.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp , nvp1 , nvep , nivep
integer          ntersl , nvlsta , nvisbr
integer          itepa(nbpmax,nivep) , ibord(nbpmax)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision volume(ncelet)
double precision dt(ncelet) , rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision taup(nbpmax) , tlag(nbpmax,3) , tempct(nbpmax,2)
double precision tsvar(nbpmax,nvp1)
double precision auxl1(nbpmax) , auxl2(nbpmax)
double precision tempf(ncelet)
double precision ra(*)

! VARIABLES LOCALES

integer          npt   , iel   , iphas  , mode
double precision ct    , aux1  , aux2   , ter1   , ter2
double precision energ , dissip

!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

iphas = ilphas

ct = 1.d0

mode = 1

!===============================================================================
! 2. Temperature moyenne Fluide en degres Celsius
!===============================================================================

if ( ippmod(icp3pl).ge.0 .or.                                     &
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

else if ( iscsth(iscalt(iphas)).eq.-1 ) then
  do iel = 1,ncel
    tempf(iel) = rtp(iel,isca(iscalt(iphas)))
  enddo

else if ( iscsth(iscalt(iphas)).eq.1 ) then
  do iel = 1,ncel
    tempf(iel) = rtp(iel,isca(iscalt(iphas))) - tkelvi
  enddo

else if ( iscsth(iscalt(iphas)).eq.2 ) then
  do iel = 1,ncel
    call usthht (mode, rtp(iel,isca(iscalt(iphas))), tempf(iel))
    !==========
  enddo
endif

!===============================================================================
! 3. INTEGRATION DE L'EDS SUR LES PARTICULES
!===============================================================================

do npt = 1,nbpart

  if ( itepa(npt,jisor).gt.0 ) then

    iel = itepa(npt,jisor)

    if (itytur(iphas).eq.2 .or. itytur(iphas).eq.3 .or.           &
         iturb(iphas).eq.50 .or. iturb(iphas).eq.60 ) then

      if ( itytur(iphas).eq.2 .or. iturb(iphas).eq.50 ) then

        energ  = rtp(iel,ik(iphas))
        dissip = rtp(iel,iep(iphas))

      else if ( itytur(iphas).eq.3 ) then

        energ  = 0.5d0 * ( rtp(iel,ir11(iphas))                   &
                         + rtp(iel,ir22(iphas))                   &
                         + rtp(iel,ir33(iphas)) )
        dissip = rtp(iel,iep(iphas))

      else if (iturb(iphas).eq.60) then

        energ  = rtp(iel,ik(iphas))
        dissip = cmu*rtp(iel,ik(iphas))*rtp(iel,iomg(iphas))

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

!===============================================================================

!=======
! FORMAT
!=======

!----
! FIN
!----

end subroutine
