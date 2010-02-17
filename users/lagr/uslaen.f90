!-------------------------------------------------------------------------------

!VERS


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

subroutine uslaen &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , nvlsta ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivarl  , ivarl1 , ivarlm , iflu   , ilpd1  , icla   ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statis , stativ , tracel ,                   &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!    SOUS-PROGRAMME UTILISATEUR (INTERVENTION NON OBLIGATOIRE)

!    POUR LES ECRITURES LISTING ET LE POST-PROCESSING :
!    MOYENNE DES VARIABLES STATISTIQUES VOLUMIQUES LAGRANGIENNES
!    INTERVENTIONS UTILISATEUR POSSIBLES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nvlsta           ! e  ! <-- ! nombre de variables stat. lagrangien           !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
!  (nfml,nprfml    !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac                      !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !                                                !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr                      !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !                                                !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ivarl            !  e ! <-- ! num de la stat (entre 1 et nvlsta)             !
! ivarl1           !  e ! <-- ! num de la stat globale+groupe                  !
!                  !    !     !   (moyenne ou variance)                        !
! ivarlm           !  e ! <-- ! num de la moyenne stat globale+groupe          !
! iflu             !  e ! <-- ! 0 : moyenne de la stat ivarl/ivarl1            !
!                  !    !     ! 1 : variance de la stat ivarl/ivarl1           !
! ilpd1            !  e ! <-- ! "pointeur" sur le poids stat globale           !
!                  !    !     !    ou sur le poids stat d'un groupe            !
! icla             !  e ! <-- !   0 : statistique globale                      !
                   !    ! <-- ! !=0 : statistique pour le groupe icla          !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! tr ! <-- ! coordonnes des noeuds                          !
! (ndim,nnod)      !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! statis(ncelet    ! tr ! <-- ! cumul des statistiques volumiques              !
!   nvlsta)        !    !     !                                                !
! stativ           ! tr ! <-- ! cumul pour les variances des                   !
!(ncelet,          !    !     !    statistiques volumiques                     !
!   nvlsta-1)      !    !     !                                                !
! tracel(ncelet    ! tr ! <-- ! tab reel valeurs cellules post                 !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! rtuser(nrtuse    ! tr ! <-- ! tableau des statistques sur les                !
!                  !    !     ! particules                                     !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "pointe.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "entsor.h"
include "cstnum.h"
include "lagpar.h"
include "lagran.h"
include "parall.h"
include "period.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas  , nvlsta
integer          nideve , nrdeve , nituse , nrtuse
integer          ivarl , ivarl1 , ivarlm , iflu , ilpd1 , icla

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1) , nodfac(lndfac)
integer          ipnfbr(nfabor+1) , nodfbr(lndfbr)
integer          idevel(nideve) , ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision tracel(ncelet)
double precision statis(ncelet,nvlsta)
double precision stativ(ncelet,nvlsta-1)
double precision rdevel(nrdeve) , rtuser(nrtuse)
double precision ra(*)

! Local variables

integer          idebia , idebra , iel
double precision aa

!===============================================================================
!===============================================================================
! 0.  PAR DEFAUT, ON CONSIDERE QUE LE SOUS PROGRAMME CI-DESSOUS CONVIENT
!       A L'UTILISATEUR, C'EST-A-DIRE QUE LA MISE EN OEUVRE DU LAGRANGIEN
!       DECLENCHE LA PRODUCTION DE STATISTIQUES STANDARD.

!     L'UTILISATEUR N'A PAS A MODIFIER LE PRESENT SOUS-PROGRAMME DANS
!       LES CONDITIONS D'UTILISATION STANDARD.
!     DANS LE CAS OU IL SOUHAITE SORTIR DES STATISTIQUES NON PREVUES,
!       IL DOIT INTERVENIR DANS LA RUBRIQUE NUMERO 2.
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1 . ZONE DE STATISTIQUES STANDARD
!===============================================================================

!--> Cas general :
!      Vitesse X des particules : IVARL=ILVX
!      Vitesse Y des particules : IVARL=ILVY
!      Vitesse Z des particules : IVARL=ILVZ
!      Temperature des particules : IVARL=ILTP
!      Diametre des particules : IVARL=ILDP
!      Masse des particules : IVARL=ILMP
!      Temperature des grains de charbon : IVARL=ILHP
!      Masse de Charbon reactif des grains de charbon : IVARL=ILMCH
!      Masse de Coke des grains de charbon : IVARL=ILMCK
!      Diametre du coeur retrecissant des grains de charbon: IVARL=ILMCK
!    sauf fraction volumique (IVARL = ILFV) et
!    somme des poids statistiques (IVARL = ILPD)


if (ivarl.ne.ilfv .and. ivarl.ne.ilpd) then


!-----> Moyenne


  if (iflu.eq.0) then

    do iel = 1, ncel
      if (statis(iel,ilpd1).gt.seuil ) then
        tracel(iel) = statis(iel,ivarl1) / statis(iel,ilpd1)
      else
        tracel(iel) = zero
      endif
    enddo

!-----> Variance

  else

    do iel = 1, ncel
      if ( statis(iel,ilpd1).gt.seuil ) then
        aa = statis(iel,ivarlm)/statis(iel,ilpd1)
        tracel(iel) =  stativ(iel,ivarl1)/statis(iel,ilpd1)       &
                    -( aa * aa )
        tracel(iel) = sqrt( max(zero,tracel(iel)))
      else
        tracel(iel) = zero
      endif
    enddo

  endif

!--> Fraction volumique (ILFV)

else if (ivarl.eq.ilfv) then

!-----> Moyenne

  if (iflu.eq.0) then

    do iel = 1, ncel
      if (statis(iel,ilpd1).gt.seuil .and. npst.gt.0) then
        tracel(iel) = statis(iel,ilfv)                            &
                      / (dble(npst) * volume(iel))
      else if (statis(iel,ilpd1).gt.seuil .and.                   &
                iplas.ge.idstnt                  ) then
        tracel(iel) = statis(iel,ilfv) / volume(iel)
      else
        tracel(iel) = zero
      endif
    enddo

  else

!-----> Variance

    do iel = 1, ncel

      if (statis(iel,ilpd1).gt.seuil .and. npst.gt.0) then

        aa = statis(iel,ivarlm) / (dble(npst) * volume(iel))
        tracel(iel) = stativ(iel,ivarl1)                          &
                 / ( dble(npst) * volume(iel) * volume(iel))      &
                 - aa*aa
        tracel(iel) = sqrt( max(zero,tracel(iel)) )

      else if ( statis(iel,ilpd1).gt.seuil .and.                  &
                iplas.ge.idstnt                  ) then

        aa =  statis(iel,ivarlm) / volume(iel)
        tracel(iel) = stativ(iel,ivarl1) / volume(iel)            &
                         - aa*aa
        tracel(iel) = sqrt( max(zero,tracel(iel)))

      else
        tracel(iel) = zero
      endif

    enddo
  endif

!--> Somme des poids statistiques

else if (ivarl.eq.ilpd) then

  if (iflu .eq.0) then
    do iel = 1, ncel
      tracel(iel) = statis(iel,ivarl1)
    enddo
  else
    write(nfecra,9000) iflu
    do iel = 1, ncel
      tracel(iel) = zero
    enddo
  endif

endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR DANS LE MODULE LAGRANGIEN (USLAEN)   ',/,&
'@    =========                                               ',/,&
'@  IL N''EST PAS POSSIBLE DE CALCULER LA VARIANCE DU POIDS   ',/,&
'@     STATISTIQUE                                            ',/,&
'@                                                            ',/,&
'@  La variance du poids statistique a ete demande            ',/,&
'@    dans USLAEN (IVARL=',   I10,' et IFLU=',  I10,').       ',/,&
'@                                                            ',/,&
'@  L''appel au sous-programme USLAEN doit etre verifie.      ',/,&
'@                                                            ',/,&
'@  Le calcul continue.                                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 2. ZONE D'INTERVENTION UTILISATEUR
!===============================================================================

!    --------------------------------------------------
!    EXEMPLE 1 : STATISTIQUE CALCULEE DANS USLAST.F ET
!                STOCKEE DANS LE TABLEAU STATIS
!    --------------------------------------------------

  if (nvlsts.gt.0) then

    if (ivarl.eq.ilvu(1)) then

!-----> Moyenne pour la concentration massique

      if (iflu.eq.0) then

        do iel = 1, ncel
          if ( statis(iel,ilpd1).gt.seuil .and. npst.gt.0 ) then
            tracel(iel) = statis(iel,ivarl1)                      &
                        / ( dble(npst) *ro0(ilphas) *volume(iel) )
          else if ( statis(iel,ilpd1).gt.seuil .and.              &
                  iplas.ge.idstnt                  ) then
            tracel(iel) = statis(iel,ivarl1)                      &
                        / ( ro0(ilphas) *volume(iel) )
          else
            tracel(iel) = zero
          endif
        enddo

      else

!-----> Variance de la concentration massique

        do iel = 1, ncel
          if (statis(iel,ilpd1).gt.seuil) then
            aa = statis(iel,ivarlm)/statis(iel,ilpd1)
            tracel(iel) = stativ(iel,ivarl1)/statis(iel,ilpd1)    &
                        -( aa * aa)
            tracel(iel) = sqrt( max(zero,tracel(iel)))
          else
            tracel(iel) = zero
          endif

        enddo

      endif

    endif

  endif

!===============================================================================


end subroutine
