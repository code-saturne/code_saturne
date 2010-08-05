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

subroutine lagsta &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   itepa  ,                                                       &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   ettp   , tepa   , statis , stativ ,                            &
   w1     ,                                                       &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

!       SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!       -----------------------------------

!    CALCUL DES STATISTIQUES SUR LES PARTICULES


!  ISTTIO = 0 : calcul instationnaire pour le lagrangien
!         = 1 : calcul stationnaire   pour le lagrangien

!  ISTALA : calcul statistiques       si  >= 1 sinon pas de stat

!  ISUIST : suite calcul statistiques si  >= 1 sinon pas de stat

!  IDSTNT : Numero du pas de temps pour debut statistque

!  NSTIST : iteration Lagrangienne du debut calcul stationnaire

!  NPST   : nombre d'iterations de calcul de stats stationnaires

!  NPSTT  : nombre d'iterations total des stats depuis le debut
!             du calcul, partie instationnaire comprise
!             (Attention : uniquement pour affichage listing,
!              ne pas faire de test dessus, preferer IDSTNT)

!  TSTAT  : temps physique de calcul des statistiques stationnaires

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
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
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
! volume(ncelet    ! tr ! <-- ! volume d'un des ncelet elements                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! statis(ncelet    ! tr ! --> ! cumul des statistiques volumiques              !
!    nvlsta)       !    !     !                                                !
! stativ           ! tr ! <-- ! cumul pour les variances des                   !
!(ncelet,          !    !     !    statistiques volumiques                     !
!   nvlsta-1)      !    !     !                                                !
! w1(ncelet)       ! tr ! --- ! tableau de travail                             !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
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
use numvar
use cstnum
use optcal
use entsor
use lagpar
use lagran
use cstphy
use ppppar
use ppthch

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nideve , nrdeve , nituse , nrtuse
integer          itepa(nbpmax,nivep)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision statis(ncelet,nvlsta)
double precision stativ(ncelet,nvlsta-1)
double precision w1(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse)
double precision ra(*)

! Local variables

integer          idebia , idebra
integer          npt , nv , iel, nv1, izcl
integer          ilvx1 , ilvy1  , ilvz1  , ilpd1  , ilfv1 , ilts1
integer          iltp1 , ildp1  , ilmp1
integer          ilhp1 , ilmch1 , ilmck1 , ildck1
double precision pis6 , concen

!===============================================================================


!===============================================================================

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

! Initialize variables to avoid compiler warnings

ildck1 = 0
ilmck1 = 0
ilmch1 = 0

! Memoire

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. INCREMENTATION DES COMPTEURS ET INITIALISATION
!===============================================================================

!      ISTTIO = 0 : calcul instationnaire pour le lagrangien
!             = 1 : calcul stationnaire   pour le lagrangien
!      NSTIST : iter de depart pour le debut calcul stationnaire
!      NPST   : nombre de pas de temps pour le cumul des stats
!               stationnaires
!      NPSTT  : nombre de pas de temps total des stats depuis le debut
!               du calcul, partie instationnaire comprise (pour listing)
!      TSTAT  : temps physique de calcul des statistiques stationnaires

!-> Si on est en instationnaire, ou si le debut des stat stationnaires
!   n'est pas encore enclenchee, on remet les stat a zero a chaque
!   pas de temps

if (isttio.eq.0 .or. (isttio.eq.1 .and. iplas.le.nstist)) then

  npst  = 0
  tstat = 0.d0

! Statistiques globales
  do nv = 1,nvlsta
    do iel = 1,ncel
      statis(iel,nv) = 0.d0
    enddo
  enddo

  do nv = 1,nvlsta-1
    do iel = 1,ncel
      stativ(iel,nv) = 0.d0
    enddo
  enddo

! Statistiques par groupe
  if (nbclst.gt.0) then

    do nv = 1, nvlsta*nbclst
      nv1 = nvlsta + nv
      do iel = 1,ncel
        statis(iel,nv1) = 0.d0
      enddo
    enddo

    do nv = 1, (nvlsta-1)*nbclst
      nv1 = (nvlsta-1) + nv
      do iel = 1,ncel
        stativ(iel,nv1) = 0.d0
      enddo
    enddo

  endif

endif

npst  = npst  + 1
tstat = tstat + dtp

npstt = iplas - idstnt + 1

!===============================================================================
! 2 - CALCUL DES STATISTIQUES PARTICULAIRES
!===============================================================================

!     * Moyenne et variance des composantes de la vitesse
!     * Moyenne et variance du taux de presence
!           (i.e. concentration volumique)
!     * Moyenne et variance du temps de séjour
!     * Somme du poids statistiques associé aux particules
!          (i.e    . nombre de particules par cellules)

!     * Moyenne et variance de la temperature
!     * Moyenne et variance du diametre
!     * Moyenne et variance de la masse

!     * Moyenne et variance de la temperature
!     * Moyenne et variance de la masse de charbon reactif
!     * Moyenne et variance de la masse de coke
!     * Moyenne et variance du diametre du coeur retrecissant


pis6 = pi / 6.d0

!===============================================================================
! 2.1 -  STATISTIQUES GLOBALES
!===============================================================================

do npt = 1,nbpart

  if (itepa(npt,jisor).gt.0) then

!     * Moyenne et variance des composantes de la vitesse

    iel = itepa(npt,jisor)

    statis(iel,ilvx) = statis(iel,ilvx)                           &
                     + tepa(npt,jrpoi) * ettp(npt,jup)
    statis(iel,ilvy) = statis(iel,ilvy)                           &
                     + tepa(npt,jrpoi) * ettp(npt,jvp)
    statis(iel,ilvz) = statis(iel,ilvz)                           &
                     + tepa(npt,jrpoi) * ettp(npt,jwp)

    stativ(iel,ilvx) = stativ(iel,ilvx)                           &
     + tepa(npt,jrpoi) * ettp(npt,jup) * ettp(npt,jup)

    stativ(iel,ilvy) = stativ(iel,ilvy)                           &
     + tepa(npt,jrpoi) * ettp(npt,jvp) * ettp(npt,jvp)
    stativ(iel,ilvz) = stativ(iel,ilvz)                           &
     + tepa(npt,jrpoi) * ettp(npt,jwp) * ettp(npt,jwp)

!     * Moyenne et variance du taux de presence

    concen = (ettp(npt,jdp)**3) *pis6

    statis(iel,ilfv) = statis(iel,ilfv)                           &
      + tepa(npt,jrpoi) * concen

    stativ(iel,ilfv) = stativ(iel,ilfv)                           &
      + tepa(npt,jrpoi) * concen * concen

!     * Moyenne et variance du temps de séjour

    statis(iel,ilts) = statis(iel,ilts)                           &
      + tepa(npt,jrpoi) *  tepa(npt,jrtsp)

    stativ(iel,ilts) = stativ(iel,ilts)                           &
      + tepa(npt,jrpoi) *  tepa(npt,jrtsp) *  tepa(npt,jrtsp)

!     * Somme du poids statistiques associé aux particules

    statis(iel,ilpd) = statis(iel,ilpd) + tepa(npt,jrpoi)

    if (iphyla.eq.1) then

!     * Moyenne et variance de la temperature

      if ( itpvar .eq. 1 ) then

        statis(iel,iltp) = statis(iel,iltp)                       &
                         + tepa(npt,jrpoi) * ettp(npt,jtp)

        stativ(iel,iltp) = stativ(iel,iltp)                       &
             + tepa(npt,jrpoi)*ettp(npt,jtp)*ettp(npt,jtp)

      endif

!     * Moyenne et variance du diametre

      if ( idpvar .eq. 1 ) then

        statis(iel,ildp) = statis(iel,ildp)                       &
                         + tepa(npt,jrpoi) * ettp(npt,jdp)

        stativ(iel,ildp) = stativ(iel,ildp)                       &
             + tepa(npt,jrpoi)*ettp(npt,jdp)*ettp(npt,jdp)

      endif

!     * Moyenne et variance de la masse

      if ( impvar .eq. 1 ) then

        statis(iel,ilmp) = statis(iel,ilmp)                       &
                         + tepa(npt,jrpoi) * ettp(npt,jmp)

        stativ(iel,ilmp) = stativ(iel,ilmp)                       &
             + tepa(npt,jrpoi)*ettp(npt,jmp)*ettp(npt,jmp)

      endif

    else if ( iphyla .eq. 2 ) then

!     * Moyenne et variance de la temperature
!     * Moyenne et variance de la masse de charbon reactif
!     * Moyenne et variance de la masse de coke
!     * Moyenne et variance du diametre du coeur retrecissant

      statis(iel,ilhp)  = statis(iel,ilhp)                        &
                        + tepa(npt,jrpoi) * ettp(npt,jhp)
      statis(iel,ilmch) = statis(iel,ilmch)                       &
                        + tepa(npt,jrpoi) * ettp(npt,jmch)
      statis(iel,ilmck) = statis(iel,ilmck)                       &
                        + tepa(npt,jrpoi) * ettp(npt,jmck)
      statis(iel,ildck) = statis(iel,ildck)                       &
                        + tepa(npt,jrpoi) * tepa(npt,jrdck)

      stativ(iel,ilhp)  = stativ(iel,ilhp)                        &
          + tepa(npt,jrpoi)*ettp(npt,jhp)*ettp(npt,jhp)
      stativ(iel,ilmch) = stativ(iel,ilmch)                       &
          + tepa(npt,jrpoi)*ettp(npt,jmch)*ettp(npt,jmch)
      stativ(iel,ilmck) = stativ(iel,ilmck)                       &
          + tepa(npt,jrpoi)*ettp(npt,jmck)*ettp(npt,jmck)
      stativ(iel,ildck) = stativ(iel,ildck)                       &
          + tepa(npt,jrpoi)*tepa(npt,jrdck)*tepa(npt,jrdck)

    endif

  endif

enddo

!===============================================================================
! 2.2 -  STATISTIQUES PAR GROUPE
!===============================================================================

if (nbclst.gt.0) then

  do izcl = 1,nbclst

    do npt = 1,nbpart

      if (itepa(npt,jisor).gt.0 .and.                             &
                         itepa(npt,jclst).eq. izcl ) then

!     * Moyenne et variance des composantes de la vitesse

        iel = itepa(npt,jisor)

        ilvx1 = ilvx + izcl*nvlsta
        ilvy1 = ilvy + izcl*nvlsta
        ilvz1 = ilvz + izcl*nvlsta

        statis(iel,ilvx1) = statis(iel,ilvx1)                     &
                       + tepa(npt,jrpoi) * ettp(npt,jup)
        statis(iel,ilvy1) = statis(iel,ilvy1)                     &
                       + tepa(npt,jrpoi) * ettp(npt,jvp)
        statis(iel,ilvz1) = statis(iel,ilvz1)                     &
                       + tepa(npt,jrpoi)*ettp(npt,jwp)

        ilvx1 = ilvx + izcl*(nvlsta-1)
        ilvy1 = ilvy + izcl*(nvlsta-1)
        ilvz1 = ilvz + izcl*(nvlsta-1)

        stativ(iel,ilvx1) = stativ(iel,ilvx1)                     &
       + tepa(npt,jrpoi) * ettp(npt,jup) * ettp(npt,jup)
        stativ(iel,ilvy1) = stativ(iel,ilvy1)                     &
       + tepa(npt,jrpoi) * ettp(npt,jvp) * ettp(npt,jvp)
        stativ(iel,ilvz1) = stativ(iel,ilvz1)                     &
       + tepa(npt,jrpoi) * ettp(npt,jwp) * ettp(npt,jwp)

!     * Moyenne et variance du taux de presence


        concen = (ettp(npt,jdp)**3) *pis6

        ilfv1 = ilfv+izcl*nvlsta
        statis(iel,ilfv1) = statis(iel,ilfv1)                     &
                          +tepa(npt,jrpoi)*concen

        ilfv1 = ilfv+izcl*(nvlsta-1)
        stativ(iel,ilfv1) = stativ(iel,ilfv1)                     &
                          +tepa(npt,jrpoi)*concen*concen

!     * Moyenne et variance du temps de séjour

        ilts1 = ilts+izcl*nvlsta
        statis(iel,ilts1) = statis(iel,ilts1)                     &
          + tepa(npt,jrpoi) *  tepa(npt,jrtsp)

        ilts1 = ilts+izcl*(nvlsta-1)
        stativ(iel,ilts1) = stativ(iel,ilts1)                     &
          + tepa(npt,jrpoi)*tepa(npt,jrtsp)*tepa(npt,jrtsp)

!     * Somme du poids statistiques associé aux particules

        ilpd1 = ilpd+izcl*nvlsta
        statis(iel,ilpd1) = statis(iel,ilpd1) + tepa(npt,jrpoi)

       if ( iphyla .eq. 1 ) then

!     * Moyenne et variance de la temperature

          if ( itpvar .eq. 1 ) then

            iltp1 = iltp+izcl*nvlsta
            statis(iel,iltp1) = statis(iel,iltp1)                 &
                         + tepa(npt,jrpoi) * ettp(npt,jtp)

            iltp1 = iltp+izcl*(nvlsta-1)
            stativ(iel,iltp1) = stativ(iel,iltp1)                 &
             + tepa(npt,jrpoi)*ettp(npt,jtp)*ettp(npt,jtp)

          endif

!     * Moyenne et variance du diametre

          if ( idpvar .eq. 1 ) then

            ildp1 = ildp+izcl*nvlsta
            statis(iel,ildp1) = statis(iel,ildp1)                 &
                         + tepa(npt,jrpoi) * ettp(npt,jdp)

            ildp1 = ildp+izcl*(nvlsta-1)
            stativ(iel,ildp1) = stativ(iel,ildp1)                 &
             + tepa(npt,jrpoi)*ettp(npt,jdp)*ettp(npt,jdp)

          endif

!     * Moyenne et variance de la masse

          if ( impvar .eq. 1 ) then

            ilmp1 = ilmp+izcl*nvlsta
            statis(iel,ilmp1) = statis(iel,ilmp1)                 &
                         + tepa(npt,jrpoi) * ettp(npt,jmp)

            ilmp1 = ilmp+izcl*(nvlsta-1)
            stativ(iel,ilmp1) = stativ(iel,ilmp1)                 &
             + tepa(npt,jrpoi)*ettp(npt,jmp)*ettp(npt,jmp)

          endif

        else if ( iphyla .eq. 2 ) then

!     * Moyenne et variance de la temperature
!     * Moyenne et variance de la masse de charbon reactif
!     * Moyenne et variance de la masse de coke
!     * Moyenne et variance du diametre du coeur retrecissant

          ilhp1  = ilhp  +izcl*nvlsta
          ilmch1 = ilmch1+izcl*nvlsta
          ilmck1 = ilmck1+izcl*nvlsta
          ildck1 = ildck1+izcl*nvlsta

          statis(iel,ilhp1)  = statis(iel,ilhp1)                  &
                        + tepa(npt,jrpoi) * ettp(npt,jhp)
          statis(iel,ilmch1) = statis(iel,ilmch1)                 &
                        + tepa(npt,jrpoi) * ettp(npt,jmch)
          statis(iel,ilmck1) = statis(iel,ilmck1)                 &
                        + tepa(npt,jrpoi) * ettp(npt,jmck)
          statis(iel,ildck1) = statis(iel,ildck1)                 &
                        + tepa(npt,jrpoi) * tepa(npt,jrdck)

          ilhp1  = ilhp  +izcl*(nvlsta-1)
          ilmch1 = ilmch1+izcl*(nvlsta-1)
          ilmck1 = ilmck1+izcl*(nvlsta-1)
          ildck1 = ildck1+izcl*(nvlsta-1)

          stativ(iel,ilhp1)  = stativ(iel,ilhp1)                  &
          + tepa(npt,jrpoi)*ettp(npt,jhp)*ettp(npt,jhp)
          stativ(iel,ilmch1) = stativ(iel,ilmch1)                 &
          + tepa(npt,jrpoi)*ettp(npt,jmch)*ettp(npt,jmch)
          stativ(iel,ilmck1) = stativ(iel,ilmck1)                 &
          + tepa(npt,jrpoi)*ettp(npt,jmck)*ettp(npt,jmck)
          stativ(iel,ildck1) = stativ(iel,ildck1)                 &
          + tepa(npt,jrpoi)*tepa(npt,jrdck)*tepa(npt,jrdck)

        endif

      endif

    enddo

  enddo

endif

!===============================================================================

!====
! FIN
!====

end subroutine
