!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine lagsta &
!================

 ( nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  ,                                                       &
   dt     , rtpa   , rtp    , propce , propfb ,                   &
   ettp   , tepa   , statis , stativ ,                            &
   w1     )

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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
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
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itepa(nbpmax,nivep)

double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision statis(ncelet,nvlsta)
double precision stativ(ncelet,nvlsta-1)
double precision w1(ncelet)

! Local variables

integer          npt , nv , iel, nv1, izcl
integer          ilvx1 , ilvy1  , ilvz1  , ilpd1  , ilfv1 , ilts1
integer          iltp1 , ildp1  , ilmp1
integer          ilhp1 , ilmwat1 , ilmch1 , ilmck1 , ildck1
double precision pis6 , concen

!===============================================================================


!===============================================================================

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

! Initialize variables to avoid compiler warnings

ildck1 = 0
ilmck1 = 0
ilmwat1 = 0
ilmch1 = 0

! Memoire


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

if ((npst.lt.2).and.(ilfv.gt.0)) then
   do iel = 1,ncel
      stativ(iel,ilfv) = 0
   enddo
endif

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
!     * Moyenne et variance de la masse d eau
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

    if (iactvx.eq.1) then

       statis(iel,ilvx) = statis(iel,ilvx)                           &
            + tepa(npt,jrpoi) * ettp(npt,jup)
       stativ(iel,ilvx) = stativ(iel,ilvx)                           &
            + tepa(npt,jrpoi) * ettp(npt,jup) * ettp(npt,jup)

    endif

    if (iactvy.eq.1) then

       statis(iel,ilvy) = statis(iel,ilvy)                           &
            + tepa(npt,jrpoi) * ettp(npt,jvp)
       stativ(iel,ilvy) = stativ(iel,ilvy)                           &
            + tepa(npt,jrpoi) * ettp(npt,jvp) * ettp(npt,jvp)

    endif

    if (iactvz.eq.1) then

       statis(iel,ilvz) = statis(iel,ilvz)                           &
                     + tepa(npt,jrpoi) * ettp(npt,jwp)
       stativ(iel,ilvz) = stativ(iel,ilvz)                           &
            + tepa(npt,jrpoi) * ettp(npt,jwp) * ettp(npt,jwp)

    endif

!     * Moyenne et variance du taux de presence

    if (iactfv.eq.1) then

       concen = (ettp(npt,jdp)**3) *pis6

       statis(iel,ilfv) = statis(iel,ilfv)                           &
            + tepa(npt,jrpoi) * concen

!    * La variance du taux de presence n'a de sens qu'en stationnaire

       if (npst.gt.1) then

          stativ(iel,ilfv) = stativ(iel,ilfv)                           &
               + tepa(npt,jrpoi) * concen * concen

       endif

    endif

!     * Moyenne et variance du temps de séjour

    if (iactts.eq.1) then

       statis(iel,ilts) = statis(iel,ilts)                           &
            + tepa(npt,jrpoi) *  tepa(npt,jrtsp)

       stativ(iel,ilts) = stativ(iel,ilts)                           &
            + tepa(npt,jrpoi) *  tepa(npt,jrtsp) *  tepa(npt,jrtsp)

    endif

!     * Somme du poids statistiques associé aux particules (toujours calcule)

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

!     * Moyenne et variance de la masse des particules
!     * Moyenne et variance de la temperature
!     * Moyenne et variance de la masse d eau
!     * Moyenne et variance de la masse de charbon reactif
!     * Moyenne et variance de la masse de coke
!     * Moyenne et variance du diametre du coeur retrecissant

      statis(iel,ilmp) = statis(iel,ilmp)                       &
                         + tepa(npt,jrpoi) * ettp(npt,jmp)
      statis(iel,ilhp)  = statis(iel,ilhp)                        &
                        + tepa(npt,jrpoi) * ettp(npt,jhp)
      statis(iel,ilmwat) = statis(iel,ilmwat)                       &
                        + tepa(npt,jrpoi) * ettp(npt,jmwat)
      statis(iel,ilmch) = statis(iel,ilmch)                       &
                        + tepa(npt,jrpoi) * ettp(npt,jmch)
      statis(iel,ilmck) = statis(iel,ilmck)                       &
                        + tepa(npt,jrpoi) * ettp(npt,jmck)
      statis(iel,ildck) = statis(iel,ildck)                       &
                        + tepa(npt,jrpoi) * tepa(npt,jrdck)

     stativ(iel,ilmp) = stativ(iel,ilmp)                       &
             + tepa(npt,jrpoi)*ettp(npt,jmp)*ettp(npt,jmp)
      stativ(iel,ilhp)  = stativ(iel,ilhp)                        &
          + tepa(npt,jrpoi)*ettp(npt,jhp)*ettp(npt,jhp)
      stativ(iel,ilmwat) = stativ(iel,ilmwat)                       &
          + tepa(npt,jrpoi)*ettp(npt,jmwat)*ettp(npt,jmwat)
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

!     * Moyenne et variance de la masse
!     * Moyenne et variance de la temperature
!     * Moyenne et variance de la masse d eau
!     * Moyenne et variance de la masse de charbon reactif
!     * Moyenne et variance de la masse de coke
!     * Moyenne et variance du diametre du coeur retrecissant

          ilmp1   = ilmp   + izcl*nvlsta
          ilhp1   = ilhp   + izcl*nvlsta
          ilmwat1 = ilmwat + izcl*nvlsta
          ilmch1  = ilmch  + izcl*nvlsta
          ilmck1  = ilmck  + izcl*nvlsta
          ildck1  = ildck  + izcl*nvlsta

          statis(iel,ilmp1)  = statis(iel,ilmp1)                  &
                        + tepa(npt,jrpoi) * ettp(npt,jmp)
          statis(iel,ilhp1)  = statis(iel,ilhp1)                  &
                        + tepa(npt,jrpoi) * ettp(npt,jhp)
          statis(iel,ilmwat1) = statis(iel,ilmwat1)                 &
                        + tepa(npt,jrpoi) * ettp(npt,jmwat)
          statis(iel,ilmch1) = statis(iel,ilmch1)                 &
                        + tepa(npt,jrpoi) * ettp(npt,jmch)
          statis(iel,ilmck1) = statis(iel,ilmck1)                 &
                        + tepa(npt,jrpoi) * ettp(npt,jmck)
          statis(iel,ildck1) = statis(iel,ildck1)                 &
                        + tepa(npt,jrpoi) * tepa(npt,jrdck)

          ilmp1   = ilmp   + izcl*(nvlsta-1)
          ilhp1   = ilhp   + izcl*(nvlsta-1)
          ilmwat1 = ilmwat + izcl*(nvlsta-1)
          ilmch1  = ilmch  + izcl*(nvlsta-1)
          ilmck1  = ilmck  + izcl*(nvlsta-1)
          ildck1  = ildck  + izcl*(nvlsta-1)

          stativ(iel,ilmp1)  = stativ(iel,ilmp1)                  &
          + tepa(npt,jrpoi)*ettp(npt,jmp)*ettp(npt,jmp)
          stativ(iel,ilhp1)  = stativ(iel,ilhp1)                  &
          + tepa(npt,jrpoi)*ettp(npt,jhp)*ettp(npt,jhp)
          stativ(iel,ilmwat1) = stativ(iel,ilmwat1)                 &
          + tepa(npt,jrpoi)*ettp(npt,jmwat)*ettp(npt,jmwat)
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
