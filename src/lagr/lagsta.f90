!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

 ( nvlsta ,                                                       &
   statis , stativ )

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
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! statis(ncelet    ! tr ! --> ! cumul des statistiques volumiques              !
!    nvlsta)       !    !     !                                                !
! stativ           ! tr ! <-- ! cumul pour les variances des                   !
!(ncelet,          !    !     !    statistiques volumiques                     !
!   nvlsta-1)      !    !     !                                                !
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

integer          nvlsta

double precision statis(ncelet,nvlsta)
double precision stativ(ncelet,nvlsta-1)

! Local variables

integer          npt , nv , iel, ilayer, nv1, izcl
integer          ilvx1 , ilvy1  , ilvz1  , ilpd1  , ilfv1 , ilts1
integer          iltp1 , ildp1  , ilmp1
integer          ilhp1(nlayer) , ilmwat1 , ilmch1(nlayer) , ilmck1(nlayer) , ildck1
double precision pis6 , concen

!===============================================================================


!===============================================================================

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================

! Initialize variables to avoid compiler warnings
do ilayer=1,nlayer
  ilhp1(ilayer) = 0
  ilmch1(ilayer) = 0
  ilmck1(ilayer) = 0
enddo
ilmwat1 = 0
ildck1 = 0

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

  if (ipepa(jisor,npt).gt.0) then

!     * Moyenne et variance des composantes de la vitesse

    iel = ipepa(jisor,npt)

    if (iactvx.eq.1) then

       statis(iel,ilvx) = statis(iel,ilvx)                           &
            + pepa(jrpoi,npt) * eptp(jup,npt)
       stativ(iel,ilvx) = stativ(iel,ilvx)                           &
            + pepa(jrpoi,npt) * eptp(jup,npt) * eptp(jup,npt)

    endif

    if (iactvy.eq.1) then

       statis(iel,ilvy) = statis(iel,ilvy)                           &
            + pepa(jrpoi,npt) * eptp(jvp,npt)
       stativ(iel,ilvy) = stativ(iel,ilvy)                           &
            + pepa(jrpoi,npt) * eptp(jvp,npt) * eptp(jvp,npt)

    endif

    if (iactvz.eq.1) then

       statis(iel,ilvz) = statis(iel,ilvz)                           &
                     + pepa(jrpoi,npt) * eptp(jwp,npt)
       stativ(iel,ilvz) = stativ(iel,ilvz)                           &
            + pepa(jrpoi,npt) * eptp(jwp,npt) * eptp(jwp,npt)

    endif

!     * Moyenne et variance du taux de presence

    if (iactfv.eq.1) then

       concen = (eptp(jdp,npt)**3) *pis6

       statis(iel,ilfv) = statis(iel,ilfv)                           &
            + pepa(jrpoi,npt) * concen

!    * La variance du taux de presence n'a de sens qu'en stationnaire

       if (npst.gt.1) then

          stativ(iel,ilfv) = stativ(iel,ilfv)                           &
               + pepa(jrpoi,npt) * concen * concen

       endif

    endif

!     * Moyenne et variance du temps de séjour

    if (iactts.eq.1) then

       statis(iel,ilts) = statis(iel,ilts)                           &
            + pepa(jrpoi,npt) *  pepa(jrtsp,npt)

       stativ(iel,ilts) = stativ(iel,ilts)                           &
            + pepa(jrpoi,npt) *  pepa(jrtsp,npt) *  pepa(jrtsp,npt)

    endif

!     * Somme du poids statistiques associé aux particules (toujours calcule)

    statis(iel,ilpd) = statis(iel,ilpd) + pepa(jrpoi,npt)

    if (iphyla.eq.1) then

!     * Moyenne et variance de la temperature

      if ( itpvar .eq. 1 ) then

        statis(iel,iltp) = statis(iel,iltp)                       &
                         + pepa(jrpoi,npt) * eptp(jtp,npt)

        stativ(iel,iltp) = stativ(iel,iltp)                       &
             + pepa(jrpoi,npt)*eptp(jtp,npt)*eptp(jtp,npt)

      endif

!     * Moyenne et variance du diametre

      if ( idpvar .eq. 1 ) then

        statis(iel,ildp) = statis(iel,ildp)                       &
                         + pepa(jrpoi,npt) * eptp(jdp,npt)

        stativ(iel,ildp) = stativ(iel,ildp)                       &
             + pepa(jrpoi,npt)*eptp(jdp,npt)*eptp(jdp,npt)

      endif

!     * Moyenne et variance de la masse

      if ( impvar .eq. 1 ) then

        statis(iel,ilmp) = statis(iel,ilmp)                       &
                         + pepa(jrpoi,npt) * eptp(jmp,npt)

        stativ(iel,ilmp) = stativ(iel,ilmp)                       &
             + pepa(jrpoi,npt)*eptp(jmp,npt)*eptp(jmp,npt)

      endif

    else if ( iphyla .eq. 2 ) then

!     * Moyenne et variance de la masse des particules
!     * Moyenne et variance de la temperature
!     * Moyenne et variance de la masse d eau
!     * Moyenne et variance de la masse de charbon reactif
!     * Moyenne et variance de la masse de coke
!     * Moyenne et variance du diametre du coeur retrecissant

      statis(iel,ilmp) = statis(iel,ilmp)                         &
                         + pepa(jrpoi,npt) * eptp(jmp,npt)
      do ilayer=1,nlayer
        statis(iel,ilhp(ilayer))  = statis(iel,ilhp(ilayer))      &
              + pepa(jrpoi,npt) * eptp(jhp(ilayer),npt)
        statis(iel,ilmch(ilayer)) = statis(iel,ilmch(ilayer))     &
              + pepa(jrpoi,npt) * eptp(jmch(ilayer),npt)
        statis(iel,ilmck(ilayer)) = statis(iel,ilmck(ilayer))     &
              + pepa(jrpoi,npt) * eptp(jmck(ilayer),npt)
      enddo
      statis(iel,ilmwat) = statis(iel,ilmwat)                     &
                        + pepa(jrpoi,npt) * eptp(jmwat,npt)
      statis(iel,ildck) = statis(iel,ildck)                       &
                        + pepa(jrpoi,npt) * pepa(jrdck,npt)

     stativ(iel,ilmp) = stativ(iel,ilmp)                          &
             + pepa(jrpoi,npt)*eptp(jmp,npt)*eptp(jmp,npt)
      do ilayer=1,nlayer
        stativ(iel,ilhp(ilayer))  = stativ(iel,ilhp(ilayer))      &
        + pepa(jrpoi,npt)*eptp(jhp(ilayer),npt)*eptp(jhp(ilayer),npt)
        stativ(iel,ilmch(ilayer)) = stativ(iel,ilmch(ilayer))     &
        + pepa(jrpoi,npt)*eptp(jmch(ilayer),npt)*eptp(jmch(ilayer),npt)
        stativ(iel,ilmck(ilayer)) = stativ(iel,ilmck(ilayer))     &
        + pepa(jrpoi,npt)*eptp(jmck(ilayer),npt)*eptp(jmck(ilayer),npt)
      enddo
      stativ(iel,ilmwat) = stativ(iel,ilmwat)                     &
          + pepa(jrpoi,npt)*eptp(jmwat,npt)*eptp(jmwat,npt)
      stativ(iel,ildck) = stativ(iel,ildck)                       &
          + pepa(jrpoi,npt)*pepa(jrdck,npt)*pepa(jrdck,npt)

    endif

  endif

enddo

!===============================================================================
! 2.2 -  STATISTIQUES PAR GROUPE
!===============================================================================

if (nbclst.gt.0) then

  do izcl = 1,nbclst

    do npt = 1,nbpart

      if (ipepa(jisor,npt).gt.0 .and.                             &
                         ipepa(jclst,npt).eq. izcl ) then

!     * Moyenne et variance des composantes de la vitesse

        iel = ipepa(jisor,npt)

        ilvx1 = ilvx + izcl*nvlsta
        ilvy1 = ilvy + izcl*nvlsta
        ilvz1 = ilvz + izcl*nvlsta

        statis(iel,ilvx1) = statis(iel,ilvx1)                     &
                       + pepa(jrpoi,npt) * eptp(jup,npt)
        statis(iel,ilvy1) = statis(iel,ilvy1)                     &
                       + pepa(jrpoi,npt) * eptp(jvp,npt)
        statis(iel,ilvz1) = statis(iel,ilvz1)                     &
                       + pepa(jrpoi,npt)*eptp(jwp,npt)

        ilvx1 = ilvx + izcl*(nvlsta-1)
        ilvy1 = ilvy + izcl*(nvlsta-1)
        ilvz1 = ilvz + izcl*(nvlsta-1)

        stativ(iel,ilvx1) = stativ(iel,ilvx1)                     &
       + pepa(jrpoi,npt) * eptp(jup,npt) * eptp(jup,npt)
        stativ(iel,ilvy1) = stativ(iel,ilvy1)                     &
       + pepa(jrpoi,npt) * eptp(jvp,npt) * eptp(jvp,npt)
        stativ(iel,ilvz1) = stativ(iel,ilvz1)                     &
       + pepa(jrpoi,npt) * eptp(jwp,npt) * eptp(jwp,npt)

!     * Moyenne et variance du taux de presence


        concen = (eptp(jdp,npt)**3) *pis6

        ilfv1 = ilfv+izcl*nvlsta
        statis(iel,ilfv1) = statis(iel,ilfv1)                     &
                          +pepa(jrpoi,npt)*concen

        ilfv1 = ilfv+izcl*(nvlsta-1)
        stativ(iel,ilfv1) = stativ(iel,ilfv1)                     &
                          +pepa(jrpoi,npt)*concen*concen

!     * Moyenne et variance du temps de séjour

        ilts1 = ilts+izcl*nvlsta
        statis(iel,ilts1) = statis(iel,ilts1)                     &
          + pepa(jrpoi,npt) *  pepa(jrtsp,npt)

        ilts1 = ilts+izcl*(nvlsta-1)
        stativ(iel,ilts1) = stativ(iel,ilts1)                     &
          + pepa(jrpoi,npt)*pepa(jrtsp,npt)*pepa(jrtsp,npt)

!     * Somme du poids statistiques associé aux particules

        ilpd1 = ilpd+izcl*nvlsta
        statis(iel,ilpd1) = statis(iel,ilpd1) + pepa(jrpoi,npt)

       if ( iphyla .eq. 1 ) then

!     * Moyenne et variance de la temperature

          if ( itpvar .eq. 1 ) then

            iltp1 = iltp+izcl*nvlsta
            statis(iel,iltp1) = statis(iel,iltp1)                 &
                         + pepa(jrpoi,npt) * eptp(jtp,npt)

            iltp1 = iltp+izcl*(nvlsta-1)
            stativ(iel,iltp1) = stativ(iel,iltp1)                 &
             + pepa(jrpoi,npt)*eptp(jtp,npt)*eptp(jtp,npt)

          endif

!     * Moyenne et variance du diametre

          if ( idpvar .eq. 1 ) then

            ildp1 = ildp+izcl*nvlsta
            statis(iel,ildp1) = statis(iel,ildp1)                 &
                         + pepa(jrpoi,npt) * eptp(jdp,npt)

            ildp1 = ildp+izcl*(nvlsta-1)
            stativ(iel,ildp1) = stativ(iel,ildp1)                 &
             + pepa(jrpoi,npt)*eptp(jdp,npt)*eptp(jdp,npt)

          endif

!     * Moyenne et variance de la masse

          if ( impvar .eq. 1 ) then

            ilmp1 = ilmp+izcl*nvlsta
            statis(iel,ilmp1) = statis(iel,ilmp1)                 &
                         + pepa(jrpoi,npt) * eptp(jmp,npt)

            ilmp1 = ilmp+izcl*(nvlsta-1)
            stativ(iel,ilmp1) = stativ(iel,ilmp1)                 &
             + pepa(jrpoi,npt)*eptp(jmp,npt)*eptp(jmp,npt)

          endif

        else if ( iphyla .eq. 2 ) then

!     * Moyenne et variance de la masse
!     * Moyenne et variance de la temperature
!     * Moyenne et variance de la masse d eau
!     * Moyenne et variance de la masse de charbon reactif
!     * Moyenne et variance de la masse de coke
!     * Moyenne et variance du diametre du coeur retrecissant

          ilmp1   = ilmp   + izcl*nvlsta
          do ilayer=1,nlayer
            ilhp1(ilayer)   = ilhp(ilayer)  + izcl*nvlsta
            ilmch1(ilayer)  = ilmch(ilayer) + izcl*nvlsta
            ilmck1(ilayer)  = ilmck(ilayer) + izcl*nvlsta
          enddo
          ilmwat1 = ilmwat + izcl*nvlsta
          ildck1  = ildck  + izcl*nvlsta

          statis(iel,ilmp1)  = statis(iel,ilmp1)                   &
                        + pepa(jrpoi,npt) * eptp(jmp,npt)
          do ilayer=1,nlayer
            statis(iel,ilhp1(ilayer))  = statis(iel,ilhp1(ilayer)) &
                  + pepa(jrpoi,npt) * eptp(jhp(ilayer),npt)
            statis(iel,ilmch1(ilayer)) = statis(iel,ilmch1(ilayer))&
                  + pepa(jrpoi,npt) * eptp(jmch(ilayer),npt)
            statis(iel,ilmck1(ilayer)) = statis(iel,ilmck1(ilayer))&
                  + pepa(jrpoi,npt) * eptp(jmck(ilayer),npt)
          enddo
          statis(iel,ilmwat1) = statis(iel,ilmwat1)               &
                        + pepa(jrpoi,npt) * eptp(jmwat,npt)
          statis(iel,ildck1) = statis(iel,ildck1)                 &
                        + pepa(jrpoi,npt) * pepa(jrdck,npt)

          ilmp1   = ilmp   + izcl*(nvlsta-1)
          do ilayer=1,nlayer
            ilhp1(ilayer)  = ilhp(ilayer)  + izcl*(nvlsta-1)
            ilmch1(ilayer) = ilmch(ilayer) + izcl*(nvlsta-1)
            ilmck1(ilayer) = ilmck(ilayer) + izcl*(nvlsta-1)

          enddo
          ilmwat1 = ilmwat + izcl*(nvlsta-1)
          ildck1  = ildck  + izcl*(nvlsta-1)

          stativ(iel,ilmp1)  = stativ(iel,ilmp1)                  &
          + pepa(jrpoi,npt)*eptp(jmp,npt)*eptp(jmp,npt)
          do ilayer=1,nlayer
            stativ(iel,ilhp1(ilayer))  = stativ(iel,ilhp1(ilayer)) &
            + pepa(jrpoi,npt)*eptp(jhp(ilayer),npt) *eptp(jhp(ilayer),npt)
            stativ(iel,ilmch1(ilayer)) = stativ(iel,ilmch1(ilayer))&
            + pepa(jrpoi,npt)*eptp(jmch(ilayer),npt)*eptp(jmch(ilayer),npt)
            stativ(iel,ilmck1(ilayer)) = stativ(iel,ilmck1(ilayer))&
            + pepa(jrpoi,npt)*eptp(jmck(ilayer),npt)*eptp(jmck(ilayer),npt)
          enddo
          stativ(iel,ilmwat1) = stativ(iel,ilmwat1)               &
          + pepa(jrpoi,npt)*eptp(jmwat,npt)*eptp(jmwat,npt)
          stativ(iel,ildck1) = stativ(iel,ildck1)                 &
          + pepa(jrpoi,npt)*pepa(jrdck,npt)*pepa(jrdck,npt)

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
