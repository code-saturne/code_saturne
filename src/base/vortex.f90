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

subroutine vortex &
!================

 ( ivorce , visco  , xyz    ,                                     &
   yzcel  , xu     , xv     , xw     ,                            &
   yzvor  , yzvora , signv  ,                                     &
   sigma  , gamma  , temps  , tpslim )

!===============================================================================
!  FONCTION  :
!  ----------

! GESTION DES ENTREES L.E.S. PAR LA METHODE DES VORTEX

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ivorce           ! te ! <-- ! numero du vortex le plus proche d'un           !
!     (nvomax)     !    !     ! vortex donne                                   !
! visco            ! tr ! <-- ! viscosite cinematique sur les faces            !
!(icvmax,nnent)    !    !     ! d'entree                                       !
! xyz(icvmax,3)    !    ! <-- ! coordonnees des cellules d'entree              !
!                  !    !     ! dans le referentiel global                     !
! yzcel            ! tr ! <-- ! coordonnees des faces d'entree dans            !
!   (nelvmx ,2)    !    !     ! le referentiel local                           !
! xu(nelvmx)       ! tr ! <-- ! composante de vitesse principale               !
! xv(nelvmx)       ! tr ! <-- ! composantes de vitesse transverses             !
! xw(nelvmx)       ! tr ! <-- !                                                !
! yzvor            ! tr ! <-- ! coordonnees du centre des vortex               !
!   (nvomax,2)     !    !     !                                                !
! yzvora           ! tr ! <-- ! anciennes coordonnees du centre                !
!   (nvomax,2)     !    !     ! des vortex                                     !
! signv(nvomax)    ! tr ! <-- ! sens de rotation des vortex                    !
! sigma            ! tr ! <-- ! taille des vortex                              !
!(nvomax,nnent)    !    !     !                                                !
! gamma            ! tr ! <-- ! intensite des vortex                           !
!(nvomax,2,nnen    !    !     ! (dans les deux directions du plan)             !
! temps            ! tr ! <-- ! temps ecoule depuis la creation                !
!     (nvomax)     !    !     ! du vortex                                      !
! tpslim           ! tr ! <-- ! duree de vie du vortex                         !
!(nvomax,nnent)    !    !     !                                                !
!__________________.____._____.________________________________________________.

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "entsor.h"
include "optcal.h"
include "vortex.h"

!===============================================================================

! Arguments

integer          ivorce(nvomax,nnent)

double precision yzcel(icvmax,2,nnent)    , visco(icvmax,nnent)
double precision xyz(icvmax,3,nnent)      , xu(icvmax,nnent)
double precision xv(icvmax,nnent)         , xw(icvmax,nnent)
double precision yzvor(nvomax,2,nnent)    , yzvora(nvomax,2,nnent)
double precision signv(nvomax,nnent)
double precision sigma(nvomax,nnent)      , gamma(nvomax,2,nnent)
double precision temps(nvomax,nnent)      , tpslim(nvomax,nnent)

!     VARIABLES LOCALES

character        ficsui*32

integer          ii, ient
integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
! 1. INITIALISATION
!===============================================================================


! L'EQUATION DE LANGEVIN RESTE A TRAVAILLER POUR UN CAS 3D QUELCONQUE
! OU A MODIFIER    . VU LE PEU D'IMPORTANCE QU'ELLE A SUR CE QUI SE PASSE
! EN AVAL DE L'ENTREE (L'IMPORTANT ETANT D'IMPOSER V' ET W'), ON NE VA
! PAS PLUS LOIN (ON ANNULE CES CONTRIBUTION POUR LE MOMENT POUR LES
! CAS 3 ET 4)

ipass = ipass + 1

do ient = 1, nnent

  if (ipass.eq.1)then

    call vorini                                                   &
    !==========
 ( icvor(ient)     , nvort(ient)     ,                            &
   ient   , ivorce(1,ient)  ,                                     &
   xyz(1,1,ient)   , yzcel(1,1,ient) ,                            &
   xu(1,ient)      , xv(1,ient)      , xw(1,ient)      ,          &
   yzvor(1,1,ient) , signv(1,ient)   , temps(1,ient)   ,          &
   tpslim(1,ient)  )

  endif

!===============================================================================
! 2. DEPLACEMENT DU VORTEX
!===============================================================================

  call vordep                                                     &
  !==========
 ( icvor(ient)     , nvort(ient)     , ient   , dtref  ,          &
   ivorce(1,ient)  , yzcel(1,1,ient) ,                            &
   xu(1,ient)      , xv(1,ient)      , xw(1,ient)      ,          &
   yzvor(1,1,ient) , yzvora(1,1,ient), signv(1,ient)   ,          &
   temps(1,ient)   , tpslim(1,ient)  )

!===============================================================================
! 3. CALCUL DE LA VITESSE
!===============================================================================

  call vorvit                                                     &
  !==========
 ( icvor(ient)     , nvort(ient)     , ient   ,                   &
   ivorce(1,ient)  , visco(1,ient)   ,                            &
   yzcel(1,1,ient) , xu(1,ient)      , xv(1,ient)       ,         &
   xw(1,ient)      , yzvor(1,1,ient) , signv(1,ient)    ,         &
   sigma(1,ient)   , gamma(1,1,ient) , temps(1,ient)    )

!===============================================================================
! 4. CALCUL DES FLUCTUATIONS DANS LE SENS DE L'ECOULEMENT
!===============================================================================

  call vorlgv                                                     &
  !==========
 ( icvor(ient)     , ient   , dtref  ,                            &
   yzcel(1,1,ient) , xu(1,ient)      ,                            &
   xv(1,ient)      , xw(1,ient)      )

enddo

!===============================================================================
! 5. ECRITURE DU FICHIER SUITE
!===============================================================================

! on ecrit a tous les pas de temps pour eviter
! les mauvaises surprises en cas de fin prematuree.
! Il ne faut pas mettre cette partie dans la boucle sur IENT
! car on accede deja a l'unite IMPMVO(=IMPDVO) dans VORINI.
! Seul le premier processeur ecrit (test avant l'appel à VORTEX)

ficsui = 'checkpoint/vortex'
open(unit=impvvo,file=ficsui)
rewind(impvvo)
do ient = 1, nnent
  write(impvvo,100) ient
  write(impvvo,100) nvort(ient)
  do ii = 1, nvort(ient)
    write(impvvo,200) yzvor(ii,1,ient),yzvor(ii,2,ient),          &
         temps(ii,ient), tpslim(ii,ient), signv(ii,ient)
  enddo
enddo
close(impvvo)

!===============================================================================
! 6. FIN
!===============================================================================

! FORMATS

 100  format(i10)
 200  format(5e13.5)

return

end subroutine
