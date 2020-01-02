!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

subroutine lwcurl &
!================

 ( ampen1 , valmoy , valvar , valmin , valmax ,                   &
   exit01 , exit02  , ampl01 , ampl02  )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DES PARAMETRES DE LA PDF
! PDF LIBBY - WILLIAMS 2 POINTS EN UTILISANT LE
! MOMENT D'ORDRE 3 DEDUIT DES RECURRENCES SUR
! LES FONCTIONS BETA

! COMMENTAIRES : HYPOTHESE DE CURL MODIFIEE
! ------------
!                a partir de la valeur moyenne d'une variable,
!                des extremas et de la variance de cette variable
!                on en deduit 2 etat autour de l'etat moyen
!                et une amplitude pour chaque etat

! LE RESULTAT EST :
! ---------------

!    CALCUL DES PARAMETRES ASSOCIES AUX FONCTIONS DIRAC

!      Les Diracs sont en position [F(.,1),Y(.,1)] et [F(.,2),Y(.,2)]
!      Leurs amplitudes respectives sont D(.,1) et D(.,2)


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ampen1           ! r  ! <-- !  amplitude totale des pics                     !
! valmoy           ! r  ! <-- !  valeur moyenne de la variable                 !
! valvar           ! r  ! <-- !  variance de la variable                       !
! valmin           ! r  ! <-- !  min de la variable                            !
! valmax           ! r  ! <-- !  max de la variable                            !
! exit01           ! r  ! --> !  etat 1                                        !
! exit02           ! r  ! --> !  etat 2                                        !
! ampl01           ! r  ! --> !  amplitude 1                                   !
! ampl02           ! r  ! --> !  amplitude 2                                   !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!==============================================================================
! Module files
!==============================================================================

use paramx
use pointe
use entsor
use cstnum
use cstphy
use ppppar
use ppthch
use coincl
use cpincl
use ppincl

!===============================================================================

implicit none

! Arguments

double precision valmoy, valvar
double precision valmin, valmax
double precision exit01, exit02
double precision ampl01, ampl02
double precision ampen1

! Local variables

double precision tvv, c, d
double precision moyadm, varadm, tvvadm
double precision epsi
!===============================================================================
! 0.  CALCULS PRELIMINAIRES
!===============================================================================

! ---> Test sur l'amplitude totale des pics
!       si elle est trop faible on positionne
!        les deux etat sur l'etat moyen

epsi = 1.d-6

if (ampen1.gt.epsi) then

! ---> Test sur la variance,
!       si elle est trop faible on positionne
!        les deux etat sur l'etat moyen

  if ((valvar.gt.epsi)) then

! ---> on travaille en variable adimentionnelle pour ce calcul

    moyadm = (valmoy-valmin)/(valmax-valmin)
    varadm = valvar/((valmax-valmin)**2)

! ---> calcul du moment d'ordre 3 adim

    tvvadm = 2.d0*varadm**2                                       &
         *((1.d0-2.d0*moyadm)/((moyadm*(1-moyadm))+varadm))

! ---> calcul du moment d'ordre 3 non adim

    tvv = ((valmax-valmin)**3)*tvvadm

! ---> calcul du termedu polynome du moment 3

    c=(4.d0+tvv**2/valvar**3)

! ---> determination du signe de la racine

    if ((1.d0-moyadm).gt.moyadm) then
      d = 0.5d0+sqrt((c-4.d0)/(4.d0*c))
    else
      d = 0.5d0-sqrt((c-4.d0)/(4.d0*c))
    endif

! ---> calcul des amplitudes des 2 pics

    ampl01 = ampen1 * d
    ampl02 = ampen1 - ampl01

! -----> Calcul des positions des deux pics

    exit01 = valmoy - sqrt((1.d0-d)/(d)*valvar)
    exit02 = valmoy + sqrt((d)/(1.d0-d)*valvar)

    exit01 = max(valmin,min(exit01,valmax))
    exit02 = max(valmin,min(exit02,valmax))

  else
! variance faible

    ampl01 = ampen1/2.d0
    ampl02 = ampl01

    exit01 = valmoy
    exit02 = valmoy

  endif

else
! amplitude totale faible

  ampl01 = ampen1/2.d0
  ampl02 = ampl01

  exit01 = valmoy
  exit02 = valmoy

endif

end subroutine
