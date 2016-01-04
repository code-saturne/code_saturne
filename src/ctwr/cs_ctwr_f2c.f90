!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

subroutine defct &
!================

 ( idimze, nomze,   imzech,  ntypze, neleze,          &
   deltat, teaueze, qeaueze, xap, xnp, surface, dgout )

!===============================================================================
! FONCTION :
! ----------

!     DEFINITION DE COUPLAGE(S) AVEC LE CODE SYRTHES

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idimze           ! e  ! <-- ! dimension de la zone d'echange                 !
! nomze            ! a  ! <-- ! caracterisation de zone d'echange              !
! imzech           ! e  ! <-- ! modele de zone d'echange (impose par usppmo)   !
! ntypze           ! e  ! <-- ! types de zone d'echanges                       !
! neleze           ! e  ! <-- ! nombre d'element pour extrusion de la zone     !
! deltat           ! r  ! <-- ! ecart de temperature impose                    !
! teaueze          ! r  ! <-- ! Temperature eau de la zone d'echange           !
! qeaueze          ! r  ! <-- ! debit d'eau traversant la zone d'echange       !
! xap              ! r  ! <-- ! coefficent                                     !
! xnp              ! r  ! <-- ! coefficient                                    !
! surface          ! r  ! <-- ! surface de la zone d'echange                   !
! dgout            ! r  ! <-- ! diametre de goutte                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

character*(*) nomze

integer          imzech,ntypze,idimze,neleze

double precision teaueze,qeaueze,deltat
double precision xap,xnp,surface
double precision dgout

! Variables locales

integer       lnomze

!===============================================================================

lnomze = len(nomze)


call defct1( idimze, nomze, lnomze, imzech, ntypze, neleze,      &
!==========
              deltat, teaueze, qeaueze, xap, xnp, surface, dgout )


return

end subroutine
