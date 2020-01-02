!-------------------------------------------------------------------------------

!VERS

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

!===============================================================================
!
!> \file usvort.f90
!>
!> \brief Unsteady inlet boundary condition for LES with the vortex method.
!>
!> See \subpage us_vort for examples.
!
!-------------------------------------------------------------------------------

subroutine usvort &
!================

 ( nvar   , nscal  ,                                              &
   iappel ,                                                       &
   dt     )

!===============================================================================
!>
!> \brief User subroutine
!>
!> METHODE DES VORTEX POUR LES CONDITIONS AUX LIMITES D'ENTREE
!> EN L.E.S. :
!> DEFINITION DES ENTREES AVEC VORTEX
!> DEFINITION DES CARACTERISTIQUES DES VORTEX
!>
!> Boundary faces identification
!>
!>
!> Boundary faces may be identified using the \ref getfbr subroutine.
!> The syntax of this subroutine is described in the
!> \ref cs_user_boundary_conditions subroutine,
!> but a more thorough description can be found in the user guide.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     iappel        indique les donnes a renvoyer
!> \param[in]     dt            time step (per cell)
!______________________________________________________________________________!

!===============================================================================
! Module files
!===============================================================================

use paramx
use optcal
use entsor
use vorinc
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          iappel

double precision dt(ncelet)

! Local variables

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine usvort

function phidat &
!==============

 ( nfecra , icas   , ndat   ,                                     &
   yy     , zz     , ydat   , zdat   ,                            &
   vardat , iii    )

!===============================================================================
!> \brief User subroutine
!>
!>
!> FONCTION PERMETTANT D'INTERPOLER LES DONNEES D'ENTREE FOURNIES
!> PAR L'UTILISATEUR AU CENTRE DES FACES D'ENTREE POUR LESQUELLES
!> EST UTILISEE LA METHODE DES VORTEX

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in] nfecra            unite
!> \param[in] icas              type de geometrie du cas
!> \param[in] ndat              nbr de lignes du fichier de donnees
!> \param[in] yy                coordoonnes dans le repere local du
!> \param[in] zz                point ou l'on cherche a connaitre la
!>                              variable vardat
!> \param[in] ydat              coordoonnes ou est connue la variable
!> \param[in] zdat              vardat dans le fichier de donnees
!> \param[in] vardat            valeur de la variable vardat
!> \param[out] iii              ligne ou a ete trouvee la donnee la
!                               plus proche du point (yy,zz)
!______________________________________________________________________________!

implicit none

integer          nfecra, icas, ndat, iii
double precision zz, yy
double precision zdat(ndat), ydat(ndat)
double precision vardat(ndat)

integer          ii
double precision phidat, dist1

! Initialize variables to avoid compiler warnings

phidat = 0.d0

return
end function phidat