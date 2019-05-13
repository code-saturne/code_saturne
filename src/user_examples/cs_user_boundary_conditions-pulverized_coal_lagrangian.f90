!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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
! Function:
! ---------
!> \file  cs_user_boundary_conditions-pulverized_coal_lagrangian.f90
!> \brief Example of cs_f_user_boundary_conditions subroutine.f90 for
!>        lagrangian pulverized coal
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[out]    icodcl        boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!> \param[in]     itrifb        indirection for boundary faces ordering
!> \param[in,out] itypfb        boundary face types
!> \param[out]    izfppp        boundary face zone number
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradt \, \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!_______________________________________________________________________________

subroutine cs_f_user_boundary_conditions &
 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     ,                                                       &
   rcodcl )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ihmpre
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use atincl
use ctincl
use cs_fuel_incl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet)
double precision rcodcl(nfabor,nvar,3)

! Local variables

!< [loc_var_dec]
integer          ifac, ii
integer          izone

integer          ilelt, nlelt

double precision uref2, d2s3
double precision xkent, xeent

integer, allocatable, dimension(:) :: lstelt
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

d2s3 = 2.d0/3.d0

!===============================================================================
! Assign boundary conditions to boundary faces here

! For each subset:
! - use selection criteria to filter boundary faces of a given subset
! - loop on faces from a subset
!   - set the boundary condition for each face
!===============================================================================

! ---- Face de type entree correspondant a une entree d'air
!        Par exemple : Air primaire , secondaire ou Air tertiaire

!< [example_1]
call getfbr('12', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!   Type de condition aux limites pour les variables standard
  itypfb(ifac) = ientre

!   Numero de zone (on les numerote de 1 a n)
  izone = 1

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

! ------ Pour ces faces d'entree , on est a debit impose

  ientat(izone) = 1
  iqimp(izone)  = 1
!      - Debit en kg/s pour l'air
  qimpat(izone) = 1.46d-03
!      - Temperature en K pour l'air
  timpat(izone) = 400.d0 + tkelvi


! ------ On impose en couleur 12 une entree a debit impose
!        L'utilisateur donne donc ici uniquement
!          la direction du vecteur vitesse

  rcodcl(ifac,iu,1) = 0.d0
  rcodcl(ifac,iv,1) = 0.d0
  rcodcl(ifac,iw,1) = 5.d0

! ------ Traitement de la turbulence

!        La turbulence est calculee par defaut si ICALKE different de 0
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference adaptes a l'entree courante si ICALKE = 1
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference et de l'intensite turvulente
!            adaptes a l'entree courante si ICALKE = 2

!      Choix pour le calcul automatique ICALKE = 1 ou 2
  icalke(izone) = 1
!      Saisie des donnees
  dh(izone)     = 0.1d0
  xintur(izone) = 0.1d0



! Exemple de cas ou ICALKE(IZONE) = 0 : DEBUT
!    Eliminer ces lignes pour la clarte si on a fait le choix ICALKE(IZONE) = 1

  if(icalke(izone).eq.0) then

!         Calcul de k et epsilon en entree (XKENT et XEENT) a partir
!           l'intensite turbulente et de lois standards en conduite
!           circulaire (leur initialisation est inutile mais plus
!           propre)
    uref2 = rcodcl(ifac,iu,1)**2                           &
           +rcodcl(ifac,iv,1)**2                           &
           +rcodcl(ifac,iw,1)**2
    uref2 = max(uref2,1.d-12)
    xkent  = epzero
    xeent  = epzero

!     (ITYTUR est un indicateur qui vaut ITURB/10)
    if    (itytur.eq.2) then

      rcodcl(ifac,ik,1)  = xkent
      rcodcl(ifac,iep,1) = xeent

    elseif(itytur.eq.3) then

      rcodcl(ifac,ir11,1) = d2s3*xkent
      rcodcl(ifac,ir22,1) = d2s3*xkent
      rcodcl(ifac,ir33,1) = d2s3*xkent
      rcodcl(ifac,ir12,1) = 0.d0
      rcodcl(ifac,ir13,1) = 0.d0
      rcodcl(ifac,ir23,1) = 0.d0
      rcodcl(ifac,iep,1)  = xeent

    elseif (iturb.eq.50) then

      rcodcl(ifac,ik,1)   = xkent
      rcodcl(ifac,iep,1)  = xeent
      rcodcl(ifac,iphi,1) = d2s3
      rcodcl(ifac,ifb,1)  = 0.d0

    elseif (iturb.eq.60) then

      rcodcl(ifac,ik,1)   = xkent
      rcodcl(ifac,iomg,1) = xeent/cmu/xkent

    elseif (iturb.eq.70) then

      rcodcl(ifac,inusa,1) = cmu*xkent**2/xeent

    endif

  endif
! Exemple de cas ou ICALKE(IZONE) = 0 : FIN

! ------ Traitement des scalaires physiques particulieres
!        Ils sont traites automatiquement


! ------ Traitement des scalaires utilisateurs

  if ( (nscal-nscapp).gt.0 ) then
    do ii = 1, (nscal-nscapp)
      rcodcl(ifac,isca(ii),1) = 1.d0
    enddo
  endif

enddo
!< [example_1]

! --- On impose en couleur 15 une paroi laterale

!< [example_2]
call getfbr('15', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          PAROI : DEBIT NUL (FLUX NUL POUR LA PRESSION)
!                  FROTTEMENT POUR LES VITESSES (+GRANDEURS TURB)
!                  FLUX NUL SUR LES SCALAIRES

!   Type de condition aux limites pour les variables standard
  itypfb(ifac)   = iparoi


!   Numero de zone (on les numerote de 1 a n)
  izone = 2

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

enddo
!< [example_2]

! --- On impose en couleur 19 une sortie

!< [example_3]
call getfbr('19', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          SORTIE : FLUX NUL VITESSE ET TEMPERATURE, PRESSION IMPOSEE
!            Noter que la pression sera recalee a P0
!                sur la premiere face de sortie libre (ISOLIB)

!   Type de condition aux limites pour les variables standard
  itypfb(ifac)   = isolib

!   Numero de zone (on les numerote de 1 a n)
  izone = 3

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

enddo
!< [example_3]


! --- On impose en couleur 14 et 4 une symetrie

!< [example_4]
call getfbr('14 or 4', nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          SYMETRIES

!   Type de condition aux limites pour les variables standard
  itypfb(ifac)   = isymet

!   Numero de zone (on les numerote de 1 a n)
  izone = 4

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

enddo
!< [example_4]

!--------
! Formats
!--------

!----
! End
!----

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine cs_f_user_boundary_conditions
