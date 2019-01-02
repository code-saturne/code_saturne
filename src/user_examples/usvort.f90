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

!< [arg]

! Arguments

integer          nvar   , nscal
integer          iappel

double precision dt(ncelet)

!< [arg]

!< [loc_var_dec]

! Local variables

integer          ifac, ient
integer          ilelt, nlelt

integer, allocatable, dimension(:) :: lstelt

!< [loc_var_dec]

!===============================================================================

!< [allocate]

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))

!< [allocate]

!< [glob_param]

!===============================================================================
! 1. PARAMETRES GLOBAUX
!===============================================================================

! --- Nombre d'entrees avec la methode des vortex

nnent = 2

! --- Nombre de vortex a mettre dans chaque entree

!   NVORT min ~ Surface d'entree/(pi*SIGMA**2)

nvort(1) = 500
nvort(2) = 500

!< [glob_param]

!< [inlet]

if (iappel.eq.1) then

!===============================================================================
! 2. DEFINITION DES ZONES D'ENTREE (AU PREMIER PASSAGE)
!===============================================================================

  do ifac = 1, nfabor
    irepvo(ifac) = 0
  enddo

! ------------------
!   ENTREE 1
! ------------------
  CALL GETFBR('3',NLELT,LSTELT)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)

    ient = 1
    irepvo(ifac) = ient

  enddo

! ------------------
!   ENTREE 2
! ------------------
  CALL GETFBR('1',NLELT,LSTELT)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)

    ient = 2
    irepvo(ifac) = ient

  enddo

!< [inlet]

!< [bc]

elseif (iappel.eq.2) then

!===============================================================================
! 3. PARAMETRES GEOMETRIQUES ET CONDITIONS LIMITES
!===============================================================================

! --- Cas traite

! ICAS = 1...Conduite rectangulaire
!        2...Conduite circulaire
!        3...Geometrie quelconque sans traitement specifique des conditions aux limites
!        4...Geometrie quelconque sans traitement specifique des conditions aux limites
!            ni fichier de donnees (la vitesse moyenne, le niveau de k et de epsilon
!            sont fournis par l'utilisateur)

  ient = 1
  icas(ient) = 1

  ient = 2
  icas(ient) = 2


! --- Repere definissant le plan d'entree

!     Si ICAS = 4, le code se charge de ces donnees
!     Sinon il faut preciser les vecteurs DIR1 et DIR2 definissant
!     un repere directe tel que DIR3 soit un vecteur entrant normal
!     a la face d'entree.

  ient = 1
  if(icas(ient).eq.1.or.icas(ient).eq.2.or.icas(ient).eq.3) then
    dir1(1,ient) = 1.d0
    dir1(2,ient) = 0.d0
    dir1(3,ient) = 0.d0

    dir2(1,ient) = 0.d0
    dir2(2,ient) = 1.d0
    dir2(3,ient) = 0.d0
  endif

  ient = 2
  if(icas(ient).eq.1.or.icas(ient).eq.2.or.icas(ient).eq.3) then
    dir1(1,ient) = 0.d0
    dir1(2,ient) = 1.d0
    dir1(3,ient) = 0.d0

    dir2(1,ient) = 0.d0
    dir2(2,ient) = 0.d0
    dir2(3,ient) = 1.d0
  endif

! --- Centre du repere local dans le plan d'entree

!     Si ICAS = 1 ou ICAS = 2, le centre du repere doit correspondre
!                 au centre de gravite de la zone d'entree (rectangle ou cercle)

  ient = 1

  cen(1,ient) = 0.d0
  cen(2,ient) = 0.d0
  cen(3,ient) = -6.05d-1

  ient = 2

  cen(1,ient) = -3.664d-1
  cen(2,ient) = 0.d0
  cen(3,ient) = 0.d0

! --- Condition aux limites

! -> Si ICAS = 1...Il faut specifier le type de condition aux limite ICLVOR
!               dans les directions DIR1, DIR2, - DIR1, -DIR2

!               Ces conditions peuvent etre de 3 types :

! ICLVOR = 1...Condition de paroi
!          2...Condition de symetrie
!          3...Condition de periodicite

!                    y = LLY/2
!                    (ICLVOR 1)
!           +-----------------------+
!           |           ^ DIR1      |
!           |           |           |
!           |           |           |
! z=- LLZ/2 |           +----> DIR2 | z = LLZ/2
! (ICLVOR 4)|                       | (ICLVOR 2)
!           |                       |
!           |                       |
!           +-----------------------+
!                    y = -LLY/2
!                    (ICLVOR 3)


! -> Si ICAS = 2, les conditions sont necessairement de type paroi
! -> Si ICAS = 3 ou 4, pas de traitement particulier

  ient = 1

  if(icas(ient).eq.1) then
    iclvor(1,ient) = 1
    iclvor(2,ient) = 2
    iclvor(3,ient) = 1
    iclvor(4,ient) = 2
  endif

! LLY et LLZ sont les dimensions de l'entree dans les directions DIR1 et DIR2
! LDD est le diametre de la conduite


  ient = 1
  lly(ient) = 0.2d0
  llz(ient) = 0.1d0

  ient = 2
  lld(2) = 0.154d0

!< [bc]

!< [param]

!===============================================================================
! 5. PARAMETRES PHYSIQUES ET MARCHE EN TEMPS
!===============================================================================

! --- " Temps de vie " limite du vortex

! ITLIVO = 1...Les vortex sont retire au bout du temps TLIMVO
!                donne par l'utilisateur
!                ( par exemple TLIMVO = 10*DTREF)

!          2...Chaque vortex a un temps d'exitence limite valant
!                5.Cmu.k^(3/2).U/epsilon
!               ( ou U est la vitesse principale suivant DIR3)

  ient = 1
  itlivo(ient) = 1

  if(itlivo(ient).eq.1) then
    tlimvo(ient) = 10.d0*dtref
  endif

  ient = 2
  itlivo(ient) = 2


! --- " Diametre " des vortex

! ISGMVO = 1...diametre constant XSGMVO donne par l'utilisateur
!          2...basee sur la formule sigma = Cmu^(3/4).k^(3/2)/epsilon
!          3...basee sur la formule sigma = max(Lt, Lk) avec
!                 Lt = (5 nu.k/epsilon)^(1/2)
!             et  Lk = 200.(nu^3/epsilon)^(1/4)

  ient = 1
  isgmvo(ient) = 1

  if(isgmvo(ient).eq.1) then
    xsgmvo(ient) = 0.01d0
  endif

  ient = 2
  isgmvo(ient) = 2


! --- Mode de deplacement des vortex

! IDEPVO = 1...Deplacement en r*UD (r aleatoire dans [0,1])
!              UD a fournir par l'utilisateur
!          2...Convection par les vortex

  ient = 1
  idepvo(ient) = 2

  ient = 2
  idepvo(ient) = 1

  if(idepvo(ient).eq.1) then
    ud(ient) = 0.7d0
  endif

!< [param]

!< [input]

!===============================================================================
! 6. PARAMETRES D'ENTREE / SORTIES ET DONNEES UTILISATEUR
!===============================================================================

! --- Fichier de donnees utilisateur

! NDAT ...Nombre de lignes du fichier de donnees contenant les donnees :
!          x | y | z | U | V | W | Grad[u.DIR3].n | k | epsilon

!         dans le plan d'entree du calcul

!         Grad[u.DIR3].n est le gradient dans la direction normale
!         a la paroi, de la vitesse principale dans le plan d'entree.
!         Cette donnees n'est utilisee qu'avec ICAS=2

! FICVOR...Nom du fichier de donnees utilisateur

  ient = 1
  ndat(ient) = 2080

  ient = 2
  ndat(ient) = 2080

! Par les defaut les fichiers sont nommes "vordat" affecte de l'indice d'entree

  ient = 1
  FICVOR(IENT) = 'entree_1.dat'

  ient = 2
  FICVOR(IENT) = 'entree_2.dat'

! Pour ICAS = 4, on precise juste la valeur moyenne de U, k et de espilon
! a l'entree

  if(icas(ient).eq.4) then
    udebit(ient) = 10.d0
    kdebit(ient) = 1.d0
    edebit(ient) = 1.d0
  endif

! --- Relecture d'un fichier suite eventuel

! ISUIVO = 0...Pas de relecture (reinitialisation des vortex)
!          1...Relecture du fichier suite de methode des vortex

  isuivo = isuite


endif

!< [input]

!< [deallocate]

! Deallocate the temporary array
deallocate(lstelt)

!< [deallocate]

return
end subroutine usvort

!===============================================================================
! 7. DEFINTION DE LA FONCTION PERMETTANT D'IMPOSER LES DONNEES D'ENTREE
!===============================================================================

function phidat &
!==============

 ( nfecra , icas   , ndat   ,                                     &
   yy     , zz     , ydat   , zdat   ,                            &
   vardat , iii    )

!===============================================================================
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

!< [phidat]

implicit none

integer          nfecra, icas, ndat, iii
double precision zz, yy
double precision zdat(ndat), ydat(ndat)
double precision vardat(ndat)

integer          ii
double precision phidat, dist1

! Initialize variables to avoid compiler warnings

phidat = 0.d0

! Dans l'exemple suivant, on se contente de retourne la valeur situee
! dans le fichier de donnee a l'abscisse la plus proche du point de
! coordonnee (Y,Z) ou l'on cherche a connaitre la valeur de la
! variable numero VARDAT.


if(icas.eq.1.or.icas.eq.2.or.icas.eq.3) then

  if(iii.eq.0) then
    dist1 = 1.d20
    do ii = 1,ndat
      if(sqrt((yy-ydat(ii))**2+(zz-zdat(ii))**2).lt.dist1) then
        dist1 = sqrt((zz-zdat(ii))**2+(yy-ydat(ii))**2)
        iii   = ii
        phidat = vardat(ii)
      endif
    enddo
  elseif(iii.ne.0) then
    phidat =  vardat(iii)
  endif

elseif(icas.eq.4) then
  phidat = vardat(1)
endif

!< [phidat]

return
end function phidat