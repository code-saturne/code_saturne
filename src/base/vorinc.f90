!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

! Module for vortex method for LES boundary conditions

module vorinc

  !=============================================================================

  ! --------------
  ! Parametres max
  ! --------------
  integer    nentmx, ndatmx
  parameter (nentmx = 10)
  parameter (ndatmx = 10000)

  ! nentmx    : nombre d'entree max
  ! ndatmx    : nombre de points max pour le fichier des donnees

  ! ----------
  ! Dimensions
  ! ----------
  integer, save :: icvor(nentmx), icvor2(nentmx), icvmax , nvomax

  ! icvor  : nombre de faces (global) utilisant des vortex
  !          pour chaque entree
  ! icvor2 : compteur du nombre local de faces utilisant des vortex
  ! icvmax : nombre max de faces utilisant des vortex (sur toutes entrees
  !          confondues)
  ! nvomax : nombre max de vortex utilise (toutes entrees confondues)

  ! ---------
  ! pointeurs
  ! ---------

  integer, save ::  iirepv , iifagl , iivrce ,                    &
                    ixyzv  , ivisv  ,                             &
                    iyzcel , iuvort , ivvort , iwvort ,           &
                    iyzvor , iyzvoa , isignv , ixsigm ,           &
                    ixgamm , ixtmp  , ixtmpl ,                    &
                    iw1x   , iw1y   , iw1z   , iw1v   ,           &
                    iw2x   , iw2y   , iw2z   , iw2v

  ! iirepv    : debut du tableau associant aux faces de bord
  !             le numero d'une entree
  ! iifagl    : debut du tableau de connectivite
  ! iivrce    : debut du tableau reperant la cellule la plus voisine
  !             de chaque vortex
  ! ixyzv     : debut du tableaux contenant les coordonnees de
  !             toutes les faces d'entree
  ! ivisv     : debut du tableau contenant la viscosite sur
  !             toutes les faces d'entree
  ! iyzcel    : debut du tableau contenant les coordonnees des
  !             faces d'entree dans le repere local
  ! iuvort,...: debuts des tableaux contenant les composantes de vitesse
  ! iyzvor    : debut du tableau contenant la position des vortex
  !             dans le repere local
  ! iyzvoa    : debut du tableau contenant la position des vortex
  !             dans le repere local au pas de temps precedent
  ! isignv    : debut du tableau contenant le sens de rotation des
  !             vortex
  ! ixsigm    : debut du tableau contenant la taille des vortex
  ! ixgamm    : debut du tableau contenant l'intensite des vortex
  ! ixtmp     : debut du tableau contenant le temps cumule
  ! ixtmpl    : debut du tableau contenant le temps de vie des vortex
  ! iw1x,...  : debut des tableaux de travails servant a communiquer
  !             les donnees aux entrees a tous les processeurs
  !             (plus utilise apres vorpre)

  ! -----------------
  ! Options de calcul
  ! -----------------

  integer, save :: nnent, nvort(nentmx), initvo(nentmx),          &
                   icas(nentmx), itlivo(nentmx),                  &
                   isgmvo(nentmx), idepvo(nentmx),                &
                   iclvor(4,nentmx), ndat(nentmx)

  ! nnent  : nombre d entrees utilisees
  ! nvort  : nombre de vortex
  ! initvo : indicateur de reinitialisation
  ! icas   : type de geometrie pour l'entree
  ! itlivo : type de modele pour la duree de vie
  ! isgmvo : type de modele pour la taille des vortex
  ! idepvo : type de modele pour la marche en temps
  ! iclvor : type de condition aux limites
  ! ndat   : nombre de lignes du fichier de donnees

  ! -------
  ! Donnees
  ! -------

  double precision, save :: tlimvo(nentmx), xsgmvo(nentmx), ud(nentmx),      &
                            xdat(ndatmx,nentmx),                             &
                            ydat(ndatmx,nentmx), zdat(ndatmx,nentmx),        &
                            udat(ndatmx,nentmx),                             &
                            vdat(ndatmx,nentmx), wdat(ndatmx,nentmx),        &
                            dudat(ndatmx,nentmx),                            &
                            kdat(ndatmx,nentmx), epsdat(ndatmx,nentmx),      &
                            udebit(nentmx), kdebit(nentmx), edebit(nentmx),  &
                            dir1(3,nentmx), dir2(3,nentmx), dir3(3,nentmx),  &
                            cen(3,nentmx) , surf(3,nentmx),                  &
                            ymax(nentmx)  , ymin(nentmx),                    &
                            zmax(nentmx)  , zmin(nentmx),                    &
                            xsurfv(nentmx), llz(nentmx),                     &
                            lly(nentmx)   , lld(nentmx)

  ! tlimvo      : temps de vie max des vortex impose par l'utilisateur
  ! xsgmvo      : diametre des vortex impose par l'utilisateur
  ! ud          : vitesse de deplacement (max) imposee par l'utilisateur
  ! xdat, ...   : coordonnees des points ou sont connues les donnees
  ! udat        : vitesse moyenne principale (fichier de donnees)
  ! vdat,wdat   : vitesse moyenne transverse (fichier de donnees)
  ! dudat       : derive normale de la vitesse principale (fichier d'entree)
  ! kdat        : ec moyenne (fichier d'entree)
  ! epsdat      : dissipation (fichier d'entree)
  ! udebit      : vitesse moyenne imposee par l'utilisateur en entree
  ! kdebit      : ec imposee par l'utilisateur en entree
  ! edebit      : dissipation imposee par l'utilisateur en entree
  ! dir1,...    : vecteurs definissant le repere local dans le plan d'entree
  ! cen         : coordonnees du centre de l'entree
  ! surf        : vecteur surface du plan d'entree (supposee plane)
  ! xmax,...    : dimensions max de l'entree dans le repere local
  ! llz,lly,lld : dimensions de l'entree dans le calcul

  character*50, save :: ficvor(nentmx)

  ! ficvor : nom du fichier de donnee
  !=============================================================================

end module vorinc

