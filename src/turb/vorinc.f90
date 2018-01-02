!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!> \file vorinc.f90
!> Module for vortex method for LES boundary conditions

module vorinc

  !=============================================================================

  implicit none

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

  integer, allocatable, dimension(:) :: irepvo
  integer, allocatable, dimension(:,:) :: ifacgl, ivorce

  double precision, allocatable, dimension(:,:,:) :: xyzv, yzcel
  double precision, allocatable, dimension(:,:,:) :: yzvor , yzvora
  double precision, allocatable, dimension(:,:) :: uvort, vvort, wvort
  double precision, allocatable, dimension(:,:) :: visv
  double precision, allocatable, dimension(:,:) :: signv, sigma
  double precision, allocatable, dimension(:,:,:) :: gamma
  double precision, allocatable, dimension(:,:) :: temps, tpslim

  ! irepv    : tableau associant aux faces de bord le numero d'une entree
  ! ifagl    : tableau de connectivite
  ! ivrce    : tableau reperant la cellule la plus voisine
  !             de chaque vortex
  ! xyzv     : tableaux contenant les coordonnees de
  !             toutes les faces d'entree
  ! visv     : tableau contenant la viscosite sur
  !             toutes les faces d'entree
  ! yzcel    : tableau contenant les coordonnees des
  !             faces d'entree dans le repere local
  ! uvort,...: tableaux contenant les composantes de vitesse
  ! yzvor    : du tableau contenant la position des vortex dans le repere local
  ! yzvoa    : debut du tableau contenant la position des vortex
  !             dans le repere local au pas de temps precedent
  ! signv    : debut du tableau contenant le sens de rotation des
  !             vortex
  ! xsigm    : debut du tableau contenant la taille des vortex
  ! xgamm    : debut du tableau contenant l'intensite des vortex
  ! xtmp     : debut du tableau contenant le temps cumule
  ! xtmpl    : debut du tableau contenant le temps de vie des vortex

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

  character(len=50), save :: ficvor(nentmx)

  ! ficvor : nom du fichier de donnee
  !=============================================================================

contains

  !=============================================================================

  subroutine init_vortex

    allocate(ivorce(nvomax,nnent))
    allocate(yzcel(icvmax,2,nnent))
    allocate(visv(icvmax,nnent))
    allocate(xyzv(icvmax,3,nnent))
    allocate(uvort(icvmax,nnent))
    allocate(vvort(icvmax,nnent))
    allocate(wvort(icvmax,nnent))
    allocate(yzvor(nvomax,2,nnent))
    allocate(yzvora(nvomax,2,nnent))
    allocate(signv(nvomax,nnent))
    allocate(sigma(nvomax,nnent))
    allocate(gamma(nvomax,2,nnent))
    allocate(temps(nvomax,nnent))
    allocate(tpslim(nvomax,nnent))

  end subroutine init_vortex

  !=============================================================================

  subroutine finalize_vortex

    deallocate(ivorce)
    deallocate(yzcel)
    deallocate(visv)
    deallocate(xyzv)
    deallocate(uvort)
    deallocate(vvort)
    deallocate(wvort)
    deallocate(yzvor)
    deallocate(yzvora)
    deallocate(signv)
    deallocate(sigma)
    deallocate(gamma)
    deallocate(temps)
    deallocate(tpslim)

  end subroutine finalize_vortex

end module vorinc

