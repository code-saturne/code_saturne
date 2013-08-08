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

!> \file entsor.f90
!> Module for input/output

module entsor

  !=============================================================================

  use paramx

  implicit none

  !=============================================================================

  ! --- sortie std
  !     nfecra : unite sortie std
  !     iwarni : niveau d'impression

  integer, save :: nfecra, iwarni(nvarmx)

  ! --- fichier methode des vortex

  !     impmvo --> fichier suite methode des vortex (ASCII
  !                obligatoirement, structure specifique)
  !     impvvo --> fichier suite methode des vortex (ASCII
  !                obligatoirement, structure specifique)
  !     impdvo --> fichier de donnees de la methode des vortex
  !                (nom FICDAT laisse a l'utilisateur dans usvort)

  integer, save :: impmvo, impvvo, impdvo

  character*13, save :: ficdat

  ! ---  frequence de sauvegarde

  integer, save :: ntsuit

  ! --- Sorties post traitement

  !     ichrvr : on sort la variable (1) ou non (0) ou non initialise

  integer, save :: ichrvr(nvppmx)

  ! --- Fichier thermophysique specifique physique particuliere

  !     imp    --> unite logique du fichier
  !     fpp    --> fichier utilisateur lorsqu'on utilise Janaf
  !     jon    --> Utilisation de Janaf ou non

  character*32, save :: ficfpp

  integer, save :: impfpp, indjon

  ! --- Input files for the atmospheric specific physics
  !        impmet --> logical unit of the meteo profile file
  !        ficmet --> name of the meteo profile file

  integer, save :: impmet

  character*32, save :: ficmet

  ! --- Fichiers historiques

  !     nushmx : nombre max de fichiers utilisateur pour historiques
  !     emphis : emplacement
  !     prehis : prefixe des fichiers
  !     impush : unite fichiers specifiques ushist
  !     ficush : nom   fichiers specifiques ushist
  !     impsth : fichier stock + unite d'ecriture des variables
  !              des structures mobiles

  integer    nushmx
  parameter(nushmx=16)

  character*80, save :: emphis, prehis
  character*13, save :: ficush(nushmx)

  integer, save :: impush(nushmx), impsth(2)

  ! ncaptm : nombre max de sondes (pour historiques)
  !          voir le format associe dans ecrhis

  integer    ncaptm
  parameter(ncaptm=100)

  ! tplfmt : time plot format (1: .dat, 2: .csv, 3: both)
  ! ncapt  : nombre de sondes total (limite a ncaptm)
  ! nthist : frequence de sortie: > 0 ou -1 (jamais) ou non initialise -999
  ! frhist : frequence de sortie en secondes
  ! nthsav : periode de sauvegarde (> 0 (fichiers ouverts et refermes) ou -1 )
  ! ihisvr : nb de sonde et numero par variable (-999 non initialise)
  ! ihistr : indicateur d'ecriture des historiques des structures
  !          mobiles internes (=0 ou 1)
  ! ncapt  : nombre de sondes total (limite a ncaptm)
  ! nodcap : element correspondant aux sondes
  ! ndrcap : rang processus contenant nodcap (parallelisme)
  ! xyzcap : position demandee des sondes
  ! tplflw : time plot flush wall-time interval (none if <= 0)

  integer, save :: tplfmt, ncapt, nthist, nthsav,          &
                   ihisvr(nvppmx,ncaptm+1), ihistr,        &
                   nodcap(ncaptm), ndrcap(ncaptm)

  double precision, save :: tplflw, frhist, xyzcap(3,ncaptm)

  ! --- Fichiers Lagrangiens

  !   - Fichier listing Lagrangien

  !     ficlal : nom du fichier
  !     implal : unite du fichier
  !     ntlal  : periode de sortie du listing

  character*6, save :: ficlal

  integer, save :: implal, ntlal

  !   - Fichiers Lagrangien

  integer, save :: impla1 , impla2 , impla3 , impla4 , impla5(15)

  ! --- Fichiers utilisateurs

  ! nusrmx : nombre max de fichiers utilisateur

  integer    nusrmx
  parameter(nusrmx=10)

  character*13, save :: ficusr(nusrmx)
  integer, save ::      impusr(nusrmx)

  ! --- Sorties listing

  !   Communes
  !     ipp*   : pointeurs de reperage des variables pour les sorties
  !     nomvar : nom des variables
  !     ilisvr : on suit la variable (1) ou non (0) ou non initialise
  !     itrsvr : numero de variable si ipp correspond a une variable resolue
  !              (p,u,k...)
  !              0 si ipp correspond a une variable annexe (cp, mut...)
  !              ou a rien
  !     ntlist : periode  ecriture
  !       ( -1 : dernier pas de temps : > 0 : periode)

  integer, save :: ipprtp(nvarmx), ipppro(npromx),                &
                   ippdt, ipptx, ippty, ipptz,                    &
                   ipp2ra(nvppmx)

  character*80, save :: nomvar(nvppmx)

  integer, save :: ilisvr(nvppmx),itrsvr(nvppmx)

  integer, save :: ntlist

  ! Parametres de convergence, norme du second membre, nombre iterations
  ! residu norme, derive.

  integer, save :: nbivar(nvppmx)

  double precision, save :: rnsmbr(nvppmx) ,                      &
                            resvar(nvppmx) , dervar(nvppmx)

  ! Parametres des sorties au bord =

  ! ipstdv = post-processed properties
  !   ipstfo: Efforts (1: all; 2: tangent; 4: normal)
  !   ipstyp: yplus
  !   ipsttp: Tplus
  !   ipstft: flux thermique reconstruit
  !   ipsttb: boundary temperature
  !   ipstnu: Nusselt

  integer    ipstfo  , ipstyp  , ipsttp   , ipstft  , ipsttb  , ipstnu
  parameter (ipstfo=1, ipstyp=2, ipsttp= 3, ipstft=4, ipsttb=5, ipstnu=6)

  !   si ipstdv = 1 = pas de sortie

  integer, save :: ipstdv(6)

  ! --- Temps CPU

  double precision, save :: tmarus

  !=============================================================================

end module entsor
