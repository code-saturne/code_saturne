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

! Module for Lagrangian: parameters

module lagpar

  !=============================================================================

  implicit none

  !=============================================================================

  !         Trois modules complementaires
  !                            lagran qui porte les non dimensions
  !                            lagdim qui porte les dimensions variables
  !                            lagpar qui porte les parametres

  !=============================================================================
  ! 1. Classes et particules

  !     NCHARM2 : NOMBRE MAXIMAL DE CLASSES DE CHARBON (voir cpincl.f90)

  integer         ncharm2
  parameter      (ncharm2 = 5)

  !     NCLSTM : NOMBRE MAXIMUM DE STATISTIQUES VOLUMIQUE PAR GROUPE

  integer         nclstm
  parameter      (nclstm = 100)

  !=============================================================================
  ! 2. Conditions aux limites

  !     NFLAGM : NOMBRE MAXIMAL DE ZONE FRONTIERES

  integer         nflagm
  parameter      (nflagm = 100)

  !=============================================================================
  ! 3. Conditions aux limites

  !     NDLAGM : NOMBRE MAXIMAL DE DONNEES SUR LES PARTICULES (REELS)

  integer         ndlagm
  parameter      (ndlagm = 50)

  !     NDLAIM : NOMBRE MAXIMAL DE DONNEES SUR LES PARTICULES (ENTIERS)

  integer         ndlaim
  parameter      (ndlaim = 10)

  !=============================================================================
  ! 4. Schema en temps

  !     NVGAUS : NOMBRE DE VARIABLES ALEATOIRES GAUSSIENNES PAR PARTICULES

  integer         nvgaus
  parameter      (nvgaus = 9)


  !=============================================================================
  ! 5. Mouvement Brownien

  !     NVGAUS : NOMBRE DE VARIABLES ALEATOIRES GAUSSIENNES PAR PARTICULES

  integer         nbrgau
  parameter      (nbrgau = 6)


  !=============================================================================
  ! 6. Variables utilisateurs supplementaires

  !     NUSVAR : Nombre maximum de variables utilisateur supplementaires

  integer         nusvar
  parameter      (nusvar = 10)

!     NUSSTA : Nombre maximum de stats utilisateur supplementaires

  integer         nussta
  parameter      (nussta = 20)

  !     NUSBRD : Nombre maximum interactions particules/frontieres
  !              utilisateur supplementaires

  integer         nusbrd
  parameter      (nusbrd = 10)

  !============================================================================
  ! 7. Affichages et fichiers suite

  !     NVPLMX : Nombre maximum de variables

  integer         nvplmx
  parameter      (nvplmx = 50)

  !=============================================================================
  ! 8. Types d'interaction au bord

  integer         ientrl     , isortl     , irebol
  integer         idepo1     , idepo2
  integer         iencrl     , jbord1     , jbord2
  integer         jbord3     , jbord4     , jbord5
  integer         idepfa     , isymtl

  parameter      (ientrl =  1, isortl =  2, irebol =  3)
  parameter      (idepo1 =  4, idepo2 =  5)
  parameter      (iencrl =  7, jbord1 =  8, jbord2 =  9)
  parameter      (jbord3 = 10, jbord4 = 11, jbord5 = 12)
  parameter      (idepfa = 13, isymtl = 14)

  !=============================================================================

end module lagpar
