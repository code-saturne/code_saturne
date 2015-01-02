!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

!> \file lagpar.f90
!> \brief Module for Lagrangian model (parameters)

module lagpar

  !=============================================================================

  implicit none

  !=============================================================================

  !> \defgroup lagpar Module for Lagrangian (parameters)

  !> \addtogroup lagpar
  !> \{

  !         Trois modules complementaires
  !                            lagran qui porte les non dimensions
  !                            lagdim qui porte les dimensions variables
  !                            lagpar qui porte les parametres

  !=============================================================================

  !> \defgroup classes_particles Classes and particles

  !> \addtogroup classes_particles
  !> \{

  !> maximal number of coal classes
  integer         ncharm2
  parameter      (ncharm2 = 5)

  !> maximal number of volumetric statistics
  integer         nclstm
  parameter      (nclstm = 100)

  !> maximal number of layer per coal particle
  integer         nlayer
  parameter      (nlayer = 5)

  !> \}

  !=============================================================================

  !> \defgroup lag_bcs Boundary conditions

  !> \addtogroup lag_bcs
  !> \{

  !> maximal number of boundary zones
  integer         nflagm
  parameter      (nflagm = 100)

  !=============================================================================

  !> maximal number of particle real data
  integer         ndlagm
  parameter      (ndlagm = 50+4*nlayer)

  !> maximal number of particle integer data
  integer         ndlaim
  parameter      (ndlaim = 10)

  !> \}

  !=============================================================================

  !> \defgroup lag_brownian_motion Brownian motion

  !> \addtogroup lag_brownian_motion
  !> \{

  !> number of gaussian random variables by particle
  integer         nvgaus
  parameter      (nvgaus = 9)

  ! TODO
  integer         nbrgau
  parameter      (nbrgau = 6)

  !> \}

  !=============================================================================

  !> \defgroup lag_user_variable Additional user variables

  !> \addtogroup lag_user_variable
  !> \{

  !> maximal number of additional user variables
  integer         nusvar
  parameter      (nusvar = 10)

  !> maximal number of additional user volume statistics
  integer         nussta
  parameter      (nussta = 20)

  !> maximal number of additional user particles/boundary interactions
  integer         nusbrd
  parameter      (nusbrd = 10)

  !> \}

  !============================================================================

  !> \defgroup lag_printing Printing

  !> \addtogroup lag_printing
  !> \{

  !> maximal number of variables
  integer         nvplmx
  parameter      (nvplmx = 50+4*nlayer)

  !> \}

  !=============================================================================
  ! 8. Types d'interaction au bord

  !> \defgroup lag_typ_bnd_interaction Boundary interation type
  !> (value of \c iusclb(izone))

  !> \addtogroup lag_typ_bnd_interaction
  !> \{


  !> particle injection zone. For each particle class associated with this zone,
  !> information must be provided. If a particle trajectory may cross an
  !> injection zone, then this particle leaves the calculation domain.
  integer         ientrl
  !> constant = 2 !TODO
  integer         isortl
  !> constant = 3
  integer         irebol
  !> constant = 4
  integer         idepo1
  !> constant = 5
  integer         idepo2
  !> constant = 6
  integer         iencrl
  !> constant = 7
  integer         jbord1
  !> constant = 8
  integer         jbord2
  !> constant = 9
  integer         jbord3
  !> constant = 10
  integer         jbord4
  !> constant = 11
  integer         jbord5
  !> constant = 12
  integer         idepfa
  !> constant = 13
  integer         isymtl
  !> constant = 14

  parameter      (ientrl =  1, isortl =  2, irebol =  3)
  parameter      (idepo1 =  4, idepo2 =  5)
  parameter      (iencrl =  7, jbord1 =  8, jbord2 =  9)
  parameter      (jbord3 = 10, jbord4 = 11, jbord5 = 12)
  parameter      (idepfa = 13, isymtl = 14)

  !> \}
  !> \}
  !=============================================================================

end module lagpar
