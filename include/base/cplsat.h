!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

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

!                              cplsat.h
!===============================================================================

!  COUPLAGE CODE / CODE - GESTION DES PARAMETRES PRINCIPAUX

! NBRCPL : NOMBRE DE COUPLAGE CODE_SATURNE / CODE_SATURNE
! IFACCP : INDICATEUR DE COUPLAGE FACE/FACE UNIQUEMENT
! IMOBIL : INDICATEUR DE MAILLAGE MOBILE POUR LES TURBOMACHINES

integer           nbrcpl, ifaccp, imobil

common / icplcs / nbrcpl, ifaccp, imobil

! NBCPMX : NOMBRE DE COUPLAGE MAX ADMISSIBLE

integer   nbcpmx
parameter(nbcpmx=10)

! ITURCP(NBCPMX,NPHSMX) : MODELE DE TURBULENCE DE L'INSTANCE DISTANTE
! IMAJCP(NBCPMX)        : INDICE DE MISE A JOUR DE LA LOCALISATION DU COUPLAGE
! NVARCP(NBCPMX)        : NOMBRE DE VARIABLE A ENVOYER/RECEVOIR
! NVARTO(NBCPMX)        : TAILLE DES TABLEAUX D'ECHANGE


integer           iturcp(nbcpmx,nphsmx), imajcp(nbcpmx)
integer           nvarcp(nbcpmx), nvarto(nbcpmx)
common / icplcs / iturcp, imajcp, nvarcp, nvarto

