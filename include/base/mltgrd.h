!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

!                              mltgrd.h
!===============================================================================

! MULTIGRILLE
! -----------
!   NCEGRM : NOMBRE MAX DE CELLULES SUR MAILLAGE LE PLUS GROSSIER
!   NGRMAX : NOMBRE MAX DE NIVEAUX DE MAILLAGES
!   NAGMX0 : PARAMETRE CONSTRUCTION DE  MAILLAGE AUTOMATIQUE
!   IAGMX0 : PARAMETRE CONSTRUCTION DE  MAILLAGE AUTOMATIQUE
!   NCPMGR : Si > 0, active le post traitement de l'agglomeration, en
!            projetant les numeros de cellules grossieres sur le
!            maillage fin (modulo NCPMGR(IVAR))

integer           ncegrm, ngrmax,                                 &
                  nagmx0(nvarmx), iagmx0(nvarmx), ncpmgr(nvarmx)
common / ioptmg / ncegrm, ngrmax,                                 &
                  nagmx0, iagmx0, ncpmgr

! FIN
