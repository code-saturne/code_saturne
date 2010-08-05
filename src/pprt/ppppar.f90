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

!                              ppppar.h

!===============================================================================

!            INCLUDE GENERAL PROPRE A LA PHYSIQUE PARTICULIERE
!                    CONTENANT DES PARAMETRES COMMUNS
!                        (A PLUSIEURS INCLUDES)
!-------------------------------------------------------------------------------

! --> NB DE ZONES DE BORD MAXIMAL
integer    nbzppm
parameter (nbzppm=2000)
! --> NUMERO DE ZONE DE BORD MAXIMAL
integer    nozppm
parameter (nozppm=2000)


!--> POINTEURS VARIABLES COMBUSTION CHARBON PULVERISE cpincl, ppincl

!       NCHARM        --> Nombre maximal de charbons
!       NCPCMX        --> Nombre maximal de classes par charbon
!       NCLCPM        --> Nombre total de classes

integer    ncharm  , ncpcmx   , nclcpm
parameter (ncharm=3, ncpcmx=10, nclcpm=ncharm*ncpcmx)
! -->
! FIN
