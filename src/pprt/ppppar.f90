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

! General module for specific physics containing common parameters

module ppppar

  !===========================================================================

  ! --> Nb de zones de bord maximal
  integer    nbzppm
  parameter (nbzppm=2000)

  ! --> Numero de zone de bord maximal
  integer    nozppm
  parameter (nozppm=2000)

  !--> Pointeurs variables combustion charbon pulverise cpincl, ppincl

  !    ncharm --> nombre maximal de charbons
  !    ncpcmx --> nombre maximal de classes par charbon
  !    nclcpm --> Nombre total de classes

  integer    ncharm  , ncpcmx   , nclcpm
  parameter (ncharm=3, ncpcmx=10, nclcpm=ncharm*ncpcmx)

  !=============================================================================

end module ppppar
