/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*-----------------------------------------------------------------------------*/

/*!
  \page us_hist Non-standard monitoring point definition (ushist.f90)

  \section ushist_h_intro Introduction

  The \ref ushist function allows one to define non-standard monitoring points.

  The following example shows the definition of 4 monitoring points and
  for each variable, a file is written.

  \section ushist_h_ushist_arg Arguments of ushist

  The \ref ushist subroutine has the following arguments :

  \snippet ushist.f90 arg

  \section ushist_h_variable Local variable declaration

  Local variables are declared hereafter :

  \snippet ushist.f90 loc_var_dec

  \section ushist_h_initialization Initialization

  The \c ipass variable, which is the current number of \ref ushist calls, is initialized.

  \snippet ushist.f90 init

  \section ushist_h_search Searching of the monitoring points

  The routine \ref findpt is used to find the number of the closest cell center
  to a (x,y,z) point.

  \snippet ushist.f90 search

  \section ushist_h_open Opening files

  For each variable, a file is opened.

  \snippet ushist.f90 open

  \section ushist_h_write Writing files

  For each variable, the file opened in \ref shist_h_open is then written. It
  contains the time step value, the physical time value and the variable value
  at each monitoring points.

  \snippet ushist.f90 write

  \section ushist_h_close Closing files

  Each file is closed.

  \snippet ushist.f90 close

*/
