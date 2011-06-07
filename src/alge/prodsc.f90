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

subroutine prodsc &
!================

 ( ncelet , ncel   , isqrt  , va     , vb     , vavb   )

!===============================================================================
! FONCTION :
! ----------
!                                    ______
! PRODUIT SCALAIRE VAPVB = VA.VB OU \/ VA.VB  SI ISQRT=1

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! isqrt            ! e  ! <-- ! indicateur = 1 pour prendre la racine          !
! va, vb(ncelet    ! tr ! <-- ! vecteurs a multiplier                          !
! vavb             ! r  ! --> ! produit scalaire                               !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use parall

!===============================================================================

implicit none

! Arguments

integer          ncelet,ncel,isqrt
double precision vavb
double precision va(ncelet),vb(ncelet)

! Local variables

integer incx, incy
double precision ddot
external         ddot

!===============================================================================

incx = 1
incy = 1
vavb = ddot(ncel, va, incx, vb, incy)

if (irangp.ge.0) call parsom (vavb)
                 !==========
if( isqrt.eq.1 ) vavb= sqrt(vavb)
!----
! FIN
!----

return

end subroutine

