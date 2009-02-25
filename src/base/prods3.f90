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

                  subroutine prods3                               &
!================

 ( ncelet , ncel   , isqrt  ,                                     &
   va1    , vb1    , va2    , vb2    , va3    , vb3    ,          &
   vavb1  , vavb2  , vavb3  )

!===============================================================================
! FONCTION :
! ----------

! CALCUL SIMULTANE DE 3 PRODUITS SCALAIRES
!                   ______
! VAPVB = VA.VB OU \/ VA.VB  SI ISQRT=1

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! isqrt            ! e  ! <-- ! indicateur = 1 pour prendre la racine          !
! va1(), vb1()     ! tr ! <-- ! premiers   vecteurs a multiplier               !
! va2(), vb2()     ! tr ! <-- ! seconds    vecteurs a multiplier               !
! va3(), vb3()     ! tr ! <-- ! troisiemes vecteurs a multiplier               !
! vavb1            ! r  ! --> ! premier   produit scalaire                     !
! vavb2            ! r  ! --> ! second    produit scalaire                     !
! vavb3            ! r  ! --> ! troisieme produit scalaire                     !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "parall.h"

!===============================================================================

! Arguments

integer          ncelet,ncel,isqrt
double precision vavb1, vavb2, vavb3
double precision va1(ncelet),vb1(ncelet)
double precision va2(ncelet),vb2(ncelet)
double precision va3(ncelet),vb3(ncelet)

! VARIABLES LOCALES

integer nvavb
double precision vavb(3)

integer incx, incy
double precision ddot
external         ddot

!===============================================================================

incx = 1
incy = 1
vavb(1) = ddot(ncel, va1, incx, vb1, incy)
vavb(2) = ddot(ncel, va2, incx, vb2, incy)
vavb(3) = ddot(ncel, va3, incx, vb3, incy)

if (irangp.ge.0) then
  nvavb = 3
  call parrsm (nvavb, vavb)
  !==========
endif

vavb1 = vavb(1)
vavb2 = vavb(2)
vavb3 = vavb(3)

if (isqrt.eq.1) then
  vavb1 = sqrt(vavb1)
  vavb2 = sqrt(vavb2)
  vavb3 = sqrt(vavb3)
endif

!----
! FIN
!----

return

end

