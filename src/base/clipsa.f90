!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2011 EDF S.A., France

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

subroutine clipsa &
!================

 ( ncelet , ncel   , nvar   , nphas  ,                            &
   iphas  , iclip  , iwarnu ,                                     &
   propce , rtp    )

!===============================================================================
! Purpose:
! --------

! Clipping of nusa for the Spalart-Allmaras model

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! e  ! <-- ! nombre de variables                            !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nphas            ! i  ! <-- ! number of phases                               !
! iphas            ! i  ! <-- ! phase number                                   !
! iclip            ! e  ! <-- ! indicateur = 0 on utilise viscl0               !
!                  !    !     !            sinon on utilise viscl              !
! iwarnu           ! e  ! <-- ! niveau d'impression                            !
! propce           ! tr ! <-- ! tableaux des variables au pdt courant          !
!(ncelet,*         !    !     !                                                !
! rtp              ! tr ! <-- ! tableaux des variables au pdt courant          !
! (ncelet     )    !    !     !                                                !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use cstphy
use cstnum
use entsor
use optcal
use parall

!===============================================================================

implicit none

! Arguments

integer          nvar, ncelet, ncel, nphas
integer          iphas, iclip, iwarnu
double precision propce(ncelet,*)
double precision rtp(ncelet,nvar)

! Local variables

integer          iclpnu,iel
integer          ivar,ipp,inuiph
double precision xnu, vmin, vmax, var
double precision epz2

!===============================================================================

inuiph = inusa(iphas)


! Une petite valeur pour eviter des valeurs exactement nulles.

epz2 = epzero**2

!===============================================================================
! ---> Stockage Min et Max pour listing
!===============================================================================

ipp  = ipprtp(inuiph)

vmin =  grand
vmax = -grand
do iel = 1, ncel
  var = rtp(iel,inuiph)
  vmin = min(vmin,var)
  vmax = max(vmax,var)
enddo
if (irangp.ge.0) then
  call parmax (vmax)
!==========
  call parmin (vmin)
!==========
endif
varmna(ipp) = vmin
varmxa(ipp) = vmax

!===============================================================================
! ---> Clipping "standard" NUSA>0
!===============================================================================


iclpnu = 0
do iel = 1, ncel
  xnu = rtp(iel,inuiph)
  if (xnu.lt.0.D0) then
    iclpnu = iclpnu + 1
    rtp(iel,inuiph) = 0.d0
  endif
enddo

if (irangp.ge.0) then
  call parcpt (iclpnu)
!==========
endif

! ---  Stockage nb de clippings pour listing

iclpmn(ipprtp(inuiph)) = iclpnu


!===============================================================================
! ---> Formats
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS clipsa                           ',/,&
'@    =========                                               ',/,&
'@     APPEL DE clipsa              AVEC OPTION = ',I10        ,/,&
'@     Phase : ',I10                                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut pas etre execute.                       ',/,&
'@                                                            ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1010 format(                                                           &
 I10,' VALEURS DE NUSA INCORRECTE                            ')

#else

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN clipsa                                ',/,&
'@    ========                                                ',/,&
'@     CALL OF clipsa               WITH OPTION = ',I10        ,/,&
'@     Phase : ',I10                                           ,/,&
'@                                                            ',/,&
'@  The calulation will not be run.                           ',/,&
'@                                                            ',/,&
'@  Contact the support.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1010 format(                                                           &
 I10,' SA    VALUES INCORRECT                       ')

#endif

return

end subroutine
