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

subroutine clpv2f &
!================

 ( ncelet , ncel   , nvar   , nphas  ,                            &
   iphas  , iwaphi ,                                              &
   propce , rtp    )

!===============================================================================
! FONCTION :
! ----------

! CLIPPING DE PHI EN V2F (PAS DE CLIPPING SUR F_BARRE)

!-------------------------------------------------------------------------------
! Arguments
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! e  ! <-- ! nombre de cellules                             !
! nvar             ! e  ! <-- ! nombre de variables                            !
! nphas            ! i  ! <-- ! number of phases                               !
! iphas            ! i  ! <-- ! phase number                                   !
! iwaphi           ! e  ! <-- ! niveau d'impression                            !
! propce           ! tr ! <-- ! tableaux des variables au pdt courant          !
!(ncelet,*         !    !     !                                                !
! rtp              ! tr ! <-- ! tableaux des variables au pdt courant          !
! (ncelet,nvar)    !    !     !                                                !
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
use entsor
use numvar
use cstnum
use parall

!===============================================================================

implicit none

! Arguments

integer          nvar, ncelet, ncel, nphas
integer          iphas, iwaphi
double precision propce(ncelet,*)
double precision rtp(ncelet,nvar)

! Local variables

integer          iel, ipp
integer          iphiph
integer          nclpmx, nclpmn
double precision xphi, vmin, vmax, var

!===============================================================================


iphiph = iphi(iphas)

!===============================================================================
!  ---> Stockage Min et Max pour listing
!===============================================================================

ipp = ipprtp(iphiph)

vmin =  grand
vmax = -grand
do iel = 1, ncel
  var = rtp(iel,iphiph)
  vmin = min(vmin,var)
  vmax = max(vmax,var)
enddo
if (irangp.ge.0) then
  call parmin(vmin)
  !==========
  call parmax(vmax)
  !==========
endif
varmna(ipp) = vmin
varmxa(ipp) = vmax

!==============================================================================
!  ---> Reperage des valeurs superieures a 2, pour affichage seulement
!==============================================================================

if (iwaphi.ge.2) then
  nclpmx = 0
  do iel = 1, ncel
    if (rtp(iel,iphiph).gt.2.d0) nclpmx = nclpmx+1
  enddo
  if(irangp.ge.0) call parcpt(nclpmx)
                  !==========
  if (nclpmx.gt.0) write(nfecra,1000) nclpmx
endif

!==============================================================================
!  ---> Clipping en valeur absolue pour les valeurs negatives
!==============================================================================

nclpmn = 0
do iel = 1, ncel
  xphi = rtp(iel,iphiph)
  if (xphi.lt.0.d0) then
    rtp(iel,iphiph) = -xphi
    nclpmn = nclpmn + 1
  endif
enddo
if(irangp.ge.0) call parcpt(nclpmn)
                !==========
iclpmn(ipp) = nclpmn

!===============================================================================
! ---> Formats
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format('ATTENTION VARIABLE PHI'                             &
     'VALEUR MAXIMALE PHYSIQUE DE 2 DEPASSEE SUR ',I10,           &
     ' CELLULES')

#else

 1000 format('WARNING VARIABLE PHI'                               &
     'MAXIMUM PHYSICAL VALUE OF 2 EXCEEDED FOR ',I10,             &
     ' CELLS')

#endif

return

end subroutine
