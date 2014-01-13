!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine clpv2f &
!================

 ( ncelet , ncel   , nvar   ,                                     &
   iwaphi ,                                                       &
   rtp    )

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
! iwaphi           ! e  ! <-- ! niveau d'impression                            !
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
use optcal
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar, ncelet, ncel
integer          iwaphi
double precision rtp(ncelet,nvar)

! Local variables

integer          iel, ipp
integer          nclpmx, nclpmn
double precision xphi, xal, vmin(1), vmax(1), var

!===============================================================================

!===============================================================================
!  1. Pour le phi-fbar et BL-v2/k model, reperage des valeurs de phi
!     superieures a 2 et clipping de phi en valeur absolue
!     pour les valeurs negatives
!===============================================================================

!===============================================================================
!     1.a Stockage Min et Max pour listing
!===============================================================================

ipp = ipprtp(iphi)

vmin(1) =  grand
vmax(1) = -grand
do iel = 1, ncel
  var = rtp(iel,iphi)
  vmin(1) = min(vmin(1),var)
  vmax(1) = max(vmax(1),var)
enddo

!==============================================================================
!     1.b Reperage des valeurs superieures a 2, pour affichage seulement
!==============================================================================

if (iwaphi.ge.2) then
  nclpmx = 0
  do iel = 1, ncel
    if (rtp(iel,iphi).gt.2.d0) nclpmx = nclpmx+1
  enddo
  if(irangp.ge.0) call parcpt(nclpmx)
                  !==========
  if (nclpmx.gt.0) write(nfecra,1000) nclpmx
endif

!==============================================================================
!     1.c Clipping en valeur absolue pour les valeurs negatives
!==============================================================================

nclpmn = 0
do iel = 1, ncel
  xphi = rtp(iel,iphi)
  if (xphi.lt.0.d0) then
    rtp(iel,iphi) = -xphi
    nclpmn = nclpmn + 1
  endif
enddo

call log_iteration_clipping_field(ivarfl(iphi), nclpmn, 0, vmin, vmax)

!===============================================================================
!  2. Pour le BL-v2/k model, clipping de alpha a 0 pour les valeurs negatives
!     et a 1 pour les valeurs superieurs a 1
!===============================================================================

if (iturb.eq.51) then

!===============================================================================
!     2.a Stockage Min et Max pour listing
!===============================================================================

  ipp = ipprtp(ial)

  vmin(1) =  grand
  vmax(1) = -grand
  do iel = 1, ncel
    var = rtp(iel,ial)
    vmin(1) = min(vmin(1),var)
    vmax(1) = max(vmax(1),var)
  enddo

!==============================================================================
!     2.b Clipping a 0 pour les valeurs negatives et a 1 pour les valeurs
!         superieures a 1
!==============================================================================

  nclpmn = 0
  nclpmx = 0
  do iel = 1, ncel
    xal = rtp(iel,ial)
    if (xal.lt.0.d0) then
      rtp(iel,ial) = 0.d0
      nclpmn = nclpmn + 1
    endif
    if (xal.gt.1.d0) then
      rtp(iel,ial) = 1.d0
      nclpmx = nclpmx + 1
    endif
  enddo

  call log_iteration_clipping_field(ivarfl(ial), nclpmn, nclpmx, vmin, vmax)

endif


!===============================================================================
! ---> Formats
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format('ATTENTION VARIABLE PHI',                            &
     'VALEUR MAXIMALE PHYSIQUE DE 2 DEPASSEE SUR ',I10,           &
     ' CELLULES')

#else

 1000 format('WARNING VARIABLE PHI',                              &
     'MAXIMUM PHYSICAL VALUE OF 2 EXCEEDED FOR ',I10,             &
     ' CELLS')

#endif

return

end subroutine
