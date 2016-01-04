!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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
 ( ncel   ,                                                       &
   iwaphi )

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
! ncel             ! e  ! <-- ! nombre de cellules                             !
! iwaphi           ! e  ! <-- ! niveau d'impression                            !
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
use field

!===============================================================================

implicit none

! Arguments

integer          ncel
integer          iwaphi

! Local variables

integer          iel
integer          nclpmx(1), nclpmn(1)
double precision xphi, xal, vmin(1), vmax(1), var

double precision, dimension(:), pointer :: cvar_al, cvar_phi

!===============================================================================

call field_get_val_s(ivarfl(iphi), cvar_phi)
if (iturb.eq.51) call field_get_val_s(ivarfl(ial), cvar_al)

!===============================================================================
!  1. Pour le phi-fbar et BL-v2/k model, reperage des valeurs de phi
!     superieures a 2 et clipping de phi en valeur absolue
!     pour les valeurs negatives
!===============================================================================

!===============================================================================
!     1.a Stockage Min et Max pour listing
!===============================================================================

vmin(1) =  grand
vmax(1) = -grand
do iel = 1, ncel
  var = cvar_phi(iel)
  vmin(1) = min(vmin(1),var)
  vmax(1) = max(vmax(1),var)
enddo

!==============================================================================
!     1.b Reperage des valeurs superieures a 2, pour affichage seulement
!==============================================================================

if (iwaphi.ge.2) then
  nclpmx(1) = 0
  do iel = 1, ncel
    if (cvar_phi(iel).gt.2.d0) nclpmx(1) = nclpmx(1)+1
  enddo
  if (irangp.ge.0) call parcpt(nclpmx(1))
  if (nclpmx(1).gt.0) write(nfecra,1000) nclpmx(1)
endif

!==============================================================================
!     1.c Clipping en valeur absolue pour les valeurs negatives
!==============================================================================

nclpmn(1) = 0
do iel = 1, ncel
  xphi = cvar_phi(iel)
  if (xphi.lt.0.d0) then
    cvar_phi(iel) = -xphi
    nclpmn(1) = nclpmn(1) + 1
  endif
enddo

call log_iteration_clipping_field(ivarfl(iphi), nclpmn(1), 0, vmin, vmax,nclpmn(1), nclpmx(1))

!===============================================================================
!  2. Pour le BL-v2/k model, clipping de alpha a 0 pour les valeurs negatives
!     et a 1 pour les valeurs superieurs a 1
!===============================================================================

if (iturb.eq.51) then

!===============================================================================
!     2.a Stockage Min et Max pour listing
!===============================================================================

  vmin(1) =  grand
  vmax(1) = -grand
  do iel = 1, ncel
    var = cvar_al(iel)
    vmin(1) = min(vmin(1),var)
    vmax(1) = max(vmax(1),var)
  enddo

!==============================================================================
!     2.b Clipping a 0 pour les valeurs negatives et a 1 pour les valeurs
!         superieures a 1
!==============================================================================

  nclpmn(1) = 0
  nclpmx(1) = 0
  do iel = 1, ncel
    xal = cvar_al(iel)
    if (xal.lt.0.d0) then
      cvar_al(iel) = 0.d0
      nclpmn(1) = nclpmn(1) + 1
    endif
    if (xal.gt.1.d0) then
      cvar_al(iel) = 1.d0
      nclpmx(1) = nclpmx(1) + 1
    endif
  enddo

  call log_iteration_clipping_field(ivarfl(ial), nclpmn(1), nclpmx(1), vmin,vmax,nclpmn(1), nclpmx(1))

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
